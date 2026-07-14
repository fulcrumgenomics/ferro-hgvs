//! In-memory, JSON-backed reference provider.
//!
//! [`JsonProvider`] loads transcripts (and optional protein and genomic
//! sequences) from a `transcripts.json` document — the output of
//! `ferro convert-gff` / `ferro build-transcript`, or a hand-authored file — and
//! serves them for normalization and projection. It also backs the Python
//! `Normalizer(reference_json=...)` path, so despite living in the historically
//! `mock`-named module it is a real, production-facing provider (transcript-level
//! normalization is genuine and correct; add a `genomic_sequences` map to make it
//! genome-capable — see `docs/transcripts_json_schema.md`). Its
//! [`JsonProvider::with_test_data`] constructor is the only genuinely test-only
//! part. The former name `MockProvider` remains as a compatibility alias.

use crate::error::FerroError;
use crate::hgvs::variant::Accession;
use crate::reference::provider::{GenomicPlacement, ReferenceProvider};
use crate::reference::transcript::{Strand, Transcript};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::{Arc, OnceLock};

/// In-memory reference provider that loads transcripts, proteins, and optional
/// genomic sequences from a `transcripts.json` document. Backs
/// `Normalizer(reference_json=...)`. See the module documentation.
#[derive(Clone)]
pub struct JsonProvider {
    transcripts: HashMap<String, Arc<Transcript>>,
    proteins: HashMap<String, String>,
    /// Genomic sequences keyed by contig name
    genomic_sequences: HashMap<String, String>,
    /// Transcript ids this provider should report as *not* version-exact from
    /// [`ReferenceProvider::has_transcript_version_exact`]. Lets tests exercise
    /// the protein path's cross-version gate (#505) without a full FASTA/cdot
    /// manifest. Empty by default, so an unconfigured mock reports every
    /// transcript as version-exact (matching the trait default).
    non_version_exact: HashSet<String>,
    /// Chromosomal placements of genomic parent references (`NG_`/`LRG_`),
    /// keyed by the parent's full versioned accession (e.g. `"NG_007485.1"`).
    /// Lets tests exercise the #480 NC→parent-frame re-anchoring without a real
    /// RefSeqGene/LRG placement source.
    /// One entry per genome build (the build is encoded in each placement's
    /// `nc` accession), so tests can exercise per-build selection (#653).
    genomic_placements: HashMap<String, Vec<GenomicPlacement>>,
    /// Legacy gene-model selector resolution: gene Symbol (upper-case) →
    /// reference-standard transcript accession (e.g. `"TIMM8B"` →
    /// `"NM_012459.4"`). Lets tests exercise the #500/#637 selector rewrite
    /// without ingesting NCBI's `LRG_RefSeqGene` table.
    legacy_gene_models: HashMap<String, String>,
    /// Silent version-substitution map: requested transcript accession (e.g.
    /// `"NM_002001.2"`) → the *different*-version accession the provider should
    /// actually serve for it (e.g. `"NM_002001.4"`). Lets tests model the #785
    /// silent version-strip substitution that a real `MultiFastaProvider`
    /// performs via `resolve_name`'s base→latest fallback, which `JsonProvider`
    /// otherwise refuses by design. Empty by default, so an unconfigured mock
    /// never substitutes a version.
    version_substitutions: HashMap<String, String>,
    /// Build-keyed transcripts: `(transcript_id, build)` → the `Transcript` whose
    /// bases/CDS the provider should serve for an input addressed on that genome
    /// build. Lets a test model a transcript whose *sequence differs between
    /// GRCh37 and GRCh38* — the real divergence a genome-reconstructed transcript
    /// exhibits, which the single build-agnostic `transcripts` map cannot capture
    /// (#843). Consulted by [`Self::get_transcript_for_accession`] when the input
    /// accession carries a build-bearing `genomic_context` (an `NC_*.10`/`.11`
    /// parent); falls back to the build-agnostic `transcripts` map otherwise.
    /// Empty by default, so an unconfigured mock is build-agnostic as before.
    build_transcripts: HashMap<(String, &'static str), Arc<Transcript>>,
}

/// Backwards-compatible alias for the former name of [`JsonProvider`].
///
/// The provider was renamed from `MockProvider` because it is not a test mock — it
/// backs the production `Normalizer(reference_json=...)` path. This alias keeps
/// existing code compiling; prefer [`JsonProvider`] in new code. (Not marked
/// `#[deprecated]` yet, to avoid warning-storming the many existing call sites; a
/// follow-up can migrate them and then deprecate.)
pub type MockProvider = JsonProvider;

/// Validate that a JSON reference declaring genomic capability actually backs
/// every placed transcript with genomic sequence bytes that are the *right length*
/// **and** the *right bases*.
///
/// Only enforced when `genomic_sequences` is non-empty — i.e. the reference claims
/// `has_genomic_data() == true`, which arms the genome-aware normalization rules.
/// Note this is a **per-file** requirement: once any `genomic_sequences` are
/// present, every placed transcript must be backed, so a partially-genomic
/// reference is rejected rather than silently granting genome capability to only
/// some transcripts. Two checks per placed transcript:
///
/// 1. **Structure** — the transcript's contig must be present in `genomic_sequences`
///    and long enough to contain its maximal genomic coordinate. Catches a missing
///    or too-short contig (#1012 comment 1).
/// 2. **Content** — when every exon carries genomic coordinates and the transcript
///    has a `sequence`, the transcript sequence is *reconstructed* from the genomic
///    bytes at the exon placement (concatenating the forward exon slices, then
///    reverse-complementing on the minus strand — exactly how
///    `convert-gff`/`build-transcript` build it) and must equal the stated
///    `sequence`. This catches genomic bytes that are the right length but the
///    **wrong bases** — a padded junk blob, a different assembly, or shifted
///    coordinates — which would otherwise silently corrupt genome-aware
///    normalization of del/dup/inv/intronic edits (whose recommended forms carry no
///    stated reference bases, so nothing downstream flags the wrong genome).
///
/// Two residual cases are the emitter's responsibility and cannot be caught here:
/// bytes in the *intronic* gaps between exons (no independent copy to check
/// against), and a coordinate frame that is *self-consistently* wrong (the
/// transcript sequence and genomic bytes agree because both were sliced from the
/// same mismatched FASTA). `convert-gff --emit-genomic-sequences` avoids both by
/// construction.
fn validate_genomic_placement_backing(
    transcripts: &[Transcript],
    genomic_sequences: &HashMap<String, String>,
) -> Result<(), FerroError> {
    if genomic_sequences.is_empty() {
        return Ok(());
    }

    for tx in transcripts {
        let Some(chromosome) = tx.chromosome.as_deref() else {
            continue;
        };
        // Highest 1-based genomic coordinate this transcript's placement references.
        let required_len = tx
            .exons
            .iter()
            .filter_map(|exon| exon.genomic_end)
            .chain(tx.genomic_end)
            .max();
        let Some(required_len) = required_len else {
            continue;
        };

        // (1) Structure.
        let contig = match genomic_sequences.get(chromosome) {
            None => {
                return Err(FerroError::Json {
                    msg: format!(
                        "reference declares genomic sequences but transcript '{}' is placed on \
                         contig '{}', which has no entry in `genomic_sequences` (once any \
                         genomic_sequences are present, every placed transcript must be backed)",
                        tx.id, chromosome
                    ),
                });
            }
            Some(sequence) if (sequence.len() as u64) < required_len => {
                return Err(FerroError::Json {
                    msg: format!(
                        "`genomic_sequences` for contig '{}' has length {}, too short for the \
                         placement of transcript '{}' (needs at least {})",
                        chromosome,
                        sequence.len(),
                        tx.id,
                        required_len
                    ),
                });
            }
            Some(sequence) => sequence,
        };

        // (2) Content.
        validate_reconstructed_transcript_sequence(tx, chromosome, contig)?;
    }

    Ok(())
}

/// Reconstruct a transcript's sequence from the genomic contig bytes at its exon
/// placement and confirm it equals the stated `sequence`. See
/// [`validate_genomic_placement_backing`] for the rationale. Skips (returns `Ok`)
/// when the transcript has no `sequence`, or any exon lacks genomic coordinates, or
/// an exon's coordinates are malformed or out of the contig's range — there is
/// nothing to reconstruct against in those cases, and the structural check has
/// already bounded the maximal coordinate.
fn validate_reconstructed_transcript_sequence(
    tx: &Transcript,
    chromosome: &str,
    contig: &str,
) -> Result<(), FerroError> {
    let Some(tx_sequence) = tx.sequence.as_deref() else {
        return Ok(());
    };

    let mut reconstructed = String::with_capacity(tx_sequence.len());
    for exon in &tx.exons {
        let (Some(g_start), Some(g_end)) = (exon.genomic_start, exon.genomic_end) else {
            return Ok(()); // partial placement — cannot fully reconstruct.
        };
        if g_start < 1 || g_end < g_start || g_end as usize > contig.len() {
            return Ok(()); // malformed / out-of-range coords — not this check's job.
        }
        reconstructed.push_str(&contig[(g_start - 1) as usize..g_end as usize]);
    }

    if tx.strand == Strand::Minus {
        reconstructed = crate::sequence::reverse_complement(&reconstructed);
    }

    if !reconstructed.eq_ignore_ascii_case(tx_sequence) {
        return Err(FerroError::Json {
            msg: format!(
                "`genomic_sequences` for contig '{}' do not reconstruct transcript '{}'s \
                 sequence from its exon placement: the genomic bytes are not the reference the \
                 transcript sequence was built from (wrong assembly, shifted coordinates, or \
                 placeholder/junk bytes)",
                chromosome, tx.id
            ),
        });
    }

    Ok(())
}

/// Highest `transcripts.json` container **major** schema version this build can
/// interpret. The container `version` is `"MAJOR.MINOR"` (e.g. `"1.0"`). A newer
/// **major** is rejected because it may reinterpret existing fields. A newer
/// *minor* is only accepted if it adds no fields: the object form uses
/// `deny_unknown_fields`, so a minor that introduces a new key is still rejected
/// (with an "unknown field" error) by an older build — minor-forward-compatibility
/// is therefore limited to value-only changes, not new fields.
const SUPPORTED_TRANSCRIPTS_JSON_MAJOR: u32 = 1;

/// Validate the container `version` from a `transcripts.json` object.
///
/// Accepts a JSON string (`"1.0"`) or number (`1`) — both are natural hand-written
/// forms. `None` (a bare-array reference, or an object with no `version`) is
/// accepted for backward compatibility. A present version must be `MAJOR[.MINOR…]`
/// with every dot-separated component numeric, and its major must not exceed
/// [`SUPPORTED_TRANSCRIPTS_JSON_MAJOR`]; anything else is a clear load-time error.
///
/// This is validated against the **raw** JSON value before the `deny_unknown_fields`
/// deserialize, so a newer-major file gets this friendly "upgrade" message even when
/// it also carries fields this build does not recognize.
fn validate_transcripts_json_version(
    version: Option<&serde_json::Value>,
) -> Result<(), FerroError> {
    let Some(version) = version else {
        return Ok(());
    };

    let unrecognized = || FerroError::Json {
        msg: format!(
            "unrecognized transcripts.json schema version {} (expected e.g. \"1.0\")",
            version
        ),
    };

    // Accept a string ("1.0") or a bare number (1 / 1.0).
    let version_str = match version {
        serde_json::Value::String(s) => s.trim().to_string(),
        serde_json::Value::Number(n) => n.to_string(),
        _ => return Err(unrecognized()),
    };

    // Every dot-separated component must be a non-empty run of ASCII digits, so
    // suffixes like "1.0-draft", "v1", "1a" are rejected consistently.
    let mut components = version_str.split('.');
    let major_str = components.next().unwrap_or("");
    let all_numeric = std::iter::once(major_str)
        .chain(components)
        .all(|c| !c.is_empty() && c.bytes().all(|b| b.is_ascii_digit()));
    if !all_numeric {
        return Err(unrecognized());
    }
    let major: u32 = major_str.parse().map_err(|_| unrecognized())?;

    if major > SUPPORTED_TRANSCRIPTS_JSON_MAJOR {
        return Err(FerroError::Json {
            msg: format!(
                "transcripts.json schema version {} is newer than this build supports \
                 (max major version {}); upgrade ferro-hgvs",
                version_str, SUPPORTED_TRANSCRIPTS_JSON_MAJOR
            ),
        });
    }

    Ok(())
}

impl JsonProvider {
    /// Create an empty mock provider
    pub fn new() -> Self {
        Self {
            transcripts: HashMap::new(),
            proteins: HashMap::new(),
            genomic_sequences: HashMap::new(),
            non_version_exact: HashSet::new(),
            genomic_placements: HashMap::new(),
            legacy_gene_models: HashMap::new(),
            version_substitutions: HashMap::new(),
            build_transcripts: HashMap::new(),
        }
    }

    /// Register a legacy gene-model → reference-standard transcript mapping
    /// (e.g. `"TIMM8B"` → `"NM_012459.4"`), keyed case-insensitively by gene
    /// Symbol, for the #500/#637 selector rewrite.
    pub fn add_legacy_gene_model(
        &mut self,
        gene_symbol: impl AsRef<str>,
        transcript: impl Into<String>,
    ) {
        self.legacy_gene_models
            .insert(gene_symbol.as_ref().to_ascii_uppercase(), transcript.into());
    }

    /// Register the chromosomal placement of a genomic parent reference
    /// (`NG_`/`LRG_`), keyed by its full versioned accession (#480). May be
    /// called more than once per accession to register a per-build placement
    /// (e.g. a GRCh37 and a GRCh38 placement), for #653 build-selection tests.
    pub fn add_genomic_placement(
        &mut self,
        parent_accession: impl Into<String>,
        placement: GenomicPlacement,
    ) {
        self.genomic_placements
            .entry(parent_accession.into())
            .or_default()
            .push(placement);
    }

    /// Mark `id` as not available at its exact requested version, so
    /// [`ReferenceProvider::has_transcript_version_exact`] returns `false` for
    /// it. Simulates a reference that only carries a sibling version's bases.
    pub fn mark_non_version_exact(&mut self, id: impl Into<String>) {
        self.non_version_exact.insert(id.into());
    }

    /// Model a silent version-strip substitution (#785): a lookup of `requested`
    /// resolves to (and serves the bases of) the *different*-version `served`
    /// accession, which must already be registered. Mirrors a real
    /// `MultiFastaProvider` that, lacking the requested version, falls back to a
    /// sibling. Used to exercise the normalization version-substitution gate.
    ///
    /// A substitution inherently means the exact requested version is not served,
    /// so this also marks `requested` as not version-exact (as
    /// [`Self::mark_non_version_exact`] would) — keeping
    /// [`ReferenceProvider::has_transcript_version_exact`] and
    /// [`ReferenceProvider::get_transcript`] mutually consistent for the gate.
    pub fn mark_version_substitution(
        &mut self,
        requested: impl Into<String>,
        served: impl Into<String>,
    ) {
        let requested = requested.into();
        self.non_version_exact.insert(requested.clone());
        self.version_substitutions.insert(requested, served.into());
    }

    /// Load reference data from a JSON file.
    ///
    /// Accepts either a bare array of `Transcript` records or an object
    /// of the form `{ transcripts, proteins, genomic_sequences }`.
    ///
    /// When `genomic_sequences` is present (the reference declares genomic
    /// capability), the placement of every placed transcript must be backed by
    /// bytes: its contig must appear in `genomic_sequences` and be long enough to
    /// contain the transcript's maximal genomic coordinate. A missing or too-short
    /// contig is rejected here rather than producing a silently-degraded
    /// genome-aware normalization (#1012 comment 1). See
    /// [`validate_genomic_placement_backing`].
    ///
    /// For the object form, an incompatible container `version` (a major newer
    /// than this build supports) is rejected at load — see
    /// [`validate_transcripts_json_version`].
    pub fn from_json(path: &Path) -> Result<Self, FerroError> {
        // `deny_unknown_fields` so a typo'd key (e.g. `transripts`) produces
        // a clear error rather than silently defaulting to an empty provider.
        #[derive(serde::Deserialize)]
        #[serde(deny_unknown_fields)]
        struct ObjectForm {
            #[serde(default)]
            transcripts: Vec<Transcript>,
            #[serde(default)]
            proteins: HashMap<String, String>,
            #[serde(default)]
            genomic_sequences: HashMap<String, String>,
            /// Container metadata emitted by `ferro convert-gff`; accepted but
            /// ignored here (validated from the raw value before this deserialize).
            /// Typed as `Value` so a string ("1.0") or number (1) both parse.
            #[serde(default, rename = "version")]
            _version: Option<serde_json::Value>,
            /// Container metadata emitted by `ferro convert-gff`; per-transcript
            /// `genome_build` on each `Transcript` is authoritative.
            #[serde(default, rename = "genome_build")]
            _genome_build: Option<serde_json::Value>,
        }

        // Parse straight from a buffered reader rather than reading the whole
        // file into a `String` first: a genome-capable transcripts.json can carry
        // whole-chromosome sequences, and `read_to_string` + `from_str` would hold
        // the raw text, the `Value` DOM, AND the final typed structures alive
        // simultaneously (~3x the file size at peak).
        let file = std::fs::File::open(path)?;
        let value: serde_json::Value = serde_json::from_reader(std::io::BufReader::new(file))?;

        // Fail loud on an incompatible container schema version rather than silently
        // loading a file this build cannot correctly interpret (#1012 comment 2).
        // The `MockProvider`/`reference_json` path previously did only
        // `deny_unknown_fields` and never checked the version, unlike the manifest
        // path (`validate_loaded_manifest`). Validated against the RAW value (before
        // the `deny_unknown_fields` deserialize) so a newer-major file gets the
        // friendly "upgrade" message even if it also carries unknown fields.
        // `Value::get` returns `None` for a bare array, which is accepted.
        validate_transcripts_json_version(value.get("version"))?;

        let (transcripts, proteins, genomic_sequences) = match value {
            serde_json::Value::Array(_) => {
                let transcripts: Vec<Transcript> = serde_json::from_value(value)?;
                (transcripts, HashMap::new(), HashMap::new())
            }
            serde_json::Value::Object(_) => {
                let obj: ObjectForm = serde_json::from_value(value)?;
                (obj.transcripts, obj.proteins, obj.genomic_sequences)
            }
            _ => {
                return Err(FerroError::Json {
                    msg: "JsonProvider JSON root must be an array or object".to_string(),
                })
            }
        };

        // A reference with no transcripts, proteins, or genomic sequences can
        // resolve nothing; fail loud rather than silently construct an empty,
        // useless provider (#1012 item 5). A genomic-only reference (zero
        // transcripts but populated `genomic_sequences`) is still valid for
        // `g.` normalization, so only a wholly-empty file is rejected.
        if transcripts.is_empty() && proteins.is_empty() && genomic_sequences.is_empty() {
            return Err(FerroError::Json {
                msg: "reference JSON has no usable reference data: it defines no \
                      transcripts, proteins, or genomic sequences"
                    .to_string(),
            });
        }

        // When the reference declares genomic capability (a non-empty
        // `genomic_sequences` map, so `has_genomic_data()` is true and the
        // genome-aware normalization rules will run), fail loud if that capability
        // is not actually backed by bytes for every placed transcript. Without this,
        // a reference that claims a genome but whose `genomic_sequences` does not
        // cover a transcript's placement produces a genome-aware result silently
        // computed against absent/short bytes (#1012 comment 1). We can only check
        // structure, not content, but a missing or too-short contig is a real,
        // detectable footgun.
        validate_genomic_placement_backing(&transcripts, &genomic_sequences)?;

        let map: HashMap<String, Arc<Transcript>> = transcripts
            .into_iter()
            .map(|tx| (tx.id.clone(), Arc::new(tx)))
            .collect();

        Ok(Self {
            transcripts: map,
            proteins,
            genomic_sequences,
            non_version_exact: HashSet::new(),
            genomic_placements: HashMap::new(),
            legacy_gene_models: HashMap::new(),
            version_substitutions: HashMap::new(),
            build_transcripts: HashMap::new(),
        })
    }

    /// Register a build-specific transcript so the provider serves
    /// build-divergent bases/CDS for the same `transcript_id` (#843). `build`
    /// is `"GRCh37"` / `"GRCh38"`. Consulted by
    /// [`ReferenceProvider::get_transcript_for_accession`] when the input
    /// accession carries a build-bearing `genomic_context` whose inferred build
    /// matches; the build-agnostic [`Self::add_transcript`] map is the fallback.
    pub fn add_transcript_on_build(&mut self, build: &'static str, transcript: Transcript) {
        self.build_transcripts
            .insert((transcript.id.clone(), build), Arc::new(transcript));
    }

    /// Add a transcript to the provider
    pub fn add_transcript(&mut self, transcript: Transcript) {
        self.transcripts
            .insert(transcript.id.clone(), Arc::new(transcript));
    }

    /// Add a protein sequence to the provider
    pub fn add_protein(&mut self, accession: impl Into<String>, sequence: impl Into<String>) {
        self.proteins.insert(accession.into(), sequence.into());
    }

    /// Add a genomic sequence for a contig/chromosome
    pub fn add_genomic_sequence(&mut self, contig: impl Into<String>, sequence: impl Into<String>) {
        self.genomic_sequences
            .insert(contig.into(), sequence.into());
    }

    /// Create a provider with some test transcripts
    pub fn with_test_data() -> Self {
        use crate::reference::transcript::{Exon, ManeStatus, Strand};

        let mut provider = Self::new();

        // Add a typical coding transcript (MANE Select)
        provider.add_transcript(Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some(
                "ATGCCCAAGGTGCTGCCCCAGATGCTGCCAGTGCTGCTGCTGCTGCTGCTGCTGCTGCTG".to_string(),
            ),
            cds_start: Some(1),
            cds_end: Some(60),
            exons: vec![
                Exon::new(1, 1, 20),
                Exon::new(2, 21, 40),
                Exon::new(3, 41, 60),
            ],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            cds_start_incomplete: false,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        // Add a transcript with UTRs
        provider.add_transcript(Transcript {
            id: "NM_001234.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("AAAAATGCCCAAGGGGGGGGGGGGGGGGGGGGGGGGGTAAAAAA".to_string()),
            cds_start: Some(5),
            cds_end: Some(38),
            exons: vec![
                Exon::new(1, 1, 15),
                Exon::new(2, 16, 30),
                Exon::new(3, 31, 44),
            ],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            cds_start_incomplete: false,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        // Add a minus strand transcript
        provider.add_transcript(Transcript {
            id: "NM_999999.1".to_string(),
            gene_symbol: Some("MINUS".to_string()),
            strand: Strand::Minus,
            sequence: Some("ATGCATGCATGCATGCATGCATGCATGCATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(30),
            exons: vec![Exon::new(1, 1, 32)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            cds_start_incomplete: false,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        // Add a transcript for testing 3' normalization of duplications
        // Sequence has "GAA" at c.8-c.9 to test that c.8dup -> c.9dup (2 A's)
        // and "GTTT" at c.21-c.24 to test that c.22dup -> c.24dup (3 T's)
        // Positions: 1234567890123456789012345678901234567890
        // Sequence:  ATGCCCGAAGCCCCCCCCCGTTTGCATGCATGCATGCAT
        //                  ^^           ^^^
        //                  c.8-9 (AA)   c.22-24 (TTT)
        provider.add_transcript(Transcript {
            id: "NM_888888.1".to_string(),
            gene_symbol: Some("DUPTEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGCCCGAAGCCCCCCCCCGTTTGCATGCATGCATGCAT".to_string()),
            cds_start: Some(1),
            cds_end: Some(39),
            exons: vec![Exon::new(1, 1, 39)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            cds_start_incomplete: false,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        // Add a non-coding transcript (NR_) for r.-axis cross-reference tests.
        // With no CDS, r. numbering is transcript-relative (== n.); positions
        // map directly to transcript bases (#773).
        provider.add_transcript(Transcript {
            id: "NR_000123.1".to_string(),
            gene_symbol: Some("NONCODING".to_string()),
            strand: Strand::Plus,
            sequence: Some("ACGTACGTACGT".to_string()),
            cds_start: None,
            cds_end: None,
            exons: vec![Exon::new(1, 1, 12)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            cds_start_incomplete: false,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        // Add test protein sequences for validation testing
        // NP_000079.2 is used for BRAF V600E-style testing
        // We need Val at position 600, Arg at 97 for frameshift tests
        // Building a sequence with known positions:
        // - Position 1: M (Met)
        // - Position 97: R (Arg) for frameshift test
        // - Position 600: V (Val) for V600E test
        // - Position 23-25: K, A, E for range tests
        let mut seq = String::new();
        // Positions 1-96 (96 chars): just padding with Ala
        seq.push('M'); // Position 1
        for _ in 2..23 {
            seq.push('A'); // Positions 2-22
        }
        seq.push('K'); // Position 23 - for Lys23
        seq.push('A'); // Position 24 - for middle of range
        seq.push('E'); // Position 25 - for Glu25
        for _ in 26..97 {
            seq.push('A'); // Positions 26-96
        }
        seq.push('R'); // Position 97 - for Arg97 frameshift
        for _ in 98..600 {
            seq.push('A'); // Positions 98-599
        }
        seq.push('V'); // Position 600 - for Val600
        for _ in 601..=700 {
            seq.push('A'); // Positions 601-700
        }
        provider.add_protein("NP_000079.2", seq);

        // A simpler test protein for easier position testing
        // Position 1=M, 2=V, 3=L, 4=S, 5=P, etc.
        provider.add_protein("NP_TEST.1", "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFL");

        provider
    }

    /// Get the number of transcripts
    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    /// Total number of genomic sequence bases across all contigs. Zero when the
    /// reference carries no genomic sequence, or only empty contig strings — so a
    /// caller can distinguish real genomic capability from a merely-present but
    /// empty `genomic_sequences` map (which `has_genomic_data()` still reports true).
    pub fn total_genomic_bases(&self) -> u64 {
        self.genomic_sequences
            .values()
            .map(|s| s.len() as u64)
            .sum()
    }

    /// Check if provider is empty
    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }

    /// Get all transcript IDs
    pub fn transcript_ids(&self) -> Vec<&str> {
        self.transcripts.keys().map(|s| s.as_str()).collect()
    }

    /// Iterate over all transcripts registered in the provider.
    pub fn all_transcripts(&self) -> impl Iterator<Item = &Transcript> {
        self.transcripts.values().map(Arc::as_ref)
    }

    /// True when `id` is declared as a chromosome/contig name — either an
    /// entry in `genomic_sequences` or the `chromosome` field on any
    /// transcript. Used by `get_sequence` to disambiguate genomic lookups
    /// from transcript-id lookups when the two name spaces collide
    /// (e.g. a transcript whose id is `chr1-gene.1`; see #311).
    fn is_known_contig(&self, id: &str) -> bool {
        if self.genomic_sequences.contains_key(id) {
            return true;
        }
        self.transcripts
            .values()
            .any(|tx| tx.chromosome.as_deref() == Some(id))
    }
}

impl Default for JsonProvider {
    fn default() -> Self {
        Self::new()
    }
}

impl ReferenceProvider for JsonProvider {
    fn genomic_placement(&self, parent: &Accession) -> Option<GenomicPlacement> {
        self.genomic_placement_on_build(parent, None)
    }

    fn genomic_placement_on_build(
        &self,
        parent: &Accession,
        build: Option<&str>,
    ) -> Option<GenomicPlacement> {
        let list = self.genomic_placements.get(&parent.full())?;
        crate::reference::provider::select_placement_for_build(
            list,
            build,
            crate::liftover::aliases::infer_genome_build_from_accession,
        )
    }

    fn resolve_legacy_gene_selector(
        &self,
        selector: &str,
        _ng_parent: Option<&crate::hgvs::variant::Accession>,
    ) -> Option<String> {
        crate::reference::legacy_selector::resolve_legacy_selector_in(selector, |g| {
            self.legacy_gene_models.get(g).cloned()
        })
    }

    /// Build-aware transcript resolution (#843). When the input accession
    /// carries a build-bearing `genomic_context` (an `NC_*.10`/`.11` parent)
    /// and a build-specific transcript was registered via
    /// [`Self::add_transcript_on_build`], serve that build's record so a test
    /// can model build-divergent transcript bases. Otherwise fall back to the
    /// build-agnostic [`Self::get_transcript`] — preserving existing behavior
    /// for every input that does not opt into build-keyed transcripts.
    ///
    /// To avoid silently masking an incomplete test setup, the fallback is
    /// **refused for a partially build-keyed transcript**: if `tx_id` has a
    /// build-keyed record for *some* build but not the requested one, return
    /// [`FerroError::ReferenceNotFound`] rather than serving the build-agnostic
    /// record (which would hide a wrong/missing `add_transcript_on_build` call).
    /// A transcript with **no** build-keyed records keeps the build-agnostic
    /// path, so cross-id fallback and the deliberate primary-build fixture
    /// entry are unaffected.
    fn get_transcript_for_accession(
        &self,
        accession: &Accession,
    ) -> Result<Arc<Transcript>, FerroError> {
        if !self.build_transcripts.is_empty() {
            if let Some(parent) = accession.genomic_context.as_deref() {
                if let Some(build) = self.infer_genome_build(parent) {
                    let tx_id = accession.transcript_accession();
                    // Apply the same unversioned → versioned base-accession rule
                    // [`Self::get_transcript`] uses (`NM_123` resolves to a
                    // stored `NM_123.1`): an exact-key lookup would miss a
                    // build-keyed versioned record for a bare/unversioned input
                    // accession and silently fall back to the build-agnostic
                    // (wrong-build) bases — the very bug #843 fixes. Only an
                    // unversioned query bridges to a versioned key; a versioned
                    // miss never resolves a different version.
                    let matches_tx_id = |keyed_id: &str| {
                        keyed_id == tx_id.as_str()
                            || (!tx_id.contains('.')
                                && keyed_id.split('.').next().unwrap_or(keyed_id) == tx_id.as_str())
                    };
                    if let Some(tx) =
                        self.build_transcripts
                            .iter()
                            .find_map(|((keyed_id, keyed_build), tx)| {
                                (matches_tx_id(keyed_id) && *keyed_build == build)
                                    .then(|| Arc::clone(tx))
                            })
                    {
                        return Ok(tx);
                    }
                    // `tx_id` is build-keyed on a different build but not this
                    // one: a partial setup. Surface it rather than silently
                    // serving the build-agnostic record.
                    let partially_keyed = self
                        .build_transcripts
                        .keys()
                        .any(|(keyed_id, _)| matches_tx_id(keyed_id));
                    if partially_keyed {
                        return Err(FerroError::ReferenceNotFound {
                            id: format!("{tx_id} (no build-keyed record for {build})"),
                        });
                    }
                }
            }
        }
        self.get_transcript(&accession.transcript_accession())
    }

    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError> {
        // Handle versioned and unversioned lookups
        if let Some(tx) = self.transcripts.get(id) {
            return Ok(Arc::clone(tx));
        }

        // Modeled silent version substitution (#785): serve a configured
        // sibling version's record (whose own `id` differs from `id`), mirroring
        // a real `MultiFastaProvider`'s base→latest fallback. Only fires when a
        // test explicitly opts in via `mark_version_substitution`.
        if let Some(served) = self.version_substitutions.get(id) {
            if let Some(tx) = self.transcripts.get(served) {
                return Ok(Arc::clone(tx));
            }
        }

        // Try without version: only match a stored key whose base id (the
        // segment before the version dot) equals the requested id. A bare
        // `starts_with` would also match unrelated keys that merely share
        // a prefix (e.g. `chr1-gene.1` for a `chr1` lookup), causing a
        // genomic accession to be resolved as a transcript — see #311.
        //
        // Gated on an unversioned query: a versioned miss (`NM_000088.5`
        // against a stored `NM_000088.3`) must NOT silently return a
        // different version. The fallback only bridges `NM_000088` →
        // `NM_000088.3`.
        if !id.contains('.') {
            for (key, tx) in &self.transcripts {
                if key.split('.').next().unwrap_or(key) == id {
                    return Ok(Arc::clone(tx));
                }
            }
        }

        Err(FerroError::ReferenceNotFound { id: id.to_string() })
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        // If the requested id is registered as a contig — either as a key
        // in `genomic_sequences` or as the `chromosome` field on any
        // transcript — resolve it strictly against the genomic path.
        // Without this dispatch, a transcript whose own id (or, before
        // #311, whose id merely shared a prefix with the contig name)
        // would short-circuit the lookup and silently interpret genomic
        // coordinates as transcript-relative indices.
        if self.is_known_contig(id) {
            return self.get_genomic_sequence(id, start, end);
        }

        // Otherwise treat the id as a transcript accession (matches
        // FastaProvider's transcript-first ordering for transcript ids).
        if let Ok(transcript) = self.get_transcript(id) {
            return transcript
                .get_sequence(start, end)
                .map(|s| s.to_string())
                .ok_or_else(|| FerroError::InvalidCoordinates {
                    msg: format!("Position {}-{} out of range for {}", start, end, id),
                });
        }

        // Fall through to contig/chromosome lookup so genomic accessions
        // resolve against `genomic_sequences` instead of returning
        // ReferenceNotFound.
        self.get_genomic_sequence(id, start, end)
    }

    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        // Handle versioned and unversioned lookups
        let protein_seq = if let Some(seq) = self.proteins.get(accession) {
            seq.clone()
        } else if !accession.contains('.') {
            // Try without version: match a stored key whose base accession
            // (the segment before the version dot) equals the requested
            // accession. Strict base-id equality avoids the prefix-collision
            // failure mode that `get_transcript` had pre-#311.
            //
            // Gated on an unversioned query — a versioned miss must not
            // cross-version fall back (mirrors `get_transcript` above).
            let mut found = None;
            for (key, seq) in &self.proteins {
                if key.split('.').next().unwrap_or(key) == accession {
                    found = Some(seq.clone());
                    break;
                }
            }
            found.ok_or_else(|| FerroError::ProteinReferenceNotAvailable {
                accession: accession.to_string(),
                start,
                end,
            })?
        } else {
            return Err(FerroError::ProteinReferenceNotAvailable {
                accession: accession.to_string(),
                start,
                end,
            });
        };

        let start = start as usize;
        let end = end as usize;

        if start >= protein_seq.len() || end > protein_seq.len() || start > end {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "Position {}-{} out of range for {} (length {})",
                    start,
                    end,
                    accession,
                    protein_seq.len()
                ),
            });
        }

        Ok(protein_seq[start..end].to_string())
    }

    fn get_protein_length(&self, accession: &str) -> Result<u64, FerroError> {
        // Mirror the accession resolution in `get_protein_sequence`:
        // exact (versioned or stored) match first, then an unversioned
        // base-accession fallback, returning the stored length directly.
        if let Some(seq) = self.proteins.get(accession) {
            return Ok(seq.len() as u64);
        }
        if !accession.contains('.') {
            for (key, seq) in &self.proteins {
                if key.split('.').next().unwrap_or(key) == accession {
                    return Ok(seq.len() as u64);
                }
            }
        }
        // An accession this provider cannot resolve maps to length `0`,
        // matching the trait contract (the default probe leaves `lo = 0`
        // for an unresolvable protein). Returning an error here would
        // diverge from that contract.
        Ok(0)
    }

    fn has_protein_data(&self) -> bool {
        !self.proteins.is_empty()
    }

    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        let genomic_seq = self.genomic_sequences.get(contig).ok_or_else(|| {
            FerroError::GenomicReferenceNotAvailable {
                contig: contig.to_string(),
                start,
                end,
            }
        })?;

        let start = start as usize;
        let end = end as usize;

        if start >= genomic_seq.len() || end > genomic_seq.len() || start > end {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "Position {}-{} out of range for {} (length {})",
                    start,
                    end,
                    contig,
                    genomic_seq.len()
                ),
            });
        }

        Ok(genomic_seq[start..end].to_string())
    }

    fn has_genomic_data(&self) -> bool {
        !self.genomic_sequences.is_empty()
    }

    fn get_sequence_length(&self, id: &str) -> Result<u64, FerroError> {
        // Mirror the routing of `get_sequence` (and the real
        // `MultiFastaProvider`, which serves a length for any indexed sequence —
        // genomic *or* transcript). Tests that need a transcript's true length
        // (e.g. the #797 poly-A walk) rely on the transcript fallback; without
        // it Mock returns `ReferenceNotFound` for a transcript id and diverges
        // from the production contract.
        //
        // 1. A known contig resolves strictly against `genomic_sequences`,
        //    mirroring `get_sequence`, which returns from the contig branch
        //    immediately. Without the early error a known contig with no loaded
        //    genomic sequence would fall through to the transcript lookup below
        //    and could return a transcript length for an ambiguous contig id.
        if self.is_known_contig(id) {
            if let Some(s) = self.genomic_sequences.get(id) {
                return Ok(s.len() as u64);
            }
            return Err(FerroError::ReferenceNotFound { id: id.to_string() });
        }
        // 2. Otherwise treat it as a transcript accession, via `get_transcript`
        //    so the unversioned-fallback / version-substitution paths apply
        //    exactly as in `get_sequence`.
        if let Ok(tx) = self.get_transcript(id) {
            if let Some(seq) = &tx.sequence {
                return Ok(seq.len() as u64);
            }
        }
        // 3. Fall through to a direct genomic lookup (a genomic accession that
        //    is not flagged as a known contig).
        if let Some(s) = self.genomic_sequences.get(id) {
            return Ok(s.len() as u64);
        }
        Err(FerroError::ReferenceNotFound { id: id.to_string() })
    }

    fn has_transcript_version_exact(&self, id: &str) -> bool {
        !self.non_version_exact.contains(id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mock_provider_with_test_data() {
        let provider = JsonProvider::with_test_data();
        assert!(!provider.is_empty());
        assert!(provider.len() >= 2);
    }

    /// The former name `MockProvider` remains a working alias for [`JsonProvider`]
    /// so existing code keeps compiling after the rename.
    #[test]
    fn mock_provider_alias_is_json_provider() {
        let via_alias: MockProvider = MockProvider::with_test_data();
        let via_new: JsonProvider = via_alias; // same type — assignment compiles
        assert!(!via_new.is_empty());
    }

    /// #653: an NG_ parent with both a GRCh37 and a GRCh38 RefSeqGene placement
    /// resolves per build — GRCh38 by default, GRCh37 on request — and an
    /// explicit build with no matching placement declines (build-match guard).
    #[test]
    fn genomic_placement_selects_per_build() {
        use crate::reference::provider::GenomicPlacement;
        use crate::reference::transcript::Strand;
        let mk = |chrom_version: u32| GenomicPlacement {
            nc: Accession::new("NC", "000001", Some(chrom_version)),
            parent_start: 1,
            nc_start: 1000,
            nc_end: 1008,
            strand: Strand::Plus,
        };
        let mut provider = JsonProvider::new();
        provider.add_genomic_placement("NG_900.1", mk(11)); // GRCh38
        provider.add_genomic_placement("NG_900.1", mk(10)); // GRCh37
        let ng = Accession::new("NG", "900", Some(1));

        // Default prefers GRCh38; explicit builds select their own placement.
        assert_eq!(
            provider.genomic_placement(&ng).unwrap().nc.to_string(),
            "NC_000001.11"
        );
        assert_eq!(
            provider
                .genomic_placement_on_build(&ng, Some("GRCh37"))
                .unwrap()
                .nc
                .to_string(),
            "NC_000001.10"
        );
        assert_eq!(
            provider
                .genomic_placement_on_build(&ng, Some("GRCh38"))
                .unwrap()
                .nc
                .to_string(),
            "NC_000001.11"
        );

        // An NG_ with only a GRCh37 placement declines a GRCh38 request (guard).
        let mut grch37_only = JsonProvider::new();
        grch37_only.add_genomic_placement("NG_901.1", mk(10));
        let ng2 = Accession::new("NG", "901", Some(1));
        assert!(grch37_only
            .genomic_placement_on_build(&ng2, Some("GRCh38"))
            .is_none());
        // ...but resolves under GRCh37 and as the default-when-sole-coverage.
        assert!(grch37_only
            .genomic_placement_on_build(&ng2, Some("GRCh37"))
            .is_some());
        assert!(grch37_only.genomic_placement(&ng2).is_some());
    }

    /// #843: the build-keyed transcript lookup must honor the same unversioned
    /// → versioned base-accession rule that [`JsonProvider::get_transcript`]
    /// applies (`NM_123` resolves to a stored `NM_123.1`). Without it, a bare
    /// `g.` allele whose transcript accession is unversioned would miss the
    /// build-specific record and silently fall back to the build-agnostic
    /// (primary) bases — exactly the wrong-build serving this fix exists to
    /// prevent. Register a versioned record on GRCh37, query with the
    /// unversioned accession carrying a GRCh37 `genomic_context`, and assert it
    /// resolves to that build's record rather than the build-agnostic one.
    #[test]
    fn build_keyed_lookup_resolves_unversioned_accession() {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand};

        let mk = |build: GenomeBuild, seq: &str| Transcript {
            id: "NM_555555.1".to_string(),
            gene_symbol: Some("BUILDDIV".to_string()),
            strand: Strand::Plus,
            sequence: Some(seq.to_string()),
            cds_start: Some(1),
            cds_end: Some(seq.len() as u64),
            exons: vec![Exon::new(1, 1, seq.len() as u64)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: build,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            cds_start_incomplete: false,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let mut provider = JsonProvider::new();
        // Build-keyed records (versioned ids) with build-divergent bases.
        provider.add_transcript_on_build("GRCh37", mk(GenomeBuild::GRCh37, "AAAA"));
        provider.add_transcript_on_build("GRCh38", mk(GenomeBuild::GRCh38, "CCCC"));
        // Build-agnostic fallback = the primary (GRCh38) bases, mirroring how a
        // real provider serves the primary record when no build is keyed.
        provider.add_transcript(mk(GenomeBuild::GRCh38, "CCCC"));

        // Unversioned transcript accession (`NM_555555`, version `None`) carrying
        // a build-bearing GRCh37 `genomic_context` (`NC_000001.10`).
        let unversioned = Accession::new("NM", "555555", None)
            .with_genomic_context(Accession::new("NC", "000001", Some(10)));
        let tx = provider
            .get_transcript_for_accession(&unversioned)
            .expect("unversioned accession must resolve to the GRCh37 build record");
        assert_eq!(
            tx.sequence.as_deref(),
            Some("AAAA"),
            "must serve the GRCh37 build-keyed bases, not the build-agnostic (GRCh38) fallback"
        );

        // GRCh38 context resolves to its own build's bases.
        let unversioned_38 = Accession::new("NM", "555555", None)
            .with_genomic_context(Accession::new("NC", "000001", Some(11)));
        let tx38 = provider
            .get_transcript_for_accession(&unversioned_38)
            .expect("unversioned accession must resolve to the GRCh38 build record");
        assert_eq!(tx38.sequence.as_deref(), Some("CCCC"));
    }

    #[test]
    fn test_get_transcript() {
        let provider = JsonProvider::with_test_data();
        let tx = provider.get_transcript("NM_000088.3").unwrap();
        assert_eq!(tx.gene_symbol, Some("COL1A1".to_string()));
    }

    #[test]
    fn test_get_transcript_not_found() {
        let provider = JsonProvider::with_test_data();
        let result = provider.get_transcript("NM_NONEXISTENT.1");
        assert!(result.is_err());
    }

    #[test]
    fn with_test_data_includes_a_noncoding_transcript() {
        let provider = JsonProvider::with_test_data();
        let tx = provider
            .get_transcript("NR_000123.1")
            .expect("NR_000123.1 should be present in test data");
        assert!(!tx.is_coding(), "NR_ transcript must be non-coding");
        assert_eq!(tx.sequence.as_deref(), Some("ACGTACGTACGT"));
    }

    #[test]
    fn test_get_sequence() {
        let provider = JsonProvider::with_test_data();
        let seq = provider.get_sequence("NM_000088.3", 0, 3).unwrap();
        assert_eq!(seq, "ATG");
    }

    #[test]
    fn test_get_sequence_falls_through_to_contig() {
        // Regression: get_sequence should fall through to contig lookup
        // when the id is not a transcript, matching FastaProvider behavior.
        let mut provider = JsonProvider::new();
        provider.add_genomic_sequence("chr1", "ACGTACGT");
        let seq = provider.get_sequence("chr1", 0, 4).unwrap();
        assert_eq!(seq, "ACGT");
    }

    #[test]
    fn test_has_transcript() {
        let provider = JsonProvider::with_test_data();
        assert!(provider.has_transcript("NM_000088.3"));
        assert!(!provider.has_transcript("NONEXISTENT"));
    }

    #[test]
    fn boxed_provider_forwards_has_transcript_version_exact() {
        // The protein gate (#505) is consulted through the provider type the
        // production CLI/benchmark paths use: `Box<dyn ReferenceProvider>`. If
        // the blanket boxed impl does not forward `has_transcript_version_exact`,
        // it silently falls through to the trait default `true` and the gate is
        // inert in production. Pin the forwarding here.
        let mut inner = JsonProvider::new();
        inner.mark_non_version_exact("NM_TEST.1");
        let boxed: Box<dyn ReferenceProvider> = Box::new(inner);
        assert!(
            !boxed.has_transcript_version_exact("NM_TEST.1"),
            "boxed provider must forward to the inner provider's verdict (false), \
             not the trait default (true)"
        );
        assert!(boxed.has_transcript_version_exact("NM_OTHER.1"));
    }

    #[test]
    fn test_from_json_object_form_with_proteins_and_genomic() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"{
      "transcripts": [{
        "id": "NM_TEST.1",
        "gene_symbol": "TEST",
        "strand": "+",
        "sequence": "ATGCATGCAT",
        "cds_start": 1,
        "cds_end": 10,
        "exons": [{"number": 1, "start": 1, "end": 10}]
      }],
      "proteins": {
        "NP_TEST.1": "MAPLE"
      },
      "genomic_sequences": {
        "chr1": "ACGTACGTACGT"
      }
    }"#;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider = JsonProvider::from_json(file.path()).unwrap();

        assert!(provider.has_transcript("NM_TEST.1"));
        assert!(provider.has_protein_data());
        assert_eq!(
            provider.get_protein_sequence("NP_TEST.1", 0, 5).unwrap(),
            "MAPLE"
        );
        assert!(provider.has_genomic_data());
        assert_eq!(provider.get_genomic_sequence("chr1", 0, 4).unwrap(), "ACGT");
    }

    /// A genomic reference whose `genomic_sequences` fully backs the placement of
    /// every placed transcript loads and reports genomic capability (#1026).
    #[test]
    fn from_json_accepts_consistent_genomic_reference() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // chr1 genomic sequence of length 12; the transcript exon spans genomic 3..8,
        // so its sequence must be the contig bytes there: "ACGTACGTACGT"[2..8] = "GTACGT".
        let json = r#"{
      "transcripts": [{
        "id": "NM_TEST.1",
        "strand": "+",
        "sequence": "GTACGT",
        "chromosome": "chr1",
        "genomic_start": 3,
        "genomic_end": 8,
        "exons": [{"number": 1, "start": 1, "end": 6, "genomic_start": 3, "genomic_end": 8}]
      }],
      "genomic_sequences": { "chr1": "ACGTACGTACGT" }
    }"#;
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider = JsonProvider::from_json(file.path()).expect("consistent reference loads");
        assert!(provider.has_genomic_data());
    }

    /// A reference whose `genomic_sequences` are long enough but do NOT reconstruct
    /// the transcript's stated sequence (wrong/junk bases) is rejected at load — the
    /// wrong-content half of the #1012 comment 1 footgun. Without this, genome-aware
    /// normalization of del/dup/inv/intronic edits (which carry no stated reference
    /// bases) would silently run against the wrong genome.
    #[test]
    fn from_json_rejects_genomic_content_that_does_not_match_transcript() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // The exon at genomic 3..8 reconstructs to "GTACGT", but the transcript
        // claims "ACGTAC" — a plausible-but-wrong (junk) genome.
        let json = r#"{
      "transcripts": [{
        "id": "NM_TEST.1",
        "strand": "+",
        "sequence": "ACGTAC",
        "chromosome": "chr1",
        "genomic_start": 3,
        "genomic_end": 8,
        "exons": [{"number": 1, "start": 1, "end": 6, "genomic_start": 3, "genomic_end": 8}]
      }],
      "genomic_sequences": { "chr1": "ACGTACGTACGT" }
    }"#;
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let err = match JsonProvider::from_json(file.path()) {
            Ok(_) => panic!("expected rejection of genomic bytes that don't match the transcript"),
            Err(e) => e,
        };
        assert!(
            format!("{err}").contains("do not reconstruct transcript"),
            "error should explain the content mismatch: {err}"
        );
    }

    /// The content check reverse-complements the genomic bytes for a minus-strand
    /// transcript, so a correctly-placed minus-strand reference loads.
    #[test]
    fn from_json_accepts_minus_strand_reconstruction() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // contig[0..6] = "ACGTAC"; on the minus strand the transcript sequence is its
        // reverse-complement: revcomp("ACGTAC") = "GTACGT".
        let json = r#"{
      "transcripts": [{
        "id": "NM_TEST.1",
        "strand": "-",
        "sequence": "GTACGT",
        "chromosome": "chr1",
        "genomic_start": 1,
        "genomic_end": 6,
        "exons": [{"number": 1, "start": 1, "end": 6, "genomic_start": 1, "genomic_end": 6}]
      }],
      "genomic_sequences": { "chr1": "ACGTACGTACGT" }
    }"#;
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider =
            JsonProvider::from_json(file.path()).expect("consistent minus-strand reference loads");
        assert!(provider.has_genomic_data());
    }

    /// A reference that declares genomic capability but has no contig entry backing a
    /// placed transcript is rejected at load, rather than silently normalizing against
    /// absent bytes (#1012 comment 1).
    #[test]
    fn from_json_rejects_placement_without_backing_contig() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // Transcript is placed on chr1, but genomic_sequences only carries chr2.
        let json = r#"{
      "transcripts": [{
        "id": "NM_TEST.1",
        "strand": "+",
        "sequence": "ATGCAT",
        "chromosome": "chr1",
        "genomic_start": 3,
        "genomic_end": 8,
        "exons": [{"number": 1, "start": 1, "end": 6, "genomic_start": 3, "genomic_end": 8}]
      }],
      "genomic_sequences": { "chr2": "ACGTACGTACGT" }
    }"#;
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let err = match JsonProvider::from_json(file.path()) {
            Ok(_) => panic!("expected an error for an unbacked placement"),
            Err(e) => e,
        };
        assert!(
            format!("{err}").contains("chr1"),
            "error should name the unbacked contig: {err}"
        );
    }

    /// A contig present but too short to contain the placement is rejected.
    #[test]
    fn from_json_rejects_too_short_genomic_sequence() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // Exon needs genomic coordinate 50 but chr1 is only 12 bases long.
        let json = r#"{
      "transcripts": [{
        "id": "NM_TEST.1",
        "strand": "+",
        "sequence": "ATGCAT",
        "chromosome": "chr1",
        "genomic_start": 45,
        "genomic_end": 50,
        "exons": [{"number": 1, "start": 1, "end": 6, "genomic_start": 45, "genomic_end": 50}]
      }],
      "genomic_sequences": { "chr1": "ACGTACGTACGT" }
    }"#;
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let err = match JsonProvider::from_json(file.path()) {
            Ok(_) => panic!("expected an error for a too-short contig"),
            Err(e) => e,
        };
        assert!(
            format!("{err}").contains("too short"),
            "error should explain the contig is too short: {err}"
        );
    }

    /// A transcripts-only reference (no `genomic_sequences`) is not subject to the
    /// placement-backing check even when transcripts carry genomic coordinates, so
    /// the pre-#1026 transcript-only references keep loading unchanged.
    #[test]
    fn from_json_skips_backing_check_without_genomic_sequences() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"{
      "transcripts": [{
        "id": "NM_TEST.1",
        "strand": "+",
        "sequence": "ATGCAT",
        "chromosome": "chr1",
        "genomic_start": 3,
        "genomic_end": 8,
        "exons": [{"number": 1, "start": 1, "end": 6, "genomic_start": 3, "genomic_end": 8}]
      }]
    }"#;
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider = JsonProvider::from_json(file.path()).expect("transcripts-only loads");
        assert!(!provider.has_genomic_data());
    }

    #[test]
    fn test_get_protein_length_resolves_and_falls_back() {
        let mut provider = JsonProvider::new();
        provider.add_protein("NP_TEST.1", "MAPLE");

        // Exact (versioned) match returns the stored length.
        assert_eq!(provider.get_protein_length("NP_TEST.1").unwrap(), 5);
        // Unversioned base accession falls back to the versioned entry.
        assert_eq!(provider.get_protein_length("NP_TEST").unwrap(), 5);
    }

    #[test]
    fn test_get_protein_length_unresolvable_is_zero() {
        // Contract: an accession the provider cannot resolve maps to a
        // length of `0` (matching the trait's default probe semantics),
        // NOT an error — otherwise `normalize()` is not behavior-preserving
        // when a provider switches from the probe loop to this API.
        let provider = JsonProvider::with_test_data();
        assert_eq!(provider.get_protein_length("NP_NONEXISTENT.1").unwrap(), 0);
    }

    #[test]
    fn test_from_json_bare_array_form_still_works() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"[{
      "id": "NM_TEST.1",
      "gene_symbol": "TEST",
      "strand": "+",
      "sequence": "ATGCATGCAT",
      "cds_start": 1,
      "cds_end": 10,
      "exons": [{"number": 1, "start": 1, "end": 10}]
    }]"#;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider = JsonProvider::from_json(file.path()).unwrap();

        assert!(provider.has_transcript("NM_TEST.1"));
        assert!(!provider.has_protein_data());
        assert!(!provider.has_genomic_data());
    }

    #[test]
    fn test_from_json_empty_object_is_rejected() {
        // An empty object has no usable reference data, so it must be rejected
        // rather than loaded as a useless empty provider (#1012 item 5).
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(b"{}").unwrap();

        match JsonProvider::from_json(file.path()) {
            Err(FerroError::Json { msg }) => assert!(
                msg.contains("no usable reference data"),
                "expected error to mention no usable reference data, got: {msg}",
            ),
            Err(e) => panic!("expected FerroError::Json, got {e}"),
            Ok(_) => panic!("expected a wholly-empty reference file to be rejected"),
        }
    }

    #[test]
    fn from_json_rejects_wholly_empty_forms() {
        // A bare empty array and an explicit `"transcripts": []` with no other
        // data are equally unusable and must be rejected (#1012 item 5).
        use std::io::Write;
        use tempfile::NamedTempFile;

        for body in [&b"[]"[..], br#"{"transcripts": []}"#] {
            let mut file = NamedTempFile::new().unwrap();
            file.write_all(body).unwrap();
            match JsonProvider::from_json(file.path()) {
                Err(FerroError::Json { msg }) => assert!(
                    msg.contains("no usable reference data"),
                    "expected error to mention no usable reference data, got: {msg}",
                ),
                Err(e) => panic!("expected FerroError::Json, got {e}"),
                Ok(_) => panic!("expected a wholly-empty reference file to be rejected"),
            }
        }
    }

    #[test]
    fn from_json_accepts_genomic_only_reference() {
        // A reference with zero transcripts but populated `genomic_sequences` is
        // still valid for `g.` normalization and must load (#1012 item 5): the
        // empty-file guard only rejects a wholly-empty reference.
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"{
            "transcripts": [],
            "genomic_sequences": {"NC_000001.11": "ACGTACGTACGT"}
        }"#;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider =
            JsonProvider::from_json(file.path()).expect("a genomic-only reference should load");
        assert!(provider.is_empty()); // no transcripts
        assert!(provider.has_genomic_data());
        assert_eq!(
            provider.get_genomic_sequence("NC_000001.11", 0, 4).unwrap(),
            "ACGT"
        );
    }

    #[test]
    fn from_json_accepts_protein_only_reference() {
        // A reference with zero transcripts and no genomic sequence but populated
        // `proteins` is still usable (protein queries) and must load — the
        // empty-file guard only rejects a wholly-empty reference (#1012 item 5).
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"{
            "transcripts": [],
            "proteins": {"NP_000001.1": "MAPLE"}
        }"#;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider =
            MockProvider::from_json(file.path()).expect("a protein-only reference should load");
        assert!(provider.is_empty()); // no transcripts
        assert!(!provider.has_genomic_data());
        assert!(provider.has_protein_data());
    }

    #[test]
    fn test_from_json_rejects_scalar_root() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(b"42").unwrap();

        match JsonProvider::from_json(file.path()) {
            Err(FerroError::Json { msg }) => {
                assert!(
                    msg.contains("array or object"),
                    "expected error message to mention 'array or object', got: {msg}",
                );
            }
            Err(e) => panic!("expected FerroError::Json, got {e}"),
            Ok(_) => panic!("expected scalar JSON root to be rejected"),
        }
    }

    #[test]
    fn test_from_json_rejects_unknown_field() {
        // A typo'd top-level key must error rather than silently parse as an
        // empty provider.
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(br#"{"transripts": []}"#).unwrap();

        match JsonProvider::from_json(file.path()) {
            Err(e) => assert!(
                format!("{e}").contains("transripts"),
                "expected error to mention the unknown field, got {e}",
            ),
            Ok(_) => panic!("expected unknown field to be rejected"),
        }
    }

    #[test]
    fn from_json_accepts_convert_gff_output_with_metadata_keys() {
        // The container shape produced by `ferro convert-gff`: the `version` and
        // `genome_build` container-metadata keys must be accepted (not rejected
        // by `deny_unknown_fields`) alongside transcript records.
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"{
            "version": "1.0",
            "genome_build": "GRCh38",
            "transcripts": [
                {
                    "id": "NM_000001.1",
                    "strand": "+",
                    "sequence": "ATGC",
                    "exons": [{"number": 1, "start": 1, "end": 4}]
                }
            ]
        }"#;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider = JsonProvider::from_json(file.path())
            .expect("convert-gff JSON with version/genome_build metadata should load");
        assert!(provider.has_transcript("NM_000001.1"));
    }

    #[test]
    fn validate_transcripts_json_version_accepts_supported_and_missing() {
        use serde_json::json;
        // Missing (bare array or object without `version`) is accepted.
        validate_transcripts_json_version(None).expect("missing version accepted");
        // Current major, any (value-only) minor. String or number form.
        validate_transcripts_json_version(Some(&json!("1.0"))).expect("1.0 accepted");
        validate_transcripts_json_version(Some(&json!("1.5"))).expect("newer minor accepted");
        validate_transcripts_json_version(Some(&json!("1"))).expect("bare major accepted");
        validate_transcripts_json_version(Some(&json!(1))).expect("numeric 1 accepted");
        validate_transcripts_json_version(Some(&json!(1.0))).expect("numeric 1.0 accepted");
    }

    #[test]
    fn validate_transcripts_json_version_rejects_newer_major_and_garbage() {
        use serde_json::json;
        let newer = validate_transcripts_json_version(Some(&json!("2.0"))).unwrap_err();
        assert!(
            format!("{newer}").contains("newer than this build supports"),
            "got: {newer}"
        );
        // A newer major as a number is caught too.
        let newer_num = validate_transcripts_json_version(Some(&json!(2))).unwrap_err();
        assert!(format!("{newer_num}").contains("newer than this build supports"));
        // Non-numeric components / suffixes are rejected consistently.
        for bad in [
            json!("banana"),
            json!("v1"),
            json!("1.0-draft"),
            json!("1a"),
            json!(""),
        ] {
            let err = validate_transcripts_json_version(Some(&bad)).unwrap_err();
            assert!(
                format!("{err}").contains("unrecognized transcripts.json schema version"),
                "expected unrecognized for {bad}, got: {err}"
            );
        }
        // A non-string, non-number version is rejected.
        let wrong_type = validate_transcripts_json_version(Some(&json!(["1.0"]))).unwrap_err();
        assert!(format!("{wrong_type}").contains("unrecognized"));
    }

    #[test]
    fn from_json_rejects_incompatible_container_version() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"{
            "version": "2.0",
            "genome_build": "GRCh38",
            "transcripts": []
        }"#;
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let err = match MockProvider::from_json(file.path()) {
            Ok(_) => panic!("expected an error for an incompatible schema version"),
            Err(e) => e,
        };
        assert!(
            format!("{err}").contains("newer than this build supports"),
            "got: {err}"
        );
    }

    #[test]
    fn from_json_reports_version_before_unknown_field_for_newer_major() {
        // A newer-major file that ALSO carries a field this build doesn't know must
        // surface the friendly "upgrade" message (version validated on the raw value
        // before `deny_unknown_fields`), not a cryptic "unknown field" error.
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"{ "version": "2.0", "brand_new_v2_field": 42, "transcripts": [] }"#;
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let err = match MockProvider::from_json(file.path()) {
            Ok(_) => panic!("expected rejection of a newer-major schema"),
            Err(e) => e,
        };
        assert!(
            format!("{err}").contains("newer than this build supports"),
            "newer-major must report the version error, not unknown-field: {err}"
        );
    }

    #[test]
    fn from_json_accepts_numeric_version() {
        // `"version": 1` (a JSON number, the natural hand-written form) is accepted.
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"{ "version": 1, "transcripts": [
            {"id": "NM_1.1", "strand": "+", "sequence": "ACGT",
             "exons": [{"number": 1, "start": 1, "end": 4}]}
        ] }"#;
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider = MockProvider::from_json(file.path()).expect("numeric version 1 loads");
        assert!(provider.has_transcript("NM_1.1"));
    }

    #[test]
    fn test_mock_provider_returns_sequence_length_for_added_sequence() {
        let mut provider = JsonProvider::new();
        provider.add_genomic_sequence("NC_012920.1", "A".repeat(16569));
        assert_eq!(provider.get_sequence_length("NC_012920.1").unwrap(), 16569);
    }

    #[test]
    fn test_mock_provider_sequence_length_errors_for_unknown_id() {
        let provider = JsonProvider::new();
        let err = provider.get_sequence_length("missing").unwrap_err();
        assert!(matches!(err, FerroError::ReferenceNotFound { .. }));
    }

    #[test]
    fn from_json_loads_convert_gff_output_with_transcript() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"{
            "version": "1.0",
            "genome_build": "GRCh38",
            "transcripts": [
                {
                    "id": "NM_000001.1",
                    "strand": "+",
                    "sequence": "ATGC",
                    "exons": [{"number": 1, "start": 1, "end": 4}]
                }
            ]
        }"#;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider = JsonProvider::from_json(file.path())
            .expect("convert-gff JSON with transcript should load");
        let tx = provider
            .get_transcript("NM_000001.1")
            .expect("tx not found");
        assert_eq!(tx.sequence.as_deref(), Some("ATGC"));
    }
}
