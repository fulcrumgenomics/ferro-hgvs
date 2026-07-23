//! Derive an `NG_`'s chromosomal [`GenomicPlacement`] from its GenBank record +
//! cdot, for versions absent from the RefSeqGene alignment snapshot (#728).
//!
//! `genomic_placement` normally matches an `NG_` **version-exact** against the
//! ingested RefSeqGene alignment GFF3, and **declines** when the pinned version
//! is absent (#655/#702) rather than mis-anchor across versions. The precomputed
//! alignment for old versions is gone from NCBI, leaving a coverage gap. This
//! module reconstructs the placement from data we *can* obtain by exact version:
//! the `NG_` GenBank record (which annotates each hosted transcript's exons in
//! `NG_` coordinates) plus cdot (which has the same transcript's exons in `NC_`
//! coordinates). Matching the two coordinate systems pins the single contiguous
//! `NG_→NC_` affine.
//!
//! **The invariant is decline-rather-than-mis-anchor.** Every step that cannot
//! be proven correct returns `None`. The decisive guard is a full base-by-base
//! comparison of the entire `NG_` sequence against the derived genome slice
//! ([`validate_sequence`]): a clean placement is ~100% identical, while any
//! indel between the `NG_` version and the assembly shifts the downstream frame
//! into a long mismatch run — so a gapped `NG_` (which cannot be one ungapped
//! affine anyway) is declined. Exon-coordinate matching only supplies the cheap
//! *candidate* offset; the sequence comparison is what makes it safe.
//!
//! All derivation runs at `ferro prepare` time; the derived placements are baked
//! into the manifest and served by the unchanged version-exact runtime path.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::hgvs::parser::accession::parse_accession;
use crate::hgvs::variant::Accession;
use crate::reference::provider::GenomicPlacement;
use crate::reference::Strand;
use crate::FerroError;

/// Minimum exon count required to trust a derived affine. A single exon
/// under-constrains it (one length match + one trivial offset is the weakest
/// possible evidence), so require at least two.
const MIN_EXONS: usize = 2;

/// Maximum fraction of mismatched bases tolerated by [`validate_sequence`].
/// Clean RefSeqGene-vs-assembly differences are isolated SNPs (well under 0.1%);
/// any indel produces a frame-shift mismatch run far above this, so the threshold
/// cleanly separates accept (clean) from decline (gapped/mis-anchored).
const MAX_MISMATCH_FRACTION: f64 = 0.005;

/// One transcript annotated in an `NG_` GenBank record: its accession and its
/// exon intervals in **`NG_` coordinates** (1-based inclusive, ascending).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HostedTranscript {
    /// The `/transcript_id` (e.g. `"NM_003002.4"`).
    pub transcript_id: String,
    /// Exon intervals `(start, end)` in `NG_` coordinates, ascending.
    pub exons: Vec<(u64, u64)>,
}

/// A candidate `NG_→NC_` affine derived from matched exon coordinates, before
/// sequence validation. `parent_start` is always 1 (an `NG_` is numbered from 1).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AffineCandidate {
    /// Chromosome coordinate of `NG_` position 1 (1-based, inclusive).
    pub nc_start: u64,
    /// Chromosome coordinate of `NG_` position `ng_len` (1-based, inclusive,
    /// `>= nc_start`).
    pub nc_end: u64,
    /// Orientation of the `NG_` relative to the chromosome.
    pub strand: Strand,
}

/// Parse the `mRNA` features of an `NG_` GenBank record into ordered
/// `(gene_upper, transcript_id)` pairs for the hosting map (#792).
///
/// Unlike [`parse_ng_hosted_transcripts`] (which skips partial `<`/`>` features
/// because their exon coordinates can't anchor an affine), this **includes**
/// partial features: a partial mRNA still authoritatively states "this `NG_`
/// hosts this transcript for this gene." Only features carrying BOTH a `/gene`
/// and a `/transcript_id` are kept. Order is the GenBank feature order.
pub fn parse_ng_gene_transcripts(genbank: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let bytes = genbank.as_bytes();
    let mut search_from = 0;
    while let Some(rel) = genbank[search_from..].find("     mRNA ") {
        let feat_start = search_from + rel;
        let block_end = next_feature_offset(bytes, feat_start + 1);
        let block = &genbank[feat_start..block_end];
        search_from = block_end;
        if let (Some(gene), Some(tx)) = (
            extract_qualifier(block, "gene"),
            extract_qualifier(block, "transcript_id"),
        ) {
            out.push((gene.to_ascii_uppercase(), tx));
        }
    }
    out
}

/// Parse the `mRNA` features of an `NG_` GenBank record into hosted transcripts.
///
/// Each `mRNA` feature carries an exon join in `NG_` coordinates and a
/// `/transcript_id`. A feature whose location contains a partial marker (`<` or
/// `>`) is **skipped**: the `NG_` holds only part of that transcript, so its
/// exon count/coordinates cannot be matched against cdot's full structure.
pub fn parse_ng_hosted_transcripts(genbank: &str) -> Vec<HostedTranscript> {
    let mut out = Vec::new();
    let bytes = genbank.as_bytes();
    let mut search_from = 0;
    while let Some(rel) = genbank[search_from..].find("     mRNA ") {
        let feat_start = search_from + rel;
        // The feature block runs until the next feature line (5-space indent +
        // a non-space) — capture the whole block so the join and the qualifiers
        // (which include /transcript_id) are together.
        let block_end = next_feature_offset(bytes, feat_start + 1);
        let block = &genbank[feat_start..block_end];
        search_from = block_end;

        if let Some(ht) = parse_mrna_block(block) {
            out.push(ht);
        }
    }
    out
}

/// Offset of the next GenBank feature line at or after `from`, i.e. a line that
/// begins with exactly five spaces followed by a non-space (a feature key), or
/// a line with a smaller indent (end of the FEATURES table). Returns the input
/// length if none is found.
fn next_feature_offset(bytes: &[u8], from: usize) -> usize {
    let text = match std::str::from_utf8(bytes) {
        Ok(t) => t,
        Err(_) => return bytes.len(),
    };
    let mut idx = from;
    while let Some(rel) = text[idx..].find('\n') {
        let line_start = idx + rel + 1;
        let line = &text[line_start..];
        // A feature key line: 5 spaces, then a non-space. A line indented fewer
        // than 5 spaces (e.g. "ORIGIN", "CONTIG") ends the table. Qualifier lines
        // (indent > 5) are part of the current feature, not a boundary.
        let indent = line.len() - line.trim_start_matches(' ').len();
        let rest = &line[indent..];
        if !rest.is_empty() && !rest.starts_with('\n') && indent <= 5 {
            return line_start;
        }
        idx = line_start;
        if line_start >= text.len() {
            break;
        }
    }
    bytes.len()
}

/// Parse a single `mRNA` feature block. Returns `None` if it is partial (`<`/`>`),
/// lacks a `/transcript_id`, or has no parseable exon intervals.
fn parse_mrna_block(block: &str) -> Option<HostedTranscript> {
    // The location spans from after "mRNA" to the first qualifier ("/..."). It
    // may wrap across lines; join and strip whitespace.
    let after_key = block.find("mRNA").map(|i| i + "mRNA".len())?;
    let qual_rel = block[after_key..].find("/")?;
    let location: String = block[after_key..after_key + qual_rel]
        .split_whitespace()
        .collect();
    if location.contains('<') || location.contains('>') {
        return None; // partial annotation — the NG holds only part of this tx
    }

    let transcript_id = extract_qualifier(block, "transcript_id")?;

    let mut exons: Vec<(u64, u64)> = Vec::new();
    // Scan for "<digits>..<digits>" pairs anywhere in the (possibly
    // complement(join(...))) location string.
    let loc = location.as_str();
    let lb = loc.as_bytes();
    let mut i = 0;
    while i < lb.len() {
        if lb[i].is_ascii_digit() {
            let start = i;
            while i < lb.len() && lb[i].is_ascii_digit() {
                i += 1;
            }
            let a: u64 = loc[start..i].parse().ok()?;
            // expect ".."
            if loc[i..].starts_with("..") {
                i += 2;
                let s2 = i;
                while i < lb.len() && lb[i].is_ascii_digit() {
                    i += 1;
                }
                if i > s2 {
                    let b: u64 = loc[s2..i].parse().ok()?;
                    exons.push((a.min(b), a.max(b)));
                }
            }
        } else {
            i += 1;
        }
    }
    if exons.is_empty() {
        return None;
    }
    exons.sort_unstable();
    Some(HostedTranscript {
        transcript_id,
        exons,
    })
}

/// Extract a single-line `/name="value"` qualifier value from a feature block.
fn extract_qualifier(block: &str, name: &str) -> Option<String> {
    let needle = format!("/{name}=\"");
    let start = block.find(&needle)? + needle.len();
    let end = block[start..].find('"')? + start;
    Some(block[start..end].to_string())
}

/// Derive the contiguous `NG_→NC_` affine from a hosted transcript's `NG_` exon
/// intervals and the same transcript's `NC_` exon intervals (from cdot, 1-based
/// inclusive). Returns `None` unless the two sets are mutually consistent:
///
/// - equal exon count and at least [`MIN_EXONS`];
/// - exon **lengths** match in exactly one orientation (ascending `NG_` against
///   ascending `NC_` ⇒ plus; against descending `NC_` ⇒ minus) — this both
///   detects strand and rejects structurally different transcripts;
/// - a **single** affine constant holds across *every* exon (any disagreement
///   means an internal indel ⇒ decline).
///
/// `ng_len` is the total `NG_` length (from its FASTA), used to place `nc_end`.
pub fn derive_affine(
    ng_exons: &[(u64, u64)],
    nc_exons: &[(u64, u64)],
    ng_len: u64,
) -> Option<AffineCandidate> {
    if ng_exons.len() != nc_exons.len() || ng_exons.len() < MIN_EXONS {
        return None;
    }
    let mut ng = ng_exons.to_vec();
    ng.sort_unstable();
    let mut nc_asc = nc_exons.to_vec();
    nc_asc.sort_unstable();

    let ng_lens: Vec<u64> = ng.iter().map(|(a, b)| b - a + 1).collect();
    let nc_asc_lens: Vec<u64> = nc_asc.iter().map(|(a, b)| b - a + 1).collect();
    let nc_desc: Vec<(u64, u64)> = nc_asc.iter().rev().copied().collect();
    let nc_desc_lens: Vec<u64> = nc_desc.iter().map(|(a, b)| b - a + 1).collect();

    if ng_lens == nc_asc_lens {
        // plus: NG ascending aligns with NC ascending.
        // NC_start_i = nc_start + (ng_start_i - 1)  =>  k = nc_start_i - ng_start_i constant.
        let mut k = None;
        for (&(ng_s, _), &(nc_s, _)) in ng.iter().zip(nc_asc.iter()) {
            let ki = nc_s.checked_sub(ng_s)?; // nc_s >= ng_s for a sane placement
            if *k.get_or_insert(ki) != ki {
                return None;
            }
        }
        let nc_start = k? + 1;
        Some(AffineCandidate {
            nc_start,
            nc_end: nc_start + ng_len - 1,
            strand: Strand::Plus,
        })
    } else if ng_lens == nc_desc_lens {
        // minus: NG ascending aligns with NC descending.
        // NG position p maps to NC = nc_end - (p - 1); pairing NG exon (ascending)
        // with NC exon (descending), the NG exon START (low p) maps to the NC
        // exon END (high coordinate): nc_end = nc_e_i + (ng_s_i - 1), constant.
        let mut e = None;
        for (&(ng_s, _), &(_, nc_e)) in ng.iter().zip(nc_desc.iter()) {
            let ei = nc_e + ng_s - 1;
            if *e.get_or_insert(ei) != ei {
                return None;
            }
        }
        let nc_end = e?;
        let nc_start = nc_end.checked_sub(ng_len - 1)?;
        Some(AffineCandidate {
            nc_start,
            nc_end,
            strand: Strand::Minus,
        })
    } else {
        None
    }
}

/// Validate a candidate affine by comparing the entire `NG_` sequence against the
/// genome slice it claims to map to. `genome_slice` must be the chromosome bases
/// for `[nc_start, nc_end]` (1-based inclusive), in chromosome `+` orientation.
///
/// Returns `true` only if the lengths are equal **and** the mismatch fraction is
/// below [`MAX_MISMATCH_FRACTION`]. For a minus-strand candidate the `NG_` is the
/// reverse complement of the slice. A length difference or a frame-shifting indel
/// drives the mismatch far above the threshold ⇒ `false` ⇒ decline. Comparison is
/// case-insensitive; `N`/ambiguity bases count as mismatches (conservative).
///
/// Note: ambiguity bases on **either** side count against the threshold, so a
/// `genome_slice` carrying a run of `N`s (e.g. an assembly gap) inflates the
/// mismatch fraction toward [`MAX_MISMATCH_FRACTION`] and can push an otherwise
/// valid placement to decline.
pub fn validate_sequence(ng_seq: &[u8], genome_slice: &[u8], strand: Strand) -> bool {
    if ng_seq.len() != genome_slice.len() || ng_seq.is_empty() {
        return false;
    }
    let mismatches = match strand {
        Strand::Plus => ng_seq
            .iter()
            .zip(genome_slice.iter())
            .filter(|(a, b)| !bases_match(**a, **b))
            .count(),
        Strand::Minus => {
            // NG[i] should equal complement(genome_slice[len-1-i]).
            ng_seq
                .iter()
                .zip(genome_slice.iter().rev())
                .filter(|(a, b)| !bases_match(**a, complement(**b)))
                .count()
        }
        Strand::Unknown => return false,
    };
    (mismatches as f64) / (ng_seq.len() as f64) < MAX_MISMATCH_FRACTION
}

/// Case-insensitive base equality (`a`==`A`).
fn bases_match(a: u8, b: u8) -> bool {
    a.eq_ignore_ascii_case(&b)
}

/// Watson-Crick complement (case-insensitive, returns uppercase). Non-ACGT maps
/// to `b'N'` so it can never spuriously match.
fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// Assemble a [`GenomicPlacement`] for `ng` on chromosome `nc` from a validated
/// candidate. (Construction helper so callers don't duplicate the field mapping;
/// `parent_start` is always 1.)
pub fn placement_from_candidate(nc: Accession, candidate: &AffineCandidate) -> GenomicPlacement {
    GenomicPlacement {
        nc,
        parent_start: 1,
        nc_start: candidate.nc_start,
        nc_end: candidate.nc_end,
        strand: candidate.strand,
    }
}

/// cdot's view of a transcript's chromosomal exons on one genome build. Supplied
/// by the orchestrator so [`derive_ng_placement`] stays decoupled from
/// `CdotMapper` (and unit-testable).
pub trait NcExonSource {
    /// The transcript's exon intervals on its chromosome (1-based inclusive),
    /// or `None` if cdot does not carry this **exact** transcript on this build.
    fn nc_exons(&self, transcript_id: &str) -> Option<Vec<(u64, u64)>>;
    /// The chromosome (`NC_`) accession the transcript maps to.
    fn nc_contig(&self, transcript_id: &str) -> Option<Accession>;
}

/// Genome base access for the sequence-validation guard (1-based inclusive,
/// chromosome `+` orientation). Supplied by the orchestrator.
pub trait GenomeSlice {
    fn slice(&self, nc: &Accession, start: u64, end: u64) -> Option<Vec<u8>>;
}

/// Derive a validated [`GenomicPlacement`] for an `NG_` from its GenBank record +
/// sequence, cdot, and the genome — or `None` (decline) if it cannot be proven.
///
/// Upholds the must-not-mis-anchor invariant via layered gates:
/// 1. parse hosted transcripts, skipping partial (`<`/`>`) annotations;
/// 2. for each transcript cdot carries on this build, derive a candidate affine
///    ([`derive_affine`]: ≥2 exons, exon-length-match orientation, single
///    consistent offset);
/// 3. **cross-transcript agreement** — every transcript that derives must yield
///    the *identical* (`nc`, candidate); any disagreement ⇒ decline;
/// 4. **full-sequence validation** — the entire `NG_` must match the derived
///    genome slice ([`validate_sequence`]); any indel/length change ⇒ decline.
pub fn derive_ng_placement(
    ng_genbank: &str,
    ng_seq: &[u8],
    nc_source: &impl NcExonSource,
    genome: &impl GenomeSlice,
) -> Option<GenomicPlacement> {
    let ng_len = ng_seq.len() as u64;
    if ng_len == 0 {
        return None;
    }

    let mut agreed: Option<(Accession, AffineCandidate)> = None;
    for ht in parse_ng_hosted_transcripts(ng_genbank) {
        let Some(nc_exons) = nc_source.nc_exons(&ht.transcript_id) else {
            continue; // cdot lacks this exact transcript on this build
        };
        let Some(nc) = nc_source.nc_contig(&ht.transcript_id) else {
            continue;
        };
        let Some(cand) = derive_affine(&ht.exons, &nc_exons, ng_len) else {
            continue; // this transcript can't anchor (count/length/offset)
        };
        match &agreed {
            None => agreed = Some((nc, cand)),
            Some((prev_nc, prev_cand)) => {
                // Two transcripts on one NG must place it identically.
                if *prev_nc != nc || *prev_cand != cand {
                    return None;
                }
            }
        }
    }

    let (nc, cand) = agreed?;
    let slice = genome.slice(&nc, cand.nc_start, cand.nc_end)?;
    if !validate_sequence(ng_seq, &slice, cand.strand) {
        return None;
    }
    Some(placement_from_candidate(nc, &cand))
}

// ----------------------------------------------------------------------------
// On-disk artifact (produced at prepare time, consumed at load)
// ----------------------------------------------------------------------------

/// Current on-disk schema version for the `derived_refseqgene_placements`
/// artifact. Bump when the shape changes incompatibly.
pub const DERIVED_PLACEMENTS_SCHEMA_VERSION: u32 = 1;

/// The committed/manifest artifact of derived `NG_`/`LRG_` placements: the
/// output of the prepare-time derivation, loaded by `MultiFastaProvider` and
/// merged into its `refseqgene_placements` map. Self-describing (strand as
/// `"+"`/`"-"`, per-entry provenance) so it is auditable independent of the
/// runtime [`GenomicPlacement`] (which is not serializable).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DerivedPlacements {
    /// On-disk schema version (#1001). Absent on artifacts written before
    /// versioning — those deserialize to `0` and are accepted as legacy. A
    /// value newer than [`DERIVED_PLACEMENTS_SCHEMA_VERSION`] is rejected at
    /// load, so a forward-incompatible artifact fails loudly instead of being
    /// silently misread.
    #[serde(default)]
    pub schema_version: u32,
    /// Human-facing provenance note (generator command, source manifest).
    #[serde(default)]
    pub description: String,
    /// One entry per derived placement.
    pub placements: Vec<DerivedPlacement>,
}

/// One derived placement entry, keyed to its versioned parent accession.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DerivedPlacement {
    /// Versioned `NG_`/`LRG_` accession this places (e.g. `"NG_012337.1"`).
    pub parent: String,
    /// Chromosome (`NC_`) accession (e.g. `"NC_000011.10"`).
    pub nc: String,
    /// Inclusive 1-based chromosome span.
    pub nc_start: u64,
    pub nc_end: u64,
    /// Parent orientation relative to the chromosome: `"+"` or `"-"`.
    pub strand: String,
    /// The transcript whose exons anchored the affine (provenance). Optional:
    /// empty when the producer does not surface it (the `dev` example producer
    /// currently records `""`); not consumed at load.
    pub anchored_by: String,
    /// The full-sequence validation mismatch fraction (provenance; ~0 for a
    /// clean placement). Optional: the `dev` example producer derives placements
    /// only after validation passes and records `0.0` rather than the exact
    /// fraction; not consumed at load.
    pub mismatch_fraction: f64,
}

impl DerivedPlacements {
    /// Load from a JSON file.
    pub fn from_json_path<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        let content = std::fs::read_to_string(path.as_ref())?;
        let parsed: Self = serde_json::from_str(&content)?;
        if parsed.schema_version > DERIVED_PLACEMENTS_SCHEMA_VERSION {
            return Err(FerroError::Io {
                msg: format!(
                    "derived_refseqgene_placements artifact at {} has schema_version {}, \
                     which is newer than this build supports (maximum {}). Upgrade ferro, or \
                     re-derive the artifact.",
                    path.as_ref().display(),
                    parsed.schema_version,
                    DERIVED_PLACEMENTS_SCHEMA_VERSION
                ),
            });
        }
        Ok(parsed)
    }

    /// Serialize to pretty JSON with a trailing newline (stable for `--check`).
    pub fn to_json(&self) -> Result<String, FerroError> {
        let mut s = serde_json::to_string_pretty(self)?;
        s.push('\n');
        Ok(s)
    }

    /// Convert to `(parent_key, GenomicPlacement)` entries for merging into the
    /// `refseqgene_placements` map. An entry with an unparseable `nc` accession
    /// or a strand other than `"+"`/`"-"` is **skipped** (never silently
    /// mis-placed). `parent_start` is always 1.
    pub fn to_placements(&self) -> Vec<(String, GenomicPlacement)> {
        self.placements
            .iter()
            .filter_map(|p| {
                let strand = match p.strand.as_str() {
                    "+" => Strand::Plus,
                    "-" => Strand::Minus,
                    _ => return None,
                };
                if p.nc_start > p.nc_end {
                    return None;
                }
                // Decline on an unparseable / trailing-garbage NC accession
                // rather than place against a malformed one.
                let nc = match parse_accession(&p.nc) {
                    Ok(("", acc)) => acc, // fully consumed — no trailing garbage
                    _ => return None,
                };
                Some((
                    p.parent.clone(),
                    GenomicPlacement {
                        nc,
                        parent_start: 1,
                        nc_start: p.nc_start,
                        nc_end: p.nc_end,
                        strand,
                    },
                ))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legacy_artifact_without_schema_version_is_accepted() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("a.json");
        // No `schema_version` key at all — the pre-versioning on-disk shape.
        std::fs::write(&p, br#"{"description":"legacy","placements":[]}"#).unwrap();
        let loaded = DerivedPlacements::from_json_path(&p).expect("legacy artifact must load");
        assert_eq!(loaded.schema_version, 0, "absent version reads as 0");
    }

    #[test]
    fn current_schema_version_is_accepted() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("a.json");
        std::fs::write(
            &p,
            format!(
                r#"{{"schema_version":{},"description":"","placements":[]}}"#,
                DERIVED_PLACEMENTS_SCHEMA_VERSION
            ),
        )
        .unwrap();
        assert!(DerivedPlacements::from_json_path(&p).is_ok());
    }

    #[test]
    fn newer_schema_version_is_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("a.json");
        std::fs::write(
            &p,
            format!(
                r#"{{"schema_version":{},"description":"","placements":[]}}"#,
                DERIVED_PLACEMENTS_SCHEMA_VERSION + 1
            ),
        )
        .unwrap();
        let err = DerivedPlacements::from_json_path(&p).unwrap_err();
        assert!(
            format!("{err}").contains("newer than this build"),
            "expected an actionable version error, got: {err}"
        );
    }

    // ---- parse_ng_gene_transcripts ----

    const NG_GENE_SAMPLE: &str = concat!(
        "     mRNA            join(1..100,200..300)\n",
        "                     /gene=\"TIMM8B\"\n",
        "                     /transcript_id=\"NM_012459.2\"\n",
        "     mRNA            join(<136..200,400..500)\n",
        "                     /gene=\"C11orf57\"\n",
        "                     /transcript_id=\"NM_001321072.1\"\n",
        "     CDS             join(50..100,200..280)\n",
        "                     /gene=\"TIMM8B\"\n",
    );

    #[test]
    fn parse_ng_gene_transcripts_captures_gene_and_includes_partials() {
        let got = parse_ng_gene_transcripts(NG_GENE_SAMPLE);
        assert_eq!(
            got,
            vec![
                ("TIMM8B".to_string(), "NM_012459.2".to_string()),
                // partial (`<136`) feature is INCLUDED for hosting membership:
                ("C11ORF57".to_string(), "NM_001321072.1".to_string()),
            ]
        );
    }

    // ---- GenBank mRNA parsing ----

    // GenBank feature keys are indented exactly 5 spaces, qualifiers more — use
    // `concat!` with explicit "\n" so the leading spaces are preserved (a `"\`
    // line-continuation would strip them).
    const SDHD_MRNA: &str = concat!(
        "     mRNA            join(5027..5113,6011..6127,7021..7165,12959..13948)\n",
        "                     /gene=\"SDHD\"\n",
        "                     /product=\"succinate dehydrogenase complex subunit D\"\n",
        "                     /transcript_id=\"NM_003002.4\"\n",
        "     CDS             join(5089..5113,6011..6127,7021..7165,12959..13251)\n",
        "                     /gene=\"SDHD\"\n",
        "                     /protein_id=\"NP_002993.1\"\n",
    );

    const PARTIAL_MRNA: &str = concat!(
        "     mRNA            join(<136..189,619..3304)\n",
        "                     /gene=\"C11orf57\"\n",
        "                     /transcript_id=\"NM_018195.3\"\n",
        "     CDS             join(700..3304)\n",
        "                     /protein_id=\"NP_060665.3\"\n",
    );

    const COMPLEMENT_MRNA: &str = concat!(
        "     mRNA            complement(join(355610..359344,359722..359915,381412..381574))\n",
        "                     /gene=\"PCCB\"\n",
        "                     /transcript_id=\"NM_022153.2\"\n",
        "     CDS             complement(join(355700..359344))\n",
        "                     /protein_id=\"NP_071436.1\"\n",
    );

    #[test]
    fn parses_plus_mrna_exons() {
        let hts = parse_ng_hosted_transcripts(SDHD_MRNA);
        assert_eq!(hts.len(), 1);
        assert_eq!(hts[0].transcript_id, "NM_003002.4");
        assert_eq!(
            hts[0].exons,
            vec![(5027, 5113), (6011, 6127), (7021, 7165), (12959, 13948)]
        );
    }

    #[test]
    fn skips_partial_mrna() {
        // `<136` marks a 5'-truncated annotation: the NG holds only part of the
        // transcript, so it must not be used to anchor.
        assert!(parse_ng_hosted_transcripts(PARTIAL_MRNA).is_empty());
    }

    #[test]
    fn parses_complement_join_exons_ascending() {
        let hts = parse_ng_hosted_transcripts(COMPLEMENT_MRNA);
        assert_eq!(hts.len(), 1);
        assert_eq!(hts[0].transcript_id, "NM_022153.2");
        // Exons are returned in ascending NG order regardless of complement().
        assert_eq!(
            hts[0].exons,
            vec![(355610, 359344), (359722, 359915), (381412, 381574)]
        );
    }

    #[test]
    fn excludes_following_cds_join_from_mrna_exons() {
        // A CDS feature immediately follows the mRNA, with its own join(...) whose
        // coordinate ranges differ from the mRNA's. `next_feature_offset` must bound
        // the mRNA block at the CDS key line so the CDS ranges are never parsed as
        // extra exons. The CDS uses coordinates (90000..90100, 91000..91100) that
        // do not appear among the mRNA exons, so a regression that folded them in
        // would surface here.
        let record = concat!(
            "     mRNA            join(5027..5113,6011..6127,7021..7165,12959..13948)\n",
            "                     /gene=\"SDHD\"\n",
            "                     /transcript_id=\"NM_003002.4\"\n",
            "     CDS             join(90000..90100,91000..91100)\n",
            "                     /gene=\"SDHD\"\n",
            "                     /protein_id=\"NP_002993.1\"\n",
        );
        let hts = parse_ng_hosted_transcripts(record);
        assert_eq!(hts.len(), 1);
        // Only the mRNA's four exons — the CDS join ranges must be excluded.
        assert_eq!(
            hts[0].exons,
            vec![(5027, 5113), (6011, 6127), (7021, 7165), (12959, 13948)]
        );
    }

    // ---- affine derivation (real ground-truth coordinates) ----

    #[test]
    fn derives_plus_affine_matching_ground_truth() {
        // NG_012337.3 / NM_003002.4 (SDHD, + on NC). Ground truth from the
        // RefSeqGene GFF3: NC_000011.10, NG pos 1 -> NC 112081847, len 39784.
        let ng = [(5027, 5113), (6011, 6127), (7021, 7165), (12959, 13948)];
        let nc = [
            (112086873, 112086959),
            (112087857, 112087973),
            (112088867, 112089011),
            (112094805, 112095794),
        ];
        let cand = derive_affine(&ng, &nc, 39784).expect("derives");
        assert_eq!(cand.strand, Strand::Plus);
        assert_eq!(cand.nc_start, 112081847);
        assert_eq!(cand.nc_end, 112081847 + 39784 - 1);
    }

    #[test]
    fn derives_minus_affine_matching_ground_truth() {
        // NG_033963.1 / NM_002839.4 (PTPRD, - on NC), first 3 exons. Ground
        // truth: NC_000009.12, minus, nc_end 10617723, len 2305478.
        // cdot NC exons (1-based) descending-paired with NG ascending.
        let ng = [(1, 100), (5000, 5100), (20000, 20050)];
        // Construct NC so that nc_end = 10617723: for NG start s, NC end = 10617723 - (s-1).
        // exon (1,100): NC end = 10617723, NC start = 10617723-99 = 10617624.
        // exon (5000,5100): NC end = 10617723-4999=10612724, start=10612724-100=10612624.
        // exon (20000,20050): NC end = 10617723-19999=10597724, start=10597724-50=10597674.
        let nc = [
            (10617624, 10617723),
            (10612624, 10612724),
            (10597674, 10597724),
        ];
        let cand = derive_affine(&ng, &nc, 2305478).expect("derives");
        assert_eq!(cand.strand, Strand::Minus);
        assert_eq!(cand.nc_end, 10617723);
        assert_eq!(cand.nc_start, 10617723 - 2305478 + 1);
    }

    #[test]
    fn declines_when_offset_inconsistent() {
        // An internal indel: exon 2 shifted by +1 relative to a clean affine.
        let ng = [(100, 200), (300, 400), (500, 600)];
        let nc = [(1100, 1200), (1301, 1401), (1500, 1600)];
        assert!(derive_affine(&ng, &nc, 1000).is_none());
    }

    #[test]
    fn declines_on_exon_count_mismatch() {
        let ng = [(100, 200), (300, 400)];
        let nc = [(1100, 1200), (1300, 1400), (1500, 1600)];
        assert!(derive_affine(&ng, &nc, 1000).is_none());
    }

    #[test]
    fn declines_single_exon_underconstrained() {
        let ng = [(100, 200)];
        let nc = [(1100, 1200)];
        assert!(derive_affine(&ng, &nc, 1000).is_none());
    }

    #[test]
    fn declines_on_length_mismatch_neither_orientation() {
        // Same count, but exon lengths differ in both orientations.
        let ng = [(100, 200), (300, 450)]; // lens 101, 151
        let nc = [(1100, 1200), (1300, 1400)]; // lens 101, 101
        assert!(derive_affine(&ng, &nc, 1000).is_none());
    }

    // ---- sequence validation ----

    #[test]
    fn validates_clean_plus() {
        let ng = b"ACGTACGT";
        assert!(validate_sequence(ng, b"ACGTACGT", Strand::Plus));
        assert!(validate_sequence(ng, b"acgtacgt", Strand::Plus)); // case-insensitive
    }

    #[test]
    fn validates_clean_minus_revcomp() {
        // NG == revcomp(genome_slice).
        let genome = b"ACGTACGT";
        let ng = b"ACGTACGT"; // revcomp of ACGTACGT is ACGTACGT
        assert!(validate_sequence(ng, genome, Strand::Minus));
        let genome2 = b"AAAACCCC";
        let ng2 = b"GGGGTTTT"; // revcomp(AAAACCCC)
        assert!(validate_sequence(ng2, genome2, Strand::Minus));
    }

    #[test]
    fn declines_frame_shifted_sequence() {
        // A 1bp indel shifts the frame: ~half mismatch -> well above threshold.
        let ng = b"ACGTACGTACGTACGTACGT";
        let shifted = b"CGTACGTACGTACGTACGTA"; // rotated by 1
        assert!(!validate_sequence(ng, shifted, Strand::Plus));
    }

    #[test]
    fn declines_length_mismatch() {
        assert!(!validate_sequence(b"ACGT", b"ACGTA", Strand::Plus));
        assert!(!validate_sequence(b"", b"", Strand::Plus));
    }

    #[test]
    fn tolerates_isolated_snps_below_threshold() {
        // One SNP in 1000 bases is well under 0.5% -> still accepts.
        let mut ng = vec![b'A'; 1000];
        let genome = vec![b'A'; 1000];
        ng[500] = b'C'; // single SNP
        assert!(validate_sequence(&ng, &genome, Strand::Plus));
    }

    // ---- end-to-end orchestration (derive_ng_placement) ----

    fn acc(prefix: &str, number: &str, version: u32) -> Accession {
        Accession::new(prefix, number, Some(version))
    }

    /// Build a minimal NG GenBank fragment with one mRNA (+ a bounding CDS) for
    /// `tid` whose exons are `exons` (NG coords).
    fn ng_record(tid: &str, exons: &[(u64, u64)]) -> String {
        let join = exons
            .iter()
            .map(|(a, b)| format!("{a}..{b}"))
            .collect::<Vec<_>>()
            .join(",");
        // Explicit "\n" + literal indentation (5 spaces for feature keys, 21 for
        // qualifiers); a `\`-continuation would strip the leading spaces.
        let q = "                     "; // 21 spaces
        format!(
            "     mRNA            join({join})\n{q}/transcript_id=\"{tid}\"\n     CDS             join({join})\n{q}/protein_id=\"NP_000000.1\"\n"
        )
    }

    struct MockNc {
        tid: &'static str,
        exons: Vec<(u64, u64)>,
        nc: Accession,
    }
    impl NcExonSource for MockNc {
        fn nc_exons(&self, t: &str) -> Option<Vec<(u64, u64)>> {
            (t == self.tid).then(|| self.exons.clone())
        }
        fn nc_contig(&self, t: &str) -> Option<Accession> {
            (t == self.tid).then(|| self.nc.clone())
        }
    }
    struct MockGenome {
        nc: Accession,
        start: u64,
        bases: Vec<u8>,
    }
    impl GenomeSlice for MockGenome {
        fn slice(&self, nc: &Accession, start: u64, end: u64) -> Option<Vec<u8>> {
            (nc == &self.nc
                && start == self.start
                && end == self.start + self.bases.len() as u64 - 1)
                .then(|| self.bases.clone())
        }
    }

    #[test]
    fn derive_ng_placement_plus_success() {
        let ng_seq = b"ACGTACGTACGTACGTACGT"; // len 20
                                              // NG exons (3,6),(11,14); NG pos 1 -> NC 1001 => NC exons (1003,1006),(1011,1014).
        let gb = ng_record("NM_X.1", &[(3, 6), (11, 14)]);
        let nc = acc("NC", "000001", 11);
        let src = MockNc {
            tid: "NM_X.1",
            exons: vec![(1003, 1006), (1011, 1014)],
            nc: nc.clone(),
        };
        let genome = MockGenome {
            nc: nc.clone(),
            start: 1001,
            bases: ng_seq.to_vec(),
        };
        let p = derive_ng_placement(&gb, ng_seq, &src, &genome).expect("derives");
        assert_eq!(p.nc, nc);
        assert_eq!(p.parent_start, 1);
        assert_eq!(p.nc_start, 1001);
        assert_eq!(p.nc_end, 1020);
        assert_eq!(p.strand, Strand::Plus);
    }

    #[test]
    fn derive_ng_placement_minus_success() {
        // Minus-strand end-to-end: NG ascending exons align with NC *descending*
        // exons, and the NG sequence is the reverse complement of the genome slice.
        //
        // NG exons (3,6) len 4 and (11,16) len 6 -> ng_lens [4, 6].
        // NC exons (1001,1006) len 6 and (1011,1014) len 4 -> nc_asc_lens [6, 4];
        // reversed nc_desc_lens [4, 6] == ng_lens, so only the minus orientation
        // matches (the plus branch declines on the [4,6] != [6,4] length check).
        // The minus affine then yields nc_end = nc1_e + (ng_s - 1) = 1006 + 10 = 1016
        // and nc_start = nc_end - (ng_len - 1) = 1016 - 19 = 997.
        let ng_seq = b"AAAACCCCGGGGTTTTACGT"; // len 20, not revcomp-palindromic
        let gb = ng_record("NM_X.1", &[(3, 6), (11, 16)]);
        let nc = acc("NC", "000001", 11);
        let src = MockNc {
            tid: "NM_X.1",
            exons: vec![(1001, 1006), (1011, 1014)],
            nc: nc.clone(),
        };
        // Genome slice [997, 1016] is the reverse complement of ng_seq.
        let genome = MockGenome {
            nc: nc.clone(),
            start: 997,
            bases: b"ACGTAAAACCCCGGGGTTTT".to_vec(),
        };
        let p = derive_ng_placement(gb.as_str(), ng_seq, &src, &genome).expect("derives");
        assert_eq!(p.nc, nc);
        assert_eq!(p.parent_start, 1);
        assert_eq!(p.nc_start, 997);
        assert_eq!(p.nc_end, 1016);
        assert_eq!(p.strand, Strand::Minus);
    }

    #[test]
    fn derive_ng_placement_declines_on_sequence_mismatch() {
        // Genome slice is frame-shifted relative to the NG -> validation fails.
        let ng_seq = b"ACGTACGTACGTACGTACGT";
        let gb = ng_record("NM_X.1", &[(3, 6), (11, 14)]);
        let nc = acc("NC", "000001", 11);
        let src = MockNc {
            tid: "NM_X.1",
            exons: vec![(1003, 1006), (1011, 1014)],
            nc: nc.clone(),
        };
        let genome = MockGenome {
            nc: nc.clone(),
            start: 1001,
            bases: b"CGTACGTACGTACGTACGTA".to_vec(), // rotated by 1
        };
        assert!(derive_ng_placement(&gb, ng_seq, &src, &genome).is_none());
    }

    #[test]
    fn derive_ng_placement_declines_when_genome_absent() {
        struct NoGenome;
        impl GenomeSlice for NoGenome {
            fn slice(&self, _: &Accession, _: u64, _: u64) -> Option<Vec<u8>> {
                None
            }
        }
        let ng_seq = b"ACGTACGTACGTACGTACGT";
        let gb = ng_record("NM_X.1", &[(3, 6), (11, 14)]);
        let nc = acc("NC", "000001", 11);
        let src = MockNc {
            tid: "NM_X.1",
            exons: vec![(1003, 1006), (1011, 1014)],
            nc,
        };
        assert!(derive_ng_placement(&gb, ng_seq, &src, &NoGenome).is_none());
    }

    #[test]
    fn derive_ng_placement_declines_on_cross_transcript_disagreement() {
        // Two transcripts on the NG that derive different placements -> decline.
        struct TwoNc {
            nc: Accession,
        }
        impl NcExonSource for TwoNc {
            fn nc_exons(&self, t: &str) -> Option<Vec<(u64, u64)>> {
                match t {
                    "NM_A.1" => Some(vec![(1003, 1006), (1011, 1014)]), // nc_start 1001
                    "NM_B.1" => Some(vec![(2003, 2006), (2011, 2014)]), // nc_start 2001
                    _ => None,
                }
            }
            fn nc_contig(&self, t: &str) -> Option<Accession> {
                t.starts_with("NM_").then(|| self.nc.clone())
            }
        }
        let ng_seq = b"ACGTACGTACGTACGTACGT";
        let gb = format!(
            "{}{}",
            ng_record("NM_A.1", &[(3, 6), (11, 14)]),
            ng_record("NM_B.1", &[(3, 6), (11, 14)])
        );
        let nc = acc("NC", "000001", 11);
        let genome = MockGenome {
            nc: nc.clone(),
            start: 1001,
            bases: ng_seq.to_vec(),
        };
        assert!(derive_ng_placement(&gb, ng_seq, &TwoNc { nc }, &genome).is_none());
    }

    // ---- DerivedPlacements artifact ----

    fn sample_entry() -> DerivedPlacement {
        DerivedPlacement {
            parent: "NG_012337.1".to_string(),
            nc: "NC_000011.10".to_string(),
            nc_start: 112081847,
            nc_end: 112097794,
            strand: "+".to_string(),
            anchored_by: "NM_003002.2".to_string(),
            mismatch_fraction: 0.0,
        }
    }

    #[test]
    fn derived_placements_json_round_trips() {
        let dp = DerivedPlacements {
            schema_version: DERIVED_PLACEMENTS_SCHEMA_VERSION,
            description: "test".to_string(),
            placements: vec![sample_entry()],
        };
        let json = dp.to_json().expect("serializes");
        assert!(json.ends_with('\n'));
        let parsed: DerivedPlacements = serde_json::from_str(&json).expect("round-trips");
        assert_eq!(parsed, dp);
    }

    #[test]
    fn derived_placements_rejects_unknown_field() {
        let json = r#"{"placements":[],"bogus":1}"#;
        assert!(serde_json::from_str::<DerivedPlacements>(json).is_err());
    }

    #[test]
    fn to_placements_builds_genomic_placement() {
        let dp = DerivedPlacements {
            description: String::new(),
            placements: vec![sample_entry()],
            ..Default::default()
        };
        let out = dp.to_placements();
        assert_eq!(out.len(), 1);
        let (key, p) = &out[0];
        assert_eq!(key, "NG_012337.1");
        assert_eq!(p.nc.full(), "NC_000011.10");
        assert_eq!(p.parent_start, 1);
        assert_eq!(p.nc_start, 112081847);
        assert_eq!(p.nc_end, 112097794);
        assert_eq!(p.strand, Strand::Plus);
    }

    #[test]
    fn to_placements_minus_strand() {
        let mut e = sample_entry();
        e.strand = "-".to_string();
        let dp = DerivedPlacements {
            description: String::new(),
            placements: vec![e],
            ..Default::default()
        };
        assert_eq!(dp.to_placements()[0].1.strand, Strand::Minus);
    }

    #[test]
    fn to_placements_skips_bad_strand_and_inverted_range() {
        let mut bad_strand = sample_entry();
        bad_strand.strand = "?".to_string();
        let mut inverted = sample_entry();
        inverted.nc_start = 200;
        inverted.nc_end = 100;
        let dp = DerivedPlacements {
            description: String::new(),
            placements: vec![bad_strand, inverted],
            ..Default::default()
        };
        // Both invalid entries are skipped rather than mis-placed.
        assert!(dp.to_placements().is_empty());
    }
}
