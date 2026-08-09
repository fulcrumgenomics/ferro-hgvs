//! Synthetic `MockProvider` builders for nucleotide-shift coverage matrices.
//!
//! See `docs/superpowers/specs/2026-05-01-A7-ins-shift-coverage-matrix-design.md`.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, Normalizer};

/// 256 bases of padding to ensure the normalizer's 100bp window stays in
/// bounds on either side of any test variant.
const PAD: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
                   ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
                   ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
                   ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

/// The 0-based offset added by padding — used by cds/noncoding/rna builders.
pub const PAD_OFFSET: u64 = 256;
const GENOMIC_CONTIG: &str = "NC_TEST.1";
const TX_ACCESSION: &str = "NM_TEST.1";
const NR_ACCESSION: &str = "NR_TEST.1";
const TX_CONTIG: &str = "chr_synth";

/// Wrap `core` in [`PAD_OFFSET`] bases of `ACGT…` padding on each side.
///
/// Used internally by every builder to keep the normalizer's 100bp window in
/// bounds. It is also public so a caller can build a **transcript** out of the
/// padded string and thereby give that transcript a real 5'UTR: passing
/// `padded(core)` to [`SyntheticBuilder::cds`] with `cds_start = PAD_OFFSET + 1`
/// puts core base `p` at CDS position `p`, while the CDS translation itself is
/// non-trivial (offset by [`PAD_OFFSET`]) rather than the identity you get from
/// `cds_start == 1`.
///
/// The padding is a perfect period-4 `ACGT` tandem: the base immediately 5' of
/// the core is `T` and the base immediately 3' of it is `A`. A core whose first
/// base is not `A` and whose last base is not `T` therefore cannot extend the
/// pad's own rotation, which is what bounds a repeat tract to the core.
pub fn padded(core: &str) -> String {
    format!("{}{}{}", PAD, core, PAD)
}

/// Builder for synthetic `MockProvider` fixtures across HGVS coord systems.
///
/// All builders pad the input core sequence with 256 bases on each side so
/// the normalizer's 100bp lookahead window stays in range. The first 1-based
/// HGVS position of the core region is `PAD_OFFSET + 1` (= 257).
pub struct SyntheticBuilder {
    inner: CoordSystem,
}

// Per-binary dead-code allowance: not every test binary that uses this
// helper exercises every coord system (e.g., the B1 matrix has no r.
// modules because r. is not in the spec's accepted reference types for
// repeats), and clippy --tests evaluates each binary independently.
#[allow(dead_code)]
enum CoordSystem {
    Genomic {
        core: String,
    },
    Cds {
        core: String,
        cds_start: u64,
        cds_end: u64,
        strand: Strand,
    },
    Noncoding {
        core: String,
        strand: Strand,
    },
    Rna {
        core: String,
        strand: Strand,
    },
}

impl SyntheticBuilder {
    /// Build a genomic provider over the contig `NC_TEST.1` with the given
    /// core sequence. The core's first base is at 1-based HGVS position 257.
    pub fn genomic(core: &str) -> Self {
        Self {
            inner: CoordSystem::Genomic {
                core: core.to_string(),
            },
        }
    }

    /// Build a CDS provider with a transcript named `NM_TEST.1`.
    ///
    /// `core` is the transcript sequence (no padding needed — the transcript
    /// is treated as the full reference for c. coordinates). `cds_start` and
    /// `cds_end` are 1-based inclusive transcript positions delimiting the
    /// CDS region. The transcript is mapped to a genomic contig `chr_synth`
    /// at offset `PAD_OFFSET` for plus strand, or reverse-complemented at
    /// the same offset for minus strand.
    pub fn cds(core: &str, cds_start: u64, cds_end: u64, strand: Strand) -> Self {
        Self {
            inner: CoordSystem::Cds {
                core: core.to_string(),
                cds_start,
                cds_end,
                strand,
            },
        }
    }

    /// Build a non-coding RNA provider with a transcript named `NR_TEST.1`.
    pub fn noncoding(core: &str, strand: Strand) -> Self {
        Self {
            inner: CoordSystem::Noncoding {
                core: core.to_string(),
                strand,
            },
        }
    }

    /// Build an RNA provider with a transcript named `NR_TEST.1`. The
    /// transcript sequence should be lowercase RNA bases (a/c/g/u). The
    /// underlying genomic mapping is built from the DNA-equivalent.
    #[allow(dead_code)] // not used by every test binary; see CoordSystem allow above
    pub fn rna(core: &str, strand: Strand) -> Self {
        Self {
            inner: CoordSystem::Rna {
                core: core.to_string(),
                strand,
            },
        }
    }

    /// Materialize the configured `MockProvider` ready for use in a test.
    pub fn build(self) -> MockProvider {
        let mut provider = MockProvider::new();
        match self.inner {
            CoordSystem::Genomic { core } => {
                provider.add_genomic_sequence(GENOMIC_CONTIG, padded(&core));
            }
            CoordSystem::Cds {
                core,
                cds_start,
                cds_end,
                strand,
            } => {
                let tx_len = core.len() as u64;
                let genomic = padded(&match strand {
                    Strand::Plus => core.clone(),
                    Strand::Minus => reverse_complement(&core),
                    Strand::Unknown => {
                        unreachable!("synthetic transcripts never have unknown strand")
                    }
                    // #[non_exhaustive] wildcard — no other variants exist today.
                    _ => unreachable!("unexpected Strand variant"),
                });
                let g_start = PAD_OFFSET + 1;
                let g_end = PAD_OFFSET + tx_len;
                let exon = Exon::with_genomic(1, 1, tx_len, g_start, g_end);
                let transcript = Transcript::new(
                    TX_ACCESSION.to_string(),
                    Some("SYNTH".to_string()),
                    strand,
                    core,
                    Some(cds_start),
                    Some(cds_end),
                    vec![exon],
                    Some(TX_CONTIG.to_string()),
                    Some(g_start),
                    Some(g_end),
                    GenomeBuild::GRCh38,
                    ManeStatus::None,
                    None,
                    None,
                );
                provider.add_genomic_sequence(TX_CONTIG, genomic);
                provider.add_transcript(transcript);
            }
            CoordSystem::Noncoding { core, strand } => {
                let tx_len = core.len() as u64;
                let genomic = padded(&match strand {
                    Strand::Minus => reverse_complement(&core),
                    _ => core.clone(),
                });
                let g_start = PAD_OFFSET + 1;
                let g_end = PAD_OFFSET + tx_len;
                let exon = Exon::with_genomic(1, 1, tx_len, g_start, g_end);
                let transcript = Transcript::new(
                    NR_ACCESSION.to_string(),
                    Some("SYNTH_NR".to_string()),
                    strand,
                    core,
                    None,
                    None,
                    vec![exon],
                    Some(TX_CONTIG.to_string()),
                    Some(g_start),
                    Some(g_end),
                    GenomeBuild::GRCh38,
                    ManeStatus::None,
                    None,
                    None,
                );
                provider.add_genomic_sequence(TX_CONTIG, genomic);
                provider.add_transcript(transcript);
            }
            CoordSystem::Rna { core, strand } => {
                let tx_len = core.len() as u64;
                let dna_core = rna_to_dna(&core);
                let genomic = padded(&match strand {
                    Strand::Minus => reverse_complement(&dna_core),
                    _ => dna_core.clone(),
                });
                let g_start = PAD_OFFSET + 1;
                let g_end = PAD_OFFSET + tx_len;
                let exon = Exon::with_genomic(1, 1, tx_len, g_start, g_end);
                let transcript = Transcript::new(
                    NR_ACCESSION.to_string(),
                    Some("SYNTH_RNA".to_string()),
                    strand,
                    dna_core,
                    None,
                    None,
                    vec![exon],
                    Some(TX_CONTIG.to_string()),
                    Some(g_start),
                    Some(g_end),
                    GenomeBuild::GRCh38,
                    ManeStatus::None,
                    None,
                    None,
                );
                provider.add_genomic_sequence(TX_CONTIG, genomic);
                provider.add_transcript(transcript);
            }
        }
        provider
    }
}

/// Builder for `MockProvider` fixtures that deliberately put a transcript
/// id and a chromosome name into the same key space — the failure mode
/// behind #311 / #312. The fixture mirrors the issue's reproducer:
///
/// - one transcript with caller-supplied `id` whose `chromosome` field
///   is the contig name (no separate `genomic_sequences` entry).
/// - the transcript sequence is laid out so that the genomic-coordinate
///   bases at `g.1000..=1002` (= `GAC`) are observably different from
///   the bases a wrong-frame lookup would return (`AAG`); a wrong-frame
///   trim therefore collapses an `ACG` delins to a single residual
///   `1001A>C`. The numerical positions match the issue's reproducer
///   exactly so downstream cross-PR tests can share assertions.
///
/// Sequence layout, all 1102 bp:
///   - 909 bp poly-A leader
///   - 93 bp coding segment (`SEGMENT` below)
///   - 100 bp poly-T trailer
///
/// With `genomic_start = 91`, genomic position 1000..=1002 maps to
/// transcript positions 910..=912 (= `SEGMENT[0..3] = "GAC"`). The
/// wrong-frame read at transcript indices 999..=1001 returns
/// `SEGMENT[90..93] = "AAG"`, which is the buggy reference span that
/// drove #311's collapse to `1001A>C`.
pub struct AdversarialNamespaceBuilder {
    transcript_id: String,
    chromosome: String,
}

const ADVERSARIAL_SEGMENT: &str =
    "GACCGGCGGTCTGCAAGTCTCCACCTTCCGAAGCTATCCATAACCGGAACCTATGACCTAAAATCAGTTCTGGGTCAGCTCGGTATCACGAAG";

impl AdversarialNamespaceBuilder {
    pub fn new(transcript_id: &str, chromosome: &str) -> Self {
        Self {
            transcript_id: transcript_id.to_string(),
            chromosome: chromosome.to_string(),
        }
    }

    pub fn build(self) -> MockProvider {
        let mut s = String::with_capacity(1102);
        s.extend(std::iter::repeat_n('A', 909));
        s.push_str(ADVERSARIAL_SEGMENT);
        s.extend(std::iter::repeat_n('T', 100));
        debug_assert_eq!(s.len(), 1102);
        debug_assert_eq!(&s[909..912], "GAC");

        let exons = vec![Exon::with_genomic(1, 1, 1102, 91, 1192)];
        let transcript = Transcript::new(
            self.transcript_id,
            None,
            Strand::Plus,
            s,
            Some(1),
            Some(1102),
            exons,
            Some(self.chromosome),
            Some(91),
            Some(1192),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );

        let mut provider = MockProvider::new();
        provider.add_transcript(transcript);
        provider
    }
}

/// Substitute 1-based numeric positions into an HGVS template by index.
///
/// Replaces `{0}`, `{1}`, ... in `template` with `args[0]`, `args[1]`, ...
/// Used by the c./n./r. matrix tests (cds, noncoding, rna) where the
/// padding offset does NOT shift positions — the transcript itself is the
/// reference, so positions are literal. The genomic matrix uses a
/// per-module helper that adds `PAD_OFFSET` to each position.
#[allow(dead_code)] // not used by every test binary; see CoordSystem allow above
pub fn hgvs(template: &str, args: &[u64]) -> String {
    let mut s = template.to_string();
    for (i, p) in args.iter().enumerate() {
        s = s.replace(&format!("{{{}}}", i), &p.to_string());
    }
    s
}

/// Normalize an HGVS variant string against a provider and return the
/// formatted output. Panics on parse or normalize failure.
pub fn normalize_to_string(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant =
        parse_hgvs(input).unwrap_or_else(|e| panic!("Failed to parse '{}': {}", input, e));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("Normalization failed for '{}': {}", input, e));
    format!("{}", normalized)
}

/// Normalize `input` against the padded synthetic genomic contig built from
/// `core`, assert the result denotes the same bases the input did, and return
/// the rendered output.
///
/// The comparison is projected through [`hgvs_to_spdi`] rather than the
/// normalizer's own applier, so a pass that silently changed the sequence shows
/// up here even when every encoding agrees with every other. Members that
/// overlap have no resulting sequence at all, which this reports as a failure
/// rather than as equality — that is the shape most of the cis-allele issue
/// regressions are about.
pub fn assert_padded_preserving(core: &str, input: &str) -> String {
    let provider = SyntheticBuilder::genomic(core).build();
    let normalizer = Normalizer::new(provider.clone());
    let reference = padded(core);
    let output = normalizer
        .normalize(&parse_hgvs(input).expect("parse"))
        .expect("normalize")
        .to_string();

    let denote = |descriptor: &str| -> Option<String> {
        let members: Vec<HgvsVariant> = match parse_hgvs(descriptor).ok()? {
            HgvsVariant::Allele(allele) => allele.variants.clone(),
            single => vec![single],
        };
        let mut triples = Vec::new();
        for member in &members {
            triples.push(hgvs_to_spdi(member, &provider).ok()?);
        }
        // 3'-to-5', and at a tied position the *longer* deletion first. A
        // deletion and an insertion can both touch one interbase (`262_263insA`
        // and `263del`); applying the zero-width one first leaves `claimed` at
        // that interbase and makes the deletion look like it overruns, so this
        // applier would report a legal input as having no resulting sequence.
        triples.sort_by_key(|t| std::cmp::Reverse((t.position, t.deletion.len())));
        let mut edited = reference.as_bytes().to_vec();
        let mut claimed = reference.len();
        for triple in &triples {
            let start = usize::try_from(triple.position).ok()?;
            let end = start.checked_add(triple.deletion.len())?;
            if end > claimed {
                return None; // overlapping members
            }
            edited.splice(start..end, triple.insertion.bytes());
            claimed = start;
        }
        String::from_utf8(edited).ok()
    };

    let from_input = denote(input).expect("input applies");
    let from_output = denote(&output)
        .unwrap_or_else(|| panic!("`{input}` -> `{output}` has no resulting sequence"));
    assert_eq!(
        from_output, from_input,
        "`{input}` -> `{output}` changed the sequence"
    );

    // Denoting the same bases is not enough: the members must also render
    // disjoint and ascending. An overlapping pair that happens to apply to the
    // right sequence is still malformed, and that is the shape most of the
    // cis-allele regressions are about.
    let members: Vec<HgvsVariant> = match parse_hgvs(&output).expect("output parses") {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    };
    let spans: Vec<(u64, u64)> = members
        .iter()
        .map(|member| {
            let triple = hgvs_to_spdi(member, &provider).expect("member has an SPDI");
            (
                triple.position,
                triple.position + triple.deletion.len() as u64,
            )
        })
        .collect();
    for (index, pair) in spans.windows(2).enumerate() {
        assert!(
            pair[0] != pair[1] && pair[1].0 >= pair[0].1,
            "`{input}` -> `{output}`: member {} at interbase {} overlaps or does not \
             follow the previous member ending at {} (spans {spans:?})",
            index + 1,
            pair[1].0,
            pair[0].1
        );
    }

    // The output is a fixed point. Several call sites' comments state this as
    // something this helper measured, and until now it measured nothing of the
    // kind: denotation equality and member disjointness both look only at where
    // the *first* pass landed. Preservation that is not a fixed point means a
    // second pass re-partitions the member the first pass kept, which is a
    // rule change these tests are meant to catch and which would otherwise stay
    // green while the comments claimed otherwise.
    let again = normalizer
        .normalize(&parse_hgvs(&output).expect("output parses"))
        .expect("normalize")
        .to_string();
    assert_eq!(
        again, output,
        "`{input}` -> `{output}` is not a fixed point: a second pass reaches `{again}`"
    );

    output
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'a' => 't',
            't' => 'a',
            'g' => 'c',
            'c' => 'g',
            other => other,
        })
        .collect()
}

fn rna_to_dna(seq: &str) -> String {
    seq.chars()
        .map(|c| match c {
            'u' => 't',
            'U' => 'T',
            other => other,
        })
        .collect::<String>()
        .to_uppercase()
}
