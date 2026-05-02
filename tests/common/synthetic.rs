//! Synthetic `MockProvider` builders for nucleotide-shift coverage matrices.
//!
//! See `docs/superpowers/specs/2026-05-01-A7-ins-shift-coverage-matrix-design.md`.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

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

fn padded(core: &str) -> String {
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
