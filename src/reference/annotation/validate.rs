//! Optional FASTA-aware validation. See spec §6 Stage 5.
//!
//! Checks performed (both are warnings — never reject, never drop):
//!
//! - `W-LOAD-200 CdsLengthNotMod3` — the CDS length is not divisible by 3.
//! - `W-LOAD-201 NonCanonicalStartCodon` — the first codon is not ATG/CTG/GTG/TTG.

use std::path::Path;

use crate::reference::provider::ReferenceProvider;
use crate::reference::transcript::{Strand, Transcript};

use super::diagnostics::{DiagnosticPayload, LoaderDiagnostic, LoaderReport, SourceLocation};

/// Validate a single transcript against a FASTA reference.
///
/// Both diagnostics are warnings; this function never drops a transcript or
/// returns an error.
pub(crate) fn validate_transcript(
    tx: &Transcript,
    fasta: &dyn ReferenceProvider,
    source_path: &Path,
    report: &mut LoaderReport,
) {
    let (cds_start, cds_end) = match (tx.cds_start, tx.cds_end) {
        (Some(s), Some(e)) => (s, e),
        _ => return, // non-coding transcript — nothing to check
    };

    // CDS length in transcript coordinates (1-based inclusive).
    let cds_len = cds_end.saturating_sub(cds_start).saturating_add(1);

    if cds_len % 3 != 0 {
        report.record(LoaderDiagnostic::warning(
            "W-LOAD-200",
            format!("CdsLengthNotMod3: tx {} cds_length={}", tx.id, cds_len),
            SourceLocation {
                path: source_path.to_path_buf(),
                line: 0,
            },
            Some(tx.id.clone()),
            DiagnosticPayload::CdsLengthNotMod3 {
                transcript_id: tx.id.clone(),
                length: cds_len,
            },
        ));
    }

    // Locate the chromosome — required for genomic sequence lookup.
    let chrom = match tx.chromosome.as_deref() {
        Some(c) => c,
        None => return,
    };

    // Build the start codon from spliced CDS bases. A start codon can be split
    // 1+2 or 2+1 across an exon junction, so walking only the exon that holds
    // `cds_start` would fetch intronic bases and produce false W-LOAD-201
    // warnings. `extract_first_codon` collects the first three CDS bases in
    // transcript order across exon boundaries and strand-corrects.
    let codon = extract_first_codon(tx, chrom, fasta, cds_start);

    if let Some(c) = codon {
        let codon_uc = c.to_uppercase();
        let canonical = matches!(codon_uc.as_str(), "ATG" | "CTG" | "GTG" | "TTG");
        if !canonical {
            report.record(LoaderDiagnostic::warning(
                "W-LOAD-201",
                format!("NonCanonicalStartCodon: tx {} codon={}", tx.id, codon_uc),
                SourceLocation {
                    path: source_path.to_path_buf(),
                    line: 0,
                },
                Some(tx.id.clone()),
                DiagnosticPayload::NonCanonicalStartCodon {
                    transcript_id: tx.id.clone(),
                    codon: codon_uc,
                },
            ));
        }
    }
}

/// Assemble the first three CDS bases (the start codon) in transcript order,
/// strand-corrected. Returns `None` when the transcript has no genomic exon
/// coordinates, the strand is unknown, the FASTA fetch fails, or fewer than
/// three CDS bases are available (e.g. a CDS that runs off the end of the
/// transcript).
///
/// Iterates exons in transcript order from the one containing `cds_start`,
/// taking up to the remaining-bases-needed from each exon. This handles
/// start codons split 1+2 or 2+1 across an exon junction without ever
/// fetching intronic bases.
fn extract_first_codon(
    tx: &Transcript,
    chrom: &str,
    fasta: &dyn ReferenceProvider,
    cds_start: u64,
) -> Option<String> {
    const TARGET: u64 = 3;
    let mut codon = String::with_capacity(3);
    let mut tx_pos = cds_start;

    for ex in &tx.exons {
        if codon.len() as u64 >= TARGET {
            break;
        }
        // Skip exons entirely upstream of the current transcript position.
        if ex.end < tx_pos {
            continue;
        }
        let (gs, ge) = match (ex.genomic_start, ex.genomic_end) {
            (Some(gs), Some(ge)) => (gs, ge),
            _ => return None,
        };
        let local_start_tx = tx_pos.max(ex.start);
        if local_start_tx > ex.end {
            continue;
        }
        let remaining = TARGET - codon.len() as u64;
        let avail = ex.end - local_start_tx + 1;
        let take = remaining.min(avail);
        let local_end_tx = local_start_tx + take - 1;

        // Map the transcript-coordinate slice [local_start_tx, local_end_tx]
        // to a contiguous genomic range and fetch. `get_sequence` is 0-based
        // half-open: 1-based genomic g → start = g - 1.
        let chunk = match tx.strand {
            Strand::Plus => {
                let g_lo = gs + (local_start_tx - ex.start);
                let g_hi = gs + (local_end_tx - ex.start);
                fasta
                    .get_sequence(chrom, g_lo.saturating_sub(1), g_hi)
                    .ok()?
            }
            Strand::Minus => {
                // tx position increases ⇒ genomic decreases on minus strand.
                // local_start_tx ↦ higher genomic; local_end_tx ↦ lower genomic.
                let g_hi = ge - (local_start_tx - ex.start);
                let g_lo = ge - (local_end_tx - ex.start);
                let raw = fasta
                    .get_sequence(chrom, g_lo.saturating_sub(1), g_hi)
                    .ok()?;
                reverse_complement(&raw)
            }
            Strand::Unknown => return None,
        };
        codon.push_str(&chunk);
        tx_pos = local_end_tx + 1;
    }

    (codon.len() as u64 == TARGET).then_some(codon)
}

/// Return the reverse complement of a DNA string.
fn reverse_complement(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'G' | 'g' => 'C',
            'C' | 'c' => 'G',
            'N' | 'n' => 'N',
            other => other,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::provider::ReferenceProvider;
    use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Transcript};
    use std::sync::OnceLock;

    /// Minimal stub that serves a single in-memory sequence for any chromosome name.
    struct StubFasta {
        seq: String,
    }

    impl ReferenceProvider for StubFasta {
        fn get_transcript(&self, id: &str) -> Result<Transcript, crate::error::FerroError> {
            Err(crate::error::FerroError::ReferenceNotFound { id: id.to_string() })
        }

        fn get_sequence(
            &self,
            _chrom: &str,
            start: u64,
            end: u64,
        ) -> Result<String, crate::error::FerroError> {
            Ok(self.seq[(start as usize)..(end as usize)].to_string())
        }
    }

    /// Build a minimal plus-strand transcript with transcript-space CDS bounds.
    ///
    /// The single exon spans positions 1–100 in both transcript and genomic space
    /// (i.e., genomic_start = tx.start = 1, genomic_end = tx.end = 100).
    fn tx_with_cds(cds_start: u64, cds_end: u64) -> Transcript {
        Transcript {
            id: "tx1".into(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: None,
            cds_start: Some(cds_start),
            cds_end: Some(cds_end),
            exons: vec![Exon::with_genomic(1, 1, 100, 1, 100)],
            chromosome: Some("chr1".into()),
            genomic_start: Some(1),
            genomic_end: Some(100),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: vec![],
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn warns_when_cds_length_not_mod3() {
        // cds_start=1, cds_end=5 → cds_len = 5, not divisible by 3
        let tx = tx_with_cds(1, 5);
        let fa = StubFasta {
            seq: "ATGCCCAAATAGNNNNNNN".to_string(),
        };
        let mut report = LoaderReport::default();
        validate_transcript(&tx, &fa, Path::new("t"), &mut report);
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-200"), Some(&1));
    }

    #[test]
    fn does_not_warn_when_cds_length_mod3_and_atg_start() {
        // cds_start=1, cds_end=9 → cds_len = 9 (divisible by 3)
        // StubFasta seq: "ATGCCCAAATAGNNNNNNN"
        // Plus strand, gs = genomic_start_of_cds = 1 + (1-1) = 1
        // fetch(0-based): start = gs-1 = 0, end = gs+2 = 3 → "ATG" ✓
        let tx = tx_with_cds(1, 9);
        let fa = StubFasta {
            seq: "ATGCCCAAATAGNNNNNNN".to_string(),
        };
        let mut report = LoaderReport::default();
        validate_transcript(&tx, &fa, Path::new("t"), &mut report);
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-200"), None);
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-201"), None);
    }

    /// Build a FASTA seq of length `len` with specific bases planted at the
    /// given 1-based genomic positions (everything else is 'N').
    fn seq_with(planted: &[(u64, char)], len: usize) -> String {
        let mut bytes = vec![b'N'; len];
        for &(pos1, base) in planted {
            bytes[(pos1 - 1) as usize] = base as u8;
        }
        String::from_utf8(bytes).unwrap()
    }

    /// Build a two-exon transcript on the requested strand.
    fn two_exon_tx(
        strand: Strand,
        ex1: (u64, u64, u64, u64),
        ex2: (u64, u64, u64, u64),
        cds_start: u64,
        cds_end: u64,
    ) -> Transcript {
        Transcript {
            id: "tx_split".into(),
            gene_symbol: None,
            strand,
            sequence: None,
            cds_start: Some(cds_start),
            cds_end: Some(cds_end),
            exons: vec![
                Exon::with_genomic(1, ex1.0, ex1.1, ex1.2, ex1.3),
                Exon::with_genomic(2, ex2.0, ex2.1, ex2.2, ex2.3),
            ],
            chromosome: Some("chr1".into()),
            genomic_start: Some(ex1.2.min(ex2.2)),
            genomic_end: Some(ex1.3.max(ex2.3)),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: vec![],
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn plus_strand_start_codon_split_1_2_across_exon_junction() {
        // Exon 1: tx [1..1]  → genomic [10..10]   (1 base of codon here)
        // Exon 2: tx [2..10] → genomic [50..58]   (next 2 bases of codon here)
        // Without splice-aware extraction, a single-window fetch would read
        // genomic [10,11,12] = 'A','N','N' and falsely warn.
        let tx = two_exon_tx(Strand::Plus, (1, 1, 10, 10), (2, 10, 50, 58), 1, 9);
        let fa = StubFasta {
            seq: seq_with(&[(10, 'A'), (50, 'T'), (51, 'G')], 100),
        };
        let mut report = LoaderReport::default();
        validate_transcript(&tx, &fa, Path::new("t"), &mut report);
        // mod3 ok (cds_len=9), start codon = ATG → no W-LOAD-201
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-200"), None);
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-201"), None);
    }

    #[test]
    fn plus_strand_start_codon_split_2_1_across_exon_junction() {
        // Exon 1: tx [1..2]  → genomic [10..11]   (first 2 bases of codon)
        // Exon 2: tx [3..9]  → genomic [50..56]   (3rd codon base here)
        let tx = two_exon_tx(Strand::Plus, (1, 2, 10, 11), (3, 9, 50, 56), 1, 9);
        let fa = StubFasta {
            seq: seq_with(&[(10, 'A'), (11, 'T'), (50, 'G')], 100),
        };
        let mut report = LoaderReport::default();
        validate_transcript(&tx, &fa, Path::new("t"), &mut report);
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-201"), None);
    }

    #[test]
    fn minus_strand_start_codon_split_1_2_across_exon_junction() {
        // Exon 1 (transcript order): tx [1..1] → genomic [80..80]
        // Exon 2 (transcript order): tx [2..10] → genomic [20..28]
        // Minus strand: tx pos n maps to (ex.ge - (n - ex.start)).
        //   tx 1 → g 80 → tx base = complement(FASTA[80]) = A iff FASTA[80]='T'
        //   tx 2 → g 28 → tx base = T iff FASTA[28]='A'
        //   tx 3 → g 27 → tx base = G iff FASTA[27]='C'
        let tx = two_exon_tx(Strand::Minus, (1, 1, 80, 80), (2, 10, 20, 28), 1, 9);
        let fa = StubFasta {
            seq: seq_with(&[(80, 'T'), (28, 'A'), (27, 'C')], 100),
        };
        let mut report = LoaderReport::default();
        validate_transcript(&tx, &fa, Path::new("t"), &mut report);
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-201"), None);
    }

    #[test]
    fn minus_strand_start_codon_split_2_1_across_exon_junction() {
        // Exon 1: tx [1..2] → genomic [80..81] (minus: tx 1 → g 81, tx 2 → g 80)
        // Exon 2: tx [3..9] → genomic [20..26] (tx 3 → g 26)
        // Codon = (complement of FASTA[81], FASTA[80], FASTA[26])
        // For ATG: FASTA[81]='T', FASTA[80]='A', FASTA[26]='C'.
        let tx = two_exon_tx(Strand::Minus, (1, 2, 80, 81), (3, 9, 20, 26), 1, 9);
        let fa = StubFasta {
            seq: seq_with(&[(81, 'T'), (80, 'A'), (26, 'C')], 100),
        };
        let mut report = LoaderReport::default();
        validate_transcript(&tx, &fa, Path::new("t"), &mut report);
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-201"), None);
    }

    #[test]
    fn split_start_codon_warns_when_not_canonical() {
        // Confirms the splice-aware path still emits W-LOAD-201 for true
        // non-canonical start codons. Same 1+2 layout as the plus-strand
        // test above but planting GAA instead of ATG.
        let tx = two_exon_tx(Strand::Plus, (1, 1, 10, 10), (2, 10, 50, 58), 1, 9);
        let fa = StubFasta {
            seq: seq_with(&[(10, 'G'), (50, 'A'), (51, 'A')], 100),
        };
        let mut report = LoaderReport::default();
        validate_transcript(&tx, &fa, Path::new("t"), &mut report);
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-201"), Some(&1));
    }

    #[test]
    fn warns_on_non_canonical_start_codon() {
        // cds_start=1, cds_end=9 → cds_len = 9 (mod 3 == 0, no W-LOAD-200)
        // StubFasta seq: "AAATTTGGGNNNNNNNNNN"
        // Plus strand fetch at 0-based [0,3) → "AAA" → non-canonical → W-LOAD-201
        let tx = tx_with_cds(1, 9);
        let fa = StubFasta {
            seq: "AAATTTGGGNNNNNNNNNN".to_string(),
        };
        let mut report = LoaderReport::default();
        validate_transcript(&tx, &fa, Path::new("t"), &mut report);
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-201"), Some(&1));
    }
}
