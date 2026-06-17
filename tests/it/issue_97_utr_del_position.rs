//! Regression test for issue #97: 5'UTR deletions on minus-strand
//! transcripts collapsed to `c.?del`, discarding positional information.
//!
//! Root cause is an off-by-one in the CDS↔tx coordinate conversions
//! (`cds_to_tx_pos` and `tx_to_cds_pos`) for 5'UTR positions: HGVS
//! coordinates skip from `c.-1` directly to `c.1` (no `c.0`), but the
//! conversion formulas treated the gap as if `c.0` existed. As a
//! result, after 3'-shifting moved a single-base UTR deletion to the
//! rightmost position in its homopolymer tract, the inverse mapping
//! produced `base = 0` — which `CdsPos::Display` renders as `?`.
//!
//! These tests fail before the fix and pass after. Independent of the
//! symptom-locked tests in `tests/coverage_gap_tests.rs`.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

const PAD: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
const PAD_OFFSET: u64 = 256;

/// Single-exon minus-strand transcript whose first 5 transcript bases
/// are `AAAAA` (5'UTR; cds_start = 6 means `c.-5..c.-1` map to
/// `tx 1..5`). The genomic sequence is the reverse complement of the
/// transcript view: `TTTTT` at PAD..PAD+4 becomes `AAAAA` in
/// transcript view at c.-5..c.-1.
fn make_minus_strand_utr_fixture() -> MockProvider {
    let core = "TTTTTGCATGCATGCATGCATTTTT";
    let genomic_seq = format!("{}{}{}", PAD, core, PAD);

    let tx_seq = "AAAAATGCATGCATGCATGCAAAAA";

    let p = PAD_OFFSET + 1;
    let transcript = Transcript::new(
        "NM_MUTR.1".to_string(),
        Some("MUTRGENE".to_string()),
        Strand::Minus,
        tx_seq.to_string(),
        Some(6),
        Some(20),
        vec![Exon::with_genomic(1, 1, 25, p, p + 24)],
        Some("chr_mutr".to_string()),
        Some(p),
        Some(p + 24),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );

    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr_mutr", genomic_seq);
    provider.add_transcript(transcript);
    provider
}

fn normalize(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).expect("parse should succeed");
    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize should succeed");
    format!("{}", normalized)
}

/// **Diagnostic 1 — 5'UTR del 3'-shift on minus strand.**
///
/// `c.-3del` (a single A in the c.-5..c.-1 homopolymer) must shift to
/// the 3'-most A in the run, c.-1. Before the fix it collapsed to
/// `c.?del`.
#[test]
fn five_prime_utr_del_shifts_to_minus_one() {
    let provider = make_minus_strand_utr_fixture();
    let result = normalize(provider, "NM_MUTR.1:c.-3del");
    assert_eq!(result, "NM_MUTR.1:c.-1del");
}

/// **Diagnostic 2 — round-trip of c.-1del.**
///
/// Re-normalizing the canonical form must be a no-op. If the inverse
/// `tx_to_cds_pos` mapping is correct, a second pass through normalize
/// produces the same output.
#[test]
fn five_prime_utr_del_round_trip() {
    let provider = make_minus_strand_utr_fixture();
    let once = normalize(provider, "NM_MUTR.1:c.-1del");
    assert_eq!(once, "NM_MUTR.1:c.-1del");
}

/// **Diagnostic 3 — substitution at c.-3 round-trip.**
///
/// Substitutions do not shift, so the round-trip must preserve the
/// position exactly. The forward-direction off-by-one is invisible to
/// substitutions on a homopolymer (every UTR base is the same letter),
/// so this test explicitly pins both the forward and inverse mappings.
#[test]
fn five_prime_utr_substitution_round_trip() {
    let provider = make_minus_strand_utr_fixture();
    let result = normalize(provider, "NM_MUTR.1:c.-3A>G");
    assert_eq!(result, "NM_MUTR.1:c.-3A>G");
}
