//! Regression test for issue #98: minus-strand intronic reference-base
//! lookup returned genomic-strand letters, defeating every rule that
//! compared the input edit against the local reference window.
//!
//! These tests fail before the fix in `normalize_intronic_cds` /
//! `normalize_intronic_tx` and pass after. Each one is a single,
//! diagnostic assertion: it does not depend on the symptom-locked tests
//! in `tests/coverage_gap_tests.rs`.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

const PAD: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
const PAD_OFFSET: u64 = 256;

/// Minus-strand transcript with intron 1 chosen so that the transcript-view
/// bases at `c.30+1..c.30+10` are `C, A, A, A, G, G, G, C, C, C` — the
/// reverse complement of the genomic-strand intron `GGGCCCTTTG`. Every
/// orientation rule below asserts a behavior that is only achievable if
/// the normalizer reads the *transcript-view* letters (A,A,A; G,G,G) when
/// canonicalizing — not the genomic-strand letters (T,T,T; C,C,C).
fn make_minus_strand_fixture() -> MockProvider {
    let core = "CATGCATGCATGCATGCATGCATGCATGCA\
                AAATTTAAAT\
                GCATGCATGCATGCATGCATGCATGCATGC\
                GGGCCCTTTG\
                ATGCATGCATGCATGCATGCATGCATGCAT";
    let genomic_seq = format!("{}{}{}", PAD, core, PAD);

    let tx_seq = "ATGCATGCATGCATGCATGCATGCATGCAT\
                   GCATGCATGCATGCATGCATGCATGCATGC\
                   TGCATGCATGCATGCATGCATGCATGCATG";

    let p = PAD_OFFSET + 1;
    let transcript = Transcript::new(
        "NM_MINUS.1".to_string(),
        Some("MINUSGENE".to_string()),
        Strand::Minus,
        tx_seq.to_string(),
        Some(1),
        Some(90),
        vec![
            Exon::with_genomic(1, 1, 30, p + 80, p + 109),
            Exon::with_genomic(2, 31, 60, p + 40, p + 69),
            Exon::with_genomic(3, 61, 90, p, p + 29),
        ],
        Some("chr_minus".to_string()),
        Some(p),
        Some(p + 109),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );

    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr_minus", genomic_seq);
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

/// **Diagnostic 1 — ins → dup canonicalization (#81 A1, intronic).**
///
/// Inserting `A` into the transcript-view `AAA` tract at `c.30+2..c.30+4`
/// must canonicalize to `dup` at the 3'-most A (c.30+4). If the
/// normalizer compares the inserted `A` against the *genomic-strand*
/// letters at this region (`T,T,T`), the rule does not fire and the
/// output is unchanged.
#[test]
fn ins_into_transcript_view_aaa_tract_canonicalizes_to_dup() {
    let provider = make_minus_strand_fixture();
    let result = normalize(provider, "NM_MINUS.1:c.30+2_30+3insA");
    assert_eq!(result, "NM_MINUS.1:c.30+4dup");
}

/// **Diagnostic 2 — single-base dup 3'-shift (#81 A6, intronic).**
///
/// `c.30+2dup` (A in the transcript-view AAA tract) must shift to the
/// 3'-most equivalent position, `c.30+4dup`.
#[test]
fn dup_in_transcript_view_aaa_tract_shifts_three_prime() {
    let provider = make_minus_strand_fixture();
    let result = normalize(provider, "NM_MINUS.1:c.30+2dup");
    assert_eq!(result, "NM_MINUS.1:c.30+4dup");
}

/// **Diagnostic 3 — single-base del 3'-shift (#81 A5, intronic).**
///
/// `c.30+5del` (G in the transcript-view GGG tract at c.30+5..c.30+7)
/// must shift to the 3'-most G in the run. For this 10-base intron
/// ferro renders the resulting position with acceptor-relative notation
/// (positions past the midpoint use `c.31-N`); `c.30+7` and `c.31-4`
/// designate the same intronic base.
#[test]
fn del_in_transcript_view_ggg_tract_shifts_three_prime() {
    let provider = make_minus_strand_fixture();
    let result = normalize(provider, "NM_MINUS.1:c.30+5del");
    assert_eq!(result, "NM_MINUS.1:c.31-4del");
}

/// **Diagnostic 4 — delins identity rewrite (PR #78, intronic).**
///
/// Replacing the transcript-view `AAA` at c.30+2..c.30+4 with `AAA` is
/// an identity edit and PR #78 should rewrite it. If the engine sees
/// the genomic-strand `TTT` as the reference, the identity check fails
/// and the rewrite is skipped.
#[test]
fn delins_identity_on_transcript_view_aaa_tract_rewrites() {
    let provider = make_minus_strand_fixture();
    let result = normalize(provider, "NM_MINUS.1:c.30+2_30+4delinsAAA");
    assert_eq!(result, "NM_MINUS.1:c.30+2_30+4=");
}

/// **Diagnostic 5 — `n.` parity for ins → dup canonicalization.**
///
/// Mirror of `ins_into_transcript_view_aaa_tract_canonicalizes_to_dup` but
/// using `n.` notation, which routes through `normalize_intronic_tx()`
/// instead of `normalize_intronic_cds()`. Both functions share the
/// `flip_intronic_for_strand` helper for the issue #98 orientation fix,
/// so this test guards against regressions on the `n.` (TxVariant) path.
#[test]
fn ins_into_transcript_view_aaa_tract_canonicalizes_to_dup_n() {
    let provider = make_minus_strand_fixture();
    let result = normalize(provider, "NM_MINUS.1:n.30+2_30+3insA");
    assert_eq!(result, "NM_MINUS.1:n.30+4dup");
}
