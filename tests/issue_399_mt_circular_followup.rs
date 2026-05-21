//! End-to-end coverage for issue #399 — the F1 follow-up to PR #380.
//!
//! These tests exercise wraparound-aware span math through the public
//! `_with_provider` APIs. The no-provider `get_indel_length` path
//! returns `None` for wraparound on `m.`/`o.` (pinned by
//! `tests/mito_circular_audit.rs`); these tests pin the spec-correct
//! values when a provider can supply contig length.

use ferro_hgvs::parse_hgvs;
use ferro_hgvs::python_helpers::get_indel_length_with_provider;
use ferro_hgvs::reference::mock::MockProvider;

fn mt_provider() -> MockProvider {
    let mut p = MockProvider::new();
    // NC_012920.1 is 16569 bp; for span math we only need the length,
    // not the bases. Use a placeholder string of the right length.
    p.add_genomic_sequence("NC_012920.1", "A".repeat(16569));
    p
}

#[test]
fn wraparound_del_indel_length_is_spec_correct_with_provider() {
    let v = parse_hgvs("NC_012920.1:m.16569_1del").unwrap();
    let p = mt_provider();
    // 2-nt deletion (positions 16569 and 1); del reports -span.
    assert_eq!(get_indel_length_with_provider(&v, &p), Some(-2));
}

#[test]
fn wraparound_dup_indel_length_is_spec_correct_with_provider() {
    let v = parse_hgvs("NC_012920.1:m.16560_5dup").unwrap();
    let p = mt_provider();
    // Wraparound span = (16569 − 16560 + 1) + 5 = 15.
    assert_eq!(get_indel_length_with_provider(&v, &p), Some(15));
}

#[test]
fn wraparound_delins_indel_length_is_spec_correct_with_provider() {
    let v = parse_hgvs("NC_012920.1:m.16569_1delinsT").unwrap();
    let p = mt_provider();
    // 2-nt window replaced with 1-nt insert; net = 1 − 2 = −1.
    assert_eq!(get_indel_length_with_provider(&v, &p), Some(-1));
}
