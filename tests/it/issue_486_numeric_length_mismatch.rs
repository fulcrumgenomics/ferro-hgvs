//! Issue #486 — numeric length-mismatch (`del<N>`/`dup<N>`/`del<N>ins…`) where
//! the stated count disagrees with the position span must reject under strict
//! normalize as `ReferenceMismatch` (mutalyzer ELENGTHMISMATCH). Confirms all
//! three edit kinds reach `validate_reference` with the span intact.
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

fn tx_provider() -> MockProvider {
    let mut p = MockProvider::new();
    let seq = "ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT".to_string();
    let len = seq.len() as u64;
    let tx = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        seq,
        Some(1),
        Some(len),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    p.add_transcript(tx);
    p
}

fn assert_rejects(input: &str) {
    let normalizer = Normalizer::with_config(tx_provider(), NormalizeConfig::strict());
    let v = parse_hgvs(input).expect("parse");
    let err = normalizer
        .normalize(&v)
        .expect_err("strict normalize must reject numeric length mismatch");
    assert!(
        format!("{err:?}").contains("ReferenceMismatch"),
        "{input}: expected ReferenceMismatch, got {err:?}"
    );
}

#[test]
fn del_numeric_length_mismatch_rejects() {
    assert_rejects("NM_TEST.1:c.5_6del4"); // span 2 != 4
}

#[test]
fn dup_numeric_length_mismatch_rejects() {
    assert_rejects("NM_TEST.1:c.5_6dup4"); // span 2 != 4
}

#[test]
fn delins_numeric_deleted_length_mismatch_rejects() {
    assert_rejects("NM_TEST.1:c.5del4insATC"); // span 1 != 4
}

#[test]
fn matching_numeric_length_does_not_reject() {
    // c.5_8del4: span 4 == 4 — must normalize cleanly even in strict mode.
    let normalizer = Normalizer::with_config(tx_provider(), NormalizeConfig::strict());
    let v = parse_hgvs("NM_TEST.1:c.5_8del4").expect("parse");
    assert!(
        normalizer.normalize(&v).is_ok(),
        "c.5_8del4 (span 4 == 4) must not reject"
    );
}
