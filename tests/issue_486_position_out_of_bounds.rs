//! Issue #486 — a position past the transcript/CDS end (mutalyzer
//! EOUTOFBOUNDARY) must reject under strict normalize as
//! `FerroError::InvalidCoordinates`. The bounds check (`check_cds_pos_past_end`
//! et al.) already exists; this pins the strict-mode rejection that the errors
//! axis relies on, independent of the manifest.
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

fn tx_provider() -> MockProvider {
    let mut p = MockProvider::new();
    // 60 bp transcript, CDS spans the whole thing (1..=60).
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

#[test]
fn cds_position_past_end_rejects_strict() {
    // c.100 is past the 60 bp CDS-end → EOUTOFBOUNDARY.
    let normalizer = Normalizer::with_config(tx_provider(), NormalizeConfig::strict());
    let v = parse_hgvs("NM_TEST.1:c.100A>C").expect("parse");
    let err = normalizer
        .normalize(&v)
        .expect_err("strict normalize must reject a position past CDS-end");
    assert!(
        format!("{err:?}").contains("InvalidCoordinates"),
        "expected InvalidCoordinates, got {err:?}"
    );
}

#[test]
fn in_bounds_position_does_not_reject_strict() {
    // c.10 is well within the 60 bp CDS (and c.10 == 'C', so the stated ref
    // matches — isolates the bounds check from a ref mismatch). Must pass.
    let normalizer = Normalizer::with_config(tx_provider(), NormalizeConfig::strict());
    let v = parse_hgvs("NM_TEST.1:c.10C>A").expect("parse");
    assert!(
        normalizer.normalize(&v).is_ok(),
        "in-bounds c.10 must not reject"
    );
}
