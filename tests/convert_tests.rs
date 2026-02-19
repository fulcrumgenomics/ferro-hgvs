//! Coordinate conversion tests

use ferro_hgvs::convert::CoordinateMapper;
use ferro_hgvs::hgvs::location::{CdsPos, TxPos};
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};

fn make_test_transcript() -> Transcript {
    // Structure: 5bp 5'UTR + 30bp CDS + 5bp 3'UTR = 40bp
    Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        "AAAAATGCCCAAAGGGTTTAGGCCCAAAGGGTTATAAA".to_string(),
        Some(6),
        Some(35),
        vec![Exon::new(1, 1, 38)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    )
}

#[test]
fn test_cds_to_tx_basic() {
    let tx = make_test_transcript();
    let mapper = CoordinateMapper::new(&tx);

    // c.1 should map to tx position 6 (first CDS base)
    let result = mapper.cds_to_tx(&CdsPos::new(1)).unwrap();
    assert_eq!(result.base, 6);

    // c.30 should map to tx position 35 (last CDS base)
    let result = mapper.cds_to_tx(&CdsPos::new(30)).unwrap();
    assert_eq!(result.base, 35);
}

#[test]
fn test_cds_to_tx_5utr() {
    let tx = make_test_transcript();
    let mapper = CoordinateMapper::new(&tx);

    // c.-1 should map to tx position 5 (last 5'UTR base)
    let result = mapper.cds_to_tx(&CdsPos::new(-1)).unwrap();
    assert_eq!(result.base, 5);

    // c.-5 should map to tx position 1 (first base)
    let result = mapper.cds_to_tx(&CdsPos::new(-5)).unwrap();
    assert_eq!(result.base, 1);
}

#[test]
fn test_cds_to_tx_3utr() {
    let tx = make_test_transcript();
    let mapper = CoordinateMapper::new(&tx);

    // c.*1 should map to tx position 36 (first 3'UTR base)
    let result = mapper.cds_to_tx(&CdsPos::utr3(1)).unwrap();
    assert_eq!(result.base, 36);
}

#[test]
fn test_tx_to_cds_basic() {
    let tx = make_test_transcript();
    let mapper = CoordinateMapper::new(&tx);

    // tx.10 (in CDS) should map to c.5
    let result = mapper.tx_to_cds(&TxPos::new(10)).unwrap();
    assert_eq!(result.base, 5);
    assert!(!result.utr3);
}

#[test]
fn test_tx_to_cds_5utr() {
    let tx = make_test_transcript();
    let mapper = CoordinateMapper::new(&tx);

    // tx.3 (in 5'UTR) should map to c.-3
    let result = mapper.tx_to_cds(&TxPos::new(3)).unwrap();
    assert_eq!(result.base, -3);
    assert!(!result.utr3);
}

#[test]
fn test_tx_to_cds_3utr() {
    let tx = make_test_transcript();
    let mapper = CoordinateMapper::new(&tx);

    // tx.37 (in 3'UTR) should map to c.*2
    let result = mapper.tx_to_cds(&TxPos::new(37)).unwrap();
    assert_eq!(result.base, 2);
    assert!(result.utr3);
}

#[test]
fn test_cds_to_protein() {
    let tx = make_test_transcript();
    let mapper = CoordinateMapper::new(&tx);

    // c.1-3 should map to p.1 (first codon)
    for pos in 1..=3 {
        let result = mapper.cds_to_protein(&CdsPos::new(pos)).unwrap();
        assert_eq!(result.number, 1, "c.{} should map to p.1", pos);
    }

    // c.4-6 should map to p.2 (second codon)
    for pos in 4..=6 {
        let result = mapper.cds_to_protein(&CdsPos::new(pos)).unwrap();
        assert_eq!(result.number, 2, "c.{} should map to p.2", pos);
    }
}

#[test]
fn test_cds_to_protein_utr_error() {
    let tx = make_test_transcript();
    let mapper = CoordinateMapper::new(&tx);

    // UTR positions should not convert to protein
    assert!(mapper.cds_to_protein(&CdsPos::new(-1)).is_err());
    assert!(mapper.cds_to_protein(&CdsPos::utr3(1)).is_err());
}

#[test]
fn test_roundtrip_conversion() {
    let tx = make_test_transcript();
    let mapper = CoordinateMapper::new(&tx);

    // CDS position roundtrip: cds -> tx -> cds
    for cds_pos in 1..=30 {
        let original = CdsPos::new(cds_pos);
        let tx_pos = mapper.cds_to_tx(&original).unwrap();
        let back = mapper.tx_to_cds(&tx_pos).unwrap();
        assert_eq!(
            original.base, back.base,
            "Roundtrip failed for c.{}",
            cds_pos
        );
    }
}
