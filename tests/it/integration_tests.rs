//! Integration tests for the full pipeline

use ferro_hgvs::{parse_hgvs, FerroError, MockProvider, Normalizer};

#[test]
fn test_full_pipeline_parse_and_normalize() {
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variants = vec![
        "NC_000001.11:g.12345A>G",
        "NM_000088.3:c.10del", // Position within mock transcript bounds
        "NM_000088.3:c.10A>G",
    ];

    for input in variants {
        // Parse
        let parsed = parse_hgvs(input);
        assert!(parsed.is_ok(), "Failed to parse: {}", input);

        // Normalize
        let normalized = normalizer.normalize(&parsed.unwrap());
        assert!(normalized.is_ok(), "Failed to normalize: {}", input);
    }
}

#[test]
fn test_error_handling() {
    // Test that errors are properly typed
    let result = parse_hgvs("invalid:::");
    assert!(matches!(result, Err(FerroError::Parse { .. })));

    let result = parse_hgvs("NM_NONEXISTENT.1:c.100del");
    // Parsing should succeed even for non-existent transcripts
    assert!(result.is_ok());
}

#[test]
fn test_display_formatting() {
    let test_cases = vec![
        ("NC_000001.11:g.12345A>G", "NC_000001.11:g.12345A>G"),
        ("NM_000088.3:c.10del", "NM_000088.3:c.10del"),
        ("NM_000088.3:c.10_11insATG", "NM_000088.3:c.10_11insATG"),
        ("NM_000088.3:c.10+5G>A", "NM_000088.3:c.10+5G>A"),
        ("NM_000088.3:c.-20G>A", "NM_000088.3:c.-20G>A"),
        ("NM_000088.3:c.*50G>A", "NM_000088.3:c.*50G>A"),
    ];

    for (input, expected) in test_cases {
        let parsed = parse_hgvs(input).unwrap();
        assert_eq!(format!("{}", parsed), expected);
    }
}

#[test]
fn test_variant_type_identification() {
    use ferro_hgvs::HgvsVariant;

    assert!(matches!(
        parse_hgvs("NC_000001.11:g.12345A>G").unwrap(),
        HgvsVariant::Genome(_)
    ));

    assert!(matches!(
        parse_hgvs("NM_000088.3:c.459A>G").unwrap(),
        HgvsVariant::Cds(_)
    ));

    assert!(matches!(
        parse_hgvs("NR_000001.1:n.100A>G").unwrap(),
        HgvsVariant::Tx(_)
    ));

    assert!(matches!(
        parse_hgvs("NP_000079.2:p.Val600Glu").unwrap(),
        HgvsVariant::Protein(_)
    ));

    assert!(matches!(
        parse_hgvs("NC_012920.1:m.3243A>G").unwrap(),
        HgvsVariant::Mt(_)
    ));
}

#[test]
fn test_mock_provider() {
    use ferro_hgvs::ReferenceProvider;

    let provider = MockProvider::with_test_data();

    // Test that we can get a transcript
    let tx = provider.get_transcript("NM_000088.3");
    assert!(tx.is_ok());

    // Test that we can get sequence
    let seq = provider.get_sequence("NM_000088.3", 0, 3);
    assert!(seq.is_ok());

    // Test non-existent transcript
    let missing = provider.get_transcript("NM_MISSING.1");
    assert!(missing.is_err());
}

#[test]
fn test_serde_roundtrip() {
    use ferro_hgvs::hgvs::variant::Accession;
    use serde_json;

    // Test accession serialization
    let acc = Accession::new("NM", "000088", Some(3));
    let json = serde_json::to_string(&acc).unwrap();
    let back: Accession = serde_json::from_str(&json).unwrap();
    assert_eq!(acc, back);
}
