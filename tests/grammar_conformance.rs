//! Grammar conformance tests against Mutalyzer grammar rules
//!
//! Each test category corresponds to a grammar production from
//! external-repos/mutalyzer-hgvs-parser/mutalyzer_hgvs_parser/ebnf/

use ferro_hgvs::parse_hgvs;

/// Test cases derived from dna.g: substitution rule
#[test]
fn test_grammar_substitution() {
    let cases = vec![
        // Basic substitution
        "NC_000001.11:g.12345A>G",
        // With sequence context
        "NC_000001.11:g.12345CAT>G",
        // Multi-base replacement
        "NC_000001.11:g.12345_12347CAT>GGGG",
    ];

    for case in cases {
        assert!(parse_hgvs(case).is_ok(), "Failed to parse: {}", case);
    }
}

/// Test cases derived from dna.g: deletion rule
#[test]
fn test_grammar_deletion() {
    let cases = vec![
        "NC_000001.11:g.12345del",
        "NC_000001.11:g.12345_12350del",
        "NC_000001.11:g.12345delA",
        "NC_000001.11:g.12345_12347delCAT",
    ];

    for case in cases {
        assert!(parse_hgvs(case).is_ok(), "Failed to parse: {}", case);
    }
}

/// Test cases derived from dna.g: insertion rule
#[test]
fn test_grammar_insertion() {
    let cases = vec![
        "NC_000001.11:g.12345_12346insA",
        "NC_000001.11:g.12345_12346insACGT",
        // Insertion with count
        "NC_000001.11:g.12345_12346ins10",
        // Insertion with uncertain length
        "NC_000001.11:g.12345_12346ins(5_10)",
    ];

    for case in cases {
        assert!(parse_hgvs(case).is_ok(), "Failed to parse: {}", case);
    }
}

/// Test cases derived from dna.g: duplication rule
#[test]
fn test_grammar_duplication() {
    let cases = vec![
        "NC_000001.11:g.12345dup",
        "NC_000001.11:g.12345_12350dup",
        "NC_000001.11:g.12345dupA",
    ];

    for case in cases {
        assert!(parse_hgvs(case).is_ok(), "Failed to parse: {}", case);
    }
}

/// Test cases derived from dna.g: deletion_insertion (delins) rule
#[test]
fn test_grammar_delins() {
    let cases = vec![
        "NC_000001.11:g.12345delinsGG",
        "NC_000001.11:g.12345_12347delinsA",
        "NC_000001.11:g.12345_12347delCATinsGGGG",
    ];

    for case in cases {
        assert!(parse_hgvs(case).is_ok(), "Failed to parse: {}", case);
    }
}

/// Test cases derived from dna.g: inversion rule
#[test]
fn test_grammar_inversion() {
    let cases = vec!["NC_000001.11:g.12345_12350inv"];

    for case in cases {
        assert!(parse_hgvs(case).is_ok(), "Failed to parse: {}", case);
    }
}

/// Test cases derived from dna.g: repeat rule
#[test]
fn test_grammar_repeat() {
    let cases = vec![
        "NC_000001.11:g.12345CAG[20]",
        "NC_000001.11:g.12345_12350[15]",
    ];

    for case in cases {
        assert!(parse_hgvs(case).is_ok(), "Failed to parse: {}", case);
    }
}

/// Test cases derived from common.g: uncertain locations
#[test]
fn test_grammar_uncertain_locations() {
    let cases = vec![
        // Uncertain point
        "NC_000001.11:g.(12345_12350)del",
        // Uncertain range
        "NC_000001.11:g.(12340_12345)_(12350_12355)del",
    ];

    for case in cases {
        assert!(parse_hgvs(case).is_ok(), "Failed to parse: {}", case);
    }
}

/// Test cases derived from protein.g rules
#[test]
fn test_grammar_protein() {
    let cases = vec![
        // Substitution
        "NP_000001.1:p.Ala123Gly",
        "NP_000001.1:p.A123G",
        // Deletion
        "NP_000001.1:p.Ala123del",
        "NP_000001.1:p.Ala123_Gly125del",
        // Insertion
        "NP_000001.1:p.Ala123_Gly124insLeu",
        // Duplication
        "NP_000001.1:p.Ala123dup",
        // Frameshift
        "NP_000001.1:p.Ala123fs",
        "NP_000001.1:p.Ala123Glyfs*10",
        // Extension
        "NP_000001.1:p.Met1ext-5",
        "NP_000001.1:p.*123Glnext*50",
        // Synonymous
        "NP_000001.1:p.Ala123=",
        // Unknown
        "NP_000001.1:p.Ala123?",
    ];

    for case in cases {
        assert!(parse_hgvs(case).is_ok(), "Failed to parse: {}", case);
    }
}
