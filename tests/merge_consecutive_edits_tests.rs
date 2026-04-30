//! Tests for HGVS-spec consecutive-edit merging in alleles.
//! See docs/superpowers/specs/2026-04-30-merge-consecutive-allele-edits-design.md.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

fn normalize_to_string(input: &str) -> String {
    let normalizer = Normalizer::new(MockProvider::new());
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

#[test]
fn test_merge_consecutive_subs_genome() {
    // Issue #72 example: two adjacent SNVs collapse to a single delins.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000G>A;1001A>C]"),
        "NC_000001.11:g.1000_1001delinsAC",
    );
}

#[test]
fn test_merge_consecutive_dels_genome() {
    // Issue #72 example: two adjacent single-nt deletions become a single ranged del.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000del;1001del]"),
        "NC_000001.11:g.1000_1001del",
    );
}

#[test]
fn test_merge_sub_then_del() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000G>A;1001del]"),
        "NC_000001.11:g.1000_1001delinsA",
    );
}

#[test]
fn test_merge_del_then_sub() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000del;1001A>C]"),
        "NC_000001.11:g.1000_1001delinsC",
    );
}

#[test]
fn test_merge_dels_drops_explicit_sequence() {
    // Per design doc: del+del with explicit ref sequences emits the no-sequence form.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000delA;1001delC]"),
        "NC_000001.11:g.1000_1001del",
    );
}

#[test]
fn test_merge_dels_drops_length() {
    // Per design doc: del+del with length specifiers emits no-sequence/no-length form.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000_1002del3;1003_1004del2]"),
        "NC_000001.11:g.1000_1004del",
    );
}

#[test]
fn test_merge_sub_then_delins() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000G>A;1001_1002delinsTT]"),
        "NC_000001.11:g.1000_1002delinsATT",
    );
}

#[test]
fn test_merge_delins_then_sub() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000_1001delinsTT;1002A>C]"),
        "NC_000001.11:g.1000_1002delinsTTC",
    );
}

#[test]
fn test_merge_delins_then_delins() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000_1001delinsTT;1002_1003delinsAA]"),
        "NC_000001.11:g.1000_1003delinsTTAA",
    );
}

#[test]
fn test_merge_skips_non_literal_delins() {
    // delins with a non-Literal payload (e.g., delins10) is not safely
    // concatenable; the design doc requires the pair to pass through.
    let input = "NC_000001.11:g.[1000G>A;1001_1010delins10]";
    let result = normalize_to_string(input);
    // The output must still contain both edits separately (unchanged).
    assert!(result.contains("1000G>A"), "expected 1000G>A in {}", result);
    assert!(
        result.contains("1001_1010delins"),
        "expected 1001_1010delins in {}",
        result
    );
    assert!(result.contains(';'), "expected separator in {}", result);
}
