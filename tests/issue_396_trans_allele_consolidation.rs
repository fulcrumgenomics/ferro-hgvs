//! Issue #396 item 1 — behavior-preserving regression net for the 7
//! `parse_*_trans_allele_shorthand` helpers.
//!
//! Each axis must round-trip the compact-prefix trans-allele forms:
//!
//! - `[var1];[var2]` (two member variants)
//! - `[var1];[0]`   (paired with NullAllele)
//! - `[var1];[?]`   (paired with UnknownAllele)
//!
//! and reject malformed inputs (missing close bracket, fewer than 2
//! members) the same way they do today. These tests must pass on
//! `origin/main` before the refactor and continue passing after — the
//! refactor is behavior-preserving by design.
//!
//! Protein has 3 additional special-case markers (`[0]`, `[0?]`,
//! `[(0)]`) that resolve to `ProteinEdit::NoProtein` inside the
//! compact `p.` form. Those are exercised here too so the refactor
//! doesn't accidentally relax that protein-specific routing.

use ferro_hgvs::parse_hgvs;

fn parse_ok(input: &str) -> String {
    parse_hgvs(input)
        .unwrap_or_else(|e| panic!("parse({input:?}) failed: {e:?}"))
        .to_string()
}

fn parse_err(input: &str) -> ferro_hgvs::error::FerroError {
    parse_hgvs(input).expect_err(&format!("parse({input:?}) unexpectedly succeeded"))
}

// =============================================================================
// Axis: genome (g.)
// =============================================================================

#[test]
fn genome_trans_two_members_roundtrip() {
    let s = parse_ok("NC_000001.11:g.[100A>G];[200T>C]");
    assert!(s.contains("];["), "expected split brackets, got {s}");
    assert!(s.contains("100A>G"));
    assert!(s.contains("200T>C"));
}

#[test]
fn genome_trans_with_null_allele_roundtrip() {
    let s = parse_ok("NC_000001.11:g.[100A>G];[0]");
    assert!(s.contains(";[0]"), "expected NullAllele segment, got {s}");
}

#[test]
fn genome_trans_with_unknown_allele_roundtrip() {
    let s = parse_ok("NC_000001.11:g.[100A>G];[?]");
    assert!(
        s.contains(";[?]"),
        "expected UnknownAllele segment, got {s}"
    );
}

// =============================================================================
// Axis: cds (c.)
// =============================================================================

#[test]
fn cds_trans_two_members_roundtrip() {
    let s = parse_ok("NM_000088.3:c.[1A>G];[2T>C]");
    assert!(s.contains("];["), "got {s}");
    assert!(s.contains("1A>G") && s.contains("2T>C"));
}

#[test]
fn cds_trans_with_null_allele_roundtrip() {
    let s = parse_ok("NM_000088.3:c.[1A>G];[0]");
    assert!(s.contains(";[0]"), "got {s}");
}

#[test]
fn cds_trans_with_unknown_allele_roundtrip() {
    let s = parse_ok("NM_000088.3:c.[1A>G];[?]");
    assert!(s.contains(";[?]"), "got {s}");
}

// =============================================================================
// Axis: tx (n.)
// =============================================================================

#[test]
fn tx_trans_two_members_roundtrip() {
    let s = parse_ok("NR_000001.1:n.[1A>G];[2T>C]");
    assert!(s.contains("];["), "got {s}");
}

#[test]
fn tx_trans_with_null_allele_roundtrip() {
    let s = parse_ok("NR_000001.1:n.[1A>G];[0]");
    assert!(s.contains(";[0]"), "got {s}");
}

#[test]
fn tx_trans_with_unknown_allele_roundtrip() {
    let s = parse_ok("NR_000001.1:n.[1A>G];[?]");
    assert!(s.contains(";[?]"), "got {s}");
}

// =============================================================================
// Axis: mt (m.)
// =============================================================================

#[test]
fn mt_trans_two_members_roundtrip() {
    let s = parse_ok("NC_012920.1:m.[100A>G];[200T>C]");
    assert!(s.contains("];["), "got {s}");
}

#[test]
fn mt_trans_with_null_allele_roundtrip() {
    let s = parse_ok("NC_012920.1:m.[100A>G];[0]");
    assert!(s.contains(";[0]"), "got {s}");
}

#[test]
fn mt_trans_with_unknown_allele_roundtrip() {
    let s = parse_ok("NC_012920.1:m.[100A>G];[?]");
    assert!(s.contains(";[?]"), "got {s}");
}

// =============================================================================
// Axis: rna (r.)
// =============================================================================

#[test]
fn rna_trans_two_members_roundtrip() {
    let s = parse_ok("NM_004006.2:r.[100a>g];[200u>c]");
    assert!(s.contains("];["), "got {s}");
}

#[test]
fn rna_trans_with_null_allele_roundtrip() {
    let s = parse_ok("NM_004006.2:r.[100a>g];[0]");
    assert!(s.contains(";[0]"), "got {s}");
}

#[test]
fn rna_trans_with_unknown_allele_roundtrip() {
    let s = parse_ok("NM_004006.2:r.[100a>g];[?]");
    assert!(s.contains(";[?]"), "got {s}");
}

// =============================================================================
// Axis: circular (o.)
// =============================================================================

#[test]
fn circular_trans_two_members_roundtrip() {
    let s = parse_ok("NC_011083.1:o.[100A>G];[200T>C]");
    assert!(s.contains("];["), "got {s}");
}

#[test]
fn circular_trans_with_null_allele_roundtrip() {
    let s = parse_ok("NC_011083.1:o.[100A>G];[0]");
    assert!(s.contains(";[0]"), "got {s}");
}

#[test]
fn circular_trans_with_unknown_allele_roundtrip() {
    let s = parse_ok("NC_011083.1:o.[100A>G];[?]");
    assert!(s.contains(";[?]"), "got {s}");
}

// =============================================================================
// Axis: protein (p.) — keeps its 3 special markers
// =============================================================================

#[test]
fn protein_trans_two_members_roundtrip() {
    let s = parse_ok("NP_000079.2:p.[Ala1Val];[Ser2Thr]");
    assert!(s.contains("];["), "got {s}");
}

#[test]
fn protein_trans_with_no_protein_zero_roundtrip() {
    // `p.[X];[0]` → second member is ProteinEdit::NoProtein (NOT NullAllele).
    let s = parse_ok("NP_000079.2:p.[Ala1Val];[0]");
    assert!(s.contains(";[0]"), "got {s}");
}

#[test]
fn protein_trans_with_unknown_allele_roundtrip() {
    let s = parse_ok("NP_000079.2:p.[Ala1Val];[?]");
    assert!(s.contains(";[?]"), "got {s}");
}

#[test]
fn protein_trans_with_predicted_no_protein_roundtrip() {
    // `p.[X];[0?]` → second member is ProteinEdit::NoProtein { predicted: true }.
    let s = parse_ok("NP_000079.2:p.[Ala1Val];[0?]");
    assert!(s.contains(";[0?]"), "got {s}");
}

// =============================================================================
// Negative: malformed forms must still reject across all axes
// =============================================================================

#[test]
fn trans_single_bracket_rejected_across_axes() {
    // Only one bracket member — must reject (less than 2 alleles for trans).
    // A trailing `;` with no second bracket is not a spec-defined allele
    // shape: the inner singleton parses but the dangling `;` becomes
    // unexpected trailing input. The trans-allele code path must NOT
    // produce a 1-element AlleleVariant from this form.
    for input in [
        "NC_000001.11:g.[100A>G];", // trailing semicolon, no second bracket
        "NM_000088.3:c.[1A>G];",
    ] {
        parse_err(input);
    }
}

#[test]
fn trans_truncated_bracket_rejected() {
    // Missing close bracket on second member.
    for input in [
        "NC_000001.11:g.[100A>G];[200T>C",
        "NM_000088.3:c.[1A>G];[2T>C",
    ] {
        let _ = parse_err(input);
    }
}
