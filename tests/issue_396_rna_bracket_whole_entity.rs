//! Issue #396 item 3 — `parse_rna_bracket_member` only knew about
//! position-based edits and the predicted-wrapper `(pos edit)`. Inside
//! an RNA allele bracket (`r.[X];[Y]` or `r.[X;Y]` or `r.[X(;)Y]`),
//! whole-entity edits like `=`, `?`, `0`, `spl`, `spl?`, `(spl)`,
//! `(spl?)` were rejected because they have no positional content.
//!
//! Item 3 dispatches to the whole-entity RNA parsers
//! (`parse_whole_rna_identity`, `parse_whole_rna_unknown`,
//! `parse_rna_no_product`, `parse_rna_splice`,
//! `parse_whole_rna_predicted_splice`) inside the bracket member,
//! attaching a dummy interval at position 1 (mirroring the
//! `parse_rna_variant` dispatcher's handling of the same forms).
//!
//! The two existing trans-allele special cases `[0]` → `NullAllele`
//! and `[?]` → `UnknownAllele` are preserved at the
//! `parse_rna_trans_allele_shorthand` layer (above this bracket
//! member), so a bare `r.[X];[0]` still routes to the
//! cross-coordinate `NullAllele`. Only the cis/unknown-flat paths and
//! members that look like `(0)` / `spl` / `spl?` / `=` / `?` (which
//! never short-circuited in the trans helper) flow into this new
//! bracket-member dispatch.
//!
//! When item 3 lands, `tests/rna_spl_marker.rs:246`
//! (`rna_spl_allele_compound_currently_rejected`) is flipped from a
//! "rejects today" pin to an explicit success assertion.

use ferro_hgvs::parse_hgvs;

fn parse_ok(input: &str) -> String {
    parse_hgvs(input)
        .unwrap_or_else(|e| panic!("parse({input:?}) must succeed; got {e:?}"))
        .to_string()
}

// =============================================================================
// Trans form: `r.[A];[B]` with whole-entity edits as one member
// =============================================================================

#[test]
fn rna_trans_spl_question_then_spl_roundtrip() {
    // The exact case pinned by `rna_spl_marker.rs:246`.
    let s = parse_ok("NM_004006.2:r.[spl?];[spl]");
    assert!(s.contains("spl?"), "expected 'spl?' in {s}");
    assert!(s.contains("spl"), "expected 'spl' in {s}");
    assert!(s.contains("];["), "expected trans-form display, got {s}");
}

#[test]
fn rna_trans_spl_pair_roundtrip() {
    let s = parse_ok("NM_004006.2:r.[spl];[spl]");
    assert!(s.contains("];["));
}

#[test]
fn rna_trans_identity_with_substitution_roundtrip() {
    // `[=];[edit]` — identity (no change) on one strand, edit on the other.
    let s = parse_ok("NM_004006.2:r.[=];[100a>g]");
    assert!(s.contains("];["));
    assert!(s.contains("=") || s.contains("=]"));
}

#[test]
fn rna_trans_predicted_spl_with_substitution_roundtrip() {
    // `[(spl?)];[edit]` — predicted-uncertain splice on one strand.
    let s = parse_ok("NM_004006.2:r.[(spl?)];[100a>g]");
    assert!(s.contains("(spl?)"), "expected '(spl?)' in {s}");
}

#[test]
fn rna_trans_predicted_no_product_with_substitution_roundtrip() {
    // `[(0)];[edit]` — predicted no product (parens form).
    // Bare `[0]` is the cross-coord NullAllele short-circuit upstream;
    // `[(0)]` flows into the bracket member.
    let s = parse_ok("NM_004006.2:r.[(0)];[100a>g]");
    assert!(s.contains("(0)"), "expected '(0)' in {s}");
}

// =============================================================================
// Cis form: `r.[A;B]` with whole-entity edits as a member
// =============================================================================

#[test]
fn rna_cis_spl_question_with_substitution_roundtrip() {
    let s = parse_ok("NM_004006.2:r.[100a>g;spl?]");
    assert!(s.contains("spl?"));
    assert!(s.contains("100a>g"));
}

#[test]
fn rna_cis_spl_then_substitution_roundtrip() {
    let s = parse_ok("NM_004006.2:r.[spl;100a>g]");
    assert!(s.contains("spl"));
    assert!(s.contains("100a>g"));
}

// =============================================================================
// Unknown-phase form: `r.[A(;)B]` with whole-entity edits as a member
// =============================================================================

#[test]
fn rna_unknown_phase_spl_with_substitution_roundtrip() {
    let s = parse_ok("NM_004006.2:r.[100a>g(;)spl?]");
    assert!(s.contains("spl?"));
    assert!(s.contains("(;)"), "expected unknown-phase '(;)' in {s}");
}

// =============================================================================
// Negative controls: invalid whole-entity-like content must still reject
// =============================================================================

#[test]
fn rna_bracket_garbage_still_rejected() {
    // `spline` is not `spl` or `spl?` — trailing chars must fail.
    assert!(parse_hgvs("NM_004006.2:r.[spline];[100a>g]").is_err());
}

#[test]
fn rna_bracket_uppercase_spl_rejected() {
    // RNA is lowercase per spec.
    assert!(parse_hgvs("NM_004006.2:r.[SPL];[100a>g]").is_err());
}
