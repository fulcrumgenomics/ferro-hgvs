//! Issue #300 — Predicted-edit wrapper Display parity between `c.` and `r.`.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/300>.
//!
//! Background (carried over from PR for #234, which pinned `c.`/`r.` Display
//! parity, and the #283 audit of `c.` → `r.` conversion):
//!
//! Both `c.(123A>G)` and `r.(123a>g)` parse as predicted-edit wrappers. Per
//! HGVS v21.0 §G3 (`uncertain.md`), the canonical Display form wraps the
//! **whole** position+edit expression in parens (i.e. `c.(123A>G)` /
//! `r.(123a>g)`), not just the edit body and not the bare edit form.
//!
//! Historically the two coordinate systems diverged:
//!   - `NM_000088.3:c.(123A>G)` → `NM_000088.3:c.123(A>G)` (parens migrated)
//!   - `NM_000088.3:r.(123a>g)` → `NM_000088.3:r.123a>g`   (parens dropped)
//!
//! The asymmetry was resolved by PR #242
//! ([`242a21d`](../../commit/242a21d)), which replaced the c.-only
//! substitution gate with `delimited('(', (parse_interval, parse_na_edit),
//! ')')` on every nucleic-acid axis and unified the Display path to emit
//! `<acc>:<axis>.(<pos><edit>)`. PR #244 extended that to compound
//! brackets; PR #246 added the whole-entity `(=)/(?)/(0)` forms.
//!
//! This file pins the canonical outer-paren round-trip shape on both axes
//! across the full edit-kind matrix so the invariant cannot regress under
//! a future refactor. Coverage is symmetrical because the parser fix made
//! the two axes share their dispatch: any edit kind that parses on one
//! must parse on the other with the corresponding alphabet (uppercase DNA
//! vs. lowercase RNA).

use ferro_hgvs::parse_hgvs;

const ACC: &str = "NM_000088.3";

/// Parse `input`, format the result, and require the output equals `expected`.
/// Then re-parse the output and assert idempotency (round-trip stability).
fn assert_canonical(input: &str, expected: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse({input:?}) failed: {e}"));
    let out = format!("{v}");
    assert_eq!(
        out, expected,
        "Display mismatch for {input:?}: expected {expected:?}, got {out:?}"
    );
    let v2 =
        parse_hgvs(&out).unwrap_or_else(|e| panic!("re-parse of canonical {out:?} failed: {e}"));
    let out2 = format!("{v2}");
    assert_eq!(
        out2, expected,
        "idempotency broken for {expected:?}: second pass became {out2:?}"
    );
}

// ============================================================================
// c.  — full edit-kind matrix.
// ============================================================================

#[test]
fn cds_predicted_substitution_preserves_outer_parens() {
    assert_canonical(&format!("{ACC}:c.(123A>G)"), &format!("{ACC}:c.(123A>G)"));
}

#[test]
fn cds_predicted_deletion_preserves_outer_parens() {
    assert_canonical(
        &format!("{ACC}:c.(123_125del)"),
        &format!("{ACC}:c.(123_125del)"),
    );
}

#[test]
fn cds_predicted_insertion_preserves_outer_parens() {
    assert_canonical(
        &format!("{ACC}:c.(123_124insATG)"),
        &format!("{ACC}:c.(123_124insATG)"),
    );
}

#[test]
fn cds_predicted_duplication_preserves_outer_parens() {
    assert_canonical(
        &format!("{ACC}:c.(123_125dup)"),
        &format!("{ACC}:c.(123_125dup)"),
    );
}

#[test]
fn cds_predicted_inversion_preserves_outer_parens() {
    assert_canonical(
        &format!("{ACC}:c.(123_125inv)"),
        &format!("{ACC}:c.(123_125inv)"),
    );
}

#[test]
fn cds_predicted_delins_preserves_outer_parens() {
    assert_canonical(
        &format!("{ACC}:c.(123_125delinsATG)"),
        &format!("{ACC}:c.(123_125delinsATG)"),
    );
}

#[test]
fn cds_predicted_identity_preserves_outer_parens() {
    assert_canonical(&format!("{ACC}:c.(123A=)"), &format!("{ACC}:c.(123A=)"));
}

#[test]
fn cds_predicted_repeat_preserves_outer_parens() {
    assert_canonical(
        &format!("{ACC}:c.(123_125AC[4])"),
        &format!("{ACC}:c.(123_125AC[4])"),
    );
}

// ============================================================================
// r.  — full edit-kind matrix.
// ============================================================================

#[test]
fn rna_predicted_substitution_preserves_outer_parens() {
    assert_canonical(&format!("{ACC}:r.(123a>g)"), &format!("{ACC}:r.(123a>g)"));
}

#[test]
fn rna_predicted_deletion_preserves_outer_parens() {
    assert_canonical(
        &format!("{ACC}:r.(123_125del)"),
        &format!("{ACC}:r.(123_125del)"),
    );
}

#[test]
fn rna_predicted_insertion_preserves_outer_parens() {
    assert_canonical(
        &format!("{ACC}:r.(123_124insaug)"),
        &format!("{ACC}:r.(123_124insaug)"),
    );
}

#[test]
fn rna_predicted_duplication_preserves_outer_parens() {
    assert_canonical(
        &format!("{ACC}:r.(123_125dup)"),
        &format!("{ACC}:r.(123_125dup)"),
    );
}

#[test]
fn rna_predicted_inversion_preserves_outer_parens() {
    assert_canonical(
        &format!("{ACC}:r.(123_125inv)"),
        &format!("{ACC}:r.(123_125inv)"),
    );
}

#[test]
fn rna_predicted_delins_preserves_outer_parens() {
    assert_canonical(
        &format!("{ACC}:r.(123_125delinsaug)"),
        &format!("{ACC}:r.(123_125delinsaug)"),
    );
}

#[test]
fn rna_predicted_identity_preserves_outer_parens() {
    assert_canonical(&format!("{ACC}:r.(123a=)"), &format!("{ACC}:r.(123a=)"));
}

#[test]
fn rna_predicted_repeat_preserves_outer_parens() {
    assert_canonical(
        &format!("{ACC}:r.(123_125ac[4])"),
        &format!("{ACC}:r.(123_125ac[4])"),
    );
}

// ============================================================================
// Compound brackets: PR #244 extended the wrapper to predicted edits inside
// a cis allele. Pin the outer-paren shape on the wrapped member.
// ============================================================================

#[test]
fn cds_predicted_substitution_in_compound_bracket_preserves_parens() {
    assert_canonical(
        &format!("{ACC}:c.[(123A>G);125C>T]"),
        &format!("{ACC}:c.[(123A>G);125C>T]"),
    );
}

#[test]
fn rna_predicted_substitution_in_compound_bracket_preserves_parens() {
    assert_canonical(
        &format!("{ACC}:r.[(123a>g);125c>u]"),
        &format!("{ACC}:r.[(123a>g);125c>u]"),
    );
}

// ============================================================================
// c./r. structural parity (PR #234 invariant — alphabet differs, structure
// does not). Cover every edit kind to lock in cross-axis symmetry.
// ============================================================================

/// Assert that `cds_input` and `rna_input` Display to the same shape modulo
/// alphabet (case-fold + DNA-to-RNA `t` → `u`) and axis marker (`:c.` →
/// `:r.`).
///
/// Normalisation removes only the spec-defined axis differences between
/// DNA and RNA descriptions: uppercase-vs-lowercase and `T` vs `U`. Any
/// other structural divergence between the two Display shapes — wrapper
/// position, separator, paren placement — will fail the assertion.
fn assert_axis_parity(cds_input: &str, rna_input: &str) {
    fn normalize(s: &str) -> String {
        s.replace(":c.", ":r.").to_lowercase().replace('t', "u")
    }
    let cds_out = format!("{}", parse_hgvs(cds_input).unwrap());
    let rna_out = format!("{}", parse_hgvs(rna_input).unwrap());
    let normalized_cds = normalize(&cds_out);
    let normalized_rna = normalize(&rna_out);
    assert_eq!(
        normalized_cds, normalized_rna,
        "c./r. predicted-wrapper Display shapes diverge: {cds_out:?} vs {rna_out:?}"
    );
}

#[test]
fn cds_and_rna_predicted_substitution_share_structure() {
    assert_axis_parity(&format!("{ACC}:c.(123A>G)"), &format!("{ACC}:r.(123a>g)"));
}

#[test]
fn cds_and_rna_predicted_deletion_share_structure() {
    assert_axis_parity(
        &format!("{ACC}:c.(123_125del)"),
        &format!("{ACC}:r.(123_125del)"),
    );
}

#[test]
fn cds_and_rna_predicted_insertion_share_structure() {
    assert_axis_parity(
        &format!("{ACC}:c.(123_124insATG)"),
        &format!("{ACC}:r.(123_124insaug)"),
    );
}

#[test]
fn cds_and_rna_predicted_duplication_share_structure() {
    assert_axis_parity(
        &format!("{ACC}:c.(123_125dup)"),
        &format!("{ACC}:r.(123_125dup)"),
    );
}

#[test]
fn cds_and_rna_predicted_inversion_share_structure() {
    assert_axis_parity(
        &format!("{ACC}:c.(123_125inv)"),
        &format!("{ACC}:r.(123_125inv)"),
    );
}

#[test]
fn cds_and_rna_predicted_delins_share_structure() {
    assert_axis_parity(
        &format!("{ACC}:c.(123_125delinsATG)"),
        &format!("{ACC}:r.(123_125delinsaug)"),
    );
}

#[test]
fn cds_and_rna_predicted_identity_share_structure() {
    assert_axis_parity(&format!("{ACC}:c.(123A=)"), &format!("{ACC}:r.(123a=)"));
}

#[test]
fn cds_and_rna_predicted_repeat_share_structure() {
    assert_axis_parity(
        &format!("{ACC}:c.(123_125AC[4])"),
        &format!("{ACC}:r.(123_125ac[4])"),
    );
}

#[test]
fn cds_and_rna_predicted_substitution_in_compound_bracket_share_structure() {
    assert_axis_parity(
        &format!("{ACC}:c.[(123A>G);125C>T]"),
        &format!("{ACC}:r.[(123a>g);125c>u]"),
    );
}

// ============================================================================
// Negative: the spec-illegal inner-paren shape `<axis>.<pos>(<edit>)` (parens
// around the edit only, not the whole position+edit expression) must be
// rejected by the parser on both axes. Pinning this on `c.` and `r.`
// prevents a future lenient parser from silently accepting the non-canonical
// input on either dispatch path.
// ============================================================================

#[test]
fn cds_inner_paren_only_is_rejected() {
    let result = parse_hgvs(&format!("{ACC}:c.123(A>G)"));
    assert!(
        result.is_err(),
        "spec-illegal `c.<pos>(<edit>)` form must be rejected, got {result:?}"
    );
}

#[test]
fn rna_inner_paren_only_is_rejected() {
    let result = parse_hgvs(&format!("{ACC}:r.123(a>g)"));
    assert!(
        result.is_err(),
        "spec-illegal `r.<pos>(<edit>)` form must be rejected, got {result:?}"
    );
}
