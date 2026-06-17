//! Issue #423 — `parse_cds_bracket_member`, `parse_tx_bracket_member`,
//! and `parse_genome_bracket_member` only knew about position-based
//! edits and the predicted-wrapper `(pos edit)`. Inside a c./n./g./m./o.
//! allele bracket (`c.[X];[Y]` or `c.[X;Y]` or `c.[X(;)Y]`), whole-entity
//! edits like `=`, `?` (and `0` for `n.`) were rejected because they
//! have no positional content.
//!
//! This issue closes the symmetry gap with PR #396 item 3 (which added
//! the same dispatch for `r.`). Each dispatcher now probes the
//! axis-specific whole-entity parsers (`parse_whole_cds_*`,
//! `parse_whole_rna_*`, `parse_whole_genome_*`, plus `parse_rna_no_product`
//! for the `n.0` tx form) before the position-based parsers, attaching
//! a dummy interval at position 1 (mirroring the top-level dispatchers'
//! handling of the same forms).
//!
//! Bare `[0]` / `[?]` in the trans-allele form `[X];[Y]` continue to
//! cross-coord-short-circuit to `NullAllele` / `UnknownAllele` via the
//! shared `parse_trans_allele_shorthand_generic`. Only cis/unknown-flat
//! paths and predicted-form members (`[(?)]`, `[(=)]`, `[(0)]`) plus
//! bare whole-entity tokens inside cis brackets flow into the new
//! probes.
//!
//! **Display-shape note:** for some inputs (e.g. `c.[100A>G(;)=]`,
//! `c.[(?)];[100A>G]`) the existing pre-#423 Display logic re-emits the
//! variant in the expanded `[ACC:c.X];[ACC:c.Y]` form rather than the
//! compact-prefix `c.[X];[Y]` form. That re-emission is pre-existing
//! behavior on this axis pair, not introduced here, so the assertions
//! below pin the actual canonical Display output (re-parse-stable) for
//! each input rather than expecting a byte-identical echo. Re-parse of
//! every canonical Display is asserted as a fixed point at the end of
//! each test.

use ferro_hgvs::parse_hgvs;

/// Parse, Display, and re-parse the input. Asserts:
/// 1. parse(`input`) succeeds and Displays to `expected_canonical`.
/// 2. parse(`expected_canonical`) succeeds.
/// 3. Display of the re-parsed variant equals `expected_canonical`
///    (Display is a fixed point on the canonical form).
fn assert_canonical(input: &str, expected_canonical: &str) {
    let parsed =
        parse_hgvs(input).unwrap_or_else(|e| panic!("parse({input:?}) must succeed; got {e:?}"));
    let displayed = parsed.to_string();
    assert_eq!(
        displayed, expected_canonical,
        "Display({input:?}) must equal {expected_canonical:?}; got {displayed:?}",
    );
    let reparsed = parse_hgvs(&displayed).unwrap_or_else(|e| {
        panic!("re-parse of canonical Display {displayed:?} must succeed; got {e:?}")
    });
    let redisplayed = reparsed.to_string();
    assert_eq!(
        redisplayed, expected_canonical,
        "Display must be a fixed point on the canonical form for {input:?}; got {redisplayed:?}",
    );
}

// =============================================================================
// CDS (c.) axis
// =============================================================================

#[test]
fn cds_cis_substitution_then_identity_roundtrip() {
    assert_canonical("NM_000088.3:c.[100A>G;=]", "NM_000088.3:c.[100A>G;=]");
}

#[test]
fn cds_unknown_phase_substitution_with_identity_roundtrip() {
    // Unknown-phase brackets (`(;)`) get re-emitted in the bracketless
    // form `c.X(;)Y` per the pre-#423 Display logic. Re-parse stable.
    assert_canonical("NM_000088.3:c.[100A>G(;)=]", "NM_000088.3:c.100A>G(;)=");
}

#[test]
fn cds_trans_identity_then_substitution_roundtrip() {
    // `c.[=];[100A>G]` — whole-CDS identity on one arm, substitution
    // on the other. The `[=]` arm flows through the bracket-member
    // dispatcher (NOT the cross-coord short-circuit, which only fires
    // for bare `[0]` / `[?]`).
    assert_canonical("NM_000088.3:c.[=];[100A>G]", "NM_000088.3:c.[=];[100A>G]");
}

#[test]
fn cds_trans_predicted_unknown_with_substitution_roundtrip() {
    // `c.[(?)];[100A>G]` — predicted-uncertain whole-CDS unknown on
    // one arm. Display re-emits in expanded form per pre-#423 behavior.
    assert_canonical(
        "NM_000088.3:c.[(?)];[100A>G]",
        "[NM_000088.3:c.(?)];[NM_000088.3:c.100A>G]",
    );
}

#[test]
fn cds_cis_unknown_with_substitution_roundtrip() {
    // `c.[100A>G;?]` — bare `?` inside a cis bracket flows through
    // the bracket member (not short-circuited because the
    // short-circuit only fires in the trans `[X];[Y]` shape).
    // Display re-emits in expanded form per pre-#423 behavior.
    assert_canonical(
        "NM_000088.3:c.[100A>G;?]",
        "[NM_000088.3:c.100A>G;NM_000088.3:c.?]",
    );
}

// =============================================================================
// Non-coding tx (n.) axis — includes the `0` / `(0)` no-product forms
// =============================================================================

#[test]
fn tx_cis_substitution_then_identity_roundtrip() {
    assert_canonical("NR_001234.1:n.[100A>G;=]", "NR_001234.1:n.[100A>G;=]");
}

#[test]
fn tx_trans_identity_with_substitution_roundtrip() {
    assert_canonical("NR_001234.1:n.[=];[100A>G]", "NR_001234.1:n.[=];[100A>G]");
}

#[test]
fn tx_trans_predicted_no_product_with_substitution_roundtrip() {
    // `n.[(0)];[100A>G]` — predicted no transcript on one arm.
    // Bare `[0]` short-circuits to `NullAllele`; `[(0)]` is the
    // tx-specific predicted no-product form (mirrors `r.(0)`).
    assert_canonical(
        "NR_001234.1:n.[(0)];[100A>G]",
        "NR_001234.1:n.[(0)];[100A>G]",
    );
}

#[test]
fn tx_trans_pure_whole_entity_pair_roundtrip() {
    // Both arms whole-entity: `n.[=];[?]`. Bare `[?]` cross-coord-short-
    // circuits to `UnknownAllele` upstream of the bracket member, so
    // this exercises the new `=` probe on the first arm AND the
    // upstream short-circuit on the second.
    let v = parse_hgvs("NR_001234.1:n.[=];[?]")
        .expect("NR_001234.1:n.[=];[?] must parse — pure whole-entity trans pair");
    // The displayed form depends on how UnknownAllele Displays inside
    // an Allele — pin re-parse stability rather than the exact bytes.
    let displayed = v.to_string();
    let reparsed = parse_hgvs(&displayed)
        .unwrap_or_else(|e| panic!("re-parse of {displayed:?} must succeed; got {e:?}"));
    assert_eq!(
        reparsed.to_string(),
        displayed,
        "Display must be a fixed point for the canonical form",
    );
}

// =============================================================================
// Genome (g./m./o.) axis
// =============================================================================

#[test]
fn genome_cis_substitution_then_identity_roundtrip() {
    assert_canonical("NC_000001.11:g.[100A>G;=]", "NC_000001.11:g.[100A>G;=]");
}

#[test]
fn genome_trans_identity_with_substitution_roundtrip() {
    assert_canonical("NC_000001.11:g.[=];[100A>G]", "NC_000001.11:g.[=];[100A>G]");
}

#[test]
fn mito_cis_substitution_then_identity_roundtrip() {
    assert_canonical("NC_012920.1:m.[100A>G;=]", "NC_012920.1:m.[100A>G;=]");
}

#[test]
fn mito_trans_predicted_unknown_with_substitution_roundtrip() {
    // `m.[(?)];[100A>G]` — predicted-uncertain whole-mito unknown on
    // one arm. Same shape as `g.[(?)]`, routed through the same
    // dispatcher with `kind = Mt`. Display expansion mirrors c.[(?)]
    // (pre-#423 trans-allele Display behavior).
    assert_canonical(
        "NC_012920.1:m.[(?)];[100A>G]",
        "[NC_012920.1:m.(?)];[NC_012920.1:m.100A>G]",
    );
}

// =============================================================================
// Protein bracket-member dispatch is intentionally out of scope here.
// `p.[=];[X]` and the protein cis-form whole-entity members
// (`p.[Arg97Cys;=]`, etc.) are not currently a supported protein
// bracket-member shape — the issue body's claim that the protein
// dispatcher already accepts `[=]` is true only for the
// position-bearing form `[pos=]` (e.g. `[Arg97=];[X]`), which has
// always worked. Closing the "whole-protein identity in a bracket" gap
// is its own follow-up because protein needs a different dummy-interval
// shape (no `1`-base dummy makes sense for `p.`) and a different
// round-trip story.
//
// =============================================================================
// Negative controls: spec-malformed inputs must reject
// =============================================================================

#[test]
fn cds_bracket_trailing_junk_after_identity_rejected() {
    // `=garbage` is not a whole-entity edit; trailing chars must fail
    // the `final_remaining.trim().is_empty()` check.
    assert!(parse_hgvs("NM_000088.3:c.[=garbage];[100A>G]").is_err());
}

#[test]
fn cds_bracket_double_identity_rejected() {
    // `==` is not a spec form. The first `=` is consumed by the
    // whole-entity-identity probe, leaving `=` as the remainder —
    // which fails the bracket member's
    // `final_remaining.trim().is_empty()` guard.
    assert!(parse_hgvs("NM_000088.3:c.[==];[100A>G]").is_err());
}

#[test]
fn cds_bracket_identity_then_unknown_concat_rejected() {
    // `=?` is not a spec form — same shape as `==` but exercising the
    // "first probe wins, second token leaks" failure path.
    assert!(parse_hgvs("NM_000088.3:c.[=?];[100A>G]").is_err());
}

#[test]
fn tx_bracket_double_marked_no_product_rejected() {
    // `(0?)` is the spec-malformed double-marked form for `n.` no
    // product (the protein-axis analog `p.(0?)` is pinned in
    // `tests/issue_289_protein_zero_predicted.rs::double_marked_p_paren_zero_question_rejects`).
    // Must reject in trans-allele context too.
    assert!(parse_hgvs("NR_001234.1:n.[(0?)];[100A>G]").is_err());
}

#[test]
fn genome_bracket_uppercase_position_letter_rejected() {
    // `X` is not a base or whole-entity marker; the position parser
    // refuses it and the bracket member returns an error.
    assert!(parse_hgvs("NC_000001.11:g.[100A>G;X]").is_err());
}
