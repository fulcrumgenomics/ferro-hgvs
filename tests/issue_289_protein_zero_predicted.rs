//! Audit for issue #289 — `p.(0)` whole-entity predicted protein null-allele.
//!
//! Follow-up to PR #246 (#245), which added `c.(=)`, `c.(?)`, `r.(=)`, `r.(?)`,
//! and `r.(0)` predicted whole-entity forms via the `Mu::Uncertain` wrapper.
//! The protein analogue `p.(0)` was flagged out-of-scope there because the
//! protein code path is structurally different: `ProteinEdit::NoProtein`
//! carries its own `predicted: bool` flag rather than relying on an outer
//! `Mu::Uncertain` wrapper.
//!
//! On main, `p.(=)` and `p.(?)` already parse correctly (routed through
//! `whole_protein_identity_predicted()` / `whole_protein_unknown_predicted()`),
//! but `p.(0)` fails with a parse error because the dispatcher only accepts
//! the bare `0` / `0?` forms.
//!
//! This audit pins:
//!
//! 1. `p.(0)` parses (alongside the spec-canonical `p.0?`) and routes to
//!    `ProteinEdit::NoProtein { predicted: true }`.
//! 2. Parse → Display → parse → Display is idempotent (after at most one
//!    canonicalization pass: the input form `p.(0)` Displays as the
//!    spec-canonical `p.0?`, which then round-trips).
//! 3. LRG accessions (`LRG_<n>p<m>`) accept `p.(0)` identically.
//! 4. Adjacent forms continue to work: `p.0` (certain no-protein), `p.0?`
//!    (predicted no-protein, spec-canonical Display), `p.(=)` (predicted
//!    identity), `p.(?)` (predicted unknown).
//!
//! Display canonicalization choice: the HGVS spec
//! (`assets/hgvs-nomenclature/docs/recommendations/protein/`, cited at
//! `tests/protein_no_protein_roundtrip.rs:1-9`) uses `p.0?` and `p.0` for
//! the predicted / certain no-protein forms — never the parenthesised
//! `p.(0)`. We therefore *accept* `p.(0)` as a tolerated alternate input
//! (mirroring the RNA `r.(0)` parse path's intent) but continue to emit
//! `p.0?` as the canonical Display. This keeps the long-standing D6
//! contract pinned by `protein_no_protein_roundtrip.rs` intact.
//!
//! Asymmetry note (resist a "make all three symmetric" refactor): `p.(=)`
//! and `p.(?)` Display preserves the parens (spec form for predicted
//! identity / unknown), but `p.(0)` canonicalizes to `p.0?` (spec form for
//! predicted no-protein). The three predicted whole-entity forms share a
//! parser but diverge on Display per the protein spec — keep them that way.

use ferro_hgvs::hgvs::edit::ProteinEdit;
use ferro_hgvs::hgvs::variant::ProteinVariant;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

/// Borrow the inner `ProteinVariant` from an `HgvsVariant`.
fn protein_variant(v: &HgvsVariant) -> &ProteinVariant {
    match v {
        HgvsVariant::Protein(p) => p,
        other => panic!("expected Protein variant, got {:?}", other),
    }
}

/// Borrow the inner `ProteinEdit` (panics if no concrete inner edit).
fn protein_edit(v: &HgvsVariant) -> &ProteinEdit {
    protein_variant(v)
        .loc_edit
        .edit
        .inner()
        .expect("p.(0) / p.0? / p.0 must carry a concrete (non-uncertain) edit")
}

// =============================================================================
// SECTION 1 — `p.(0)` parses and routes to NoProtein { predicted: true }
// =============================================================================

/// `p.(0)` must parse as a protein variant on a standard `NP_` accession.
#[test]
fn p_paren_zero_parses_on_np_accession() {
    let v = parse_hgvs("NP_000088.3:p.(0)")
        .expect("NP_000088.3:p.(0) must parse — issue #289 whole-entity predicted form");
    assert!(matches!(v, HgvsVariant::Protein(_)));
}

/// The parsed `p.(0)` must route to `ProteinEdit::NoProtein { predicted: true }`,
/// not to a positional edit, not to a different ProteinEdit variant, and not
/// silently collapse to the certain form.
#[test]
fn p_paren_zero_routes_to_no_protein_predicted_true() {
    let v = parse_hgvs("NP_000088.3:p.(0)").expect("parse p.(0)");
    match protein_edit(&v) {
        ProteinEdit::NoProtein { predicted } => assert!(
            *predicted,
            "p.(0) must carry predicted=true; got NoProtein {{ predicted: false }}"
        ),
        other => panic!(
            "p.(0) did not parse as ProteinEdit::NoProtein; got {:?}",
            other
        ),
    }
}

/// `p.(0)` and the spec-canonical `p.0?` must be structurally equal: same
/// accession, same loc_edit, same `predicted` flag. A consumer pattern-matching
/// on `ProteinEdit::NoProtein { predicted: true }` must accept both inputs
/// identically.
#[test]
fn p_paren_zero_is_structurally_equal_to_p_zero_question() {
    let paren = parse_hgvs("NP_000088.3:p.(0)").expect("parse p.(0)");
    let canonical = parse_hgvs("NP_000088.3:p.0?").expect("parse p.0?");
    assert_eq!(
        paren, canonical,
        "p.(0) and p.0? must parse to structurally equal variants"
    );
}

// =============================================================================
// SECTION 2 — Round-trip idempotency
// =============================================================================

/// Parse `p.(0)`, Display, then re-parse, re-Display: the canonical Display
/// form is `p.0?` (spec). The first Display already canonicalizes (`p.(0)` →
/// `p.0?`); subsequent round-trips are fixed-points.
#[test]
fn p_paren_zero_round_trip_is_idempotent_after_canonicalization() {
    let parsed = parse_hgvs("NP_000088.3:p.(0)").expect("parse p.(0)");
    let displayed = format!("{}", parsed);

    // Display canonicalizes to the spec form.
    assert_eq!(
        displayed, "NP_000088.3:p.0?",
        "p.(0) Display must canonicalize to spec form p.0?"
    );

    // Re-parsing the canonical form yields a structurally equal variant.
    let reparsed = parse_hgvs(&displayed).expect("re-parse canonical Display");
    assert_eq!(
        parsed, reparsed,
        "re-parse of canonical Display must equal the original parse"
    );

    // And re-Display is a fixed point.
    assert_eq!(
        format!("{}", reparsed),
        displayed,
        "Display is a fixed point on the canonical form"
    );
}

// =============================================================================
// SECTION 3 — LRG accession path
// =============================================================================

/// `LRG_1p1:p.(0)` must parse just like the `NP_` form. LRG protein accessions
/// use the `LRG_<n>p<m>` shape (see `protein_no_protein_roundtrip.rs::CASES`
/// for `LRG_199p1:p.0?`).
#[test]
fn p_paren_zero_parses_on_lrg_accession() {
    let v = parse_hgvs("LRG_1p1:p.(0)").expect("LRG_1p1:p.(0) must parse");
    match protein_edit(&v) {
        ProteinEdit::NoProtein { predicted } => assert!(*predicted),
        other => panic!("LRG p.(0) did not route to NoProtein, got {:?}", other),
    }
    // Same canonicalization on the LRG path.
    assert_eq!(format!("{}", v), "LRG_1p1:p.0?");
}

// =============================================================================
// SECTION 4 — Regression guards on adjacent whole-protein forms
// =============================================================================

/// `p.0` (certain no-protein) must continue to parse and Display as `p.0`.
#[test]
fn certain_p_zero_still_works() {
    let v = parse_hgvs("NP_000088.3:p.0").expect("parse p.0");
    match protein_edit(&v) {
        ProteinEdit::NoProtein { predicted } => {
            assert!(!*predicted, "p.0 must carry predicted=false")
        }
        other => panic!("p.0 did not route to NoProtein, got {:?}", other),
    }
    assert_eq!(format!("{}", v), "NP_000088.3:p.0");
}

/// `p.0?` (predicted no-protein, spec form) must continue to parse and
/// Display as `p.0?` — this is the canonical Display form.
#[test]
fn predicted_p_zero_question_still_works_and_is_canonical_display() {
    let v = parse_hgvs("NP_000088.3:p.0?").expect("parse p.0?");
    match protein_edit(&v) {
        ProteinEdit::NoProtein { predicted } => assert!(*predicted),
        other => panic!("p.0? did not route to NoProtein, got {:?}", other),
    }
    assert_eq!(format!("{}", v), "NP_000088.3:p.0?");
}

/// `p.(=)` (predicted whole-protein identity) must continue to round-trip.
/// Pinned here as a boundary guard so a refactor of the `(0)` arm cannot
/// inadvertently break the `(=)` arm.
#[test]
fn predicted_identity_paren_still_round_trips() {
    let v = parse_hgvs("NP_000088.3:p.(=)").expect("parse p.(=)");
    assert!(matches!(v, HgvsVariant::Protein(_)));
    assert_eq!(format!("{}", v), "NP_000088.3:p.(=)");
}

/// `p.(?)` (predicted whole-protein unknown) must continue to round-trip.
/// Pinned here as a boundary guard for the same reason as `p.(=)`.
#[test]
fn predicted_unknown_paren_still_round_trips() {
    let v = parse_hgvs("NP_000088.3:p.(?)").expect("parse p.(?)");
    assert!(matches!(v, HgvsVariant::Protein(_)));
    assert_eq!(format!("{}", v), "NP_000088.3:p.(?)");
}

/// Bare-form predicted edits must continue to work: `p.=`, `p.?`.
#[test]
fn bare_predicted_edits_still_round_trip() {
    for (input, expected) in &[
        ("NP_000088.3:p.=", "NP_000088.3:p.="),
        ("NP_000088.3:p.?", "NP_000088.3:p.?"),
    ] {
        let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e}"));
        assert_eq!(format!("{}", v), *expected);
    }
}

// =============================================================================
// SECTION 5 — Negative guards: spec-malformed variants stay rejected
// =============================================================================

/// `p.(0?)` (double-marked predicted: inner `?` plus outer parens) is not in
/// the spec and must continue to reject. Pinned in
/// `protein_no_protein_roundtrip.rs::malformed_adjacent_inputs_reject`; we
/// re-pin it here as a focused boundary guard for issue #289 so a future
/// expansion of the `(0)` arm cannot loosen acceptance to `(0?)`.
#[test]
fn double_marked_p_paren_zero_question_rejects() {
    assert!(
        parse_hgvs("NP_000088.3:p.(0?)").is_err(),
        "p.(0?) is not a spec form and must reject"
    );
}

/// Trailing junk after `p.(0)` must reject. The `(0)` arm must consume
/// exactly the three characters `(0)` and nothing more.
#[test]
fn p_paren_zero_with_trailing_junk_rejects() {
    for bad in &[
        "NP_000088.3:p.(0)X",
        "NP_000088.3:p.(0)a",
        "NP_000088.3:p.(0)0",
        "NP_000088.3:p.(0)?",
    ] {
        assert!(
            parse_hgvs(bad).is_err(),
            "{bad:?} must reject (trailing junk after p.(0))"
        );
    }
}

// =============================================================================
// SECTION 6 — Trans-allele compact-form (neighborhood of #277)
// =============================================================================

/// `p.[(0)];[Arg97Trp]` parses and Displays correctly. The `(0)` arm is the
/// tolerated alternate input for predicted no-protein inside a compact-form
/// trans-allele bracket; it routes through
/// `parse_protein_trans_allele_shorthand` to the same
/// `ProteinEdit::NoProtein { predicted: true }` value as `[0?]`. Display
/// canonicalises the `(0)` arm to `[0?]` (spec form), mirroring the bare
/// dispatcher's `p.(0)` → `p.0?` canonicalisation. This sits in the
/// neighborhood of #277's trans-allele `[0]` / `[0?]` routing work.
#[test]
fn p_paren_zero_in_trans_allele_parses_and_displays() {
    use ferro_hgvs::hgvs::variant::AllelePhase;

    let parsed = parse_hgvs("NP_000088.3:p.[(0)];[Arg97Trp]")
        .expect("NP_000088.3:p.[(0)];[Arg97Trp] must parse — issue #289 trans-allele form");

    // Outer shape: a trans Allele with two arms.
    let allele = match &parsed {
        HgvsVariant::Allele(a) => a,
        other => panic!("expected Allele, got {:?}", other),
    };
    assert_eq!(allele.phase, AllelePhase::Trans);
    assert_eq!(allele.variants.len(), 2);

    // First arm: predicted no-protein.
    match protein_edit(&allele.variants[0]) {
        ProteinEdit::NoProtein { predicted } => assert!(
            *predicted,
            "trans-allele [(0)] arm must carry predicted=true"
        ),
        other => panic!(
            "trans-allele [(0)] arm did not route to NoProtein; got {:?}",
            other
        ),
    }

    // Second arm: ordinary substitution — still parses as Protein.
    assert!(
        matches!(&allele.variants[1], HgvsVariant::Protein(_)),
        "second arm must be Protein(Arg97Trp), got {:?}",
        allele.variants[1]
    );

    // Display canonicalises the `(0)` arm to `[0?]` (spec form), matching the
    // bare-dispatcher canonicalisation pinned in
    // `p_paren_zero_round_trip_is_idempotent_after_canonicalization`.
    let displayed = format!("{}", parsed);
    assert_eq!(
        displayed, "NP_000088.3:p.[0?];[Arg97Trp]",
        "trans-allele Display must canonicalise [(0)] to [0?]"
    );

    // Re-parse → Display is a fixed point on the canonical form.
    let reparsed = parse_hgvs(&displayed).expect("re-parse canonical trans Display");
    assert_eq!(
        parsed, reparsed,
        "re-parse of canonical trans Display must equal the original parse"
    );
    assert_eq!(
        format!("{}", reparsed),
        displayed,
        "Display is a fixed point on the canonical trans form"
    );
}

/// Regression for issue #424: the canonical spec form `p.[0?];[X]` must route
/// to `ProteinEdit::NoProtein { predicted: true }`. This is the first-arm
/// path through `parse_protein_trans_allele_shorthand` (the arm that handles
/// both `[0]` and `[0?]`, distinguishing them by `predicted`). The
/// `[(0)]` form is covered by `p_paren_zero_in_trans_allele_parses_and_displays`
/// above and routes through the third arm.
///
/// Pinning both arms as separate tests guards against a future refactor that
/// collapses them and loses one of the input shapes — the third arm used to
/// carry a dead `|| content == "0?"` clause that masked exactly this kind of
/// loss.
#[test]
fn p_zero_question_in_trans_allele_parses_and_displays() {
    use ferro_hgvs::hgvs::variant::AllelePhase;

    let parsed = parse_hgvs("NP_000088.3:p.[0?];[Arg97Trp]")
        .expect("NP_000088.3:p.[0?];[Arg97Trp] must parse — issue #424 first-arm regression");

    let allele = match &parsed {
        HgvsVariant::Allele(a) => a,
        other => panic!("expected Allele, got {:?}", other),
    };
    assert_eq!(allele.phase, AllelePhase::Trans);
    assert_eq!(allele.variants.len(), 2);

    match protein_edit(&allele.variants[0]) {
        ProteinEdit::NoProtein { predicted } => assert!(
            *predicted,
            "trans-allele [0?] arm must carry predicted=true"
        ),
        other => panic!(
            "trans-allele [0?] arm did not route to NoProtein; got {:?}",
            other
        ),
    }

    // Spec form Displays back to itself byte-for-byte.
    assert_eq!(format!("{}", parsed), "NP_000088.3:p.[0?];[Arg97Trp]");
}

/// Negative guard for issue #424: with the third arm now tightened to
/// `content == "(0)"` (exact match), the double-marked shape `p.[(0?)];[X]`
/// must reject — `(0?)` matches neither the first arm (exact `0` / `0?`),
/// the second arm (`?`), nor the third arm (`(0)`), and falls through to
/// the generic `parse_prot_interval` path which rejects the parens. The
/// bare-form `p.(0?)` rejection is already pinned at
/// `double_marked_p_paren_zero_question_rejects` above; this is the
/// trans-allele bracket sibling.
#[test]
fn p_paren_zero_question_in_trans_allele_rejects() {
    for bad in &[
        "NP_000088.3:p.[(0?)];[Arg97Trp]",
        "NP_000088.3:p.[Arg97Trp];[(0?)]",
    ] {
        assert!(
            parse_hgvs(bad).is_err(),
            "{bad:?} must reject (double-marked predicted in trans-allele bracket is not a spec form)"
        );
    }
}
