//! Issue #468 — protein bracket-member dispatch for whole-entity edits.
//!
//! PR #433 (closes #423) extended whole-entity bracket-member dispatch
//! (`=`, `?`, `0`, splice markers) to the nucleic-acid axes
//! (`c.`/`n.`/`g.`/`m.`/`o.`), mirroring PR #396 item 3 for `r.`. The
//! protein axis (`p.`) was intentionally left out of #433 because `p.`
//! does not share the `NaEdit` type and `p.[=];[X]` needs a different
//! dummy-interval/edit shape (a `ProtInterval` at position `Met1` carrying
//! a `ProteinEdit` rather than a `NaEdit`).
//!
//! This file pins the protein-axis behavior. The spec-clear whole-entity
//! protein bracket member is `[=]` → whole-protein identity (`p.=`):
//!
//! - HGVS v21 `recommendations/protein/alleles.md` line 58-60:
//!   `NP_003997.1:p.[Ser68Arg];[=]` is a valid description and is
//!   explicitly *different* from `p.[Ser68Arg];[Ser68=]` — the `[=]`
//!   member means "the entire protein reference sequence was analysed".
//! - HGVS v21 `recommendations/protein/substitution.md` line 43:
//!   `p.=` means the **entire** protein coding region was analysed and no
//!   sequence-changing variant was found.
//!
//! For input tolerance (mirroring the established `p.(=)` / `p.(0)`
//! handling from PR #246 / #289 and PR #433's bracket probes), the
//! parenthesised predicted forms `[(=)]` and `[(?)]` are also accepted in
//! bracket-member position. The certain `[(=)]` Displays back as `[(=)]`.
//!
//! Pre-existing, unchanged behavior (NOT introduced here):
//!   - bare `[0]` / `[0?]` / `[(0)]` in a `p.`-prefixed trans-allele
//!     bracket → `ProteinEdit::NoProtein` (issue #277 / #289).
//!   - bare `[?]` in a `p.`-prefixed trans-allele bracket →
//!     `HgvsVariant::UnknownAllele` (cross-coord short-circuit).
//!   - bare `[?]` inside a cis bracket → position-unknown interval.

use ferro_hgvs::hgvs::edit::ProteinEdit;
use ferro_hgvs::hgvs::uncertainty::Mu;
use ferro_hgvs::hgvs::variant::{AllelePhase, ProteinVariant};
use ferro_hgvs::{parse_hgvs, HgvsVariant};

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

/// Extract the `(AllelePhase, &[HgvsVariant])` from a parsed allele.
fn allele_parts(variant: &HgvsVariant) -> (AllelePhase, &[HgvsVariant]) {
    match variant {
        HgvsVariant::Allele(a) => (a.phase, a.variants.as_slice()),
        other => panic!("expected HgvsVariant::Allele, got {other:?}"),
    }
}

/// Assert that the given allele arm is a protein variant carrying
/// `ProteinEdit::Identity { whole_protein: true, predicted }`.
fn assert_whole_protein_identity_arm(arm: &HgvsVariant, predicted: bool) {
    let pv: &ProteinVariant = match arm {
        HgvsVariant::Protein(pv) => pv,
        other => panic!("expected Protein arm, got {other:?}"),
    };
    let edit = match &pv.loc_edit.edit {
        Mu::Certain(e) => e,
        Mu::Uncertain(e) => e,
        other => panic!("expected Certain/Uncertain protein edit, got {other:?}"),
    };
    match edit {
        ProteinEdit::Identity {
            whole_protein: true,
            predicted: p,
        } => assert_eq!(
            *p, predicted,
            "predicted flag mismatch on whole-protein identity arm"
        ),
        other => panic!("expected whole-protein Identity, got {other:?}"),
    }
}

// =============================================================================
// Trans phase: `p.[X];[=]`
// =============================================================================

/// The canonical spec example (alleles.md:58-60): one allele carries a
/// substitution, the other carries whole-protein identity (`p.=`).
#[test]
fn trans_substitution_then_whole_protein_identity_roundtrips() {
    assert_canonical(
        "NP_003997.1:p.[Ser68Arg];[=]",
        "NP_003997.1:p.[Ser68Arg];[=]",
    );
}

/// Whole-protein identity on the *first* arm parses symmetrically.
#[test]
fn trans_whole_protein_identity_first_arm_roundtrips() {
    assert_canonical(
        "NP_003997.1:p.[=];[Ser68Arg]",
        "NP_003997.1:p.[=];[Ser68Arg]",
    );
}

/// Structural check: the `[=]` arm is a Protein variant carrying
/// whole-protein identity, NOT a `NullAllele`/`UnknownAllele`.
#[test]
fn trans_whole_protein_identity_arm_is_protein_identity() {
    let parsed = parse_hgvs("NP_003997.1:p.[Ser68Arg];[=]")
        .expect("parse trans allele with whole-protein identity");
    let (phase, arms) = allele_parts(&parsed);
    assert_eq!(phase, AllelePhase::Trans);
    assert_eq!(arms.len(), 2);
    assert!(
        matches!(&arms[0], HgvsVariant::Protein(_)),
        "first arm must be Protein, got {:?}",
        arms[0]
    );
    assert_whole_protein_identity_arm(&arms[1], false);
}

// =============================================================================
// Cis phase: `p.[X;=]`
// =============================================================================

/// Whole-protein identity as a trailing cis bracket member.
#[test]
fn cis_substitution_then_whole_protein_identity_roundtrips() {
    assert_canonical("NP_003997.1:p.[Ser68Arg;=]", "NP_003997.1:p.[Ser68Arg;=]");
}

/// Whole-protein identity as a leading cis bracket member.
#[test]
fn cis_whole_protein_identity_then_substitution_roundtrips() {
    assert_canonical("NP_003997.1:p.[=;Ser68Arg]", "NP_003997.1:p.[=;Ser68Arg]");
}

// =============================================================================
// Unknown phase: `p.[X(;)=]`
// =============================================================================

/// Whole-protein identity as an unknown-phase *bracket* member.
///
/// In-scope behavior for #468 is that the bracketed unknown-phase form
/// `p.[X(;)=]` parses symmetrically with the nucleic-acid axes and yields
/// an `AllelePhase::Unknown` allele whose second arm is whole-protein
/// identity.
///
/// We deliberately do NOT assert a Display round-trip fixed point here.
/// The protein unknown-phase Display canonicalises the bracketed form to
/// the bare compact form (`p.Ser68Arg(;)=`), and re-parsing that compact
/// form is a *pre-existing* protein-axis gap unrelated to whole-entity
/// members: even `p.[Ser68Arg(;)Asn594del]` → `p.Ser68Arg(;)Asn594del`
/// fails to re-parse on `main` (the bracketless unknown-phase protein
/// parser only admits fully parenthesised members `(X)(;)(Y)`). Fixing
/// that broad Display/parse asymmetry is out of scope for #468.
#[test]
fn unknown_phase_substitution_then_whole_protein_identity_parses() {
    let parsed = parse_hgvs("NP_003997.1:p.[Ser68Arg(;)=]")
        .expect("parse unknown-phase allele with whole-protein identity");
    let (phase, arms) = allele_parts(&parsed);
    assert_eq!(phase, AllelePhase::Unknown);
    assert_eq!(arms.len(), 2);
    assert!(
        matches!(&arms[0], HgvsVariant::Protein(_)),
        "first arm must be Protein, got {:?}",
        arms[0]
    );
    assert_whole_protein_identity_arm(&arms[1], false);
}

// =============================================================================
// Predicted-wrapper whole-entity members: `[(=)]`, `[(?)]`
// =============================================================================

/// Parenthesised predicted whole-protein identity in trans bracket
/// position. Tolerated input mirroring `p.(=)` (#246) and PR #433's
/// bracket probes. The certain Display of `(=)` keeps the parens.
#[test]
fn trans_predicted_whole_protein_identity_parses() {
    let parsed = parse_hgvs("NP_003997.1:p.[Ser68Arg];[(=)]")
        .expect("parse trans allele with predicted whole-protein identity");
    let (phase, arms) = allele_parts(&parsed);
    assert_eq!(phase, AllelePhase::Trans);
    assert_eq!(arms.len(), 2);
    assert_whole_protein_identity_arm(&arms[1], true);

    // Display is a re-parse fixed point.
    let displayed = parsed.to_string();
    let reparsed = parse_hgvs(&displayed)
        .unwrap_or_else(|e| panic!("re-parse of {displayed:?} must succeed; got {e:?}"));
    assert_eq!(reparsed.to_string(), displayed, "Display must be stable");
}

/// Parenthesised predicted whole-protein unknown (`p.(?)`) in trans
/// bracket position. Spec example alleles.md:55-56 uses bare `[?]` for
/// this; the explicit `[(?)]` is the tolerated parenthesised input.
#[test]
fn trans_predicted_whole_protein_unknown_parses() {
    let parsed = parse_hgvs("NP_003997.1:p.[Ser68Arg];[(?)]")
        .expect("parse trans allele with predicted whole-protein unknown");
    let (phase, arms) = allele_parts(&parsed);
    assert_eq!(phase, AllelePhase::Trans);
    assert_eq!(arms.len(), 2);
    let pv = match &arms[1] {
        HgvsVariant::Protein(pv) => pv,
        other => panic!("expected Protein arm, got {other:?}"),
    };
    let edit = match &pv.loc_edit.edit {
        Mu::Certain(e) | Mu::Uncertain(e) => e,
        other => panic!("expected protein edit, got {other:?}"),
    };
    assert!(
        matches!(
            edit,
            ProteinEdit::Unknown {
                whole_protein: true,
                predicted: true
            }
        ),
        "expected predicted whole-protein Unknown, got {edit:?}"
    );

    let displayed = parsed.to_string();
    let reparsed = parse_hgvs(&displayed)
        .unwrap_or_else(|e| panic!("re-parse of {displayed:?} must succeed; got {e:?}"));
    assert_eq!(reparsed.to_string(), displayed, "Display must be stable");
}

// =============================================================================
// Singleton bracket: `p.[=]`
// =============================================================================

/// A singleton bracketed whole-protein identity `p.[=]` unwraps to the
/// spec-sanctioned `p.=`, exactly as every other singleton bracket
/// unwraps its sole member (`p.[Ser68Arg]` → `p.Ser68Arg`). This is a
/// consequence of admitting whole-entity members in bracket position;
/// the prior `tests/protein_silent_eq.rs` pin that rejected `p.[=]` is
/// updated accordingly.
#[test]
fn singleton_bracket_whole_protein_identity_unwraps() {
    assert_canonical("NP_003997.1:p.[=]", "NP_003997.1:p.=");
}

// =============================================================================
// Pre-existing behavior preserved (regression guards)
// =============================================================================

/// `[0]` in a protein trans-allele bracket still routes to
/// `ProteinEdit::NoProtein` (issue #277), NOT changed by this work.
#[test]
fn trans_bracket_zero_still_no_protein() {
    let parsed =
        parse_hgvs("NP_003997.1:p.[Ser68Arg];[0]").expect("parse trans allele with bracketed zero");
    let (_, arms) = allele_parts(&parsed);
    let pv = match &arms[1] {
        HgvsVariant::Protein(pv) => pv,
        other => panic!("expected Protein arm for [0], got {other:?}"),
    };
    assert!(
        matches!(
            pv.loc_edit.edit,
            Mu::Certain(ProteinEdit::NoProtein { predicted: false })
        ),
        "[0] must remain NoProtein, got {:?}",
        pv.loc_edit.edit
    );
}

/// `[?]` in a protein trans-allele bracket still routes to
/// `HgvsVariant::UnknownAllele` (cross-coord short-circuit), unchanged.
#[test]
fn trans_bracket_question_still_unknown_allele() {
    let parsed = parse_hgvs("NP_003997.1:p.[Ser68Arg];[?]")
        .expect("parse trans allele with bracketed question");
    let (_, arms) = allele_parts(&parsed);
    assert!(
        matches!(&arms[1], HgvsVariant::UnknownAllele),
        "[?] must remain UnknownAllele, got {:?}",
        arms[1]
    );
}

/// Position-bound identity `[Ser68=]` is distinct from whole-protein
/// `[=]` and continues to parse as a position-bound identity (alleles.md
/// explicitly contrasts the two).
#[test]
fn trans_position_bound_identity_still_distinct() {
    assert_canonical(
        "NP_003997.1:p.[Ser68Arg];[Ser68=]",
        "NP_003997.1:p.[Ser68Arg];[Ser68=]",
    );
}

// =============================================================================
// Negative controls (spec-disallowed forms must still reject)
// =============================================================================

/// `p.[=;0]` mixes whole-protein identity with no-protein inside a single
/// cis bracket. `[0]` is not a valid cis bracket member (issue #277 pins
/// `p.[X;0]` as forbidden), so this must reject.
#[test]
fn cis_bracket_with_zero_member_rejects() {
    for input in [
        "NP_003997.1:p.[=;0]",
        "NP_003997.1:p.[0;=]",
        "NP_003997.1:p.[Ser68Arg;0]",
    ] {
        assert!(
            parse_hgvs(input).is_err(),
            "{input:?} must reject ([0] is not a valid cis bracket member)"
        );
    }
}

/// A whole-entity edit followed by trailing garbage in bracket position
/// must reject rather than silently consuming the prefix.
#[test]
fn whole_entity_member_with_trailing_garbage_rejects() {
    for input in [
        "NP_003997.1:p.[=X];[Ser68Arg]",
        "NP_003997.1:p.[Ser68Arg];[=foo]",
        "NP_003997.1:p.[Ser68Arg;=bar]",
    ] {
        assert!(
            parse_hgvs(input).is_err(),
            "{input:?} must reject (trailing chars after whole-entity edit)"
        );
    }
}
