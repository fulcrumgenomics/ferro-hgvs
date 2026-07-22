//! Reference-free minimization of a protein `delins` whose affix-trim leaves a
//! **pure insertion between two named endpoints** — #1126, extending #1119.
//!
//! # Scope
//!
//! `try_protein_delins_canonicalize` trims a ≤2-residue `delins` against the
//! residues its own location names. When the trim consumes the whole deleted
//! window the residual is a pure insertion, which previously required a
//! protein-backed provider (the ins→dup search needs the reference). But when
//! the trim leaves a named residue on *both* sides of the insertion point, the
//! minimal form is derivable from the AST alone:
//!
//! `p.Val559_Glu560delinsValSerGlu` — leading `Val` matches `Val559`, trailing
//! `Glu` matches `Glu560` — is `p.Val559_Glu560insSer`.
//!
//! The ins→dup refinement stays reference-gated: with a protein-backed provider
//! the same input may reduce further to a `dup`, which only the reference can
//! establish. Reference-free we emit the minimal `ins`.
//!
//! When either flank is **unnamed** — the trim ran off the start of the named
//! window, or consumed it entirely — the insertion point cannot be anchored
//! without the reference, and the input `delins` is kept.
//!
//! These rewrites use an empty `MockProvider` (no protein data) to prove the
//! reference-free path.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Pin the canonical reference-free form of `input`, and require that form to
/// be a normalization fixed point (idempotent) — a `delins`→`ins` rewrite must
/// not then be re-trimmed or oscillate on a second pass.
fn canonicalizes_to(input: &str, expected: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
    let normalized = Normalizer::new(MockProvider::new())
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("normalize {input:?}: {e}"));
    assert_eq!(normalized.to_string(), expected, "for {input:?}");

    let reparsed =
        parse_hgvs(expected).unwrap_or_else(|e| panic!("parse expected {expected:?}: {e}"));
    let renormalized = Normalizer::new(MockProvider::new())
        .normalize(&reparsed)
        .unwrap_or_else(|e| panic!("re-normalize {expected:?}: {e}"));
    assert_eq!(
        renormalized.to_string(),
        expected,
        "canonical form {expected:?} is not a normalization fixed point",
    );
}

/// The #1126 repro: both flanks are named, so the residual single-residue
/// insertion is expressible without a reference.
#[test]
fn both_flanks_named_trims_to_insertion() {
    canonicalizes_to(
        "NP_003997.1:p.Val559_Glu560delinsValSerGlu",
        "NP_003997.1:p.Val559_Glu560insSer",
    );
}

/// A multi-residue residual insertion trims the same way — the whole residual
/// becomes the inserted sequence.
#[test]
fn multi_residue_residual_trims_to_insertion() {
    canonicalizes_to(
        "NP_003997.1:p.Val559_Glu560delinsValSerTrpGlu",
        "NP_003997.1:p.Val559_Glu560insSerTrp",
    );
}

/// The predicted `( )` wrapper survives the rewrite.
#[test]
fn predicted_delins_trims_to_predicted_insertion() {
    canonicalizes_to(
        "NP_003997.1:p.(Val559_Glu560delinsValSerGlu)",
        "NP_003997.1:p.(Val559_Glu560insSer)",
    );
}

/// No leading match: the trim comes entirely off the trailing end, so the
/// insertion point sits 5′ of the named window and its left flank residue is
/// unknown. Reference-free this cannot be anchored — keep the input.
#[test]
fn unnamed_left_flank_keeps_the_delins() {
    canonicalizes_to(
        "NP_003997.1:p.Val559_Glu560delinsSerValGlu",
        "NP_003997.1:p.Val559_Glu560delinsSerValGlu",
    );
}

/// The leading match consumes the ENTIRE named window, so the insertion point
/// sits 3′ of it and its right flank residue is unknown — keep the input.
#[test]
fn unnamed_right_flank_keeps_the_delins() {
    canonicalizes_to(
        "NP_003997.1:p.Val559_Val560delinsValValVal",
        "NP_003997.1:p.Val559_Val560delinsValValVal",
    );
}

/// A 1-residue `delins` can never leave a pure insertion between two *named*
/// flanks — one side is always outside the single named residue — so both
/// trim directions keep the input.
#[test]
fn single_residue_delins_keeps_the_delins_in_both_trim_directions() {
    canonicalizes_to(
        "NP_003997.1:p.Val559delinsValSer",
        "NP_003997.1:p.Val559delinsValSer",
    );
    canonicalizes_to(
        "NP_003997.1:p.Val559delinsSerVal",
        "NP_003997.1:p.Val559delinsSerVal",
    );
}
