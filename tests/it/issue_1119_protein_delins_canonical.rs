//! Reference-free canonicalization of a protein `delins` to its minimal form
//! (#1119). A `delins` that is structurally a smaller edit is rewritten per the
//! HGVS minimal-form / edit-priority rule (`protein/delins.md`), even when no
//! protein reference sequence is available — for the cases derivable from the
//! AST alone.
//!
//! # Scope
//!
//! Without a protein sequence, the deleted residues are known only where the
//! `delins` **location names them** — its start and end residues. That covers a
//! range spanning **≤ 2 residues** (every deleted residue is a named endpoint):
//!
//! - a 1-residue `delins` with a 1-residue insert is a **substitution**
//!   (`p.Val559delinsGly` → `p.Val559Gly`), regardless of shared affix
//!   (`delins.md`: sub > delins);
//! - a 2-residue `delins` whose inserted sequence shares a boundary residue
//!   with a named endpoint trims to a substitution or a pure deletion.
//!
//! A range spanning ≥ 3 residues has unnamed interior residues, so it cannot be
//! trimmed without the protein sequence and is left unchanged (the
//! reference-gated path in `try_protein_delins_canonicalize` handles it when a
//! protein-backed provider is present).
//!
//! These rewrites use an empty `MockProvider` (no protein data) to prove the
//! reference-free path.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

fn canonicalizes_to(input: &str, expected: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
    let normalized = Normalizer::new(MockProvider::new())
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("normalize {input:?}: {e}"));
    assert_eq!(normalized.to_string(), expected, "for {input:?}");

    // The canonical form must be a normalization fixed point: re-normalizing it
    // (reference-free) reproduces it exactly. This guards against a rewrite that
    // is not idempotent — e.g. one that keeps trimming or oscillates — across
    // every case the suite pins.
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

/// A 1→1 delins is a substitution even with no shared affix.
#[test]
fn single_residue_delins_becomes_substitution() {
    canonicalizes_to("NP_003997.1:p.Val559delinsGly", "NP_003997.1:p.Val559Gly");
}

/// A predicted `( )` delins is canonicalized inside its wrapper on the
/// reference-free path too: the minimized edit stays predicted, so the wrapper
/// is preserved rather than dropped or duplicated.
#[test]
fn predicted_delins_canonicalizes_inside_parens_reference_free() {
    canonicalizes_to(
        "NP_003997.1:p.(Val559delinsGly)",
        "NP_003997.1:p.(Val559Gly)",
    );
}

/// A 2-residue delins whose trailing inserted residue matches the named end
/// residue trims to a pure deletion of the remaining named residue.
#[test]
fn two_residue_delins_trailing_match_becomes_deletion() {
    canonicalizes_to(
        "NP_003997.1:p.Val559_Glu560delinsGlu",
        "NP_003997.1:p.Val559del",
    );
}

/// A 2-residue delins whose trailing inserted residue matches the named end
/// residue trims to a substitution at the remaining named residue.
#[test]
fn two_residue_delins_trailing_match_becomes_substitution() {
    canonicalizes_to(
        "NP_003997.1:p.Val559_Glu560delinsSerGlu",
        "NP_003997.1:p.Val559Ser",
    );
}

/// A 2-residue delins whose leading inserted residue matches the named start
/// residue trims to a substitution at the remaining named residue.
#[test]
fn two_residue_delins_leading_match_becomes_substitution() {
    canonicalizes_to(
        "NP_003997.1:p.Val559_Glu560delinsValSer",
        "NP_003997.1:p.Glu560Ser",
    );
}

/// A 2-residue delins whose only inserted residue matches the named start
/// residue trims the leading match away, leaving a pure deletion of the
/// remaining named residue (the leading-trim counterpart of the trailing-trim
/// deletion above).
#[test]
fn two_residue_delins_leading_match_becomes_deletion() {
    canonicalizes_to(
        "NP_003997.1:p.Val559_Glu560delinsVal",
        "NP_003997.1:p.Glu560del",
    );
}

/// A 2-residue delins with no shared boundary residue cannot be reduced
/// (reference-free) and is left unchanged.
#[test]
fn two_residue_delins_no_affix_is_unchanged() {
    canonicalizes_to(
        "NP_003997.1:p.Val559_Glu560delinsSerTrp",
        "NP_003997.1:p.Val559_Glu560delinsSerTrp",
    );
}

/// A ≥3-residue delins has unnamed interior residues, so without a protein
/// sequence it is left unchanged even when a boundary residue matches.
#[test]
fn three_residue_delins_is_unchanged_without_reference() {
    canonicalizes_to(
        "NP_003997.1:p.Val559_Trp561delinsAlaGlyTrp",
        "NP_003997.1:p.Val559_Trp561delinsAlaGlyTrp",
    );
}
