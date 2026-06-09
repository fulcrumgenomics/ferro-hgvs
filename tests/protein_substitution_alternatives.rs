//! Residue-level and/or `^` on a protein substitution (#544).
//!
//! `p.(Gly56Ala^Ser^Cys)` means residue Gly56 is changed to Ala, Ser, OR Cys
//! — one reference + position, multiple alternative residues
//! (protein/substitution.md).
//!
//! This form shares its outer `p.(...^...)` shape with the whole-edit and/or
//! operator from #547 (`p.(Gly23GlufsTer7^Gly23CysfsTer26)`, uncertain.md).
//! The two are disambiguated by the operands: residue-alternatives have bare
//! amino-acid trailing operands (no position), whole-edit and/or operands are
//! full edits (each carries a position). These tests pin that boundary so the
//! two `^` readings stay distinct.

use ferro_hgvs::hgvs::variant::AllelePhase;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

#[test]
fn protein_substitution_alternatives_round_trip() {
    for s in [
        "NP_003997.1:p.(Gly56Ala^Ser^Cys)",
        "NP_003997.1:p.(Gly719Ala^Ser)",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        assert!(
            matches!(v, HgvsVariant::Protein(_)),
            "expected Protein (residue-alternatives) for `{s}`, got {v:?}"
        );
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}

/// The whole-edit and/or form (trailing operands are full edits carrying a
/// position) must still route to #547's and/or allele, NOT the residue path.
#[test]
fn whole_edit_and_or_still_routes_to_allele() {
    for s in [
        "NP_003997.1:p.(Gly23GlufsTer7^Gly23CysfsTer26)",
        // Two full substitutions at the same position — positioned trailing
        // operand, so this is and/or, not residue-alternatives.
        "NP_003997.1:p.(Gly56Ala^Gly56Ser)",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        match &v {
            HgvsVariant::Allele(a) => {
                assert_eq!(
                    a.phase,
                    AllelePhase::AndOr,
                    "expected and/or phase for `{s}`"
                )
            }
            other => panic!("expected and/or Allele for `{s}`, got {other:?}"),
        }
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}

/// The compact-mosaic somatic form (`<acc>:p.<pos>=/<bare-rhs>`) accepts a
/// residue-alternatives RHS (`Ala^Ser`), but ONLY when the inherited identity
/// LHS is a point. A range identity LHS would build a `SubstitutionAlternatives`
/// spanning the whole interval, which breaks the single-residue contract, so it
/// must be rejected rather than silently mis-parsed.
#[test]
fn compact_mosaic_residue_alternatives_rejects_range_lhs() {
    // Point LHS: residue-alternatives RHS is accepted (mosaic somatic allele).
    let point = "NP_003997.1:p.Trp24=/Ala^Ser";
    let v = parse_hgvs(point).unwrap_or_else(|e| panic!("point LHS must parse `{point}`: {e}"));
    match &v {
        HgvsVariant::Allele(a) => {
            assert_eq!(
                a.phase,
                AllelePhase::Mosaic,
                "expected mosaic for `{point}`"
            );
            assert_eq!(a.variants.len(), 2, "expected 2 members for `{point}`");
        }
        other => panic!("expected mosaic Allele for `{point}`, got {other:?}"),
    }

    // Range LHS: must NOT produce a residue-alternatives edit on a range.
    let range = "NP_003997.1:p.Trp24_Cys26=/Ala^Ser";
    assert!(
        parse_hgvs(range).is_err(),
        "range identity LHS must be rejected for residue-alternatives RHS, got {:?}",
        parse_hgvs(range)
    );
}

/// Fail-safe pin: a trailing whole-protein token (`=`) is NOT an amino acid,
/// so it is not mistaken for a residue alternative — the input routes to the
/// and/or path rather than being swallowed (and mis-rejected) by the residue
/// shorthand.
#[test]
fn trailing_non_residue_token_is_not_residue_alternatives() {
    let s = "NP_003997.1:p.(Trp24Cys^=)";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert!(
        matches!(v, HgvsVariant::Allele(_)),
        "`{s}` must route to the and/or allele path, got {v:?}"
    );
}
