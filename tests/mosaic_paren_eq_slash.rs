//! Parenthesized (predicted) and whole-entity mosaic `=/` forms (#544).
//!
//! - `p.(Val7=/del)` — a predicted (`(...)`) mosaic `=/`.
//! - `r.=/6_8del` and `r.(=/6_8del)` — mosaic with a whole-entity `r.=` LHS
//!   and a position-bearing RHS.
//!
//! Builds on the protein `=/` (#548) and the `uncertain` wrapper (#547).

use ferro_hgvs::hgvs::AllelePhase;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

#[test]
fn parenthesized_predicted_mosaic_eq_slash_round_trips() {
    for s in [
        "LRG_199t1:p.(Trp24=/Cys)",
        "NP_003997.1:p.(Val7=/del)",
        "NP_003997.1:p.(Val7=/dup)",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        let HgvsVariant::Allele(a) = &v else {
            panic!("expected Allele for `{s}`, got {v:?}");
        };
        assert_eq!(a.phase, AllelePhase::Mosaic, "phase for `{s}`");
        assert!(a.uncertain, "uncertain flag for `{s}`");
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}

#[test]
fn whole_entity_rna_mosaic_eq_slash_round_trips() {
    for s in ["LRG_199t1:r.=/6_8del", "LRG_199t1:r.(=/6_8del)"] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}

/// A whole-entity identity on *both* sides of the mosaic (`=/=`), including
/// the predicted (`(=/=)`) wrapper, must round-trip. The predicted form is the
/// load-bearing case: the Display "Branch 1b" path is the only branch that
/// emits the `(...)` wrapper for a whole-entity-identity RHS, so it must keep
/// handling this shape rather than deferring to the `var/=` cleanup branch
/// (which drops the `uncertain` parens).
#[test]
fn whole_entity_identity_both_sides_mosaic_round_trips() {
    for s in [
        "NC_000001.11:g.=/=",
        "NC_000001.11:g.(=/=)",
        "NM_004006.2:c.=/=",
        "NM_004006.2:c.(=/=)",
        "LRG_199t1:r.=/=",
        "LRG_199t1:r.(=/=)",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}

/// A malformed RHS following a whole-entity identity LHS must surface the
/// RHS's *own* parse diagnostic, not the generic "Unknown variant type prefix"
/// fallback. `g.=/6_8delQ` rebuilds to `…:g.6_8delQ`, whose genuine failure is
/// the trailing `Q`; the old `.ok()` swallowed it and reported a misleading
/// "found '6_8'" prefix error instead.
#[test]
fn whole_entity_mosaic_malformed_rhs_surfaces_specific_error() {
    let bad = "NC_000001.11:g.=/6_8delQ";
    let err = parse_hgvs(bad).expect_err("malformed RHS must fail to parse");
    let msg = err.to_string();
    assert!(
        msg.contains("trailing") && msg.contains('Q'),
        "expected the specific trailing-`Q` diagnostic for `{bad}`, got: {msg}"
    );
    assert!(
        !msg.contains("Unknown variant type prefix"),
        "must not fall through to the generic prefix error for `{bad}`, got: {msg}"
    );

    // Other malformed whole-entity-mosaic RHS shapes must also be rejected.
    for bad in ["NC_000001.11:g.=/6_8xyz", "NM_004006.2:c.=/6_8delZ"] {
        assert!(
            parse_hgvs(bad).is_err(),
            "malformed RHS `{bad}` must fail to parse"
        );
    }
}
