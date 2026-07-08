//! Predicted-wrapper and whole-entity-LHS mosaic `=/` forms (issue #955).
//!
//! Two gaps beyond the already-supported position-bound compact mosaic
//! (`p.Trp24=/Cys`, `r.85=/u>c`):
//!   - Gap A: the predicted `(...)` wrapper — `p.(Trp24=/Cys)` (protein
//!     substitution/deletion/duplication), `r.(=/6_8del)` (RNA). The
//!     slash-splitters are `[]`-aware but not `()`-aware, so the wrapped `/`
//!     was mis-routed.
//!   - Gap B: a whole-entity `=` LHS with a positioned RHS — `r.=/6_8del`.

use ferro_hgvs::hgvs::AllelePhase;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

fn assert_round_trips(input: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("failed to parse {input:?}: {e}"));
    assert_eq!(v.to_string(), input, "round-trip mismatch for {input:?}");
}

/// Independent structural check (not just round-trip): confirm a parsed
/// predicted `(...)` mosaic/chimeric built the expected AST — an `Allele` with
/// the right phase, the predicted-uncertain wrapper flag set, and the expected
/// member count — verifying the parse itself is correct, not merely that
/// Display inverts it.
fn assert_predicted_allele(input: &str, phase: AllelePhase, members: usize) {
    match parse_hgvs(input).unwrap_or_else(|e| panic!("failed to parse {input:?}: {e}")) {
        HgvsVariant::Allele(a) => {
            assert_eq!(a.phase, phase, "phase for {input:?}");
            assert!(
                a.uncertain,
                "predicted `(...)` wrapper must set the uncertain flag for {input:?}"
            );
            assert_eq!(a.variants.len(), members, "member count for {input:?}");
        }
        other => panic!("{input:?} must parse as an Allele, got {other:?}"),
    }
}

#[test]
fn predicted_mosaic_parses_to_expected_ast() {
    assert_predicted_allele("LRG_199t1:p.(Trp24=/Cys)", AllelePhase::Mosaic, 2);
    assert_predicted_allele("LRG_199t1:r.(=/6_8del)", AllelePhase::Mosaic, 2);
    assert_predicted_allele("NM_004006.3:r.(85=//u>c)", AllelePhase::Chimeric, 2);
}

#[test]
fn predicted_protein_mosaic_round_trips() {
    // Gap A — protein substitution / deletion / duplication under `p.(...)`.
    assert_round_trips("LRG_199t1:p.(Trp24=/Cys)");
    assert_round_trips("NP_003997.1:p.(Val7=/del)");
    assert_round_trips("NP_003997.1:p.(Val7=/dup)");
}

#[test]
fn whole_entity_lhs_rna_mosaic_round_trips() {
    // Gap B — whole-entity `r.=` LHS with a positioned RHS.
    assert_round_trips("LRG_199t1:r.=/6_8del");
}

#[test]
fn predicted_rna_whole_entity_mosaic_round_trips() {
    // Gap A + Gap B — predicted wrapper around a whole-entity-LHS RNA mosaic.
    assert_round_trips("LRG_199t1:r.(=/6_8del)");
}

#[test]
fn existing_compact_mosaic_still_round_trips() {
    // Regression: position-bound compact mosaic (already supported).
    assert_round_trips("LRG_199t1:p.Trp24=/Cys");
    assert_round_trips("NM_004006.3:r.85=/u>c");
}

#[test]
fn predicted_chimeric_round_trips() {
    // Gap A for the chimeric `//` separator (was untested): the predicted
    // wrapper around a 2-member chimeric substitution must round-trip.
    assert_round_trips("NM_004006.3:r.(85=//u>c)");
}

#[test]
fn predicted_wrapper_rejects_more_than_two_members() {
    // The predicted `(...)` wrapper is only defined for a 2-member mosaic. A
    // 3-member predicted wrapper must be rejected, not silently accepted and
    // then rendered WITHOUT the parens (which would change observed-vs-predicted
    // meaning).
    for bad in [
        "LRG_199t1:r.(=/6_8del/10_12del)",
        "LRG_199t1:p.(Trp24=/Cys/Ser)",
    ] {
        assert!(
            parse_hgvs(bad).is_err(),
            ">2-member predicted mosaic wrapper must be rejected: {bad:?}"
        );
    }
}
