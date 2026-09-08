//! The R2 winner family as library constructors.
//!
//! The four arms that tied at 119/120 on the blessed corpus and were validated
//! at scale (0-invalid, ~93–98% agreement with production) share one L2 typer
//! (`OperatorExtractL2::v2`) and one Dial-B config (`DialBConfig::ledger_r5`),
//! differing only in their L1 segmenter. They previously existed only as
//! string-keyed closures inside `examples/bakeoff_production_pick.rs`
//! (`l1s()`/`l2s()`/`build_all_arms()`); this module lifts them to `pub fn`s so
//! the confluence / idempotency harness (`tests/it/bakeoff_winner_confluence.rs`)
//! and any driver can target the exact winner without reconstructing it, the same
//! way `sweep::winning_arm`/`candidate_arm` are library constructors.
//!
//! Zero new logic: each function is a direct transcription of the corresponding
//! `l1s()`/`l2s()` cell. If those cells change, these must be kept in step —
//! `tests/it/bakeoff_winner_confluence.rs::the_r2_winner_names_are_the_graded_arms`
//! pins the names against `build_all_arms()`.

use crate::partition::arm::DialBConfig;
use crate::partition::arms::anchor_peel::AnchorPeelL1;
use crate::partition::arms::arms::{AxisAware, TrimOnlyL1, WallPolicy};
use crate::partition::arms::op_extract::OperatorExtractL2;
use crate::partition::strategy::Strategy;

/// `trim-only/op-extract-v2/ledger-r5` — the best ceiling-fidelity L1 (bare
/// prefix/suffix trim).
pub fn trim_only() -> Strategy {
    Strategy::Grid {
        name: "trim-only/op-extract-v2/ledger-r5".into(),
        l1: Box::new(TrimOnlyL1),
        l2: Box::new(OperatorExtractL2::v2(4, 2, "op-extract-v2")),
        dial_b: DialBConfig::ledger_r5(),
    }
}

/// `axis-peel/op-extract-v2/ledger-r5` — the primary/TSV-witness arm: anchor-peel
/// L1 under the DNA coincidence-collapse `AxisAware` wrapper.
pub fn axis_peel() -> Strategy {
    Strategy::Grid {
        name: "axis-peel/op-extract-v2/ledger-r5".into(),
        l1: Box::new(AxisAware::new(
            Box::new(AnchorPeelL1::new(0.01, "peel-strict")),
            WallPolicy::all_on(),
        )),
        l2: Box::new(OperatorExtractL2::v2(4, 2, "op-extract-v2")),
        dial_b: DialBConfig::ledger_r5(),
    }
}

/// `peel-comp/op-extract-v2/ledger-r5` — composition-aware anchor-peel L1.
pub fn peel_comp() -> Strategy {
    Strategy::Grid {
        name: "peel-comp/op-extract-v2/ledger-r5".into(),
        l1: Box::new(AnchorPeelL1::new_composition(0.01, "peel-comp")),
        l2: Box::new(OperatorExtractL2::v2(4, 2, "op-extract-v2")),
        dial_b: DialBConfig::ledger_r5(),
    }
}

/// `peel-strict/op-extract-v2/ledger-r5` — strict-floor anchor-peel L1, unwrapped.
pub fn peel_strict() -> Strategy {
    Strategy::Grid {
        name: "peel-strict/op-extract-v2/ledger-r5".into(),
        l1: Box::new(AnchorPeelL1::new(0.01, "peel-strict")),
        l2: Box::new(OperatorExtractL2::v2(4, 2, "op-extract-v2")),
        dial_b: DialBConfig::ledger_r5(),
    }
}

/// All four R2 winner arms, in `(name, strategy)` pairs.
pub fn all() -> Vec<Strategy> {
    vec![trim_only(), axis_peel(), peel_comp(), peel_strict()]
}
