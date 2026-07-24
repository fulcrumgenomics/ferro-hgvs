//! Regression test for a 5'-direction insertion bug found by the #1157
//! idempotency proptest campaign.
//!
//! When a multi-base literal insertion is 5'-shifted through a repeat region,
//! the inserted sequence must ROTATE as it crosses each repeat base — inserting
//! `CA` after position p equals inserting `AC` after position p-1 in an A-run.
//! The emit-path rotation used `result.start.saturating_sub(shuffle_start)`,
//! which clamps to 0 for a leftward (5') shift, so the alt was emitted UNROTATED
//! at the shifted position — a DIFFERENT haplotype than the input, and
//! non-idempotent (it slid one more base on each subsequent normalize). The fix
//! makes the rotation direction-aware.
//!
//! NB: this pins only the one shape the fix corrects. The 5' direction has other
//! (separate, unfixed) non-idempotency / non-confluence bugs in the shuffle +
//! boundary-clamp interaction; those are a dedicated follow-up and are NOT
//! covered here.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

fn norm5(core: &str, input: &str) -> String {
    let normalizer = Normalizer::with_config(
        SyntheticBuilder::genomic(core).build(),
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    normalizer
        .normalize(&parse_hgvs(input).expect("parse"))
        .expect("normalize")
        .to_string()
}

/// Core is 8 A's at g.257..=264. `g.260_261insCA` 5'-shifts left one base; the
/// alt must rotate `CA -> AC`, giving `g.259_260insAC` — a fixed point. The bug
/// emitted the unrotated `g.259_260insCA` (a different sequence) that then kept
/// sliding.
#[test]
fn five_prime_multibase_insertion_rotates_and_is_idempotent() {
    assert_eq!(
        norm5("AAAAAAAA", "NC_TEST.1:g.260_261insCA"),
        "NC_TEST.1:g.259_260insAC",
    );
}

#[test]
fn five_prime_rotated_insertion_is_a_fixed_point() {
    let once = norm5("AAAAAAAA", "NC_TEST.1:g.260_261insCA");
    assert_eq!(
        norm5("AAAAAAAA", &once),
        once,
        "norm(norm(x)) must equal norm(x)"
    );
}
