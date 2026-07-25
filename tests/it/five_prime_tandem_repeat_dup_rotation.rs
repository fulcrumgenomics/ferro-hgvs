//! Regression tests for the 5'-direction tandem-repeat `dup` idempotency
//! (idempotency-campaign follow-up to PR #1169/#1170).
//!
//! A duplication over an ambiguous alternating tandem tract routes to `unit[N]`
//! repeat notation. The `dup`→`unit[N]` emission is anchored at the shuffle
//! stage's direction-appropriate position (5'-most under `FivePrime`), but the
//! repeat *re-normalization* (`normalize_repeat`) previously always applied the
//! 3' phase alignment (`three_prime_align_tract`) regardless of direction — so a
//! `FivePrime` `GA[4]` anchored at the tract's 5' phase drifted one base 3' to
//! `AG[4]` and rotated its unit on the second pass, breaking idempotency.
//!
//! The fix makes `normalize_repeat` direction-aware: `FivePrime` uses the mirror
//! `five_prime_align_tract`, so the 5'-most repeat is a fixed point AND both
//! spellings of the haplotype converge to it. `ThreePrime` is unchanged.
//!
//! Reference core `GGAGAGGGC` places an ambiguous `GAGAG` tract at g.258..262
//! (as `GA` it is 2 clean copies at g.258; as `AG` it is 2 clean copies at
//! g.259). `g.258_261dup` extends it to 4 copies either way — the same
//! haplotype, so the two representations `g.258_261GA[4]` and `g.259_262AG[4]`
//! describe identical sequence and differ only in canonical anchor/phase.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

const CORE: &str = "GGAGAGGGC";

fn norm(input: &str, dir: ShuffleDirection) -> String {
    Normalizer::with_config(
        SyntheticBuilder::genomic(CORE).build(),
        NormalizeConfig::default().with_direction(dir),
    )
    .normalize(&parse_hgvs(input).expect("parse"))
    .expect("normalize")
    .to_string()
}

/// 5' normalizes the dup to the 5'-most repeat phase `GA[4]` at g.258.
#[test]
fn five_prime_tandem_dup_anchors_at_five_prime_phase() {
    assert_eq!(
        norm("NC_TEST.1:g.258_261dup", ShuffleDirection::FivePrime),
        "NC_TEST.1:g.258_261GA[4]",
    );
}

/// That 5' repeat is a fixed point — it previously drifted to `g.259_262AG[4]`.
#[test]
fn five_prime_tandem_dup_is_idempotent() {
    let once = norm("NC_TEST.1:g.258_261dup", ShuffleDirection::FivePrime);
    assert_eq!(
        norm(&once, ShuffleDirection::FivePrime),
        once,
        "norm(norm(x)) must equal norm(x)"
    );
}

/// Both spellings of the one haplotype converge to the 5'-most `GA[4]` under 5'
/// (idempotency alone would not catch a per-spelling fixed point).
#[test]
fn five_prime_tandem_repeat_spellings_are_confluent() {
    for input in ["NC_TEST.1:g.258_261GA[4]", "NC_TEST.1:g.259_262AG[4]"] {
        assert_eq!(
            norm(input, ShuffleDirection::FivePrime),
            "NC_TEST.1:g.258_261GA[4]",
            "5' confluence for {input}",
        );
    }
}

/// 3' is unchanged: it phases to the 3'-most `AG[4]` at g.259 and is idempotent.
#[test]
fn three_prime_tandem_dup_unchanged_and_idempotent() {
    let once = norm("NC_TEST.1:g.258_261dup", ShuffleDirection::ThreePrime);
    assert_eq!(once, "NC_TEST.1:g.259_262AG[4]");
    assert_eq!(
        norm(&once, ShuffleDirection::ThreePrime),
        once,
        "3' must remain idempotent"
    );
}
