//! A single-position repeat anchor must find the *maximal* tandem tract, not
//! whichever phase the literal unit happens to match first.
//!
//! `normalize_repeat` tries `count_tandem_repeats` with the caller-spelled unit
//! and only falls back to the rotation/offset search when that returns `None`.
//! So when the literal unit matches a SHORT tract while a rotation of it matches
//! a LONGER tract covering the same anchor, the short one wins and the longer
//! one is never considered — the search is a fallback rather than a
//! maximization.
//!
//! That makes normalization non-idempotent: pass 1 emits a window covering only
//! the short tract, which pass 2 reads as an explicit *range*, discovers the
//! real (longer) tract around it, and absorbs the extra reference copies into
//! the count. Window and count both grow.
//!
//! Found by the projector-path idempotency oracle on a real spec input,
//! `NM_004006.3:r.-124ug[14]`, whose reference reads `…ACTGAAGTGTTACTTT…`:
//! literal `TG` matches 1 copy at one offset while `GT` matches 2 copies
//! spanning the same anchor.
//!
//! ```text
//! r.-124ug[14]      -> r.-125_-124gu[14]    (1-copy TG tract)
//! r.-125_-124gu[14] -> r.-127_-124gu[15]    (2-copy GT tract + absorption)
//! ```
//!
//! The core below reproduces that shape on the genomic axis: `GTGT` at
//! g.259..262 is 2 copies of `GT`, and the anchor g.262 is the tract's last
//! base, where literal `TG` tiles only 1 copy.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

/// `AA` + `GTGT` (g.259..262) + `TA` — the `GT` tract is 2 copies, cleanly
/// bounded (`A` before, `T` after breaks the `GT` phase on both sides).
const CORE: &str = "AAGTGTTA";

fn norm(input: &str, dir: ShuffleDirection) -> String {
    Normalizer::with_config(
        SyntheticBuilder::genomic(CORE).build(),
        NormalizeConfig::default().with_direction(dir),
    )
    .normalize(&parse_hgvs(input).expect("parse"))
    .expect("normalize")
    .to_string()
}

/// The one canonical form every spelling and direction must reach.
///
/// The maximal tract is the 2-copy `GT` at g.259..262. Neither direction can
/// slide it — the flanking `A` (g.258) and `T` (g.263) break the `GT` phase on
/// both sides — so the window is identical under 3' and 5', and the count is the
/// requested 6. Asserting this value, and not merely that the spellings agree,
/// is what distinguishes "both found the maximal tract" from "both found the
/// same NON-maximal tract" (the pre-fix 1-copy `TG` window would satisfy a bare
/// equality check).
const CANONICAL: &str = "NC_TEST.1:g.259_262GT[6]";

/// An expansion spelled on a single anchor whose literal unit tiles the shorter
/// phase must still be a fixed point — and must land on the maximal tract.
#[test]
fn single_anchor_repeat_expansion_is_idempotent_three_prime() {
    let once = norm("NC_TEST.1:g.262TG[6]", ShuffleDirection::ThreePrime);
    assert_eq!(once, CANONICAL, "must maximize to the 2-copy GT tract");
    let twice = norm(&once, ShuffleDirection::ThreePrime);
    assert_eq!(
        once, twice,
        "NOT IDEMPOTENT\n  once ={once}\n  twice={twice}"
    );
}

#[test]
fn single_anchor_repeat_expansion_is_idempotent_five_prime() {
    let once = norm("NC_TEST.1:g.262TG[6]", ShuffleDirection::FivePrime);
    assert_eq!(
        once, CANONICAL,
        "the tract cannot slide 5' (g.258 is `A`), so 5' matches 3'"
    );
    let twice = norm(&once, ShuffleDirection::FivePrime);
    assert_eq!(
        once, twice,
        "NOT IDEMPOTENT\n  once ={once}\n  twice={twice}"
    );
}

/// The anchored and range spellings of the same tract must agree, which is the
/// property whose violation produced the drift: the single anchor found 1 copy
/// where the explicit range found 2.
#[test]
fn anchored_and_range_spellings_agree_on_the_maximal_tract() {
    let anchored = norm("NC_TEST.1:g.262TG[6]", ShuffleDirection::ThreePrime);
    let ranged = norm("NC_TEST.1:g.259_262GT[6]", ShuffleDirection::ThreePrime);
    assert_eq!(anchored, CANONICAL, "anchored spelling");
    assert_eq!(ranged, CANONICAL, "ranged spelling");
    assert_eq!(
        anchored, ranged,
        "a single anchor and an explicit range over the same tract must \
         normalize identically\n  anchored={anchored}\n  ranged  ={ranged}"
    );
}
