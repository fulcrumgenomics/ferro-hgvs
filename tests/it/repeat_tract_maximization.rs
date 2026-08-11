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
//!
//! # #1618: reaching the maximal tract is only half of it
//!
//! Maximizing in one pass fixed idempotency and **broke the denoted sequence**.
//! `[6]` counts copies of the unit the caller spelled; widening from the 1-copy
//! `TG` tract to the 2-copy `GT` tract swallows a reference copy, so holding the
//! count fixed deletes two bases. Three of the four seam oracles pass on that
//! output — it is well-formed, in bounds, re-parses and is a fixed point — which
//! is why it survived; only `FERRO_ASSERT_SEQUENCE` sees it.
//!
//! The absorption the two-pass behaviour above performed by accident is now done
//! deliberately in one pass, so the count carries the copies the window gains:
//! `g.262TG[6]` -> `g.259_262GT[7]`, denoting the same 14 bases as its input. A
//! pure re-phase absorbs nothing and keeps its count, which is what
//! `DNA/repeated.md:97`'s 11 `TG` -> 11 `GT` requires.
//!
//! `tests/it/issue_1618_anchored_repeat_semantics.rs` holds the sequence-level
//! guards and the spec's own worked cases.

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

/// The canonical form the **anchored** spelling must reach.
///
/// The maximal tract is the 2-copy `GT` at g.259..262. Neither direction can
/// slide it — the flanking `A` (g.258) and `T` (g.263) break the `GT` phase on
/// both sides — so the window is identical under 3' and 5'. Asserting this
/// value, and not merely that a spelling is stable, is what distinguishes "found
/// the maximal tract" from "found the same NON-maximal tract" (the pre-fix
/// 1-copy `TG` window would satisfy a bare equality check).
///
/// **The count is 7, not the requested 6 (#1618).** `[6]` counts copies of the
/// unit the caller spelled, `TG`, whose tract here is a single copy; widening to
/// the 2-copy `GT` tract swallows one reference copy, and it has to be carried
/// or it is silently deleted. `g.262TG[6]` denotes `GT`×7 = 14 bases, and
/// `GT[7]` is that sequence on the maximal tract. `DNA/repeated.md:97-99` is the
/// spec's own instance of the same conservation: re-spelling 11 `TG` copies as
/// `GT` slides the tract one base 3' and the spec decrements the neighbouring
/// run `T[7]` -> `T[6]` to keep the total intact.
const CANONICAL: &str = "NC_TEST.1:g.259_262GT[7]";

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

/// Both spellings must land on the maximal tract — and must NOT converge,
/// because they do not denote the same variant.
///
/// **This assertion was inverted until #1618.** It required `g.262TG[6]` and
/// `g.259_262GT[6]` to normalize identically, which reads as confluence but is
/// not: `TG[6]` against the 1-copy `TG` tract denotes 14 bases, `GT[6]` against
/// the 2-copy `GT` tract denotes 12. Converging them is not two spellings of one
/// variant agreeing — it is one of them losing two bases. Confluence is a claim
/// about descriptions of the SAME variant, and these are descriptions of
/// different ones.
///
/// What both must share is the canonical *window*: the maximal 2-copy tract.
/// They differ only in the count each denotes.
#[test]
fn both_spellings_reach_the_maximal_tract_without_converging() {
    let anchored = norm("NC_TEST.1:g.262TG[6]", ShuffleDirection::ThreePrime);
    let ranged = norm("NC_TEST.1:g.259_262GT[6]", ShuffleDirection::ThreePrime);
    assert_eq!(anchored, CANONICAL, "anchored spelling");
    assert_eq!(
        ranged, "NC_TEST.1:g.259_262GT[6]",
        "the ranged spelling states its own tract, so it is already canonical"
    );
    assert!(
        anchored.starts_with("NC_TEST.1:g.259_262GT[")
            && ranged.starts_with("NC_TEST.1:g.259_262GT["),
        "both must land on the maximal 2-copy tract\n  anchored={anchored}\n  ranged  ={ranged}"
    );
    assert_ne!(
        anchored, ranged,
        "these denote 14 and 12 bases respectively; converging them destroys one"
    );
}
