//! A `c.`-axis cis allele must not normalize to an insertion spanning the
//! CDS/3'UTR junction (#1398).
//!
//! ferro's design forbids re-axing a CDS-interior edit onto the `c.*N` axis: the
//! CDS-end clamp turns a saturating insertion into `c.<cds_end>delins…` rather
//! than the spanning `c.<cds_end>_*1ins…`. That rule is stated and locked for
//! single variants in `five_prime_boundary_delins_unification.rs`.
//!
//! The merge path reached the same junction without it. `c.[11_12dup;16_17insC]`
//! normalized to `c.18_*1insCAT`, and re-normalizing *that* applied the clamp
//! after all — so the first answer was not a fixed point, and the two shuffle
//! directions did not even agree on the second pass:
//!
//! ```text
//! once   c.18_*1insCAT
//! twice  c.18delinsTCAT   (3')
//! twice  c.16_17insATC    (5')
//! ```
//!
//! Found by the transcript-axis sweep added for #1283, and only at the full
//! corpus — the two triggering sequences are corpus indices 26 and 36, outside
//! the 4-seed default prefix (#1295).

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::normalize::{NormalizeConfig, ShuffleDirection};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// The transcript from the issue: 20 bases, CDS `1..=18`, so `c.18` is the last
/// coding base and `c.*1` the first 3'UTR base.
const CORE: &str = "TAAATAATAAATATATATTA";
const CDS_START: u64 = 1;
const CDS_END: u64 = 18;

/// The second sequence in the same corpus that reaches the same first pass.
const CORE_SIBLING: &str = "TAATTTATTAATATATATAT";

/// A range starting on a CDS coordinate and ending on a `*` one, e.g. `18_*1`.
///
/// The leading `[.\[;]` is load-bearing. A bare `\d+_\*\d+` also matches the
/// `1_*2` *inside* `*1_*2` — a range wholly within the 3'UTR, which is a
/// legitimate answer and not what this rule forbids — so the left endpoint is
/// pinned to a position that actually starts one: after the `c.`, after the
/// allele's `[`, or after a member separator.
static SPANS_JUNCTION: std::sync::LazyLock<regex::Regex> =
    std::sync::LazyLock::new(|| regex::Regex::new(r"[.\[;]\d+_\*\d+").expect("valid regex"));

fn normalize_once(core: &str, input: &str, dir: ShuffleDirection) -> String {
    let normalizer = Normalizer::with_config(
        SyntheticBuilder::cds(core, CDS_START, CDS_END, Strand::Plus).build(),
        NormalizeConfig::default()
            .with_direction(dir)
            .allow_crossing_boundaries(),
    );
    normalizer
        .normalize(&parse_hgvs(input).expect("parse input"))
        .expect("normalize")
        .to_string()
}

/// Normalizing twice must reproduce the first answer.
///
/// This is what `FERRO_ASSERT_IDEMPOTENT` asserts globally; pinning it here
/// keeps the case named and reproducible without depending on that flag being
/// set, and covers both shuffle directions, which disagreed.
/// The answer each `(core, direction)` settles on, pinned literally.
///
/// The two tests below were property-only — "it is stable" and "it does not
/// span the junction" — and a regression producing a *stable, non-spanning,
/// wrong* string satisfies both. These are what make them value guards. The
/// four differ from one another, which is itself the point: the clamp's answer
/// depends on both the sibling and the direction, so a fix that collapsed them
/// to one form would be caught here rather than passing quietly.
const PINNED: &[(&str, ShuffleDirection, &str)] = &[
    (
        CORE,
        ShuffleDirection::ThreePrime,
        "NM_TEST.1:c.18delinsTCAT",
    ),
    (
        CORE_SIBLING,
        ShuffleDirection::ThreePrime,
        "NM_TEST.1:c.[16_17insC;*1_*2dup]",
    ),
    (CORE, ShuffleDirection::FivePrime, "NM_TEST.1:c.16_17insATC"),
    (
        CORE_SIBLING,
        ShuffleDirection::FivePrime,
        "NM_TEST.1:c.16_17insATC",
    ),
];

/// Every pinned answer must be the one the normalizer actually produces.
#[test]
fn the_clamped_answers_are_unchanged() {
    for (core, dir, expected) in PINNED {
        let once = normalize_once(core, "NM_TEST.1:c.[11_12dup;16_17insC]", *dir);
        assert_eq!(
            &once, expected,
            "the clamped answer moved for {core} under {dir:?}"
        );
    }
}

#[test]
fn a_cis_allele_across_the_cds_utr3_junction_is_a_fixed_point() {
    for dir in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        for core in [CORE, CORE_SIBLING] {
            let once = normalize_once(core, "NM_TEST.1:c.[11_12dup;16_17insC]", dir);
            let twice = normalize_once(
                core,
                &format!("NM_TEST.1:{}", once.split(':').nth(1).unwrap()),
                dir,
            );
            assert_eq!(
                once, twice,
                "normalizing twice changed the answer for {core} under {dir:?}: \
                 {once} -> {twice}"
            );
        }
    }
}

/// The output must not span the CDS/3'UTR junction.
///
/// Stated separately from idempotency because it is the *design* rule, and a
/// fix that merely made the spanning form stable would satisfy the test above
/// while still contradicting the CDS-end clamp every other path applies.
#[test]
fn a_cis_allele_does_not_re_axis_onto_the_utr3() {
    for dir in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        for core in [CORE, CORE_SIBLING] {
            let once = normalize_once(core, "NM_TEST.1:c.[11_12dup;16_17insC]", dir);
            // A range is only *spanning* when it starts on a CDS coordinate and
            // ends on a `*` one. A range wholly inside the 3'UTR (`*1_*2`) is a
            // legitimate answer, so match the junction shape rather than any
            // `_*`, which would flag it too.
            assert!(
                !SPANS_JUNCTION.is_match(&once),
                "output spans the CDS/3'UTR junction for {core} under {dir:?}: {once}"
            );
        }
    }
}
