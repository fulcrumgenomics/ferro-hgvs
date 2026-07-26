//! Issue #1202 — a shuffled insertion that saturates at a **transcript**
//! boundary must be clamped on every axis numbered against the transcript.
//!
//! `normalize_cds` has carried a boundary rewrite since #1170 (5', keyed on
//! `cds_start`) and #387 (3', keyed on `cds_end`), but both are gated on
//! `AxisRegion::Cds`. Nothing clamped the *transcript* bounds, so a saturated
//! insertion escaped in two shapes, both of which break the round-trip:
//!
//! | direction | escaped as | why it is broken |
//! |---|---|---|
//! | 5' | `n.1ins<A'>` | single-position insertion — ferro's own parser rejects it (DNA/insertion.md:95-101) |
//! | 3' | `n.<len>_<len+1>ins<A'>` | second anchor is past the transcript end — strict re-normalization raises W4004 `PositionPastEnd` |
//!
//! Three axes reach it, and none of them is "non-coding transcripts":
//!
//! - `n.` — numbers the whole transcript with no CDS sub-axis (#334), so an
//!   `n.` variant on a **coding** `NM_` fails identically to one on an `NR_`.
//! - `r.` — same, via `normalize_rna`.
//! - `c.` in a **UTR** region — outside the `AxisRegion::Cds` gate, so neither
//!   #1170 nor #387 fires; reachable once `cds_start > 1`.
//!
//! Every expectation is the value derived from the fixture by hand (see the
//! per-test comments), not merely `norm(norm(x)) == norm(x)` — an
//! idempotency-only assertion is satisfied by a normalizer that returns its
//! input unchanged, which is exactly the bug class here.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Core used for the 5' cases. Positions 1-4 are `C`, 5-11 `A`, 12-15 `G`.
const CORE_5P: &str = "CCCCAAAAAAAGGGG";

/// Core used for the 3' cases. Positions 13, 14, 15 are `T`, `G`, `C` — chosen
/// so a `GCA` insertion phase-walks to the transcript end while the rotated
/// `AGC` still differs from the preceding `TGC`, i.e. it does NOT promote to a
/// `dup` on the way and really does saturate as a plain insertion.
const CORE_3P: &str = "CCCCAAAAAAAATGC";

fn normalize(provider: MockProvider, input: &str, dir: ShuffleDirection) -> String {
    let normalizer =
        Normalizer::with_config(provider, NormalizeConfig::default().with_direction(dir));
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize {input}: {e}"))
        .to_string()
}

/// The RNA spelling of a DNA core: the `r.` axis alphabet is `a/c/g/u`.
fn as_rna(core: &str) -> String {
    core.to_lowercase().replace('t', "u")
}

/// Assert the exact canonical output, that it re-parses, and that a second
/// pass is a fixed point.
fn assert_canonical(provider: MockProvider, input: &str, dir: ShuffleDirection, expected: &str) {
    let once = normalize(provider.clone(), input, dir);
    assert_eq!(once, expected, "{input} ({dir:?}) canonical form");
    parse_hgvs(&once)
        .unwrap_or_else(|e| panic!("normalized output {once} must re-parse, got: {e}"));
    let twice = normalize(provider, &once, dir);
    assert_eq!(twice, once, "{once} ({dir:?}) must be a fixed point");
}

// ---------------------------------------------------------------------------
// 5' saturation at transcript position 1
// ---------------------------------------------------------------------------

/// `n.1_2insGC` 5'-shuffles to rest immediately before base 1 with the rotated
/// payload `A' = CG`. Identity: insert `CG` before base 1 ≡ delete `ref[1]`
/// (= `C`) and insert `CG ++ C`.
///
///   input   `C|GC|CCCAAAAAAAGGGG` -> `CGCCCCAAAAAAAGGGG`
///   output  `[C->CGC]CCCAAAAAAAGGGG` -> `CGCCCCAAAAAAAGGGG`   (identical)
#[test]
fn noncoding_five_prime_saturated_insertion_clamps_to_transcript_start() {
    assert_canonical(
        SyntheticBuilder::noncoding(CORE_5P, Strand::Plus).build(),
        "NR_TEST.1:n.1_2insGC",
        ShuffleDirection::FivePrime,
        "NR_TEST.1:n.1delinsCGC",
    );
}

/// The `n.` axis numbers the whole transcript with no CDS sub-axis, so an `n.`
/// variant on a **coding** transcript takes the same path and must clamp the
/// same way. This is the case that shows the defect was never about
/// non-coding transcripts.
#[test]
fn coding_transcript_n_axis_five_prime_saturated_insertion_clamps() {
    assert_canonical(
        SyntheticBuilder::cds(CORE_5P, 1, CORE_5P.len() as u64, Strand::Plus).build(),
        "NM_TEST.1:n.1_2insGC",
        ShuffleDirection::FivePrime,
        "NM_TEST.1:n.1delinsCGC",
    );
}

/// Same identity on the `r.` axis; the RNA alphabet renders lowercase.
#[test]
fn rna_five_prime_saturated_insertion_clamps_to_transcript_start() {
    assert_canonical(
        SyntheticBuilder::rna(&as_rna(CORE_5P), Strand::Plus).build(),
        "NR_TEST.1:r.1_2insgc",
        ShuffleDirection::FivePrime,
        "NR_TEST.1:r.1delinscgc",
    );
}

/// Minus-strand non-coding transcript: the `n.` axis is transcript-relative, so
/// the clamp is strand-independent. Paired with
/// `noncoding_minus_strand_three_prime_saturated_insertion_clamps` below so
/// that claim is tested at *both* bounds, not only this one.
#[test]
fn noncoding_minus_strand_five_prime_saturated_insertion_clamps() {
    assert_canonical(
        SyntheticBuilder::noncoding(CORE_5P, Strand::Minus).build(),
        "NR_TEST.1:n.1_2insGC",
        ShuffleDirection::FivePrime,
        "NR_TEST.1:n.1delinsCGC",
    );
}

// ---------------------------------------------------------------------------
// 3' saturation at the transcript end
// ---------------------------------------------------------------------------

/// `n.13_14insGCA` 3'-shuffles to rest immediately after base 15 with the
/// rotated payload `A' = AGC`. Identity: insert `AGC` after base 15 ≡ delete
/// `ref[15]` (= `C`) and insert `C ++ AGC`.
///
///   input   `CCCCAAAAAAAAT|GCA|GC` -> `CCCCAAAAAAAATGCAGC`
///   output  `CCCCAAAAAAAATG[C->CAGC]` -> `CCCCAAAAAAAATGCAGC`  (identical)
///
/// Before the fix this emitted `n.15_16insAGC`, whose `n.16` anchor is past the
/// 15-base transcript end.
#[test]
fn noncoding_three_prime_saturated_insertion_clamps_to_transcript_end() {
    assert_canonical(
        SyntheticBuilder::noncoding(CORE_3P, Strand::Plus).build(),
        "NR_TEST.1:n.13_14insGCA",
        ShuffleDirection::ThreePrime,
        "NR_TEST.1:n.15delinsCAGC",
    );
}

#[test]
fn coding_transcript_n_axis_three_prime_saturated_insertion_clamps() {
    assert_canonical(
        SyntheticBuilder::cds(CORE_3P, 1, CORE_3P.len() as u64, Strand::Plus).build(),
        "NM_TEST.1:n.13_14insGCA",
        ShuffleDirection::ThreePrime,
        "NM_TEST.1:n.15delinsCAGC",
    );
}

#[test]
fn rna_three_prime_saturated_insertion_clamps_to_transcript_end() {
    assert_canonical(
        SyntheticBuilder::rna(&as_rna(CORE_3P), Strand::Plus).build(),
        "NR_TEST.1:r.13_14insgca",
        ShuffleDirection::ThreePrime,
        "NR_TEST.1:r.15delinscagc",
    );
}

/// 3'-bound mirror of the minus-strand 5' case above.
#[test]
fn noncoding_minus_strand_three_prime_saturated_insertion_clamps() {
    assert_canonical(
        SyntheticBuilder::noncoding(CORE_3P, Strand::Minus).build(),
        "NR_TEST.1:n.13_14insGCA",
        ShuffleDirection::ThreePrime,
        "NR_TEST.1:n.15delinsCAGC",
    );
}

// ---------------------------------------------------------------------------
// The `c.` axis: UTR regions saturate at the transcript bounds too
// ---------------------------------------------------------------------------

/// `normalize_cds`'s own clamps are keyed on `cds_start`/`cds_end` and gated on
/// `AxisRegion::Cds`, so a **UTR-region** `c.` insertion is outside them and
/// saturates against the *transcript* bounds instead. Reachable only when a
/// 5'UTR exists (`cds_start > 1`); with `cds_start == 1` the #1170 clamp sits on
/// the transcript start and masks it.
///
/// With `cds_start = 5`, `c.-4` is transcript position 1. `c.-4_-3insGC`
/// 5'-shuffles to rest before base 1 with `A' = CG`, so the identity gives
/// `delete ref[1]` (= `C`), `insert CG ++ C`.
///
///   input   `C|GC|CCCAAAAAAAGGGG` -> `CGCCCCAAAAAAAGGGG`
///   output  `[C->CGC]CCCAAAAAAAGGGG` -> `CGCCCCAAAAAAAGGGG`   (identical)
///
/// Before the fix this emitted the degenerate `c.-4insCG`.
#[test]
fn cds_five_utr_saturated_insertion_clamps_to_transcript_start() {
    assert_canonical(
        SyntheticBuilder::cds(CORE_5P, 5, CORE_5P.len() as u64, Strand::Plus).build(),
        "NM_TEST.1:c.-4_-3insGC",
        ShuffleDirection::FivePrime,
        "NM_TEST.1:c.-4delinsCGC",
    );
}

/// 3'UTR mirror. With `cds_end = 11` on a 15-base transcript, `c.*4` is
/// transcript position 15 and `c.*5` does not exist. `c.*2_*3insGCA`
/// 3'-shuffles to rest after base 15 with `A' = AGC`, so the identity gives
/// `delete ref[15]` (= `C`), `insert C ++ AGC`.
///
/// Before the fix this emitted `c.*4_*5insAGC`, whose `c.*5` is past the end.
#[test]
fn cds_three_utr_saturated_insertion_clamps_to_transcript_end() {
    assert_canonical(
        SyntheticBuilder::cds(CORE_3P, 1, 11, Strand::Plus).build(),
        "NM_TEST.1:c.*2_*3insGCA",
        ShuffleDirection::ThreePrime,
        "NM_TEST.1:c.*4delinsCAGC",
    );
}

/// The CDS-region clamps (#1170 / #387) must keep their own boundary — adding
/// the transcript-bound clamp to `normalize_cds` must not move them. With
/// `cds_start == 1` the two boundaries coincide and #1170 owns the rewrite.
#[test]
fn cds_region_clamps_keep_their_own_boundary() {
    assert_canonical(
        SyntheticBuilder::cds(CORE_5P, 1, CORE_5P.len() as u64, Strand::Plus).build(),
        "NM_TEST.1:c.1_2insGC",
        ShuffleDirection::FivePrime,
        "NM_TEST.1:c.1delinsCGC",
    );
}

/// The clamped output must survive a **strict** re-normalization. This is the
/// assertion that fails loudest for the 3' shape: the pre-fix
/// `n.15_16insAGC` is rejected with W4004 `PositionPastEnd`.
#[test]
fn clamped_output_survives_strict_renormalization() {
    // Covers BOTH call sites: the `n.` cases go through `normalize_tx`, the `c.`
    // UTR cases through `normalize_cds`'s own invocation of the clamp. The `c.`
    // rows matter most for the 3' shape — a `c.*<n>_*<n+1>ins` that escapes the
    // clamp is exactly what strict mode rejects with W4004 `PositionPastEnd`.
    for (input, dir, provider) in [
        (
            "NR_TEST.1:n.13_14insGCA",
            ShuffleDirection::ThreePrime,
            SyntheticBuilder::noncoding(CORE_3P, Strand::Plus).build(),
        ),
        (
            "NR_TEST.1:n.1_2insGC",
            ShuffleDirection::FivePrime,
            SyntheticBuilder::noncoding(CORE_5P, Strand::Plus).build(),
        ),
        (
            "NM_TEST.1:c.-4_-3insGC",
            ShuffleDirection::FivePrime,
            SyntheticBuilder::cds(CORE_5P, 5, CORE_5P.len() as u64, Strand::Plus).build(),
        ),
        (
            "NM_TEST.1:c.*2_*3insGCA",
            ShuffleDirection::ThreePrime,
            SyntheticBuilder::cds(CORE_3P, 1, 11, Strand::Plus).build(),
        ),
    ] {
        let out = normalize(provider.clone(), input, dir);
        let strict = Normalizer::with_config(
            provider.clone(),
            NormalizeConfig::strict().with_direction(dir),
        );
        let variant = parse_hgvs(&out).unwrap_or_else(|e| panic!("re-parse {out}: {e}"));
        strict
            .normalize(&variant)
            .unwrap_or_else(|e| panic!("strict re-normalize of {out} (from {input}): {e}"));
    }
}

// ---------------------------------------------------------------------------
// Non-regression: unsaturated insertions are untouched
// ---------------------------------------------------------------------------

/// An insertion resting on the **last interior anchor pair** — `n.14_15` on a
/// 15-base transcript, one base short of saturation — must keep its `ins` form.
///
/// This is the off-by-one guard for the 3' condition: the clamp tests
/// `end == tx_len + 1`, and if it were written `end == tx_len` this case would
/// be wrongly rewritten to `n.14delins…`. None of the other interior cases sit
/// close enough to the 3' bound to catch that.
#[test]
fn insertion_one_base_short_of_the_three_prime_bound_is_not_clamped() {
    // Core positions 13, 14, 15 are `T`, `G`, `C`; inserting `T` between 14 and
    // 15 matches neither neighbour, so it neither shifts nor promotes to a dup.
    assert_canonical(
        SyntheticBuilder::noncoding(CORE_3P, Strand::Plus).build(),
        "NR_TEST.1:n.14_15insT",
        ShuffleDirection::ThreePrime,
        "NR_TEST.1:n.14_15insT",
    );
}

/// An insertion that comes to rest strictly inside the transcript keeps its
/// two-anchor `ins` (or `dup`/repeat) form — the clamp must not fire early.
#[test]
fn interior_insertions_are_not_clamped() {
    for (input, dir, expected) in [
        (
            "NR_TEST.1:n.1_2insCG",
            ShuffleDirection::FivePrime,
            "NR_TEST.1:n.1_2insCG",
        ),
        (
            "NR_TEST.1:n.4_5insA",
            ShuffleDirection::FivePrime,
            "NR_TEST.1:n.5dup",
        ),
        (
            "NR_TEST.1:n.1_2insA",
            ShuffleDirection::ThreePrime,
            "NR_TEST.1:n.1_2insA",
        ),
    ] {
        assert_canonical(
            SyntheticBuilder::noncoding(CORE_5P, Strand::Plus).build(),
            input,
            dir,
            expected,
        );
    }
}
