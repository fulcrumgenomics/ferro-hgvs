//! Regression tests for the 5'/3' homopolymer-boundary insertion/dup
//! unification (idempotency-campaign follow-up to PR #1169; resolves the #387
//! CDS-end saturation latent and the dup-path non-confluence).
//!
//! At a CDS-interior homopolymer that saturates against a coding boundary, every
//! spelling of the same haplotype — dup, direct insertion, `delins`, `unit[N]`,
//! at any start position — must normalize to ONE valid, idempotent form. That
//! convergence is what this file locks. **Which** form it converges to depends on
//! whether a duplication describes the change, and both regimes are covered here:
//!
//! - **The added bases copy an adjacent same-length reference tract** → a `dup`,
//!   which needs no clamp: it sits wholly inside the CDS with two valid anchors.
//!   `repeated.md` L22 prescribes exactly this ("use `NM_024312.4:c.2692_2693dup`
//!   and **not** `c.2686A[10]`") and biocommons agrees, rewriting even a
//!   boundary `delins` into it (`NM_212556.2:c.1delinsAA` → `c.1dup`).
//!   Re-blessed from `c.1delinsAAA` / `c.7delinsAAA` by #1204.
//! - **The added bases are longer than the tract**, so no duplication describes
//!   them → the insertion saturates the boundary and is clamped to
//!   `c.<cds_boundary>delins…`. ferro's design forbids re-axing a CDS-interior
//!   edit onto the 3'UTR (`c.*N`) axis (see the CDS-start/-end clamps in
//!   `src/normalize/mod.rs`), so this remains the canonical boundary form rather
//!   than a `c.<cds_end>_*1ins…` spanning insertion.
//!
//! Spanning *duplications* (`c.7_*1dup`) remain the spec-canonical form and are
//! preserved in both regimes (#401).
//!
//! Before the original fix the dup path returned the 3'-anchored `c.7_*1insAA`
//! verbatim (bypassing both boundary clamps), the direct-insertion path was
//! non-confluent (`c.1delinsAAA` / `c.1_3delinsAAAAA` / `c.1_4delinsAAAAAA`), and
//! every 5' form re-normalized to the INVALID single-position `c.1insAA` (which
//! ferro's own parser rejects).

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

/// Normalize `input` against a plus-strand CDS whose transcript is `core`, with
/// the CDS spanning `cds_start..=cds_end` (1-based), in shuffle `dir`.
fn norm(core: &str, cds_start: u64, cds_end: u64, input: &str, dir: ShuffleDirection) -> String {
    let normalizer = Normalizer::with_config(
        SyntheticBuilder::cds(core, cds_start, cds_end, Strand::Plus).build(),
        NormalizeConfig::default().with_direction(dir),
    );
    normalizer
        .normalize(&parse_hgvs(input).expect("parse input"))
        .expect("normalize")
        .to_string()
}

/// Assert `s` is valid HGVS that ferro can re-parse (guards against emitting the
/// invalid single-position `c.1insAA`).
fn assert_reparses(s: &str) {
    parse_hgvs(s).unwrap_or_else(|e| panic!("normalizer emitted un-parseable HGVS {s:?}: {e}"));
}

/// Assert every spelling in `inputs` normalizes to `expected`, that each output
/// re-parses, and that `expected` is a fixed point.
///
/// Convergence and the fixed point are separate properties and this file needs
/// both: one target every spelling reaches, which re-normalizes to itself.
fn assert_all_converge(
    core: &str,
    cds_start: u64,
    cds_end: u64,
    inputs: &[&str],
    dir: ShuffleDirection,
    expected: &str,
) {
    for input in inputs {
        let out = norm(core, cds_start, cds_end, &format!("NM_TEST.1:{input}"), dir);
        assert_eq!(out, expected, "{dir:?} {input}");
        assert_reparses(&out);
    }
    assert_eq!(
        norm(core, cds_start, cds_end, expected, dir),
        expected,
        "{expected} must be a {dir:?} fixed point",
    );
}

const A7: &str = "AAAAAAA"; // 7-A CDS, cds 1..=7, no UTR either side.

/// Every spelling of "+2 `A` into the 7-`A` tract". The `A[9]` row is what
/// #1204's second half added: before it, the `unit[N]` spelling short-circuited
/// to a raw `ins` and was the one spelling that did NOT converge with the others.
const PLUS_TWO_SPELLINGS: &[&str] = &[
    "c.1_2dup",
    "c.3_4dup",
    "c.6_7dup",
    "c.1_2insAA",
    "c.3_4insAA",
    "c.5_6insAA",
    "c.6_7insAA",
    "c.1delinsAAA",
    "c.1A[9]",
];

/// Eight added `A`s, written from three different start positions. No `A` tract
/// eight bases long exists to copy, so no duplication describes this and it must
/// stay an insertion — which then saturates the coding boundary and clamps.
const PLUS_EIGHT_SPELLINGS: &[&str] = &["c.1_2insAAAAAAAA", "c.3_4insAAAAAAAA", "c.6_7insAAAAAAAA"];

// ---- The duplication regime: all spellings converge to one dup --------------

/// The 5'-most duplicated pair in the tract is `c.1_2`.
#[test]
fn five_prime_boundary_spellings_converge_to_start_dup() {
    assert_all_converge(
        A7,
        1,
        7,
        PLUS_TWO_SPELLINGS,
        ShuffleDirection::FivePrime,
        "NM_TEST.1:c.1_2dup",
    );
}

/// …and the 3'-most is `c.6_7`, which is the shape `repeated.md` L22 shows:
/// `c.2692_2693dup` names the two 3'-most bases of the `A[8]` tract at
/// `c.2686_2693`.
#[test]
fn three_prime_boundary_spellings_converge_to_end_dup() {
    assert_all_converge(
        A7,
        1,
        7,
        PLUS_TWO_SPELLINGS,
        ShuffleDirection::ThreePrime,
        "NM_TEST.1:c.6_7dup",
    );
}

// ---- The clamp regime: an alt longer than the tract has no dup --------------

/// Eight added `A`s against a seven-base tract clamp to the minimal
/// single-base-anchored `delins`: `delete c.1` (`A`), `insert AAAAAAAAA`, net +8.
///
/// This is the case that keeps the #1170 / #387 / #1208 boundary clamps covered
/// now that the dup-able shapes above resolve before reaching them. Without a row
/// of this shape the clamp would be left with no direct test.
#[test]
fn five_prime_over_long_insertion_still_clamps_to_start_delins() {
    assert_all_converge(
        A7,
        1,
        7,
        PLUS_EIGHT_SPELLINGS,
        ShuffleDirection::FivePrime,
        "NM_TEST.1:c.1delinsAAAAAAAAA",
    );
}

/// The 3' mirror: `delete c.7`, `insert AAAAAAAAA`.
#[test]
fn three_prime_over_long_insertion_still_clamps_to_end_delins() {
    assert_all_converge(
        A7,
        1,
        7,
        PLUS_EIGHT_SPELLINGS,
        ShuffleDirection::ThreePrime,
        "NM_TEST.1:c.7delinsAAAAAAAAA",
    );
}

// ---- WITH a real 3'UTR: spanning dups preserved, not over-clamped -----------

/// `AAAAAAAGCG`: A-run CDS c.1..=7, 3'UTR = `GCG` (so c.*1 = 'G' EXISTS).
/// A dup-able homopolymer addition resolves to the CDS-interior dup, and the
/// existence of a representable far side (`c.*1`) does not tempt it onto the UTR
/// axis.
#[test]
fn with_3utr_homopolymer_insertion_becomes_the_end_dup() {
    let out = norm(
        "AAAAAAAGCG",
        1,
        7,
        "NM_TEST.1:c.6_7insAA",
        ShuffleDirection::ThreePrime,
    );
    assert_eq!(out, "NM_TEST.1:c.6_7dup");
}

/// …and the over-long alt on the same transcript still clamps at the CDS end
/// rather than escaping onto `c.*1`.
#[test]
fn with_3utr_over_long_insertion_still_clamps_to_end_delins() {
    let out = norm(
        "AAAAAAAGCG",
        1,
        7,
        "NM_TEST.1:c.6_7insAAAAAAAA",
        ShuffleDirection::ThreePrime,
    );
    assert_eq!(out, "NM_TEST.1:c.7delinsAAAAAAAAA");
}

/// …but a genuine boundary-spanning *duplication* (`AG` = dup of c.7_*1) is the
/// spec-canonical form and must be PRESERVED, never clamped to a delins.
#[test]
fn with_3utr_spanning_duplication_is_preserved() {
    let out = norm(
        "AAAAAAAGCG",
        1,
        7,
        "NM_TEST.1:c.6_7insAG",
        ShuffleDirection::ThreePrime,
    );
    assert_eq!(out, "NM_TEST.1:c.7_*1dup");
}
