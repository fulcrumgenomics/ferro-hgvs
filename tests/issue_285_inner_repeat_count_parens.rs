//! Audit (closes #285): inner repeat-count range `[(150_180)]` must preserve
//! its parentheses on Display.
//!
//! Background: PR #238 (closes #237) fixed the **outer** single-uncertain
//! position range form `c.(a_b)<edit>`. PR #238 explicitly called out the
//! **inner** repeat-count range form as out of scope:
//!
//!   > Inner repeat-count range form `insN[(150_180)]` — the outer `(a_b)`
//!   > is now correctly preserved, but the inner `[(150_180)]` count-range
//!   > drops parens to `[150_180]`. Three fixture rows pin this as a
//!   > remaining divergence.
//!
//! HGVS distinguishes two related shapes:
//!
//! - `[m_n]`   — an **explicit fixed range** (e.g. tandem repeats observed
//!   anywhere from m to n copies).
//! - `[(m_n)]` — an **uncertain count** that lies somewhere in `[m, n]`. The
//!   parentheses are the spec's way of marking the range as uncertain rather
//!   than as a discrete observation interval.
//!
//! Pre-fix, the parser collapsed both shapes into `RepeatCount::Range(m, n)`
//! and Display always emitted `[m_n]` (no parens). This audit pins:
//!
//! 1. `insN[(150_180)]` parses and Display preserves the parens.
//! 2. Round-trip parse → Display → parse → Display is idempotent on the
//!    second pass.
//! 3. The fix touches only the `ins` shape that today accepts the
//!    `N[(m_n)]` form via `parse_repeated_base_insertion`; other shapes
//!    (`dup`, `del`, `inv`) do not parse the inner-parens form on main and
//!    are intentionally left out (they would need orthogonal parser work).
//!    The `delins` shape with a raw `[(m_n)]` body goes through a separate
//!    `InsertedSequence::Range` path that is not the `RepeatCount::Range`
//!    code being fixed here.
//! 4. The pre-existing **explicit fixed-range** forms `[150_180]` and exact
//!    `[150]` are still accepted unchanged and Display unchanged.

use ferro_hgvs::parse_hgvs;

/// Helper: parse, Display, parse again, Display again. Asserts that the
/// second-pass Display equals the first-pass Display (idempotent round trip)
/// and returns the canonical Display string for further assertions.
fn round_trip(input: &str) -> String {
    let v1 = parse_hgvs(input).unwrap_or_else(|e| panic!("first parse of {input:?} failed: {e}"));
    let s1 = v1.to_string();
    let v2 = parse_hgvs(&s1).unwrap_or_else(|e| panic!("second parse of {s1:?} failed: {e}"));
    let s2 = v2.to_string();
    assert_eq!(
        s1, s2,
        "round-trip not idempotent: first display {s1:?} second display {s2:?}",
    );
    s1
}

// ---------------------------------------------------------------------------
// Primary fix: insN[(150_180)] preserves parens on Display.
// ---------------------------------------------------------------------------

#[test]
fn ins_n_inner_uncertain_range_preserves_parens() {
    let input = "NM_000088.3:c.123_124insN[(150_180)]";
    let displayed = round_trip(input);
    assert_eq!(
        displayed, input,
        "inner repeat-count uncertain range must preserve parens"
    );
    assert!(
        displayed.contains("[(150_180)]"),
        "expected '[(150_180)]' substring in {displayed:?}",
    );
}

// ---------------------------------------------------------------------------
// Negative controls: the explicit fixed-range and exact-count forms must
// still work, with no parens added.
// ---------------------------------------------------------------------------

#[test]
fn ins_n_explicit_fixed_range_unchanged() {
    // [m_n] without parens is the explicit fixed-range form and must NOT
    // gain parens after the fix.
    let input = "NM_000088.3:c.123_124insN[150_180]";
    let displayed = round_trip(input);
    assert_eq!(displayed, input);
    assert!(
        !displayed.contains("[(150_180)]"),
        "explicit fixed range must not be re-emitted with parens: {displayed:?}",
    );
}

#[test]
fn ins_n_exact_count_unchanged() {
    let input = "NM_000088.3:c.123_124insN[150]";
    let displayed = round_trip(input);
    assert_eq!(displayed, input);
}

// ---------------------------------------------------------------------------
// Adjacent boundary-uncertain forms inside [] must keep working unchanged.
// These were already supported by `parse_repeat_count` and are pinned here
// as guardrails so the fix doesn't regress them.
// ---------------------------------------------------------------------------

#[test]
fn ins_n_min_uncertain_unchanged() {
    // [(150_?)] form — already supported and rendered as [150_?].
    let input = "NM_000088.3:c.123_124insN[150_?]";
    let displayed = round_trip(input);
    assert_eq!(displayed, input);
}

#[test]
fn ins_n_max_uncertain_unchanged() {
    let input = "NM_000088.3:c.123_124insN[?_180]";
    let displayed = round_trip(input);
    assert_eq!(displayed, input);
}

#[test]
fn ins_n_unknown_count_unchanged() {
    let input = "NM_000088.3:c.123_124insN[?]";
    let displayed = round_trip(input);
    assert_eq!(displayed, input);
}

// ---------------------------------------------------------------------------
// Same fix applies regardless of base identity (the codepath is shared).
// ---------------------------------------------------------------------------

#[test]
fn ins_a_inner_uncertain_range_preserves_parens() {
    let input = "NM_000088.3:c.123_124insA[(5_10)]";
    let displayed = round_trip(input);
    assert_eq!(displayed, input);
}

#[test]
fn ins_g_inner_uncertain_range_preserves_parens() {
    let input = "NM_000088.3:c.123_124insG[(2_4)]";
    let displayed = round_trip(input);
    assert_eq!(displayed, input);
}

// ---------------------------------------------------------------------------
// Boundary cases: pin observed behavior for degenerate / inverted bounds so
// that any future tightening (rejecting `m>n`, rejecting `(0_0)`, collapsing
// `(m_m)` to `[m]`, …) is an intentional change with a touchpoint here.
//
// Today the parser accepts these shapes; no semantic validation runs on
// `RepeatCount::UncertainRange(min, max)` beyond u64 parsing. These tests do
// NOT assert that the shapes are *correct* — they assert only that the
// parser+Display round-trips them verbatim, so we'll notice on the day we
// decide to add validation.
// ---------------------------------------------------------------------------

#[test]
fn ins_n_inner_uncertain_range_equal_bounds_round_trips() {
    // `(m_m)` — degenerate "uncertain" range where both bounds are equal.
    // Today: parsed and rendered verbatim. A future tightening might collapse
    // to `RepeatCount::Exact(m)` and Display as `[m]`.
    let input = "NM_000088.3:c.123_124insN[(5_5)]";
    let displayed = round_trip(input);
    assert_eq!(displayed, input);
}

#[test]
fn ins_n_inner_uncertain_range_zero_bounds_round_trips() {
    // `(0_0)` — zero-copy uncertain range. Semantically meaningless (no
    // insertion at all), but the parser accepts it today. Pin so future
    // validation (e.g. rejecting count == 0) is a deliberate decision.
    let input = "NM_000088.3:c.123_124insN[(0_0)]";
    let displayed = round_trip(input);
    assert_eq!(displayed, input);
}

#[test]
fn ins_n_inner_uncertain_range_inverted_bounds_round_trips() {
    // `(m_n)` with `m > n`. The HGVS spec implies `min <= max`, but the
    // parser does not enforce ordering today. Pin so a future
    // `min <= max` check is a deliberate, test-visible change.
    let input = "NM_000088.3:c.123_124insN[(180_150)]";
    let displayed = round_trip(input);
    assert_eq!(displayed, input);
}

// ---------------------------------------------------------------------------
// Normalize no-op: `insN[(150_180)]` flows through `Normalizer::normalize`
// unchanged. Verifies `normalize/mod.rs`'s early-return for non-`Exact`
// repeat counts (see the `let RepeatCount::Exact(..) = count else { ... }`
// guard in the `NaEdit::Repeat` branch).
// ---------------------------------------------------------------------------

#[test]
fn ins_n_inner_uncertain_range_normalize_is_noop() {
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
    use ferro_hgvs::{MockProvider, Normalizer};

    // 300 bp of `N` so the c.123_124 ins position has reference bytes to
    // hand to the normalizer (the actual bytes don't matter — non-`Exact`
    // counts early-return before touching the reference window).
    let seq: String = "N".repeat(300);
    let len = seq.len() as u64;
    let tx = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        seq,
        Some(1),
        Some(len),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    let mut provider = MockProvider::new();
    provider.add_transcript(tx);

    let input = "NM_TEST.1:c.123_124insN[(150_180)]";
    let variant = parse_hgvs(input).expect("parse");
    let normalizer = Normalizer::new(provider);
    let out = normalizer.normalize(&variant).expect("normalize");
    assert_eq!(
        format!("{}", out),
        input,
        "non-Exact repeat count must early-return unchanged from normalize"
    );
}
