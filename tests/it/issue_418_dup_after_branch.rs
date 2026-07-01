//! Regression test for issue #418 subgroup (b): NM_001166478.1:c.36_37insTC
//! mis-emitted as `c.33_34dup` (semantic mismatch) where biocommons
//! emits the haplotype-correct `c.35_36dup` under 5'+cross.
//!
//! GenBank `NM_001166478.1` has `c.30..c.40 = "ATTTTTCTTTC"`. Inserting
//! `TC` between `c.36 = C` and `c.37 = T` is haplotype-equivalent to:
//!
//!   * `c.35_36dup` (canonical biocommons form — duplicate the existing
//!     c.35_36 = "TC" tract).
//!   * `c.34_35insTC` (the 5'-most insertion form ferro's shuffle lands
//!     on after walking left through c.35..c.36).
//!
//! and is NOT equivalent to `c.33_34dup` (which would duplicate c.33_34 =
//! "TT", producing a different haplotype).
//!
//! ## Root cause
//!
//! After the 5'-shuffle stops at `result.start = 34` (0-based), the
//! ins→dup recognizer in `normalize_na_edit` calls
//! `insertion_is_duplication`, which returns `true` because *either* the
//! preceding ref tract (`ref[pos-L..pos]`) OR the following one
//! (`ref[pos..pos+L]`) equals the (rotated) alt. The post-shuffle code
//! then unconditionally computed the dup position assuming the BEFORE
//! tract matched (`c.{X-L+1}..c.{X}`), which is correct for the 3'-shift
//! stopping shape but WRONG for the 5'-shift one — at the 5'-most
//! stopping point the alt aligns with the tract just AHEAD (the one
//! shuffle walked through), so the dup region must be `c.{X+1}..c.{X+L}`.
//! Issue #418's `c.36_37insTC` example shifts by exactly 2 (`alt.len()`),
//! and the alt's cyclic period happens to bring the rotation back to its
//! original ("TC") at the stopping position, so the BEFORE-vs-AFTER
//! choice is the *only* signal that distinguishes the two dup tracts.
//!
//! ## Fix
//!
//! `normalize_na_edit`'s ins→dup branch now distinguishes which side
//! matched and computes dup positions accordingly. Tests pin both
//! directions to catch a future regression that conflates them.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Synthetic transcript with `c.30..c.45 = "ATTTTTCTTTCTGGTC"`, matching
/// the NM_001166478.1 geometry around c.36_37 (the issue's spec).
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let mut seq = String::new();
    for _ in 0..29 {
        seq.push('A');
    } // c.1..c.29 filler
    seq.push_str("ATTTTTCTTTCTGGTC"); // c.30..c.45 (16 chars)
    while seq.len() < 1500 {
        seq.push('G');
    }
    let len = seq.len() as u64;
    let transcript = Transcript::new(
        "NM_TEST418B.1".to_string(),
        Some("TEST418B".to_string()),
        Strand::Plus,
        seq,
        Some(1),    // cds_start: c.1 at tx byte 1 (1-based)
        Some(1300), // cds_end (synthetic)
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

fn normalize_with(direction: ShuffleDirection, input: &str) -> String {
    let normalizer = Normalizer::with_config(
        provider(),
        NormalizeConfig::default()
            .with_direction(direction)
            .allow_crossing_boundaries(),
    );
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{}", normalized)
}

#[test]
fn ins_tc_five_prime_emits_after_branch_dup() {
    // The off-by-2 bug: pre-fix this emitted `c.33_34dup` (duplicating
    // the WRONG tract — c.33_34 = "TT", not the alt's "TC"). Fixed to
    // emit `c.35_36dup` (the canonical after-tract dup matching
    // biocommons).
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST418B.1:c.36_37insTC"),
        "NM_TEST418B.1:c.35_36dup",
    );
}

#[test]
fn ins_tc_three_prime_unaffected() {
    // 3'-direction must not change behavior. At the 3'-shift stopping
    // point the BEFORE tract matches, so the existing
    // `(X - L + 1, X)` formula stays canonical.
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, "NM_TEST418B.1:c.36_37insTC"),
        "NM_TEST418B.1:c.36_37dup",
    );
}

#[test]
fn ins_tc_at_c35_five_prime_already_canonical() {
    // #882: `c.35_36insTC` sits at the OUT-OF-PHASE cut relative to the
    // adjacent tract — the candidate dup `c.35_36dup` would decode to a
    // different sequence than inserting "TC" here, so the phase gate rejects
    // the conversion and the variant stays a plain `ins` (duplication.md:19).
    // (The in-phase sibling `c.36_37insTC` still becomes a dup — see
    // `ins_tc_five_prime_emits_after_branch_dup`.)
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST418B.1:c.35_36insTC"),
        "NM_TEST418B.1:c.35_36insTC",
    );
}

#[test]
fn ins_tc_at_c35_three_prime_already_canonical() {
    // #882: same out-of-phase input under 3'-direction. The candidate dup
    // would decode to a different sequence than the input insertion, so the
    // phase gate rejects it and the variant stays a plain `ins`
    // (duplication.md:19). Contrast the in-phase `c.36_37insTC`, which still
    // becomes `c.36_37dup` (see `ins_tc_three_prime_unaffected`).
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, "NM_TEST418B.1:c.35_36insTC"),
        "NM_TEST418B.1:c.35_36insTC",
    );
}

#[test]
fn ins_outside_tract_unaffected() {
    // Negative: inserting OUTSIDE the TC tract (here at c.42_43 which
    // sits in the GG region) must not trigger the after-branch fix.
    // The alt "TC" doesn't match adjacent ref on either side; the
    // recognizer leaves it as a plain ins. Pin the result so a stray
    // after-branch fire on non-dup inputs is caught.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST418B.1:c.42_43insTC"),
        "NM_TEST418B.1:c.42_43insTC",
    );
}

#[test]
fn ins_single_base_after_branch_five_prime_emits_following_position() {
    // Single-base AFTER-match regression (CodeRabbit on PR #420). On
    // the synthetic transcript above, `c.40..c.45 = "CTGGTC"`, so
    // `c.40 = C` and `c.41 = T`. Inserting `T` at the `c.40_41`
    // boundary creates an AFTER-only single-copy match: the alt `"T"`
    // equals `ref[c.41]` (AFTER) but not `ref[c.40] = "C"` (BEFORE),
    // and the T-tract has `ref_count == 1` so the path (a) tract-
    // aligned early return declines (it only fires for `ref_count
    // >= 2`). Under 5'-direction the shuffle cannot walk further left
    // (`c.39 = T` gives a different haplotype than `c.40 = C` after
    // insertion), so `result.start` stays at the original position
    // and path (b)'s `insertion_is_duplication` shortcut fires with
    // `before_matches = false` and `dup_len == 1`.
    //
    // Pre-fix the `dup_len == 1` shortcut emitted `(new_start,
    // new_start)` = `c.40dup` regardless of which side matched —
    // duplicating the wrong base (the preceding C, not the alt's T).
    // Post-fix the shortcut respects `before_matches` symmetrically
    // with the multi-base case and emits `c.41dup`, the position of
    // the actual duplicated base.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST418B.1:c.40_41insT"),
        "NM_TEST418B.1:c.41dup",
    );
}

#[test]
fn ins_single_base_before_branch_three_prime_unchanged() {
    // Single-base BEFORE-match negative control: same input
    // `c.40_41insT` under 3'-direction. The 3'-shuffle walks right
    // by one to `c.41_42` (the haplotype is preserved across that
    // step), at which point the BEFORE tract `ref[c.41]` = "T"
    // matches the alt. `before_matches = true` selects the original
    // `(new_start, new_start)` formula, which must continue to emit
    // `c.41dup`. Pins that the single-base BEFORE branch keeps
    // anchoring at `new_start`.
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, "NM_TEST418B.1:c.40_41insT"),
        "NM_TEST418B.1:c.41dup",
    );
}
