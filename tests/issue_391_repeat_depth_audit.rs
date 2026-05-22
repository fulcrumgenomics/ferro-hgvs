//! Audit regression file for issue #391 item 1 (the deferred half of #334):
//! "repeat-region depth-by-1" off-by-one in the 3'-rule shuffle.
//!
//! # Audit hypothesis (from #334)
//!
//! > investigate whether ferro's repeat termination uses `<` vs `<=`;
//! > verify against `SVD-WG009` semantics
//!
//! # Audit verdict
//!
//! **No off-by-one exists.** The strict-less-than predicates in
//! `src/normalize/shuffle.rs` (`while new_end < boundaries.right` and the
//! symmetric `while new_start > boundaries.left`) are the correct
//! half-open exclusive bounds — the documented contract for
//! `Boundaries` (see `shuffle.rs` module docs) is that `right` is
//! 0-based exclusive, so `<` terminates one step before the bound is
//! breached. Switching to `<=` would index past the array or violate
//! the boundary contract.
//!
//! The tract-extension predicates in `rules.rs` (`extend_tandem_tract`,
//! `count_tandem_repeats`, `find_homopolymer_at`) walk in `unit_len`
//! steps and terminate exactly when the candidate slot doesn't match
//! the unit. The right-walk bound is `end + u <= ref_seq.len()`
//! (`rules.rs:483`), the left-walk bound is `start >= u`
//! (`rules.rs:478`) — both half-open shapes correct.
//!
//! The SVD-WG009 citation in #334 / `failure-patterns.md` row 9 was a
//! mis-attribution: SVD-WG009 retires the `con` (conversion) variant
//! type (see `assets/hgvs-nomenclature/docs/consultation/SVD-WG009.md`)
//! and has no bearing on repeat depth.
//!
//! # What this file pins
//!
//! Tests that the depth-by-1 hypothesis would have flipped:
//!
//! 1. 3'-direction homopolymer shuffle terminates at the deepest
//!    (most-3') position within the tract.
//! 2. 5'-direction homopolymer shuffle terminates at the most-5'
//!    position — symmetric mirror.
//! 3. 3'-direction multi-base tandem shuffle terminates unit-aligned
//!    on the deepest unit slot.
//! 4. 5'-direction multi-base tandem shuffle terminates unit-aligned
//!    on the most-5' unit slot (5'/3' symmetry for tandems).
//! 5. Insertion-to-repeat rewrite anchors on the full reference tract
//!    and counts `ref_count + added_copies` (not `ref_count + k - 1`).
//!    Uses a 2-base `insAA` so the input actually routes through
//!    `insertion_to_repeat` (single-base `insA` short-circuits to
//!    `insertion_to_duplication` and never invokes the repeat-count
//!    arithmetic this test is meant to probe).
//! 6. Duplication-to-repeat rewrite reports `total_count = ref_count +
//!    dup_len` for a multi-base dup at the tract's 3' edge. Uses a
//!    2-base `g.134_135dup` so the input clears `duplication_to_repeat`'s
//!    `dup_len >= 2` gate (single-base dup falls through and exits as
//!    a plain `dup`, not exercising the count math).
//!
//! All six are g.-axis fixtures so the CDS↔UTR axis clamp (issue #337)
//! and the exon-junction exception (#334's first half) are out of
//! scope here. Any of these failing would indicate a real off-by-one.

use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// 256 bp of `ACGT` padding — matches the convention in
/// `tests/issue_209_repeat_3prime_remaining.rs::PAD` (also 256 bp).
/// Wider than the default 100 bp normalizer window so the shuffle
/// window never hits an array edge, with comfortable headroom in case
/// the default `window_size` grows in a future config bump. Flat
/// `ACGTACGT…` so the pad never accidentally overlaps a homopolymer or
/// tandem unit in the test core.
const PAD: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
                   ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
                   ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
                   ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
const PAD_LEN_HGVS: u64 = 256; // pad lives at HGVS 1..256, core starts at HGVS 257

fn padded(core: &str) -> String {
    format!("{PAD}{core}{PAD}")
}

fn normalize(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant =
        parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for {:?}: {:?}", input, e));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize failed for {:?}: {:?}", input, e));
    format!("{}", normalized)
}

fn normalize_with_direction(
    provider: MockProvider,
    direction: ShuffleDirection,
    input: &str,
) -> String {
    let config = NormalizeConfig::default().with_direction(direction);
    let normalizer = Normalizer::with_config(provider, config);
    let variant =
        parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for {:?}: {:?}", input, e));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize failed for {:?}: {:?}", input, e));
    format!("{}", normalized)
}

// --- Item 1.1: 3' homopolymer reaches the deepest position --------------

/// Padded sequence — HGVS positions relative to the start:
/// ```text
///   PAD (256 bp, HGVS 1..256) | "CCCC" 257..260 | "AAAAAAAA" 261..268 | "GT" 269..270 | PAD
/// ```
///
/// A 3'-direction shuffle of `g.261del` (delete the leftmost A) must
/// walk through all 8 A's and land at `g.268del` — the rightmost A,
/// immediately before the G terminator. If the `<` in shuffle.rs were
/// secretly off-by-one, the result would be `g.267del`.
#[test]
fn three_prime_homopolymer_shuffle_reaches_deepest_position() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr391a", padded("CCCCAAAAAAAAGT"));
    let first_a = PAD_LEN_HGVS + 5; // 261
    let last_a = PAD_LEN_HGVS + 12; // 268
    let out = normalize(provider, &format!("chr391a:g.{first_a}del"));
    assert_eq!(
        out,
        format!("chr391a:g.{last_a}del"),
        "3'-direction homopolymer shuffle must land at the deepest A \
         (g.{last_a}del)",
    );
}

// --- Item 1.2: 5' homopolymer reaches the most-5' position --------------

/// Same fixture as 1.1. From the rightmost A, a 5' shuffle must walk
/// back to the leftmost A at `g.{first_a}del`. An off-by-one in the 5'
/// loop (`while new_start > boundaries.left`) would stop one position
/// short.
#[test]
fn five_prime_homopolymer_shuffle_reaches_most_5prime_position() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr391a", padded("CCCCAAAAAAAAGT"));
    let first_a = PAD_LEN_HGVS + 5; // 261
    let last_a = PAD_LEN_HGVS + 12; // 268
    let out_5p = normalize_with_direction(
        provider,
        ShuffleDirection::FivePrime,
        &format!("chr391a:g.{last_a}del"),
    );
    assert_eq!(
        out_5p,
        format!("chr391a:g.{first_a}del"),
        "5'-direction homopolymer shuffle must land at the most-5' A \
         (g.{first_a}del)",
    );

    // Re-asserts the 3'-direction fixed-point invariant on an already-
    // canonical input. Uses a fresh contig name so the two providers
    // can't accidentally share state.
    let mut provider_3p = MockProvider::new();
    provider_3p.add_genomic_sequence("chr391a_3p", padded("CCCCAAAAAAAAGT"));
    let out_3p = normalize_with_direction(
        provider_3p,
        ShuffleDirection::ThreePrime,
        &format!("chr391a_3p:g.{last_a}del"),
    );
    assert_eq!(
        out_3p,
        format!("chr391a_3p:g.{last_a}del"),
        "3'-direction shuffle of already-3'-canonical input must be a \
         fixed point",
    );
}

// --- Item 1.3: 3' multi-base tandem terminates on the deepest unit ------

/// Padded sequence — HGVS positions:
/// ```text
///   PAD | "ACACACACAC" 257..266 (five AC) | "GT" 267..268 | PAD
/// ```
///
/// `g.257_258del` (delete the leftmost AC) must 3'-shuffle to
/// `g.265_266del` (the rightmost AC). ferro's shuffle walks
/// base-by-base; the rotation predicate (`ref[del_start] == ref[del_end]`)
/// preserves unit alignment across each step. A depth-by-1
/// termination would land one unit short of the canonical depth.
#[test]
fn three_prime_tandem_repeat_shuffle_terminates_on_deepest_unit() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr391b", padded("ACACACACACGT"));
    let first_unit_start = PAD_LEN_HGVS + 1; // 257
    let last_unit_start = PAD_LEN_HGVS + 9; // 265
    let last_unit_end = PAD_LEN_HGVS + 10; // 266
    let out = normalize(
        provider,
        &format!(
            "chr391b:g.{first_unit_start}_{end}del",
            end = first_unit_start + 1
        ),
    );
    // ferro may surface the canonicalized form as either a plain
    // `g.{last_unit_start}_{last_unit_end}del` or as an AC[k] tract
    // anchored at the same unit boundary. Both shapes pin the same
    // depth.
    let plain_del = format!("chr391b:g.{last_unit_start}_{last_unit_end}del");
    let tract_repeat = format!("chr391b:g.{first_unit_start}_{last_unit_end}AC[");
    assert!(
        out == plain_del || out.starts_with(&tract_repeat),
        "3'-direction tandem AC shuffle must terminate on the deepest \
         unit ({plain_del} or {tract_repeat}…); got {out:?}",
    );
}

// --- Item 1.4: 5' multi-base tandem terminates on the most-5' unit ------

/// Same fixture as 1.3, but starting from the rightmost AC. The
/// 5'-direction shuffle must walk back to the leftmost AC at
/// `g.257_258del`. Pins the 5'/3' symmetry for the multi-base tandem
/// case — exercises the `while new_start > boundaries.left` loop in
/// `shuffle.rs:120`.
#[test]
fn five_prime_tandem_repeat_shuffle_terminates_on_most_5prime_unit() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr391b", padded("ACACACACACGT"));
    let first_unit_start = PAD_LEN_HGVS + 1; // 257
    let first_unit_end = PAD_LEN_HGVS + 2; // 258
    let last_unit_start = PAD_LEN_HGVS + 9; // 265
    let last_unit_end = PAD_LEN_HGVS + 10; // 266
    let out = normalize_with_direction(
        provider,
        ShuffleDirection::FivePrime,
        &format!("chr391b:g.{last_unit_start}_{last_unit_end}del"),
    );
    let plain_del = format!("chr391b:g.{first_unit_start}_{first_unit_end}del");
    let tract_repeat = format!("chr391b:g.{first_unit_start}_{last_unit_end}AC[");
    assert!(
        out == plain_del || out.starts_with(&tract_repeat),
        "5'-direction tandem AC shuffle must terminate on the most-5' \
         unit ({plain_del} or {tract_repeat}…); got {out:?}",
    );
}

// --- Item 1.5: insertion-to-repeat counts ref_count + added_copies ------

/// Padded sequence — HGVS positions:
/// ```text
///   PAD | "CCC" 257..259 | "AAAA" 260..263 | "GT" 264..265 | PAD
/// ```
///
/// `g.263_264insAA` inserts two A's between the last tract A (HGVS
/// 263) and the G (HGVS 264). A 2-base insertion clears
/// `insertion_to_repeat`'s `added_copies >= 2` gate
/// (`rules.rs:204`) — a single-base `insA` would short-circuit to
/// `insertion_to_duplication` and never invoke the repeat-count
/// arithmetic this test is meant to probe.
///
/// After canonicalization the result is `g.260_263A[6]` — the full
/// tract span (HGVS 260..263, 4 reference A's) with `total_count =
/// ref_count + added_copies = 4 + 2 = 6`. A depth-by-1 bug in the
/// count math would surface as `A[5]` (forgot one added copy) or
/// `A[7]` (double-counted). A shift-by-one in `find_tandem_extent`'s
/// anchor selection would shift the position window by ±1.
#[test]
fn insertion_to_repeat_counts_added_copies_at_full_tract() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr391c", padded("CCCAAAAGT"));
    let tract_start = PAD_LEN_HGVS + 4; // 260
    let last_tract_a = PAD_LEN_HGVS + 7; // 263
    let after_tract = PAD_LEN_HGVS + 8; // 264
    let out = normalize(
        provider,
        &format!("chr391c:g.{last_tract_a}_{after_tract}insAA"),
    );

    let canonical = format!("chr391c:g.{tract_start}_{last_tract_a}A[6]");
    assert_eq!(
        out, canonical,
        "2-base insertion at the 3' edge of an AAAA tract must canonicalize \
         to A[6] anchored on the full 4-A tract",
    );

    // Negative assertions: depth-by-1 / shift-by-one shapes.
    assert!(
        !out.contains("A[5]"),
        "depth-by-1 (under-count): A[5] indicates total_count = ref_count + \
         k - 1 instead of ref_count + k; got {out:?}",
    );
    assert!(
        !out.contains("A[7]"),
        "depth-by-1 (over-count): A[7] indicates total_count = ref_count + \
         k + 1 instead of ref_count + k; got {out:?}",
    );
    let off_by_one_start = tract_start + 1;
    let off_by_one_end = last_tract_a + 1;
    assert!(
        !out.contains(&format!("g.{off_by_one_start}_{off_by_one_end}A[")),
        "shift-by-one: anchor at g.{off_by_one_start}_{off_by_one_end} \
         indicates the tract-finder picked an offset-by-one start position; \
         got {out:?}",
    );
}

// --- Item 1.6: duplication-to-repeat sums ref_count + dup_len -----------

/// Same fixture as 1.5. `g.262_263dup` duplicates the last two A's of
/// the tract (HGVS 262 and 263). A 2-base dup clears
/// `duplication_to_repeat`'s `dup_len >= 2` gate
/// (`rules.rs:1095`) — a single-base `dup` falls through and exits as
/// a plain `dup` without invoking the homopolymer-to-repeat conversion
/// that owns the count arithmetic.
///
/// Canonical result: `g.260_263A[6]` (4 reference + 2 duplicated = 6
/// total). A depth-by-1 bug in the count math (`total_count =
/// analysis.ref_count + dup_len` at `rules.rs:1111`) would surface as
/// `A[5]` (loses one duplicated base) or `A[7]` (double-counts).
#[test]
fn duplication_to_repeat_sums_ref_count_and_dup_len() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr391d", padded("CCCAAAAGT"));
    let tract_start = PAD_LEN_HGVS + 4; // 260
    let last_tract_a = PAD_LEN_HGVS + 7; // 263
    let dup_start = PAD_LEN_HGVS + 6; // 262
    let out = normalize(
        provider,
        &format!("chr391d:g.{dup_start}_{last_tract_a}dup"),
    );

    let canonical = format!("chr391d:g.{tract_start}_{last_tract_a}A[6]");
    assert_eq!(
        out, canonical,
        "2-base dup of the tract's 3'-edge A's must canonicalize to A[6] \
         (4 reference + 2 duplicated)",
    );

    assert!(
        !out.contains("A[5]"),
        "depth-by-1 (under-count): A[5] indicates total_count = ref_count + \
         dup_len - 1 instead of ref_count + dup_len; got {out:?}",
    );
    assert!(
        !out.contains("A[7]"),
        "depth-by-1 (over-count): A[7] indicates total_count = ref_count + \
         dup_len + 1 instead of ref_count + dup_len; got {out:?}",
    );
}
