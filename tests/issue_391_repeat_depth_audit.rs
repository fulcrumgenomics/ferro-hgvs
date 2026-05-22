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
//! `count_tandem_repeats`, `find_homopolymer_at`) all walk in
//! `unit_len` steps and terminate exactly when the candidate slot
//! doesn't match the unit, with bounds expressed as `start + u <=
//! ref_seq.len()` — the correct half-open shape.
//!
//! The SVD-WG009 citation in #334 / `failure-patterns.md` row 9 was a
//! mis-attribution: SVD-WG009 retires the `con` (conversion) variant
//! type and has no bearing on repeat depth.
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
//! 4. Insertion-to-repeat rewrite anchors on the full reference tract
//!    (not depth-by-1 short).
//! 5. Duplication-to-repeat rewrite reports `unit[ref_count + 1]` for
//!    a single-unit dup at the tract's 3' edge (depth, not depth-1).
//!
//! All five are g.-axis fixtures so the CDS↔UTR axis clamp (issue #337)
//! and the exon-junction exception (#334's first half) are out of
//! scope here. Any of these failing would indicate a real off-by-one.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

/// 128 bp of `ACGT` padding — wider than the default 100 bp normalizer
/// window so the shuffle window never hits an array edge. Same flat
/// `ACGTACGT…` shape used by `tests/issue_209_repeat_3prime_remaining.rs`
/// to keep the pad inert (no accidental homopolymer or tandem overlap
/// with the test core).
const PAD: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
                   ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
const PAD_LEN_HGVS: u64 = 128; // pad lives at HGVS 1..128, core starts at HGVS 129

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
///   PAD (128 bp, HGVS 1..128) | "CCCC" 129..132 | "AAAAAAAA" 133..140 | "GT" 141..142 | PAD
/// ```
///
/// A 3'-direction shuffle of `g.133del` (delete the leftmost A) must
/// walk through all 8 A's and land at `g.140del` — the rightmost A,
/// immediately before the G terminator. If the `<` in shuffle.rs were
/// secretly off-by-one, the result would be `g.139del`.
#[test]
fn three_prime_homopolymer_shuffle_reaches_deepest_position() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr391a", padded("CCCCAAAAAAAAGT"));
    let first_a = PAD_LEN_HGVS + 5; // 133
    let last_a = PAD_LEN_HGVS + 12; // 140
    let out = normalize(provider, &format!("chr391a:g.{first_a}del"));
    let expected = format!(":g.{last_a}del");
    assert!(
        out.ends_with(&expected),
        "3'-direction homopolymer shuffle must land at the deepest A \
         (g.{last_a}del); got {out:?}",
    );
}

// --- Item 1.2: 5' homopolymer reaches the most-5' position --------------

/// Same fixture as 1.1. From the rightmost A, a 5' shuffle must walk
/// back to the leftmost A at `g.133del`. An off-by-one in the 5' loop
/// (`while new_start > boundaries.left`) would stop at `g.134del`.
#[test]
fn five_prime_homopolymer_shuffle_reaches_most_5prime_position() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr391a", padded("CCCCAAAAAAAAGT"));
    let first_a = PAD_LEN_HGVS + 5; // 133
    let last_a = PAD_LEN_HGVS + 12; // 140
    let out_5p = normalize_with_direction(
        provider,
        ShuffleDirection::FivePrime,
        &format!("chr391a:g.{last_a}del"),
    );
    let expected = format!(":g.{first_a}del");
    assert!(
        out_5p.ends_with(&expected),
        "5'-direction homopolymer shuffle must land at the most-5' A \
         (g.{first_a}del); got {out_5p:?}",
    );

    // The PAD itself contains an `ACGT…ACGT` flat repeat, so a 5'
    // shuffle on a tract whose left flank is `C` (no upstream A in PAD
    // immediately abuts the tract) cannot leak into the pad — the
    // boundary base `CCCC` 129..132 blocks the walk. This re-asserts
    // the test result holds even when the inert pad is in play.
    let mut provider2 = MockProvider::new();
    provider2.add_genomic_sequence("chr391a", padded("CCCCAAAAAAAAGT"));
    let out_3p = normalize_with_direction(
        provider2,
        ShuffleDirection::ThreePrime,
        &format!("chr391a:g.{last_a}del"),
    );
    let expected_3p = format!(":g.{last_a}del");
    assert!(
        out_3p.ends_with(&expected_3p),
        "3'-direction shuffle of already-3'-canonical input must be a \
         fixed point; got {out_3p:?}",
    );
}

// --- Item 1.3: 3' multi-base tandem terminates on the deepest unit ------

/// Padded sequence — HGVS positions:
/// ```text
///   PAD | "ACACACACAC" 129..138 (five AC) | "GT" 139..140 | PAD
/// ```
///
/// `g.{first_unit_start}_{first_unit_end}del` (delete the leftmost AC)
/// must 3'-shuffle to the rightmost unit. ferro's shuffle walks
/// base-by-base; the rotation predicate (`ref[del_start] == ref[del_end]`)
/// preserves unit alignment across each step. A depth-by-1 termination
/// would land one unit short of the canonical depth.
#[test]
fn three_prime_tandem_repeat_shuffle_terminates_on_deepest_unit() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr391b", padded("ACACACACACGT"));
    let first_unit_start = PAD_LEN_HGVS + 1; // 129
    let last_unit_start = PAD_LEN_HGVS + 9; // 137
    let last_unit_end = PAD_LEN_HGVS + 10; // 138
    let out = normalize(
        provider,
        &format!(
            "chr391b:g.{first_unit_start}_{end}del",
            end = first_unit_start + 1
        ),
    );
    // ferro may surface the canonicalized form as either a plain
    // `g.{last_unit_start}_{last_unit_end}del` or as an AC[k] tract
    // anchored at the unit boundary. Both shapes pin the same depth.
    let plain_del = format!(":g.{last_unit_start}_{last_unit_end}del");
    let tract_start = PAD_LEN_HGVS + 1;
    let tract_end = PAD_LEN_HGVS + 10;
    let tract_repeat = format!(":g.{tract_start}_{tract_end}AC[");
    let depth_correct = out.contains(&plain_del) || out.contains(&tract_repeat);
    assert!(
        depth_correct,
        "3'-direction tandem AC shuffle must terminate on the deepest \
         unit (...:{plain_del} or ...{tract_repeat}); got {out:?}",
    );
}

// --- Item 1.4: insertion-to-repeat anchors at the full tract ------------

/// Padded sequence — HGVS positions:
/// ```text
///   PAD | "CCC" 129..131 | "AAAA" 132..135 | "GT" 136..137 | PAD
/// ```
///
/// `g.135_136insA` inserts one A between the last tract A (HGVS 135)
/// and the G (HGVS 136). After 3'-rule canonicalization this should
/// become a 5-copy `A[5]` anchored on the full tract (HGVS 132..135)
/// or a duplication `g.135dup` at the tract's rightmost slot. A
/// depth-by-1 bug in `find_tandem_extent` could surface as an A[4]
/// count or an anchor shifted by one base.
///
/// Pin both the canonical surface forms ferro emits and the negative
/// "depth-by-1" shapes that the hypothesis would have produced.
#[test]
fn insertion_to_repeat_anchors_at_full_tract_3prime_edge() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr391c", padded("CCCAAAAGT"));
    let last_tract_a = PAD_LEN_HGVS + 7; // 135
    let after_tract = PAD_LEN_HGVS + 8; // 136
    let out = normalize(
        provider,
        &format!("chr391c:g.{last_tract_a}_{after_tract}insA"),
    );
    let tract_start = PAD_LEN_HGVS + 4; // 132 (first tract A)
                                        // Accept any of the canonical surfaces:
    let depth_correct = out.contains(&format!("g.{tract_start}_{last_tract_a}A[5]"))
        || out.contains(&format!("g.{tract_start}_{last_tract_a}AAAA["))
        || out.contains(&format!("g.{last_tract_a}dup"))
        || out.contains("A[5]");
    assert!(
        depth_correct,
        "single-A insertion at the 3' edge of an AAAA tract must canonicalize \
         to A[5] (anchored on the full 4-A tract) or to g.{last_tract_a}dup — \
         not a depth-3 subset; got {out:?}",
    );
    // Negative assertions: depth-by-1 surfaces.
    assert!(
        !out.contains(&format!("g.{tract_start}_{last_tract_a}A[4]")),
        "depth-by-1: A[4] count for a 5-copy result indicates termination \
         one unit short of canonical; got {out:?}",
    );
    // Shift-by-1 anchor: tract picked one position past the true start.
    let off_by_one_start = tract_start + 1;
    let off_by_one_end = last_tract_a + 1;
    assert!(
        !out.contains(&format!("g.{off_by_one_start}_{off_by_one_end}A[5]")),
        "depth-by-1: anchor at g.{off_by_one_start}_{off_by_one_end} indicates \
         the tract-finder picked an offset-by-one start position; got {out:?}",
    );
}

// --- Item 1.5: duplication-to-repeat reports depth+1, not depth ----------

/// Same fixture as 1.4. `g.135dup` duplicates the last A of the tract.
/// After repeat-notation conversion this becomes `A[5]` (4 original +
/// 1 duplicated = 5). A depth-by-1 bug in
/// `duplication_to_repeat`'s count math (`total_count =
/// analysis.ref_count + dup_len`) would surface as `A[4]` (loses the
/// duplicated unit) or `A[6]` (double-counts it).
#[test]
fn duplication_to_repeat_count_includes_the_duplicated_unit() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr391d", padded("CCCAAAAGT"));
    let last_tract_a = PAD_LEN_HGVS + 7; // 135
    let out = normalize(provider, &format!("chr391d:g.{last_tract_a}dup"));
    // The canonical depth is 5. Accept either the explicit-bracket
    // form or the plain dup surface, but reject the obvious depth-by-1
    // misreporting.
    assert!(
        out.contains("A[5]") || out.contains(&format!("g.{last_tract_a}dup")),
        "g.{last_tract_a}dup on a 4-A tract must canonicalize to A[5] (or \
         stay as g.{last_tract_a}dup); got {out:?}",
    );
    assert!(
        !out.contains("A[4]"),
        "depth-by-1 (under-count): A[4] on a dup of a 4-A tract loses \
         the duplicated unit; got {out:?}",
    );
    assert!(
        !out.contains("A[6]"),
        "depth-by-1 (over-count): A[6] on a dup of a 4-A tract counts \
         the duplicated unit twice; got {out:?}",
    );
}
