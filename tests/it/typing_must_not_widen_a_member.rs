//! Typing may widen a member's rendered extent — but not onto a sibling.
//!
//! The partition model's invariant is that **a member's boundaries come from the
//! input's asserted partition and its bases come from the derived sequence**.
//! `partition_block_preserving` honours it structurally: it takes each member's
//! span from the description and reads that member's payload out of the applied
//! window.
//!
//! There is a third stage after partitioning, and it is not covered by that
//! guarantee: **typing**, which decides whether a `(span, payload)` pair renders
//! as `dup`, `ins`, `delins`, `inv` or a repeat. `repeated.md` B2 spells an edit
//! inside a tandem tract as a copy count over the **whole tract**, so typing can
//! *widen* a member's rendered extent — the asserted member is `g.2_3dup`, the
//! rendered member is `g.1_9T[11]`. Same bases, same asserted member, wider
//! rendered span. Under the ruling that the partition is the unit of
//! normalization, that widening is a **boundary change**, and it needs the same
//! licence any other boundary move needs.
//!
//! Widening as such is legal and must stay legal — a lone member has no sibling
//! to reach, and `a_widening_with_no_sibling_on_the_tract_is_still_licensed`
//! pins that. What is not licensed is widening **into territory a sibling
//! claims**:
//!
//! ```text
//! TTTTTTTTTAATATATTTTA
//!   g.[2_3dup;5_6insA]   ->  g.[1_9T[11];5_6insA]   junction 5|6 inside 1_9
//!   g.[2dup;4dup;6T>A]   ->  g.[1_9T[11];6T>A]      base 6 inside 1_9
//! ```
//!
//! The pair claims one interbase point twice, so it denotes no single sequence
//! at all: `parse_hgvs` accepts it, the SPDI apply oracle declines it. That is a
//! defect and not a reading — an overlapping allele has no clause to weigh —
//! so nothing here is re-blessable.
//!
//! # Where the widening was caught, and why the repair is where it is
//!
//! Measured by instrumenting every per-member correction pass in
//! `Normalizer::normalize_allele` and printing each member's rendered extent
//! after each stage, for `TEMPLATE:g.[2_3dup;5_6insA]` under a 5' shuffle. Two
//! distinct widenings, one repaired and one not:
//!
//! | stage | member before | member after | repaired? |
//! |---|---|---|---|
//! | `normalize_core` (`rules::duplication_to_repeat`) | `g.2_3dup` | `g.1_9T[11]` | yes |
//! | `canonical_split_for_variant` → `normalize_core` (`rules::insertion_to_repeat`) | `g.1delinsTTT` | `g.1_9T[11]` | **no** |
//!
//! The second reaches the allele a second time, from
//! `Normalizer::canonicalize_from_sequence`: the re-derivation hands the member
//! back as `g.1delinsTTT`, and `merge::demote_repeats_spanning_siblings` is
//! *differential* — it reads the length from the pre-normalization member's edit
//! type (`Deletion` / `Duplication` / `Insertion`) — so a `delins` yields
//! nothing to re-spell to and the pass declines. It declines a second time on
//! the fixed point, where `before == after == g.1_9T[11]` and there is no growth
//! to difference at all. Both declines were observed directly, not inferred.
//!
//! Typing itself is sibling-blind by construction — `normalize_core` sees one
//! member — so the licence check has to sit where siblings are visible, which is
//! that same pass. It now falls back to reading the growth off the repeat and
//! re-spells it at **the asserted member's own junction**.
//!
//! # The near-miss, and why the negative half names it
//!
//! The first cut re-spelled the growth at the tract's 3' junction, which is what
//! the pass's differential routes do — safe there only because the clamps that
//! run immediately after pull the member back off the sibling. Both clamps refuse
//! a member whose span *grew* relative to its snapshot, and a one-base `delins`
//! re-spelled as a two-base `dup` is exactly that geometry, so nothing bounded
//! it. The result was `g.[5_6insA;8_9dup]`: disjoint, re-parsing, warning-free,
//! in bounds, and denoting a different sequence.
//!
//! That cut turned every overlapping row in the three exhaustive sweeps in
//! `cis_junction_crossing_shift.rs` into a *sequence-changing* row — strictly
//! worse, because a consumer can reject an overlap and cannot see this. The
//! sweeps' `overlapping.is_empty()` assertion went green on it. So the negative
//! half below names that string too: a guard that only asserts the right answer
//! would have reported the wrong fix as progress.

use crate::common::cis_apply_oracle::{apply, member_interbase_spans, normalize_in};
use ferro_hgvs::ShuffleDirection;

/// Nine `T` at positions 1-9 — a tract long enough to canonicalise to a repeat,
/// with non-`T` bases after it so the tract has a definite 3' end. The same
/// constant the three exhaustive sweeps use.
const TRACT: &str = "TTTTTTTTTAATATATTTTA";

#[test]
fn a_typed_repeat_does_not_widen_over_a_siblings_junction() {
    // The asserted `2_3dup` types as `1_9T[11]`, whose tract contains the
    // sibling's `5|6` junction. Re-spelled at the boundary the partition
    // asserted and then clamped, the pair settles two members apart.
    assert_eq!(
        normalize_in(
            TRACT,
            "TEMPLATE:g.[2_3dup;5_6insA]",
            ShuffleDirection::FivePrime
        ),
        "TEMPLATE:g.[3_4dup;5_6insA]"
    );
}

#[test]
fn a_typed_repeat_does_not_widen_over_a_siblings_bases() {
    // Same widening, reached from two duplications that merge; the sibling
    // claims base 6 rather than a junction.
    assert_eq!(
        normalize_in(
            TRACT,
            "TEMPLATE:g.[2dup;4dup;6T>A]",
            ShuffleDirection::FivePrime
        ),
        "TEMPLATE:g.6_7insTA"
    );
}

/// Both answers denote what the input denoted, checked through the SPDI applier
/// rather than through the normalizer.
///
/// The two rows above pin strings, which is a *canonical form* claim. This is
/// the claim that actually cannot be re-blessed: an overlapping allele denotes no
/// sequence, and a disjoint allele denoting the wrong one is the silent failure
/// the near-miss produced.
#[test]
fn the_repaired_answers_denote_the_input_sequence() {
    for input in ["TEMPLATE:g.[2_3dup;5_6insA]", "TEMPLATE:g.[2dup;4dup;6T>A]"] {
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            let output = normalize_in(TRACT, input, direction);
            let want = apply(TRACT, input).expect("the input denotes a sequence");
            let got = apply(TRACT, &output).unwrap_or_else(|| {
                panic!("{input} [{direction:?}] -> {output} denotes no single sequence")
            });
            assert_eq!(got, want, "{input} [{direction:?}] -> {output}");
        }
    }
}

/// The forbidden outputs, named exactly.
///
/// Two shapes, and both must be named. The first two strings are the swallowing
/// form this defect is: a tract-wide copy count containing a sibling. The third
/// is the near-miss described in the module docs — disjoint, so every oracle in
/// the tree is blind to it except an apply comparison, and produced by re-spelling
/// the growth at the tract's 3' end instead of at the asserted boundary.
#[test]
fn the_swallowing_and_near_miss_forms_are_never_emitted() {
    for (input, forbidden) in [
        (
            "TEMPLATE:g.[2_3dup;5_6insA]",
            "TEMPLATE:g.[1_9T[11];5_6insA]",
        ),
        ("TEMPLATE:g.[2_3dup;5_6insA]", "TEMPLATE:g.[5_6insA;8_9dup]"),
        ("TEMPLATE:g.[2dup;4dup;6T>A]", "TEMPLATE:g.[1_9T[11];6T>A]"),
    ] {
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            let output = normalize_in(TRACT, input, direction);
            assert_ne!(
                output, forbidden,
                "{input} [{direction:?}] emitted the forbidden form"
            );
        }
    }
}

/// The negative control the repair must not break: widening is what
/// `repeated.md` B2 asks for, and with no sibling on the tract there is nothing
/// to swallow.
///
/// Both entries are lone members, which is the only configuration in which the
/// full tract-wide extent is unambiguously right. If either of these starts
/// reporting a narrower spelling, the licence check has become a ban.
#[test]
fn a_widening_with_no_sibling_on_the_tract_is_still_licensed() {
    assert_eq!(
        normalize_in(TRACT, "TEMPLATE:g.1_2dup", ShuffleDirection::FivePrime),
        "TEMPLATE:g.1_9T[11]"
    );
    assert_eq!(
        normalize_in(TRACT, "TEMPLATE:g.2_3insTT", ShuffleDirection::FivePrime),
        "TEMPLATE:g.1_9T[11]"
    );
}

/// The scope boundary, and the reason the repair carries two gates rather than
/// one.
///
/// An **authored** overlap is a different question with a settled and opposite
/// answer: ferro reports the conflict and hands the description back rather than
/// picking an order the spec does not rank (#395, #486, #1004). So the swallowing
/// form entered *as input* must survive verbatim — typing did not widen it, it
/// arrived that way — and the repair must not launder it into a tidy single
/// member that strict mode would then accept.
///
/// This is what the repair's `widened` and "the asserted member did not already
/// span a sibling" gates buy. Without them the same fallback also rewrites
/// `g.[263_265A[7];264_265insC]` and `g.[1005_1009inv;1005_1006A[4]]`, which are
/// pinned as preserved by `repeat_lowering_sibling_junction` and
/// `idempotency_tests` respectively.
#[test]
fn an_authored_swallowing_allele_is_still_preserved_verbatim() {
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_eq!(
            normalize_in(TRACT, "TEMPLATE:g.[1_9T[11];5_6insA]", direction),
            "TEMPLATE:g.[1_9T[11];5_6insA]",
            "an authored conflict must be reported, not repaired ({direction:?})"
        );
    }
}

/// The invariant itself, over the whole family the defect lives in, stated
/// without naming a single output string.
///
/// For every pairing of a tract-growing first member with a sibling — a junction
/// member, a base-claiming member, in both shuffle directions — the rendered
/// members must occupy **pairwise disjoint interbase territory**, and the allele
/// must denote what it denoted. Reading extents from `member_interbase_spans`
/// rather than from the printed endpoints is what makes "rendered extent" mean
/// the same thing here as in the diagnosis: `a_b ins` names two positions while
/// occupying one junction, and a repeat names its tract.
///
/// Deliberately a property and not a table. The three exhaustive sweeps in
/// `cis_junction_crossing_shift.rs` already enumerate far more shapes; what this
/// adds is the *extent* reading, so a future widening that denotes the right
/// sequence while claiming a sibling's territory fails here even where an apply
/// comparison would pass.
#[test]
fn no_rendered_member_claims_a_siblings_territory() {
    let mut checked = 0usize;
    let mut violations: Vec<String> = Vec::new();

    for first_start in 1..=6usize {
        for first_len in 1..=2usize {
            let first_end = first_start + first_len - 1;
            let span = if first_len == 1 {
                format!("{first_start}")
            } else {
                format!("{first_start}_{first_end}")
            };
            // `dup` and `ins` both reach the tract through a junction, and both
            // type as a copy count over it. `del` is the shape the differential
            // route already covered.
            for first in [format!("{span}dup"), format!("{span}del")] {
                for second_pos in first_end + 2..=9usize {
                    for second in [
                        format!("{second_pos}_{}insA", second_pos + 1),
                        format!("{second_pos}T>A"),
                        format!("{second_pos}del"),
                        format!("{second_pos}dup"),
                    ] {
                        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime]
                        {
                            let input = format!("TEMPLATE:g.[{first};{second}]");
                            // An input that denotes no sequence has nothing to
                            // preserve; the generator can emit one at the tract
                            // edges.
                            let Some(want) = apply(TRACT, &input) else {
                                continue;
                            };
                            let output = normalize_in(TRACT, &input, direction);
                            checked += 1;

                            let spans = member_interbase_spans(TRACT, &output);
                            for (i, first) in spans.iter().enumerate() {
                                for second in &spans[i + 1..] {
                                    if claims_same_territory(*first, *second) {
                                        violations.push(format!(
                                            "{input} [{direction:?}] -> {output} \
                                             claims {first:?} and {second:?} at once"
                                        ));
                                    }
                                }
                            }
                            match apply(TRACT, &output) {
                                Some(got) if got == want => {}
                                Some(got) => violations.push(format!(
                                    "{input} [{direction:?}] -> {output} \
                                     (want {want}, got {got})"
                                )),
                                None => violations.push(format!(
                                    "{input} [{direction:?}] -> {output} denotes no sequence"
                                )),
                            }
                        }
                    }
                }
            }
        }
    }

    // A floor, so the enumeration cannot go hollow while the assertion below
    // stays green. The loop bounds are fixed, so this is a structural count and
    // not a measured one.
    assert!(
        checked > 400,
        "the family went hollow: only {checked} cases enumerated"
    );
    assert!(
        violations.is_empty(),
        "rendered members claim one another's territory in {checked} cases: {violations:#?}"
    );
}

/// Do two half-open interbase spans claim the same territory?
///
/// Flush adjacency (`hi == other_lo`) is **not** an overlap — that is the #999
/// adjacency the collapse pass exists to catch, and treating it as one would make
/// this guard reject the correct answers. Three things are:
///
/// * two *zero-length* spans at one interbase, since nothing orders the two
///   payloads (#1286);
/// * a zero-length span strictly inside a non-empty one, which is this defect —
///   a copy count over a tract cannot express another member adding sequence
///   among those bases (#1287);
/// * two non-empty spans that intersect at all.
fn claims_same_territory((lo, hi): (u64, u64), (other_lo, other_hi): (u64, u64)) -> bool {
    match (lo == hi, other_lo == other_hi) {
        (true, true) => lo == other_lo,
        (true, false) => lo > other_lo && lo < other_hi,
        (false, true) => other_lo > lo && other_lo < hi,
        (false, false) => lo < other_hi && other_lo < hi,
    }
}
