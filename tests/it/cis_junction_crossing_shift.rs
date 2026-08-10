//! An insertion or duplication must not carry its junction over a sibling.
//!
//! The sibling-crossing clamp governs edits that consume the bases under their
//! span. An insertion or duplication consumes none — it adds sequence at the
//! junction between two bases — so it is excluded there, deliberately: a member
//! landing flush against one is the adjacency the collapse pass exists to catch
//! (#999, #1135). A duplication is the one partial exception, and only in the
//! 5' direction: the bases under its span are what it copies, so a sibling
//! shifting 5' onto them is bounded (`blocks_sibling_shift`, #1286). Reaching
//! it from the 3' side is still the adjacency above.
//!
//! It can still cross. The junction 3'-shifts through a tract like anything
//! else, and moving it past a base a sibling edits changes what the allele
//! denotes:
//!
//! ```text
//! reference    ACAAAAAAAACGTACGTACG        A-run at 3-10
//! g.[2_3insA;5A>G]  ->  g.[5A>G;10dup]     junction moved 2|3 -> 10|11
//! input applied     ACAAAGAAAAACGTACGTACG
//! output applied    ACAAGAAAAAACGTACGTACG  ← the substituted base moved
//! ```
//!
//! Well-formed, disjoint, warning-free, and wrong — the silent shape.
//!
//! The same defect reaches the repeat path. A `dup` or `ins` inside a tandem
//! tract canonicalises to a copy count over the **whole tract**, and that span
//! can swallow a sibling outright:
//!
//! ```text
//! reference    TTTTTTTTTAATATATTTTA
//! g.[1_2dup;4T>A]  ->  g.[1_9T[11];4T>A]   `4T>A` sits inside `1_9`
//! ```
//!
//! The demotion pass re-spelled only repeats that grew from a *deletion*, so
//! these were left. Both are fixed here: the demotion covers `dup`/`ins`-grown
//! repeats too, and the junction is then bounded at the sibling's 5' edge —
//! still flush against it, so the #999 collapse keeps firing.

use crate::common::cis_apply_oracle::{
    apply, apply_parsed_with, apply_with, assert_normalizes_preserving,
    assert_normalizes_preserving_in, normalize, normalize_in, normalizers_for, provider,
    sweep_seeds, sweep_sequences,
};
use ferro_hgvs::{parse_hgvs, ShuffleDirection};
use rayon::prelude::*;

/// How many offenders a failing sweep names.
///
/// Applied per drawn sequence and again to the concatenation, which is what
/// keeps a parallel sweep's report byte-identical to a serial one: each task
/// keeps its own first `REPORTED`, and truncating the in-order concatenation
/// yields the same set the serial loop's global cap kept.
const REPORTED: usize = 10;

/// What one drawn sequence contributed to a sweep.
///
/// The counters are counts over a partition of the cases, so summing them is
/// exact; the two lists are capped samples for the failure message. Returning
/// one of these rather than mutating shared state is what makes the sequences
/// independent, and [`SweepTally::absorb`] is the one place the fold lives.
///
/// It is the **union** of what the three sweeps in this file tally, so two
/// fields are used by one sweep each and stay zero elsewhere: `residual` by the
/// genomic two-member sweep (its `dup`+`del`+5' shape), and
/// `input_declined`/`by_cause` by the three-member sweep. One struct with two
/// idle fields beats three near-identical ones with three copies of the fold.
#[derive(Default)]
struct SweepTally {
    checked: usize,
    residual: usize,
    skipped: usize,
    input_declined: usize,
    by_cause: std::collections::BTreeMap<&'static str, usize>,
    overlapping: Vec<String>,
    changed: Vec<String>,
}

impl SweepTally {
    /// Fold one sequence's tally in, preserving sweep order.
    ///
    /// The caps are re-applied by the caller after the whole concatenation, so
    /// this appends rather than truncating: each task already kept only its own
    /// first [`REPORTED`], and truncating the in-order concatenation yields the
    /// set the serial loop's global cap kept.
    fn absorb(&mut self, other: SweepTally) {
        self.checked += other.checked;
        self.residual += other.residual;
        self.skipped += other.skipped;
        self.input_declined += other.input_declined;
        for (cause, count) in other.by_cause {
            *self.by_cause.entry(cause).or_default() += count;
        }
        self.overlapping.extend(other.overlapping);
        self.changed.extend(other.changed);
    }
}

/// Nine `T` at positions 1-9 — a tract long enough to canonicalise to a repeat.
const TRACT: &str = "TTTTTTTTTAATATATTTTA";

/// An `A`-run at positions 3-10 — long enough for a junction to travel, short
/// enough that a two-copy insertion does not reach repeat notation.
const RUN: &str = "ACAAAAAAAACGTACGTACG";

/// Sequence-changing cases the exhaustive sweep below finds in the `dup` +
/// `del` + 5'-shuffle shape. **Zero**, and pinned as an exact count so the
/// shape stays measured separately from the rest of the sweep.
///
/// This was 74 — the residual tracked as #1266 — until `blocks_sibling_shift`
/// made a duplication a barrier to a sibling's shift *and* made a duplication's
/// own shift subject to the clamp. A rise means a new sequence-changing defect
/// in that shape specifically.
const FIVE_PRIME_DUP_DEL_SEQUENCE_CHANGES: usize = 0;

// The 5'-shuffle insertion-junction shape — the insertion half of #1267 — used
// to be excluded from this file's sequence-preservation assertion and counted
// instead, as `TRACT_FIVE_PRIME_INSERTION_SEQUENCE_CHANGES_BY_POSITION =
// [(15, 100), (36, 16)]`: 116 cases pinned per sibling position so the shape
// stayed measured while it was broken.
//
// It is fixed, so the exclusion is gone and those cases are asserted like every
// other — a regression there now reports the offending descriptions rather than
// a bare count. The two halves of that 116 needed separate bounds, and measuring
// per-position is what showed the second half existed: bounding the junction
// against base-claiming siblings alone took the map to `{15: 50, 36: 8}`,
// exactly half, because the other half's sibling was a `dup` or an `ins` —
// junction-occupying, so invisible to `claims_bases`. Mirroring the
// junction-vs-junction bound (#1290) closed the rest.

#[test]
fn an_insertion_does_not_shift_past_a_substituted_base() {
    // The junction may travel the `A`-run only as far as the base before the
    // substitution; clamped there it is flush against it, and the two coalesce.
    assert_normalizes_preserving(RUN, "TEMPLATE:g.[2_3insA;5A>G]", "TEMPLATE:g.5_6insG");
}

#[test]
fn a_duplication_does_not_shift_past_a_substituted_base() {
    assert_normalizes_preserving(RUN, "TEMPLATE:g.[4dup;7A>G]", "TEMPLATE:g.7_8insG");
}

#[test]
fn a_repeat_grown_from_a_duplication_does_not_span_a_sibling() {
    // `1_2dup` canonicalised to `1_9T[11]`, whose span covers the substituted
    // base at 4. Demoted back to a duplication, clamped, then coalesced.
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[1_2dup;4T>A]", "TEMPLATE:g.5_6insAT");
}

#[test]
fn a_repeat_grown_from_a_duplication_does_not_span_a_deletion() {
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[1_2dup;4del]", "TEMPLATE:g.9dup");
}

#[test]
fn a_repeat_grown_from_an_insertion_does_not_span_a_sibling() {
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[2_3insTT;5T>A]", "TEMPLATE:g.6_7insAT");
}

#[test]
fn a_lone_duplication_still_reaches_repeat_notation() {
    // Negative control: no sibling, nothing to span, so the tract-wide repeat
    // is exactly right and must survive.
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.1_2dup", "TEMPLATE:g.1_9T[11]");
}

#[test]
fn an_insertion_may_still_land_flush_against_its_sibling() {
    // Negative control for #999, the rule this clamp must not break. The
    // inserted `C` 3'-shifts to a dup at 306, landing at the junction
    // `306|307` — the last position the clamp permits, since `307G>T` claims
    // base 307 and the junction has not passed it. The two then coalesce.
    let mut seq = vec![b'A'; 1000];
    for (i, b) in "CATCCTCGCTCCT".bytes().enumerate() {
        seq[299 + i] = b;
    }
    let seq = String::from_utf8(seq).unwrap();
    assert_eq!(
        normalize(&seq, "TEMPLATE:g.[305_306insC;307G>T]"),
        "TEMPLATE:g.307delinsCT"
    );
}

#[test]
fn an_insertion_reaching_the_end_of_its_run_still_collapses() {
    // Negative control: the junction travels the whole `A`-run and lands flush
    // against `11C>G`, which is outside the run — no crossing, so the two
    // coalesce exactly as before.
    assert_normalizes_preserving(RUN, "TEMPLATE:g.[2_3insA;11C>G]", "TEMPLATE:g.11delinsAG");
}

/// A 5'-shuffled insertion merging with a deletion sibling must not move the
/// base it contributes (#1342).
///
/// Surfaced by giving the sweep's first-member `ins` a payload that is not the
/// span's own reference base. In the `T`-run at 1-9, inserting an `A` at the
/// `2|3` junction and deleting a `T` at 4 is net "one `T` becomes an `A`" — and
/// the `A` lands at position 3. It cannot shuffle 5': an `A` in a `T`-run is
/// unique, so moving it is a different sequence, not a different spelling of
/// the same one. The 3' direction already answered `g.3T>A`.
#[test]
fn an_insertion_merging_with_a_deletion_keeps_its_base_in_place() {
    assert_normalizes_preserving_in(
        TRACT,
        "TEMPLATE:g.[2_3insA;4del]",
        "TEMPLATE:g.3T>A",
        ShuffleDirection::FivePrime,
    );
}

/// Negative controls for the 5' junction barrier: it must bound only the shifts
/// that actually reorder the two edits, never a member that may legitimately
/// travel to the 5' end of its tract.
///
/// This is the failure mode the barrier's own assertion cannot see. A member
/// held back too far still *preserves sequence*, so the exhaustive sweep below
/// would pass while every answer quietly under-shifted. These pin the shifts
/// that must still complete.
#[test]
fn the_five_prime_junction_barrier_does_not_over_clamp() {
    // Payload equals the run base, so the insertion and the deletion cancel:
    // the pair must still reduce all the way to an identity rather than being
    // stranded as two members either side of the junction.
    assert_normalizes_preserving_in(
        TRACT,
        "TEMPLATE:g.[2_3insT;4del]",
        "TEMPLATE:g.2_3=",
        ShuffleDirection::FivePrime,
    );
    // Net +1 `T` in the run. The 5'-most spelling is a duplication at 1, which
    // is 5' of the `2|3` junction — reachable because the merged member no
    // longer straddles it.
    assert_normalizes_preserving_in(
        TRACT,
        "TEMPLATE:g.[2_3insTT;5del]",
        "TEMPLATE:g.1dup",
        ShuffleDirection::FivePrime,
    );
    // The deletion is outside the tract entirely, so there is no sweep and
    // nothing for the barrier to bound.
    //
    // **Re-blessed from `TEMPLATE:g.[2_3insA;12del]`, and it is a SPEC-CONFORMANT
    // move rather than a loss.** With the input-relative weight bound deleted
    // (`rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]`) this
    // allele is re-derived from the sequence instead of being handed back
    // verbatim. The block is `TTTTTTTAAT` -> `ATTTTTTTAA`, **equal length 10/10**,
    // so the column correspondence is unique and columns {0, 7, 9} — positions
    // 3, 10 and 12 — change as a *fact* rather than a choice.
    // `DNA/delins.md:15` types each single changed column as a substitution
    // ("by definition, when **one** nucleotide is replaced by **one** other
    // nucleotide, the change is a substitution") and `:17` describes variants
    // separated by one or more nucleotides individually. `g.[3T>A;10A>T;12T>A]`
    // is that answer verbatim; the old `[2_3insA;12del]` is the input's own
    // spelling, which no clause derives from the sequence.
    //
    // **This row no longer bounds the barrier**, and that is measured rather
    // than assumed: deleting `.filter(|_| a.junction.is_none())` from the 5'
    // branch of `clamp_sibling_crossing_shifts` leaves this answer
    // byte-identical, because the sequence-first derivation answers the allele
    // and never consults the per-member clamp. The two rows above still bound
    // it. See the note on
    // `a_third_member_clear_of_the_tract_is_re_derived_from_the_sequence`,
    // which lost the same bite at the same time; the junction filter's guard is
    // now
    // `a_third_member_past_the_derivation_window_keeps_the_duplication_reaching_its_five_prime_most_position`
    // (#1603).
    //
    // Both spellings are asserted, and the authored one is the load-bearing
    // half. The paragraph above claims a MOVE — `[2_3insA;12del]` is re-derived
    // onto the three substitutions — and pinning only the output's
    // fixed-pointedness leaves that claim unobservable: a change that stopped
    // re-deriving the authored form would hand it straight back and this row
    // would stay green while its own doc had become false.
    assert_normalizes_preserving_in(
        TRACT,
        "TEMPLATE:g.[2_3insA;12del]",
        "TEMPLATE:g.[3T>A;10A>T;12T>A]",
        ShuffleDirection::FivePrime,
    );
    assert_normalizes_preserving_in(
        TRACT,
        "TEMPLATE:g.[3T>A;10A>T;12T>A]",
        "TEMPLATE:g.[3T>A;10A>T;12T>A]",
        ShuffleDirection::FivePrime,
    );
}

/// The tract `[4_5insC;5_6dup]` shuffles through: an `AT`-alternating run at
/// 3-11, long enough for a two-base duplication to travel and short enough that
/// a third member placed past 11 sits clear of it.
const DUP_RUN: &str = "TTAATATATAATAATATTAT";

/// The pair alone. Both members land inside one run of change, so #1235's
/// sequence-first derivation answers the whole allele from the bases it denotes
/// and emits a single insertion.
///
/// **Re-blessed from `TEMPLATE:g.[4_5insC;4_5dup]`.** The merged form is derived
/// from the sequence rather than assembled per member; `[4_5insC;4_5dup]` was
/// only ever the per-member repair's spelling of the same bases. Verified
/// sequence-preserving: `assert_normalizes_preserving_in` asserts the string
/// *before* it applies both descriptions, so the failing string assertion left
/// the preservation check unrun — with the string re-blessed it executes, and
/// `apply(DUP_RUN, …)` gives `TTAACTATATATAATAATATTAT` for input and output
/// alike.
///
/// **This row no longer bounds anything.** It was #1357's negative control for
/// the over-clamp guard on the 5' junction barrier, and it passes now whether
/// that guard runs or not — the output contains no duplication at all, so the
/// property in its old name (`…_still_reaches_its_five_prime_most_position`) is
/// not exercised here. That coverage moved to
/// [`a_third_member_past_the_derivation_window_keeps_the_duplication_reaching_its_five_prime_most_position`],
/// which is the row that still goes red when the guard is removed (#1603).
///
/// It is **not**
/// [`a_third_member_clear_of_the_tract_is_re_derived_from_the_sequence`], which
/// this line used to name: that row lost its bite to the same deletion, and its
/// own doc records the mutation leaving it byte-identical. Pointing here at a
/// row that no longer discriminates is how a filter comes to look guarded when
/// nothing guards it.
#[test]
fn a_multi_base_duplication_beside_an_insertion_merges_from_the_sequence() {
    assert_normalizes_preserving_in(
        DUP_RUN,
        "TEMPLATE:g.[4_5insC;5_6dup]",
        "TEMPLATE:g.3_4insACT",
        ShuffleDirection::FivePrime,
    );
}

/// #1357's discriminator, restored on the nucleotide axis.
///
/// The barrier bounds a member by its **span start**, which is the right edge to
/// bound only for a member that claims bases. A `dup` carries its own junction at
/// its *end*, so bounding its start would hold it back by `len - 1` more than the
/// invariant requires — which is why the 5' barrier carries
/// `.filter(|_| a.junction.is_none())`.
///
/// `3_4dup` really is illegal here — it would share interbase 4 with the
/// insertion — so `4_5dup` is the 5'-most legal spelling, and the member must
/// still reach it. Without the guard it is stranded at `5_6dup`, its input
/// position: sequence-preserving, so neither exhaustive sweep can see it, and
/// silently non-canonical. The sweeps are blind to this class twice over — their
/// only oracle is `apply`, and their second member is always a single-base `dup`,
/// exactly the length at which `len - 1` is zero and the over-clamp disappears.
///
/// The third member is what keeps the property visible. `15del` stops the
/// *merge* without touching the shift: the block is no longer one run of change,
/// so the derivation declines to collapse the allele and the per-member pipeline
/// decides the output again, which is where the barrier lives. The `dup` then
/// reaches `4_5dup` and the spelling is observable.
///
/// **Distance measured, not assumed.** Sweeping every third-member position on
/// `DUP_RUN`: at `10del` and `11del` the member joins the block and the allele
/// merges again (to `g.[4_5insC;10_11insT]`) — a test written there would be
/// vacuous in a new way, pinning a merged form while claiming to pin a barrier.
/// `12del` and `13del` are clear but only 1-2 nt past the tract; at `14del` and
/// `18del` the third member itself 5'-shifts, which muddies the row. `15del` is
/// 4 nt clear, stays exactly where it is written, and reaches the tie.
///
/// # It stopped bounding anything, and that is measured
///
/// It used to **fail both ways**: deleting `.filter(|_| a.junction.is_none())`
/// from the 5' branch of `clamp_sibling_crossing_shifts` turned the output into
/// `g.[4_5insC;5_6dup;15del]`, and restoring it turned it back. That is no
/// longer true. With the input-relative weight bound deleted
/// (`rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]`) the third
/// member no longer keeps the allele out of the sequence-first derivation: the
/// whole thing is re-derived, the per-member pipeline never runs, and the clamp
/// is not consulted. Applying that same mutation now leaves this row
/// **byte-identical**.
///
/// The coverage that cost was real: running the mutation against the whole
/// suite added **zero** failures, so this row and the `[2_3insA;12del]` row in
/// [`the_five_prime_junction_barrier_does_not_over_clamp`] were the filter's
/// only two guards and both lost their bite for the same reason. It is restored
/// by
/// [`a_third_member_past_the_derivation_window_keeps_the_duplication_reaching_its_five_prime_most_position`],
/// which reaches the same property through a decline the deleted bound had
/// nothing to do with (#1603). This row is kept as the *characterisation* of how
/// the shape reads once the derivation answers it. (The *other* half of the
/// barrier — removing the whole 5' `across_junctions` branch — was never guarded
/// here either; that is
/// `an_insertion_merging_with_a_deletion_keeps_its_base_in_place`, which still
/// bites.)
///
/// The re-blessed output itself is **SPEC-SILENT**: the block is
/// `TATATAATAAT` -> `CTATATATAATAA`, unequal 11/13, so no column correspondence
/// exists and `delins.md:15`/`:16`/`:17` have no defined input. The `dup` that
/// disappears is not a `duplication.md:18` violation — see the ruling record for
/// why `:18` ranks type labels for one span rather than requiring a partition
/// that manufactures a `dup` member.
///
/// Renamed with the re-bless, on the precedent
/// `issue_1304_junction_barrier_snapshot.rs` sets at its own line 50: the old
/// name promised a `dup` reaching its 5'-most position and the form this pins
/// carries no `dup` at all, so the name stated the opposite of the assertion. A
/// failure message and a `grep` both have to land on the guard that actually
/// broke.
#[test]
fn a_third_member_clear_of_the_tract_is_re_derived_from_the_sequence() {
    assert_normalizes_preserving_in(
        DUP_RUN,
        "TEMPLATE:g.[4_5insC;5_6dup;15del]",
        "TEMPLATE:g.[5_10delinsCTATAT;16A[3]]",
        ShuffleDirection::FivePrime,
    );
}

/// ADJUDICATED 2026-08-09 — one variant, two stable normalized strings.
///
/// # What is measured
///
/// The row above and the row here denote the **same 22-mer**,
/// `TTAACTATATATAATAAATTAT`, verified by applying each spelling to `DUP_RUN`
/// through `hgvs_to_spdi` independently of the normalizer (the `apply`
/// assertions below). Yet each direction lands on a stable answer and does not
/// move off it:
///
/// ```text
/// 5'  g.[4_5insC;5_6dup;15del] -> g.[4_5insC;4_5dup;15del]   (output is stable)
/// 5'  g.[4_5insCTA;15del]      -> g.[3_4insACT;15del]        (output is stable)
/// 5'  g.[3_4insACT;15del]      -> g.[3_4insACT;15del]        (fixed point)
/// 3'  g.[4_5insC;5_6dup;15del] -> g.[4_5insC;9_10dup;15del]  (output is stable)
/// 3'  g.[3_4insACT;15del]      -> g.[4_5insCTA;15del]        (output is stable)
/// 3'  g.[4_5insCTA;15del]      -> g.[4_5insCTA;15del]        (fixed point)
/// ```
///
/// That is the #1235 shape — two stable fixed points for one variant — inside
/// this module's own fixture, in **both** directions. It is not caused by the
/// junction barrier: the barrier decides *where* the `dup` lands, not whether
/// the allele stays three members.
///
/// # Why it is three members at all
///
/// Only because of the distant `15del`. Without it the pair merges:
/// `g.[4_5insC;5_6dup]` -> `g.3_4insACT` at 5' and `g.4_5insCTA` at 3', which
/// `a_multi_base_duplication_beside_an_insertion_merges_from_the_sequence`
/// already pins. The third member breaks the block into more than one run of
/// change, the sequence-first derivation declines, and the per-member pipeline
/// hands back the input's own partition.
///
/// # The ruling, and the clause it turns on
///
/// The locus carries **one** change: aligning `DUP_RUN` against the denoted
/// 22-mer leaves a single contiguous 3 nt insertion (`CTA` at junction 4|5,
/// which 5'-rolls to `ACT` at 3|4). `general.md:34` —
///
/// > two variants separated by one or more nucleotides should be described
/// > individually and **not** as a "delins"
///
/// is stated over *two variants*, so it does not license splitting that one
/// insertion into an `ins` at junction 4|5 and a `dup` inserting at 5|6. The
/// separation of 1 that the three-member form presents is a property of the
/// spelling, not of the variant — the case
/// `separation-is-a-property-of-the-spelling-not-of-the-variant` records — and
/// it exists here only because the `AT` tract lets the three inserted bases be
/// dealt out on either side of an unchanged base.
///
/// `DNA/duplication.md:18` ("when a variant can be described as a duplication,
/// it **must** be described as a duplication") does not rescue the `dup`
/// either: the variant is the 3 nt insertion, and neither `ACT` nor `CTA` is a
/// copy of the reference bases it abuts.
///
/// So the adjudicated-canonical form is the **re-derived one-member** insertion
/// — `g.[3_4insACT;15del]` at 5', `g.[4_5insCTA;15del]` at 3' — which is what
/// the operator ruling `canonical-form-choice-when-both-legal` (decided
/// 2026-08-07: derive from the resulting sequence, do not preserve the input's
/// spelling) selects, and which ferro already emits for the equal-denoting
/// one-member spelling.
///
/// # Recorded as a gap, not asserted
///
/// Ferro does **not** do this today, so this row pins both fixed points and the
/// fact that they are one variant. It is deliberately not written as
/// `assert_eq!(three_member_output, one_member_output)`: that would be a red
/// test, and closing it is a normalizer change with a representation-change
/// declaration attached, not a test edit. What this row makes impossible is
/// quietly re-blessing either string as *the* canonical one — see the ruling
/// record `contiguous-insertion-split-by-a-blocked-derivation`.
///
/// A sweep document dated 2026-08-08 reached the opposite verdict, deriving
/// `g.[4_5insC;4_5dup;15del]` as correct from the observation that the two
/// members are separation 1 rather than separation 0. The arithmetic is right
/// and the conclusion does not follow: it evaluates `general.md:34` on the
/// spelling, which is the very thing the sibling record says cannot be done.
///
/// # RED ON THIS BRANCH, DELIBERATELY — DO NOT TAKE THE RE-BLESS THE MESSAGE OFFERS
///
/// Deleting the input-relative weight bound
/// (`rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]`) removed the
/// decline that kept the three-member spelling out of the sequence-first
/// derivation, so the first pin below now reads
/// `g.[5_10delinsCTATAT;16A[3]]` rather than `g.[4_5insC;4_5dup;15del]`.
///
/// The `assert_ne!` further down says "this gap is closed, so re-bless the four
/// pins above". **It is not closed.** `contiguous-insertion-split-by-a-blocked-derivation`
/// decides for the re-derived ONE-MEMBER insertion — `g.[3_4insACT;15del]` at 5'
/// and `g.[4_5insCTA;15del]` at 3' — and the branch reaches neither. The output
/// moved off the recorded gap onto a THIRD form, which is an adjudication rather
/// than a re-bless, and it is the same frameless merge the ruling record's
/// `THE FRAMELESS RESIDUAL IS A BUG` paragraph says must be closed rather than
/// declared. Leave this red until it is.
#[test]
fn the_three_member_spelling_and_its_one_member_form_are_two_fixed_points() {
    const THREE_MEMBER: &str = "TEMPLATE:g.[4_5insC;5_6dup;15del]";
    const ONE_MEMBER_FIVE_PRIME: &str = "TEMPLATE:g.[3_4insACT;15del]";
    const ONE_MEMBER_THREE_PRIME: &str = "TEMPLATE:g.[4_5insCTA;15del]";

    // The premise. If these ever stopped denoting one sequence the rest of this
    // row would be comparing two different variants, which are *supposed* to
    // normalize apart.
    let denoted = apply(DUP_RUN, THREE_MEMBER).expect("three-member spelling applies");
    for spelling in [ONE_MEMBER_FIVE_PRIME, ONE_MEMBER_THREE_PRIME] {
        assert_eq!(
            apply(DUP_RUN, spelling).as_deref(),
            Some(denoted.as_str()),
            "`{spelling}` and `{THREE_MEMBER}` must denote one sequence",
        );
    }

    // Both are fixed points, in the direction that is theirs.
    //
    // The two `ONE_MEMBER_*` rows assert that directly — the output IS the
    // input. The two `THREE_MEMBER` rows cannot: their output is a different
    // string from their input, so pinning it says only where the first pass
    // landed. The claim this row is making is that the landing place is
    // *stable*, so each pinned output is normalized a second time in the same
    // direction and required to be byte-identical. Without that, a
    // non-idempotent regression on either output keeps this row green while the
    // "two stable normalized strings" headline above becomes false.
    //
    // All SIX rows of the table above are here, including the two cross-direction
    // ones — each one-member spelling read in the *other* direction's arm. Those
    // two are not decoration: the divergence loop below compares against
    // `ONE_MEMBER_FIVE_PRIME` in both directions, which is only the one-member
    // answer under 3' because that spelling rolls onto `ONE_MEMBER_THREE_PRIME`
    // there. An unasserted row that another assertion leans on is the shape this
    // whole branch is about.
    for (input, expected, direction) in [
        (
            THREE_MEMBER,
            "TEMPLATE:g.[4_5insC;4_5dup;15del]",
            ShuffleDirection::FivePrime,
        ),
        (
            ONE_MEMBER_FIVE_PRIME,
            ONE_MEMBER_FIVE_PRIME,
            ShuffleDirection::FivePrime,
        ),
        (
            ONE_MEMBER_THREE_PRIME,
            ONE_MEMBER_FIVE_PRIME,
            ShuffleDirection::FivePrime,
        ),
        (
            THREE_MEMBER,
            "TEMPLATE:g.[4_5insC;9_10dup;15del]",
            ShuffleDirection::ThreePrime,
        ),
        (
            ONE_MEMBER_THREE_PRIME,
            ONE_MEMBER_THREE_PRIME,
            ShuffleDirection::ThreePrime,
        ),
        (
            ONE_MEMBER_FIVE_PRIME,
            ONE_MEMBER_THREE_PRIME,
            ShuffleDirection::ThreePrime,
        ),
    ] {
        let once = normalize_in(DUP_RUN, input, direction);
        assert_eq!(
            once, expected,
            "{direction:?}: `{input}` no longer normalizes to `{expected}`",
        );
        let twice = normalize_in(DUP_RUN, &once, direction);
        assert_eq!(
            twice, once,
            "{direction:?}: `{input}` -> `{once}` is not a fixed point (it \
             normalizes on to `{twice}`), so this row's claim that the variant \
             has two STABLE forms does not hold",
        );
    }

    // And they are two, not one. Stated as its own assertion so the divergence
    // is what fails if a fix lands, rather than four string pins failing for
    // four reasons.
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        let from_three = normalize_in(DUP_RUN, THREE_MEMBER, direction);
        let from_one = normalize_in(DUP_RUN, ONE_MEMBER_FIVE_PRIME, direction);
        assert_ne!(
            from_three, from_one,
            "{direction:?}: the two spellings converged — this gap is closed, so \
             re-bless the four pins above and say which form moved \
             (`contiguous-insertion-split-by-a-blocked-derivation`)",
        );
    }

    // The one thing that must never be true of either answer.
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        for spelling in [THREE_MEMBER, ONE_MEMBER_FIVE_PRIME, ONE_MEMBER_THREE_PRIME] {
            let output = normalize_in(DUP_RUN, spelling, direction);
            assert_eq!(
                apply(DUP_RUN, &output).as_deref(),
                Some(denoted.as_str()),
                "{direction:?}: `{spelling}` -> `{output}` changed the sequence",
            );
        }
    }
}

/// Where a third member has to sit for the sequence-first derivation to decline
/// on width alone.
///
/// `canonicalize_from_sequence` pads the member-span union by `CANONICAL_PAD`
/// (128) at each end and refuses a window wider than `MAX_CANONICAL_WINDOW`
/// (4096). The **pair alone** pads to `1..=134`, comfortably inside the limit —
/// so the pair is not what declines. A member here widens the span union to
/// 4-4200, whose padded window is `1..=4328` and therefore over the limit. That
/// is what takes the allele out of the derivation, without depending on anything
/// the derivation *decides* — which is the property the guard below needs. Both
/// figures are stated because the next reader will use this arithmetic to pick a
/// new constant.
const PAST_THE_DERIVATION_WINDOW: usize = 4200;

/// [`DUP_RUN`] followed by a plain `ACGT` filler long enough to carry a member
/// at [`PAST_THE_DERIVATION_WINDOW`].
///
/// The filler is 4-periodic on purpose: it holds no homopolymer or dinucleotide
/// tract a member could shuffle through, so the third member stays exactly where
/// it is written and contributes nothing to the case but its distance.
fn dup_run_past_the_derivation_window() -> String {
    let mut seq = String::from(DUP_RUN);
    while seq.len() < PAST_THE_DERIVATION_WINDOW + 100 {
        seq.push_str("ACGT");
    }
    seq
}

/// #1603 — the 5' junction filter, guarded again.
///
/// [`a_third_member_clear_of_the_tract_is_re_derived_from_the_sequence`]
/// used to be the discriminator for `.filter(|_| a.junction.is_none())` in the 5'
/// branch of `clamp_sibling_crossing_shifts`, and stopped being one when the
/// input-relative weight bound was deleted: its third member no longer keeps the
/// allele out of the sequence-first derivation, so the per-member pipeline — the
/// only place the filter runs — never decides the output. Measured at the time:
/// deleting that line left the **whole suite** green.
///
/// The property is unchanged and is worth restating, because it is invisible to
/// every other check here. The barrier bounds a member by its **span start**,
/// which is the right edge only for a member that claims bases. A `dup` carries
/// its junction at its *end* (`junction_of`), so bounding its start holds it back
/// by `len - 1` more than the invariant requires. `3_4dup` really is illegal —
/// it would share interbase 4 with the insertion — so `4_5dup` is the 5'-most
/// legal spelling and the member must reach it. Without the filter it is
/// stranded at `5_6dup`, its input position: sequence-preserving, so the
/// exhaustive sweeps below cannot see it, and silently non-canonical. The sweeps
/// are blind to the class twice over — their only oracle is `apply`, and their
/// second member is always a single-base `dup`, exactly the length at which
/// `len - 1` is zero and the over-clamp disappears.
///
/// What is new is the *decline*. A third member past
/// [`PAST_THE_DERIVATION_WINDOW`] refuses the derivation on window width, which
/// is a property of the members' geometry rather than of what the derivation
/// would have produced — so, unlike the old `15del`, it survives every change to
/// the rules governing the derived pieces, which is what re-admitted that row.
///
/// State that claim at its actual width, because the wider one is false and was
/// written here once: `CANONICAL_PAD` and `MAX_CANONICAL_WINDOW` **are** the
/// derivation's own bounds, so raising either re-admits this row exactly as
/// raising the weight bound re-admitted the old `15del`. The decline is durable
/// against what the derivation *produces*, not against how far it is allowed to
/// look. The 3' row is the control: it pins
/// the same allele reaching `9_10dup`, which is what shows the derivation really
/// declined (the pair alone merges to a single member — see
/// [`a_multi_base_duplication_beside_an_insertion_merges_from_the_sequence`])
/// and that the filter is 5'-only.
///
/// Verified both ways on this tree: with the filter deleted the 5' row emits
/// `g.[4_5insC;5_6dup;4200T>A]` and this test fails, while the 3' row is
/// byte-identical.
#[test]
fn a_third_member_past_the_derivation_window_keeps_the_duplication_reaching_its_five_prime_most_position(
) {
    let seq = dup_run_past_the_derivation_window();
    // Fixture sanity: the third member's stated reference base is part of the
    // case, and a mis-stated one is refused upstream by
    // `stated_reference_bases_match` — which would decline the derivation for
    // the wrong reason and make this row pass whatever the clamp does.
    assert_eq!(
        seq.as_bytes()[PAST_THE_DERIVATION_WINDOW - 1],
        b'T',
        "the ACGT filler must put a T at {PAST_THE_DERIVATION_WINDOW}"
    );
    let input = format!("TEMPLATE:g.[4_5insC;5_6dup;{PAST_THE_DERIVATION_WINDOW}T>A]");

    assert_normalizes_preserving_in(
        &seq,
        &input,
        &format!("TEMPLATE:g.[4_5insC;4_5dup;{PAST_THE_DERIVATION_WINDOW}T>A]"),
        ShuffleDirection::FivePrime,
    );
    assert_normalizes_preserving_in(
        &seq,
        &input,
        &format!("TEMPLATE:g.[4_5insC;9_10dup;{PAST_THE_DERIVATION_WINDOW}T>A]"),
        ShuffleDirection::ThreePrime,
    );
}

#[test]
fn no_two_member_allele_normalizes_to_a_different_sequence() {
    // The widened class. The sweep in `repeat_span_sibling_overlap.rs` uses a
    // deletion as the first member; this one adds duplications and insertions,
    // which reach the same tracts through the junction rather than through a
    // consumed span.
    //
    // Two failure modes are counted separately: an *overlapping* output, which
    // a consumer can at least reject, and a *sequence-changing* one, which is
    // silent. Before this fix the same sweep reported 1,318 of the first and
    // 1,033 of the second; it now reports 0 and 0.
    //
    // The `dup` + `del` + 5' shape is still *tallied* on its own, pinned at
    // zero by `FIVE_PRIME_DUP_DEL_SEQUENCE_CHANGES` and characterised in
    // `a_five_prime_duplication_beside_a_deletion_keeps_its_sequence`. It was
    // the last residual class here — 74 cases, closed by `blocks_sibling_shift`
    // — and keeping the separate tally says which shape regressed if one does.

    // 24 seeds, not 48. Widening the second member from 2 shapes to 4 (#1283)
    // doubles the work per sequence, and this is already the suite's slowest
    // test — measured 45s before, 117s at 48 seeds. Halving the seed count
    // holds the runtime and the *case count* at roughly what they were while
    // covering twice the shapes, which is the trade the issue asks for: its
    // five gaps are all shape gaps, not sequence-diversity gaps, and every
    // blocking defect found so far lived in a shape the generator could not
    // emit rather than in a sequence it did not draw.
    //
    // One cost to keep on the record: `FIVE_PRIME_DUP_DEL_SEQUENCE_CHANGES`
    // below is an exact count over *this* corpus, so halving the seeds halves
    // the sequence diversity backing that pinned zero. It stays zero here, and
    // a defect would have to live only in a dropped sequence *and* in the
    // dup+del+5' shape to escape — but the guard is weaker than it was, and
    // #1295's seed knob is what restores it without paying the runtime.
    //
    // #1295 delivered that knob, so the full count goes back to 48 and the
    // diversity behind the pinned zero is restored. `sweep_seeds` returns the
    // full 48 when CI asks (`FERRO_SWEEP_SEEDS=full`) and a 4-seed prefix
    // otherwise, so a local run no longer pays for it: this test was 79.6s of an
    // 86.6s local suite at 24 seeds, and the prefix is a strict subset of the
    // corpus, not a different one.
    let seeds = sweep_seeds(48);

    // One rayon task per drawn sequence. Each already builds its own `template`
    // provider and its own normalizers — that is what the comment inside is
    // about — so the sequences share nothing, and this changes the ORDER cases
    // run in, never which cases run.
    //
    // The totals are byte-identical to the serial sweep, which the exactly
    // pinned `FIVE_PRIME_DUP_DEL_SEQUENCE_CHANGES` and the two bounds below
    // require. `checked`, `residual` and `skipped` count over a partition of the
    // cases, so summing is exact; `overlapping` and `changed` are "the first
    // `REPORTED` offenders in sweep order", reproduced by capping per sequence
    // and truncating the in-order concatenation — `collect` preserves input
    // order however the tasks finish.
    let per_sequence: Vec<SweepTally> = sweep_sequences(seeds)
        .into_par_iter()
        .map(|seq| {
            let mut checked = 0usize;
            let mut residual = 0usize;
            let mut skipped = 0usize;
            let mut overlapping: Vec<String> = Vec::new();
            let mut changed: Vec<String> = Vec::new();
        // Built once per sequence rather than once per case. `normalize_in` and
        // `apply` each construct a `TEMPLATE` provider internally, so the previous
        // shape built three per case — the normalizer's, the input apply's and the
        // output apply's — over half a million cases.
        let template = provider(&seq);
        let normalizers = normalizers_for(&seq);
        for first_start in 2..=13usize {
            for first_len in 1..=2usize {
                let first_end = first_start + first_len - 1;
                let span = if first_len == 1 {
                    format!("{first_start}")
                } else {
                    format!("{first_start}_{first_end}")
                };
                // #1342: the payload must not equal either reference base
                // flanking the insertion point, for the same reason the second
                // member's payload excludes both of its neighbours below.
                // This entry used to insert `seq[first_start - 1..first_end]`,
                // the span's own reference bases — so for `first_len == 1` it
                // denoted exactly `{first_start}dup`, the entry directly above
                // it. Half the `first_len` values contributed two distinct
                // first-member shapes while the array read as three, and the
                // single-base span is the case where a junction has the least
                // context to disambiguate it.
                //
                // `first_start <= 13` and `sweep_sequences` yields 20-base
                // sequences, so index `first_start` — the 3' neighbour — is
                // always in bounds.
                let first_base = seq.as_bytes()[first_start - 1] as char;
                let first_next_base = seq.as_bytes()[first_start] as char;
                let first_ins_payload = ['A', 'C', 'G', 'T']
                    .into_iter()
                    .find(|b| *b != first_base && *b != first_next_base)
                    .expect("four bases, at most two excluded");
                let inserted: String = std::iter::repeat_n(first_ins_payload, first_len).collect();
                let firsts = [
                    format!("{span}del"),
                    format!("{span}dup"),
                    format!("{first_start}_{}ins{inserted}", first_start + 1),
                ];
                for first in firsts {
                    for second_start in first_end + 2..=19usize {
                        let base = seq.as_bytes()[second_start - 1] as char;
                        let alt = if base == 'A' { 'G' } else { 'A' };
                        // #1283 gap 1: the second member used to be only
                        // `{del, sub}`, so no case in ~276,000 ever placed a
                        // `dup` or an `ins` downstream of the first member —
                        // a whole class of sibling interactions (#1279's
                        // shape) was unreachable while the file read as
                        // exhaustive.
                        //
                        // The payload is not a reference base adjacent to the
                        // junction, which the first attempt at this list got
                        // wrong twice over. Inserting the base immediately 5' of
                        // the junction denotes the same edit as
                        // `{second_start}dup` — the entry right above it — and
                        // inserting the base immediately 3' of it denotes
                        // `{second_start + 1}dup`. Either way the list reads as
                        // four shapes while carrying three, which is precisely
                        // the gap-1 defect this block exists to remove. So the
                        // payload excludes *both* neighbours rather than only
                        // the 5' one.
                        //
                        // (`firsts` builds its own payload the dup-equivalent
                        // way for `first_len == 1`. That predates this PR and is
                        // tracked as #1342 rather than widened into it.)
                        //
                        // `second_start <= 19` and `sweep_sequences` yields
                        // 20-base sequences (its argument is the *seed* count,
                        // not a length), so index `second_start` — the 3'
                        // neighbour — is at worst the last valid byte.
                        let next_base = seq.as_bytes()[second_start] as char;
                        let second_ins_payload = ['A', 'C', 'G', 'T']
                            .into_iter()
                            .find(|b| *b != base && *b != next_base)
                            .expect("four bases, at most two excluded");
                        for second in [
                            format!("{second_start}del"),
                            format!("{second_start}{base}>{alt}"),
                            format!("{second_start}dup"),
                            format!("{second_start}_{}ins{second_ins_payload}", second_start + 1),
                        ] {
                            // Direction-invariant, so hoisted out of the loop
                            // below: the description, its parse, and the bases it
                            // denotes are all properties of the input rather than
                            // of the shuffle direction.
                            let input = format!("TEMPLATE:g.[{first};{second}]");
                            let variant = parse_hgvs(&input).expect("generated input parses");
                            let want = apply_parsed_with(&template, &seq, &variant);
                            for (direction, normalizer) in &normalizers {
                                let output = normalizer
                                    .normalize(&variant)
                                    .expect("normalize")
                                    .to_string();
                                checked += 1;
                                let Some(want) = want.as_deref() else {
                                    skipped += 1;
                                    continue;
                                };
                                // The dup + del + 5' shape is *counted*, not
                                // excluded. Skipping it would drop roughly
                                // a twelfth of the sweep — some 23,000
                                // cases — out of the sequence-preservation
                                // check, and a new defect anywhere in that
                                // shape would then pass silently. Counting
                                // keeps it measured: the tally is pinned at
                                // zero, so a regression fails just as the
                                // main assertion would, while still naming
                                // the shape it landed in.
                                let residual_shape = first.ends_with("dup")
                                    && second.ends_with("del")
                                    && *direction == ShuffleDirection::FivePrime;
                                // The *output* keeps going through the string-taking
                                // oracle: re-parsing what the normalizer produced is
                                // part of what this sweep asserts.
                                match apply_with(&template, &seq, &output) {
                                    None if overlapping.len() < REPORTED => {
                                        overlapping.push(format!("{seq}: {input} -> {output}"))
                                    }
                                    None => {}
                                    Some(got) if got != want && residual_shape => {
                                        residual += 1;
                                    }
                                    Some(got) if got != want && changed.len() < REPORTED => {
                                        changed.push(format!(
                                            "{seq}: {input} [{direction:?}] -> {output} (want {want}, got {got})"
                                        ));
                                    }
                                    Some(_) => {}
                                }
                            }
                        }
                    }
                }
            }
        }
            SweepTally {
                checked,
                residual,
                skipped,
                overlapping,
                changed,
                ..SweepTally::default()
            }
        })
        .collect();

    let mut total = SweepTally::default();
    for tally in per_sequence {
        total.absorb(tally);
    }
    total.overlapping.truncate(REPORTED);
    total.changed.truncate(REPORTED);
    let SweepTally {
        checked,
        residual,
        overlapping,
        changed,
        skipped,
        ..
    } = total;

    // Per seed, not absolute (#1295): the seed count is now a knob, so a fixed
    // floor would either fail at the default prefix or be vacuous at the full
    // corpus. Measured at 11,520 cases per seed — the loop bounds below are what
    // fix that, and they do not depend on the sequences drawn — so this floor
    // sits deliberately loose, at roughly a third, and guards against the sweep
    // being hollowed out by a lost loop rather than against exact drift.
    const CASES_PER_SEED_FLOOR: usize = 4_000;
    let floor = CASES_PER_SEED_FLOOR * seeds as usize;
    assert!(
        checked > floor,
        "sweep covered too little: {checked} cases over {seeds} seeds (floor {floor})"
    );
    // `checked` counts enumerated cases, but a case whose *input* will not
    // apply contributes nothing to either property below. Bound that share, so
    // the sweep cannot go quietly hollow while still clearing the floor.
    assert!(
        skipped * 10 < checked,
        "too many cases skipped as unconvertible: {skipped} of {checked}"
    );
    assert!(
        overlapping.is_empty(),
        "overlapping or unconvertible output in {checked} cases: {overlapping:#?}"
    );
    assert!(
        changed.is_empty(),
        "sequence-changing normalizations in {checked} cases: {changed:#?}"
    );
    // Pinned exactly at zero. Any rise is a new sequence-changing defect in the
    // dup + del + 5' shape; there is no longer a tolerated residual class here.
    assert_eq!(
        residual, FIVE_PRIME_DUP_DEL_SEQUENCE_CHANGES,
        "the dup + del + 5' shape changed"
    );
}

#[test]
fn a_five_prime_duplication_beside_a_deletion_keeps_its_sequence() {
    // #1266. Under 5' shuffle a duplication beside a deletion used to cancel
    // wrongly, emitting a shorter duplication in the wrong place:
    //
    //   reference        TTATTTAAATAAAAATAAAA
    //   g.[6_7dup;9del]  applies to  TTATTTATAATAAAAATAAAA
    //   normalized to    g.4dup      TTATTTTAAATAAAAATAAAA   <- different bases
    //
    // This was the whole 74-case residual of the sweep above. The deletion was
    // free to shift 5' through the duplicated span because a `Duplication`
    // reports claiming no reference bases — yet the bases under its span are
    // exactly what it copies. `blocks_sibling_shift` makes it a barrier.
    let seq = "TTATTTAAATAAAAATAAAA";
    assert_normalizes_preserving_in(
        seq,
        "TEMPLATE:g.[6_7dup;9del]",
        "TEMPLATE:g.7_8insT",
        ShuffleDirection::FivePrime,
    );
}

#[test]
fn a_five_prime_duplication_does_not_cross_an_upstream_sibling() {
    // #1267, the part `blocks_sibling_shift` reaches. Under 5' shuffle a
    // duplication's span travelled toward an **upstream** sibling and crossed
    // it, moving the base the sibling had edited:
    //
    //   reference        ACAAAAAAAACGTACGTACG        A-run at 3-10
    //   g.[4A>G;9dup]    applies to  ACAGAAAAAAACGTACGTACG
    //   normalized to    g.4_5insG   ACAAGAAAAAACGTACGTACG   <- moved base
    //
    // The sweep above cannot reach this: it always places the sibling
    // *downstream*, so under 5' shuffle the member travels away from it. These
    // were found by hand and stand in for the sweep until its generator is
    // widened to upstream siblings (#1283).
    //
    // The expectations are the *merged* forms since #1235: the pair is now
    // derived from the sequence it denotes rather than assembled per member, and
    // a substitution abutting an A-run duplication is one inserted `G`. The
    // property under test is unchanged — the sequence assertion below is what
    // catches a crossing — and `g.[4A>G;5dup]` was only ever the per-member
    // repair's spelling of the same bases.
    let seq = "ACAAAAAAAACGTACGTACG";
    for (input, expected) in [
        ("TEMPLATE:g.[4A>G;9dup]", "TEMPLATE:g.3_4insG"),
        ("TEMPLATE:g.[5A>G;10dup]", "TEMPLATE:g.4_5insG"),
    ] {
        assert_normalizes_preserving_in(seq, input, expected, ShuffleDirection::FivePrime);
    }
}

#[test]
fn an_insertion_junction_does_not_cross_an_upstream_sibling() {
    // The rest of #1267: the member written as an insertion. Under 5' shuffle
    // its payload travelled 5' past an upstream sibling, so the allele denoted
    // different bases while staying well-formed, disjoint and warning-free:
    //
    //   reference             ACAAAAAAAACGTACGTACG        A-run at 3-10
    //   g.[4A>G;9_10insA]     applies to  ACAGAAAAAAACGTACGTACG
    //   emitted  g.4_5insG                ACAAGAAAAAACGTACGTACG   <- moved base
    //
    // `blocks_sibling_shift` does not reach it, and neither does the span clamp:
    // the member is an `Insertion` on the way in but canonicalises to a
    // `Duplication` (`g.9_10insA` alone -> `g.3dup`), so it *grew*, and a member
    // that grew is refused there for reasons #1266/#1279 measured. What bounds
    // it is the junction clamp's 5' half, which governs **every** junction mover
    // regardless of how the input spelled it — restricting it to inputs written
    // as insertions was tried and refuted by measurement (see
    // `clamp_sibling_crossing_junctions`: over 73,376 duplication-mover cases the
    // unrestricted bound leaves zero sequence changes against the restricted
    // form's twelve).
    //
    // Member order is an input the normalizer must be indifferent to, so both
    // orders are asserted.
    let seq = "ACAAAAAAAACGTACGTACG";
    for input in ["TEMPLATE:g.[4A>G;9_10insA]", "TEMPLATE:g.[9_10insA;4A>G]"] {
        let actual = normalize_in(seq, input, ShuffleDirection::FivePrime);
        // Pinned, not merely sequence-preserving. The `apply` comparison below
        // is satisfied by an output identical to the input — which is exactly
        // what a clamp that failed to fire would produce — so on its own it
        // cannot tell "bounded correctly" from "nothing happened". The junction
        // must land at `3|4`, one short of the substituted base at 4, where it
        // coalesces with the sibling.
        assert_eq!(
            actual, "TEMPLATE:g.3_4insG",
            "{input} must bound its junction at the sibling's 5' edge"
        );
        assert_ne!(actual, input, "the clamp must actually move the junction");
        assert_eq!(
            apply(seq, &actual).expect("output applies"),
            apply(seq, input).expect("input applies"),
            "{input} -> {actual} changed the denoted sequence"
        );
    }
}

/// `(CAG)x10` at 1-based 11..=40, flanked by `T`. The tract the tandem sweep
/// uses; a junction can travel it without the payload being in phase everywhere,
/// which is what makes it the corpus for the 5' insertion-junction shape.
fn cag_tract() -> String {
    format!("{}{}{}", "T".repeat(10), "CAG".repeat(10), "T".repeat(20))
}

#[test]
fn an_insertion_junction_does_not_cross_an_upstream_junction_sibling() {
    // The other half of #1267's 5' shape. Here the upstream sibling adds
    // sequence at a junction instead of claiming bases, so `claims_bases` — the
    // bound that stops the substitution case above — does not see it at all:
    //
    //   g.[16_17insCAG;15dup]  ->  g.[11_13dup;15dup]
    //
    // The two members stay disjoint (11-13 against 15), so nothing flags it, but
    // the payload has moved from 3' of the sibling to 5' of it. Two junctions'
    // relative order is the only thing fixing the order of their payloads, which
    // is the same reasoning the 3' half already applies (#1290) — mirrored.
    let seq = cag_tract();
    for input in [
        "TEMPLATE:g.[16_17insCAG;15dup]",
        "TEMPLATE:g.[15dup;16_17insCAG]",
        "TEMPLATE:g.[37_38insCAG;36_37insTT]",
        "TEMPLATE:g.[36_37insTT;37_38insCAG]",
    ] {
        let actual = normalize_in(&seq, input, ShuffleDirection::FivePrime);
        assert_eq!(
            apply(&seq, &actual).expect("output applies"),
            apply(&seq, input).expect("input applies"),
            "{input} -> {actual} changed the denoted sequence"
        );
    }
}

#[test]
fn a_duplication_junction_does_not_cross_an_upstream_junction_sibling() {
    // The 5' junction bound applies to a **duplication** mover too, not only one
    // the input spelled as an insertion.
    //
    // Restricting it to input insertions looks like the conservative reading of
    // #1259 — a duplication's span and junction move together under a 5' shift,
    // and a blanket mirror of the 3' rule was measured there to turn 80
    // previously-correct outputs silently wrong. But that mirror bounded the
    // junction against a sibling the member was moving *away* from, which is a
    // different rule from this one. Restricting by input spelling instead leaves
    // the shape below broken: a `dup` whose junction sweeps past an upstream
    // sibling's junction, reordering the two payloads.
    //
    //   reference   GGGGGGGG + (A x 24) + GGGGGGGG
    //   g.[3dup;2_3insTT]   applies to  GGTTGGGGGGGAAA…
    //   unbounded           g.[1dup;2_3insTT]   GGGTTGGGGGG AAA…   <- payloads swapped
    //
    // Measured over 73,376 duplication-mover cases with an upstream sibling:
    // bounding every junction mover leaves **zero** sequence changes, while
    // restricting the bound to input insertions leaves **12** — all of this
    // shape. So the unrestricted rule is both simpler and strictly more correct,
    // and it does not decide the output from the input's encoding.
    let seq = format!("{}{}{}", "G".repeat(8), "A".repeat(24), "G".repeat(8));
    for input in [
        "TEMPLATE:g.[3dup;2_3insTT]",
        "TEMPLATE:g.[2_3insTT;3dup]",
        "TEMPLATE:g.[4dup;2_3insTT]",
        "TEMPLATE:g.[2_3insTT;4dup]",
    ] {
        let actual = normalize_in(&seq, input, ShuffleDirection::FivePrime);
        assert_eq!(
            apply(&seq, &actual).expect("output applies"),
            apply(&seq, input).expect("input applies"),
            "{input} -> {actual} changed the denoted sequence"
        );
    }
}

#[test]
fn a_deletion_does_not_shift_across_an_insertion_junction() {
    // The converse of the junction clamp above, and the case it did not cover.
    //
    // An insertion-like member is not a *barrier* for a span-claiming sibling —
    // it consumes no base, so `claims_reference_bases` says it claims nothing —
    // but its junction is still a position the sibling must not slide past.
    // Sliding a deletion across the point where a sibling adds sequence
    // reorders the two edits:
    //
    //   reference (A-run at 258-263)
    //   g.[258del;259_260insC]  ->  g.[259_260insC;263del]
    //   input applied   …T A C A A A A T…
    //   output applied  …T A A C A A A T…   ← the inserted base moved
    //
    // Clamped, the deletion stops at the junction, becomes adjacent to the
    // insertion, and the two coalesce.
    let mut bases = vec![b'T'; 300];
    for position in 258..=263 {
        bases[position - 1] = b'A';
    }
    let seq = String::from_utf8(bases).unwrap();

    // All three spellings of the same variant reach one form.
    for input in [
        "TEMPLATE:g.[258del;259_260insC]",
        "TEMPLATE:g.[259del;259_260insC]",
        "TEMPLATE:g.259delinsC",
    ] {
        assert_normalizes_preserving(&seq, input, "TEMPLATE:g.259A>C");
    }

    // A junction further into the tract clamps to that junction instead.
    assert_normalizes_preserving(&seq, "TEMPLATE:g.[258del;261_262insC]", "TEMPLATE:g.261A>C");

    // Negative control: with no sibling the deletion still shifts fully.
    assert_normalizes_preserving(&seq, "TEMPLATE:g.258del", "TEMPLATE:g.263del");
}

/// A non-homopolymer tandem tract, swept at every rotation and junction
/// (#1283 gap 2).
///
/// The sweep above cannot reach this class at all, and the reason is structural
/// rather than a matter of drawing more sequences:
///
/// * its corpus is 20 bp of `{A,T}` or `{A,C,G,T}` noise, which contains no
///   tandem tract long enough to canonicalise to a **non-homopolymer**
///   `Repeat`; and
/// * its insertion payloads are always a single copy, and a single-copy
///   insertion inside a tract canonicalises to a `dup` — so the multi-copy
///   repeat path is never entered.
///
/// #1280 lives exactly there: an out-of-phase payload in a `(CAG)` tract needs
/// a >=6-unit tract *and* a 2+-copy insertion. This sweeps a 10-unit `(CAG)`
/// tract, inserting 1..=3 copies of the unit at **every rotation** (`CAG`,
/// `AGC`, `GCA`) and at **every junction** across the tract, each against a
/// sibling on either side and in both shuffle directions.
///
/// The oracle is the existing SPDI apply: whatever ferro chooses to call the
/// result, it must denote the same bases. That is why this is a generator
/// change and not a new assertion — correctness is already decided.
#[test]
fn no_tandem_tract_allele_normalizes_to_a_different_sequence() {
    // 10 copies of CAG at 1-based 11..=40, flanked so members can sit outside
    // the tract on either side.
    const UNIT: &str = "CAG";
    const COPIES: usize = 10;
    let tract_start = 11usize;
    let tract: String = UNIT.repeat(COPIES);
    let tract_end = tract_start + tract.len() - 1; // 40
    let seq = format!("{}{tract}{}", "T".repeat(tract_start - 1), "T".repeat(20));

    let mut checked = 0usize;
    let mut skipped = 0usize;
    let mut overlapping: Vec<String> = Vec::new();
    let mut changed: Vec<String> = Vec::new();

    // One provider and one normalizer per direction for the whole sweep — this
    // one draws a single fixed reference, so there is nothing per-case about
    // either. `normalize_in` / `apply` built a fresh `TEMPLATE` provider per call.
    let template = provider(&seq);
    let normalizers = normalizers_for(&seq);

    // Every rotation of the unit: an insertion whose payload is out of phase
    // with the adjacent reference is the shape #1280 is about.
    let rotations = [UNIT.to_string(), "AGC".to_string(), "GCA".to_string()];

    // `tract_start - 1 ..= tract_end`, not `tract_start..tract_end`: the two
    // boundary junctions (`10_11` and `40_41`) are where an insertion can shift
    // out of the tract entirely, which is the shape closest to #1267. Excluding
    // them left the sweep claiming "every junction across the tract" while
    // covering only the strictly-interior ones.
    for junction in tract_start - 1..=tract_end {
        for rotation in &rotations {
            for copies in 1..=3usize {
                let payload = rotation.repeat(copies);
                let first = format!("{junction}_{}ins{payload}", junction + 1);
                // A sibling on each side of the tract, and two *inside* it.
                // The inside pair is the point: an insertion in a tract
                // canonicalises to a copy count over the **whole** tract, and
                // that span can swallow a sibling sitting within it — the
                // shape `a_repeat_grown_from_an_insertion_does_not_span_a_sibling`
                // pins by hand for a homopolymer, here swept over a
                // non-homopolymer one at every rotation.
                //
                // All four sibling shapes, not just `{sub, del}`. Restricting
                // the sibling to those two is exactly the gap-1 defect this PR
                // removes from the two-member sweep above, and it would be
                // reintroduced here otherwise: a `dup` or an `ins` sibling
                // never appearing means the tract class reads as exhaustive
                // over sibling *placement* while covering half the sibling
                // *shapes*. `dup` matters most — it is the one edit whose span
                // the 5' clamp treats as claimed bases (`blocks_sibling_shift`,
                // #1286), so it is the shape most likely to interact with a
                // repeat span.
                let inside_5p = tract_start + 4;
                let inside_3p = tract_end - 4;
                let base_at = |pos: usize| seq.as_bytes()[pos - 1] as char;
                // Each sibling names its own position, so the shapes are grouped
                // by comment rather than carried alongside a redundant index.
                let siblings = [
                    format!("{}T>A", tract_start - 2), // 5' of tract
                    format!("{}del", tract_start - 2),
                    format!("{}dup", tract_start - 2),
                    format!("{}_{}insA", tract_start - 2, tract_start - 1),
                    format!("{}T>A", tract_end + 2), // 3' of tract
                    format!("{}del", tract_end + 2),
                    format!("{}dup", tract_end + 2),
                    format!("{}_{}insA", tract_end + 2, tract_end + 3),
                    format!("{inside_5p}{}>T", base_at(inside_5p)), // inside
                    format!("{inside_5p}dup"),
                    format!("{inside_3p}del"),
                    // `T` is absent from a `(CAG)` tract, so this stays an
                    // `Insertion`. Inserting `base_at(inside_3p)` here — the
                    // reference base immediately 5' of the junction — denotes the
                    // same edit as `{inside_3p}dup`, so the array would have
                    // claimed four sibling shapes while carrying three.
                    format!("{inside_3p}_{}insTT", inside_3p + 1),
                ];
                for sibling in siblings {
                    // Author the sibling first or second: member order is
                    // an input the normalizer must be indifferent to.
                    //
                    // The member-order and direction loops are **nested the other
                    // way round** from how they were written, so that the input's
                    // parse and its applied sequence — neither of which depends on
                    // the shuffle direction — are computed once per description
                    // instead of once per direction. The case set is identical; only
                    // the order of enumeration changes, and nothing here is pinned
                    // to it (every bucket is asserted empty, and the samplers are
                    // capped at ten purely to keep a failure readable).
                    for input in [
                        format!("TEMPLATE:g.[{first};{sibling}]"),
                        format!("TEMPLATE:g.[{sibling};{first}]"),
                    ] {
                        let variant = parse_hgvs(&input).expect("generated input parses");
                        let want = apply_parsed_with(&template, &seq, &variant);
                        for (direction, normalizer) in &normalizers {
                            let output = normalizer
                                .normalize(&variant)
                                .expect("normalize")
                                .to_string();
                            checked += 1;
                            let Some(want) = want.as_deref() else {
                                skipped += 1;
                                continue;
                            };
                            // Every shape is asserted, including the 5'
                            // insertion-junction one that used to be excluded and
                            // counted (see the note above the sweep). A sibling
                            // upstream of the junction is the #1267 shape; it is
                            // bounded now, so a sequence change there is a
                            // regression like any other and reports the
                            // description rather than incrementing a count.
                            match apply_with(&template, &seq, &output) {
                                None if overlapping.len() < 10 => {
                                    overlapping.push(format!("{input} -> {output}"));
                                }
                                None => {}
                                Some(got) if got != want && changed.len() < 10 => {
                                    changed.push(format!(
                                        "{input} [{direction:?}] -> {output} \
                                         (want {want}, got {got})"
                                    ));
                                }
                                Some(_) => {}
                            }
                        }
                    }
                }
            }
        }
    }

    // 13,392 cases as written (31 junctions x 3 rotations x 3 copies x 12
    // siblings x 2 directions x 2 member orders). The floor sits ~20% under
    // that: low enough not to be brittle, high enough that dropping a rotation
    // or a copy count breaches it rather than passing quietly on a thinner
    // sweep.
    assert!(
        checked > 11_000,
        "tract sweep covered too little: {checked}"
    );
    // Same hollowness bound the two-member sweep carries: a case whose *input*
    // does not apply proves nothing, so cap that share rather than letting the
    // sweep pass by skipping.
    assert!(
        skipped * 10 < checked,
        "too many tract cases skipped as unconvertible: {skipped} of {checked}"
    );
    assert!(
        overlapping.is_empty(),
        "overlapping or unconvertible output in {checked} tract cases: {overlapping:#?}"
    );
    // No excluded shape: this covers the 5' insertion-junction cases too, which
    // is where the 116 pinned residuals used to live.
    assert!(
        changed.is_empty(),
        "sequence-changing normalizations in {checked} tract cases: {changed:#?}"
    );
}

/// A smoke slice of the same class on the **transcript axes** (#1283 gap 4).
///
/// Every sweep above runs on `g.` only, so the file's claim — that no two-member
/// cis allele normalizes to a different sequence — was never tested where the
/// axis itself changes the answer. It does change it: the CDS/UTR and
/// transcript-end insertion clamps, the exon-junction exception to the 3' rule,
/// the coding one-amino-acid exception to the separation rule and the `c.` codon
/// gate on repeat notation all fire here and nowhere above. #1284, which had to
/// land first because a transcript-axis member could abort the run outright under
/// `FERRO_ASSERT_REPARSE`, is closed.
///
/// Deliberately a **smoke slice**, not a second exhaustive sweep. It reuses the
/// same corpus, member shapes and oracle, so what it adds is the axis rather than
/// new shapes — and the genomic sweeps above already carry the shape coverage.
/// Sizing it that way is also what keeps it affordable: it multiplies the case
/// count by the number of axes, and the seed knob (#1295) is what makes that a
/// prefix locally and the full corpus in CI.
///
/// `cds_start = 1` so `c.N` is transcript position `N` and the coordinate
/// arithmetic matches the genomic sweeps exactly. The CDS is 18 of the 20 bases —
/// six whole codons, so the codon gate is not tripped by a ragged reading frame —
/// leaving `c.*1_*2` as a 3'UTR the members deliberately stay out of, since UTR
/// boundary behaviour has its own dedicated tests and would only add noise to a
/// smoke slice.
///
/// # No longer `#[ignore]`d — it took five fixes to get here
///
/// It was committed `#[ignore]`d rather than narrowed, deliberately: trimming the
/// position range until the boundary shape was unreachable would have been the
/// exact move #1283 exists to complain about — a generator that reads wider than
/// it is — and the gap it closes was worth the wait. What it found, in the order
/// each was unmasked by the one before it:
///
/// 1. **#1398** (#1417) — the original fire. `c.[11_12dup;16_17insC]` over
///    `TAAATAATAAATATATATTA` normalized to `c.18_*1insCAT`, an insertion across
///    the CDS/3'UTR boundary, and re-normalizing that did not return it
///    (`c.18delinsTCAT` under 3', `c.16_17insATC` under 5').
/// 2. **#1418** (#1425) — a sibling-crossing shift unbounded across a region
///    boundary.
/// 3. **#1426, 3' half** (#1427) — a junction insertion needing two passes to
///    finish its shift.
/// 4. **#1426, 5' half** (#1428) — the same defect arriving from the other side,
///    via the axis clamp rather than the #918 relaxation.
/// 5. **#1429** (this PR) — the sequence-first derivation shifting 3' whatever
///    direction was asked for.
///
/// Each was a genuine defect, and each was invisible until its predecessor was
/// fixed. That is the argument for committing a sweep `#[ignore]`d instead of
/// narrowing it: a narrowed generator would have found none of the five.
///
/// It runs in the `sweeps` job only — `cis_junction_crossing_shift` is inside
/// `SWEEP_FILTER`, which `test` and `test-oracle` negate — so it is exercised at
/// the full corpus with the oracles set, which is the configuration that found
/// all five. Measured at 280s for 829,440 cases over two axes (`c.`, `n.`).
///
/// #1471 added a third (non-coding `r.`), taking it to **1,244,160 cases**. That
/// count is exact and deterministic — three axes over the same corpus, shapes and
/// positions — so it is quoted where a wall-clock is not: the machine this was
/// developed on was running three concurrent copies of this sweep at a load
/// average of 81, which makes any timing taken there contention noise rather than
/// a measurement. Expect roughly 1.5x the 280s above on a quiet machine, and
/// treat the `sweeps` job's own duration as the figure of record.
///
/// Run it locally in the configuration `sweeps` uses — the oracles are not
/// optional here, since four of the five defects above were *reported* by
/// `FERRO_ASSERT_IDEMPOTENT` rather than by this sweep's own sequence comparison:
///
/// ```text
/// FERRO_SWEEP_SEEDS=full \
///   FERRO_ASSERT_IDEMPOTENT=1 FERRO_ASSERT_REPARSE=1 FERRO_ASSERT_IN_BOUNDS=1 \
///   cargo nextest run --features dev \
///   -E 'test(no_two_member_transcript_axis_allele)'
/// ```
/// The three transcript axes this sweep covers, each with the accession it is
/// drawn against and the name of the test that runs it.
///
/// One row per `#[test]` below, and [`every_transcript_axis_has_its_own_sweep`]
/// checks the two stay in step. That guard replaces the `per_axis` key assertion
/// this sweep used to carry: when all three axes ran inside one test, asserting
/// the key set was what made a deleted arm fail loudly instead of silently
/// shrinking a total that still cleared the floor. Split one test per axis, each
/// axis now has its own floor — stronger — but deleting a whole `#[test]` would
/// delete its own guard along with it, so the guard has to live outside all three.
const TRANSCRIPT_AXIS_SWEEPS: [(&str, &str, &str); 3] = [
    (
        "NM_TEST.1",
        "c",
        "no_two_member_cds_axis_allele_normalizes_to_a_different_sequence",
    ),
    (
        "NR_TEST.1",
        "n",
        "no_two_member_noncoding_axis_allele_normalizes_to_a_different_sequence",
    ),
    (
        "NR_TEST.1",
        "r",
        "no_two_member_rna_axis_allele_normalizes_to_a_different_sequence",
    ),
];

/// One axis of the transcript-axis sweep. See the module-level notes above the
/// three `#[test]` wrappers for what it covers and why each axis is swept.
///
/// **Split one test per axis to shorten CI's critical path.** All three ran in one
/// test, and that test alone set the `sweeps` job's duration — a job cannot finish
/// faster than its longest single test, so sharding the job could not help while
/// this was one 114s test. Three tests of a third the work run concurrently under
/// nextest instead. Nothing about the enumeration changes: same accessions, same
/// axes, same shapes, same positions, same directions, and each axis keeps the
/// per-seed floor it already had to clear individually.
fn sweep_one_transcript_axis(accession: &str, axis: &str) {
    use crate::common::cis_apply_oracle::{apply_parsed_with, apply_with};
    use crate::common::synthetic::SyntheticBuilder;
    use ferro_hgvs::reference::transcript::Strand;
    use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

    /// 1-based inclusive CDS end within the 20-base transcript: six codons.
    const CDS_END: u64 = 18;

    let seeds = sweep_seeds(48);

    // One rayon task per drawn sequence — see the genomic two-member sweep
    // above for why this is byte-identical rather than merely equivalent.
    // Parallelizing only ONE sweep in this file was measured as a wash: the
    // `Exhaustive sweeps` job runs several at once, so the parallel one simply
    // took cores from its serial siblings (job 1m42s -> 1m46s, with two
    // siblings newly crossing the SLOW line). They have to go together.
    let per_sequence: Vec<SweepTally> = sweep_sequences(seeds)
        .into_par_iter()
        .map(|seq| {
            let mut checked = 0usize;
            let mut skipped = 0usize;
            let mut overlapping: Vec<String> = Vec::new();
            let mut changed: Vec<String> = Vec::new();
            // All three transcript axes share the generator.
            //
            // `r.` was excluded until #1471 on two grounds, both of which were
            // measured and neither of which held:
            //
            // * *"its `U`/`T` alphabet makes the payload construction below a
            //   different problem"* — it does not. The parser accepts the
            //   DNA-alphabet payloads this generator already builds and renders
            //   them in the RNA alphabet on output, so not one line below changes:
            //   `NR_TEST.1:r.5_6insAA` parses as `r.5_6insaa` and `r.5T>A` as
            //   `r.5u>a`. The oracle needs no adjustment either — `apply_with`
            //   applies an `r.` description to the same bases as its `n.` twin
            //   (verified byte-identical), because `hgvs_to_spdi` folds the
            //   alphabet.
            // * *"`issue_1235_transcript_axes.rs` pins that axis by hand"* — it
            //   does, with 13 `r.` cases including a non-coding one, but they are
            //   substitutions and `delins`. Those reach `canonicalize_from_sequence`
            //   and never the *repair* passes, so none of them produces two
            //   junction-occupying members shifting onto one junction.
            //
            // That second gap is what #1453 fell through: on a non-coding `r.`
            // record every repair routed through `respell_at_gap` was a silent
            // no-op, so `r.[9dup;10dup;11dup]` normalized to `r.[9dup;11dup;11dup]`
            // — one interbase point claimed twice, denoting no sequence. It was
            // invisible to this suite and to two of the three seam oracles: the
            // output is well-formed, so `FERRO_ASSERT_REPARSE` accepts it, and every
            // coordinate is in range, so `FERRO_ASSERT_IN_BOUNDS` does too. Only the
            // apply oracle distinguishes it, as `CoincidentInsertions`.
            //
            // Verified against the defect rather than assumed: with #1465's fix
            // reverted, this arm fails at **four** seeds on the earliest shapes it
            // emits — `r.[2dup;4dup] -> r.[9dup;9dup]` and nine more, the collector's
            // cap. (The 849-output figure quoted in #1453 and #1471 comes from a
            // purpose-built 63 nt census in #1465, not from this sweep's 20-mers;
            // the two corpora are not comparable and only the *direction* carries
            // over.)
            //
            // **Non-coding `r.` only, and coding `r.` is a deliberate remaining
            // gap** — not a covered one. It would be convenient to say `c.` already
            // sweeps it, since `axis_frame` gives both `reading_frame: true`, but
            // that is only true of the frame: the two still take different code
            // paths (`build_rna_merged` vs `build_cds_merged`, and
            // `cds_relative_gap(Region::Rna)` vs `(Region::Cds)`), so a coding `r.`
            // defect could hide exactly the way the non-coding one did.
            //
            // The reason to stop here is cost, stated plainly: each axis is +50% on
            // what is already the slowest test in the suite. Non-coding `r.` buys a
            // *measured* gap — 849 corrupt outputs — and coding `r.` buys an
            // unmeasured one. Worth adding if a defect ever turns up there, or if
            // this sweep is made cheaper.
            // Built once per (sequence, axis) and cloned per case. `SyntheticBuilder`
            // pads, maps to a genomic contig and reverse-complements on demand, so
            // constructing one inside the innermost loop dominated the test:
            // measured at 39.8s for the 4-seed prefix before hoisting, against
            // 18.3s for the far larger genomic sweep beside it.
            // Keyed on the *accession*, not the axis, because that is what the
            // builders actually decide: `cds` names its transcript `NM_TEST.1` and
            // `noncoding` names its own `NR_TEST.1`, and every description below is
            // built as `{accession}:{axis}.…`. Testing `axis == "c"` instead would
            // route the coding `r.` row named in the gap note above — which is
            // `("NM_TEST.1", "r")` — to an `NR_TEST.1` provider, so every lookup
            // would fail to resolve and normalization would hand the input straight
            // back: a sweep that is green because it checked nothing, with the
            // axis label in every assertion message still reading `r.`.
            // `every_transcript_axis_has_its_own_sweep` does not reach this choice,
            // so an unrecognized accession is loud rather than defaulted.
            let axis_provider = match accession {
                "NM_TEST.1" => SyntheticBuilder::cds(&seq, 1, CDS_END, Strand::Plus).build(),
                "NR_TEST.1" => SyntheticBuilder::noncoding(&seq, Strand::Plus).build(),
                other => panic!(
                    "no reference shape is wired for `{other}` (sweeping the `{axis}.` axis); add \
                 one to this match rather than letting the row fall through to a transcript \
                 whose accession the descriptions never name"
                ),
            };
            for first_start in 2..=11usize {
                for first_len in 1..=2usize {
                    let first_end = first_start + first_len - 1;
                    let span = if first_len == 1 {
                        format!("{first_start}")
                    } else {
                        format!("{first_start}_{first_end}")
                    };
                    // Same payload rule as the genomic sweeps: excluding both
                    // flanking reference bases keeps the `ins` entry from
                    // denoting the `dup` beside it.
                    let first_base = seq.as_bytes()[first_start - 1] as char;
                    let first_next = seq.as_bytes()[first_start] as char;
                    let first_payload = ['A', 'C', 'G', 'T']
                        .into_iter()
                        .find(|b| *b != first_base && *b != first_next)
                        .expect("four bases, at most two excluded");
                    let inserted: String = std::iter::repeat_n(first_payload, first_len).collect();
                    let firsts = [
                        format!("{span}del"),
                        format!("{span}dup"),
                        format!("{first_start}_{}ins{inserted}", first_start + 1),
                    ];
                    for first in firsts {
                        for second_start in first_end + 2..=(CDS_END as usize - 1) {
                            let base = seq.as_bytes()[second_start - 1] as char;
                            let alt = if base == 'A' { 'G' } else { 'A' };
                            let next_base = seq.as_bytes()[second_start] as char;
                            let second_payload = ['A', 'C', 'G', 'T']
                                .into_iter()
                                .find(|b| *b != base && *b != next_base)
                                .expect("four bases, at most two excluded");
                            for second in [
                                format!("{second_start}del"),
                                format!("{second_start}{base}>{alt}"),
                                format!("{second_start}dup"),
                                format!("{second_start}_{}ins{second_payload}", second_start + 1),
                            ] {
                                // Everything that depends only on the *description*
                                // is hoisted out of the direction loop below. The
                                // input string, its parse and the bases it denotes
                                // are all properties of the description, not of the
                                // shuffle direction, and computing them per
                                // direction did each of them exactly twice.
                                let input = format!("{accession}:{axis}.[{first};{second}]");
                                let variant = parse_hgvs(&input).expect("generated input parses");
                                // Lazily, and at most once: a case both directions
                                // decline must stay as cheap as it was.
                                let mut input_applied: Option<Option<String>> = None;
                                for direction in
                                    [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime]
                                {
                                    // Deliberately built per case rather than once
                                    // per (sequence, axis). Hoisting it is strictly
                                    // less work and IS faster in the unoptimized
                                    // test profile, but it measured **14% slower**
                                    // on this test under the `soak` profile CI's
                                    // `sweeps` job actually uses — 11.54s -> 13.16s,
                                    // three reps a side at CV <= 0.5%. The
                                    // direction-invariant hoisting above (the
                                    // description, its parse, its applied sequence)
                                    // is kept: that is where the duplicated work was.
                                    let normalizer = Normalizer::with_config(
                                        axis_provider.clone(),
                                        NormalizeConfig::default()
                                            .with_direction(direction)
                                            .allow_crossing_boundaries(),
                                    );
                                    let Ok(normalized) = normalizer.normalize(&variant) else {
                                        // A refusal is a legitimate answer here —
                                        // the transcript axes decline shapes the
                                        // genomic axis accepts — and says nothing
                                        // about sequence preservation.
                                        skipped += 1;
                                        continue;
                                    };
                                    let output = format!("{normalized}");
                                    checked += 1;

                                    let Some(want) = input_applied
                                        .get_or_insert_with(|| {
                                            apply_parsed_with(&axis_provider, &seq, &variant)
                                        })
                                        .as_deref()
                                    else {
                                        // An *input* that does not apply is a case
                                        // this sweep cannot speak for. Counted
                                        // separately from an unconvertible output,
                                        // which is a failure.
                                        skipped += 1;
                                        checked -= 1;
                                        continue;
                                    };
                                    // Both collections are capped at ten, matching
                                    // the two sweeps above. The assertions below
                                    // fire on `is_empty()`, so every failure past
                                    // the tenth changes nothing about the verdict
                                    // and only makes the panic output harder to
                                    // read — and a broad regression on the full
                                    // corpus would otherwise accumulate a string
                                    // per case.
                                    match apply_with(&axis_provider, &seq, &output) {
                                        None if overlapping.len() < REPORTED => overlapping
                                            .push(format!("{input} [{direction:?}] -> {output}")),
                                        None => {}
                                        Some(got) if got != want && changed.len() < REPORTED => {
                                            changed.push(format!(
                                                "{input} [{direction:?}] -> {output} \
                                             (want {want}, got {got})"
                                            ))
                                        }
                                        Some(_) => {}
                                    }
                                }
                            }
                        }
                    }
                }
            }
            SweepTally {
                checked,
                skipped,
                overlapping,
                changed,
                ..SweepTally::default()
            }
        })
        .collect();

    let mut total = SweepTally::default();
    for tally in per_sequence {
        total.absorb(tally);
    }
    total.overlapping.truncate(REPORTED);
    total.changed.truncate(REPORTED);
    let SweepTally {
        checked,
        skipped,
        overlapping,
        changed,
        ..
    } = total;

    // Per seed, for the reason recorded on the floor in the first sweep above.
    // This is the SAME per-axis floor the combined test applied to each axis, so
    // splitting the test did not weaken it — it removed the aggregate bound, which
    // was the weaker of the two (three axes summed clear a per-axis floor even
    // with one arm deleted, which is exactly why the combined test also asserted
    // per-axis).
    const CASES_PER_SEED_FLOOR: usize = 1_000;
    let floor = CASES_PER_SEED_FLOOR * seeds as usize;
    assert!(
        checked > floor,
        "the `{axis}.` transcript-axis sweep covered too little: {checked} over \
         {seeds} seeds (floor {floor})"
    );
    // Now a per-axis share rather than an aggregate one, which is strictly
    // stricter: an axis with a high skip rate can no longer hide behind two with
    // low ones.
    //
    // **`checked` here EXCLUDES the skipped cases**, unlike the two genomic
    // sweeps above, where `checked` is incremented before the skip test and
    // `skipped` is therefore a subset of it. Both of this sweep's skip paths are
    // disjoint from `checked`: a normalize refusal returns before the increment,
    // and an input that does not apply decrements it back out. So the two floors
    // bound different ratios and must not be read as the same quantity — a share
    // of the enumerated total here, a share of the *accepted* cases there. The
    // enumerated total is stated in the message for that reason, matching how
    // the three-member sweep reports `input_declined`.
    assert!(
        skipped * 2 < checked,
        "too many `{axis}.` transcript-axis cases skipped: {skipped} of {} enumerated \
         ({checked} checked; `checked` excludes skipped cases on this axis)",
        skipped + checked
    );
    assert!(
        overlapping.is_empty(),
        "overlapping or unconvertible output in {checked} `{axis}.` transcript-axis \
         cases: {overlapping:#?}"
    );
    assert!(
        changed.is_empty(),
        "sequence-changing normalizations in {checked} `{axis}.` transcript-axis \
         cases: {changed:#?}"
    );
}

/// Sweep the [`TRANSCRIPT_AXIS_SWEEPS`] row whose third field is `test_name`.
///
/// **By name and not by index, because an index is silent when it is wrong.** The
/// wrappers below first selected their row as `TRANSCRIPT_AXIS_SWEEPS[0..=2]`, and
/// nothing in the module reached the subscript: write `[0]` in both the `cds` and
/// the `rna` wrapper and the `r.` axis is never swept, all three tests pass,
/// [`every_transcript_axis_has_its_own_sweep`] passes, and the only trace is two
/// `c.`-axis runs quoting different axis labels in their assertion messages. That
/// is the exact failure — "quietly reducing what this module covers" — the guard
/// exists to prevent, so the split must not reintroduce it one layer down.
///
/// A name cannot be wrong quietly: there is no row to run, so this panics.
fn sweep_axis_named(test_name: &str) {
    let (accession, axis, _) = TRANSCRIPT_AXIS_SWEEPS
        .iter()
        .find(|(_, _, name)| *name == test_name)
        .unwrap_or_else(|| {
            panic!("no TRANSCRIPT_AXIS_SWEEPS row names `{test_name}`, so it sweeps no axis")
        });
    sweep_one_transcript_axis(accession, axis);
}

#[test]
fn no_two_member_cds_axis_allele_normalizes_to_a_different_sequence() {
    sweep_axis_named("no_two_member_cds_axis_allele_normalizes_to_a_different_sequence");
}

#[test]
fn no_two_member_noncoding_axis_allele_normalizes_to_a_different_sequence() {
    sweep_axis_named("no_two_member_noncoding_axis_allele_normalizes_to_a_different_sequence");
}

#[test]
fn no_two_member_rna_axis_allele_normalizes_to_a_different_sequence() {
    sweep_axis_named("no_two_member_rna_axis_allele_normalizes_to_a_different_sequence");
}

/// Every axis in [`TRANSCRIPT_AXIS_SWEEPS`] has a `#[test]` that runs it.
///
/// The guard the split owes the suite. When all three axes ran inside one test,
/// `assert_eq!(per_axis.keys(), ["c","n","r"])` was what made deleting an arm fail
/// loudly; per-axis tests cannot carry that, because deleting the test deletes the
/// assertion with it. So this reads the module's own source — the same technique
/// `msto_regression_corpus::every_cataloged_test_exists_in_its_source_file` uses —
/// and fails when a row has no test, or a test has been renamed out from under its
/// row.
///
/// It is deliberately a source scan and not a runtime registry: nextest gives each
/// test its own process, so no test can observe whether its siblings ran.
///
/// # It checks the wiring too, not only the name
///
/// A test's *existence* is the weaker half. [`sweep_axis_named`] records why the
/// wrappers select their row by name rather than by subscript; the hole a name
/// leaves is two wrappers naming *one* row, which panics nowhere — both lookups
/// succeed, and an axis goes unswept while three tests pass. So each wrapper's body
/// is checked against the whole row set: between `fn <name>()` and the brace that
/// closes it, the literal `"<name>"` must appear and no *other* row's name may. The
/// second half is what actually rules out the duplicate — naming itself is
/// satisfiable by a body that also names something else.
///
/// Both halves are matched on the bare literal rather than on
/// `sweep_axis_named("<name>")`, because `rustfmt` breaks a long call across lines
/// and a contiguous match would then fail on formatting alone. No row name is a
/// substring of another, so `contains` cannot confuse two rows.
#[test]
fn every_transcript_axis_has_its_own_sweep() {
    let source = include_str!("cis_junction_crossing_shift.rs");
    for (accession, axis, test_name) in TRANSCRIPT_AXIS_SWEEPS {
        let signature = format!("fn {test_name}()");
        let start = source.find(&signature).unwrap_or_else(|| {
            panic!(
                "the `{axis}.` axis (drawn against {accession}) has no sweep: expected a \
                 test named `{test_name}`. Deleting an axis must fail here rather than \
                 quietly reducing what this module covers."
            )
        });
        // The wrapper is a single statement, so the first line that closes a block
        // at column zero ends it.
        let body = &source[start..];
        let body = &body[..body.find("\n}").map_or(body.len(), |end| end + 2)];
        assert!(
            body.contains(&format!("\"{test_name}\"")),
            "`{test_name}` does not look up its own row, so the `{axis}.` axis \
             (drawn against {accession}) may be unswept while this test passes. Its \
             body must name itself:\n{body}"
        );
        for (_, other_axis, other) in TRANSCRIPT_AXIS_SWEEPS {
            assert!(
                other == test_name || !body.contains(&format!("\"{other}\"")),
                "`{test_name}` names the `{other_axis}.` row as well as its own, so one \
                 of the two axes is swept twice and the other not at all:\n{body}"
            );
        }
    }
    // And the rows are distinct, so three rows cannot name one test.
    let names: std::collections::BTreeSet<&str> =
        TRANSCRIPT_AXIS_SWEEPS.iter().map(|(_, _, n)| *n).collect();
    assert_eq!(
        names.len(),
        TRANSCRIPT_AXIS_SWEEPS.len(),
        "two axes name the same test, so one axis is unswept"
    );
}

/// Three-member cis alleles, with the failure causes separated (#1268).
///
/// Both committed sweeps above generate exactly **two** members, so the ask
/// CodeRabbit made on #1257 — "extend the randomized overlap/convergence coverage
/// to alleles with three or more members" — was reported as addressed while the
/// artifact never covered it. Three members is not a marginal widening: the
/// measurement #1268 records found 0 overlapping outputs at two members and tens
/// of thousands at three.
///
/// # The causes are counted apart, which is the point
///
/// That measurement's other half was a 61,765-case "unapplicable" bucket that
/// conflated genuine overlap, `hgvs_to_spdi` conversion failure and parse
/// rejection — three quite different severities. #1268 asks that no three-member
/// number be quoted until they are separated, so this sweep reports
/// [`ApplyFailure`] per cause rather than a single count, and asserts on the two
/// that are defects in a normalizer's *output*:
///
/// * `Overlapping` / `CoincidentInsertions` — the output denotes no single
///   sequence. Both are failures here.
/// * `Unconvertible` / `Unparseable` / `OutOfBounds` / `StatedBasesMismatch` on
///   an **input** — a case this sweep cannot speak for, counted and bounded as a
///   share rather than asserted at zero.
///
/// `FERRO_ASSERT_REPARSE` (#1259) already covers parse rejection of normalizer
/// output suite-wide, so this does not re-derive it.
///
/// # Sizing
///
/// Three nested member loops over the same positions would be roughly the
/// two-member sweep cubed. The first and third members therefore **stride** their
/// position range by two rather than walking it: striding keeps both ends of the
/// range reachable, where truncating the range would silently drop the cases
/// nearest the sequence ends — which is where the shifts this file is about
/// actually terminate. Measured, at the 4-seed default: 131,328 cases and 60.8s
/// walking, 38,016 and 17.2s striding, against 18.3s for the two-member sweep
/// beside it. The seed knob (#1295) does the rest.
///
/// # What it measures today
///
/// All four buckets are **empty** over this corpus — 0 overlapping, 0
/// coincident-insertion, 0 sequence-changing outputs, and 0 inputs declined for
/// any cause. That is the re-measurement #1268 asks for before any three-member
/// number is quoted again, and it is a long way from what that issue recorded
/// from the throwaway probe: 10,881 sequence-changing and a 61,765-case
/// undifferentiated "unapplicable" bucket. The confluence work between then and
/// now is what closed it; this sweep is what keeps it closed.
///
/// Because every bucket is zero, there is nothing to pin as a known residual —
/// the `assert!(…is_empty())` form reports the offending descriptions directly,
/// which is strictly better than a count, and #1268's request for named
/// exclusions cross-referenced to issues has no subject.
#[test]
fn no_three_member_allele_normalizes_to_a_different_sequence() {
    use crate::common::cis_apply_oracle::{
        apply_parsed_reason, apply_reason, provider as genomic_provider, ApplyFailure,
    };
    use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

    /// The tally key for one decline cause. Named so the per-direction tallying
    /// below stays a one-liner after the apply itself was hoisted out of that loop.
    fn cause_name(cause: ApplyFailure) -> &'static str {
        match cause {
            ApplyFailure::Unparseable => "unparseable",
            ApplyFailure::Unconvertible => "unconvertible",
            ApplyFailure::OutOfBounds => "out-of-bounds",
            ApplyFailure::Overlapping => "overlapping",
            ApplyFailure::CoincidentInsertions => "coincident-insertions",
            ApplyFailure::StatedBasesMismatch => "stated-bases-mismatch",
        }
    }

    // 16 rather than the 48 its neighbours use. Each seed here carries roughly
    // three times their cases, so 16 seeds is comparable *total* coverage — and
    // at 48 this one test was 411s in CI's sweeps job, against 183s for all three
    // existing sweeps together. It is a deliberate starting point for a sweep
    // whose buckets are all currently zero, not a ceiling: raise it if a
    // three-member defect ever turns up in a sequence beyond the sixteenth.
    let seeds = sweep_seeds(16);

    // One rayon task per drawn sequence — see the genomic two-member sweep
    // above for why this is byte-identical rather than merely equivalent.
    // Parallelizing only ONE sweep in this file was measured as a wash: the
    // `Exhaustive sweeps` job runs several at once, so the parallel one simply
    // took cores from its serial siblings (job 1m42s -> 1m46s, with two
    // siblings newly crossing the SLOW line). They have to go together.
    let per_sequence: Vec<SweepTally> = sweep_sequences(seeds)
        .into_par_iter()
        .map(|seq| {
            let mut checked = 0usize;
            let mut input_declined = 0usize;
            let mut overlapping: Vec<String> = Vec::new();
            let mut changed: Vec<String> = Vec::new();
            // Per-cause tallies for the *input* declines, so the share below is
            // decomposable rather than another single opaque number.
            let mut by_cause: std::collections::BTreeMap<&'static str, usize> =
                std::collections::BTreeMap::new();
            // One provider and one normalizer per direction per sequence. The
            // innermost loop used to build `genomic_provider(&seq)` three times per
            // case — for the input apply, the normalizer and the output apply.
            let template = genomic_provider(&seq);
            let normalizers =
                [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime].map(|direction| {
                    (
                        direction,
                        Normalizer::with_config(
                            template.clone(),
                            NormalizeConfig::default()
                                .with_direction(direction)
                                .allow_crossing_boundaries(),
                        ),
                    )
                });
            // Even positions only, and the bound says 8 rather than 9 because
            // `step_by(2)` from 2 can never reach an odd endpoint — writing `..=9`
            // named a position this loop does not visit. Sampling every other
            // position is deliberate (the third member's loop below multiplies this
            // one), but the range should not claim coverage it does not provide.
            for first_start in (2..=8usize).step_by(2) {
                let first_base = seq.as_bytes()[first_start - 1] as char;
                let first_next = seq.as_bytes()[first_start] as char;
                let first_payload = ['A', 'C', 'G', 'T']
                    .into_iter()
                    .find(|b| *b != first_base && *b != first_next)
                    .expect("four bases, at most two excluded");
                let firsts = [
                    format!("{first_start}del"),
                    format!("{first_start}dup"),
                    format!("{first_start}_{}ins{first_payload}", first_start + 1),
                ];
                for first in &firsts {
                    // Second and third members walk the remaining positions, each
                    // separated from its predecessor by at least one unchanged base
                    // so every case starts as three genuinely separate members.
                    for second_start in first_start + 2..=14usize {
                        let second_base = seq.as_bytes()[second_start - 1] as char;
                        let second_alt = if second_base == 'A' { 'G' } else { 'A' };
                        for second in [
                            format!("{second_start}del"),
                            format!("{second_start}{second_base}>{second_alt}"),
                            format!("{second_start}dup"),
                        ] {
                            // `..=19` is right here, unlike the `first_start` loop
                            // above: `second_start` walks every integer, so when it
                            // is odd this range is odd-valued and 19 IS visited.
                            for third_start in (second_start + 2..=19usize).step_by(2) {
                                let third_base = seq.as_bytes()[third_start - 1] as char;
                                let third_alt = if third_base == 'A' { 'G' } else { 'A' };
                                for third in [
                                    format!("{third_start}del"),
                                    format!("{third_start}{third_base}>{third_alt}"),
                                ] {
                                    // The description, its parse and the bases it
                                    // denotes do not depend on the shuffle direction,
                                    // so each is computed once per case rather than
                                    // once per direction. A parse failure is folded
                                    // into the `ApplyFailure` the string-taking oracle
                                    // would have reported for it, so the tallying below
                                    // is unchanged — including that it counts a
                                    // declined input once *per direction*, which is the
                                    // quantity the share assertion is written against.
                                    let input = format!("TEMPLATE:g.[{first};{second};{third}]");
                                    let parsed =
                                        parse_hgvs(&input).map_err(|_| ApplyFailure::Unparseable);
                                    let input_applied =
                                        parsed.as_ref().map_err(|c| *c).and_then(|variant| {
                                            apply_parsed_reason(&template, &seq, variant)
                                        });
                                    for (direction, normalizer) in &normalizers {
                                        let want = match &input_applied {
                                            Ok(want) => want.as_str(),
                                            Err(cause) => {
                                                // An input this sweep cannot speak
                                                // for. Tallied by cause so the share
                                                // asserted below can be read.
                                                *by_cause.entry(cause_name(*cause)).or_default() +=
                                                    1;
                                                input_declined += 1;
                                                continue;
                                            }
                                        };
                                        let variant = parsed.as_ref().expect(
                                            "a parse failure is reported as a decline above",
                                        );
                                        let Ok(normalized) = normalizer.normalize(variant) else {
                                            input_declined += 1;
                                            *by_cause.entry("normalize-declined").or_default() += 1;
                                            continue;
                                        };
                                        let output = format!("{normalized}");
                                        checked += 1;

                                        match apply_reason(&template, &seq, &output) {
                                            Ok(got) if got == want => {}
                                            Ok(got) => changed.push(format!(
                                                "{input} [{direction:?}] -> {output} \
                                             (want {want}, got {got})"
                                            )),
                                            Err(cause) => overlapping.push(format!(
                                                "{input} [{direction:?}] -> {output} ({cause:?})"
                                            )),
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            SweepTally {
                checked,
                input_declined,
                by_cause,
                overlapping,
                changed,
                ..SweepTally::default()
            }
        })
        .collect();

    let mut total = SweepTally::default();
    for tally in per_sequence {
        total.absorb(tally);
    }
    total.overlapping.truncate(REPORTED);
    total.changed.truncate(REPORTED);
    let SweepTally {
        checked,
        input_declined,
        by_cause,
        overlapping,
        changed,
        ..
    } = total;

    const CASES_PER_SEED_FLOOR: usize = 2_000;
    let floor = CASES_PER_SEED_FLOOR * seeds as usize;
    assert!(
        checked > floor,
        "three-member sweep covered too little: {checked} over {seeds} seeds (floor {floor})"
    );
    // The share this sweep cannot speak for, bounded so it cannot go hollow
    // while still clearing the floor. Decomposed in the message, per #1268 —
    // a single "unapplicable" number is what that issue exists to reject.
    assert!(
        input_declined < checked,
        "more three-member inputs declined than checked: {input_declined} of \
         {} enumerated, by cause: {by_cause:#?}",
        input_declined + checked
    );
    assert!(
        overlapping.is_empty(),
        "output denotes no single sequence in {checked} three-member cases \
         (inputs declined: {input_declined}, by cause {by_cause:#?}): {overlapping:#?}"
    );
    assert!(
        changed.is_empty(),
        "sequence-changing normalizations in {checked} three-member cases: {changed:#?}"
    );
}
