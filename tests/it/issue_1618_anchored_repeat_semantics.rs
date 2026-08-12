//! An anchored repeat's count is counted against the tract of the unit it
//! spells, and re-phasing that unit must carry the copies it absorbs (#1618).
//!
//! # What the spec fixes, and what it leaves to us
//!
//! `RNA/repeated.md:31` states the semantics: a repeat description names a unit
//! and a **total copy count**, "with the first unit located from position …, is
//! present in 14 copies". `:33` then pins it against an independent yardstick —
//! with a reference tract of **15** units, `ug[14]` is *preferred over*
//! `-97_-96del`, i.e. it denotes a net **−2** bases, and `ug[17]` a net **+4**.
//!
//! That equivalence is what makes this testable without picking a side in
//! #1618: it judges the reading from outside ferro. `hgvs_to_spdi` passes it
//! (measured: `g.258GT[14]` over a 15-copy tract resolves to
//! `257:GT×15:GT×14`, net −2), so its reading of an anchored repeat — the
//! maximal tract **of the spelled unit** — is the spec's, and it is sound to
//! use as the baseline for what an input denotes.
//!
//! # The defect
//!
//! The reference `GTGT` at `g.259..262` is 2 copies of `GT` but only 1 copy of
//! `TG`. So `g.262TG[6]` denotes 6 copies of `TG` replacing the 1-copy `TG`
//! tract — 14 bases — while `g.259_262GT[6]` denotes 6 copies of `GT` replacing
//! the 2-copy `GT` tract — 12. **Both are correct for their own unit.**
//!
//! `normalize_repeat` re-phases `TG` to `GT` to land on the longer tract and
//! keeps the count, so the extra reference copy it absorbs is silently dropped.
//! Re-phasing is right (it reaches a canonical, maximal tract); keeping the
//! count is not. The absorbed copy has to be carried: `GT[7]`, not `GT[6]`.
//!
//! Three of the four seam oracles pass on the wrong output — it is well-formed,
//! in bounds, re-parses and is a fixed point — which is why this survived. Only
//! `FERRO_ASSERT_SEQUENCE` sees it, and these tests assert the same property
//! unconditionally so it is guarded whether or not that flag is set.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::spdi::{compare_denoted_sequences, DenotedSequenceComparison};
use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer, ShuffleDirection};

/// `AA` + `GTGT` (g.259..262) + `TA`. Core base 1 is g.257. The flanking `A`
/// (g.258) and `T` (g.263) break the `GT` phase on both sides, so the tract
/// cannot slide in either direction and both shuffle directions must agree.
const CORE: &str = "AAGTGTTA";

fn norm_in(core: &str, input: &str, direction: ShuffleDirection) -> String {
    Normalizer::with_config(
        SyntheticBuilder::genomic(core).build(),
        NormalizeConfig::default().with_direction(direction),
    )
    .normalize(&parse_hgvs(input).expect("parse"))
    .expect("normalize")
    .to_string()
}

fn norm(core: &str, input: &str) -> String {
    norm_in(core, input, ShuffleDirection::ThreePrime)
}

fn parsed(input: &str) -> HgvsVariant {
    parse_hgvs(input).expect("parse")
}

/// Fail unless normalizing `input` leaves the denoted bases untouched.
///
/// Asserts on the bases rather than on the output string, so it keeps its
/// meaning if the canonical spelling is ever revised — the property is that
/// normalization does not change what the description means.
fn assert_normalization_preserves_bases(core: &str, input: &str, direction: ShuffleDirection) {
    let output = norm_in(core, input, direction);
    let provider = SyntheticBuilder::genomic(core).build();
    match compare_denoted_sequences(&parsed(input), &parsed(&output), &provider) {
        DenotedSequenceComparison::Agree => {}
        DenotedSequenceComparison::Differ {
            reference,
            from_input,
            from_output,
            ..
        } => panic!(
            "normalization changed the denoted sequence\n  input     {input}\n  output    \
             {output}\n  reference {reference}\n  from input  {from_input} ({} bases)\n  from \
             output {from_output} ({} bases)",
            from_input.len(),
            from_output.len(),
        ),
        other => panic!("expected a comparison, got {other:?}\n  input {input}\n  output {output}"),
    }
}

#[test]
fn re_phasing_a_repeat_unit_preserves_the_denoted_sequence_three_prime() {
    assert_normalization_preserves_bases(
        CORE,
        "NC_TEST.1:g.262TG[6]",
        ShuffleDirection::ThreePrime,
    );
}

#[test]
fn re_phasing_a_repeat_unit_preserves_the_denoted_sequence_five_prime() {
    assert_normalization_preserves_bases(CORE, "NC_TEST.1:g.262TG[6]", ShuffleDirection::FivePrime);
}

/// The exact canonical spelling, so a fix that preserves the bases by declining
/// to re-phase at all is distinguishable from one that re-phases and carries the
/// absorbed copy. The maximal tract is the 2-copy `GT`, and the input's 14 bases
/// are 7 copies of it.
#[test]
fn a_re_phased_repeat_carries_the_copies_it_absorbs() {
    assert_eq!(
        norm(CORE, "NC_TEST.1:g.262TG[6]"),
        "NC_TEST.1:g.259_262GT[7]",
        "re-phasing TG -> GT absorbs one reference copy, so the count must grow 6 -> 7"
    );
}

/// The two spellings denote *different* variants (14 bases against 12), so they
/// must NOT converge. Pinned because the pre-fix behaviour converged them, and a
/// bare "they agree" assertion reads as confluence while actually asserting that
/// one of them lost bases.
#[test]
fn spellings_that_denote_different_sequences_do_not_converge() {
    let anchored = norm(CORE, "NC_TEST.1:g.262TG[6]");
    let ranged = norm(CORE, "NC_TEST.1:g.259_262GT[6]");
    assert_ne!(
        anchored, ranged,
        "TG[6] denotes 14 bases and GT[6] denotes 12; converging them destroys one"
    );
    assert_eq!(
        ranged, "NC_TEST.1:g.259_262GT[6]",
        "the ranged spelling is already canonical"
    );
}

/// Re-phasing must still be idempotent — the property the maximization was
/// written for, which the fix must not cost.
#[test]
fn the_re_phased_output_is_a_fixed_point() {
    let once = norm(CORE, "NC_TEST.1:g.262TG[6]");
    let twice = norm(CORE, &once);
    assert_eq!(
        once, twice,
        "NOT IDEMPOTENT\n  once ={once}\n  twice={twice}"
    );
}

// ---------------------------------------------------------------------------
// The spec's own equivalence, and the preference it states alongside it.
// ---------------------------------------------------------------------------

/// `A` + `GT`×15 + `A`: core base 1 is g.257, so the 15-copy tract is
/// g.258..287 — the shape `RNA/repeated.md:33` states its equivalences on.
fn fifteen_copy_core() -> String {
    format!("A{}A", "GT".repeat(15))
}

/// The expansion direction keeps the repeat form.
///
/// Pinned here because it is half of an **asymmetry** filed as #1675: a
/// contraction by one unit is rewritten to `del` (`g.258GT[14]` ->
/// `g.286_287del`) although `RNA/repeated.md:33` states a preference for the
/// repeat description in both directions. This test is what shows that rewrite
/// is not a uniform policy of preferring `del`/`dup`.
///
/// The contraction assertion is deliberately **not** here. Its only citation is
/// RNA-axis while the reproducer is genomic, and applying it to `g.`/`c.` is a
/// cross-axis generalization that needs a ruling first — the same class of move
/// #1664 was about. Pinning it as correct before that ruling would be smuggling
/// in the answer.
#[test]
fn an_expanding_repeat_keeps_the_repeat_form() {
    let core = fifteen_copy_core();
    assert_eq!(
        norm(&core, "NC_TEST.1:g.258GT[17]"),
        "NC_TEST.1:g.258_287GT[17]"
    );
}

/// Both directions must preserve the bases they denote, which is the invariant
/// the rule 2 preference above must not be bought at the expense of.
#[test]
fn the_spec_equivalence_shapes_preserve_their_denoted_bases() {
    let core = fifteen_copy_core();
    for input in ["NC_TEST.1:g.258GT[14]", "NC_TEST.1:g.258GT[17]"] {
        assert_normalization_preserves_bases(&core, input, ShuffleDirection::ThreePrime);
    }
}

// ---------------------------------------------------------------------------
// `DNA/repeated.md:91-99` — the CFTR intron 9 discussion, which is the spec's
// own worked case on this exact `GT`/`TG` pair.
// ---------------------------------------------------------------------------

/// The CFTR shape: `A` + `TG`×11 + `T`×7 + `A`. Core base 1 is g.257, so the
/// `TG` tract runs g.258..279 and the `T` run g.280..286.
///
/// The spec's question spells the tract as `TG`; its answer is:
///
/// > First, note that by applying the **3'rule** it is a **variable GT and not
/// > a TG stretch**. When the coding DNA reference sequence has 11 TG copies
/// > followed by 7 T copies, the reference allele is described as
/// > `c.1210-33_1210-6GT[11]T[6]`.
///
/// Two things are pinned by that answer, and both are load-bearing here.
/// **Re-phasing preserves the copy count** — 11 `TG` becomes 11 `GT`, not 10 or
/// 12. And **when the re-phased tract covers different bases, the counts move
/// with it**: sliding one base 3' swallows the first `T`, so the spec decrements
/// that run 7 -> 6 and the total is conserved (1 + 22 + 6 = 29 = 22 + 7). It
/// does not leave `T[7]` standing and silently gain a base.
///
/// That conservation is the principle this file's other tests extend from a
/// slide to a widening.
const CFTR_CORE_TG: usize = 11;
const CFTR_CORE_T: usize = 7;

fn cftr_core() -> String {
    format!("A{}{}A", "TG".repeat(CFTR_CORE_TG), "T".repeat(CFTR_CORE_T))
}

/// `repeated.md:97`: "by applying the **3'rule** it is a **variable GT and not a
/// TG stretch**". A `TG`-spelled tract must come back spelled `GT`.
#[test]
fn a_tg_stretch_is_re_spelled_as_a_gt_stretch() {
    let out = norm(&cftr_core(), "NC_TEST.1:g.258TG[13]");
    assert!(
        out.contains("GT["),
        "the 3' rule makes this a variable GT stretch, not TG; got {out}"
    );
    assert!(
        !out.contains("TG["),
        "the TG phase must not survive the 3' rule; got {out}"
    );
}

/// The conservation the spec demonstrates by moving `T[7]` to `T[6]`: re-phasing
/// must not change how many bases the description denotes.
#[test]
fn re_phasing_the_cftr_tract_conserves_the_denoted_bases() {
    let core = cftr_core();
    for input in [
        "NC_TEST.1:g.258TG[13]", // the questioner's expansion, TG-spelled
        "NC_TEST.1:g.258TG[9]",  // a contraction, the other direction
    ] {
        assert_normalization_preserves_bases(&core, input, ShuffleDirection::ThreePrime);
    }
}

/// A pure re-phase preserves the count — 11 `TG` is 11 `GT`, per the spec's own
/// `GT[11]`. Pinned as the boundary of the fix: carrying absorbed copies must
/// not turn into inflating the count on a tract that merely slid.
#[test]
fn a_pure_re_phase_preserves_the_copy_count() {
    let core = cftr_core();
    let out = norm(&core, "NC_TEST.1:g.258TG[13]");
    assert!(
        out.contains("GT[13]"),
        "sliding phase must not change the count the caller asked for; got {out}"
    );
}

// ---------------------------------------------------------------------------
// The other direction: a re-phase that lands on a SMALLER tract.
//
// The maximization seeds `best` with the literal tract only when that tract
// spans the anchor, so in the #1618 shape — where it does not — the winning
// rotation is not constrained to hold more copies than the spelled unit's own
// tract. The count therefore has to be re-based in both directions, and an
// arithmetic that measures the window's two edges independently
// (`(ls - s) + (e - le)`, each saturating) can only ever report a gain.
// ---------------------------------------------------------------------------

/// `AA` + `TGTGTG` (g.259..264) + `GT` (g.265..266) + `AA`.
///
/// Anchored at g.265, the spelled `TG` tiles a **3**-copy tract that ends at the
/// anchor, while the `GT` rotation tiles a **1**-copy tract that spans it. So
/// the maximization re-phases onto a tract two copies SHORTER than the one the
/// count was written against.
const NARROWING_CORE: &str = "AATGTGTGGTAA";

/// The tract the re-phase lands on holds two copies fewer than the tract the
/// spelled unit names, so the count must fall by two — `TG[6]` against a 3-copy
/// `TG` tract is a net +6 bases, and only `GT[4]` against a 1-copy `GT` tract
/// carries that net through.
///
/// The per-edge form shipped +1 here (`GT[7]`, a net +12): the tract's 5' edge
/// moved 3' by three units, which `ls.saturating_sub(s)` reads as zero, while
/// its 3' edge moved 3' by one, which is reported as a gain.
///
/// **And re-basing is not enough here.** The two tracts are *disjoint* — `TG`
/// runs g.259..264 and `GT` g.265..266 — so the description's window moves to a
/// different stretch of reference and the re-based count conserves the net
/// length and nothing else: `g.265_266GT[4]` denotes `TGTGTGGTGTGTGT` where the
/// input denotes `TGTGTGTGTGTGGT`, measured by `compare_denoted_sequences`. The
/// re-phase is therefore **declined**, and the count stands against the tract it
/// was written for.
#[test]
fn a_re_phase_onto_a_narrower_tract_keeps_the_spelled_units_tract() {
    assert_eq!(
        norm(NARROWING_CORE, "NC_TEST.1:g.265TG[6]"),
        "NC_TEST.1:g.259_264TG[6]",
        "the GT tract is disjoint from the TG tract the count was written \
         against, so the re-phase is declined rather than re-based"
    );
    assert_normalization_preserves_bases(
        NARROWING_CORE,
        "NC_TEST.1:g.265TG[6]",
        ShuffleDirection::ThreePrime,
    );
}

/// The saturating region, which no test reached and which cost two bases.
///
/// `g.265TG[1]` re-bases as `1 + 1 - 3`. The shipped `saturating_add`/
/// `saturating_sub` pair clamped that to `0`, and `specified_count >= 1` then
/// routed it to a two-base deletion:
///
/// ```text
/// NC_TEST.1:g.265TG[1]  ->  NC_TEST.1:g.265_266del
///   reference    TGTGTGGT
///   from input   TGGT      (4 bases — the input denotes a net -4)
///   from output  TGTGTG    (6 bases — the output denotes a net -2)
/// ```
///
/// `[0]` is never emitted, so nothing downstream could catch it. The saturation
/// was justified in comment by the *disjoint* case, where `literal == winner`
/// and nothing can saturate at all; it fires on `literal > specified + winner`,
/// which is a different shape entirely.
#[test]
fn a_re_phase_the_count_cannot_survive_keeps_the_spelled_units_tract() {
    assert_eq!(
        norm(NARROWING_CORE, "NC_TEST.1:g.265TG[1]"),
        "NC_TEST.1:g.259_264TG[1]",
        "re-basing 1 + 1 - 3 underflows, so the re-phase must be declined \
         rather than clamped to a count that denotes different bases"
    );
    assert_normalization_preserves_bases(
        NARROWING_CORE,
        "NC_TEST.1:g.265TG[1]",
        ShuffleDirection::ThreePrime,
    );
}

/// The whole narrowing family, on the bases rather than on any one spelling.
///
/// `[2]` is the boundary the saturation reached but did not clamp — `2 + 1 - 3`
/// is exactly `0` — and it shipped `g.265_266del`, which denotes `TGTGTG` where
/// the input denotes `TGTGGT`: the same *length*, different bases. A test that
/// only checked the net change would have called it correct, which is why this
/// asserts through `compare_denoted_sequences`.
#[test]
fn every_count_over_the_narrowing_tract_preserves_its_denoted_bases() {
    for count in 1..=8 {
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            assert_normalization_preserves_bases(
                NARROWING_CORE,
                &format!("NC_TEST.1:g.265TG[{count}]"),
                direction,
            );
        }
    }
}

/// `GG` + `AC` (g.259..260) + `CA` (g.261..262) + `GG`.
///
/// The narrowing case's boundary: anchored at g.261 the spelled `AC` tract
/// (g.259..260) and the winning `CA` rotation (g.261..262) hold **one copy
/// each**, so the tract has moved without changing size. The per-edge form read
/// the 3' edge's one-unit advance as an absorbed copy and shipped `CA[6]`, and
/// re-basing corrected the count to `CA[5]` — while leaving the window on a
/// different stretch of reference.
///
/// Equal counts are exactly the case that shows re-basing is the wrong lever
/// here: there is nothing for it to correct, and the output still denotes
/// different bases. The two tracts are disjoint, so the re-phase is declined.
const DISJOINT_CORE: &str = "GGACCAGG";

#[test]
fn a_disjoint_re_phase_keeps_the_spelled_units_tract() {
    assert_eq!(
        norm(DISJOINT_CORE, "NC_TEST.1:g.261AC[5]"),
        "NC_TEST.1:g.259_260AC[5]",
        "the AC and CA tracts are disjoint, so the count stays on the AC tract \
         it was written against rather than being re-phased onto CA"
    );
    assert_normalization_preserves_bases(
        DISJOINT_CORE,
        "NC_TEST.1:g.261AC[5]",
        ShuffleDirection::ThreePrime,
    );
}

// ---------------------------------------------------------------------------
// A2: the `None` literal arm, which takes no adjustment at all.
// ---------------------------------------------------------------------------

/// `count_tandem_repeats` returns `None` when the spelled unit tiles nowhere at
/// or immediately 5' of the anchor — the #852 shape the rotation search exists
/// for. There are then no copies of the spelled unit to re-base against, so the
/// winning rotation's count is emitted as spelled and the two anchors into one
/// and the same `GT` tract disagree by one:
///
/// ```text
/// g.259TG[6]  ->  g.259_262GT[6]    literal TG tract: None,     no adjustment
/// g.262TG[6]  ->  g.259_262GT[7]    literal TG tract: 1 copy,   +1 carried
/// ```
///
/// **This pins a limit, not a fix**, and the reason it is left as one is
/// measured rather than assumed: `compare_denoted_sequences` reports
/// `InputDenotesNoSequence` for `g.259TG[6]`, i.e. `hgvs_to_spdi` — the
/// independent reading this file's baseline rests on — finds no bases for the
/// input to denote. There is no baseline to conserve, so any count chosen here
/// would be a ruling on what a repeat spelled against a tract that does not
/// exist means, not a rescue of one. The `g.262TG[6]` sibling, which *does*
/// have a baseline, is asserted alongside so the asymmetry is visible rather
/// than inferred.
#[test]
fn the_none_literal_arm_takes_no_adjustment() {
    assert_eq!(
        norm(CORE, "NC_TEST.1:g.259TG[6]"),
        "NC_TEST.1:g.259_262GT[6]",
        "no copies of the spelled TG at or 5' of g.259, so nothing is re-based"
    );
    assert_eq!(
        norm(CORE, "NC_TEST.1:g.262TG[6]"),
        "NC_TEST.1:g.259_262GT[7]",
        "one copy of the spelled TG ends at g.262, so the re-phase carries it"
    );

    let provider = SyntheticBuilder::genomic(CORE).build();
    let output = norm(CORE, "NC_TEST.1:g.259TG[6]");
    assert!(
        matches!(
            compare_denoted_sequences(&parsed("NC_TEST.1:g.259TG[6]"), &parsed(&output), &provider),
            DenotedSequenceComparison::NotComparable(_)
        ),
        "the None arm is left alone because the input denotes no sequence to \
         conserve; if this input ever gains a baseline, the arm needs a ruling"
    );
}

// ---------------------------------------------------------------------------
// The SEEDED overlapping re-phase — the shape where the re-basing moves the
// emitted **edit kind**, not merely the count.
//
// Every shape above is either non-seeded (the spelled unit's tract does not
// span the anchor, so `best` starts empty) or disjoint (the re-phase is
// declined outright). This one is neither: the spelled unit's tract *spans* the
// anchor and therefore seeds `best`, and a strictly longer rotation of it
// displaces the seed while still overlapping it. The re-phase is carried, the
// count is re-based by exactly the copies the wider window swallows — and
// because a repeat whose re-based count lands one unit either side of the
// reference count is rendered as `del`/`dup` rather than as a repeat, shifting
// the count by one also shifts *which* count renders as which edit.
//
// A consumer reading "the count now carries the copies the wider window
// swallows" would not anticipate a stored `repeat` coming back as a `del` or a
// `dup`, so it is pinned here and called out in the representation-change
// declaration.
// ---------------------------------------------------------------------------

/// `A` + `TGTGTGTG` (g.258..265) + `A`. Core base 1 is g.257.
///
/// Anchored at **g.259** the spelled `GT` tract is g.259..264 — **3 copies**,
/// and it *spans* the anchor, so it seeds `best`. The `TG` rotation taken from
/// g.258 tiles g.258..265 — **4 copies** — and being strictly longer it
/// displaces the seed. The two tracts overlap, so the re-phase is carried and
/// the count is re-based `N + 4 - 3 = N + 1`.
const SEEDED_OVERLAP_CORE: &str = "ATGTGTGTGA";

/// The count carries the one copy the wider `TG` window swallows.
///
/// Measured against `origin/main`, which emitted the count unchanged and so
/// denoted two bases fewer than the input at every count:
///
/// ```text
/// g.259GT[1]   TG[1] -> TG[2]      main denoted TG      (2) for an input of TGTG (4)
/// g.259GT[6]   TG[6] -> TG[7]      main denoted 12 bases for an input of 14
/// ```
#[test]
fn a_seeded_overlapping_re_phase_carries_one_copy() {
    assert_eq!(
        norm(SEEDED_OVERLAP_CORE, "NC_TEST.1:g.259GT[1]"),
        "NC_TEST.1:g.258_265TG[2]",
        "the TG tract holds one copy more than the GT tract the count was \
         written against, so re-phasing must carry it"
    );
    assert_eq!(
        norm(SEEDED_OVERLAP_CORE, "NC_TEST.1:g.259GT[6]"),
        "NC_TEST.1:g.258_265TG[7]",
        "GT[6] at g.259 denotes 16 bases; only TG[7] over g.258_265 denotes 16"
    );
}

/// **The emitted edit KIND moves.** A repeat whose count sits one unit either
/// side of the reference tract's count renders as `del` or `dup` rather than as
/// a repeat, so re-basing the count by one also re-bases which count renders as
/// which edit. The reference tract here is 4 copies of `TG`, so the `del`/
/// unchanged/`dup` triple sits at re-based counts 3/4/5, i.e. at spelled counts
/// **2/3/4** — where `origin/main`, which did not re-base, put it at 3/4/5.
///
/// Measured, `origin/main` -> this branch:
///
/// ```text
/// g.259GT[2]   Repeat g.258_265TG[2]  ->  Deletion    g.264_265del
/// g.259GT[3]   Deletion g.264_265del  ->  unchanged   g.259GT[3]
/// g.259GT[4]   unchanged g.259GT[4]   ->  Duplication g.264_265dup
/// ```
///
/// Every new value is the one that denotes the input's bases — `origin/main`
/// denoted two bases too few on all three — but a consumer holding the old
/// strings sees a stored `repeat` become a `del` and a `dup`, and a stored `del`
/// become a repeat left alone.
#[test]
fn a_seeded_overlapping_re_phase_can_change_the_emitted_edit_kind() {
    assert_eq!(
        norm(SEEDED_OVERLAP_CORE, "NC_TEST.1:g.259GT[2]"),
        "NC_TEST.1:g.264_265del",
        "re-based to 3 against a 4-copy tract, which renders as a one-unit del"
    );
    assert_eq!(
        norm(SEEDED_OVERLAP_CORE, "NC_TEST.1:g.259GT[3]"),
        "NC_TEST.1:g.259GT[3]",
        "re-based to 4 against a 4-copy tract, which denotes the reference and \
         so leaves the description alone"
    );
    assert_eq!(
        norm(SEEDED_OVERLAP_CORE, "NC_TEST.1:g.259GT[4]"),
        "NC_TEST.1:g.264_265dup",
        "re-based to 5 against a 4-copy tract, which renders as a one-unit dup"
    );
}

/// The property behind all three kind changes, asserted on the **bases** rather
/// than on any spelling, through `hgvs_to_spdi` rather than through the
/// normalizer agreeing with itself.
///
/// Every count in the range moved between `origin/main` and this branch, and on
/// `origin/main` seven of these eight inputs denoted a different sequence after
/// normalization than before it. All eight agree here.
#[test]
fn every_count_over_the_seeded_overlap_tract_preserves_its_denoted_bases() {
    for count in 1..=8 {
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            assert_normalization_preserves_bases(
                SEEDED_OVERLAP_CORE,
                &format!("NC_TEST.1:g.259GT[{count}]"),
                direction,
            );
        }
    }
}

/// The seeded overlap is idempotent, including across the three counts whose
/// edit kind changes — a `del` or `dup` output must not re-enter the repeat
/// path and move again.
#[test]
fn the_seeded_overlap_outputs_are_fixed_points() {
    for count in 1..=8 {
        let input = format!("NC_TEST.1:g.259GT[{count}]");
        let once = norm(SEEDED_OVERLAP_CORE, &input);
        let twice = norm(SEEDED_OVERLAP_CORE, &once);
        assert_eq!(
            once, twice,
            "NOT IDEMPOTENT for {input}\n  once ={once}\n  twice={twice}"
        );
    }
}
