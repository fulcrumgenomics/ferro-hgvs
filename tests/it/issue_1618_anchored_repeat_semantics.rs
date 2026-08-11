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
