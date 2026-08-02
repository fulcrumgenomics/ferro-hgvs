//! Acceptance corpus for the cis-allele normalizer's sequence-first rewrite (v1, genomic
//! axis).
//!
//! The rewrite (#1235) replaces "normalize each member independently, then repair the
//! damage with seven passes" with "apply the allele, align the result to the reference,
//! partition at the alignment's dominators, render." This module pins the open defects
//! that rewrite is meant to fix, and states in one place which open defects it is **not**
//! meant to fix.
//!
//! ## PARTITION targets asserted here
//!
//! Every assertion below is sequence-level — never a hardcoded "correct" output string,
//! since these are open defects and nobody has stated what the fixed rendering should be.
//! Two forms:
//!
//! - *Confluence* (#1260, #1262): two spellings of one variant must reach the same
//!   normalized string. They don't, today. `the_two_confluence_targets_still_diverge`
//!   asserts they still don't; `every_confluence_target_denotes_one_variant` guards that
//!   both spellings really are the same variant, checked by `cis_apply_oracle::apply`,
//!   independent of the normalizer under test.
//! - *Denotation* (#1267, #1325): normalizing must not change what the description
//!   denotes. It does, today. `assert_denotation_currently_broken` asserts the normalized
//!   output either denotes a different sequence than the input (#1267), or denotes no
//!   sequence at all (#1325) because its members now overlap and the independent applier
//!   declines to splice them.
//!
//! Every test here is a **pinned failure**, not a `#[ignore]`d one: it currently passes
//! *because* it asserts the defect reproduces. The moment the rewrite fixes a case, the
//! assertion flips and the test goes red — that's the signal to delete it (or move it to
//! a "fixed" note), never something to silence.
//!
//! ## PARTITION issues not asserted
//!
//! - **#1235** — the tracking/design issue for the whole rewrite. It supplies no
//!   reproducer of its own (a table of one-line sketches, not full HGVS strings); its
//!   acceptance criteria (confluence, no overlapping/out-of-order members, idempotence)
//!   are exactly what `cis_allele_confluence_proptest.rs` and
//!   `issue_1235_cis_allele_confluence.rs` already exercise. Nothing new to assert here.
//! - **#1271** — the issue's own worked example (the HGVS spec's `LRG_199` delins) is
//!   **not a live defect** on this branch: the regression that used to split it four ways
//!   was already fixed (a regime-aware `MAX_UNGUARDED_SPLIT_BLOCK`, per the issue's own
//!   account). What #1271 tracks is a *principled* follow-up (extending
//!   `separations_are_meaningful` to net deletions) with no constructed failing case yet.
//!   Its illustrative example is also stated as a coding-axis (`c.`) description, and the
//!   split spelling's coordinates are relative/illustrative rather than full genomic
//!   positions — there is no sequence content to build a synthetic genomic reproducer from
//!   without inventing it. Excluded rather than forced into a frame it doesn't fit.
//!
//! ## BOUNDARY — explicitly not this rewrite's targets
//!
//! These four are contig-bounds, arithmetic-underflow, or circular-wraparound bugs, not
//! member-partitioning bugs; the rewrite's apply/align/partition/render pipeline has no
//! reason to touch any of them. Documented here so the exclusion is a decision, not an
//! oversight. None of the four are asserted below.
//!
//! - **#1274** — a cis insertion+deletion that cancels at a contig's last base emits a
//!   coordinate one past the end (`g.10_11=` on a 10-base contig). The sequence itself is
//!   right; only the stated span is out of bounds. A clamp/validation fix, not a
//!   partition fix.
//! - **#1282** — 5'-shifting a member to position 1 underflows `hgvs_pos_to_index`'s
//!   `(pos - 1)` and panics in debug builds (silently wraps to a garbage index in
//!   release). An input-validation/guard fix at the coordinate layer, not a partition fix.
//! - **#1307** — respelling a duplication that ends at a contig's last base as an
//!   insertion places it at a gap that does not exist (`24_25insC` on a 24-base contig). A
//!   bounds check on `respell_at_gap`, not a partition fix.
//! - **#1327** — the same `respell_at_gap` gap-placement defect as #1307, but on the
//!   mitochondrial/circular (`m.`) axis: a junction at the contig's last base should wrap
//!   to `16569_1` instead of running off the end. No concrete reproducer exists (a
//!   code-level report only); circular wraparound is also out of a genomic-only v1's scope
//!   regardless of reproducibility.
//!
//! ## OTHER
//!
//! - **#1270** — an unconfirmed coding-axis (`c.`) codon-frame-gate asymmetry between the
//!   deletion-to-repeat and insertion-to-repeat paths. The issue's own author states no
//!   case has been constructed that reaches it ("an asymmetry with a plausible consequence
//!   rather than a demonstrated defect"). Not genomic, not reproduced, nothing to assert.
//!
//! ## An unexpected pass, found while building this corpus
//!
//! #1267's issue body gives three lines over reference `ACAAAAAAAACGTACGTACG`. Only the
//! insertion form asserted below (`g.[4A>G;9_10insA]` -> `g.4_5insG`) still reproduces on
//! this branch. The two duplication-form lines the issue also lists —
//! `g.[4A>G;9dup]` -> `g.[4A>G;5dup]` and `g.[5A>G;10dup]` -> `g.[5A>G;6dup]` — were
//! independently checked here (via this module's own applier) and already normalize
//! correctly, preserving the denoted sequence; `blocks_sibling_shift` fixed them, and
//! `cis_junction_crossing_shift.rs`'s `a_five_prime_duplication_does_not_cross_an_upstream_sibling`
//! already pins them as passing. They are not asserted here as failures — reported here so
//! the unexpected pass isn't silently dropped, per the corpus's own rules.

use crate::common::cis_apply_oracle::{apply, normalize, normalize_in};
use crate::common::synthetic::padded;
use ferro_hgvs::ShuffleDirection;

/// `(issue, core, split spelling, spanning/merged spelling)` — the two spellings that
/// currently reach two different fixed points instead of one.
const CONFLUENCE_TARGETS: &[(&str, &str, &str, &str)] = &[
    (
        "#1260",
        "AAAAAA",
        "TEMPLATE:g.[258_259insC;259_260insC]",
        "TEMPLATE:g.258_260delinsACACA",
    ),
    (
        "#1262",
        "AAAAAA",
        "TEMPLATE:g.[258A>C;260del]",
        "TEMPLATE:g.258_260delinsCA",
    ),
];

/// Both spellings in every row must denote the same sequence, or the row is a broken pin
/// rather than evidence about the normalizer (this has happened before in this campaign).
#[test]
fn every_confluence_target_denotes_one_variant() {
    for (issue, core, a, b) in CONFLUENCE_TARGETS {
        let seq = padded(core);
        let from_a = apply(&seq, a).unwrap_or_else(|| panic!("{issue}: `{a}` does not apply"));
        let from_b = apply(&seq, b).unwrap_or_else(|| panic!("{issue}: `{b}` does not apply"));
        assert_eq!(
            from_a, from_b,
            "{issue}: `{a}` and `{b}` are not the same variant"
        );
    }
}

/// #1260 and #1262: pinned as *currently diverging*. Fixing either fails this test loudly.
#[test]
fn the_two_confluence_targets_still_diverge() {
    for (issue, core, a, b) in CONFLUENCE_TARGETS {
        let seq = padded(core);
        let (norm_a, norm_b) = (normalize(&seq, a), normalize(&seq, b));
        assert_ne!(
            norm_a, norm_b,
            "{issue} appears fixed — `{a}` and `{b}` now both normalize to `{norm_a}`. \
             Move this case out of the PARTITION target list and delete this test."
        );
    }
}

/// Assert that normalizing `input` in `direction` currently changes what it denotes: the
/// output either splices to a different sequence than the input, or (when its members now
/// overlap) does not splice to a well-defined sequence at all. Flips to a failing
/// assertion, on purpose, the moment the rewrite fixes the case — that is the point of
/// this corpus.
fn assert_denotation_currently_broken(
    seq: &str,
    input: &str,
    direction: ShuffleDirection,
    issue: &str,
) {
    let actual = normalize_in(seq, input, direction);
    let from_input = apply(seq, input)
        .unwrap_or_else(|| panic!("{issue}: `{input}` has no well-defined resulting sequence"));
    let from_output = apply(seq, &actual);
    assert_ne!(
        from_output,
        Some(from_input),
        "{issue} appears fixed — `{input}` now normalizes to `{actual}` (under {direction:?} \
         shuffle), which denotes the same sequence as the input. Move this case out of the \
         PARTITION target list and delete this test."
    );
}

/// #1267: a 5'-shifting insertion's junction crosses an upstream sibling substitution,
/// moving the base the sibling edited — well-formed, disjoint, warning-free, and wrong.
#[test]
fn a_five_prime_insertion_junction_still_crosses_an_upstream_sibling() {
    let seq = "ACAAAAAAAACGTACGTACG";
    assert_denotation_currently_broken(
        seq,
        "TEMPLATE:g.[4A>G;9_10insA]",
        ShuffleDirection::FivePrime,
        "#1267",
    );
}

/// #1325: two insertions collapse into a repeat correctly, but adding a third pushes a
/// junction inside the tract; the tract grew by more bases than it can express as a
/// duplication, so the demotion pass declines and the repeat is left spanning (swallowing)
/// the sibling's junction. The independent applier refuses to splice the resulting
/// members at all, since they overlap — there is no single well-defined resulting sequence
/// for them, which is itself the defect this pins (#1235's criterion 2).
#[test]
fn a_repeat_growth_exceeding_its_tract_still_swallows_a_sibling_junction() {
    let seq = padded("GATCATAAATTCAGC");
    assert_denotation_currently_broken(
        &seq,
        "TEMPLATE:g.[262_263insAA;263_264insAA;264_265insC]",
        ShuffleDirection::ThreePrime,
        "#1325",
    );
}
