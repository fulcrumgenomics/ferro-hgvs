//! Regression coverage for issue #1254 — a cis member's 3'-shift must not
//! carry it over a sibling's bases.
//!
//! The per-member 3'-shift is sibling-unaware: each member goes to its
//! *standalone* most-3' position. When a sibling edits a base inside the tract
//! the member rotates through, the member leapfrogs it and the pair describes a
//! different sequence:
//!
//! ```text
//! reference   TAATATATATAATATATATT
//! input       TEMPLATE:g.[3_4del;9del]
//! was      -> TEMPLATE:g.12_14del      applied: TAATATATATATATATT
//! input applied                                 TAATATTAATATATATT   ← different
//! now      -> TEMPLATE:g.7_9del        applied: TAATATTAATATATATT
//! ```
//!
//! `3_4del` shifted to `10_11del`, straight past the `9del`; the two were then
//! adjacent and merged to `9_11del`, which shifted on to `12_14del`. Nothing
//! about that output invites suspicion — one well-formed member, no overlap, no
//! warning — so the corruption is silent.
//!
//! Every assertion here is against an **independently applied** reference
//! sequence, reconstructed from SPDI triples rather than from the normalizer.
//! Asserting `normalize(input) == normalize(output)` cannot catch this class: it
//! normalizes both sides, so it passed on the broken behavior.

use crate::common::cis_apply_oracle::{
    apply, assert_normalizes_preserving, assert_normalizes_preserving_in, normalize,
    sweep_sequences,
};
use ferro_hgvs::ShuffleDirection;

/// The reported reference: an `(AT)` tract that ends at position 11, where
/// `11=A, 12=A` breaks the repeat.
const TEMPLATE: &str = "TAATATATATAATATATATT";

#[test]
fn separated_deletions_do_not_merge_into_a_different_deletion() {
    // The reported case. `3_4del` may 3'-shift only as far as `7_8del`, which
    // stops one base short of the `9del`; the two are then adjacent (no
    // intervening nucleotide) and merge into one three-base deletion.
    assert_normalizes_preserving(TEMPLATE, "TEMPLATE:g.[3_4del;9del]", "TEMPLATE:g.7_9del");
}

#[test]
fn separated_deletions_with_a_downstream_substitution_stay_well_formed() {
    // The three-member variant of the same input. It used to expose the
    // corruption *visibly*, as the overlapping `g.[12_14del;13T>A]`; with the
    // shift clamped, the deletion lands where it belongs and the substitution
    // is untouched.
    assert_normalizes_preserving(
        TEMPLATE,
        "TEMPLATE:g.[3_4del;9del;13T>A]",
        "TEMPLATE:g.[7_9del;13T>A]",
    );
}

#[test]
fn every_spelling_of_the_variant_reaches_one_canonical_form() {
    // Confluence: the pre-shifted, the shifted-and-adjacent, the already-merged
    // and the reordered spellings are one variant and must share a key.
    for input in [
        "TEMPLATE:g.[3_4del;9del]",
        "TEMPLATE:g.[7_8del;9del]",
        "TEMPLATE:g.[9del;3_4del]",
        "TEMPLATE:g.7_9del",
    ] {
        assert_normalizes_preserving(TEMPLATE, input, "TEMPLATE:g.7_9del");
    }
}

#[test]
fn a_deletion_does_not_shift_onto_a_substituted_base() {
    // A deletion inside a homopolymer whose 3'-most standalone position is the
    // substituted base itself. Clamped to `8del`, it is adjacent to `9A>G` and
    // the two render as one delins (`delins.md:16`). This previously emitted the
    // overlapping `g.[9A>G;9del]`.
    assert_normalizes_preserving(
        "GGGAAAAAAGGG",
        "TEMPLATE:g.[4del;9A>G]",
        "TEMPLATE:g.8_9delinsG",
    );
}

#[test]
fn a_member_that_cannot_reach_its_sibling_still_shifts_fully() {
    // Negative control against over-clamping: `4del` 3'-shifts through the A-run
    // to `9del`, and the substitution at `11` is two bases beyond it — never
    // inside the swept window — so the shift must complete and the two must stay
    // separate members, separated by the unchanged `10G` (`general.md:34`).
    assert_normalizes_preserving(
        "GGGAAAAAAGCG",
        "TEMPLATE:g.[4del;11C>T]",
        "TEMPLATE:g.[9del;11C>T]",
    );
}

#[test]
fn a_five_prime_shift_does_not_cross_a_sibling_either() {
    // The mirror image. Under `ShuffleDirection::FivePrime`, `9_10del` rotates
    // 5' through the same `(AT)` tract and its standalone 5'-most position is
    // past the `5del`. Clamped to `6_7del` it is adjacent to the `5del`, and the
    // two merge into `g.5_7del`. This previously emitted `g.1_3del`, which
    // denotes `TATATATAATATATATT` where the input denotes `TAATTATAATATATATT`.
    assert_normalizes_preserving_in(
        TEMPLATE,
        "TEMPLATE:g.[5del;9_10del]",
        "TEMPLATE:g.5_7del",
        ShuffleDirection::FivePrime,
    );
}

// The #999 negative control for this clamp — an insertion claims no reference
// base, so a member landing flush against one is adjacency, not a crossing —
// lives once, as `cis_junction_crossing_shift::
// an_insertion_may_still_land_flush_against_its_sibling`. Both clamps run on
// that one input, so a second identical copy here would only be a second thing
// to keep in step.

#[test]
fn shifts_never_change_the_sequence_across_an_exhaustive_two_member_sweep() {
    // The class, not the instance. Every two-member cis allele of a deletion
    // plus a downstream deletion or substitution, over 128 deterministic
    // 20-mers (64 xorshift seeds x the {A,T} and {A,C,G,T} alphabets), checked
    // by applying both sides to the reference. Before the clamp this sweep
    // reported silent sequence-changing normalizations — well-formed,
    // warning-free, disjoint output — in the low thousands; over-clamping would
    // show up as a canonical-form change in the pinned cases above.
    //
    // Two failure modes are counted separately, matching the sibling sweep in
    // `cis_junction_crossing_shift.rs`: an *input* that does not apply is a case
    // this sweep cannot speak for and is skipped, while an *output* that does
    // not apply is a failure — an overlapping or unconvertible normalization is
    // precisely the defect the sweep exists to catch, and folding it into
    // `skipped` let it pass silently.
    let mut enumerated = 0usize;
    let mut skipped = 0usize;
    let mut overlapping: Vec<String> = Vec::new();
    let mut mismatched: Vec<String> = Vec::new();

    for seq in sweep_sequences(64) {
        for first_start in 1..=14usize {
            for first_len in 1..=2usize {
                let first_end = first_start + first_len - 1;
                let first = if first_len == 1 {
                    format!("{first_start}del")
                } else {
                    format!("{first_start}_{first_end}del")
                };
                // Leave at least one unchanged base between the members, so
                // every case starts out as a genuinely separated pair.
                for second_start in first_end + 2..=19usize {
                    let base = seq.as_bytes()[second_start - 1] as char;
                    let alt = if base == 'A' { 'G' } else { 'A' };
                    for second in [
                        format!("{second_start}del"),
                        format!("{second_start}{base}>{alt}"),
                    ] {
                        let input = format!("TEMPLATE:g.[{first};{second}]");
                        let output = normalize(&seq, &input);
                        // Counted before the skip, so the floor below
                        // measures enumeration and not how many cases
                        // happened to be convertible — otherwise a small
                        // change in the skip rate trips a coverage floor
                        // for a reason unrelated to coverage.
                        enumerated += 1;
                        let Some(want) = apply(&seq, &input) else {
                            skipped += 1; // the input itself is unconvertible
                            continue;
                        };
                        let Some(got) = apply(&seq, &output) else {
                            if overlapping.len() < 10 {
                                overlapping.push(format!("{seq}: {input} -> {output}"));
                            }
                            continue;
                        };
                        if want != got && mismatched.len() < 10 {
                            mismatched.push(format!(
                                "{seq}: {input} -> {output} (want {want}, got {got})"
                            ));
                        }
                    }
                }
            }
        }
    }

    // Floor on the enumeration, so it measures the generator rather than the
    // convertible fraction. The skip count is asserted separately: it is the
    // number this sweep cannot speak for, and letting it grow silently would
    // hollow the property out without moving `enumerated` at all.
    assert!(
        enumerated > 70_000,
        "sweep covered too little: {enumerated} enumerated"
    );
    assert!(
        skipped * 10 < enumerated,
        "too many cases skipped as unconvertible: {skipped} of {enumerated}"
    );
    assert!(
        overlapping.is_empty(),
        "overlapping or unconvertible output in {enumerated} enumerated cases: \
         {overlapping:#?}"
    );
    assert!(
        mismatched.is_empty(),
        "sequence-changing normalizations in {enumerated} enumerated cases \
         ({skipped} skipped): {mismatched:#?}"
    );
}
