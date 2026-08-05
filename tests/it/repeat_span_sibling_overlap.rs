//! A repeat's tract-wide span must not swallow a cis sibling.
//!
//! `deletion_to_repeat` (`repeated.md` B2) re-expresses a deletion inside a
//! tandem tract as a copy count over the **whole tract**, not over the deleted
//! bases. On a nine-`T` run, `g.1_2del` becomes `g.1_9T[7]` — correct for a lone
//! variant, wrong the moment a cis sibling lives inside that tract:
//!
//! ```text
//! reference   TTTTTTTTTAATATATTTTA        (positions 1-9 are T)
//! input       TEMPLATE:g.[1_2del;4del]
//! was      -> TEMPLATE:g.[1_9T[7];9del]   `9del` sits inside `1_9`
//! now      -> TEMPLATE:g.1_9T[6]
//! ```
//!
//! The conversion runs per member, deep inside `normalize_na_edit`, with no
//! sibling context reaching it, so the widened span is undone in the allele loop
//! instead: the repeat is re-spelled as the equivalent 3'-most deletion, the
//! sibling-crossing clamp pulls it back, and the merge finishes the job. Repeat
//! notation is restored once the members have coalesced and there is no sibling
//! left to span.
//!
//! This is the residual left by #1254 and #1234, which fixed the *shift* into a
//! sibling but not the *span expansion* over one — `NaEdit::Repeat` claims no
//! reference bases as far as the clamp is concerned, so it skipped these. It is
//! pre-existing on `main` and independent of both.

use crate::common::cis_apply_oracle::{
    apply, assert_normalizes_preserving, normalize_in, sweep_seeds, sweep_sequences,
};
use ferro_hgvs::ShuffleDirection;

/// Nine `T` at positions 1-9, then `AA` breaking the run.
const TRACT: &str = "TTTTTTTTTAATATATTTTA";

#[test]
fn two_deletions_in_one_tract_reach_a_single_repeat() {
    // The reported case. Three bases leave the nine-`T` tract, so the whole
    // allele is one repeat of six copies — not a seven-copy repeat with a
    // stray deletion inside it.
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[1_2del;4del]", "TEMPLATE:g.1_9T[6]");
}

#[test]
fn three_deletions_in_one_tract_reach_the_same_repeat() {
    // Same total, spelled as three single-base members.
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[1del;5del;9del]", "TEMPLATE:g.1_9T[6]");
}

#[test]
fn a_substitution_inside_the_tract_blocks_the_repeat_form() {
    // A substituted base inside the tract cannot be described by a copy count,
    // so the deletion stays a deletion, clamps against the substitution, and
    // the two coalesce.
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[1_2del;4T>A]", "TEMPLATE:g.2_4delinsA");
}

#[test]
fn a_substitution_at_the_tract_end_blocks_it_too() {
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[1_2del;9T>A]", "TEMPLATE:g.7_9delinsA");
}

#[test]
fn the_repeat_and_deletion_spellings_reach_one_canonical_form() {
    // Confluence. `[1_2del;4T>A]` used to reach the overlapping
    // `[1_9T[7];4T>A]` while `[2_3del;4T>A]` — the same variant, already at the
    // clamp position — reached `2_4delinsA`. Both now land on the latter.
    for input in ["TEMPLATE:g.[1_2del;4T>A]", "TEMPLATE:g.[2_3del;4T>A]"] {
        assert_normalizes_preserving(TRACT, input, "TEMPLATE:g.2_4delinsA");
    }
}

#[test]
fn a_lone_deletion_still_becomes_a_repeat() {
    // Negative control: with no sibling there is nothing to span, so B2 applies
    // exactly as before. Guards against "fixing" this by disabling the repeat
    // form for deletions.
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.1_2del", "TEMPLATE:g.1_9T[7]");
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.1_3del", "TEMPLATE:g.1_9T[6]");
}

#[test]
fn a_sibling_outside_the_tract_leaves_the_repeat_alone() {
    // Negative control: `10A>G` is past the tract's last base (9), so the
    // repeat span claims nothing the sibling claims and must survive.
    assert_normalizes_preserving(
        TRACT,
        "TEMPLATE:g.[1_2del;10A>G]",
        "TEMPLATE:g.[1_9T[7];10A>G]",
    );
}

#[test]
fn a_short_tract_with_a_sibling_against_its_edge_merges() {
    // Was `a_short_tract_with_a_clear_sibling_keeps_its_repeat`, pinning
    // `g.[4_7T[2];8A>G]`. The sibling is not in fact clear of the tract: the
    // deletion 3'-shifts to 6_7 and the substitution at 8 abuts it, so the two
    // are consecutive changes and `delins.md:16` makes them one `delins`.
    //
    // The repeat spelling used to survive only because the sequence-first pass
    // refused any group containing one; with repeats lowered for derivation
    // (#1296) it re-derives the merged form. `TTA` -> `G` over 6_8.
    //
    // The negative control this row was written to be lives on in
    // `a_short_tract_with_a_genuinely_clear_sibling_keeps_its_repeat` below.
    assert_normalizes_preserving(
        "ACGTTTTACGTACGTACGTA",
        "TEMPLATE:g.[5_6del;8A>G]",
        "TEMPLATE:g.6_8delinsG",
    );
}

#[test]
fn a_short_tract_with_a_genuinely_clear_sibling_keeps_its_repeat() {
    // Negative control on the same four-base tract (`T` at 4-7): `5_6del`
    // becomes `4_7T[2]`, and the substitution at 9 is separated from the
    // shifted deletion by the unchanged base at 8 — so there is nothing for the
    // derivation to merge, and the repeat survives.
    //
    // A control, not a discriminator, and measured as such: this row passes
    // whether repeat lowering runs or not, because its input carries no repeat
    // member and its expected output *keeps* one. It guards the opposite risk
    // from the rows above — the derivation over-reaching onto a genuinely
    // separated sibling — which is exactly what it is here to do. The rows that
    // redden when lowering is removed are the five in `issue_1296_*`,
    // `issue_1135_*`, `cis_spelling_confluence_gap` and
    // `a_short_tract_with_a_sibling_against_its_edge_merges`.
    assert_normalizes_preserving(
        "ACGTTTTACGTACGTACGTA",
        "TEMPLATE:g.[5_6del;9C>G]",
        "TEMPLATE:g.[4_7T[2];9C>G]",
    );
}

#[test]
fn a_five_prime_shuffle_is_unaffected() {
    // B2 is gated on 3' shuffle, so under FivePrime no repeat is produced and
    // this pass has nothing to undo: the two deletions shuffle 5' into one
    // three-base deletion at the tract's start.
    //
    // Pinned exactly rather than asserted as "contains no `[`". That weaker
    // form cannot tell repeat notation from allele brackets — a well-formed
    // two-member output renders as `g.[...;...]` and would fail it, while any
    // single-member output passes it however far its span had drifted.
    let out = normalize_in(
        TRACT,
        "TEMPLATE:g.[1_2del;4del]",
        ShuffleDirection::FivePrime,
    );
    assert_eq!(out, "TEMPLATE:g.1_3del");
    assert_eq!(
        apply(TRACT, &out).expect("output applies"),
        apply(TRACT, "TEMPLATE:g.[1_2del;4del]").expect("input applies"),
    );
}

#[test]
fn no_two_member_allele_normalizes_to_overlapping_members() {
    // The class. Every two-member cis allele of a deletion plus a downstream
    // deletion or substitution, over deterministic 20-mers, in both shuffle
    // directions — checked for overlapping output *and* for sequence
    // preservation against an independently applied reference.
    //
    // Overlapping outputs across this sweep: 2,578 before the sibling-crossing
    // clamp, 526 after it (all of them this repeat-span class), 0 now.
    let mut checked = 0usize;
    let mut skipped = 0usize;
    let mut overlapping: Vec<String> = Vec::new();
    let mut changed: Vec<String> = Vec::new();

    // #1295: full 64 when CI asks (`FERRO_SWEEP_SEEDS=full`), a 4-seed prefix
    // otherwise. This sweep was 33.1s of an 86.6s local suite.
    let seeds = sweep_seeds(64);
    for seq in sweep_sequences(seeds) {
        for first_start in 1..=14usize {
            for first_len in 1..=2usize {
                let first_end = first_start + first_len - 1;
                let first = if first_len == 1 {
                    format!("{first_start}del")
                } else {
                    format!("{first_start}_{first_end}del")
                };
                for second_start in first_end + 2..=19usize {
                    let base = seq.as_bytes()[second_start - 1] as char;
                    let alt = if base == 'A' { 'G' } else { 'A' };
                    for second in [
                        format!("{second_start}del"),
                        format!("{second_start}{base}>{alt}"),
                    ] {
                        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime]
                        {
                            let input = format!("TEMPLATE:g.[{first};{second}]");
                            let output = normalize_in(&seq, &input, direction);
                            checked += 1;
                            let Some(want) = apply(&seq, &input) else {
                                skipped += 1; // unconvertible input
                                continue;
                            };
                            match apply(&seq, &output) {
                                // `apply` declines an overlapping output, so
                                // `None` here is the overlap signal.
                                None if overlapping.len() < 10 => {
                                    overlapping.push(format!("{seq}: {input} -> {output}"));
                                }
                                None => {}
                                Some(got) if got != want && changed.len() < 10 => {
                                    changed.push(format!(
                                        "{seq}: {input} -> {output} (want {want}, got {got})"
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

    // Per seed since #1295, for the reason recorded on the equivalent floor in
    // `cis_junction_crossing_shift.rs`: the seed count is a knob, so an absolute
    // floor would either fail at the default prefix or go vacuous at the full
    // corpus.
    const CASES_PER_SEED_FLOOR: usize = 1_500;
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
}
