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
    apply, assert_normalizes_preserving, assert_normalizes_preserving_in, normalize, normalize_in,
    sweep_sequences,
};
use ferro_hgvs::ShuffleDirection;

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
    assert_eq!(normalize(TRACT, "TEMPLATE:g.1_2dup"), "TEMPLATE:g.1_9T[11]");
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
    let mut checked = 0usize;
    let mut residual = 0usize;
    let mut skipped = 0usize;
    let mut overlapping: Vec<String> = Vec::new();
    let mut changed: Vec<String> = Vec::new();

    for seq in sweep_sequences(48) {
        for first_start in 2..=13usize {
            for first_len in 1..=2usize {
                let first_end = first_start + first_len - 1;
                let span = if first_len == 1 {
                    format!("{first_start}")
                } else {
                    format!("{first_start}_{first_end}")
                };
                let inserted = &seq[first_start - 1..first_end];
                let firsts = [
                    format!("{span}del"),
                    format!("{span}dup"),
                    format!("{first_start}_{}ins{inserted}", first_start + 1),
                ];
                for first in firsts {
                    for second_start in first_end + 2..=19usize {
                        let base = seq.as_bytes()[second_start - 1] as char;
                        let alt = if base == 'A' { 'G' } else { 'A' };
                        for second in [
                            format!("{second_start}del"),
                            format!("{second_start}{base}>{alt}"),
                        ] {
                            for direction in
                                [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime]
                            {
                                let input = format!("TEMPLATE:g.[{first};{second}]");
                                let output = normalize_in(&seq, &input, direction);
                                checked += 1;
                                let Some(want) = apply(&seq, &input) else {
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
                                    && direction == ShuffleDirection::FivePrime;
                                match apply(&seq, &output) {
                                    None if overlapping.len() < 10 => {
                                        overlapping.push(format!("{seq}: {input} -> {output}"))
                                    }
                                    None => {}
                                    Some(got) if got != want && residual_shape => {
                                        residual += 1;
                                    }
                                    Some(got) if got != want && changed.len() < 10 => {
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
    }

    assert!(checked > 100_000, "sweep covered too little: {checked}");
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
    let seq = "ACAAAAAAAACGTACGTACG";
    for (input, expected) in [
        ("TEMPLATE:g.[4A>G;9dup]", "TEMPLATE:g.[4A>G;5dup]"),
        ("TEMPLATE:g.[5A>G;10dup]", "TEMPLATE:g.[5A>G;6dup]"),
    ] {
        assert_normalizes_preserving_in(seq, input, expected, ShuffleDirection::FivePrime);
    }
}

#[test]
fn a_five_prime_insertion_junction_with_an_upstream_sibling_is_a_known_gap() {
    // The rest of #1267, which `blocks_sibling_shift` deliberately does not
    // reach. That predicate makes a member a barrier because the bases under
    // its span are what it copies — true of a `Duplication`, false of a pure
    // `Insertion`, whose payload is a literal and whose span is zero-width.
    // There is nothing for a sibling to land *on*.
    //
    // Bounding the insertion's **junction** in the 5' direction is the fix this
    // shape needs, and #1259 measured a 5' mirror of the junction clamp turning
    // 80 previously-correct outputs into silently wrong ones. So this stays
    // open, narrowed to the insertion case.
    //
    // Pinned as *currently wrong*: fixing it fails this test loudly.
    let seq = "ACAAAAAAAACGTACGTACG";
    let input = "TEMPLATE:g.[4A>G;9_10insA]";
    let actual = normalize_in(seq, input, ShuffleDirection::FivePrime);
    assert_eq!(actual, "TEMPLATE:g.4_5insG", "known-gap shape changed");
    assert_ne!(
        apply(seq, &actual).expect("output applies"),
        apply(seq, input).expect("input applies"),
        "{input} is fixed — update this test and close the insertion half of #1267"
    );
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
        let actual = normalize(&seq, input);
        assert_eq!(actual, "TEMPLATE:g.259A>C", "normalizing {input}");
        assert_eq!(
            apply(&seq, &actual).expect("output applies"),
            apply(&seq, input).expect("input applies"),
            "{input} -> {actual} changed the sequence"
        );
    }

    // A junction further into the tract clamps to that junction instead.
    let input = "TEMPLATE:g.[258del;261_262insC]";
    let actual = normalize(&seq, input);
    assert_eq!(actual, "TEMPLATE:g.261A>C");
    assert_eq!(
        apply(&seq, &actual).expect("output applies"),
        apply(&seq, input).expect("input applies"),
    );

    // Negative control: with no sibling the deletion still shifts fully.
    assert_eq!(normalize(&seq, "TEMPLATE:g.258del"), "TEMPLATE:g.263del");
}
