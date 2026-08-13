//! The confluence pairs reported from a downstream pipeline (#1419/#1420/#1421).
//!
//! # Why these are committed rather than left in the issues
//!
//! Each row is two spellings of **one** variant that ferro normalizes to two
//! different strings. They are the only externally-reported instances of the
//! #1235 defect class, they were found by running a real pipeline rather than by
//! a generator, and until this module existed **not one of them was guarded by a
//! test**. Filed is not the same as guarded: an issue records that something was
//! once broken, but nothing re-checks it, so a change that silently reintroduces
//! the defect passes CI and the report has to be made again.
//!
//! They also cover a shape the generated corpus in `cis_confluence_axis` does
//! not reach — #1421's net insertions, where the payload is longer than the span
//! it replaces.
//!
//! # What this pins, and what it deliberately does not
//!
//! It pins **convergence** — the two spellings must produce one output — and
//! **sequence preservation**, checked through the SPDI applier in
//! [`crate::common::cis_apply_oracle`] rather than through the normalizer, so a
//! pass that converged on a wrong sequence fails here instead of looking like a
//! success.
//!
//! It does **not** pin *which* string a pair converges on. That is an open
//! product decision (#1235), and pinning a winner here would freeze it by
//! accident. The census reports the current answer per row when it moves, so a
//! change is visible in review without being asserted as correct.
//!
//! # The per-spelling answers are pinned next door, in both directions
//!
//! Declining to pin a *winner* is not the same as declining to pin an *answer*,
//! and the two were conflated while this was the only module in the family.
//! `reported_partition_verdicts` now pins what each of the eighteen spellings
//! prints — under 3' *and* under 5' — together with the form its issue asks for
//! and what backs that form under the project's precedence policy
//! (**spec-explicit > Mutalyzer > our judgement**). Six of those eighteen
//! spellings — three pairs — have a spec line that answers them outright, which
//! is the split `the_spec_authority_census_holds` pins next door as `(6, 12)`
//! rows; #1419's six and #1421's six do not, because `delins.md:44-47`
//! recommends the merged delins for the same alignment coincidence each pair's
//! split relies on (`general.md:56` for #1419, the unchanged interior
//! reappearing in the payload for #1421) — so those stay recorded rather than
//! targeted.
//!
//! So this module keeps exactly one job: the **ratchet**.
//! [`CONVERGING_PAIRS_THREE_PRIME`] and [`CONVERGING_PAIRS_FIVE_PRIME`] are the
//! numbers that must only ever go up, and each is deliberately a count and not
//! nine strings, because the goal is total convergence and partial progress
//! should read as one number moving. They are tracked separately per direction
//! rather than as one shared constant, since the two directions are not
//! guaranteed to converge the same rows.
//!
//! # Current status
//!
//! [`CONVERGING_PAIRS_THREE_PRIME`] and [`CONVERGING_PAIRS_FIVE_PRIME`] record
//! how many converge today (currently 0 under both directions). Both spellings
//! of every row are well-formed and sequence-preserving, so a split row is a pure
//! representation difference rather than a correctness bug — which is precisely
//! why it is damaging downstream: key-based aggregation silently files one
//! variant under two keys and halves its count.

use crate::common::cis_apply_oracle::{apply, normalize_in};
use ferro_hgvs::ShuffleDirection;

/// A 125 nt contig. Position 1 is offset 0 — the oracle's provider serves the
/// core unpadded under the accession `TEMPLATE`.
///
/// Reproduced verbatim from the reports so these rows run against the same bases
/// they were observed on. Synthetic: it names no real sequence.
pub(crate) const TEMPLATE: &str = concat!(
    "ATGCACCAGTCACCAGTCTGATGCGGATCACGTGCAATTGCACGTGCAATTGGATCCGATCG",
    "TACGTACGATCGGCATGCATGCTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGA"
);

/// `(label, spelling A, spelling B)`.
///
/// A is the multi-member form as authored downstream; B is the spanning form of
/// the same variant. The two are kept distinguishable rather than collapsed into
/// an unordered pair because *which one a consumer has stored* decides what a
/// convergence fix costs them.
pub(crate) const REPORTED_PAIRS: &[(&str, &str, &str)] = &[
    // #1419 — two deletions versus the spanning delins.
    (
        "1419-r1",
        "TEMPLATE:g.[19_23del;27_33del]",
        "TEMPLATE:g.19_33delinsCGG",
    ),
    (
        "1419-r2",
        "TEMPLATE:g.[19_22del;26_36del]",
        "TEMPLATE:g.19_36delinsGCG",
    ),
    (
        "1419-r3",
        "TEMPLATE:g.[19_24del;28_33del]",
        "TEMPLATE:g.19_33delinsGGA",
    ),
    // #1420 — a dup, an insertion and a delins, each beside a deletion.
    (
        "1420-v2",
        "TEMPLATE:g.[37dup;41del]",
        "TEMPLATE:g.38_41delinsATTG",
    ),
    (
        "1420-v3",
        "TEMPLATE:g.[36_37insC;40del]",
        "TEMPLATE:g.37_40delinsCATT",
    ),
    (
        "1420-v4",
        "TEMPLATE:g.[21delinsGC;24del]",
        "TEMPLATE:g.21_24delinsGCTG",
    ),
    // #1421 — net insertions: the payload is longer than the span it replaces.
    (
        "1421-n1",
        "TEMPLATE:g.[29C>A;32_33delinsACATACTG]",
        "TEMPLATE:g.29_33delinsAACACATACTG",
    ),
    (
        "1421-n2",
        "TEMPLATE:g.[32G>T;35_36delinsGAATCGAC]",
        "TEMPLATE:g.32_36delinsTTGGAATCGAC",
    ),
    (
        "1421-n3",
        "TEMPLATE:g.[34G>T;37_38delinsCCTTTACG]",
        "TEMPLATE:g.34_38delinsTCACCTTTACG",
    ),
];

/// How many reported pairs converge today under `ShuffleDirection::ThreePrime`.
///
/// **This may only ever go up.** A drop means a change has regressed a defect
/// somebody outside the project reported, and must not be re-blessed silently.
///
/// Kept as a separate constant from [`CONVERGING_PAIRS_FIVE_PRIME`] — the two
/// directions reach the same partitioner but are not guaranteed to converge the
/// same rows, and one shared constant cannot state two different true counts:
/// whichever direction differed would have to be re-blessed against the other's
/// number. Not because sharing would *hide* a single-direction regression — it
/// would not. [`the_reported_pair_census_is_unchanged`] asserts inside the
/// per-direction loop, so each direction is already compared against the
/// constant on its own and no sum is ever taken. Each direction's current value
/// is its own constant; the change sections below record how it got there, and
/// no figure is restated here to drift away from them.
///
/// # A rise is not automatically progress — read this before raising either
///
/// The number going up looks like the goal and can be the opposite of it. A
/// pair converges when both spellings reach one string; that string can be one
/// **neither spelling asserted**, arrived at by re-partitioning across bases the
/// input left unchanged. `1420-v2` is the specific hazard: its cis spelling
/// `g.[37dup;41del]` has three unchanged nucleotides between its members (38,
/// 39, 40), and a re-derivation that merges them lands on
/// `g.[38T>A;40_41delinsTG]` — which happens to be exactly the form #1420 asks
/// for. That coincidence is how a defect gets banked as a fix, and it has been
/// observed on an experimental partitioner.
///
/// A guard below, `the_1420_v2_pair_does_not_converge_by_re_derivation`, used to
/// name that string so the coincidence could not happen quietly. **It was
/// deleted by operator ruling of 2026-08-13** — see the `3 -> 9 (#1616)` section
/// on [`CONVERGING_PAIRS_THREE_PRIME`]. Its ground was that `1420-v2`'s members
/// "are separated by three unchanged nucleotides", which reads separation off the
/// *input's spelling*; the decided
/// `rulings[separation-is-a-property-of-the-spelling-not-of-the-variant]` reads
/// it off the partition re-derived from the sequence, where the output is two
/// members at separation one, described individually. **The hazard this
/// paragraph describes is real and unretracted** — what changed is that this
/// particular row is not an instance of it.
///
/// So a rise must name **which** pair moved and **which clause** carried it. A
/// rise that cannot name a clause is a re-derivation, not a convergence fix.
///
/// # 0 -> 3 (#1649), the pairs and the clause
///
/// The three that converged are **`1419-r1`, `1419-r2` and `1419-r3`**, in both
/// directions — so `CONVERGING_PAIRS_FIVE_PRIME` moves by the same 3, and the
/// hazard row `1420-v2` is **not** among them.
///
/// In each pair it is the `/span` spelling that moved, onto the two-deletion
/// form its `/cis` sibling already printed; the `/cis` spellings are unchanged,
/// which is checkable from `reported_partition_verdicts`' per-row pins rather
/// than taken on trust. The clause is `general.md:34` read through
/// `canonical-form-choice-when-both-legal` (`decided`): the partition is
/// re-derived from the **resulting sequence** and what falls out is emitted.
/// #1649 lets `merge`'s payload splitter express *deletion, retained reference,
/// deletion* — the mirror of the two-insertion shape it already expressed — so
/// the form the sequence supports is now spellable and both spellings reach it.
///
/// **This is a convergence, not a closure, and the distinction is the whole
/// reason that census lives in the other module.** The string the three pairs
/// converge on is still not the spanning `delins` the decided chain wants, so
/// all three stay `PairState::NeitherReaches` in `reported_partition_verdicts`'
/// `PAIR_STATES` and its `OPEN_GAPS` stays at twelve. What moved is that the
/// family stopped having two canonical forms.
///
/// # 3 -> 9 (#1616), and the two clauses that carry it
///
/// **All nine reported pairs now converge; `still split` is 0.** This counter's
/// own failure message requires a rise to name the clause behind it, because "a
/// rise with no clause behind it is a re-derivation onto a form neither spelling
/// asserted". Two carry this one:
///
/// * `rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]` (decided)
///   deleted the input-relative weight bound, so a multi-member spelling is
///   re-derived from the resulting sequence instead of handed back. That is what
///   lets the two spellings of each pair meet at all.
/// * `rulings[canonical-form-choice-when-both-legal]` (decided) says which form
///   they meet on: derive from the resulting sequence and emit what falls out.
///
/// **`1420-v2` is among them, and it was previously guarded against.** The guard
/// `the_1420_v2_pair_does_not_converge_by_re_derivation` was deleted by operator
/// ruling of 2026-08-13 — see the block comment where it stood. Its ground was
/// that the pair's members "are separated by three unchanged nucleotides", which
/// reads separation off the INPUT'S SPELLING; the decided
/// `rulings[separation-is-a-property-of-the-spelling-not-of-the-variant]` says
/// it is read off the partition re-derived from the sequence, and read that way
/// the output `g.[38T>A;40_41delinsTG]` is two members at separation one,
/// described individually. **Nothing merged.**
///
/// Measured, both directions, not composed: 9 converged / 0 split.
const CONVERGING_PAIRS_THREE_PRIME: usize = 9;

/// How many reported pairs converge today under `ShuffleDirection::FivePrime`.
///
/// See [`CONVERGING_PAIRS_THREE_PRIME`] for why this is tracked separately
/// rather than shared, and for why a rise in either direction has to name the
/// clause that carried it. This constant's current value is its own; the change
/// sections below record how it got there, and no figure is restated here to
/// drift away from them. #1649 raised it from 0 to 3 — the three `1419-r*`
/// pairs, same as 3'.
///
/// The two directions agreeing here is worth reading rather than skipping: the
/// pairs converge under 5' on a *different* string in `1419-r3`'s case, because
/// the 5' pass shifts that pair's leading deletion one base left. Convergence is
/// a property of the partition the splitter can express, and that is what #1649
/// changed; which end of the run the members then shuffle to is downstream of
/// it and does not decide whether the two spellings meet.
///
/// # 3 -> 9 (#1616)
///
/// Raised with 3' and by the same two clauses — see
/// [`CONVERGING_PAIRS_THREE_PRIME`], which carries the reasoning. **Measured
/// here rather than assumed to match 3'**, because this constant exists
/// precisely so the two directions cannot be inferred from one another: run
/// separately, 5' also reports 9 converged and 0 still split.
const CONVERGING_PAIRS_FIVE_PRIME: usize = 9;

/// Both directions, because a rule that converges only under the default 3'
/// direction has not solved the problem: 5' is a supported option and reaches
/// the same partitioner.
const DIRECTIONS: [ShuffleDirection; 2] =
    [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime];

/// The premise every row rests on: the two spellings really are one variant.
///
/// Asserted first and separately. If a pair ever stopped denoting one sequence,
/// the census below would be measuring nothing — two genuinely different
/// variants are *supposed* to normalize differently.
#[test]
fn every_reported_pair_denotes_one_sequence() {
    for (label, a, b) in REPORTED_PAIRS {
        let left = apply(TEMPLATE, a);
        let right = apply(TEMPLATE, b);
        assert!(
            left.is_some(),
            "{label}: {a} does not apply to the reference"
        );
        assert!(
            right.is_some(),
            "{label}: {b} does not apply to the reference"
        );
        assert_eq!(
            left, right,
            "{label}: the two spellings denote different sequences, so this row \
             is not a confluence pair and the census is measuring nothing"
        );
    }
}

/// Whatever each spelling normalizes to must still denote the sequence it
/// started from.
///
/// Asserted outright rather than counted: converging on a *wrong* sequence would
/// be worse than the divergence this module tracks, so there is no acceptable
/// non-zero value.
#[test]
fn no_reported_pair_normalizes_to_a_different_sequence() {
    for (label, a, b) in REPORTED_PAIRS {
        let expected = apply(TEMPLATE, a).expect("pair applies");
        for direction in DIRECTIONS {
            for input in [a, b] {
                let output = normalize_in(TEMPLATE, input, direction);
                assert_eq!(
                    apply(TEMPLATE, &output).as_deref(),
                    Some(expected.as_str()),
                    "{label}: under {direction:?}, normalizing {input} to {output} \
                     changed the sequence it denotes"
                );
            }
        }
    }
}

/// The census, and the ratchet.
///
/// A count rather than nine separate assertions: the goal is total convergence,
/// so partial progress should show up as one number moving rather than as nine
/// test-file edits. The failure message prints every row's current pair of
/// outputs, which is what a reader needs to decide whether a move is progress or
/// a regression.
#[test]
fn the_reported_pair_census_is_unchanged() {
    for direction in DIRECTIONS {
        // The catch-all arm is required, not defensive padding, and deleting it
        // for "compile-time exhaustiveness" does not compile: `ShuffleDirection`
        // is `#[non_exhaustive]` (see its own doc comment), so from outside the
        // defining crate — which `tests/it` is — a `match` on it must carry a
        // wildcard however many variants it has today. A runtime panic is
        // therefore the strongest check available here, and it is reachable only
        // via `DIRECTIONS`, so a new variant must be added there too before this
        // fires. Prefer widening both together over relying on the panic.
        let expected = match direction {
            ShuffleDirection::ThreePrime => CONVERGING_PAIRS_THREE_PRIME,
            ShuffleDirection::FivePrime => CONVERGING_PAIRS_FIVE_PRIME,
            _ => panic!(
                "unhandled shuffle direction {direction:?} — add a per-direction census constant"
            ),
        };
        let mut converged: Vec<String> = Vec::new();
        let mut split: Vec<String> = Vec::new();
        for (label, a, b) in REPORTED_PAIRS {
            let left = normalize_in(TEMPLATE, a, direction);
            let right = normalize_in(TEMPLATE, b, direction);
            if left == right {
                converged.push(format!("{label} -> {left}"));
            } else {
                split.push(format!(
                    "{label}\n      A {a}\n        -> {left}\n      B {b}\n        -> {right}"
                ));
            }
        }
        assert_eq!(
            converged.len(),
            expected,
            "reported-pair convergence moved under {direction:?}.\n\n\
             converged ({}):\n    {}\n\nstill split ({}):\n    {}\n\n\
             If this went UP, raise CONVERGING_PAIRS_THREE_PRIME/CONVERGING_PAIRS_FIVE_PRIME \
             (whichever direction moved) and say so in the PR — it is a \
             representation change for anyone storing the losing spelling — AND \
             name the clause that carried the move. A rise with no clause behind \
             it is a re-derivation onto a form neither spelling asserted, which \
             is not progress even when it lands on the form the issue asks for. \
             If it went DOWN, a change has regressed an externally-reported \
             defect.",
            converged.len(),
            converged.join("\n    "),
            split.len(),
            split.join("\n    "),
        );
    }
}

// ---------------------------------------------------------------------------
// DELETED: `the_1420_v2_pair_does_not_converge_by_re_derivation`
// ---------------------------------------------------------------------------
//
// Removed by OPERATOR RULING, 2026-08-13, under
// `rulings[separation-is-a-property-of-the-spelling-not-of-the-variant]`
// (decided). The trace is kept here because a deleted guard leaves none
// otherwise, and because the guard's own closing line asked for exactly this:
// "If a change genuinely licenses the merge, this assertion is the one to argue
// with — cite the clause, then delete it in the same PR that raises the census."
//
// WHAT IT ASSERTED. That `TEMPLATE:g.[37dup;41del]` must not normalize to
// `TEMPLATE:g.[38T>A;40_41delinsTG]`, the string its span sibling
// `g.38_41delinsATTG` prints, on the ground that the cis members "are separated
// by three unchanged nucleotides — 38, 39 and 40".
//
// WHY IT IS OBSOLETE, AND IT IS NOT THAT THE RULE CHANGED.
//
// 1. That ground reads separation off the INPUT'S SPELLING — the members sit at
//    37 and 41, so 38/39/40 are the bases *that spelling* leaves between them.
//
// 2. `separation-is-a-property-of-the-spelling-not-of-the-variant` holds that
//    the separation `general.md:34` keys on is read off the partition
//    RE-DERIVED FROM THE RESULTING SEQUENCE, never off the input's. The guard
//    does the one thing that record forbids.
//
// 3. Read off the re-derived partition, the output is TWO MEMBERS AT SEPARATION
//    ONE, DESCRIBED INDIVIDUALLY — which is what `:34` asks for. **NOTHING
//    MERGED.** The guard's phrase "the merge-across-unchanged-bases veto
//    stopped firing" describes an event that did not occur, and that phrase is
//    what made this guard look load-bearing. There is no `delins` here spanning
//    an unchanged base: `40_41delinsTG` is two consecutive changed columns.
//
// 4. The destination is independently spec-correct, and this file already
//    commits the argument in the sibling assertion: `general.md:56` ranks
//    substitution above duplication, so 38 must not be spelled `dup`;
//    `delins.md:17` keeps the runs individual across the unchanged 39; and
//    Mutalyzer merges but gets no vote under `adjudication-precedence-order`'s
//    spec-explicit > Mutalyzer.
//
// 5. The clause that carried the move is
//    `rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]`. Its scope
//    paragraph warns the deletion "does not license the merges the comparand
//    happened to be blocking" — RESPECTED, because this is not a merge.
//
// Consequences landed in the same change: `CONVERGING_PAIRS_*` raised with the
// clause named, `1420-v2`'s verdict flipped to `Canonical` in
// `reported_partition_verdicts` with `OPEN_GAPS` lowered, and a MEASURED
// addendum on the record above.

/// The three #1421 pairs converge — **without pinning the string they converge
/// on**.
///
/// # Why this is not a string pin
///
/// Measured on this branch, all three land on a form that is neither spelling's
/// `wanted`, and under 5' all three do:
///
/// | pair | 3' | 5' |
/// |---|---|---|
/// | `1421-n1` | `g.[29C>A;32_33delinsACATACTG]` (its `wanted`) | `g.[29delinsAACACAT;32_33delinsTG]` |
/// | `1421-n2` | `g.[32G>T;35C>G;38_39insCGACAT]` | `g.[32G>T;35C>G;36_37insATCGAC]` |
/// | `1421-n3` | `g.[34G>T;37_38delinsCCTTTACG]` (its `wanted`) | `g.[34G>T;35_36insACCTTT;37_38delinsCG]` |
///
/// `reported_partition_verdicts` pins those strings per row, which is what makes
/// the movement disclosed rather than silent. **This test asserts the durable
/// half instead**: that the two spellings of each pair reach ONE string,
/// whatever that string is. A literal pin stops meaning anything the moment the
/// form moves again; this keeps meaning.
///
/// # It rules on nothing
///
/// `rulings[canonical-form-choice-when-both-legal]` (decided) licenses deriving
/// from the resulting sequence and emitting what falls out, so a form neither
/// spelling authored is not wrong on its face. Whether THESE forms are right is
/// a `:47`/`:34` question, and it is deliberately not answered here or by this
/// PR.
///
/// **And this is not the `OPEN_GAPS` hazard.** That hazard is a row reaching
/// *its wanted string* by re-derivation — the shape that banks a defect as a
/// fix. A pair landing on a form that is not `wanted` banks nothing.
#[test]
fn the_1421_pairs_converge_whatever_form_they_land_on() {
    for pair in ["1421-n1", "1421-n2", "1421-n3"] {
        for direction in DIRECTIONS {
            let (label, a, b) = REPORTED_PAIRS
                .iter()
                .find(|(label, _, _)| *label == pair)
                .unwrap_or_else(|| panic!("{pair} is not a reported pair"));
            let (first, second) = (
                normalize_in(TEMPLATE, a, direction),
                normalize_in(TEMPLATE, b, direction),
            );
            assert_eq!(
                first, second,
                "{label} under {direction:?}: the two spellings must reach ONE string. \
                 The string itself is deliberately NOT pinned here — see the doc \
                 above; `reported_partition_verdicts` carries the per-row pins. \
                 If this fails the pair has split apart again, which is a \
                 regression and not a re-bless.",
            );
        }
    }
}

/// Idempotence, which convergence alone does not imply.
///
/// A pass could converge two spellings onto a string that itself normalizes to
/// something else. That would satisfy the census and still be unusable as a
/// stable key, which is the whole point of the exercise.
#[test]
fn every_reported_output_is_a_fixed_point() {
    for (label, a, b) in REPORTED_PAIRS {
        for direction in DIRECTIONS {
            for input in [a, b] {
                let once = normalize_in(TEMPLATE, input, direction);
                let twice = normalize_in(TEMPLATE, &once, direction);
                assert_eq!(
                    once, twice,
                    "{label}: under {direction:?}, {input} normalized to {once}, \
                     which is not a fixed point (it normalizes on to {twice})"
                );
            }
        }
    }
}
