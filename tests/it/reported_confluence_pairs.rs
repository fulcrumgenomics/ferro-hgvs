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
//! (**spec-explicit > Mutalyzer > our judgement**). Twelve of those eighteen
//! spellings — six pairs — have a spec line that answers them outright, which is
//! the split `the_spec_authority_census_holds` pins next door as `(12, 6)` rows:
//! #1420's six on `delins.md:16`/`:17`, and #1421's six on `:17`, whose only
//! competitor `delins.md:47` is scoped away from the `g.` axis and away from net
//! insertions by two decided ledger records. #1419's six do not, because
//! `delins.md:46-47` recommends the merged delins for the alignment coincidence
//! their split relies on, and their `wanted` form is that recommendation's own
//! answer — so those stay recorded rather than targeted.
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
/// constant on its own and no sum is ever taken. Measured: both are currently 3.
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
/// [`the_1420_v2_pair_does_not_converge_by_re_derivation`] below names the
/// string so it cannot happen quietly here.
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
/// **3 -> 9 under the partition default flip: ALL NINE reported pairs now
/// converge, and the clause that carried it is
/// `rulings[canonical-form-choice-when-both-legal]`** (decided) — ferro derives
/// the description from the resulting sequence and emits what falls out, so two
/// spellings of one variant derive the same pieces and reach one answer. The
/// six that joined are `1420-v2`, `1420-v3`, `1420-v4`, `1421-n1`, `1421-n2` and
/// `1421-n3`.
///
/// **Read what it did NOT buy, because the test's own warning is about exactly
/// this.** A rise "onto a form neither spelling asserted ... is not progress
/// even when it lands on the form the issue asks for" — and none of the six
/// lands on the form its issue asks for. Every one of them converges while
/// moving OFF `wanted`, which is why all six pairs also move to
/// `PairState::NeitherReaches` in `reported_partition_verdicts.rs` and
/// `OPEN_GAPS` rises 12 -> 18 in the same change. Convergence and conformance
/// moved in opposite directions here; the census going up is not a fix.
///
/// Two of the six have a further, separate finding: `1420-v2` and `1420-v3`
/// converge on a form `general.md:56` ranks BELOW the one they left (a `dup`
/// and an `ins` where a substitution was available). That is a README rule-2
/// preference miss, filed as issue #1878. `1420-v4` is licensed instead — the
/// `g.` axis is outside `delins.md:47` by
/// `rulings[delins-payload-coincidence-carve-out-is-coding-dna-scoped]`.
const CONVERGING_PAIRS_THREE_PRIME: usize = 9;

/// How many reported pairs converge today under `ShuffleDirection::FivePrime`.
///
/// See [`CONVERGING_PAIRS_THREE_PRIME`] for why this is tracked separately
/// rather than shared, and for why a rise in either direction has to name the
/// clause that carried it. Measured: currently 3, same as 3' — the same three
/// `1419-r*` pairs, raised from 0 by #1649.
///
/// The two directions agreeing here is worth reading rather than skipping: the
/// pairs converge under 5' on a *different* string in `1419-r3`'s case, because
/// the 5' pass shifts that pair's leading deletion one base left. Convergence is
/// a property of the partition the splitter can express, and that is what #1649
/// changed; which end of the run the members then shuffle to is downstream of
/// it and does not decide whether the two spellings meet.
/// **3 -> 9 under the partition default flip, the same six pairs and the same
/// clause as 3' — see [`CONVERGING_PAIRS_THREE_PRIME`], including why this
/// rise is not by itself progress.** The two directions agree on the count
/// but not always on the string: `1421-n3` converges on
/// `g.[34G>T;37delinsCC;39_40insTACGT]` at 3' and
/// `g.[34G>T;37delinsCC;38_39insTTACG]` at 5', which is why both its rows
/// joined `FIVE_PRIME_MOVERS` next door.
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

/// `1420-v2` must not converge by re-derivation, in either direction.
///
/// The one row in this family where a defect coincides with an issue's ask, and
/// therefore the one that can be banked as a fix without anybody noticing.
///
/// `g.[37dup;41del]`'s members are separated by three unchanged nucleotides —
/// 38, 39 and 40 — so `general.md:34` ("two variants separated by one or more
/// nucleotides should be described individually and **not** as a \"delins\"")
/// keeps them individual. Re-partitioning the block from the sequence instead
/// yields `g.[38T>A;40_41delinsTG]`, which is both what the *span* spelling
/// `g.38_41delinsATTG` prints and what #1420 asks for. So the census would show
/// this pair converging, on the issue's own wanted string, purely because the
/// merge-across-unchanged-bases veto stopped firing.
///
/// Observed exactly that way on an experimental partitioner arm, where the 5'
/// count went 0 -> 1 on this row alone. Named here so the string, not the
/// count, is what fails.
///
/// If a change genuinely licenses the merge, this assertion is the one to argue
/// with — cite the clause, then delete it in the same PR that raises the census.
#[test]
fn the_1420_v2_pair_does_not_converge_by_re_derivation() {
    // RE-EXPRESSED FOR THE PARTITION DEFAULT FLIP, AND THE HAZARD IT NAMES DID
    // NOT OCCUR. This guard was written against one mechanism: the CIS spelling
    // merging across the unchanged 38-40 to reach its span sibling's answer,
    // which `general.md:34` forbids and which would have coincided with #1420's
    // wanted form — a defect that reads as a fix.
    //
    // Under `canonical-coalesced` the pair does converge, but by the OPPOSITE
    // route and without that merge. The cis is returned UNCHANGED at
    // `g.[37dup;41del]` — its members still individual across the unchanged
    // nucleotides, so `:34` is satisfied — and it is the SPAN that re-derives
    // down onto it, off `g.[38T>A;40_41delinsTG]`.
    //
    // So the original `assert_ne!` would now fail on a convergence that is not
    // the thing it was protecting against, while the real invariant — no merge
    // across the unchanged interior — is untested. Both are fixed here: the cis
    // is pinned positively, and the merged spelling is named as forbidden
    // directly rather than inferred from equality with the span.
    //
    // What the convergence DOES cost is that neither spelling now reaches
    // #1420's wanted form, because the surviving form spells the change at 38 as
    // a `dup` where `general.md:56` ranks a substitution above it. That is a
    // README rule-2 preference miss, disclosed in the PR trailer, pinned in
    // `reported_partition_verdicts.rs` (`1420-v2/span` -> `Gap`, the pair ->
    // `NeitherReaches`, `OPEN_GAPS` 12 -> 13) and filed as issue #1878. Closing
    // #1878 should move the span back onto `MERGED_ACROSS_UNCHANGED` below —
    // which is why that string stays named here rather than being deleted.
    const SPAN: &str = "TEMPLATE:g.38_41delinsATTG";
    const CIS: &str = "TEMPLATE:g.[37dup;41del]";
    const CIS_OUTPUT: &str = "TEMPLATE:g.[37dup;41del]";
    const CIS_OUTPUT_FIVE_PRIME: &str = "TEMPLATE:g.[36dup;41del]";
    // The form #1420 asks for, and the form a merge across the unchanged 39
    // would produce. Not currently reached by either spelling; see #1878.
    const MERGED_ACROSS_UNCHANGED: &str = "TEMPLATE:g.[38T>A;40_41delinsTG]";

    for direction in DIRECTIONS {
        let cis = normalize_in(TEMPLATE, CIS, direction);
        let expected_cis = match direction {
            ShuffleDirection::FivePrime => CIS_OUTPUT_FIVE_PRIME,
            _ => CIS_OUTPUT,
        };
        assert_eq!(
            cis, expected_cis,
            "{direction:?}: the `1420-v2` cis spelling moved. Its members are \
             three unchanged nucleotides apart (38, 39, 40), so `general.md:34` \
             keeps them individual — a cis output that is NOT this string is \
             the thing this guard exists to catch, whatever it converged with",
        );
        assert_ne!(
            cis, MERGED_ACROSS_UNCHANGED,
            "{direction:?}: the `1420-v2` cis spelling merged across three \
             unchanged nucleotides (38, 39, 40), which `general.md:34` forbids. \
             It coincides with #1420's wanted form; that is not a licence to \
             raise the census",
        );
        let span = normalize_in(TEMPLATE, SPAN, direction);
        assert_ne!(
            span, MERGED_ACROSS_UNCHANGED,
            "{direction:?}: the `1420-v2` span spelling now reaches #1420's \
             wanted form again. If that is issue #1878 being fixed, this is the \
             good news — re-pin `1420-v2/span` to `Canonical`, the pair to \
             `OneReaches`, and `OPEN_GAPS` back to 12, and delete this \
             assertion in the same change",
        );
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
