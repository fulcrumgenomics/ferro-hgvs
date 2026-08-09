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
//! (**spec-explicit > Mutalyzer > our judgement**). Twelve of the eighteen rows
//! have a spec line that answers them outright; #1419's six do not, because
//! `delins.md:44-47` recommends the merged delins for the very alignment
//! coincidence `general.md:56` is being used to split — so those stay recorded
//! rather than targeted.
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
/// constant on its own and no sum is ever taken. Measured: both are currently 0.
const CONVERGING_PAIRS_THREE_PRIME: usize = 0;

/// How many reported pairs converge today under `ShuffleDirection::FivePrime`.
///
/// See [`CONVERGING_PAIRS_THREE_PRIME`] for why this is tracked separately
/// rather than shared. Measured: currently 0, same as 3'.
const CONVERGING_PAIRS_FIVE_PRIME: usize = 0;

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
             (whichever direction moved) and say so in the PR — it is a representation \
             change for anyone storing the losing spelling. If it went DOWN, a change has \
             regressed an externally-reported defect.",
            converged.len(),
            converged.join("\n    "),
            split.len(),
            split.join("\n    "),
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
