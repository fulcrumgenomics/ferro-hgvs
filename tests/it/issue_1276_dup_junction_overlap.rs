//! Issue #1276 — a duplication occupies an insertion junction, so it must be
//! subject to the same overlap detection an insertion is.
//!
//! The fix itself landed with #1243, which registers a `dup` (and a `repeat`)
//! as a junction occupant in `junction_writing_kind`. This file is the
//! end-to-end lock for it: #1243 covers the detector directly, in
//! `overlap.rs`'s own unit tests, and asserts on warnings; what is pinned here
//! is the *caller-visible* contract those warnings exist to produce — what
//! strict mode rejects, and what lenient mode settles on. Six of the seven
//! below were confirmed to fail with the `Duplication` arm of
//! `junction_writing_kind` deleted; the seventh is the negative case and
//! correctly does not move.
//!
//! `x_ydup` places its copy at the junction after `y`. `overlap.rs` registered
//! it only as a *span*, never as a junction occupant — and normalization is
//! precisely what rewrites an interior insertion **into** a duplication. So the
//! conflict `detect_insertion_overlaps` reported on the input went undetected on
//! its own output, with two consequences pinned below:
//!
//!   1. lenient normalize turned a strict-**rejected** allele into one strict
//!      mode **accepted**; and
//!   2. `normalize` was not idempotent, because
//!      `sort_cis_members_by_genomic_order` is skipped only for conflicting
//!      alleles — so the member order flipped on the second pass.
//!
//! Exercised on the `c.` axis. The `g.` spelling of the same shape is currently
//! collapsed by the sequence-first canonicalization before the dup is ever
//! rendered, which would mask the behaviour under test; the transcript axes are
//! not reached by that pass, so the dup survives and the detector is what
//! decides. Reference-free (`MockProvider` via `SyntheticBuilder`), so this
//! holds independent of the manifest.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::error::FerroError;
use ferro_hgvs::normalize::NormalizationWarning;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

/// `c.4..10` is `CAATTGG` — no tandem repeat, so a `dup` here has no room to
/// 3'-shift and these assertions are about overlap detection, not shuffling.
const CORE: &str = "GGGCAATTGGGCCCAAATTTGGG";

fn provider() -> MockProvider {
    SyntheticBuilder::cds(CORE, 1, 21, Strand::Plus).build()
}

fn normalize_lenient(input: &str) -> String {
    let variant = parse_hgvs(input).expect("parse");
    Normalizer::new(provider())
        .normalize(&variant)
        .expect("lenient mode must not reject")
        .to_string()
}

/// Whether strict mode rejects `input` specifically as `W5002`. Panics on any
/// other error variant — and on an `InvalidCoordinates` that is not the overlap
/// rejection — so a rejection for an unrelated reason cannot be mistaken for the
/// contract under test.
///
/// The promotion site emits both halves of the tag in one message
/// (`… (OverlapConflictingEdits / W5002)`), so both are required: accepting
/// either alone would let the message keep passing after the code half of the
/// contract regressed.
fn strict_rejects_as_overlap(input: &str) -> bool {
    let variant = parse_hgvs(input).expect("parse");
    match Normalizer::with_config(provider(), NormalizeConfig::strict()).normalize(&variant) {
        Ok(_) => false,
        Err(FerroError::InvalidCoordinates { msg }) => {
            assert!(
                msg.contains("W5002") && msg.contains("OverlapConflictingEdits"),
                "expected `OverlapConflictingEdits / W5002`; got: {msg}"
            );
            true
        }
        Err(other) => panic!("unexpected error variant: {other:?}"),
    }
}

/// A duplication whose copy lands strictly inside an inversion has no defined
/// resulting sequence — invert-then-duplicate and duplicate-then-invert differ.
///
/// Accepted before #1243.
#[test]
fn strict_rejects_a_duplication_interior_to_an_inversion() {
    assert!(strict_rejects_as_overlap("NM_TEST.1:c.[4_10inv;5_6dup]"));
}

/// Same for a delins: the copy would land in territory the sibling replaces.
///
/// (There is no `del` counterpart — the parser rejects an overlapping
/// `del`+`dup` earlier still, as a self-cancelling allele.)
#[test]
fn strict_rejects_a_duplication_interior_to_a_delins() {
    assert!(strict_rejects_as_overlap(
        "NM_TEST.1:c.[4_10delinsTT;5_6dup]"
    ));
}

/// A duplication and an insertion competing for one junction are ambiguous in
/// exactly the way two insertions at one junction are: nothing orders them.
#[test]
fn strict_rejects_a_duplication_and_an_insertion_at_one_junction() {
    assert!(strict_rejects_as_overlap("NM_TEST.1:c.[5_6dup;6_7insA]"));
}

/// The regression proper: lenient normalization must not emit a form that its
/// own strict mode accepts.
///
/// `[4_10inv;5_6insAA]` is an insertion strictly interior to an inversion —
/// rejected under strict since #486. Normalizing it leniently rewrites that
/// insertion into a `dup`, and before #1243 the result passed strict
/// validation. A caller who normalized leniently and re-validated strictly got a
/// clean answer for an allele with no defined resulting sequence.
#[test]
fn lenient_output_of_a_conflict_is_still_rejected_by_strict() {
    const INPUT: &str = "NM_TEST.1:c.[4_10inv;5_6insAA]";
    assert!(
        strict_rejects_as_overlap(INPUT),
        "the #486 contract: strict rejects an insertion interior to a span"
    );

    let normalized = normalize_lenient(INPUT);
    assert_eq!(
        normalized, "NM_TEST.1:c.[5_9inv;5_6dup]",
        "precondition: normalization rewrites the interior insertion as a dup, \
         which is what used to escape detection"
    );
    assert!(
        strict_rejects_as_overlap(&normalized),
        "lenient output `{normalized}` must not pass the strict check its own \
         input failed"
    );
}

/// The second symptom. `sort_cis_members_by_genomic_order` is deliberately
/// skipped for overlap-conflicting alleles, so a conflict that vanished on the
/// second pass let the member order flip: `[4_10inv;5_6insAA]` settled as
/// `[5_9inv;5_6dup]` on pass one and `[5_6dup;5_9inv]` on pass two.
#[test]
fn normalizing_a_dup_conflict_is_idempotent() {
    for input in [
        "NM_TEST.1:c.[4_10inv;5_6insAA]",
        "NM_TEST.1:c.[4_10inv;5_6dup]",
        "NM_TEST.1:c.[5_6dup;4_10inv]",
        "NM_TEST.1:c.[5_6dup;6_7insA]",
    ] {
        let once = normalize_lenient(input);
        let twice = normalize_lenient(&once);
        assert_eq!(
            once, twice,
            "`{input}` must settle in one pass; got `{once}` then `{twice}`"
        );
    }
}

/// The interior bound is strict at the 3' end (`gap < end`), so a duplication
/// never flags *itself* and the shapes that must stay valid still do.
#[test]
fn non_overlapping_duplications_are_still_accepted() {
    for input in [
        // A lone dup has no sibling to conflict with — the self-flag check.
        "NM_TEST.1:c.4_6dup",
        // Siblings that do not touch the dup's junction at c.6.
        "NM_TEST.1:c.[4_6dup;8T>A]",
        "NM_TEST.1:c.[4_6dup;7_8insG]",
        "NM_TEST.1:c.[4_6dup;2_3insG]",
        // Two disjoint duplications.
        "NM_TEST.1:c.[4_6dup;10_12dup]",
    ] {
        assert!(
            !strict_rejects_as_overlap(input),
            "`{input}` has no overlap and must stay accepted"
        );
    }
}

/// The warning names each conflicting member by its real edit kind. Before
/// #1243 a same-junction group was assumed to be all insertions and hardcoded
/// `"ins"`, which would now mislabel the duplication.
#[test]
fn the_warning_names_the_duplication_as_a_dup() {
    let variant = parse_hgvs("NM_TEST.1:c.[5_6dup;6_7insA]").expect("parse");
    let result = Normalizer::new(provider())
        .normalize_with_diagnostics(&variant)
        .expect("lenient mode must not reject");

    // Select the overlap warning *first*, then assert on its own member labels:
    // matching `"dup"` anywhere in the rendered warning list would also be
    // satisfied by an unrelated warning that merely mentions a duplication.
    let overlap_kinds: Vec<&Vec<String>> = result
        .warnings
        .iter()
        .filter_map(|w| match w {
            NormalizationWarning::OverlapConflict { edit_kinds, .. } => Some(edit_kinds),
            _ => None,
        })
        .collect();
    assert!(
        !overlap_kinds.is_empty(),
        "expected an OverlapConflict warning; got {:?}",
        result.warnings
    );
    // The contract is that *each* member is named by its real edit kind, so
    // pin both labels: asserting only `dup` would still pass if the sibling
    // insertion were mislabelled.
    assert!(
        overlap_kinds.iter().any(|kinds| {
            let mut sorted: Vec<&str> = kinds.iter().map(String::as_str).collect();
            sorted.sort_unstable();
            sorted == ["dup", "ins"]
        }),
        "the overlap warning must name its members as `dup` and `ins`; got {overlap_kinds:?}"
    );
}
