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
//! Exercised on the `c.` axis. Reference-free (`MockProvider` via
//! `SyntheticBuilder`), so this holds independent of the manifest.
//!
//! **Axis rationale corrected by #1284.** This header used to say the `c.` axis
//! was chosen because "the transcript axes are not reached by that pass, so the
//! dup survives". That was true and is no longer: #1284 lifted the collision
//! repair's genomic-only gate, so a normalization-produced dup that collides
//! with a sibling is now respelled on `c.` exactly as on `g.`. The tests below
//! that author a `dup` **directly** are unaffected — they hand the detector a
//! dup rather than relying on normalization to make one — and they are what
//! carries the #1243 coverage.

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

/// A duplication and an insertion at one junction are **accepted**.
///
/// This test asserted the opposite until the dup-plus-insertion ruling, and its
/// old rationale — "ambiguous in exactly the way two insertions at one junction
/// are: nothing orders them" — is the conflation that ruling removed. Two
/// insertions at one junction genuinely have nothing to order them, because both
/// write at the same interbase point. A duplication and an insertion do not
/// compete: the insertion consumes no reference base, so only one member claims
/// a base and the pair has a defined resulting sequence.
///
/// Ruled in `delins-adjacent-members-when-both-consume-reference`, whose
/// duplication half decides the *accept* question (yes — `duplication.md:90`
/// publishes this geometry, `:91` glosses the member order).
///
/// The **form** question that ruling leaves open is deliberately not blessed
/// here. The settled string below is recorded as measured behaviour, not as a
/// verdict: the pair currently coalesces into a single member, and whether that
/// or the split spelling is canonical cannot be decided until the choice is made
/// after the member set is final (#1946). If a change moves this string, that is
/// a fact to weigh against that open question — not automatically a regression.
#[test]
fn strict_accepts_a_duplication_and_an_insertion_at_one_junction() {
    const INPUT: &str = "NM_TEST.1:c.[5_6dup;6_7insA]";
    assert!(
        !strict_rejects_as_overlap(INPUT),
        "the pair writes disjoint content and must be accepted"
    );
    // Measured, not blessed. The duplication is absorbed into the insertion:
    // its two copied bases and the inserted `A` become one three-base payload
    // at the shared junction, and the `dup` label does not survive.
    assert_eq!(
        normalize_lenient(INPUT),
        "NM_TEST.1:c.6_7insAAA",
        "records the coalesced form; the form question is open, see the doc comment"
    );
}

/// The regression proper: lenient normalization must not emit a form that its
/// own strict mode accepts.
///
/// `[4_10inv;5_6insAA]` is an insertion strictly interior to an inversion —
/// rejected under strict since #486. Normalizing it leniently rewrote that
/// insertion into a `dup`, and before #1243 the result passed strict
/// validation. A caller who normalized leniently and re-validated strictly got a
/// clean answer for an allele with no defined resulting sequence.
///
/// **Settled form updated by #1284**, which is why the pinned string below is
/// an `ins` and not the `dup` this test was written around: the dup the
/// pipeline used to leave here now collides with the inversion's bases and is
/// repaired into the equivalent insertion, as it already was on `g.`. The
/// contract this test exists for is untouched and is asserted twice — strict
/// rejects the input, and strict still rejects whatever lenient mode emits for
/// it. The dup-specific half of #1243 is carried by the direct-authoring tests
/// above, which do not depend on normalization producing a dup.
#[test]
fn lenient_output_of_a_conflict_is_still_rejected_by_strict() {
    const INPUT: &str = "NM_TEST.1:c.[4_10inv;5_6insAA]";
    assert!(
        strict_rejects_as_overlap(INPUT),
        "the #486 contract: strict rejects an insertion interior to a span"
    );

    let normalized = normalize_lenient(INPUT);
    assert_eq!(
        normalized, "NM_TEST.1:c.[5_9inv;6_7insAA]",
        "precondition: normalization settles the interior insertion at its own \
         junction, still strictly inside the inversion"
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
///
/// The input moved when the dup-plus-insertion ruling made the original
/// `[5_6dup;6_7insA]` a **non**-conflict — it now coalesces and emits no
/// `OverlapConflict` at all, so it can no longer witness a labelling contract
/// about conflicting members. #1243's contract is untouched and is asserted
/// here on a shape that still conflicts: a duplication interior to an
/// inversion, which the sibling test above pins as rejected under strict.
///
/// Keeping the dup in the input is the point. A hardcoded `"ins"` is exactly
/// what this test exists to catch, and only a conflicting group *containing a
/// duplication* can catch it.
#[test]
fn the_warning_names_the_duplication_as_a_dup() {
    let variant = parse_hgvs("NM_TEST.1:c.[4_10inv;5_6dup]").expect("parse");
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
    // were mislabelled.
    assert!(
        overlap_kinds.iter().any(|kinds| {
            let mut sorted: Vec<&str> = kinds.iter().map(String::as_str).collect();
            sorted.sort_unstable();
            sorted == ["dup", "inv"]
        }),
        "the overlap warning must name its members by their real kinds; got {overlap_kinds:?}"
    );
}
