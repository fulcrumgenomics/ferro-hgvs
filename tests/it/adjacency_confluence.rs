//! Confluence: two legal spellings of one variant must reach one normalized
//! form.
//!
//! # Why this is the property that matters downstream
//!
//! Consumers key equality on the normalized string. If two descriptions of the
//! same edit settle apart, every such consumer sees two variants where there is
//! one — and the failure is silent, because both outputs are individually valid
//! and both denote the right sequence. Sequence-preservation checks cannot see
//! it; only comparing the two settled forms can.
//!
//! This has regressed before in exactly this family: #1301 is
//! `g.[263_264insAC;264_265insAA]` settling apart from `g.264_265insCAAA`, and
//! `cis_spelling_confluence_gap` exists because a strict-mode change reopened
//! it. So confluence is asserted directly rather than inferred from the pieces.
//!
//! # Deterministic here, explored in `adjacency_confluence_proptest`
//!
//! The pairs below are the ones whose convergence is *load-bearing* — a spelling
//! a caller plausibly writes against the spelling the pipeline produces. The
//! proptest generates the same property over random geometry; this module is
//! what survives when a generator is later narrowed and stops reaching a case.
//!
//! # One axis this module deliberately does not sweep
//!
//! `FERRO_PARTITION` selects the member-splitting rule and is read **once per
//! process and cached** (see `ferro normalize --help`), so it cannot be varied
//! from inside a test binary. Cross-arm confluence is a CLI-level question and
//! belongs to the bake-off harness, not here. Asserting it in-process would
//! silently test one arm four times.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer, ShuffleDirection};

use crate::common::hg38_window::{local_desc, provider, LOCAL_CONTIG};

fn normalize_in(body: &str, direction: ShuffleDirection) -> Result<String, String> {
    let input = local_desc(body);
    let variant: HgvsVariant = parse_hgvs(&input).map_err(|e| format!("parse: {e}"))?;
    Normalizer::with_config(
        provider(),
        NormalizeConfig::default()
            .with_error_mode(ErrorMode::Strict)
            .with_direction(direction),
    )
    .normalize(&variant)
    .map(|v| v.to_string())
    .map_err(|e| format!("{e}"))
}

fn normalize(body: &str) -> Result<String, String> {
    normalize_in(body, ShuffleDirection::ThreePrime)
}

/// Assert two spellings converge, in **both** shuffle directions.
///
/// Both directions, because a confluence that holds only 3'-ward is a
/// coincidence of the shift rather than a property of the merge — and the 5'
/// direction is a supported configuration, not a curiosity.
///
/// # Each side is unwrapped before the comparison, deliberately
///
/// Comparing the two `Result`s directly would make this helper satisfiable by a
/// **paired refusal**: two `Err`s carrying the same message are equal, so a
/// shape that stopped being accepted in both spellings would converge trivially
/// and every test routed through here would stay green. That is verbatim the
/// defect this change fixes in `issue_1749` ("two spellings that both refuse
/// agree trivially"), and it would be the same defect written a second time.
/// Convergence is a claim about two *accepted* forms, so a refusal on either
/// side is a failure of the property, not evidence for it.
#[track_caller]
fn assert_converges(left: &str, right: &str, why: &str) {
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        let l = normalize_in(left, direction)
            .unwrap_or_else(|e| panic!("{left} was refused in {direction:?}: {e}\n  ({why})"));
        let r = normalize_in(right, direction)
            .unwrap_or_else(|e| panic!("{right} was refused in {direction:?}: {e}\n  ({why})"));
        assert_eq!(
            l, r,
            "\n  {left}\n  and\n  {right}\n  settled apart in {direction:?}: {l:?} vs {r:?}\n  \
             ({why})"
        );
    }
}

// ---------------------------------------------------------------------------
// Spelling confluence: the same edit written different ways.
// ---------------------------------------------------------------------------

/// A tandem copy written as an insertion and as a duplication.
///
/// The pipeline itself re-spells between these two (`insertion.md:17` requires
/// the `dup`), so a caller who writes either must land in the same place.
#[test]
fn an_insertion_and_the_duplication_it_denotes_converge() {
    assert_converges(
        "304_305insTCT",
        "302_304dup",
        "g.302_304 is TCT, so inserting TCT after 304 is a tandem duplication",
    );
}

/// A merged pair and its split spelling.
///
/// `delins.md:86-89` makes the merged form the correct one for an abutting
/// substitution and insertion; a caller who writes the split form must reach it.
#[test]
fn a_split_abutting_pair_and_its_merged_form_converge() {
    assert_converges(
        "[301_302insG;302T>C]",
        "302delinsGC",
        "delins.md:86-89 — the abutting pair IS the delins",
    );
}

/// Two insertions at one junction, written split and written as one payload.
///
/// `general.md:78` gives `ins[A;B]` as the spelling for multiple inserted
/// sequences at one position, so the two-member form and the bracketed payload
/// describe the same thing and must not diverge. This is #1301's shape.
#[test]
fn adjacent_insertions_and_their_combined_payload_converge() {
    assert_converges(
        "[301_302insG;302_303insC]",
        "[301dup;303dup]",
        "#1301 — adjacent-gap insertions must not settle apart from their merged form",
    );
}

/// Deletion-plus-insertion written both ways round, and as a delins.
#[test]
fn a_deletion_with_an_abutting_insertion_converges_from_either_side() {
    assert_converges(
        "[301_302insG;302_304del]",
        "[302_304del;304_305insG]",
        "the junction side is not recoverable from the result",
    );
    assert_converges(
        "[302_304del;304_305insG]",
        "302_304delinsG",
        "both denote replacing 302_304 with G",
    );
}

// ---------------------------------------------------------------------------
// Order confluence: member order is not part of what an allele denotes.
// ---------------------------------------------------------------------------

/// Every permutation of a three-member allele settles identically.
///
/// Three members rather than two because #1103 and #1261 are both order leaks
/// that a two-member test cannot distinguish from a stable sort.
///
/// Every permutation is required to be **accepted** before the forms are
/// compared, for the reason spelled out on [`assert_converges`]: comparing
/// `Result`s would let a uniformly-refused allele satisfy this test, and a
/// uniform refusal is exactly what an order leak looks like once it starts
/// refusing every spelling.
#[test]
fn every_permutation_of_a_three_member_allele_converges() {
    let members = ["301_302insG", "302_306inv", "308del"];
    let permutations = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let mut settled: Option<String> = None;
    for perm in permutations {
        let body = format!(
            "[{};{};{}]",
            members[perm[0]], members[perm[1]], members[perm[2]]
        );
        let got = normalize(&body)
            .unwrap_or_else(|e| panic!("permutation {perm:?} ({body}) was refused: {e}"));
        match &settled {
            None => settled = Some(got),
            Some(first) => assert_eq!(
                &got, first,
                "permutation {perm:?} settled differently: {got:?} vs {first:?}"
            ),
        }
    }
}

// ---------------------------------------------------------------------------
// Idempotence: normalizing a normalized form is a no-op.
// ---------------------------------------------------------------------------

/// Re-normalizing an output reproduces it exactly.
///
/// A form that moves on a second pass is not canonical, and it breaks any
/// consumer that normalizes defensively. `cds_utr3_crossing_shift_idempotency`
/// covers this for the seam; this is the adjacency family's own check.
#[test]
fn normalized_adjacency_forms_are_idempotent() {
    let bodies = [
        "[301_302insG;302T>C]",
        "[302T>C;302_303insG]",
        "[301_302insG;302_304del]",
        "[302_304del;304_305insG]",
        "[301_302insG;302_306inv]",
        "[302_306inv;306_307insG]",
        "[301_302insG;302_304dup]",
        "[301_302insG;302_303insC]",
        "[301_302insG;302_306inv;308del]",
        "304_305insTCT",
    ];
    for body in bodies {
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            let once = normalize_in(body, direction)
                .unwrap_or_else(|e| panic!("{body} was refused in {direction:?}: {e}"));
            let stripped = once
                .strip_prefix(&format!("{LOCAL_CONTIG}:g."))
                .expect("accession survives");
            let twice = normalize_in(stripped, direction)
                .unwrap_or_else(|e| panic!("re-normalizing {once} failed: {e}"));
            assert_eq!(
                once, twice,
                "{body} is not idempotent in {direction:?}: {once} -> {twice}"
            );
        }
    }
}

/// A refusal is a property of the allele, not of how it was written.
///
/// If a conflicting shape is refused in one member order and accepted in
/// another, strict mode's verdict depends on authoring accident — which is the
/// #1508 defect class, where a member became invisible to the detector and the
/// same conflict was accepted or rejected depending on spelling.
#[test]
fn refusals_do_not_depend_on_member_order() {
    for (a, b) in [
        (
            "[302_304delinsGG;302_303insG]",
            "[302_303insG;302_304delinsGG]",
        ),
        ("[302_306inv;303_304insG]", "[303_304insG;302_306inv]"),
        ("[302_304dup;304_305insG]", "[304_305insG;302_304dup]"),
    ] {
        let ra = normalize(a);
        let rb = normalize(b);
        assert_eq!(
            ra.is_err(),
            rb.is_err(),
            "verdict depends on member order: {a} -> {ra:?} but {b} -> {rb:?}"
        );
    }
}
