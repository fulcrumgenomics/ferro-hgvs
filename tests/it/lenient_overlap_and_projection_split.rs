//! Two validity items: a conflict lenient mode drops, and a projection split.
//!
//! Both are descriptions ferro should refuse and does not, and both are rank-1
//! validity failures — "the output re-parses, denotes a sequence, and violates no
//! absolute prohibition" is rank 1 of the operator's precedence order
//! (`rulings[adjudication-precedence-order]`, tabulated at
//! `tests/it/spec_conformance_axis.rs:20`). They are kept in one module because
//! the *shape* of each defect is the same: a stage detects something and the
//! finding does not survive to the caller.
//!
//! # Item 1 — the conflict is detected and then dropped at the API seam
//!
//! An allele whose members claim the same reference territory denotes no
//! sequence. Ferro knows this: `detect_overlap_conflicts` fires on the raw
//! members, `ErrorMode::Strict` promotes the W5002 warning to an error, and
//! lenient mode hands the description back exactly as authored (#395/#486/#1004,
//! and #1406's gate at `src/normalize/mod.rs:3450`). All of that is correct and
//! is pinned below as such.
//!
//! What is **not** correct is that the finding never reaches the caller.
//! `Normalizer::normalize` returns `Result<HgvsVariant, _>` — one value, no
//! warning channel (`src/normalize/mod.rs:1764`) — so a lenient caller receives
//! `Ok(<a description denoting no sequence>)` and nothing else. Measured through
//! the CLI, `--error-mode lenient` and `--error-mode silent` are byte-identical:
//! the description on stdout, empty stderr, exit 0. Yet `ErrorMode` itself says
//! they differ — `emits_warnings()` is `true` for `Lenient` and `false` for
//! `Silent` (`src/error_handling/types.rs:46-48`), and `Lenient`'s own doc says
//! "warnings will be generated" (`:20-24`).
//!
//! So the laundering is not in the normalizer's algorithm. It is at the seam:
//! **the mode that promises a diagnostic has no channel to deliver one through
//! the entry point that returns the output.**
//!
//! ## The tension, stated rather than glossed
//!
//! Validity is rank 1 and is never traded. Lenient mode exists precisely to keep
//! going on bad input. Refusing honours the first and abolishes the second;
//! silently accepting honours the second and abolishes the first.
//!
//! **The adjudication: pass the input through, and report.** The precedence order
//! ranks candidate outputs a normalizer *chose*; preservation is a decline to
//! choose, not a choice. Handing back an invalid description is not ferro
//! asserting the description is valid — it becomes exactly that once the
//! invalidity goes unreported, because then no caller can tell the two apart.
//! The diagnostic is therefore not a nicety attached to the pass-through; it is
//! what makes the pass-through a decline rather than an endorsement. Refusal
//! stays strict mode's job, and strict mode already does it.
//!
//! # Item 2 — a single-member description projects into several members
//!
//! See the section header further down. Same shape: one variant in, a bracketed
//! allele out, on an axis whose input had no allele to split.

use ferro_hgvs::conformance::spec_corpus::{
    corpus_cores, denotation_of, Denotation, Frame, RefShape, DENSE_CORE_LEN,
};
use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

/// The `AT`-alphabet core the corpus draws first, as `spec_corpus_regressions`
/// uses it. Long homopolymer runs, which is where an ambiguous shift lives.
fn at_core() -> String {
    corpus_cores(1, DENSE_CORE_LEN).remove(0)
}

/// The three conflicting geometries `spec_corpus_regressions::
/// conflicting_member_geometries_are_normalized_instead_of_refused` names, each
/// of which the corpus verifies denotes no sequence.
///
/// Repeated here rather than imported so this module states its own premise: a
/// row that quietly stopped being a conflict would otherwise satisfy every
/// assertion below while pinning nothing.
const CONFLICTING_ALLELES: [(&str, &str); 3] = [
    ("nested", "NM_TEST.1:c.[9_12del;10_11del]"),
    ("overlapping", "NM_TEST.1:c.[9_12del;11_14del]"),
    ("coincident insertions", "NM_TEST.1:c.[8_9insAC;8_9insCC]"),
];

/// A normalizer over `frame` in the given error mode.
fn normalizer(frame: &Frame, mode: ErrorMode) -> Normalizer<MockProvider> {
    Normalizer::with_config(
        frame.provider().clone(),
        NormalizeConfig::default().with_error_mode(mode),
    )
}

/// Normalize through the plain entry point — the one that returns only a value.
fn normalize_in(frame: &Frame, mode: ErrorMode, input: &str) -> Result<String, String> {
    normalizer(frame, mode)
        .normalize(&parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}")))
        .map(|v| v.to_string())
        .map_err(|e| e.to_string())
}

/// The W5002 diagnostics `input` produces in `mode`, via the entry point that
/// has a channel for them.
fn conflict_diagnostics(frame: &Frame, mode: ErrorMode, input: &str) -> usize {
    normalizer(frame, mode)
        .normalize_with_diagnostics(
            &parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}")),
        )
        .expect("a non-strict mode must not reject")
        .warnings
        .iter()
        .filter(|w| {
            matches!(
                w,
                ferro_hgvs::normalize::NormalizationWarning::OverlapConflict { .. }
            )
        })
        .count()
}

// ---------------------------------------------------------------------------
// Item 1: what lenient mode does with an allele that denotes no sequence
// ---------------------------------------------------------------------------

/// **Question.** What does lenient mode emit for an allele whose members claim
/// one another's territory, and does the output denote anything?
///
/// **It emits the input, verbatim, and the output denotes nothing** — the same
/// nothing the input denoted. Measured on all three geometries the corpus builds:
/// input and output are byte-identical and both classify `Denotation::NoSequence`.
///
/// The relevant clause is `general.md:58` — "descriptions removing part of a
/// reference sequence and replacing it with part of the same sequence are not
/// allowed (e.g., `NM_004006.2:c.[762_768del;767_774dup]`)" — and it is worth
/// being exact about how far it reaches: it names the `del`+`dup` spelling, which
/// ferro already rejects at parse. It does not name `del`+`del` or two insertions
/// at one interbase. The plainer statement, the one that does reach all three, is
/// that a description must denote one sequence and these denote none.
///
/// **ADJUDICATED CORRECT, in part.** Emitting the input rather than a repaired
/// form is ferro's standing answer since #395/#486/#1004 and is right: the
/// independent applier declines an overlapping description outright, so there are
/// no "input's bases" a repair could preserve — any output would come from
/// *choosing* an apply order, which is what W5002 refuses to do. The half that is
/// not right is pinned in the next test.
#[test]
fn lenient_hands_back_a_conflicting_allele_that_denotes_no_sequence() {
    let core = at_core();
    let frame = Frame::build(RefShape::CodingSingleExon, &core);
    let served = frame.served();

    for (geometry, input) in CONFLICTING_ALLELES {
        // The premise. Without it this test asserts a preference, not a defect.
        assert_eq!(
            denotation_of(frame.provider(), served, input),
            Denotation::NoSequence,
            "{geometry}: `{input}` must denote no sequence for this test to mean anything"
        );

        let output = normalize_in(&frame, ErrorMode::Lenient, input)
            .unwrap_or_else(|e| panic!("{geometry}: lenient must not reject `{input}`: {e}"));
        assert_eq!(
            output, input,
            "{geometry}: the conflicting allele must come back as authored"
        );
        assert_eq!(
            denotation_of(frame.provider(), served, &output),
            Denotation::NoSequence,
            "{geometry}: and the output denotes no sequence — a rank-1 validity failure \
             that lenient mode returns as `Ok`"
        );
        // Stable in one pass, so preserving is an answer rather than a pause.
        assert_eq!(
            normalize_in(&frame, ErrorMode::Lenient, &output).expect("lenient"),
            output,
            "{geometry}: the preserved form must be a fixed point"
        );
    }
}

/// **Question.** Is the conflict reported to a lenient caller?
///
/// **Only to one that asks a different question.** The W5002 `OverlapConflict`
/// warning is produced — `normalize_with_diagnostics` carries exactly one per
/// input — but `Normalizer::normalize`, the entry point the CLI and every
/// in-repo caller use, returns `Result<HgvsVariant, _>` and has nowhere to put
/// it (`src/normalize/mod.rs:1764`).
///
/// **PINNED DEFECT.** The consequence, asserted rather than described: through
/// `normalize`, `ErrorMode::Lenient` and `ErrorMode::Silent` are
/// indistinguishable — same `Ok`, same string — while `ErrorMode` itself
/// declares they differ (`emits_warnings()`, `src/error_handling/types.rs:46-48`)
/// and `Lenient`'s doc promises "warnings will be generated" (`:20-24`).
/// Confirmed end-to-end at the CLI: `ferro normalize --error-mode lenient` and
/// `--error-mode silent` both print the description on stdout with empty stderr
/// and exit 0.
///
/// **Correct behaviour:** a lenient normalization that declined to canonicalize
/// must say so through the same call that returned the description — the input
/// still passes through unchanged, refusal stays strict mode's job. Until then a
/// caller cannot distinguish "ferro normalized this" from "ferro refused to, and
/// this denotes no sequence".
#[test]
fn the_conflict_is_dropped_by_the_entry_point_that_returns_the_output() {
    let core = at_core();
    let frame = Frame::build(RefShape::CodingSingleExon, &core);

    for (geometry, input) in CONFLICTING_ALLELES {
        // The finding exists. Asserted first, so the pin below is about the
        // channel rather than about detection having quietly stopped.
        assert_eq!(
            conflict_diagnostics(&frame, ErrorMode::Lenient, input),
            1,
            "{geometry}: the conflict must be detected for this test to be about reporting"
        );
        // And strict mode reaches it, which is what makes the drop a seam
        // problem rather than a detection one.
        let strict = normalize_in(&frame, ErrorMode::Strict, input)
            .expect_err("strict mode must reject a conflicting allele");
        assert!(
            strict.contains("W5002") || strict.contains("OverlapConflictingEdits"),
            "{geometry}: strict must reject `{input}` AS a conflict; got: {strict}"
        );

        // The defect. `ErrorMode` says these two modes differ...
        assert!(
            ErrorMode::Lenient.emits_warnings() && !ErrorMode::Silent.emits_warnings(),
            "the modes must claim to differ for the pin below to be a contradiction"
        );
        // ...and through `normalize` they do not.
        assert_eq!(
            normalize_in(&frame, ErrorMode::Lenient, input),
            normalize_in(&frame, ErrorMode::Silent, input),
            "PINNED DEFECT ({geometry}) — lenient and silent are indistinguishable through \
             `Normalizer::normalize`, which has no warning channel. Correct behaviour: a \
             lenient caller must learn from this call that `{input}` was declined and denotes \
             no sequence."
        );
    }
}

/// **Question.** Does the pass-through at least keep the conflict *visible in the
/// description*, so that re-reading the output reaches the same verdict?
///
/// **Yes, and that is the property to protect.** #1406 row 3 was the case where
/// it did not: lenient mode repaired the members one at a time until the output
/// no longer looked conflicting, so strict mode accepted a description it had
/// just rejected. The gate at `src/normalize/mod.rs:3450-3481` closed that, and
/// `issue_1406_lenient_output_keeps_the_conflict.rs` pins it for the
/// `[inv;ins]` family.
///
/// **ADJUDICATED CORRECT**, asserted here for the three *corpus* geometries,
/// which that file does not cover. It is the reason refusal is not needed in
/// lenient mode: the invalidity is still recoverable from the output by anyone
/// who re-reads it strictly. What it does not fix is that the caller who got the
/// output was never told to.
#[test]
fn the_lenient_output_is_still_refused_by_strict() {
    let core = at_core();
    let frame = Frame::build(RefShape::CodingSingleExon, &core);

    for (geometry, input) in CONFLICTING_ALLELES {
        let output = normalize_in(&frame, ErrorMode::Lenient, input).expect("lenient");
        let strict = normalize_in(&frame, ErrorMode::Strict, &output)
            .expect_err("the lenient output must still be a conflict");
        assert!(
            strict.contains("W5002") || strict.contains("OverlapConflictingEdits"),
            "{geometry}: lenient turned `{input}` into `{output}`, which strict does not \
             reject as a conflict — the conflict was laundered out of the description"
        );
    }
}
