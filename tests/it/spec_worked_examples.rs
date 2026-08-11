//! The HGVS recommendations' own worked normalization examples, run against
//! real reference bases (GATE 5).
//!
//! # Why this file exists
//!
//! Most conformance rows are somebody's *reading* of the spec. These twelve are
//! not: the spec names an input and prints, in the same sentence, the
//! description to use instead — "use `NM_024312.4:c.2692_2693dup` and **not**
//! `NM_024312.4:c.2686A[10]`". There is nothing to interpret, so a divergence
//! is a defect and a match is worth pinning.
//!
//! They had nonetheless never been executed. Both rules they cover — the
//! coding-sequence restriction on repeat notation and the 3'rule's exon-junction
//! exception — are sequence-dependent, and every existing home for them is
//! sequence-blind:
//!
//! - the harvested spec fixture (`hgvs_spec_normalization_tests.rs`) normalizes
//!   through `MockProvider::new()`, the empty constructor, so no rule that reads
//!   bases can fire. It carries `LRG_199t1:c.3921del`, but that is a *different
//!   edit* against an empty provider — a parser/renderer check, not a
//!   normalization one;
//! - every family in `examples/dump_normalized_corpus.rs` is single-exon, so no
//!   generated row can exhibit an exon junction at all (#1478);
//! - the manifest-backed axes skip when the prepared reference is absent, which
//!   in CI is always.
//!
//! So the flagship worked example of `general.md:44` had never run. This module
//! is its blocking home: a `WindowProvider` over the committed slice beside
//! `cases.json` serves the real bases with no `FERRO_MANIFEST`, so the gate
//! cannot silently skip.
//!
//! # Why three modes per row
//!
//! Pinning one output would have hidden both of the interesting results. The
//! shipped `--error-mode` modes disagree on these inputs, and *which* mode
//! produces which answer is the finding:
//!
//! - `NM_024312.4:r.-6_-3g[6]` — the plain library path leaves the spec's
//!   published form alone, lenient strips the redundant unit label (W3013) and
//!   emits `r.-6_-3[6]`, strict rejects the input. The divergence therefore
//!   belongs to the **preprocessor's repair**, not to the normalizer. It is a
//!   spec self-conflict either way: `RNA/repeated.md:22` calls range-plus-unit
//!   redundant and invalid, `:27` publishes exactly that shape as valid. All
//!   three answers are pinned rather than reconciled — see the undecided ruling
//!   `rna-repeat-range-plus-unit-redundancy`.
//! - `NM_024312.4:c.1738TA[6]` / `r.1738ua[6]` — **strict mode rejects two of
//!   the spec's own published inputs**, because the reference tract is not a
//!   whole number of `TA` copies. Recorded, not endorsed.
//!
//! # The exon-junction pair must converge, and does not yet
//!
//! `LRG_199t1:c.3921dup` and `c.3922dup` denote one transcript sequence (both
//! positions are `T`) but are two fixed points, projecting 2,790 bp apart. The
//! ruling `exon-junction-dup-converge-from-the-far-side` is **`decided`:
//! CONVERGE** — `c.3922dup` normalizes to `c.3921dup`, because the 3'rule's
//! exon-junction clamp binds from **both** sides. Ferro does not do that yet
//! (#1621), so this file pins each row's own output as **pre-fix behaviour**
//! and `the_decided_target_is_convergence_on_the_near_side` asserts the decided
//! answer under `#[ignore]`.
//!
//! That also makes `LRG_199t1:c.3922dup` a **recorded divergence** rather than
//! a row the spec is silent on: its `spec_expected` was null on the reasoning
//! that the spec published no answer, and the spec publishes one three times
//! (`DNA/duplication.md:26`, `:60`, `:148`).
//!
//! # Regenerating the slice
//!
//! `cargo run --features dev --example extract_spec_worked_example_windows -- \
//!  --manifest <manifest>`, and never by hand.

use std::path::{Path, PathBuf};

use ferro_hgvs::conformance::reference_window::{WindowFixture, WindowProvider};
use ferro_hgvs::conformance::spec_worked_examples::{
    Case, Fixture, CASES_PATH, SPEC_DIR, WINDOWS_PATH,
};
use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

/// Which shipped entry point a row is being run through.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Mode {
    /// The plain library API: `parse_hgvs` + `Normalizer::new`. No
    /// input-hygiene preprocessing; the normalizer's own default error config
    /// is lenient.
    Default,
    /// The shipped `ferro normalize --error-mode lenient` path.
    Lenient,
    /// The shipped `ferro normalize --error-mode strict` path.
    Strict,
}

impl Mode {
    fn name(self) -> &'static str {
        match self {
            Mode::Default => "default",
            Mode::Lenient => "lenient",
            Mode::Strict => "strict",
        }
    }

    /// The pinned output for this mode, `None` meaning "errors".
    fn pinned(self, case: &Case) -> Option<&str> {
        match self {
            Mode::Default => case.default_output.as_deref(),
            Mode::Lenient => case.lenient_output.as_deref(),
            Mode::Strict => case.strict_output.as_deref(),
        }
    }
}

/// Anchored on `CARGO_MANIFEST_DIR` so the paths resolve regardless of the
/// test's working directory.
fn repo_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

fn cases() -> Fixture {
    Fixture::from_json_path(&repo_path(CASES_PATH)).expect("load spec worked-example cases")
}

/// The hermetic provider: the committed slice, served through
/// `ReferenceProvider`. No manifest, no environment variable, no skip path.
fn provider() -> WindowProvider {
    WindowFixture::from_json_path(&repo_path(WINDOWS_PATH))
        .expect("load the committed spec worked-example reference slice")
        .to_provider()
}

/// Run one input through one shipped entry point, mirroring `src/bin/ferro.rs`:
/// the error configuration reaches both the preprocessor and the normalizer
/// (#1181/#1197), which is exactly what makes the modes differ here.
fn run(provider: &WindowProvider, input: &str, mode: Mode) -> Result<String, String> {
    let error_config = match mode {
        Mode::Default => None,
        Mode::Lenient => Some(ErrorConfig::lenient()),
        Mode::Strict => Some(ErrorConfig::strict()),
    };
    let (variant, normalize_config) = match error_config {
        None => (
            parse_hgvs(input).map_err(|e| format!("parse error: {e}"))?,
            NormalizeConfig::default(),
        ),
        Some(config) => {
            let parsed = parse_hgvs_with_config(input, config.clone())
                .map_err(|e| format!("preprocess error: {e}"))?;
            (
                parsed.result,
                NormalizeConfig::for_entry_point(ShuffleDirection::ThreePrime, config),
            )
        }
    };
    Normalizer::with_config(provider.clone(), normalize_config)
        .normalize(&variant)
        .map(|n| n.to_string())
        .map_err(|e| format!("normalize error: {e}"))
}

/// Every worked example produces its recorded answer, in every shipped mode.
///
/// Reported as one batch rather than row-by-row: when a normalizer change moves
/// these, which rows moved together — and in which modes — is the diagnosis.
#[test]
fn the_spec_worked_examples_produce_their_recorded_answers() {
    let cases = cases();
    let provider = provider();
    let mut failures = Vec::new();
    let mut checked = 0usize;
    for case in &cases.cases {
        for mode in [Mode::Default, Mode::Lenient, Mode::Strict] {
            checked += 1;
            let actual = run(&provider, &case.input, mode);
            let matches = match (mode.pinned(case), &actual) {
                (Some(pinned), Ok(got)) => pinned == got,
                (None, Err(_)) => true,
                _ => false,
            };
            if !matches {
                failures.push(format!(
                    "  {} [{}]\n    expected: {}\n    actual:   {}\n    {} — {}",
                    case.input,
                    mode.name(),
                    mode.pinned(case).unwrap_or("<error>"),
                    actual.as_deref().unwrap_or_else(|e| e),
                    case.citation.clause,
                    case.note
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {checked} pinned (row, mode) results moved:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// The spec's published answer is what ferro produces, wherever the spec
/// publishes one and the row is not the recorded self-conflict.
///
/// The test above pins ferro against itself, which catches movement but proves
/// nothing about correctness. This one is the correctness half: it compares
/// against the spec's own printed string, so it can only be satisfied by being
/// right.
#[test]
fn the_lenient_path_reproduces_the_spec_where_the_spec_is_self_consistent() {
    let cases = cases();
    let provider = provider();
    let mut mismatched = Vec::new();
    let mut compared = 0usize;
    for case in &cases.cases {
        let Some(spec) = case.spec_expected.as_deref() else {
            continue;
        };
        if !case.lenient_matches_spec() {
            // The recorded divergence; the row below pins that it stays the
            // only one.
            continue;
        }
        compared += 1;
        let actual = run(&provider, &case.input, Mode::Lenient);
        if actual.as_deref() != Ok(spec) {
            mismatched.push(format!(
                "  {} → {:?}, spec says {spec} ({})",
                case.input, actual, case.citation.clause
            ));
        }
    }
    assert_eq!(compared, 10, "the set of spec-conformant rows changed size");
    assert!(
        mismatched.is_empty(),
        "{} spec worked example(s) no longer produce the spec's published answer:\n{}",
        mismatched.len(),
        mismatched.join("\n")
    );
}

/// Normalization is idempotent on every answer these rows produce.
///
/// Cheap, and worth having: each output is itself a legal description, so
/// re-normalizing it must be a no-op. A rule that reaches the right answer by
/// one shift too many and one back would satisfy the tests above and fail this.
#[test]
fn every_worked_example_answer_is_a_fixed_point() {
    let cases = cases();
    let provider = provider();
    for case in &cases.cases {
        for mode in [Mode::Default, Mode::Lenient, Mode::Strict] {
            let Ok(once) = run(&provider, &case.input, mode) else {
                continue;
            };
            let twice = run(&provider, &once, mode);
            assert_eq!(
                twice.as_deref(),
                Ok(once.as_str()),
                "{} [{}] normalized to {once}, which is not a fixed point",
                case.input,
                mode.name()
            );
        }
    }
}

/// Every row cites the spec, and the citation still quotes it.
///
/// A pinned output with no authority behind it is a change detector, not a
/// record: it tells a future reader that something moved, not whether the old
/// answer or the new one was right. Verifying the quote verbatim is what keeps
/// the citation honest across a spec-submodule bump — a bare line number
/// resolves against any file long enough to have that line.
#[test]
fn every_row_cites_the_spec_verbatim() {
    let spec_dir = repo_path(SPEC_DIR);
    assert!(
        spec_dir.join("docs/recommendations/general.md").is_file(),
        "the vendored HGVS spec checkout at {SPEC_DIR} is empty. Initialise it:\n    \
         git submodule update --init {SPEC_DIR}"
    );
    let cases = cases();
    for case in &cases.cases {
        assert!(
            !case.note.trim().is_empty(),
            "{} records no reason for being here",
            case.input
        );
        if let Err(e) = case.citation.verify(&spec_dir) {
            panic!("{}: {e}", case.input);
        }
    }
}

/// The census: twelve rows, exactly one divergence from the spec, exactly one
/// row the spec does not answer.
///
/// Without this, a future edit could quietly demote a conforming row to a
/// divergence — recording the regression as though it had always been one — and
/// every other test here would stay green.
#[test]
fn the_divergence_set_is_exactly_the_one_recorded_conflict() {
    let cases = cases();
    assert_eq!(cases.cases.len(), 12, "the corpus changed size");

    let unanswered: Vec<&str> = cases
        .cases
        .iter()
        .filter(|c| c.spec_expected.is_none())
        .map(|c| c.input.as_str())
        .collect();
    assert!(
        unanswered.is_empty(),
        "the set of rows the spec publishes no answer for changed: {unanswered:?}. Adding one \
         means claiming the spec is silent on an input it may well settle — an adjudication, not \
         a fixture edit. It was `[LRG_199t1:c.3922dup]` until the \
         `exon-junction-dup-converge-from-the-far-side` ruling found the spec answering that \
         input three times over (`DNA/duplication.md:26`, `:60`, `:148`)."
    );

    // Sorted, and the expected literal sorted to match: this set has two members
    // as of the exon-junction ruling, so an ordered comparison would key on row
    // order in the fixture. Inserting an unrelated row between these two would
    // then fail with a message about the *set* having changed when only the
    // order did — the misleading diagnostic every other message here avoids.
    let mut diverging: Vec<&str> = cases
        .cases
        .iter()
        .filter(|c| c.spec_expected.is_some() && !c.lenient_matches_spec())
        .map(|c| c.input.as_str())
        .collect();
    diverging.sort_unstable();
    assert_eq!(
        diverging,
        ["LRG_199t1:c.3922dup", "NM_024312.4:r.-6_-3g[6]"],
        "the set of rows diverging from the spec changed. Two are recorded: the \
         RNA/repeated.md:22-vs-:27 self-conflict, and the exon-junction row that the \
         `exon-junction-dup-converge-from-the-far-side` ruling decided against ferro's current \
         output (#1621). A third is a defect, not a fixture edit."
    );
}

/// The two exon-junction rows are two fixed points for one transcript sequence.
///
/// **Pre-fix behaviour against a decided target**, not an observation any more.
/// The `exon-junction-dup-converge-from-the-far-side` ruling is `decided`:
/// CONVERGE, so the `assert_ne!` below records a defect rather than a choice.
/// It stays because the premise it also checks — both positions carry the same
/// base, so the pair really is two spellings of one sequence — is what the
/// ruling rests on, and because a test that goes green on the day the fix lands
/// is the signal to move the fixture rows with it.
///
/// The decided answer is asserted by
/// [`the_decided_target_is_convergence_on_the_near_side`], which is `#[ignore]`d
/// until #1621. When that lands, **delete this assertion rather than inverting
/// it** — the ignored guard already says the thing worth saying.
#[test]
fn the_exon_junction_pair_does_not_converge_and_that_is_recorded() {
    let fixture = WindowFixture::from_json_path(&repo_path(WINDOWS_PATH))
        .expect("load the committed spec worked-example reference slice");
    let tx = fixture
        .transcripts
        .iter()
        .find(|t| t.id == "LRG_199t1")
        .expect("LRG_199t1 in the committed slice");
    let sequence = tx.sequence.as_deref().expect("LRG_199t1 bases");
    let cds_start = tx.cds_start.expect("LRG_199t1 CDS start") as usize;

    // `c.n` is transcript position `cds_start + n - 1`, 1-based; index is one
    // less again.
    let base_at = |c: usize| sequence.as_bytes()[cds_start + c - 2] as char;
    assert_eq!(
        (base_at(3920), base_at(3921), base_at(3922)),
        ('A', 'T', 'T'),
        "the premise of this record is that c.3921 and c.3922 are the same base, so the two \
         spellings denote one transcript sequence"
    );

    let provider = fixture.to_provider();
    let first = run(&provider, "LRG_199t1:c.3921dup", Mode::Lenient);
    let second = run(&provider, "LRG_199t1:c.3922dup", Mode::Lenient);
    assert_ne!(
        first, second,
        "the exon-junction pair converged — which is the DECIDED target (ruling \
         `exon-junction-dup-converge-from-the-far-side`), so this failure is good news. Delete \
         this assertion, un-ignore `the_decided_target_is_convergence_on_the_near_side`, and \
         move the `LRG_199t1:c.3922dup` rows in `tests/fixtures/spec-worked-examples/cases.json` \
         and `tests/fixtures/case-harvest/cases.json` to `c.3921dup`. See #1621."
    );
}

/// **The decided target**: `LRG_199t1:c.3922dup` normalizes to
/// `LRG_199t1:c.3921dup`, in every mode.
///
/// **Authority.** The `OPERATOR RULING, 2026-08-10` paragraph of ruling record
/// `exon-junction-dup-converge-from-the-far-side`. The canonical position is
/// the most 3' position that does not cross an exon/exon junction, reached from
/// **either** side: `c.3921dup` is at the clamp and stays, `c.3922dup` is past
/// it and is pulled back. `DNA/duplication.md` says so three times — `:26`
/// names `c.3921dup` as the description "and not as `c.3922dup`", `:60`
/// restates it on the GRCh38 coordinates, and `:148` gives the reason, that
/// `c.3922dup` translated back to a genomic position lands "at the wrong
/// nucleotide, in the wrong exon".
///
/// **`#[ignore]`d because ferro does not do this yet**, not because the answer
/// is in doubt. #1617's sibling; the implementation is **#1621**, and this test
/// is its acceptance criterion.
#[test]
#[ignore = "decided target, not yet implemented — see #1621"]
fn the_decided_target_is_convergence_on_the_near_side() {
    let fixture = WindowFixture::from_json_path(&repo_path(WINDOWS_PATH))
        .expect("load the committed spec worked-example reference slice");
    let provider = fixture.to_provider();

    for mode in [Mode::Default, Mode::Lenient, Mode::Strict] {
        assert_eq!(
            run(&provider, "LRG_199t1:c.3922dup", mode),
            Ok("LRG_199t1:c.3921dup".to_string()),
            "{mode:?}: the exon-junction clamp binds from the far side too"
        );
    }
}
