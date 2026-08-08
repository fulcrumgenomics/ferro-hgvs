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
//! # The exon-junction pair is deliberately not asserted equal
//!
//! `LRG_199t1:c.3921dup` and `c.3922dup` denote one transcript sequence (both
//! positions are `T`) but are two fixed points, projecting 2,790 bp apart. This
//! file pins each row's own output and does **not** assert they converge —
//! whether they should is an open interpretive question, recorded as the
//! undecided ruling `exon-junction-dup-converge-from-the-far-side`.
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
    assert_eq!(
        unanswered,
        ["LRG_199t1:c.3922dup"],
        "the set of rows the spec publishes no answer for changed. Adding one means claiming the \
         spec is silent on an input it may well settle — an adjudication, not a fixture edit."
    );

    let diverging: Vec<&str> = cases
        .cases
        .iter()
        .filter(|c| c.spec_expected.is_some() && !c.lenient_matches_spec())
        .map(|c| c.input.as_str())
        .collect();
    assert_eq!(
        diverging,
        ["NM_024312.4:r.-6_-3g[6]"],
        "the set of rows diverging from the spec changed. The only recorded divergence is the \
         RNA/repeated.md:22-vs-:27 self-conflict; a second one is a defect, not a fixture edit."
    );
}

/// The two exon-junction rows are two fixed points for one transcript sequence.
///
/// Stated as an observation, not a requirement — see the module docs and the
/// `exon-junction-dup-converge-from-the-far-side` ruling. What is asserted is
/// only that the situation is what the record says it is: both positions carry
/// the same base, so the pair really is two spellings of one sequence, and
/// neither row is shifting.
///
/// If a future change *does* converge them, this test fails and points at the
/// ruling — which is the correct outcome, because converging them is a
/// representation change that needs a decision behind it, not a side effect.
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
        "the exon-junction pair converged. That may well be right — but it is a representation \
         change and an open question (ruling `exon-junction-dup-converge-from-the-far-side`), so \
         it must be decided deliberately rather than land as a side effect."
    );
}
