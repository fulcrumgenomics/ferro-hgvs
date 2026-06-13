//! Biocommons normalizer test cases.
//!
//! Imports tests from biocommons/hgvs's `tests/test_hgvs_normalizer.py` (via
//! the hgvs-rs Rust port at `src/normalizer.rs`) plus six normalize-bearing
//! regressions from biocommons/hgvs `tests/issues/test_02xx.py` and
//! `test_03xx.py`. Pinned via `scripts/refresh-biocommons-fixtures.py`.
//!
//! Two-layer test pattern per `tests/fixtures/CORPUS_LAYOUT.md`:
//!
//! 1. **`regression_under_mock_normalized`** (CI-always): runs every case
//!    under `MockProvider` and diffs ferro's output against the pinned
//!    `mock-pin/normalized.txt` snapshot. Catches refactor-side regressions
//!    in mock-mode behavior. The pin records *ferro's current behavior*
//!    (often "input unchanged" — MockProvider has no reference bases),
//!    not upstream truth. Regenerate with:
//!    `BLESS_MOCK_PIN=1 cargo nextest run --features dev -E 'test(regression_under_mock_normalized) & binary(biocommons_normalize_tests)'`
//!
//! 2. **`axis_normalized`** (manifest-or-skip): runs every case with a
//!    real `MultiFastaProvider` loaded via the same path convention used
//!    by `tests/mutalyzer_normalize_tests.rs` and
//!    `tests/real_data_normalization_tests.rs`. Strict-asserts ferro's
//!    output equals the biocommons-expected `normalized` field; FAIL inputs
//!    are written to `/tmp/ferro-xfail/biocommons-normalized.{txt,tsv}` for
//!    burn-down. Tests skip when the manifest is absent (e.g. CI).
//!
//! See `cases.json` (and the generated `failure-patterns.md` summary, produced
//! by `examples/generate_conformance_summary`, never hand-maintained — #509)
//! for the disposition of currently-known divergences.

use ferro_hgvs::conformance::biocommons::{Case, Disposition, Fixture};
use ferro_hgvs::conformance::reference_window::WindowFixture;
use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::{
    parse_hgvs, FerroError, MultiFastaProvider, NormalizeConfig, Normalizer, ReferenceProvider,
    ShuffleDirection,
};
use std::fs;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

const FIXTURE_PATH: &str = "tests/fixtures/biocommons-normalize/cases.json";
const WINDOWS_FIXTURE_PATH: &str = "tests/fixtures/biocommons-normalize/reference-windows.json";
const MOCK_PIN_PATH: &str = "tests/fixtures/biocommons-normalize/mock-pin/normalized.txt";
const XFAIL_REPORT_DIR: &str = "/tmp/ferro-xfail";
const FAIL_PRINT_LIMIT: usize = 10;

fn direction_from_str(s: &str) -> ShuffleDirection {
    match s {
        "3prime" => ShuffleDirection::ThreePrime,
        "5prime" => ShuffleDirection::FivePrime,
        other => panic!("unknown shuffle_direction in fixture: {other}"),
    }
}

fn build_config(case: &Case) -> NormalizeConfig {
    let mut cfg =
        NormalizeConfig::default().with_direction(direction_from_str(&case.shuffle_direction));
    if case.cross_boundaries {
        cfg = cfg.allow_crossing_boundaries();
    }
    cfg
}

// ----------------------------------------------------------------------------
// Shared state
// ----------------------------------------------------------------------------

fn fixture() -> &'static Fixture {
    static FIXTURE: OnceLock<Fixture> = OnceLock::new();
    FIXTURE.get_or_init(|| {
        let content = fs::read_to_string(FIXTURE_PATH)
            .unwrap_or_else(|e| panic!("failed to read {FIXTURE_PATH}: {e}"));
        serde_json::from_str(&content)
            .unwrap_or_else(|e| panic!("failed to parse {FIXTURE_PATH}: {e}"))
    })
}

fn manifest_path() -> Option<PathBuf> {
    // Same convention as tests/mutalyzer_normalize_tests.rs and
    // tests/real_data_normalization_tests.rs. FERRO_MANIFEST is authoritative
    // when set; otherwise fall back to well-known paths.
    if let Ok(path) = std::env::var("FERRO_MANIFEST") {
        let p = PathBuf::from(path);
        return if p.exists() { Some(p) } else { None };
    }
    for candidate in [
        "/Volumes/scratch-00001/work/clients/fulcrum/ferro-hgvs/data/ferro/manifest.json",
        "benchmark-output/manifest.json",
    ] {
        let p = Path::new(candidate);
        if p.exists() {
            return Some(p.to_path_buf());
        }
    }
    None
}

fn provider() -> Option<Arc<MultiFastaProvider>> {
    static PROVIDER: OnceLock<Option<Arc<MultiFastaProvider>>> = OnceLock::new();
    PROVIDER
        .get_or_init(|| {
            // `manifest_path()` returning `None` is the legitimate "fixtures not
            // generated yet" skip path. But if the manifest IS present, loading
            // it must succeed — otherwise we'd silently mask a misconfigured
            // manifest and report broken setups as no-op skips.
            let path = manifest_path()?;
            Some(Arc::new(
                MultiFastaProvider::from_manifest(&path)
                    .unwrap_or_else(|e| panic!("from_manifest({}) failed: {e}", path.display())),
            ))
        })
        .clone()
}

/// `Arc<MultiFastaProvider>` newtype that implements `ReferenceProvider` so it
/// can be plugged into `Normalizer::with_config` without cloning the
/// underlying index for each case.
#[derive(Clone)]
struct ArcProvider(Arc<MultiFastaProvider>);

impl ReferenceProvider for ArcProvider {
    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError> {
        self.0.get_transcript(id)
    }
    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        self.0.get_sequence(id, start, end)
    }
    fn has_transcript(&self, id: &str) -> bool {
        self.0.has_transcript(id)
    }
    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        self.0.get_genomic_sequence(contig, start, end)
    }
    fn has_genomic_data(&self) -> bool {
        self.0.has_genomic_data()
    }
    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        self.0.get_protein_sequence(accession, start, end)
    }
    fn has_protein_data(&self) -> bool {
        self.0.has_protein_data()
    }
}

/// Catch panics inside a per-case closure so a single bad case can't abort
/// the whole tally.
fn catch_panics(body: impl FnOnce() -> Result<String, String>) -> Result<String, String> {
    match catch_unwind(AssertUnwindSafe(body)) {
        Ok(r) => r,
        Err(payload) => {
            let msg = payload
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| {
                    payload
                        .downcast_ref::<&'static str>()
                        .map(|s| s.to_string())
                })
                .unwrap_or_else(|| "<non-string panic payload>".to_string());
            Err(format!("panic: {msg}"))
        }
    }
}

// ----------------------------------------------------------------------------
// Smoke tests
// ----------------------------------------------------------------------------

#[test]
fn loads_fixture() {
    let f = fixture();
    let active = f.cases.iter().filter(|c| c.to_test).count();
    println!(
        "biocommons-normalize: loaded {} total cases ({} to_test) @ commit {}",
        f.cases.len(),
        active,
        f.source_commit
    );
    assert!(active > 0, "fixture should have at least one to_test case");
}

#[test]
fn manifest_or_skip() {
    match manifest_path() {
        Some(p) => println!("biocommons-normalize: using manifest {}", p.display()),
        None => println!(
            "biocommons-normalize: skipping — no manifest at FERRO_MANIFEST, \
             /Volumes/scratch-00001/.../manifest.json, or benchmark-output/manifest.json"
        ),
    }
}

// ----------------------------------------------------------------------------
// Layer 1: regression_under_mock_normalized — MockProvider pin (CI-always)
// ----------------------------------------------------------------------------

fn render_mock_normalized() -> String {
    let mut lines: Vec<String> = Vec::new();
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let normalizer = Normalizer::with_config(MockProvider::new(), build_config(case));
        let output = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse error: {e}"))?;
            let n = normalizer
                .normalize(&v)
                .map_err(|e| format!("normalize error: {e}"))?;
            Ok(format!("{n}"))
        })
        .unwrap_or_else(|e| e);
        lines.push(format!(
            "{}\t{}\t{}\t{}",
            case.input, case.shuffle_direction, case.cross_boundaries, output
        ));
    }
    let mut out = lines.join("\n");
    out.push('\n');
    out
}

#[test]
fn regression_under_mock_normalized() {
    let observed = render_mock_normalized();

    if std::env::var("BLESS_MOCK_PIN").is_ok() {
        fs::write(MOCK_PIN_PATH, &observed)
            .unwrap_or_else(|e| panic!("write {MOCK_PIN_PATH}: {e}"));
        eprintln!("blessed {MOCK_PIN_PATH} ({} bytes)", observed.len());
        return;
    }

    let expected = fs::read_to_string(MOCK_PIN_PATH).unwrap_or_else(|e| {
        panic!(
            "read {MOCK_PIN_PATH}: {e}; if this is a fresh corpus, run \
             `BLESS_MOCK_PIN=1 cargo nextest run --features dev \
              -E 'test(regression_under_mock_normalized) & binary(biocommons_normalize_tests)'`"
        )
    });

    if observed == expected {
        let total = observed.lines().count();
        println!(
            "regression_under_mock_normalized: {} cases pinned, all match",
            total
        );
        return;
    }

    let mut diffs = Vec::new();
    for (i, (o, e)) in observed.lines().zip(expected.lines()).enumerate() {
        if o != e {
            diffs.push(format!("  line {i}:\n    pinned:   {e}\n    observed: {o}"));
            if diffs.len() >= FAIL_PRINT_LIMIT {
                break;
            }
        }
    }
    let observed_total = observed.lines().count();
    let expected_total = expected.lines().count();
    if observed_total != expected_total {
        diffs.push(format!(
            "  line count differs: observed={observed_total}, expected={expected_total} \
             (cases.json may have changed without a re-bless)"
        ));
    }

    panic!(
        "mock pin drifted ({} divergence(s) shown).\n\n{}\n\n\
         If this drift is intentional, regenerate the pin:\n  \
         BLESS_MOCK_PIN=1 cargo nextest run --features dev \
         -E 'test(regression_under_mock_normalized) & binary(biocommons_normalize_tests)'",
        diffs.len(),
        diffs.join("\n"),
    );
}

// ----------------------------------------------------------------------------
// Layer 2: axis_normalized — manifest-mode strict assertion
// ----------------------------------------------------------------------------

/// Outcome of comparing ferro's output for one row against biocommons and the
/// row's [`Disposition`]. See [`classify_outcome`].
#[derive(Debug, Clone, PartialEq, Eq)]
enum RecordOutcome {
    /// ferro matches biocommons; no annotation needed.
    Match,
    /// ferro diverges exactly as an annotated, spec-correct accepted divergence.
    AcceptedDivergence,
    /// ferro diverges exactly as an annotated known bug (xfail).
    KnownBug,
    /// ferro diverges exactly as an annotated tracked improvement (xfail).
    Improvement,
    /// A failure that must break the suite.
    Fail(FailKind),
}

/// The distinct ways a row can fail the conformance suite.
#[derive(Debug, Clone, PartialEq, Eq)]
enum FailKind {
    /// Divergence from biocommons with no annotation explaining it.
    Unannotated,
    /// A `known_bug` row now matches biocommons — the bug appears fixed (XPASS);
    /// remove the annotation and demote the row.
    KnownBugFixed,
    /// An `improvement` row now matches biocommons — ferro converged to the
    /// spec-preferred form (XPASS); remove the annotation and demote the row.
    ImprovementConverged,
    /// An `accepted_divergence` row now matches biocommons — the divergence is
    /// gone; remove the annotation.
    AcceptedDivergenceGone,
    /// ferro's output no longer equals the annotation's recorded `ferro_output`.
    AnnotationDrift,
}

impl FailKind {
    fn label(&self) -> &'static str {
        match self {
            FailKind::Unannotated => "UNANNOTATED",
            FailKind::KnownBugFixed => "XPASS",
            FailKind::ImprovementConverged => "XPASS",
            FailKind::AcceptedDivergenceGone => "XPASS",
            FailKind::AnnotationDrift => "DRIFT",
        }
    }
}

/// Pure classification of one row, given biocommons' `expected` value, ferro's
/// `actual` result, and the row's optional `disposition`. The heart of #478
/// pillars 1-2: unannotated divergences fail; an annotation that still holds is
/// accounted for; an annotation that no longer holds (XPASS or drift) fails so
/// it gets corrected.
fn classify_outcome(
    expected: &str,
    actual: &Result<String, String>,
    disposition: Option<&Disposition>,
) -> RecordOutcome {
    let actual_str = actual.as_deref().ok();
    let matches_upstream = actual_str == Some(expected);
    match disposition {
        None => {
            if matches_upstream {
                RecordOutcome::Match
            } else {
                RecordOutcome::Fail(FailKind::Unannotated)
            }
        }
        Some(Disposition::KnownBug { ferro_output, .. }) => {
            if matches_upstream {
                RecordOutcome::Fail(FailKind::KnownBugFixed)
            } else if actual_str == Some(ferro_output.as_str()) {
                RecordOutcome::KnownBug
            } else {
                RecordOutcome::Fail(FailKind::AnnotationDrift)
            }
        }
        Some(Disposition::Improvement { ferro_output, .. }) => {
            if matches_upstream {
                RecordOutcome::Fail(FailKind::ImprovementConverged)
            } else if actual_str == Some(ferro_output.as_str()) {
                RecordOutcome::Improvement
            } else {
                RecordOutcome::Fail(FailKind::AnnotationDrift)
            }
        }
        Some(Disposition::AcceptedDivergence { ferro_output, .. }) => {
            if matches_upstream {
                RecordOutcome::Fail(FailKind::AcceptedDivergenceGone)
            } else if actual_str == Some(ferro_output.as_str()) {
                RecordOutcome::AcceptedDivergence
            } else {
                RecordOutcome::Fail(FailKind::AnnotationDrift)
            }
        }
    }
}

#[derive(Debug)]
struct AxisTally {
    axis: &'static str,
    pass: usize,
    /// Annotated, spec-correct divergences whose `ferro_output` still holds.
    accepted: usize,
    /// Annotated known bugs (xfail) still producing their recorded output.
    known_bug: usize,
    /// Annotated tracked improvements (xfail) still producing their recorded
    /// spec-non-preferred output.
    improvement: usize,
    fail: Vec<(String, String)>,
    skipped: usize,
}

impl AxisTally {
    fn new(axis: &'static str) -> Self {
        Self {
            axis,
            pass: 0,
            accepted: 0,
            known_bug: 0,
            improvement: 0,
            fail: Vec::new(),
            skipped: 0,
        }
    }

    fn record(
        &mut self,
        case_input: &str,
        expected: &str,
        actual: Result<String, String>,
        disposition: Option<&Disposition>,
    ) {
        match classify_outcome(expected, &actual, disposition) {
            RecordOutcome::Match => self.pass += 1,
            RecordOutcome::AcceptedDivergence => self.accepted += 1,
            RecordOutcome::KnownBug => self.known_bug += 1,
            RecordOutcome::Improvement => self.improvement += 1,
            RecordOutcome::Fail(kind) => {
                let got = match &actual {
                    Ok(s) => format!("got={s:?}"),
                    Err(e) => format!("err={e}"),
                };
                // Use the disposition's own fields to make the failure
                // actionable (and to read them, satisfying dead-code lints).
                let detail = match (&kind, disposition) {
                    (
                        FailKind::KnownBugFixed,
                        Some(Disposition::KnownBug { tracking_issue, .. }),
                    ) => format!(
                        " — known_bug #{tracking_issue} now matches biocommons; the fix appears \
                         to have landed: remove the annotation and demote the row"
                    ),
                    (
                        FailKind::ImprovementConverged,
                        Some(Disposition::Improvement { tracking_issue, .. }),
                    ) => format!(
                        " — improvement #{tracking_issue} now matches biocommons; ferro converged \
                         to the spec-preferred form: remove the annotation and demote the row"
                    ),
                    (
                        FailKind::AcceptedDivergenceGone,
                        Some(Disposition::AcceptedDivergence { reason, .. }),
                    ) => format!(
                        " — accepted_divergence no longer diverges (recorded reason: {reason}); \
                         remove the annotation"
                    ),
                    (FailKind::AnnotationDrift, _) => {
                        " — output no longer matches the recorded ferro_output; update or remove \
                         the annotation"
                            .to_string()
                    }
                    _ => String::new(),
                };
                self.fail.push((
                    case_input.to_string(),
                    format!("[{}] expected={expected:?} {got}{detail}", kind.label()),
                ));
            }
        }
    }

    fn finish(self) {
        let dir = Path::new(XFAIL_REPORT_DIR);
        let _ = fs::create_dir_all(dir);
        let report_path = dir.join(format!("biocommons-{}.txt", self.axis));
        let tsv_path = dir.join(format!("biocommons-{}.tsv", self.axis));
        let body: String = self.fail.iter().map(|(i, _)| format!("{i}\n")).collect();
        let _ = fs::write(&report_path, &body);
        let tsv: String = self
            .fail
            .iter()
            .map(|(i, d)| format!("{i}\t{d}\n"))
            .collect();
        let _ = fs::write(&tsv_path, &tsv);

        println!(
            "{}: {} pass / {} accepted-divergence / {} known-bug / {} improvement / {} FAIL / \
             {} skipped (FAIL inputs -> {})",
            self.axis,
            self.pass,
            self.accepted,
            self.known_bug,
            self.improvement,
            self.fail.len(),
            self.skipped,
            report_path.display()
        );

        for (input, diag) in self.fail.iter().take(FAIL_PRINT_LIMIT) {
            eprintln!("  FAIL  [{}] {} | {}", self.axis, input, diag);
        }
        if self.fail.len() > FAIL_PRINT_LIMIT {
            eprintln!(
                "  ... {} more FAIL — full list in {}",
                self.fail.len() - FAIL_PRINT_LIMIT,
                report_path.display()
            );
        }

        assert!(
            self.fail.is_empty(),
            "{}: {} divergence(s) from biocommons — see {} and tests/fixtures/biocommons-normalize/failure-patterns.md",
            self.axis,
            self.fail.len(),
            report_path.display()
        );
    }
}

/// Run the `normalized` axis for every `to_test` case against `provider`,
/// returning the tally (the caller asserts via [`AxisTally::finish`]). Shared by
/// the per-PR hermetic gate ([`axis_normalized_hermetic`]) and the manifest
/// nightly tier ([`axis_normalized`]) so both apply identical config,
/// expects-error handling, and disposition classification — the only difference
/// is the reference provider behind them.
fn run_normalized_axis<P: ReferenceProvider + Clone>(provider: P) -> AxisTally {
    let mut t = AxisTally::new("normalized");
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        // biocommons expects-error cases: ferro should also error. Run in strict
        // mode so the W4xxx / W5xxx rejectable diagnostics (`PositionPastEnd`,
        // `VariantExceedsReference`, etc.) fire as typed errors rather than
        // lenient-mode warnings — otherwise every rejectable W-code ferro emits
        // would round-trip as `Ok(canonicalized_input)` and fail the harness.
        if case.expects_error {
            let cfg = build_config(case).with_error_mode(ErrorMode::Strict);
            let normalizer = Normalizer::with_config(provider.clone(), cfg);
            // Map ferro's outcome onto the `<expects error>` sentinel so
            // `classify_outcome` treats this row uniformly with the value-axis
            // rows. ferro erroring is the match (`Ok("<expects error>")`);
            // ferro *accepting* is the divergence, and we surface the accepted
            // HGVS string as the `Ok` value so an `accepted_divergence`
            // annotation can pin it (e.g. #253: biocommons raises, ferro
            // spec-correctly accepts `NM_001166478.1:c.59_61del`). Returning the
            // accepted form as `Ok` — not `Err` — is what lets the disposition's
            // `ferro_output` compare equal; an unannotated accept still fails as
            // `Unannotated` (it can never equal the `<expects error>` sentinel).
            let actual = catch_panics(|| -> Result<String, String> {
                let v = parse_hgvs(&case.input).map_err(|e| format!("parse error: {e}"))?;
                match normalizer.normalize(&v) {
                    Err(_) => Ok("<expects error>".to_string()),
                    Ok(n) => Ok(n.to_string()),
                }
            });
            t.record(
                &case.input,
                "<expects error>",
                actual,
                case.disposition.as_ref(),
            );
            continue;
        }
        let Some(expected) = case.normalized.as_deref() else {
            t.skipped += 1;
            continue;
        };

        let normalizer = Normalizer::with_config(provider.clone(), build_config(case));
        let actual = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse error: {e}"))?;
            let n = normalizer
                .normalize(&v)
                .map_err(|e| format!("normalize error: {e}"))?;
            Ok(format!("{n}"))
        });

        t.record(&case.input, expected, actual, case.disposition.as_ref());
    }
    t
}

/// Per-PR hermetic merge gate (CI-always): runs the `normalized` axis against a
/// [`WindowProvider`] built from the committed `reference-windows.json` — the
/// exact reference bases the manifest pass touches (transcripts whole, padded
/// genomic windows), with zero out-of-band data (#478 pillar 4). This is the
/// gate that blocks PRs. Regenerate the fixture from the manifest with
/// `cargo run --features dev --example extract_biocommons_windows` whenever
/// cases.json or normalize behavior changes the reference access set.
#[test]
fn axis_normalized_hermetic() {
    let window_fixture = WindowFixture::from_json_path(Path::new(WINDOWS_FIXTURE_PATH))
        .unwrap_or_else(|e| panic!("load {WINDOWS_FIXTURE_PATH}: {e}"));
    run_normalized_axis(window_fixture.to_provider()).finish();
}

/// Manifest-mode nightly tier: the same axis against the full
/// `MultiFastaProvider`. Skips when the manifest is absent (e.g. CI) — the
/// hermetic gate above is the per-PR merge gate. This tier catches anything the
/// committed windows miss and is the source the extraction tool regenerates
/// from; it must agree with the hermetic gate row-for-row.
#[test]
fn axis_normalized() {
    let Some(provider) = provider() else {
        eprintln!(
            "axis_normalized: skipping — no manifest (per-PR gate is axis_normalized_hermetic)"
        );
        return;
    };
    run_normalized_axis(ArcProvider(provider)).finish();
}

// ----------------------------------------------------------------------------
// #478 pillars 1-2: classify_outcome — annotation + XPASS unit tests.
// These are hermetic (no manifest) and run in CI on every PR.
// ----------------------------------------------------------------------------

mod classify_outcome_tests {
    use super::*;

    fn ok(s: &str) -> Result<String, String> {
        Ok(s.to_string())
    }

    fn accepted(ferro_output: &str) -> Disposition {
        Disposition::AcceptedDivergence {
            reason: "spec-correct per HGVS edit-type priority".to_string(),
            spec_citation: None,
            ferro_output: ferro_output.to_string(),
            cluster: None,
        }
    }

    fn known_bug(ferro_output: &str) -> Disposition {
        Disposition::KnownBug {
            tracking_issue: 999,
            ferro_output: ferro_output.to_string(),
            cluster: None,
        }
    }

    fn improvement(ferro_output: &str) -> Disposition {
        Disposition::Improvement {
            tracking_issue: 500,
            section: "HGVS §RefSeqGene transcript selection".to_string(),
            ferro_output: ferro_output.to_string(),
            cluster: None,
        }
    }

    #[test]
    fn no_disposition_match_is_match() {
        assert_eq!(
            classify_outcome("c.5dup", &ok("c.5dup"), None),
            RecordOutcome::Match
        );
    }

    #[test]
    fn no_disposition_divergence_is_unannotated_fail() {
        assert_eq!(
            classify_outcome("c.5dup", &ok("c.4_5dup"), None),
            RecordOutcome::Fail(FailKind::Unannotated)
        );
    }

    #[test]
    fn known_bug_producing_recorded_output_is_known_bug() {
        // ferro still produces the recorded wrong answer → accounted xfail.
        assert_eq!(
            classify_outcome(
                "c.36_37dup",
                &ok("c.35_36dup"),
                Some(&known_bug("c.35_36dup"))
            ),
            RecordOutcome::KnownBug
        );
    }

    #[test]
    fn known_bug_now_matching_upstream_is_xpass_fail() {
        // THE KEYSTONE: a fixed known_bug must fail loudly so it can't linger.
        assert_eq!(
            classify_outcome(
                "c.36_37dup",
                &ok("c.36_37dup"),
                Some(&known_bug("c.35_36dup"))
            ),
            RecordOutcome::Fail(FailKind::KnownBugFixed)
        );
    }

    #[test]
    fn known_bug_drifted_output_is_drift_fail() {
        // ferro's wrong answer changed to a third value → annotation is stale.
        assert_eq!(
            classify_outcome(
                "c.36_37dup",
                &ok("c.34_35dup"),
                Some(&known_bug("c.35_36dup"))
            ),
            RecordOutcome::Fail(FailKind::AnnotationDrift)
        );
    }

    #[test]
    fn improvement_producing_recorded_output_is_improvement() {
        // ferro still produces the recorded spec-non-preferred output →
        // accounted xfail (tracked to convergence).
        assert_eq!(
            classify_outcome(
                "c.12_13dup",
                &ok("c.13_14dup"),
                Some(&improvement("c.13_14dup"))
            ),
            RecordOutcome::Improvement
        );
    }

    #[test]
    fn improvement_now_matching_upstream_is_xpass_fail() {
        // A converged improvement must fail loudly so it can't linger — the
        // same XPASS keystone as known_bug.
        assert_eq!(
            classify_outcome(
                "c.12_13dup",
                &ok("c.12_13dup"),
                Some(&improvement("c.13_14dup"))
            ),
            RecordOutcome::Fail(FailKind::ImprovementConverged)
        );
    }

    #[test]
    fn improvement_drifted_output_is_drift_fail() {
        // ferro's output changed to a third value → annotation is stale.
        assert_eq!(
            classify_outcome(
                "c.12_13dup",
                &ok("c.11_12dup"),
                Some(&improvement("c.13_14dup"))
            ),
            RecordOutcome::Fail(FailKind::AnnotationDrift)
        );
    }

    #[test]
    fn accepted_divergence_producing_recorded_output_is_accepted() {
        assert_eq!(
            classify_outcome(
                "c.-1_1insAC",
                &ok("c.-1_1dup"),
                Some(&accepted("c.-1_1dup"))
            ),
            RecordOutcome::AcceptedDivergence
        );
    }

    #[test]
    fn accepted_divergence_now_matching_upstream_is_xpass_fail() {
        assert_eq!(
            classify_outcome(
                "c.-1_1insAC",
                &ok("c.-1_1insAC"),
                Some(&accepted("c.-1_1dup"))
            ),
            RecordOutcome::Fail(FailKind::AcceptedDivergenceGone)
        );
    }

    #[test]
    fn accepted_divergence_drifted_output_is_drift_fail() {
        assert_eq!(
            classify_outcome(
                "c.-1_1insAC",
                &ok("c.1delinsACA"),
                Some(&accepted("c.-1_1dup"))
            ),
            RecordOutcome::Fail(FailKind::AnnotationDrift)
        );
    }

    #[test]
    fn expects_error_accept_with_accepted_divergence_is_accepted() {
        // The expects-error axis encodes "ferro errored" as the `<expects error>`
        // sentinel and "ferro accepted" as the accepted HGVS string (see the
        // expects_error branch of `axis_normalized`). An `accepted_divergence`
        // pinning that accepted form must classify as AcceptedDivergence — this
        // is #253 (biocommons raises, ferro spec-correctly accepts).
        assert_eq!(
            classify_outcome(
                "<expects error>",
                &ok("NM_001166478.1:c.59_61del"),
                Some(&accepted("NM_001166478.1:c.59_61del"))
            ),
            RecordOutcome::AcceptedDivergence
        );
    }

    #[test]
    fn expects_error_when_ferro_errors_with_accepted_divergence_is_xpass() {
        // If an accepted_divergence expects ferro to accept but ferro now errors,
        // the recorded divergence is gone — XPASS, so the annotation gets removed.
        assert_eq!(
            classify_outcome(
                "<expects error>",
                &ok("<expects error>"),
                Some(&accepted("NM_001166478.1:c.59_61del"))
            ),
            RecordOutcome::Fail(FailKind::AcceptedDivergenceGone)
        );
    }

    #[test]
    fn annotated_row_that_errors_is_drift_fail() {
        // An annotation recording a concrete ferro_output, but ferro now errors,
        // is drift worth surfacing.
        let err: Result<String, String> = Err("normalize error: boom".to_string());
        assert_eq!(
            classify_outcome("c.5dup", &err, Some(&known_bug("c.4_5dup"))),
            RecordOutcome::Fail(FailKind::AnnotationDrift)
        );
    }

    #[test]
    fn improvement_error_is_drift_fail() {
        // An annotated improvement that records a concrete ferro_output but now
        // errors is drift — a normalize failure must not be silently absorbed
        // into the improvement xfail bucket.
        let err: Result<String, String> = Err("normalize error: boom".to_string());
        assert_eq!(
            classify_outcome("c.12_13dup", &err, Some(&improvement("c.13_14dup"))),
            RecordOutcome::Fail(FailKind::AnnotationDrift)
        );
    }
}
