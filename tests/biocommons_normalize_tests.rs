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
//! See `tests/fixtures/biocommons-normalize/failure-patterns.md` for the
//! disposition of currently-known divergences.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::{
    parse_hgvs, FerroError, MultiFastaProvider, NormalizeConfig, Normalizer, ReferenceProvider,
    ShuffleDirection,
};
use serde::Deserialize;
use std::fs;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

const FIXTURE_PATH: &str = "tests/fixtures/biocommons-normalize/cases.json";
const MOCK_PIN_PATH: &str = "tests/fixtures/biocommons-normalize/mock-pin/normalized.txt";
const XFAIL_REPORT_DIR: &str = "/tmp/ferro-xfail";
const FAIL_PRINT_LIMIT: usize = 10;

// ----------------------------------------------------------------------------
// Fixture model
// ----------------------------------------------------------------------------

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct Fixture {
    description: String,
    source: String,
    source_commit: String,
    biocommons_upstream: String,
    license: String,
    refreshed_at: String,
    cases: Vec<Case>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct Case {
    input: String,
    /// `None` means biocommons expected `normalize()` to raise an error.
    /// Mutually informed by `expects_error: true`.
    normalized: Option<String>,
    #[serde(default)]
    expects_error: bool,
    shuffle_direction: String,
    cross_boundaries: bool,
    source_function: String,
    #[serde(default = "default_true")]
    to_test: bool,
}

fn default_true() -> bool {
    true
}

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
    fn get_transcript(&self, id: &str) -> Result<Transcript, FerroError> {
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

#[derive(Debug)]
struct AxisTally {
    axis: &'static str,
    pass: usize,
    fail: Vec<(String, String)>,
    skipped: usize,
}

impl AxisTally {
    fn new(axis: &'static str) -> Self {
        Self {
            axis,
            pass: 0,
            fail: Vec::new(),
            skipped: 0,
        }
    }

    fn record(&mut self, case_input: &str, expected: &str, actual: Result<String, String>) {
        if matches!(&actual, Ok(s) if s == expected) {
            self.pass += 1;
        } else {
            let diag = match actual {
                Ok(got) => format!("expected={expected:?} got={got:?}"),
                Err(e) => format!("expected={expected:?} err={e}"),
            };
            self.fail.push((case_input.to_string(), diag));
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
            "{}: {} pass / {} FAIL / {} skipped (FAIL inputs -> {})",
            self.axis,
            self.pass,
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

#[test]
fn axis_normalized() {
    let Some(provider) = provider() else {
        eprintln!("axis_normalized: skipping — no manifest");
        return;
    };

    let mut t = AxisTally::new("normalized");
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        // biocommons expects-error cases: ferro should also error.
        if case.expects_error {
            let normalizer =
                Normalizer::with_config(ArcProvider(provider.clone()), build_config(case));
            let actual = catch_panics(|| -> Result<String, String> {
                let v = parse_hgvs(&case.input).map_err(|e| format!("parse error: {e}"))?;
                match normalizer.normalize(&v) {
                    Err(_) => Ok("<expects error>".to_string()),
                    Ok(n) => Err(format!("accepted: {n}")),
                }
            });
            t.record(&case.input, "<expects error>", actual);
            continue;
        }
        let Some(expected) = case.normalized.as_deref() else {
            t.skipped += 1;
            continue;
        };

        let normalizer = Normalizer::with_config(ArcProvider(provider.clone()), build_config(case));
        let actual = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse error: {e}"))?;
            let n = normalizer
                .normalize(&v)
                .map_err(|e| format!("normalize error: {e}"))?;
            Ok(format!("{n}"))
        });

        t.record(&case.input, expected, actual);
    }
    t.finish();
}
