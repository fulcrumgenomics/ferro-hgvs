//! Mutalyzer normalizer test cases.
//!
//! Imports tests from mutalyzer/mutalyzer's `tests/variants_set.py` (pinned
//! to a specific SHA via `scripts/refresh-mutalyzer-fixtures.py`) and runs
//! them against ferro-hgvs's normalizer + projector.
//!
//! Strict mode: any divergence between ferro-hgvs output and the upstream
//! expected output fails the test loudly. There is no xfail / xpass
//! mechanism — burn-down happens in follow-up PRs that fix ferro-hgvs and
//! demote inputs from
//! `tests/fixtures/mutalyzer-normalize/baseline-failures/<axis>.txt`.
//!
//! Backed by the ferro-prepared reference manifest. When the manifest is
//! absent (CI by default), every axis test reports `skipping — no manifest`
//! and exits 0. Same convention as `tests/real_data_normalization_tests.rs`.
//!
//! Each axis test also writes its current FAIL list to
//! `/tmp/ferro-xfail/<axis>.{txt,tsv}` so contributors can regenerate the
//! committed `baseline-failures/` snapshot or feed
//! `tests/fixtures/mutalyzer-normalize/failure-patterns.md`.
//!
//! **Regression layer (CI-always):** `regression_under_mock_normalized` runs
//! every `to_test` case with a `normalized` expectation through ferro under
//! `MockProvider` and diffs against `mock-pin/normalized.txt`. This catches
//! refactor-side regressions even on CI without a manifest. The pin file is
//! intentionally distinct from `cases.json` — `cases.json` holds upstream
//! truth (durable), `mock-pin/normalized.txt` holds ferro's current
//! MockProvider behavior (ephemeral, regenerable). See
//! `tests/fixtures/CORPUS_LAYOUT.md` for the unified shape.
//!
//! Regenerate the mock pin (after intentional changes):
//!     `BLESS_MOCK_PIN=1 cargo nextest run --features dev -E 'test(regression_under_mock_normalized)'`

use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::{
    parse_hgvs, FerroError, HgvsVariant, MultiFastaProvider, Normalizer, ReferenceProvider,
    VariantProjector,
};
use serde::Deserialize;
use std::fs;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

/// Run `body` and convert any panic into a `Result::Err("panic: …")` so a
/// single bad case (e.g. an arithmetic overflow inside the library) cannot
/// abort the whole axis test before it has a chance to tally.
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

const FIXTURE_PATH: &str = "tests/fixtures/mutalyzer-normalize/cases.json";
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
    license: String,
    refreshed_at: String,
    cases: Vec<Case>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct Case {
    #[serde(default)]
    keywords: Vec<String>,
    input: String,
    #[serde(default)]
    normalized: Option<String>,
    #[serde(default)]
    genomic: Option<String>,
    #[serde(default)]
    coding_protein_descriptions: Option<Vec<Vec<String>>>,
    #[serde(default)]
    protein_description: Option<String>,
    #[serde(default)]
    rna_description: Option<String>,
    #[serde(default)]
    errors: Option<Vec<String>>,
    #[serde(default)]
    infos: Option<Vec<String>>,
    #[serde(default)]
    noncoding: Option<Vec<String>>,
    #[serde(default = "default_true")]
    to_test: bool,
}

fn default_true() -> bool {
    true
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
    // `FERRO_MANIFEST`, when set, is authoritative — no fallback to the
    // well-known paths. This lets CI explicitly disable the runner via
    // `FERRO_MANIFEST=/nonexistent` even on a host that happens to have
    // one of the well-known paths mounted.
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
            let path = manifest_path()?;
            MultiFastaProvider::from_manifest(&path)
                .map_err(|e| eprintln!("from_manifest({}) failed: {e}", path.display()))
                .ok()
                .map(Arc::new)
        })
        .clone()
}

/// `Arc<MultiFastaProvider>` newtype that implements `ReferenceProvider + Clone`
/// so it can be plugged into [`VariantProjector`], which requires `Clone`.
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

fn normalizer() -> Option<Normalizer<ArcProvider>> {
    Some(Normalizer::new(ArcProvider(provider()?)))
}

fn variant_projector() -> Option<VariantProjector<ArcProvider>> {
    let p = provider()?;
    let cdot = p.cdot_mapper()?.clone();
    let projector = Projector::new(cdot);
    Some(VariantProjector::new(projector, ArcProvider(p)))
}

// ----------------------------------------------------------------------------
// Per-axis tally
// ----------------------------------------------------------------------------

#[derive(Debug)]
struct AxisTally {
    axis: &'static str,
    pass: usize,
    fail: Vec<(String, String)>, // (input, diagnostic)
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
        let matches = matches!(&actual, Ok(s) if s == expected);
        if matches {
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
        // Write FAIL inputs (one per line) to /tmp/ferro-xfail/<axis>.txt and
        // (input \t diagnostic) pairs to /tmp/ferro-xfail/<axis>.tsv so the
        // committed baseline-failures/<axis>.txt and failure-patterns.md can
        // be regenerated from a single run.
        let dir = Path::new(XFAIL_REPORT_DIR);
        let _ = fs::create_dir_all(dir);
        let report_path = dir.join(format!("{}.txt", self.axis));
        let tsv_path = dir.join(format!("{}.tsv", self.axis));
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
            "{}: {} divergence(s) from mutalyzer — see {} and tests/fixtures/mutalyzer-normalize/failure-patterns.md",
            self.axis,
            self.fail.len(),
            report_path.display()
        );
    }
}

// ----------------------------------------------------------------------------
// Helpers
// ----------------------------------------------------------------------------

fn transcript_of(v: &HgvsVariant) -> Option<String> {
    match v {
        HgvsVariant::Cds(c) => Some(c.accession.full()),
        HgvsVariant::Tx(t) => Some(t.accession.full()),
        HgvsVariant::Rna(r) => Some(r.accession.full()),
        _ => None,
    }
}

fn format_pairs(pairs: &[Vec<String>]) -> String {
    let inner: Vec<String> = pairs
        .iter()
        .filter(|p| p.len() == 2)
        .map(|p| format!("({:?}, {:?})", p[0], p[1]))
        .collect();
    format!("[{}]", inner.join(", "))
}

/// Map a mutalyzer error/info code to a substring that should appear in the
/// `Debug` representation of the corresponding ferro-hgvs error. Test-side
/// (not in `src/`) so changes don't ripple through production code. Unmapped
/// codes count as FAIL with a `no mapping for X` diagnostic — extend the
/// table to fix.
fn map_mutalyzer_code(code: &str) -> Option<&'static str> {
    Some(match code {
        "ESEQUENCEMISMATCH" => "SequenceMismatch",
        // Extend as mappings are discovered. Unmapped codes default to FAIL.
        _ => return None,
    })
}

// ----------------------------------------------------------------------------
// Smoke tests
// ----------------------------------------------------------------------------

#[test]
fn loads_fixture() {
    let f = fixture();
    let active = f.cases.iter().filter(|c| c.to_test).count();
    println!(
        "mutalyzer-normalize: loaded {} total cases ({} to_test) @ commit {}",
        f.cases.len(),
        active,
        f.source_commit
    );
    assert!(active > 0, "fixture should have at least one to_test case");
}

#[test]
fn manifest_or_skip() {
    match manifest_path() {
        Some(p) => println!("mutalyzer-normalize: using manifest {}", p.display()),
        None => println!(
            "mutalyzer-normalize: skipping — no manifest at FERRO_MANIFEST, \
             /Volumes/scratch-00001/.../manifest.json, or benchmark-output/manifest.json"
        ),
    }
}

// ----------------------------------------------------------------------------
// Axis: normalized
// ----------------------------------------------------------------------------

#[test]
fn axis_normalized() {
    let Some(normalizer) = normalizer() else {
        eprintln!("axis_normalized: skipping — no manifest");
        return;
    };

    let mut t = AxisTally::new("normalized");
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected) = case.normalized.as_deref() else {
            t.skipped += 1;
            continue;
        };

        let actual = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse: {e}"))?;
            let n = normalizer
                .normalize(&v)
                .map_err(|e| format!("normalize: {e}"))?;
            Ok(format!("{n}"))
        });

        t.record(&case.input, expected, actual);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: genomic
// ----------------------------------------------------------------------------

#[test]
fn axis_genomic() {
    let Some(normalizer) = normalizer() else {
        eprintln!("axis_genomic: skipping — no manifest");
        return;
    };

    let mut t = AxisTally::new("genomic");
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected) = case.genomic.as_deref() else {
            t.skipped += 1;
            continue;
        };

        let actual = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse: {e}"))?;
            let n = normalizer
                .normalize(&v)
                .map_err(|e| format!("normalize: {e}"))?;
            match &n {
                HgvsVariant::Genome(_) => Ok(format!("{n}")),
                _ => Err("c→g projection API not yet exposed".to_string()),
            }
        });

        t.record(&case.input, expected, actual);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: protein_description
// ----------------------------------------------------------------------------

#[test]
fn axis_protein_description() {
    let Some(vp) = variant_projector() else {
        eprintln!("axis_protein_description: skipping — no manifest");
        return;
    };

    let mut t = AxisTally::new("protein_description");
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected) = case.protein_description.as_deref() else {
            t.skipped += 1;
            continue;
        };

        let actual = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse: {e}"))?;
            let tx_id =
                transcript_of(&v).ok_or_else(|| "could not infer transcript_id".to_string())?;
            let result = vp
                .project_variant(&v, &tx_id)
                .map_err(|e| format!("project: {e}"))?;
            result
                .protein
                .as_ref()
                .map(|p| format!("{p}"))
                .ok_or_else(|| "no protein predicted".to_string())
        });

        t.record(&case.input, expected, actual);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: coding_protein_descriptions
// ----------------------------------------------------------------------------

#[test]
fn axis_coding_protein_descriptions() {
    let Some(vp) = variant_projector() else {
        eprintln!("axis_coding_protein_descriptions: skipping — no manifest");
        return;
    };

    let mut t = AxisTally::new("coding_protein_descriptions");
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(pairs) = case.coding_protein_descriptions.as_ref() else {
            t.skipped += 1;
            continue;
        };
        let expected_repr = format_pairs(pairs);

        let actual = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse: {e}"))?;
            let results = vp
                .project_variant_all(&v)
                .map_err(|e| format!("project_all: {e}"))?;
            let actual_pairs: Vec<(String, String)> = results
                .iter()
                .filter_map(|r| {
                    let c = r.coding.as_ref()?.to_string();
                    let p = r.protein.as_ref()?.to_string();
                    Some((c, p))
                })
                .collect();

            for pair in pairs {
                if pair.len() != 2 {
                    continue;
                }
                let want_c = &pair[0];
                let want_p = &pair[1];
                if !actual_pairs.iter().any(|(c, p)| c == want_c && p == want_p) {
                    return Err(format!(
                        "missing pair ({want_c:?}, {want_p:?}); got {actual_pairs:?}"
                    ));
                }
            }
            Ok(expected_repr.clone())
        });

        t.record(&case.input, &expected_repr, actual);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: rna_description
// ----------------------------------------------------------------------------

#[test]
fn axis_rna_description() {
    // Manifest-gate even though the current body doesn't consume the
    // manifest — keeps CI green and matches the other axes' shape.
    if manifest_path().is_none() {
        eprintln!("axis_rna_description: skipping — no manifest");
        return;
    }

    let mut t = AxisTally::new("rna_description");
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected) = case.rna_description.as_deref() else {
            t.skipped += 1;
            continue;
        };
        let actual: Result<String, String> =
            Err("ferro-hgvs r. prediction surface not yet wired into this runner".to_string());
        t.record(&case.input, expected, actual);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: noncoding
// ----------------------------------------------------------------------------

#[test]
fn axis_noncoding() {
    if manifest_path().is_none() {
        eprintln!("axis_noncoding: skipping — no manifest");
        return;
    }

    let mut t = AxisTally::new("noncoding");
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected_list) = case.noncoding.as_ref() else {
            t.skipped += 1;
            continue;
        };
        let expected = expected_list.join(" | ");
        let actual: Result<String, String> =
            Err("ferro-hgvs n. projection surface not yet wired into this runner".to_string());
        t.record(&case.input, &expected, actual);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: errors
// ----------------------------------------------------------------------------

#[test]
fn axis_errors() {
    let Some(normalizer) = normalizer() else {
        eprintln!("axis_errors: skipping — no manifest");
        return;
    };

    let mut t = AxisTally::new("errors");
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected) = case.errors.as_ref() else {
            t.skipped += 1;
            continue;
        };
        let expected_repr = expected.join(",");

        let actual = catch_panics(|| -> Result<String, String> {
            let mapped: Vec<&str> = expected
                .iter()
                .map(|c| map_mutalyzer_code(c).unwrap_or("<unmapped>"))
                .collect();
            if mapped.iter().any(|s| *s == "<unmapped>") {
                return Err(format!("no mapping for one of {expected:?}"));
            }

            let v_result = parse_hgvs(&case.input);
            let normalize_err = match v_result {
                Ok(v) => normalizer.normalize(&v).err().map(|e| format!("{e:?}")),
                Err(e) => Some(format!("{e:?}")),
            };

            match normalize_err {
                None => Err(format!(
                    "ferro produced no error; mutalyzer expected {expected:?}"
                )),
                Some(actual_err) => {
                    for tag in &mapped {
                        if !actual_err.contains(tag) {
                            return Err(format!(
                                "ferro error {actual_err:?} does not contain expected tag {tag:?}"
                            ));
                        }
                    }
                    Ok(expected_repr.clone())
                }
            }
        });

        t.record(&case.input, &expected_repr, actual);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: infos
// ----------------------------------------------------------------------------

#[test]
fn axis_infos() {
    if manifest_path().is_none() {
        eprintln!("axis_infos: skipping — no manifest");
        return;
    }

    let mut t = AxisTally::new("infos");
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected_list) = case.infos.as_ref() else {
            t.skipped += 1;
            continue;
        };
        let expected = expected_list.join(",");
        let actual: Result<String, String> =
            Err("ferro-hgvs info-code surface not yet wired into this runner".to_string());
        t.record(&case.input, &expected, actual);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Regression layer: MockProvider-pinned `normalized` axis (CI-always)
// ----------------------------------------------------------------------------
//
// The axis_* tests above gate on a real reference manifest, so CI cannot
// detect refactor-side regressions in ferro's mock-mode behavior. This test
// fills that gap by pinning ferro's output under `MockProvider::new()` for
// every `normalized`-bearing case, then diffing against the committed
// `mock-pin/normalized.txt` snapshot.
//
// Crucially, the pin records *ferro's current behavior*, not upstream truth.
// A line reading `INPUT\tINPUT` means ferro returned the input verbatim (no
// shift possible without reference bases) — that's the regression baseline,
// not a correctness claim. Correctness is asserted by `axis_normalized`
// against `cases.json`.
//
// The pin format is one line per case: `<input>\t<ferro_output_or_error>`.
// Inputs are emitted in cases.json order for byte-stable diffs.

const MOCK_PIN_PATH: &str = "tests/fixtures/mutalyzer-normalize/mock-pin/normalized.txt";

fn render_mock_normalized() -> String {
    let normalizer = Normalizer::new(MockProvider::new());
    let mut lines: Vec<String> = Vec::new();
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        if case.normalized.is_none() {
            continue;
        }
        let output = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse error: {e}"))?;
            let n = normalizer
                .normalize(&v)
                .map_err(|e| format!("normalize error: {e}"))?;
            Ok(format!("{n}"))
        })
        .unwrap_or_else(|e| e);
        lines.push(format!("{}\t{}", case.input, output));
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

    let expected = fs::read_to_string(MOCK_PIN_PATH)
        .unwrap_or_else(|e| panic!("read {MOCK_PIN_PATH}: {e}; if this is a fresh corpus, run `BLESS_MOCK_PIN=1 cargo nextest run --features dev -E 'test(regression_under_mock_normalized)'`"));

    if observed == expected {
        let total = observed.lines().count();
        println!(
            "regression_under_mock_normalized: {} cases pinned, all match",
            total
        );
        return;
    }

    // Find the first few diverging lines for a tight error message.
    let mut diffs = Vec::new();
    for (i, (o, e)) in observed.lines().zip(expected.lines()).enumerate() {
        if o != e {
            diffs.push(format!("  line {i}:\n    pinned:   {e}\n    observed: {o}"));
            if diffs.len() >= 10 {
                break;
            }
        }
    }
    let observed_total = observed.lines().count();
    let expected_total = expected.lines().count();
    if observed_total != expected_total {
        diffs.push(format!(
            "  line count differs: observed={observed_total}, expected={expected_total} (cases.json may have changed without a re-bless)"
        ));
    }

    panic!(
        "mock pin drifted ({} divergence(s) shown).\n\n{}\n\n\
         If this drift is intentional, regenerate the pin:\n  \
         BLESS_MOCK_PIN=1 cargo nextest run --features dev -E 'test(regression_under_mock_normalized)'",
        diffs.len(),
        diffs.join("\n"),
    );
}
