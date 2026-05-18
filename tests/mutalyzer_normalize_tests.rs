//! Mutalyzer normalizer test cases.
//!
//! Imports tests from mutalyzer/mutalyzer's `tests/variants_set.py` (pinned
//! to a specific SHA via `scripts/refresh-mutalyzer-fixtures.py`) and runs
//! them against ferro-hgvs's normalizer + projector. Per-axis `xfail` markers
//! in `cases.json` track capability gaps; XPASS (xfailed axis that starts
//! passing) fails the test so it gets promoted.
//!
//! Backed by the ferro-prepared reference manifest. When the manifest is
//! absent the tests print a `skipping — no manifest` line and exit 0 — this
//! matches the convention in `tests/real_data_normalization_tests.rs`. CI
//! that wants these tests enforced must provision the manifest.
//!
//! Each per-axis test writes its FAIL inputs to
//! `/tmp/ferro-xfail/<axis>.txt` so they can be batch-annotated via:
//!
//!     pixi run python scripts/refresh-mutalyzer-fixtures.py \
//!       annotate-xfails --axis <axis> --from /tmp/ferro-xfail/<axis>.txt

use ferro_hgvs::data::projection::Projector;
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
    #[serde(default)]
    xfail: Vec<String>,
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
    if let Ok(path) = std::env::var("FERRO_MANIFEST") {
        let p = PathBuf::from(path);
        if p.exists() {
            return Some(p);
        }
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
    xfail: usize,
    fail: Vec<(String, String)>, // (input, diagnostic)
    xpass: Vec<String>,          // inputs that XPASSed
    skipped: usize,
}

impl AxisTally {
    fn new(axis: &'static str) -> Self {
        Self {
            axis,
            pass: 0,
            xfail: 0,
            fail: Vec::new(),
            xpass: Vec::new(),
            skipped: 0,
        }
    }

    fn record(
        &mut self,
        case_input: &str,
        expected: &str,
        actual: Result<String, String>,
        xfailed: bool,
    ) {
        let matches = matches!(&actual, Ok(s) if s == expected);
        match (matches, xfailed) {
            (true, false) => self.pass += 1,
            (false, false) => {
                let diag = match actual {
                    Ok(got) => format!("expected={expected:?} got={got:?}"),
                    Err(e) => format!("expected={expected:?} err={e}"),
                };
                self.fail.push((case_input.to_string(), diag));
            }
            (false, true) => self.xfail += 1,
            (true, true) => self.xpass.push(case_input.to_string()),
        }
    }

    fn finish(self) {
        // Write FAIL inputs to /tmp/ferro-xfail/<axis>.txt for annotation.
        let dir = Path::new(XFAIL_REPORT_DIR);
        let _ = fs::create_dir_all(dir);
        let report_path = dir.join(format!("{}.txt", self.axis));
        let body: String = self.fail.iter().map(|(i, _)| format!("{i}\n")).collect();
        let _ = fs::write(&report_path, &body);

        println!(
            "{}: {} pass / {} xfail / {} FAIL / {} XPASS / {} skipped (FAIL inputs -> {})",
            self.axis,
            self.pass,
            self.xfail,
            self.fail.len(),
            self.xpass.len(),
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
        for input in &self.xpass {
            eprintln!(
                "  XPASS [{}] {} | promote (remove from xfail in cases.json)",
                self.axis, input
            );
        }

        assert!(
            self.fail.is_empty() && self.xpass.is_empty(),
            "{}: {} unexpected FAIL + {} XPASS",
            self.axis,
            self.fail.len(),
            self.xpass.len()
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

fn axis_xfailed(case: &Case, axis: &str) -> bool {
    case.xfail.iter().any(|a| a == axis || a == "*")
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
/// `Debug` representation of the corresponding ferro-hgvs error. Lives here
/// (test-side) so changes don't ripple through production code. Unmapped
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
        let xfailed = axis_xfailed(case, "normalized");

        let actual = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse: {e}"))?;
            let n = normalizer
                .normalize(&v)
                .map_err(|e| format!("normalize: {e}"))?;
            Ok(format!("{n}"))
        });

        t.record(&case.input, expected, actual, xfailed);
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
        let xfailed = axis_xfailed(case, "genomic");

        // ferro-hgvs has no exposed c./n. → g. API right now (would require
        // a `VariantProjector::coding_to_genomic` we don't ship yet). For now
        // we only assert when normalization yields a g. variant directly.
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

        t.record(&case.input, expected, actual, xfailed);
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
        let xfailed = axis_xfailed(case, "protein_description");

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

        t.record(&case.input, expected, actual, xfailed);
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
        let xfailed = axis_xfailed(case, "coding_protein_descriptions");
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

        t.record(&case.input, &expected_repr, actual, xfailed);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: rna_description
// ----------------------------------------------------------------------------

#[test]
fn axis_rna_description() {
    let mut t = AxisTally::new("rna_description");
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected) = case.rna_description.as_deref() else {
            t.skipped += 1;
            continue;
        };
        let xfailed = axis_xfailed(case, "rna_description");
        let actual: Result<String, String> =
            Err("ferro-hgvs r. prediction surface not yet wired into this runner".to_string());
        t.record(&case.input, expected, actual, xfailed);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: noncoding
// ----------------------------------------------------------------------------

#[test]
fn axis_noncoding() {
    let mut t = AxisTally::new("noncoding");
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected_list) = case.noncoding.as_ref() else {
            t.skipped += 1;
            continue;
        };
        let xfailed = axis_xfailed(case, "noncoding");
        let expected = expected_list.join(" | ");
        let actual: Result<String, String> =
            Err("ferro-hgvs n. projection surface not yet wired into this runner".to_string());
        t.record(&case.input, &expected, actual, xfailed);
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
        let xfailed = axis_xfailed(case, "errors");
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

        t.record(&case.input, &expected_repr, actual, xfailed);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: infos
// ----------------------------------------------------------------------------

#[test]
fn axis_infos() {
    let mut t = AxisTally::new("infos");
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected_list) = case.infos.as_ref() else {
            t.skipped += 1;
            continue;
        };
        let xfailed = axis_xfailed(case, "infos");
        let expected = expected_list.join(",");
        let actual: Result<String, String> =
            Err("ferro-hgvs info-code surface not yet wired into this runner".to_string());
        t.record(&case.input, &expected, actual, xfailed);
    }
    t.finish();
}
