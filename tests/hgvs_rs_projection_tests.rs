//! hgvs-rs projection test cases.
//!
//! Imports the biocommons hgvs-rs mapper corpus (curated subset, pinned via
//! `scripts/refresh-hgvs-rs-projection-fixtures.py`) and runs it against
//! ferro-hgvs's projector. This is the third consumer of the corpus harness,
//! mirroring `tests/mutalyzer_normalize_tests.rs` (the multi-axis, projection-
//! focused template) rather than the single-axis biocommons normalize corpus.
//!
//! The case `input` is a genome-anchored HGVSg (or a bare HGVSc for the direct
//! `c→p` rows); the expected HGVSc/HGVSp/`n.` values drive the projection axes:
//!
//! - `axis_coding_protein_descriptions` — `project_variant_all(input)`; asserts
//!   each expected `[c, p]` pair is among the `(coding, protein)` results
//!   (`g→c→p`).
//! - `axis_protein_description` — `project_variant(input, transcript_of(input))`
//!   → `result.protein` (direct `c→p` from a bare `c.` input).
//! - `axis_coding` — `project_variant_all(input)`; asserts the expected `c.`
//!   is among the `coding` results (`g→c`-only UTR rows).
//! - `axis_noncoding` — manifest-gated wiring placeholder (NR_ rows).
//!
//! Backed by the ferro-prepared reference manifest. When the manifest is absent
//! (CI by default), every axis test reports `skipping — no manifest` and exits
//! 0. Same convention as `tests/mutalyzer_normalize_tests.rs`.
//!
//! Each axis test writes its current FAIL list to `/tmp/ferro-xfail/<axis>.{txt,
//! tsv}` for burn-down tracking. The `failure-patterns.md` summary is generated
//! from `cases.json` (see `examples/generate_conformance_summary`), never hand-
//! maintained (#509).
//!
//! **Regression layer (CI-always):** `regression_under_mock_<axis>` runs every
//! applicable case through ferro under `MockProvider` and diffs against
//! `mock-pin/<axis>.txt`. Under Mock there is no reference, so most projections
//! error / return None — the pin captures whatever ferro does, asserting only
//! stability. Regenerate after intentional changes:
//!     `BLESS_MOCK_PIN=1 cargo nextest run --features dev -E 'test(regression_under_mock) & binary(hgvs_rs_projection_tests)'`

use ferro_hgvs::conformance::hgvs_rs_projection::{
    AcceptedDivergence, Axis, Case, Fixture, Improvement, KnownBug, Policy, SpecCitation,
    SpecSection,
};
use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::{
    parse_hgvs, FerroError, HgvsVariant, MultiFastaProvider, ReferenceProvider, VariantProjector,
};
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

const FIXTURE_PATH: &str = "tests/fixtures/hgvs-rs-projection/cases.json";
const XFAIL_REPORT_DIR: &str = "/tmp/ferro-xfail";
const FAIL_PRINT_LIMIT: usize = 10;

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
    // `FERRO_MANIFEST=/nonexistent` even on a host that happens to have one of
    // the well-known paths mounted.
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

/// `Arc<MultiFastaProvider>` newtype that implements `ReferenceProvider + Clone`
/// so it can be plugged into [`VariantProjector`], which requires `Clone`.
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

fn variant_projector() -> Option<VariantProjector<ArcProvider>> {
    let p = provider()?;
    let cdot = p.cdot_mapper()?.clone();
    let projector = Projector::new(cdot);
    Some(VariantProjector::new(projector, ArcProvider(p)))
}

// ----------------------------------------------------------------------------
// Per-axis tally (mirrors tests/mutalyzer_normalize_tests.rs)
// ----------------------------------------------------------------------------

#[derive(Debug)]
struct AxisTally {
    axis: Axis,
    pass: usize,
    divergence_accepted: Vec<(String, String)>,
    spec_overridden: Vec<(String, String)>,
    known_bug: Vec<(String, String)>,
    improvement: Vec<(String, String)>,
    fail: Vec<(String, String)>, // (input, diagnostic)
    skipped: usize,
}

impl AxisTally {
    fn new(axis: Axis) -> Self {
        Self {
            axis,
            pass: 0,
            divergence_accepted: Vec::new(),
            spec_overridden: Vec::new(),
            known_bug: Vec::new(),
            improvement: Vec::new(),
            fail: Vec::new(),
            skipped: 0,
        }
    }

    /// Bucket a single case against `expected`/`actual`. See the mutalyzer
    /// harness for the full contract: a match XPASS-fails when a disposition is
    /// stale; a string mismatch on the disposition's axis routes to the
    /// matching bucket; an `Err` (panic / parse / project failure) is always a
    /// real FAIL even when a disposition is present.
    fn record(&mut self, case: &Case, expected: &str, actual: Result<String, String>) {
        let matches = matches!(&actual, Ok(s) if s == expected);
        if matches {
            if let Some(ad) = &case.accepted_divergence {
                if ad.axis == self.axis {
                    self.fail.push((
                        case.input.clone(),
                        format!(
                            "XPASS: accepted_divergence (policy {}) now matches hgvs-rs; the divergence is gone — remove the annotation",
                            ad.policy
                        ),
                    ));
                    return;
                }
            }
            if let Some(kb) = &case.known_bug {
                if kb.axis == self.axis {
                    self.fail.push((
                        case.input.clone(),
                        format!(
                            "XPASS: known_bug #{} now matches hgvs-rs; the fix appears to have landed — remove the annotation and demote the row",
                            kb.tracking_issue
                        ),
                    ));
                    return;
                }
            }
            if let Some(imp) = &case.improvement {
                if imp.axis == self.axis {
                    self.fail.push((
                        case.input.clone(),
                        format!(
                            "XPASS: improvement #{} now matches hgvs-rs; ferro reached the spec-preferred form — remove the annotation and demote the row",
                            imp.tracking_issue
                        ),
                    ));
                    return;
                }
            }
            self.pass += 1;
            return;
        }
        if actual.is_ok() {
            if let Some(ad) = &case.accepted_divergence {
                if ad.axis == self.axis {
                    self.divergence_accepted
                        .push((case.input.clone(), ad.policy.to_string()));
                    return;
                }
            }
            if let Some(kb) = &case.known_bug {
                if kb.axis == self.axis {
                    self.known_bug
                        .push((case.input.clone(), format!("#{}", kb.tracking_issue)));
                    return;
                }
            }
            if let Some(imp) = &case.improvement {
                if imp.axis == self.axis {
                    self.improvement
                        .push((case.input.clone(), format!("#{}", imp.tracking_issue)));
                    return;
                }
            }
            if let Some(sc) = &case.spec_citation {
                if sc.axis == self.axis {
                    self.spec_overridden
                        .push((case.input.clone(), sc.section.to_string()));
                    return;
                }
            }
        }
        let diag = match actual {
            Ok(got) => format!("expected={expected:?} got={got:?}"),
            Err(e) => format!("expected={expected:?} err={e}"),
        };
        self.fail.push((case.input.clone(), diag));
    }

    /// Single-line, grep-friendly summary of the tally state.
    fn summary_line(&self, report_path: &Path) -> String {
        format!(
            "{}: {} pass / {} divergence_accepted / {} known_bug / {} improvement / {} spec_overridden / {} FAIL / {} skipped (FAIL inputs -> {})",
            self.axis,
            self.pass,
            self.divergence_accepted.len(),
            self.known_bug.len(),
            self.improvement.len(),
            self.spec_overridden.len(),
            self.fail.len(),
            self.skipped,
            report_path.display()
        )
    }

    fn finish(self) {
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

        println!("{}", self.summary_line(&report_path));

        for (input, policy) in self.divergence_accepted.iter().take(FAIL_PRINT_LIMIT) {
            eprintln!(
                "  DIVERGENCE_ACCEPTED  [{}] {} | policy={}",
                self.axis, input, policy
            );
        }
        if self.divergence_accepted.len() > FAIL_PRINT_LIMIT {
            eprintln!(
                "  ... {} more divergence_accepted",
                self.divergence_accepted.len() - FAIL_PRINT_LIMIT,
            );
        }

        for (input, issue) in self.known_bug.iter().take(FAIL_PRINT_LIMIT) {
            eprintln!(
                "  KNOWN_BUG  [{}] {} | tracking_issue={}",
                self.axis, input, issue
            );
        }
        if self.known_bug.len() > FAIL_PRINT_LIMIT {
            eprintln!(
                "  ... {} more known_bug",
                self.known_bug.len() - FAIL_PRINT_LIMIT,
            );
        }

        for (input, issue) in self.improvement.iter().take(FAIL_PRINT_LIMIT) {
            eprintln!(
                "  IMPROVEMENT  [{}] {} | tracking_issue={}",
                self.axis, input, issue
            );
        }
        if self.improvement.len() > FAIL_PRINT_LIMIT {
            eprintln!(
                "  ... {} more improvement",
                self.improvement.len() - FAIL_PRINT_LIMIT,
            );
        }

        for (input, section) in self.spec_overridden.iter().take(FAIL_PRINT_LIMIT) {
            eprintln!(
                "  SPEC_OVERRIDDEN  [{}] {} | section={}",
                self.axis, input, section
            );
        }
        if self.spec_overridden.len() > FAIL_PRINT_LIMIT {
            eprintln!(
                "  ... {} more spec_overridden",
                self.spec_overridden.len() - FAIL_PRINT_LIMIT,
            );
        }

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
            "{}: {} divergence(s) from hgvs-rs — see {} and tests/fixtures/hgvs-rs-projection/failure-patterns.md",
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
    // Return the *bare* transcript accession (e.g. `NM_000077.4`), not the
    // compound `genomic_context(transcript)` form `full()` produces.
    match v {
        HgvsVariant::Cds(c) => Some(c.accession.transcript_accession()),
        HgvsVariant::Tx(t) => Some(t.accession.transcript_accession()),
        HgvsVariant::Rna(r) => Some(r.accession.transcript_accession()),
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

// ----------------------------------------------------------------------------
// Smoke tests
// ----------------------------------------------------------------------------

#[test]
fn loads_fixture() {
    let f = fixture();
    let active = f.cases.iter().filter(|c| c.to_test).count();
    println!(
        "hgvs-rs-projection: loaded {} total cases ({} to_test) @ commit {}",
        f.cases.len(),
        active,
        f.source_commit
    );
    assert!(active > 0, "fixture should have at least one to_test case");
}

#[test]
fn manifest_or_skip() {
    match manifest_path() {
        Some(p) => println!("hgvs-rs-projection: using manifest {}", p.display()),
        None => println!(
            "hgvs-rs-projection: skipping — no manifest at FERRO_MANIFEST, \
             /Volumes/scratch-00001/.../manifest.json, or benchmark-output/manifest.json"
        ),
    }
}

// ----------------------------------------------------------------------------
// Axis: coding_protein_descriptions (g→c→p)
// ----------------------------------------------------------------------------

#[test]
fn axis_coding_protein_descriptions() {
    let Some(vp) = variant_projector() else {
        eprintln!("axis_coding_protein_descriptions: skipping — no manifest");
        return;
    };

    let mut t = AxisTally::new(Axis::CodingProteinDescriptions);
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

        t.record(case, &expected_repr, actual);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: protein_description (direct c→p)
// ----------------------------------------------------------------------------

#[test]
fn axis_protein_description() {
    let Some(vp) = variant_projector() else {
        eprintln!("axis_protein_description: skipping — no manifest");
        return;
    };

    let mut t = AxisTally::new(Axis::ProteinDescription);
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

        t.record(case, expected, actual);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: coding (g→c only — UTR rows with no protein)
// ----------------------------------------------------------------------------

#[test]
fn axis_coding() {
    let Some(vp) = variant_projector() else {
        eprintln!("axis_coding: skipping — no manifest");
        return;
    };

    let mut t = AxisTally::new(Axis::Coding);
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected) = case.coding.as_deref() else {
            t.skipped += 1;
            continue;
        };

        let actual = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse: {e}"))?;
            let results = vp
                .project_variant_all(&v)
                .map_err(|e| format!("project_all: {e}"))?;
            let coding: Vec<String> = results
                .iter()
                .filter_map(|r| r.coding.as_ref().map(|c| c.to_string()))
                .collect();
            if coding.iter().any(|c| c == expected) {
                Ok(expected.to_string())
            } else {
                Err(format!("missing coding {expected:?}; got {coding:?}"))
            }
        });

        t.record(case, expected, actual);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: noncoding (g→n — NR_ rows)
// ----------------------------------------------------------------------------

#[test]
fn axis_noncoding() {
    let Some(vp) = variant_projector() else {
        eprintln!("axis_noncoding: skipping — no manifest");
        return;
    };

    let mut t = AxisTally::new(Axis::Noncoding);
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected_list) = case.noncoding.as_ref() else {
            t.skipped += 1;
            continue;
        };
        let expected = expected_list.join(" | ");

        let actual = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse: {e}"))?;
            let results = vp
                .project_variant_all(&v)
                .map_err(|e| format!("project_all: {e}"))?;
            let noncoding: Vec<String> = results
                .iter()
                .filter_map(|r| r.noncoding.as_ref().map(|n| n.to_string()))
                .collect();
            for want in expected_list {
                if !noncoding.iter().any(|n| n == want) {
                    return Err(format!("missing noncoding {want:?}; got {noncoding:?}"));
                }
            }
            Ok(expected.clone())
        });

        t.record(case, &expected, actual);
    }
    t.finish();
}

// ----------------------------------------------------------------------------
// Regression layer: MockProvider-pinned axes (CI-always)
// ----------------------------------------------------------------------------
//
// The axis_* tests above gate on a real reference manifest, so CI cannot detect
// refactor-side regressions in ferro's mock-mode behavior. These tests fill that
// gap by pinning ferro's output under `MockProvider::new()` for every applicable
// case, then diffing against the committed `mock-pin/<axis>.txt` snapshot.
//
// Under Mock there is no reference data, so most projections error / return
// None. The pin records *ferro's current behavior*, not upstream truth — it
// asserts stability, not correctness. Correctness is asserted by the axis_*
// tests against `cases.json` with a manifest.
//
// The pin format is one line per case: `<input>\t<ferro_output_or_error>`.

const MOCK_PIN_DIR: &str = "tests/fixtures/hgvs-rs-projection/mock-pin";

/// Build a Mock-backed projector with an empty cdot mapper. `MockProvider` has
/// no genome→transcript alignment data, so projections error / return None —
/// exactly the "stability, not correctness" baseline the mock pin records.
fn mock_projector() -> VariantProjector<MockProvider> {
    let projector = Projector::new(CdotMapper::new());
    VariantProjector::new(projector, MockProvider::new())
}

/// Render one Mock-pin file for the cases selected by `applies`, projecting via
/// `project`. The output string for a case is whatever ferro produces (a value
/// or an error message), captured verbatim for stability pinning.
fn render_mock_pin(
    applies: impl Fn(&Case) -> bool,
    project: impl Fn(&VariantProjector<MockProvider>, &HgvsVariant) -> Result<String, String>,
) -> String {
    let mut lines: Vec<String> = Vec::new();
    let projector = mock_projector();
    for case in &fixture().cases {
        if !case.to_test || !applies(case) {
            continue;
        }
        let output = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse error: {e}"))?;
            project(&projector, &v)
        })
        .unwrap_or_else(|e| e);
        lines.push(format!("{}\t{}", case.input, output));
    }
    let mut out = lines.join("\n");
    out.push('\n');
    out
}

/// Bless-or-diff a rendered pin against the committed `mock-pin/<axis>.txt`.
fn check_mock_pin(axis: &str, observed: &str) {
    let path = format!("{MOCK_PIN_DIR}/{axis}.txt");

    if std::env::var("BLESS_MOCK_PIN").is_ok() {
        if let Some(parent) = Path::new(&path).parent() {
            let _ = fs::create_dir_all(parent);
        }
        fs::write(&path, observed).unwrap_or_else(|e| panic!("write {path}: {e}"));
        eprintln!("blessed {path} ({} bytes)", observed.len());
        return;
    }

    let expected = fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!(
            "read {path}: {e}; if this is a fresh corpus, run \
             `BLESS_MOCK_PIN=1 cargo nextest run --features dev -E 'test(regression_under_mock) & binary(hgvs_rs_projection_tests)'`"
        )
    });

    if observed == expected {
        let total = observed.lines().count();
        println!("regression_under_mock_{axis}: {total} cases pinned, all match");
        return;
    }

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
        "mock pin {axis} drifted ({} divergence(s) shown).\n\n{}\n\n\
         If this drift is intentional, regenerate the pin:\n  \
         BLESS_MOCK_PIN=1 cargo nextest run --features dev -E 'test(regression_under_mock) & binary(hgvs_rs_projection_tests)'",
        diffs.len(),
        diffs.join("\n"),
    );
}

#[test]
fn regression_under_mock_coding_protein_descriptions() {
    let observed = render_mock_pin(
        |c| c.coding_protein_descriptions.is_some(),
        |vp, v| {
            let results = vp
                .project_variant_all(v)
                .map_err(|e| format!("project_all error: {e}"))?;
            let pairs: Vec<String> = results
                .iter()
                .filter_map(|r| {
                    let c = r.coding.as_ref()?.to_string();
                    let p = r.protein.as_ref()?.to_string();
                    Some(format!("({c:?}, {p:?})"))
                })
                .collect();
            Ok(format!("[{}]", pairs.join(", ")))
        },
    );
    check_mock_pin("coding_protein_descriptions", &observed);
}

#[test]
fn regression_under_mock_protein_description() {
    let observed = render_mock_pin(
        |c| c.protein_description.is_some(),
        |vp, v| {
            let tx_id =
                transcript_of(v).ok_or_else(|| "could not infer transcript_id".to_string())?;
            let result = vp
                .project_variant(v, &tx_id)
                .map_err(|e| format!("project error: {e}"))?;
            result
                .protein
                .as_ref()
                .map(|p| format!("{p}"))
                .ok_or_else(|| "no protein predicted".to_string())
        },
    );
    check_mock_pin("protein_description", &observed);
}

#[test]
fn regression_under_mock_coding() {
    let observed = render_mock_pin(
        |c| c.coding.is_some(),
        |vp, v| {
            let results = vp
                .project_variant_all(v)
                .map_err(|e| format!("project_all error: {e}"))?;
            let coding: Vec<String> = results
                .iter()
                .filter_map(|r| r.coding.as_ref().map(|c| c.to_string()))
                .collect();
            Ok(format!("[{}]", coding.join(", ")))
        },
    );
    check_mock_pin("coding", &observed);
}

// ----------------------------------------------------------------------------
// Comparator unit tests (hermetic — no manifest, no fixture, no I/O)
// ----------------------------------------------------------------------------

#[cfg(test)]
mod comparator_tests {
    use super::*;
    use std::path::PathBuf;

    /// Build a minimal `Case` with the given input. All optional fields default.
    fn make_case(input: &str) -> Case {
        Case {
            keywords: Vec::new(),
            input: input.to_string(),
            coding: None,
            protein_description: None,
            coding_protein_descriptions: None,
            noncoding: None,
            to_test: true,
            accepted_divergence: None,
            known_bug: None,
            improvement: None,
            spec_citation: None,
        }
    }

    // A match with no annotation is a PASS.
    #[test]
    fn tally_pass_when_output_matches() {
        let mut t = AxisTally::new(Axis::Coding);
        let case = make_case("g");
        t.record(&case, "X", Ok("X".to_string()));
        assert_eq!(t.pass, 1);
        assert!(t.fail.is_empty());
    }

    // A mismatch with no annotation is a FAIL.
    #[test]
    fn tally_fail_when_no_annotation() {
        let mut t = AxisTally::new(Axis::Coding);
        let case = make_case("g");
        t.record(&case, "X", Ok("Y".to_string()));
        assert_eq!(t.pass, 0);
        assert_eq!(t.fail.len(), 1);
        assert!(t.fail[0].1.contains("expected=\"X\" got=\"Y\""));
    }

    // A string mismatch on the annotated `coding` axis routes to
    // divergence_accepted, not FAIL.
    #[test]
    fn tally_divergence_accepted_on_coding_axis() {
        let mut t = AxisTally::new(Axis::Coding);
        let mut case = make_case("g");
        case.accepted_divergence = Some(AcceptedDivergence {
            axis: Axis::Coding,
            policy: Policy::GeneSymbolSelector121,
            note: None,
            cluster: None,
        });
        t.record(&case, "X", Ok("Y".to_string()));
        assert_eq!(t.pass, 0);
        assert_eq!(t.divergence_accepted.len(), 1);
        assert!(t.fail.is_empty());
    }

    // An annotation on a different axis does not catch the mismatch.
    #[test]
    fn tally_fail_when_annotation_axis_mismatches() {
        let mut t = AxisTally::new(Axis::Coding);
        let mut case = make_case("g");
        case.accepted_divergence = Some(AcceptedDivergence {
            axis: Axis::ProteinDescription,
            policy: Policy::GeneSymbolSelector121,
            note: None,
            cluster: None,
        });
        t.record(&case, "X", Ok("Y".to_string()));
        assert!(t.divergence_accepted.is_empty());
        assert_eq!(t.fail.len(), 1);
    }

    // An `Err` with an annotation on the matching axis is STILL a FAIL — the
    // annotation accepts a string difference, not a catastrophic failure.
    #[test]
    fn tally_err_with_annotation_still_fails() {
        let mut t = AxisTally::new(Axis::Coding);
        let mut case = make_case("g");
        case.known_bug = Some(KnownBug {
            axis: Axis::Coding,
            tracking_issue: 999,
            note: None,
            cluster: None,
        });
        t.record(&case, "X", Err("project: boom".to_string()));
        assert!(t.known_bug.is_empty());
        assert_eq!(t.fail.len(), 1);
        assert!(t.fail[0].1.contains("err=project: boom"));
    }

    // XPASS on an annotated axis fails loudly so the stale annotation is removed.
    #[test]
    fn tally_known_bug_xpass_fails() {
        let mut t = AxisTally::new(Axis::Coding);
        let mut case = make_case("g");
        case.known_bug = Some(KnownBug {
            axis: Axis::Coding,
            tracking_issue: 999,
            note: None,
            cluster: None,
        });
        t.record(&case, "X", Ok("X".to_string()));
        assert_eq!(t.pass, 0);
        assert_eq!(t.fail.len(), 1);
        assert!(t.fail[0].1.contains("XPASS"));
    }

    // An `improvement` (with required section) routes a mismatch on its axis.
    #[test]
    fn tally_improvement_on_matching_axis() {
        let mut t = AxisTally::new(Axis::ProteinDescription);
        let mut case = make_case("c");
        case.improvement = Some(Improvement {
            axis: Axis::ProteinDescription,
            tracking_issue: 500,
            section: SpecSection::ProteinReference,
            note: None,
            cluster: None,
        });
        t.record(&case, "X", Ok("Y".to_string()));
        assert_eq!(t.improvement, vec![("c".to_string(), "#500".to_string())]);
        assert!(t.fail.is_empty());
    }

    // A `spec_citation` routes a string mismatch on the cited axis.
    #[test]
    fn tally_spec_citation_routes_to_spec_overridden() {
        let mut t = AxisTally::new(Axis::ProteinDescription);
        let mut case = make_case("c");
        case.spec_citation = Some(SpecCitation {
            axis: Axis::ProteinDescription,
            section: SpecSection::ProteinReference,
            note: None,
            cluster: None,
        });
        t.record(&case, "X", Ok("Y".to_string()));
        assert_eq!(t.spec_overridden.len(), 1);
        assert!(t.fail.is_empty());
    }

    // The summary line carries the documented stable format including the new
    // `coding` axis name.
    #[test]
    fn tally_summary_line_format() {
        let mut t = AxisTally::new(Axis::Coding);
        t.record(&make_case("a"), "X", Ok("X".to_string()));
        t.record(&make_case("b"), "X", Ok("Y".to_string()));
        let path = PathBuf::from("/tmp/example.txt");
        assert_eq!(
            t.summary_line(&path),
            "coding: 1 pass / 0 divergence_accepted / 0 known_bug / 0 improvement / 0 spec_overridden / 1 FAIL / 0 skipped \
             (FAIL inputs -> /tmp/example.txt)"
        );
    }
}
