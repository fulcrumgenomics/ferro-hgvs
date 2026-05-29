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
//! The per-axis `axis_*` tests are pure pass/fail gates and do NOT touch the
//! filesystem. To (re)generate the FAIL lists under `/tmp/ferro-xfail/<axis>.{txt,tsv}`
//! — used to seed the committed `baseline-failures/` snapshot and
//! `failure-patterns.md` — run the single-process `regenerate_baselines`
//! generator (it is `#[ignore]`d in the normal suite):
//!
//!     FERRO_MANIFEST=/path/to/manifest.json \
//!       cargo nextest run --features dev -E 'test(regenerate_baselines)' \
//!       --run-ignored only --no-capture
//!
//! Running it in one process (rather than relying on best-effort writes from
//! the parallel, panicking per-axis test processes) avoids intermittently
//! dropped artifact files.
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

/// Policy label recorded for `infos`-axis accepted divergences. The upstream
/// mutalyzer code is one ferro intentionally does not model
/// (`error_handling::info_map::NO_FERRO_INFO_EQUIV`): a mutalyzer *internal*
/// diagnostic (e.g. `ISORTEDVARIANTS`, `ICORRECTEDCOORDINATESYSTEM`) that is
/// not part of the HGVS nomenclature spec, so ferro is under no obligation to
/// mirror it. See #326 (spec arbitration).
const INFO_NO_SPEC_EQUIVALENT: &str = "mutalyzer-info-no-spec-equivalent";

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

#[derive(Debug, Deserialize, Clone)]
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
    /// Marks a mismatch on a specific axis as an accepted (non-bug)
    /// divergence — e.g. a ferro policy decision that intentionally
    /// disagrees with mutalyzer while both outputs remain HGVS-spec-allowed.
    /// Tallied into `AxisTally::divergence_accepted` instead of `fail`.
    #[serde(default)]
    accepted_divergence: Option<AcceptedDivergence>,
    /// Documentary annotation for cases where mutalyzer's expected output
    /// has been corrected in `cases.json` to ferro's spec-correct value.
    /// Has NO effect on tally bucketing — the corrected expected string
    /// makes the case PASS by the usual match check; this field only
    /// records the spec citation for review.
    #[serde(default)]
    spec_citation: Option<SpecCitation>,
}

/// Closed enum of corpus-runner axes. Replaces the previous free-form
/// `String` field on `AcceptedDivergence` and `SpecCitation` — serde
/// deserialization rejects unknown variants, so a typo in `cases.json`
/// surfaces as a hard parse error rather than silently no-op'ing the
/// annotation (#397 item 2). The variants match the 8 `AxisTally::new`
/// call sites; the lowercase serde rename keeps the JSON wire format
/// stable.
#[derive(Debug, Deserialize, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
#[allow(dead_code)]
enum Axis {
    Normalized,
    Genomic,
    ProteinDescription,
    CodingProteinDescriptions,
    RnaDescription,
    Noncoding,
    Errors,
    Infos,
}

impl Axis {
    /// Stable kebab-case string used in summary lines, file paths, and
    /// log messages. Matches the lowercase strings the previous
    /// `&'static str` axis carried.
    fn as_str(self) -> &'static str {
        match self {
            Axis::Normalized => "normalized",
            Axis::Genomic => "genomic",
            Axis::ProteinDescription => "protein_description",
            Axis::CodingProteinDescriptions => "coding_protein_descriptions",
            Axis::RnaDescription => "rna_description",
            Axis::Noncoding => "noncoding",
            Axis::Errors => "errors",
            Axis::Infos => "infos",
        }
    }
}

impl std::fmt::Display for Axis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Closed enum of accepted-divergence policies. Each policy identifier
/// names a specific reason a mismatch is non-failing (e.g. ferro's spec
/// arbitration on gene-symbol selectors). Adding a new policy requires
/// a code change here, which is intentional: it forces deliberate
/// review of every accepted divergence rather than allowing ad-hoc
/// string typing in `cases.json` (#397 item 2).
#[derive(Debug, Deserialize, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
enum Policy {
    /// ferro's gene-symbol-selector handling per issue #121 (the
    /// `NM_X(GENE):c.X` rendering and the related dispatcher
    /// arbitration). See PR #146 for the spec arbitration trail.
    #[serde(rename = "ferro-policy-121-gene-symbol-selector")]
    GeneSymbolSelector121,
}

impl Policy {
    /// Stable identifier used in summary lines so reviewers can grep
    /// across runs.
    fn as_str(self) -> &'static str {
        match self {
            Policy::GeneSymbolSelector121 => "ferro-policy-121-gene-symbol-selector",
        }
    }
}

impl std::fmt::Display for Policy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Accepted-divergence annotation. When attached to a `Case`, the
/// corpus runner treats a mismatch on `axis` as a tracked-but-non-failing
/// divergence (counted in `divergence_accepted`) rather than a hard FAIL.
///
/// `axis` and `policy` are now closed enums — `serde` rejects unknown
/// variants at fixture-parse time, so misspellings surface as a hard
/// parse error rather than silently no-op'ing the annotation.
#[derive(Debug, Deserialize, Clone)]
#[allow(dead_code)]
struct AcceptedDivergence {
    /// Axis this annotation applies to.
    axis: Axis,
    /// Closed policy identifier explaining *why* the divergence is
    /// acceptable. Surfaced in the per-axis tally summary so reviewers
    /// can grep for the policy.
    policy: Policy,
    /// Optional human-readable note expanding on the policy.
    #[serde(default)]
    note: Option<String>,
}

/// Closed enum of HGVS spec section identifiers cited from `cases.json`.
/// Mirrors [`Policy`]'s shape: adding a new section requires a code
/// change here, which forces deliberate review of every new citation
/// rather than allowing ad-hoc string typing in `cases.json` (#430,
/// follow-up to #397 item 2). The `#[serde(rename = "...")]` annotation
/// keeps the JSON wire format stable.
#[derive(Debug, Deserialize, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
enum SpecSection {
    /// HGVS §Prioritization — the spec's edit-prioritization rules (dup
    /// over ins, etc.) that ferro applies during normalization. This is
    /// a ferro-internal label, not a stable spec heading: the upstream
    /// spec submodule (`assets/hgvs-nomenclature/docs/recommendations/`)
    /// has no `Prioritization` anchor today, so binding to a path here
    /// would create a brittle pointer. The human label is the stable
    /// catalog key; the rename string preserves the byte sequence in
    /// `cases.json`.
    #[serde(rename = "HGVS §Prioritization")]
    HgvsPrioritization,
}

impl SpecSection {
    /// Stable identifier used in summary lines so reviewers can grep
    /// across runs. Matches the `#[serde(rename)]` strings byte-for-byte
    /// to preserve the pre-enum summary-line format.
    fn as_str(self) -> &'static str {
        match self {
            SpecSection::HgvsPrioritization => "HGVS §Prioritization",
        }
    }
}

impl std::fmt::Display for SpecSection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Spec-citation annotation. Documentary only — does NOT affect tally
/// bucketing. Used to record that the `cases.json` expected output for
/// `axis` has been corrected to ferro's spec-correct value (versus
/// mutalyzer's spec-incorrect output), with a citation to the relevant
/// HGVS spec section.
///
/// `axis` and `section` are now closed enums (see [`Axis`], [`SpecSection`]);
/// serde rejects unknown variants at fixture-parse time, so a misspelled
/// section identifier surfaces as a hard parse error rather than silently
/// losing cataloguability.
#[derive(Debug, Deserialize, Clone)]
#[allow(dead_code)]
struct SpecCitation {
    /// Axis this citation applies to.
    axis: Axis,
    /// HGVS spec section the citation references (e.g.
    /// `SpecSection::HgvsPrioritization`).
    section: SpecSection,
    /// Optional human-readable note expanding on the citation.
    #[serde(default)]
    note: Option<String>,
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
    axis: Axis,
    pass: usize,
    /// Mismatches gated by an `accepted_divergence` annotation whose
    /// `axis` matches this tally's axis. Tracked-but-non-failing; reported
    /// in the summary line as `divergence_accepted` and NOT written to
    /// `baseline-failures/<axis>.txt`. Stores `(input, policy)`.
    divergence_accepted: Vec<(String, String)>,
    /// Mismatches routed to the spec-overridden bucket via a
    /// `spec_citation` annotation whose `axis` matches this tally's axis.
    /// Tracked-but-non-failing; reported in the summary line as
    /// `spec_overridden` and NOT written to `baseline-failures/<axis>.txt`.
    /// Stores `(input, section)`. The `cases.json` `expected` field is
    /// expected to carry ferro's spec-correct value, but for axes where
    /// the upstream mutalyzer expectation still differs, the citation
    /// keeps the case visible without surfacing as a FAIL.
    spec_overridden: Vec<(String, String)>,
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
            fail: Vec::new(),
            skipped: 0,
        }
    }

    /// Bucket a single case against `expected`/`actual`:
    /// 1. PASS when `actual == expected`.
    /// 2. `divergence_accepted` when the case carries
    ///    `accepted_divergence` matching this tally's axis.
    /// 3. `spec_overridden` when the case carries a `spec_citation`
    ///    matching this tally's axis.
    /// 4. FAIL otherwise.
    ///
    /// `accepted_divergence` and `spec_citation` only catch a *string
    /// mismatch* from a successful run. An `Err` means ferro failed to
    /// produce any result (panic, parse error, normalize failure) —
    /// those are real bugs and must surface as FAIL even when an
    /// annotation is present.
    fn record(&mut self, case: &Case, expected: &str, actual: Result<String, String>) {
        let matches = matches!(&actual, Ok(s) if s == expected);
        if matches {
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

    /// Single-line, grep-friendly summary of the tally state. When
    /// `report_path` is `Some`, the FAIL-list path is appended (used by the
    /// regenerator, which writes that file); the gate path passes `None`.
    fn summary_line(&self, report_path: Option<&Path>) -> String {
        let dest = match report_path {
            Some(p) => format!(" (FAIL inputs -> {})", p.display()),
            None => String::new(),
        };
        format!(
            "{}: {} pass / {} divergence_accepted / {} spec_overridden / {} FAIL / {} skipped{}",
            self.axis,
            self.pass,
            self.divergence_accepted.len(),
            self.spec_overridden.len(),
            self.fail.len(),
            self.skipped,
            dest
        )
    }

    /// Print the per-bucket diagnostics (accepted / spec-overridden / fail),
    /// truncated to `FAIL_PRINT_LIMIT` per bucket. Shared by the gate
    /// (`assert_clean`) and the regenerator (`write_artifacts`).
    fn print_diagnostics(&self) {
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
                "  ... {} more FAIL (full list via `regenerate_baselines`)",
                self.fail.len() - FAIL_PRINT_LIMIT,
            );
        }
    }

    /// Gate path (used by every `axis_*` test): print a summary + diagnostics
    /// and assert the axis has zero divergences. Performs **no** filesystem
    /// I/O — baseline / `.tsv` artifacts are written only by
    /// [`regenerate_baselines`], in a single process, so we never depend on
    /// best-effort writes from parallel, panicking per-axis test processes.
    fn assert_clean(self) {
        println!("{}", self.summary_line(None));
        self.print_diagnostics();
        assert!(
            self.fail.is_empty(),
            "{}: {} divergence(s) from mutalyzer — run `regenerate_baselines` (see module docs) and consult tests/fixtures/mutalyzer-normalize/failure-patterns.md",
            self.axis,
            self.fail.len(),
        );
    }

    /// Regeneration path (used only by [`regenerate_baselines`]): write the
    /// FAIL inputs to `/tmp/ferro-xfail/<axis>.txt` (one per line) and
    /// `(input \t diagnostic)` pairs to `<axis>.tsv`, then print the summary +
    /// diagnostics. Never asserts. Because the regenerator runs all axes in a
    /// single process, these writes are race-free.
    fn write_artifacts(&self) {
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
        println!("{}", self.summary_line(Some(report_path.as_path())));
        self.print_diagnostics();
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
/// `Debug` representation of the corresponding ferro-hgvs error. Delegates to
/// `ferro_hgvs::error_handling::mutalyzer_to_ferro` so the canonical mapping
/// lives next to the production error taxonomy. Unmapped codes (returned as
/// `None`) count as FAIL with a `no mapping for X` diagnostic.
fn map_mutalyzer_code(code: &str) -> Option<&'static str> {
    ferro_hgvs::error_handling::mutalyzer_to_ferro(code).map(|t| t.debug_tag())
}

/// Map a mutalyzer info code (`I*`) to the `NormalizationInfo::code()`
/// string emitted by ferro for the equivalent signal. Delegates to
/// `ferro_hgvs::error_handling::mutalyzer_info_to_ferro`. Unmapped codes
/// (`None`) name signals ferro intentionally does not emit (see
/// `info_map::NO_FERRO_INFO_EQUIV`).
fn map_mutalyzer_info_code(code: &str) -> Option<&'static str> {
    ferro_hgvs::error_handling::mutalyzer_info_to_ferro(code).map(|t| t.code())
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

fn run_normalized() -> Option<AxisTally> {
    let normalizer = normalizer()?;
    let mut t = AxisTally::new(Axis::Normalized);
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

        t.record(case, expected, actual);
    }
    Some(t)
}

#[test]
fn axis_normalized() {
    match run_normalized() {
        Some(t) => t.assert_clean(),
        None => eprintln!("axis_normalized: skipping — no manifest"),
    }
}

// ----------------------------------------------------------------------------
// Axis: genomic
// ----------------------------------------------------------------------------

fn run_genomic() -> Option<AxisTally> {
    let normalizer = normalizer()?;
    let mut t = AxisTally::new(Axis::Genomic);
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

        t.record(case, expected, actual);
    }
    Some(t)
}

#[test]
fn axis_genomic() {
    match run_genomic() {
        Some(t) => t.assert_clean(),
        None => eprintln!("axis_genomic: skipping — no manifest"),
    }
}

// ----------------------------------------------------------------------------
// Axis: protein_description
// ----------------------------------------------------------------------------

fn run_protein_description() -> Option<AxisTally> {
    let vp = variant_projector()?;
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
    Some(t)
}

#[test]
fn axis_protein_description() {
    match run_protein_description() {
        Some(t) => t.assert_clean(),
        None => eprintln!("axis_protein_description: skipping — no manifest"),
    }
}

// ----------------------------------------------------------------------------
// Axis: coding_protein_descriptions
// ----------------------------------------------------------------------------

fn run_coding_protein_descriptions() -> Option<AxisTally> {
    let vp = variant_projector()?;
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
    Some(t)
}

#[test]
fn axis_coding_protein_descriptions() {
    match run_coding_protein_descriptions() {
        Some(t) => t.assert_clean(),
        None => eprintln!("axis_coding_protein_descriptions: skipping — no manifest"),
    }
}

// ----------------------------------------------------------------------------
// Axis: rna_description
// ----------------------------------------------------------------------------

fn run_rna_description() -> Option<AxisTally> {
    // Manifest-gate even though the current body doesn't consume the
    // manifest — keeps CI green and matches the other axes' shape.
    manifest_path()?;

    let mut t = AxisTally::new(Axis::RnaDescription);
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
        t.record(case, expected, actual);
    }
    Some(t)
}

#[test]
fn axis_rna_description() {
    match run_rna_description() {
        Some(t) => t.assert_clean(),
        None => eprintln!("axis_rna_description: skipping — no manifest"),
    }
}

// ----------------------------------------------------------------------------
// Axis: noncoding
// ----------------------------------------------------------------------------

fn run_noncoding() -> Option<AxisTally> {
    manifest_path()?;

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
        let actual: Result<String, String> =
            Err("ferro-hgvs n. projection surface not yet wired into this runner".to_string());
        t.record(case, &expected, actual);
    }
    Some(t)
}

#[test]
fn axis_noncoding() {
    match run_noncoding() {
        Some(t) => t.assert_clean(),
        None => eprintln!("axis_noncoding: skipping — no manifest"),
    }
}

// ----------------------------------------------------------------------------
// Axis: errors
// ----------------------------------------------------------------------------

fn run_errors() -> Option<AxisTally> {
    let normalizer = normalizer()?;
    let mut t = AxisTally::new(Axis::Errors);
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
            if mapped.contains(&"<unmapped>") {
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

        t.record(case, &expected_repr, actual);
    }
    Some(t)
}

#[test]
fn axis_errors() {
    match run_errors() {
        Some(t) => t.assert_clean(),
        None => eprintln!("axis_errors: skipping — no manifest"),
    }
}

// ----------------------------------------------------------------------------
// Axis: infos
// ----------------------------------------------------------------------------

fn run_infos() -> Option<AxisTally> {
    let normalizer = normalizer()?;
    let mut t = AxisTally::new(Axis::Infos);
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected_list) = case.infos.as_ref() else {
            t.skipped += 1;
            continue;
        };
        let expected_repr = expected_list.join(",");

        // Partition the expected mutalyzer codes: those ferro models (a
        // concrete `NormalizationInfo` tag) vs those ferro intentionally does
        // NOT model (`info_map::NO_FERRO_INFO_EQUIV`). The latter are mutalyzer
        // *internal* diagnostics (e.g. `ISORTEDVARIANTS`, `ICORRECTED*`) that
        // are absent from the HGVS nomenclature spec, so an expected code with
        // no ferro equivalent is an accepted divergence, not a failure (see
        // #326). A *mapped* code ferro fails to emit remains a hard FAIL.
        let has_no_equiv = expected_list
            .iter()
            .any(|c| map_mutalyzer_info_code(c).is_none());

        let actual = catch_panics(|| -> Result<String, String> {
            // Only the modelled codes are required; drop the no-equivalent
            // ones from the comparison entirely.
            let mapped: Vec<&str> = expected_list
                .iter()
                .filter_map(|c| map_mutalyzer_info_code(c))
                .collect();

            let v = parse_hgvs(&case.input).map_err(|e| format!("parse error: {e:?}"))?;
            let result = normalizer
                .normalize_with_diagnostics(&v)
                .map_err(|e| format!("normalize error: {e:?}"))?;
            let emitted: Vec<&str> = result.infos.iter().map(|i| i.code()).collect();

            // Compare per-code multiplicities over the modelled codes. ferro
            // must emit every modelled code the upstream expected (under-
            // emission = real bug) and must not emit extras the upstream did
            // not (over-emission = divergence). No-equivalent expected codes
            // are excluded above, so they never require a ferro emission; an
            // empty `mapped` therefore requires ferro to emit nothing, which
            // the over-emission check below enforces.
            let mut expected_counts = std::collections::HashMap::<&str, usize>::new();
            for tag in &mapped {
                *expected_counts.entry(*tag).or_insert(0) += 1;
            }
            let mut emitted_counts = std::collections::HashMap::<&str, usize>::new();
            for tag in &emitted {
                *emitted_counts.entry(*tag).or_insert(0) += 1;
            }
            for (tag, need) in &expected_counts {
                let got = emitted_counts.get(tag).copied().unwrap_or(0);
                if got < *need {
                    return Err(format!(
                        "ferro infos {emitted:?} missing expected code {tag:?} x{need} (got {got})"
                    ));
                }
            }
            // Symmetric over-emission check: any ferro info beyond what the
            // upstream corpus expected is a divergence (catches `[A, A]` vs
            // `[A]` and unknown extras like `[A, B]` vs `[A]`).
            for (tag, got) in &emitted_counts {
                let need = expected_counts.get(tag).copied().unwrap_or(0);
                if *got > need {
                    return Err(format!(
                        "ferro emitted extra info {tag:?} x{got} (expected x{need}) — full emit list {emitted:?}, expected codes {expected_list:?}"
                    ));
                }
            }
            Ok(expected_repr.clone())
        });

        // A clean match where the upstream expected at least one code ferro
        // intentionally does not model is an accepted divergence (ferro does
        // not mirror mutalyzer's non-spec internal info codes), tallied
        // separately so it stays visible. Errors (under/over-emission, parse,
        // normalize, panic) are real failures routed through `record`.
        if actual.is_ok() && has_no_equiv {
            t.divergence_accepted
                .push((case.input.clone(), INFO_NO_SPEC_EQUIVALENT.to_string()));
        } else {
            t.record(case, &expected_repr, actual);
        }
    }
    Some(t)
}

#[test]
fn axis_infos() {
    match run_infos() {
        Some(t) => t.assert_clean(),
        None => eprintln!("axis_infos: skipping — no manifest"),
    }
}

// ----------------------------------------------------------------------------
// Baseline regeneration (single-process; run explicitly)
// ----------------------------------------------------------------------------

/// Regenerate the per-axis FAIL artifacts under `/tmp/ferro-xfail/` in a
/// **single process**, so the writes are race-free. The per-axis `axis_*`
/// gate tests never touch the filesystem (see `AxisTally::assert_clean`); they
/// run as separate, parallel, panic-on-divergence processes, and relying on
/// best-effort writes from those occasionally dropped a file. This test is the
/// one supported way to (re)generate the `.txt`/`.tsv` snapshots that seed
/// `tests/fixtures/mutalyzer-normalize/baseline-failures/<axis>.txt` and
/// `failure-patterns.md`.
///
/// `#[ignore]` so it does not run in the normal suite (it would always
/// "fail-as-skip" without a manifest and is a generator, not a gate). Run it:
///
///     FERRO_MANIFEST=/path/to/manifest.json \
///       cargo nextest run --features dev -E 'test(regenerate_baselines)' \
///       --run-ignored only --no-capture
fn run_all_axes() -> Vec<AxisTally> {
    [
        run_normalized(),
        run_genomic(),
        run_protein_description(),
        run_coding_protein_descriptions(),
        run_rna_description(),
        run_noncoding(),
        run_errors(),
        run_infos(),
    ]
    .into_iter()
    .flatten()
    .collect()
}

#[test]
#[ignore = "generator, not a gate: run explicitly to (re)write /tmp/ferro-xfail artifacts"]
fn regenerate_baselines() {
    assert!(
        manifest_path().is_some(),
        "regenerate_baselines requires a manifest (set FERRO_MANIFEST or provide a well-known one)"
    );
    let tallies = run_all_axes();
    assert!(
        !tallies.is_empty(),
        "manifest present but no axis produced a tally"
    );
    for t in &tallies {
        t.write_artifacts();
    }
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

// ----------------------------------------------------------------------------
// Comparator unit tests
// ----------------------------------------------------------------------------
//
// Hermetic — no manifest, no fixture, no I/O. Each test builds a synthetic
// `Case` and drives `AxisTally::record` directly to assert the bucketing
// rules for `accepted_divergence` and `spec_citation`.

#[cfg(test)]
mod comparator_tests {
    use super::*;
    use std::path::PathBuf;

    /// Build a minimal `Case` with the given input and optional
    /// annotations. All other fields default to None/empty/true.
    fn make_case(
        input: &str,
        accepted_divergence: Option<AcceptedDivergence>,
        spec_citation: Option<SpecCitation>,
    ) -> Case {
        Case {
            keywords: Vec::new(),
            input: input.to_string(),
            normalized: None,
            genomic: None,
            coding_protein_descriptions: None,
            protein_description: None,
            rna_description: None,
            errors: None,
            infos: None,
            noncoding: None,
            to_test: true,
            accepted_divergence,
            spec_citation,
        }
    }

    fn parse_case(json: &str) -> Case {
        serde_json::from_str(json).expect("Case should deserialize")
    }

    // (1) Backward compat: a case without the new optional fields parses.
    #[test]
    fn case_parses_without_new_fields() {
        let case = parse_case(
            r#"{"input": "NM_000088.3:c.459A>G", "normalized": "NM_000088.3:c.459A>G"}"#,
        );
        assert_eq!(case.input, "NM_000088.3:c.459A>G");
        assert!(case.accepted_divergence.is_none());
        assert!(case.spec_citation.is_none());
    }

    // (1b) Typo in `axis` or `policy` is rejected at deserialization
    // (#397 item 2): the closed enums turn what used to be a silent
    // no-op into a hard parse error during fixture load.
    #[test]
    fn case_axis_typo_rejected_at_deserialization() {
        let result: Result<Case, _> = serde_json::from_str(
            r#"{
                "input": "x",
                "accepted_divergence": {
                    "axis": "normalised",
                    "policy": "ferro-policy-121-gene-symbol-selector"
                }
            }"#,
        );
        assert!(
            result.is_err(),
            "expected typo'd axis (`normalised` vs `normalized`) to be rejected; got {:?}",
            result.ok(),
        );
    }

    #[test]
    fn case_policy_typo_rejected_at_deserialization() {
        let result: Result<Case, _> = serde_json::from_str(
            r#"{
                "input": "x",
                "accepted_divergence": {
                    "axis": "normalized",
                    "policy": "ferro-policy-999-bogus"
                }
            }"#,
        );
        assert!(
            result.is_err(),
            "expected unknown policy to be rejected; got {:?}",
            result.ok(),
        );
    }

    /// Issue #430 closes `SpecCitation.section` to a [`SpecSection`] enum so
    /// a misspelled section identifier surfaces as a hard parse error
    /// rather than silently losing cataloguability in `cases.json`. Mirrors
    /// `case_axis_typo_rejected_at_deserialization` and
    /// `case_policy_typo_rejected_at_deserialization`.
    #[test]
    fn case_section_typo_rejected_at_deserialization() {
        // Missing the `§` sigil — looks plausible but is not the
        // canonical identifier.
        let result: Result<Case, _> = serde_json::from_str(
            r#"{
                "input": "x",
                "spec_citation": {
                    "axis": "normalized",
                    "section": "HGVS Prioritization"
                }
            }"#,
        );
        assert!(
            result.is_err(),
            "expected typo'd section (`HGVS Prioritization` vs `HGVS §Prioritization`) to be rejected; got {:?}",
            result.ok(),
        );
    }

    // (2) `accepted_divergence` deserializes including the optional `note`.
    #[test]
    fn case_parses_with_accepted_divergence() {
        let case = parse_case(
            r#"{
                "input": "NM_000088.3:c.459A>G",
                "accepted_divergence": {
                    "axis": "normalized",
                    "policy": "ferro-policy-121-gene-symbol-selector",
                    "note": "ferro emits (COL1A1) per #121"
                }
            }"#,
        );
        let ad = case
            .accepted_divergence
            .as_ref()
            .expect("accepted_divergence present");
        assert_eq!(ad.axis, Axis::Normalized);
        assert_eq!(ad.policy, Policy::GeneSymbolSelector121);
        assert_eq!(ad.note.as_deref(), Some("ferro emits (COL1A1) per #121"));
    }

    // (3) `spec_citation` deserializes including the optional `note`.
    #[test]
    fn case_parses_with_spec_citation() {
        let case = parse_case(
            r#"{
                "input": "NM_000143.3:c.1_2insCAT",
                "spec_citation": {
                    "axis": "normalized",
                    "section": "HGVS §Prioritization",
                    "note": "dup over ins"
                }
            }"#,
        );
        let sc = case.spec_citation.as_ref().expect("spec_citation present");
        assert_eq!(sc.axis, Axis::Normalized);
        assert_eq!(sc.section, SpecSection::HgvsPrioritization);
        assert_eq!(sc.note.as_deref(), Some("dup over ins"));
    }

    // (4) `record` PASSes when actual == expected, regardless of annotations.
    #[test]
    fn tally_pass_when_output_matches() {
        let mut t = AxisTally::new(Axis::Normalized);
        let case = make_case("NM_000088.3:c.459A>G", None, None);
        t.record(&case, "X", Ok("X".to_string()));
        assert_eq!(t.pass, 1);
        assert!(t.divergence_accepted.is_empty());
        assert!(t.fail.is_empty());
    }

    // (5) Mismatch on the annotated axis goes into `divergence_accepted`.
    #[test]
    fn tally_divergence_accepted_on_matching_axis() {
        let mut t = AxisTally::new(Axis::Normalized);
        let case = make_case(
            "in",
            Some(AcceptedDivergence {
                axis: Axis::Normalized,
                policy: Policy::GeneSymbolSelector121,
                note: None,
            }),
            None,
        );
        t.record(&case, "X", Ok("Y".to_string()));
        assert_eq!(t.pass, 0);
        assert_eq!(
            t.divergence_accepted,
            vec![(
                "in".to_string(),
                "ferro-policy-121-gene-symbol-selector".to_string()
            )]
        );
        assert!(t.fail.is_empty());
    }

    // (6) Mismatch when the annotation's axis does NOT match the tally's
    // axis is still a FAIL — annotations are scoped strictly per-axis.
    #[test]
    fn tally_fail_when_annotation_axis_mismatches() {
        let mut t = AxisTally::new(Axis::Genomic);
        let case = make_case(
            "in",
            Some(AcceptedDivergence {
                axis: Axis::Normalized,
                policy: Policy::GeneSymbolSelector121,
                note: None,
            }),
            None,
        );
        t.record(&case, "X", Ok("Y".to_string()));
        assert_eq!(t.pass, 0);
        assert!(t.divergence_accepted.is_empty());
        assert_eq!(t.fail.len(), 1);
        assert_eq!(t.fail[0].0, "in");
    }

    // (7a) `Err` result with `accepted_divergence` on the matching axis is
    // STILL FAIL — the annotation accepts a string difference, not a
    // catastrophic failure (panic / parse / normalize Err). A regression
    // that turned a passing-with-divergence case into an Err must not be
    // silenced by the annotation.
    #[test]
    fn tally_err_with_accepted_divergence_still_fails() {
        let mut t = AxisTally::new(Axis::Normalized);
        let case = make_case(
            "in",
            Some(AcceptedDivergence {
                axis: Axis::Normalized,
                policy: Policy::GeneSymbolSelector121,
                note: None,
            }),
            None,
        );
        t.record(&case, "X", Err("normalize: panic boom".to_string()));
        assert_eq!(t.pass, 0);
        assert!(
            t.divergence_accepted.is_empty(),
            "Err must not silence into divergence_accepted bucket"
        );
        assert_eq!(t.fail.len(), 1);
        assert!(t.fail[0].1.contains("err=normalize: panic boom"));
    }

    // (7) Mismatch with no annotation is FAIL.
    #[test]
    fn tally_fail_when_no_annotation() {
        let mut t = AxisTally::new(Axis::Normalized);
        let case = make_case("in", None, None);
        t.record(&case, "X", Ok("Y".to_string()));
        assert_eq!(t.pass, 0);
        assert!(t.divergence_accepted.is_empty());
        assert_eq!(t.fail.len(), 1);
        assert!(t.fail[0].1.contains("expected=\"X\" got=\"Y\""));
    }

    // (8) `spec_citation` routes a string mismatch on the cited axis into
    // the `spec_overridden` bucket — non-failing, separately tracked.
    #[test]
    fn tally_spec_citation_routes_to_spec_overridden() {
        let mut t = AxisTally::new(Axis::Normalized);
        let case = make_case(
            "in",
            None,
            Some(SpecCitation {
                axis: Axis::Normalized,
                section: SpecSection::HgvsPrioritization,
                note: None,
            }),
        );
        t.record(&case, "X", Ok("Y".to_string()));
        assert_eq!(t.pass, 0);
        assert!(t.divergence_accepted.is_empty());
        assert_eq!(t.spec_overridden.len(), 1);
        assert_eq!(t.spec_overridden[0].0, "in");
        assert_eq!(t.spec_overridden[0].1, "HGVS §Prioritization");
        assert!(t.fail.is_empty());
    }

    // (8b) `spec_citation` only catches *string mismatches*, not `Err`
    // results. A panic / parse error on a cited row is still a real bug
    // and must surface as FAIL — same contract as `accepted_divergence`.
    #[test]
    fn tally_err_with_spec_citation_still_fails() {
        let mut t = AxisTally::new(Axis::Normalized);
        let case = make_case(
            "in",
            None,
            Some(SpecCitation {
                axis: Axis::Normalized,
                section: SpecSection::HgvsPrioritization,
                note: None,
            }),
        );
        t.record(&case, "X", Err("normalize: panic boom".to_string()));
        assert_eq!(t.pass, 0);
        assert!(
            t.spec_overridden.is_empty(),
            "Err must not silence into spec_overridden bucket"
        );
        assert_eq!(t.fail.len(), 1);
    }

    // (8c) `spec_citation` on a different axis from the tally's must NOT
    // route to `spec_overridden` — citations are axis-scoped.
    #[test]
    fn tally_spec_citation_other_axis_is_fail() {
        let mut t = AxisTally::new(Axis::Normalized);
        let case = make_case(
            "in",
            None,
            Some(SpecCitation {
                axis: Axis::Errors,
                section: SpecSection::HgvsPrioritization,
                note: None,
            }),
        );
        t.record(&case, "X", Ok("Y".to_string()));
        assert!(t.spec_overridden.is_empty());
        assert_eq!(t.fail.len(), 1);
    }

    // (9) `summary_line` includes the `divergence_accepted` AND
    // `spec_overridden` buckets in the documented stable format.
    #[test]
    fn tally_finish_summary_string_includes_new_bucket() {
        let mut t = AxisTally::new(Axis::Normalized);
        // 1 pass, 1 divergence_accepted, 1 spec_overridden, 1 fail,
        // 0 skipped.
        t.record(&make_case("a", None, None), "X", Ok("X".to_string()));
        t.record(
            &make_case(
                "b",
                Some(AcceptedDivergence {
                    axis: Axis::Normalized,
                    policy: Policy::GeneSymbolSelector121,
                    note: None,
                }),
                None,
            ),
            "X",
            Ok("Y".to_string()),
        );
        t.record(
            &make_case(
                "d",
                None,
                Some(SpecCitation {
                    axis: Axis::Normalized,
                    section: SpecSection::HgvsPrioritization,
                    note: None,
                }),
            ),
            "X",
            Ok("Y".to_string()),
        );
        t.record(&make_case("c", None, None), "X", Ok("Y".to_string()));

        let path = PathBuf::from("/tmp/example.txt");
        let line = t.summary_line(Some(path.as_path()));
        assert_eq!(
            line,
            "normalized: 1 pass / 1 divergence_accepted / 1 spec_overridden / 1 FAIL / 0 skipped \
             (FAIL inputs -> /tmp/example.txt)"
        );
        // `None` omits the path suffix (the gate path).
        let line_no_path = t.summary_line(None);
        assert_eq!(
            line_no_path,
            "normalized: 1 pass / 1 divergence_accepted / 1 spec_overridden / 1 FAIL / 0 skipped"
        );
    }

    // (10) A case may carry BOTH annotations and the bucketing rules
    // remain orthogonal: a mismatch on the matching `accepted_divergence`
    // axis still routes to `divergence_accepted`, regardless of any
    // co-attached `spec_citation`.
    #[test]
    fn tally_coexisting_annotations_route_by_accepted_divergence() {
        let mut t = AxisTally::new(Axis::Normalized);
        let case = make_case(
            "in",
            Some(AcceptedDivergence {
                axis: Axis::Normalized,
                policy: Policy::GeneSymbolSelector121,
                note: None,
            }),
            Some(SpecCitation {
                axis: Axis::Normalized,
                section: SpecSection::HgvsPrioritization,
                note: None,
            }),
        );
        t.record(&case, "X", Ok("Y".to_string()));
        assert_eq!(t.divergence_accepted.len(), 1);
        assert!(t.fail.is_empty());
    }
}
