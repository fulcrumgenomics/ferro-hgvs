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
//! committed `baseline-failures/` snapshot. The
//! `tests/fixtures/mutalyzer-normalize/failure-patterns.md` summary is
//! *generated* from `cases.json` (see `examples/generate_conformance_summary`),
//! never hand-maintained (#509).
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

use ferro_hgvs::conformance::mutalyzer::{
    AcceptedDivergence, Axis, Case, Fixture, Improvement, KnownBug, Policy,
    ReferenceUnavailableReason, SpecCitation, SpecSection,
};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::error_handling::{ErrorOverride, ErrorType};
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::{
    parse_hgvs, FerroError, HgvsVariant, MultiFastaProvider, NormalizeConfig, Normalizer,
    ReferenceProvider, VariantProjector,
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

const FIXTURE_PATH: &str = "tests/fixtures/mutalyzer-normalize/cases.json";
const XFAIL_REPORT_DIR: &str = "/tmp/ferro-xfail";
const FAIL_PRINT_LIMIT: usize = 10;

/// Prefix an axis runner stamps on the `Err` it returns when ferro produced an
/// **empty / degenerate projection** for a row — no protein predicted, or
/// `project_variant_all` yielded zero coding/protein pairs (issue #651).
///
/// This is an *output-quality* signal distinct from bucket membership: the
/// burn-down compares which inputs fail, not the *quality* of each row's `got`,
/// so a populated projection silently degrading to empty (#498's regression
/// shape: a framed `NG_(NM_):c…/NG_(NP_):p…` pair collapsing to `[]`) stays in
/// the same FAIL/annotated bucket and the membership diff reports "0
/// regressions". [`AxisTally::record`] recognizes this prefix and counts it in
/// `empty_got` so a rise is caught even when membership is unchanged. The prefix
/// is stripped from the human-facing diagnostic.
const EMPTY_PROJECTION_SENTINEL: &str = "<empty-projection> ";

/// Directory of committed per-axis empty-projection count baselines (one
/// integer per `<axis>.count` file). A run with a manifest fails if an axis's
/// `empty_got` rises above its baseline; regenerate with
/// `BLESS_EMPTY_PROJECTION=1`. Absent file ⇒ baseline 0.
const EMPTY_PROJECTION_BASELINE_DIR: &str = "tests/fixtures/mutalyzer-normalize/empty-projection";

/// Policy label recorded for `infos`-axis accepted divergences. The upstream
/// mutalyzer code is one ferro intentionally does not model
/// (`error_handling::info_map::NO_FERRO_INFO_EQUIV`): a mutalyzer *internal*
/// diagnostic (e.g. `ISORTEDVARIANTS`, `ICORRECTEDCOORDINATESYSTEM`) that is
/// not part of the HGVS nomenclature spec, so ferro is under no obligation to
/// mirror it. See #326 (spec arbitration).
const INFO_NO_SPEC_EQUIVALENT: &str = "mutalyzer-info-no-spec-equivalent";

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
    let p = Path::new("benchmark-output/manifest.json");
    if p.exists() {
        return Some(p.to_path_buf());
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
    fn genomic_placement(
        &self,
        parent: &ferro_hgvs::hgvs::variant::Accession,
    ) -> Option<ferro_hgvs::reference::GenomicPlacement> {
        // Must delegate, or NG_/LRG_ re-anchoring (#480) silently no-ops behind
        // the wrapper via the trait's `None` default.
        self.0.genomic_placement(parent)
    }
    // The remaining overridable methods must be forwarded too — the same
    // wrapper-drops-overrides bug #726 fixed for the blanket `Arc<T>`/`Box<T>`
    // impls also affected this hand-written newtype, which #726's audit missed.
    // Without these, the conformance harness silently used the trait defaults
    // (`None` / the `genomic_context`-blind `get_transcript` fallback), so rows
    // that ferro resolves in production (legacy gene-model selectors per
    // #500/#709; intronic normalization via the build-aware transcript probe)
    // surfaced here as spurious `Err` divergences rather than real ferro
    // behavior. Mirror the blanket impl exactly so this wrapper measures ferro,
    // not the wrapper.
    fn get_transcript_for_variant(
        &self,
        variant: &ferro_hgvs::hgvs::variant::HgvsVariant,
    ) -> Result<Arc<Transcript>, FerroError> {
        self.0.get_transcript_for_variant(variant)
    }
    fn get_transcript_for_accession(
        &self,
        accession: &ferro_hgvs::hgvs::variant::Accession,
    ) -> Result<Arc<Transcript>, FerroError> {
        self.0.get_transcript_for_accession(accession)
    }
    fn has_transcript_version_exact(&self, id: &str) -> bool {
        self.0.has_transcript_version_exact(id)
    }
    fn genomic_placement_on_build(
        &self,
        parent: &ferro_hgvs::hgvs::variant::Accession,
        build: Option<&str>,
    ) -> Option<ferro_hgvs::reference::GenomicPlacement> {
        self.0.genomic_placement_on_build(parent, build)
    }
    fn resolve_legacy_gene_selector(&self, selector: &str) -> Option<String> {
        self.0.resolve_legacy_gene_selector(selector)
    }
    fn get_protein_length(&self, accession: &str) -> Result<u64, FerroError> {
        self.0.get_protein_length(accession)
    }
    fn get_sequence_length(&self, id: &str) -> Result<u64, FerroError> {
        self.0.get_sequence_length(id)
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
    /// Mismatches gated by a `known_bug` annotation whose `axis` matches
    /// this tally's axis. Tracked-but-non-failing (xfail); reported in the
    /// summary line as `known_bug` and NOT written to
    /// `baseline-failures/<axis>.txt`. Stores `(input, "#<issue>")`.
    known_bug: Vec<(String, String)>,
    /// Mismatches gated by an `improvement` annotation whose `axis` matches
    /// this tally's axis. Tracked-but-non-failing (xfail); reported in the
    /// summary line as `improvement` and NOT written to
    /// `baseline-failures/<axis>.txt`. Stores `(input, "#<issue>")`.
    improvement: Vec<(String, String)>,
    /// Mismatches gated by a `reference_unavailable` annotation whose `axis`
    /// matches this tally's axis: the prepared manifest lacks the reference, so
    /// ferro no-ops — a reference-availability artifact, NOT a ferro defect.
    /// Tracked-but-non-failing; reported separately from `known_bug`. Stores
    /// `(input, "<reason> #<issue>")`.
    reference_unavailable: Vec<(String, String)>,
    fail: Vec<(String, String)>, // (input, diagnostic)
    skipped: usize,
    /// Count of rows for which ferro produced an **empty / degenerate
    /// projection** (tagged by [`EMPTY_PROJECTION_SENTINEL`]) — an
    /// output-quality signal independent of which bucket the row lands in
    /// (#651). Gated against a committed baseline in [`finish`](Self::finish).
    empty_got: usize,
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
            reference_unavailable: Vec::new(),
            fail: Vec::new(),
            skipped: 0,
            empty_got: 0,
        }
    }

    /// Bucket a single case against `expected`/`actual`:
    /// 1. On a match (`actual == expected`):
    ///    - XPASS FAIL when the case carries an `accepted_divergence`,
    ///      `known_bug`, or `improvement` annotation for this tally's axis —
    ///      the annotation is now stale (the divergence/bug/non-preferred form
    ///      is gone) and must be removed.
    ///    - PASS otherwise.
    /// 2. `divergence_accepted` when the case carries
    ///    `accepted_divergence` matching this tally's axis.
    /// 3. `known_bug` (xfail) when the case carries a `known_bug`
    ///    matching this tally's axis.
    /// 4. `improvement` (xfail) when the case carries an `improvement`
    ///    matching this tally's axis.
    /// 5. `reference_unavailable` when the case carries a `reference_unavailable`
    ///    annotation matching this tally's axis (the prepared manifest lacks the
    ///    reference, so ferro no-ops — not a ferro defect).
    /// 6. `spec_overridden` when the case carries a `spec_citation`
    ///    matching this tally's axis.
    /// 7. FAIL otherwise.
    ///
    /// `accepted_divergence`, `known_bug`, `improvement`, `reference_unavailable`,
    /// and `spec_citation` catch a *string mismatch* from a successful run on
    /// the projection axes. On those axes an `Err` means ferro failed to produce
    /// any result (panic, parse error, normalize failure) — a real bug that must
    /// surface as FAIL even when an annotation is present.
    ///
    /// The **errors axis is the exception**: there a non-matching `Err` is the
    /// *expected* "ferro diverges" signal — ferro emitted no error, or a
    /// different error class than the mutalyzer code the row expects — not a
    /// crash. So on the errors axis the annotation buckets also apply to an
    /// `Err` actual, mirroring the OK-path on every other axis. A genuine panic
    /// (`"panic: …"`) is still a hard FAIL there, so real crashes are never
    /// masked by an annotation.
    fn record(&mut self, case: &Case, expected: &str, actual: Result<String, String>) {
        // Output-quality signal (#651): the axis runner stamps
        // `EMPTY_PROJECTION_SENTINEL` on its `Err` when ferro produced nothing
        // (no protein / zero pairs). Count it independently of the bucket the
        // row falls into below, so a populated→empty degradation is visible
        // even when bucket membership is unchanged. (An empty projection is
        // always an `Err`, so it still falls through to `fail` as before.)
        if matches!(&actual, Err(e) if e.starts_with(EMPTY_PROJECTION_SENTINEL)) {
            self.empty_got += 1;
        }

        let matches = matches!(&actual, Ok(s) if s == expected);
        if matches {
            // XPASS detection: a match on an axis that carries an
            // `accepted_divergence` or `known_bug` annotation means the
            // annotation is stale — the divergence/bug the row documents is
            // gone. Fail loudly so the annotation (and, for known_bug, the
            // fixed row) is cleaned up rather than silently rotting.
            if let Some(ad) = &case.accepted_divergence {
                if ad.axis == self.axis {
                    self.fail.push((
                        case.input.clone(),
                        format!(
                            "XPASS: accepted_divergence (policy {}) now matches mutalyzer; the divergence is gone — remove the annotation",
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
                            "XPASS: known_bug #{} now matches mutalyzer; the fix appears to have landed — remove the annotation and demote the row",
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
                            "XPASS: improvement #{} now matches mutalyzer; ferro reached the spec-preferred form — remove the annotation and demote the row",
                            imp.tracking_issue
                        ),
                    ));
                    return;
                }
            }
            if let Some(ru) = &case.reference_unavailable {
                if ru.axis == self.axis {
                    self.fail.push((
                        case.input.clone(),
                        format!(
                            "XPASS: reference_unavailable ({}) #{} now matches mutalyzer; the reference became available — remove the annotation and demote the row",
                            ru.reason, ru.tracking_issue
                        ),
                    ));
                    return;
                }
            }
            self.pass += 1;
            return;
        }
        // On the errors axis a non-panic `Err` is the expected "ferro diverges"
        // outcome, so annotations bucket it the same way an `Ok` mismatch is
        // bucketed elsewhere; a genuine panic stays a hard FAIL. See the
        // `record` doc comment.
        let errors_axis_divergence =
            self.axis == Axis::Errors && matches!(&actual, Err(e) if !e.starts_with("panic:"));
        if actual.is_ok() || errors_axis_divergence {
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
            if let Some(ru) = &case.reference_unavailable {
                if ru.axis == self.axis {
                    self.reference_unavailable.push((
                        case.input.clone(),
                        format!("{} #{}", ru.reason, ru.tracking_issue),
                    ));
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
            Err(e) => {
                // Strip the internal empty-projection sentinel from the
                // human-facing diagnostic (it's already accounted in
                // `empty_got`); tag it so the FAIL line still reads clearly.
                match e.strip_prefix(EMPTY_PROJECTION_SENTINEL) {
                    Some(rest) => format!("expected={expected:?} err=[empty projection] {rest}"),
                    None => format!("expected={expected:?} err={e}"),
                }
            }
        };
        self.fail.push((case.input.clone(), diag));
    }

    /// Single-line, grep-friendly summary of the tally state. Stable
    /// format used by both `finish()` and unit tests.
    fn summary_line(&self, report_path: &Path) -> String {
        format!(
            "{}: {} pass / {} divergence_accepted / {} known_bug / {} improvement / {} reference_unavailable / {} spec_overridden / {} FAIL / {} empty / {} skipped (FAIL inputs -> {})",
            self.axis,
            self.pass,
            self.divergence_accepted.len(),
            self.known_bug.len(),
            self.improvement.len(),
            self.reference_unavailable.len(),
            self.spec_overridden.len(),
            self.fail.len(),
            self.empty_got,
            self.skipped,
            report_path.display()
        )
    }

    fn finish(self) {
        // Write FAIL inputs (one per line) to /tmp/ferro-xfail/<axis>.txt and
        // (input \t diagnostic) pairs to /tmp/ferro-xfail/<axis>.tsv so the
        // committed baseline-failures/<axis>.txt can be regenerated from a
        // single run. (The failure-patterns.md summary is generated from
        // cases.json instead — see examples/generate_conformance_summary.)
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

        for (input, reason) in self.reference_unavailable.iter().take(FAIL_PRINT_LIMIT) {
            eprintln!(
                "  REFERENCE_UNAVAILABLE  [{}] {} | {}",
                self.axis, input, reason
            );
        }
        if self.reference_unavailable.len() > FAIL_PRINT_LIMIT {
            eprintln!(
                "  ... {} more reference_unavailable",
                self.reference_unavailable.len() - FAIL_PRINT_LIMIT,
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

        // Output-quality gate (#651): an empty/degenerate-projection count that
        // rose above the committed baseline is a regression even if the FAIL
        // *set* is unchanged. `finish` only runs with a manifest (the axis tests
        // early-return otherwise), so this never affects manifest-less CI.
        let baseline_path = empty_projection_baseline_path(self.axis);
        if std::env::var("BLESS_EMPTY_PROJECTION").is_ok() {
            if self.empty_got > 0 {
                let _ = fs::create_dir_all(Path::new(EMPTY_PROJECTION_BASELINE_DIR));
                fs::write(&baseline_path, format!("{}\n", self.empty_got))
                    .unwrap_or_else(|e| panic!("write {}: {e}", baseline_path.display()));
            } else {
                // An axis with no empty projections needs no baseline file
                // (load treats a missing file as 0); drop any stale one.
                let _ = fs::remove_file(&baseline_path);
            }
            eprintln!(
                "blessed empty-projection baseline {} = {}",
                self.axis, self.empty_got
            );
            return;
        }
        let baseline = load_empty_projection_baseline(self.axis);
        if let Err(msg) = empty_budget_verdict(self.empty_got, baseline) {
            panic!(
                "{}: {msg} — a populated projection degraded to empty (output quality regressed \
                 even though the FAIL set may be unchanged). Investigate, or re-bless with \
                 BLESS_EMPTY_PROJECTION=1 if the increase is intentional.",
                self.axis
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

/// Path of the committed empty-projection count baseline for `axis`.
fn empty_projection_baseline_path(axis: Axis) -> PathBuf {
    Path::new(EMPTY_PROJECTION_BASELINE_DIR).join(format!("{axis}.count"))
}

/// Load the committed empty-projection baseline count for `axis`. A missing file
/// is treated as baseline `0` (no empty projections expected); an existing but
/// unparseable file panics so corruption surfaces immediately rather than being
/// silently read as `0` (which would lose the true baseline and misreport the
/// regression as `0 -> N`).
fn load_empty_projection_baseline(axis: Axis) -> usize {
    let path = empty_projection_baseline_path(axis);
    match fs::read_to_string(&path) {
        Ok(s) => s.trim().parse::<usize>().unwrap_or_else(|e| {
            panic!("baseline file {} is not a valid usize: {e}", path.display())
        }),
        Err(_) => 0, // missing file -> baseline 0
    }
}

/// Verdict for the #651 empty-projection budget: `Ok` while the observed count
/// is at or below the committed baseline (equal = unchanged, below =
/// improvement), `Err` with a message when it rises (a regression).
fn empty_budget_verdict(observed: usize, baseline: usize) -> Result<(), String> {
    if observed > baseline {
        Err(format!(
            "empty-projection count rose {baseline} -> {observed}"
        ))
    } else {
        Ok(())
    }
}

// ----------------------------------------------------------------------------
// Helpers
// ----------------------------------------------------------------------------

fn transcript_of(v: &HgvsVariant) -> Option<String> {
    // Return the *bare* transcript accession (e.g. `NM_000077.4`), not the
    // compound `genomic_context(transcript)` form `full()` produces (e.g.
    // `NG_007485.1(NM_000077.4)`). `VariantProjector::project_variant`
    // compares the requested transcript id against `transcript_accession()`,
    // so the compound form spuriously trips its transcript_id-mismatch guard
    // for every `NG_(NM_)` input (#326).
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
            "mutalyzer-normalize: skipping — no manifest at FERRO_MANIFEST \
             or benchmark-output/manifest.json"
        ),
    }
}

/// Issue #506: the direct c.→p. path must extend to bare `n.` (Tx) inputs.
///
/// `NM_003002.4:n.206_210del` carries no `genomic_context` parent, so the
/// genome-pivot path cannot run; before #506 it hard-errored with "input
/// variant has no parent reference". The direct path converts the 1-based
/// transcript position to a CDS position via the transcript's `cds_start`,
/// then predicts protein identically to the coding path. Pinned against a
/// live manifest so the worked example resolves end-to-end.
///
/// The protein *consequence* matches mutalyzer's `p.(Gly58GlnfsTer9)`. The
/// protein *accession* renders in ferro's spec-correct bare-`NP_` form
/// (`NP_002993.1:p.…`) rather than the `NM_(NP_)` wrapper the upstream
/// `protein_description` axis carries — that wrapper-vs-bare-NP difference is
/// the pre-existing, axis-wide format-only divergence tracked under #498
/// (#483/#495), out of scope here. This test pins #506's guarantee: a bare
/// `n.` input now *produces* the correct prediction instead of erroring.
#[test]
fn issue_506_bare_n_transcript_predicts_protein() {
    let Some(vp) = variant_projector() else {
        eprintln!("issue_506_bare_n_transcript_predicts_protein: skipping — no manifest");
        return;
    };
    let v = ferro_hgvs::parse_hgvs("NM_003002.4:n.206_210del").expect("parse n. input");
    let proj = vp
        .project_variant(&v, "NM_003002.4")
        .expect("bare n. input should project via the direct n.→CDS→p. path");
    let protein = proj
        .protein
        .as_ref()
        .expect("protein consequence should be predicted for a bare n. input");
    assert_eq!(
        protein.to_string(),
        "NP_002993.1:p.(Gly58GlnfsTer9)",
        "worked example n.206_210del must predict the spec frameshift consequence"
    );
    // The coding axis must reframe n.206_210del to the spec-normalized c. form.
    let coding = proj
        .coding
        .as_ref()
        .expect("coding (c.) axis should be derived");
    assert_eq!(
        coding.to_string(),
        "NM_003002.4(SDHD):c.171_175del",
        "n.206_210del should reframe to c.171_175del"
    );
    // A bare-NM_ n. input has no genome alignment, so no genomic form.
    assert!(
        proj.genomic.is_none(),
        "bare n. input must report genomic = None, got {:?}",
        proj.genomic
    );
}

/// Issue #506 regression: a 3'-downstream (`*N`) position on a bare `n.` / `r.`
/// input must NOT be mis-projected through the direct n./r.→CDS→p. path.
///
/// The parser stores `n.*5` as `base = 5, downstream = true`, but the exon-aware
/// `tx_to_cds` mapper the direct path uses classifies purely by `base` vs the CDS
/// bounds and never reads the `downstream` flag. Without an explicit guard,
/// `n.*5` (3'UTR / transcript-tail) is silently read as in-transcript position 5,
/// which for a CDS-at-1 transcript lands inside the CDS — producing a bogus `c.`
/// form and running protein prediction on it. The direct path has no exon-aware
/// 3'UTR translation (that lives on the genome-pivot path), so it must decline
/// with `UnsupportedProjection` rather than emit a wrong answer.
#[test]
fn issue_506_bare_downstream_position_is_rejected() {
    let Some(vp) = variant_projector() else {
        eprintln!("issue_506_bare_downstream_position_is_rejected: skipping — no manifest");
        return;
    };

    // Both a bare `n.*N` and a bare `r.*N` input must be declined, not
    // mis-projected into a spurious `c.` form / protein consequence.
    for input in ["NM_003002.4:n.*5del", "NM_003002.4:r.*5del"] {
        let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
        match vp.project_variant(&v, "NM_003002.4") {
            Err(FerroError::UnsupportedProjection { .. }) => {} // expected
            Err(other) => panic!(
                "{input}: expected UnsupportedProjection for a 3'-downstream bare position, \
                 got a different error: {other}"
            ),
            Ok(proj) => panic!(
                "{input}: a 3'-downstream bare position must be rejected, but it projected to \
                 coding={:?} protein={:?}",
                proj.coding, proj.protein
            ),
        }
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
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: genomic
// ----------------------------------------------------------------------------

#[test]
fn axis_genomic() {
    let (Some(projector), Some(normalizer)) = (variant_projector(), normalizer()) else {
        eprintln!("axis_genomic: skipping — no manifest");
        return;
    };

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
            // Project transcript-coordinate (c./n./r.) inputs onto their genomic
            // reference frame; pass g./allele/etc. inputs through unchanged
            // (project_to_genomic only handles single c./n./r. variants — #328).
            // Then normalize the genomic form so the output matches mutalyzer's
            // normalized g. (3'-shift, ins→dup, …).
            let g = if matches!(
                v,
                HgvsVariant::Cds(_) | HgvsVariant::Tx(_) | HgvsVariant::Rna(_)
            ) {
                projector
                    .project_to_genomic(&v)
                    .map_err(|e| format!("project: {e}"))?
            } else {
                v
            };
            let n = normalizer
                .normalize(&g)
                .map_err(|e| format!("normalize: {e}"))?;
            Ok(format!("{n}"))
        });

        t.record(case, expected, actual);
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
                // Empty/degenerate projection (#651): expected a protein but
                // ferro predicted none. Tag it so the harness counts the
                // output-quality degradation, not just the bucket membership.
                .ok_or_else(|| format!("{EMPTY_PROJECTION_SENTINEL}no protein predicted"))
        });

        t.record(case, expected, actual);
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

            // Empty/degenerate projection (#651): mutalyzer expects at least
            // one coding/protein pair but `project_variant_all` yielded none.
            // Tag it as an empty projection (distinct from "produced some pairs
            // but missing a specific one") so a populated→empty degradation is
            // counted as an output-quality regression, not just a set member.
            if actual_pairs.is_empty() && pairs.iter().any(|p| p.len() == 2) {
                return Err(format!(
                    "{EMPTY_PROJECTION_SENTINEL}project_variant_all produced no coding/protein pairs; expected {expected_repr}"
                ));
            }

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
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: errors
// ----------------------------------------------------------------------------

#[test]
fn axis_errors() {
    let Some(p) = provider() else {
        eprintln!("axis_errors: skipping — no manifest");
        return;
    };
    // The errors axis asserts STRICT-mode behavior for reference mismatches:
    // ferro's library default intentionally warn-corrects them (an accepted
    // divergence from mutalyzer), so the default lenient `normalizer()` used by
    // the other axes can never reproduce mutalyzer's ESEQUENCEMISMATCH /
    // EREPEATMISMATCH rejections. Strict mode does. See #486 and the design doc.
    let cfg = NormalizeConfig::new()
        .with_error_override(ErrorType::RefSeqMismatch, ErrorOverride::Reject)
        // Position-past-end (#486 EOUTOFBOUNDARY) — same strict-mode-on-the-
        // errors-axis rationale as RefSeqMismatch above; ferro's default is
        // lenient warn-only. The bounds check itself is pre-existing.
        .with_error_override(ErrorType::PositionPastEnd, ErrorOverride::Reject)
        // Intronic-on-bare-transcript (#486 EINTRONIC) — same strict-mode-on-
        // the-errors-axis rationale: ferro's default is a lenient warn; the
        // spec-compliant reject is asserted here. See the #486 design doc.
        .with_error_override(ErrorType::IntronicOnBareTranscript, ErrorOverride::Reject)
        // Overlapping cis-allele edits (#486 EOVERLAP) — ferro's default is a
        // lenient W5002 warn; the spec-compliant reject is asserted on the
        // errors axis, same rationale as the categories above. Covers both
        // coincident-bounds and insertion-junction overlaps.
        .with_error_override(ErrorType::OverlapConflictingEdits, ErrorOverride::Reject);
    let normalizer = Normalizer::with_config(ArcProvider(p), cfg);

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
    t.finish();
}

// ----------------------------------------------------------------------------
// Axis: infos
// ----------------------------------------------------------------------------

#[test]
fn axis_infos() {
    let Some(normalizer) = normalizer() else {
        eprintln!("axis_infos: skipping — no manifest");
        return;
    };

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
        // separately so it stays visible.
        //
        // An over-emission `Err` where the case carries an `accepted_divergence`
        // on `Axis::Infos` is also a tracked divergence — ferro emits a
        // spec-valid info code that mutalyzer does not (e.g. per-member
        // `SHUFFLE_APPLIED` for a genuine 3' shift in a compound allele,
        // where mutalyzer only reports the allele-level `ISORTEDVARIANTS`).
        // See #499 (precedent #481). All other errors (under-emission, parse,
        // normalize, panic) are real failures routed through `record`.
        if actual.is_ok() && has_no_equiv {
            t.divergence_accepted
                .push((case.input.clone(), INFO_NO_SPEC_EQUIVALENT.to_string()));
        } else if let Err(ref e) = actual {
            if e.starts_with("ferro emitted extra info") {
                if let Some(ad) = &case.accepted_divergence {
                    if ad.axis == Axis::Infos {
                        t.divergence_accepted
                            .push((case.input.clone(), ad.policy.to_string()));
                        continue;
                    }
                }
            }
            t.record(case, &expected_repr, actual);
        } else {
            t.record(case, &expected_repr, actual);
        }
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

// ----------------------------------------------------------------------------
// Hermetic gate: snapshot-backed `normalized` axis (CI-always)
// ----------------------------------------------------------------------------
//
// Unlike `regression_under_mock_normalized` (a `MockProvider` with NO bases, so
// it pins ferro's no-op behavior), this gate serves the committed transcript
// snapshot's REAL version-exact bases + CDS (issue #719, increment I4). It runs
// every snapshot-covered `normalized` case against `cases.json` truth with NO
// manifest, so the normalized axis becomes a real per-PR correctness gate, not
// just a refactor pin.
//
// Coverage is the snapshot's transcript set. A row whose bare transcript the
// snapshot does not carry — a genomic-anchored parent needing NG bases, or a
// transcript absent from the local reference — is skipped here; it stays covered
// by the manifest-backed `axis_normalized` (and, for genomic, the deferred
// genomic snapshot). A skip is not a failure.
//
// Annotation-aware via the shared `AxisTally`: a snapshot-covered row carrying a
// `known_bug`/`accepted_divergence`/`improvement`/`spec_citation`/
// `reference_unavailable` annotation for the normalized axis is bucketed, not
// failed — so this gate asserts only genuine, unannotated divergences on real
// bases.

const SNAPSHOT_DIR: &str = "tests/fixtures/mutalyzer-normalize/reference-snapshot";

/// Whether a row is soundly coverable by the **transcript** snapshot — used by
/// both the normalized and protein hermetic gates. The snapshot models each
/// transcript as a single-exon, transcript-space sequence + CDS with no cdot or
/// genomic context, so snapshot-derived `c.` normalization / `c.→p.` projection
/// equal the manifest-backed axes for bare `c.` exonic rows but NOT for:
/// - a genomic / compound context `(...)` — needs the NG bases (deferred genomic
///   snapshot);
/// - intronic offset positions (`c.52+1`, `c.378-17`, `c.*824+10`) — need the
///   multi-exon structure;
/// - `n.` input on a coding transcript — mutalyzer's `n.→c.` conversion needs
///   cdot; ferro keeps `n.` without it (so restrict to `c.` / [`HgvsVariant::Cds`]);
/// - a `delins` against a *coordinate range* (`delinsN_M…`) — rendering the
///   referenced/inverted segment needs reference context the snapshot lacks.
///
/// Excluded rows are not failures here; they stay covered by the manifest-backed
/// axes (and, for genomic, the deferred genomic snapshot).
fn snapshot_coverable(input: &str, v: &HgvsVariant) -> bool {
    if !matches!(v, HgvsVariant::Cds(_)) {
        return false;
    }
    if input.contains('(') {
        return false;
    }
    let coord = match input.split_once(":c.") {
        Some((_, c)) => c,
        None => return false,
    };
    !has_intron_offset(coord) && !delins_coordinate_range(coord)
}

/// True if `coord` has an intron offset: a digit immediately followed by `+`/`-`
/// then a digit (e.g. `52+1`, `378-17`, `*824+10`). The `+`/`-` in a UTR marker
/// (`-1`, `*1`) is not between two digits, so those exonic positions are kept.
fn has_intron_offset(coord: &str) -> bool {
    let b = coord.as_bytes();
    (1..b.len().saturating_sub(1)).any(|i| {
        (b[i] == b'+' || b[i] == b'-') && b[i - 1].is_ascii_digit() && b[i + 1].is_ascii_digit()
    })
}

/// True if `coord` is a `delins` against a coordinate range — `delins` directly
/// followed by a digit (e.g. `delins190_220inv`), as opposed to literal bases
/// (`delinsATC`).
fn delins_coordinate_range(coord: &str) -> bool {
    coord
        .find("delins")
        .and_then(|idx| coord.as_bytes().get(idx + "delins".len()))
        .is_some_and(u8::is_ascii_digit)
}

#[test]
fn gate_normalized_snapshot() {
    use ferro_hgvs::conformance::reference_snapshot::{load_provider, TranscriptSnapshot};
    use std::collections::HashSet;

    let provider = load_provider(SNAPSHOT_DIR)
        .unwrap_or_else(|e| panic!("load snapshot provider from {SNAPSHOT_DIR}: {e}"));
    // Derive coverage from transcripts the provider actually serves, not raw
    // metadata keys: `to_provider` skips any metadata entry whose bases are
    // absent from the FASTA, so a metadata-only entry would otherwise count as
    // covered yet no-op (a divergent FAIL) instead of being skipped. Filtering
    // through `has_transcript` keeps `covered` and the served set in lockstep.
    let covered: HashSet<String> = TranscriptSnapshot::from_json_path(
        Path::new(SNAPSHOT_DIR).join("transcripts.metadata.json"),
    )
    .expect("snapshot metadata loads")
    .transcripts
    .keys()
    .filter(|accession| provider.has_transcript(accession))
    .cloned()
    .collect();
    let normalizer = Normalizer::new(provider);

    let mut t = AxisTally::new(Axis::Normalized);
    let mut covered_count = 0usize;
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected) = case.normalized.as_deref() else {
            continue;
        };
        // A parse error or non-transcript input is out of this gate's scope
        // (still covered by the manifest-backed axis), not a failure here.
        let Ok(v) = parse_hgvs(&case.input) else {
            continue;
        };
        let Some(tx) = transcript_of(&v) else {
            continue;
        };
        // Coverage: the snapshot must carry the transcript, and the row must be
        // soundly coverable by a transcript-only (single-exon, no-cdot) model.
        if !covered.contains(&tx) || !snapshot_coverable(&case.input, &v) {
            continue;
        }
        covered_count += 1;

        let actual = catch_panics(|| -> Result<String, String> {
            let n = normalizer
                .normalize(&v)
                .map_err(|e| format!("normalize: {e}"))?;
            Ok(format!("{n}"))
        });
        t.record(case, expected, actual);
    }

    assert!(
        covered_count > 0,
        "gate_normalized_snapshot: snapshot covered zero normalized rows — the \
         coverage filter or the committed snapshot is broken"
    );
    assert!(
        t.fail.is_empty(),
        "gate_normalized_snapshot: {} hermetic divergence(s) from cases.json on \
         {covered_count} snapshot-covered normalized row(s) (real version-exact \
         bases, no manifest):\n{}",
        t.fail.len(),
        t.fail
            .iter()
            .take(20)
            .map(|(i, d)| format!("  {i} | {d}"))
            .collect::<Vec<_>>()
            .join("\n"),
    );
    println!(
        "gate_normalized_snapshot: {covered_count} snapshot-covered normalized row(s) \
         pass hermetically (no manifest); {} known_bug, {} improvement, {} reference_unavailable",
        t.known_bug.len(),
        t.improvement.len(),
        t.reference_unavailable.len(),
    );
}

/// Hermetic gate for the `protein_description` axis (CI-always), the analog of
/// `gate_normalized_snapshot`. The snapshot carries each transcript's CDS + bases,
/// so a `CdotMapper` built from it (`TranscriptSnapshot::to_cdot`) drives
/// `c.→p.` projection with no manifest or network — the protein consequence is
/// transcript-frame translation. Coverage is the same transcript-only subset
/// (`snapshot_coverable`); genomic-anchored / intronic rows stay on the nightly
/// manifest axis.
#[test]
fn gate_protein_snapshot() {
    use ferro_hgvs::conformance::reference_snapshot::{load_sequences, TranscriptSnapshot};
    use ferro_hgvs::data::projection::Projector;
    use std::collections::HashSet;

    let dir = Path::new(SNAPSHOT_DIR);
    let snapshot = TranscriptSnapshot::from_json_path(dir.join("transcripts.metadata.json"))
        .expect("snapshot metadata loads");
    let sequences = load_sequences(dir).expect("snapshot FASTA loads");
    let covered: HashSet<String> = snapshot.transcripts.keys().cloned().collect();
    let projector = VariantProjector::new(
        Projector::new(snapshot.to_cdot(&sequences)),
        snapshot.to_provider(&sequences),
    );

    let mut t = AxisTally::new(Axis::ProteinDescription);
    let mut covered_count = 0usize;
    let mut empty_skipped = 0usize;
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected) = case.protein_description.as_deref() else {
            continue;
        };
        let Ok(v) = parse_hgvs(&case.input) else {
            continue;
        };
        let Some(tx) = transcript_of(&v) else {
            continue;
        };
        if !covered.contains(&tx) || !snapshot_coverable(&case.input, &v) {
            continue;
        }

        let actual = catch_panics(|| -> Result<String, String> {
            let result = projector
                .project_variant(&v, &tx)
                .map_err(|e| format!("project: {e}"))?;
            result
                .protein
                .as_ref()
                .map(|p| format!("{p}"))
                .ok_or_else(|| "<empty>".to_string())
        });
        // Empty protein predictions (start-loss `p.?`, synonymous `p.(=)`) are a
        // separate output-quality concern (#651), and occur on transcripts the
        // manifest axis treats as reference_unavailable (so it never validated
        // them). Report, don't gate.
        if matches!(&actual, Err(e) if e == "<empty>") {
            empty_skipped += 1;
            continue;
        }
        covered_count += 1;
        // Compare the spec-correct **bare** `NP_:p.` consequence: mutalyzer wraps
        // it as `NM_x(NP_y):p.…`, which ferro intentionally does not emit (#498,
        // the corpus's `bare-np-protein` spec_citation cluster).
        t.record(case, &to_bare_protein(expected), actual);
    }

    assert!(
        covered_count > 0,
        "gate_protein_snapshot: snapshot covered zero protein rows — the coverage \
         filter or the committed snapshot is broken"
    );
    assert!(
        t.fail.is_empty(),
        "gate_protein_snapshot: {} hermetic divergence(s) from cases.json on \
         {covered_count} snapshot-covered protein row(s) (real bases + CDS, no manifest):\n{}",
        t.fail.len(),
        t.fail
            .iter()
            .take(20)
            .map(|(i, d)| format!("  {i} | {d}"))
            .collect::<Vec<_>>()
            .join("\n"),
    );
    println!(
        "gate_protein_snapshot: {covered_count} snapshot-covered protein row(s) pass \
         hermetically (no manifest); {empty_skipped} empty-projection skipped (#651); \
         {} spec_overridden, {} known_bug",
        t.spec_overridden.len(),
        t.known_bug.len(),
    );
}

/// Committed hermetic genomic-window fixture for the genomic-axis gate (#737):
/// real transcripts (multi-exon), padded chromosome windows, and the
/// version-independent `NG_`/`LRG_` placements the projection re-anchors through.
/// Regenerated by `examples/extract_mutalyzer_genomic_windows.rs`.
const GENOMIC_WINDOWS_PATH: &str = "tests/fixtures/mutalyzer-normalize/genomic-windows.json";

/// Whether a row is soundly coverable by the **genomic** window fixture. The
/// genomic gate projects `c./n./r.` → `g.` (re-anchored into the input's
/// `NG_`/`LRG_`/`NC_` parent frame) and normalizes the result, so coverage
/// requires the reference data the fixture commits — a real transcript the cdot
/// can map and (for an `NG_`/`LRG_` parent) its captured placement. It covers a
/// single exonic `c.` variant whose parent carries an **`NM_`/`NR_` transcript
/// selector**, excluding the classes the fixture/gate deliberately do not gate:
/// - a gene-symbol / legacy locus selector (`(TIMM8B_v001)`, `(BRCA2)`) — the
///   #500-class selector-resolution divergence (ferro keeps the input selector);
/// - intronic offset positions — need flanking the gate does cover, but the
///   exon-junction coordinate convention is the deferred #742 fix;
/// - compound (`[...;...]`) and coordinate-range / repeat / special (`pter`)
///   forms — rendering needs context out of this gate's scope.
///
/// Excluded rows stay on the manifest-backed `axis_genomic`.
fn genomic_coverable(input: &str) -> bool {
    // Must carry an explicit NM_/NR_ transcript selector on a genomic parent.
    let Some((parent, rest)) = input.split_once('(') else {
        return false;
    };
    if !(parent.starts_with("NG_") || parent.starts_with("LRG_") || parent.starts_with("NC_")) {
        return false;
    }
    // NG_008835.1 is a large (>400 KB) RefSeqGene with homopolymer runs where the
    // genomic 3' shift / repeat detection needs more flanking sequence than the
    // committed padded windows carry; serving it via windows truncates the shift
    // (and a deletion that the full sequence renders as the spec-valid repeat
    // contraction — an accepted divergence from mutalyzer's del, #745 — collapses
    // to del). Those rows stay on the manifest-backed axis_genomic rather than
    // being faithfully gateable here.
    if parent.starts_with("NG_008835.") {
        return false;
    }
    let Some((selector, _)) = rest.split_once(')') else {
        return false;
    };
    if !(selector.starts_with("NM_") || selector.starts_with("NR_")) {
        return false;
    }
    // Single exonic c. variant only (no compound, intronic, range-delins, repeat,
    // or special position).
    let Some((_, coord)) = input.split_once(":c.") else {
        return false;
    };
    if input.contains('[')
        || coord.contains("pter")
        || coord.contains("qter")
        || coord.contains('[')
    {
        return false;
    }
    !has_intron_offset(coord) && !delins_coordinate_range(coord)
}

/// True if the *expected* genomic output renders the edit against a coordinate
/// range / cross-reference (`ins[…]`, `ins<digit>…`) rather than literal bases.
/// Reproducing that rendering needs cross-reference context the hermetic windows
/// do not commit (ferro emits the equivalent literal `ins`), so such rows stay
/// on the manifest-backed `axis_genomic`.
fn expected_uses_range_reference(expected: &str) -> bool {
    expected.contains("ins[")
        || expected
            .find("ins")
            .and_then(|i| expected.as_bytes().get(i + 3))
            .is_some_and(u8::is_ascii_digit)
}

/// Hermetic gate for the `genomic` axis (#737), the analog of
/// `gate_normalized_snapshot`. Builds a [`WindowProvider`] from the committed
/// `genomic-windows.json` (real transcripts + chromosome windows + `NG_`/`LRG_`
/// placements) and a [`CdotMapper`] from the same transcripts, then projects
/// each covered `c.` row to `g.` and normalizes it — exercising the #737
/// genomic-frame re-normalization with no manifest or network. Divergences that
/// remain a tracked bug (`known_bug`) or spec-non-preferred (`improvement`) are
/// bucketed by `AxisTally` from the cases.json annotations, not failed.
#[test]
fn gate_genomic_snapshot() {
    use ferro_hgvs::conformance::reference_window::WindowFixture;
    use ferro_hgvs::data::cdot::CdotMapper;
    use std::collections::HashSet;

    let windows = WindowFixture::from_json_path(Path::new(GENOMIC_WINDOWS_PATH))
        .unwrap_or_else(|e| panic!("load genomic windows from {GENOMIC_WINDOWS_PATH}: {e}"));
    let cdot = CdotMapper::from_transcripts(windows.transcripts.iter());
    let provider = windows.to_provider();
    // Coverage is the set of transcripts the rebuilt cdot can actually map.
    let covered: HashSet<String> = windows
        .transcripts
        .iter()
        .filter(|tx| cdot.get_transcript(&tx.id).is_some())
        .map(|tx| tx.id.clone())
        .collect();
    let projector = VariantProjector::new(Projector::new(cdot), provider.clone());
    let normalizer = Normalizer::new(provider);

    let mut t = AxisTally::new(Axis::Genomic);
    let mut covered_count = 0usize;
    let mut declined = 0usize;
    for case in &fixture().cases {
        if !case.to_test {
            continue;
        }
        let Some(expected) = case.genomic.as_deref() else {
            continue;
        };
        let Ok(v) = parse_hgvs(&case.input) else {
            continue;
        };
        let Some(tx) = transcript_of(&v) else {
            continue;
        };
        if !covered.contains(&tx)
            || !genomic_coverable(&case.input)
            || expected_uses_range_reference(expected)
        {
            continue;
        }

        // A `project_to_genomic` decline (#655 — the variant cannot be expressed
        // in the parent frame: a transcript version absent from cdot, an extent
        // exceeding the placed span, a position past the annotated UTR) is out of
        // this hermetic *rendering* gate's scope; the manifest-backed
        // `axis_genomic` retains full coverage of those. Skip rather than fail.
        let g = match projector.project_to_genomic(&v) {
            Ok(g) => g,
            Err(_) => {
                declined += 1;
                continue;
            }
        };
        covered_count += 1;

        let actual = catch_panics(|| -> Result<String, String> {
            let n = normalizer
                .normalize(&g)
                .map_err(|e| format!("normalize: {e}"))?;
            Ok(format!("{n}"))
        });
        t.record(case, expected, actual);
    }

    assert!(
        covered_count > 0,
        "gate_genomic_snapshot: fixture covered zero genomic rows — the coverage \
         filter or the committed window fixture is broken"
    );
    assert!(
        t.fail.is_empty(),
        "gate_genomic_snapshot: {} hermetic divergence(s) from cases.json on \
         {covered_count} window-covered genomic row(s) (real bases + placements, no \
         manifest):\n{}",
        t.fail.len(),
        t.fail
            .iter()
            .take(20)
            .map(|(i, d)| format!("  {i} | {d}"))
            .collect::<Vec<_>>()
            .join("\n"),
    );
    println!(
        "gate_genomic_snapshot: {covered_count} window-covered genomic row(s) pass \
         hermetically (no manifest); {} known_bug, {} improvement; {declined} \
         declined-out-of-scope (#655, retained on the manifest axis)",
        t.known_bug.len(),
        t.improvement.len(),
    );
}

/// Strip mutalyzer's non-standard `NM_x(NP_y):p.…` genomic-context wrapper to the
/// spec-correct bare `NP_y:p.…` ferro emits (#498). A value that is not wrapped
/// (the inner token before `:` is not an `NP_`/`XP_` accession) is returned
/// unchanged — including a bare `NP_y:p.(…)` whose first `(` belongs to the
/// consequence.
fn to_bare_protein(desc: &str) -> String {
    if let Some(open) = desc.find('(') {
        if let Some(close_rel) = desc[open..].find(')') {
            let inner = &desc[open + 1..open + close_rel];
            let rest = &desc[open + close_rel + 1..];
            if inner.starts_with("NP_") || inner.starts_with("XP_") {
                return format!("{inner}{rest}");
            }
        }
    }
    desc.to_string()
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
            known_bug: None,
            improvement: None,
            reference_unavailable: None,
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
        assert_eq!(sc.section, SpecSection::Prioritization);
        assert_eq!(sc.note.as_deref(), Some("dup over ins"));
    }

    // (3b) The `SubstitutionConsequence` section deserializes from its wire
    // string and routes the cited row to `spec_overridden` when ferro's value
    // differs (the c.41A>C substitution-not-frameshift arbitration). Pinning
    // the rename here keeps the `cases.json` byte sequence stable.
    #[test]
    fn case_parses_with_substitution_consequence_section() {
        let case = parse_case(
            r#"{
                "input": "NG_009299.1(NM_017668.3):c.41A>C",
                "protein_description": "NP_060138.1:p.(Glu14Ala)",
                "spec_citation": {
                    "axis": "protein_description",
                    "section": "HGVS §Substitution (no frameshift)",
                    "note": "a substitution cannot be a frameshift"
                }
            }"#,
        );
        let sc = case.spec_citation.as_ref().expect("spec_citation present");
        assert_eq!(sc.axis, Axis::ProteinDescription);
        assert_eq!(sc.section, SpecSection::SubstitutionConsequence);

        // A string mismatch on the cited axis routes to spec_overridden, not FAIL.
        let mut t = AxisTally::new(Axis::ProteinDescription);
        t.record(
            &case,
            "NP_060138.1:p.(Glu14Ala)",
            Ok("NP_060138.1:p.(Asp14Ala)".to_string()),
        );
        assert_eq!(t.spec_overridden.len(), 1);
        assert!(t.fail.is_empty());
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
                cluster: None,
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
                cluster: None,
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
                cluster: None,
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

    // (7b) Under-emission `Err` with `accepted_divergence` on `Axis::Infos` is
    // STILL FAIL. The infos-axis accepted-divergence path in `axis_infos` only
    // silences over-emission errors (those starting with `"ferro emitted extra
    // info"`); under-emission, parse, and normalize errors must still surface
    // as real failures even when the case carries the annotation.
    //
    // On every axis EXCEPT errors, `AxisTally::record` never silences `Err`
    // results — this test mirrors `tally_err_with_accepted_divergence_still_fails`
    // for `Axis::Infos`. (The errors axis is the deliberate exception, where a
    // non-panic `Err` is the expected divergence signal; see the errors-axis
    // tests below.)
    #[test]
    fn tally_infos_under_emission_err_with_accepted_divergence_still_fails() {
        let mut t = AxisTally::new(Axis::Infos);
        let case = make_case(
            "in",
            Some(AcceptedDivergence {
                axis: Axis::Infos,
                policy: Policy::ShuffleAppliedCompoundAllele499,
                note: None,
                cluster: None,
            }),
            None,
        );
        t.record(
            &case,
            "SHUFFLE_APPLIED",
            Err("ferro infos [] missing expected code \"SHUFFLE_APPLIED\" x1 (got 0)".to_string()),
        );
        assert_eq!(t.pass, 0);
        assert!(
            t.divergence_accepted.is_empty(),
            "under-emission Err must not silence into divergence_accepted bucket"
        );
        assert_eq!(t.fail.len(), 1);
    }

    // (7d) Errors axis: a non-panic `Err` (ferro emitted no error, or a
    // different error class than the expected mutalyzer code) WITH an
    // `accepted_divergence` on the errors axis buckets as divergence_accepted —
    // not FAIL. This is the errors-axis exception to the Err-is-always-FAIL rule
    // the projection axes follow (#486).
    #[test]
    fn tally_errors_axis_err_with_accepted_divergence_buckets() {
        let mut t = AxisTally::new(Axis::Errors);
        let case = make_case(
            "in",
            Some(AcceptedDivergence {
                axis: Axis::Errors,
                policy: Policy::ParseTimeRejectionTaxonomy486,
                note: None,
                cluster: None,
            }),
            None,
        );
        t.record(
            &case,
            "EINTRONIC",
            Err("ferro produced no error; mutalyzer expected [\"EINTRONIC\"]".to_string()),
        );
        assert_eq!(t.pass, 0);
        assert_eq!(t.fail.len(), 0, "errors-axis divergence must not FAIL");
        assert_eq!(
            t.divergence_accepted,
            vec![(
                "in".to_string(),
                "ferro-policy-486-parse-time-rejection-taxonomy".to_string()
            )]
        );
    }

    // (7e) Errors axis: a genuine PANIC (`"panic: …"`) is STILL a hard FAIL even
    // with an errors-axis `accepted_divergence` — real crashes are never masked.
    #[test]
    fn tally_errors_axis_panic_still_fails_despite_annotation() {
        let mut t = AxisTally::new(Axis::Errors);
        let case = make_case(
            "in",
            Some(AcceptedDivergence {
                axis: Axis::Errors,
                policy: Policy::ParseTimeRejectionTaxonomy486,
                note: None,
                cluster: None,
            }),
            None,
        );
        t.record(&case, "EINTRONIC", Err("panic: boom".to_string()));
        assert_eq!(t.pass, 0);
        assert!(
            t.divergence_accepted.is_empty(),
            "a panic must not bucket into divergence_accepted"
        );
        assert_eq!(t.fail.len(), 1);
    }

    // (7f) Errors axis: if ferro NOW emits the expected tag (the annotated
    // divergence is gone), the row is an XPASS and FAILs loudly so the stale
    // annotation gets cleaned up — same guard as every other axis.
    #[test]
    fn tally_errors_axis_xpass_fails_to_flag_stale_annotation() {
        let mut t = AxisTally::new(Axis::Errors);
        let case = make_case(
            "in",
            Some(AcceptedDivergence {
                axis: Axis::Errors,
                policy: Policy::ParseTimeRejectionTaxonomy486,
                note: None,
                cluster: None,
            }),
            None,
        );
        t.record(&case, "EINTRONIC", Ok("EINTRONIC".to_string()));
        assert_eq!(t.pass, 0);
        assert!(t.divergence_accepted.is_empty());
        assert_eq!(
            t.fail.len(),
            1,
            "stale errors-axis annotation must XPASS-FAIL"
        );
        assert!(t.fail[0].1.contains("XPASS"));
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

    // (7c) `reference_unavailable` deserializes including its closed `reason`
    // enum and optional `note`.
    #[test]
    fn case_parses_with_reference_unavailable() {
        let case = parse_case(
            r#"{
                "input": "NM_003002.2:c.273del",
                "normalized": "NM_003002.2:c.274del",
                "reference_unavailable": {
                    "axis": "normalized",
                    "reason": "accession-version-absent",
                    "tracking_issue": 672,
                    "note": "NM_003002.2 absent from the prepared manifest; ferro no-ops"
                }
            }"#,
        );
        let ru = case
            .reference_unavailable
            .as_ref()
            .expect("reference_unavailable present");
        assert_eq!(ru.axis, Axis::Normalized);
        assert_eq!(
            ru.reason,
            ReferenceUnavailableReason::AccessionVersionAbsent
        );
        assert_eq!(ru.tracking_issue, 672);
        assert!(ru.note.as_deref().unwrap().contains("no-ops"));
    }

    // (7d) A typo'd `reason` is rejected at deserialization — the closed enum
    // turns an ad-hoc string into a hard parse error (mirrors the axis/policy/
    // section typo guards, #397 item 2).
    #[test]
    fn case_reference_unavailable_reason_typo_rejected_at_deserialization() {
        let result: Result<Case, _> = serde_json::from_str(
            r#"{
                "input": "x",
                "reference_unavailable": {
                    "axis": "normalized",
                    "reason": "version-missing",
                    "tracking_issue": 672
                }
            }"#,
        );
        assert!(
            result.is_err(),
            "expected unknown reason (`version-missing`) to be rejected; got {:?}",
            result.ok(),
        );
    }

    // (7e) A no-op mismatch on the annotated axis routes to the
    // `reference_unavailable` bucket — NOT `fail` and NOT `known_bug` (ferro
    // is not wrong; the oracle simply lacks the reference). ferro echoes the
    // input (`Ok`), so the bucket is reached via the `actual.is_ok()` path.
    #[test]
    fn tally_reference_unavailable_routes_to_its_bucket() {
        let mut t = AxisTally::new(Axis::Normalized);
        let case = parse_case(
            r#"{
                "input": "NM_003002.2:c.273del",
                "reference_unavailable": {
                    "axis": "normalized",
                    "reason": "accession-version-absent",
                    "tracking_issue": 672
                }
            }"#,
        );
        // expected = mutalyzer's normalized output; actual = ferro echoing input.
        t.record(
            &case,
            "NM_003002.2:c.274del",
            Ok("NM_003002.2:c.273del".to_string()),
        );
        assert_eq!(t.pass, 0);
        assert!(t.known_bug.is_empty(), "must not land in known_bug");
        assert!(t.fail.is_empty(), "must not land in fail");
        assert_eq!(
            t.reference_unavailable,
            vec![(
                "NM_003002.2:c.273del".to_string(),
                "accession-version-absent #672".to_string()
            )]
        );
    }

    // (7f) XPASS guard: if a `reference_unavailable` row now MATCHES, the
    // reference became available (the provisioning landed) — fail loudly so the
    // annotation is removed and the row demoted.
    #[test]
    fn tally_reference_unavailable_xpass_fails() {
        let mut t = AxisTally::new(Axis::Normalized);
        let case = parse_case(
            r#"{
                "input": "NM_003002.2:c.273del",
                "reference_unavailable": {
                    "axis": "normalized",
                    "reason": "accession-version-absent",
                    "tracking_issue": 672
                }
            }"#,
        );
        t.record(
            &case,
            "NM_003002.2:c.274del",
            Ok("NM_003002.2:c.274del".to_string()),
        );
        assert_eq!(t.pass, 0);
        assert!(t.reference_unavailable.is_empty());
        assert_eq!(t.fail.len(), 1);
        assert!(
            t.fail[0].1.contains("XPASS: reference_unavailable"),
            "diagnostic must flag the stale annotation: {}",
            t.fail[0].1
        );
    }

    // (7g) Axis scoping: a `reference_unavailable` whose axis does NOT match the
    // tally's axis is still a FAIL — dispositions are strictly per-axis.
    #[test]
    fn tally_reference_unavailable_axis_mismatch_is_fail() {
        let mut t = AxisTally::new(Axis::Genomic);
        let case = parse_case(
            r#"{
                "input": "in",
                "reference_unavailable": {
                    "axis": "normalized",
                    "reason": "ensembl-absent-from-refseq-manifest",
                    "tracking_issue": 671
                }
            }"#,
        );
        t.record(&case, "X", Ok("Y".to_string()));
        assert!(t.reference_unavailable.is_empty());
        assert_eq!(t.fail.len(), 1);
    }

    // (7h) An `Err` result with `reference_unavailable` on the matching axis is
    // STILL FAIL — the annotation accepts a no-op string difference, not a
    // catastrophic failure (panic / parse / normalize Err). The bucket is
    // reached only inside the `actual.is_ok()` block, so an `Err` must fall
    // through to `fail`; this pins that contract (mirrors the sibling
    // `tally_err_with_{accepted_divergence,known_bug,improvement,spec_citation}_still_fails`).
    #[test]
    fn tally_err_with_reference_unavailable_still_fails() {
        let mut t = AxisTally::new(Axis::Normalized);
        let case = parse_case(
            r#"{
                "input": "in",
                "reference_unavailable": {
                    "axis": "normalized",
                    "reason": "accession-version-absent",
                    "tracking_issue": 672
                }
            }"#,
        );
        t.record(&case, "X", Err("normalize: panic boom".to_string()));
        assert_eq!(t.pass, 0);
        assert!(
            t.reference_unavailable.is_empty(),
            "Err must not silence into the reference_unavailable bucket"
        );
        assert_eq!(t.fail.len(), 1);
        assert!(t.fail[0].1.contains("err=normalize: panic boom"));
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
                section: SpecSection::Prioritization,
                note: None,
                cluster: None,
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
                section: SpecSection::Prioritization,
                note: None,
                cluster: None,
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
                section: SpecSection::Prioritization,
                note: None,
                cluster: None,
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
        // 1 pass, 1 divergence_accepted, 1 reference_unavailable,
        // 1 spec_overridden, 1 fail, 0 skipped.
        t.record(&make_case("a", None, None), "X", Ok("X".to_string()));
        t.record(
            &make_case(
                "b",
                Some(AcceptedDivergence {
                    axis: Axis::Normalized,
                    policy: Policy::GeneSymbolSelector121,
                    note: None,
                    cluster: None,
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
                    section: SpecSection::Prioritization,
                    note: None,
                    cluster: None,
                }),
            ),
            "X",
            Ok("Y".to_string()),
        );
        t.record(
            &parse_case(
                r#"{"input":"e","reference_unavailable":{"axis":"normalized",
                   "reason":"accession-version-absent","tracking_issue":672}}"#,
            ),
            "X",
            Ok("e".to_string()),
        );
        t.record(&make_case("c", None, None), "X", Ok("Y".to_string()));

        let path = PathBuf::from("/tmp/example.txt");
        let line = t.summary_line(&path);
        assert_eq!(
            line,
            "normalized: 1 pass / 1 divergence_accepted / 0 known_bug / 0 improvement / 1 reference_unavailable / 1 spec_overridden / 1 FAIL / 0 empty / 0 skipped \
             (FAIL inputs -> /tmp/example.txt)"
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
                cluster: None,
            }),
            Some(SpecCitation {
                axis: Axis::Normalized,
                section: SpecSection::Prioritization,
                note: None,
                cluster: None,
            }),
        );
        t.record(&case, "X", Ok("Y".to_string()));
        assert_eq!(t.divergence_accepted.len(), 1);
        assert!(t.fail.is_empty());
    }

    // (11) `known_bug` deserializes including the optional `note`. The
    // `axis` field is the closed [`Axis`] enum, so a typo would be rejected
    // at parse time (same contract as `accepted_divergence`).
    #[test]
    fn case_parses_with_known_bug() {
        let case = parse_case(
            r#"{
                "input": "NM_000088.3:c.459A>G",
                "known_bug": {
                    "axis": "normalized",
                    "tracking_issue": 325,
                    "note": "ferro mis-shifts this duplication"
                }
            }"#,
        );
        let kb = case.known_bug.as_ref().expect("known_bug present");
        assert_eq!(kb.axis, Axis::Normalized);
        assert_eq!(kb.tracking_issue, 325);
        assert_eq!(
            kb.note.as_deref(),
            Some("ferro mis-shifts this duplication")
        );
    }

    // (12) A string mismatch on an axis with a matching `known_bug` routes
    // into the `known_bug` (xfail) bucket — non-failing, separately tracked,
    // and never counted as a pass.
    #[test]
    fn tally_known_bug_on_matching_axis() {
        let mut t = AxisTally::new(Axis::Normalized);
        let mut case = make_case("in", None, None);
        case.known_bug = Some(KnownBug {
            axis: Axis::Normalized,
            tracking_issue: 325,
            note: None,
            cluster: None,
        });
        t.record(&case, "X", Ok("Y".to_string()));
        assert_eq!(t.pass, 0);
        assert_eq!(t.known_bug, vec![("in".to_string(), "#325".to_string())]);
        assert!(t.fail.is_empty());
    }

    // (13) XPASS on a `known_bug` axis: ferro now matches mutalyzer, so the
    // bug appears fixed. This must FAIL loudly so the stale annotation and
    // fixed row are cleaned up — not silently counted as a pass.
    #[test]
    fn tally_known_bug_xpass_fails() {
        let mut t = AxisTally::new(Axis::Normalized);
        let mut case = make_case("in", None, None);
        case.known_bug = Some(KnownBug {
            axis: Axis::Normalized,
            tracking_issue: 325,
            note: None,
            cluster: None,
        });
        t.record(&case, "X", Ok("X".to_string()));
        assert_eq!(t.pass, 0);
        assert!(t.known_bug.is_empty());
        assert_eq!(t.fail.len(), 1);
        assert!(t.fail[0].1.contains("XPASS"));
        assert!(t.fail[0].1.contains("now matches"));
    }

    // (14) XPASS on an `accepted_divergence` axis: ferro now matches
    // mutalyzer, so the documented divergence is gone. This must FAIL so the
    // stale annotation is removed (mirrors the known_bug XPASS contract).
    #[test]
    fn tally_accepted_divergence_xpass_fails() {
        let mut t = AxisTally::new(Axis::Normalized);
        let case = make_case(
            "in",
            Some(AcceptedDivergence {
                axis: Axis::Normalized,
                policy: Policy::GeneSymbolSelector121,
                note: None,
                cluster: None,
            }),
            None,
        );
        t.record(&case, "X", Ok("X".to_string()));
        assert_eq!(t.pass, 0);
        assert!(t.divergence_accepted.is_empty());
        assert_eq!(t.fail.len(), 1);
        assert!(t.fail[0].1.contains("XPASS"));
    }

    // (15) `Err` result with `known_bug` on the matching axis is STILL FAIL
    // — the annotation accepts a string difference (an xfail), not a
    // catastrophic failure (panic / parse / normalize Err). Mirrors
    // `tally_err_with_accepted_divergence_still_fails`.
    #[test]
    fn tally_err_with_known_bug_still_fails() {
        let mut t = AxisTally::new(Axis::Normalized);
        let mut case = make_case("in", None, None);
        case.known_bug = Some(KnownBug {
            axis: Axis::Normalized,
            tracking_issue: 325,
            note: None,
            cluster: None,
        });
        t.record(&case, "X", Err("normalize: panic boom".to_string()));
        assert_eq!(t.pass, 0);
        assert!(
            t.known_bug.is_empty(),
            "Err must not silence into known_bug bucket"
        );
        assert_eq!(t.fail.len(), 1);
        assert!(t.fail[0].1.contains("err=normalize: panic boom"));
    }

    // (16) `improvement` deserializes including the optional `note`. Both
    // `axis` and `section` are closed enums, so a typo in either is rejected
    // at parse time (same contract as `known_bug` / `spec_citation`).
    #[test]
    fn case_parses_with_improvement() {
        let case = parse_case(
            r#"{
                "input": "NG_012337.1(TIMM8B_v001):c.12_13insGATC",
                "known_bug": null,
                "improvement": {
                    "axis": "normalized",
                    "tracking_issue": 500,
                    "section": "HGVS §RefSeqGene transcript selection",
                    "note": "ferro preserves the gene-symbol selector; spec prefers the NM_ accession"
                }
            }"#,
        );
        let imp = case.improvement.as_ref().expect("improvement present");
        assert_eq!(imp.axis, Axis::Normalized);
        assert_eq!(imp.tracking_issue, 500);
        assert_eq!(imp.section, SpecSection::RefSeqGeneSelector);
        assert_eq!(
            imp.note.as_deref(),
            Some("ferro preserves the gene-symbol selector; spec prefers the NM_ accession")
        );
    }

    // (16b) `improvement` requires a `section`: omitting it is a hard parse
    // error, so every tracked improvement is substantiated by a spec citation.
    #[test]
    fn case_improvement_without_section_rejected() {
        let result: Result<Case, _> = serde_json::from_str(
            r#"{
                "input": "x",
                "improvement": { "axis": "normalized", "tracking_issue": 500 }
            }"#,
        );
        assert!(
            result.is_err(),
            "expected improvement without a section to be rejected; got {:?}",
            result.ok(),
        );
    }

    // (17) A string mismatch on an axis with a matching `improvement` routes
    // into the `improvement` (xfail) bucket — non-failing, separately tracked,
    // and never counted as a pass.
    #[test]
    fn tally_improvement_on_matching_axis() {
        let mut t = AxisTally::new(Axis::Normalized);
        let mut case = make_case("in", None, None);
        case.improvement = Some(Improvement {
            axis: Axis::Normalized,
            tracking_issue: 500,
            section: SpecSection::RefSeqGeneSelector,
            note: None,
            cluster: None,
        });
        t.record(&case, "X", Ok("Y".to_string()));
        assert_eq!(t.pass, 0);
        assert_eq!(t.improvement, vec![("in".to_string(), "#500".to_string())]);
        assert!(t.fail.is_empty());
    }

    // (18) XPASS on an `improvement` axis: ferro now matches mutalyzer, so it
    // reached the spec-preferred form. This must FAIL loudly so the stale
    // annotation and converged row are cleaned up (mirrors the known_bug /
    // accepted_divergence XPASS contract).
    #[test]
    fn tally_improvement_xpass_fails() {
        let mut t = AxisTally::new(Axis::Normalized);
        let mut case = make_case("in", None, None);
        case.improvement = Some(Improvement {
            axis: Axis::Normalized,
            tracking_issue: 500,
            section: SpecSection::RefSeqGeneSelector,
            note: None,
            cluster: None,
        });
        t.record(&case, "X", Ok("X".to_string()));
        assert_eq!(t.pass, 0);
        assert!(t.improvement.is_empty());
        assert_eq!(t.fail.len(), 1);
        assert!(t.fail[0].1.contains("XPASS"));
        assert!(t.fail[0].1.contains("now matches"));
    }

    // (19) `improvement` on a different axis from the tally's must NOT route
    // to the improvement bucket — annotations are axis-scoped.
    #[test]
    fn tally_improvement_other_axis_is_fail() {
        let mut t = AxisTally::new(Axis::Genomic);
        let mut case = make_case("in", None, None);
        case.improvement = Some(Improvement {
            axis: Axis::Normalized,
            tracking_issue: 500,
            section: SpecSection::RefSeqGeneSelector,
            note: None,
            cluster: None,
        });
        t.record(&case, "X", Ok("Y".to_string()));
        assert_eq!(t.pass, 0);
        assert!(t.improvement.is_empty());
        assert_eq!(t.fail.len(), 1);
        assert_eq!(t.fail[0].0, "in");
    }

    // (20) `Err` result with `improvement` on the matching axis is STILL FAIL
    // — the annotation accepts a string difference (an xfail), not a
    // catastrophic failure (panic / parse / normalize Err). Mirrors
    // `tally_err_with_known_bug_still_fails`.
    #[test]
    fn tally_err_with_improvement_still_fails() {
        let mut t = AxisTally::new(Axis::Normalized);
        let mut case = make_case("in", None, None);
        case.improvement = Some(Improvement {
            axis: Axis::Normalized,
            tracking_issue: 500,
            section: SpecSection::RefSeqGeneSelector,
            note: None,
            cluster: None,
        });
        t.record(&case, "X", Err("normalize: panic boom".to_string()));
        assert_eq!(t.pass, 0);
        assert!(
            t.improvement.is_empty(),
            "Err must not silence into improvement bucket"
        );
        assert_eq!(t.fail.len(), 1);
        assert!(t.fail[0].1.contains("err=normalize: panic boom"));
    }

    // ---- Output-quality / empty-projection signal (#651) -------------------

    // (21) An empty-projection-sentinel `Err` increments `empty_got`, and the
    // sentinel is stripped (and tagged) in the FAIL diagnostic.
    #[test]
    fn tally_empty_projection_sentinel_counts_and_strips() {
        let mut t = AxisTally::new(Axis::ProteinDescription);
        let case = make_case("in", None, None);
        t.record(
            &case,
            "NP_x:p.(Arg1Cys)",
            Err(format!("{EMPTY_PROJECTION_SENTINEL}no protein predicted")),
        );
        assert_eq!(t.empty_got, 1);
        // Still a FAIL (an empty projection is an Err), with a clean diagnostic.
        assert_eq!(t.fail.len(), 1);
        assert!(
            t.fail[0]
                .1
                .contains("err=[empty projection] no protein predicted"),
            "diagnostic should strip the sentinel and tag it: {}",
            t.fail[0].1
        );
        assert!(!t.fail[0].1.contains(EMPTY_PROJECTION_SENTINEL));
    }

    // (22) A non-empty (ordinary) `Err` mismatch does NOT count as empty.
    #[test]
    fn tally_ordinary_err_is_not_counted_empty() {
        let mut t = AxisTally::new(Axis::CodingProteinDescriptions);
        let case = make_case("in", None, None);
        t.record(
            &case,
            "[...]",
            Err("missing pair (\"c\", \"p\"); got []".to_string()),
        );
        assert_eq!(t.empty_got, 0);
        assert_eq!(t.fail.len(), 1);
    }

    // (23) An empty projection on an *annotated* axis is counted as empty even
    // though the annotation would otherwise silence a string mismatch — and an
    // `Err` is never silenced into the annotation bucket, so it still FAILs.
    // This is the core #651 guarantee: a populated→empty degradation surfaces.
    #[test]
    fn tally_empty_projection_on_annotated_row_still_surfaces() {
        let mut t = AxisTally::new(Axis::CodingProteinDescriptions);
        let mut case = make_case("in", None, None);
        case.known_bug = Some(KnownBug {
            axis: Axis::CodingProteinDescriptions,
            tracking_issue: 326,
            note: None,
            cluster: None,
        });
        t.record(
            &case,
            "[(\"c\", \"p\")]",
            Err(format!(
                "{EMPTY_PROJECTION_SENTINEL}project_variant_all produced no coding/protein pairs"
            )),
        );
        assert_eq!(t.empty_got, 1);
        assert!(
            t.known_bug.is_empty(),
            "empty-projection Err must not be silenced into the known_bug bucket"
        );
        assert_eq!(t.fail.len(), 1);
    }

    // (24) The budget verdict: at-or-below baseline passes, above fails.
    #[test]
    fn empty_budget_verdict_gates_only_on_rise() {
        assert!(empty_budget_verdict(0, 0).is_ok());
        assert!(empty_budget_verdict(5, 5).is_ok(), "equal = unchanged");
        assert!(empty_budget_verdict(3, 5).is_ok(), "below = improvement");
        let err = empty_budget_verdict(6, 5).expect_err("rise must fail");
        assert!(err.contains("5 -> 6"), "message names the rise: {err}");
    }

    // ------------------------------------------------------------------------
    // Snapshot-coverage string-parser unit tests
    // ------------------------------------------------------------------------
    //
    // The snapshot gates exercise these byte-level helpers only indirectly via
    // the corpus. These tests pin the named edge cases from each helper's doc
    // comment directly against its behavior.

    // (25) Intron-offset detection: a `+`/`-` between two digits flags an
    // intronic position; a leading `-` (UTR) or `*` (3'UTR) marker does not.
    #[test]
    fn has_intron_offset_flags_intronic_and_keeps_utr() {
        // Intronic: digit `+`/`-` digit.
        assert!(has_intron_offset("52+1"), "c.52+1 is intronic");
        assert!(has_intron_offset("378-17"), "c.378-17 is intronic");
        assert!(has_intron_offset("*824+10"), "c.*824+10 is intronic");
        // UTR / exonic: the `+`/`-` is not flanked by two digits.
        assert!(
            !has_intron_offset("-1"),
            "c.-1 is a 5'UTR position, not intronic"
        );
        assert!(
            !has_intron_offset("*1"),
            "c.*1 is a 3'UTR position, not intronic"
        );
        assert!(!has_intron_offset("273del"), "c.273del has no offset");
    }

    // (26) Delins coordinate-range detection: `delins` immediately followed by
    // a digit is a range delins; `delins` followed by literal bases is not.
    #[test]
    fn delins_coordinate_range_distinguishes_range_from_literal() {
        assert!(
            delins_coordinate_range("206_210delins190_220inv"),
            "delins followed by a digit references a coordinate range"
        );
        assert!(
            !delins_coordinate_range("100delinsATC"),
            "delins followed by literal bases is not a coordinate range"
        );
        assert!(
            !delins_coordinate_range("273del"),
            "plain del has no delins"
        );
    }

    // (27) `snapshot_coverable` composes the above filters: a bare exonic `c.`
    // [`HgvsVariant::Cds`] row is coverable; non-`Cds`, parenthesized,
    // intronic, and range-delins rows are excluded.
    #[test]
    fn snapshot_coverable_covers_bare_exonic_cds_only() {
        let coverable = parse_hgvs("NM_003002.2:c.273del").expect("parse c.273del");
        assert!(
            snapshot_coverable("NM_003002.2:c.273del", &coverable),
            "bare exonic c. on a coding transcript is snapshot-coverable"
        );

        // Intronic offset → excluded (needs multi-exon structure).
        let intronic = parse_hgvs("NM_003002.2:c.274+20C>T").expect("parse intronic");
        assert!(
            !snapshot_coverable("NM_003002.2:c.274+20C>T", &intronic),
            "intronic offset rows are excluded"
        );

        // Range delins → excluded (needs reference context to render the segment).
        let range_delins =
            parse_hgvs("NM_003002.4:c.206_210delins190_220inv").expect("parse range delins");
        assert!(
            !snapshot_coverable("NM_003002.4:c.206_210delins190_220inv", &range_delins),
            "coordinate-range delins rows are excluded"
        );

        // Non-`Cds` (genomic) → excluded regardless of the coordinate.
        let genomic = parse_hgvs("NC_012920.1:g.3243A>G").expect("parse genomic");
        assert!(
            !matches!(genomic, HgvsVariant::Cds(_)),
            "fixture genomic input must not be a Cds variant"
        );
        assert!(
            !snapshot_coverable("NC_012920.1:g.3243A>G", &genomic),
            "non-Cds variants are excluded"
        );
    }

    // (28) `to_bare_protein` strips mutalyzer's `NM_x(NP_y):p.…` wrapper to the
    // bare `NP_y:p.…` and leaves an already-bare `NP_y:p.(…)` unchanged.
    #[test]
    fn to_bare_protein_unwraps_wrapper_and_preserves_bare() {
        // Wrapped: first `(` encloses the NP_ accession.
        assert_eq!(
            to_bare_protein("NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGlu)"),
            "NP_478102.2:p.(Asp68_Gly69delinsGlu)",
        );
        // Already bare: the first `(` belongs to the `p.(…)` consequence, whose
        // inner token does not start with `NP_`/`XP_`, so it is returned as-is.
        assert_eq!(
            to_bare_protein("NP_002993.1:p.(Asp92ThrfsTer43)"),
            "NP_002993.1:p.(Asp92ThrfsTer43)",
        );
        // No parentheses at all → returned unchanged.
        assert_eq!(to_bare_protein("NP_002993.1:p.="), "NP_002993.1:p.=");
    }

    // (29) `genomic_coverable` gates a row onto the hermetic genomic axis only
    // when it is a bare exonic `c.` variant on an `NM_`/`NR_` selector under a
    // genomic parent. Pins the load-bearing inclusion/exclusion edges from the
    // helper's doc comment directly.
    #[test]
    fn genomic_coverable_covers_bare_exonic_nm_selector_only() {
        // Positive: NG_ parent, NM_ selector, bare exonic c. del.
        assert!(
            genomic_coverable("NG_012337.1(NM_003002.2):c.1del"),
            "bare exonic c. with an NM_ selector on a genomic parent is coverable"
        );

        // #745 carve-out: NG_008835.* is explicitly excluded.
        assert!(
            !genomic_coverable("NG_008835.1(NM_022124.6):c.3304_3305del"),
            "NG_008835.* is the #745 homopolymer carve-out and is excluded"
        );

        // Gene-symbol selector (not NM_/NR_) → excluded.
        assert!(
            !genomic_coverable("NG_012772.1(BRCA2_v001):c.1del"),
            "a gene-symbol selector lacks the required NM_/NR_ prefix"
        );

        // Intronic offset → excluded.
        assert!(
            !genomic_coverable("NG_012337.1(NM_003002.2):c.52+1del"),
            "intronic offset positions are excluded"
        );

        // Compound `[...;...]` → excluded.
        assert!(
            !genomic_coverable("NG_012772.1(NM_000059.3):c.[632-?_681+?del;681+4del]"),
            "compound allele rows are excluded"
        );
    }

    // (30) `expected_uses_range_reference` detects an expected genomic output that
    // renders the inserted segment as a coordinate range / cross-reference
    // (`ins[…]` or `ins<digit>…`) rather than literal bases.
    #[test]
    fn expected_uses_range_reference_detects_coordinate_reference() {
        // Bracketed cross-reference insertion.
        assert!(
            expected_uses_range_reference("NG_008939.1:g.5207_5212delins[4300_4320]"),
            "ins[…] is a bracketed coordinate cross-reference"
        );
        // `ins` immediately followed by a digit (coordinate range).
        assert!(
            expected_uses_range_reference("NG_008939.1:g.5207_5212delins4300_4320"),
            "ins followed by a digit references a coordinate range"
        );
        // Literal inserted bases → not a coordinate reference.
        assert!(
            !expected_uses_range_reference("NG_008939.1:g.100_110delinsAAA"),
            "ins followed by literal bases is not a coordinate reference"
        );
    }
}
