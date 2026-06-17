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
    // Selection coverage (reported metric, not a pass/fail gate): of all
    // applicable cases, how many had ferro return the expected *base*
    // transcript accession at all — independent of whether the consequence
    // matched. Quantifies the different-base ("ferro picked a different
    // transcript") portion of the divergences as a tracked number.
    selection_hits: usize,
    selection_total: usize,
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
            selection_hits: 0,
            selection_total: 0,
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

    /// Record a discovery-axis case with structural data-source quarantine.
    ///
    /// `expected_repr`/`actual` are the same display/disposition pair the
    /// [`AxisTally::record`] path uses (for PASS/XPASS and per-case annotation
    /// routing). `expected_entries` are the individual expected rendered
    /// strings (one for single-value axes, several for the pair/list axes) and
    /// `returned` is the full set of rendered strings ferro produced for this
    /// case.
    ///
    /// When the case is a genuine divergence (`actual` is `Ok` but did not
    /// match, or `Err`), every *diverging* expected entry is classified via
    /// [`classify_divergence`] and the case is routed **structurally** — never
    /// by a per-case annotation:
    ///
    /// - All diverging entries are [`Divergence::SelectionMiss`] →
    ///   `divergence_accepted` (cluster `transcript-selection-vs-uta`).
    /// - All diverging entries are quarantinable and at least one is a
    ///   [`Divergence::CoordinateSkew`] on a [`SOURCE_SKEW_TRANSCRIPTS`]
    ///   transcript → `divergence_accepted` (cluster `alignment-source-skew`).
    /// - Otherwise (a [`Divergence::FormOnly`], or a `CoordinateSkew` on a
    ///   transcript NOT on the gate list, or an `Err` projection failure) →
    ///   FAIL, so genuine convention/algorithm signal and gate-list gaps are
    ///   surfaced rather than silently accepted.
    ///
    /// A clean match (or a stale per-case annotation XPASS) delegates to
    /// [`AxisTally::record`] unchanged.
    fn record_classified(
        &mut self,
        case: &Case,
        expected_repr: &str,
        actual: Result<String, String>,
        expected_entries: &[String],
        returned: &[String],
    ) {
        // The discovery-axis closures signal a clean match by returning
        // `Ok(expected_repr)` and a divergence by returning `Err(...)` (a
        // missing-entry diagnostic). A clean match (or a stale-annotation
        // XPASS) is handled by the existing buckets; anything else is a
        // candidate for structural quarantine.
        let is_match = matches!(&actual, Ok(s) if s == expected_repr);
        if is_match {
            self.record(case, expected_repr, actual);
            return;
        }

        // A hard projection failure (parse / project / panic) is always a real
        // FAIL — it is not a data-source divergence. The discovery-axis closures
        // signal a value divergence with a `missing …` diagnostic; any other
        // `Err` is a hard failure that must surface via the normal path.
        if let Err(e) = &actual {
            if !e.starts_with("missing ") {
                self.record(case, expected_repr, actual);
                return;
            }
        }

        // Classify each expected entry that did not match against ferro's set.
        let diverging: Vec<Divergence> = expected_entries
            .iter()
            .map(|e| classify_divergence(e, returned))
            .filter(|d| *d != Divergence::Match)
            .collect();

        // Defensive: if nothing classified as diverging (shouldn't happen given
        // `is_divergence`), fall back to the standard path.
        if diverging.is_empty() {
            self.record(case, expected_repr, actual);
            return;
        }

        let all_selection_miss = diverging.iter().all(|d| *d == Divergence::SelectionMiss);
        // A coordinate skew is quarantinable only on a gate-listed transcript.
        let listed_coordinate_skew = |entry: &str| -> bool {
            base_accession(entry).is_some_and(|b| SOURCE_SKEW_TRANSCRIPTS.contains(&b.as_str()))
        };
        let all_quarantinable = expected_entries
            .iter()
            .map(|e| (e, classify_divergence(e, returned)))
            .all(|(e, d)| match d {
                Divergence::Match | Divergence::SelectionMiss => true,
                Divergence::CoordinateSkew => listed_coordinate_skew(e),
                Divergence::FormOnly => false,
            });
        let has_listed_skew = expected_entries.iter().any(|e| {
            classify_divergence(e, returned) == Divergence::CoordinateSkew
                && listed_coordinate_skew(e)
        });

        if all_selection_miss {
            self.divergence_accepted
                .push((case.input.clone(), CLUSTER_TRANSCRIPT_SELECTION.to_string()));
        } else if all_quarantinable && has_listed_skew {
            self.divergence_accepted.push((
                case.input.clone(),
                CLUSTER_ALIGNMENT_SOURCE_SKEW.to_string(),
            ));
        } else {
            // FormOnly, or CoordinateSkew on a non-listed transcript (gate gap),
            // or a mix that isn't cleanly one source-skew category: surface it.
            let diag = match actual {
                Ok(got) => format!("expected={expected_repr:?} got={got:?}"),
                Err(e) => format!("expected={expected_repr:?} err={e}"),
            };
            self.fail.push((case.input.clone(), diag));
        }
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
        println!(
            "{}: selection coverage {}/{} (ferro returned the expected base transcript)",
            self.axis, self.selection_hits, self.selection_total
        );

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

/// Extract the *bare* base accession from a rendered variant string.
///
/// Takes the substring before the first `:`, drops a trailing `(...)`
/// gene-symbol group if present (ferro renders e.g. `NM_000682.7(ADRA2B)`),
/// then drops the `.N` version suffix. Returns the bare base accession, e.g.
/// both `NM_000682.7(ADRA2B):c.5A>T` and `NM_000682.5:c.5A>T` map to
/// `Some("NM_000682")`.
///
/// Returns `None` when the input has no `:` (malformed — no accession/
/// consequence boundary).
fn base_accession(rendered: &str) -> Option<String> {
    let (accession, _consequence) = rendered.split_once(':')?;
    // Drop a trailing `(GENE)` gene-symbol group if present.
    let accession = match accession.split_once('(') {
        Some((before_paren, _gene)) => before_paren,
        None => accession,
    };
    // Drop the `.N` version suffix if present.
    let base = match accession.rsplit_once('.') {
        Some((base, _version)) => base,
        None => accession,
    };
    Some(base.to_string())
}

/// Return everything after the first `:` in a rendered variant string — the
/// `c.`/`n.`/`p.` consequence body. Returns `""` when there is no `:`.
fn consequence_part(rendered: &str) -> &str {
    match rendered.split_once(':') {
        Some((_accession, consequence)) => consequence,
        None => "",
    }
}

/// Version-insensitive consequence comparison: `expected` and `actual` match
/// iff they share the same base accession (version digit and `(GENE)` suffix
/// stripped) AND the same consequence body. Only the version digit and gene
/// symbol are forgiven — any coordinate or edit difference is a mismatch.
fn consequence_matches(expected: &str, actual: &str) -> bool {
    match (base_accession(expected), base_accession(actual)) {
        (Some(expected_base), Some(actual_base)) => {
            expected_base == actual_base && consequence_part(expected) == consequence_part(actual)
        }
        _ => false,
    }
}

/// Structural detector for the **predicted-parenthesization** improvement: is
/// `actual` exactly `expected` with the `p.`-token wrapped in `(...)`?
///
/// Returns `true` iff `expected` and `actual` share the same base accession
/// (version digit and `(GENE)` suffix forgiven, like [`consequence_matches`])
/// and `actual`'s consequence is `expected`'s consequence with the body after
/// the `p.` datum prefix wrapped in parentheses — e.g. `p.Arg268Trp` →
/// `p.(Arg268Trp)`. The inner token must be byte-identical; a different inner
/// token (`p.(Arg268Cys)`), an already-wrapped `expected`, or a base-accession
/// mismatch all return `false`.
///
/// This is ferro's spec-correct form: a DNA-predicted protein consequence must
/// be parenthesized (`docs/recommendations/protein/substitution.md` — predicted
/// consequences are uncertain and rendered `p.(…)`). hgvs-rs emits the bare
/// legacy form. The classifier routes these to `improvement` rather than FAIL.
fn is_predicted_parenthesization(expected: &str, actual: &str) -> bool {
    let (Some(expected_base), Some(actual_base)) =
        (base_accession(expected), base_accession(actual))
    else {
        return false;
    };
    if expected_base != actual_base {
        return false;
    }
    let expected_consequence = consequence_part(expected);
    let actual_consequence = consequence_part(actual);
    // Split each consequence into its `p.` datum prefix and the body token.
    let (Some((expected_datum, expected_body)), Some((actual_datum, actual_body))) = (
        expected_consequence.split_once('.'),
        actual_consequence.split_once('.'),
    ) else {
        return false;
    };
    // The datum type (`p`) must agree; only the body wrapping may differ.
    if expected_datum != actual_datum {
        return false;
    }
    // The bare body must not already be wrapped (no wrap to add), and the
    // wrapped body must be exactly the bare body inside one `(...)` group.
    if expected_body.starts_with('(') {
        return false;
    }
    actual_body == format!("({expected_body})")
}

/// Why a case diverges from the corpus expectation, used to route structural
/// data-source skew into `divergence_accepted` while keeping genuine
/// algorithm/convention deltas as FAILs.
///
/// - [`Divergence::SelectionMiss`] — ferro never returned the expected base
///   transcript accession (transcript selection/absence: ferro's cdot set vs
///   the 2021 UTA snapshot).
/// - [`Divergence::CoordinateSkew`] — ferro returned the expected base
///   transcript but at a *different coordinate* (the `c.`/`n.` position prefix
///   differs). This is the alignment-source skew (RefSeq-GFF vs UTA splign).
/// - [`Divergence::FormOnly`] — ferro returned the expected base AND the same
///   coordinate position, but a different *edit form* (e.g. `dup` vs `ins`).
///   This is NOT data-source — it is a genuine convention/algorithm delta and
///   stays a FAIL.
/// - [`Divergence::Match`] — some actual consequence matched the expectation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Divergence {
    SelectionMiss,
    CoordinateSkew,
    FormOnly,
    Match,
}

/// Extract the *position prefix* of a `c.`/`n.` (or `p.`) consequence body: the
/// leading coordinate token(s) before the edit operation. For a consequence
/// body like `c.*1865G>T` this returns `*1865`; for `c.*1353+1856G>T` it
/// returns `*1353+1856`; for an `c.10+830_10+831insG` range it returns
/// `10+830_10+831`.
///
/// The position prefix is the maximal leading run of *coordinate* characters:
/// digits, the UTR/offset sign markers `*`/`-`/`+`, and the range separator
/// `_`. The first character that is none of those (an edit letter `A`–`Z`, a
/// `>` substitution arrow, or an operation keyword like `del`/`ins`/`dup`)
/// terminates the prefix. The leading datum-type letter and `.` (the `c.`/`n.`
/// in `c.*1865`) are skipped first.
fn position_prefix(consequence: &str) -> &str {
    // Drop the `c.`/`n.`/`p.`/`r.`/`m.`/`g.` datum-type prefix (`<letter>.`).
    let body = match consequence.split_once('.') {
        Some((_datum, rest)) => rest,
        None => consequence,
    };
    let end = body
        .char_indices()
        .find(|(_, c)| !(c.is_ascii_digit() || matches!(c, '*' | '-' | '+' | '_')))
        .map(|(i, _)| i)
        .unwrap_or(body.len());
    &body[..end]
}

/// Classify why `expected` diverges from `actual_set` (the rendered strings
/// ferro returned for this case), for structural quarantine routing.
///
/// See [`Divergence`] for the meaning of each variant. The ordering is
/// significant: `Match` short-circuits; then a missing base accession is a
/// `SelectionMiss`; then a position-prefix difference is `CoordinateSkew`;
/// otherwise the base + position matched and only the edit form differs, which
/// is `FormOnly`.
fn classify_divergence(expected: &str, actual_set: &[String]) -> Divergence {
    if actual_set.iter().any(|a| consequence_matches(expected, a)) {
        return Divergence::Match;
    }
    if !set_has_base(actual_set, expected) {
        return Divergence::SelectionMiss;
    }
    // Base accession is present. Compare the position prefix of the expected
    // consequence against every actual that shares the expected base accession.
    let Some(want_base) = base_accession(expected) else {
        return Divergence::SelectionMiss;
    };
    let want_prefix = position_prefix(consequence_part(expected));
    let same_position = actual_set.iter().any(|a| {
        base_accession(a).as_deref() == Some(want_base.as_str())
            && position_prefix(consequence_part(a)) == want_prefix
    });
    if same_position {
        Divergence::FormOnly
    } else {
        Divergence::CoordinateSkew
    }
}

/// Base accessions where ferro's RefSeq-GFF alignment is known (per the gate;
/// see `tests/fixtures/hgvs-rs-projection/ALIGNMENT_SOURCE_DIVERGENCE.md`) to
/// disagree with the 2021 UTA splign snapshot the corpus expectations come
/// from. A [`Divergence::CoordinateSkew`] on one of these is structural
/// data-source skew and is auto-quarantined as `divergence_accepted`; a
/// coordinate skew on any *other* transcript stays a FAIL (the gate list is
/// incomplete — surface it, don't silently accept).
///
/// Two sub-kinds, both alignment-source skew (see the provenance doc):
/// - **gapped-CIGAR skew**: the UTA splign CIGAR carries a within-exon `D`/`I`
///   gap or an exon-boundary off-by-one, so the `c.`/`n.` coordinate shifts
///   downstream of the gap.
/// - **genomic-anchor skew**: the per-exon CIGARs are ungapped *and* ferro's
///   CDS start agrees with UTA's `cds_start_i`, yet a uniform `+2` `c.`-offset
///   remains — a 2 bp difference in where the two sources anchor the transcript
///   on the genome. Confirmed for the last three entries below.
const SOURCE_SKEW_TRANSCRIPTS: &[&str] = &[
    // Gapped-CIGAR skew (within-exon indel or boundary off-by-one).
    "NM_001080519",
    "NM_001077527",
    "NM_000804",
    "NM_000682",
    "NM_003777",
    "NM_006158",
    "NM_001277115",
    // Genomic-anchor skew (ungapped CIGAR, CDS start agrees with UTA, uniform
    // +2 c.-offset from a 2 bp genomic alignment-anchor difference).
    "NM_020451",
    "NM_007199",
    "NM_002386",
];

/// Cluster id for the alignment-source skew (coordinate differs, base matches,
/// transcript is on the gate list).
const CLUSTER_ALIGNMENT_SOURCE_SKEW: &str = "alignment-source-skew";

/// Cluster id for transcript selection/absence (expected base never returned).
const CLUSTER_TRANSCRIPT_SELECTION: &str = "transcript-selection-vs-uta";

/// Cluster id for the predicted-parenthesization improvement (ferro emits the
/// spec-preferred `p.(…)` predicted form; hgvs-rs emits the bare legacy form).
const CLUSTER_PREDICTED_PARENTHESIZATION: &str = "predicted-parenthesization";

/// Detect the **predicted-parenthesization** improvement on a `(c., p.)`-pair
/// axis: for every expected `[c, p]` pair, ferro returned a pair whose coding
/// matches version-insensitively ([`consequence_matches`]) AND whose protein is
/// exactly the expected protein with the predicted-consequence parentheses
/// added ([`is_predicted_parenthesization`]).
///
/// When this holds, the only divergence is that ferro emits the spec-preferred
/// parenthesized predicted form while the corpus carries the bare legacy form —
/// a tracked `improvement`, not a FAIL. Returns `false` if any expected pair is
/// unmatched on the coding component or differs by anything other than the
/// protein parenthesization (so genuine coordinate/selection/form deltas still
/// route through the normal classifier).
fn is_predicted_parenthesization_pairs(
    expected_pairs: &[Vec<String>],
    returned_pairs: &[(String, String)],
) -> bool {
    let applicable: Vec<&Vec<String>> = expected_pairs.iter().filter(|p| p.len() == 2).collect();
    if applicable.is_empty() {
        return false;
    }
    applicable.iter().all(|pair| {
        let want_c = &pair[0];
        let want_p = &pair[1];
        returned_pairs.iter().any(|(c, p)| {
            consequence_matches(want_c, c) && is_predicted_parenthesization(want_p, p)
        })
    })
}

/// Version-insensitive set membership on base accession: does any rendered
/// string in `rendered_set` share the same base accession (version digit and
/// `(GENE)` suffix stripped) as `expected`? Used for the selection-coverage
/// metric — "did ferro return the expected base transcript at all", regardless
/// of whether its consequence matched.
fn set_has_base(rendered_set: &[String], expected: &str) -> bool {
    let Some(want) = base_accession(expected) else {
        return false;
    };
    rendered_set
        .iter()
        .any(|r| base_accession(r).as_deref() == Some(want.as_str()))
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
            "hgvs-rs-projection: skipping — no manifest at FERRO_MANIFEST \
             or benchmark-output/manifest.json"
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

        // Capture ferro's returned coding accessions for the selection-coverage
        // metric — keyed on the coding (c.) transcript, the projection target.
        // `returned_pairs` carries the full (coding, protein) set so the
        // predicted-parenthesization improvement can be detected on the protein
        // component (which `record_classified` does not inspect).
        let mut returned_coding: Vec<String> = Vec::new();
        let mut returned_pairs: Vec<(String, String)> = Vec::new();
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
            returned_coding = actual_pairs.iter().map(|(c, _)| c.clone()).collect();
            returned_pairs = actual_pairs.clone();

            for pair in pairs {
                if pair.len() != 2 {
                    continue;
                }
                let want_c = &pair[0];
                let want_p = &pair[1];
                if !actual_pairs
                    .iter()
                    .any(|(c, p)| consequence_matches(want_c, c) && consequence_matches(want_p, p))
                {
                    return Err(format!(
                        "missing pair ({want_c:?}, {want_p:?}); got {actual_pairs:?}"
                    ));
                }
            }
            Ok(expected_repr.clone())
        });

        // Selection coverage: ferro returned the expected base coding
        // transcript for every expected (coding, protein) pair.
        t.selection_total += 1;
        if pairs
            .iter()
            .filter(|p| p.len() == 2)
            .all(|p| set_has_base(&returned_coding, &p[0]))
        {
            t.selection_hits += 1;
        }
        // Predicted-parenthesization improvement: when the case diverges only
        // because ferro emits the spec-preferred parenthesized predicted
        // protein form (`p.(X)`) over the corpus's bare legacy `p.X` — with the
        // coding component matching — route to `improvement` (cluster
        // `predicted-parenthesization`) rather than the data-source classifier.
        // A clean match (or stale-annotation XPASS) still goes through the
        // normal path so the XPASS guard fires.
        let is_match = matches!(&actual, Ok(s) if s == &expected_repr);
        if !is_match && is_predicted_parenthesization_pairs(pairs, &returned_pairs) {
            t.improvement.push((
                case.input.clone(),
                CLUSTER_PREDICTED_PARENTHESIZATION.to_string(),
            ));
            continue;
        }

        // Classify on the coding (c.) component against ferro's returned coding
        // set — transcript selection and alignment-source skew manifest on the
        // coding coordinate; the protein is derived downstream.
        let expected_coding: Vec<String> = pairs
            .iter()
            .filter(|p| p.len() == 2)
            .map(|p| p[0].clone())
            .collect();
        t.record_classified(
            case,
            &expected_repr,
            actual,
            &expected_coding,
            &returned_coding,
        );
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

        // Capture ferro's returned coding set so we can tally selection
        // coverage (did ferro return the expected base transcript at all)
        // independently of whether the consequence matched.
        let mut returned: Vec<String> = Vec::new();
        let actual = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse: {e}"))?;
            let results = vp
                .project_variant_all(&v)
                .map_err(|e| format!("project_all: {e}"))?;
            let coding: Vec<String> = results
                .iter()
                .filter_map(|r| r.coding.as_ref().map(|c| c.to_string()))
                .collect();
            returned = coding.clone();
            if coding.iter().any(|c| consequence_matches(expected, c)) {
                Ok(expected.to_string())
            } else {
                Err(format!("missing coding {expected:?}; got {coding:?}"))
            }
        });

        t.selection_total += 1;
        if set_has_base(&returned, expected) {
            t.selection_hits += 1;
        }
        t.record_classified(case, expected, actual, &[expected.to_string()], &returned);
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

        let mut returned: Vec<String> = Vec::new();
        let actual = catch_panics(|| -> Result<String, String> {
            let v = parse_hgvs(&case.input).map_err(|e| format!("parse: {e}"))?;
            let results = vp
                .project_variant_all(&v)
                .map_err(|e| format!("project_all: {e}"))?;
            let noncoding: Vec<String> = results
                .iter()
                .filter_map(|r| r.noncoding.as_ref().map(|n| n.to_string()))
                .collect();
            returned = noncoding.clone();
            for want in expected_list {
                if !noncoding.iter().any(|n| consequence_matches(want, n)) {
                    return Err(format!("missing noncoding {want:?}; got {noncoding:?}"));
                }
            }
            Ok(expected.clone())
        });

        // Selection coverage: ferro returned the expected base transcript for
        // every expected noncoding entry (mirrors the all-must-match contract).
        t.selection_total += 1;
        if expected_list
            .iter()
            .all(|want| set_has_base(&returned, want))
        {
            t.selection_hits += 1;
        }
        t.record_classified(case, &expected, actual, expected_list, &returned);
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

    // `base_accession` strips both the version suffix and ferro's `(GENE)`
    // gene-symbol group, returning the bare base accession.
    #[test]
    fn base_accession_strips_version_and_gene_suffix() {
        assert_eq!(
            base_accession("NM_000682.7(ADRA2B):c.5A>T"),
            Some("NM_000682".to_string())
        );
        assert_eq!(
            base_accession("NM_000682.5:c.5A>T"),
            Some("NM_000682".to_string())
        );
    }

    // Version-only and gene-symbol-suffix differences are forgiven when the base
    // accession and consequence are identical.
    #[test]
    fn consequence_matches_is_version_insensitive_and_ignores_gene_suffix() {
        assert!(consequence_matches(
            "NM_000682.5:c.*1865G>T",
            "NM_000682.7(ADRA2B):c.*1865G>T"
        ));
        assert!(consequence_matches(
            "NM_000682.5:c.*1865G>T",
            "NM_000682.5:c.*1865G>T"
        ));
    }

    // The matcher must NOT forgive a coordinate/model-skew difference: same base
    // accession, different consequence coordinate → not a match.
    #[test]
    fn consequence_matches_rejects_coordinate_model_skew() {
        assert!(!consequence_matches(
            "NM_000682.5:c.*1865G>T",
            "NM_000682.7(ADRA2B):c.*1353+1856G>T"
        ));
    }

    // A different base accession is never a match, even with an identical
    // consequence.
    #[test]
    fn consequence_matches_rejects_different_base_accession() {
        assert!(!consequence_matches(
            "NR_111984.1:n.44G>A",
            "NM_001103170.2(AADACL3):n.51G>A"
        ));
        assert!(!consequence_matches(
            "NM_000001.1:c.5A>T",
            "NM_000002.1:c.5A>T"
        ));
    }

    // `set_has_base` is version-insensitive set membership on the base
    // accession: present iff some element shares the expected base accession,
    // regardless of version, gene suffix, or consequence.
    #[test]
    fn set_has_base_is_version_insensitive_membership() {
        let returned = vec![
            "NM_000682.7(ADRA2B):c.*1353+1856G>T".to_string(),
            "NM_000682.6(ADRA2B):c.*1353+1856G>T".to_string(),
        ];
        // Same base accession (different version + different consequence) → hit.
        assert!(set_has_base(&returned, "NM_000682.5:c.*1865G>T"));
        // Different base accession → miss, even though the set is non-empty.
        assert!(!set_has_base(&returned, "NR_111984.1:n.44G>A"));
        // Empty set → miss.
        assert!(!set_has_base(&[], "NM_000682.5:c.5A>T"));
    }

    // classify_divergence: a selection miss — the expected NR_ base accession
    // is absent from ferro's set (ferro returned a different transcript).
    #[test]
    fn classify_selection_miss() {
        assert!(matches!(
            classify_divergence(
                "NR_111984.1:n.44G>A",
                &["NM_001103170.2(AADACL3):n.51G>A".into()]
            ),
            Divergence::SelectionMiss
        ));
    }

    // classify_divergence: coordinate skew — same base accession, different
    // position prefix (`*1865` vs `*1353+1856`).
    #[test]
    fn classify_coordinate_skew() {
        assert!(matches!(
            classify_divergence(
                "NM_000682.5:c.*1865G>T",
                &["NM_000682.7(ADRA2B):c.*1353+1856G>T".into()]
            ),
            Divergence::CoordinateSkew
        ));
    }

    // classify_divergence: form-only — same base, same position prefix, the
    // edit form differs (`insG` over a range vs `dup`). NOT data-source.
    #[test]
    fn classify_form_only() {
        assert!(matches!(
            classify_divergence(
                "NM_X.1:c.10+830_10+831insG",
                &["NM_X.2(G):c.10+830_10+831dup".into()]
            ),
            Divergence::FormOnly
        ));
    }

    // classify_divergence: a match short-circuits to `Match`.
    #[test]
    fn classify_match() {
        assert!(matches!(
            classify_divergence("NM_X.1:c.5A>T", &["NM_X.2(G):c.5A>T".into()]),
            Divergence::Match
        ));
    }

    // position_prefix extracts the leading coordinate token before the edit,
    // handling `*` UTR, `+`/`-` intronic offsets, and `_` ranges.
    #[test]
    fn position_prefix_handles_utr_offset_and_range() {
        assert_eq!(position_prefix("c.*1865G>T"), "*1865");
        assert_eq!(position_prefix("c.*1353+1856G>T"), "*1353+1856");
        assert_eq!(position_prefix("c.10+830_10+831insG"), "10+830_10+831");
        assert_eq!(position_prefix("n.51G>A"), "51");
        assert_eq!(position_prefix("c.-12del"), "-12");
        assert_eq!(position_prefix("c.5_7dup"), "5_7");
    }

    // is_predicted_parenthesization: true iff `actual` is `expected` with the
    // `p.`-token wrapped in `(...)` — same base accession, same inner token.
    #[test]
    fn predicted_parenthesization_wraps_bare_protein() {
        assert!(is_predicted_parenthesization(
            "NP_000673.2:p.Arg268Trp",
            "NP_000673.2:p.(Arg268Trp)"
        ));
    }

    // Version digit and `(GENE)` suffix differences on the base accession are
    // forgiven (the corpus is version-insensitive), as long as the inner
    // protein token is the bare → wrapped pair.
    #[test]
    fn predicted_parenthesization_is_version_and_gene_insensitive() {
        assert!(is_predicted_parenthesization(
            "NP_000673.2:p.Lys233Ile",
            "NP_000673.2(ADRA2B):p.(Lys233Ile)"
        ));
    }

    // A different inner protein token is NOT a parenthesization-only divergence.
    #[test]
    fn predicted_parenthesization_rejects_inner_token_difference() {
        assert!(!is_predicted_parenthesization(
            "NP_000673.2:p.Arg268Trp",
            "NP_000673.2:p.(Arg268Cys)"
        ));
    }

    // Already-equal (both wrapped) adds no wrap → not a parenthesization
    // divergence (the harness handles equality as a match/XPASS elsewhere).
    #[test]
    fn predicted_parenthesization_rejects_already_equal() {
        assert!(!is_predicted_parenthesization(
            "NP_000673.2:p.(Arg268Trp)",
            "NP_000673.2:p.(Arg268Trp)"
        ));
    }

    // A bare→bare pair (no wrap added) is not a parenthesization divergence.
    #[test]
    fn predicted_parenthesization_rejects_no_wrap() {
        assert!(!is_predicted_parenthesization(
            "NP_000673.2:p.Arg268Trp",
            "NP_000673.2:p.Arg268Trp"
        ));
    }

    // A different base accession is never a parenthesization-only divergence,
    // even if the inner token would wrap.
    #[test]
    fn predicted_parenthesization_rejects_different_base() {
        assert!(!is_predicted_parenthesization(
            "NP_000673.2:p.Arg268Trp",
            "NP_000042.3:p.(Arg268Trp)"
        ));
    }

    // is_predicted_parenthesization_pairs: the pairs-axis wrapper that drives
    // the `improvement` routing in `axis_coding_protein_descriptions`. It holds
    // iff every applicable (c., p.) pair has a returned pair whose coding
    // matches version-insensitively AND whose protein is the bare→wrapped
    // predicted form. A single clean matching pair → true.
    #[test]
    fn predicted_parenthesization_pairs_holds_for_matching_pair() {
        let expected = vec![vec![
            "NM_000673.2:c.803G>A".to_string(),
            "NP_000673.2:p.Arg268Trp".to_string(),
        ]];
        // Coding matches version- and gene-insensitively; protein is bare→wrapped.
        let returned = vec![(
            "NM_000673.5(GENE):c.803G>A".to_string(),
            "NP_000673.2:p.(Arg268Trp)".to_string(),
        )];
        assert!(is_predicted_parenthesization_pairs(&expected, &returned));
    }

    // A coding-component mismatch (coordinate skew on the c. position) is not a
    // parenthesization-only divergence, even when the protein wraps cleanly.
    #[test]
    fn predicted_parenthesization_pairs_rejects_coding_mismatch() {
        let expected = vec![vec![
            "NM_000673.2:c.803G>A".to_string(),
            "NP_000673.2:p.Arg268Trp".to_string(),
        ]];
        let returned = vec![(
            "NM_000673.5(GENE):c.999G>A".to_string(),
            "NP_000673.2:p.(Arg268Trp)".to_string(),
        )];
        assert!(!is_predicted_parenthesization_pairs(&expected, &returned));
    }

    // An inner protein-token difference (Trp vs Cys) is a genuine delta, not a
    // parenthesization-only divergence — false even though the coding matches.
    #[test]
    fn predicted_parenthesization_pairs_rejects_inner_token_difference() {
        let expected = vec![vec![
            "NM_000673.2:c.803G>A".to_string(),
            "NP_000673.2:p.Arg268Trp".to_string(),
        ]];
        let returned = vec![(
            "NM_000673.5(GENE):c.803G>A".to_string(),
            "NP_000673.2:p.(Arg268Cys)".to_string(),
        )];
        assert!(!is_predicted_parenthesization_pairs(&expected, &returned));
    }

    // No applicable (length-2) expected pair → false: the `applicable.is_empty()`
    // guard must not route a pair-less or malformed case to `improvement`.
    #[test]
    fn predicted_parenthesization_pairs_rejects_empty_applicable() {
        let returned = vec![(
            "NM_000673.5(GENE):c.803G>A".to_string(),
            "NP_000673.2:p.(Arg268Trp)".to_string(),
        )];
        // A malformed singleton (len != 2) leaves `applicable` empty → false.
        let malformed = vec![vec!["NM_000673.2:c.803G>A".to_string()]];
        assert!(!is_predicted_parenthesization_pairs(&malformed, &returned));
        // No expected pairs at all → also false.
        assert!(!is_predicted_parenthesization_pairs(&[], &returned));
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

// ----------------------------------------------------------------------------
// Hermetic #615 stop-loss extension regression guards (no manifest, no corpus)
// ----------------------------------------------------------------------------
//
// The 4 ex-`stoploss-xaa` rows in cases.json (now folded into the
// `extter-vs-legacy-glyph` improvement cluster) exercise the same stop-loss path
// against the corpus manifest. This module pins the behavior hermetically: a
// hand-built `MockProvider` transcript (CDS ending in a stop codon + a 3'UTR)
// plus a single-base indel that disrupts the terminator codon, projected through
// the real `VariantProjector`. #615 is fixed (PR #621 routes stop-disrupting
// indels to a C-terminal extension), so these assert ferro now emits the spec
// `p.(Ter…ext…)` extension — not the old `p.(Xaa…)` placeholder — independent of
// the manifest.
#[cfg(test)]
mod stoploss_615_reproducer {
    use super::*;
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand};

    /// Build a single-exon plus-strand coding transcript from a literal mRNA
    /// sequence, CDS bounds (1-based, inclusive), and an explicit protein
    /// accession so the direct `c.→p.` path can name the predicted consequence.
    fn stoploss_transcript() -> Transcript {
        // mRNA:  ATG CGC AAA TAA  GCATAAGGGCCC
        //        Met Arg Lys Ter  └─── 3'UTR ───┘
        // CDS = c.1..12 (ATGCGCAAATAA → Met-Arg-Lys-Ter). The terminator codon
        // `TAA` sits at c.10..12; a deletion inside it is a stop-loss.
        let seq = "ATGCGCAAATAAGCATAAGGGCCC";
        let len = seq.len() as u64;
        Transcript::new(
            "NM_700000.1".to_string(),
            Some("STOPLOSS".to_string()),
            Strand::Plus,
            seq.to_string(),
            Some(1),
            Some(12),
            vec![Exon::new(1, 1, len)],
            None,
            None,
            None,
            Default::default(),
            ManeStatus::Select,
            None,
            None,
        )
        .with_protein_id("NP_700000.1".to_string())
    }

    /// A Mock-backed projector seeded with the stop-loss transcript. The cdot
    /// mapper is empty; the direct bare-`NM_` `c.→p.` path reads the CDS and
    /// protein straight from the transcript record, so no genome alignment is
    /// needed.
    fn stoploss_projector() -> VariantProjector<MockProvider> {
        let mut provider = MockProvider::new();
        provider.add_transcript(stoploss_transcript());
        let projector = Projector::new(CdotMapper::new());
        VariantProjector::new(projector, provider)
    }

    fn project_protein(input: &str) -> String {
        let vp = stoploss_projector();
        let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse `{input}`: {e}"));
        let result = vp
            .project_variant(&v, "NM_700000.1")
            .unwrap_or_else(|e| panic!("project `{input}`: {e}"));
        result
            .protein
            .as_ref()
            .map(|p| p.to_string())
            .unwrap_or_else(|| panic!("no protein predicted for `{input}`"))
    }

    // #615 FIXED (PR #621): a deletion that disrupts the terminator codon is now
    // rendered as a C-terminal stop-loss extension (`p.(Ter…ext…)`), not the old
    // `p.(Xaa…)` placeholder. Regression guard for the fix.
    #[test]
    fn stoploss_deletion_in_terminator_emits_extension() {
        // c.11delA deletes the middle base of the `TAA` terminator (c.10..12),
        // a stop-loss → spec C-terminal extension.
        let got = project_protein("NM_700000.1:c.11delA");
        assert_eq!(
            got, "NP_700000.1:p.(Ter4extTer0)",
            "#615 fixed: stop-loss deletion should yield a Ter…ext extension; got {got}"
        );
        assert!(
            got.contains("ext") && !got.contains("Xaa"),
            "#615 fixed: expected a spec extension, not an Xaa placeholder; got {got}"
        );
    }

    // #615 FIXED (PR #621, duplication variant): a single-base duplication that
    // frameshifts through the terminator is now rendered as the spec C-terminal
    // extension, not an `Xaa` placeholder. Mirrors the corpus `dup` rows (e.g.
    // NM_000425.3:c.3772dupT, now an extter-vs-legacy-glyph improvement).
    #[test]
    fn stoploss_duplication_through_terminator_emits_extension() {
        let got = project_protein("NM_700000.1:c.9dupA");
        assert_eq!(
            got, "NP_700000.1:p.(Ter4IleextTer?)",
            "#615 fixed: stop-disrupting dup should yield a Ter…ext extension; got {got}"
        );
        assert!(
            got.contains("ext") && !got.contains("Xaa"),
            "#615 fixed: expected a spec extension, not an Xaa placeholder; got {got}"
        );
    }

    // Control: a duplication that cleanly converts the terminator to a sense
    // codon DOES produce the spec extension form (`extTer…`), proving the
    // predictor is capable of the correct output and #615 is a specific
    // disruption-path bug, not a wholesale failure. This must stay GREEN both
    // before and after the #615 fix.
    #[test]
    fn stoploss_clean_readthrough_emits_extension_not_xaa() {
        let got = project_protein("NM_700000.1:c.10dupT");
        assert!(
            got.contains("ext") && !got.contains("Xaa"),
            "a clean stop-loss readthrough should yield an extension, not Xaa; got {got}"
        );
    }
}
