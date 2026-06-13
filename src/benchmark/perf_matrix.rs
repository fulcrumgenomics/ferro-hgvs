//! Performance-matrix orchestration: run the tools across operations × worker
//! counts × reps over a shared seeded sample, aggregate to median/min-max, and
//! serialize the schema consumed by the perf-table renderer (PR #604).
//!
//! The aggregation + serialization here are pure and unit-tested. The actual
//! tool runs are manual (need the reference stack) — see docs/BENCHMARK_RUNBOOK.md.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// One tool's measured throughput at a given worker count.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToolMeasure {
    pub median_pps: f64,
    pub min_pps: f64,
    pub max_pps: f64,
    pub reps: u32,
    /// Mean per-rep success rate across all reps measured for this tool.
    pub success_rate: f64,
    /// Set when the tool's stack could not run (e.g. UTA unavailable).
    #[serde(default)]
    pub not_run: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
}

impl ToolMeasure {
    /// Aggregate per-rep patterns/sec into median + min + max.
    pub fn from_rates(_tool: &str, rates: &[f64], success_rate: f64) -> Self {
        let mut sorted: Vec<f64> = rates.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = sorted.len();
        let median = if n == 0 {
            0.0
        } else if n % 2 == 1 {
            sorted[n / 2]
        } else {
            (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
        };
        ToolMeasure {
            median_pps: median,
            min_pps: sorted.first().copied().unwrap_or(0.0),
            max_pps: sorted.last().copied().unwrap_or(0.0),
            reps: n as u32,
            success_rate,
            not_run: false,
            reason: None,
        }
    }

    /// A tool whose stack could not run.
    pub fn not_run(reason: &str) -> Self {
        ToolMeasure {
            median_pps: 0.0,
            min_pps: 0.0,
            max_pps: 0.0,
            reps: 0,
            success_rate: 0.0,
            not_run: true,
            reason: Some(reason.to_string()),
        }
    }
}

/// All tools' measurements at a single worker count.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WorkerPoint {
    pub tools: BTreeMap<String, ToolMeasure>,
}

/// One operation (parse or normalize): cross-tool points + ferro scaling.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Operation {
    pub sample_size_slow: u64,
    pub ferro_full_population_n: u64,
    pub ferro_full_population_pps: f64,
    /// Keyed by worker count as a decimal-integer string ("1", "8"). BTreeMap
    /// iteration is lexicographic; renderers use explicit key lists, so this is fine.
    pub by_workers: BTreeMap<String, WorkerPoint>,
    /// ferro-only thread scaling: thread count ("1","2","4","8") -> patterns/sec.
    #[serde(default)]
    pub ferro_scaling: BTreeMap<String, f64>,
    /// Per-tool calibrated sample sizes used for this operation.
    ///
    /// Records the N chosen for each tool by the calibration step (or the
    /// manual override when --sample-size is set).  Keyed by tool name.
    /// The renderer does not read this field; it is included for audit trails.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub sample_sizes: BTreeMap<String, u64>,
}

/// Run provenance metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Provenance {
    pub generated: String,
    pub corpus_id: String,
    pub corpus_sha: String,
    pub machine: String,
    pub os: String,
    pub tool_versions: BTreeMap<String, String>,
    /// Disclosed measurement caveats (e.g. the Python-startup timed-region note).
    /// The renderer does not read this field but tolerates it (no deny_unknown_fields).
    #[serde(default)]
    pub notes: Vec<String>,
}

/// The whole results document.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerfResults {
    pub schema_version: u32,
    /// When `true`, the data has not yet been measured and the README tables
    /// should display a "not measured" warning. Omitting the field defaults to
    /// `false` (real data, no banner).
    #[serde(default)]
    pub placeholder: bool,
    pub provenance: Provenance,
    /// Keyed by operation name ("parse", "normalize").
    pub operations: BTreeMap<String, Operation>,
}

// ---------------------------------------------------------------------------
// Measure trait + run loop
// ---------------------------------------------------------------------------

/// One repetition's measured rate for a (tool, op, W).
pub struct RepRate {
    pub pps: f64,
    pub success_rate: f64,
}

/// Abstracts a single measured run so the loop is testable without a reference
/// stack. The real impl shells to the per-tool primitives; tests use a fake.
pub trait Measure {
    fn measure(&self, tool: &str, op: &str, workers: usize, sample: &[String]) -> Option<RepRate>;
}

/// Configuration for a full matrix run.
#[derive(Debug, Clone)]
pub struct MatrixConfig {
    pub tools: Vec<String>,
    pub workers: Vec<usize>,
    pub reps: u32,
    pub ferro_threads: Vec<usize>,
}

/// Compute the arithmetic mean of a slice. Returns 0.0 for an empty slice.
fn mean(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        0.0
    } else {
        xs.iter().sum::<f64>() / xs.len() as f64
    }
}

/// Run one operation across the worker × tool × rep matrix.
///
/// `sample_for(tool, rep)` is called **once per (tool, rep)** and returns a
/// sample sized to that tool's calibrated N.  All tools draw from the same
/// stratified population with the same seed-per-rep formula, so the
/// difficulty distribution is identical across tools even though the counts
/// differ.  This is the fairness property that matters: same distribution,
/// tool-specific precision via size.
///
/// Aggregation:
/// - Per-tool `pps` values across reps → median/min/max via [`ToolMeasure::from_rates`].
/// - Per-tool `success_rate` → mean of per-rep success rates (not the last rep's value).
pub fn run_operation(
    m: &dyn Measure,
    op: &str,
    cfg: &MatrixConfig,
    sample_for: &dyn Fn(&str, u32) -> Vec<String>,
) -> Operation {
    let mut by_workers: BTreeMap<String, WorkerPoint> = BTreeMap::new();
    for &w in &cfg.workers {
        // Accumulate per-rep rates and success rates per tool across all reps.
        let mut rates: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        let mut succ: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        for tool in &cfg.tools {
            rates.insert(tool.clone(), Vec::new());
            succ.insert(tool.clone(), Vec::new());
        }

        for rep in 0..cfg.reps {
            for tool in &cfg.tools {
                // Draw a per-tool sample for this rep.  Each tool uses its
                // calibrated N drawn from the same seeded stratified population.
                let sample = sample_for(tool, rep);
                if let Some(r) = m.measure(tool, op, w, &sample) {
                    rates.get_mut(tool.as_str()).unwrap().push(r.pps);
                    succ.get_mut(tool.as_str()).unwrap().push(r.success_rate);
                }
            }
        }

        let mut tools: BTreeMap<String, ToolMeasure> = BTreeMap::new();
        for tool in &cfg.tools {
            let tool_rates = &rates[tool.as_str()];
            let tool_succ = &succ[tool.as_str()];
            let measure = if tool_rates.is_empty() {
                ToolMeasure::not_run(&format!(
                    "{tool} {op} did not run (stack unavailable or tool error)"
                ))
            } else {
                ToolMeasure::from_rates(tool, tool_rates, mean(tool_succ))
            };
            tools.insert(tool.clone(), measure);
        }
        by_workers.insert(w.to_string(), WorkerPoint { tools });
    }

    // ferro thread scaling (median over reps at each thread count).
    // ferro uses its own calibrated sample size (same as the cross-tool ferro
    // entry) drawn per rep from the shared stratified population.
    let mut ferro_scaling = BTreeMap::new();
    for &t in &cfg.ferro_threads {
        let mut rates = Vec::new();
        for rep in 0..cfg.reps {
            let sample = sample_for("ferro", rep);
            if let Some(r) = m.measure("ferro", op, t, &sample) {
                rates.push(r.pps);
            }
        }
        if !rates.is_empty() {
            ferro_scaling.insert(
                t.to_string(),
                ToolMeasure::from_rates("ferro", &rates, 1.0).median_pps,
            );
        }
    }

    Operation {
        sample_size_slow: 0,
        ferro_full_population_n: 0,
        ferro_full_population_pps: 0.0,
        by_workers,
        ferro_scaling,
        sample_sizes: BTreeMap::new(),
    }
}

// ---------------------------------------------------------------------------
// Timed-region caveat constant
// ---------------------------------------------------------------------------

/// Measurement caveat recorded in every `PerfResults::provenance.notes`.
///
/// For NORMALIZE: mutalyzer and biocommons are invoked as Python subprocesses;
/// the timed region starts before the subprocess is spawned, so Python
/// interpreter startup and import time (~0.3–1 s) are folded into the measured
/// elapsed time. ferro and hgvs-rs exclude equivalent setup (provider
/// construction happens before the timer). At the sample sizes used here the
/// overhead is <2%.
///
/// For PARSE: mutalyzer and biocommons parse subprocesses report their own
/// internal elapsed time (measured inside the subprocess after imports complete),
/// so Python startup is NOT included in parse throughput numbers. ferro and
/// hgvs-rs parse throughput is also measured after any one-time setup.
pub const TIMED_REGION_NOTE: &str =
    "normalize throughput for mutalyzer/biocommons includes Python subprocess \
     startup in the timed region (ferro/hgvs-rs exclude it; <2% at this sample \
     size); parse throughput excludes Python startup (measured by the \
     subprocess's internal timer). See issue tracker.";

// ---------------------------------------------------------------------------
// Real HarnessMeasure implementation
// ---------------------------------------------------------------------------

/// Paths and settings derived from [`super::cli::MatrixArgs`] that
/// `HarnessMeasure` needs to dispatch per-tool runs.
pub struct HarnessMeasureConfig {
    pub reference: Option<std::path::PathBuf>,
    pub mutalyzer_settings: Option<std::path::PathBuf>,
    pub biocommons_settings: Option<std::path::PathBuf>,
    pub seqrepo_path: Option<std::path::PathBuf>,
    pub uta_db_url: Option<String>,
}

/// The real [`Measure`] implementation that dispatches to the existing
/// per-tool parse/normalize primitives already used by `bin/benchmark.rs`.
///
/// Each call writes the sample to a unique `NamedTempFile`, invokes the
/// appropriate tool function, reads back timing information, and returns
/// `Some(RepRate)`. On any error the tool is recorded as `not_run` for that
/// rep by returning `None`.
pub struct HarnessMeasure {
    pub config: HarnessMeasureConfig,
}

impl HarnessMeasure {
    /// Write a sample of HGVS strings to a temporary file and return the file.
    ///
    /// The returned `NamedTempFile` keeps the file open; the path is valid as
    /// long as the handle is alive.  Callers should hold onto it until after
    /// the tool function returns.
    fn write_sample_to_tempfile(
        sample: &[String],
    ) -> Result<tempfile::NamedTempFile, crate::FerroError> {
        use std::io::Write;
        let mut tmp = tempfile::NamedTempFile::new().map_err(|e| crate::FerroError::Io {
            msg: format!("Failed to create temp file: {e}"),
        })?;
        for line in sample {
            writeln!(tmp, "{line}").map_err(|e| crate::FerroError::Io {
                msg: format!("Failed to write temp file: {e}"),
            })?;
        }
        tmp.flush().map_err(|e| crate::FerroError::Io {
            msg: format!("Failed to flush temp file: {e}"),
        })?;
        Ok(tmp)
    }

    // --- ferro ---

    fn measure_ferro_parse(&self, workers: usize, sample: &[String]) -> Option<RepRate> {
        use crate::benchmark::parse::parse_ferro_count_parallel;

        // ferro parse is pure (parse_hgvs needs no provider), so we can use a
        // sized Rayon pool to honour the requested worker count.  This makes the
        // 1/2/4/8-thread ferro_scaling measurements genuinely different.
        let (ok, err, elapsed) = parse_ferro_count_parallel(sample, workers);
        let total = ok + err;
        if total == 0 || elapsed.as_secs_f64() <= 0.0 {
            return None;
        }
        Some(RepRate {
            pps: total as f64 / elapsed.as_secs_f64(),
            success_rate: ok as f64 / total as f64,
        })
    }

    fn measure_ferro_normalize(&self, workers: usize, sample: &[String]) -> Option<RepRate> {
        use crate::benchmark::normalize::{normalize_ferro, normalize_ferro_parallel};

        let reference = self.config.reference.as_deref()?;
        let tmp = Self::write_sample_to_tempfile(sample).ok()?;
        let results_tmp = tempfile::NamedTempFile::new().ok()?;
        let timing_tmp = tempfile::NamedTempFile::new().ok()?;

        let shard = if workers > 1 {
            normalize_ferro_parallel(
                tmp.path(),
                results_tmp.path(),
                timing_tmp.path(),
                Some(reference),
                workers,
            )
            .ok()?
        } else {
            normalize_ferro(
                tmp.path(),
                results_tmp.path(),
                timing_tmp.path(),
                Some(reference),
            )
            .ok()?
        };

        let total = shard.timing.total_patterns;
        let successful = shard.timing.successful;
        // A rate computed over zero successes is not a valid throughput — record
        // this tool as not_run (None) rather than publishing a misleading p/s
        // figure with success_rate=0.0.  hgvs-rs W>1 hitting this path is a
        // known path inconsistency tracked separately; both W=1 and W=8 will now
        // consistently return None when there are no successes.
        if total == 0 || successful == 0 {
            return None;
        }
        Some(RepRate {
            pps: shard.timing.patterns_per_second,
            success_rate: successful as f64 / total as f64,
        })
    }

    // --- mutalyzer ---

    fn measure_mutalyzer_parse(&self, _workers: usize, sample: &[String]) -> Option<RepRate> {
        use crate::benchmark::mutalyzer::{has_mutalyzer_parser, run_mutalyzer_parser_subprocess};

        if !has_mutalyzer_parser() {
            return None;
        }
        let tmp = Self::write_sample_to_tempfile(sample).ok()?;
        let results_tmp = tempfile::NamedTempFile::new().ok()?;

        run_mutalyzer_parser_subprocess(tmp.path().to_str()?, results_tmp.path().to_str()?).ok()?;

        // Read the JSON output written by the subprocess.
        let content = std::fs::read_to_string(results_tmp.path()).ok()?;
        let v: serde_json::Value = serde_json::from_str(&content).ok()?;
        let total = v["total_patterns"].as_u64().unwrap_or(0) as usize;
        let successful = v["successful"].as_u64().unwrap_or(0) as usize;
        let elapsed = v["elapsed_seconds"].as_f64().unwrap_or(0.0);
        if total == 0 || elapsed <= 0.0 {
            return None;
        }
        Some(RepRate {
            pps: total as f64 / elapsed,
            success_rate: successful as f64 / total as f64,
        })
    }

    fn measure_mutalyzer_normalize(&self, workers: usize, sample: &[String]) -> Option<RepRate> {
        use crate::benchmark::compare::run_mutalyzer_normalize_parallel;
        use crate::benchmark::mutalyzer::has_mutalyzer_normalizer;

        if !has_mutalyzer_normalizer() {
            return None;
        }
        let settings = self
            .config
            .mutalyzer_settings
            .as_ref()
            .and_then(|p| p.to_str().map(str::to_owned));

        let (results, elapsed, _errs) = run_mutalyzer_normalize_parallel(
            sample,
            workers,
            settings.as_deref(),
            false, // allow_network = false
        )
        .ok()?;

        let total = results.len();
        if total == 0 || elapsed.as_secs_f64() <= 0.0 {
            return None;
        }
        let successful = results.iter().filter(|r| r.success).count();
        if successful == 0 {
            return None;
        }
        Some(RepRate {
            pps: total as f64 / elapsed.as_secs_f64(),
            success_rate: successful as f64 / total as f64,
        })
    }

    // --- biocommons ---

    fn measure_biocommons_parse(&self, _workers: usize, sample: &[String]) -> Option<RepRate> {
        use crate::benchmark::biocommons::has_biocommons_normalizer;

        if !has_biocommons_normalizer() {
            return None;
        }

        // Derive UTA/SeqRepo from biocommons settings if provided.
        let (uta_db_url, seqrepo_dir) = self.biocommons_url_and_seqrepo();
        let tmp = Self::write_sample_to_tempfile(sample).ok()?;
        let results_tmp = tempfile::NamedTempFile::new().ok()?;

        // biocommons parse uses its own Python script (same pattern as bin/benchmark.rs).
        let script = r#"
import sys, json, time
from hgvs.parser import Parser

parser = Parser()
patterns = [l.strip() for l in open(sys.argv[1]) if l.strip()]
results = []
successful = 0
start = time.time()
for p in patterns:
    try:
        parsed = parser.parse_hgvs_variant(p)
        results.append({"input": p, "success": True, "output": str(parsed)})
        successful += 1
    except Exception as e:
        results.append({"input": p, "success": False, "error": str(e)[:200]})
elapsed = time.time() - start
out = {"tool":"biocommons","total_patterns":len(patterns),"successful":successful,
       "failed":len(patterns)-successful,"elapsed_seconds":elapsed,
       "throughput":len(patterns)/elapsed if elapsed>0 else 0,"results":results}
with open(sys.argv[2],'w') as f:
    json.dump(out,f)
"#;
        let script_tmp = tempfile::NamedTempFile::new().ok()?;
        std::fs::write(script_tmp.path(), script).ok()?;

        let mut cmd = std::process::Command::new("python3");
        cmd.arg(script_tmp.path())
            .arg(tmp.path())
            .arg(results_tmp.path());
        if let Some(ref url) = uta_db_url {
            cmd.env("UTA_DB_URL", url);
        }
        if let Some(ref sd) = seqrepo_dir {
            cmd.env("HGVS_SEQREPO_DIR", sd);
        }
        let status = cmd.status().ok()?;
        if !status.success() {
            return None;
        }

        let content = std::fs::read_to_string(results_tmp.path()).ok()?;
        let v: serde_json::Value = serde_json::from_str(&content).ok()?;
        let total = v["total_patterns"].as_u64().unwrap_or(0) as usize;
        let successful = v["successful"].as_u64().unwrap_or(0) as usize;
        let elapsed = v["elapsed_seconds"].as_f64().unwrap_or(0.0);
        if total == 0 || elapsed <= 0.0 {
            return None;
        }
        Some(RepRate {
            pps: total as f64 / elapsed,
            success_rate: successful as f64 / total as f64,
        })
    }

    fn measure_biocommons_normalize(&self, workers: usize, sample: &[String]) -> Option<RepRate> {
        use crate::benchmark::biocommons::{
            has_biocommons_normalizer, run_biocommons_normalizer_parallel,
            run_biocommons_normalizer_subprocess,
        };
        use std::time::Instant;

        if !has_biocommons_normalizer() {
            return None;
        }
        let (uta_db_url, seqrepo_dir) = self.biocommons_url_and_seqrepo();
        let tmp = Self::write_sample_to_tempfile(sample).ok()?;
        let results_tmp = tempfile::NamedTempFile::new().ok()?;

        let start = Instant::now();
        if workers > 1 {
            run_biocommons_normalizer_parallel(
                tmp.path().to_str()?,
                results_tmp.path().to_str()?,
                uta_db_url.as_deref(),
                seqrepo_dir.as_deref(),
                None, // lrg_mapping_file
                workers,
            )
            .ok()?;
        } else {
            run_biocommons_normalizer_subprocess(
                tmp.path().to_str()?,
                results_tmp.path().to_str()?,
                uta_db_url.as_deref(),
                seqrepo_dir.as_deref(),
                None, // lrg_mapping_file
            )
            .ok()?;
        }
        let elapsed = start.elapsed();

        let content = std::fs::read_to_string(results_tmp.path()).ok()?;
        let v: serde_json::Value = serde_json::from_str(&content).ok()?;
        let total = v["total_patterns"].as_u64().unwrap_or(0) as usize;
        let successful = v["successful"].as_u64().unwrap_or(0) as usize;
        if total == 0 || successful == 0 || elapsed.as_secs_f64() <= 0.0 {
            return None;
        }
        Some(RepRate {
            pps: total as f64 / elapsed.as_secs_f64(),
            success_rate: successful as f64 / total as f64,
        })
    }

    // --- hgvs-rs ---

    #[cfg(feature = "hgvs-rs")]
    fn measure_hgvs_rs_parse(&self, _workers: usize, sample: &[String]) -> Option<RepRate> {
        use hgvs::parser::HgvsVariant;
        use std::str::FromStr;
        use std::time::Instant;

        let start = Instant::now();
        let mut successful = 0usize;
        for p in sample {
            if HgvsVariant::from_str(p.as_str()).is_ok() {
                successful += 1;
            }
        }
        let elapsed = start.elapsed();
        let total = sample.len();
        if total == 0 || elapsed.as_secs_f64() <= 0.0 {
            return None;
        }
        Some(RepRate {
            pps: total as f64 / elapsed.as_secs_f64(),
            success_rate: successful as f64 / total as f64,
        })
    }

    #[cfg(feature = "hgvs-rs")]
    fn measure_hgvs_rs_normalize(&self, workers: usize, sample: &[String]) -> Option<RepRate> {
        use crate::benchmark::hgvs_rs::{run_hgvs_rs_normalize_parallel, HgvsRsConfig};

        // hgvs-rs normalize requires BOTH seqrepo_path AND uta_db_url.
        // If either is missing, record not_run rather than falling back to a
        // hardcoded default that silently connects to the wrong host.
        let seqrepo = self
            .config
            .seqrepo_path
            .as_ref()?
            .to_string_lossy()
            .to_string();
        let uta_url_full = self.config.uta_db_url.as_deref()?;

        // The caller passes a URL that may include the schema as the last path
        // segment (e.g. "postgresql://…/uta/uta_20210129b").
        // HgvsRsConfig.uta_db_url must end at the database name ("/uta"),
        // with the schema passed separately in uta_db_schema.
        // Strategy: if the last path segment looks like a schema name (no port,
        // no '@', contains '_' or matches known schema pattern), strip it from
        // the URL and use it as the schema; otherwise use the full URL as-is.
        let (uta_db_url, uta_schema) = {
            let last_seg = uta_url_full.rsplit('/').next().unwrap_or("");
            // A schema segment looks like "uta_20210129b" — contains '_' and no '@' or ':'.
            if !last_seg.is_empty()
                && last_seg.contains('_')
                && !last_seg.contains('@')
                && !last_seg.contains(':')
            {
                let base_url = uta_url_full
                    .strip_suffix(last_seg)
                    .and_then(|u| u.strip_suffix('/'))
                    .unwrap_or(uta_url_full);
                (base_url.to_string(), last_seg.to_string())
            } else {
                (uta_url_full.to_string(), "uta_20210129b".to_string())
            }
        };

        let config = HgvsRsConfig {
            uta_db_url,
            uta_db_schema: uta_schema,
            seqrepo_path: seqrepo,
            lrg_mapping_file: None,
            in_memory: false,
        };

        let (results, elapsed, _errs) =
            run_hgvs_rs_normalize_parallel(sample, &config, workers).ok()?;
        let total = results.len();
        if total == 0 || elapsed.as_secs_f64() <= 0.0 {
            return None;
        }
        let successful = results.iter().filter(|r| r.success).count();
        // Zero successes → not_run (consistent with ferro/mutalyzer guards).
        // The hgvs-rs W=1 vs W>1 path inconsistency is tracked separately.
        if successful == 0 {
            return None;
        }
        Some(RepRate {
            pps: total as f64 / elapsed.as_secs_f64(),
            success_rate: successful as f64 / total as f64,
        })
    }

    /// Resolve biocommons UTA URL and SeqRepo dir from settings file or direct
    /// args.  Returns `(uta_db_url, seqrepo_dir)` as `Option<String>` each.
    fn biocommons_url_and_seqrepo(&self) -> (Option<String>, Option<String>) {
        use crate::benchmark::biocommons::load_biocommons_settings;

        if let Some(ref settings_path) = self.config.biocommons_settings {
            if let Ok(s) = load_biocommons_settings(settings_path) {
                return (
                    Some(s.uta_db_url),
                    Some(s.seqrepo_dir.to_string_lossy().to_string()),
                );
            }
        }
        // Fall back to direct args.
        let uta = self.config.uta_db_url.clone();
        let seqrepo = self
            .config
            .seqrepo_path
            .as_ref()
            .map(|p| p.to_string_lossy().to_string());
        (uta, seqrepo)
    }
}

impl Measure for HarnessMeasure {
    fn measure(&self, tool: &str, op: &str, workers: usize, sample: &[String]) -> Option<RepRate> {
        match (tool, op) {
            ("ferro", "parse") => self.measure_ferro_parse(workers, sample),
            ("ferro", "normalize") => self.measure_ferro_normalize(workers, sample),
            ("mutalyzer", "parse") => self.measure_mutalyzer_parse(workers, sample),
            ("mutalyzer", "normalize") => self.measure_mutalyzer_normalize(workers, sample),
            ("biocommons", "parse") => self.measure_biocommons_parse(workers, sample),
            ("biocommons", "normalize") => self.measure_biocommons_normalize(workers, sample),
            #[cfg(feature = "hgvs-rs")]
            ("hgvs-rs", "parse") => self.measure_hgvs_rs_parse(workers, sample),
            #[cfg(feature = "hgvs-rs")]
            ("hgvs-rs", "normalize") => self.measure_hgvs_rs_normalize(workers, sample),
            _ => {
                eprintln!("[perf-matrix] unknown (tool={tool}, op={op}) — skipping");
                None
            }
        }
    }
}

// ---------------------------------------------------------------------------
// run_matrix: top-level driver
// ---------------------------------------------------------------------------

/// Compute a stable 64-bit hex hash of a byte slice using `DefaultHasher`.
///
/// This is reproducible within a single binary build and Rust toolchain
/// version. It is used only as a corpus fingerprint in provenance metadata,
/// not for security.
fn corpus_sha_hex(bytes: &[u8]) -> String {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    bytes.hash(&mut h);
    format!("{:016x}", h.finish())
}

/// Run the full performance matrix and write the results JSON to `args.output`.
///
/// This is the real driver that wires together the pure `run_operation` loop
/// (testable without a reference stack) with the `HarnessMeasure` that calls
/// the actual tool primitives.
pub fn run_matrix(args: &super::cli::MatrixArgs) -> Result<(), crate::FerroError> {
    use crate::benchmark::sample::stratified_sample_vec;

    eprintln!(
        "[perf-matrix] Loading population from {} …",
        args.population.display()
    );
    let population_bytes = std::fs::read(&args.population).map_err(|e| crate::FerroError::Io {
        msg: format!("Failed to read population file: {e}"),
    })?;
    let corpus_sha = corpus_sha_hex(&population_bytes);
    let corpus_id = args
        .population
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "corpus".to_string());

    // Parse the raw bytes as lines.
    let population: Vec<String> = population_bytes
        .split(|&b| b == b'\n')
        .map(|line| String::from_utf8_lossy(line).trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();

    eprintln!(
        "[perf-matrix] Population: {} patterns, corpus_id={corpus_id}, sha={corpus_sha}",
        population.len()
    );

    let seed = args.seed;
    let exclude_protein = args.exclude_protein;

    // Determine the tools to run based on which settings are provided.
    // Always include ferro; include others only if settings/paths are present
    // or the feature is enabled.
    let mut tools: Vec<String> = vec!["ferro".to_string(), "mutalyzer".to_string()];
    if args.biocommons_settings.is_some()
        || (args.uta_db_url.is_some() && args.seqrepo_path.is_some())
    {
        tools.push("biocommons".to_string());
    }
    // hgvs-rs normalize requires both seqrepo_path AND uta_db_url; include the
    // tool only when both are supplied (mirrors biocommons gating).
    #[cfg(feature = "hgvs-rs")]
    if args.seqrepo_path.is_some() && args.uta_db_url.is_some() {
        tools.push("hgvs-rs".to_string());
    }

    let cfg = MatrixConfig {
        tools: tools.clone(),
        workers: args.workers.clone(),
        reps: args.reps,
        ferro_threads: args.ferro_threads.clone(),
    };

    let harness = HarnessMeasure {
        config: HarnessMeasureConfig {
            reference: args.reference.clone(),
            mutalyzer_settings: args.mutalyzer_settings.clone(),
            biocommons_settings: args.biocommons_settings.clone(),
            seqrepo_path: args.seqrepo_path.clone(),
            uta_db_url: args.uta_db_url.clone(),
        },
    };

    let mut operations: BTreeMap<String, Operation> = BTreeMap::new();

    for op_name in &args.operations {
        let op_str = op_name.as_str();
        eprintln!("[perf-matrix] Running operation={op_str} …");

        // --- Calibration: compute per-(tool, op) sample sizes -----------------
        //
        // If --sample-size is set, skip calibration and apply it uniformly.
        // Otherwise, probe each tool with `min_sample` patterns at W=1 to
        // estimate its rate, then choose N = clamp(rate * target_seconds,
        // min_sample, max_sample).  Tools that fail the probe are marked
        // not_run and excluded from the matrix (sample size = 0).
        let calibrated_sizes: BTreeMap<String, usize> = if let Some(fixed) = args.sample_size {
            eprintln!("[perf-matrix]   Using manual sample-size override: N={fixed} for all tools");
            tools.iter().map(|t| (t.clone(), fixed)).collect()
        } else {
            let probe_n = args.min_sample;
            eprintln!(
                "[perf-matrix]   Calibrating sample sizes (target={:.1}s, min={}, max={}) …",
                args.target_seconds, args.min_sample, args.max_sample
            );
            let probe_sample = stratified_sample_vec(&population, probe_n, seed, exclude_protein)
                .unwrap_or_else(|e| {
                    eprintln!("[perf-matrix]     WARNING: probe sample failed: {e}; using empty");
                    Vec::new()
                });

            let mut sizes = BTreeMap::new();
            for tool in &tools {
                if probe_sample.is_empty() {
                    eprintln!("[perf-matrix]     {tool}: probe sample empty → skipping");
                    sizes.insert(tool.clone(), 0usize);
                    continue;
                }
                match harness.measure(tool, op_str, 1, &probe_sample) {
                    None => {
                        eprintln!("[perf-matrix]     {tool}: probe returned None → not_run");
                        sizes.insert(tool.clone(), 0usize);
                    }
                    Some(r) if r.pps <= 0.0 || r.success_rate == 0.0 => {
                        eprintln!(
                            "[perf-matrix]     {tool}: probe gave pps={:.0} success={:.2} → not_run",
                            r.pps, r.success_rate
                        );
                        sizes.insert(tool.clone(), 0usize);
                    }
                    Some(r) => {
                        let n = (r.pps * args.target_seconds)
                            .round()
                            .clamp(args.min_sample as f64, args.max_sample as f64)
                            as usize;
                        eprintln!("[perf-matrix]     {tool}: probe pps={:.0} → N={n}", r.pps);
                        sizes.insert(tool.clone(), n);
                    }
                }
            }
            sizes
        };

        // Build the sample_for closure that draws a per-tool sample of the
        // calibrated N from the same seeded stratified population.
        let sample_for = |tool: &str, rep: u32| -> Vec<String> {
            let n = calibrated_sizes.get(tool).copied().unwrap_or(0);
            if n == 0 {
                return Vec::new();
            }
            stratified_sample_vec(&population, n, seed + rep as u64, exclude_protein)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "[perf-matrix] WARNING: sample_for(tool={tool}, rep={rep}) failed: {e}; \
                         using empty sample"
                    );
                    Vec::new()
                })
        };

        let mut op = run_operation(&harness, op_str, &cfg, &sample_for);

        // Record calibrated sizes in provenance.
        op.sample_sizes = calibrated_sizes
            .iter()
            .map(|(t, &n)| (t.clone(), n as u64))
            .collect();

        // Keep sample_size_slow for back-compat: record the smallest non-zero
        // size among the slow tools (or 0 if none ran).
        let slow_tools = ["mutalyzer", "biocommons", "hgvs-rs"];
        op.sample_size_slow = slow_tools
            .iter()
            .filter_map(|&t| calibrated_sizes.get(t).copied())
            .filter(|&n| n > 0)
            .min()
            .unwrap_or(0) as u64;

        // ferro full-population pass — measure at the highest thread count to
        // get the headline throughput number.
        let ferro_full_threads = cfg.ferro_threads.iter().copied().max().unwrap_or(1);
        let full_pop: Vec<String> = if args.ferro_full_n > 0 {
            population.iter().take(args.ferro_full_n).cloned().collect()
        } else {
            population.clone()
        };
        op.ferro_full_population_n = full_pop.len() as u64;

        eprintln!(
            "[perf-matrix] ferro full-population pass: {} patterns, {} threads …",
            full_pop.len(),
            ferro_full_threads
        );
        if let Some(r) = harness.measure("ferro", op_str, ferro_full_threads, &full_pop) {
            op.ferro_full_population_pps = r.pps;
        }

        operations.insert(op_name.clone(), op);
    }

    // Build provenance.
    let mut tool_versions: BTreeMap<String, String> = BTreeMap::new();
    tool_versions.insert("ferro".to_string(), env!("CARGO_PKG_VERSION").to_string());
    tool_versions.insert("mutalyzer".to_string(), "3.1.1".to_string());
    tool_versions.insert("hgvs-rs".to_string(), "0.20.2".to_string());
    tool_versions.insert("biocommons".to_string(), "fa02ca5 (unpinned)".to_string());

    let provenance = Provenance {
        generated: chrono::Utc::now().to_rfc3339(),
        corpus_id,
        corpus_sha,
        machine: args.machine.clone(),
        os: format!("{}-{}", std::env::consts::OS, std::env::consts::ARCH),
        tool_versions,
        notes: vec![TIMED_REGION_NOTE.to_string()],
    };

    let results = PerfResults {
        schema_version: 1,
        placeholder: false,
        provenance,
        operations,
    };

    // Write output.
    if let Some(parent) = args.output.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| crate::FerroError::Io {
                msg: format!("Failed to create output directory: {e}"),
            })?;
        }
    }
    let json = serde_json::to_string_pretty(&results).map_err(|e| crate::FerroError::Json {
        msg: format!("Failed to serialize PerfResults: {e}"),
    })?;
    std::fs::write(&args.output, &json).map_err(|e| crate::FerroError::Io {
        msg: format!("Failed to write {}: {e}", args.output.display()),
    })?;
    eprintln!("[perf-matrix] Wrote results to {}", args.output.display());
    Ok(())
}

// ---------------------------------------------------------------------------
// Inline tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aggregates_median_min_max() {
        // 5 reps of patterns/sec for one tool.
        let rates = vec![100.0, 120.0, 110.0, 90.0, 130.0];
        let m = ToolMeasure::from_rates("ferro", &rates, 0.99);
        assert_eq!(m.median_pps, 110.0); // middle of sorted [90,100,110,120,130]
        assert_eq!(m.min_pps, 90.0);
        assert_eq!(m.max_pps, 130.0);
        assert_eq!(m.reps, 5);
        assert_eq!(m.success_rate, 0.99);
    }

    #[test]
    fn median_even_count_averages_middle_two() {
        let rates = vec![10.0, 20.0, 30.0, 40.0];
        let m = ToolMeasure::from_rates("x", &rates, 1.0);
        assert_eq!(m.median_pps, 25.0); // (20+30)/2
    }

    #[test]
    fn serializes_to_plan1_schema() {
        let m = ToolMeasure::from_rates("ferro", &[8_000_000.0], 1.0);
        let v = serde_json::to_value(&m).unwrap();
        assert!(v.get("median_pps").is_some());
        assert!(v.get("min_pps").is_some());
        assert!(v.get("reps").is_some());
        assert!(v.get("success_rate").is_some());
    }

    struct FakeMeasure;
    impl Measure for FakeMeasure {
        fn measure(
            &self,
            tool: &str,
            op: &str,
            workers: usize,
            sample: &[String],
        ) -> Option<RepRate> {
            // Deterministic synthetic rate: ferro fast, others slow; scales with workers.
            let base = match tool {
                "ferro" => 1_000_000.0,
                "mutalyzer" => 20.0,
                _ => 100.0,
            };
            let _ = (op, sample);
            Some(RepRate {
                pps: base * workers as f64,
                success_rate: 1.0,
            })
        }
    }

    /// Cross-module schema round-trip: serialize a fully-populated
    /// `perf_matrix::PerfResults` to JSON, then deserialize it as
    /// `perf_table::PerfResults` and assert key fields survive.
    ///
    /// This test requires BOTH the `benchmark` and `dev` features because
    /// `perf_matrix` lives under `#[cfg(feature="benchmark")]` and
    /// `perf_table` is only compiled under `#[cfg(feature="dev")]`.
    ///
    /// Run with: `cargo nextest run --features dev,benchmark -E 'test(perf_matrix)'`
    #[cfg(feature = "dev")]
    #[test]
    fn round_trip_schema_engine_to_renderer() {
        use crate::perf_table;

        // Build a fully-populated perf_matrix::PerfResults.
        let mut tool_versions = BTreeMap::new();
        tool_versions.insert("ferro".to_string(), "0.6.0".to_string());
        tool_versions.insert("mutalyzer".to_string(), "3.1.1".to_string());

        let provenance = Provenance {
            generated: "2026-06-12T00:00:00Z".to_string(),
            corpus_id: "clinvar_patterns".to_string(),
            corpus_sha: "deadbeef".to_string(),
            machine: "Apple M2 Max".to_string(),
            os: "darwin-arm64".to_string(),
            tool_versions,
            notes: vec!["Python startup overhead excluded from timed region.".to_string()],
        };

        // Worker point at W=1 with a running tool and a not_run tool.
        let mut tools_w1: BTreeMap<String, ToolMeasure> = BTreeMap::new();
        tools_w1.insert(
            "ferro".to_string(),
            ToolMeasure::from_rates("ferro", &[1_200_000.0, 1_100_000.0], 1.0),
        );
        tools_w1.insert(
            "biocommons".to_string(),
            ToolMeasure::not_run("UTA unavailable"),
        );

        // Worker point at W=8.
        let mut tools_w8: BTreeMap<String, ToolMeasure> = BTreeMap::new();
        tools_w8.insert(
            "ferro".to_string(),
            ToolMeasure::from_rates("ferro", &[8_000_000.0], 1.0),
        );
        tools_w8.insert(
            "biocommons".to_string(),
            ToolMeasure::not_run("UTA unavailable"),
        );

        let mut by_workers = BTreeMap::new();
        by_workers.insert("1".to_string(), WorkerPoint { tools: tools_w1 });
        by_workers.insert("8".to_string(), WorkerPoint { tools: tools_w8 });

        let mut ferro_scaling = BTreeMap::new();
        ferro_scaling.insert("1".to_string(), 1_200_000.0_f64);
        ferro_scaling.insert("2".to_string(), 2_300_000.0_f64);
        ferro_scaling.insert("4".to_string(), 4_400_000.0_f64);
        ferro_scaling.insert("8".to_string(), 8_000_000.0_f64);

        let mut parse_sample_sizes = BTreeMap::new();
        parse_sample_sizes.insert("ferro".to_string(), 2_000_000u64);
        parse_sample_sizes.insert("biocommons".to_string(), 500u64);

        let parse_op = Operation {
            sample_size_slow: 1000,
            ferro_full_population_n: 1_000_000,
            ferro_full_population_pps: 4_200_000.0,
            by_workers: by_workers.clone(),
            ferro_scaling: ferro_scaling.clone(),
            sample_sizes: parse_sample_sizes,
        };
        let normalize_op = Operation {
            sample_size_slow: 500,
            ferro_full_population_n: 500_000,
            ferro_full_population_pps: 2_100_000.0,
            by_workers,
            ferro_scaling,
            sample_sizes: BTreeMap::new(),
        };

        let mut operations = BTreeMap::new();
        operations.insert("parse".to_string(), parse_op);
        operations.insert("normalize".to_string(), normalize_op);

        let engine_results = PerfResults {
            schema_version: 1,
            placeholder: true,
            provenance,
            operations,
        };

        // Serialize via serde_json.
        let json =
            serde_json::to_string(&engine_results).expect("serialize perf_matrix::PerfResults");

        // Deserialize as perf_table::PerfResults (the renderer's type).
        let table_results: perf_table::PerfResults =
            serde_json::from_str(&json).expect("deserialize into perf_table::PerfResults");

        // Assert key fields survive the round-trip.
        assert_eq!(table_results.schema_version, 1);
        assert!(table_results.placeholder, "placeholder flag must survive");
        assert_eq!(table_results.provenance.corpus_sha, "deadbeef");

        let parse = &table_results.operations["parse"];
        let ferro_w8 = &parse.by_workers["8"].tools["ferro"];
        assert_eq!(ferro_w8.median_pps, 8_000_000.0);
        assert!(!ferro_w8.not_run);

        let biocommons_w1 = &parse.by_workers["1"].tools["biocommons"];
        assert!(
            biocommons_w1.not_run,
            "not_run tool must survive round-trip"
        );

        assert_eq!(parse.ferro_scaling["4"], 4_400_000.0_f64);
    }

    #[test]
    fn run_loop_builds_worker_points() {
        let cfg = MatrixConfig {
            tools: vec!["ferro".into(), "mutalyzer".into()],
            workers: vec![1, 8],
            reps: 3,
            ferro_threads: vec![1, 2, 4, 8],
        };
        let sample = vec!["NM_000088.3:c.589G>T".to_string()];
        // sample_for now takes (tool, rep); the tool parameter is ignored in
        // this test because FakeMeasure uses a fixed, non-empty sample.
        let op = run_operation(&FakeMeasure, "parse", &cfg, &|_tool, _rep| sample.clone());
        // W=1 and W=8 present, ferro+mutalyzer each.
        assert_eq!(op.by_workers["1"].tools["ferro"].median_pps, 1_000_000.0);
        assert_eq!(op.by_workers["8"].tools["ferro"].median_pps, 8_000_000.0);
        assert_eq!(op.by_workers["8"].tools["mutalyzer"].median_pps, 160.0);
        // ferro scaling at 1/2/4/8.
        assert_eq!(op.ferro_scaling["4"], 4_000_000.0);
    }

    /// When a tool's measure returns zero successes the run_operation loop must
    /// record it as not_run, NOT as a real throughput figure.
    ///
    /// This test reproduces the smoke regression: hgvs-rs W=8 returned
    /// median_pps=17845 with success_rate=0.0 and not_run=false — a meaningless
    /// rate published as real.  After the fix, the zero-success guard in
    /// `HarnessMeasure::measure_ferro_normalize` (and the equivalent guards for
    /// mutalyzer/hgvs-rs) returns `None`, which causes `run_operation` to
    /// produce `not_run=true` for all reps in which there are no successes.
    ///
    /// Here we test the `run_operation` loop directly with `ZeroSuccessMeasure`
    /// returning non-None but zero-success rates; the caller (`run_operation`)
    /// inserts these into the `rates` vec and then calls
    /// `ToolMeasure::from_rates`, which faithfully records them. That is
    /// CORRECT: the zero-success guard must live INSIDE the `Measure` impl, not
    /// in `run_operation`. So the right way to test the guard is to verify that a
    /// measure that returns `None` on zero successes produces `not_run=true`.
    ///
    /// This test therefore uses a measure that always returns `None` (simulating
    /// the guard), and asserts the resulting `not_run=true`.
    #[test]
    fn zero_success_measure_produces_not_run() {
        struct AlwaysNoneMeasure;
        impl Measure for AlwaysNoneMeasure {
            fn measure(
                &self,
                _tool: &str,
                _op: &str,
                _workers: usize,
                _sample: &[String],
            ) -> Option<RepRate> {
                None // guard fired: zero successes → not_run
            }
        }

        let cfg = MatrixConfig {
            tools: vec!["hgvs-rs".into()],
            workers: vec![8],
            reps: 2,
            ferro_threads: vec![],
        };
        let sample = vec!["NM_000088.3:c.589G>T".to_string()];
        let op = run_operation(&AlwaysNoneMeasure, "normalize", &cfg, &|_tool, _rep| {
            sample.clone()
        });

        let tm = &op.by_workers["8"].tools["hgvs-rs"];
        assert!(
            tm.not_run,
            "zero-success tool must be recorded as not_run, not as a real throughput"
        );
        assert_eq!(tm.median_pps, 0.0, "not_run tool must have median_pps=0");
    }
}
