//! Performance-comparison tables: types + rendering.
//!
//! Loads `data/benchmark/perf_results.json` (produced by the benchmark
//! measurement engine) and renders the README perf tables. See
//! `docs/superpowers/specs/2026-06-12-benchmark-perf-tables-design.md`.

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

/// Canonical column order for every table.
pub const TOOL_ORDER: [&str; 4] = ["ferro", "mutalyzer", "biocommons", "hgvs-rs"];

/// Provenance for one benchmark run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Provenance {
    pub generated: String,
    pub corpus_id: String,
    pub corpus_sha: String,
    pub machine: String,
    pub os: String,
    pub tool_versions: BTreeMap<String, String>,
}

/// One tool's measured throughput at a given worker count.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToolMeasure {
    pub median_pps: f64,
    pub min_pps: f64,
    pub max_pps: f64,
    pub reps: u32,
    pub success_rate: f64,
    /// Set when the tool's stack could not run (e.g. UTA unavailable).
    #[serde(default)]
    pub not_run: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
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
    /// Keys are decimal-integer strings; renderers iterate explicit `["1","2","4","8"]`,
    /// so lexicographic BTreeMap order does not affect output.
    #[serde(default)]
    pub ferro_scaling: BTreeMap<String, f64>,
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

impl PerfResults {
    /// A markdown warning to render above the tables when the data has not yet
    /// been measured (`placeholder == true`). `None` once real measurements
    /// replace the placeholder, so the warning auto-disappears.
    pub fn placeholder_banner(&self) -> Option<&'static str> {
        if self.placeholder {
            Some("> ⚠️ **Placeholder data — not yet measured.** These figures are illustrative pending a real benchmark run; do not cite them.")
        } else {
            None
        }
    }

    /// Parse from a JSON string.
    pub fn from_json_str(s: &str) -> Result<Self, String> {
        serde_json::from_str(s).map_err(|e| format!("parse perf_results: {e}"))
    }

    /// Load from a file path.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, String> {
        let p = path.as_ref();
        let s = std::fs::read_to_string(p).map_err(|e| format!("read {}: {e}", p.display()))?;
        Self::from_json_str(&s)
    }

    /// Render the cross-tool throughput table for an operation ("parse" or
    /// "normalize"): one row per tool present in the JSON (in `TOOL_ORDER`),
    /// columns for W=1 and W=8 throughput and ferro's speedup at W=8.
    pub fn render_cross_tool(&self, op: &str) -> Result<String, String> {
        let operation = self
            .operations
            .get(op)
            .ok_or_else(|| format!("no operation {op}"))?;
        let w1 = operation
            .by_workers
            .get("1")
            .ok_or_else(|| format!("op {op} missing worker point 1"))?;
        let w8 = operation
            .by_workers
            .get("8")
            .ok_or_else(|| format!("op {op} missing worker point 8"))?;

        let cell = |wp: &WorkerPoint, tool: &str| -> String {
            match wp.tools.get(tool) {
                None => "—".to_string(),
                Some(m) if m.not_run => "—".to_string(),
                Some(m) => fmt_pps(m.median_pps),
            }
        };

        let ferro8 = w8.tools.get("ferro").filter(|m| !m.not_run);

        let mut out = String::new();
        out.push_str(
            "| Tool | Throughput @ 1 worker | Throughput @ 8 workers | ferro speedup @ 8w |\n",
        );
        out.push_str(
            "|------|----------------------:|-----------------------:|-------------------:|\n",
        );
        for tool in TOOL_ORDER {
            // Skip tools absent from BOTH worker points entirely.
            if !w1.tools.contains_key(tool) && !w8.tools.contains_key(tool) {
                continue;
            }
            let speedup = if tool == "ferro" {
                "—".to_string()
            } else {
                match (ferro8, w8.tools.get(tool).filter(|m| !m.not_run)) {
                    (Some(f), Some(m)) => fmt_speedup(f.median_pps, m.median_pps),
                    _ => "—".to_string(),
                }
            };
            out.push_str(&format!(
                "| {} | {} | {} | {} |\n",
                tool,
                cell(w1, tool),
                cell(w8, tool),
                speedup
            ));
        }
        Ok(out)
    }

    /// Render ferro's thread-scaling row (1/2/4/8 threads) for an operation.
    /// Returns `Err` if the operation has no `ferro_scaling` data.
    pub fn render_ferro_scaling(&self, op: &str) -> Result<String, String> {
        let operation = self
            .operations
            .get(op)
            .ok_or_else(|| format!("no operation {op}"))?;
        if operation.ferro_scaling.is_empty() {
            return Err(format!("op {op} has no ferro_scaling data"));
        }
        let mut out = String::new();
        out.push_str("| Threads | 1 | 2 | 4 | 8 |\n");
        out.push_str("|---------|--:|--:|--:|--:|\n");
        out.push_str(&format!("| ferro {op} |"));
        for t in ["1", "2", "4", "8"] {
            let v = operation
                .ferro_scaling
                .get(t)
                .copied()
                .map(fmt_pps)
                .unwrap_or_else(|| "—".to_string());
            out.push_str(&format!(" {v} |"));
        }
        out.push('\n');
        Ok(out)
    }
}

/// Human-readable throughput: `>=1e6` as `X.YM/s`, `>=1e3` as `X.Yk/s`,
/// `>=1` as integer `/s`, else one decimal `/s`.
pub fn fmt_pps(pps: f64) -> String {
    if pps >= 1_000_000.0 {
        format!("{:.1}M/s", pps / 1_000_000.0)
    } else if pps >= 1_000.0 {
        format!("{:.1}k/s", pps / 1_000.0)
    } else if pps >= 1.0 {
        format!("{:.0}/s", pps)
    } else {
        format!("{:.1}/s", pps)
    }
}

/// Speedup `ferro/tool`, rounded to two significant figures and grouped with
/// thousands separators, suffixed `×`. Returns `1×` for equal rates.
pub fn fmt_speedup(ferro_pps: f64, tool_pps: f64) -> String {
    if tool_pps <= f64::EPSILON {
        return "—".to_string();
    }
    let ratio = ferro_pps / tool_pps;
    if !ratio.is_finite() || ratio <= 0.0 {
        return "—".to_string();
    }
    let rounded = round_2sig(ratio);
    format!("{}×", group_thousands(rounded))
}

/// Round to 2 significant figures. Returns 0 for non-finite or non-positive input.
fn round_2sig(x: f64) -> u64 {
    if x <= 0.0 || !x.is_finite() {
        return 0;
    }
    let mag = x.log10().floor();
    let scale = 10f64.powf(mag - 1.0);
    let rounded = ((x / scale).round() * scale).round();
    rounded.min(u64::MAX as f64) as u64
}

/// Group an integer with `,` thousands separators.
fn group_thousands(n: u64) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let mut out = String::new();
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(*b as char);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = r#"{
      "schema_version": 1,
      "provenance": {
        "generated": "2026-06-12T00:00:00Z", "corpus_id": "clinvar_patterns",
        "corpus_sha": "abc123", "machine": "Apple M2 Max (12 cores)", "os": "darwin-arm64",
        "tool_versions": {"ferro":"0.6.0","mutalyzer":"3.1.1","biocommons":"fa02ca5","hgvs-rs":"0.20.2"}
      },
      "operations": {
        "parse": {
          "sample_size_slow": 1000, "ferro_full_population_n": 1000000,
          "ferro_full_population_pps": 4200000.0,
          "by_workers": {
            "1": {"tools": {
              "ferro": {"median_pps": 1200000.0, "min_pps": 1150000.0, "max_pps": 1250000.0, "reps": 10, "success_rate": 1.0},
              "mutalyzer": {"median_pps": 21.0, "min_pps": 20.0, "max_pps": 22.0, "reps": 5, "success_rate": 0.994}
            }},
            "8": {"tools": {
              "ferro": {"median_pps": 8000000.0, "min_pps": 7800000.0, "max_pps": 8200000.0, "reps": 10, "success_rate": 1.0},
              "mutalyzer": {"median_pps": 150.0, "min_pps": 140.0, "max_pps": 160.0, "reps": 5, "success_rate": 0.994}
            }}
          },
          "ferro_scaling": {"1": 1200000.0, "2": 2300000.0, "4": 4400000.0, "8": 8000000.0}
        }
      }
    }"#;

    #[test]
    fn loads_results() {
        let r = PerfResults::from_json_str(SAMPLE).expect("parse");
        assert_eq!(r.schema_version, 1);
        let parse = &r.operations["parse"];
        assert_eq!(parse.by_workers["8"].tools["ferro"].median_pps, 8000000.0);
        assert_eq!(parse.ferro_scaling["4"], 4400000.0);
    }

    #[test]
    fn formats_throughput() {
        assert_eq!(fmt_pps(4_200_000.0), "4.2M/s");
        assert_eq!(fmt_pps(8_000_000.0), "8.0M/s");
        assert_eq!(fmt_pps(1_200_000.0), "1.2M/s");
        assert_eq!(fmt_pps(150.0), "150/s");
        assert_eq!(fmt_pps(21.0), "21/s");
        assert_eq!(fmt_pps(0.2), "0.2/s");
    }

    #[test]
    fn formats_speedup() {
        assert_eq!(fmt_speedup(8_000_000.0, 150.0), "53,000×");
        assert_eq!(fmt_speedup(1_200_000.0, 21.0), "57,000×");
        assert_eq!(fmt_speedup(100.0, 100.0), "1×");
    }

    #[test]
    fn renders_cross_tool_table() {
        let r = PerfResults::from_json_str(SAMPLE).unwrap();
        let out = r.render_cross_tool("parse").unwrap();
        assert!(out.contains(
            "| Tool | Throughput @ 1 worker | Throughput @ 8 workers | ferro speedup @ 8w |"
        ));
        // ferro is the baseline row (no speedup vs itself).
        assert!(out.contains("| ferro | 1.2M/s | 8.0M/s | — |"));
        // mutalyzer: 8000000/150 -> ~53,000x.
        assert!(out.contains("| mutalyzer | 21/s | 150/s | 53,000× |"));
        // Tools absent from the JSON are omitted (only ferro+mutalyzer in SAMPLE).
        assert!(!out.contains("biocommons"));
    }

    #[test]
    fn renders_not_run_tool() {
        let json = SAMPLE.replace(
            r#""mutalyzer": {"median_pps": 150.0, "min_pps": 140.0, "max_pps": 160.0, "reps": 5, "success_rate": 0.994}"#,
            r#""mutalyzer": {"median_pps": 0.0, "min_pps": 0.0, "max_pps": 0.0, "reps": 0, "success_rate": 0.0, "not_run": true, "reason": "UTA unavailable"}"#,
        );
        let r = PerfResults::from_json_str(&json).unwrap();
        let out = r.render_cross_tool("parse").unwrap();
        assert!(out.contains("| mutalyzer | 21/s | — | — |")); // W=1 ran, W=8 not_run
    }

    #[test]
    fn renders_ferro_scaling() {
        let r = PerfResults::from_json_str(SAMPLE).unwrap();
        let out = r.render_ferro_scaling("parse").unwrap();
        assert!(out.contains("| Threads | 1 | 2 | 4 | 8 |"));
        assert!(out.contains("| ferro parse | 1.2M/s | 2.3M/s | 4.4M/s | 8.0M/s |"));
    }

    #[test]
    fn speedup_edge_cases() {
        // Sub-1 ratio (ferro slower): round_2sig(0.5) -> 1 (Rust rounds half away from zero).
        assert_eq!(fmt_speedup(50.0, 100.0), "1×");
        // Zero / non-finite guards render the em-dash.
        assert_eq!(fmt_speedup(100.0, 0.0), "—");
        assert_eq!(fmt_speedup(f64::INFINITY, 1.0), "—");
        assert_eq!(fmt_speedup(f64::NAN, 1.0), "—");
    }

    #[test]
    fn ferro_scaling_missing_thread_renders_dash() {
        // Drop the "2" key -> that column renders "—".
        let json = SAMPLE.replace(r#""2": 2300000.0, "#, "");
        let r = PerfResults::from_json_str(&json).unwrap();
        let out = r.render_ferro_scaling("parse").unwrap();
        assert!(out.contains("| ferro parse | 1.2M/s | — | 4.4M/s | 8.0M/s |"));
    }

    #[test]
    fn placeholder_banner_toggles() {
        let mut r = PerfResults::from_json_str(SAMPLE).unwrap();
        assert!(r.placeholder_banner().is_none()); // SAMPLE has no placeholder field -> default false
        r.placeholder = true;
        assert!(r.placeholder_banner().unwrap().contains("Placeholder data"));
    }
}
