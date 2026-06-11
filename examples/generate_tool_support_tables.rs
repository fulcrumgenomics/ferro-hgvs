//! Generator for the tool-support comparison tables.
//!
//! Reads `docs/tool_support_matrix.json` (the single source of truth) and
//! renders:
//!   - the "Normalization Capabilities" markdown table into README.md and
//!     docs/BENCHMARK_GUIDE.md (between `tool-support:<view>` markers), and
//!   - the render-ready website JSON into
//!     src/service/web/static/data/tool_support_matrix.json.
//!
//! Modes:
//!   - default:  rewrite the targets in place.
//!   - --check:  re-render in memory; exit non-zero if any target differs.
//!
//! Run:   cargo run --features dev --example generate_tool_support_tables
//! Check: cargo run --features dev --example generate_tool_support_tables -- --check
//!
//! See: docs/superpowers/specs/2026-06-09-tool-support-matrix-design.md

use std::path::Path;
use std::process::ExitCode;

use clap::Parser;
use ferro_hgvs::tool_support::Matrix;

/// `today` is injected (not read from the clock) so the staleness check is
/// reproducible; CI passes the commit date or leaves the default.
#[derive(Parser, Debug)]
#[command(about = "Generate tool-support comparison tables")]
struct Cli {
    /// Re-render in memory and fail if on-disk output differs.
    #[arg(long)]
    check: bool,
    /// Path to the source matrix.
    #[arg(long, default_value = "docs/tool_support_matrix.json")]
    matrix: String,
    /// Staleness reference date (YYYY-MM-DD); warnings only. Frozen at the
    /// matrix authoring date so output is reproducible; pass --today=$(date +%F)
    /// (or the commit date) to check staleness against the present.
    #[arg(long, default_value = "2026-06-09")]
    today: String,
    /// Staleness window in months.
    #[arg(long, default_value_t = 6)]
    stale_months: i64,
}

const VIEW: &str = "normalization_capabilities";
const README: &str = "README.md";
const GUIDE: &str = "docs/BENCHMARK_GUIDE.md";
const WEB_JSON: &str = "src/service/web/static/data/tool_support_matrix.json";

fn marker_begin() -> String {
    format!("<!-- BEGIN tool-support:{VIEW} -->")
}
fn marker_end() -> String {
    format!("<!-- END tool-support:{VIEW} -->")
}

/// Replace the text between the begin/end markers with `body`. Returns the new
/// file content, or an error if markers are missing/disordered.
fn splice(content: &str, body: &str) -> Result<String, String> {
    let begin = marker_begin();
    let end = marker_end();
    let b = content
        .find(&begin)
        .ok_or_else(|| format!("missing {begin}"))?;
    let after_begin = b + begin.len();
    let e = content[after_begin..]
        .find(&end)
        .ok_or_else(|| format!("missing {end} after {begin}"))?
        + after_begin;
    let mut out = String::new();
    out.push_str(&content[..after_begin]);
    out.push('\n');
    out.push_str(body.trim_end());
    out.push('\n');
    out.push_str(&content[e..]);
    Ok(out)
}

/// One generated target: path + desired content.
struct Target {
    path: String,
    content: String,
}

fn compute_targets(m: &Matrix) -> Result<Vec<Target>, String> {
    let body = m.render_markdown_view(VIEW)?;
    let mut targets = Vec::new();
    for path in [README, GUIDE] {
        let cur = std::fs::read_to_string(path).map_err(|e| format!("read {path}: {e}"))?;
        targets.push(Target {
            path: path.to_string(),
            content: splice(&cur, &body)?,
        });
    }
    let web = serde_json::to_string_pretty(&m.render_website_json())
        .map_err(|e| format!("serialize web json: {e}"))?;
    targets.push(Target {
        path: WEB_JSON.to_string(),
        content: format!("{web}\n"),
    });
    Ok(targets)
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let m = match Matrix::load(&cli.matrix) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };
    if let Err(e) = m.validate_invariants() {
        eprintln!("error: matrix invariants: {e}");
        return ExitCode::FAILURE;
    }

    // Staleness warnings (never fail).
    for s in m.stale_cells(&cli.today, cli.stale_months) {
        eprintln!(
            "warning: stale tool-support cell {s} (>{} months)",
            cli.stale_months
        );
    }

    let targets = match compute_targets(&m) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };

    let mut drift = false;
    for t in &targets {
        let on_disk = std::fs::read_to_string(&t.path).unwrap_or_default();
        if cli.check {
            if on_disk != t.content {
                drift = true;
                eprintln!("error: {} is out of date", t.path);
            }
        } else if on_disk != t.content {
            if let Some(parent) = Path::new(&t.path).parent() {
                if let Err(e) = std::fs::create_dir_all(parent) {
                    eprintln!("error: create dir {}: {e}", parent.display());
                    return ExitCode::FAILURE;
                }
            }
            if let Err(e) = std::fs::write(&t.path, &t.content) {
                eprintln!("error: write {}: {e}", t.path);
                return ExitCode::FAILURE;
            }
            println!("updated {}", t.path);
        }
    }

    if cli.check && drift {
        eprintln!(
            "tool-support tables are out of date; rerun: \
             cargo run --features dev --example generate_tool_support_tables \
             (edit docs/tool_support_matrix.json, not the generated files)"
        );
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}
