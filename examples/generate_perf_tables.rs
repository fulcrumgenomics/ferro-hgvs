//! Generator for the README performance tables.
//!
//! Reads `data/benchmark/perf_results.json` and renders the parse + normalize
//! cross-tool tables and the ferro thread-scaling table into README.md between
//! `perf:<id>` markers.
//!
//! Run:   cargo run --features dev --example generate_perf_tables
//! Check: cargo run --features dev --example generate_perf_tables -- --check
//!
//! See: docs/superpowers/specs/2026-06-12-benchmark-perf-tables-design.md

use std::process::ExitCode;

use clap::Parser;
use ferro_hgvs::perf_table::PerfResults;

#[derive(Parser, Debug)]
#[command(about = "Generate the README performance tables")]
struct Cli {
    /// Re-render in memory and fail if README.md differs.
    #[arg(long)]
    check: bool,
    /// Path to the results JSON.
    #[arg(long, default_value = "data/benchmark/perf_results.json")]
    results: String,
}

const README: &str = "README.md";

/// Replace the text between `<!-- BEGIN perf:<id> -->` and `<!-- END perf:<id> -->`.
fn splice(content: &str, id: &str, body: &str) -> Result<String, String> {
    let begin = format!("<!-- BEGIN perf:{id} -->");
    let end = format!("<!-- END perf:{id} -->");
    let b = content
        .find(&begin)
        .ok_or_else(|| format!("missing {begin}"))?;
    let after = b + begin.len();
    let e = content[after..]
        .find(&end)
        .ok_or_else(|| format!("missing {end} after {begin}"))?
        + after;
    let mut out = String::new();
    out.push_str(&content[..after]);
    out.push('\n');
    out.push_str(body.trim_end());
    out.push('\n');
    out.push_str(&content[e..]);
    Ok(out)
}

fn render_all(r: &PerfResults, current: &str) -> Result<String, String> {
    let banner = r.placeholder_banner();
    let with_banner = |table: String| match banner {
        Some(b) => format!("{b}\n\n{table}"),
        None => table,
    };
    let mut out = current.to_string();
    out = splice(&out, "parse", &with_banner(r.render_cross_tool("parse")?))?;
    out = splice(
        &out,
        "normalize",
        &with_banner(r.render_cross_tool("normalize")?),
    )?;
    out = splice(
        &out,
        "ferro_scaling",
        &with_banner(r.render_ferro_scaling("parse")?),
    )?;
    Ok(out)
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let r = match PerfResults::load(&cli.results) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };
    let current = match std::fs::read_to_string(README) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("error: read {README}: {e}");
            return ExitCode::FAILURE;
        }
    };
    let rendered = match render_all(&r, &current) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };

    if cli.check {
        if rendered != current {
            eprintln!(
                "README perf tables are out of date; rerun: \
                 cargo run --features dev --example generate_perf_tables \
                 (edit data/benchmark/perf_results.json, not README.md)"
            );
            return ExitCode::FAILURE;
        }
    } else if rendered != current {
        if let Err(e) = std::fs::write(README, &rendered) {
            eprintln!("error: write {README}: {e}");
            return ExitCode::FAILURE;
        }
        println!("updated {README}");
    }
    ExitCode::SUCCESS
}
