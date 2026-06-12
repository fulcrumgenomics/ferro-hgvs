//! Generator for the conformance summary docs (`failure-patterns.md`).
//!
//! Renders each normalize corpus's `failure-patterns.md` from the disposition
//! annotations + cluster taxonomy in its `cases.json` (issue #509). The doc is a
//! derived view — never hand-edited — so its counts and per-row dispositions
//! cannot drift from the schema the harness enforces.
//!
//! Two modes:
//!   - default:  regenerate each `failure-patterns.md` in place
//!   - --check:  regenerate in memory and exit non-zero if the committed file
//!     differs (used in CI to catch drift)
//!
//! Run:   `cargo run --features dev --example generate_conformance_summary`
//! Check: `cargo run --features dev --example generate_conformance_summary -- --check`

use std::path::Path;
use std::process::ExitCode;

use clap::Parser;
use ferro_hgvs::conformance::{biocommons, hgvs_rs_projection, mutalyzer, summary};

#[derive(Parser, Debug)]
#[command(about = "Generate the normalize conformance summary docs from cases.json")]
struct Cli {
    /// Check that each committed `failure-patterns.md` is byte-identical to a
    /// fresh regeneration. Exits non-zero on drift; writes nothing.
    #[arg(long)]
    check: bool,
}

/// One corpus's input `cases.json` and output `failure-patterns.md`.
struct Corpus {
    cases_json: &'static str,
    summary_md: &'static str,
}

const MUTALYZER: Corpus = Corpus {
    cases_json: "tests/fixtures/mutalyzer-normalize/cases.json",
    summary_md: "tests/fixtures/mutalyzer-normalize/failure-patterns.md",
};

const BIOCOMMONS: Corpus = Corpus {
    cases_json: "tests/fixtures/biocommons-normalize/cases.json",
    summary_md: "tests/fixtures/biocommons-normalize/failure-patterns.md",
};

const HGVS_RS_PROJECTION: Corpus = Corpus {
    cases_json: "tests/fixtures/hgvs-rs-projection/cases.json",
    summary_md: "tests/fixtures/hgvs-rs-projection/failure-patterns.md",
};

fn main() -> ExitCode {
    let cli = Cli::parse();
    let mut all_ok = true;
    for result in [
        mutalyzer_model(&MUTALYZER).and_then(|m| emit(&MUTALYZER, &m, cli.check)),
        biocommons_model(&BIOCOMMONS).and_then(|m| emit(&BIOCOMMONS, &m, cli.check)),
        hgvs_rs_projection_model(&HGVS_RS_PROJECTION)
            .and_then(|m| emit(&HGVS_RS_PROJECTION, &m, cli.check)),
    ] {
        if let Err(e) = result {
            eprintln!("{e}");
            all_ok = false;
        }
    }
    if all_ok {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

/// Load + validate the mutalyzer corpus into a summary model.
fn mutalyzer_model(corpus: &Corpus) -> Result<summary::SummaryModel, String> {
    let content = std::fs::read_to_string(corpus.cases_json)
        .map_err(|e| format!("read {}: {e}", corpus.cases_json))?;
    let fixture: mutalyzer::Fixture =
        serde_json::from_str(&content).map_err(|e| format!("parse {}: {e}", corpus.cases_json))?;
    fixture.validate_clusters()?;
    Ok(fixture.to_summary())
}

/// Load + validate the biocommons corpus into a summary model.
fn biocommons_model(corpus: &Corpus) -> Result<summary::SummaryModel, String> {
    let content = std::fs::read_to_string(corpus.cases_json)
        .map_err(|e| format!("read {}: {e}", corpus.cases_json))?;
    let fixture: biocommons::Fixture =
        serde_json::from_str(&content).map_err(|e| format!("parse {}: {e}", corpus.cases_json))?;
    fixture.validate_clusters()?;
    Ok(fixture.to_summary())
}

/// Load + validate the hgvs-rs-projection corpus into a summary model.
fn hgvs_rs_projection_model(corpus: &Corpus) -> Result<summary::SummaryModel, String> {
    let content = std::fs::read_to_string(corpus.cases_json)
        .map_err(|e| format!("read {}: {e}", corpus.cases_json))?;
    let fixture: hgvs_rs_projection::Fixture =
        serde_json::from_str(&content).map_err(|e| format!("parse {}: {e}", corpus.cases_json))?;
    fixture.validate_clusters()?;
    Ok(fixture.to_summary())
}

/// Write or `--check` the rendered summary for one corpus.
fn emit(corpus: &Corpus, model: &summary::SummaryModel, check: bool) -> Result<(), String> {
    let path = Path::new(corpus.summary_md);
    if check {
        summary::check(path, model)
    } else {
        std::fs::write(path, summary::render(model))
            .map_err(|e| format!("write {}: {e}", corpus.summary_md))
    }
}
