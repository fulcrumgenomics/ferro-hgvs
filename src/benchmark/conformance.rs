// Copyright (c) 2024-2025 Fulcrum Genomics LLC
// SPDX-License-Identifier: MIT

//! The `ferro-benchmark conformance census` command (#1890).
//!
//! Thin glue over [`crate::conformance::census`]: it maps the CLI enums to the
//! library relations, runs the census, prints a human summary, and writes the
//! machine-readable artifact, the findings rows, and the burn-down. All of the
//! measurement — and both refusals, the armed seam oracle and the incomplete
//! corpus — live in the library module; this file does I/O and turns a
//! `CensusError` into a message and a non-zero exit.

use std::io::Write;
use std::path::Path;

use crate::benchmark::cli::{CensusCorpus, CensusDirection, CensusEquivalence};
use crate::conformance::census::{
    compare_reports, report, run_census, CensusReport, CorpusSource, CounterMovement,
    DirectionOfVirtue, Equivalence, Movement,
};
use crate::FerroError;
use crate::ShuffleDirection;

impl From<CensusEquivalence> for Equivalence {
    fn from(value: CensusEquivalence) -> Self {
        match value {
            CensusEquivalence::Sequence => Equivalence::Sequence,
            CensusEquivalence::Partition => Equivalence::Partition,
            CensusEquivalence::Spdi => Equivalence::Spdi,
        }
    }
}

impl From<CensusCorpus> for CorpusSource {
    fn from(value: CensusCorpus) -> Self {
        match value {
            CensusCorpus::Spec => CorpusSource::Spec,
        }
    }
}

impl From<CensusDirection> for ShuffleDirection {
    fn from(value: CensusDirection) -> Self {
        match value {
            CensusDirection::ThreePrime => ShuffleDirection::ThreePrime,
        }
    }
}

/// Run `ferro-benchmark conformance census`.
///
/// Returns `Ok(())` on a successful census. A seam oracle in the environment, or a
/// corpus that dropped more than its allowance permits, is a refusal rather than an
/// error to propagate: it prints a clear message and exits non-zero, the same shape
/// the other refusing commands in this binary use, so a flattering artifact is
/// never written.
///
/// **Both refusals are the library's**, not this function's. The oracle check used
/// to be duplicated here — an `eprintln!` restating `CensusError::OracleArmed`'s
/// own `Display` — which left the variant unconstructible and, worse, left the
/// protection attached to the CLI rather than to the census: an in-process caller
/// of `run_census` got none of it. Both now arrive as an `Err` from
/// [`run_census`], and this function's only job is to print it and exit.
pub fn run_census_command(
    equivalence: CensusEquivalence,
    corpus: CensusCorpus,
    direction: CensusDirection,
    out: Option<&Path>,
    findings: Option<&Path>,
    compare: Option<&Path>,
) -> Result<(), FerroError> {
    let relation: Equivalence = equivalence.into();
    let run = match run_census(corpus.into(), direction.into(), relation) {
        Ok(run) => run,
        Err(err) => {
            eprintln!("Error: {err}");
            std::process::exit(2);
        }
    };

    // The human summary goes to stderr so `--out -`-style piping of the JSON stays
    // possible and the machine artifact is never mixed with prose.
    eprint!("{}", report(&run.report.direction, &run.measured));

    if let Some(path) = out {
        let json = serde_json::to_string_pretty(&run.report).map_err(|e| FerroError::Json {
            msg: format!("serializing the census to {}: {e}", path.display()),
        })?;
        std::fs::write(path, json).map_err(|e| FerroError::Io {
            msg: format!("writing the census to {}: {e}", path.display()),
        })?;
        eprintln!("wrote census to {}", path.display());
    }

    if let Some(path) = findings {
        write_findings(path, &run.measured)?;
        eprintln!("wrote findings to {}", path.display());
    }

    if let Some(path) = compare {
        let baseline = load_report(path)?;
        print_burn_down(&baseline, &run.report);
    }

    Ok(())
}

/// Write every finding and divergence as JSONL, one row per line — the row-level
/// detail the totals in `--out` cannot carry.
///
/// **Every** is load-bearing and is what the flag's own help promises. It used to
/// be a lie in one direction: `Measured::divergences` was truncated to the
/// summary's `MAX_DIVERGENCES = 12`, so an `--equivalence partition` run that
/// found 721 divergent families exported 12 of them, with nothing in the file or
/// in `--out` saying so. #1891 diffs two releases' findings files to name *which*
/// families moved, and two truncated files hold the same corpus-order-first 12
/// whatever happened behind them — a release in which 400 families changed read as
/// no change. The cap is now applied where it was sized to be applied, on the
/// human summary, and this export carries the whole list.
///
/// Each row carries a `record` discriminator (`"finding"` / `"divergence"`) so a
/// consumer reads a tagged union instead of sniffing which keys are present.
fn write_findings(
    path: &Path,
    measured: &crate::conformance::census::Measured,
) -> Result<(), FerroError> {
    let file = std::fs::File::create(path).map_err(|e| FerroError::Io {
        msg: format!("creating the findings file {}: {e}", path.display()),
    })?;
    let mut writer = std::io::BufWriter::new(file);
    let io = |e: std::io::Error| FerroError::Io {
        msg: format!("writing findings to {}: {e}", path.display()),
    };
    let json = |e: serde_json::Error| FerroError::Json {
        msg: format!("serializing a finding for {}: {e}", path.display()),
    };
    for finding in &measured.findings {
        let line = serde_json::to_string(finding).map_err(json)?;
        writeln!(writer, "{line}").map_err(io)?;
    }
    for divergence in &measured.divergences {
        let line = serde_json::to_string(divergence).map_err(json)?;
        writeln!(writer, "{line}").map_err(io)?;
    }
    writer.flush().map_err(io)?;
    Ok(())
}

/// Load a previously written census artifact.
fn load_report(path: &Path) -> Result<CensusReport, FerroError> {
    let raw = std::fs::read_to_string(path).map_err(|e| FerroError::Io {
        msg: format!("reading the comparison census {}: {e}", path.display()),
    })?;
    serde_json::from_str(&raw).map_err(|e| FerroError::Json {
        msg: format!("parsing the comparison census {}: {e}", path.display()),
    })
}

/// Print a per-counter burn-down with the direction of virtue attached. No overall
/// verdict and no percentage — rank-1 and rank-2 are not commensurable.
fn print_burn_down(before: &CensusReport, after: &CensusReport) {
    eprintln!(
        "\nburn-down: {} ({}) -> {} ({})",
        before.equivalence, before.direction, after.equivalence, after.direction
    );
    if before.equivalence != after.equivalence {
        eprintln!(
            "  NOTE: the two runs used different equivalence relations, so a confluence movement \
             is a change of question, not of behaviour — the confluence counters below are \
             therefore reported as `moved`, ungraded"
        );
    }
    for movement in compare_reports(before, after) {
        eprintln!("  {}", format_movement(&movement));
    }
}

/// One movement line, e.g.
/// `converged                       9140 -> 3509  WORSE     (higher is better)`.
///
/// The `moved` / `unchanged` distinction is the whole reason this is a named
/// function with a test of its own — see
/// `format_movement_distinguishes_a_neutral_move_from_no_movement` below. A
/// neutral counter that actually moved must not read as "unchanged", because that
/// describes the *value* and the value changed. A denominator that moved means the
/// corpus changed under the comparison, and that is the one thing a reader must
/// not miss while scanning a column of verdicts.
pub(crate) fn format_movement(movement: &CounterMovement) -> String {
    let verdict = match movement.movement {
        Movement::Improved => "IMPROVED",
        Movement::Worse => "WORSE",
        Movement::Moved => "moved",
        Movement::Unchanged => "unchanged",
    };
    let virtue = match movement.virtue {
        DirectionOfVirtue::LowerIsBetter => "lower is better",
        DirectionOfVirtue::HigherIsBetter => "higher is better",
        DirectionOfVirtue::Neutral => "neutral: a population, not a verdict",
    };
    format!(
        "{:<42} {:>7} -> {:<7} {:<10} ({virtue})",
        movement.name, movement.before, movement.after, verdict
    )
}

#[cfg(test)]
mod tests {
    use super::format_movement;
    use crate::conformance::census::{classify_movement, CounterMovement, DirectionOfVirtue};

    fn line(name: &'static str, before: usize, after: usize, virtue: DirectionOfVirtue) -> String {
        format_movement(&CounterMovement {
            name,
            before,
            after,
            virtue,
            movement: classify_movement(before, after, virtue),
        })
    }

    /// The one distinction this function exists to carry: a **neutral** counter
    /// that moved says `moved`, and only a counter that did not move says
    /// `unchanged`.
    ///
    /// It had no test. A neutral counter is a denominator or an instrument, so its
    /// movement is exactly the "the corpus changed under this comparison" signal a
    /// reader must not scan past — and the two words differ by one column in a
    /// wall of output, which is the kind of thing that regresses silently.
    #[test]
    fn format_movement_distinguishes_a_neutral_move_from_no_movement() {
        let moved = line("outputs", 5, 7, DirectionOfVirtue::Neutral);
        assert!(
            moved.contains("moved") && !moved.contains("unchanged"),
            "a neutral counter that moved must not read as unchanged: {moved}"
        );
        assert!(
            moved.contains("neutral: a population, not a verdict"),
            "the movement must carry its (absent) direction of virtue: {moved}"
        );

        let still = line("outputs", 5, 5, DirectionOfVirtue::Neutral);
        assert!(
            still.contains("unchanged"),
            "a counter that did not move must read as unchanged: {still}"
        );
    }

    /// A graded counter still reads IMPROVED/WORSE, in both directions of virtue.
    #[test]
    fn format_movement_grades_a_counter_that_has_a_direction_of_virtue() {
        assert!(line("declined", 9, 3, DirectionOfVirtue::LowerIsBetter).contains("IMPROVED"));
        assert!(line("declined", 3, 9, DirectionOfVirtue::LowerIsBetter).contains("WORSE"));
        assert!(line("converged", 3, 9, DirectionOfVirtue::HigherIsBetter).contains("IMPROVED"));
        assert!(line("converged", 9, 3, DirectionOfVirtue::HigherIsBetter).contains("WORSE"));
    }
}
