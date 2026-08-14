//! Measure the HGVS spec corpus against the spec's own stated form, **once per
//! `FERRO_PARTITION` arm**, over a prepared reference.
//!
//! # Why this exists
//!
//! `generate_spec_fixture` normalizes with an empty `MockProvider`, so no row
//! can fetch a block and `merge::canonicalize_from_sequence` — the only consumer
//! of the `FERRO_PARTITION` switch — is never reached. The fixture is therefore
//! **byte-identical on all four arms**, and every assertion the spec suite makes
//! is constant across them by construction.
//!
//! **That blindness is not yet pinned by anything in this tree.** The guard is
//! `it::hgvs_spec_normalization_tests::the_spec_fixture_is_blind_to_the_partition_switch`,
//! which lands with #1822 — open, not merged — so until it does, the byte-identity
//! above rests on the measurement recorded in this PR's description rather than on
//! a test. Do not cite it as coverage before checking it exists:
//!
//! ```text
//! grep -rn the_spec_fixture_is_blind_to_the_partition_switch tests/
//! ```
//!
//! This generator is the measurement that closes the blindness, and it
//! deliberately does not change the fixture or that suite.
//!
//! # What it measures, stated so nobody reads it as something else
//!
//! **Conformance = agreement with `spec_expected`, per arm.** `spec_expected` is
//! the *spec's* form and does not move with the arm, so four runs measure ferro
//! against the spec four times. The delta between two runs is the flip's
//! conformance effect. This is **not** an arm-versus-arm comparison and must not
//! be reported as one — two arms can disagree with each other while agreeing
//! with the spec equally often.
//!
//! # `spec_expected`'s provenance
//!
//! Read off `build_rows` in `generate_spec_fixture`, in precedence order:
//!
//! 1. an `hgvs_spec_normalization_overrides.json` `by_input` entry — curated;
//! 2. `None` when the spec wraps the string in `<code class="invalid">` —
//!    structural, extracted from the spec text;
//! 3. otherwise the harvested spec string itself, with a mechanical default
//!    accession prepended for the bare fragments the spec writes as shorthand.
//!
//! Case 3 carries the overwhelming majority, and none of the three consults
//! ferro's normalizer — so comparing ferro to this column is not circular. The
//! **one** exception is curated and named: `NP_060250.2:p.Gln746_Lys747ins*63`,
//! whose override pins ferro's own `insTer63` rendering with a note saying so.
//! It is reported separately (`ferro_derived_expectations`) and excluded from the
//! headline numerator, so a reader can take the measurement with or without it.
//!
//! # Effective denominator
//!
//! A conformance delta of zero is only evidence if rows actually reached the
//! switch. `normalize::partition_blocks_cut()` counts blocks cut on **every**
//! arm (`partition_decline_counts` cannot — it is structurally zero on `live`),
//! and this generator reads it either side of each row, so `reached_partitioner`
//! is a measured count and not a replica of the gate predicate. A run whose
//! `reached_partitioner` is zero has measured nothing about the arms, whatever
//! its census says.
//!
//! Run (one process per arm — the arm is read once into a `OnceLock`):
//!
//! ```text
//! FERRO_MANIFEST=<manifest> FERRO_PARTITION=canonical-coalesced \
//!   cargo run --features dev --example measure_spec_conformance_per_arm -- \
//!   --fixture <fixture.json> --output <arm.json>
//! ```

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use serde::{Deserialize, Serialize};

use ferro_hgvs::conformance::completeness::{CaptureCounts, CaptureLedger};
use ferro_hgvs::conformance::error_mode_stamp::ErrorModeStamp;
use ferro_hgvs::hgvs::variant::AllelePhase;
use ferro_hgvs::normalize::ShuffleDirection;
use ferro_hgvs::reference::multi_fasta::MultiFastaProvider;
use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer};

#[path = "common/spec_harvest.rs"]
mod spec_harvest;
use spec_harvest::classify;

/// The one row whose `spec_expected` is ferro's own rendering rather than the
/// spec's string. See the module header; the override carries the admission.
const FERRO_DERIVED_EXPECTATION_INPUTS: [&str; 1] = ["NP_060250.2:p.Gln746_Lys747ins*63"];

#[derive(Parser, Debug)]
#[command(about = "Measure spec conformance under one FERRO_PARTITION arm, over a real reference")]
struct Cli {
    /// The generated spec fixture. Only `input`, `input_prefixed` and
    /// `spec_expected` are read, all three of which are provider-independent.
    #[arg(
        long,
        default_value = "tests/fixtures/grammar/hgvs_spec_normalization.json"
    )]
    fixture: PathBuf,

    /// Prepared-reference manifest. Falls back to `$FERRO_MANIFEST`.
    #[arg(long)]
    manifest: Option<PathBuf>,

    /// Where to write this arm's measurement.
    #[arg(long)]
    output: PathBuf,

    /// Shuffle direction to normalize under.
    #[arg(long, default_value = "three-prime")]
    direction: Direction,

    /// Measure only the rows named in this file, one `input` per line.
    ///
    /// The cross-arm shortcut, and it is exact rather than a sample. Every bail
    /// on the way to `partition_block_for_rule` is arm-independent — the switch
    /// is first consulted *at* that call — so a row that cuts no block on one
    /// arm cuts none on any, and its output is identical across all four. Run
    /// the full corpus once, take the rows whose `blocks_cut` is non-zero, and
    /// the remaining arms need only those. Roughly 900 of 934 rows are then
    /// skipped, which is the difference between a half-hour arm and a brisk one.
    ///
    /// Reported in the artifact as `restricted_to`, so a census taken this way
    /// can never be mistaken for a full one.
    #[arg(long)]
    only_inputs: Option<PathBuf>,

    /// Print a progress line every N rows. Zero disables.
    ///
    /// Not a nicety: the first full run of this generator spent 29 minutes of
    /// CPU without finishing, and with no per-row output there was no way to
    /// tell a slow corpus from a stalled one. The line carries the elapsed time
    /// and the row being measured, so a pathological row (a whole-chromosome
    /// inversion, a `pter`-anchored delins) names itself.
    #[arg(long, default_value_t = 25)]
    progress_every: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
enum Direction {
    ThreePrime,
    FivePrime,
}

impl Direction {
    fn as_str(self) -> &'static str {
        match self {
            Direction::ThreePrime => "three-prime",
            Direction::FivePrime => "five-prime",
        }
    }

    fn shuffle(self) -> ShuffleDirection {
        match self {
            Direction::ThreePrime => ShuffleDirection::ThreePrime,
            Direction::FivePrime => ShuffleDirection::FivePrime,
        }
    }
}

/// The subset of a fixture row this generator reads.
#[derive(Debug, Deserialize)]
struct FixtureRow {
    input: String,
    #[serde(default)]
    input_prefixed: Option<String>,
    spec_expected: Option<String>,
}

#[derive(Debug, Deserialize)]
struct Fixture {
    rows: Vec<FixtureRow>,
}

/// One row, measured.
#[derive(Debug, Serialize)]
struct MeasuredRow {
    input: String,
    /// What ferro was run against.
    target: String,
    /// The spec's form, or `null` where the spec rejects the input.
    spec_expected: Option<String>,
    /// Ferro's rendering, or its refusal text.
    current: String,
    /// The taxonomy bucket, from the same `classify` the fixture uses.
    status: String,
    /// Blocks this row cut. Non-zero means the row reached the partitioner and
    /// so could see the arm at all.
    blocks_cut: u64,
    /// Members of the top-level cis allele, when the *normalized* form is one.
    /// `sequence_first_pass` admits a cis allele only at two or more.
    #[serde(skip_serializing_if = "Option::is_none")]
    normalized_cis_members: Option<usize>,
}

#[derive(Debug, Serialize)]
struct Measurement {
    description: &'static str,
    /// The raw `FERRO_PARTITION` value this process ran under; `live` when unset.
    arm: String,
    direction: &'static str,
    /// Manifest file name only — an absolute path would leak the local layout.
    manifest: String,
    /// The normalizer's error-handling precondition these numbers were taken
    /// under. A census compared against one from another mode is as wrong as one
    /// built from a partial pass, and nothing else in the file would say so.
    normalized_under: ErrorModeStamp,
    /// One row in, one row out. A refusal is *recorded* rather than dropped, so
    /// the only drop this pass can produce is a panic.
    ///
    /// **`dropped` is always zero here, and that is a property of the control
    /// flow rather than a measurement**: `ledger.finish()?` refuses before the
    /// artifact is written, so a run with a drop produces no file at all. Read
    /// this field as "attempted == succeeded == the row count below", which is
    /// the claim it can actually support — not as evidence that drops were
    /// looked for and found to be none. A dropped row is reported on stderr by
    /// the refusal, never here.
    completeness: CaptureCounts,
    /// `Some` when `--only-inputs` restricted the run. A census taken this way
    /// covers only the named rows and is not comparable to a full one.
    #[serde(skip_serializing_if = "Option::is_none")]
    restricted_to: Option<String>,
    /// Rows the restriction skipped. Zero on a full run.
    skipped_by_restriction: usize,
    /// Wall seconds this measurement took, for the record. Not a benchmark —
    /// these runs share a heavily loaded machine.
    elapsed_seconds: u64,
    total: usize,
    by_status: BTreeMap<String, usize>,
    /// `preserved` — ferro renders the spec's stated form. The conformance
    /// numerator, and the number a release note wants.
    conforms: usize,
    /// `conforms` minus the rows whose expectation is ferro-derived.
    conforms_excluding_ferro_derived: usize,
    ferro_derived_expectations: Vec<String>,
    /// Rows that cut at least one block, i.e. actually reached the switch.
    reached_partitioner: usize,
    /// Blocks cut over the whole run.
    blocks_cut: u64,
    /// Rows whose normalized form is a cis allele with two or more members —
    /// `sequence_first_pass`'s own admission gate, for context beside the
    /// measured `reached_partitioner`.
    multi_member_cis_after_normalize: usize,
    /// Rows that failed to evaluate, i.e. the reference did not serve them.
    unevaluable: usize,
    rows: Vec<MeasuredRow>,
}

/// Blocks this process has cut, or `0` in a release build, where the counter
/// does not exist. [`main`] refuses before any release run reaches this.
#[cfg(debug_assertions)]
fn blocks_cut() -> u64 {
    ferro_hgvs::normalize::partition_blocks_cut()
}

#[cfg(not(debug_assertions))]
fn blocks_cut() -> u64 {
    0
}

fn main() -> ExitCode {
    // A release build has no `partition_blocks_cut`, so it cannot report the
    // effective denominator. Refuse rather than emit a census whose
    // `reached_partitioner` is zero for a reason that has nothing to do with the
    // corpus — that zero is the exact misreading this generator exists to
    // prevent. `cfg!` rather than `#[cfg]` so the rest of `main` stays live code
    // in both profiles.
    if !cfg!(debug_assertions) {
        eprintln!(
            "error: this measurement needs `normalize::partition_blocks_cut`, which is \
             compiled out in release builds. Run it in a debug build."
        );
        return ExitCode::FAILURE;
    }
    let cli = Cli::parse();
    // The positive control the whole measurement rests on: a misspelled
    // `FERRO_PARTITION` silently falls back to `live` in a release build, which
    // would file `live`'s numbers under a candidate's name.
    if let Some(message) = ferro_hgvs::normalize::partition_switch_startup_error() {
        eprintln!("error: {message}");
        return ExitCode::FAILURE;
    }
    match run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// Members of a top-level **cis** allele, or `None` for any other shape.
fn cis_member_count(variant: &HgvsVariant) -> Option<usize> {
    match variant {
        HgvsVariant::Allele(allele) if allele.phase == AllelePhase::Cis && !allele.uncertain => {
            Some(allele.variants.len())
        }
        _ => None,
    }
}

fn run(cli: &Cli) -> anyhow::Result<()> {
    let manifest = cli
        .manifest
        .clone()
        .or_else(|| std::env::var("FERRO_MANIFEST").ok().map(PathBuf::from))
        .ok_or_else(|| {
            anyhow::anyhow!(
                "no manifest (pass --manifest or set FERRO_MANIFEST). Without a real reference \
                 this measurement is the vacuous one it exists to replace."
            )
        })?;

    let text = std::fs::read_to_string(&cli.fixture).map_err(|e| {
        anyhow::anyhow!(
            "read {}: {e}. Generate it first: \
             `cargo run --features dev --bin generate_spec_fixture`",
            cli.fixture.display()
        )
    })?;
    let fixture: Fixture = serde_json::from_str(&text)
        .map_err(|e| anyhow::anyhow!("parse {}: {e}", cli.fixture.display()))?;
    if fixture.rows.is_empty() {
        anyhow::bail!("the fixture holds no rows; there is nothing to measure");
    }

    let provider = MultiFastaProvider::from_manifest(&manifest)
        .map_err(|e| anyhow::anyhow!("load manifest {}: {e}", manifest.display()))?;
    // The same configuration `generate_spec_fixture` measures under — which is
    // `ErrorConfig::lenient()`, not strict, whatever "default" suggests — with
    // the direction as the only difference. Kept identical on purpose: a census
    // taken under another mode is not comparable to the fixture's.
    let config = NormalizeConfig::default().with_direction(cli.direction.shuffle());
    let normalized_under = ErrorModeStamp::of(&config.error_config);
    let normalizer = Normalizer::with_config(provider, config);

    let arm = std::env::var("FERRO_PARTITION").unwrap_or_else(|_| "live".to_string());
    let ferro_derived: Vec<String> = FERRO_DERIVED_EXPECTATION_INPUTS
        .iter()
        .map(|s| (*s).to_string())
        .collect();

    // Every fixture row must reach the artifact. A refusal is a *result* here —
    // `parse-error` and `needs-reference` are buckets in the taxonomy — so this
    // pass has no legitimate drop, and the ledger is what makes that visible
    // rather than assumed. The one way a row could vanish is a panic, which the
    // ledger records as a shortfall at `finish` instead of leaving the artifact
    // one row shorter than the corpus with nothing to say so.
    let only: Option<std::collections::BTreeSet<String>> = match &cli.only_inputs {
        None => None,
        Some(path) => {
            let text = std::fs::read_to_string(path)
                .map_err(|e| anyhow::anyhow!("read {}: {e}", path.display()))?;
            let set: std::collections::BTreeSet<String> = text
                .lines()
                .map(str::trim)
                .filter(|line| !line.is_empty())
                .map(str::to_string)
                .collect();
            if set.is_empty() {
                anyhow::bail!(
                    "{} names no rows; a restriction to nothing would report an empty census \
                     that reads like a measurement",
                    path.display()
                );
            }
            Some(set)
        }
    };

    let started = std::time::Instant::now();
    let mut skipped_by_restriction = 0usize;
    let mut ledger = CaptureLedger::new("spec fixture rows");
    let mut rows = Vec::with_capacity(fixture.rows.len());
    let mut by_status: BTreeMap<String, usize> = BTreeMap::new();
    // Rows this pass will actually measure — the progress line's denominator.
    // Counting the *fixture* instead would report `[7/934]` on a 35-row
    // restricted run, and, worse, would step the counter on rows the restriction
    // skips: with ~900 of 934 skipped, `index % progress_every` then fires on
    // roughly one measured row in 25 by coincidence of position, so the
    // recommended `--only-inputs` path — the one this flag exists to make
    // survivable — could print no progress at all.
    let to_measure = match &only {
        None => fixture.rows.len(),
        Some(set) => fixture
            .rows
            .iter()
            .filter(|row| set.contains(&row.input))
            .count(),
    };
    let mut measured = 0usize;
    for row in &fixture.rows {
        if only.as_ref().is_some_and(|set| !set.contains(&row.input)) {
            skipped_by_restriction += 1;
            continue;
        }
        if cli.progress_every > 0 && measured.is_multiple_of(cli.progress_every) {
            eprintln!(
                "  [{}/{}] {:>5}s  {}",
                measured,
                to_measure,
                started.elapsed().as_secs(),
                row.input
            );
        }
        measured += 1;
        let target = row
            .input_prefixed
            .clone()
            .unwrap_or_else(|| row.input.clone());
        let before = blocks_cut();
        // Contained, because a real reference opens code paths the hermetic run
        // never took and one panic would cost the whole 934-row run. A panic is
        // *recorded* as its own status rather than swallowed — a skip that reads
        // as a pass is the failure mode this whole generator exists to remove —
        // and `panicked` rows are excluded from `conforms` by construction,
        // since `classify` never returns that string.
        let outcome =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| match parse_hgvs(&target) {
                Err(e) => (format!("parse error: {e}"), false, false, Vec::new(), None),
                Ok(variant) => match normalizer.normalize_with_diagnostics(&variant) {
                    Err(e) => (
                        format!("normalize error: {e}"),
                        true,
                        false,
                        Vec::new(),
                        None,
                    ),
                    Ok(n) => (
                        format!("{}", n.result),
                        true,
                        true,
                        n.result.spec_equivalent_renderings(),
                        cis_member_count(&n.result),
                    ),
                },
            }));
        let cut = blocks_cut() - before;

        let (current, status, members) = match outcome {
            Ok((current, parse_ok, normalize_ok, equivalents, members)) => {
                let status = classify::classify(
                    parse_ok,
                    normalize_ok,
                    &current,
                    &equivalents,
                    row.spec_expected.as_deref(),
                );
                ledger.record_success();
                (current, status, members)
            }
            Err(_) => {
                // A row this pass could not measure. Counted as a drop, so
                // `finish` refuses the run rather than letting the artifact
                // claim a complete census it does not have.
                ledger.record_drop(row.input.clone(), "normalization panicked");
                ("panicked".to_string(), "panicked", None)
            }
        };
        *by_status.entry(status.to_string()).or_default() += 1;
        rows.push(MeasuredRow {
            input: row.input.clone(),
            target,
            spec_expected: row.spec_expected.clone(),
            current,
            status: status.to_string(),
            blocks_cut: cut,
            normalized_cis_members: members,
        });
    }
    // Panics are the only drop this pass can produce, and one invalidates the
    // census, so the ceiling is zero and `finish` is the right closer. If a
    // future corpus legitimately carries a panicking row, waive it by name with
    // `finish_with(Allowance::at_most(n, "why"))` — never by widening this.
    let completeness = ledger.finish()?;

    let conforms = by_status.get("preserved").copied().unwrap_or(0);
    let conforms_excluding_ferro_derived = rows
        .iter()
        .filter(|r| r.status == "preserved" && !ferro_derived.contains(&r.input))
        .count();
    let measurement = Measurement {
        description:
            "Per-`FERRO_PARTITION`-arm conformance of the HGVS spec corpus against the SPEC's \
             stated form, over a prepared reference. `spec_expected` does not move with the arm, \
             so this measures ferro against the spec once per arm — never the arms against each \
             other. A zero delta is evidence only where `reached_partitioner` is non-zero.",
        arm,
        direction: cli.direction.as_str(),
        manifest: manifest
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| "manifest.json".to_string()),
        normalized_under,
        completeness,
        restricted_to: cli.only_inputs.as_ref().map(|p| p.display().to_string()),
        skipped_by_restriction,
        elapsed_seconds: started.elapsed().as_secs(),
        total: rows.len(),
        conforms,
        conforms_excluding_ferro_derived,
        ferro_derived_expectations: ferro_derived,
        reached_partitioner: rows.iter().filter(|r| r.blocks_cut > 0).count(),
        blocks_cut: rows.iter().map(|r| r.blocks_cut).sum(),
        multi_member_cis_after_normalize: rows
            .iter()
            .filter(|r| r.normalized_cis_members.is_some_and(|n| n >= 2))
            .count(),
        unevaluable: rows
            .iter()
            .filter(|r| r.current.starts_with("normalize error:"))
            .count(),
        by_status,
        rows,
    };

    // A zero here is a fact about the corpus, not about the arms, and it is the
    // one misreading this generator exists to prevent. Say so on the run rather
    // than leaving it to be inferred from a field in the artifact.
    if measurement.reached_partitioner == 0 {
        eprintln!(
            "WARNING: no row cut a block, so this run measured NOTHING about the \
             FERRO_PARTITION arm. Any cross-arm delta computed from it is structural. \
             Check that the manifest actually serves this corpus's accessions."
        );
    }

    let mut out = serde_json::to_string_pretty(&measurement)?;
    out.push('\n');
    std::fs::write(&cli.output, out)?;
    eprintln!(
        "arm={} direction={} conforms={}/{} reached_partitioner={} blocks_cut={}",
        measurement.arm,
        measurement.direction,
        measurement.conforms,
        measurement.total,
        measurement.reached_partitioner,
        measurement.blocks_cut
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The control this whole measurement is a change *from*: with an empty
    /// provider — exactly what `generate_spec_fixture` uses — a two-member cis
    /// allele cuts **no** block, so `reached_partitioner` would be zero for
    /// every row and the run would be as vacuous as the fixture it replaces.
    ///
    /// The positive half cannot be asserted here: moving the counter needs real
    /// reference bases, and CI has no prepared manifest. It is asserted at
    /// report time instead, by the generator refusing to let a zero
    /// `reached_partitioner` pass unremarked — see `run`.
    ///
    /// Deliberately not a fixed value — the count is process-wide and monotone,
    /// so only its *movement* is a stable property.
    #[test]
    fn an_empty_provider_cuts_no_block() {
        use ferro_hgvs::reference::mock::MockProvider;

        // An empty provider cannot serve a block, so nothing is cut: this is the
        // control that gives the assertion below something to be a change from,
        // and it is the exact condition that makes the spec fixture blind.
        let hermetic = Normalizer::new(MockProvider::new());
        let variant = parse_hgvs("NM_004006.2:c.[79G>T;80C>T]").expect("parse");
        let before = blocks_cut();
        let _ = hermetic.normalize_with_diagnostics(&variant);
        assert_eq!(
            blocks_cut(),
            before,
            "an empty provider cut a block; the counter is not measuring what this \
             generator reads it for"
        );
    }

    /// `cis_member_count` must answer only for the shape `sequence_first_pass`
    /// admits, so the context figure beside `reached_partitioner` cannot be read
    /// as a wider population than the gate.
    #[test]
    fn only_a_certain_cis_allele_is_counted_as_multi_member() {
        let cis = parse_hgvs("NM_004006.2:c.[79G>T;80C>T]").expect("parse cis");
        assert_eq!(cis_member_count(&cis), Some(2));

        let single = parse_hgvs("NM_004006.2:c.79G>T").expect("parse single");
        assert_eq!(cis_member_count(&single), None);

        let trans = parse_hgvs("NM_004006.2:c.[79G>T];[80C>T]").expect("parse trans");
        assert_eq!(
            cis_member_count(&trans),
            None,
            "a trans allele is refused by `sequence_first_pass` and must not be counted"
        );
    }
}
