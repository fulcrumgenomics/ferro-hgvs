//! The `FERRO_PARTITION` refusal must be delivered by every entry point.
//!
//! `partition_rule()` used to panic on a value naming no arm, in every build,
//! from inside `canonicalize_from_sequence` — infallible, under every public
//! normalization entry point, and reached across the FFI boundary by the Python
//! bindings. A development-only bake-off switch must not be able to abort a
//! process that merely happens to have the variable set, so a release build of
//! the library now falls safe to the shipped `live` rule and retains the refusal
//! in `partition_switch_startup_error()` for a caller to report.
//!
//! That trade only holds if the callers actually report it. Falling safe without
//! reporting is the *original* defect this switch was given a refusal to prevent:
//! a bake-off served the shipped rule under a candidate's name, coming back with
//! a clean empty diff that reads as "the candidate changes nothing". So the
//! wiring is the contract, and these two tests pin both halves of it:
//!
//! - the behaviour, end to end, through the `ferro` binary — a misspelling exits
//!   non-zero with the message, a real arm still runs, and `--help` keeps working;
//! - the completeness, over `Cargo.toml`'s own `[[bin]]` table — so a binary
//!   added later cannot skip the call and leave `README.md`'s "every binary in
//!   this repository reports it" quietly false.
//!
//! The second is a source scan, so it proves the call is *written*, not that it
//! is reached. The first is what proves it is reached, and it is written against
//! the binary a user actually runs.

use std::collections::BTreeSet;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

/// A value no build has an arm for. Deliberately the near-miss `README.md` uses,
/// since a plausible misspelling is the input the refusal exists for.
const MISSPELLED: &str = "canonicl";

/// The repository root, from this test target's manifest directory.
fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// A deletion inside the `CCCC` homopolymer at c.16-19 of the mock `NM_000088.3`,
/// so the normalizer 3'-shifts it and the positive control below exercises the
/// partitioner rather than only the parser. Within the mock CDS (60 bases), which
/// the more obvious `c.459A>G` is not — that one is rejected as past the CDS end
/// under the default strict mode.
const SHIFTED: &str = "NM_000088.3:c.16del";

/// Run `ferro normalize` over one variant with `FERRO_PARTITION` set to `value`.
///
/// The mock provider serves `NM_000088.3`, so no `--reference` and no prepared
/// reference directory is needed.
fn normalize_one(value: &str) -> std::process::Output {
    let mut input = tempfile::Builder::new()
        .suffix(".txt")
        .tempfile()
        .expect("create input file");
    writeln!(input, "{SHIFTED}").expect("write input");
    input.flush().expect("flush input");

    Command::new(env!("CARGO_BIN_EXE_ferro"))
        .arg("normalize")
        .arg("-i")
        .arg(input.path())
        .env("FERRO_PARTITION", value)
        .output()
        .expect("run ferro normalize")
}

/// A misspelled arm name stops the CLI before it reads input, and a real one does
/// not.
///
/// Both directions are asserted in one test on purpose: a binary that refused
/// *every* value would pass the first assertion while making the switch useless,
/// and that failure would look identical in a log to the fix working.
#[test]
fn the_ferro_cli_refuses_a_misspelled_partition_switch_and_still_serves_a_real_arm() {
    let refused = normalize_one(MISSPELLED);
    assert!(
        !refused.status.success(),
        "a value naming no arm must fail the process; stdout was:\n{}",
        String::from_utf8_lossy(&refused.stdout)
    );
    let stderr = String::from_utf8_lossy(&refused.stderr);
    assert!(
        stderr.contains(MISSPELLED) && stderr.contains("is not a partitioner this build has"),
        "the refusal must name the value given and say what it is; stderr was:\n{stderr}"
    );
    assert!(
        stderr.contains("canonical-coalesced"),
        "the refusal must name THIS build's arms, so a value that exists on \
         another branch is reported as absent here; stderr was:\n{stderr}"
    );
    assert!(
        refused.stdout.is_empty(),
        "the refusal must precede any output that could be mistaken for a \
         measurement; stdout was:\n{}",
        String::from_utf8_lossy(&refused.stdout)
    );

    let accepted = normalize_one("canonical");
    assert!(
        accepted.status.success(),
        "a value that IS an arm must still run; stderr was:\n{}",
        String::from_utf8_lossy(&accepted.stderr)
    );
}

/// `--help` still answers under a stale `FERRO_PARTITION`.
///
/// Reading the help is how somebody finds out what the arms are, so refusing it
/// would leave them holding a rejected value and no way to look up the right one.
/// This is why the check sits after `Cli::parse()` rather than before it — clap
/// answers `--help` by exiting from inside `parse`, and nothing has been read or
/// normalized at that point either way.
#[test]
fn help_still_works_when_the_partition_switch_names_no_arm() {
    let out = Command::new(env!("CARGO_BIN_EXE_ferro"))
        .arg("--help")
        .env("FERRO_PARTITION", MISSPELLED)
        .output()
        .expect("run ferro --help");
    assert!(
        out.status.success(),
        "`--help` must not be refused; stderr was:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        !out.stdout.is_empty(),
        "`--help` must still print the help text"
    );
}

/// Every `[[bin]]` target's source calls `partition_switch_startup_error`.
///
/// Parsed as TOML rather than substring-matched, for the reason
/// `generator_completeness.rs` gives: a `path` inside a comment, or under a
/// `[[example]]` table, must not count.
///
/// The denominator is asserted non-zero, so "0 of 0 binaries are missing it"
/// cannot pass as a result — the failure mode this repository has been bitten by
/// often enough to make it a house rule.
#[test]
fn every_binary_target_reports_a_refused_partition_switch() {
    let manifest_text =
        std::fs::read_to_string(repo_root().join("Cargo.toml")).expect("read Cargo.toml");
    let manifest: toml::Value = manifest_text.parse().expect("parse Cargo.toml as TOML");

    let mut sources = BTreeSet::new();
    for entry in manifest
        .get("bin")
        .and_then(toml::Value::as_array)
        .expect("Cargo.toml declares at least one [[bin]] target")
    {
        let name = entry
            .get("name")
            .and_then(toml::Value::as_str)
            .expect("every [[bin]] has a name");
        let path = entry
            .get("path")
            .and_then(toml::Value::as_str)
            .unwrap_or_else(|| panic!("[[bin]] {name} declares no path"));
        sources.insert((name.to_string(), path.to_string()));
    }
    assert!(
        sources.len() >= 5,
        "expected the five known binaries (ferro, ferro-web, ferro-benchmark and \
         the two spec generators); found {sources:?}"
    );

    // The bake-off harness is not a `[[bin]]`, but it is the entry point a
    // blast-radius measurement runs through, so it carries the same obligation:
    // a refused value there yields an empty diff AND a silent decline census.
    let mut expected: Vec<(String, String)> = sources.into_iter().collect();
    expected.push((
        "dump_normalized_corpus".to_string(),
        "examples/dump_normalized_corpus.rs".to_string(),
    ));

    let missing: Vec<&str> = expected
        .iter()
        .filter(|(_, path)| {
            let text = std::fs::read_to_string(repo_root().join(path))
                .unwrap_or_else(|e| panic!("read {path}: {e}"));
            !text.contains("partition_switch_startup_error()")
        })
        .map(|(name, _)| name.as_str())
        .collect();

    assert!(
        missing.is_empty(),
        "these entry points normalize without reporting a refused \
         FERRO_PARTITION, so they would serve the shipped rule under a \
         candidate's name in a release build: {missing:?}. Call \
         `ferro_hgvs::normalize::partition_switch_startup_error()` after \
         argument parsing and fail on `Some`."
    );
}
