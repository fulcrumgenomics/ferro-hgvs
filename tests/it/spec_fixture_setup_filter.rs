//! The spec-fixture setup script must be bound to its consumers — all of them,
//! and only them.
//!
//! # The two defects this sits between
//!
//! `.config/nextest.toml` wires `scripts/pregenerate-spec-fixtures.sh` as a
//! nextest setup script so the two gitignored spec fixtures are generated once,
//! serially, before any test process starts (#1608). Which tests that binding is
//! attached to is the whole question, and it can be got wrong in two opposite
//! directions.
//!
//! **Too narrow** is the failure the original `filter = 'all()'` was written to
//! avoid: a module that reads a fixture but is not named here runs with the
//! artifact absent, falls back to `common::fixture_gen`'s nested `cargo run`,
//! and several such processes race the cargo build lock — which is the flake
//! #1608 reports.
//!
//! **Too wide** is what `all()` actually cost, and it is the more expensive of
//! the two because it fails jobs that never asked for a fixture. `soak` and
//! `sweeps` check out with no submodule and run pre-built binaries from a
//! nextest archive, so the spec checkout they would have to harvest is empty and
//! the generator cannot succeed there at any price:
//!
//! ```text
//! error: no HGVS strings harvested from assets/hgvs-nomenclature — the spec checkout is empty.
//!   SETUP FAIL [  88.702s] pregenerate-spec-fixtures
//!   Cancelling due to setup script failure
//!      Summary [  88.703s] 0/54 tests run: 0 passed, 80 skipped
//! ```
//!
//! Measured on run 31931487583: 8 soak shards and all 3 sweeps ran **zero**
//! tests, having first spent 68.7-88.7s compiling in jobs designed to do no
//! building.
//!
//! # What is checked, and what it is worth
//!
//! Both directions, from **derived** sets rather than from a roster kept here —
//! a hardcoded list of consumers would rot in exactly the way the `all()`
//! argument predicted:
//!
//! 1. the setup script's filter names precisely the `tests/it` modules that call
//!    `ensure_spec_fixture()` or `ensure_spec_enumeration()`; and
//! 2. no CI job that runs from a nextest archive **without** downloading the
//!    `spec-fixtures` artifact selects any of those modules — i.e. the filter
//!    cannot fire where the fixtures are neither present nor buildable.
//!
//! The workflow half is **parsed**, for the reason `generated_fixture_ci_wiring`
//! gives: which step sits in which job is structural. The source half is a
//! **token scan** and is a floor, not a proof — it finds the call spelled with
//! its parentheses, so a consumer reaching the fixture through a new helper, or
//! aliasing the call, is invisible to it. What it buys is that a consumer
//! written the way all four current ones are written cannot be added in silence,
//! and the question gets asked.

use std::collections::BTreeSet;
use std::path::PathBuf;

use serde_yaml::Value;

/// This file. Excluded from the scan for the same reason
/// `generated_fixture_ci_wiring.rs` excludes itself: it carries [`CALLS`]'s
/// matcher literals, so it would match itself and demand its own name in a
/// filter it only describes.
const SELF: &str = "spec_fixture_setup_filter.rs";

/// The calls that mark a module as reading a generated spec fixture, spelled
/// with the parenthesis so a prose mention in a doc comment does not count.
/// Both helpers live in `tests/it/common/`, which the scan skips.
const CALLS: [&str; 2] = ["ensure_spec_fixture()", "ensure_spec_enumeration()"];

/// A floor on the derived consumer set, so a matcher that has silently stopped
/// matching cannot pass as "nothing to check". Four modules consume today:
/// `hgvs_spec_normalization_tests`, `idempotency_tests`,
/// `normalize_axis_preserving` and `spec_enumeration_tests`.
const MINIMUM_CONSUMERS: usize = 4;

/// The setup script the filter under test is attached to.
const SETUP_SCRIPT: &str = "pregenerate-spec-fixtures";

/// The artifact carrying the generated spec fixtures.
const FIXTURE_ARTIFACT: &str = "spec-fixtures";

/// The workflows that run tests from a nextest archive.
///
/// Both are scanned rather than only `ci.yml`. `nightly-soak.yml` runs the same
/// `test(proptest)` selection, off the same `nextest-archive-soak` artifact,
/// with the same uninitialised submodule — so a consumer swept into that
/// selection breaks there identically. It is non-gating, which makes it a more
/// likely place for such a break to sit unnoticed, not a less likely one.
///
/// This list is NAMED, which is a floor: a third workflow running the archive
/// would be missed until somebody added it here. That is what
/// [`the_named_archive_workflows_are_the_ones_that_run_archives`] re-derives, so
/// the constant cannot silently fall behind the `.github/workflows` directory.
const ARCHIVE_WORKFLOWS: [&str; 2] = ["ci.yml", "nightly-soak.yml"];

/// The wrapper every archive-running job invokes, and so the marker that
/// identifies such a workflow from its source text.
const ARCHIVE_RUNNER: &str = "run_nextest_archive.sh";

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Every `tests/it` module that reads a generated spec fixture, by file stem —
/// which is also the leading segment of every test name nextest gives it, and so
/// the thing a `test(...)` filter matches.
fn fixture_consuming_modules() -> BTreeSet<String> {
    let dir = repo_root().join("tests/it");
    let entries = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("read {}: {e}", dir.display()))
        .filter_map(Result::ok);

    let mut modules = BTreeSet::new();
    for entry in entries {
        let path = entry.path();
        if path.extension().is_none_or(|extension| extension != "rs") {
            continue;
        }
        let name = path
            .file_name()
            .expect("a file has a name")
            .to_string_lossy()
            .into_owned();
        if name == SELF {
            continue;
        }
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
        if CALLS.iter().any(|call| text.contains(call)) {
            modules.insert(name.trim_end_matches(".rs").to_string());
        }
    }

    assert!(
        modules.len() >= MINIMUM_CONSUMERS,
        "found only {} spec-fixture consumer(s) ({:?}); at least {MINIMUM_CONSUMERS} read one \
         today, so the scan for {CALLS:?} has gone blind rather than the repository having shrunk",
        modules.len(),
        modules
    );
    modules
}

/// The setup script's filterset, read from `.config/nextest.toml`.
fn setup_script_filter() -> String {
    let path = repo_root().join(".config/nextest.toml");
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let config: toml::Value =
        toml::from_str(&text).expect(".config/nextest.toml is not valid TOML");

    // The script must be DEFINED as well as bound; a binding naming a script
    // that does not exist is a config error nextest reports, but a definition
    // with no binding is silent and would leave the fixtures ungenerated.
    config
        .get("scripts")
        .and_then(|scripts| scripts.get("setup"))
        .and_then(|setup| setup.get(SETUP_SCRIPT))
        .and_then(|script| script.get("command"))
        .unwrap_or_else(|| {
            panic!("`.config/nextest.toml` defines no `[scripts.setup.{SETUP_SCRIPT}]` command")
        });

    let bindings = config
        .get("profile")
        .and_then(|profile| profile.get("default"))
        .and_then(|default| default.get("scripts"))
        .and_then(toml::Value::as_array)
        .unwrap_or_else(|| panic!("`.config/nextest.toml` has no `[[profile.default.scripts]]`"));

    let binding = bindings
        .iter()
        .find(|binding| binding.get("setup").and_then(toml::Value::as_str) == Some(SETUP_SCRIPT))
        .unwrap_or_else(|| panic!("no `[[profile.default.scripts]]` entry binds `{SETUP_SCRIPT}`"));

    binding
        .get("filter")
        .and_then(toml::Value::as_str)
        .unwrap_or_else(|| {
            panic!(
                "the `{SETUP_SCRIPT}` binding has no `filter`; without one nextest runs it for \
                 every test, which is what broke the archive jobs"
            )
        })
        .to_string()
}

/// The arguments of every `test(...)` predicate in a filterset.
///
/// A token scan over a small, fixed grammar rather than a filterset parser: the
/// selections in play here are all `test(a) + test(b)` unions, and the failure
/// mode of a scan that stops matching is covered by the emptiness assertions at
/// both call sites.
fn test_predicate_arguments(filterset: &str) -> BTreeSet<String> {
    let mut arguments = BTreeSet::new();
    let mut rest = filterset;
    while let Some(at) = rest.find("test(") {
        let after = &rest[at + "test(".len()..];
        match after.find(')') {
            Some(end) => {
                arguments.insert(after[..end].trim().to_string());
                rest = &after[end..];
            }
            None => break,
        }
    }
    arguments
}

fn workflow_dir() -> PathBuf {
    repo_root().join(".github/workflows")
}

fn workflow(file: &str) -> Value {
    let path = workflow_dir().join(file);
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    serde_yaml::from_str(&text).unwrap_or_else(|e| panic!("{file} is not valid YAML: {e}"))
}

fn steps(job: &Value) -> &[Value] {
    job.get("steps")
        .and_then(Value::as_sequence)
        .map(Vec::as_slice)
        .unwrap_or(&[])
}

/// The workflow-level `env:` block, whose values the job selections interpolate.
fn workflow_env(workflow: &Value) -> Vec<(String, String)> {
    workflow
        .get("env")
        .and_then(Value::as_mapping)
        .map(|env| {
            env.iter()
                .filter_map(|(key, value)| {
                    Some((key.as_str()?.to_string(), value.as_str()?.to_string()))
                })
                .collect()
        })
        .unwrap_or_default()
}

/// A job's `run:` text with `$VAR` / `${VAR}` expanded from the workflow `env:`.
///
/// The selections this file cares about are written as `-E "$SWEEP_FILTER + …"`,
/// so an unexpanded scan would read the variable's NAME and find no module in
/// it — passing vacuously on exactly the job the check exists for.
fn expanded_run_text(job: &Value, env: &[(String, String)]) -> String {
    let mut text: String = steps(job)
        .iter()
        .filter_map(|step| step.get("run").and_then(Value::as_str))
        .collect::<Vec<_>>()
        .join("\n");
    for (key, value) in env {
        text = text
            .replace(&format!("${{{key}}}"), value)
            .replace(&format!("${key}"), value);
    }
    text
}

/// The artifact names a job downloads.
fn downloaded_artifacts(job: &Value) -> BTreeSet<String> {
    steps(job)
        .iter()
        .filter(|step| {
            step.get("uses")
                .and_then(Value::as_str)
                .unwrap_or_default()
                .starts_with("actions/download-artifact")
        })
        .filter_map(|step| {
            step.get("with")
                .and_then(|with| with.get("name"))
                .and_then(Value::as_str)
                .map(str::to_string)
        })
        .collect()
}

/// The setup script fires for every module that reads a spec fixture, and for no
/// other module.
///
/// Equality, not containment, and both halves have a bug behind them: a missing
/// name reintroduces #1608's build-lock race, and an extra one is how `all()`
/// emptied `soak` and `sweeps`.
#[test]
fn the_setup_script_is_bound_to_exactly_the_fixture_consumers() {
    let filterset = setup_script_filter();
    assert!(
        !filterset.contains("all()"),
        "the `{SETUP_SCRIPT}` binding is `{filterset}`. `all()` fires the script in jobs that \
         run pre-built binaries with no spec submodule, where the generator cannot succeed — \
         it emptied 8 soak shards and 3 sweeps on run 31931487583."
    );

    let selected = test_predicate_arguments(&filterset);
    assert!(
        !selected.is_empty(),
        "no `test(...)` predicate found in the `{SETUP_SCRIPT}` filter `{filterset}`; the scan \
         has stopped matching and this guard is checking nothing"
    );

    let consumers = fixture_consuming_modules();
    assert_eq!(
        selected,
        consumers,
        "the `{SETUP_SCRIPT}` filter and the set of modules calling {CALLS:?} disagree.\n  \
         named in the filter but reads no fixture: {:?}\n  \
         reads a fixture but is not in the filter: {:?}\n\
         Add or remove the `test(<module>)` term in `.config/nextest.toml`. A module missing \
         from the filter falls back to the nested `cargo run` this script exists to remove \
         (#1608); a module present but not a reader fires the script in jobs that may be \
         unable to run it.",
        selected.difference(&consumers).collect::<Vec<_>>(),
        consumers.difference(&selected).collect::<Vec<_>>(),
    );
}

/// No archive job lacking the `spec-fixtures` artifact selects a fixture
/// consumer.
///
/// The sibling of the check above, closing the half it cannot see. That one
/// keeps the filter equal to the consumer set; this one keeps the consumer set
/// out of the selections where the script could not do its job — a new
/// `*_proptest` module reading the spec fixture would be picked up by `soak`'s
/// `test(proptest)` without anybody naming it there.
#[test]
fn no_fixtureless_archive_job_selects_a_spec_fixture_consumer() {
    let consumers = fixture_consuming_modules();
    let mut checked = 0usize;

    for file in ARCHIVE_WORKFLOWS {
        let parsed = workflow(file);
        let env = workflow_env(&parsed);
        let jobs = parsed
            .get("jobs")
            .and_then(Value::as_mapping)
            .unwrap_or_else(|| panic!("{file} has no jobs mapping"));

        for (name, job) in jobs {
            let job_name = name.as_str().unwrap_or("<unnamed>");
            let artifacts = downloaded_artifacts(job);
            // Only jobs that run from an archive are in scope: a job building
            // from source can always generate a missing fixture.
            if !artifacts.iter().any(|a| a.starts_with("nextest-archive")) {
                continue;
            }
            if artifacts.contains(FIXTURE_ARTIFACT) {
                continue;
            }
            checked += 1;

            let run = expanded_run_text(job, &env);
            let selection: BTreeSet<String> = test_predicate_arguments(&run);
            assert!(
                !selection.is_empty(),
                "job `{job_name}` in {file} runs from a nextest archive but no `test(...)` \
                 predicate was found in its expanded `run:` text; the scan has gone blind \
                 rather than the job having no selection"
            );

            for consumer in &consumers {
                for selected in &selection {
                    // `test(x)` is a SUBSTRING match on the test name, so a
                    // selection term contained in a consumer's module name
                    // selects it just as surely as an exact one.
                    assert!(
                        !consumer.contains(selected.as_str()),
                        "job `{job_name}` in {file} selects `test({selected})`, which matches \
                         the spec-fixture consumer `{consumer}`, but the job does not download \
                         the `{FIXTURE_ARTIFACT}` artifact. The setup script would fire there \
                         with the fixtures absent and no spec submodule to build them from — \
                         the shape that ran 0 tests across 8 soak shards and 3 sweeps on run \
                         31931487583. Either download `{FIXTURE_ARTIFACT}` in that job, or keep \
                         the module out of its selection."
                    );
                }
            }
        }
    }

    assert!(
        checked >= 3,
        "found {checked} archive job(s) without the `{FIXTURE_ARTIFACT}` artifact; `ci.yml`'s \
         `soak` and `sweeps` and `nightly-soak.yml`'s shard job are all such jobs, so the \
         derivation has gone blind"
    );
}

/// [`ARCHIVE_WORKFLOWS`] names every workflow that runs tests from an archive.
///
/// The check above is only as wide as that constant, and a named list is exactly
/// the thing this whole file exists to distrust — the `all()` binding it
/// replaces was chosen to avoid maintaining one. So the list is re-derived from
/// the workflow directory rather than trusted.
#[test]
fn the_named_archive_workflows_are_the_ones_that_run_archives() {
    let dir = workflow_dir();
    let entries = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("read {}: {e}", dir.display()))
        .filter_map(Result::ok);

    let mut found = BTreeSet::new();
    for entry in entries {
        let path = entry.path();
        if path.extension().is_none_or(|e| e != "yml" && e != "yaml") {
            continue;
        }
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
        if text.contains(ARCHIVE_RUNNER) {
            found.insert(
                path.file_name()
                    .expect("a file has a name")
                    .to_string_lossy()
                    .into_owned(),
            );
        }
    }

    let named: BTreeSet<String> = ARCHIVE_WORKFLOWS.iter().map(|w| w.to_string()).collect();
    assert_eq!(
        found,
        named,
        "ARCHIVE_WORKFLOWS and the workflows invoking `{ARCHIVE_RUNNER}` disagree.\n  \
         named but does not run an archive: {:?}\n  \
         runs an archive but is not named: {:?}\n\
         A workflow running archived binaries with no spec submodule is exactly where an \
         over-wide setup-script filter empties a job, so it must be in scope for the check \
         above.",
        named.difference(&found).collect::<Vec<_>>(),
        found.difference(&named).collect::<Vec<_>>(),
    );
}
