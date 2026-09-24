//! Every repository script a workflow runs must exist, and a step that fetches
//! test data must be able to fail its job.
//!
//! # The failure this exists for
//!
//! `external-validation.yml` ran `scripts/fetch_mutalyzer_normalized.py`,
//! `fetch_lovd_hgvs.py`, `fetch_mavedb_hgvs.py` and `fetch_civic_hgvs.py` every
//! week. None of the four has ever existed in this repository. Each step died
//! with `can't open file`, `continue-on-error: true` rewrote the failure to
//! success, and the tests downstream of the Mutalyzer fetch passed having
//! checked nothing (#2257). `actionlint` cannot see this: a `run:` block naming a
//! missing file is valid YAML and valid shell.
//!
//! So two checks, both over the whole of `.github/workflows/`:
//!
//! - every `scripts/…` or `.github/scripts/…` path a `run:` block names exists;
//! - no step whose `run:` block runs a `scripts/fetch…` script carries
//!   `continue-on-error: true`. A fetch that cannot fail its step is how the
//!   missing scripts went unnoticed, and the tests that read its output then
//!   test nothing, or test committed data in place of fetched data.

use std::path::PathBuf;
use std::sync::LazyLock;

use regex::Regex;
use serde_yaml::Value;

/// Where the workflows live. The set is derived from this directory rather than
/// named, so a new workflow is covered without editing this file.
const WORKFLOW_DIR: &str = ".github/workflows";

/// A floor on the number of script references found, so a traversal or pattern
/// that has gone blind fails instead of reporting a clean tree. Set well below
/// today's count: its subject is the scanner, not the repository's size.
const MINIMUM_SCRIPT_REFERENCES: usize = 8;

/// A floor on the number of fetch steps found, for the same reason: the
/// `continue-on-error` check passes vacuously over an empty set.
const MINIMUM_FETCH_STEPS: usize = 3;

/// A repository script path as it appears in a `run:` block: `scripts/x.py`,
/// `./scripts/x.sh`, `.github/scripts/x.sh`.
static SCRIPT_PATH: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"(?:^|[\s'\x22=(])(?:\./)?((?:\.github/)?scripts/[A-Za-z0-9_./-]+\.(?:py|sh))\b")
        .expect("the script-path pattern compiles")
});

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// One workflow step that has a `run:` block.
struct RunStep {
    locus: String,
    script: String,
    continue_on_error: bool,
}

/// Every step with a `run:` block, across every workflow, sorted by locus.
fn run_steps() -> Vec<RunStep> {
    let dir = repo_root().join(WORKFLOW_DIR);
    let mut steps = Vec::new();
    for entry in std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("read {}: {e}", dir.display()))
        .filter_map(Result::ok)
    {
        let path = entry.path();
        if path
            .extension()
            .is_none_or(|extension| extension != "yml" && extension != "yaml")
        {
            continue;
        }
        let name = path
            .file_name()
            .expect("a file has a name")
            .to_string_lossy()
            .into_owned();
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
        let document: Value = serde_yaml::from_str(&text)
            .unwrap_or_else(|e| panic!("{WORKFLOW_DIR}/{name} is not valid YAML: {e}"));
        let Some(jobs) = document.get("jobs").and_then(Value::as_mapping) else {
            continue;
        };
        for (job_name, job) in jobs {
            let job_name = job_name.as_str().unwrap_or("?");
            let Some(job_steps) = job.get("steps").and_then(Value::as_sequence) else {
                continue;
            };
            for (index, step) in job_steps.iter().enumerate() {
                let Some(script) = step.get("run").and_then(Value::as_str) else {
                    continue;
                };
                steps.push(RunStep {
                    locus: format!("{WORKFLOW_DIR}/{name} job `{job_name}` step {}", index + 1),
                    script: script.to_string(),
                    continue_on_error: masks_failure(step.get("continue-on-error")),
                });
            }
        }
    }
    steps.sort_by(|left, right| left.locus.cmp(&right.locus));
    steps
}

/// Whether a step's `continue-on-error` value can let a failed step pass. Only
/// an absent value or a literal `false` cannot: GitHub evaluates an expression
/// such as `${{ true }}` there, so any other value counts as masking.
fn masks_failure(value: Option<&Value>) -> bool {
    value.is_some_and(|value| value.as_bool() != Some(false))
}

/// The repository script paths a `run:` block names.
fn script_paths(script: &str) -> Vec<String> {
    SCRIPT_PATH
        .captures_iter(script)
        .map(|captures| captures[1].to_string())
        .collect()
}

#[test]
fn every_script_a_workflow_runs_exists() {
    let mut references = 0;
    let mut missing = Vec::new();
    for step in run_steps() {
        for path in script_paths(&step.script) {
            references += 1;
            if !repo_root().join(&path).is_file() {
                missing.push(format!("{}: {path}", step.locus));
            }
        }
    }
    assert!(
        references >= MINIMUM_SCRIPT_REFERENCES,
        "found {references} script reference(s) in {WORKFLOW_DIR}; at least \
         {MINIMUM_SCRIPT_REFERENCES} are there today, so the scan has gone blind"
    );
    assert!(
        missing.is_empty(),
        "workflow steps run scripts that do not exist:\n  {}",
        missing.join("\n  ")
    );
}

#[test]
fn a_fetch_step_can_fail_its_job() {
    let fetch_steps: Vec<RunStep> = run_steps()
        .into_iter()
        .filter(|step| {
            script_paths(&step.script).iter().any(|path| {
                path.rsplit('/')
                    .next()
                    .is_some_and(|file| file.starts_with("fetch"))
            })
        })
        .collect();
    assert!(
        fetch_steps.len() >= MINIMUM_FETCH_STEPS,
        "found {} fetch step(s); at least {MINIMUM_FETCH_STEPS} are there today, so the scan \
         has gone blind",
        fetch_steps.len()
    );
    let masked: Vec<&str> = fetch_steps
        .iter()
        .filter(|step| step.continue_on_error)
        .map(|step| step.locus.as_str())
        .collect();
    assert!(
        masked.is_empty(),
        "fetch steps carry `continue-on-error: true`, so a broken fetch passes and the tests \
         reading its output check nothing:\n  {}",
        masked.join("\n  ")
    );
}

#[test]
fn the_script_path_pattern_finds_each_invocation_form() {
    assert_eq!(
        script_paths(
            "python scripts/fetch_x.py --out a\n./scripts/run.sh\nsource .github/scripts/v.sh"
        ),
        [
            "scripts/fetch_x.py",
            "scripts/run.sh",
            ".github/scripts/v.sh"
        ]
    );
    // A longer name that merely ends in `scripts/` is not a repository script.
    assert!(script_paths("python myscripts/x.py").is_empty());
}

#[test]
fn any_continue_on_error_but_a_literal_false_masks_a_failure() {
    let parse = |text: &str| serde_yaml::from_str::<Value>(text).expect("valid YAML");
    assert!(!masks_failure(None));
    assert!(!masks_failure(Some(&parse("false"))));
    assert!(masks_failure(Some(&parse("true"))));
    assert!(masks_failure(Some(&parse("'${{ true }}'"))));
    assert!(masks_failure(Some(&parse("'${{ matrix.experimental }}'"))));
}
