//! Every on-demand generated fixture must be built and shipped by CI.
//!
//! # The defect this closes, which shipped and was invisible for weeks
//!
//! `common::fixture_gen` gives a test module a fallback: if its generated
//! fixture is absent, shell out to `cargo run --features dev --example <gen>`
//! and build it. That fallback exists so a fresh `cargo test` just works, and
//! locally it costs a few seconds against a warm `target/`.
//!
//! In CI it is a trap. `test` and `test-oracle` run from a nextest **archive**
//! and have no `target/` at all, so the nested `cargo run` compiles the whole
//! crate *inside the test* — and nextest gives every test its own process, so
//! **every test in the module pays it**, not the first one.
//!
//! `cis_confluence_nr_axis` landed in #1646 with no wiring in `spec-fixtures`,
//! and nothing went red. Measured on run 31670057963, `Test oracle (1/3)`:
//! `the_rna_axis_is_spelled_in_the_rna_alphabet` took **129.2s** against 0.15s
//! once the file exists — a test that normalizes nothing and only compares
//! strings — while its `g,c` sibling `three_prime_confluence_census` ran
//! **twice** the spellings in 88.0s, because *its* corpus was downloaded.
//!
//! **The failure mode is silence, and the flattering kind of silence.** Nothing
//! fails, no artifact is missing, no assertion moves; the suite is simply slower
//! in a way that reads as "those censuses are expensive". That is the same shape
//! as `oracle_exclude_invariant.rs`'s: a module added later, not named in a CI
//! list, that degrades quietly rather than loudly.
//!
//! # What is checked, and what that is worth
//!
//! Three properties, each derived from the workflow rather than from a list
//! kept here — a hardcoded roster of fixtures would rot in exactly the way the
//! defect did:
//!
//! 1. every generated fixture path a `tests/it` module names is in **every**
//!    `actions/cache` step's `path:` in `spec-fixtures` (restore and save must
//!    agree, or the entry saved is not the entry restored);
//! 2. it is published on some `actions/upload-artifact` from that job;
//! 3. every job consuming the `nextest-archive` artifact — that is, every job
//!    that runs the `it` binary with no warm `target/` — downloads each artifact
//!    carrying a generated fixture.
//!
//! The workflow half is **parsed**, for the reason `workflow_gh_repo_context.rs`
//! gives: which step sits in which job is structural, and a substring scan
//! cannot answer it. The source half is a **token scan** and is a floor, not a
//! proof: it finds a `"tests/fixtures/….json"` literal in a file that also
//! mentions `ensure_generated`, so a module holding its path somewhere else, or
//! composing it at run time, is invisible to it. What it does buy is that a
//! fixture declared the way all four current ones are declared cannot be added
//! in silence, and the question gets asked.

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

use serde_yaml::Value;

/// This file. Excluded from the scan for the same reason
/// `oracle_exclude_invariant.rs` excludes itself: it carries [`CALLS`]'s matcher
/// literals, so it matches its own consumer test and would be read for fixture
/// paths it has no business declaring. It happens to declare none today —
/// [`FIXTURE_PREFIX`] is not `.json`-terminated, so [`fixture_literals`] finds
/// nothing here — but that is a property of the current wording, not a
/// guarantee, and the exclusion is what makes it one.
const SELF: &str = "generated_fixture_ci_wiring.rs";

/// The module that *defines* the fallback rather than using it. It names no
/// fixture, but it is the one file guaranteed to mention the call.
const DEFINITION: &str = "fixture_gen.rs";

/// The call that marks a module as depending on a generated fixture, in both its
/// spellings — the `[[bin]]` door and the `[[example]]` one.
const CALLS: [&str; 2] = [
    "ensure_generated_fixture(",
    "ensure_generated_example_fixture(",
];

/// The prefix of a fixture path literal, relative to the manifest directory —
/// the form `fixture_gen::fixture_path` takes.
const FIXTURE_PREFIX: &str = "tests/fixtures/";

/// The artifact carrying the test binaries that `test` and `test-oracle` run.
///
/// Consumers are derived from this rather than named, so a new job that runs the
/// suite from the archive inherits the requirement. `soak` and `sweeps` consume
/// `nextest-archive-soak` and a narrow selection, and are deliberately not in
/// scope: they do not run these modules.
const TEST_ARCHIVE: &str = "nextest-archive";

/// The job that generates and publishes the fixtures.
const PRODUCER_JOB: &str = "spec-fixtures";

/// A floor on what the scan must find, so a matcher that silently stops matching
/// cannot pass as "nothing to check". Four fixtures are wired today: the two
/// spec artifacts and the two cis confluence corpora.
const MINIMUM_FIXTURES: usize = 4;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Every `tests/it` source that reaches for a generated fixture, as
/// `(file name, contents)`.
fn fixture_consuming_sources() -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut stack = vec![repo_root().join("tests/it")];
    while let Some(dir) = stack.pop() {
        let entries = std::fs::read_dir(&dir)
            .unwrap_or_else(|e| panic!("read {}: {e}", dir.display()))
            .filter_map(Result::ok);
        for entry in entries {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
                continue;
            }
            if path.extension().is_none_or(|e| e != "rs") {
                continue;
            }
            let name = path
                .file_name()
                .expect("a file has a name")
                .to_string_lossy()
                .into_owned();
            if name == SELF || name == DEFINITION {
                continue;
            }
            let text = std::fs::read_to_string(&path)
                .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
            if CALLS.iter().any(|call| text.contains(call)) {
                out.push((name, text));
            }
        }
    }
    out.sort();
    out
}

/// The `tests/fixtures/….json` literals in `text`.
///
/// A deliberate token scan, not a parse. It reads from the opening quote of the
/// prefix to the next quote, so it picks up the `const … : &str = "…"` form
/// every current caller uses and ignores the same path mentioned in prose,
/// which has no quote around it.
fn fixture_literals(text: &str) -> BTreeSet<String> {
    let mut found = BTreeSet::new();
    let needle = format!("\"{FIXTURE_PREFIX}");
    let mut rest = text;
    while let Some(at) = rest.find(&needle) {
        let after = &rest[at + 1..];
        if let Some(end) = after.find('"') {
            let literal = &after[..end];
            if literal.ends_with(".json") {
                found.insert(literal.to_string());
            }
            rest = &after[end..];
        } else {
            break;
        }
    }
    found
}

/// Every generated fixture path a `tests/it` module names, and the module that
/// names it.
fn generated_fixtures() -> BTreeMap<String, String> {
    let sources = fixture_consuming_sources();
    assert!(
        !sources.is_empty(),
        "no tests/it source mentions {CALLS:?}; the matcher has stopped matching \
         and this guard is checking nothing"
    );
    let mut out = BTreeMap::new();
    for (name, text) in sources {
        for literal in fixture_literals(&text) {
            out.insert(literal, name.clone());
        }
    }
    assert!(
        out.len() >= MINIMUM_FIXTURES,
        "found only {} generated fixture path(s) ({:?}); at least {MINIMUM_FIXTURES} are wired \
         today, so the scan has gone blind rather than the repository having shrunk",
        out.len(),
        out.keys().collect::<Vec<_>>()
    );
    out
}

fn ci_workflow() -> Value {
    let path = repo_root().join(".github/workflows/ci.yml");
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    serde_yaml::from_str(&text).expect("ci.yml is not valid YAML")
}

fn job<'a>(workflow: &'a Value, name: &str) -> &'a Value {
    workflow
        .get("jobs")
        .and_then(|jobs| jobs.get(name))
        .unwrap_or_else(|| panic!("ci.yml has no `{name}` job"))
}

fn steps(job: &Value) -> &[Value] {
    job.get("steps")
        .and_then(Value::as_sequence)
        .map(Vec::as_slice)
        .unwrap_or(&[])
}

fn uses(step: &Value) -> &str {
    step.get("uses").and_then(Value::as_str).unwrap_or_default()
}

/// A step's `with.path`, split into lines — the multi-line block scalar form the
/// cache and artifact steps use, and the single-path form too.
fn with_paths(step: &Value) -> Vec<String> {
    step.get("with")
        .and_then(|with| with.get("path"))
        .and_then(Value::as_str)
        .map(|paths| {
            paths
                .lines()
                .map(str::trim)
                .filter(|line| !line.is_empty())
                .map(str::to_string)
                .collect()
        })
        .unwrap_or_default()
}

fn with_name(step: &Value) -> Option<&str> {
    step.get("with")
        .and_then(|with| with.get("name"))
        .and_then(Value::as_str)
}

/// Where a `download-artifact` step unpacks, with any trailing slash trimmed.
///
/// `actions/download-artifact` defaults to the workspace root when `path` is
/// unset, which is why the fallback is `.` rather than a panic — an omitted
/// `path:` is a legal step and a wrong destination, and [`landing_path`] is what
/// reports it.
fn with_dest(step: &Value) -> String {
    step.get("with")
        .and_then(|with| with.get("path"))
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|dest| !dest.is_empty())
        .unwrap_or(".")
        .trim_end_matches('/')
        .to_string()
}

/// Where `fixture` ends up when its artifact is unpacked at `dest`.
///
/// `upload-artifact` roots a multi-path artifact at its paths' common ancestor,
/// and both fixture artifacts here publish files from a single directory — so
/// each entry is a bare file name and the landing path is `dest/<file name>`.
fn landing_path(dest: &str, fixture: &str) -> String {
    let file_name = fixture.rsplit('/').next().unwrap_or(fixture);
    if dest == "." {
        file_name.to_string()
    } else {
        format!("{dest}/{file_name}")
    }
}

/// Every generated fixture is cached by `spec-fixtures` — by **every** cache
/// step in it, restore and save alike.
///
/// Split from the artifact check because the two fail for different reasons and
/// the fix differs: a path missing from the cache is re-generated on every run
/// (wasteful but correct), while a path missing from the artifact is not
/// generated in the consuming job at all (the defect this file records).
#[test]
fn every_generated_fixture_is_cached_by_the_producer_job() {
    let workflow = ci_workflow();
    let producer = job(&workflow, PRODUCER_JOB);

    let cache_steps: Vec<&Value> = steps(producer)
        .iter()
        .filter(|step| uses(step).starts_with("actions/cache"))
        .collect();
    assert!(
        !cache_steps.is_empty(),
        "`{PRODUCER_JOB}` has no actions/cache step; this check would pass vacuously"
    );

    for (fixture, module) in generated_fixtures() {
        for step in &cache_steps {
            assert!(
                with_paths(step).contains(&fixture),
                "{module} regenerates {fixture} on demand, but `{PRODUCER_JOB}`'s \
                 `{}` step does not list it. Every cache step in that job must \
                 carry the same paths, or the entry saved is not the entry \
                 restored.",
                uses(step)
            );
        }
    }
}

/// Every generated fixture is published from `spec-fixtures`, and every job that
/// runs the test archive downloads the artifact carrying it.
///
/// This is the one that would have failed on #1646.
#[test]
fn every_generated_fixture_reaches_the_jobs_that_run_the_tests() {
    let workflow = ci_workflow();
    let producer = job(&workflow, PRODUCER_JOB);

    // artifact name -> the paths it publishes
    let mut published: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for step in steps(producer) {
        if !uses(step).starts_with("actions/upload-artifact") {
            continue;
        }
        let name = with_name(step)
            .unwrap_or_else(|| panic!("`{PRODUCER_JOB}` uploads an artifact with no name"));
        published
            .entry(name.to_string())
            .or_default()
            .extend(with_paths(step));
    }
    assert!(
        !published.is_empty(),
        "`{PRODUCER_JOB}` uploads no artifacts; this check would pass vacuously"
    );

    // Which artifacts the consumers must therefore download, and onto which path
    // each fixture has to land.
    let mut required: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for (fixture, module) in generated_fixtures() {
        let carrier = published
            .iter()
            .find(|(_, paths)| paths.contains(&fixture))
            .map(|(name, _)| name.clone())
            .unwrap_or_else(|| {
                panic!(
                    "{module} regenerates {fixture} on demand, and `{PRODUCER_JOB}` \
                     publishes it on no artifact — so every job running from the \
                     nextest archive rebuilds the crate inside the test to produce \
                     it. Add it to an upload-artifact step's `path:`."
                )
            });
        required.entry(carrier).or_default().push(fixture);
    }

    // Consumers are derived, not named: any job downloading the test archive
    // runs the `it` binary with no warm `target/`.
    let mut consumers = 0usize;
    for (name, definition) in workflow
        .get("jobs")
        .and_then(Value::as_mapping)
        .expect("ci.yml has a jobs mapping")
    {
        // artifact name -> where that job unpacks it.
        let downloads: BTreeMap<&str, String> = steps(definition)
            .iter()
            .filter(|step| uses(step).starts_with("actions/download-artifact"))
            .filter_map(|step| with_name(step).map(|name| (name, with_dest(step))))
            .collect();
        if !downloads.contains_key(TEST_ARCHIVE) {
            continue;
        }
        consumers += 1;
        let job_name = name.as_str().unwrap_or("<unnamed>");
        for (artifact, fixtures) in &required {
            let dest = downloads.get(artifact.as_str()).unwrap_or_else(|| {
                panic!(
                    "job `{job_name}` runs the `{TEST_ARCHIVE}` test binaries but does \
                     not download `{artifact}`, which carries a generated fixture. \
                     Every test needing that fixture will rebuild the crate inside \
                     its own process to produce it."
                )
            });
            // Downloading it is not enough — it has to unpack where the test
            // reads it. `upload-artifact` roots a multi-path artifact at its
            // paths' common ancestor, so `cis-corpus` holds bare basenames and
            // only lands correctly because the consumer sets `path:
            // tests/fixtures/cis/`. A `path: .` would put every corpus in the
            // repository root and send every test straight back to the nested
            // build, with the name-only check above still green.
            for fixture in fixtures {
                let landed = landing_path(dest, fixture);
                assert_eq!(
                    &landed, fixture,
                    "job `{job_name}` downloads `{artifact}` to `{dest}`, so \
                     {fixture} lands at {landed} — not where the test reads it. \
                     Set the download `path:` to the fixture's own directory."
                );
            }
        }
    }
    assert!(
        consumers >= 2,
        "found {consumers} job(s) downloading `{TEST_ARCHIVE}`; `test` and \
         `test-oracle` both do, so the derivation has gone blind"
    );
}
