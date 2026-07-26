//! The spec generators must be usable on a **fresh worktree**, and must blame
//! the right thing when they cannot run.
//!
//! Two failure modes made `git push` impossible on a newly created worktree,
//! and neither named its actual cause:
//!
//! 1. With the `assets/hgvs-nomenclature` submodule uninitialised, the harvest
//!    yielded only the compiled-in `ported_legacy_probes`, so `build_rows`
//!    validated the curated overrides against a near-empty input set and
//!    reported `overrides reference unknown inputs (typo or spec drift)` — a
//!    committed, correct file accused of being stale because a directory was
//!    empty.
//! 2. `--check` then aborted with a bare `No such file or directory` for the
//!    gitignored output artifact, which on a fresh worktree has simply never
//!    been generated.
//!
//! These tests pin both: the harvest fails with an actionable message naming the
//! submodule, and `--check` treats an absent artifact as a cold cache rather
//! than drift.

use std::path::PathBuf;
use std::process::{Command, Output};

fn manifest_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// A unique scratch directory under `target/`, removed and recreated per test so
/// reruns are deterministic. Kept inside `target/` so it is gitignored and never
/// pollutes the worktree.
fn scratch(name: &str) -> PathBuf {
    let dir = manifest_dir()
        .join("target")
        .join("spec-precond")
        .join(name);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    dir
}

/// Build the two generator examples exactly once per test binary, and return
/// the path to the requested one.
///
/// # Why not just `cargo run --example`
///
/// That is what this module did originally, and it was pathologically slow —
/// but not for the reason it looks like. Cargo takes a **file lock on
/// `target/`**, so the four tests here, which nextest runs concurrently,
/// serialised behind one another: each reported ~25s even against a fully warm
/// build tree, and in CI they blew past nextest's 60s SLOW threshold and made
/// whichever shard held them the critical path of the whole test job.
///
/// # Why not just invoke `target/<profile>/examples/<name>` directly
///
/// Tried, and it is **wrong**. `cargo test` / `cargo nextest run` do *not*
/// rebuild examples, so the binary sitting in `examples/` can be arbitrarily
/// old. Observed directly: a `generate_spec_enumeration` binary from the
/// previous day produced the pre-#1201 error message, and this module's test
/// for #1201's guard failed against source that unmistakably contained it. A
/// stale binary makes these tests assert things about code that is not the code
/// under test — a far worse failure than being slow.
///
/// # What this does
///
/// One `cargo build` for both examples, behind a `Once`, then direct execution.
/// The build is current by construction, and it happens exactly once no matter
/// how many tests run concurrently — so the cargo lock is taken once rather
/// than per test, which is where the speedup comes from.
///
/// `FERRO_PREBUILT_EXAMPLES` skips even that single build, for the one context
/// where the binaries are known-current and cargo is unavailable: CI's sharded
/// test jobs, which run from a nextest archive after `test-build` has already
/// built and executed both generators.
fn built_example(example: &str) -> PathBuf {
    static BUILD: std::sync::Once = std::sync::Once::new();
    static BUILD_OK: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

    // CI builds both examples in the `test-build` job (it runs them to generate
    // the spec artifacts) and sets this variable, so the sharded test jobs —
    // which execute from a nextest archive and have no cargo build state at all
    // — must not shell out to cargo. A cargo invocation there would rebuild the
    // world per shard.
    let prebuilt = std::env::var_os("FERRO_PREBUILT_EXAMPLES").is_some();

    BUILD.call_once(|| {
        if prebuilt {
            BUILD_OK.store(true, std::sync::atomic::Ordering::SeqCst);
            return;
        }
        let status = Command::new(env!("CARGO"))
            .current_dir(manifest_dir())
            .args([
                "build",
                "--quiet",
                "--features",
                "dev",
                "--example",
                "generate_spec_fixture",
                "--example",
                "generate_spec_enumeration",
            ])
            .status();
        let ok = matches!(status, Ok(s) if s.success());
        BUILD_OK.store(ok, std::sync::atomic::Ordering::SeqCst);
    });
    assert!(
        BUILD_OK.load(std::sync::atomic::Ordering::SeqCst),
        "failed to build the generator examples"
    );

    // The test binary lives at `…/target/<profile>/deps/it-<hash>`, so the
    // examples directory is two levels up. Deriving it from `current_exe`
    // rather than hardcoding `debug/` keeps this correct under `--release` and
    // under a nextest archive extracted into the workspace.
    let exe = std::env::current_exe().expect("current_exe");
    let path = exe
        .parent()
        .and_then(|p| p.parent())
        .expect("target/<profile>")
        .join("examples")
        .join(example);
    assert!(
        path.is_file(),
        "expected a freshly built example at {}",
        path.display()
    );
    path
}

/// Run a `--features dev` example with `args`, returning its captured output.
fn run_example(example: &str, args: &[&str]) -> Output {
    Command::new(built_example(example))
        .current_dir(manifest_dir())
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("failed to run `{example}`: {e}"))
}

fn stderr_of(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).to_string()
}

/// One pre-commit hook, flattened out of the nested `repos[].hooks[]` shape.
struct Hook {
    id: String,
    entry: String,
    stages: Vec<String>,
}

/// Every hook in `.pre-commit-config.yaml`, in declared order.
///
/// Parsed as YAML rather than substring-matched: the file's comments discuss
/// `--check` on purpose, so a naive `contains` can match prose, and a hook that
/// was commented out would still "exist".
fn pre_commit_hooks() -> Vec<Hook> {
    let text = std::fs::read_to_string(manifest_dir().join(".pre-commit-config.yaml"))
        .expect("read .pre-commit-config.yaml");
    let doc: serde_yaml::Value =
        serde_yaml::from_str(&text).expect("parse .pre-commit-config.yaml");

    doc.get("repos")
        .and_then(serde_yaml::Value::as_sequence)
        .expect("`repos` sequence")
        .iter()
        .filter_map(|repo| repo.get("hooks").and_then(serde_yaml::Value::as_sequence))
        .flatten()
        .map(|hook| Hook {
            id: hook
                .get("id")
                .and_then(serde_yaml::Value::as_str)
                .unwrap_or_default()
                .to_string(),
            entry: hook
                .get("entry")
                .and_then(serde_yaml::Value::as_str)
                .unwrap_or_default()
                .to_string(),
            stages: hook
                .get("stages")
                .and_then(serde_yaml::Value::as_sequence)
                .map(|s| {
                    s.iter()
                        .filter_map(serde_yaml::Value::as_str)
                        .map(str::to_string)
                        .collect()
                })
                .unwrap_or_default(),
        })
        .collect()
}

/// Every `run:` command across every step of every job in `ci.yml`.
fn ci_run_commands() -> Vec<String> {
    let text = std::fs::read_to_string(manifest_dir().join(".github/workflows/ci.yml"))
        .expect("read .github/workflows/ci.yml");
    let doc: serde_yaml::Value = serde_yaml::from_str(&text).expect("parse ci.yml");

    doc.get("jobs")
        .and_then(serde_yaml::Value::as_mapping)
        .expect("`jobs` mapping")
        .values()
        .filter_map(|job| job.get("steps").and_then(serde_yaml::Value::as_sequence))
        .flatten()
        .filter_map(|step| step.get("run").and_then(serde_yaml::Value::as_str))
        .map(str::to_string)
        .collect()
}

/// The real spec checkout, or `None` when the submodule is not initialised —
/// in which case the tests that need real inputs are skipped rather than failed,
/// matching how the rest of the suite treats absent optional inputs.
fn spec_dir_if_initialised() -> Option<PathBuf> {
    let dir = manifest_dir().join("assets/hgvs-nomenclature");
    dir.join("docs/recommendations").is_dir().then_some(dir)
}

// ---------------------------------------------------------------------------
// 1. An empty spec checkout blames the submodule, not the overrides file
// ---------------------------------------------------------------------------

/// An empty `--spec-dir` must fail with a message that names the empty checkout
/// and the command that fixes it — and must NOT report the overrides file as
/// carrying unknown inputs, which is what sent the last reader after phantom
/// spec drift.
#[test]
fn empty_spec_checkout_names_the_submodule_not_the_overrides() {
    let empty = scratch("empty-spec-dir");
    let out = run_example(
        "generate_spec_fixture",
        &[
            "--spec-dir",
            empty.to_str().expect("utf-8 scratch path"),
            "--output",
            empty
                .join("unused.json")
                .to_str()
                .expect("utf-8 scratch path"),
        ],
    );
    let stderr = stderr_of(&out);

    assert!(
        !out.status.success(),
        "an empty spec checkout must fail, got success. stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("no HGVS strings harvested"),
        "error must say the harvest came up empty. stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("git submodule update --init assets/hgvs-nomenclature"),
        "error must name the command that fixes it. stderr:\n{stderr}"
    );
    assert!(
        !stderr.contains("overrides reference unknown inputs"),
        "an empty checkout must not be reported as stale overrides — that \
         misattribution is the bug this pins. stderr:\n{stderr}"
    );
}

/// The same guard protects the enumeration generator, which harvests the same
/// checkout through the same `sources::discover`.
#[test]
fn empty_spec_checkout_is_rejected_by_the_enumeration_generator_too() {
    let empty = scratch("empty-spec-dir-enum");
    let out = run_example(
        "generate_spec_enumeration",
        &[
            "--spec-dir",
            empty.to_str().expect("utf-8 scratch path"),
            "--output",
            empty
                .join("unused.json")
                .to_str()
                .expect("utf-8 scratch path"),
        ],
    );
    let stderr = stderr_of(&out);

    assert!(
        !out.status.success(),
        "an empty spec checkout must fail, got success. stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("no HGVS strings harvested"),
        "error must say the harvest came up empty. stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("git submodule update --init assets/hgvs-nomenclature"),
        "error must name the command that fixes it. stderr:\n{stderr}"
    );
}

// ---------------------------------------------------------------------------
// 2. `--check` treats an absent artifact as a cold cache, not drift
// ---------------------------------------------------------------------------

/// This is the fresh-worktree case: `--check` pointed at an artifact that has
/// never been generated must succeed and create it, not abort on ENOENT. A
/// gitignored build artifact has no committed baseline to have drifted from.
#[test]
fn check_generates_an_absent_artifact_instead_of_failing() {
    let Some(spec_dir) = spec_dir_if_initialised() else {
        eprintln!("skipping: assets/hgvs-nomenclature not initialised");
        return;
    };
    let dir = scratch("check-absent");
    let output = dir.join("hgvs_spec_normalization.json");
    assert!(!output.exists(), "scratch output must start absent");

    let out = run_example(
        "generate_spec_fixture",
        &[
            "--spec-dir",
            spec_dir.to_str().expect("utf-8 spec path"),
            "--output",
            output.to_str().expect("utf-8 scratch path"),
            "--check",
        ],
    );
    let stderr = stderr_of(&out);

    assert!(
        out.status.success(),
        "--check on an absent artifact must succeed. stderr:\n{stderr}"
    );
    assert!(
        output.exists(),
        "--check must materialise the absent artifact. stderr:\n{stderr}"
    );
    assert!(
        !stderr.contains("No such file or directory"),
        "the raw ENOENT is the bug this pins. stderr:\n{stderr}"
    );
}

/// The staleness verdict itself must survive: an artifact that exists and
/// differs is still a failure, and is left on disk for inspection.
#[test]
fn check_still_reports_a_stale_artifact_and_leaves_it_alone() {
    let Some(spec_dir) = spec_dir_if_initialised() else {
        eprintln!("skipping: assets/hgvs-nomenclature not initialised");
        return;
    };
    let dir = scratch("check-stale");
    let output = dir.join("hgvs_spec_normalization.json");
    let stale = "{\"not\": \"the real fixture\"}\n";
    std::fs::write(&output, stale).expect("seed a stale artifact");

    let out = run_example(
        "generate_spec_fixture",
        &[
            "--spec-dir",
            spec_dir.to_str().expect("utf-8 spec path"),
            "--output",
            output.to_str().expect("utf-8 scratch path"),
            "--check",
        ],
    );
    let stderr = stderr_of(&out);

    assert!(
        !out.status.success(),
        "--check must fail on a stale artifact. stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("is out of date"),
        "error must say the artifact is out of date. stderr:\n{stderr}"
    );
    assert_eq!(
        std::fs::read_to_string(&output).expect("read seeded artifact"),
        stale,
        "--check must not rewrite a stale artifact; the old content stays \
         available for inspection",
    );
}

// ---------------------------------------------------------------------------
// 3. The hooks and CI must agree
// ---------------------------------------------------------------------------

/// The pre-push hooks exist to catch a broken committed input before a CI cycle
/// is spent, so they must run the same commands CI runs. `--check` answers a
/// different question (local artifact freshness) about a gitignored file, and
/// gating a push on it is what made a fresh worktree unpushable.
#[test]
fn pre_push_hooks_run_plain_generation_like_ci() {
    let hooks = pre_commit_hooks();
    let ci_runs = ci_run_commands();

    for example in ["generate_spec_fixture", "generate_spec_enumeration"] {
        let plain = format!("cargo run --features dev --example {example}");

        let hook_id = example.replace('_', "-");
        let hook = hooks
            .iter()
            .find(|h| h.id == hook_id)
            .unwrap_or_else(|| panic!("no pre-commit hook with id `{hook_id}`"));
        assert_eq!(
            hook.entry, plain,
            "hook `{}` must run plain generation, not a variant",
            hook.id
        );
        assert!(
            hook.stages.iter().any(|s| s == "pre-push"),
            "hook `{}` must be scoped to the pre-push stage, got {:?}",
            hook.id,
            hook.stages,
        );

        // The agreement has two sides. Asserting only the hook would let CI
        // drift back to `--check` without this guard noticing — and "the hooks
        // run what CI runs" is precisely the claim being pinned.
        assert!(
            ci_runs.iter().any(|r| r.trim() == plain),
            "expected CI to run plain `{example}`, matching the pre-push hook",
        );
    }

    // Inspect parsed invocations only — the comments in both files discuss
    // `--check` on purpose, explaining why it is not used here.
    let gated: Vec<String> = hooks
        .iter()
        .map(|h| (format!("hook {}", h.id), h.entry.clone()))
        .chain(ci_runs.iter().map(|r| ("ci.yml".to_string(), r.clone())))
        .filter(|(_, cmd)| cmd.contains("generate_spec_") && cmd.contains("--check"))
        .map(|(origin, cmd)| format!("{origin}: {cmd}"))
        .collect();
    assert!(
        gated.is_empty(),
        "neither the pre-push hooks nor CI may gate on `--check` for the \
         gitignored spec artifacts: {gated:?}",
    );
}

/// Guard against the dependency ordering silently regressing: the enumeration
/// generator reads the normalization fixture, and that input is satisfied only
/// because pre-commit runs hooks in declared order.
#[test]
fn fixture_hook_is_declared_before_the_enumeration_hook() {
    let hooks = pre_commit_hooks();
    let position = |id: &str| {
        hooks
            .iter()
            .position(|h| h.id == id)
            .unwrap_or_else(|| panic!("hook `{id}` is declared"))
    };
    assert!(
        position("generate-spec-fixture") < position("generate-spec-enumeration"),
        "the fixture hook must come first — it produces the input the \
         enumeration generator reads",
    );
}
