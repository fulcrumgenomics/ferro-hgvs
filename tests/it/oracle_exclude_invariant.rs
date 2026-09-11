//! Liveness check for the pairing between the spec corpus and the armed profile.
//!
//! The `oracle` profile in `.config/nextest.toml` is the policy half of `ci.yml`'s
//! `test-oracle` job. Its `default-filter` carries two exclusions, and its setup
//! script `scripts/arm-oracles.sh` arms the flags. CI and
//! `scripts/run_oracle_suite.sh` both select that profile, so this file reads the
//! profile and nothing reads `ci.yml` to recover it.
//!
//! Every module that measures over the spec corpus must be named in the first
//! exclusion, and every name there must still measure over the corpus. The second
//! exclusion is the debt list withheld from the denoted-sequence oracle; the
//! `oracle-rerun` profile must select exactly those rows and arm the other three
//! oracles. The tests below check each pairing in both directions.
//!
//! The failure mode is the flattering kind. A corpus module added later and not
//! named in the first exclusion does not go red; it reports a better census, which
//! reads as progress rather than as the lost evidence it is.
//!
//! This is not a coverage exemption. The corpus modules run unarmed in the plain
//! `test` job, and the corpus measures idempotency itself.
//!
//! See docs/ORACLES.md, section "What the oracle profile excludes".

use std::path::PathBuf;

/// This file. Excluded from the consumer scan below for the same reason
/// `sweep_filter_invariant.rs` excludes itself: the scan's **matcher literal**
/// is the import path it searches for, and that literal lives here. Without the
/// exclusion the guard would demand it be named in a CI filter that has no
/// reason to run it.
const SELF: &str = "oracle_exclude_invariant.rs";

/// The modules under `ferro_hgvs::conformance` whose import marks a test module
/// as built on the spec corpus.
///
/// `spec_corpus` is the corpus itself. **`census` is the corpus's measurement**,
/// added by #2063: the census was extracted out of `spec_conformance_axis` into a
/// library module so it could be run rather than only pinned, and
/// `conformance_census_instrument` therefore measures over the corpus while
/// importing nothing named `spec_corpus`. That is the scan's blind spot in its
/// flattering direction — an unmatched module is never demanded in the
/// exclusion, so its census is taken with the seam oracles armed and reads
/// better than the truth — and it is the exact shape this file's module doc
/// predicts ("a corpus module added later and not named in the first exclusion
/// does not go red, it reports a better census").
///
/// Path forms, not bare words, so a mention in prose — of which these modules have
/// many — does not count as consumption.
const CORPUS_MODULES: [&str; 2] = ["spec_corpus", "census"];

/// The brace-grouped form of the same imports, `conformance::{…spec_corpus…}`.
const CORPUS_IMPORT_GROUP: &str = "conformance::{";

/// Whether `text` imports any of [`CORPUS_MODULES`], in **either** spelling.
///
/// A single `contains(path)` misses `use
/// ferro_hgvs::conformance::{spec_corpus, summary};`, and it misses it in the
/// flattering direction: a module written that way is never demanded in the
/// exclusion, so its census is taken with the oracle armed and reads
/// *better* than the truth. That is the exact failure this file exists to close,
/// so the matcher may not be blind to a spelling rustfmt will happily produce.
///
/// Shared by both scans deliberately. Widening one and not the other is how the
/// same rule kept in two copies drifts — the reason [`CORPUS_MODULES`] is a
/// constant rather than an inline literal in the first place.
fn imports_corpus(text: &str) -> bool {
    if CORPUS_MODULES
        .iter()
        .any(|module| text.contains(&format!("conformance::{module}")))
    {
        return true;
    }
    // `conformance::{a, spec_corpus, b}` — read the group and look for the
    // module as a group item, so `spec_corpus_regressions` alongside it does not
    // count.
    //
    // Each item is reduced to its LEADING IDENTIFIER before comparing, which is
    // one rule covering three spellings that each escaped a narrower one, all in
    // the flattering direction — an unmatched module is never demanded in the
    // exclusion, so its census is taken with the seam oracles armed and reads
    // better than the truth:
    //
    // | item as written              | leading identifier |
    // |------------------------------|--------------------|
    // | `spec_corpus`                | `spec_corpus`      |
    // | `spec_corpus as corpus`      | `spec_corpus`      |
    // | `spec_corpus::{Frame, Row}`  | `spec_corpus`      |
    // | `spec_corpus_regressions`    | (no match)         |
    //
    // Whole-item equality missed the alias; splitting on whitespace alone still
    // missed the nested use tree. Stripping ` as …` and `::…` together is what
    // makes the three agree, and `spec_corpus_regressions` stays rejected
    // because the reduction never splits an identifier mid-word.
    //
    // Only GROUPED forms need any of this: `use …::conformance::spec_corpus as
    // corpus;` and `use …::conformance::spec_corpus::{…};` both contain the
    // path literal, so the check above already catches them.
    text.match_indices(CORPUS_IMPORT_GROUP).any(|(at, _)| {
        let rest = &text[at + CORPUS_IMPORT_GROUP.len()..];
        rest.find('}').is_some_and(|end| {
            rest[..end]
                .split(',')
                .any(|item| CORPUS_MODULES.contains(&leading_identifier(item)))
        })
    })
}

/// The module name a `use`-group item names, with any alias or nested path
/// removed: `spec_corpus as c` and `spec_corpus::{Frame}` both reduce to
/// `spec_corpus`.
fn leading_identifier(item: &str) -> &str {
    // No `.trim()` first: `split_whitespace` already skips leading whitespace,
    // and clippy rejects the pair.
    item.split_whitespace()
        .next()
        .unwrap_or("")
        .split("::")
        .next()
        .unwrap_or("")
        .trim_end_matches('{')
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

const NEXTEST_CONFIG: &str = ".config/nextest.toml";
const ORACLE_PROFILE: &str = "oracle";
const RERUN_PROFILE: &str = "oracle-rerun";

/// The setup script both profiles bind. It writes the `FERRO_ASSERT_*` flags to
/// `$NEXTEST_ENV`, so it is the only place the flag names live.
const ARM_SCRIPT: &str = "arm-oracles";

const SEQUENCE_FLAG: &str = "FERRO_ASSERT_SEQUENCE";

fn nextest_config() -> toml::Value {
    let path = repo_root().join(NEXTEST_CONFIG);
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    toml::from_str(&text).unwrap_or_else(|e| panic!("{NEXTEST_CONFIG} is not valid TOML: {e}"))
}

fn profile<'a>(config: &'a toml::Value, name: &str) -> &'a toml::Value {
    config
        .get("profile")
        .and_then(|profiles| profiles.get(name))
        .unwrap_or_else(|| panic!("`{NEXTEST_CONFIG}` defines no `[profile.{name}]`"))
}

/// A profile's `default-filter`, whitespace-normalised to one line.
fn default_filter(name: &str) -> String {
    profile(&nextest_config(), name)
        .get("default-filter")
        .and_then(toml::Value::as_str)
        .unwrap_or_else(|| {
            panic!(
                "`[profile.{name}]` has no `default-filter`. Without one the profile selects \
                 everything, and an armed run over everything is the known-red wall the \
                 profile exists to avoid."
            )
        })
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

/// The bodies of a filterset's `not (…)` clauses, in order.
///
/// A paren-depth scan rather than a filterset parser: the filter in play is
/// `not (A) and not (B)` where `A` and `B` are `test(…)` unions, and the two
/// call sites assert the clause they want exists.
fn not_clauses(filter: &str) -> Vec<String> {
    let mut clauses = Vec::new();
    let mut rest = filter;
    while let Some(at) = rest.find("not (") {
        let body = &rest[at + "not (".len()..];
        let mut depth = 1usize;
        let end = body
            .char_indices()
            .find_map(|(i, c)| {
                match c {
                    '(' => depth += 1,
                    ')' => depth -= 1,
                    _ => {}
                }
                (depth == 0).then_some(i)
            })
            .unwrap_or_else(|| panic!("unbalanced parentheses in filterset: {filter}"));
        clauses.push(body[..end].trim().to_string());
        rest = &body[end..];
    }
    clauses
}

/// One of the `oracle` profile's two exclusions, by position.
fn exclusion(index: usize, what: &str) -> String {
    let filter = default_filter(ORACLE_PROFILE);
    let clauses = not_clauses(&filter);
    clauses.get(index).cloned().unwrap_or_else(|| {
        panic!(
            "`[profile.{ORACLE_PROFILE}]`'s default-filter has {} `not (…)` clause(s), so the \
             {what} exclusion is not at position {index}: {filter}",
            clauses.len()
        )
    })
}

/// The spec-corpus modules: the `oracle` profile's first exclusion.
fn oracle_exclude() -> String {
    exclusion(0, "spec-corpus")
}

/// The rows withheld from the denoted-sequence oracle: the `oracle` profile's
/// second exclusion.
fn sequence_oracle_exclude() -> String {
    exclusion(1, "denoted-sequence")
}

/// Whether a profile binds [`ARM_SCRIPT`].
///
/// A profile with the exclusions and no binding runs the selection unarmed and
/// reports it as an oracle pass, which is the worse of the two ways to be wrong.
fn binds_arm_script(name: &str) -> bool {
    profile(&nextest_config(), name)
        .get("scripts")
        .and_then(toml::Value::as_array)
        .is_some_and(|bindings| {
            bindings.iter().any(|binding| {
                binding.get("setup").and_then(toml::Value::as_str) == Some(ARM_SCRIPT)
            })
        })
}

/// The `FERRO_ASSERT_*` keys [`ARM_SCRIPT`] arms for a profile, sorted.
///
/// Read by running the script the way nextest does — `NEXTEST_PROFILE` set and
/// `NEXTEST_ENV` pointing at a file it appends to — rather than by parsing its
/// text, so a rewrite of the script cannot leave this guard reading prose.
fn armed_flags(name: &str) -> Vec<String> {
    let command = nextest_config()
        .get("scripts")
        .and_then(|scripts| scripts.get("setup"))
        .and_then(|setup| setup.get(ARM_SCRIPT))
        .and_then(|script| script.get("command"))
        .and_then(toml::Value::as_str)
        .unwrap_or_else(|| {
            panic!("`{NEXTEST_CONFIG}` defines no `[scripts.setup.{ARM_SCRIPT}]` command")
        })
        .to_string();
    let env_file = std::env::temp_dir().join(format!("{ARM_SCRIPT}-{name}-{}", std::process::id()));
    let _ = std::fs::remove_file(&env_file);

    let output = std::process::Command::new(repo_root().join(&command))
        .env("NEXTEST_PROFILE", name)
        .env("NEXTEST_ENV", &env_file)
        .current_dir(repo_root())
        .output()
        .unwrap_or_else(|e| panic!("run {command}: {e}"));
    assert!(
        output.status.success(),
        "{command} failed under NEXTEST_PROFILE={name} ({}):\n{}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(&env_file)
        .unwrap_or_else(|e| panic!("{command} wrote nothing to $NEXTEST_ENV: {e}"));
    let _ = std::fs::remove_file(&env_file);

    let mut flags: Vec<String> = text
        .lines()
        .filter_map(|line| line.strip_suffix("=1"))
        .filter(|key| key.starts_with("FERRO_ASSERT_"))
        .map(str::to_string)
        .collect();
    assert!(
        !flags.is_empty(),
        "{command} armed no FERRO_ASSERT_* flag under NEXTEST_PROFILE={name}; the profile would \
         run its whole selection unarmed and report an oracle pass:\n{text}"
    );
    flags.sort();
    flags
}

/// `ci.yml`, parsed. Shared by the env-filter reader and the step guard.
fn ci_workflow() -> serde_yaml::Value {
    let path = repo_root().join(".github/workflows/ci.yml");
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    serde_yaml::from_str(&text).unwrap_or_else(|e| panic!("ci.yml is not valid YAML: {e}"))
}

/// A filter from `ci.yml`'s top-level `env:`.
///
/// `SWEEP_FILTER` and `CENSUS_FILTER` are scheduling, which CI owns; the profile
/// is policy, which this repo owns. The two must stay disjoint, so the guards
/// below read both.
fn ci_env_filter(key: &str) -> String {
    ci_workflow()
        .get("env")
        .and_then(|env| env.get(key))
        .and_then(serde_yaml::Value::as_str)
        .unwrap_or_else(|| panic!("ci.yml's top-level `env:` defines no `{key}`"))
        .to_string()
}

/// One step of the `test-oracle` job, found by its `name:`.
fn test_oracle_step(name: &str) -> serde_yaml::Value {
    ci_workflow()["jobs"]["test-oracle"]["steps"]
        .as_sequence()
        .unwrap_or_else(|| panic!("ci.yml's `test-oracle` job has no steps"))
        .iter()
        .find(|step| step["name"].as_str() == Some(name))
        .cloned()
        .unwrap_or_else(|| panic!("ci.yml's `test-oracle` job has no step named {name:?}"))
}

/// The profile a step's `run:` selects: the token after `--profile`.
fn profile_selected_by(step: &serde_yaml::Value) -> Option<String> {
    let mut tokens = step["run"].as_str()?.split_whitespace();
    while let Some(token) = tokens.next() {
        if token == "--profile" {
            return tokens.next().map(str::to_string);
        }
    }
    None
}

/// The profiles arm the oracles, so each `test-oracle` step must select its
/// profile by name and set no `FERRO_ASSERT_*` flag of its own. A run line that
/// dropped `--profile oracle` would run unarmed over an unfiltered selection and
/// go green: the vacuous pass this repository keeps guarding against.
#[test]
fn each_test_oracle_step_selects_its_profile_and_arms_nothing_itself() {
    let armed = test_oracle_step("Run Rust tests with the normalization self-checks");
    let rerun =
        test_oracle_step("Run the sequence-oracle exclusions under the other three oracles");

    assert_eq!(
        profile_selected_by(&armed).as_deref(),
        Some("oracle"),
        "the armed `test-oracle` step must pass `--profile oracle` on its `run:` line in \
         .github/workflows/ci.yml; the profile carries the exclusions and arms the flags, so a \
         step without it runs unarmed and green"
    );
    assert_eq!(
        profile_selected_by(&rerun).as_deref(),
        Some("oracle-rerun"),
        "the re-run `test-oracle` step must pass `--profile oracle-rerun` on its `run:` line in \
         .github/workflows/ci.yml; the profile carries the exclusions and arms the flags, so a \
         step without it runs unarmed and green"
    );

    let rerun_run = rerun["run"].as_str().expect("the re-run step has a `run:`");
    assert!(
        rerun_run
            .split_whitespace()
            .any(|token| token == "--no-tests=fail"),
        "the re-run step must carry --no-tests=fail, so a renamed row empties it loudly"
    );

    for (label, step) in [("armed", &armed), ("re-run", &rerun)] {
        let stray: Vec<&str> = step["env"]
            .as_mapping()
            .map(|env| {
                env.keys()
                    .filter_map(serde_yaml::Value::as_str)
                    .filter(|key| key.starts_with("FERRO_ASSERT_"))
                    .collect()
            })
            .unwrap_or_default();
        assert!(
            stray.is_empty(),
            "the {label} step sets {stray:?} in its own env; the profile owns the flags"
        );
    }
}

/// The local runner must select the `oracle` profile too.
///
/// `scripts/run_oracle_suite.sh` is the other consumer of the `oracle` profile,
/// beside `ci.yml`'s `test-oracle` job. It runs the suite twice, on a `list` line
/// and a `run` line, and each must pass `--profile oracle`, or the local run is
/// the vacuous kind: unarmed over an unfiltered selection, green for the wrong
/// reason. Read only the script here; the CI steps are the other test's job.
#[test]
fn the_local_runner_selects_the_oracle_profile() {
    let path = repo_root().join("scripts/run_oracle_suite.sh");
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let nextest_lines: Vec<&str> = text
        .lines()
        .filter(|line| line.contains("cargo nextest"))
        .collect();
    assert!(
        nextest_lines.len() >= 2,
        "scripts/run_oracle_suite.sh has {} `cargo nextest` line(s); it must have at least the \
         `list` line and the `run` line",
        nextest_lines.len()
    );
    for line in nextest_lines {
        assert!(
            line.contains("--profile oracle"),
            "scripts/run_oracle_suite.sh has a `cargo nextest` line without `--profile oracle`: \
             {line:?}\nThe profile carries the exclusions and arms the flags, so without it a \
             local run is unarmed and green."
        );
    }
}

/// Integration-test modules built on the spec corpus, by module name.
fn corpus_modules() -> Vec<String> {
    let dir = repo_root().join("tests/it");
    let mut modules: Vec<String> = std::fs::read_dir(&dir)
        .expect("tests/it is readable")
        .filter_map(|entry| {
            let path = entry.expect("readable dir entry").path();
            let name = path.file_name()?.to_str()?.to_string();
            if !name.ends_with(".rs") || name == SELF {
                return None;
            }
            let text = std::fs::read_to_string(&path).ok()?;
            imports_corpus(&text).then(|| name.trim_end_matches(".rs").to_string())
        })
        .collect();
    modules.sort();
    modules
}

/// Module names a nextest filter expression names via `test(...)`.
fn modules_named_in(filter: &str) -> Vec<String> {
    filter
        .match_indices("test(")
        .filter_map(|(at, _)| {
            let rest = &filter[at + "test(".len()..];
            rest.find(')').map(|end| rest[..end].trim().to_string())
        })
        .collect()
}

/// Every module that measures over the spec corpus must be named in the
/// `oracle` profile's first exclusion, or its census is taken with an oracle
/// armed and reads better than the truth.
#[test]
fn every_spec_corpus_module_is_named_in_the_oracle_exclude() {
    let filter = oracle_exclude();
    let modules = corpus_modules();
    assert!(
        !modules.is_empty(),
        "no module is built on the spec corpus; either the scan broke or the \
         corpus was removed — both should fail loudly rather than pass vacuously"
    );

    let missing: Vec<&String> = modules
        .iter()
        .filter(|module| !filter.contains(&format!("test({module})")))
        .collect();
    assert!(
        missing.is_empty(),
        "these modules measure over the spec corpus but are not named in the `oracle` \
         profile's first exclusion in {NEXTEST_CONFIG}, so `test-oracle` runs them with \
         the seam oracles armed. A panicking row contributes no output, which does not \
         redden the job — it makes confluence read HIGHER than it is: {missing:#?}\n\
         The exclusion is: {filter}\n\
         Add `+ test(<module>)` to it, or stop measuring over the corpus."
    );
}

/// The converse: a name in the first exclusion that no longer measures over the
/// corpus is withholding a module from the armed job for no reason.
#[test]
fn every_module_named_in_the_oracle_exclude_measures_over_the_corpus() {
    let filter = oracle_exclude();
    let modules = corpus_modules();

    let named = modules_named_in(&filter);
    assert!(
        !named.is_empty(),
        "the `oracle` profile's first exclusion names no modules; its formatting changed: {filter}"
    );

    let stale: Vec<&String> = named
        .iter()
        .filter(|module| !modules.contains(module))
        .collect();
    assert!(
        stale.is_empty(),
        "the `oracle` profile's first exclusion in {NEXTEST_CONFIG} names these modules, \
         but they do not measure over the spec corpus — so they are being withheld from \
         the armed job for no reason: {stale:#?}\n\
         Remove them from the exclusion."
    );
}

/// The scan's one structural blind spot, closed by forbidding the route rather
/// than by widening the matcher.
///
/// [`corpus_modules`] recognises a consumer by a literal import path from
/// [`CORPUS_MODULES`], which is right for avoiding prose false positives but
/// cannot see a module that reaches the corpus **indirectly** — through a shared
/// helper in `tests/it/common/` that does the importing. Such a module would
/// carry no matching literal, would not be demanded in the exclusion, and
/// would then measure with the seam oracles armed. Per this file's module doc
/// that does not go red; it reports a better census.
///
/// There is no such helper today: every corpus consumer defines its own `built()`
/// and imports the corpus directly. So the gap is closed at the only place it
/// could open — a `common/` helper importing the corpus is refused outright,
/// which is a mechanical guarantee rather than a list of symbol names that would
/// have to be maintained alongside the helpers it names.
#[test]
fn no_shared_helper_hides_a_corpus_consumer_from_the_scan() {
    let dir = repo_root().join("tests/it/common");
    let mut importers: Vec<String> = std::fs::read_dir(&dir)
        .expect("tests/it/common is readable")
        .filter_map(|entry| {
            let path = entry.expect("readable dir entry").path();
            let name = path.file_name()?.to_str()?.to_string();
            if !name.ends_with(".rs") {
                return None;
            }
            imports_corpus(&std::fs::read_to_string(&path).ok()?).then_some(name)
        })
        .collect();
    importers.sort();

    assert!(
        importers.is_empty(),
        "these shared helpers import the spec corpus: {importers:#?}\n\
         A module consuming the corpus THROUGH one of them carries no \
         `conformance::…` literal of its own for any of {CORPUS_MODULES:?}, so \
         `every_spec_corpus_module_is_named_in_the_oracle_exclude` would not \
         demand it be named in the exclusion — and it would then measure with \
         the seam oracles armed, which reports a better census rather than \
         going red.\n\
         Either import the corpus directly in each consuming module, or teach \
         `corpus_modules` to follow this helper."
    );
}

/// Every row withheld from the denoted-sequence oracle must name something that
/// exists.
///
/// `test()` is a substring predicate, so a typo does not error — it selects
/// nothing. In the `oracle` profile that fails in the **loud** direction (the
/// armed run stops withholding the row and goes red on it), which is the good
/// case and is why this is a cheap guard rather than a critical one. The
/// expensive half is `oracle-rerun`: `--no-tests=fail` turns a filter that
/// selects nothing into a red step, so a typo there is loud too, but a filter
/// with *one* good row and one typo would keep the step green while silently
/// dropping half of what it exists to re-run.
///
/// Checked against the **source tree** rather than by listing tests, so the guard
/// does not need a nextest subprocess: each name must be either an integration
/// module (`tests/it/<name>.rs`) or a `fn <name>(` somewhere under `tests/` or
/// `examples/`. That is a floor and not a proof — it cannot tell a `#[test]` from
/// a helper — but it catches the failure that actually happens, which is a
/// misspelling or a rename.
#[test]
fn every_row_withheld_from_the_sequence_oracle_names_something_that_exists() {
    let named = modules_named_in(&sequence_oracle_exclude());
    assert!(
        !named.is_empty(),
        "the `oracle` profile's second exclusion names nothing; either its formatting \
         changed (in which case fix it) or its last row retired — in which case delete \
         that clause, the `oracle-rerun` profile, its step in test-oracle, and these \
         guards, in one change."
    );

    let mut haystack = String::new();
    for dir in ["tests", "examples"] {
        collect_rust_sources(&repo_root().join(dir), &mut haystack);
    }
    assert!(
        haystack.len() > 100_000,
        "the source scan read only {} bytes, which cannot be the whole of tests/ and \
         examples/ — this guard has gone vacuous",
        haystack.len()
    );

    let missing: Vec<&String> = named
        .iter()
        .filter(|name| {
            !repo_root().join(format!("tests/it/{name}.rs")).is_file()
                && !haystack.contains(&format!("fn {name}("))
        })
        .collect();
    assert!(
        missing.is_empty(),
        "the `oracle` profile's second exclusion names these, and nothing in tests/ or \
         examples/ defines them: {missing:#?}\n\
         A `test()` term that matches nothing selects nothing, so `oracle-rerun` would \
         silently stop re-running the row it is there to protect."
    );
}

/// Append every `.rs` file under `dir`, recursively, to `into`.
fn collect_rust_sources(dir: &std::path::Path, into: &mut String) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            collect_rust_sources(&path, into);
        } else if path.extension().is_some_and(|e| e == "rs") {
            if let Ok(text) = std::fs::read_to_string(&path) {
                into.push_str(&text);
                into.push('\n');
            }
        }
    }
}

/// The rows withheld from the denoted-sequence oracle must still run under the
/// other three.
///
/// This is the guard that keeps arming the fourth oracle a **superset**, not a
/// trade. The `oracle` profile's second exclusion withdraws those rows from all
/// four oracles at once; the `oracle-rerun` profile puts three of them back
/// (#1815). Delete that profile, or let its filter drift from the exclusion, and
/// the change quietly becomes "three oracles surrendered to gain one", with
/// nothing red, which is why it is asserted rather than left to the comment
/// beside it.
///
/// Three things are checked, and the third is the one a reader would omit:
///
/// 1. `oracle-rerun` selects **exactly** the second exclusion. The two spellings
///    cannot reference each other, so they are compared here.
/// 2. both profiles bind the arm script. A profile with the exclusions and no
///    binding runs its selection unarmed and reports an oracle pass.
/// 3. the script arms `FERRO_ASSERT_SEQUENCE` for `oracle` and exactly the other
///    flags for `oracle-rerun`. Arming the fourth there would fire on the very
///    rows the profile exists to run, and the failure would read as a fresh
///    normalizer defect rather than as this wiring mistake.
#[test]
fn the_sequence_oracle_exclusions_still_run_under_the_other_three_oracles() {
    let withheld = sequence_oracle_exclude();
    assert_eq!(
        default_filter(RERUN_PROFILE),
        withheld,
        "`[profile.{RERUN_PROFILE}]`'s default-filter must equal the `oracle` profile's second \
         exclusion, row for row. Filtersets cannot reference each other, so the list is \
         spelled twice in {NEXTEST_CONFIG}; a drift here re-runs the wrong rows under the \
         other three oracles while the armed run withholds the right ones from all four."
    );

    for name in [ORACLE_PROFILE, RERUN_PROFILE] {
        assert!(
            binds_arm_script(name),
            "`[profile.{name}]` has no `[[profile.{name}.scripts]]` entry binding `{ARM_SCRIPT}`, \
             so it would run its selection with no oracle armed and report an oracle pass."
        );
    }

    let armed = armed_flags(ORACLE_PROFILE);
    assert!(
        armed.iter().any(|flag| flag == SEQUENCE_FLAG),
        "scripts/arm-oracles.sh no longer arms {SEQUENCE_FLAG} for the `{ORACLE_PROFILE}` \
         profile, so there is nothing for the second exclusion to withhold it from. If the \
         flag is being un-armed, remove that clause, the `{RERUN_PROFILE}` profile and this \
         guard in the same change rather than leaving a re-run for an exclusion that excludes \
         nothing. Armed: {armed:?}"
    );
    let expected: Vec<String> = armed
        .iter()
        .filter(|flag| *flag != SEQUENCE_FLAG)
        .cloned()
        .collect();
    assert_eq!(
        armed_flags(RERUN_PROFILE),
        expected,
        "the `{RERUN_PROFILE}` profile must arm exactly the flags `{ORACLE_PROFILE}` arms MINUS \
         {SEQUENCE_FLAG}. Deriving the expectation from the `{ORACLE_PROFILE}` profile rather \
         than from a literal list is deliberate — a fifth oracle added there must reach these \
         rows too, and a hardcoded three would not notice."
    );
}

/// The debt list must stay disjoint from the permanent filters.
///
/// An overlap is not merely untidy here, unlike the sweep/oracle pair below. A
/// row in both of the `oracle` profile's exclusions is withheld from the armed
/// run twice — harmless — but `oracle-rerun` would then run it under three
/// oracles that the armed run deliberately never runs it under, because the
/// first exclusion's whole point is that those instruments destroy each other
/// on those modules. So the overlap would *create* the red that the first
/// exclusion exists to prevent, in a new step, for a reason nothing states.
///
/// The `SWEEP_FILTER` half is the milder version: those rows run in `sweeps` with
/// all four oracles armed already, so re-running them here would be redundant
/// rather than wrong — but it would also mean the debt list silently governs a
/// job it says nothing about.
#[test]
fn the_sequence_oracle_exclude_is_disjoint_from_the_permanent_filters() {
    let debt = modules_named_in(&sequence_oracle_exclude());
    assert!(
        !debt.is_empty(),
        "the `oracle` profile's second exclusion names nothing; this guard has gone vacuous"
    );
    let permanent = [
        ("the `oracle` profile's first exclusion", oracle_exclude()),
        ("SWEEP_FILTER", ci_env_filter("SWEEP_FILTER")),
        ("CENSUS_FILTER", ci_env_filter("CENSUS_FILTER")),
    ];
    for (label, filter) in permanent {
        let other = modules_named_in(&filter);
        assert!(
            !other.is_empty(),
            "{label} names no module; this guard has gone vacuous"
        );
        let both: Vec<&String> = debt.iter().filter(|m| other.contains(m)).collect();
        assert!(
            both.is_empty(),
            "these are named in BOTH the second exclusion and {label}: {both:#?}\n\
             The debt list is temporary and carries issue numbers; {label} is a standing \
             statement about what the armed run must never run. A row in both means \
             `oracle-rerun` re-runs, under three armed oracles, a module {label} withholds \
             from them."
        );
    }
}

/// The two filters must stay disjoint.
///
/// They are negated together in one expression, so an overlap is harmless today
/// — but it means one job's list silently governs the other's, and the next
/// person editing either would reasonably read them as independent.
#[test]
fn the_sweep_filter_and_the_oracle_exclude_are_disjoint() {
    let sweeps = modules_named_in(&ci_env_filter("SWEEP_FILTER"));
    let excluded = modules_named_in(&oracle_exclude());

    let both: Vec<&String> = sweeps.iter().filter(|m| excluded.contains(m)).collect();
    assert!(
        both.is_empty(),
        "these modules are named in BOTH ci.yml's SWEEP_FILTER and the `oracle` profile's \
         first exclusion, which are meant to select disjoint sets for different reasons: \
         {both:#?}"
    );
}

/// [`imports_corpus`] recognises both spellings, and neither a prose mention nor
/// a neighbouring module whose name merely starts the same way.
///
/// Written because the grouped form was the gap: the scan matched a single
/// literal, so `use ferro_hgvs::conformance::{spec_corpus, summary};` read as
/// "does not consume the corpus" and the module was never demanded in the
/// exclusion.
///
/// The **aliased grouped** form was the same gap a second time. Widening the
/// matcher to read groups was not enough, because it then compared the whole
/// item for equality and `spec_corpus as corpus` is not equal to `spec_corpus`.
/// Both misses point the flattering way, which is why each spelling is pinned
/// here rather than left to the reader of `imports_corpus` to imagine.
#[test]
fn the_corpus_import_scan_sees_every_spelling() {
    assert!(imports_corpus(
        "use ferro_hgvs::conformance::spec_corpus::{denotation_of, Frame};"
    ));
    assert!(imports_corpus(
        "use ferro_hgvs::conformance::{spec_corpus, summary};"
    ));
    assert!(imports_corpus(
        "use ferro_hgvs::conformance::{\n    summary,\n    spec_corpus,\n};"
    ));
    // Aliased, grouped — the form an equality test on the whole item misses.
    assert!(imports_corpus(
        "use ferro_hgvs::conformance::{spec_corpus as corpus, summary};"
    ));
    assert!(imports_corpus(
        "use ferro_hgvs::conformance::{\n    summary,\n    spec_corpus as corpus,\n};"
    ));
    // Aliased, un-grouped — caught by the path literal, not by the group scan.
    assert!(imports_corpus(
        "use ferro_hgvs::conformance::spec_corpus as corpus;"
    ));
    // NESTED use trees — the third spelling, and the one a whitespace-only
    // split still missed: the item reads `spec_corpus::{Frame`, whose first
    // whitespace segment is not `spec_corpus`.
    assert!(imports_corpus(
        "use ferro_hgvs::conformance::{spec_corpus::{Frame, Row}, summary};"
    ));
    assert!(imports_corpus(
        "use ferro_hgvs::conformance::{summary, spec_corpus::{Frame, Row}};"
    ));
    assert!(imports_corpus(
        "use ferro_hgvs::conformance::{\n    summary,\n    spec_corpus::{Frame, Row},\n};"
    ));

    // `census` is the corpus's MEASUREMENT, and importing it is importing the
    // corpus transitively — the route `conformance_census_instrument` takes, and
    // the reason `CORPUS_MODULES` is a list rather than one path. Every spelling
    // above is pinned for it too, because the widening is worth nothing if it
    // only covers the spelling that module happens to use today.
    assert!(imports_corpus(
        "use ferro_hgvs::conformance::census::{measure, Equivalence};"
    ));
    assert!(imports_corpus(
        "use ferro_hgvs::conformance::{census, summary};"
    ));
    assert!(imports_corpus(
        "use ferro_hgvs::conformance::{census as c, summary};"
    ));
    assert!(imports_corpus(
        "use ferro_hgvs::conformance::{summary, census::{measure, Census}};"
    ));
    assert!(imports_corpus("use ferro_hgvs::conformance::census as c;"));

    // Prose, not consumption — the reason the matcher is a path and not the
    // bare word.
    assert!(!imports_corpus(
        "//! The spec_corpus generator builds these rows."
    ));
    // A different module that merely shares a prefix inside a group. Still
    // rejected after the widening: the alias split is on whitespace, and
    // `spec_corpus_regressions` is a single segment.
    assert!(!imports_corpus(
        "use ferro_hgvs::conformance::{spec_corpus_regressions, summary};"
    ));
    assert!(!imports_corpus(
        "use ferro_hgvs::conformance::{spec_corpus_regressions as regressions, summary};"
    ));
    // …and nested, which is where a reduction that split identifiers on prefix
    // rather than on a separator would go wrong.
    assert!(!imports_corpus(
        "use ferro_hgvs::conformance::{spec_corpus_regressions::{Frame}, summary};"
    ));
    // The same prefix hazard for the second entry, so widening the list did not
    // widen what counts as a match.
    assert!(!imports_corpus(
        "use ferro_hgvs::conformance::{census_filter_invariant, summary};"
    ));
    assert!(!imports_corpus(
        "use ferro_hgvs::conformance::{census_filter_invariant as c, summary};"
    ));
}
