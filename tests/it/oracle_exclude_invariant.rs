//! Liveness check for the pairing between the spec corpus and CI's armed job.
//!
//! `ci.yml`'s `test-oracle` job runs the suite with `FERRO_ASSERT_IDEMPOTENT`,
//! `FERRO_ASSERT_REPARSE` and `FERRO_ASSERT_IN_BOUNDS` set. Those oracles
//! **panic** at the normalization seam. The spec-corpus modules **count** the
//! same defects as pinned figures, and a panic returns no value — so a row the
//! oracle fires on never reaches the family's output set, and the census reads
//! better than the truth.
//!
//! Measured on `main` at 35de96c8, both directions, all 11 affected families:
//! the non-idempotent output is in every case the one that disagreed with its
//! siblings, so dropping it makes the family look unanimous. Confluence reads
//! 9147 converged / 2428 split-2 armed, against a true 9140 / 2435.
//!
//! **The failure mode is the flattering kind.** A corpus module added later and
//! not named in `ORACLE_EXCLUDE` does not go red — it reports a *better* census,
//! which reads as progress rather than as the lost evidence it is. That is the
//! same shape as #1460 and #1478, where a corpus that could not build a thing
//! reported its absence as a zero.
//!
//! This is not a coverage exemption. Both modules run in full in the plain
//! `test` job, and the corpus measures idempotency itself — a count and a
//! classification, where the oracle offers only a panic.
//!
//! Modelled on `sweep_filter_invariant.rs`, which makes the same class of silent
//! config rot loud for the sweep-seed knob.

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
/// flattering direction — an unmatched module is never demanded in
/// `ORACLE_EXCLUDE`, so its census is taken with the seam oracles armed and reads
/// better than the truth — and it is the exact shape this file's module doc
/// predicts ("a corpus module added later and not named in ORACLE_EXCLUDE does not
/// go red, it reports a better census").
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
/// flattering direction: a module written that way is never demanded in
/// `ORACLE_EXCLUDE`, so its census is taken with the oracle armed and reads
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
    // the flattering direction — an unmatched module is never demanded in
    // `ORACLE_EXCLUDE`, so its census is taken with the seam oracles armed and
    // reads better than the truth:
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

/// A `>-` folded scalar from `ci.yml`, folded back to one line.
///
/// Reading only the first line would silently exempt whichever module happens to
/// sit on the second, which is exactly the blind spot this file exists to close.
fn ci_filter(key: &str) -> String {
    let ci = std::fs::read_to_string(repo_root().join(".github/workflows/ci.yml"))
        .expect("ci.yml is readable");
    let mut lines = ci.lines();
    let header = lines
        .by_ref()
        .find(|l| l.trim_start().starts_with(&format!("{key}:")))
        .unwrap_or_else(|| panic!("ci.yml defines {key}"));

    // An inline value is returned, not merely permitted: accepting a shape the
    // parser cannot read is its own small version of the rot this file catches.
    let inline = header
        .split_once(':')
        .map(|(_, value)| value.trim())
        .unwrap_or_default();
    if !inline.is_empty() && !inline.starts_with(['>', '|']) {
        return inline.to_string();
    }
    assert!(
        inline.starts_with(['>', '|']),
        "{key} is neither a block scalar nor an inline value: {header:?}"
    );

    let key_indent = header.len() - header.trim_start().len();
    let mut value = String::new();
    for line in lines {
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            break;
        }
        let indent = line.len() - line.trim_start().len();
        if indent <= key_indent {
            break;
        }
        value.push(' ');
        value.push_str(line.trim());
    }
    assert!(
        !value.trim().is_empty(),
        "{key} parsed as empty; ci.yml's formatting changed"
    );
    value
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

/// Every module that measures over the spec corpus must be named in
/// `ORACLE_EXCLUDE`, or its census is taken with an oracle armed and reads
/// better than the truth.
#[test]
fn every_spec_corpus_module_is_named_in_the_oracle_exclude() {
    let filter = ci_filter("ORACLE_EXCLUDE");
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
        "these modules measure over the spec corpus but are not named in \
         ci.yml's ORACLE_EXCLUDE, so `test-oracle` runs them with the seam \
         oracles armed. A panicking row contributes no output, which does not \
         redden the job — it makes confluence read HIGHER than it is: {missing:#?}\n\
         ORACLE_EXCLUDE is:{filter}\n\
         Add `+ test(<module>)` to it, or stop measuring over the corpus."
    );
}

/// The converse: a name in `ORACLE_EXCLUDE` that no longer measures over the
/// corpus is withholding a module from the armed job for no reason.
#[test]
fn every_module_named_in_the_oracle_exclude_measures_over_the_corpus() {
    let filter = ci_filter("ORACLE_EXCLUDE");
    let modules = corpus_modules();

    let named = modules_named_in(&filter);
    assert!(
        !named.is_empty(),
        "ORACLE_EXCLUDE names no modules; its formatting changed: {filter}"
    );

    let stale: Vec<&String> = named
        .iter()
        .filter(|module| !modules.contains(module))
        .collect();
    assert!(
        stale.is_empty(),
        "ci.yml's ORACLE_EXCLUDE names these modules, but they do not measure \
         over the spec corpus — so they are being withheld from the armed job \
         for no reason: {stale:#?}\n\
         Remove them from ORACLE_EXCLUDE."
    );
}

/// The scan's one structural blind spot, closed by forbidding the route rather
/// than by widening the matcher.
///
/// [`corpus_modules`] recognises a consumer by a literal import path from
/// [`CORPUS_MODULES`], which is right for avoiding prose false positives but
/// cannot see a module that reaches the corpus **indirectly** — through a shared
/// helper in `tests/it/common/` that does the importing. Such a module would
/// carry no matching literal, would not be demanded in `ORACLE_EXCLUDE`, and
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
         demand it be named in ORACLE_EXCLUDE — and it would then measure with \
         the seam oracles armed, which reports a better census rather than \
         going red.\n\
         Either import the corpus directly in each consuming module, or teach \
         `corpus_modules` to follow this helper."
    );
}

/// The local oracle runner, which mirrors `test-oracle` outside CI.
const LOCAL_RUNNER: &str = "scripts/run_oracle_suite.sh";

/// The `FERRO_ASSERT_*` keys the `test-oracle` job sets, derived in Rust.
///
/// Deliberately anchored **differently** from the awk in [`LOCAL_RUNNER`],
/// which scopes to the window between that job's step `name:` and its `run:`.
/// This one bounds the whole `test-oracle:` job by indentation and keys on the
/// quoted `: "1"` value. A second opinion that shares the other's derivation is
/// not a second opinion — it is the same reading written twice, and this repo
/// has already shipped a defect that way (`check_changelog_grouping.py`'s
/// rationale).
///
/// Comment lines are skipped for the reason the awk skips them: the job's
/// comment block *mentions* every flag in prose, including the history of why
/// `FERRO_ASSERT_SEQUENCE` used to be absent. A scan that read prose as a
/// setting would demand the runner arm oracles CI does not, and the rows that
/// prose names would then be red locally and green in CI.
///
/// **Scoped to the ARMED step, and that scoping is load-bearing as of #1815.**
/// That change gave `test-oracle` a *second* step — the compensating run that
/// re-executes `SEQUENCE_ORACLE_EXCLUDE`'s rows under the other three oracles —
/// which sets three of the same keys. A whole-job scan therefore returns seven
/// entries with three duplicated, against the awk's four, and this guard fails
/// on a correctly-wired file. (It did, when the step was first added; that is
/// how this scoping came to exist.)
///
/// The discriminator is deliberately **not** the step's `name:`, which is what
/// [`LOCAL_RUNNER`]'s awk anchors on: the armed step is the one that
/// `--partition`s the suite, and the compensating step is un-partitioned by
/// design. Keying on the `run:` body rather than on the label keeps the two
/// derivations independent — a renamed step breaks one of them and not the
/// other, which is the whole point of having two.
fn test_oracle_job_flags() -> Vec<String> {
    let armed = test_oracle_steps()
        .into_iter()
        .find(|step| step.runs.contains("--partition"))
        .expect(
            "ci.yml's test-oracle job has a step whose `run:` partitions the suite; \
             its shape changed and this guard would otherwise read the wrong step's flags",
        );
    armed.flags
}

/// One step of the `test-oracle:` job: the `FERRO_ASSERT_*` keys its `env:`
/// sets, and the text of its `run:`.
struct OracleStep {
    /// `FERRO_ASSERT_*` keys set to `"1"`, in file order, comments skipped.
    flags: Vec<String>,
    /// Every line of the step at or after its `run:`, joined.
    runs: String,
    /// The `-E` expression the step hands to nextest, if it passes one.
    selection: Option<String>,
}

/// The `test-oracle:` job's steps, split on the `- name:` at step indent.
///
/// Exists because #1815 made that job carry two nextest steps with overlapping
/// flag sets, so "the flags of `test-oracle`" stopped being a well-formed
/// question — every guard below has to say *which* step it means.
fn test_oracle_steps() -> Vec<OracleStep> {
    let mut steps: Vec<OracleStep> = Vec::new();
    let mut in_run = false;
    for line in test_oracle_job_lines() {
        let trimmed = line.trim();
        // A step boundary: `- name:` at the job's step indent.
        if line.starts_with("      - name:") {
            steps.push(OracleStep {
                flags: Vec::new(),
                runs: String::new(),
                selection: None,
            });
            in_run = false;
            continue;
        }
        let Some(step) = steps.last_mut() else {
            continue;
        };
        if trimmed.starts_with('#') {
            continue;
        }
        if trimmed == "run: |" || trimmed.starts_with("run:") {
            in_run = true;
            continue;
        }
        if in_run {
            step.runs.push('\n');
            step.runs.push_str(&line);
            if let Some(at) = line.find("-E \"") {
                let rest = &line[at + "-E \"".len()..];
                if let Some(end) = rest.find('"') {
                    step.selection = Some(rest[..end].to_string());
                }
            }
            continue;
        }
        if let Some(key) = trimmed
            .strip_suffix(": \"1\"")
            .filter(|key| key.starts_with("FERRO_ASSERT_"))
        {
            step.flags.push(key.to_string());
        }
    }
    steps
}

/// The lines of `ci.yml`'s `test-oracle:` job, bounded by indentation.
fn test_oracle_job_lines() -> Vec<String> {
    let ci = std::fs::read_to_string(repo_root().join(".github/workflows/ci.yml"))
        .expect("ci.yml is readable");
    let mut in_job = false;
    let mut lines = Vec::new();
    for line in ci.lines() {
        if line.starts_with("  test-oracle:") {
            in_job = true;
            continue;
        }
        if in_job {
            // A non-blank, non-comment line at job-key indent ends the job.
            let indent = line.len() - line.trim_start().len();
            if !line.trim().is_empty() && !line.trim_start().starts_with('#') && indent <= 2 {
                break;
            }
            lines.push(line.to_string());
        }
    }
    lines
}

/// The complete `-E` expression `test-oracle` hands to nextest, with `ci.yml`'s
/// two variable references expanded.
///
/// The job negates four things — the proptest modules, `SWEEP_FILTER`,
/// `ORACLE_EXCLUDE` and `CENSUS_FILTER` — and the runner first shipped negating
/// only the second-to-last of them, so a local "as CI runs it" also executed the
/// proptest modules and the three exhaustive sweeps. Comparing only the
/// exclusion could not see that, which is why the whole expression is compared
/// here.
fn ci_oracle_selection() -> String {
    // The ARMED step's selection, identified by `--partition` rather than by
    // position. #1815 gave this job a second nextest step, and "the first `-E` in
    // the job" is a positional accident rather than a statement about which step
    // is meant — reordering them would silently retarget this whole comparison.
    let template = test_oracle_steps()
        .into_iter()
        .find(|step| step.runs.contains("--partition"))
        .and_then(|step| step.selection)
        .expect("ci.yml's test-oracle job passes a -E selection to nextest in its armed step");

    let expanded = template
        .replace("$SWEEP_FILTER", ci_filter("SWEEP_FILTER").trim())
        // Before `$ORACLE_EXCLUDE` only for readability — the order does not
        // matter, because the pattern carries the leading `$` and
        // `$SEQUENCE_ORACLE_EXCLUDE` holds no second one. Measured both ways.
        // `scripts/run_oracle_suite.sh` carries the same note; the resemblance
        // between the two names is the trap, and it is not a real one.
        .replace(
            "$SEQUENCE_ORACLE_EXCLUDE",
            ci_filter("SEQUENCE_ORACLE_EXCLUDE").trim(),
        )
        .replace("$ORACLE_EXCLUDE", ci_filter("ORACLE_EXCLUDE").trim())
        .replace("$CENSUS_FILTER", ci_filter("CENSUS_FILTER").trim());
    assert!(
        !expanded.contains('$'),
        "test-oracle's -E selection references a variable this test does not expand: {expanded}"
    );
    expanded
}

/// What [`LOCAL_RUNNER`] says it would run.
struct RunnerSelection {
    /// The `ORACLE_EXCLUDE` value it extracted from `ci.yml`.
    exclude: String,
    /// The complete `-E` expression it would hand to nextest.
    selection: String,
    /// The `FERRO_ASSERT_*` flags it would arm.
    flags: Vec<String>,
}

/// What [`LOCAL_RUNNER`] says it would run, read from its own `--print-selection`.
fn local_runner_selection() -> RunnerSelection {
    let script = repo_root().join(LOCAL_RUNNER);
    assert!(
        script.is_file(),
        "{LOCAL_RUNNER} is missing. It is the only local command that runs the seam \
         oracles over the same selection CI does; without it the documented bare \
         `FERRO_ASSERT_IDEMPOTENT=1 cargo nextest run` is the only recipe, and that \
         one cannot pass."
    );
    let output = std::process::Command::new("bash")
        .arg(&script)
        .arg("--print-selection")
        .current_dir(repo_root())
        .output()
        .expect("bash can run the local oracle runner");
    assert!(
        output.status.success(),
        "{LOCAL_RUNNER} --print-selection failed ({}):\n{}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).expect("the runner prints UTF-8");

    let mut exclude = String::new();
    let mut selection = String::new();
    let mut flags = Vec::new();
    for line in stdout.lines() {
        if let Some(value) = line.strip_prefix("ORACLE_EXCLUDE=") {
            exclude = value.trim().to_string();
        } else if let Some(value) = line.strip_prefix("SELECTION=") {
            selection = value.trim().to_string();
        } else if let Some(value) = line.strip_prefix("FLAG=") {
            flags.push(value.trim().to_string());
        }
    }
    assert!(
        !selection.is_empty(),
        "{LOCAL_RUNNER} --print-selection printed no SELECTION= line; \
         its output shape changed and the comparison below would go vacuous"
    );
    RunnerSelection {
        exclude,
        selection,
        flags,
    }
}

/// The local runner must exclude exactly what `test-oracle` excludes.
///
/// A drift here fails in the **flattering** direction, which is why it is worth
/// a test rather than a comment: a runner excluding more than CI goes green
/// locally on a defect CI is red on, and the operator's evidence for "my change
/// is clean" is then a run that never touched the modules in question.
#[test]
fn the_local_oracle_runner_excludes_exactly_what_ci_excludes() {
    assert_eq!(
        local_runner_selection().exclude,
        ci_filter("ORACLE_EXCLUDE").trim(),
        "{LOCAL_RUNNER} and ci.yml disagree about ORACLE_EXCLUDE. The runner reads it \
         from ci.yml with awk; if that extraction broke, fix the awk — do not inline a \
         copy of the list, which is the drift this test exists to catch."
    );
}

/// The local runner must run the **whole** selection `test-oracle` runs, not
/// just its exclusion.
///
/// `the_local_oracle_runner_excludes_exactly_what_ci_excludes` compares one of
/// the expression's three negated terms, so it passed throughout the period the
/// runner also executed the proptest modules and the three exhaustive sweeps —
/// tests `test-oracle` does not run. A guard over a subset of the contract reads
/// as a guard over the contract, which is the failure this file is otherwise
/// about.
#[test]
fn the_local_oracle_runner_selects_exactly_what_ci_selects() {
    assert_eq!(
        local_runner_selection().selection,
        ci_oracle_selection(),
        "{LOCAL_RUNNER} and ci.yml's test-oracle job disagree about the nextest selection. \
         The runner reads the whole `-E` expression out of ci.yml and expands its variable \
         references; if that extraction broke, fix the awk — do not inline a copy of the \
         expression, which is the drift this test exists to catch."
    );
}

/// …and the agreed selection must actually **negate** all three things CI
/// negates.
///
/// The equality test above cannot tell a correct expression from a matching
/// pair of wrong ones: both sides read the same line of `ci.yml`, so a `-E`
/// that stopped negating the sweeps would satisfy it. This one keys on the
/// `not (…)` wrapper instead.
///
/// **Asserting that each module is merely NAMED in the selection would be
/// vacuous**, and that was this test's first form. The selection is built *by
/// substituting those two filters' values into* the template, so every module
/// they name is guaranteed to appear as a `test(<module>)` substring however
/// the template is spelled — including a template that had dropped the `not`.
/// Dropping it is not hypothetical: `and ($SWEEP_FILTER)` in place of
/// `and not ($SWEEP_FILTER)` inverts the job from "the suite minus the three
/// sweeps" to "only the three sweeps", i.e. ~8,900 tests down to ~30, which is
/// the quiet-narrowing failure the runner's own header warns about.
/// **Its name may not contain `proptest`, and that is not cosmetic.** `test()`
/// is a substring predicate over the whole test name, so a function carrying
/// that token is selected by the `soak` job's `-E 'test(proptest)'` and negated
/// by `test` and `test-oracle` — which used to mean this assertion about
/// `ci.yml` ran only inside the 1M-case soak, and now that the soak archive
/// holds only the modules `tests-soak/tests/soak/main.rs` compiles, would mean
/// it ran nowhere at all. `tests/it/soak_package_membership.rs` fails on any
/// such name.
#[test]
fn the_ci_oracle_selection_negates_the_property_tests_the_sweeps_and_the_corpus_modules() {
    let selection = ci_oracle_selection();
    assert!(
        selection.contains("not test(proptest)"),
        "test-oracle's selection no longer negates the proptest modules, which the `soak` \
         job owns at 125k cases per shard: {selection}"
    );

    for key in ["SWEEP_FILTER", "ORACLE_EXCLUDE", "SEQUENCE_ORACLE_EXCLUDE"] {
        let value = ci_filter(key);
        let value = value.trim();
        assert!(
            !modules_named_in(value).is_empty(),
            "ci.yml's {key} names no module; its formatting changed and this guard has \
             gone vacuous"
        );
        assert!(
            selection.contains(&format!("not ({value})")),
            "test-oracle's -E selection does not negate {key}, so the armed job RUNS those \
             modules instead of withholding them. Note a bare `and (${key})` would still \
             mention every one of them, which is why this asserts the `not (…)` wrapper \
             rather than the module names.\nSelection is: {selection}"
        );
    }
}

/// …and arm exactly the oracles `test-oracle` arms.
///
/// Both directions matter and neither is symmetric with the other. Arming
/// **fewer** makes a local run weaker than CI while reading as an oracle pass.
/// Arming **more** makes the runner red on rows no PR caused, which teaches the
/// operator to ignore it.
///
/// `FERRO_ASSERT_SEQUENCE` was the live candidate for the second failure until
/// #1815, and its history is the argument for keeping this guard: armed over this
/// selection it read red at 5 tests (`674e9c8b`), 2 (`c9207d7e`, after #1990
/// closed #1690) and 3 (`1aecc93a`, after #2051 added a gate that fires) — so a
/// runner that had armed it ahead of the job would have been red on every PR for
/// months, at a count that moved in both directions. It is armed in both places
/// now, and the runner did not have to be changed to follow: it reads the flag
/// set out of the job.
///
/// The candidate has not gone away, it has moved. `censuses`' armed step now
/// arms three where this job arms four, deliberately and unmeasured — see that
/// job's header — so the next plausible "restore parity" edit is there.
#[test]
fn the_local_oracle_runner_arms_exactly_the_flags_test_oracle_arms() {
    let runner_flags = local_runner_selection().flags;
    let ci_flags = test_oracle_job_flags();
    assert!(
        !ci_flags.is_empty(),
        "no FERRO_ASSERT_* flag was found in ci.yml's test-oracle job; its formatting \
         changed and this guard has gone vacuous"
    );
    assert_eq!(
        runner_flags, ci_flags,
        "{LOCAL_RUNNER} and ci.yml's test-oracle job disagree about which oracles to arm"
    );
}

/// The runner may not carry its own copy of the exclusion list.
///
/// The equality tests above compare *values*, so they would still pass against a
/// hardcoded copy that happens to be current today. This one forbids the copy
/// itself, because the value only has to be right at the moment someone runs the
/// test — and the failure mode of a stale copy is silent and flattering.
#[test]
fn the_local_oracle_runner_does_not_hardcode_the_exclusion() {
    let text = std::fs::read_to_string(repo_root().join(LOCAL_RUNNER))
        .expect("the local oracle runner is readable");
    // Strip comments before scanning: the script's header explains the failure
    // by NAMING the modules, which is exactly what makes the header useful and
    // must not be mistaken for a hardcoded filter.
    let code: String = text
        .lines()
        .filter(|line| !line.trim_start().starts_with('#'))
        .collect::<Vec<_>>()
        .join("\n");
    let named = modules_named_in(&ci_filter("ORACLE_EXCLUDE"));
    // The same non-vacuity guard its siblings carry: with no module names to
    // look for, the filter below is trivially empty and this test passes
    // without having checked anything.
    assert!(
        !named.is_empty(),
        "ORACLE_EXCLUDE names no modules; its formatting changed and this guard has gone \
         vacuous"
    );
    let hardcoded: Vec<String> = named
        .into_iter()
        .filter(|module| code.contains(&format!("test({module})")))
        .collect();
    assert!(
        hardcoded.is_empty(),
        "{LOCAL_RUNNER} hardcodes these modules instead of reading ORACLE_EXCLUDE from \
         ci.yml: {hardcoded:#?}\n\
         Two copies of one list drift, and this one drifts flatteringly — a stale copy \
         excludes a module CI still runs armed."
    );
}

/// Every row withheld from the denoted-sequence oracle must name something that
/// exists.
///
/// `test()` is a substring predicate, so a typo does not error — it selects
/// nothing. In `SEQUENCE_ORACLE_EXCLUDE` that fails in the **loud** direction
/// (the armed job stops withholding the row and goes red on it), which is the
/// good case and is why this is a cheap guard rather than a critical one. The
/// expensive half is the compensating step: `--no-tests=fail` turns a filter that
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
    let named = modules_named_in(&ci_filter("SEQUENCE_ORACLE_EXCLUDE"));
    assert!(
        !named.is_empty(),
        "ci.yml's SEQUENCE_ORACLE_EXCLUDE names nothing; either its formatting changed \
         (in which case fix it) or its last row retired — in which case delete the \
         variable, the `and not (…)` term in test-oracle's -E, the compensating step, \
         the block in scripts/run_oracle_suite.sh, and these guards, in one change."
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
        "SEQUENCE_ORACLE_EXCLUDE names these, and nothing in tests/ or examples/ \
         defines them: {missing:#?}\n\
         A `test()` term that matches nothing selects nothing, so the compensating \
         step would silently stop re-running the row it is there to protect."
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
/// This is the guard that makes arming the fourth oracle a **superset** of what
/// `test-oracle` ran before #1815 rather than a trade. A nextest `-E` is one
/// expression, so the armed step's `and not ($SEQUENCE_ORACLE_EXCLUDE)` withdraws
/// those rows from all four oracles at once; the compensating step is what puts
/// three of them back. Delete that step and the change quietly becomes "three
/// oracles surrendered to gain one" — with nothing red, which is why it is
/// asserted rather than left to the comment beside it.
///
/// Three things are checked, and the third is the one a reader would omit:
///
/// 1. the step exists and selects **exactly** `$SEQUENCE_ORACLE_EXCLUDE`, by
///    reference and not by a copy of its value;
/// 2. it arms exactly the other three oracles;
/// 3. it does **not** arm `FERRO_ASSERT_SEQUENCE`. Without this, a copy-paste of
///    the armed step's `env:` block would satisfy (2) as a subset check and make
///    the step red on the very rows it exists to run — turning a green job red
///    for a reason that looks like a real defect.
#[test]
fn the_sequence_oracle_exclusions_still_run_under_the_other_three_oracles() {
    let steps = test_oracle_steps();
    let armed = steps
        .iter()
        .find(|step| step.runs.contains("--partition"))
        .expect("test-oracle has a partitioned, armed step");
    assert!(
        armed.flags.iter().any(|f| f == "FERRO_ASSERT_SEQUENCE"),
        "test-oracle's armed step no longer sets FERRO_ASSERT_SEQUENCE, so there is \
         nothing for SEQUENCE_ORACLE_EXCLUDE to withhold it from. If the flag is being \
         un-armed, remove the filter and this guard in the same change rather than \
         leaving a compensating step for an exclusion that excludes nothing."
    );

    let compensating: Vec<&OracleStep> = steps
        .iter()
        .filter(|step| {
            step.selection.as_deref() == Some("$SEQUENCE_ORACLE_EXCLUDE")
                && !step.runs.contains("--partition")
        })
        .collect();
    assert_eq!(
        compensating.len(),
        1,
        "expected exactly one un-partitioned step in test-oracle selecting \
         `$SEQUENCE_ORACLE_EXCLUDE`, found {}.\n\
         Without it, the armed step's `and not ($SEQUENCE_ORACLE_EXCLUDE)` withdraws \
         those rows from FERRO_ASSERT_IDEMPOTENT, _REPARSE and _IN_BOUNDS as well — \
         which they pass — so arming the fourth oracle would cost three.\n\
         Note the selection must be the VARIABLE REFERENCE, not a copy of its value: \
         two copies of one list drift, and this one drifts flatteringly.",
        compensating.len()
    );
    let compensating = compensating[0];

    let mut expected: Vec<String> = armed
        .flags
        .iter()
        .filter(|f| *f != "FERRO_ASSERT_SEQUENCE")
        .cloned()
        .collect();
    let mut actual = compensating.flags.clone();
    expected.sort();
    actual.sort();
    assert_eq!(
        actual, expected,
        "the compensating step must arm exactly the oracles the armed step arms MINUS \
         FERRO_ASSERT_SEQUENCE.\n\
         Armed step: {:?}\nCompensating step: {:?}\n\
         Deriving the expectation from the armed step rather than from a literal list is \
         deliberate — a fifth oracle added above must reach these rows too, and a hardcoded \
         three would not notice.",
        armed.flags, compensating.flags
    );
    assert!(
        !compensating
            .flags
            .iter()
            .any(|f| f == "FERRO_ASSERT_SEQUENCE"),
        "the compensating step arms FERRO_ASSERT_SEQUENCE, which is the one oracle these \
         rows are withheld from. It would fire on exactly the rows the step exists to run, \
         reddening a job that is otherwise green — and the failure would read as a fresh \
         normalizer defect rather than as this wiring mistake."
    );
    assert!(
        compensating.runs.contains("--no-tests=fail"),
        "the compensating step must pass --no-tests=fail. Its selection is five tests \
         named by substring, so a rename that empties the filter would otherwise be a \
         step that goes green having run nothing — the vacuous pass this repository \
         keeps meeting."
    );
}

/// The debt list must stay disjoint from the two permanent filters.
///
/// An overlap is not merely untidy here, unlike the sweep/oracle pair below. A
/// row in both `ORACLE_EXCLUDE` and `SEQUENCE_ORACLE_EXCLUDE` is withheld from
/// the armed step twice — harmless — but the compensating step would then run it
/// under three oracles that the armed job deliberately never runs it under,
/// because `ORACLE_EXCLUDE`'s whole point is that those instruments destroy each
/// other on those modules. So the overlap would *create* the red that
/// `ORACLE_EXCLUDE` exists to prevent, in a new step, for a reason nothing states.
///
/// The `SWEEP_FILTER` half is the milder version: those rows run in `sweeps` with
/// all four oracles armed already, so re-running them here would be redundant
/// rather than wrong — but it would also mean the debt list silently governs a
/// job it says nothing about.
#[test]
fn the_sequence_oracle_exclude_is_disjoint_from_the_permanent_filters() {
    let debt = modules_named_in(&ci_filter("SEQUENCE_ORACLE_EXCLUDE"));
    assert!(
        !debt.is_empty(),
        "SEQUENCE_ORACLE_EXCLUDE names nothing; this guard has gone vacuous"
    );
    for key in ["ORACLE_EXCLUDE", "SWEEP_FILTER", "CENSUS_FILTER"] {
        let other = modules_named_in(&ci_filter(key));
        assert!(
            !other.is_empty(),
            "ci.yml's {key} names no module; this guard has gone vacuous"
        );
        let both: Vec<&String> = debt.iter().filter(|m| other.contains(m)).collect();
        assert!(
            both.is_empty(),
            "these are named in BOTH SEQUENCE_ORACLE_EXCLUDE and {key}: {both:#?}\n\
             The debt list is temporary and carries issue numbers; {key} is a standing \
             statement about what the armed job must never run. A row in both means the \
             compensating step re-runs, under three armed oracles, a module {key} \
             withholds from them."
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
    let sweeps = modules_named_in(&ci_filter("SWEEP_FILTER"));
    let excluded = modules_named_in(&ci_filter("ORACLE_EXCLUDE"));

    let both: Vec<&String> = sweeps.iter().filter(|m| excluded.contains(m)).collect();
    assert!(
        both.is_empty(),
        "these modules are named in BOTH ci.yml filters, which are meant to \
         select disjoint sets for different reasons: {both:#?}"
    );
}

/// [`imports_corpus`] recognises both spellings, and neither a prose mention nor
/// a neighbouring module whose name merely starts the same way.
///
/// Written because the grouped form was the gap: the scan matched a single
/// literal, so `use ferro_hgvs::conformance::{spec_corpus, summary};` read as
/// "does not consume the corpus" and the module was never demanded in
/// `ORACLE_EXCLUDE`.
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
