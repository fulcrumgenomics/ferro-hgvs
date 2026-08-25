//! Liveness check for the pairing between the optimized archive and the jobs
//! that run off it.
//!
//! `test-build-soak` no longer archives the whole `it` binary. It archives one
//! test target owned by the `ferro-hgvs-soak-tests` package, whose driver
//! (`tests-soak/tests/soak/main.rs`) `#[path]`-includes exactly the modules the
//! three consuming jobs — `soak`, `sweeps`, `censuses` — select. That buys the
//! compile-time saving the split exists for, and it introduces one list that can
//! drift: the driver's.
//!
//! **The drift fails silently and in the flattering direction.** A module named
//! in `SWEEP_FILTER` or `CENSUS_FILTER`, or matching the `soak` job's
//! `test(proptest)`, is NEGATED by `test` and `test-oracle`. If it is not also
//! compiled into the soak driver, it is in no archive that any job selects, so
//! it runs **nowhere** — and nothing goes red. `--no-tests=fail` on each of the
//! three jobs catches a selection that matches *nothing*, never one module of
//! four going missing. This is the same shape as `sweep_filter_invariant`'s
//! unnamed-consumer case and as `census_filter_invariant`'s unrun-module case,
//! one level further out: those two ask whether the FILTERS agree with the jobs,
//! and this one asks whether the BINARY agrees with the filters.
//!
//! The converse is cheaper but still worth failing on: a module compiled into
//! the driver that no job selects is paid for on the critical path — the exact
//! cost the split was made to remove — while contributing nothing. With one
//! earned exception, which the driver has to demonstrate rather than assert: a
//! module another compiled module names as `crate::<module>::` is a link-time
//! dependency, not dead weight. See
//! [`the_soak_driver_compiles_nothing_the_archive_jobs_do_not_select`].
//!
//! # Why this reads `ci.yml` itself rather than sharing a helper
//!
//! The fourth file in `tests/it/` to parse `ci.yml`, and deliberately so, for
//! the reason `census_filter_invariant.rs` sets out at length: these are
//! liveness guards, a bug in one shared extractor would neuter all of them at
//! once, and it would do so by returning an empty selection — which makes every
//! assertion vacuous rather than red. Each guard asserts its own read is
//! non-empty, so four independent readers are four independent chances to
//! notice.
//!
//! # The one thing it does not see
//!
//! `test(X)` is a substring predicate over the whole test NAME, so it can match
//! a test function rather than a module — `test(proptest)` would select a
//! `fn …_proptest()` sitting in a module this driver does not compile. The scan
//! below expands each term over module names only. There is no such function
//! today (asserted by [`no_test_function_hides_behind_a_module_level_term`]),
//! which is what keeps the module-level reading complete.

use std::collections::BTreeSet;
use std::path::PathBuf;

/// The workflow this file is entirely about.
const CI_YML: &str = ".github/workflows/ci.yml";

/// The soak driver, whose `#[path]` list is the membership under test.
const SOAK_DRIVER: &str = "tests-soak/tests/soak/main.rs";

/// The jobs that run off the archive `test-build-soak` produces. Every module
/// any of them selects has to be compiled into it. `censuses-plain` is the
/// un-armed half of the census job, split out so it runs in parallel with the
/// armed half; both consume the soak archive.
const ARCHIVE_CONSUMERS: &[&str] = &["soak", "sweeps", "censuses", "censuses-plain"];

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn ci_text() -> String {
    std::fs::read_to_string(repo_root().join(CI_YML)).expect("ci.yml is readable")
}

/// A `>-` folded scalar from `ci.yml`'s top-level `env:`, folded back to one
/// line, or `None` when no such key exists.
fn ci_filter(key: &str) -> Option<String> {
    let ci = ci_text();
    let mut lines = ci.lines();
    let header = lines
        .by_ref()
        .find(|line| line.trim_start().starts_with(&format!("{key}:")))?;

    let inline = header
        .split_once(':')
        .map(|(_, value)| value.trim())
        .unwrap_or_default();
    if !inline.is_empty() && !inline.starts_with(['>', '|']) {
        return Some(inline.to_string());
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
    Some(value)
}

/// Every `-E "…"` selection a top-level job hands to nextest, in file order,
/// with each `$VAR` reference expanded from `ci.yml`'s own `env:` block.
///
/// Expanding rather than reading the variables separately is what makes the
/// scan agree with what the runner actually selects: `sweeps` passes
/// `"$SWEEP_FILTER + test(issue_1615_denoted_sequence_oracle)"`, and only the
/// expansion carries both halves.
fn selections_in(job: &str) -> Vec<String> {
    let ci = ci_text();
    let header = format!("  {job}:");
    let mut lines = ci.lines().skip_while(|line| *line != header);
    assert!(
        lines.next().is_some(),
        "ci.yml defines no `{job}` job; this guard cannot see what it runs"
    );
    let body: Vec<&str> = lines
        // A non-blank, non-comment line back at job-key indent ends the block.
        .take_while(|line| {
            line.starts_with("    ") || line.trim().is_empty() || line.trim_start().starts_with('#')
        })
        .collect();

    body.iter()
        .filter(|line| !line.trim_start().starts_with('#'))
        .filter_map(|line| {
            let at = line.find("-E \"").or_else(|| line.find("-E '"))?;
            let quote = line.as_bytes()[at + 3] as char;
            let rest = &line[at + 4..];
            rest.find(quote).map(|end| expand_vars(&rest[..end]))
        })
        .collect()
}

/// Substitute every `$NAME` in a selection with that key's value from `ci.yml`'s
/// `env:` block.
fn expand_vars(selection: &str) -> String {
    let mut out = selection.to_string();
    while let Some(at) = out.find('$') {
        let name: String = out[at + 1..]
            .chars()
            .take_while(|c| c.is_ascii_uppercase() || *c == '_')
            .collect();
        assert!(
            !name.is_empty(),
            "a bare `$` in a ci.yml selection: {selection}"
        );
        let value = ci_filter(&name).unwrap_or_else(|| {
            panic!("selection references ${name}, which ci.yml does not define")
        });
        out.replace_range(at..at + 1 + name.len(), &value);
    }
    out
}

/// The `test(...)` terms of a nextest filter expression.
fn terms_in(filter: &str) -> Vec<String> {
    filter
        .match_indices("test(")
        .filter_map(|(at, _)| {
            let rest = &filter[at + "test(".len()..];
            rest.find(')').map(|end| rest[..end].trim().to_string())
        })
        .collect()
}

/// Every top-level module of `tests/it/`, by name.
fn it_modules() -> BTreeSet<String> {
    std::fs::read_dir(repo_root().join("tests/it"))
        .expect("tests/it is readable")
        .filter_map(|entry| {
            let path = entry.expect("readable dir entry").path();
            let name = path.file_name()?.to_str()?;
            name.strip_suffix(".rs").map(str::to_string)
        })
        .collect()
}

/// The modules the three archive-consuming jobs select, by expanding each
/// `test(X)` term over the module names it matches.
///
/// `test(X)` is a substring predicate, which is why `test(proptest)` — the whole
/// of the `soak` job's selection — resolves to three modules without being
/// special-cased here.
fn modules_selected_from_the_archive() -> BTreeSet<String> {
    let modules = it_modules();
    let mut selected = BTreeSet::new();
    for job in ARCHIVE_CONSUMERS {
        let selections = selections_in(job);
        assert!(
            !selections.is_empty(),
            "the `{job}` job passes no -E selection to nextest; its shape changed \
             and this guard cannot see what it runs"
        );
        for selection in selections {
            for term in terms_in(&selection) {
                let matched: Vec<&String> = modules.iter().filter(|m| m.contains(&term)).collect();
                assert!(
                    !matched.is_empty(),
                    "the `{job}` job selects `test({term})`, which matches no module \
                     under tests/it/ — a rename, or a filter this guard can no \
                     longer read"
                );
                selected.extend(matched.into_iter().cloned());
            }
        }
    }
    selected
}

/// The modules `tests-soak/tests/soak/main.rs` compiles, read from its `#[path]`
/// attributes.
///
/// The helpers under `tests/it/common/` are excluded: they are infrastructure
/// the modules reach through, not things a job selects.
fn soak_driver_modules() -> BTreeSet<String> {
    let text = std::fs::read_to_string(repo_root().join(SOAK_DRIVER))
        .unwrap_or_else(|e| panic!("read {SOAK_DRIVER}: {e}"));
    let needle = "#[path = \"../../../tests/it/";
    let mut found = BTreeSet::new();
    let mut rest = text.as_str();
    while let Some(at) = rest.find(needle) {
        let after = &rest[at + needle.len()..];
        let end = after
            .find('"')
            .expect("a #[path] attribute closes its string");
        let relative = &after[..end];
        rest = &after[end..];
        if let Some(name) = relative.strip_suffix(".rs") {
            if !name.contains('/') {
                found.insert(name.to_string());
            }
        }
    }
    assert!(
        !found.is_empty(),
        "{SOAK_DRIVER} declares no tests/it module; its `#[path]` form changed and \
         every assertion in this file would be vacuous"
    );
    found
}

/// Every module the three jobs select must be compiled into the archive they
/// select it from.
#[test]
fn every_module_the_archive_jobs_select_is_compiled_into_the_soak_driver() {
    let selected = modules_selected_from_the_archive();
    let compiled = soak_driver_modules();

    let missing: Vec<&String> = selected.difference(&compiled).collect();
    assert!(
        missing.is_empty(),
        "ci.yml's `soak` / `sweeps` / `censuses` jobs select these modules from the \
         optimized archive, but {SOAK_DRIVER} does not compile them into it: \
         {missing:?}\n\
         `test` and `test-oracle` negate every one of them, so they run in NO job \
         — CI stays green and gets faster, which is the coverage deletion this \
         guard exists to make loud.\n\
         Add a `#[path = \"../../../tests/it/<module>.rs\"] mod <module>;` to that \
         file."
    );
}

/// …and nothing else it cannot account for, or the critical path pays to compile
/// a module no job runs.
///
/// **A module no job selects is allowed exactly one justification: something
/// else the driver compiles reaches into it.** `minimal_alignment_enumeration_proptest`
/// cross-checks the enumerator against
/// `crate::issue_1539_split_member_separation::forced_unchanged_columns`, so that
/// module has to be in the binary for the selected one to link — while nothing
/// selects its own ~40 tests, which go on running from `it`.
///
/// The exemption is **derived, not declared**. A marker comment in the driver
/// would be a claim this test then trusts; requiring a `crate::<module>::`
/// reference from a compiled sibling makes the driver *demonstrate* the need,
/// and makes the exemption lapse by itself the moment the reference is deleted.
/// That is what stops "support module" becoming a place to park dead compile
/// cost on the workflow's critical path, which is the cost this package was
/// split out to remove.
#[test]
fn the_soak_driver_compiles_nothing_the_archive_jobs_do_not_select() {
    let selected = modules_selected_from_the_archive();
    let compiled = soak_driver_modules();

    // What every compiled module's source says, so a reference can be looked for
    // across the whole binary rather than only in the selected half — a support
    // module may legitimately be reached by another support module.
    let sources: Vec<(String, String)> = compiled
        .iter()
        .map(|module| {
            let path = repo_root().join(format!("tests/it/{module}.rs"));
            let text = std::fs::read_to_string(&path)
                .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
            (module.clone(), text)
        })
        .collect();

    let mut unaccounted = Vec::new();
    for module in compiled.difference(&selected) {
        let needle = format!("crate::{module}::");
        let reached_by: Vec<&str> = sources
            .iter()
            .filter(|(name, text)| name != module && text.contains(&needle))
            .map(|(name, _)| name.as_str())
            .collect();
        if reached_by.is_empty() {
            unaccounted.push(module.clone());
        }
    }

    assert!(
        unaccounted.is_empty(),
        "{SOAK_DRIVER} compiles these modules, no `soak` / `sweeps` / `censuses` \
         selection reaches them, and no module it compiles names them as \
         `crate::<module>::` either: {unaccounted:?}\n\
         They are pure compile cost on the workflow's critical path, which is the \
         cost this package was split out to remove. Drop them, or select them."
    );
}

/// The scan reads `test(X)` over module names, so a test FUNCTION whose name
/// carries one of those terms would be selected by a job and invisible here.
///
/// Closed by forbidding the shape rather than by widening the matcher: a
/// function named for one of these terms outside the modules the driver
/// compiles would be negated by `test` / `test-oracle` and absent from the
/// archive, i.e. run nowhere, which is the same silent loss the two tests above
/// exist for.
#[test]
fn no_test_function_hides_behind_a_module_level_term() {
    let terms: BTreeSet<String> = ARCHIVE_CONSUMERS
        .iter()
        .flat_map(|job| selections_in(job))
        .flat_map(|selection| terms_in(&selection))
        .collect();
    assert!(!terms.is_empty(), "no job selects on a test() term");

    let compiled = soak_driver_modules();
    let mut offenders = Vec::new();
    for module in it_modules() {
        if compiled.contains(&module) {
            continue;
        }
        let path = repo_root().join(format!("tests/it/{module}.rs"));
        let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {module}: {e}"));
        for line in text.lines() {
            let Some(at) = line.find("fn ") else { continue };
            if !line.trim_start().starts_with("fn ") && !line.trim_start().starts_with("pub fn ") {
                continue;
            }
            let name: String = line[at + 3..]
                .chars()
                .take_while(|c| c.is_alphanumeric() || *c == '_')
                .collect();
            for term in &terms {
                if name.contains(term.as_str()) {
                    offenders.push(format!("{module}::{name} matches test({term})"));
                }
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "these functions live outside the soak driver but carry a name one of the \
         archive jobs selects on, so they would be negated by `test` / \
         `test-oracle` and absent from the archive — run nowhere: {offenders:#?}\n\
         Rename them, or compile their module into {SOAK_DRIVER}."
    );
}
