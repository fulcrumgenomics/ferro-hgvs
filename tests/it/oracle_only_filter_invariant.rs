//! Liveness check for `ORACLE_ONLY_FILTER` — the filter negated from `test` and
//! from `test` alone.
//!
//! The claim it encodes is narrow and load-bearing: these tests do not need an
//! un-armed copy, **because they still run armed in `test-oracle`**. The whole
//! saving rests on that second clause, and nothing else in the repository
//! checks it. Negating the filter from `test-oracle` as well is a one-token
//! edit that reads like tidying, turns both shards green, makes them *faster*,
//! and deletes the census — the failure direction this repository keeps hitting
//! (see `tests/it/soak_package_membership.rs` and the bulk-fixture skip-green
//! note in `CONTRIBUTING.md`, "Bulk corpora: a skip that reads as a pass").
//!
//! So the assertions are, in order of what they defend:
//!
//! 1. The filter parses non-empty. Every other assertion here is vacuous
//!    against an empty filter, and an empty one would ALSO make `test`'s
//!    `and not ()` select everything — silently undoing the change while
//!    looking correct.
//! 2. `test` negates it. Otherwise nothing was saved.
//! 3. `test-oracle` does **not** negate it, and does not exclude it by any
//!    other route either. This is the coverage guarantee.
//! 4. It is not also named in `CENSUS_FILTER`. That combination is the specific
//!    trap the `ci.yml` comment describes: `CENSUS_FILTER` is negated from both
//!    shards and selected off the soak archive, which contains no `src/` unit
//!    tests, so a unit test named there would run nowhere at all.
//! 5. Every term it names still resolves to modules that exist. A rename fails
//!    safe — the filter matches nothing, `test` negates nothing, and the
//!    duplicate quietly returns to the debug shard — but that is the same
//!    silent config rot the file exists for, and it is what the `WHY NOT
//!    CENSUS_FILTER` paragraph is guarding against elsewhere.
//!
//! # Both assertions 3 and 4 compare terms in BOTH directions, deliberately
//!
//! `test(P)` is a **substring predicate on the test's name**: nextest matches
//! when `P` is contained in the name, not the reverse. So asking only whether
//! the foreign filter's *text* contains this filter's *module path* — the
//! obvious spelling, and the one this file shipped with — fires only on a
//! verbatim duplication of the full path, and misses the natural spelling
//! entirely.
//!
//! Concretely: `test(direction_symmetry)` added to `ORACLE_EXCLUDE` is the
//! spelling every entry there already uses. nextest then excludes
//! `normalize::merge::tests::direction_symmetry::…` from `test-oracle`, `test`
//! already negates it here, and the census runs **nowhere** — while
//! `"…test(direction_symmetry)…".contains("merge::tests::direction_symmetry")`
//! is `false`, so a one-directional guard passes. That is precisely the failure
//! this file was written to catch.
//!
//! Hence: extract the `test(...)` terms from **both** sides and require that
//! neither contains the other. The wider term wins whichever way round it is
//! written, so both directions are needed. `soak_package_membership.rs` uses
//! the same containment direction when it expands a term over module names.
//!
//! This file reads `ci.yml` as text rather than parsing YAML, for the reason
//! `census_filter_invariant.rs` gives: the repository has no YAML dependency in
//! its test tree, and the helpers below are shaped to that file's folded
//! scalars. It duplicates those helpers deliberately — two independent readers
//! of the same file are two chances to notice a formatting change, and sharing
//! them would mean a single parse bug silently disarming both guards.

use std::path::PathBuf;

/// The workflow this file is entirely about.
const CI_YML: &str = ".github/workflows/ci.yml";

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn ci_text() -> String {
    std::fs::read_to_string(repo_root().join(CI_YML)).expect("ci.yml is readable")
}

/// A `>-` folded scalar from `ci.yml`'s top-level `env:`, folded back to one line.
fn ci_filter(key: &str) -> String {
    let ci = ci_text();
    let mut lines = ci.lines();
    let header = lines
        .by_ref()
        .find(|line| line.trim_start().starts_with(&format!("{key}:")))
        .unwrap_or_else(|| panic!("ci.yml defines {key}"));

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

/// The body lines of one top-level job block, by its YAML key.
fn job_lines(job: &str) -> Vec<String> {
    let ci = ci_text();
    let header = format!("  {job}:");
    let mut lines = ci.lines().skip_while(|line| *line != header);
    assert!(
        lines.next().is_some(),
        "ci.yml defines no `{job}` job; this guard cannot see what it runs"
    );
    lines
        .take_while(|line| {
            line.starts_with("    ") || line.trim().is_empty() || line.trim_start().starts_with('#')
        })
        .map(str::to_string)
        .collect()
}

/// Every `-E "…"` selection a job hands to nextest, in file order.
fn selections_in(job: &str) -> Vec<String> {
    job_lines(job)
        .iter()
        .filter(|line| !line.trim_start().starts_with('#'))
        .filter_map(|line| {
            let at = line.find("-E \"")? + "-E \"".len();
            let rest = &line[at..];
            rest.find('"').map(|end| rest[..end].to_string())
        })
        .collect()
}

/// The `test(...)` terms of a nextest filter expression, in file order.
///
/// Named `terms_in` rather than `modules_named_in` (its `census_filter_invariant.rs`
/// counterpart) because what a term names is a **substring of a test path**, not
/// necessarily a module: `merge::tests::direction_symmetry` is three segments of
/// one, and `direction_symmetry` on its own would match the same tests. Treating
/// these as opaque strings and comparing them both ways round is the whole point
/// of the checks below.
fn terms_in(filter: &str) -> Vec<String> {
    filter
        .match_indices("test(")
        .filter_map(|(at, _)| {
            let rest = &filter[at + "test(".len()..];
            rest.find(')').map(|end| rest[..end].trim().to_string())
        })
        .filter(|term| !term.is_empty())
        .collect()
}

/// The terms of one `ci.yml` filter, asserted non-empty.
///
/// An empty list would make every loop below vacuous — the flattering direction,
/// and the one this whole file exists to keep loud.
fn required_terms_in(key: &str) -> Vec<String> {
    let terms = terms_in(&ci_filter(key));
    assert!(
        !terms.is_empty(),
        "{key} yields no test(...) terms; ci.yml's formatting changed and the \
         checks reading it would be vacuous"
    );
    terms
}

/// Every `.rs` file under a directory, recursively.
fn rust_sources_under(dir: &std::path::Path) -> Vec<PathBuf> {
    let mut found = Vec::new();
    let Ok(entries) = std::fs::read_dir(dir) else {
        return found;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            found.extend(rust_sources_under(&path));
        } else if path.extension().and_then(|e| e.to_str()) == Some("rs") {
            found.push(path);
        }
    }
    found
}

/// Assertion 1: the filter is non-empty.
///
/// Split out from the assertions that use it so a formatting change in `ci.yml`
/// reports as "the filter stopped parsing" rather than as "the shards stopped
/// negating it", which would send the reader to the wrong file.
#[test]
fn the_oracle_only_filter_parses_to_something() {
    let filter = ci_filter("ORACLE_ONLY_FILTER");
    assert!(
        filter.contains("test("),
        "ORACLE_ONLY_FILTER names no test: {filter:?}"
    );
}

/// Assertion 2: `test` negates it, or the change saved nothing.
#[test]
fn the_test_job_negates_the_oracle_only_filter() {
    let selections = selections_in("test");
    assert!(
        !selections.is_empty(),
        "the `test` job hands nextest no -E selection; this guard is blind"
    );
    assert!(
        selections
            .iter()
            .all(|s| s.contains("not ($ORACLE_ONLY_FILTER)")),
        "the `test` job does not negate ORACLE_ONLY_FILTER, so nothing was \
         moved off that shard: {selections:?}"
    );
}

/// Assertion 3, the one that matters: `test-oracle` still runs them.
///
/// The saving is only legitimate because an armed run subsumes an un-armed one.
/// If `test-oracle` ever negates this filter too — or excludes it through
/// `ORACLE_EXCLUDE`, which is the same thing by another spelling — the census
/// runs nowhere, both shards go green, and both get *faster*.
#[test]
fn the_test_oracle_job_still_runs_what_the_test_job_dropped() {
    let selections = selections_in("test-oracle");
    assert!(
        !selections.is_empty(),
        "the `test-oracle` job hands nextest no -E selection; this guard is blind"
    );
    assert!(
        selections.iter().all(|s| !s.contains("ORACLE_ONLY_FILTER")),
        "`test-oracle` mentions ORACLE_ONLY_FILTER. It must not: the `test` job \
         drops these tests only because this job still runs them armed. \
         Negating here deletes the census while turning both shards green and \
         faster: {selections:?}"
    );

    // The same deletion, spelled through the other exclusion mechanism.
    //
    // Compared TERM against TERM and in both directions, because `test(P)` is a
    // substring predicate on the test's name: an `ORACLE_EXCLUDE` entry spelled
    // `test(direction_symmetry)` — the spelling every entry there already uses —
    // excludes this module without ever containing its full path. See the module
    // doc comment.
    assert_terms_do_not_overlap("ORACLE_EXCLUDE", |ours, theirs| {
        format!(
            "`{ours}` (ORACLE_ONLY_FILTER) and `{theirs}` (ORACLE_EXCLUDE) select \
             overlapping tests: `test(P)` matches any test whose name CONTAINS P, \
             so one of these terms subsumes the other. Those tests are negated \
             from `test` by ORACLE_ONLY_FILTER and excluded from `test-oracle` by \
             ORACLE_EXCLUDE — they run NOWHERE, both shards go green, and both get \
             faster. That is the failure this file exists to catch."
        )
    });
}

/// Assertions 3b and 4, which differ only in the filter they read.
///
/// For every term of `ORACLE_ONLY_FILTER` and every term of the foreign filter,
/// neither may contain the other. Containment either way means nextest's
/// substring predicate selects overlapping tests, which is what the two callers'
/// messages then explain the consequence of.
fn assert_terms_do_not_overlap(foreign_key: &str, message: impl Fn(&str, &str) -> String) {
    let ours = required_terms_in("ORACLE_ONLY_FILTER");
    let theirs = required_terms_in(foreign_key);
    for our_term in &ours {
        for their_term in &theirs {
            assert!(
                !our_term.contains(their_term.as_str()) && !their_term.contains(our_term.as_str()),
                "{}",
                message(our_term, their_term)
            );
        }
    }
}

/// Assertion 4: not also in `CENSUS_FILTER`, which would route it to an archive
/// that cannot contain it.
///
/// `CENSUS_FILTER` is negated from both shards and selected by `censuses` off
/// the soak archive, built `-p ferro-hgvs-soak-tests`. That archive holds the
/// soak driver's test target and nothing else — in particular no `#[cfg(test)]`
/// module from under `src/`. A unit test named in both filters is therefore
/// negated everywhere and selected nowhere, and `--no-tests=fail` does not
/// notice because the step's other terms still match.
#[test]
fn no_oracle_only_test_is_also_claimed_by_the_censuses_job() {
    assert_terms_do_not_overlap("CENSUS_FILTER", |ours, theirs| {
        format!(
            "`{ours}` (ORACLE_ONLY_FILTER) and `{theirs}` (CENSUS_FILTER) select \
             overlapping tests — `test(P)` matches any test whose name CONTAINS P. \
             CENSUS_FILTER is negated from both shards and selected off the soak \
             archive, which contains no `src/` unit tests, so those tests are \
             negated everywhere and selected nowhere. `--no-tests=fail` does not \
             notice, because the step's other terms still match."
        )
    });
}

/// Assertion 5: the filter still names something that exists.
///
/// `the_oracle_only_filter_parses_to_something` asserts only that the string
/// contains `test(`. Rename the module and the filter matches nothing: `test`
/// negates nothing, `test-oracle` runs it as before, and the duplicate quietly
/// returns to the debug shard. That direction is *safe* — no coverage is lost —
/// which is why it is a separate, weaker assertion rather than folded into the
/// three above. It is still silent config rot, and it silently un-does the
/// change this filter was added to make.
///
/// The check is deliberately cheap and structural: every `::`-separated segment
/// of each term must be declared as a module somewhere under `src/`. It does not
/// shell out to `cargo nextest list` — nothing else in `tests/it/` does, and a
/// guard that needs a built archive to answer is a guard that gets disabled.
#[test]
fn every_oracle_only_term_still_names_a_module_that_exists() {
    let sources: Vec<String> = rust_sources_under(&repo_root().join("src"))
        .iter()
        .filter_map(|path| std::fs::read_to_string(path).ok())
        .collect();
    assert!(
        !sources.is_empty(),
        "no sources read from src/; this guard cannot see what exists"
    );

    for term in required_terms_in("ORACLE_ONLY_FILTER") {
        for segment in term.split("::") {
            let declaration = format!("mod {segment}");
            assert!(
                sources.iter().any(|text| text.contains(&declaration)),
                "ORACLE_ONLY_FILTER names `test({term})`, but no `{declaration}` \
                 is declared anywhere under src/. A rename fails SAFE — the filter \
                 matches nothing, so `test` negates nothing and the test runs in \
                 both shards again — but the carve-out this variable exists for is \
                 silently gone. Update the filter, or drop it."
            );
        }
    }
}
