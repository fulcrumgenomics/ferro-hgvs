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

/// The import that marks a module as built on the spec corpus.
///
/// The path form, not the bare word `spec_corpus`, so a mention in prose — of
/// which the two modules have many — does not count as consumption.
const CORPUS_IMPORT: &str = "conformance::spec_corpus";

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
            text.contains(CORPUS_IMPORT)
                .then(|| name.trim_end_matches(".rs").to_string())
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
/// [`corpus_modules`] recognises a consumer by the literal import path
/// [`CORPUS_IMPORT`], which is right for avoiding prose false positives but
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
            std::fs::read_to_string(&path)
                .ok()?
                .contains(CORPUS_IMPORT)
                .then_some(name)
        })
        .collect();
    importers.sort();

    assert!(
        importers.is_empty(),
        "these shared helpers import the spec corpus: {importers:#?}\n\
         A module consuming the corpus THROUGH one of them carries no \
         `{CORPUS_IMPORT}` literal of its own, so \
         `every_spec_corpus_module_is_named_in_the_oracle_exclude` would not \
         demand it be named in ORACLE_EXCLUDE — and it would then measure with \
         the seam oracles armed, which reports a better census rather than \
         going red.\n\
         Either import the corpus directly in each consuming module, or teach \
         `corpus_modules` to follow this helper."
    );
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
