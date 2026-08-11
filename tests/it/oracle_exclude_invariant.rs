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

/// The brace-grouped form of the same import, `conformance::{…spec_corpus…}`.
const CORPUS_IMPORT_GROUP: &str = "conformance::{";

/// Whether `text` imports the spec corpus, in **either** spelling.
///
/// A single `contains(CORPUS_IMPORT)` misses `use
/// ferro_hgvs::conformance::{spec_corpus, summary};`, and it misses it in the
/// flattering direction: a module written that way is never demanded in
/// `ORACLE_EXCLUDE`, so its census is taken with the oracle armed and reads
/// *better* than the truth. That is the exact failure this file exists to close,
/// so the matcher may not be blind to a spelling rustfmt will happily produce.
///
/// Shared by both scans deliberately. Widening one and not the other is how the
/// same rule kept in two copies drifts — the reason [`CORPUS_IMPORT`] is a
/// constant rather than an inline literal in the first place.
fn imports_corpus(text: &str) -> bool {
    if text.contains(CORPUS_IMPORT) {
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
    // `CORPUS_IMPORT` path literal, so the check above already catches them.
    text.match_indices(CORPUS_IMPORT_GROUP).any(|(at, _)| {
        let rest = &text[at + CORPUS_IMPORT_GROUP.len()..];
        rest.find('}').is_some_and(|end| {
            rest[..end]
                .split(',')
                .any(|item| leading_identifier(item) == "spec_corpus")
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
            imports_corpus(&std::fs::read_to_string(&path).ok()?).then_some(name)
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
}
