//! Liveness check for the pairing between the sweep seed knob and CI.
//!
//! `sweep_seeds` (#1295) draws a small prefix of the corpus by default, and
//! `ci.yml`'s `sweeps` job restores the exhaustive run by setting
//! `FERRO_SWEEP_SEEDS: full`. That job selects on `SWEEP_FILTER`, so the pairing
//! only holds while every module consulting the knob is named in that filter.
//!
//! Nothing enforced it. #1295's own comment said "if a sweep is ever moved out
//! of `SWEEP_FILTER`, it silently drops to the prefix here — check both when
//! editing that filter", which is a human instruction, and #1397 asks for the
//! invariant instead.
//!
//! **The failure mode is the expensive kind.** A new sweep module that consults
//! the knob but is not named in `SWEEP_FILTER` runs at the default prefix in
//! `test` / `test-oracle` and never receives `full` anywhere. CI stays green and
//! gets *faster*, which reads as success rather than as the coverage loss it is.
//!
//! The existing guards do not cover it. `SWEEP_FILTER`'s own comment records a
//! partition — the negated selection plus the filter's own selection account for
//! the whole suite — and `sweeps` passes `--no-tests=fail`. Both catch a filter
//! that matches *nothing*, i.e. a rename. Neither catches a module missing *from*
//! it: an unnamed consumer simply lands on the negated side, where it still runs
//! and still counts, just never at the full corpus.
//!
//! Modelled on `coderabbit_config_paths.rs`, which makes the same class of
//! silent config rot loud.

use std::path::PathBuf;

/// This file. Excluded from the consumer scan below, which would otherwise
/// match it — not because of the prose (that spells the knob in backticks, and
/// the scan looks for the call form `sweep_seeds(`), but because the scan's own
/// **matcher literal** is that call form, and it lives here. So this file
/// contains the string it searches for, and without this exclusion the guard
/// would demand it be named in a CI filter that has no reason to run it.
///
/// Worth stating precisely: the self-match is a property of how the scan is
/// written, so it would survive a rewrite of every doc comment above and
/// disappear the moment the matcher changed shape.
const SELF: &str = "sweep_filter_invariant.rs";

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The `SWEEP_FILTER` value from `ci.yml`, folded to one line.
///
/// It is a `>-` folded scalar spanning two lines, so the continuation is read
/// and joined rather than the first line being taken alone — reading only the
/// first line would silently exempt whichever module happens to sit on the
/// second, which is exactly the blind spot this file exists to close.
fn sweep_filter() -> String {
    let ci = std::fs::read_to_string(repo_root().join(".github/workflows/ci.yml"))
        .expect("ci.yml is readable");
    let mut lines = ci.lines();
    let header = lines
        .by_ref()
        .find(|l| l.trim_start().starts_with("SWEEP_FILTER:"))
        .expect("ci.yml defines SWEEP_FILTER");
    // An inline value is returned here, not merely permitted. The header check
    // used to accept inline syntax and then fall through to the continuation
    // scan, which reads nothing for a one-line value — so a perfectly valid
    // `SWEEP_FILTER: test(a) + test(b)` parsed as empty and tripped the "changed
    // formatting" assertion below. Accepting a shape the parser cannot read is
    // its own small version of the rot this file exists to catch.
    let inline = header
        .split_once(':')
        .map(|(_, value)| value.trim())
        .unwrap_or_default();
    if !inline.is_empty() && !inline.starts_with(['>', '|']) {
        return inline.to_string();
    }
    assert!(
        inline.starts_with(['>', '|']),
        "SWEEP_FILTER is neither a block scalar nor an inline value: {header:?}"
    );

    // A block scalar's body is every following line indented deeper than the
    // key, up to the first line that is not.
    let key_indent = header.len() - header.trim_start().len();
    let mut value = String::new();
    for line in lines {
        if line.trim().is_empty() {
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
        "SWEEP_FILTER parsed as empty; ci.yml's formatting changed"
    );
    value
}

/// Integration-test modules that call the seed knob, by module name.
///
/// Only the top level of `tests/it/` is scanned: `common/cis_apply_oracle.rs`
/// *defines* the knob rather than consulting it, and is not a test target CI
/// could select anyway.
fn seed_knob_consumers() -> Vec<String> {
    let dir = repo_root().join("tests/it");
    let mut consumers: Vec<String> = std::fs::read_dir(&dir)
        .expect("tests/it is readable")
        .filter_map(|entry| {
            let path = entry.expect("readable dir entry").path();
            let name = path.file_name()?.to_str()?.to_string();
            if !name.ends_with(".rs") || name == SELF {
                return None;
            }
            let text = std::fs::read_to_string(&path).ok()?;
            // The call form, not the bare name, so a mention in prose does not
            // count as consumption.
            text.contains("sweep_seeds(")
                .then(|| name.trim_end_matches(".rs").to_string())
        })
        .collect();
    consumers.sort();
    consumers
}

/// Every module that draws its corpus through the seed knob must be named in
/// `SWEEP_FILTER`, or it never runs at the full corpus anywhere.
#[test]
fn every_seed_knob_consumer_is_named_in_the_sweep_filter() {
    let filter = sweep_filter();
    let consumers = seed_knob_consumers();
    assert!(
        !consumers.is_empty(),
        "no module calls the seed knob; either the scan broke or the knob was \
         removed — both should fail loudly rather than pass vacuously"
    );

    let missing: Vec<&String> = consumers
        .iter()
        .filter(|module| !filter.contains(&format!("test({module})")))
        .collect();
    assert!(
        missing.is_empty(),
        "these modules draw their corpus through the seed knob but are not named \
         in ci.yml's SWEEP_FILTER, so they run at the default prefix in every job \
         and never receive FERRO_SWEEP_SEEDS=full: {missing:#?}\n\
         SWEEP_FILTER is:{filter}\n\
         Add `+ test(<module>)` to it, or stop drawing through the knob."
    );
}

/// The converse: every module `SWEEP_FILTER` names must exist and must actually
/// consult the knob.
///
/// A filter naming a module that no longer draws through the knob is not a
/// coverage hole, but it is a claim that stopped being true — and the next person
/// editing the filter would reasonably trust it. `--no-tests=fail` on the
/// `sweeps` job catches only a filter matching nothing at all, not one stale
/// entry among several.
#[test]
fn every_module_named_in_the_sweep_filter_consumes_the_seed_knob() {
    let filter = sweep_filter();
    let consumers = seed_knob_consumers();

    let named: Vec<String> = filter
        .match_indices("test(")
        .filter_map(|(at, _)| {
            let rest = &filter[at + "test(".len()..];
            rest.find(')').map(|end| rest[..end].trim().to_string())
        })
        .collect();
    assert!(
        !named.is_empty(),
        "SWEEP_FILTER names no modules; its formatting changed: {filter}"
    );

    let stale: Vec<&String> = named
        .iter()
        .filter(|module| !consumers.contains(module))
        .collect();
    assert!(
        stale.is_empty(),
        "ci.yml's SWEEP_FILTER names these modules, but they do not draw their \
         corpus through the seed knob — the filter is claiming coverage it does \
         not provide: {stale:#?}\n\
         Remove them from SWEEP_FILTER, or have them use the knob."
    );
}
