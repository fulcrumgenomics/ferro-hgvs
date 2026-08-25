//! Liveness check for the pairing between `CENSUS_FILTER` and the job that owns it.
//!
//! `ci.yml`'s `CENSUS_FILTER` names the slow census modules that moved off the
//! `test` / `test-oracle` partitions and onto the optimized archive. `test` and
//! `test-oracle` **negate** it; the `censuses` and `censuses-plain` jobs
//! **select** on it — one with the seam oracles armed, one without, because
//! `spec_conformance_axis` may not be measured with `FERRO_ASSERT_IDEMPOTENT`
//! set (it counts the very rows that oracle panics on). The two are separate
//! jobs (not two steps of one) so their ~43s and ~50s halves run in parallel.
//!
//! **The failure mode is the silent, flattering kind, and there are two of it.**
//!
//! A module named in `CENSUS_FILTER` but selected by *neither* step is negated
//! everywhere and run **nowhere**. Nothing goes red. CI gets *faster*, which
//! reads as the improvement this filter exists to deliver rather than as the
//! coverage deletion it is. That is the same shape as `sweep_filter_invariant`'s
//! unnamed-consumer case, and the same shape as #1460 / #1478, where a corpus
//! that could not build a thing reported its absence as a zero.
//!
//! A module selected by a step but *not* named in `CENSUS_FILTER` fails the
//! other way: it is not negated, so it runs in `censuses` **and** in the shard
//! that hashes it, paying for it twice while looking like a carve-out.
//!
//! Neither is caught by the guards already in place. `--no-tests=fail` on each
//! step catches a selection matching nothing, i.e. a rename that empties a whole
//! step; it cannot see one module of three going missing. `CENSUS_FILTER`'s own
//! comment records the partition against `test` / `test-oracle`, which is a
//! claim about the negated side.
//!
//! Modelled on `sweep_filter_invariant.rs` and `oracle_exclude_invariant.rs`,
//! which make the same class of silent config rot loud for the sweep-seed knob
//! and the oracle exclusion.
//!
//! # Why this reads `ci.yml` itself rather than sharing a helper
//!
//! This is the third file in `tests/it/` to parse a folded scalar out of
//! `ci.yml`, and the duplication is deliberate. These three are *liveness*
//! guards: a bug in one shared extractor would neuter all three at once, and it
//! would do so in the flattering direction — an extractor that silently returned
//! an empty filter makes every one of them vacuous rather than red. Each guard
//! asserts its own read is non-empty, so three independent readers are three
//! independent chances to notice.

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
///
/// Reading only the first line would silently exempt whichever module happens to
/// sit on a continuation line — and `CENSUS_FILTER` spans two, so a first-line
/// read would exempt `spec_conformance_axis` specifically.
fn ci_filter(key: &str) -> String {
    let ci = ci_text();
    let mut lines = ci.lines();
    let header = lines
        .by_ref()
        .find(|line| line.trim_start().starts_with(&format!("{key}:")))
        .unwrap_or_else(|| panic!("ci.yml defines {key}"));

    // An inline value is returned rather than merely tolerated: a perfectly
    // valid `CENSUS_FILTER: test(a) + test(b)` must not parse as empty and make
    // every assertion below vacuous.
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
///
/// Blank and comment lines are kept, because they do not end a block — the
/// `censuses` job's steps are separated by both, and a `take_while` that stopped
/// at the first one would read only its first step and never see the second
/// selection.
fn job_lines(job: &str) -> Vec<String> {
    let ci = ci_text();
    let header = format!("  {job}:");
    let mut lines = ci.lines().skip_while(|line| *line != header);
    assert!(
        lines.next().is_some(),
        "ci.yml defines no `{job}` job; this guard cannot see what it runs"
    );
    lines
        // A non-blank, non-comment line back at job-key indent ends the block.
        .take_while(|line| {
            line.starts_with("    ") || line.trim().is_empty() || line.trim_start().starts_with('#')
        })
        .map(str::to_string)
        .collect()
}

/// Every `-E "…"` selection a job hands to nextest, in file order.
///
/// Comment lines are dropped here rather than in [`job_lines`], so a `-E` that
/// has been commented out is not counted as a live selection while still
/// bounding the block correctly.
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

/// The jobs that between them SELECT `CENSUS_FILTER`, armed and un-armed. The
/// armed half (`normalize_axis_preserving` &c.) and the un-armed half
/// (`spec_conformance_axis`, which may not be measured with
/// `FERRO_ASSERT_IDEMPOTENT` set) run as two SEPARATE jobs so they parallelize;
/// together their steps must cover the whole filter and nothing more.
const CENSUS_JOBS: [&str; 2] = ["censuses", "censuses-plain"];

/// Every `-E` selection across both census jobs, in job order.
fn census_selections() -> Vec<String> {
    CENSUS_JOBS.iter().flat_map(|j| selections_in(j)).collect()
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

/// The `censuses` job must actually run every module the filter reserves for it.
///
/// Without this, a module dropped from the job's two steps — or renamed in one
/// place and not the other — is negated by `test` and `test-oracle` and selected
/// by nothing, so it runs in no job at all and CI reports a faster green run.
#[test]
fn every_module_named_in_the_census_filter_is_run_by_the_censuses_job() {
    let named = modules_named_in(&ci_filter("CENSUS_FILTER"));
    assert!(
        !named.is_empty(),
        "CENSUS_FILTER names no modules; its formatting changed and every \
         assertion in this file would be vacuous"
    );

    let selections = census_selections();
    assert!(
        !selections.is_empty(),
        "the census jobs pass no -E selection to nextest; their shape changed \
         and this guard cannot see what they run"
    );
    let selected: Vec<String> = selections
        .iter()
        .flat_map(|s| modules_named_in(s))
        .collect();

    let unrun: Vec<&String> = named.iter().filter(|m| !selected.contains(m)).collect();
    assert!(
        unrun.is_empty(),
        "ci.yml's CENSUS_FILTER reserves these modules for the census jobs, but \
         none of their steps selects them: {unrun:?}\n\
         `test` and `test-oracle` negate CENSUS_FILTER, so a module here runs in NO \
         job — CI stays green and gets faster, which is the coverage loss this \
         guard exists to make loud.\n\
         Either add them to a step's -E in the `censuses` job, or remove them from \
         CENSUS_FILTER so the shards run them again."
    );
}

/// …and it must run nothing else, or that module is paid for twice.
///
/// A step selecting a module absent from `CENSUS_FILTER` does not lose coverage
/// — it duplicates it. The module still hashes into a `test` shard, so the
/// carve-out adds a runner without removing any work, which is the opposite of
/// what the job is for.
#[test]
fn the_censuses_job_runs_nothing_the_census_filter_does_not_reserve() {
    let named = modules_named_in(&ci_filter("CENSUS_FILTER"));
    let selected: Vec<String> = census_selections()
        .iter()
        .flat_map(|s| modules_named_in(s))
        .collect();
    assert!(
        !selected.is_empty(),
        "the census jobs select no modules; their shape changed"
    );

    let unreserved: Vec<&String> = selected.iter().filter(|m| !named.contains(m)).collect();
    assert!(
        unreserved.is_empty(),
        "the census jobs select these modules, but CENSUS_FILTER does not \
         name them: {unreserved:?}\n\
         `test` and `test-oracle` negate CENSUS_FILTER, so anything missing from it \
         still runs in whichever shard hashes it — these modules are being run \
         twice, and the carve-out is costing a runner without saving one."
    );
}

/// No module may be selected by both census jobs.
///
/// The two jobs differ only in whether the seam oracles are armed, so a module in
/// both runs its whole corpus twice per CI run — the duplication the split was
/// made to avoid, relocated rather than removed. Each job runs a single armed or
/// un-armed step, so each carries exactly one `-E` selection.
#[test]
fn no_module_is_selected_by_both_census_jobs() {
    let armed_sel = selections_in("censuses");
    let plain_sel = selections_in("censuses-plain");
    assert_eq!(
        (armed_sel.len(), plain_sel.len()),
        (1, 1),
        "each census job is documented as one step (armed / un-armed); found \
         {} armed and {} plain selection(s). If that is deliberate, this guard \
         needs updating along with the jobs' comments.",
        armed_sel.len(),
        plain_sel.len()
    );

    let armed = modules_named_in(&armed_sel[0]);
    let plain = modules_named_in(&plain_sel[0]);
    // `test(P)` selects any test whose name CONTAINS `P`, so two terms select
    // overlapping tests when either one contains the other — exact equality
    // misses a subsuming term (armed `test(spec_conformance)` against plain
    // `test(spec_conformance_axis)` picks the same module in both jobs). Compare
    // both directions, as `oracle_only_filter_invariant.rs::
    // assert_terms_do_not_overlap` already does.
    let both: Vec<(&String, &String)> = armed
        .iter()
        .flat_map(|a| {
            plain
                .iter()
                .filter(move |p| a.contains(p.as_str()) || p.contains(a.as_str()))
                .map(move |p| (a, p))
        })
        .collect();
    assert!(
        both.is_empty(),
        "these modules are selected by BOTH census jobs, so their corpora run \
         twice per CI run: {both:?}"
    );
}
