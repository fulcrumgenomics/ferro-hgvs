//! Generator completeness invariants (#1550).
//!
//! Generators under `examples/` build the corpora, fixtures and tables the rest
//! of the suite adjudicates against. They share one failure mode: a fallible
//! step whose failure is representable as a legitimate value — `unwrap_or_default()`,
//! `else { continue }`, a discarded `Result` — drops part of the population and
//! the generator writes its artifact anyway. A partial run and a clean run then
//! produce indistinguishable output, and the only thing separating them is
//! whether a reviewer read the right five lines.
//!
//! Two mechanical guards live here. Neither can see a *semantic* drop; that is
//! what [`CaptureLedger`](ferro_hgvs::conformance::completeness::CaptureLedger)
//! is for, and the second guard is the prompt that points generators at it.
//!
//! 1. [`examples_with_unit_tests_opt_into_running_them`] — a generator carrying
//!    `#[cfg(test)]` must set `test = true` on its cargo target, because cargo
//!    does not build a target's tests unless it opts in. `Cargo.toml` has stated
//!    this rule in prose since the two spec generators became `[[bin]]`s;
//!    nothing checked it, and `report_conformance_reference_gaps` shipped a
//!    `#[test]` that had never once executed. **Only half of this one is
//!    exact.** The manifest side is parsed as TOML, so `test = true` is read
//!    rather than guessed; the source side is a scan for the token
//!    `#[cfg(test)]`, which matches it in a doc comment or a string literal and
//!    — the direction that costs something — misses the equivalent
//!    `#[cfg(all(test, feature = "dev"))]`. A generator gating its tests that
//!    way is invisible to this guard.
//! 2. [`artifact_writers_use_the_capture_ledger_or_are_allowlisted`] — a
//!    generator whose source *looks like* it writes a file must either look like
//!    it routes its population through the ledger, or name itself in
//!    [`LEDGER_EXEMPT`] against #1550.
//!
//! **Read the second guard's reach honestly, because this file's whole subject
//! is a check that reads as stronger than it is.** It is a pair of substring
//! heuristics over source text — see [`appears_to_write_an_artifact`] and
//! [`appears_to_route_through_the_ledger`] for exactly which strings, and for
//! the write idioms it is known to miss. It is a **floor that catches the common
//! case**: a generator that writes with `std::fs::write` and never mentions the
//! ledger cannot be added in silence. It is *not* a proof that every writer is
//! accounted for, and it cannot be — deciding whether a population was routed
//! through the ledger is a semantic question, and answering it mechanically
//! would mean parsing and type-resolving Rust for a three-line rule. The
//! semantic guarantee is carried by the type and by review; this guard's job is
//! to make sure the question gets asked.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

// ---------------------------------------------------------------------------
// Reading the manifest and the generator sources
// ---------------------------------------------------------------------------

/// One cargo target whose source lives under `examples/`.
#[derive(Debug, Clone)]
struct GeneratorTarget {
    /// Cargo target name (`report_conformance_reference_gaps`).
    name: String,
    /// Which table declared it — `[[example]]` or `[[bin]]` — quoted verbatim in
    /// failures so the fix names the block to edit.
    table: &'static str,
    /// Repo-relative source path (`examples/report_conformance_reference_gaps.rs`).
    source: PathBuf,
    /// Whether the target sets `test = true`.
    builds_tests: bool,
}

/// Every `[[example]]` and `[[bin]]` target whose source lives under `examples/`.
///
/// Parsed as TOML rather than substring-matched: `Cargo.toml` discusses
/// `test = true` at length in comments, and a `[[bench]]` setting it would be a
/// false positive for the rule this file enforces.
fn generator_targets() -> Vec<GeneratorTarget> {
    let manifest_text =
        std::fs::read_to_string(repo_root().join("Cargo.toml")).expect("read Cargo.toml");
    let manifest: toml::Value = manifest_text.parse().expect("parse Cargo.toml as TOML");

    let mut targets = Vec::new();
    for (key, table) in [("example", "[[example]]"), ("bin", "[[bin]]")] {
        let Some(entries) = manifest.get(key).and_then(toml::Value::as_array) else {
            continue;
        };
        for entry in entries {
            let name = entry
                .get("name")
                .and_then(toml::Value::as_str)
                .unwrap_or_else(|| panic!("a {table} target with no `name`"))
                .to_string();
            // An `[[example]]` without an explicit `path` is `examples/<name>.rs`
            // by cargo's own convention; the two spec generators are `[[bin]]`s
            // that name their `examples/` source explicitly.
            let source = match entry.get("path").and_then(toml::Value::as_str) {
                Some(path) => PathBuf::from(path),
                None if key == "example" => PathBuf::from("examples").join(format!("{name}.rs")),
                None => continue,
            };
            if !source.starts_with("examples") {
                continue;
            }
            targets.push(GeneratorTarget {
                name,
                table,
                builds_tests: entry
                    .get("test")
                    .and_then(toml::Value::as_bool)
                    .unwrap_or(false),
                source,
            });
        }
    }
    targets
}

/// Every `examples/*.rs` that is a generator entry point (top level, so the
/// `common/` modules that are `#[path]`-included into them are excluded).
///
/// **Scope difference worth knowing:** this reads the *directory*, while
/// [`generator_targets`] reads the *manifest*. A declared target whose `path`
/// points into a subdirectory (`examples/sub/gen.rs`) would therefore be seen by
/// the `test = true` invariant and not by the artifact-writer ratchet. No such
/// target exists today, and
/// [`every_generator_source_is_declared_as_a_cargo_target`] keeps the two views
/// in agreement for everything at the top level.
fn generator_sources_on_disk() -> BTreeSet<PathBuf> {
    std::fs::read_dir(repo_root().join("examples"))
        .expect("read examples/")
        .map(|entry| entry.expect("examples/ dir entry").path())
        .filter(|path| path.extension().and_then(|e| e.to_str()) == Some("rs"))
        .map(|path| PathBuf::from("examples").join(path.file_name().expect("example file name")))
        .collect()
}

/// The `#[path = "…"]` module includes declared in `text`, resolved relative to
/// `containing_dir`.
///
/// A one-line scan rather than a parse: the includes are all of the literal form
/// `#[path = "common/x.rs"]`, and a test that needed `syn` to enforce a
/// three-line rule would not be worth its own maintenance.
///
/// **Boundary, stated rather than implied:** only `#[path]` is followed. A plain
/// `mod helper;` — cargo's own `examples/<target>/helper.rs` convention — is not
/// resolved, so `#[cfg(test)]` reached that way is invisible to the guards
/// below. There is no such target today; the repo uses `#[path]` throughout, and
/// the day one appears this comment is the thing to read.
fn path_includes(text: &str, containing_dir: &Path) -> Vec<PathBuf> {
    let mut found = Vec::new();
    for line in text.lines() {
        let Some(rest) = line.trim().strip_prefix("#[path") else {
            continue;
        };
        let Some(open) = rest.find('"') else { continue };
        let Some(close) = rest[open + 1..].find('"') else {
            continue;
        };
        found.push(containing_dir.join(&rest[open + 1..open + 1 + close]));
    }
    found
}

/// Every source file compiled into a target: its own file plus, transitively,
/// the `common/` modules it `#[path]`-includes.
///
/// The transitive walk is the point. `test = true` is a property of the *target*,
/// so a `#[cfg(test)]` in a shared `common/` module is equally dead in every
/// target that includes it without opting in — and reading only the entry point
/// would miss exactly that.
///
/// An unreadable source **panics** rather than being skipped. Skipping it would
/// be `else { continue }` on a fallible step, in the file whose module doc names
/// that idiom as the failure mode: an entry point that failed to read would
/// yield no sources, and a target with no sources reads as carrying no
/// `#[cfg(test)]` and writing nothing — so *both* guards would pass silently for
/// exactly the file nobody could read (#1550 review).
fn compiled_sources(entry_point: &Path) -> Vec<(PathBuf, String)> {
    let root = repo_root();
    let mut queue = vec![entry_point.to_path_buf()];
    let mut seen = BTreeSet::new();
    let mut sources = Vec::new();
    while let Some(relative) = queue.pop() {
        if !seen.insert(relative.clone()) {
            continue;
        }
        let absolute = root.join(&relative);
        let text = std::fs::read_to_string(&absolute).unwrap_or_else(|error| {
            panic!(
                "cannot read a generator source reached from {}: {} ({error}). \
                 Skipping it would make both guards in this file pass vacuously \
                 for that target.",
                entry_point.display(),
                absolute.display(),
            )
        });
        let containing_dir = relative.parent().unwrap_or(Path::new("")).to_path_buf();
        queue.extend(path_includes(&text, &containing_dir));
        sources.push((relative, text));
    }
    sources
}

// ---------------------------------------------------------------------------
// 1. A generator with unit tests must opt into running them
// ---------------------------------------------------------------------------

/// Cargo does not build an example's or a binary's `#[cfg(test)]` code unless
/// the target sets `test = true`, so a generator that carries unit tests without
/// it ships a test that has never executed. That is worse than no test, because
/// it reads as coverage: `Cargo.toml`'s own comment on the two spec generators
/// says so, most generator targets already set the flag, and
/// `report_conformance_reference_gaps` still shipped a `#[test]` that
/// `cargo nextest run -E 'binary(report_conformance_reference_gaps)'` reported
/// as `0 tests run` (#1550).
///
/// No count is quoted here on purpose. An earlier draft said "six targets", which
/// was true when #1550 was filed and stale by the time this branch rebased twice —
/// a comment that rots with every new generator, in a PR about claims that stop
/// being true. The test itself is the count.
///
/// The manifest side of this is parsed; the source side is a token scan for
/// `#[cfg(test)]` and misses `#[cfg(all(test, ..))]`. See the module docs — the
/// claim that this guard is exact end-to-end was itself an overclaim (#1551
/// review), which is the defect this whole file is about.
#[test]
fn examples_with_unit_tests_opt_into_running_them() {
    let mut violations = Vec::new();
    for target in generator_targets() {
        let carriers: Vec<String> = compiled_sources(&target.source)
            .into_iter()
            .filter(|(_, text)| text.contains("#[cfg(test)]"))
            .map(|(path, _)| path.display().to_string())
            .collect();
        if carriers.is_empty() || target.builds_tests {
            continue;
        }
        violations.push(format!(
            "{} `{}` compiles #[cfg(test)] code in {} but does not set `test = true`, \
             so those tests never run",
            target.table,
            target.name,
            carriers.join(", "),
        ));
    }
    assert!(
        violations.is_empty(),
        "every generator target carrying #[cfg(test)] must set `test = true` in \
         Cargo.toml — cargo does not build a target's tests unless it opts in:\n  {}",
        violations.join("\n  "),
    );
}

/// `autoexamples = false`, so a generator that is not declared is not built at
/// all — the whole-file version of the silent drop this file is about. It would
/// also make the guard above vacuous for that file.
#[test]
fn every_generator_source_is_declared_as_a_cargo_target() {
    let declared: BTreeSet<PathBuf> = generator_targets()
        .into_iter()
        .map(|target| target.source)
        .collect();
    let undeclared: Vec<String> = generator_sources_on_disk()
        .difference(&declared)
        .map(|path| path.display().to_string())
        .collect();
    assert!(
        undeclared.is_empty(),
        "`autoexamples = false`, so an undeclared generator is never built and \
         its tests never run; declare each in Cargo.toml: {undeclared:?}",
    );
}

// ---------------------------------------------------------------------------
// 2. The artifact-writer ratchet
// ---------------------------------------------------------------------------

/// Generators that write an artifact without routing their population through
/// [`CaptureLedger`], each recorded against #1550.
///
/// **This list may only shrink**, and the two shrink directions are not equally
/// trustworthy — say which you are relying on:
///
/// - A row whose generator has been **deleted**, or which **no longer matches
///   [`appears_to_write_an_artifact`]**, fails the test below. The first is a
///   fact; the second is a heuristic, so it can also fire because the generator
///   switched to a write idiom the predicate does not know.
/// - A row whose generator now **matches [`appears_to_route_through_the_ledger`]**
///   also fails, and that direction is the one that can *lose* a row. It is
///   deliberately hard to trip by accident (entry point only, import plus a
///   `finish` call — see that function), because deleting a row on the strength
///   of a substring would put a generator permanently outside this list on
///   evidence nobody checked. If it fires, confirm the generator really does
///   account for its population before deleting the row; if it does not, the
///   right fix is the generator, not the list.
///
/// Migrating these is deliberately *not* #1550's job: #1544, #1545 and #1549 own
/// several of these files, and rewriting them here would mix mechanical churn
/// into a new contract and conflict with three open PRs.
///
/// Adding a row is a legitimate act — it is the "I considered this and declined"
/// answer, exactly as `Representation-Change: none` is. What is not legitimate is
/// silence: a new generator that *the predicate can see* writing an artifact, and
/// that does neither, fails the build.
const LEDGER_EXEMPT: &[(&str, &str)] = &[
    (
        "examples/build_conformance_snapshot.rs",
        "#1550: pre-existing; harvests snapshot bases per corpus accession",
    ),
    (
        "examples/derive_ng_placements.rs",
        "#1550: pre-existing; derives NG_ placements from the prepared reference",
    ),
    (
        "examples/derive_tx_placements.rs",
        "#1550: pre-existing; derives transcript placements from cdot",
    ),
    (
        "examples/dump_confluence_divergences.rs",
        "#1550: landed with #1549 while this contract was in review; migrating it \
         belongs to that generator's own change",
    ),
    (
        "examples/dump_normalized_corpus.rs",
        "#1550: pre-existing; the representation-change measurement harness",
    ),
    (
        "examples/extract_biocommons_windows.rs",
        "#1550: pre-existing; records the reference windows one pass touched",
    ),
    (
        "examples/extract_case_harvest_windows.rs",
        "#1550: landed with #1545 while this contract was in review; migrating it \
         belongs to that generator's own change",
    ),
    (
        "examples/extract_mutalyzer_genomic_windows.rs",
        "#1550: pre-existing; records the reference windows one pass touched",
    ),
    (
        "examples/extract_spec_enumeration_windows.rs",
        "#1550: pre-existing; records the reference windows one pass touched",
    ),
    (
        "examples/extract_spec_worked_example_windows.rs",
        "#1550: landed with #1544 while this contract was in review; migrating it \
         belongs to that generator's own change",
    ),
    (
        "examples/extract_split_member_separation_windows.rs",
        "#1550: landed with #1549 while this contract was in review; it already \
         implements this refusal by hand — it is the file `CaptureLedger` was \
         promoted out of — so migrating it belongs to that generator's own change",
    ),
    (
        "examples/generate_cis_confluence_corpus.rs",
        "#1550: pre-existing; synthesises the cis-confluence corpus",
    ),
    (
        "examples/generate_conformance_summary.rs",
        "#1550: pre-existing; renders the failure-pattern views",
    ),
    (
        "examples/generate_perf_tables.rs",
        "#1550: pre-existing; renders the README performance tables",
    ),
    (
        "examples/generate_spec_enumeration.rs",
        "#1550: pre-existing; harvests the spec checkout into the enumeration",
    ),
    (
        "examples/generate_spec_fixture.rs",
        "#1550: pre-existing; harvests the spec checkout into the fixture",
    ),
    (
        "examples/generate_tool_support_tables.rs",
        "#1550: pre-existing; renders the committed tool-support tables",
    ),
    (
        "examples/harvest_multi_member_cis.rs",
        "#1550: pre-existing; harvests multi-member cis rows into a fixture",
    ),
    (
        "examples/projector-bench.rs",
        "#1550: pre-existing; writes benchmark output, not a committed corpus",
    ),
];

/// Write idioms this file can recognise. **Not exhaustive, and cannot be** —
/// there is no closed set of ways to write a file in Rust.
///
/// The rule for what belongs here is "the call **names a path**", which is what
/// separates opening an artifact from wrapping a sink somebody else opened.
/// `BufWriter::new(` and `serde_json::to_writer(` are therefore deliberately
/// *out*: they say nothing about whether the destination is a file (the only
/// `BufWriter` in `examples/` that is not already wrapping a `File::create` is
/// over `stdout`), and every file-backed use of them passes through one of the
/// four below anyway.
///
/// Measured misses, recorded so nobody rediscovers them the hard way (#1551
/// review): the first revision listed only `fs::write(` and `File::create(`, so
/// `OpenOptions::new()…open()` + `write_all` and `csv::Writer::from_path(..)`
/// were both **invisible** to the ratchet. Those two are now covered. What is
/// still not covered is anything that reaches a writer through a helper's return
/// value, a trait object, or a crate whose constructor is spelled some other way.
const WRITE_IDIOMS: &[&str] = &[
    "fs::write(",
    "File::create(",
    "OpenOptions::new(",
    "from_path(",
];

/// Whether a generator's compiled sources *look like* they write a file.
///
/// A substring scan over source text, deliberately: see [`WRITE_IDIOMS`] for the
/// list and for what it is known to miss. Two consequences worth being explicit
/// about, because the PR that introduced this check originally claimed more than
/// it delivers:
///
/// - A **false negative** (a writer using an idiom not in the list) is invisible
///   to the ratchet. So "a new artifact-writing generator cannot be added
///   without one of the two deliberate acts" is **not** what this enforces. What
///   it enforces is that the *common* way of writing a file — the way every
///   generator in this repo currently writes — cannot be added in silence.
/// - A **false positive** (a `from_path(` that opens something for reading) makes
///   a non-writer ask the question anyway, which costs a `LEDGER_EXEMPT` row and
///   nothing else. That asymmetry is the right way round for a floor.
///
/// The semantic question — "did this write get built from a complete pass?" — is
/// what no grep can answer and what the ledger exists to carry.
fn appears_to_write_an_artifact(sources: &[(PathBuf, String)]) -> bool {
    sources
        .iter()
        .any(|(_, text)| WRITE_IDIOMS.iter().any(|idiom| text.contains(idiom)))
}

/// Whether a generator's **own entry point** looks like it routes a population
/// through [`CaptureLedger`] — an import of the module *and* a call to one of
/// the closing methods, both in that one file.
///
/// This is still a substring scan and still cannot see whether the population
/// actually flows through the ledger. Two restrictions make it hard to satisfy by
/// accident, and both exist because of a demonstrated laundering path (#1551
/// review):
///
/// - **Entry point only, not the transitive `#[path]` includes.** Appending a
///   single comment mentioning `CaptureLedger` to the shared
///   `examples/common/recording.rs` previously marked *three* separate
///   generators as having adopted the ledger, and the shrink-only half then
///   demanded their [`LEDGER_EXEMPT`] rows be deleted as stale — removing three
///   generators from the list on the strength of one line in a file they merely
///   share. `recording.rs` is exactly the module the #1544/#1545 migrations will
///   touch, so that was a live hazard, not a hypothetical one.
/// - **A `use` line, plus the type name, plus a `finish` call.** A bare mention
///   of the module — in a comment, in a doc link, in a `//` note pointing a
///   reader at the type — is not evidence of adoption, and `.finish()` on its
///   own is the terminal call of half the builders in the ecosystem
///   (`csv::Writer`, `indicatif`, encoders), so either string alone is cheap to
///   produce by accident. Requiring an actual `use` declaration naming
///   `conformance::completeness`, the [`CaptureLedger`] type somewhere in the
///   file, *and* the call that is the whole point of the type (`finish` is the
///   gate; reading `counts()` instead is the printed-count pattern the type
///   replaces) means the cheapest way to satisfy this predicate is to actually
///   use it (#1551 review).
///
/// The cost of the entry-point restriction, stated: a generator that genuinely
/// delegates its accounting to a shared `common/` module reads as *not*
/// ledgered here and keeps its [`LEDGER_EXEMPT`] row. The `use`-line
/// restriction costs the same way, and it is worth naming the two shapes it
/// misses rather than leaving them to be discovered: a fully-qualified
/// `ferro_hgvs::conformance::completeness::CaptureLedger::new(..)` at the call
/// site with no import, and a nested import group broken across lines by
/// rustfmt. Both read as *not* ledgered. That is the direction to err in — a
/// row that is kept too long is a visible line of text, and a row deleted
/// wrongly is a generator that silently leaves the list.
fn appears_to_route_through_the_ledger(sources: &[(PathBuf, String)], entry_point: &Path) -> bool {
    let Some((_, text)) = sources.iter().find(|(path, _)| path == entry_point) else {
        return false;
    };
    let imports_the_module = text
        .lines()
        .map(str::trim)
        .any(|line| line.starts_with("use ") && line.contains("conformance::completeness"));
    let names_the_type = text.contains("CaptureLedger");
    let closes_a_ledger = text.contains(".finish()") || text.contains(".finish_with(");
    imports_the_module && names_the_type && closes_a_ledger
}

/// The predicate above is the half of the ratchet that *deletes* a row, so the
/// shapes it must not accept are worth pinning rather than leaving to the next
/// reader of the doc comment. The accidental-adoption case is a doc link plus an
/// unrelated builder's `.finish()` — both strings present, nothing imported.
#[test]
fn a_mention_and_an_unrelated_finish_do_not_read_as_ledger_adoption() {
    let entry_point = PathBuf::from("examples/plausible.rs");
    let route = |text: &str| {
        appears_to_route_through_the_ledger(
            &[(entry_point.clone(), text.to_string())],
            &entry_point,
        )
    };

    // A doc link naming the module and the type, and a `.finish()` belonging to
    // some other builder entirely. Every substring the old predicate looked for
    // is present.
    assert!(!route(
        "//! Accounting is described in `conformance::completeness` (`CaptureLedger`).\n\
         fn main() { let w = csv::Writer::from_path(p).unwrap(); w.finish().unwrap(); }\n"
    ));
    // An import with no closing call is a generator that reads its `counts()`
    // and writes anyway — the printed-count pattern the ledger replaces.
    assert!(!route(
        "use ferro_hgvs::conformance::completeness::CaptureLedger;\n\
         fn main() { let l = CaptureLedger::new(\"rows\"); dbg!(l.counts()); }\n"
    ));
    // Real adoption: an import, the type, and the gate.
    assert!(route(
        "use ferro_hgvs::conformance::completeness::CaptureLedger;\n\
         fn main() { let l = CaptureLedger::new(\"rows\"); let c = l.finish().unwrap(); }\n"
    ));
    // …and via the allowance path, which is the other closing call.
    assert!(route(
        "use ferro_hgvs::conformance::completeness::{Allowance, CaptureLedger};\n\
         fn main() { let l = CaptureLedger::new(\"rows\"); \
         let c = l.finish_with(Allowance::at_most(1, \"why\")).unwrap(); }\n"
    ));
}

/// A **heuristic** ratchet, not a proof. Read
/// [`appears_to_write_an_artifact`] and [`appears_to_route_through_the_ledger`]
/// before quoting what this test guarantees: both sides are substring scans over
/// source text, so this catches the common case and says nothing about the rest.
/// What it does buy is that the question cannot go unasked for a generator
/// written the way every generator in this repo is currently written.
#[test]
fn artifact_writers_use_the_capture_ledger_or_are_allowlisted() {
    let exempt: BTreeMap<&str, &str> = LEDGER_EXEMPT.iter().copied().collect();
    assert_eq!(
        exempt.len(),
        LEDGER_EXEMPT.len(),
        "LEDGER_EXEMPT has a duplicate entry",
    );

    let mut undeclared = Vec::new();
    let mut stale = Vec::new();
    let mut covered = BTreeSet::new();

    for source in generator_sources_on_disk() {
        let key = source.display().to_string();
        let sources = compiled_sources(&source);
        let writes = appears_to_write_an_artifact(&sources);
        let ledgered = appears_to_route_through_the_ledger(&sources, &source);
        match exempt.get(key.as_str()) {
            Some(_) => {
                covered.insert(key.clone());
                if !writes {
                    stale.push(format!(
                        "{key}: no longer matches any idiom in WRITE_IDIOMS, so it no \
                         longer looks like an artifact writer"
                    ));
                } else if ledgered {
                    stale.push(format!(
                        "{key}: its entry point now imports \
                         `conformance::completeness` and calls `finish`"
                    ));
                }
            }
            None if writes && !ledgered => undeclared.push(key),
            None => {}
        }
    }

    // Shrink-only, half one: an entry that has outlived its subject must go, or
    // the list stops describing the tree and starts excusing files nobody
    // checked. Reported before the additions so a migration's own PR is told to
    // delete its row rather than being told twice about the same file.
    let deleted: Vec<&str> = exempt
        .keys()
        .filter(|key| !covered.contains(**key))
        .copied()
        .collect();
    assert!(
        deleted.is_empty(),
        "LEDGER_EXEMPT names generators that no longer exist; delete these rows: {deleted:?}",
    );
    // …and half two, which is a *heuristic* observation and says so, because
    // acting on it deletes a row. Both predicates are substring scans; a row that
    // is deleted wrongly puts its generator outside this list with nothing left
    // to notice.
    assert!(
        stale.is_empty(),
        "LEDGER_EXEMPT is shrink-only, and these rows look like they are no longer \
         needed. Confirm the generator really does account for its population — \
         both checks below are substring heuristics, not usage analysis — and then \
         delete the row. If it does not, fix the generator rather than the list:\n  {}",
        stale.join("\n  "),
    );

    // A new artifact writer must make one of the two deliberate choices. Silence
    // is what this file exists to make harder — not impossible: a writer using an
    // idiom outside WRITE_IDIOMS is invisible here, and closing that would mean
    // resolving Rust rather than reading it.
    assert!(
        undeclared.is_empty(),
        "these generators look like they write an artifact without routing their \
         population through `ferro_hgvs::conformance::completeness::CaptureLedger`. \
         Either adopt the ledger, or add a row to LEDGER_EXEMPT in this file saying \
         why not: {undeclared:?}",
    );
}
