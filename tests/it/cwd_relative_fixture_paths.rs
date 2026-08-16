//! Guard against cwd-relative fixture reads (#2019).
//!
//! A test's working directory is a property of its **package**, not of the
//! repository. `cargo nextest` runs a test from the root of the package that
//! owns it, so a module compiled into the root `ferro-hgvs` package's `it`
//! target runs from the workspace root, while the *same* module `#[path]`-
//! included by the `ferro-hgvs-soak-tests` member
//! (`tests-soak/tests/soak/main.rs`) runs from `tests-soak/`. A fixture read
//! against a bare relative path therefore works in `it` and breaks the moment
//! the module is pulled into the soak driver — and it breaks by reporting the
//! file missing, which reads as an environment problem rather than a path bug.
//!
//! The class had been rediscovered three times — `common/bulk_fixtures.rs` and
//! `common/fixture_gen.rs` during #1977, and `clinvar_hgvs_tests.rs` in CI
//! *after* #2001's push, on a module that had been correct in `it` for a long
//! time and became wrong only when compiled into the driver, with no edit to
//! itself. The fix each time is the same: resolve the fixture through
//! [`crate::common::fixture_gen::fixture_path`], which anchors on
//! `CARGO_MANIFEST_DIR` and ascends to the workspace `Cargo.lock`, so it is
//! correct from either package. This guard stops instance four: it scans every
//! `tests/it` source and fails on any filesystem read whose direct argument is
//! a `"tests/fixtures/…"` string literal.
//!
//! # What it catches, and what it deliberately does not
//!
//! It matches a `"tests/fixtures/…"` literal passed **directly** to a
//! filesystem-accessing call — `fs::read_to_string("…")`, `File::open("…")`,
//! and the other accessors in [`FS_ACCESSORS`]. That is exactly the shape the
//! three prior instances took, and the shape a new test is most likely to
//! reintroduce by copy-paste.
//!
//! It is a token scan, not a parse, so it is a floor rather than a proof. A
//! path bound to a variable and read one line later (`let p =
//! Path::new("tests/fixtures/…"); fs::read_to_string(p)`), or composed at run
//! time (`format!("tests/fixtures/{name}")`), is invisible to it — the fs call
//! sees an identifier, not the literal. Anchoring those correctly is still
//! required; this guard simply cannot see them. What it does buy is that the
//! direct form — the one that keeps recurring — cannot be added in silence, and
//! the question gets asked. Widening it to chase indirection would trade a
//! precise, zero-false-positive check for a heuristic one.

use std::collections::BTreeSet;
use std::path::PathBuf;

/// This file. Excluded from the scan because it carries the very literals it
/// searches for — the `FS_ACCESSORS` names beside `"tests/fixtures/` fragments
/// in the unit test below — and would otherwise flag itself. Same reasoning as
/// `generated_fixture_ci_wiring.rs`'s and `oracle_exclude_invariant.rs`'s
/// self-exclusions.
const SELF: &str = "cwd_relative_fixture_paths.rs";

/// The prefix that marks a repository-relative fixture path literal.
const FIXTURE_PREFIX: &str = "tests/fixtures/";

/// The final path segment of a call that touches the filesystem. A
/// `"tests/fixtures/…"` literal handed directly to one of these resolves
/// against the process's working directory, which is the bug this guard exists
/// to prevent. Matching on the final segment lets one entry cover every spelling
/// of the same call — `fs::read`, `std::fs::read` and a bare imported `read` all
/// end in `read`; `File::open` and `OpenOptions::new().open` both end in `open`.
///
/// Deliberately excluded, because they do not read against the cwd: `new`
/// (`Path::new`/`PathBuf::from` build a path without touching disk), `join`,
/// `push`, `fixture_path` (the anchored resolver itself), `present_or_skip` and
/// `load` (helpers that anchor internally), and `format!` (composes a string).
const FS_ACCESSORS: &[&str] = &[
    "read_to_string",
    "read_to_end",
    "read",
    "read_dir",
    "read_link",
    "open",
    "create",
    "create_new",
    "write",
    "metadata",
    "canonicalize",
    "remove_file",
    "copy",
    "rename",
];

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// A single cwd-relative fixture read: the file it lives in, the 1-based line
/// number, and the offending literal.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Finding {
    file: String,
    line: usize,
    literal: String,
}

/// Whether `token` — the identifier immediately preceding the `(` of a call —
/// names a filesystem accessor, judged by its final `::`-separated segment.
fn is_fs_accessor(token: &str) -> bool {
    let last_segment = token.rsplit("::").next().unwrap_or(token);
    FS_ACCESSORS.contains(&last_segment)
}

/// Every cwd-relative fixture read in `text`, one [`Finding`] per site.
///
/// The scan walks each `"tests/fixtures/…"` literal, steps back over whitespace
/// to the character before it, and — when that character is the `(` of a call —
/// reads the call's identifier and asks [`is_fs_accessor`]. A literal that is
/// not the direct argument of a call (a `const`/`let` binding, or an argument to
/// `Path::new`, `fixture_path`, `format!`, …) has something other than `(`
/// there and is skipped.
///
/// The walk is over the whole text rather than line-by-line, so it sees through
/// `rustfmt`'s wrapping: `fs::read_to_string(` on one line and its
/// `"tests/fixtures/…"` argument on the next is the same finding as the
/// one-liner. It steps back over intervening whitespace *and* newlines to reach
/// the `(`.
fn cwd_relative_fixture_reads(file: &str, text: &str) -> Vec<Finding> {
    let mut findings = Vec::new();
    let needle = format!("\"{FIXTURE_PREFIX}");
    let mut search_from = 0;
    while let Some(rel) = text[search_from..].find(&needle) {
        let quote_at = search_from + rel;
        search_from = quote_at + needle.len();

        // Everything before the opening quote, with trailing whitespace (spaces,
        // newlines, indentation) stripped so a wrapped argument reads the same
        // as an inline one.
        let before = text[..quote_at].trim_end();
        let Some(open_paren) = before.strip_suffix('(') else {
            continue; // not `call(… "tests/fixtures/…"` — a binding or a non-call arg.
        };
        // The call identifier is the run of path characters immediately before
        // the `(` (again tolerating whitespace between the two).
        let ident: String = open_paren
            .trim_end()
            .chars()
            .rev()
            .take_while(|c| c.is_alphanumeric() || *c == '_' || *c == ':')
            .collect::<String>()
            .chars()
            .rev()
            .collect();
        if !is_fs_accessor(&ident) {
            continue;
        }
        // Recover the full literal (to the closing quote) for the message.
        let after_quote = &text[quote_at + 1..];
        let literal = after_quote
            .find('"')
            .map(|end| &after_quote[..end])
            .unwrap_or(after_quote);
        // 1-based line of the opening quote.
        let line = text[..quote_at].bytes().filter(|&b| b == b'\n').count() + 1;
        findings.push(Finding {
            file: file.to_string(),
            line,
            literal: literal.to_string(),
        });
    }
    findings
}

/// `(file name, contents)` for every `.rs` source under `tests/it`, excluding
/// this file.
fn integration_test_sources() -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut stack = vec![repo_root().join("tests/it")];
    while let Some(dir) = stack.pop() {
        let entries = std::fs::read_dir(&dir)
            .unwrap_or_else(|e| panic!("read {}: {e}", dir.display()))
            .filter_map(Result::ok);
        for entry in entries {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
                continue;
            }
            if path.extension().is_none_or(|e| e != "rs") {
                continue;
            }
            let name = path
                .file_name()
                .expect("a file has a name")
                .to_string_lossy()
                .into_owned();
            if name == SELF {
                continue;
            }
            let text = std::fs::read_to_string(&path)
                .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
            out.push((name, text));
        }
    }
    out.sort();
    out
}

/// No integration-test module may read a fixture against a cwd-relative literal.
///
/// The remedy for any finding is to route the literal through
/// `crate::common::fixture_gen::fixture_path`, which resolves the same path from
/// either package. See the module docs for why this matters and for what the
/// scan cannot see.
#[test]
fn no_test_reads_a_fixture_against_the_cwd() {
    let offenders: Vec<Finding> = integration_test_sources()
        .iter()
        .flat_map(|(name, text)| cwd_relative_fixture_reads(name, text))
        .collect();

    assert!(
        offenders.is_empty(),
        "these fixture reads resolve against the process working directory and \
         will break when the module is compiled into `ferro-hgvs-soak-tests` \
         (which runs from `tests-soak/`, not the workspace root). Route each \
         through `crate::common::fixture_gen::fixture_path(\"…\")`:\n{}",
        offenders
            .iter()
            .map(|f| format!("  {}:{}  {}", f.file, f.line, f.literal))
            .collect::<Vec<_>>()
            .join("\n"),
    );
}

/// The detector must recognise the shapes it is meant to catch and ignore the
/// shapes that are already correct — a guard whose matcher silently stops
/// matching would pass vacuously (see `generated_fixture_ci_wiring.rs`'s minimum
/// floor for the same concern). This pins both directions against a synthetic
/// sample rather than against the live tree, so it keeps its meaning after the
/// tree is swept clean.
#[test]
fn the_detector_matches_direct_reads_and_only_those() {
    let sample = r#"
        let a = fs::read_to_string("tests/fixtures/a.json").unwrap();
        let b = std::fs::read_to_string("tests/fixtures/b.json").unwrap();
        let c = File::open("tests/fixtures/c.gz").unwrap();
        // rustfmt-wrapped direct read — fs call and its literal on separate
        // lines. Must STILL be flagged:
        let m = fs::read_to_string(
            "tests/fixtures/m.json",
        ).unwrap();
        // correctly anchored — must NOT be flagged:
        let d = fs::read_to_string(fixture_path("tests/fixtures/d.json")).unwrap();
        let e = fs::read_to_string(crate::common::fixture_gen::fixture_path("tests/fixtures/e.json"));
        let f = Path::new("tests/fixtures/f.json");
        let g = format!("tests/fixtures/{name}");
        const H: &str = "tests/fixtures/h.json";
        buf.push("tests/fixtures/i.json");
    "#;

    let flagged: BTreeSet<String> = cwd_relative_fixture_reads("sample.rs", sample)
        .into_iter()
        .map(|f| f.literal)
        .collect();

    let expected: BTreeSet<String> = [
        "tests/fixtures/a.json",
        "tests/fixtures/b.json",
        "tests/fixtures/c.gz",
        "tests/fixtures/m.json",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();

    assert_eq!(flagged, expected);
}
