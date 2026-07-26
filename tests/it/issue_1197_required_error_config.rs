//! Issue #1197 — a dropped `ErrorConfig` must not be a silent lenient default.
//!
//! `NormalizeConfig::default().with_direction(d)` fills in every field the
//! caller does not name, and `error_config` defaults to **lenient**. That is
//! the shape that produced #1181: the CLI built an `ErrorConfig` from
//! `--error-mode`, passed it into `run_normalize`, and then constructed the
//! normalizer's config without it — so the flag was inert from the initial
//! commit through 678 commits and roughly five months, taking `--ignore` and
//! `--reject` with it. #1191 fixed that one call site; the *pattern* that
//! allowed the omission survived at every other entry point.
//!
//! The web service never had the bug, and not by discipline — by construction:
//! it builds the config with a struct literal, which must name every field, so
//! omitting `error_config` is a compile error.
//!
//! This module pins the generalisation of that property. Entry-point seams —
//! the CLI, the PyO3 bindings, the web service — must construct a
//! `NormalizeConfig` in a way that *forces* the error configuration to be
//! named: either [`NormalizeConfig::for_entry_point`], whose `error_config`
//! parameter cannot be omitted, or a struct literal.
//!
//! ## Why a source scan
//!
//! There is no live defect to regress against — every seam threads the config
//! correctly today (verified while investigating #1196). A behavioral test
//! would therefore pass before and after, proving nothing. The property being
//! added is *structural*: "the next entry point cannot get this wrong". A scan
//! is the only thing that can fail when it is violated, and it is the same
//! technique `error_code_audit.rs` already uses for the error-code registry.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

/// Entry-point seams: everything that turns an external request (a CLI
/// invocation, a Python call, an HTTP request) into a `Normalizer`. These are
/// exactly the places where an error mode arrives from outside and can be
/// dropped on the floor.
const ENTRY_POINT_PATHS: &[&str] = &["src/bin/ferro.rs", "src/python.rs", "src/service"];

/// Constructors that decide `error_config` without the caller naming it, so an
/// error mode arriving from outside is silently discarded. Using one of these
/// inside an [`ENTRY_POINT_PATHS`] file is the defect shape #1181 shipped with.
///
/// The preset builders matter as much as the defaulting ones:
/// `NormalizeConfig::lenient().with_direction(d)` reintroduces #1181 exactly,
/// while containing neither `default()` nor `new()`. `Normalizer::new` is the
/// documented idiomatic constructor and hardcodes a lenient config.
const ELIDING_CONSTRUCTORS: &[&str] = &[
    "NormalizeConfig::default()",
    "NormalizeConfig::new()",
    "NormalizeConfig::lenient()",
    "NormalizeConfig::strict()",
    "NormalizeConfig::silent()",
    "Normalizer::new(",
];

/// True when `needle` occurs in `haystack` at an identifier boundary — i.e. not
/// as the tail of a longer path. Without this, `Normalizer::new(` matches
/// `HgvsRsNormalizer::new(` (`src/service/tools/hgvs_rs.rs`), which is an
/// unrelated type and a false positive.
fn contains_at_identifier_boundary(haystack: &str, needle: &str) -> bool {
    let mut from = 0;
    while let Some(rel) = haystack[from..].find(needle) {
        let at = from + rel;
        let preceded_by_ident = haystack[..at]
            .chars()
            .next_back()
            .is_some_and(|c| c.is_alphanumeric() || c == '_' || c == ':');
        if !preceded_by_ident {
            return true;
        }
        from = at + 1;
    }
    false
}

/// Collect `path:line` for every construction in an entry-point file that
/// elides the error configuration. Whole-line comments are skipped (they may
/// legitimately name the pattern they warn about); `#[cfg(test)]` bodies are
/// deliberately left in scope, since a test that builds an entry-point config
/// the old way is also worth flagging.
fn defaulting_constructions() -> BTreeSet<String> {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let mut hits = BTreeSet::new();
    for rel in ENTRY_POINT_PATHS {
        scan(&manifest.join(rel), rel, &mut hits);
    }
    hits
}

/// Count the `.rs` files an [`ENTRY_POINT_PATHS`] entry actually resolves to.
/// A seam that has been renamed or moved would otherwise make the scan below
/// pass by scanning nothing at all.
fn rust_files_under(rel: &str) -> usize {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel);
    let mut count = 0;
    count_rust_files(&path, &mut count);
    count
}

fn count_rust_files(path: &Path, count: &mut usize) {
    if path.is_dir() {
        if let Ok(rd) = std::fs::read_dir(path) {
            for entry in rd.flatten() {
                count_rust_files(&entry.path(), count);
            }
        }
    } else if path.is_file() && path.extension().and_then(|e| e.to_str()) == Some("rs") {
        // `is_file` is the load-bearing half: a path's extension is `rs`
        // whether or not the file exists, so testing the extension alone would
        // happily "find" a seam that had been renamed out from under us.
        *count += 1;
    }
}

fn scan(path: &Path, rel: &str, hits: &mut BTreeSet<String>) {
    if !path.exists() {
        return;
    }
    if path.is_dir() {
        if let Ok(rd) = std::fs::read_dir(path) {
            for entry in rd.flatten() {
                let child = entry.path();
                let name = child.file_name().unwrap_or_default().to_string_lossy();
                scan(&child, &format!("{rel}/{name}"), hits);
            }
        }
        return;
    }
    if path.extension().and_then(|e| e.to_str()) != Some("rs") {
        return;
    }
    let Ok(text) = std::fs::read_to_string(path) else {
        return;
    };
    let lines: Vec<&str> = text.lines().collect();
    for (i, line) in lines.iter().enumerate() {
        let code = line.trim_start();
        // A comment may legitimately name the pattern it is warning about.
        if code.starts_with("//") {
            continue;
        }
        if ELIDING_CONSTRUCTORS
            .iter()
            .any(|c| contains_at_identifier_boundary(line, c))
        {
            hits.insert(format!("{rel}:{}", i + 1));
            continue;
        }
        // A `NormalizeConfig { .. }` struct literal is only safe because it
        // must name every field. `..Default::default()` re-opens the hole — and
        // it is the exact shape `for_entry_point` itself uses, so it would
        // otherwise read as idiomatic. Flag any literal whose body never
        // mentions `error_config`.
        if opens_a_struct_literal(line) && struct_literal_elides_error_config(&lines, i) {
            hits.insert(format!("{rel}:{}", i + 1));
        }
    }
}

/// True when `line` opens a `NormalizeConfig { .. }` *value*, as opposed to
/// naming the type in a declaration. `-> NormalizeConfig {` (a function return
/// type), `struct NormalizeConfig {` and `impl .. for NormalizeConfig {` all
/// end in the same two tokens and are not constructions.
fn opens_a_struct_literal(line: &str) -> bool {
    let needle = "NormalizeConfig {";
    let Some(at) = line.find(needle) else {
        return false;
    };
    if !contains_at_identifier_boundary(line, needle) {
        return false;
    }
    let before = line[..at].trim_end();
    if before.ends_with("->") {
        return false;
    }
    let head = line.trim_start();
    !(head.starts_with("struct ")
        || head.starts_with("pub struct ")
        || head.starts_with("impl ")
        || before.contains(" for "))
}

/// Walk a `NormalizeConfig {` literal opened on `start` to its balanced close,
/// reporting whether the body ever names `error_config`.
fn struct_literal_elides_error_config(lines: &[&str], start: usize) -> bool {
    let mut depth = 0i32;
    let mut names_error_config = false;
    for line in &lines[start..] {
        let code = line.trim_start();
        if !code.starts_with("//") && line.contains("error_config") {
            names_error_config = true;
        }
        for ch in line.chars() {
            match ch {
                '{' => depth += 1,
                '}' => depth -= 1,
                _ => {}
            }
        }
        if depth <= 0 {
            break;
        }
    }
    !names_error_config
}

/// An entry point must not build its `NormalizeConfig` from a defaulting
/// constructor — the lenient `error_config` it silently supplies is precisely
/// the bug #1181 shipped for five months.
#[test]
fn entry_points_do_not_build_a_defaulting_normalize_config() {
    let hits = defaulting_constructions();
    assert!(
        hits.is_empty(),
        "entry-point seams must construct `NormalizeConfig` in a way that forces the error \
         configuration to be named — `NormalizeConfig::for_entry_point(direction, error_config)` \
         or a struct literal. A defaulting constructor supplies a silent *lenient* \
         `error_config`, which is how `--error-mode` stayed inert for five months (#1181, \
         #1197). Offending sites:\n  {:?}",
        hits
    );
}

/// Every declared seam must actually resolve to Rust source. Without this, a
/// renamed or relocated entry point turns the scan above into a no-op that
/// passes because it looked at nothing — the precise failure mode a structural
/// test has to rule out before its green result means anything.
#[test]
fn every_declared_entry_point_path_resolves_to_rust_source() {
    for rel in ENTRY_POINT_PATHS {
        let found = rust_files_under(rel);
        assert!(
            found > 0,
            "entry-point seam {rel:?} matched no `.rs` files — it was renamed or moved, and the \
             scan is silently covering nothing. Update `ENTRY_POINT_PATHS` to the new location.",
        );
    }
}

/// The scanner must actually be able to see the pattern it forbids, otherwise
/// the test above could pass because the scan is broken rather than because
/// the tree is clean. Exercises the real [`scan`] function against a written
/// file rather than re-implementing the match inline — an inline check would
/// pass even if `scan` itself were broken.
#[test]
fn the_scan_recognises_a_defaulting_construction() {
    let dir = std::env::temp_dir().join(format!("ferro-1197-scan-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("create temp scan dir");
    let file = dir.join("seam.rs");
    std::fs::write(
        &file,
        "fn build(d: ShuffleDirection) -> NormalizeConfig {\n    \
         let config = NormalizeConfig::default().with_direction(d);\n    \
         // NormalizeConfig::default() in a comment must not be flagged\n    \
         config\n}\n",
    )
    .expect("write temp seam");

    let mut hits = BTreeSet::new();
    scan(&file, "seam.rs", &mut hits);
    std::fs::remove_dir_all(&dir).ok();

    assert_eq!(
        hits,
        BTreeSet::from(["seam.rs:2".to_string()]),
        "`scan` must flag the bare defaulting construction on line 2 and ignore the comment \
         naming the same pattern on line 3",
    );
}

/// Write `body` to a scratch `.rs` file, scan it, and return the hit lines.
fn scan_snippet(tag: &str, body: &str) -> BTreeSet<String> {
    let dir = std::env::temp_dir().join(format!("ferro-1197-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("create temp scan dir");
    let file = dir.join("seam.rs");
    std::fs::write(&file, body).expect("write temp seam");
    let mut hits = BTreeSet::new();
    scan(&file, "seam.rs", &mut hits);
    std::fs::remove_dir_all(&dir).ok();
    hits
}

/// The preset builders are the bypass a `default()`/`new()`-only list misses:
/// `NormalizeConfig::lenient()` reintroduces #1181 verbatim while containing
/// neither forbidden token.
#[test]
fn the_scan_flags_the_preset_builders_that_hardcode_an_error_config() {
    for preset in ["lenient", "strict", "silent"] {
        let body = format!("fn build(d: ShuffleDirection) -> NormalizeConfig {{\n    NormalizeConfig::{preset}().with_direction(d)\n}}\n");
        assert_eq!(
            scan_snippet(preset, &body),
            BTreeSet::from(["seam.rs:2".to_string()]),
            "`NormalizeConfig::{preset}()` at an entry point discards the caller's error mode \
             and must be flagged",
        );
    }
}

/// `Normalizer::new` hardcodes a lenient config, but the matcher must not trip
/// on an unrelated type whose name merely ends in `Normalizer`.
#[test]
fn the_scan_flags_normalizer_new_without_colliding_with_a_longer_type_name() {
    assert_eq!(
        scan_snippet(
            "norm-new",
            "fn a(p: P) -> Normalizer {\n    Normalizer::new(p)\n}\n"
        ),
        BTreeSet::from(["seam.rs:2".to_string()]),
        "`Normalizer::new` hardcodes a lenient error config and must be flagged",
    );
    assert!(
        scan_snippet(
            "norm-collide",
            "fn a(c: C) -> R {\n    HgvsRsNormalizer::new(&c)\n}\n"
        )
        .is_empty(),
        "`HgvsRsNormalizer::new` is an unrelated type — matching it would be a false positive \
         (it is live at src/service/tools/hgvs_rs.rs)",
    );
}

/// A struct literal is only a safeguard while it must name every field.
/// `..Default::default()` re-opens the hole silently.
#[test]
fn the_scan_flags_a_struct_literal_that_defaults_away_the_error_config() {
    let eliding = "fn a(d: ShuffleDirection) -> NormalizeConfig {\n    \
                   NormalizeConfig {\n        shuffle_direction: d,\n        \
                   ..Default::default()\n    }\n}\n";
    assert_eq!(
        scan_snippet("lit-elide", eliding),
        BTreeSet::from(["seam.rs:2".to_string()]),
        "a `NormalizeConfig` literal that defaults away `error_config` must be flagged",
    );

    let explicit = "fn a(d: ShuffleDirection, ec: ErrorConfig) -> NormalizeConfig {\n    \
                    NormalizeConfig {\n        shuffle_direction: d,\n        \
                    error_config: ec,\n        ..Default::default()\n    }\n}\n";
    assert!(
        scan_snippet("lit-explicit", explicit).is_empty(),
        "a literal that names `error_config` is exactly the safe shape and must not be flagged",
    );
}
