//! Every ruling id cited in Rust source must resolve, and must not contradict
//! the ledger's `status`.
//!
//! # Why a source scan
//!
//! A ruling record in `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`
//! is how this project records its reading of the spec where the spec conflicts
//! with itself. `hgvs_spec_normalization_tests.rs`'s `RULING_STATUSES` already
//! pins each record's `status`, so a record cannot be silently upgraded from
//! `undecided` to `decided` or deleted.
//!
//! What nothing pinned is the **prose around the citations**. A ruling gets
//! decided, the ledger and the pin both move in the same reviewable diff, and
//! the half-dozen doc comments that say "this is still open" go on saying it.
//! Those comments are the first thing a future session reads, and a stale one is
//! actively worse than no comment: it sends the reader off to re-argue a
//! question that has already been answered, which is exactly what happened to
//! `delins-merge-vs-individual-gap-two-or-more` (see
//! `delins_equal_vs_unequal_length_discriminator.rs` for the row that
//! re-litigation kept landing on).
//!
//! Only a scan can fail when that is violated, which is the same reasoning
//! `error_code_audit.rs` and `issue_1197_required_error_config.rs` use for their
//! registries.
//!
//! # What it caught on landing
//!
//! Three sites, all asserting "open"/"unsettled"/"undecided" about a record the
//! ledger records as `decided`, and all fixed in the commit that added this
//! file:
//!
//! ```text
//! src/normalize/merge.rs:1890                  adjudication-precedence-order  "itself unsettled"
//! src/normalize/merge.rs:5366                  delins-merge-vs-individual-gap-two-or-more  "is open"
//! tests/it/hgvs_spec_normalization_tests.rs:575 delins-merge-vs-individual-gap-two-or-more  "the undecided ruling"
//! ```
//!
//! The last one is the sharpest: that file pins the ruling `decided` 225 lines
//! below the comment calling it undecided.
//!
//! # Known drift this guard does NOT cover, deliberately
//!
//! The same stale prose also appears in fixture **JSON** `note` fields —
//! `tests/fixtures/case-harvest/cases.json` (three sites) and
//! `tests/fixtures/mutalyzer-normalize/cases.json` (one). They are left alone
//! here because those files are blessed conformance inputs and editing their
//! prose is a separate, reviewable change; widening this scan to `.json` is a
//! reasonable follow-up, not a silent extension.
//!
//! # Two scoping decisions, both deliberate
//!
//! **Status words are read on the citation's own line, and no further.** The
//! obvious improvement — a window of a line or two either side — is wrong here,
//! and measurably so: in `RULING_STATUSES` the records are listed one after
//! another, so a two-line window around `canonical-form-choice-when-both-legal`
//! (`decided`) picks up the word "Undecided" from the comment introducing the
//! *next* record and reports a contradiction that does not exist. Line-scoping
//! is conservative — it misses citations whose status word wrapped onto the
//! following line — but every miss is a citation nobody is misled by today,
//! whereas every false positive is a test that has to be argued with.
//!
//! **Existence is checked only where the text says "ruling".** Scanning every
//! backticked kebab-case token for ledger membership would flag ordinary prose
//! (`assert-then-flip`, `column-for-column`). Requiring the word "ruling" on the
//! line, or the explicit `rulings[...]` form, matches how the repo actually
//! cites a record and keeps the check free of judgement calls.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

/// Where the ledger lives, relative to the crate root.
const LEDGER_RELATIVE_PATH: &str = "tests/fixtures/grammar/hgvs_spec_normalization_overrides.json";

/// Roots scanned for citations, relative to the crate root.
const SCAN_ROOTS: &[&str] = &["src", "tests"];

/// This file cites every id and every status word by construction, so scanning
/// it would report itself. Excluded by path suffix.
const SELF_PATH: &str = "tests/it/ruling_citation_currency.rs";

/// Words that assert a record IS settled, and words that assert it is NOT.
///
/// Matched at identifier boundaries, which is what keeps the two lists from
/// colliding: "undecided" contains "decided" and "unsettled" contains
/// "settled", and in both cases the preceding character is alphabetic, so the
/// shorter word does not match inside the longer one.
const DECIDED_WORDS: &[&str] = &["decided", "settled"];

/// See [`DECIDED_WORDS`]. "open" earns its place here because that is the word
/// `merge.rs` actually used for a decided record; it is matched as a whole word
/// so `open-issues.md` and `opened` do not trip it.
const UNDECIDED_WORDS: &[&str] = &["undecided", "unsettled", "unresolved", "open"];

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// `id -> status`, read from the committed ledger.
fn ledger() -> BTreeMap<String, String> {
    let path = crate_root().join(LEDGER_RELATIVE_PATH);
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let value: serde_json::Value =
        serde_json::from_str(&text).unwrap_or_else(|e| panic!("parse {}: {e}", path.display()));
    let rulings = value
        .get("rulings")
        .and_then(|r| r.as_array())
        .unwrap_or_else(|| panic!("{} has no `rulings` array", path.display()));
    let map: BTreeMap<String, String> = rulings
        .iter()
        .map(|r| {
            let id = r
                .get("id")
                .and_then(|v| v.as_str())
                .unwrap_or_else(|| panic!("a ruling has no `id` in {}", path.display()));
            let status = r
                .get("status")
                .and_then(|v| v.as_str())
                .unwrap_or_else(|| panic!("ruling {id} has no `status`"));
            (id.to_string(), status.to_string())
        })
        .collect();
    assert_eq!(
        map.len(),
        rulings.len(),
        "duplicate ruling id in {}",
        path.display()
    );
    assert!(
        !map.is_empty(),
        "{} lists no rulings — the scan below would be vacuous",
        path.display()
    );
    map
}

/// Whether `needle` occurs in `haystack` with non-identifier characters on both
/// sides.
///
/// Ruling ids are kebab-case, so `-` counts as part of an identifier: without
/// that, `delins-codon-carve-out-gap-one` would be found inside a longer id.
fn contains_at_identifier_boundary(haystack: &str, needle: &str) -> bool {
    let is_word = |c: char| c.is_ascii_alphanumeric() || c == '_' || c == '-';
    let mut from = 0usize;
    while let Some(offset) = haystack[from..].find(needle) {
        let start = from + offset;
        let end = start + needle.len();
        let before_ok = haystack[..start]
            .chars()
            .next_back()
            .is_none_or(|c| !is_word(c));
        let after_ok = haystack[end..].chars().next().is_none_or(|c| !is_word(c));
        if before_ok && after_ok {
            return true;
        }
        from = start + 1;
    }
    false
}

/// One `.rs` line, with the path it came from and its 1-based number.
struct SourceLine {
    path: String,
    number: usize,
    text: String,
}

/// Every line of every `.rs` file under [`SCAN_ROOTS`], excluding this file.
fn source_lines() -> Vec<SourceLine> {
    let root = crate_root();
    let mut lines = Vec::new();
    for scan_root in SCAN_ROOTS {
        let dir = root.join(scan_root);
        assert!(
            dir.is_dir(),
            "scan root {} does not exist — the scan would be vacuous",
            dir.display()
        );
        collect(&dir, scan_root, &mut lines);
    }
    lines
}

fn collect(dir: &Path, rel: &str, out: &mut Vec<SourceLine>) {
    let entries =
        std::fs::read_dir(dir).unwrap_or_else(|e| panic!("read_dir {}: {e}", dir.display()));
    let mut paths: Vec<PathBuf> = entries
        .map(|e| {
            e.unwrap_or_else(|err| panic!("dir entry under {}: {err}", dir.display()))
                .path()
        })
        .collect();
    paths.sort();
    for path in paths {
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or_default()
            .to_string();
        let child_rel = format!("{rel}/{name}");
        if path.is_dir() {
            collect(&path, &child_rel, out);
            continue;
        }
        if path.extension().and_then(|e| e.to_str()) != Some("rs") || child_rel == SELF_PATH {
            continue;
        }
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
        for (index, line) in text.lines().enumerate() {
            out.push(SourceLine {
                path: child_rel.clone(),
                number: index + 1,
                text: line.to_string(),
            });
        }
    }
}

/// Ruling-id-shaped tokens: lowercase alphanumerics in three or more
/// hyphen-separated segments. Matches the shortest real id
/// (`adjudication-precedence-order`) without matching two-word prose.
fn looks_like_a_ruling_id(token: &str) -> bool {
    let segments: Vec<&str> = token.split('-').collect();
    segments.len() >= 3
        && segments.iter().all(|s| {
            !s.is_empty()
                && s.chars()
                    .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit())
        })
}

/// Id-shaped tokens on a line that cites a ruling — either by the explicit
/// `rulings[...]` form or by naming a backticked token near the word "ruling".
fn id_shaped_citations(line: &str) -> BTreeSet<String> {
    let mut found = BTreeSet::new();

    // `rulings[<id>]` — the explicit machine-readable form.
    let mut rest = line;
    while let Some(open) = rest.find("rulings[") {
        let after = &rest[open + "rulings[".len()..];
        match after.find(']') {
            Some(close) => {
                let token = after[..close].trim_matches('`');
                if looks_like_a_ruling_id(token) {
                    found.insert(token.to_string());
                }
                rest = &after[close..];
            }
            None => break,
        }
    }

    // A backticked token on a line that says "ruling".
    if line.to_ascii_lowercase().contains("ruling") {
        for (index, chunk) in line.split('`').enumerate() {
            if index % 2 == 1 && looks_like_a_ruling_id(chunk) {
                found.insert(chunk.to_string());
            }
        }
    }
    found
}

/// Whole-word status claims on one line.
fn status_words(line: &str) -> (bool, bool) {
    let lower = line.to_ascii_lowercase();
    let says_decided = DECIDED_WORDS
        .iter()
        .any(|w| contains_at_identifier_boundary(&lower, w));
    let says_undecided = UNDECIDED_WORDS
        .iter()
        .any(|w| contains_at_identifier_boundary(&lower, w));
    (says_decided, says_undecided)
}

/// The scan is not vacuous: it read a realistic number of files and found the
/// ledger.
///
/// Without this, every assertion below passes trivially if the roots move or the
/// extension filter breaks — the failure mode this repo's `CLAUDE.md` calls a
/// structural zero.
#[test]
fn the_scan_reads_the_tree_it_claims_to() {
    let ledger = ledger();
    assert!(
        ledger.len() >= 10,
        "expected at least the 10 committed rulings, found {}",
        ledger.len()
    );
    let lines = source_lines();
    assert!(
        lines.len() > 100_000,
        "only {} source lines scanned — the scan roots or the `.rs` filter are wrong",
        lines.len()
    );
    let cited: BTreeSet<&String> = ledger
        .keys()
        .filter(|id| {
            lines
                .iter()
                .any(|l| contains_at_identifier_boundary(&l.text, id))
        })
        .collect();
    assert!(
        cited.len() >= 5,
        "only {} of {} ruling ids are cited in Rust source; the matcher is probably broken",
        cited.len(),
        ledger.len()
    );
}

/// (a) Every id cited as a ruling resolves in the ledger.
///
/// Catches a typo, and catches a record being deleted or renamed out from under
/// a comment that still points at it.
#[test]
fn every_cited_ruling_id_exists_in_the_ledger() {
    let ledger = ledger();
    let mut unknown: Vec<String> = Vec::new();
    for line in source_lines() {
        for token in id_shaped_citations(&line.text) {
            if !ledger.contains_key(&token) {
                unknown.push(format!("{}:{} cites `{token}`", line.path, line.number));
            }
        }
    }
    assert!(
        unknown.is_empty(),
        "these citations name a ruling the ledger does not have — fix the id, or restore the \
         record in {LEDGER_RELATIVE_PATH}:\n  {}",
        unknown.join("\n  ")
    );
}

/// (b) No citation contradicts the ledger's `status`.
///
/// Only lines that make a status claim are judged; a citation that merely names
/// a record is left alone. See the module docs for why the window is one line.
#[test]
fn no_ruling_citation_contradicts_the_ledger_status() {
    let ledger = ledger();
    let mut wrong: Vec<String> = Vec::new();
    for line in source_lines() {
        for (id, status) in &ledger {
            if !contains_at_identifier_boundary(&line.text, id) {
                continue;
            }
            let (says_decided, says_undecided) = status_words(&line.text);
            let claimed = match (says_decided, says_undecided) {
                (false, false) => continue,
                (true, true) => {
                    wrong.push(format!(
                        "{}:{} cites `{id}` on a line claiming both decided and undecided:\n      {}",
                        line.path,
                        line.number,
                        line.text.trim()
                    ));
                    continue;
                }
                (true, false) => "decided",
                (false, true) => "undecided",
            };
            if claimed != status {
                wrong.push(format!(
                    "{}:{} says `{id}` is {claimed}, ledger says {status}:\n      {}",
                    line.path,
                    line.number,
                    line.text.trim()
                ));
            }
        }
    }
    assert!(
        wrong.is_empty(),
        "these citations contradict {LEDGER_RELATIVE_PATH}. Update the prose (or the ledger, if \
         the ledger is what is wrong):\n  {}",
        wrong.join("\n  ")
    );
}
