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
//! # The other direction
//!
//! This file judges prose that *already cites* a record. It is blind to the more
//! common failure: someone reads a spec clause, forms a view from it, and never
//! learns a record governs that clause — so there is no citation to judge.
//! `tests/it/clause_ruling_index.rs` covers that direction, with a clause-to-
//! record index rendered into its own module docs. Both files read the ledger
//! through `tests/it/common/rulings.rs`, so "what a record is" is defined once.
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

use crate::common::rulings::{self, looks_like_a_record_id, LEDGER_RELATIVE_PATH};

/// Roots scanned for citations, relative to the crate root.
const SCAN_ROOTS: &[&str] = &["src", "tests"];

/// This file cites every id and every status word by construction, so scanning
/// it would report itself. Excluded by its scan-relative path, which
/// [`the_scan_reads_the_tree_it_claims_to`] asserts both exists and matches.
const SELF_PATH: &str = "tests/it/ruling_citation_currency.rs";

/// Floor for how many lines must make a status claim about a cited record.
///
/// Measured at **84** on the commit that added this constant. The floor is set
/// at less than half of that, so ordinary prose edits cannot trip it while a
/// matcher that stopped recognising claims still would.
const STATUS_CLAIMS_FLOOR: usize = 40;

/// Words that assert a record IS settled, and words that assert it is NOT.
///
/// Matched at identifier boundaries, which is what keeps the two lists from
/// colliding: "undecided" contains "decided" and "unsettled" contains
/// "settled", and in both cases the preceding character is alphabetic, so the
/// shorter word does not match inside the longer one.
const DECIDED_WORDS: &[&str] = &["decided", "settled"];

/// See [`DECIDED_WORDS`]. "open" earns its place here because that is the word
/// `merge.rs` actually used for a decided record; it is matched as a whole word
/// so `open-issues.md` and `opened` do not trip it, and via [`contains_as_prose`]
/// so `File::open(` does not either.
const UNDECIDED_WORDS: &[&str] = &["undecided", "unsettled", "unresolved", "open"];

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// `id -> status`, read from the committed ledger.
///
/// The parsing lives in `common::rulings` so that this scan and
/// `clause_ruling_index.rs` share one definition of a record.
fn ledger() -> BTreeMap<String, String> {
    rulings::statuses()
}

/// Whether `needle` occurs in `haystack` with non-identifier characters on both
/// sides.
///
/// Ruling ids are kebab-case, so `-` counts as part of an identifier: without
/// that, `delins-codon-carve-out-gap-one` would be found inside a longer id.
fn identifier_boundary_matches(haystack: &str, needle: &str) -> Vec<usize> {
    let is_word = |c: char| c.is_ascii_alphanumeric() || c == '_' || c == '-';
    let mut found = Vec::new();
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
            found.push(start);
        }
        from = start + 1;
    }
    found
}

fn contains_at_identifier_boundary(haystack: &str, needle: &str) -> bool {
    !identifier_boundary_matches(haystack, needle).is_empty()
}

/// Whether `needle` occurs as prose rather than as part of a code identifier.
///
/// `open` is the hazard, being both a status word and one of the most common
/// method names there is: `File::open(`, `.open()` and `provider.open` all
/// contain it at an identifier boundary, so a line carrying such a call *and* a
/// ruling id would be read as claiming that record is unsettled. No line in the
/// tree does today, and this keeps it that way.
///
/// Deliberately narrower than "only look at comments", which was the other
/// option: several status claims in this repo live inside assertion messages and
/// `const` doc strings, and a comments-only scan would stop judging them.
fn contains_as_prose(haystack: &str, needle: &str) -> bool {
    identifier_boundary_matches(haystack, needle)
        .into_iter()
        .any(|start| {
            let before = haystack[..start].trim_end();
            let called = haystack[start + needle.len()..].starts_with('(');
            let field_or_path = before.ends_with('.') || before.ends_with("::");
            !called && !field_or_path
        })
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
                if looks_like_a_record_id(token) {
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
            if index % 2 == 1 && looks_like_a_record_id(chunk) {
                found.insert(chunk.to_string());
            }
        }
    }
    found
}

/// Whole-word status claims on one line.
fn status_words(line: &str) -> (bool, bool) {
    let lower = line.to_ascii_lowercase();
    let says_decided = DECIDED_WORDS.iter().any(|w| contains_as_prose(&lower, w));
    let says_undecided = UNDECIDED_WORDS.iter().any(|w| contains_as_prose(&lower, w));
    (says_decided, says_undecided)
}

/// A status word used as a code identifier is not a status claim.
///
/// Pins both directions of [`contains_as_prose`], because the risk runs both
/// ways: too loose and `File::open` on a citing line reports a contradiction
/// that does not exist; too tight and real claims stop being judged.
#[test]
fn a_status_word_used_as_code_is_not_a_status_claim() {
    for line in [
        "    let text = File::open(&path)?;",
        "    provider.open();",
        "        .open()",
    ] {
        assert_eq!(
            status_words(line),
            (false, false),
            "{line} uses a status word as code, not as a claim"
        );
    }
    // Prose is still judged — including inside a string literal, where several
    // of this repo's real claims live.
    assert_eq!(status_words("// the question is open"), (false, true));
    assert_eq!(
        status_words("\"...is decided, scoped to...\""),
        (true, false)
    );
    assert_eq!(status_words("// nothing to say here"), (false, false));
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

    // The self-exclusion must fire, and must not be vacuous. This file names
    // every id and every status word by construction, so if `SELF_PATH` stopped
    // matching, the two checks below would report its own examples as real
    // contradictions — a confusing failure pointing at the guard rather than at
    // the drift.
    assert!(
        crate_root().join(SELF_PATH).is_file(),
        "{SELF_PATH} is not there — this file moved, so its self-exclusion now excludes nothing"
    );
    assert!(
        !lines.iter().any(|l| l.path == SELF_PATH),
        "{SELF_PATH} was not excluded from the scan"
    );

    // The *status* half must not be vacuous either. Only lines that make a
    // status claim are judged, so a matcher that stopped recognising claims
    // would leave `no_ruling_citation_contradicts_the_ledger_status` passing by
    // never judging anything — which is exactly the risk carried by narrowing
    // the match to prose (see `contains_as_prose`).
    let judged = lines
        .iter()
        .filter(|l| {
            status_words(&l.text) != (false, false)
                && ledger
                    .keys()
                    .any(|id| contains_at_identifier_boundary(&l.text, id))
        })
        .count();
    assert!(
        judged >= STATUS_CLAIMS_FLOOR,
        "only {judged} lines make a status claim about a cited record, below the floor of \
         {STATUS_CLAIMS_FLOOR} — the status matcher is probably broken, and the contradiction \
         check would pass by judging nothing"
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
        // Hoisted out of the per-id loop below: the status claim is a property of
        // the line, not of the record, so computing it per id repeated the whole
        // word scan once per ledger entry.
        let (says_decided, says_undecided) = status_words(&line.text);
        if (says_decided, says_undecided) == (false, false) {
            continue;
        }
        for (id, status) in &ledger {
            if !contains_at_identifier_boundary(&line.text, id) {
                continue;
            }
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
