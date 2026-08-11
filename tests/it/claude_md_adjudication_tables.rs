//! `CLAUDE.md`'s two adjudication tables must partition the ledger by status.
//!
//! # Why this exists
//!
//! `CLAUDE.md` lists the ruling records in two tables — one for the open
//! questions, one for the decided rulings and the scope each was decided at.
//! Which table a record sits in is a **status claim**, and it was hand-maintained.
//!
//! It drifted three records deep (#1584): `codon-carve-out-shape-restriction`,
//! `exon-junction-dup-converge-from-the-far-side` and
//! `separation-is-a-property-of-the-spelling-not-of-the-variant` had all been
//! ruled on while the table still carried them as live questions. A fourth,
//! `ring-telomere-anchoring`, appeared in neither table at all.
//!
//! (The wording above is deliberately status-word-free next to each id.
//! `ruling_citation_currency` scans this file and reads an id sitting on a line
//! with contradicting status words as a false claim — which it did, on the first
//! draft of this comment. Describing historical drift and asserting a current
//! status look identical to a line-based scan, so keep the two apart rather than
//! adding an exclusion, which would blind the scan to this file for good.)
//!
//! That is the more expensive direction of the two. A stale *count* reads as a
//! claim about coverage; a stale *status* reads as an open question that is in
//! fact settled, and the next reader re-derives a decision that has already been
//! made. It came within one step of doing exactly that: the
//! `separation-is-a-property-of-the-spelling-not-of-the-variant` record's **id
//! states its question, and its ruling is the opposite of it**, so a reader who
//! trusts the stale table and reads no further gets the answer backwards.
//!
//! # Why the existing currency scan could not catch it
//!
//! `ruling_citation_currency` already compares status claims against the ledger,
//! and it would not have helped here for two independent reasons:
//!
//! 1. It scans `.rs` files under `SCAN_ROOTS = ["src", "tests"]`. `CLAUDE.md` is
//!    Markdown, at the repository root — outside that twice over.
//! 2. It is **line-based**. In these tables the status word ("open") sits in the
//!    introducing paragraph while the record ids sit in table rows, so no single
//!    line carries both an id and a status word. Widening the scan's roots or
//!    extensions would still have found nothing.
//!
//! So the check has to be *section-aware*, which is what this file adds. The two
//! guards are complements, not duplicates: that one judges prose anywhere, this
//! one judges membership of two specific tables.
//!
//! # What it does not claim
//!
//! Only membership. Nothing here reads the prose in the second column, so a row
//! can sit in the right table and still describe the ruling wrongly — that is
//! what review is for. What it does buy is that a record cannot be *silently* on
//! the wrong side of the line, and that a newly added record cannot be omitted.

use std::collections::BTreeSet;
use std::path::PathBuf;

use crate::common::rulings;

/// The header row that opens the open-questions table.
const OPEN_TABLE_HEADER: &str = "| record | what is open |";

/// The header row that opens the decided-rulings table.
const DECIDED_TABLE_HEADER: &str = "| record | ruling |";

fn claude_md() -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("CLAUDE.md");
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

/// The backticked ids in the first column of the table opened by `header`, in
/// document order and **with duplicates kept**.
///
/// Reads from the header to the first line that is not a table row, so adding
/// prose after a table cannot silently extend it.
///
/// Returning a `Vec` rather than a set is deliberate. Collecting straight into a
/// `BTreeSet` discards a repeated row, so one record listed twice in the same
/// table compared equal to the ledger and passed — under a test named
/// `every_ledger_record_appears_exactly_once`. `unique_ids` is what makes the
/// name true.
fn ids_in_table(markdown: &str, header: &str) -> Vec<String> {
    let mut lines = markdown.lines().skip_while(|line| line.trim() != header);
    assert!(
        lines.next().is_some(),
        "CLAUDE.md has no table opened by {header:?} — the guard would be vacuous. If the table \
         was renamed, rename it here too rather than deleting this assertion"
    );

    let mut ids = Vec::new();
    for line in lines {
        let line = line.trim();
        if !line.starts_with('|') {
            break;
        }
        if line.starts_with("|---") {
            continue;
        }
        // First cell, between the leading `|` and the next one.
        let first_cell = line
            .trim_start_matches('|')
            .split('|')
            .next()
            .unwrap_or_default();
        if let Some(id) = first_cell.split('`').nth(1) {
            ids.push(id.to_string());
        }
    }
    ids
}

/// The ids of `rows` as a set, refusing a table that lists one record twice.
///
/// A repeated row is a documentation defect on its own — two rows for one record
/// can carry two different descriptions of the ruling, and nothing says which is
/// current — and it is invisible to any set-versus-set comparison.
fn unique_ids(rows: Vec<String>, table: &str) -> BTreeSet<String> {
    let mut seen = BTreeSet::new();
    let mut repeated = Vec::new();
    for id in rows {
        if !seen.insert(id.clone()) {
            repeated.push(id);
        }
    }
    assert!(
        repeated.is_empty(),
        "the {table} table lists these records more than once, so one row of each pair is \
         unreachable prose that no comparison here can check: {repeated:?}"
    );
    seen
}

#[test]
fn the_two_tables_partition_the_ledger_by_status() {
    let markdown = claude_md();
    let open = unique_ids(ids_in_table(&markdown, OPEN_TABLE_HEADER), "open-questions");
    let decided = unique_ids(
        ids_in_table(&markdown, DECIDED_TABLE_HEADER),
        "decided-rulings",
    );

    let ledger = rulings::records();
    let mut want_open = BTreeSet::new();
    let mut want_decided = BTreeSet::new();
    for record in &ledger {
        match record.status.as_str() {
            "undecided" => want_open.insert(record.id.clone()),
            "decided" => want_decided.insert(record.id.clone()),
            other => panic!("record {} has unknown status {other:?}", record.id),
        };
    }

    // Non-vacuous in both directions: `0 of 0` must not pass as a result, and a
    // ledger that had become all-decided would silently empty one table.
    assert!(
        !want_open.is_empty() && !want_decided.is_empty(),
        "the ledger has no record on one side of the line, so this guard measures nothing"
    );

    assert_eq!(
        open,
        want_open,
        // `want_open` is built from the `undecided` records, so the second list is
        // records the ledger leaves OPEN. Calling them "decided" here would state
        // the opposite of the ledger — the precise failure this guard exists to
        // catch, reproduced in its own diagnostic, and it would send the reader to
        // move the row into the wrong table.
        "CLAUDE.md's open-questions table disagrees with the ledger.\n  \
         listed as open but decided in the ledger: {:?}\n  \
         undecided in the ledger but missing from the open table: {:?}\n\
         Moving a record between the two tables is part of deciding it.",
        open.difference(&want_open).collect::<Vec<_>>(),
        want_open.difference(&open).collect::<Vec<_>>(),
    );

    assert_eq!(
        decided,
        want_decided,
        "CLAUDE.md's decided-rulings table disagrees with the ledger.\n  \
         listed as decided but still open: {:?}\n  \
         decided in the ledger but missing from the table: {:?}",
        decided.difference(&want_decided).collect::<Vec<_>>(),
        want_decided.difference(&decided).collect::<Vec<_>>(),
    );
}

#[test]
fn every_ledger_record_appears_exactly_once() {
    let markdown = claude_md();
    let open = unique_ids(ids_in_table(&markdown, OPEN_TABLE_HEADER), "open-questions");
    let decided = unique_ids(
        ids_in_table(&markdown, DECIDED_TABLE_HEADER),
        "decided-rulings",
    );

    let overlap: Vec<_> = open.intersection(&decided).collect();
    assert!(
        overlap.is_empty(),
        "these records are in BOTH tables, so one of the two is a false status claim: {overlap:?}"
    );

    // The union check is what would have caught `ring-telomere-anchoring`, which
    // was in neither table — an omission that no per-table comparison sees if
    // you only look at the table you happen to be reading.
    let listed: BTreeSet<_> = open.union(&decided).cloned().collect();
    let ledger: BTreeSet<_> = rulings::records().iter().map(|r| r.id.clone()).collect();
    assert_eq!(
        listed,
        ledger,
        "CLAUDE.md and the ledger name different record sets.\n  \
         in CLAUDE.md but not the ledger: {:?}\n  \
         in the ledger but documented nowhere: {:?}",
        listed.difference(&ledger).collect::<Vec<_>>(),
        ledger.difference(&listed).collect::<Vec<_>>(),
    );
}

/// The parser's own guard, on synthetic input.
///
/// The two tests above read the real `CLAUDE.md`, so neither can exercise a
/// malformed table without editing a committed file. This one hands `ids_in_table`
/// a table it controls, which is the only way to pin that a repeated row is
/// *rejected* rather than quietly deduped — the defect the `BTreeSet` collect
/// hid, under a test whose name promised otherwise.
#[test]
#[should_panic(expected = "lists these records more than once")]
fn a_record_listed_twice_in_one_table_is_rejected() {
    let markdown = format!(
        "{OPEN_TABLE_HEADER}\n|---|---|\n| `alpha` | first |\n| `beta` | second |\n| `alpha` | \
         a second row for alpha, saying something else |\n"
    );
    let _ = unique_ids(ids_in_table(&markdown, OPEN_TABLE_HEADER), "open-questions");
}

/// And the negative control: a well-formed table of the same shape passes, so the
/// guard above cannot be satisfied by `ids_in_table` failing for some other reason.
#[test]
fn a_table_with_no_repeats_is_accepted() {
    let markdown =
        format!("{OPEN_TABLE_HEADER}\n|---|---|\n| `alpha` | first |\n| `beta` | second |\n");
    let ids = unique_ids(ids_in_table(&markdown, OPEN_TABLE_HEADER), "open-questions");
    assert_eq!(
        ids,
        ["alpha", "beta"].iter().map(|s| s.to_string()).collect(),
        "the parser must read both rows of a well-formed table"
    );
}
