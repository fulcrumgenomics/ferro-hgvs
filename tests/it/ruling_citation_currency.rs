//! Every ruling id cited in Rust source must resolve **and** must not contradict
//! the ledger's `status`; every id cited in the ledger's own prose must resolve.
//!
//! The asymmetry is deliberate and is stated rather than implied, because the
//! opening line used to claim both halves for both sources. A Rust doc comment
//! saying "this is still open" about a `decided` record is the failure this file
//! was built for, and the source scan checks exactly that. The ledger's own
//! rationales are discursive — they narrate what a record *used* to say, what a
//! withdrawn amendment argued, which neighbouring question a record does not
//! settle — so a status-claim matcher over that prose would fire on history
//! rather than on drift. What is checked there is resolution: an id named in a
//! rationale must be a record the ledger has, or a declared retired id.
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
//!
//! # The ledger's own cross-references
//!
//! The scan above reads `.rs` files, so it never looks at the ledger — and the
//! ledger cross-references itself constantly, because that is how one record
//! says which neighbouring question it does *not* settle.
//! [`every_record_to_record_citation_resolves`] closes that, with a looser
//! extractor ([`record_ids_named_in`]) because a rationale cites its siblings by
//! bare backticked id rather than with the word "ruling".
//!
//! # The ledger's citations of its own GUARDS
//!
//! Everything above is about record **ids**. A rationale also cites the test
//! that enforces its ruling — "Pinned by `defect_371_transcript_exit::the_exit_follows_…`"
//! — and until [`every_ruling_to_guard_citation_resolves`] nothing checked those
//! at all. Rename the guard and the record goes on naming the old one, reading
//! as enforced while enforcing nothing.
//!
//! That check is `#[ignore]`d because it is **red on two `decided` records** on
//! the base it landed on, and both are real rather than fixable by editing a
//! name: each rename inverted the claim, so substituting the new name would
//! change what the record asserts. Its doc comment carries the census and the
//! two rows.
//!
//! # What is deliberately NOT checked: whether a quoted sentence is really the
//! # record's
//!
//! The sibling failure is prose that puts words in a record's mouth — a doc
//! comment quoting a sentence and attributing it to a ruling, where the ledger
//! contains no such sentence. That happened: a module doc quoted
//! `canonical-form-choice-when-both-legal` as saying "Sequence-level equivalence
//! still needs an answer for consumers who dedupe … a groupable SPDI key, not
//! the canonical string", which is text from an unmerged branch's copy of the
//! ledger and appears nowhere in `main`'s.
//!
//! It is not checked here, for two measured reasons rather than as an oversight:
//!
//! 1. **Attribution is not mechanically decidable.** Prose in this repo quotes
//!    the *spec* far more often than it quotes a record, frequently in the same
//!    paragraph as a record citation. Deciding which of the two a quoted run
//!    belongs to is the judgement call the rest of this file is built to avoid.
//! 2. **The population is empty, so the check would be vacuous.** Measured on
//!    the tree at the commit that added this: blocks citing a record *and*
//!    carrying a quoted run of 40+ characters *and* naming no spec clause —
//!    **zero**. A guard over an empty population reads as coverage and provides
//!    none, which is the failure `the_index_is_not_vacuous` exists to name.
//!
//! If the population ever becomes non-empty, the missing ingredient is an
//! explicit attribution convention in the prose (a marker saying "this record
//! says"), not a cleverer heuristic over the text we have.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use crate::common::rulings::{self, looks_like_a_record_id, LEDGER_RELATIVE_PATH};

/// Roots scanned for citations, relative to the crate root.
const SCAN_ROOTS: &[&str] = &["src", "tests"];

/// This file cites every id and every status word by construction, so scanning
/// it would report itself. Excluded by its scan-relative path, which
/// [`the_scan_reads_the_tree_it_claims_to`] asserts both exists and matches.
const SELF_PATH: &str = "tests/it/ruling_citation_currency.rs";

/// Record ids that no longer exist, and the record that carries their question
/// now.
///
/// A rationale may name a **former** id on purpose — `inversion-vs-two-delins-76-83`
/// opens with "THE RECORD'S FORMER ID AND QUESTION WERE BOTH WRONG … It was
/// `inversion-vs-two-substitutions-76-83`", because the rename *is* the argument:
/// the old id asserted the competing members were substitutions, and they are
/// two-base `delins`.
///
/// Declaring those here rather than pattern-matching on words like "former"
/// keeps [`every_record_to_record_citation_resolves`] free of judgement calls,
/// and turns the retired id into a fact the ledger states rather than one a
/// reader has to infer. [`no_retired_id_is_also_a_live_record`] stops an entry
/// masking a live record.
const RETIRED_RECORD_IDS: &[(&str, &str)] = &[
    (
        "inversion-vs-two-substitutions-76-83",
        "renamed to `inversion-vs-two-delins-76-83`, which names its own former id \
         because the rename corrected the argument",
    ),
    // Not a rename and not a deletion: this id names a record that **never
    // existed on `main`**. `adjudication-precedence-order` cites it to record
    // that a 2026-08-08 amendment removed a rank "on the strength of
    // `partition-is-the-unit-of-normalization`, a record that has never existed
    // on `main`" — the citation IS the finding, so it must survive.
    //
    // Added after measuring the interop rather than after CI reported it: this
    // citation arrives with #1604, so on `main` today the check passes and on
    // `main` plus #1604 it fails. Neither PR's own CI can see that — each is
    // green alone. Declaring it here makes the merge order irrelevant.
    //
    // It is also why the second field is prose and not a replacement id: there
    // is no record that carries this question now, and a schema demanding one
    // could not express that.
    (
        "partition-is-the-unit-of-normalization",
        "never existed on `main`; cited by `adjudication-precedence-order` as the \
         non-existent record a withdrawn amendment was argued from",
    ),
];

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
    source_files()
        .into_iter()
        .flat_map(|(path, text)| {
            text.lines()
                .enumerate()
                .map(|(index, line)| SourceLine {
                    path: path.clone(),
                    number: index + 1,
                    text: line.to_string(),
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

/// Every `.rs` file under [`SCAN_ROOTS`] as `(scan-relative path, text)`,
/// excluding this file.
///
/// Split out from [`source_lines`] because [`guard_sites`] needs each file
/// *whole*: it brace-matches function bodies, and a flattened line list has
/// lost the file boundaries that makes possible. One walker rather than two, so
/// "what the scan reads" has a single definition — the same reasoning
/// `common::rulings` gives for parsing the ledger in one place.
fn source_files() -> Vec<(String, String)> {
    let root = crate_root();
    let mut files = Vec::new();
    for scan_root in SCAN_ROOTS {
        let dir = root.join(scan_root);
        assert!(
            dir.is_dir(),
            "scan root {} does not exist — the scan would be vacuous",
            dir.display()
        );
        collect(&dir, scan_root, &mut files);
    }
    files
}

fn collect(dir: &Path, rel: &str, out: &mut Vec<(String, String)>) {
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
        out.push((child_rel, text));
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

// ---------------------------------------------------------------------------
// The ledger's own cross-references
// ---------------------------------------------------------------------------

/// Every record id a **rationale** names, as the ledger spells them.
///
/// Deliberately not [`id_shaped_citations`]: that requires the word "ruling" on
/// the line, which is how *Rust prose* cites a record. A rationale cites its
/// siblings by bare backticked id — "see `delins-merge-vs-individual-gap-two-or-more`,
/// which remains decided" — so requiring the keyword would find almost nothing.
///
/// The looser rule is safe here because of what else appears in a rationale:
/// code identifiers are `snake_case` or `CamelCase` and fail the hyphen grammar,
/// and clause citations carry a `/` and a `.md:` and fail it too. Measured on the
/// ledger at the commit that added this: 8 tokens matched across 11 records, of
/// which 7 were live records and 1 was the declared former id below.
fn record_ids_named_in(rationale: &str) -> BTreeSet<String> {
    let mut found = BTreeSet::new();
    for (index, chunk) in rationale.split('`').enumerate() {
        if index % 2 == 1 && looks_like_a_record_id(chunk) {
            found.insert(chunk.to_string());
        }
    }
    found
}

/// (c) A record's rationale may not cite a record that does not exist.
///
/// The source scan above reads `.rs` files, so it never looks at the ledger —
/// and the ledger cross-references itself constantly, because that is how one
/// record says which neighbouring question it does *not* settle. A record
/// deleted or renamed leaves those references pointing at nothing, and nothing
/// noticed.
///
/// This is not hypothetical. PR #1569 rewrote the `decided` record
/// `adjudication-precedence-order` and removed a rank from the precedence order,
/// on the stated ground that `partition-is-the-unit-of-normalization` was
/// decided — a record that has never existed on `main`. The rewrite was caught
/// by a human reading the diff, twice, which is not a mechanism.
#[test]
fn every_record_to_record_citation_resolves() {
    let records = rulings::records();
    let live: BTreeSet<&str> = records.iter().map(|r| r.id.as_str()).collect();
    let retired: BTreeMap<&str, &str> = RETIRED_RECORD_IDS.iter().copied().collect();

    let mut unknown: Vec<String> = Vec::new();
    let mut checked = 0usize;
    for record in &records {
        for cited in record_ids_named_in(&record.rationale) {
            if cited == record.id {
                continue;
            }
            checked += 1;
            if !live.contains(cited.as_str()) && !retired.contains_key(cited.as_str()) {
                unknown.push(format!("`{}` cites `{cited}`", record.id));
            }
        }
    }

    assert!(
        unknown.is_empty(),
        "these rationales in {LEDGER_RELATIVE_PATH} name a record the ledger does not have. \
         Fix the id, restore the record, or — if the reference is to a former id on purpose — \
         add it to `RETIRED_RECORD_IDS` with the record that carries the question now:\n  {}",
        unknown.join("\n  ")
    );

    // Non-vacuity. The ledger cross-references itself by construction, so a run
    // that checked nothing means the extractor stopped recognising citations —
    // and this test would then pass forever by looking at an empty set. Measured
    // at 8 on the commit that added this; the floor is deliberately far below
    // that so ordinary prose edits cannot trip it.
    assert!(
        checked >= 4,
        "only {checked} record-to-record citations were found in {LEDGER_RELATIVE_PATH} — \
         `record_ids_named_in` is probably broken, and this check would pass by reading nothing"
    );
}

/// A retired id must not name a record that still exists.
///
/// Without this, `RETIRED_RECORD_IDS` is a way to silence the check above rather
/// than to declare a rename: adding a live id would make every stale citation of
/// it resolve through the exemption instead of through the ledger.
#[test]
fn no_retired_id_is_also_a_live_record() {
    let live: BTreeSet<String> = rulings::records().into_iter().map(|r| r.id).collect();
    let overlapping: Vec<&str> = RETIRED_RECORD_IDS
        .iter()
        .map(|(id, _)| *id)
        .filter(|id| live.contains(*id))
        .collect();
    assert!(
        overlapping.is_empty(),
        "these ids are listed as retired but are live records in {LEDGER_RELATIVE_PATH}: \
         {overlapping:?} — remove them from `RETIRED_RECORD_IDS`, which is for ids the ledger \
         no longer has"
    );

    // And the successor an entry NAMES must itself resolve.
    //
    // The second field is prose, so a successor id inside it was checked by
    // nothing — "renamed to `inversion-vs-two-delnis-76-83`" would read as a
    // complete declaration while pointing at a record that does not exist,
    // which is the same dangling-citation defect this file exists to close,
    // one level up. Every backticked record-shaped token in the note must be
    // live.
    //
    // Deliberately NOT a required field. `partition-is-the-unit-of-normalization`
    // never existed and has no successor, so a schema demanding one could only
    // be satisfied by inventing a record or by lying. An entry naming no
    // successor is legitimate; an entry naming a wrong one is not.
    let mut dangling: Vec<String> = Vec::new();
    for (retired_id, note) in RETIRED_RECORD_IDS {
        for named in record_ids_named_in(note) {
            if named != *retired_id && !live.contains(&named) {
                dangling.push(format!("`{retired_id}`'s note names `{named}`"));
            }
        }
    }
    assert!(
        dangling.is_empty(),
        "a `RETIRED_RECORD_IDS` note names a record the ledger does not have:\n  {}\n\
         Either the successor id is misspelled, or the record it names was itself removed.",
        dangling.join("\n  ")
    );
}

// --------------------------------------------------------------------------
// (d) A record's rationale may not cite a GUARD that does not exist.
//
// The three checks above are all about record *ids*. This one is about the
// other half of a record's prose: the test function it names as the thing that
// enforces its ruling.
// --------------------------------------------------------------------------

/// Words a guard name in this repo starts with.
///
/// Test names here are sentences — `a_span_outweighs_a_split_that_keeps_reference_bases`,
/// `the_exit_follows_the_intron_bases_not_the_strand` — while source functions
/// are verb or noun phrases (`coalesce_coding_frame_separation`,
/// `insertion_is_duplication`, `count_tandem_repeats`). That difference is what
/// [`looks_like_a_guard_name`] keys on, and it is the same move
/// [`looks_like_a_record_id`] makes: recognise a citation by its **shape**, so
/// the check still fires when the thing it names has been deleted.
///
/// A grammar that instead asked "is this token a test that exists" would be
/// vacuous in precisely the case it is for — delete the test and the token
/// stops matching, so nothing is reported.
const GUARD_NAME_OPENERS: &[&str] = &[
    "a_", "an_", "the_", "two_", "three_", "four_", "five_", "every_", "no_", "non_", "both_",
    "each_", "one_",
];

/// Guard names a rationale cites on purpose although the tree no longer has
/// them, with why.
///
/// The sibling of [`RETIRED_RECORD_IDS`], and it exists for the same reason: a
/// rationale may name a **former** guard because the rename is part of the
/// argument. Modelled on how `Representation-Change: none` works — declining is
/// a first-class answer, and what is rejected is *silence*.
///
/// It is deliberately EMPTY. Both citations that fail today are stale rather
/// than deliberate, and listing them here would be using the escape hatch to
/// silence the finding instead of to declare an intent — see this check's own
/// doc comment for the two, and note that
/// [`no_retired_guard_name_is_also_a_live_guard`] stops an entry masking a
/// guard that does exist.
const RETIRED_GUARD_NAMES: &[(&str, &str)] = &[];

/// Whether `token` has the shape of a guard name; see [`GUARD_NAME_OPENERS`].
///
/// Three conjuncts, each measured against the tree rather than chosen: the
/// sentence-opener, at least four words, and at least twenty characters. Over
/// `origin/main` at `d4552167` this selects 1,649 of the 9,652 test functions
/// and only **11** of the 7,689 source functions that are not also tests — and
/// none of those eleven is cited `::`-qualified by any record.
fn looks_like_a_guard_name(token: &str) -> bool {
    GUARD_NAME_OPENERS
        .iter()
        .any(|opener| token.starts_with(opener))
        && token.len() >= 20
        && token.matches('_').count() >= 3
        && token
            .chars()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_')
}

/// Every guard a **rationale** names, as a bare function name.
///
/// Only the `::`-qualified forms are read, which is the deliberate narrowing
/// that keeps this free of judgement calls. The ledger spells a guard citation
/// four ways, and all four qualify it:
///
/// ```text
/// it::cis_junction_crossing_shift::the_three_member_spelling_and_...
/// defect_371_transcript_exit::the_exit_follows_the_intron_bases_not_the_strand
/// tests/it/corpus_prohibited_inputs.rs::a_bare_transcript_intronic_position_...
/// normalize::tests::the_junction_exit_fold_declines_a_mixed_accession
/// ```
///
/// A **bare** backticked name is left alone. Rationales do use that form, but
/// so does ordinary prose about the code, and telling the two apart is the
/// judgement call this file is built to avoid. Reading only the qualified form
/// costs recall and buys a rule with no arguable cases; the qualified spelling
/// is also the one a reader can act on, since it names where to look.
fn guard_names_named_in(rationale: &str) -> BTreeSet<String> {
    let mut found = BTreeSet::new();
    for (index, chunk) in rationale.split('`').enumerate() {
        if index % 2 == 0 || !chunk.contains("::") {
            continue;
        }
        if let Some(last) = chunk.rsplit("::").next() {
            if looks_like_a_guard_name(last) {
                found.insert(last.to_string());
            }
        }
    }
    found
}

/// Where a guard is defined, and the two properties that decide whether it can
/// enforce anything.
struct GuardSite {
    path: String,
    line: usize,
    /// `#[ignore]`d, so it does not run in an ordinary CI job.
    ignored: bool,
    /// Its body contains at least one assertion.
    asserts: bool,
}

/// Whether a function body asserts.
///
/// Matches `assert` anywhere — which covers the `assert!`/`assert_eq!` macros
/// **and** this repo's asserting helpers (`assert_normalizes_preserving_in`,
/// `assert_seam_oracles`) — plus the diverging macros.
///
/// Helper calls are the reason the substring form is used rather than a macro
/// pattern. A stricter matcher that required a `!` reported
/// `a_third_member_clear_of_the_tract_keeps_the_duplication_reaching_its_five_prime_most_position`
/// as assertionless while its whole body is one `assert_normalizes_preserving_in`
/// call — a false positive on the exact shape this check exists to find, which
/// would have been reported as a defect.
fn body_asserts(body: &str) -> bool {
    body.contains("assert")
        || body.contains("panic!")
        || body.contains("unreachable!")
        || body.contains("todo!")
        || body.contains("unimplemented!")
}

/// `guard name -> where it is defined`, over [`SCAN_ROOTS`].
///
/// Indexed by bare function name because that is how a citation resolves: the
/// ledger's module qualifier is prose (`it::…`, `tests/it/….rs::…`,
/// `normalize::tests::…` are three spellings of the same idea), so matching on
/// it would reject correct citations for how they were punctuated.
fn guard_sites() -> BTreeMap<String, Vec<GuardSite>> {
    let mut sites: BTreeMap<String, Vec<GuardSite>> = BTreeMap::new();
    for (path, text) in source_files() {
        let lines: Vec<&str> = text.lines().collect();
        let mut ignored = false;
        for (index, line) in lines.iter().enumerate() {
            let trimmed = line.trim_start();
            if trimmed.starts_with("#[ignore") {
                ignored = true;
                continue;
            }
            let Some(name) = function_name_declared_on(trimmed) else {
                // Attributes and doc comments sit between `#[ignore]` and the
                // `fn`, so only a line that is neither resets the flag.
                if !trimmed.is_empty() && !trimmed.starts_with("#[") && !trimmed.starts_with("//") {
                    ignored = false;
                }
                continue;
            };
            if looks_like_a_guard_name(&name) {
                sites.entry(name).or_default().push(GuardSite {
                    path: path.clone(),
                    line: index + 1,
                    ignored,
                    asserts: body_asserts(&body_starting_at(&lines, index)),
                });
            }
            ignored = false;
        }
    }
    sites
}

/// The function name declared on `line`, if it declares one.
///
/// `pub fn` is accepted, so non-test declarations are in scope on purpose and
/// the remaining visibility and qualifier forms are listed for the same reason:
/// `src/hgvs/variant.rs`'s `pub(crate) fn non_flanking_genomic_insertion_anchor`
/// and `src/hgvs/parser/variant.rs`'s `pub(crate) fn non_spec_mosaic_form_error`
/// both satisfy [`looks_like_a_guard_name`], and omitting their prefix would
/// leave a citation of either reported as "defined nowhere" while the function
/// sits in the tree — the most misleading message this check can emit.
fn function_name_declared_on(line: &str) -> Option<String> {
    let after_fn = line
        .strip_prefix("fn ")
        .or_else(|| line.strip_prefix("pub fn "))
        .or_else(|| line.strip_prefix("pub(crate) fn "))
        .or_else(|| line.strip_prefix("pub(super) fn "))
        .or_else(|| line.strip_prefix("const fn "))
        .or_else(|| line.strip_prefix("unsafe fn "))
        .or_else(|| line.strip_prefix("async fn "))
        .or_else(|| line.strip_prefix("pub async fn "))?;
    let name: String = after_fn
        .chars()
        .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
        .collect();
    (!name.is_empty()).then_some(name)
}

/// The body of the function whose signature starts at `start`, by brace
/// matching.
fn body_starting_at(lines: &[&str], start: usize) -> String {
    let mut depth = 0usize;
    let mut opened = false;
    let mut body = String::new();
    for line in lines.iter().skip(start) {
        body.push_str(line);
        body.push('\n');
        for c in line.chars() {
            match c {
                '{' => {
                    depth += 1;
                    opened = true;
                }
                '}' => depth = depth.saturating_sub(1),
                _ => {}
            }
        }
        if opened && depth == 0 {
            break;
        }
    }
    body
}

/// (d) Every guard a **decided** record names must exist, run, and assert.
///
/// # The failure this closes
///
/// A `decided` record and the guard that enforces it are connected by nothing
/// but prose. Rename the guard and the record goes on naming the old one; the
/// ledger still reads as enforced, and nothing anywhere fails.
///
/// `ruling_citation_currency`'s other three checks cannot see it — they judge
/// record **ids**, and a guard name is not id-shaped. Nor can
/// `spec_equivalence_classes_converge`, whose reach is the four decided records
/// carrying an `equivalence_classes` array.
///
/// # Why this is `#[ignore]`d, and what un-ignoring it costs
///
/// It is **red on `origin/main` at `d4552167`**, on two `decided` records, and
/// both are real. Landing it armed would land a red suite; exempting them via
/// [`RETIRED_GUARD_NAMES`] would use the escape hatch to hide the finding.
/// Neither is fixable mechanically, because **both renames inverted the claim**
/// — substituting the new name would silently change what the record's sentence
/// asserts, which is an adjudication and not a citation repair:
///
/// | record | cites | the tree now has |
/// |---|---|---|
/// | `alignment-only-symbol-in-a-description` | `…_refused_for_dash_and_accepted_for_x` | `…_refused_in_every_mode_for_both_x_and_dash` |
/// | `derivation-may-not-be-bounded-by-the-inputs-spelling` | `a_dup_flush_against_a_del_is_re_derived_into_the_span` | `a_dup_flush_against_a_del_is_left_alone` |
///
/// The first looks like a plain rename — the guard was renamed *when the ruling
/// was implemented* (#1684), so the record cites the pre-fix name, which
/// describes the behaviour the ruling **overturned**. The second is
/// substantive: the record cites that row as reading "AGAINST THE BOUND rather
/// than for it" because it was re-derived into the span, and the guard standing
/// in its place asserts it is left alone.
///
/// This is the idiom `corpus_prohibited_inputs`'s own
/// `the_decided_target_is_a_mode_gated_refusal` already uses: `#[ignore]`d,
/// asserting the decided answer, with an issue as its acceptance criterion.
/// Here that issue is **#1881** — settling both citations and un-`#[ignore]`ing
/// this is what closes it.
///
/// # Census on `origin/main` at `d4552167`
///
/// 28 records, 25 `decided`, 3 `undecided`. **No structured guard field
/// exists** — `equivalence_classes` is the only structured pointer, on 4
/// decided records, and it is already enforced end to end (the generator bails
/// on an unknown class id; `spec_equivalence_classes_converge` asserts and
/// guards its own vacuity). Guard references are otherwise prose. Of the 14
/// qualified citations this check reads, 12 resolve; all 12 assert, none is
/// manifest-gated or skip-on-absent, and **none is `#[ignore]`d**.
///
/// That last one is worth stating rather than leaving to be inferred, because
/// the ledger *does* name an `#[ignore]`d guard and this check does not reach
/// it: `absolute-prohibition-enforcement-stage` names
/// `the_decided_target_is_a_mode_gated_refusal` — the #1630/#1627/#1628
/// deferral — by **bare** backticked name, and [`guard_names_named_in`] reads
/// only the `::`-qualified form. So the ignored-guard arm below is currently
/// exercised by no live citation; what keeps it from being dead code is
/// [`the_guard_scan_reads_the_tree_it_claims_to`], which asserts the attribute
/// is observable at all.
///
/// **Reach, which "14 citations" does not convey.** Of the 25 `decided`
/// records, **5** carry a `::`-qualified guard citation and are checked here;
/// 3 more name a guard only by bare backticked name and are skipped by the
/// narrowing above; and 17 name no guard-shaped token at all. This check
/// judges the citations that exist — it does not require a record to have one.
#[test]
#[ignore = "red on two decided records whose guards were renamed in ways that inverted the claim; \
            fixing either is an adjudication, not a citation repair — see #1881"]
fn every_ruling_to_guard_citation_resolves() {
    let records = rulings::records();
    let sites = guard_sites();
    let retired: BTreeMap<&str, &str> = RETIRED_GUARD_NAMES.iter().copied().collect();

    let mut faults: Vec<String> = Vec::new();
    let mut checked = 0usize;
    for record in &records {
        for guard in guard_names_named_in(&record.rationale) {
            checked += 1;
            if retired.contains_key(guard.as_str()) {
                continue;
            }
            let Some(defined) = sites.get(&guard) else {
                faults.push(format!(
                    "`{}` ({}) names `{guard}`, which is defined nowhere under {SCAN_ROOTS:?}",
                    record.id, record.status
                ));
                continue;
            };
            // A guard that runs but checks nothing passes every other
            // instrument there is — not pass/fail, and not duration.
            if !defined.iter().any(|site| site.asserts) {
                faults.push(format!(
                    "`{}` ({}) names `{guard}`, which asserts nothing ({})",
                    record.id,
                    record.status,
                    defined
                        .iter()
                        .map(|s| format!("{}:{}", s.path, s.line))
                        .collect::<Vec<_>>()
                        .join(", ")
                ));
                continue;
            }
            // `#[ignore]`d everywhere it is defined means it never runs, so the
            // ruling it is offered as enforcing is not enforced by it.
            if defined.iter().all(|site| site.ignored) {
                faults.push(format!(
                    "`{}` ({}) names `{guard}`, which is `#[ignore]`d ({})",
                    record.id,
                    record.status,
                    defined
                        .iter()
                        .map(|s| format!("{}:{}", s.path, s.line))
                        .collect::<Vec<_>>()
                        .join(", ")
                ));
            }
        }
    }

    assert!(
        faults.is_empty(),
        "{} guard citation(s) in {LEDGER_RELATIVE_PATH} name a guard that cannot enforce the \
         record:\n  {}\n\n\
         A decided record's guard is the only thing standing between its ruling and silent \
         rot. Fix the citation to name the guard that enforces the ruling now, restore the \
         guard, or — if the record names a former guard on purpose, because the rename is part \
         of the argument — declare it in `RETIRED_GUARD_NAMES` with the reason. What is \
         rejected is silence.",
        faults.len(),
        faults.join("\n  ")
    );

    // Non-vacuity, for the same reason `every_record_to_record_citation_resolves`
    // carries one: if `guard_names_named_in` stopped recognising citations this
    // check would pass forever by reading an empty set. Measured at 14 on
    // `origin/main` @ `d4552167`; the floor is far below that so ordinary prose
    // edits cannot trip it.
    assert!(
        checked >= 6,
        "only {checked} guard citations were found in {LEDGER_RELATIVE_PATH} — \
         `guard_names_named_in` is probably broken, and this check would pass by reading nothing"
    );
}

/// A retired guard name must not name a guard that still exists.
///
/// Without this, [`RETIRED_GUARD_NAMES`] is a way to silence the check above
/// rather than to declare a rename — the same reasoning as
/// [`no_retired_id_is_also_a_live_record`]. It runs armed although the list is
/// empty today, so the first entry added is checked on arrival.
#[test]
fn no_retired_guard_name_is_also_a_live_guard() {
    let sites = guard_sites();
    let overlapping: Vec<&str> = RETIRED_GUARD_NAMES
        .iter()
        .map(|(name, _)| *name)
        .filter(|name| sites.contains_key(*name))
        .collect();
    assert!(
        overlapping.is_empty(),
        "these guards are listed as retired but are defined in the tree: {overlapping:?} — \
         remove them from `RETIRED_GUARD_NAMES`, which is for guards the tree no longer has"
    );
}

/// The guard-name grammar recognises a citation by shape, not by resolution.
///
/// Both directions are pinned because the risk runs both ways. Too loose and
/// ordinary source functions named in a rationale are demanded to be tests; too
/// tight and the check reads an empty set and passes forever — and the second
/// failure is the one that looks like success.
#[test]
fn the_guard_name_grammar_separates_guards_from_source_functions() {
    // Guard-shaped: sentences, whether or not they exist in the tree. The
    // second is the deleted citation this check exists to report, so it MUST
    // match while resolving to nothing.
    for name in [
        "the_exit_follows_the_intron_bases_not_the_strand",
        "a_dup_flush_against_a_del_is_re_derived_into_the_span",
        "two_insertions_sharing_a_start_merge_into_one_member",
    ] {
        assert!(
            looks_like_a_guard_name(name),
            "{name} should be read as a guard name"
        );
    }

    // Not guard-shaped: source functions the ledger cites `::`-qualified in the
    // same rationales. Demanding these be asserting tests would make the check
    // unsatisfiable.
    for name in [
        "coalesce_coding_frame_separation",
        "insertion_is_duplication",
        "count_tandem_repeats",
        "hgvs_to_spdi",
        "reparent_junction_exit",
    ] {
        assert!(
            !looks_like_a_guard_name(name),
            "{name} is a source function and should not be read as a guard name"
        );
    }

    // Extraction reads the qualified form and takes the final segment, and
    // leaves the bare form alone.
    let named = guard_names_named_in(
        "Pinned by `it::cis_junction_crossing_shift::the_three_member_spelling_and_its_one_member_form_are_two_fixed_points` \
         and by `tests/it/corpus_prohibited_inputs.rs::a_bare_transcript_intronic_position_is_refused_in_strict_only`, \
         over `g.100_200del::300_400dup`, using `merge::is_tandem_duplication`, \
         see also `a_bare_unqualified_name_is_deliberately_not_read`.",
    );
    assert_eq!(
        named,
        [
            "a_bare_transcript_intronic_position_is_refused_in_strict_only",
            "the_three_member_spelling_and_its_one_member_form_are_two_fixed_points"
        ]
        .iter()
        .map(|s| s.to_string())
        .collect::<BTreeSet<_>>(),
        "extraction should read the two qualified guard citations and nothing else — \
         not the ring-junction HGVS string, not the source function, not the bare name"
    );
}

/// The assertion detector reads this repo's asserting helpers, not just macros.
///
/// Pinned because a stricter matcher reported a guard whose entire body is one
/// `assert_normalizes_preserving_in` call as assertionless — a false positive
/// on the very shape the check exists to find.
#[test]
fn the_assertion_detector_reads_helper_calls_as_assertions() {
    assert!(body_asserts("{ assert_eq!(a, b); }"));
    assert!(body_asserts(
        "{ assert_normalizes_preserving_in(DUP_RUN, a, b, ShuffleDirection::FivePrime); }"
    ));
    assert!(body_asserts("{ panic!(\"no\") }"));
    assert!(
        !body_asserts("{ let normalized = normalize(input); eprintln!(\"{normalized}\"); }"),
        "a body that only computes and prints asserts nothing — the #1858 shape"
    );
}

/// The guard scan is not vacuous, and reads `#[ignore]` through the doc
/// comments that sit between the attribute and the `fn`.
#[test]
fn the_guard_scan_reads_the_tree_it_claims_to() {
    let sites = guard_sites();
    // Measured at 1,649 distinct guard-shaped names on `origin/main` @
    // `d4552167`. The floor is an order of magnitude below that, so it cannot
    // be tripped by ordinary test churn but still fires if the walk breaks.
    assert!(
        sites.len() >= 200,
        "the guard scan found only {} guard-shaped functions under {SCAN_ROOTS:?} — \
         the walk or the grammar is broken",
        sites.len()
    );
    // Every citation-bearing property this check asserts on must be observable
    // at least once, or the corresponding arm is dead code that always passes.
    assert!(
        sites.values().flatten().any(|site| site.ignored),
        "no `#[ignore]`d guard was found — the attribute is not being read, so the \
         ignored-guard arm of the check above can never fire"
    );
    assert!(
        sites.values().flatten().any(|site| site.asserts),
        "no asserting guard was found — `body_asserts` is not seeing function bodies"
    );
}
