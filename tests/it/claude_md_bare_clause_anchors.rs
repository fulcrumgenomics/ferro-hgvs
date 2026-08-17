//! A bare `` `:N` `` shorthand in `CLAUDE.md`'s two adjudication tables must name
//! a line the row's own ledger record cites.
//!
//! # The hole this closes
//!
//! Four guards already read `CLAUDE.md` against the ledger —
//! [`claude_md_clause_anchors`], [`claude_md_adjudication_tables`],
//! [`ledger_prose_clause_anchors`] and [`adjudication_prose_regrowth`]. All 31 of
//! their assertions were **green for the whole window in which seven bare anchors
//! were wrong** (#2083, corrected by #2091). That is not a bug in any of them:
//! [`claude_md_clause_anchors`]'s own doc comment names the gap in as many words,
//! under "What it does NOT catch" — it keys on the literal `general.md:` token, so
//! a shorthand written `` `:34` `` is invisible to it, and nothing enforced that
//! such a shorthand agrees with the row it sits in.
//!
//! The class is unguarded in both directions. A shorthand can be wrong the day it
//! is written, and a spec-submodule bump can invalidate every one of them at once
//! while every number on the page still resolves to a real, non-blank line.
//!
//! The shape worth mechanising is the `coding-axis-merges-…` row, where a stale
//! `` `:34` `` sat four times beside a *correct* qualified `general.md:34` and the
//! two clauses argue opposite ways — `general.md:33` is the separation rule
//! ("described individually"), `general.md:34` the codon exception ("should be
//! described as a delins"). Nothing about the prose looks wrong; both numbers have
//! to be resolved against the spec before the contradiction is visible.
//!
//! # Read the spec at the GITLINK, never at the submodule working tree
//!
//! This is the trap the class is made of, and it inverts the verdict rather than
//! blurring it. A worktree in this layout routinely pins the submodule at an older
//! commit than the superproject's gitlink — `6f85311` against `565b9734` while
//! #2091 was being reviewed — and the two numberings are each internally
//! consistent:
//!
//! | line | at `6f85311` | at `565b9734` |
//! |---|---|---|
//! | `:25` | `c` coding DNA | `g` linear genomic |
//! | `:33` | *(blank)* | separation / describe individually |
//! | `:34` | separation / describe individually | codon exception |
//! | `:55` | *(blank)* | prioritisation |
//!
//! So a guard reading `assets/hgvs-nomenclature/docs/…` off the filesystem is
//! green on a stale checkout and red on a fresh one **for byte-identical source**,
//! and it reports #2091's corrections as regressions with high confidence. Every
//! spec byte this file reads therefore comes from [`SpecCheckout`], which resolves
//! the gitlink with `git ls-tree HEAD assets/hgvs-nomenclature` and reads content
//! with `git -C … show <sha>:<path>`. There is no `fs::read_to_string` on the
//! submodule path anywhere below, and
//! [`the_spec_text_comes_from_the_gitlink_not_the_submodule_working_tree`] is what
//! keeps it that way.
//!
//! # What the comparison is, and what it deliberately is not
//!
//! For each table row: the row's id names a ledger record, and that record's cited
//! lines — its structured `clauses`, its `governing` clause, and the `file:line`
//! tokens in its own `question` and `rationale` prose — are the comparand. Every
//! bare shorthand in the row must intersect that set.
//!
//! **The comparison is on line numbers, not on `(file, line)` pairs, and that is a
//! measurement rather than a preference.** The obvious rule — a shorthand inherits
//! the file of the nearest preceding qualified anchor — was implemented and run
//! over the tables: it reports **17** failures on a `CLAUDE.md` that is correct.
//! The prose simply does not obey it. In the `self-cancelling-across-ring-junctions`
//! row a `` `:130` `` follows `complex.md:5` and means `complex.md:130`; in the
//! `coding-axis-merges-…` row a `` `:33` `` follows `DNA/delins.md:81` and means
//! `general.md:33`. Both are legible to a human and neither is mechanically
//! recoverable, so the file is used for *scoping* only and never as an expectation.
//!
//! Two consequences, stated rather than left to be discovered:
//!
//! - **A number the record cites on some other file passes.** This is the price of
//!   the paragraph above, and it is the reason the four stale `` `:34` `` anchors
//!   in the `coding-axis-merges-…` row are **not** caught here: that record cites
//!   `general.md:34` legitimately, in the opposite role, so `34` is a member of the
//!   comparand. Measured on the pre-#2091 tree this file catches **3** of the 7,
//!   one in each of the other three affected rows. A content-based discriminator
//!   was prototyped for the remaining four — match the citing sentence against word
//!   n-grams of each cited spec line at the gitlink — and **rejected on
//!   measurement**: restricted to sentences carrying exactly one anchor it judges
//!   3 of 54 anchors and produces a false failure on the
//!   `duplication-must-ranks-the-label-not-the-partition` row at every n-gram width
//!   tried. A guard that judges 6% of the population and misfires is worse than the
//!   hole it patches.
//! - **Whether the row's *argument* is right is not checked.** Only that a number
//!   it writes is a number its record cites.
//!
//! # What is skipped, and why a skip cannot hide a fault silently
//!
//! `CLAUDE.md`'s rows cite source lines in the same shorthand — `` `:1963` `` after
//! `spec_conformance_axis.rs:1982` is a line of Rust, not a clause. Those are
//! skipped, on two rules:
//!
//! 1. the nearest preceding qualified anchor names a file that is not under
//!    `docs/` in the spec checkout (2 anchors), or
//! 2. the shorthand is past the end of every spec file its context could name
//!    **and** the row cites at least one non-spec file (1 anchor).
//!
//! Rule 2's second conjunct is load-bearing: without it a stale shorthand that
//! happens to fall past the end of the file would be skipped instead of failed. In
//! a row that cites no source file at all, a past-the-end shorthand is a **fault**,
//! not a skip.
//!
//! Both counts are reported in every message, and
//! [`JUDGED_BARE_ANCHORS_FLOOR`] fails the run if the scanner stops finding
//! anchors — the failure this repo names most often, where a parse that matches
//! nothing satisfies every assertion about what it matched.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Mutex, OnceLock};

use crate::common::rulings::{self, Record};

/// The vendored spec checkout, relative to the crate root.
const SPEC_DIR: &str = "assets/hgvs-nomenclature";

/// The directory inside the spec checkout that holds citable clauses. The ledger
/// spells every structured citation `docs/…`, and restricting resolution to it is
/// what keeps the spec's own root `README.md` from capturing `CLAUDE.md`'s
/// citations of *this repository's* `README.md` ruleset.
const SPEC_DOCS_ROOT: &str = "docs/";

/// One spec file, used only by the gitlink self-test.
const SPEC_PROBE_FILE: &str = "docs/recommendations/general.md";

/// The document under test, relative to the crate root.
const CLAUDE_MD: &str = "CLAUDE.md";

/// The header row that opens the open-questions table.
const OPEN_TABLE_HEADER: &str = "| record | what is open |";

/// The header row that opens the decided-rulings table.
const DECIDED_TABLE_HEADER: &str = "| record | ruling |";

/// Floor for how many rows the two tables yield between them.
///
/// Measured at 34 on the commit that added this file — one per ledger record. The
/// floor sits well below that so ordinary ledger growth or shrinkage cannot trip
/// it, while a header rename or a table reformat that left the scanner reading
/// nothing still would.
const TABLE_ROWS_FLOOR: usize = 20;

/// Floor for how many bare shorthands are actually judged.
///
/// Measured at 54 judged, 3 skipped, on the commit that added this file. Half is
/// the floor, for the same reason [`claude_md_clause_anchors`]'s
/// `ANCHORED_CITATIONS_FLOOR` gives: a scanner that stopped recognising the
/// `` `:N` `` form would pass every assertion below by judging nothing.
const JUDGED_BARE_ANCHORS_FLOOR: usize = 27;

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Run `git` and return its stdout, or the reason it could not be run.
///
/// The ambient `GIT_*` environment is cleared for the same reason
/// `generate_spec_fixture::read_submodule_commit` clears it (#1046): git exports
/// `GIT_DIR`/`GIT_WORK_TREE` pointing at the **outer** repository when it runs a
/// hook, and those override `-C`'s repository discovery — so under the pre-push
/// hook a submodule query would silently answer about the superproject.
pub(crate) fn git(cwd: &Path, args: &[&str]) -> Result<String, String> {
    let output = Command::new("git")
        .arg("-C")
        .arg(cwd)
        .args(args)
        .env_remove("GIT_DIR")
        .env_remove("GIT_WORK_TREE")
        .env_remove("GIT_INDEX_FILE")
        .output()
        .map_err(|e| format!("running `git {}`: {e}", args.join(" ")))?;
    if !output.status.success() {
        return Err(format!(
            "`git {}` in {} failed: {}",
            args.join(" "),
            cwd.display(),
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }
    Ok(String::from_utf8_lossy(&output.stdout).to_string())
}

/// The spec checkout as the commit under test records it — never as the submodule
/// working tree happens to be checked out.
///
/// Shared with [`claude_md_clause_anchors`], which resolves its one spec file the
/// same way (#2124) — the gitlink resolver is deliberately in one place so both
/// guards read the spec at the commit under test rather than at whatever the
/// submodule working tree happens to be checked out to.
pub(crate) struct SpecCheckout {
    /// The commit the superproject's tree records for [`SPEC_DIR`].
    pub(crate) gitlink: String,
    /// Every `.md` path under [`SPEC_DOCS_ROOT`] at that commit.
    paths: Vec<String>,
    /// `path -> lines`, filled on demand. A handful of files are ever needed.
    lines: Mutex<BTreeMap<String, Vec<String>>>,
}

impl SpecCheckout {
    fn load() -> Self {
        let root = crate_root();
        let listing = git(&root, &["ls-tree", "HEAD", SPEC_DIR]).unwrap_or_else(|e| {
            panic!(
                "cannot read the {SPEC_DIR} gitlink from the commit under test: {e}\n\
                 This guard resolves the spec from the gitlink on purpose — reading the submodule \
                 working tree would make it pin-dependent, which is the defect it exists for."
            )
        });
        let fields: Vec<&str> = listing.split_whitespace().collect();
        assert_eq!(
            fields.first().copied(),
            Some("160000"),
            "`git ls-tree HEAD {SPEC_DIR}` did not report a gitlink: {listing:?}"
        );
        let gitlink = fields
            .get(2)
            .unwrap_or_else(|| panic!("no object id in `git ls-tree HEAD {SPEC_DIR}`: {listing:?}"))
            .to_string();

        let spec_dir = root.join(SPEC_DIR);
        let names =
            git(&spec_dir, &["ls-tree", "-r", "--name-only", &gitlink]).unwrap_or_else(|e| {
                panic!(
                    "cannot list {SPEC_DIR} at gitlink {gitlink}: {e}\n\
                     The spec submodule is probably not initialised. Run\n    \
                     git -c protocol.file.allow=always submodule update --init {SPEC_DIR}"
                )
            });
        let paths: Vec<String> = names
            .lines()
            .filter(|p| p.starts_with(SPEC_DOCS_ROOT) && p.ends_with(".md"))
            .map(str::to_string)
            .collect();
        assert!(
            paths.len() > 20,
            "only {} `.md` files under {SPEC_DOCS_ROOT} at gitlink {gitlink} — the checkout is not \
             what this expects, and every resolution below would silently fall through",
            paths.len()
        );

        Self {
            gitlink,
            paths,
            lines: Mutex::new(BTreeMap::new()),
        }
    }

    pub(crate) fn get() -> &'static Self {
        static SPEC: OnceLock<SpecCheckout> = OnceLock::new();
        SPEC.get_or_init(SpecCheckout::load)
    }

    /// The lines of one spec file, read at the gitlink.
    pub(crate) fn lines_of(&self, path: &str) -> Vec<String> {
        if let Some(cached) = self.lines.lock().expect("spec cache").get(path) {
            return cached.clone();
        }
        let spec_dir = crate_root().join(SPEC_DIR);
        let text = git(&spec_dir, &["show", &format!("{}:{path}", self.gitlink)])
            .unwrap_or_else(|e| panic!("reading {path} at gitlink {}: {e}", self.gitlink));
        let lines: Vec<String> = text.lines().map(str::to_string).collect();
        self.lines
            .lock()
            .expect("spec cache")
            .insert(path.to_string(), lines.clone());
        lines
    }

    fn line_count(&self, path: &str) -> usize {
        self.lines_of(path).len()
    }

    /// The spec files a citation's file token could name.
    ///
    /// Empty means "not a spec file" — a `.rs` source path, or this repository's
    /// own `README.md`. Several means the token is a bare basename that exists in
    /// more than one molecule directory (`delins.md` is under `DNA/`, `RNA/` and
    /// `protein/`), which `CLAUDE.md` writes 22 times. That ambiguity is a
    /// jurisdiction question this guard does not adjudicate; it keeps every
    /// candidate and lets the line-number comparison do the work.
    fn candidates(&self, token: &str) -> Vec<String> {
        if !token.ends_with(".md") {
            return Vec::new();
        }
        let suffix = format!("/{token}");
        self.paths
            .iter()
            .filter(|p| p.as_str() == token || p.ends_with(&suffix))
            .cloned()
            .collect()
    }
}

/// A `file:N`, `file:N-M` or bare `` `:N` `` / `` `:N-M` `` token, as written.
#[derive(Debug, PartialEq, Eq)]
enum Token {
    /// A citation naming its own file.
    Qualified {
        file: String,
        first: usize,
        last: usize,
    },
    /// A shorthand whose file is left to the surrounding prose.
    Bare {
        first: usize,
        last: usize,
        spelling: String,
    },
}

/// Whether a byte may appear in a citation's file token.
fn is_path_byte(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'.' | b'/' | b'-')
}

/// Whether a run of path bytes reads as a filename rather than as prose.
///
/// Requires an extension of one to five ASCII letters, which is what keeps an
/// HGVS description (`c.850_901delinsTTCC…`) and a version string out.
fn looks_like_a_file(token: &str) -> bool {
    match token.rsplit_once('.') {
        Some((stem, extension)) => {
            !stem.is_empty()
                && (1..=5).contains(&extension.len())
                && extension.chars().all(|c| c.is_ascii_alphabetic())
        }
        None => false,
    }
}

/// Every citation token in `text`, in document order.
///
/// Hand-rolled, matching [`claude_md_clause_anchors`]'s `citations_in`, so this
/// guard's own measured constants cannot be perturbed by a change to a shared
/// scanner — and so [`the_scanner_reads_both_citation_forms`] can pin it against
/// literals. A bare shorthand is recognised only inside backticks (`` `:33` ``),
/// which is how `CLAUDE.md` writes them and what keeps an ordinary `1:30` out.
fn tokens_in(text: &str) -> Vec<Token> {
    let bytes = text.as_bytes();
    let mut out = Vec::new();
    let mut index = 0usize;
    while index < bytes.len() {
        if bytes[index] != b':' {
            index += 1;
            continue;
        }
        let digits_start = index + 1;
        let mut cursor = digits_start;
        while cursor < bytes.len() && bytes[cursor].is_ascii_digit() {
            cursor += 1;
        }
        if cursor == digits_start {
            index += 1;
            continue;
        }
        let first: usize = text[digits_start..cursor].parse().expect("ascii digits");
        let mut last = first;
        let mut end = cursor;
        if bytes.get(cursor) == Some(&b'-') {
            let tail_start = cursor + 1;
            let mut tail = tail_start;
            while tail < bytes.len() && bytes[tail].is_ascii_digit() {
                tail += 1;
            }
            if tail > tail_start {
                last = text[tail_start..tail].parse().expect("ascii digits");
                end = tail;
            }
        }

        let mut prefix_start = index;
        while prefix_start > 0 && is_path_byte(bytes[prefix_start - 1]) {
            prefix_start -= 1;
        }
        let prefix = &text[prefix_start..index];
        if prefix.is_empty() {
            let opened = index > 0 && bytes[index - 1] == b'`';
            let closed = bytes.get(end) == Some(&b'`');
            if opened && closed {
                out.push(Token::Bare {
                    first,
                    last,
                    spelling: text[index..end].to_string(),
                });
            }
        } else if looks_like_a_file(prefix) {
            out.push(Token::Qualified {
                file: prefix.to_string(),
                first,
                last,
            });
        }
        index = end;
    }
    out
}

/// One row of an adjudication table.
struct Row {
    id: String,
    text: String,
}

/// The rows of the table opened by `header`, in document order.
///
/// Reads from the header to the first line that is not a table row, so prose after
/// a table cannot silently extend it — the same boundary rule
/// [`claude_md_adjudication_tables`] uses. That guard is deliberately not reused:
/// its helpers are private, and #2100/#2107 are open against its file and against
/// the ledger, so a shared edit would collide with work in review.
fn rows_in_table(markdown: &str, header: &str) -> Vec<Row> {
    let mut lines = markdown.lines().skip_while(|line| line.trim() != header);
    assert!(
        lines.next().is_some(),
        "CLAUDE.md has no table opened by {header:?}, so this guard would judge nothing. If the \
         table was renamed, rename it here too rather than deleting this assertion"
    );

    let mut rows = Vec::new();
    for line in lines {
        let line = line.trim();
        if !line.starts_with('|') {
            break;
        }
        if line.starts_with("|---") {
            continue;
        }
        let first_cell = line
            .trim_start_matches('|')
            .split('|')
            .next()
            .unwrap_or_default();
        if let Some(id) = first_cell.split('`').nth(1) {
            rows.push(Row {
                id: id.to_string(),
                text: line.to_string(),
            });
        }
    }
    rows
}

/// Both tables' rows, with the non-vacuity floor applied once, here, so no test
/// below can run against an empty parse.
fn table_rows(markdown: &str) -> Vec<Row> {
    let mut rows = rows_in_table(markdown, OPEN_TABLE_HEADER);
    rows.extend(rows_in_table(markdown, DECIDED_TABLE_HEADER));
    assert!(
        rows.len() >= TABLE_ROWS_FLOOR,
        "only {} adjudication-table rows were parsed out of {CLAUDE_MD}, below the floor of \
         {TABLE_ROWS_FLOOR} — the row scanner is broken, and every assertion below would pass by \
         reading almost nothing",
        rows.len()
    );
    rows
}

fn claude_md() -> String {
    let path = crate_root().join(CLAUDE_MD);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

/// Every `(spec file, line)` the record cites, from its structured clauses, its
/// governing clause and its own prose.
///
/// Prose is included because the ledger genuinely cites clauses there that the
/// `clauses` array does not carry — `absolute-prohibition-enforcement-stage`
/// discusses `checklist.md:45` in its rationale and lists four other checklist
/// lines as clauses — and `CLAUDE.md`'s row summarises the whole record, not only
/// its structured half. Restricting the comparand to `clauses` was measured and
/// produces exactly that one false failure.
fn cited_lines(record: &Record, spec: &SpecCheckout) -> BTreeSet<(String, usize)> {
    let mut text = String::new();
    text.push_str(&record.question);
    text.push(' ');
    text.push_str(&record.rationale);
    text.push(' ');
    if let Some(governing) = &record.governing {
        text.push_str(governing);
        text.push(' ');
    }
    for citation in &record.citations {
        text.push_str(&citation.clause);
        text.push(' ');
        text.push_str(&citation.quote);
        text.push(' ');
    }

    let mut cited = BTreeSet::new();
    for token in tokens_in(&text) {
        let Token::Qualified { file, first, last } = token else {
            continue;
        };
        for path in spec.candidates(&file) {
            let count = spec.line_count(&path);
            for line in first..=last {
                if line >= 1 && line <= count {
                    cited.insert((path.clone(), line));
                }
            }
        }
    }
    cited
}

/// What the prose has most recently established about which file a bare shorthand
/// belongs to.
enum FileContext {
    /// No qualified citation yet in this row.
    Unset,
    /// The candidate spec files of the last qualified citation. Never empty.
    Spec(Vec<String>),
    /// The last qualified citation named a file that is not a spec clause file —
    /// a `.rs` source path, or this repository's own `README.md`. A shorthand
    /// after one of those is a source line, not a clause.
    Source,
}

/// The tally a run produces, so a message can say what it judged rather than only
/// what it found wrong.
#[derive(Default)]
struct Census {
    judged: usize,
    skipped_source_context: usize,
    skipped_past_end_of_file: usize,
}

/// Judge one row, appending to `faults` and `census`.
fn judge_row(
    row: &Row,
    record: &Record,
    spec: &SpecCheckout,
    faults: &mut Vec<String>,
    census: &mut Census,
) {
    let cited = cited_lines(record, spec);
    let cited_numbers: BTreeSet<usize> = cited.iter().map(|(_, line)| *line).collect();
    let tokens = tokens_in(&row.text);
    let row_cites_a_non_spec_file = tokens.iter().any(|token| match token {
        Token::Qualified { file, .. } => spec.candidates(file).is_empty(),
        Token::Bare { .. } => false,
    });

    let mut context = FileContext::Unset;
    for token in &tokens {
        match token {
            Token::Qualified { file, .. } => {
                let candidates = spec.candidates(file);
                context = if candidates.is_empty() {
                    FileContext::Source
                } else {
                    FileContext::Spec(candidates)
                };
            }
            Token::Bare {
                first,
                last,
                spelling,
            } => {
                if matches!(context, FileContext::Source) {
                    census.skipped_source_context += 1;
                    continue;
                }
                if let FileContext::Spec(candidates) = &context {
                    let longest = candidates
                        .iter()
                        .map(|path| spec.line_count(path))
                        .max()
                        .unwrap_or(0);
                    if *last > longest {
                        if row_cites_a_non_spec_file {
                            census.skipped_past_end_of_file += 1;
                            continue;
                        }
                        faults.push(format!(
                            "`{}` in the `{}` row names line {last}, past the end of every file \
                             its context could mean ({}, longest {longest} lines), and the row \
                             cites no source file that could explain it",
                            spelling,
                            row.id,
                            candidates.join(", ")
                        ));
                        continue;
                    }
                }
                census.judged += 1;
                if (*first..=*last).any(|line| cited_numbers.contains(&line)) {
                    continue;
                }
                let mut cited_here: Vec<String> = cited
                    .iter()
                    .map(|(path, line)| {
                        let name = path.rsplit('/').next().unwrap_or(path.as_str());
                        format!("{name}:{line}")
                    })
                    .collect();
                cited_here.sort();
                cited_here.dedup();
                faults.push(format!(
                    "`{}` in the `{}` row names no line the record cites. It cites: {}",
                    spelling,
                    row.id,
                    cited_here.join(", ")
                ));
            }
        }
    }
}

/// The guard.
///
/// See the module docs for what a bare shorthand is compared against, why the
/// comparison is on line numbers rather than on `(file, line)` pairs, and what it
/// therefore does not catch.
#[test]
fn every_bare_shorthand_names_a_line_its_record_cites() {
    let spec = SpecCheckout::get();
    let markdown = claude_md();
    let rows = table_rows(&markdown);
    let records: BTreeMap<String, Record> = rulings::records()
        .into_iter()
        .map(|record| (record.id.clone(), record))
        .collect();

    let mut faults = Vec::new();
    let mut census = Census::default();
    for row in &rows {
        let Some(record) = records.get(&row.id) else {
            // Not a skip. A row naming no record has no comparand, and treating
            // that as "nothing to check" is the failure mode this whole file is
            // about. `claude_md_adjudication_tables` reports the set difference
            // properly; here it is enough that the run cannot go green on it.
            faults.push(format!(
                "the `{}` row names no ledger record, so its shorthands have nothing to be \
                 checked against",
                row.id
            ));
            continue;
        };
        judge_row(row, record, spec, &mut faults, &mut census);
    }

    assert!(
        faults.is_empty(),
        "{} bare `:N` shorthand(s) in {CLAUDE_MD}'s adjudication tables disagree with the ledger \
         record their row summarises (spec gitlink {}, {} judged, {} skipped as source-line \
         references, {} skipped as past-end-of-file):\n  {}\n\n\
         A shorthand is stale when the spec submodule moves under it, and every number still \
         resolves to a real, non-blank line afterwards — so re-point it at the clause the record \
         cites. Resolve the number against the spec at the GITLINK, not against the submodule \
         working tree: this layout routinely pins it older, and under the older numbering the \
         stale anchors read as correct.",
        faults.len(),
        &spec.gitlink[..8.min(spec.gitlink.len())],
        census.judged,
        census.skipped_source_context,
        census.skipped_past_end_of_file,
        faults.join("\n  "),
    );

    assert!(
        census.judged >= JUDGED_BARE_ANCHORS_FLOOR,
        "only {} bare shorthands were judged, below the floor of {JUDGED_BARE_ANCHORS_FLOOR} \
         ({} skipped as source-line references, {} as past-end-of-file) — the scanner, the row \
         parser or the spec resolution stopped matching, and the assertion above would pass by \
         judging nothing",
        census.judged,
        census.skipped_source_context,
        census.skipped_past_end_of_file,
    );
}

/// The spec bytes must come from the gitlink, checked against a second derivation.
///
/// The equality below is true by construction today; the point is that it stops
/// being true the moment somebody replaces [`SpecCheckout::lines_of`] with a
/// filesystem read, which is the one edit that would silently make this whole file
/// pin-dependent. The two sides are derived differently on purpose —
/// `show <rev>:<path>` against a tree lookup plus `cat-file blob` — because two
/// readings of the same command are one reading.
///
/// It also *reports* the pin, so a run in a worktree whose submodule is checked out
/// elsewhere says so in its output instead of leaving the reader to wonder.
#[test]
fn the_spec_text_comes_from_the_gitlink_not_the_submodule_working_tree() {
    let spec = SpecCheckout::get();
    let spec_dir = crate_root().join(SPEC_DIR);

    let listing = git(
        &spec_dir,
        &["ls-tree", &spec.gitlink, "--", SPEC_PROBE_FILE],
    )
    .expect("list the probe file at the gitlink");
    let blob = listing
        .split_whitespace()
        .nth(2)
        .unwrap_or_else(|| panic!("no blob id for {SPEC_PROBE_FILE} in {listing:?}"))
        .to_string();
    let independent = git(&spec_dir, &["cat-file", "blob", &blob]).expect("read the blob");
    let independent_lines: Vec<String> = independent.lines().map(str::to_string).collect();

    let through_the_guard = spec.lines_of(SPEC_PROBE_FILE);
    assert!(
        through_the_guard.len() > 100,
        "{SPEC_PROBE_FILE} has only {} lines at gitlink {} — the checkout is not what this expects",
        through_the_guard.len(),
        spec.gitlink
    );
    assert_eq!(
        through_the_guard, independent_lines,
        "the text this guard reads for {SPEC_PROBE_FILE} is not the blob the commit's gitlink \
         records. Whatever it is reading instead — most likely the submodule working tree — makes \
         every line number below a property of the local checkout rather than of the commit."
    );

    // Non-vacuous only where the two pins differ, which is the interesting case
    // and the one that cannot be arranged in CI. Report it either way.
    let head = git(&spec_dir, &["rev-parse", "HEAD"])
        .map(|s| s.trim().to_string())
        .unwrap_or_default();
    let on_disk = std::fs::read_to_string(spec_dir.join(SPEC_PROBE_FILE)).ok();
    if head == spec.gitlink {
        eprintln!(
            "submodule working tree is AT the gitlink ({}); the pin-independence half of this \
             test is vacuous here and is exercised by checking out another revision",
            &spec.gitlink[..8]
        );
    } else if let Some(on_disk) = on_disk {
        let on_disk_lines: Vec<String> = on_disk.lines().map(str::to_string).collect();
        eprintln!(
            "submodule working tree is at {} while the gitlink is {}; working-tree copy {} the \
             gitlink's",
            &head[..8.min(head.len())],
            &spec.gitlink[..8],
            if on_disk_lines == through_the_guard {
                "happens to match"
            } else {
                "DIFFERS from"
            }
        );
        assert_eq!(
            through_the_guard, independent_lines,
            "with the working tree at a different revision, the guard must still read the gitlink"
        );
    }
}

/// The two scanners, pinned against literals.
///
/// Both are hand-rolled and both fail in the direction that looks like success: a
/// citation scanner that finds nothing, and a row parser that returns no rows,
/// each leave the guard green while judging nothing.
#[test]
fn the_scanner_reads_both_citation_forms() {
    let found = tokens_in(
        "`general.md:33` governs, `:34` is the exception, `DNA/delins.md:44-47` is the passage, \
         `:79-84` follows it, `merge.rs:2000` is code and 1:30 is a time.",
    );
    assert_eq!(
        found,
        vec![
            Token::Qualified {
                file: "general.md".to_string(),
                first: 33,
                last: 33
            },
            Token::Bare {
                first: 34,
                last: 34,
                spelling: ":34".to_string()
            },
            Token::Qualified {
                file: "DNA/delins.md".to_string(),
                first: 44,
                last: 47
            },
            Token::Bare {
                first: 79,
                last: 84,
                spelling: ":79-84".to_string()
            },
            Token::Qualified {
                file: "merge.rs".to_string(),
                first: 2000,
                last: 2000
            },
        ],
        "the scanner must read both forms and their ranges, and must not read an unbackticked \
         `1:30` as a shorthand"
    );

    // An HGVS description carries a `.` and a `:` and must not be read as a file.
    assert!(
        tokens_in("`LRG_199t1:c.850_901delinsTTCCTCGATGCCTG` is one variant").is_empty(),
        "an accession-qualified description is not a citation"
    );

    let rows = rows_in_table(
        &format!(
            "prose\n{OPEN_TABLE_HEADER}\n|---|---|\n| `alpha-beta-gamma` | cites `:1` |\n\
             | `delta-epsilon-zeta` | cites `:2` |\nprose after\n| `not-a-row-any-more` | x |\n"
        ),
        OPEN_TABLE_HEADER,
    );
    let ids: Vec<&str> = rows.iter().map(|row| row.id.as_str()).collect();
    assert_eq!(
        ids,
        ["alpha-beta-gamma", "delta-epsilon-zeta"],
        "the row parser must read every row and must stop at the first non-table line"
    );
}

/// The sabotage, as a committed test rather than as a claim.
///
/// The real tables are correct, so neither of the tests above can show that this
/// guard is *able* to fail. This one hands [`judge_row`] a row it controls,
/// against a real ledger record, and pins both directions: the number the record
/// cites passes, and a neighbouring number it does not cite fails by name.
///
/// The record chosen is the one whose two adjacent clauses argue opposite ways, so
/// the fixture is the defect's own shape rather than an arbitrary pair.
#[test]
fn a_shorthand_the_record_does_not_cite_is_reported_by_row_id() {
    let spec = SpecCheckout::get();
    let records = rulings::records();
    let subject = "contiguous-insertion-split-by-a-blocked-derivation";
    let record = records
        .iter()
        .find(|record| record.id == subject)
        .unwrap_or_else(|| panic!("the ledger no longer has a record named {subject}"));

    let good = Row {
        id: subject.to_string(),
        text: format!(
            "| `{subject}` | `general.md:33` is stated over two variants, and `:33` \
                       does not reach a one-variant locus |"
        ),
    };
    let mut faults = Vec::new();
    let mut census = Census::default();
    judge_row(&good, record, spec, &mut faults, &mut census);
    assert!(
        faults.is_empty(),
        "the negative control must pass, or the failure below proves nothing: {faults:?}"
    );
    assert_eq!(
        census.judged, 1,
        "the control must actually judge its shorthand"
    );

    let stale = Row {
        id: subject.to_string(),
        text: format!(
            "| `{subject}` | `general.md:33` is stated over two variants, and `:34` \
                       does not reach a one-variant locus |"
        ),
    };
    let mut faults = Vec::new();
    let mut census = Census::default();
    judge_row(&stale, record, spec, &mut faults, &mut census);
    assert_eq!(
        census.judged, 1,
        "the stale shorthand must be judged, not skipped"
    );
    assert_eq!(
        faults.len(),
        1,
        "expected exactly one fault, got {faults:?}"
    );
    assert!(
        faults[0].contains(subject) && faults[0].contains(":34"),
        "the message must name the row and the offending number: {}",
        faults[0]
    );
}

/// A shorthand after a source-file citation is a source line, and a shorthand past
/// the end of the spec file is only excused by a row that cites source at all.
///
/// Both are skips in the real tables, and a skip that cannot be told from a pass is
/// the failure this file's module docs open on — so both rules are pinned here,
/// including the conjunct that makes the second one safe.
#[test]
fn a_source_line_shorthand_is_skipped_and_an_unexplained_one_is_not() {
    let spec = SpecCheckout::get();
    let records = rulings::records();
    let subject = "contiguous-insertion-split-by-a-blocked-derivation";
    let record = records
        .iter()
        .find(|record| record.id == subject)
        .unwrap_or_else(|| panic!("the ledger no longer has a record named {subject}"));

    let source_context = Row {
        id: subject.to_string(),
        text: format!(
            "| `{subject}` | the guard is at `spec_conformance_axis.rs:1982`; `:1963` is its \
             message |"
        ),
    };
    let mut faults = Vec::new();
    let mut census = Census::default();
    judge_row(&source_context, record, spec, &mut faults, &mut census);
    assert!(
        faults.is_empty(),
        "a source-line shorthand is not a clause: {faults:?}"
    );
    assert_eq!(
        (census.judged, census.skipped_source_context),
        (0, 1),
        "it must be skipped as a source-line reference, and counted as one"
    );

    let unexplained = Row {
        id: subject.to_string(),
        text: format!("| `{subject}` | `general.md:33` governs, and so does `:9999` |"),
    };
    let mut faults = Vec::new();
    let mut census = Census::default();
    judge_row(&unexplained, record, spec, &mut faults, &mut census);
    assert_eq!(
        census.skipped_past_end_of_file, 0,
        "a row citing no source file may not have a past-the-end shorthand excused"
    );
    assert_eq!(
        faults.len(),
        1,
        "expected exactly one fault, got {faults:?}"
    );
    assert!(
        faults[0].contains("past the end"),
        "the message must say why: {}",
        faults[0]
    );
}
