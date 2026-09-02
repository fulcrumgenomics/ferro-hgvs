//! A `path.rs::name` citation in a ruling record's **prose** must resolve to an
//! item that actually exists in the file it names.
//!
//! # The gap this closes
//!
//! A ruling record narrates its reasoning in prose, and that prose pins the
//! ruling to the tests that enforce it — "Pinned by
//! `tests/it/foo.rs::a_guard_name`". Two registers can then disagree *inside one
//! record*: the structured `guard.tests` array (checked by
//! `ruling_guard_field.rs`) names the live test, while the sentence a human reads
//! names one that no longer exists. That is the more misleading direction — a
//! reader who follows the prose name concludes the ruling is unguarded when it is
//! in fact guarded.
//!
//! Measured on the ledger when this file was written: two decided records carried
//! exactly that split. `alignment-only-symbol-in-a-description`'s prose named
//! `…_refused_for_dash_and_accepted_for_x` while its `guard.tests` — and the tree
//! — carried `…_refused_in_every_mode_for_both_x_and_dash`; and
//! `delins-merge-vs-individual-gap-two-or-more`'s prose named
//! `…_converge_on_the_recommended_form` while the test that pins it is
//! `…_are_two_fixed_points` (whose body asserts the convergence after #1835). See
//! #2016.
//!
//! # Why this form, and this form only
//!
//! A prose guard reference is spelled several ways — `path.rs::name`, bare
//! `module::name`, and `it::module::name` among them. This check reads **only**
//! the `path.rs::name` form, for two reasons that reinforce each other:
//!
//! - **It is the actionable form.** The `.rs` path names *where* to look, so the
//!   citation can be resolved precisely — the leaf must be defined in **that
//!   file**, not merely exist somewhere in the tree. A leaf that resolves in a
//!   different file is the failure this form makes visible and the looser forms
//!   cannot.
//! - **It keeps the check free of judgement calls.** A `path.rs::name` token
//!   cannot be confused with ordinary prose about the code, so the extractor
//!   needs no heuristic to decide whether a backticked token is a citation.
//!
//! The parser is `common::rulings::split_guard_citation`, the same grammar
//! `ruling_guard_field.rs` resolves `guard.tests` with, so the two never drift on
//! what a `path.rs::name` citation is.
//!
//! # It composes with the two existing guards; none subsumes another
//!
//! - **`ruling_guard_field.rs`** reads the structured `guard.tests` field, not
//!   prose. A record can carry a live `guard.tests` and still describe itself
//!   wrongly in the sentence beside it — which is exactly the #2016 shape.
//! - **`ruling_citation_currency.rs`**'s `every_ruling_to_guard_citation_resolves`
//!   reads the bare `module::name` prose forms this file leaves alone, resolving
//!   each leaf against the whole tree rather than against the named file. It is
//!   `#[ignore]`d pending #1881, because one bare-`module::name` citation
//!   (`derivation-may-not-be-bounded-by-the-inputs-spelling`'s) names a guard
//!   whose rename **inverted the claim** — repairing that is an adjudication, not
//!   a citation fix. This file does not reach that citation (it is a bare
//!   `module::` form, not a `path.rs::` one), so activating this check does not
//!   depend on that adjudication.
//!
//! # A resolved leaf need not be a `#[test]`
//!
//! A `path.rs::name` citation may legitimately name a **module** rather than a
//! test — `issue_390_hgvs_vcf_provider_coverage.rs::range_payload_never_reaches_alt`
//! is a `mod`, and reporting it as a missing test would be a false positive. So
//! the check resolves the leaf against every *item* a file declares — `fn`,
//! `mod`, `const`, `static`, `struct`, `enum`, `trait`, `type`, `macro_rules!` —
//! not the test set alone. Whether the item is the *right* enforcement for the
//! ruling is a judgement this check does not make; that it exists where the prose
//! says it does is what it buys.

use std::collections::BTreeSet;
use std::path::PathBuf;

use crate::common::rulings::{self, split_guard_citation, LEDGER_RELATIVE_PATH};

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Floor for how many `path.rs::name` prose citations the ledger must carry.
///
/// Measured at 7 when this file was written. The floor sits below that so
/// ordinary ledger edits cannot trip it, while a reader that stopped seeing
/// citations — the failure that would make the check pass by resolving nothing —
/// still would.
const PROSE_CITATION_FLOOR: usize = 4;

/// The keywords that introduce a named Rust item, each as it appears at the start
/// of a (trimmed) declaration line.
///
/// `macro_rules!` is spelled with its bang because that is how the name follows
/// it (`macro_rules! foo`); the rest are followed by the name after a space.
const ITEM_KEYWORDS: &[&str] = &[
    "fn ",
    "mod ",
    "const ",
    "static ",
    "struct ",
    "enum ",
    "trait ",
    "type ",
    "union ",
    "macro_rules!",
];

/// The backticked runs in `prose`.
///
/// Odd-indexed chunks of a split on the backtick are the quoted runs; this is the
/// same extraction `ruling_guard_field.rs` uses for exemption reasons.
fn backticked_tokens(prose: &str) -> Vec<String> {
    prose
        .split('`')
        .enumerate()
        .filter(|(index, _)| index % 2 == 1)
        .map(|(_, chunk)| chunk.to_string())
        .collect()
}

/// Every `path.rs::name` citation in a record's prose, as `(path, name)`.
///
/// The prose scanned is the `question` and `rationale` fields — the record's two
/// free-prose fields. `guard.tests` is deliberately not read here: that is the
/// structured field `ruling_guard_field.rs` owns, and a `path.rs::name` in it is
/// not prose.
fn prose_citations(record: &rulings::Record) -> Vec<(String, String)> {
    let mut found = Vec::new();
    for prose in [&record.question, &record.rationale] {
        for token in backticked_tokens(prose) {
            if let Some((path, name)) = split_guard_citation(&token) {
                found.push((path.to_string(), name.to_string()));
            }
        }
    }
    found
}

/// Every item name declared in `text`.
///
/// A coarse line scanner rather than a parser: the question is only "does this
/// file declare an item of this name", and a syntax-tree dependency for that
/// would be a large amount of machinery pointed at a small question. It reads
/// items at any indentation, so a `mod tests { fn … }` unit test and a
/// module-scoped `const` are both found.
///
/// It does **not** mask raw-string literals, so an item *name* that appears only
/// inside a fixture string in the named file would resolve. The consequence is
/// bounded and one-directional: a stale citation could slip through only by
/// colliding with an item name that exists nowhere in the file except inside a
/// string literal — a far narrower hazard than the tree-wide phantom problem
/// `ruling_guard_field.rs` masks against, because resolution here is pinned to
/// the one file the citation names.
fn declared_items(text: &str) -> BTreeSet<String> {
    let mut names = BTreeSet::new();
    for raw in text.lines() {
        let mut line = raw.trim_start();
        // Step over visibility and qualifier prefixes so `pub fn`,
        // `pub(crate) fn`, `pub(super) mod`, `const fn`, `unsafe fn`,
        // `async fn`, `pub async fn` and the like all reach the keyword.
        for prefix in [
            "pub(crate) ",
            "pub(super) ",
            "pub(self) ",
            "pub ",
            "const ",
            "unsafe ",
            "async ",
            "default ",
        ] {
            if let Some(rest) = line.strip_prefix(prefix) {
                // `const NAME` is itself an item; only step over `const` when a
                // function keyword follows it (`const fn`).
                if prefix == "const " && !rest.trim_start().starts_with("fn ") {
                    break;
                }
                line = rest.trim_start();
            }
        }
        for keyword in ITEM_KEYWORDS {
            if let Some(rest) = line.strip_prefix(keyword) {
                let name: String = rest
                    .trim_start()
                    .chars()
                    .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
                    .collect();
                if !name.is_empty() {
                    names.insert(name);
                }
                break;
            }
        }
    }
    names
}

/// The item names declared in the file at `scan_relative`, or `None` if the file
/// does not exist.
fn items_in_file(scan_relative: &str) -> Option<BTreeSet<String>> {
    let path = crate_root().join(scan_relative);
    if !path.is_file() {
        return None;
    }
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    Some(declared_items(&text))
}

// ---------------------------------------------------------------------------
// The check
// ---------------------------------------------------------------------------

/// Every `path.rs::name` prose citation names a file that exists and an item
/// defined in it.
///
/// This is the check the file exists for. A guard renamed out from under a
/// record's *prose* leaves the sentence a reader follows pointing at nothing,
/// while the structured `guard.tests` — and every other instrument — stays green.
#[test]
fn every_prose_pin_citation_resolves() {
    let records = rulings::records();
    assert!(
        records.len() >= 10,
        "only {} records read from {LEDGER_RELATIVE_PATH} — this check would pass by reading \
         nothing",
        records.len()
    );

    let mut checked = 0usize;
    let mut unresolved: Vec<String> = Vec::new();
    for record in &records {
        for (path, name) in prose_citations(record) {
            checked += 1;
            match items_in_file(&path) {
                None => unresolved.push(format!(
                    "`{}` names `{path}::{name}`, but {path} does not exist under the crate root",
                    record.id
                )),
                Some(items) if !items.contains(&name) => unresolved.push(format!(
                    "`{}` names `{path}::{name}`, but {path} declares no item `{name}`",
                    record.id
                )),
                Some(_) => {}
            }
        }
    }

    // Non-vacuity: if the extractor stopped recognising citations, every branch
    // above would pass by iterating an empty set. The floor is well below the
    // measured count so ordinary ledger edits cannot trip it.
    assert!(
        checked >= PROSE_CITATION_FLOOR,
        "only {checked} `path.rs::name` prose citations found in {LEDGER_RELATIVE_PATH}, below the \
         floor of {PROSE_CITATION_FLOOR} — `prose_citations` is probably broken, and this check \
         would pass by reading nothing"
    );

    assert!(
        unresolved.is_empty(),
        "these `path.rs::name` citations in the prose of {LEDGER_RELATIVE_PATH} name an item that \
         does not exist where the prose says it does. Repair the citation to name the item that \
         pins the ruling now — or, if the guard was renamed in a way that changed what it claims, \
         that is an adjudication and belongs in the record, not a citation swap:\n  {}",
        unresolved.join("\n  ")
    );
}

// ---------------------------------------------------------------------------
// The scanner's own tests
// ---------------------------------------------------------------------------

/// Extraction reads the `path.rs::name` form and nothing else — not a bare
/// `module::name`, not an `it::module::name`, not an HGVS string, not prose.
#[test]
fn extraction_reads_only_the_file_anchored_form() {
    let record = record_with_rationale(
        "Pinned by `tests/it/foo.rs::a_guard_name` and by \
         `defect_371_transcript_exit::a_bare_module_form`, over `g.100_200del::300_400dup`, \
         see also `it::cis_junction_crossing_shift::a_qualified_but_pathless_form`.",
    );
    let cited: BTreeSet<(String, String)> = prose_citations(&record).into_iter().collect();
    assert_eq!(
        cited,
        [("tests/it/foo.rs".to_string(), "a_guard_name".to_string())]
            .into_iter()
            .collect::<BTreeSet<_>>(),
        "only the path-anchored form is a citation this check reads"
    );
}

/// The `question` field is scanned too, not only `rationale`.
#[test]
fn extraction_reads_both_prose_fields() {
    let mut record = record_with_rationale("Pinned by `tests/it/bar.rs::from_the_rationale`.");
    record.question = "See `tests/it/baz.rs::from_the_question`.".to_string();
    let names: BTreeSet<String> = prose_citations(&record)
        .into_iter()
        .map(|(_, name)| name)
        .collect();
    assert!(names.contains("from_the_rationale"));
    assert!(
        names.contains("from_the_question"),
        "a citation in the question field must be read too"
    );
}

/// The item scanner finds a test `fn`, a `mod`, and a module-scoped `const`, and
/// does not invent a name from a line that declares none.
///
/// The `mod` case is the one that matters: a `path.rs::name` citation may name a
/// module (the ledger has one), and resolving it as an item is what stops this
/// check reporting a real module as a missing test.
#[test]
fn the_item_scanner_finds_functions_modules_and_consts() {
    let source = "\
#[test]
fn a_guard_function() {
    assert!(true);
}

mod range_payload_never_reaches_alt {
    fn helper() {}
}

const A_PINNED_COUNT: usize = 3;
pub(crate) fn a_visible_helper() {}

    let not_an_item = 1;
";
    let items = declared_items(source);
    assert!(
        items.contains("a_guard_function"),
        "a #[test] fn must be found"
    );
    assert!(
        items.contains("range_payload_never_reaches_alt"),
        "a module declaration must be found, or a `path.rs::module` citation reads as missing"
    );
    assert!(items.contains("A_PINNED_COUNT"), "a const must be found");
    assert!(
        items.contains("a_visible_helper"),
        "a pub(crate) fn must be found"
    );
    assert!(items.contains("helper"), "an indented fn must be found");
    assert!(
        !items.contains("not_an_item"),
        "a `let` binding is not an item declaration"
    );
}

/// A missing file and a present-file-but-missing-item are both faults, and a real
/// citation is not. Exercises the resolver against the live tree so it cannot
/// pass by reading nothing.
#[test]
fn resolution_flags_a_dangling_citation_and_accepts_a_live_one() {
    // A real citation from the ledger — the file exists and declares the item.
    let live = record_with_rationale(
        "Pinned by \
         `tests/it/corpus_prohibited_inputs.rs::an_alignment_only_symbol_is_refused_in_every_mode_for_both_x_and_dash`.",
    );
    assert!(
        prose_citations_all_resolve(&live),
        "a citation naming a real item in a real file must resolve"
    );

    // A file that does not exist.
    let missing_file =
        record_with_rationale("Pinned by `tests/it/there_is_no_such_file.rs::whatever`.");
    assert!(
        !prose_citations_all_resolve(&missing_file),
        "a citation naming a missing file must not resolve"
    );

    // A real file, an item it does not declare.
    let missing_item = record_with_rationale(
        "Pinned by `tests/it/corpus_prohibited_inputs.rs::this_test_does_not_exist_anywhere`.",
    );
    assert!(
        !prose_citations_all_resolve(&missing_item),
        "a citation naming an absent item in a real file must not resolve"
    );
}

/// A `path.rs::module` citation resolves to the module, and is not reported as a
/// missing test — the false positive the issue warns about.
#[test]
fn a_module_citation_resolves_as_an_item() {
    let record = record_with_rationale(
        "Pinned by `tests/it/issue_390_hgvs_vcf_provider_coverage.rs::range_payload_never_reaches_alt`.",
    );
    assert!(
        prose_citations_all_resolve(&record),
        "a `path.rs::module` citation must resolve as an item, not read as a missing test"
    );
}

// ---------------------------------------------------------------------------
// Test helpers
// ---------------------------------------------------------------------------

/// Whether every `path.rs::name` prose citation in `record` resolves.
fn prose_citations_all_resolve(record: &rulings::Record) -> bool {
    prose_citations(record)
        .into_iter()
        .all(|(path, name)| items_in_file(&path).is_some_and(|items| items.contains(&name)))
}

/// A minimal record carrying only the prose fields this check reads, so a
/// scanner assertion is about behaviour rather than about a live ledger row.
fn record_with_rationale(rationale: &str) -> rulings::Record {
    rulings::Record {
        id: "a-synthetic-record".to_string(),
        status: "decided".to_string(),
        question: String::new(),
        rationale: rationale.to_string(),
        summary: None,
        applies_to: Vec::new(),
        equivalence_classes: Vec::new(),
        guard: rulings::Guard::Declined("synthetic".to_string()),
        house_choice: None,
        citations: Vec::new(),
        governing: None,
    }
}

/// The scan root a citation's path is resolved against exists.
///
/// A guard on the guard: if the crate root stopped containing `tests/it`, every
/// citation would report its file as missing and the check would be loud rather
/// than silent — but this pins the assumption directly.
#[test]
fn the_scan_root_exists() {
    assert!(
        crate_root().join("tests/it").is_dir(),
        "tests/it not found under the crate root — every citation would read as a missing file"
    );
}
