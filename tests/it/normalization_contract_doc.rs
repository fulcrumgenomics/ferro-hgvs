//! Generates and gates `docs/NORMALIZATION_CONTRACT.md` — the published
//! rendering of the adjudication ledger (#1552).
//!
//! # What the document is for
//!
//! Ferro's decisions about what to emit where the HGVS recommendations are
//! silent, ambiguous or self-contradictory live in the `rulings` section of
//! `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`. That is the
//! right home — the records are read by the spec-fixture generator and pinned by
//! `ruling_records_are_intact`, so a record that no test reads cannot rot.
//!
//! But the audience for those decisions is not only ferro's own test suite. A
//! contributor deciding a new case, and a downstream consumer trying to predict
//! what a string will normalize to, both need to scan every ruling at once. The
//! document is that index: one row per record — status, governing clause, and
//! the record's own one-sentence summary — over the full reasoning in the ledger.
//!
//! # Why it is GENERATED rather than written
//!
//! #1552 makes it a standing constraint that the document be *generated from or
//! checked against* the records rather than hand-maintained beside them: two
//! hand-maintained copies of a normative rule is the failure mode this project
//! has already hit, with the census constants and with the counts in
//! `clause_ruling_index.rs`'s own header. Generation is the stronger of the two
//! options offered, so it is the one taken — each record's status, governing
//! clause and one-sentence summary below is the record's own text, copied
//! mechanically. Nothing in the rendered body is paraphrase, and there is no
//! editorial pass in which a paraphrase could be introduced. The full reasoning,
//! the clause quotes and any scope are not rendered here; they live in the ledger.
//!
//! The hand-authored part is the preamble, and it is authored **here**, in this
//! module, so that it too has exactly one copy.
//!
//! # Why a test rather than an `[[example]]` with `--check`
//!
//! The ledger reader is `common::rulings`, which lives in this integration-test
//! crate; an example target cannot use it, so an example generator would need a
//! second parser for the same file. That is the drift this document exists to
//! avoid, one level down. Living here also means the gate needs no new CI wiring
//! — the required `Test` context is a rollup over the job that runs this binary,
//! so a stale document is already a merge blocker.
//!
//! Regenerate with:
//!
//! ```text
//! BLESS_CONTRACT_DOC=1 cargo nextest run --features dev --test it \
//!   -E 'test(normalization_contract_doc)'
//! ```
//!
//! # What this module deliberately does not do
//!
//! It states none of the normalization rules, and
//! [`the_preamble_does_not_restate_the_ruleset`] enforces that against the
//! ruleset page itself rather than against a copy of it. Single-sourcing that
//! ruleset is part of a ruling in the ledger this document renders, so a
//! document that restated it would contradict its own contents.
//!
//! It also names no record id, for the reason `common/rulings.rs` gives: this
//! tree is scanned by `ruling_citation_currency.rs`, and an id written into a
//! renderer is a second place a status claim could go stale. Every id in the
//! output arrives from the ledger at render time.

use std::fmt::Write as _;
use std::path::PathBuf;

use super::common::rulings::{records, Record, Role};

/// Where the generated document lives, relative to the crate root.
const DOC_RELATIVE_PATH: &str = "docs/NORMALIZATION_CONTRACT.md";

/// The environment variable that turns [`the_published_document_is_current`]
/// from a comparison into a write.
///
/// Named in the failure message so the fix never has to be looked up, following
/// `BLESS_MOCK_PIN` in the projection and biocommons suites.
const BLESS_VAR: &str = "BLESS_CONTRACT_DOC";

/// The source file the shipped default partition arm is read out of.
///
/// Read rather than restated: see [`shipped_default_partition_arm`].
const PARTITION_SOURCE_RELATIVE_PATH: &str = "src/normalize/merge.rs";

/// Heading that separates the hand-authored preamble from the rendered records.
///
/// [`the_preamble_does_not_restate_the_ruleset`] scopes itself to the
/// text above this line, because the text below it is the ledger's and is not
/// this module's to police.
const RECORDS_HEADING: &str = "## The records";

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn read(relative: &str) -> String {
    let path = repo_root().join(relative);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

// --------------------------------------------------------------------------
// The shipped default partition arm.
//
// The rulings are decisions about what ferro's output *should* be, and are
// independent of which block partitioner is selected. Several of them are
// nevertheless explicit about whether their ruling is implemented on the
// shipped default or only under a candidate arm, so the document has to say
// which default it was generated against — and a hand-written sentence saying
// so would be stale on the day the default flips, which is exactly the change
// in flight.
//
// So it is read out of the code. `PartitionRule` is private to the crate, so
// this is a scan of one match arm rather than a call; the scan is narrow enough
// that a rename fails it loudly instead of answering wrongly.
// --------------------------------------------------------------------------

/// The `FERRO_PARTITION` value that an unset environment selects, as spelled in
/// `PARTITION_RULE_NAMES`.
///
/// Located by finding the arm of `partition_rule_from_env` that matches `None`
/// and reading the `PartitionRule` variant it returns. Panics if the shape is
/// not found, which is the correct outcome: a silent fallback here would let the
/// document publish `live` after the default had moved, and the whole point of
/// deriving it is that it cannot.
fn shipped_default_partition_arm() -> String {
    let source = read(PARTITION_SOURCE_RELATIVE_PATH);
    let body = source
        .split_once("fn partition_rule_from_env")
        .unwrap_or_else(|| {
            panic!(
                "no `fn partition_rule_from_env` in {PARTITION_SOURCE_RELATIVE_PATH} — the \
                 function was renamed or moved, and the default arm this document publishes can \
                 no longer be derived"
            )
        })
        .1;

    let arm = body
        .lines()
        .take_while(|line| !line.starts_with('}'))
        .find(|line| line.contains("None"))
        .unwrap_or_else(|| {
            panic!(
                "no match arm in `partition_rule_from_env` handles `None` — the unset case moved, \
                 so the shipped default can no longer be read from \
                 {PARTITION_SOURCE_RELATIVE_PATH}"
            )
        });

    // The arm may name the variant directly (`None => Ok(PartitionRule::Live)`)
    // or, since #1835, name a CONSTANT that holds it
    // (`None | Some("") => Ok(DEFAULT_PARTITION_RULE)`). Both are resolved
    // here, and the indirection is followed rather than special-cased away:
    // the point of deriving this from source is that the document cannot
    // publish a default the code no longer has, and a scraper that only
    // understood the literal form would fail closed on a refactor that changed
    // nothing about the value. It still fails closed on anything it cannot
    // resolve — the `expect` below has no fallback.
    let variant = match arm.split_once("PartitionRule::") {
        Some((_, rest)) => rest
            .chars()
            .take_while(char::is_ascii_alphanumeric)
            .collect::<String>(),
        None => {
            let constant = arm
                .rsplit_once("Ok(")
                .map(|(_, rest)| {
                    rest.chars()
                        .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
                        .collect::<String>()
                })
                .filter(|name| !name.is_empty())
                .unwrap_or_else(|| {
                    panic!(
                        "the `None` arm of `partition_rule_from_env` names neither a \
                         `PartitionRule::` variant nor a constant this can resolve: {arm:?}"
                    )
                });
            let decl = format!("const {constant}: PartitionRule = PartitionRule::");
            let rest = source
                .split_once(&decl)
                .map(|(_, rest)| rest)
                .unwrap_or_else(|| {
                    panic!(
                    "the `None` arm returns `{constant}`, but no `{decl}…` declaration exists in \
                     {PARTITION_SOURCE_RELATIVE_PATH} — the shipped default cannot be derived"
                )
                });
            rest.chars()
                .take_while(char::is_ascii_alphanumeric)
                .collect::<String>()
        }
    };
    assert!(
        !variant.is_empty(),
        "could not read a variant name out of the `None` arm: {arm:?}"
    );

    let name = kebab_case(&variant);
    let names = declared_partition_rule_names(&source);
    assert!(
        names.contains(&name),
        "the unset arm of `partition_rule_from_env` returns `PartitionRule::{variant}`, whose \
         kebab-case name {name:?} is not one of the declared `PARTITION_RULE_NAMES` {names:?} — \
         either the naming convention changed or the arm is not offered by the diagnostic"
    );
    name
}

/// `CanonicalCoalesced` -> `canonical-coalesced`.
fn kebab_case(variant: &str) -> String {
    let mut out = String::new();
    for (index, ch) in variant.chars().enumerate() {
        if ch.is_ascii_uppercase() && index > 0 {
            out.push('-');
        }
        out.extend(ch.to_lowercase());
    }
    out
}

/// The arm names `PARTITION_RULE_NAMES` declares, read from the same source.
fn declared_partition_rule_names(source: &str) -> Vec<String> {
    let literal = source
        .split_once("const PARTITION_RULE_NAMES")
        .unwrap_or_else(|| panic!("no `PARTITION_RULE_NAMES` in {PARTITION_SOURCE_RELATIVE_PATH}"))
        .1;
    // `= [` rather than `[`: the declaration's type annotation (`[&str; 4]`)
    // comes first, and splitting on the bare bracket reads that instead — which
    // fails as "the arm is not offered by the diagnostic", a message pointing at
    // the wrong file entirely.
    let list = literal
        .split_once("= [")
        .and_then(|(_, rest)| rest.split_once(']'))
        .expect("`PARTITION_RULE_NAMES` is initialised from an array literal")
        .0;
    list.split(',')
        .filter_map(|entry| {
            let entry = entry.trim().trim_matches('"');
            (!entry.is_empty()).then(|| entry.to_string())
        })
        .collect()
}

// --------------------------------------------------------------------------
// Rendering.
// --------------------------------------------------------------------------

/// The whole document.
fn render() -> String {
    let records = records();
    let decided: Vec<&Record> = sorted(&records, "decided");
    let open: Vec<&Record> = sorted(&records, "undecided");
    assert_eq!(
        decided.len() + open.len(),
        records.len(),
        "a record carries a status that is neither `decided` nor `undecided`, so it would be \
         rendered under no heading and silently vanish from this document"
    );

    let mut out = String::new();
    out.push_str(&preamble(
        records.len(),
        decided.len(),
        open.len(),
        &shipped_default_partition_arm(),
    ));
    let _ = writeln!(out, "{RECORDS_HEADING}\n");
    let _ = writeln!(
        out,
        "Every ruling record on one screen: its status, the clause that governs it where one \
         applies (a decided record names its governing clause or, as a house choice, cites none; \
         an open record names no authority), a one-sentence statement of the ruling, and — where \
         ferro's shipped output lags the decision — a note on what is implemented. Read the \
         ruling, not the id \
         — an id states the record's *question*, and at least one of the questions below is \
         answered in the negative. An `open` record names a conflict ferro has **not** settled; \
         whatever ferro does with that case today is the status quo, not a ruling, and must not \
         be cited as one. The full reasoning, the clauses each record quotes, and any scope on \
         the ruling live in the ledger itself, \
         [`hgvs_spec_normalization_overrides.json`](../tests/fixtures/grammar/hgvs_spec_normalization_overrides.json); \
         this document is an index into it, not a copy of it.\n"
    );
    out.push_str(&render_index(&open, &decided));

    // Exactly one trailing newline, and no trailing whitespace on any line.
    // Not cosmetic: `end-of-file-fixer` and `trailing-whitespace` are wired as
    // pre-commit hooks, so an output that ends `\n\n` is rewritten the moment it
    // is committed — and the next run of this test compares the rewritten file
    // against the un-rewritten render and fails, on a document nobody touched.
    // Pinned by `the_render_survives_the_file_hygiene_hooks`.
    let mut out: String = out
        .lines()
        .map(str::trim_end)
        .collect::<Vec<_>>()
        .join("\n");
    out.truncate(out.trim_end().len());
    out.push('\n');
    out
}

/// Records with `status`, ordered by id so the document is stable under
/// reordering of the ledger.
fn sorted<'a>(records: &'a [Record], status: &str) -> Vec<&'a Record> {
    let mut selected: Vec<&Record> = records.iter().filter(|r| r.status == status).collect();
    selected.sort_by(|a, b| a.id.cmp(&b.id));
    selected
}

/// The hand-authored half, and the only prose in the file this module owns.
fn preamble(total: usize, decided: usize, open: usize, default_arm: &str) -> String {
    let mut out = String::new();
    let _ = write!(
        out,
        r#"<!--
GENERATED FILE — do not edit by hand.

Rendered from the `rulings` section of
tests/fixtures/grammar/hgvs_spec_normalization_overrides.json
by tests/it/normalization_contract_doc.rs. Edit the ledger, then regenerate:

    BLESS_CONTRACT_DOC=1 cargo nextest run --features dev --test it \
      -E 'test(normalization_contract_doc)'

An edit made here instead is reverted by the next regeneration, and fails CI
before that.
-->

# Ferro's normalization contract

**{total} adjudication records — {decided} decided, {open} open.**

## What this document is

The HGVS recommendations are, in places, silent, ambiguous, or self-contradictory. A
normalizer still has to emit one string. Where ferro has had to decide such a question, the
decision is recorded as a **ruling record**, and this document is an index of every one of
those records: its status, the clause that governs it where one applies (a decided record names
its governing clause or, as a house choice, cites none; an open record names no authority), the
record's own one-sentence statement of the ruling, and — where ferro's shipped output lags the
decision — a note on what is implemented.

It is generated, not written. Each ruling summary is the record's own `summary` field, copied
verbatim from
[`tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`](../tests/fixtures/grammar/hgvs_spec_normalization_overrides.json);
nothing here is paraphrased by the renderer. The full reasoning behind each ruling, the spec
clauses it quotes, the text each clause was quoted against, and any scope on the ruling are
**not** reproduced here — they live in the ledger, which is what the build enforces.

## What this document is not

**It is not the ruleset.** What ferro's output is allowed to be — which properties are
absolute and which are best effort, what happens where the spec determines no answer, and what
must be disclosed when a choice changes — is stated once, in
[the normalization rules](src/reference/normalization-rules.md). That statement is
deliberately not reproduced here, and single-sourcing it is itself part of one of the rulings
below. Where this document and the ruleset page appear to disagree, the ruleset page governs and this
document has a bug.

**It is not a substitute for the records.** It is a reading of them. The records are what the
build enforces.

**It is not a spec.** The HGVS recommendations are upstream — rendered at
[hgvs-nomenclature.org](https://hgvs-nomenclature.org/) and sourced from the
`HGVSnomenclature/hgvs-nomenclature` repository that `assets/hgvs-nomenclature` vendors. Each
governing clause below is spelled as `path:line` and links to that exact line in the **pinned**
commit of that repository, so the citation resolves to the spec version ferro actually pins
rather than to whatever the site currently shows.

## How to read the index

- **The id states the QUESTION, not the ruling, and the two can be opposites.** Read the
  ruling column, not the id. One record below is titled as the position it *rejects*.
- **`undecided` is a first-class state**, not an oversight. An open record states a conflict
  and declines to settle it; whatever ferro does with that case today is the status quo, and
  citing the behaviour as a decision is the error the record exists to prevent.
- **A one-line ruling states the decision, not how to carry it out.** Many records add an
  explicit scope — one axis, one direction, one shape — and the mechanism or tie-break that
  turns the decision into an output; this index shows neither. Read it to scan every ruling,
  and open the ledger record before *acting* on one. A summary that reads like a complete rule
  may still be narrower, or more conditional, than it looks.
- **The Implementation column names a KNOWN GAP between the ruling and shipped output**, in
  the record's own words, where one exists — a decision implemented for only one shape, one
  axis, or only under a candidate partition arm. A dash means no such gap is recorded, which
  is not a guarantee the ruling ships in full: the ruling is the decision, and its shipped
  status is the ledger record's to state. The note is a verbatim quote of the record's
  reasoning, so it cannot drift from it.

## Which build this describes

The rulings are decisions about what ferro's output *should* be, and they do not depend on
which block partitioner is selected at run time. Several records are nevertheless explicit
about whether their ruling is already live in shipped output or is implemented only under a
candidate arm, so the answer depends on the default — and that default is currently in motion.

**As generated, the shipped default — what `FERRO_PARTITION` unset selects — is
`{default_arm}`.** That sentence is not written by hand: the generator reads the arm out of
`src/normalize/merge.rs`, so a change to the default fails this document's own test until it
is regenerated, and cannot leave a stale claim behind. See
[Comparing normalization rules](src/guide/comparing-rules.md#comparing-normalization-rules-ferro_partition)
for the knob and its traps.

"#
    );
    out
}

/// The index table: one row per record, open records first, then decided,
/// each already sorted by id. The columns are the record's id, its status, the
/// clause that governs it, a one-sentence statement of the ruling, and — where
/// shipped output lags the decision — a note on what is implemented.
fn render_index(open: &[&Record], decided: &[&Record]) -> String {
    let commit = pinned_spec_commit();
    let mut out = String::new();
    let _ = writeln!(
        out,
        "| Record | Status | Governing clause | Ruling | Implementation |"
    );
    let _ = writeln!(out, "|---|---|---|---|---|");
    for record in open.iter().chain(decided.iter()) {
        let _ = writeln!(
            out,
            "| `{}` | {} | {} | {} | {} |",
            record.id,
            record.status,
            governing_cell(record, &commit),
            ruling_cell(record),
            implementation_cell(record),
        );
    }
    let _ = writeln!(out);
    out
}

/// The implementation cell. Present only where a decided ruling's shipped output
/// lags its decision, and then rendered as the record's own note — a verbatim
/// substring of its rationale (see [`Record::implemented_note`]), so nothing here
/// is paraphrase. An em dash where the record notes no gap; that is the absence
/// of a *noted* gap, not a guarantee the ruling ships in full.
fn implementation_cell(record: &Record) -> String {
    match record.implemented_note.as_deref() {
        Some(note) => cell(note),
        None => "—".to_string(),
    }
}

/// The governing-clause cell. A decided record either names a governing clause
/// or is a house choice that cites none; an open record names no authority. The
/// deviated-from and also-cited clauses are the ledger's to carry, not this
/// index's. Each named clause links to its exact line in the pinned spec
/// checkout (see [`clause_link`]).
fn governing_cell(record: &Record, commit: &str) -> String {
    if record.house_choice.is_some() {
        return "house choice (cites none)".to_string();
    }
    let clauses: Vec<String> = record
        .citations
        .iter()
        .filter(|c| c.role == Role::Governing)
        .map(|c| clause_link(&c.clause, commit))
        .collect();
    if clauses.is_empty() {
        "—".to_string()
    } else {
        clauses.join(", ")
    }
}

/// The upstream spec repository's blob tree, so a `path:line` clause citation
/// can link to the exact source line it names.
const SPEC_REPO_BLOB_BASE: &str = "https://github.com/HGVSnomenclature/hgvs-nomenclature/blob";

/// The pinned `assets/hgvs-nomenclature` commit, read from the superproject's
/// recorded gitlink so a clause link points at exactly the spec version ferro
/// pins — not at a moving branch, where the line numbers would drift.
///
/// Derived, never hardcoded: a submodule bump changes it and fails
/// [`the_published_document_is_current`] until the document is regenerated, so a
/// link can never survive pointing at a version the tree no longer pins (the
/// same discipline as the derived default-partition arm). Reading the gitlink
/// out of `HEAD`'s tree means the submodule need not even be checked out.
///
/// This assumes a git checkout: it shells `git` and reads `HEAD`, so it panics
/// rather than skips in an environment without a `.git` (a source tarball or
/// vendored tree). That is acceptable because this whole suite is a `dev`-only
/// generator gate that already resolves everything relative to the source tree,
/// and CI runs it inside the checkout.
fn pinned_spec_commit() -> String {
    let output = std::process::Command::new("git")
        .arg("-C")
        .arg(repo_root())
        .args(["rev-parse", "HEAD:assets/hgvs-nomenclature"])
        .output()
        .expect("run `git rev-parse HEAD:assets/hgvs-nomenclature`");
    assert!(
        output.status.success(),
        "`git rev-parse HEAD:assets/hgvs-nomenclature` failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let sha = String::from_utf8(output.stdout)
        .expect("git prints utf-8")
        .trim()
        .to_string();
    assert!(
        sha.len() == 40 && sha.bytes().all(|b| b.is_ascii_hexdigit()),
        "expected a 40-char submodule commit sha, got {sha:?}"
    );
    sha
}

/// A clause citation rendered as a link to its source line, or as plain
/// backticked text if it is not in `path:line` shape. The link label stays the
/// `path:line` string a reader recognises; only the target is added.
fn clause_link(clause: &str, commit: &str) -> String {
    match spec_url(clause, commit) {
        Some(url) => format!("[`{clause}`]({url})"),
        None => format!("`{clause}`"),
    }
}

/// A GitHub source URL for a `path:line` (or `path:start-end`) clause citation
/// at the pinned spec `commit`. `None` if the clause is not in that shape —
/// including when the line component is not numeric, so a malformed citation
/// falls back to plain text rather than minting a broken `#L…` anchor.
///
/// `?plain=1` forces GitHub's source view: a Markdown file otherwise renders
/// rich, where a `#L…` anchor resolves to nothing. The plain view shows line
/// numbers and highlights the cited span.
fn spec_url(clause: &str, commit: &str) -> Option<String> {
    let (path, lines) = clause.rsplit_once(':')?;
    if path.is_empty() {
        return None;
    }
    let numeric = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let anchor = match lines.split_once('-') {
        Some((start, end)) if numeric(start) && numeric(end) => format!("#L{start}-L{end}"),
        None if numeric(lines) => format!("#L{lines}"),
        _ => return None,
    };
    Some(format!(
        "{SPEC_REPO_BLOB_BASE}/{commit}/{path}?plain=1{anchor}"
    ))
}

/// The ruling cell. A decided record renders its one-sentence `summary`
/// verbatim; an open record has no ruling to state, so its question stands in,
/// marked as unsettled. Both are collapsed to one line and their pipes escaped
/// so a Markdown table row is not broken.
fn ruling_cell(record: &Record) -> String {
    match record.summary.as_deref() {
        Some(summary) => cell(summary),
        None => format!("*Undecided.* {}", cell(record.question.trim())),
    }
}

/// One table cell: whitespace collapsed to a single line, and `|` escaped so it
/// does not terminate the cell.
fn cell(text: &str) -> String {
    collapse(text).replace('|', "\\|")
}

/// A quote or paragraph on one line, so a Markdown blockquote does not break
/// across the record's own line wrapping.
fn collapse(text: &str) -> String {
    text.split_whitespace().collect::<Vec<_>>().join(" ")
}

// --------------------------------------------------------------------------
// Tests.
// --------------------------------------------------------------------------

/// The committed document matches what the ledger renders to.
///
/// With [`BLESS_VAR`] set this writes the document instead of comparing, which
/// is how it is regenerated; the failure message names the command, so the fix
/// never has to be looked up.
#[test]
fn the_published_document_is_current() {
    let rendered = render();
    let path = repo_root().join(DOC_RELATIVE_PATH);

    if std::env::var(BLESS_VAR).is_ok() {
        std::fs::write(&path, &rendered)
            .unwrap_or_else(|e| panic!("write {}: {e}", path.display()));
        return;
    }

    let committed = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!(
            "read {}: {e}\n\nThe published contract document is missing. Generate it with:\n  \
             {BLESS_VAR}=1 cargo nextest run --features dev --test it \
             -E 'test(normalization_contract_doc)'",
            path.display()
        )
    });

    if committed != rendered {
        let (committed_lines, rendered_lines) =
            (committed.lines().count(), rendered.lines().count());
        let first_difference = committed
            .lines()
            .zip(rendered.lines())
            .position(|(a, b)| a != b)
            .map(|index| index + 1);
        panic!(
            "{DOC_RELATIVE_PATH} is stale against the ledger \
             ({committed_lines} committed lines vs {rendered_lines} rendered; first differing \
             line: {first_difference:?}).\n\nRegenerate it — do not edit it by hand:\n  \
             {BLESS_VAR}=1 cargo nextest run --features dev --test it \
             -E 'test(normalization_contract_doc)'"
        );
    }
}

/// Every record reaches the index as a row that states its status.
///
/// The table has no structural grouping, so a record must appear as its own row
/// carrying its status or it silently vanishes from this document. `render`
/// asserts the status partition; this asserts the rendered result.
#[test]
fn every_record_is_published_with_its_status() {
    let rendered = render();
    let records = records();
    assert!(
        !records.is_empty(),
        "the ledger holds no records, so this document would be vacuous"
    );

    for record in &records {
        let row_prefix = format!("| `{}` | {} |", record.id, record.status);
        assert!(
            rendered.contains(&row_prefix),
            "record {} is not published as a `{}` row in {DOC_RELATIVE_PATH}",
            record.id,
            record.status
        );
    }
}

/// The ruling reaches the reader, not just the id: a decided record shows its
/// summary, an open record its question, and every record its governing clause.
///
/// A document that published ids alone would look complete and answer nothing —
/// the state #1552 was filed about. The full quotes and reasoning are the
/// ledger's; this index publishes the answer and the authority.
#[test]
fn every_record_publishes_its_ruling_and_governing_clause() {
    let rendered = render();
    for record in records() {
        match record.summary.as_deref() {
            Some(summary) => assert!(
                rendered.contains(&cell(summary)),
                "decided record {} is published without its ruling summary",
                record.id
            ),
            None => assert!(
                rendered.contains(&format!("*Undecided.* {}", cell(record.question.trim()))),
                "open record {} is published without its question",
                record.id
            ),
        }
        for citation in record
            .citations
            .iter()
            .filter(|c| c.role == Role::Governing)
        {
            assert!(
                rendered.contains(&format!("`{}`", citation.clause)),
                "record {} names governing clause {} without publishing it",
                record.id,
                citation.clause
            );
        }
    }
}

/// The index carries each decided record's full row — id, status, governing
/// clause, and its one-sentence summary — so a reader can scan every ruling on
/// one screen. An id states the question, not the answer, and the previous
/// place this scan existed was a hand-maintained agent-guidance table that
/// drifted three times.
#[test]
fn the_index_carries_every_decided_summary() {
    let rendered = render();
    for record in records().iter().filter(|r| r.status == "decided") {
        let summary = record
            .summary
            .as_deref()
            .unwrap_or_else(|| panic!("decided record {} has no summary", record.id));
        let row = format!(
            "| `{}` | decided | {} | {} | {} |",
            record.id,
            governing_cell(record, &pinned_spec_commit()),
            cell(summary),
            implementation_cell(record),
        );
        assert!(
            rendered.contains(&row),
            "index is missing the row for decided record {}",
            record.id
        );
    }
}

/// Every record carrying an `implemented_note` publishes it, and the column it
/// renders into exists. The note is the record's own text (a verbatim substring
/// of its rationale), so a reader scanning the table sees the shipped-status gap
/// without opening the ledger.
///
/// Non-vacuous by assertion: if no record carries a note the column proves
/// nothing, and the guard says so rather than passing over an empty set — the
/// same discipline as `the_index_is_not_vacuous`.
#[test]
fn the_index_publishes_every_implementation_note() {
    let rendered = render();
    assert!(
        rendered.contains("| Record | Status | Governing clause | Ruling | Implementation |"),
        "the index has no Implementation column header"
    );
    let mut noted = 0usize;
    for record in records() {
        if let Some(note) = record.implemented_note.as_deref() {
            noted += 1;
            assert!(
                rendered.contains(&cell(note)),
                "record {} carries an `implemented_note` that is not published in the index",
                record.id
            );
        }
    }
    assert!(
        noted > 0,
        "no record carries an `implemented_note`, so this guard and the Implementation column \
         prove nothing — either the field has fallen out of use (drop the column) or the ledger \
         regressed"
    );
}

/// The preamble points at the ruleset page and does not restate it.
///
/// Single-sourcing that ruleset is part of a ruling this very document renders,
/// so a preamble that restated it would contradict its own contents. The rule
/// openers are read from the page itself — a copy of them kept here would be
/// a second copy of the thing being guarded.
///
/// Scoped to the preamble deliberately. Below [`RECORDS_HEADING`] the text is
/// the ledger's own, and a record is entitled to quote the ruleset it rules
/// under; policing that would be policing the ledger from a renderer.
#[test]
fn the_preamble_does_not_restate_the_ruleset() {
    let rendered = render();
    let preamble = rendered
        .split_once(RECORDS_HEADING)
        .expect("the document has a records section")
        .0;

    assert!(
        preamble.contains("src/reference/normalization-rules.md"),
        "the preamble must send the reader to the ruleset page rather than restating it"
    );

    let restated = crate::common::ruleset_page::restated_rule_openers(preamble);
    assert!(
        restated.is_empty(),
        "the preamble restates the ruleset page's rules {restated:?}; link to it instead"
    );
}

/// The published default partition arm is the one the code actually selects.
///
/// This is what stops the document going stale across the default flip: the arm
/// is derived, so flipping it fails [`the_published_document_is_current`] until
/// the document is regenerated.
#[test]
fn the_published_default_partition_arm_is_derived_from_the_code() {
    let arm = shipped_default_partition_arm();
    let rendered = render();
    assert!(
        rendered.contains(&format!(
            "**As generated, the shipped default — what `FERRO_PARTITION` unset selects — is\n`{arm}`.**"
        )),
        "the document does not publish the derived default arm `{arm}`"
    );

    let source = read(PARTITION_SOURCE_RELATIVE_PATH);
    let names = declared_partition_rule_names(&source);
    assert_eq!(
        names.len(),
        4,
        "`PARTITION_RULE_NAMES` declares {names:?}; this document's preamble describes a knob \
         with the four arms that set has always had, so a change in arity is a prose change too"
    );
}

#[test]
fn kebab_case_converts_every_declared_arm_name() {
    assert_eq!(kebab_case("Live"), "live");
    assert_eq!(kebab_case("Shadow"), "shadow");
    assert_eq!(kebab_case("Canonical"), "canonical");
    assert_eq!(kebab_case("CanonicalCoalesced"), "canonical-coalesced");
}

#[test]
fn spec_url_links_to_the_pinned_source_line() {
    assert_eq!(
        spec_url("docs/recommendations/DNA/delins.md:81", "0123abc").as_deref(),
        Some(
            "https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/0123abc/\
             docs/recommendations/DNA/delins.md?plain=1#L81"
        )
    );
    assert_eq!(
        spec_url("docs/recommendations/DNA/delins.md:44-47", "0123abc").as_deref(),
        Some(
            "https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/0123abc/\
             docs/recommendations/DNA/delins.md?plain=1#L44-L47"
        )
    );
    assert_eq!(spec_url("not-a-citation", "0123abc"), None);
    assert_eq!(spec_url("trailing-colon:", "0123abc"), None);
    // A non-numeric line component is not a citation: fall back to plain text
    // rather than mint a broken `#L…` anchor.
    assert_eq!(spec_url("docs/foo.md:bar", "0123abc"), None);
    assert_eq!(spec_url("docs/foo.md:12-x", "0123abc"), None);
}

/// Every governing clause in the document links to its exact line in the pinned
/// spec commit, so a reader clicks through to the authority rather than
/// resolving a `path:line` by hand. Also asserts the check is non-vacuous.
#[test]
fn governing_clauses_link_to_the_pinned_spec() {
    let rendered = render();
    let commit = pinned_spec_commit();
    let mut checked = 0;
    for record in records() {
        for citation in record
            .citations
            .iter()
            .filter(|c| c.role == Role::Governing)
        {
            let url = spec_url(&citation.clause, &commit).unwrap_or_else(|| {
                panic!(
                    "governing clause {} is not a `path:line` citation",
                    citation.clause
                )
            });
            assert!(
                rendered.contains(&format!("[`{}`]({url})", citation.clause)),
                "governing clause {} is not published as a link to the pinned spec",
                citation.clause
            );
            checked += 1;
        }
    }
    assert!(
        checked > 0,
        "no governing clause was checked — the document or the role filter is wrong"
    );
}

/// The rendered document is already in the shape the file-hygiene hooks want.
///
/// `trailing-whitespace` and `end-of-file-fixer` run as pre-commit hooks over
/// every file in the tree, including generated ones. If the render disagreed
/// with them, committing the document would rewrite it and
/// [`the_published_document_is_current`] would then fail against a file nobody
/// had edited — a failure that reads as ledger drift and is not.
#[test]
fn the_render_survives_the_file_hygiene_hooks() {
    let rendered = render();
    assert!(
        rendered.ends_with('\n') && !rendered.ends_with("\n\n"),
        "the render must end with exactly one newline, or `end-of-file-fixer` rewrites it"
    );
    for (number, line) in rendered.lines().enumerate() {
        assert_eq!(
            line.trim_end(),
            line,
            "line {} of the render carries trailing whitespace, which `trailing-whitespace` \
             strips: {line:?}",
            number + 1
        );
    }
}
