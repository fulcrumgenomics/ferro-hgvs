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
//! what a string will normalize to, both need to read the ladder end to end, and
//! a JSON `rationale` blob is a poor medium for that. The document is that
//! reading.
//!
//! # Why it is GENERATED rather than written
//!
//! #1552 makes it a standing constraint that the document be *generated from or
//! checked against* the records rather than hand-maintained beside them: two
//! hand-maintained copies of a normative rule is the failure mode this project
//! has already hit, with the census constants and with the counts in
//! `clause_ruling_index.rs`'s own header. Generation is the stronger of the two
//! options offered, so it is the one taken — every question, verdict, clause and
//! quote below is the record's own text, copied mechanically. Nothing in the
//! rendered body is paraphrase, and there is no editorial pass in which a
//! paraphrase could be introduced.
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
//! It states none of `README.md`'s normalization rules, and
//! [`the_preamble_does_not_restate_the_readme_ruleset`] enforces that against
//! `README.md` itself rather than against a copy of it. Single-sourcing that
//! ruleset is part of a ruling in the ledger this document renders, so a
//! document that restated it would contradict its own contents.
//!
//! It also names no record id, for the reason `common/rulings.rs` gives: this
//! tree is scanned by `ruling_citation_currency.rs`, and an id written into a
//! renderer is a second place a status claim could go stale. Every id in the
//! output arrives from the ledger at render time.

use std::collections::BTreeMap;
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
/// [`the_preamble_does_not_restate_the_readme_ruleset`] scopes itself to the
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
    out.push_str(&contents(&decided, &open));
    let _ = writeln!(out, "{RECORDS_HEADING}\n");
    let _ = writeln!(
        out,
        "### Open questions\n\nThese are recorded conflicts that ferro has **not** settled. \
         Whatever ferro currently does with them is the status quo, not a ruling, and must not \
         be cited as one.\n"
    );
    for record in &open {
        out.push_str(&render_record(record));
    }
    let _ = writeln!(
        out,
        "### Decided\n\nEach heading is the record's id. Read the ruling rather than the id — an \
         id states the record's *question*, and at least one of the questions below is answered \
         in the negative.\n"
    );
    for record in &decided {
        out.push_str(&render_record(record));
    }

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
decision is recorded as a **ruling record**, and this document is a rendering of every one of
those records: its question, its verdict, the spec clauses it names, the text each clause was
quoted against, and the reasoning in full.

It is generated, not written. Every question, verdict, clause and quote below is the record's
own text, reproduced mechanically from
[`tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`](../tests/fixtures/grammar/hgvs_spec_normalization_overrides.json).
The one-line summary beside each decided record in the contents list is the record's own
`summary` field, copied verbatim; nothing here is paraphrased by the renderer.

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

**It is not a spec.** The HGVS recommendations are upstream, at
[hgvs-nomenclature.org](https://hgvs-nomenclature.org/); the clause citations below are into
the pinned `assets/hgvs-nomenclature` checkout, spelled as `path:line`.

## How to read a record

- **The id states the QUESTION, not the ruling, and the two can be opposites.** Read the
  ruling. One record below is titled as the position it *rejects*.
- **`undecided` is a first-class state**, not an oversight. An open record states a conflict
  and declines to settle it; whatever ferro does with that case today is the status quo, and
  citing the behaviour as a decision is the error the record exists to prevent.
- **A decided record is often narrower than it sounds.** Several carry an explicit scope — to
  one axis, to one direction, to one shape — inside the ruling text. The scope is part of the
  ruling.
- **Counter-evidence is recorded inside the ruling, not omitted from it.** Where the spec
  argues the other way, the record says so and says why it was outweighed.

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

/// The two id lists, so a reader arriving with an id can find it.
fn contents(decided: &[&Record], open: &[&Record]) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "## Contents\n");
    let _ = writeln!(out, "**Open questions**\n");
    for record in open {
        let _ = writeln!(out, "- [`{}`](#{})", record.id, anchor(&record.id));
    }
    let _ = writeln!(out, "\n**Decided**\n");
    for record in decided {
        match record.summary.as_deref() {
            Some(summary) => {
                let _ = writeln!(
                    out,
                    "- [`{}`](#{}) — {}",
                    record.id,
                    anchor(&record.id),
                    summary
                );
            }
            None => {
                let _ = writeln!(out, "- [`{}`](#{})", record.id, anchor(&record.id));
            }
        }
    }
    let _ = writeln!(out);
    out
}

/// GitHub's anchor for a heading whose whole text is a backticked id.
///
/// Backticks are dropped and the id is already lowercase kebab-case, so the
/// anchor is the id itself. Kept as a function so the assumption has one home
/// and is stated where it is relied on.
fn anchor(id: &str) -> String {
    id.to_string()
}

/// One record.
fn render_record(record: &Record) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "#### `{}`\n", record.id);
    let _ = writeln!(out, "**Status:** {}\n", record.status);
    let _ = writeln!(out, "**The question.** {}\n", record.question.trim());

    if !record.applies_to.is_empty() {
        let _ = writeln!(
            out,
            "**Applies to.** {}\n",
            record
                .applies_to
                .iter()
                .map(|d| format!("`{d}`"))
                .collect::<Vec<_>>()
                .join(", ")
        );
    }

    // A house choice and an open question both reach the loop below with no
    // governing clause, and the reason differs completely — one has chosen and
    // has no clause to choose *under*, the other has a clause and has not
    // chosen. The published note must not tell the reader the wrong one, so the
    // house-choice case is stated here and the empty-governing note is suppressed
    // for it.
    let empty_governing_note = match &record.house_choice {
        Some(choice) => {
            let _ = writeln!(
                out,
                "**House choice.** This ruling is the project's own, made under {} where the \
                 recommendations do not decide. It cites no governing clause and must never be \
                 quoted as conformance.\n",
                choice.under.label()
            );
            let _ = writeln!(
                out,
                "**Considered and rejected.** {}\n",
                collapse(&choice.considered_and_rejected)
            );
            None
        }
        None => Some(
            "**Governing clause.** None. An open record names a conflict without choosing a \
             side, and the generator refuses to build one that names an authority.",
        ),
    };

    for (role, heading, empty_note) in [
        (
            Role::Governing,
            "**Governing clause.**",
            empty_governing_note,
        ),
        (Role::DeviatesFrom, "**Deviates from.**", None),
        (Role::Cited, "**Also cited.**", None),
    ] {
        let clauses: Vec<&_> = record.citations.iter().filter(|c| c.role == role).collect();
        if clauses.is_empty() {
            if let Some(note) = empty_note {
                let _ = writeln!(out, "{note}\n");
            }
            continue;
        }
        let _ = writeln!(out, "{heading}\n");
        for citation in clauses {
            let _ = writeln!(
                out,
                "- `{}`\n  > {}",
                citation.clause,
                collapse(&citation.quote)
            );
        }
        let _ = writeln!(out);
    }

    if !record.equivalence_classes.is_empty() {
        let _ = writeln!(
            out,
            "**Convergence pinned by equivalence class.** {}\n",
            record
                .equivalence_classes
                .iter()
                .map(|c| format!("`{c}`"))
                .collect::<Vec<_>>()
                .join(", ")
        );
    }

    let _ = writeln!(out, "**The ruling.**\n");
    for paragraph in record.rationale.trim().split("\n\n") {
        let _ = writeln!(out, "{}\n", collapse(paragraph));
    }
    out
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

/// Every record reaches the document, under a heading that states its status.
///
/// The rendering groups by status, so a record with an unexpected status would
/// be dropped rather than mis-filed. `render` asserts the partition; this
/// asserts the result, which is the property a reader depends on.
#[test]
fn every_record_is_published_with_its_status() {
    let rendered = render();
    let records = records();
    assert!(
        !records.is_empty(),
        "the ledger holds no records, so this document would be vacuous"
    );

    let statuses: BTreeMap<&str, &str> = records
        .iter()
        .map(|r| (r.id.as_str(), r.status.as_str()))
        .collect();
    for (id, status) in &statuses {
        let heading = format!("#### `{id}`\n");
        assert!(
            rendered.contains(&heading),
            "record {id} has no heading in {DOC_RELATIVE_PATH}"
        );
        let body = rendered
            .split_once(&heading)
            .expect("the heading was just found")
            .1;
        assert!(
            body.starts_with(&format!("\n**Status:** {status}\n")),
            "record {id} is published without its `{status}` status"
        );
    }

    let (open_half, decided_half) = rendered
        .split_once("### Decided")
        .expect("the document has a Decided section");
    for (id, status) in &statuses {
        let half = if *status == "decided" {
            decided_half
        } else {
            open_half
        };
        assert!(
            half.contains(&format!("#### `{id}`")),
            "record {id} is `{status}` but is not rendered in the `{status}` half"
        );
    }
}

/// Every record's question and ruling reach the reader, not just its id.
///
/// A document that published ids and clause citations alone would look complete
/// and answer nothing — which is the state #1552 was filed about.
#[test]
fn every_record_publishes_its_question_and_its_ruling() {
    let rendered = render();
    for record in records() {
        assert!(
            rendered.contains(&collapse(&record.question)),
            "record {} is published without its question",
            record.id
        );
        let first_paragraph = collapse(
            record
                .rationale
                .trim()
                .split("\n\n")
                .next()
                .expect("split always yields one element"),
        );
        assert!(
            rendered.contains(&first_paragraph),
            "record {} is published without the opening of its ruling",
            record.id
        );
        for citation in &record.citations {
            assert!(
                rendered.contains(&format!("`{}`", citation.clause)),
                "record {} cites {} without publishing it",
                record.id,
                citation.clause
            );
            assert!(
                rendered.contains(&collapse(&citation.quote)),
                "record {} publishes {} without the text it was quoted against",
                record.id,
                citation.clause
            );
        }
    }
}

/// The contents list carries each decided record's one-sentence summary, so
/// a reader can scan every ruling on one screen. An id states the question,
/// not the answer, and the previous place this scan existed was a
/// hand-maintained agent-guidance table that drifted three times.
#[test]
fn the_contents_list_carries_every_decided_summary() {
    let rendered = render();
    let contents_section = rendered
        .split("## Contents")
        .nth(1)
        .and_then(|rest| rest.split(RECORDS_HEADING).next())
        .expect("the contents list sits between its heading and the records heading");
    for record in records().iter().filter(|r| r.status == "decided") {
        let summary = record
            .summary
            .as_deref()
            .unwrap_or_else(|| panic!("decided record {} has no summary", record.id));
        let line = format!("- [`{}`](#{}) — {}", record.id, anchor(&record.id), summary);
        assert!(
            contents_section.contains(&line),
            "contents list is missing the summary line for {}",
            record.id
        );
    }
}

/// The preamble points at the README ruleset and does not restate it.
///
/// Single-sourcing that ruleset is part of a ruling this very document renders,
/// so a preamble that restated it would contradict its own contents. The check
/// is against `README.md` itself — a copy of the rule openers kept here would be
/// the third copy of the thing being guarded.
///
/// Scoped to the preamble deliberately. Below [`RECORDS_HEADING`] the text is
/// the ledger's own, and a record is entitled to quote the ruleset it rules
/// under; policing that would be policing the ledger from a renderer.
#[test]
fn the_preamble_does_not_restate_the_readme_ruleset() {
    let rendered = render();
    let preamble = rendered
        .split_once(RECORDS_HEADING)
        .expect("the document has a records section")
        .0;

    assert!(
        preamble.contains("src/reference/normalization-rules.md"),
        "the preamble must send the reader to the ruleset page rather than restating it"
    );

    let page = read("docs/src/reference/normalization-rules.md");
    let section = page
        .split_once("# Normalization rules\n")
        .expect("the ruleset page has a title heading")
        .1;

    let openers: Vec<&str> = section
        .lines()
        .filter_map(|line| {
            let line = line.trim_start();
            let rest = line.strip_prefix(char::is_numeric)?.strip_prefix(". **")?;
            rest.split_once("**").map(|(opener, _)| opener)
        })
        .collect();
    assert!(
        openers.len() >= 7,
        "found {} rule openers in the README ruleset, expected the seven rules — the section's \
         shape moved and this guard is no longer reading it",
        openers.len()
    );
    for opener in openers {
        assert!(
            !preamble.contains(&format!("**{opener}**")),
            "the preamble restates the README ruleset's rule {opener:?}; link to it instead"
        );
    }
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
