//! Generates and gates `docs/NORMALIZATION_STAGE_AUDIT.md` — the stage-by-stage
//! inventory of what governs each decision the normalizer makes (#1553).
//!
//! # The question it answers
//!
//! Ferro's normalization pipeline makes a decision at every stage, and only some
//! of those decisions are backed by a governing clause or a ruling record. The
//! rest are backed by whatever the implementation happens to do. Without an
//! inventory saying which is which, a contributor cannot tell a deliberate
//! choice from an accident, and neither can a reviewer.
//!
//! The audit is `tests/fixtures/grammar/normalization_stage_audit.json`: one row
//! per stage, naming the decision the stage makes, the records and clauses that
//! govern it, and — the deliverable — the decisions that **nothing** governs,
//! each classified as a genuine gap or as a case the spec settles so plainly no
//! record is warranted.
//!
//! # What is curated and what is checked
//!
//! The assignment of a record to a stage is editorial. Everything a reader might
//! act on is not:
//!
//! * every `rulings` id must exist in the adjudication ledger
//!   ([`every_cited_ruling_exists_in_the_ledger`]);
//! * every clause must resolve to real lines of the pinned spec checkout, and
//!   the quote must be found there
//!   ([`every_cited_clause_quote_is_in_the_spec_checkout`]);
//! * the stage list must be exactly the closed set the issue enumerated, in
//!   pipeline order ([`the_audit_covers_every_stage_exactly_once`]).
//!
//! **No status is written down in the fixture.** A row cites a record; whether
//! that record is decided or open is read from the ledger when the document is
//! rendered. That is deliberate — a copied status is the drift this repository
//! has been bitten by more than once, and here it would be worse than usual,
//! because "this stage is settled" is exactly the claim the audit exists to
//! test.
//!
//! # Relationship to the sibling document
//!
//! This inventory says *which* record governs a stage. The records themselves —
//! question, verdict, clauses, quotes and reasoning — are published separately.
//! Neither document restates `README.md`'s ruleset; that is stated once, there,
//! and both link to it.
//!
//! Regenerate with:
//!
//! ```text
//! BLESS_STAGE_AUDIT_DOC=1 cargo nextest run --features dev --test it \
//!   -E 'test(normalization_stage_audit_doc)'
//! ```

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::PathBuf;

use super::common::rulings::records;

/// Where the curated audit lives, relative to the crate root.
const AUDIT_RELATIVE_PATH: &str = "tests/fixtures/grammar/normalization_stage_audit.json";

/// Where the generated document lives, relative to the crate root.
const DOC_RELATIVE_PATH: &str = "docs/NORMALIZATION_STAGE_AUDIT.md";

/// The pinned upstream spec checkout, which every clause citation resolves into.
const SPEC_DIR: &str = "assets/hgvs-nomenclature";

/// The environment variable that turns [`the_published_document_is_current`]
/// from a comparison into a write.
const BLESS_VAR: &str = "BLESS_STAGE_AUDIT_DOC";

/// The stages, in pipeline order.
///
/// A closed list, held in code rather than inferred from the fixture, so a stage
/// cannot be dropped from the audit by deleting its row — which is the one edit
/// that would make the inventory read as complete while covering less. The set
/// is the one #1553 enumerates.
const STAGES: [&str; 7] = [
    "parse-validate",
    "axis-projection",
    "sequence-derivation",
    "partition-and-typing",
    "shifting",
    "merge-split",
    "re-spelling",
];

/// The verdicts an ungoverned decision may carry.
///
/// The distinction is the whole point of recording an ungoverned decision at
/// all: #1553 asks for it explicitly, because "nothing governs this" is a
/// finding only when the spec has not in fact settled it plainly.
const UNGOVERNED_VERDICTS: [&str; 2] = ["genuine-gap", "no-record-warranted"];

/// How normalization-relevant an axis-scoped clause is.
///
/// `unmodellable` is a first-class answer, not a synonym for `open`: a clause
/// keyed on provenance or on population frequency is not a function of the
/// sequences a normalizer holds, so no amount of work closes it. Leaving such a
/// clause blank would read as a gap someone could close.
const RELEVANCES: [&str; 4] = [
    "modelled",
    "open",
    "not-a-normalization-rule",
    "unmodellable",
];

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn read(relative: &str) -> String {
    let path = repo_root().join(relative);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

// --------------------------------------------------------------------------
// The curated audit.
// --------------------------------------------------------------------------

/// One clause a row cites, with the text it is cited for.
struct Clause {
    clause: String,
    quote: String,
}

/// A decision at a stage that nothing governs.
struct Ungoverned {
    decision: String,
    verdict: String,
    why: String,
}

/// One pipeline stage.
struct Stage {
    id: String,
    title: String,
    decision: String,
    rulings: Vec<String>,
    clauses: Vec<Clause>,
    ungoverned: Vec<Ungoverned>,
    note: Option<String>,
}

/// One clause whose reach is scoped to particular axes.
struct AxisClause {
    clause: String,
    quote: String,
    scope: String,
    relevance: String,
    note: Option<String>,
}

/// The whole curated audit.
struct Audit {
    stages: Vec<Stage>,
    axis_scoped_clauses: Vec<AxisClause>,
}

/// A required string field, reported by name.
fn field<'a>(value: &'a serde_json::Value, name: &str, subject: &str) -> &'a str {
    match value.get(name) {
        None | Some(serde_json::Value::Null) => panic!("{subject} has no `{name}`"),
        Some(found) => found
            .as_str()
            .unwrap_or_else(|| panic!("{subject} has a non-string `{name}`: {found}")),
    }
}

/// An optional string field. Absent, `null` and a blank string all read as
/// absent — a note present but empty renders as an empty paragraph, which looks
/// like a dropped sentence rather than like no note.
fn optional_field(value: &serde_json::Value, name: &str, subject: &str) -> Option<String> {
    match value.get(name) {
        None | Some(serde_json::Value::Null) => None,
        Some(found) => {
            let text = found
                .as_str()
                .unwrap_or_else(|| panic!("{subject} has a non-string `{name}`: {found}"));
            (!text.trim().is_empty()).then(|| text.to_string())
        }
    }
}

fn string_list(value: &serde_json::Value, name: &str, subject: &str) -> Vec<String> {
    match value.get(name) {
        None | Some(serde_json::Value::Null) => Vec::new(),
        Some(found) => found
            .as_array()
            .unwrap_or_else(|| panic!("{subject} has a non-array `{name}`: {found}"))
            .iter()
            .map(|entry| {
                entry
                    .as_str()
                    .unwrap_or_else(|| panic!("{subject} has a non-string `{name}` entry: {entry}"))
                    .to_string()
            })
            .collect(),
    }
}

fn clauses(value: &serde_json::Value, subject: &str) -> Vec<Clause> {
    match value.get("clauses") {
        None | Some(serde_json::Value::Null) => Vec::new(),
        Some(found) => found
            .as_array()
            .unwrap_or_else(|| panic!("{subject} has a non-array `clauses`: {found}"))
            .iter()
            .map(|entry| Clause {
                clause: field(entry, "clause", subject).to_string(),
                quote: field(entry, "quote", subject).to_string(),
            })
            .collect(),
    }
}

/// Parse the curated audit. Panics on anything malformed, for the reason
/// `common::rulings` gives: a broken fixture is not a case to degrade through.
fn audit() -> Audit {
    let text = read(AUDIT_RELATIVE_PATH);
    let value: serde_json::Value =
        serde_json::from_str(&text).unwrap_or_else(|e| panic!("parse {AUDIT_RELATIVE_PATH}: {e}"));

    let stages = value
        .get("stages")
        .and_then(|s| s.as_array())
        .unwrap_or_else(|| panic!("{AUDIT_RELATIVE_PATH} has no `stages` array"))
        .iter()
        .map(|entry| {
            let id = field(entry, "id", "a stage").to_string();
            let subject = format!("stage {id}");
            Stage {
                title: field(entry, "title", &subject).to_string(),
                decision: field(entry, "decision", &subject).to_string(),
                rulings: string_list(entry, "rulings", &subject),
                clauses: clauses(entry, &subject),
                ungoverned: match entry.get("ungoverned") {
                    None | Some(serde_json::Value::Null) => Vec::new(),
                    Some(found) => found
                        .as_array()
                        .unwrap_or_else(|| panic!("{subject} has a non-array `ungoverned`"))
                        .iter()
                        .map(|u| {
                            let verdict = field(u, "verdict", &subject).to_string();
                            assert!(
                                UNGOVERNED_VERDICTS.contains(&verdict.as_str()),
                                "{subject} records an ungoverned decision with verdict \
                                 {verdict:?}, which is not one of {UNGOVERNED_VERDICTS:?} — the \
                                 distinction between a genuine gap and a case the spec settles \
                                 plainly is what makes the row a finding"
                            );
                            Ungoverned {
                                decision: field(u, "decision", &subject).to_string(),
                                verdict,
                                why: field(u, "why", &subject).to_string(),
                            }
                        })
                        .collect(),
                },
                note: optional_field(entry, "note", &subject),
                id,
            }
        })
        .collect();

    let axis_scoped_clauses = value
        .get("axis_scoped_clauses")
        .and_then(|s| s.as_array())
        .unwrap_or_else(|| panic!("{AUDIT_RELATIVE_PATH} has no `axis_scoped_clauses` array"))
        .iter()
        .map(|entry| {
            let clause = field(entry, "clause", "an axis-scoped clause").to_string();
            let subject = format!("axis-scoped clause {clause}");
            let relevance = field(entry, "relevance", &subject).to_string();
            assert!(
                RELEVANCES.contains(&relevance.as_str()),
                "{subject} has relevance {relevance:?}, which is not one of {RELEVANCES:?}"
            );
            AxisClause {
                quote: field(entry, "quote", &subject).to_string(),
                scope: field(entry, "scope", &subject).to_string(),
                note: optional_field(entry, "note", &subject),
                relevance,
                clause,
            }
        })
        .collect();

    Audit {
        stages,
        axis_scoped_clauses,
    }
}

// --------------------------------------------------------------------------
// Resolving a clause into the spec checkout.
// --------------------------------------------------------------------------

/// The cited lines of the spec checkout, joined by spaces.
///
/// `docs/recommendations/general.md:34` and `…:79-84` are both accepted; a range
/// is joined, which is what makes a quote spanning several lines resolvable.
/// This mirrors what the spec-fixture generator does for the ledger's own
/// citations, and shares its limit: the comparison is whitespace-collapsed
/// containment, not a byte-for-byte match. It catches a clause moving out from
/// under a citation; it does not certify that a quote is exact.
fn cited_text(clause: &str) -> String {
    let (path, span) = clause
        .rsplit_once(':')
        .unwrap_or_else(|| panic!("clause {clause:?} carries no `:line`"));
    let (first, last) = match span.split_once('-') {
        Some((a, b)) => (parse_line(a, clause), parse_line(b, clause)),
        None => {
            let only = parse_line(span, clause);
            (only, only)
        }
    };
    assert!(
        first >= 1 && last >= first,
        "clause {clause:?} names an empty or inverted line range"
    );

    let file = repo_root().join(SPEC_DIR).join(path);
    let text = std::fs::read_to_string(&file).unwrap_or_else(|e| {
        panic!(
            "clause {clause:?} names {}, which cannot be read: {e}. If the spec checkout is \
             empty, initialise it:\n    git submodule update --init {SPEC_DIR}",
            file.display()
        )
    });
    let lines: Vec<&str> = text.lines().collect();
    assert!(
        last <= lines.len(),
        "clause {clause:?} names line {last} of a {}-line file",
        lines.len()
    );
    collapse(&lines[first - 1..last].join(" "))
}

fn parse_line(raw: &str, clause: &str) -> usize {
    raw.parse()
        .unwrap_or_else(|_| panic!("clause {clause:?} has a non-numeric line number {raw:?}"))
}

/// Runs of whitespace collapsed to single spaces.
fn collapse(text: &str) -> String {
    text.split_whitespace().collect::<Vec<_>>().join(" ")
}

// --------------------------------------------------------------------------
// Rendering.
// --------------------------------------------------------------------------

/// `id -> status`, from the ledger. Never from the fixture.
fn ledger_statuses() -> BTreeMap<String, String> {
    records().into_iter().map(|r| (r.id, r.status)).collect()
}

fn render() -> String {
    let audit = audit();
    let statuses = ledger_statuses();

    let open_rulings: BTreeSet<&str> = audit
        .stages
        .iter()
        .flat_map(|s| s.rulings.iter())
        .filter(|id| statuses.get(*id).map(String::as_str) == Some("undecided"))
        .map(String::as_str)
        .collect();
    let gaps: Vec<(&Stage, &Ungoverned)> = audit
        .stages
        .iter()
        .flat_map(|s| s.ungoverned.iter().map(move |u| (s, u)))
        .collect();

    let mut out = String::new();
    out.push_str(&preamble(audit.stages.len(), open_rulings.len(), &gaps));

    let _ = writeln!(out, "## The stages\n");
    for stage in &audit.stages {
        out.push_str(&render_stage(stage, &statuses));
    }

    out.push_str(&render_axis_table(&audit.axis_scoped_clauses));

    let mut out: String = out
        .lines()
        .map(str::trim_end)
        .collect::<Vec<_>>()
        .join("\n");
    out.truncate(out.trim_end().len());
    out.push('\n');
    out
}

fn preamble(stages: usize, open: usize, gaps: &[(&Stage, &Ungoverned)]) -> String {
    let mut out = String::new();
    let _ = write!(
        out,
        r#"<!--
GENERATED FILE — do not edit by hand.

Rendered from tests/fixtures/grammar/normalization_stage_audit.json and the
adjudication ledger, by tests/it/normalization_stage_audit_doc.rs. Edit the
fixture, then regenerate:

    BLESS_STAGE_AUDIT_DOC=1 cargo nextest run --features dev --test it \
      -E 'test(normalization_stage_audit_doc)'
-->

# What governs each normalization stage

Ferro's normalizer makes a decision at every stage of its pipeline. Some of those decisions are
backed by a clause of the HGVS recommendations, some by an adjudication record, and some by
nothing more than what the implementation happens to do. This document is the inventory of
which is which, so that a deliberate choice can be told from an accident — by a contributor
deciding a new case, and by a reviewer reading a diff.

**{stages} stages. {open} of the records cited below are still open. {gap_count} decisions are
governed by nothing.**

## How to read it

Each stage names the decision it makes, the records that rule on it, and the clauses it is
decided under. A record's status is read from the ledger when this file is generated and is
never copied into the fixture, so a stage cannot be described as settled by a record that is
not.

The rows worth reading first are the **ungoverned** ones. Each is classified as either a
*genuine gap* — a decision with real consequences that no clause and no record settles — or as
*no record warranted*, meaning the spec settles it so plainly that writing a record would be
ceremony. Only the first kind is a finding.

## What this document is not

**It is not the ruleset.** What ferro's output is allowed to be, and what happens where the
spec determines no answer, is stated once, in
[README.md, *Normalization rules*](../README.md#normalization-rules), and is deliberately not
restated here.

**It is not the records.** It says which record governs a stage; it does not reproduce the
record. The records live in
[`tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`](../tests/fixtures/grammar/hgvs_spec_normalization_overrides.json).

**It is not exhaustive about clauses.** A stage's clause list carries the clauses that decide
it, not every clause that touches it.

## The ungoverned decisions

"#,
        gap_count = gaps.len()
    );

    for (stage, gap) in gaps {
        let _ = writeln!(
            out,
            "### {} — {}\n\n*{}*\n\n{}\n",
            stage.title,
            gap.verdict,
            collapse(&gap.decision),
            collapse(&gap.why)
        );
    }

    out
}

fn render_stage(stage: &Stage, statuses: &BTreeMap<String, String>) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "### {} (`{}`)\n", stage.title, stage.id);
    let _ = writeln!(out, "**Decides.** {}\n", collapse(&stage.decision));

    if stage.rulings.is_empty() {
        let _ = writeln!(out, "**Ruled by.** No adjudication record.\n");
    } else {
        let _ = writeln!(out, "**Ruled by.**\n");
        for id in &stage.rulings {
            let status = statuses
                .get(id)
                .unwrap_or_else(|| panic!("stage {} cites unknown record {id}", stage.id));
            let _ = writeln!(out, "- `{id}` — {status}");
        }
        let _ = writeln!(out);
    }

    if !stage.clauses.is_empty() {
        let _ = writeln!(out, "**Decided under.**\n");
        for clause in &stage.clauses {
            let _ = writeln!(
                out,
                "- `{}`\n  > {}",
                clause.clause,
                collapse(&clause.quote)
            );
        }
        let _ = writeln!(out);
    }

    if stage.ungoverned.is_empty() {
        let _ = writeln!(
            out,
            "**Governed by nothing.** Nothing recorded at this stage.\n"
        );
    } else {
        let _ = writeln!(out, "**Governed by nothing.**\n");
        for gap in &stage.ungoverned {
            let _ = writeln!(out, "- **{}** — {}", gap.verdict, collapse(&gap.decision));
        }
        let _ = writeln!(out);
    }

    if let Some(note) = &stage.note {
        let _ = writeln!(out, "{}\n", collapse(note));
    }
    out
}

fn render_axis_table(clauses: &[AxisClause]) -> String {
    let mut out = String::new();
    let _ = write!(
        out,
        r#"## The axis dimension

The recommendations are not uniform across axes, and the non-uniformity is deliberate rather
than accidental. `general.md:162` settles that outright: asked whether the `g.` and `c.`
descriptions of one variant may name different nucleotides, the spec answers **yes**, for a
gene on the minus strand inside a repeated sequence.

That matters because a normalizer's axes disagreeing has, in this repository, always turned out
to be a bug — a dropped offset, a skipped pass, a missing gate. Every such case so far was an
*accidental* divergence, and none was a deliberate axis-scoped rule. Without a table recording
which is which, the next deliberate divergence is indistinguishable from the next defect.

Two things shape the table. **Rules are keyed on the property, not on the axis label**: what
matters for the codon carve-out is having a reading frame, which correctly groups `n.` with
`g.` rather than with `c.`, a pairing an axis-name-keyed table would get wrong. And **some
rules cannot be modelled at all** — a clause keyed on provenance or on population frequency is
not a function of the sequences a normalizer holds — so those are recorded as `unmodellable`
rather than left looking unimplemented.

| clause | scope | normalization-relevant? |
|---|---|---|
"#
    );
    for clause in clauses {
        let _ = writeln!(
            out,
            "| `{}` | {} | {} |",
            clause.clause,
            collapse(&clause.scope),
            clause.relevance
        );
    }
    let _ = writeln!(out, "\n### The clauses, with their text\n");
    for clause in clauses {
        let _ = writeln!(
            out,
            "#### `{}` — {}\n\n> {}\n\n**Scope.** {}\n",
            clause.clause,
            clause.relevance,
            collapse(&clause.quote),
            collapse(&clause.scope)
        );
        if let Some(note) = &clause.note {
            let _ = writeln!(out, "{}\n", collapse(note));
        }
    }
    out
}

// --------------------------------------------------------------------------
// Tests.
// --------------------------------------------------------------------------

/// The committed document matches what the audit and the ledger render to.
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
            "read {}: {e}\n\nGenerate it with:\n  {BLESS_VAR}=1 cargo nextest run --features dev \
             --test it -E 'test(normalization_stage_audit_doc)'",
            path.display()
        )
    });
    assert!(
        committed == rendered,
        "{DOC_RELATIVE_PATH} is stale against the audit fixture or the ledger. Regenerate it — \
         do not edit it by hand:\n  {BLESS_VAR}=1 cargo nextest run --features dev --test it \
         -E 'test(normalization_stage_audit_doc)'"
    );
}

/// Every record the audit cites is in the ledger.
///
/// Without this the inventory could name a record that had been renamed or
/// retired, and read as coverage while pointing at nothing.
#[test]
fn every_cited_ruling_exists_in_the_ledger() {
    let statuses = ledger_statuses();
    let audit = audit();
    let mut cited = 0usize;
    for stage in &audit.stages {
        for id in &stage.rulings {
            assert!(
                statuses.contains_key(id),
                "stage {} cites `{id}`, which is not a record in the ledger",
                stage.id
            );
            cited += 1;
        }
    }
    assert!(
        cited > 0,
        "the audit cites no records at all, so this check is vacuous"
    );
}

/// Every clause the audit cites resolves in the spec checkout, and the quote is
/// found at the cited lines.
///
/// A `file:line` nobody checks is the failure mode the ledger's own citations
/// were given quotes to close; an inventory of what governs each stage must not
/// reintroduce it one file over.
#[test]
fn every_cited_clause_quote_is_in_the_spec_checkout() {
    let audit = audit();
    let mut checked = 0usize;

    let mut check = |clause: &str, quote: &str, subject: &str| {
        let cited = cited_text(clause);
        assert!(
            cited.contains(&collapse(quote)),
            "{subject}: the quote for `{clause}` is not at those lines.\n  quoted: {}\n  found:  \
             {cited}",
            collapse(quote)
        );
        checked += 1;
    };

    for stage in &audit.stages {
        for clause in &stage.clauses {
            check(
                &clause.clause,
                &clause.quote,
                &format!("stage {}", stage.id),
            );
        }
    }
    for clause in &audit.axis_scoped_clauses {
        check(&clause.clause, &clause.quote, "the axis table");
    }
    assert!(
        checked > 0,
        "no clause was checked, so this test would pass over an empty audit"
    );
}

/// The audit covers exactly the stages [`STAGES`] names, once each, in order.
///
/// Held in code rather than derived from the fixture: deleting a row is the one
/// edit that makes the inventory read as complete while covering less.
#[test]
fn the_audit_covers_every_stage_exactly_once() {
    let ids: Vec<String> = audit().stages.into_iter().map(|s| s.id).collect();
    assert_eq!(
        ids,
        STAGES.map(String::from).to_vec(),
        "the audit's stages differ from the closed list this pipeline has. Adding a stage is a \
         legitimate answer; dropping one silently is not"
    );
}

/// Every stage says something about every field a reader relies on.
///
/// A row that names a stage and leaves the governance blank is worse than no
/// row: it occupies the place where the answer would be.
#[test]
fn no_stage_row_is_blank() {
    for stage in audit().stages {
        assert!(
            !stage.title.trim().is_empty() && !stage.decision.trim().is_empty(),
            "stage {} has no title or no decision",
            stage.id
        );
        assert!(
            !stage.rulings.is_empty() || !stage.clauses.is_empty() || !stage.ungoverned.is_empty(),
            "stage {} records neither a record, nor a clause, nor an ungoverned decision — the \
             row states nothing",
            stage.id
        );
    }
}

/// Every ungoverned decision the audit records reaches the document, with its
/// verdict.
///
/// These rows are the deliverable. Rendering them only inside their stage would
/// bury the finding in the middle of a long file, so they are also collected at
/// the top, and this asserts both copies exist.
#[test]
fn every_ungoverned_decision_is_published_with_its_verdict() {
    let rendered = render();
    let audit = audit();
    let mut found = 0usize;
    for stage in &audit.stages {
        for gap in &stage.ungoverned {
            let decision = collapse(&gap.decision);
            assert!(
                rendered.matches(decision.as_str()).count() >= 2,
                "the ungoverned decision at stage {} is not published both in the summary and \
                 under its stage: {decision}",
                stage.id
            );
            assert!(
                rendered.contains(&collapse(&gap.why)),
                "the ungoverned decision at stage {} is published without its reasoning",
                stage.id
            );
            assert!(
                rendered.contains(&format!("### {} — {}", stage.title, gap.verdict)),
                "the ungoverned decision at stage {} is published without its verdict",
                stage.id
            );
            found += 1;
        }
    }
    assert!(
        found > 0,
        "the audit records no ungoverned decision anywhere. That is a possible answer, but it is \
         the one that would make this document pointless, so it must be arrived at deliberately \
         rather than by an empty fixture"
    );
}

/// The document does not restate the README ruleset.
///
/// Same reason as its sibling: single-sourcing that ruleset is part of a ruling
/// in the ledger this audit indexes.
#[test]
fn the_document_does_not_restate_the_readme_ruleset() {
    let rendered = render();
    assert!(
        rendered.contains("../README.md#normalization-rules"),
        "the document must send the reader to the README ruleset rather than restating it"
    );

    let readme = read("README.md");
    let section = readme
        .split_once("\n## Normalization rules\n")
        .expect("README has a Normalization rules section")
        .1
        .split_once("\n## ")
        .expect("the section is followed by another")
        .0;
    let openers: Vec<&str> = section
        .lines()
        .filter_map(|line| {
            let rest = line
                .trim_start()
                .strip_prefix(char::is_numeric)?
                .strip_prefix(". **")?;
            rest.split_once("**").map(|(opener, _)| opener)
        })
        .collect();
    assert!(
        openers.len() >= 7,
        "found {} rule openers in the README ruleset, expected seven — the section's shape moved \
         and this guard is no longer reading it",
        openers.len()
    );
    for opener in openers {
        assert!(
            !rendered.contains(&format!("**{opener}**")),
            "the document restates the README ruleset's rule {opener:?}; link to it instead"
        );
    }
}

/// The render is already in the shape the file-hygiene pre-commit hooks want, so
/// committing it cannot rewrite it into a failure that reads as drift.
#[test]
fn the_render_survives_the_file_hygiene_hooks() {
    let rendered = render();
    assert!(
        rendered.ends_with('\n') && !rendered.ends_with("\n\n"),
        "the render must end with exactly one newline"
    );
    for (number, line) in rendered.lines().enumerate() {
        assert_eq!(
            line.trim_end(),
            line,
            "line {} of the render carries trailing whitespace: {line:?}",
            number + 1
        );
    }
}
