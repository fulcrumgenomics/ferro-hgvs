//! Reads the adjudication ledger — the `rulings` section of
//! `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`.
//!
//! Two integration tests need the same file for different reasons:
//! `ruling_citation_currency.rs` wants `id -> status` so it can judge prose that
//! makes a status claim, and `clause_ruling_index.rs` wants the citations so it
//! can invert the ledger into a clause-to-record index. Parsing it twice would
//! give two definitions of "a record", so both read it from here.
//!
//! Nothing in this module names a record id, deliberately: it is scanned by
//! `ruling_citation_currency.rs` along with every other `.rs` file in the tree.
//! The parser's own tests at the bottom hold to that too. Their fixture id has
//! to *be* id-shaped now that the parser enforces the grammar, so it is kept out
//! of the scan's way instead: the scan only reads a token as a citation when the
//! line also says "ruling" and the token is backticked, or when it appears in the
//! `rulings[…]` form, and the fixture's id is a bare JSON string on neither kind
//! of line.

use std::collections::BTreeMap;
use std::path::PathBuf;

/// Where the ledger lives, relative to the crate root.
pub const LEDGER_RELATIVE_PATH: &str =
    "tests/fixtures/grammar/hgvs_spec_normalization_overrides.json";

/// How a record relates to one of the clauses it names.
///
/// A record may govern at most one clause and may deviate from several; every
/// other clause it names is context for the conflict it records.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Role {
    /// The clause the record rules governs the conflict.
    Governing,
    /// A clause the record deliberately departs from.
    DeviatesFrom,
    /// Named as part of the conflict, without a verdict attached.
    Cited,
}

impl Role {
    /// Lower-case label used in rendered output.
    pub fn label(self) -> &'static str {
        match self {
            Role::Governing => "governing",
            Role::DeviatesFrom => "deviates-from",
            Role::Cited => "cited",
        }
    }
}

/// One clause a record names, as the record spells it.
///
/// `clause` keeps the record's own spelling, including a line range
/// (`…/delins.md:79-84`). Callers that need per-line lookup expand it
/// themselves; see `clause_ruling_index.rs`.
#[derive(Clone, Debug)]
pub struct Citation {
    /// The clause as the record spells it — a repo-relative path into the spec
    /// checkout with a line or line range (`docs/…/delins.md:79-84`).
    pub clause: String,
    /// How the record relates to this clause; see [`Role`].
    pub role: Role,
}

/// One adjudication record.
#[derive(Clone, Debug)]
pub struct Record {
    /// The record's `id`, the kebab-case token prose cites it by.
    pub id: String,
    /// `decided` or `undecided`, verbatim from the ledger. No other value
    /// parses; see [`STATUSES`].
    pub status: String,
    /// The record's `rationale` prose, verbatim.
    ///
    /// Exposed because a rationale cites *other records* — that is how the
    /// ledger cross-references itself — and those citations are not covered by
    /// the source scan in `ruling_citation_currency.rs`, which reads `.rs` files
    /// only. See `every_record_to_record_citation_resolves` there.
    pub rationale: String,
    /// Every clause the record names, in the ledger's own order, each tagged
    /// with the role the record's verdict fields give it.
    pub citations: Vec<Citation>,
}

/// The only two `status` values a record may carry.
///
/// `ruling_citation_currency.rs` judges prose against exactly these two words,
/// so a third value would be a status no citation could ever agree with.
pub const STATUSES: &[&str] = &["decided", "undecided"];

/// Whether `token` has the shape of a record id: three or more hyphen-separated
/// segments of lowercase alphanumerics.
///
/// Defined here rather than in `ruling_citation_currency.rs` because both sides
/// need the *same* grammar and a rule kept in two copies drifts. That scan uses
/// it to decide whether a backticked token in prose is citing a record, which is
/// why [`records_from_value`] enforces it on the ledger side too: an id outside
/// this shape can never be *recognised* as a citation, so every stale reference
/// to it would sail past the existence check. The shortest real id
/// (three segments) is the floor; two-word prose must not match.
pub fn looks_like_a_record_id(token: &str) -> bool {
    let segments: Vec<&str> = token.split('-').collect();
    segments.len() >= 3
        && segments.iter().all(|s| {
            !s.is_empty()
                && s.chars()
                    .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit())
        })
}

/// Absolute path to the ledger.
pub fn ledger_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(LEDGER_RELATIVE_PATH)
}

/// Every record, in the ledger's own order.
///
/// Panics rather than returning a `Result`: a malformed ledger is a broken
/// fixture, and every caller is a test whose only sensible response is to fail
/// loudly at the point of the malformation.
pub fn records() -> Vec<Record> {
    let path = ledger_path();
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let value: serde_json::Value =
        serde_json::from_str(&text).unwrap_or_else(|e| panic!("parse {}: {e}", path.display()));
    records_from_value(&value, &path.display().to_string())
}

/// An optional string field.
///
/// Absent or `null` means "not set"; **any other** non-string value panics
/// rather than reading as absent. A silently-dropped verdict field is the
/// dangerous shape here: it turns a malformed claim into a missing one, so the
/// clause index renders `cited` where the ledger meant `governing` and nothing
/// fails.
fn optional_str<'a>(record: &'a serde_json::Value, field: &str, id: &str) -> Option<&'a str> {
    match record.get(field) {
        None | Some(serde_json::Value::Null) => None,
        Some(value) => Some(value.as_str().unwrap_or_else(|| {
            panic!("record {id} has a non-string `{field}`: {value}");
        })),
    }
}

/// A required string field.
///
/// `subject` names what is being parsed, for the panic message. Absent and
/// present-but-not-a-string are reported **differently**: "has no `id`" sends
/// the reader looking for a missing key, which is the wrong hunt when the key is
/// there holding a number.
fn required_str<'a>(value: &'a serde_json::Value, field: &str, subject: &str) -> &'a str {
    match value.get(field) {
        None | Some(serde_json::Value::Null) => panic!("{subject} has no `{field}`"),
        Some(found) => found
            .as_str()
            .unwrap_or_else(|| panic!("{subject} has a non-string `{field}`: {found}")),
    }
}

/// Every record in a parsed ledger document.
///
/// Split from [`records`] so the malformed-input cases below can drive it
/// without a file on disk.
fn records_from_value(value: &serde_json::Value, origin: &str) -> Vec<Record> {
    let rulings = value
        .get("rulings")
        .and_then(|r| r.as_array())
        .unwrap_or_else(|| panic!("{origin} has no `rulings` array"));

    let records: Vec<Record> = rulings
        .iter()
        .map(|record| {
            let id = required_str(record, "id", &format!("a record in {origin}")).to_string();
            assert!(
                looks_like_a_record_id(&id),
                "record id {id:?} is not id-shaped — see `looks_like_a_record_id`. \
                 `ruling_citation_currency.rs` recognises a citation *by* that shape, so an id \
                 outside it could never be matched in prose and every stale reference to this \
                 record would pass the existence check unnoticed"
            );
            let status = required_str(record, "status", &format!("record {id}")).to_string();
            assert!(
                STATUSES.contains(&status.as_str()),
                "record {id} has status {status:?}, which is not one of {STATUSES:?} — the \
                 citation-currency scan judges prose against exactly those two words, so a third \
                 value is a status no citation could agree with"
            );

            let rationale = required_str(record, "rationale", &format!("record {id}")).to_string();

            let governing = optional_str(record, "governing", &id);
            let deviates_from: Vec<&str> = match record.get("deviates_from") {
                None | Some(serde_json::Value::Null) => Vec::new(),
                Some(value) => value
                    .as_array()
                    .unwrap_or_else(|| {
                        panic!("record {id} has a non-array `deviates_from`: {value}")
                    })
                    .iter()
                    .map(|entry| {
                        entry.as_str().unwrap_or_else(|| {
                            panic!("record {id} has a non-string `deviates_from` entry: {entry}")
                        })
                    })
                    .collect(),
            };

            // A clause cannot be both the authority and the departure. Without
            // this, the role selection below resolves the contradiction silently
            // — `governing` is tested first, so it just wins.
            if let Some(governing) = governing {
                assert!(
                    !deviates_from.contains(&governing),
                    "record {id} both governs and deviates from {governing}"
                );
            }

            let clauses = record
                .get("clauses")
                .and_then(|v| v.as_array())
                .unwrap_or_else(|| panic!("record {id} has no `clauses` array"));
            let citations: Vec<Citation> = clauses
                .iter()
                .map(|clause| {
                    let clause =
                        required_str(clause, "clause", &format!("a clause of record {id}"))
                            .to_string();
                    let role = if governing == Some(clause.as_str()) {
                        Role::Governing
                    } else if deviates_from.contains(&clause.as_str()) {
                        Role::DeviatesFrom
                    } else {
                        Role::Cited
                    };
                    Citation { clause, role }
                })
                .collect();

            // The index built from `clauses` is only complete if the verdict
            // fields cannot name a clause the array omits.
            if let Some(governing) = governing {
                assert!(
                    citations.iter().any(|c| c.clause == governing),
                    "record {id} governs {governing}, which is not in its `clauses` array — any \
                     index built from `clauses` would miss it"
                );
            }
            for deviated in &deviates_from {
                assert!(
                    citations.iter().any(|c| c.clause == *deviated),
                    "record {id} deviates from {deviated}, which is not in its `clauses` array — \
                     any index built from `clauses` would miss it"
                );
            }

            Record {
                id,
                status,
                rationale,
                citations,
            }
        })
        .collect();

    let mut seen = BTreeMap::new();
    for record in &records {
        assert!(
            seen.insert(record.id.clone(), ()).is_none(),
            "duplicate record id {} in {origin}",
            record.id
        );
    }
    assert!(
        !records.is_empty(),
        "{origin} lists no records — every scan built on it would be vacuous"
    );
    records
}

/// `id -> status`, for callers that only need the verdict.
pub fn statuses() -> BTreeMap<String, String> {
    records()
        .into_iter()
        .map(|r| (r.id, r.status))
        .collect::<BTreeMap<_, _>>()
}

// --------------------------------------------------------------------------
// Parser tests.
//
// Plain `#[test]`, deliberately **not** inside a `#[cfg(test)] mod tests`: this
// tree is an integration-test binary, which compiles without `cfg(test)`, so a
// gated module would never run and would read as coverage it does not provide
// (see the repository `CLAUDE.md` on committed tests that have never executed).
//
// A well-formed document plus one mutation per field. The mutations are the
// point: every one of them used to parse as "field absent", which converts a
// malformed claim into a missing one and lets the clause index render a role
// the ledger never wrote.
// --------------------------------------------------------------------------

/// A minimal well-formed document: one record, one clause, a verdict on it.
fn one_valid_record() -> serde_json::Value {
    serde_json::json!({
        "rulings": [{
            "id": "a-test-record",
            "status": "decided",
            "rationale": "A test record, with prose that cites no other record.",
            "governing": "docs/a.md:1",
            "clauses": [{ "clause": "docs/a.md:1" }],
        }],
    })
}

/// The positive control. Without it a `should_panic` below could pass for an
/// unrelated reason — a missing field in the fixture rather than the mutation.
#[test]
fn a_well_formed_document_parses() {
    let records = records_from_value(&one_valid_record(), "<test>");
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].status, "decided");
    assert_eq!(records[0].citations[0].role, Role::Governing);
}

/// An id outside the citation grammar is rejected. Not cosmetic: the scan
/// recognises a citation by that shape, so such a record could be cited staler
/// and staler without the existence check ever firing.
#[test]
#[should_panic(expected = "is not id-shaped")]
fn an_id_outside_the_citation_grammar_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["id"] = serde_json::json!("Not An Id");
    records_from_value(&document, "<test>");
}

/// A clause cannot be both the authority and the departure. The role selection
/// would otherwise resolve the contradiction silently in favour of governing.
#[test]
#[should_panic(expected = "both governs and deviates from")]
fn a_clause_that_both_governs_and_is_deviated_from_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["deviates_from"] = serde_json::json!(["docs/a.md:1"]);
    records_from_value(&document, "<test>");
}

/// A present-but-non-string required field reports itself as non-string. Not
/// cosmetic: "has no `id`" sends the reader hunting for a key that is right
/// there holding a number.
#[test]
#[should_panic(expected = "has a non-string `id`")]
fn a_non_string_required_field_is_not_reported_as_missing() {
    let mut document = one_valid_record();
    document["rulings"][0]["id"] = serde_json::json!(17);
    records_from_value(&document, "<test>");
}

/// A missing `rationale` is rejected rather than read as empty prose.
///
/// Empty prose would silently satisfy `every_record_to_record_citation_resolves`
/// in `ruling_citation_currency.rs`: a record with no rationale cites nothing,
/// so it can never fail, and a ledger that lost its prose would still pass the
/// scan that exists to read it.
#[test]
#[should_panic(expected = "has no `rationale`")]
fn a_missing_rationale_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]
        .as_object_mut()
        .expect("record is an object")
        .remove("rationale");
    records_from_value(&document, "<test>");
}

/// A `status` outside [`STATUSES`] is rejected, not carried through.
#[test]
#[should_panic(expected = "which is not one of")]
fn an_unknown_status_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["status"] = serde_json::json!("pending");
    records_from_value(&document, "<test>");
}

/// A non-string `governing` is rejected rather than read as absent.
#[test]
#[should_panic(expected = "non-string `governing`")]
fn a_non_string_governing_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["governing"] = serde_json::json!(17);
    records_from_value(&document, "<test>");
}

/// A `deviates_from` that is not an array is rejected.
#[test]
#[should_panic(expected = "non-array `deviates_from`")]
fn a_non_array_deviates_from_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["deviates_from"] = serde_json::json!("docs/a.md:1");
    records_from_value(&document, "<test>");
}

/// A `deviates_from` array with a non-string entry is rejected, rather than
/// dropping just that entry.
#[test]
#[should_panic(expected = "non-string `deviates_from` entry")]
fn a_mixed_type_deviates_from_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["clauses"] = serde_json::json!([
        { "clause": "docs/a.md:1" },
        { "clause": "docs/b.md:2" },
    ]);
    document["rulings"][0]["deviates_from"] = serde_json::json!(["docs/b.md:2", 17]);
    records_from_value(&document, "<test>");
}

/// An explicit `null` verdict field still means "not set" — that is a
/// legitimate JSON spelling of absent, and the ledger's own records simply omit
/// the keys.
#[test]
fn a_null_verdict_field_reads_as_absent() {
    let mut document = one_valid_record();
    document["rulings"][0]["governing"] = serde_json::Value::Null;
    document["rulings"][0]["deviates_from"] = serde_json::Value::Null;
    let records = records_from_value(&document, "<test>");
    assert_eq!(records[0].citations[0].role, Role::Cited);
}
