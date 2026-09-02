//! Reads the adjudication ledger — the `rulings` section of
//! `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`.
//!
//! Four integration tests need the same file for different reasons:
//! `ruling_citation_currency.rs` wants `id -> status` so it can judge prose that
//! makes a status claim, `clause_ruling_index.rs` wants the citations so it
//! can invert the ledger into a clause-to-record index,
//! `ledger_clause_jurisdiction.rs` wants the citations, their quotes and
//! `applies_to` so it can compare what a record rules on against the molecule
//! directories it cites, and `ruling_guard_field.rs` wants the `guard` field so
//! it can resolve each record's declared enforcement against the tree. Parsing
//! it four times would give four definitions of "a record", so all four read it
//! from here.
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
    /// The clause's quoted text, verbatim from the ledger.
    ///
    /// Required rather than optional, matching what `generate_spec_fixture`
    /// already demands: a citation with no quote is a `file:line` nothing checks
    /// against the spec checkout. `ledger_clause_jurisdiction.rs` also reads it,
    /// so an absent quote would silently narrow the authority a clause supplies.
    pub quote: String,
    /// How the record relates to this clause; see [`Role`].
    pub role: Role,
}

/// What enforces a record's ruling — or the record's stated reason that nothing
/// does.
///
/// Modelled on the `Representation-Change:` trailer, where declining is a
/// first-class answer and what is rejected is *silence*. A record may name the
/// tests that fail when its ruling stops holding, or it may say in words that no
/// test enforces it and why. What it may not do is leave the question open: an
/// absent field is indistinguishable from an unconsidered one, and reads as
/// enforcement that is not there.
///
/// The precedent is [`Record::equivalence_classes`], the ledger's only other
/// structured pointer at pinned evidence. This is the same move for the much
/// larger population of records whose evidence is an ordinary test.
///
/// # Silence is REPRESENTED here, and rejected one layer out
///
/// [`Guard::Absent`] exists so that a record which never mentions `guard` still
/// *parses*. That is not a softening of the rule — the rule is enforced by
/// [`ruling_guard_field::every_record_declares_a_guard`], which names the
/// offending records and the fix. It is a change of *stage*, and the stage is
/// the whole point.
///
/// Refusing at parse made silence fatal to every consumer of the ledger rather
/// than to the check that cares. Measured on this branch, rebased onto
/// `origin/main` at `426a944b`, with the two records `main` had added carrying
/// no `guard`: **50** tests failed — 28 on this reader's own panic, 20 because
/// `generate_spec_fixture` could not deserialize the ledger at all (so the
/// generated fixture never existed and everything downstream of it died too),
/// and 2 in the generator's own precondition suite. **Five** of the fifty were
/// about guards. Seven open PRs write ledger records at any given time, so every
/// one of them met that wall on rebase, and what it told them was that the spec
/// fixture was broken.
///
/// The variant is deliberately a **third state** rather than a reuse of an
/// existing one. Mapping absence onto `Tests(vec![])` or `Declined(String::new())`
/// would make it indistinguishable from a record whose author deliberately wrote
/// an empty list or a blank reason — and those two are separately rejected, with
/// separate messages, precisely because they are different mistakes. No JSON a
/// record can carry produces `Absent`: it is reachable only by the key not being
/// there.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Guard {
    /// `<path>.rs::<function>` citations of the tests that enforce the ruling.
    /// Never empty; see [`split_guard_citation`] for the grammar.
    Tests(Vec<String>),
    /// A stated, reasoned declaration that no test enforces this record.
    ///
    /// Carries the reason rather than a bare flag, for the reason
    /// `RETIRED_RECORD_IDS`'s second field is prose: "nothing enforces this"
    /// is only a usable answer when it says *why*, and the whys differ —
    /// an undecided record has no ruling to enforce, while a decided one may
    /// rule on how a departure is classified rather than on any output.
    Declined(String),
    /// The record carries no `guard` key at all — the silence the field exists
    /// to abolish.
    ///
    /// Constructed only from a missing or `null` key, never from anything an
    /// author could write, so it cannot be confused with a deliberate blank.
    /// Every reader below treats it as "names no test and declines nothing",
    /// which is the honest reading; the state is *reported* by
    /// [`ruling_guard_field::every_record_declares_a_guard`] rather than
    /// tolerated.
    Absent,
}

impl Guard {
    /// The cited tests, or an empty slice for a declared exemption or silence.
    pub fn tests(&self) -> &[String] {
        match self {
            Guard::Tests(tests) => tests,
            Guard::Declined(_) | Guard::Absent => &[],
        }
    }

    /// The stated reason, for a declared exemption only.
    ///
    /// [`Guard::Absent`] answers `None` — it states no reason, which is what is
    /// wrong with it. Use [`Guard::is_absent`] to ask about silence; asking this
    /// would conflate "declined nothing" with "declined without saying why".
    pub fn declined_reason(&self) -> Option<&str> {
        match self {
            Guard::Tests(_) | Guard::Absent => None,
            Guard::Declined(reason) => Some(reason),
        }
    }

    /// Whether the record said nothing at all about what enforces it.
    pub fn is_absent(&self) -> bool {
        matches!(self, Guard::Absent)
    }
}

/// The `README.md` rule a clause-free house choice is made under.
///
/// A closed set of two, mirroring `generate_spec_fixture`'s `overrides::HouseRule`.
/// Free text here would let a record write a spec clause into the slot reserved
/// for the project's own authority, which is the substitution [`HouseChoice`]
/// exists to prevent.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum HouseRule {
    /// `README.md` rule 5, silent limb — the recommendations are incomplete
    /// here, so ferro decides under rule 6 and violates nothing.
    RuleFiveSilent,
    /// `README.md` rule 6 — among multiple conformant forms, the maintainers
    /// choose.
    RuleSix,
}

impl HouseRule {
    /// The ledger's own spelling, which is also what rendered output uses.
    pub fn token(self) -> &'static str {
        match self {
            HouseRule::RuleFiveSilent => "rule-five-silent",
            HouseRule::RuleSix => "rule-six",
        }
    }

    /// Human-readable form, for rendered documents.
    pub fn label(self) -> &'static str {
        match self {
            HouseRule::RuleFiveSilent => "`README.md` rule 5 (silent limb)",
            HouseRule::RuleSix => "`README.md` rule 6",
        }
    }

    /// Parse the ledger's spelling. `None` for anything outside the closed set.
    pub fn from_token(token: &str) -> Option<Self> {
        match token {
            "rule-five-silent" => Some(HouseRule::RuleFiveSilent),
            "rule-six" => Some(HouseRule::RuleSix),
            _ => None,
        }
    }
}

/// A ruling that is the project's own, made where the recommendations are
/// silent.
///
/// See `generate_spec_fixture`'s `overrides::HouseChoice` for the full
/// reasoning; the short version is that the ledger previously had no way to say
/// "no clause reaches this and we chose X", so such choices lived as prose in
/// `src/` and `tests/` and drifted apart there.
///
/// Carrying it as a *record* rather than as a `status` value is deliberate: a
/// house choice is `decided` — a choice has been made — and the two published
/// tables (`docs/NORMALIZATION_CONTRACT.md`) partition the ledger
/// by status. A third status would have re-partitioned both and made every
/// house choice read as a species of open question, which is the opposite of
/// what it is.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HouseChoice {
    /// The `README.md` rule the choice is made under.
    pub under: HouseRule,
    /// What was weighed and put aside in reaching it. Never blank.
    pub considered_and_rejected: String,
}

/// Splits a guard citation into the file it names and the function within it.
///
/// The grammar is deliberately narrow — a repo-relative path ending in `.rs`,
/// then `::`, then a single snake_case identifier — so that a citation is
/// resolvable *mechanically* rather than by judgement. That is the whole
/// difference between this field and the prose it replaces: a bare
/// `tests/….rs` path names a file, not a proposition, and a module-qualified
/// `foo::bar` does not say which of several `foo.rs` it means.
///
/// Returns `None` for anything outside the grammar; callers report that as a
/// malformed citation rather than silently skipping it.
pub fn split_guard_citation(token: &str) -> Option<(&str, &str)> {
    let (path, function) = token.split_once("::")?;
    if !path.ends_with(".rs") || path.starts_with('/') || path.contains("..") {
        return None;
    }
    if function.is_empty() || function.contains("::") {
        return None;
    }
    if !function.starts_with(|c: char| c.is_ascii_lowercase()) {
        return None;
    }
    if !function
        .chars()
        .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_')
    {
        return None;
    }
    Some((path, function))
}

/// One adjudication record.
#[derive(Clone, Debug)]
pub struct Record {
    /// The record's `id`, the kebab-case token prose cites it by.
    pub id: String,
    /// `decided` or `undecided`, verbatim from the ledger. No other value
    /// parses; see [`STATUSES`].
    pub status: String,
    /// The record's `question` prose, verbatim — what the record is *about*,
    /// as opposed to how it came out.
    ///
    /// Required rather than optional. A record with no question is a verdict
    /// with nothing to be a verdict on, and the id cannot stand in for it: an
    /// id states the question in shorthand and, for at least one record in this
    /// ledger, states it as the negation of its own ruling. Exposed because
    /// `normalization_contract_doc.rs` renders it, and a rendering that fell
    /// back to the id would publish that negation as the question.
    pub question: String,
    /// The record's `rationale` prose, verbatim.
    ///
    /// Exposed because a rationale cites *other records* — that is how the
    /// ledger cross-references itself — and those citations are not covered by
    /// the source scan in `ruling_citation_currency.rs`, which reads `.rs` files
    /// only. See `every_record_to_record_citation_resolves` there.
    pub rationale: String,
    /// The descriptions the record declares itself to apply to, verbatim, in the
    /// ledger's own order. Absent or `null` reads as none.
    ///
    /// Exposed because it is the record's own statement of what it rules on, and
    /// so the strongest available evidence of which axes it reaches;
    /// `ledger_clause_jurisdiction.rs` reads the axis prefix of each entry.
    pub applies_to: Vec<String>,
    /// The curated equivalence classes whose convergence this record's ruling
    /// decides, verbatim, in the ledger's own order. Absent or `null` reads as
    /// none.
    ///
    /// Exposed because it is the record's own pointer at the pinned evidence
    /// that its ruling is enforced — `spec_equivalence_classes_converge` is what
    /// fails if the ruling stops holding — and a published rendering of a record
    /// that omitted it would drop the only link from the prose to the guard.
    pub equivalence_classes: Vec<String>,
    /// What enforces this record's ruling, or its stated reason that nothing
    /// does. Required — see [`Guard`].
    pub guard: Guard,
    /// Present when the ruling is the project's own rather than the spec's —
    /// see [`HouseChoice`]. `None` for every record that rules on a clause.
    ///
    /// Exposed rather than folded into [`Self::governing`]'s absence because the
    /// two are different states and conflating them is what the field exists to
    /// stop: a `decided` record with no governing clause used to be
    /// unrepresentable, and `normalization_contract_doc.rs` publishes the note
    /// "an open record names a conflict without choosing a side" for a missing
    /// governing clause — which would be exactly wrong on a house choice.
    pub house_choice: Option<HouseChoice>,
    /// Every clause the record names, in the ledger's own order, each tagged
    /// with the role the record's verdict fields give it.
    pub citations: Vec<Citation>,
    /// The clause the record rules **under**, verbatim, when it names one.
    ///
    /// Derivable from `citations` by looking for [`Role::Governing`], but
    /// exposed directly because "this record has no governing clause" and "this
    /// record's governing clause is molecule-agnostic" are different states and
    /// `ledger_clause_jurisdiction.rs` distinguishes them.
    pub governing: Option<String>,
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

/// An optional array-of-strings field.
///
/// Absent or `null` reads as empty; a present value of any other shape panics,
/// for the same reason [`optional_str`] does. Silently reading a malformed
/// `applies_to` as empty would narrow a record's declared reach to nothing and
/// pass every jurisdiction check built on it.
fn string_array(record: &serde_json::Value, field: &str, id: &str) -> Vec<String> {
    match record.get(field) {
        None | Some(serde_json::Value::Null) => Vec::new(),
        Some(value) => value
            .as_array()
            .unwrap_or_else(|| panic!("record {id} has a non-array `{field}`: {value}"))
            .iter()
            .map(|entry| {
                entry
                    .as_str()
                    .unwrap_or_else(|| {
                        panic!("record {id} has a non-string `{field}` entry: {entry}")
                    })
                    .to_string()
            })
            .collect(),
    }
}

/// The `guard` field of one record.
///
/// **Absence parses, as [`Guard::Absent`], and is refused one layer out** by
/// `ruling_guard_field::every_record_declares_a_guard`. Everything below that is
/// *malformed* rather than *absent* still panics here, and the split is
/// deliberate: absence is what a record in flight legitimately looks like before
/// its author has answered the question, so it must reach a check that can name
/// the fix; a `guard` that is present and self-contradictory is a broken record
/// and there is nothing useful to do with it but refuse.
///
/// The states this has to keep apart are "enforced by these tests", "enforced by
/// nothing, for this stated reason", and "said nothing". The third reads from
/// the outside exactly like the first, which is why it needs a name of its own
/// rather than being folded into an empty list.
///
/// Exactly one of `tests` / `none` may be set, for the same reason
/// `governing` and `deviates_from` may not name the same clause: a record
/// carrying both has said two incompatible things, and any reader that tested
/// one key first would resolve the contradiction silently in its favour.
fn guard_from(record: &serde_json::Value, id: &str) -> Guard {
    let value = match record.get("guard") {
        None | Some(serde_json::Value::Null) => return Guard::Absent,
        Some(value) => value,
    };
    let object = value
        .as_object()
        .unwrap_or_else(|| panic!("record {id} has a non-object `guard`: {value}"));
    for key in object.keys() {
        assert!(
            key == "tests" || key == "none",
            "record {id} has an unknown `guard` key {key:?}; only `tests` and `none` parse"
        );
    }

    let tests = object.get("tests");
    let none = object.get("none");
    match (tests, none) {
        (Some(_), Some(_)) => panic!(
            "record {id} sets both `guard.tests` and `guard.none` — it has both named an \
             enforcing test and declared that nothing enforces it"
        ),
        (None, None) => {
            panic!("record {id} has an empty `guard` object; set exactly one of `tests` or `none`")
        }
        (Some(tests), None) => {
            let entries = tests
                .as_array()
                .unwrap_or_else(|| panic!("record {id} has a non-array `guard.tests`: {tests}"));
            assert!(
                !entries.is_empty(),
                "record {id} has an empty `guard.tests`; an empty list claims enforcement and \
                 supplies none — declare `guard.none` with a reason instead"
            );
            let citations: Vec<String> = entries
                .iter()
                .map(|entry| {
                    let token = entry.as_str().unwrap_or_else(|| {
                        panic!("record {id} has a non-string `guard.tests` entry: {entry}")
                    });
                    assert!(
                        split_guard_citation(token).is_some(),
                        "record {id} cites {token:?} as a guard, which is not of the form \
                         `<path>.rs::<function>` — a citation that does not resolve \
                         mechanically is the prose this field replaces"
                    );
                    token.to_string()
                })
                .collect();
            Guard::Tests(citations)
        }
        (None, Some(none)) => {
            let reason = none
                .as_str()
                .unwrap_or_else(|| panic!("record {id} has a non-string `guard.none`: {none}"));
            // A blank reason is the absent case wearing a string: it satisfies
            // the key and states nothing, which is precisely the silence the
            // field exists to reject.
            assert!(
                !reason.trim().is_empty(),
                "record {id} declines a guard with a blank reason; `guard.none` is a first-class \
                 answer only when it says why"
            );
            Guard::Declined(reason.to_string())
        }
    }
}

/// The `house_choice` field of one record, if it carries one.
///
/// Absence is the ordinary case and reads as `None`; anything present is
/// validated in full, for the reason [`guard_from`] gives about a present
/// self-contradictory field. The two required keys are checked explicitly
/// rather than by a permissive `get`-and-default, because a house choice that
/// silently lost its `considered_and_rejected` would publish a conclusion with
/// its reasoning deleted and nothing would fail.
fn house_choice_from(record: &serde_json::Value, id: &str) -> Option<HouseChoice> {
    let value = match record.get("house_choice") {
        None | Some(serde_json::Value::Null) => return None,
        Some(value) => value,
    };
    let object = value
        .as_object()
        .unwrap_or_else(|| panic!("record {id} has a non-object `house_choice`: {value}"));
    for key in object.keys() {
        assert!(
            key == "under" || key == "considered_and_rejected",
            "record {id} has an unknown `house_choice` key {key:?}; only `under` and \
             `considered_and_rejected` parse"
        );
    }
    let under_token = required_str(
        value,
        "under",
        &format!("the `house_choice` of record {id}"),
    );
    let under = HouseRule::from_token(under_token).unwrap_or_else(|| {
        panic!(
            "record {id} makes its house choice under {under_token:?}, which is not one of \
             {:?}. `under` names the `README.md` rule the project chose under, never a spec \
             clause — a house choice that could name a clause here would be claiming the \
             authority it exists to disclaim",
            [
                HouseRule::RuleFiveSilent.token(),
                HouseRule::RuleSix.token()
            ]
        )
    });
    let considered_and_rejected = required_str(
        value,
        "considered_and_rejected",
        &format!("the `house_choice` of record {id}"),
    );
    assert!(
        !considered_and_rejected.trim().is_empty(),
        "record {id} has a blank `house_choice.considered_and_rejected`; a house choice with no \
         rejected alternative is a conclusion with its reasoning deleted"
    );
    Some(HouseChoice {
        under,
        considered_and_rejected: considered_and_rejected.to_string(),
    })
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

            let question = required_str(record, "question", &format!("record {id}")).to_string();
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
                .map(|entry| {
                    let subject = format!("a clause of record {id}");
                    let clause = required_str(entry, "clause", &subject).to_string();
                    let quote =
                        required_str(entry, "quote", &format!("clause {clause} of record {id}"))
                            .to_string();
                    // A blank quote is the absent case wearing a string: it
                    // checks against nothing in the spec checkout and supplies
                    // no axis text, which is precisely what requiring the field
                    // was meant to prevent.
                    assert!(
                        !quote.trim().is_empty(),
                        "clause {clause} of record {id} has a blank `quote`"
                    );
                    let role = if governing == Some(clause.as_str()) {
                        Role::Governing
                    } else if deviates_from.contains(&clause.as_str()) {
                        Role::DeviatesFrom
                    } else {
                        Role::Cited
                    };
                    Citation {
                        clause,
                        quote,
                        role,
                    }
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

            let applies_to: Vec<String> = string_array(record, "applies_to", &id);

            let equivalence_classes: Vec<String> = string_array(record, "equivalence_classes", &id);

            let guard = guard_from(record, &id);

            // A house choice states that no clause reaches the point. The two
            // verdict fields state that one does. Refused here as well as in the
            // generator because this reader is what the committed tests and both
            // published documents are built from: a record that reached them
            // carrying both would be rendered as a house choice *and* indexed as
            // a governing authority for the clause it names.
            let house_choice = house_choice_from(record, &id);
            if house_choice.is_some() {
                assert!(
                    status == "decided",
                    "record {id} carries a `house_choice` but is {status:?}; a choice that has \
                     been made is decided by definition"
                );
                if let Some(governing) = governing {
                    panic!(
                        "record {id} is a `house_choice` and names {governing:?} as governing. A \
                         house choice is what the project does where the recommendations do NOT \
                         decide, so it cannot also hold a clause to govern"
                    );
                }
                assert!(
                    deviates_from.is_empty(),
                    "record {id} is a `house_choice` and names {deviates_from:?} as deviated \
                     from; there is nothing to deviate from where the recommendations are silent"
                );
            } else {
                assert!(
                    status != "decided" || governing.is_some(),
                    "record {id} is decided but names neither a governing clause nor a \
                     `house_choice`"
                );
            }

            Record {
                id,
                status,
                question,
                rationale,
                applies_to,
                equivalence_classes,
                guard,
                citations,
                governing: governing.map(str::to_string),
                house_choice,
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
// (see the repository `CONTRIBUTING.md` on committed tests that have never executed).
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
            "question": "Does a well-formed record parse?",
            "rationale": "A test record, with prose that cites no other record.",
            "governing": "docs/a.md:1",
            "guard": { "tests": ["tests/it/a_file.rs::a_guard"] },
            "clauses": [{ "clause": "docs/a.md:1", "quote": "a quoted clause" }],
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

/// A record whose ruling is the project's own, for the house-choice cases.
///
/// Built by *removing* `governing` from the valid record rather than by writing
/// a second literal, so the two fixtures cannot drift into disagreeing about
/// everything except the field under test.
fn one_house_choice_record() -> serde_json::Value {
    let mut document = one_valid_record();
    let record = &mut document["rulings"][0];
    record
        .as_object_mut()
        .expect("the fixture record is an object")
        .remove("governing");
    record["house_choice"] = serde_json::json!({
        "under": "rule-six",
        "considered_and_rejected": "the alternative form, rejected for a stated reason",
    });
    document
}

/// The positive control for the house-choice arm, and its discriminator.
///
/// The second half is what gives the first any force: with the `house_choice`
/// removed the very same record is refused, so "decided with no governing
/// clause" is legal *because of* the declaration and not by accident.
#[test]
fn a_house_choice_record_parses_and_is_the_only_clause_free_decided_shape() {
    let records = records_from_value(&one_house_choice_record(), "<test>");
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].status, "decided");
    assert!(records[0].governing.is_none());
    let choice = records[0]
        .house_choice
        .as_ref()
        .expect("the record declares a house choice");
    assert_eq!(choice.under, HouseRule::RuleSix);
    // Every clause it names is context, never an authority.
    assert!(records[0].citations.iter().all(|c| c.role == Role::Cited));
}

#[test]
#[should_panic(expected = "names neither a governing clause nor a `house_choice`")]
fn a_decided_record_with_no_governing_clause_and_no_house_choice_is_rejected() {
    let mut document = one_house_choice_record();
    document["rulings"][0]
        .as_object_mut()
        .expect("record")
        .remove("house_choice");
    records_from_value(&document, "<test>");
}

/// The refusal the field exists for: a house choice may not claim that a clause
/// governs. Overclaiming spec authority is what turned one implementer's choice
/// into six prose restatements presented as compliance.
#[test]
#[should_panic(expected = "names \"docs/a.md:1\" as governing")]
fn a_house_choice_that_claims_a_governing_clause_is_rejected() {
    let mut document = one_house_choice_record();
    document["rulings"][0]["governing"] = serde_json::json!("docs/a.md:1");
    records_from_value(&document, "<test>");
}

/// The sibling refusal: a clause that reaches far enough to be departed from is
/// a clause that reaches.
#[test]
#[should_panic(expected = "as deviated from")]
fn a_house_choice_that_deviates_from_a_clause_is_rejected() {
    let mut document = one_house_choice_record();
    document["rulings"][0]["deviates_from"] = serde_json::json!(["docs/a.md:1"]);
    records_from_value(&document, "<test>");
}

#[test]
#[should_panic(expected = "decided by definition")]
fn an_undecided_house_choice_is_rejected() {
    let mut document = one_house_choice_record();
    document["rulings"][0]["status"] = serde_json::json!("undecided");
    document["rulings"][0]["clauses"] = serde_json::json!([
        { "clause": "docs/a.md:1", "quote": "a quoted clause" },
        { "clause": "docs/b.md:2", "quote": "another quoted clause" },
    ]);
    records_from_value(&document, "<test>");
}

/// `under` names a `README.md` rule and nothing else — in particular not a spec
/// clause, which is the value a drifting record would most plausibly reach for.
#[test]
#[should_panic(expected = "which is not one of")]
fn a_house_choice_under_a_spec_clause_is_rejected() {
    let mut document = one_house_choice_record();
    document["rulings"][0]["house_choice"]["under"] = serde_json::json!("docs/a.md:1");
    records_from_value(&document, "<test>");
}

#[test]
#[should_panic(expected = "blank `house_choice.considered_and_rejected`")]
fn a_house_choice_with_nothing_rejected_is_rejected() {
    let mut document = one_house_choice_record();
    document["rulings"][0]["house_choice"]["considered_and_rejected"] = serde_json::json!("  ");
    records_from_value(&document, "<test>");
}

#[test]
#[should_panic(expected = "unknown `house_choice` key")]
fn an_unknown_house_choice_key_is_rejected() {
    let mut document = one_house_choice_record();
    document["rulings"][0]["house_choice"]["governing"] = serde_json::json!("docs/a.md:1");
    records_from_value(&document, "<test>");
}

/// Absence reads as absence, so every record predating the field is unaffected.
#[test]
fn an_absent_house_choice_reads_as_none() {
    let records = records_from_value(&one_valid_record(), "<test>");
    assert!(records[0].house_choice.is_none());
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

/// A missing `question` is rejected rather than read as empty prose.
///
/// The id is not a fallback for it. At least one record in this ledger is
/// decided the *opposite* way to what its id reads as, so a renderer that
/// substituted the id for an absent question would publish that record's
/// question — and by implication its answer — backwards, and it would look like
/// prose rather than like a missing field. (The record is named in
/// `clause_ruling_index.rs`; it is deliberately not named here, per this
/// module's header.)
#[test]
#[should_panic(expected = "has no `question`")]
fn a_missing_question_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]
        .as_object_mut()
        .expect("record is an object")
        .remove("question");
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
        { "clause": "docs/a.md:1", "quote": "a quoted clause" },
        { "clause": "docs/b.md:2", "quote": "another quoted clause" },
    ]);
    document["rulings"][0]["deviates_from"] = serde_json::json!(["docs/b.md:2", 17]);
    records_from_value(&document, "<test>");
}

/// A clause with no `quote` is rejected rather than read as an empty quote.
///
/// An empty quote is not a harmless absence: it is a `file:line` the spec-fixture
/// generator has nothing to check against the spec checkout, and it supplies no
/// axis text for `ledger_clause_jurisdiction.rs` to read. Both would go on
/// passing.
#[test]
#[should_panic(expected = "has no `quote`")]
fn a_clause_with_no_quote_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["clauses"] = serde_json::json!([{ "clause": "docs/a.md:1" }]);
    records_from_value(&document, "<test>");
}

/// A present-but-blank `quote` is rejected too.
///
/// Requiring the key buys nothing on its own: `""` satisfies it, checks against
/// no line of the spec checkout, and supplies `ledger_clause_jurisdiction.rs`
/// with no axis text — the same silent narrowing the required field exists to
/// stop, one step further in.
#[test]
#[should_panic(expected = "has a blank `quote`")]
fn a_clause_with_a_blank_quote_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["clauses"] =
        serde_json::json!([{ "clause": "docs/a.md:1", "quote": "   " }]);
    records_from_value(&document, "<test>");
}

/// A non-array `applies_to` is rejected rather than read as absent — the same
/// reasoning as `deviates_from` above, and it matters more here because
/// `ledger_clause_jurisdiction.rs` reads that array to decide which axes a
/// record rules on. Silently reading it as empty would narrow the record's
/// declared reach to nothing and pass every jurisdiction check.
#[test]
#[should_panic(expected = "non-array `applies_to`")]
fn a_non_array_applies_to_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["applies_to"] = serde_json::json!("c.235_237delinsTAT");
    records_from_value(&document, "<test>");
}

/// An absent `applies_to` reads as none, which is what most records spell.
#[test]
fn an_absent_applies_to_reads_as_none() {
    let records = records_from_value(&one_valid_record(), "<test>");
    assert!(records[0].applies_to.is_empty());
    assert_eq!(records[0].governing.as_deref(), Some("docs/a.md:1"));
}

/// An explicit `null` verdict field still means "not set" — that is a
/// legitimate JSON spelling of absent, and the ledger's own records simply omit
/// the keys.
///
/// Driven on an `undecided` record because that is the status for which "both
/// verdict fields absent" is the *correct* shape. It used to be driven on the
/// `decided` fixture, which now trips the rule that a decided record names
/// either a governing clause or a `house_choice` — a rule that did not exist
/// when this test was written, and whose absence is what made a clause-free
/// decided ruling unrepresentable in the first place.
#[test]
fn a_null_verdict_field_reads_as_absent() {
    let mut document = one_valid_record();
    document["rulings"][0]["status"] = serde_json::json!("undecided");
    document["rulings"][0]["governing"] = serde_json::Value::Null;
    document["rulings"][0]["deviates_from"] = serde_json::Value::Null;
    let records = records_from_value(&document, "<test>");
    assert_eq!(records[0].citations[0].role, Role::Cited);
}

// --------------------------------------------------------------------------
// `guard`.
//
// The field's whole purpose is that ONE of its states is unrepresentable, so
// the mutations below are the specification. Each one used to be expressible
// — as an absent field, an empty list, a blank string — and each reads from
// the outside like a record whose ruling something enforces.
// --------------------------------------------------------------------------

/// A record with no `guard` **parses**, as [`Guard::Absent`].
///
/// It is refused by `ruling_guard_field::every_record_declares_a_guard`, not
/// here. Rejecting it at parse took the whole ledger down with it — 50 tests on
/// the measurement recorded on [`Guard`] — for a fault in one record, so the
/// reader's job is to *represent* the silence faithfully and let the check that
/// is about guards be the one that fails.
#[test]
fn a_record_with_no_guard_parses_as_absent() {
    let mut document = one_valid_record();
    document["rulings"][0]
        .as_object_mut()
        .expect("record is an object")
        .remove("guard");
    let records = records_from_value(&document, "<test>");
    assert_eq!(records[0].guard, Guard::Absent);
    assert!(records[0].guard.is_absent());
}

/// An explicit `null` guard is the same silence spelled differently, and reads
/// the same way.
#[test]
fn a_null_guard_parses_as_absent() {
    let mut document = one_valid_record();
    document["rulings"][0]["guard"] = serde_json::Value::Null;
    let records = records_from_value(&document, "<test>");
    assert_eq!(records[0].guard, Guard::Absent);
}

/// Silence is a state of its own, distinct from every blank an author could
/// deliberately write.
///
/// This is the property the whole staging change rests on. If absence collapsed
/// onto an empty `tests` list or a blank `none`, the check one layer out could
/// not tell "unconsidered" from "considered and mis-stated" — and those need
/// different messages, because they are different mistakes.
///
/// The two deliberate blanks are refused at parse and stay refused; that half is
/// pinned by [`an_empty_guard_test_list_is_rejected`] and
/// [`a_blank_exemption_reason_is_rejected`]. What is asserted here is the half
/// those cannot state: that the value absence *does* produce is equal to neither
/// of them, so no reader can conflate the three.
#[test]
fn absence_is_a_state_of_its_own() {
    let mut document = one_valid_record();
    document["rulings"][0]
        .as_object_mut()
        .expect("record is an object")
        .remove("guard");
    let absent = records_from_value(&document, "<test>")[0].guard.clone();

    assert_eq!(absent, Guard::Absent);
    assert_ne!(
        absent,
        Guard::Tests(Vec::new()),
        "silence must not read as an empty citation list"
    );
    assert_ne!(
        absent,
        Guard::Declined(String::new()),
        "silence must not read as an exemption with a blank reason"
    );
    // And it answers the accessors the way the checks below read them: it names
    // no test, and it declines nothing.
    assert!(absent.tests().is_empty());
    assert_eq!(absent.declined_reason(), None);
}

/// A declared exemption parses, and carries its reason.
#[test]
fn a_declared_exemption_parses() {
    let mut document = one_valid_record();
    document["rulings"][0]["guard"] =
        serde_json::json!({ "none": "process record; it mandates no output" });
    let records = records_from_value(&document, "<test>");
    assert_eq!(
        records[0].guard,
        Guard::Declined("process record; it mandates no output".to_string())
    );
    assert!(records[0].guard.tests().is_empty());
}

/// A blank reason is the absent case wearing a string: it satisfies the key and
/// states nothing.
#[test]
#[should_panic(expected = "blank reason")]
fn a_blank_exemption_reason_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["guard"] = serde_json::json!({ "none": "   " });
    records_from_value(&document, "<test>");
}

/// An empty `tests` array claims enforcement and supplies none — the same
/// silence, one level in.
#[test]
#[should_panic(expected = "empty `guard.tests`")]
fn an_empty_guard_test_list_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["guard"] = serde_json::json!({ "tests": [] });
    records_from_value(&document, "<test>");
}

/// Setting both keys is a contradiction, and must not be resolved silently in
/// favour of whichever is tested first.
#[test]
#[should_panic(expected = "both `guard.tests` and `guard.none`")]
fn a_guard_that_both_cites_and_declines_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["guard"] =
        serde_json::json!({ "tests": ["tests/it/a.rs::b"], "none": "and also nothing" });
    records_from_value(&document, "<test>");
}

/// An unknown key is rejected rather than ignored, so a typo (`test`, `tests_`)
/// cannot read as a record that declared nothing.
#[test]
#[should_panic(expected = "unknown `guard` key")]
fn an_unknown_guard_key_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["guard"] = serde_json::json!({ "test": ["tests/it/a.rs::b"] });
    records_from_value(&document, "<test>");
}

/// A citation outside the grammar is rejected at parse, not skipped by the
/// resolver — a skipped citation is a guard nobody checks.
#[test]
#[should_panic(expected = "which is not of the form")]
fn a_guard_citation_outside_the_grammar_is_rejected() {
    let mut document = one_valid_record();
    document["rulings"][0]["guard"] = serde_json::json!({ "tests": ["tests/it/a_file.rs"] });
    records_from_value(&document, "<test>");
}

/// The grammar, both ways. Too loose and a bare module path (`merge::foo`)
/// passes while naming no file; too tight and a real `src/` unit test stops
/// being citable.
#[test]
fn the_guard_citation_grammar_admits_a_path_and_nothing_else() {
    assert_eq!(
        split_guard_citation("tests/it/foo.rs::a_guard_name"),
        Some(("tests/it/foo.rs", "a_guard_name"))
    );
    assert_eq!(
        split_guard_citation("src/normalize/merge.rs::a_unit_test_2"),
        Some(("src/normalize/merge.rs", "a_unit_test_2"))
    );
    for rejected in [
        "merge::a_guard_name",           // a module path names no file
        "tests/it/foo.rs",               // a file names no proposition
        "a_guard_name",                  // a bare function name
        "tests/it/foo.rs::mod::a_guard", // more than one `::`
        "tests/it/foo.rs::AGuard",       // not a function name
        "tests/it/foo.rs::",             // nothing after the separator
        "/abs/foo.rs::a_guard",          // not repo-relative
        "tests/../../foo.rs::a_guard",   // escapes the tree
    ] {
        assert_eq!(
            split_guard_citation(rejected),
            None,
            "{rejected} must not parse as a guard citation"
        );
    }
}
