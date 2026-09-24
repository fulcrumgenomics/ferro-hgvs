//! The ruling ledger's pointer at the **Normalization rules** page,
//! `docs/src/reference/normalization-rules.md`, must resolve.
//!
//! The `adjudication-precedence-order` record does not restate the rules; it names
//! the page that states them. A pointer can dangle, and nothing else checks this
//! one: `ruling_citation_currency.rs` checks the ledger's citations of *spec
//! clauses*, not of this repository's own documents.
//!
//! The record must also stay a pointer: "do not restate the rules" is the
//! substance of its ruling. The rule names it must not repeat are read from the
//! page, so this module keeps no copy of them.
//!
//! The page's own wording is deliberately not pinned. It is the single source of
//! the rules, so a test that required it to keep particular text would be a
//! second copy of it.

use std::path::PathBuf;

use crate::common::ruleset_page::{restated_rule_openers, RULESET_PAGE};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

#[test]
fn the_ledgers_pointer_at_the_ruleset_resolves() {
    let ledger_path =
        repo_root().join("tests/fixtures/grammar/hgvs_spec_normalization_overrides.json");
    let text = std::fs::read_to_string(&ledger_path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", ledger_path.display()));
    let ledger: serde_json::Value =
        serde_json::from_str(&text).expect("the overrides ledger is valid JSON");

    let record = ledger["rulings"]
        .as_array()
        .expect("`rulings` is an array")
        .iter()
        .find(|r| r["id"] == "adjudication-precedence-order")
        .expect(
            "the `adjudication-precedence-order` record must exist; it is what names the \
             ruleset page",
        );
    let rationale = record["rationale"].as_str().expect("a string rationale");

    assert!(
        rationale.contains(RULESET_PAGE),
        "`adjudication-precedence-order` no longer names `{RULESET_PAGE}`. If the ruleset \
         page moved, update the record and this constant together"
    );
    assert!(
        repo_root().join(RULESET_PAGE).is_file(),
        "`adjudication-precedence-order` points at `{RULESET_PAGE}`, which does not exist"
    );

    let restated = restated_rule_openers(rationale);
    assert!(
        restated.is_empty(),
        "`adjudication-precedence-order` restates the rules {restated:?} it is supposed to point \
         at. The ruling is that the rules are stated in exactly one place"
    );
}
