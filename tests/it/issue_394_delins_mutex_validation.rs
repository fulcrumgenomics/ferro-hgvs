//! Issue #394 item 2 — validation rejects `Delins { deleted: Some(_),
//! deleted_length: Some(_) }`.
//!
//! `src/hgvs/edit.rs:472-474` docstring states: "At most one of
//! `deleted` / `deleted_length` is `Some`". Previously this was a
//! documentation-only invariant; the parser doesn't emit the
//! dual-Some shape, but direct struct construction could. The
//! `Display` impl silently prefers `deleted` and drops
//! `deleted_length`, masking the inconsistency.
//!
//! This file promotes the invariant to runtime validation via a new
//! `E004` error in `validate_na_edit`.

use std::str::FromStr;

use ferro_hgvs::hgvs::edit::{InsertedSequence, NaEdit, Sequence};
use ferro_hgvs::hgvs::validation::rules::validate_na_edit;

#[test]
fn validation_rejects_delins_with_both_deleted_and_deleted_length() {
    let edit = NaEdit::Delins {
        sequence: InsertedSequence::Literal(Sequence::from_str("G").unwrap()),
        deleted: Some(Sequence::from_str("ATGC").unwrap()),
        deleted_length: Some(4),
        substitution_reference: None,
    };
    let result = validate_na_edit(&edit);
    assert!(
        !result.valid,
        "validation must reject Delins with both `deleted` and `deleted_length` set; \
         got valid result: {:?}",
        result,
    );
    let has_mutex_error = result
        .errors
        .iter()
        .any(|e| e.message.to_lowercase().contains("deleted"));
    assert!(
        has_mutex_error,
        "validation error must mention `deleted` / `deleted_length` mutual exclusion; \
         got errors: {:?}",
        result.errors,
    );
}

#[test]
fn validation_accepts_delins_with_deleted_only() {
    let edit = NaEdit::Delins {
        sequence: InsertedSequence::Literal(Sequence::from_str("G").unwrap()),
        deleted: Some(Sequence::from_str("ATGC").unwrap()),
        deleted_length: None,
        substitution_reference: None,
    };
    let result = validate_na_edit(&edit);
    assert!(
        result.valid,
        "validation must accept Delins with only `deleted` set; got errors: {:?}",
        result.errors,
    );
}

#[test]
fn validation_accepts_delins_with_deleted_length_only() {
    let edit = NaEdit::Delins {
        sequence: InsertedSequence::Literal(Sequence::from_str("G").unwrap()),
        deleted: None,
        deleted_length: Some(4),
        substitution_reference: None,
    };
    let result = validate_na_edit(&edit);
    assert!(
        result.valid,
        "validation must accept Delins with only `deleted_length` set; got errors: {:?}",
        result.errors,
    );
}

#[test]
fn validation_accepts_canonical_short_form_delins() {
    let edit = NaEdit::Delins {
        sequence: InsertedSequence::Literal(Sequence::from_str("G").unwrap()),
        deleted: None,
        deleted_length: None,
        substitution_reference: None,
    };
    let result = validate_na_edit(&edit);
    assert!(
        result.valid,
        "validation must accept canonical short-form Delins; got errors: {:?}",
        result.errors,
    );
}
