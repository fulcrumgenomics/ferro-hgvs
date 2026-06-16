//! Issue #394 item 3 — SPDI returns `UnsupportedEditType` (not
//! `MissingReferenceData`) for non-literal delins inserted sequences.
//!
//! `src/spdi/convert.rs:1000-1004` previously returned
//! `ConversionError::MissingReferenceData` for non-literal delins
//! inserted sequences (`InsertedSequence::Reference`,
//! `PositionRange`, `Empty`). Semantically that's "reference data
//! missing"; the actual category is "this edit shape cannot be encoded
//! as a single SPDI variant". Sibling arms for `Repeat` /
//! `Inversion` / `Conversion` already use `UnsupportedEditType` for
//! the same shape-not-encodable category — this file pins the
//! consistency.

use ferro_hgvs::spdi::{hgvs_to_spdi, ConversionError};
use ferro_hgvs::{parse_hgvs, MockProvider};

#[test]
fn spdi_non_literal_delins_inserted_sequence_returns_unsupported_edit_type() {
    // `g.100_105delins[NC_000002.12:g.500_510]` — a delins whose
    // inserted sequence is a cross-reference range. SPDI carries no
    // multi-reference encoding; the conversion must fail with
    // `UnsupportedEditType`, NOT `MissingReferenceData`.
    let v = parse_hgvs("NC_000001.11:g.100_105delins[NC_000002.12:g.500_510]").unwrap();
    let provider = MockProvider::new();
    let result = hgvs_to_spdi(&v, &provider);
    let err = result.expect_err("non-literal delins must error");
    match err {
        ConversionError::UnsupportedEditType { description } => {
            assert!(
                description.to_lowercase().contains("delins")
                    || description.to_lowercase().contains("literal"),
                "error description must reference delins or literal sequence; got: {description}",
            );
        }
        ConversionError::MissingReferenceData { description } => {
            panic!(
                "regressed: non-literal delins must surface UnsupportedEditType, not \
                 MissingReferenceData (which means 'reference fetch failed'). \
                 Description: {description}",
            );
        }
        other => panic!("unexpected error variant: {other:?}"),
    }
}
