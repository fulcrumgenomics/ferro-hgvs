//! Issue #486 — overlapping cis-allele edits that involve an insertion
//! (mutalyzer `EOVERLAP`) must reject under strict normalize as
//! `FerroError::InvalidCoordinates`, and warn-and-preserve under lenient.
//!
//! The pre-existing coincident-bounds detector (`overlap.rs`) excluded
//! insertions; this pins the extension covering the two insertion-overlap
//! shapes the mutalyzer corpus exercises:
//!   - two insertions at the same interspace (`[4_5insT;4_5insA]`); and
//!   - an insertion whose junction is strictly interior to a span edit
//!     (`[274_275delinsT;274_275insA]`).
//!
//! Reference-free (MockProvider) so it holds independent of the manifest,
//! mirroring `issue_395_overlap_conflict_strict_rejection.rs`.
use ferro_hgvs::error::FerroError;
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

fn provider() -> MockProvider {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_000001.11", "A".repeat(200));
    p
}

#[test]
fn strict_rejects_two_insertions_at_same_junction() {
    let normalizer = Normalizer::with_config(provider(), NormalizeConfig::strict());
    let v = parse_hgvs("NC_000001.11:g.[100_101insT;100_101insA]").expect("parse");
    let err = normalizer
        .normalize(&v)
        .expect_err("strict mode must reject two insertions at one junction");
    match err {
        FerroError::InvalidCoordinates { msg } => assert!(
            msg.contains("W5002") || msg.contains("OverlapConflictingEdits"),
            "expected W5002 / OverlapConflictingEdits; got: {msg}"
        ),
        other => panic!("unexpected error variant: {other:?}"),
    }
}

#[test]
fn strict_rejects_insertion_interior_to_delins() {
    let normalizer = Normalizer::with_config(provider(), NormalizeConfig::strict());
    let v = parse_hgvs("NC_000001.11:g.[100_103delinsTT;101_102insG]").expect("parse");
    let err = normalizer
        .normalize(&v)
        .expect_err("strict mode must reject an insertion interior to a delins span");
    match err {
        FerroError::InvalidCoordinates { msg } => assert!(
            msg.contains("W5002") || msg.contains("OverlapConflictingEdits"),
            "expected W5002 / OverlapConflictingEdits; got: {msg}"
        ),
        other => panic!("unexpected error variant: {other:?}"),
    }
}

#[test]
fn lenient_warns_and_preserves_overlapping_insertions() {
    // Default (lenient) mode must NOT reject — it preserves the input and
    // attaches the W5002 warning (matches the issue_395 contract).
    let normalizer = Normalizer::new(provider());
    let v = parse_hgvs("NC_000001.11:g.[100_101insT;100_101insA]").expect("parse");
    let result = normalizer
        .normalize_with_diagnostics(&v)
        .expect("lenient mode must accept");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "OVERLAP_CONFLICTING_EDITS"),
        "expected OVERLAP_CONFLICTING_EDITS warning, got {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>()
    );
}

#[test]
fn strict_accepts_non_overlapping_flanking_insertions() {
    // Insertions either side of a single-base substitution do not overlap it
    // (the junction abuts but is not interior). Must NOT reject — the
    // spec-valid `[99_100ins;100A>C;100_101ins]` shape.
    let normalizer = Normalizer::with_config(provider(), NormalizeConfig::strict());
    let v = parse_hgvs("NC_000001.11:g.[99_100insT;100A>C;100_101insG]").expect("parse");
    assert!(
        normalizer.normalize(&v).is_ok(),
        "non-overlapping flanking insertions must not reject"
    );
}
