//! Issue #1182 — `ferro project` emitted unparseable `n.0`, extrapolated
//! off-chromosome genomic coordinates, and discarded W4004.
//!
//! Three separable defects, one per section below:
//!
//!   1. A `c.` position 5' of the transcript's first base mapped to a
//!      non-positive transcript coordinate and was *rendered* — `n.0del`, which
//!      ferro's own parser rejects, and `n.-62del`, which parses cleanly and so
//!      propagated silently.
//!   2. `VariantProjection::normalization_warnings` was populated by the
//!      projector but read by nobody: `AxisOutcome` carried no warnings field
//!      and `output_projection` had nowhere to put one, so W4004 was
//!      structurally unreachable from `ferro project` in any configuration.
//!   3. `project` had no `--error-mode` flag at all, so even once the warning
//!      was reachable there was no way to ask for it to be a rejection.
//!
//! Only the reporting-surface half lives here. The two in-crate halves are
//! covered where their fixtures already exist, because neither type can be built
//! from a separate crate: `noncoding_from_coding`'s refusal in
//! `src/project/transcript_axis.rs` (`Transcript` has private fields) and
//! `select_axis`'s warning plumbing in `src/cli/project.rs`
//! (`VariantProjection` is `#[non_exhaustive]`).

use ferro_hgvs::cli::project::ProjectionWarning;
use ferro_hgvs::normalize::NormalizationWarning;

// =============================================================================
// 2. Warnings must survive the AxisOutcome hop
// =============================================================================

/// `ProjectionWarning` must preserve both halves of a normalizer warning: the
/// stable code (for machine consumers reading `--format json`) and the rendered
/// message (for humans).
#[test]
fn a_normalization_warning_flattens_to_code_and_message() {
    let warning = NormalizationWarning::RefSeqMismatch {
        stated_ref: "A".to_string(),
        actual_ref: "C".to_string(),
        position: "10-10".to_string(),
        corrected: true,
        details: None,
    };
    let flat = ProjectionWarning::from_normalization(&warning);
    assert_eq!(flat.code, warning.code());
    assert_eq!(flat.message, warning.to_string());
    assert!(
        !flat.message.is_empty(),
        "the rendered message must not be empty"
    );
}

/// Order is preserved, because the normalizer emits warnings in a meaningful
/// sequence (the expansion warning first, then downstream ones).
#[test]
fn a_warning_set_flattens_in_order() {
    let warnings = vec![
        NormalizationWarning::RefSeqMismatch {
            stated_ref: "A".to_string(),
            actual_ref: "C".to_string(),
            position: "10-10".to_string(),
            corrected: true,
            details: None,
        },
        NormalizationWarning::RefSeqMismatch {
            stated_ref: "G".to_string(),
            actual_ref: "T".to_string(),
            position: "20-20".to_string(),
            corrected: true,
            details: None,
        },
    ];
    let flat = ProjectionWarning::from_normalization_set(&warnings);
    assert_eq!(flat.len(), 2);
    assert!(flat[0].message.contains("10-10"), "got {:?}", flat[0]);
    assert!(flat[1].message.contains("20-20"), "got {:?}", flat[1]);
}

/// An empty warning set stays empty rather than becoming a placeholder — a
/// clean projection must not grow a spurious diagnostic.
#[test]
fn a_clean_projection_carries_no_warnings() {
    assert!(ProjectionWarning::from_normalization_set(&[]).is_empty());
}
