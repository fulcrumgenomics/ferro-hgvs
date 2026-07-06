//! Issue #951 — `o.` (circular, SVD-WG006) normalization is a deliberate
//! pass-through: the variant is returned unchanged, with no 3'-shift and no
//! warning. A genuine circular normalizer would 3'-shift with origin-wraparound
//! semantics, which is not yet implemented (tracked by the open #951; the prior
//! in-code reference to the now-closed, MT-scoped #129 was misleading).
//!
//! These tests pin the current documented behavior so a future implementer of
//! circular shuffling knows exactly what contract they are changing.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// A `del` is a classic 3'-shift candidate on the linear axes, yet the `o.`
/// normalizer returns it verbatim with no shuffle applied and no warning.
#[test]
fn circular_deletion_is_returned_unchanged_without_warnings() {
    let normalizer = Normalizer::new(MockProvider::new());
    let input = "NC_001416.1:o.5del";
    let variant = parse_hgvs(input).expect("parse o. deletion");
    let outcome = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("normalize o. deletion");
    assert_eq!(
        format!("{}", outcome.result),
        input,
        "circular o. normalization must be a pass-through (no 3'-shift) — #951",
    );
    assert!(
        outcome.warnings.is_empty(),
        "circular o. pass-through must emit no warnings; got {:?}",
        outcome.warnings,
    );
}

/// Idempotence follows trivially from the pass-through, but pin it too so a
/// future shuffling implementation cannot silently make the second pass differ.
#[test]
fn circular_normalization_is_idempotent() {
    let normalizer = Normalizer::new(MockProvider::new());
    let variant = parse_hgvs("NC_001416.1:o.100_105del").expect("parse");
    let once = normalizer.normalize(&variant).expect("first normalize");
    let twice = normalizer.normalize(&once).expect("second normalize");
    assert_eq!(
        format!("{once}"),
        format!("{twice}"),
        "o. normalize idempotent"
    );
}
