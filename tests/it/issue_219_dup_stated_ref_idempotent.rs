//! Audit for issue #219 — `c.<pos>dup<base>` with a mismatched
//! stated-ref base is not idempotent across two normalize passes.
//!
//! Surfaced during the #218 (C4 mixed-accession) audit. Reproducer:
//! `c.10dupA` against a transcript where position 10 is `C` (start of
//! a `CCCCC` homopolymer) normalizes to `c.10dup` on pass 1 (the
//! stated-ref base is stripped because it doesn't match the
//! reference, but the position is not 3'-shifted), then shifts to
//! `c.14dup` on pass 2 (the position-only `c.10dup` form correctly
//! 3'-shifts through the homopolymer).
//!
//! The normalizer invariant `normalize(normalize(v)) == normalize(v)`
//! requires a single pass to produce the fixed point. After the fix,
//! pass 1 produces `c.14dup` directly.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Build a `MockProvider` with `NM_TEST.1` containing a 60bp sequence
/// laid out as: `ATGC | AAAAA | CCCCC | GGGGG | TTTTT | AAAAA | CCCCC
/// | GGGGG | TTTTT | AAAAA | CCCCC | GGGGG | T`. A dup at any
/// position inside the `CCCCC` run starting at c.10 should 3'-shift to
/// the rightmost copy (c.14).
fn provider_with_test_transcript() -> MockProvider {
    let mut p = MockProvider::new();
    let seq = "ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT".to_string();
    let len = seq.len() as u64;
    let tx = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        seq,
        Some(1),
        Some(len),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    p.add_transcript(tx);
    p
}

fn normalize_pass(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let v = parse_hgvs(input).expect("parse");
    let n = normalizer.normalize(&v).expect("normalize");
    format!("{}", n)
}

// =============================================================================
// SECTION 1 — solo `c.<pos>dup<base>` idempotency (the core bug)
// =============================================================================

mod solo_dup_stated_ref {
    use super::*;

    /// Single pass on `c.10dupA` produces the 3'-shifted canonical
    /// form `c.14dup`, NOT the intermediate `c.10dup`.
    #[test]
    fn dup_with_mismatched_stated_ref_shifts_in_one_pass() {
        let out = normalize_pass(provider_with_test_transcript(), "NM_TEST.1:c.10dupA");
        assert_eq!(out, "NM_TEST.1:c.14dup");
    }

    /// Two-pass idempotency: pass 2 over pass 1's output is a fixed
    /// point. This is the universal invariant the bug was breaking.
    #[test]
    fn dup_with_mismatched_stated_ref_is_idempotent() {
        let normalizer = Normalizer::new(provider_with_test_transcript());
        let parsed = parse_hgvs("NM_TEST.1:c.10dupA").expect("parse");
        let p1 = normalizer.normalize(&parsed).expect("normalize pass 1");
        let p2 = normalizer.normalize(&p1).expect("normalize pass 2");
        assert_eq!(format!("{}", p1), format!("{}", p2));
    }

    /// Companion case: `c.10dup` (no stated-ref) already shifts in one
    /// pass today. Pin as regression guard so a fix doesn't accidentally
    /// regress the working path.
    #[test]
    fn dup_position_only_already_shifts() {
        let out = normalize_pass(provider_with_test_transcript(), "NM_TEST.1:c.10dup");
        assert_eq!(out, "NM_TEST.1:c.14dup");
    }

    /// `c.10dupC` (matching stated-ref): existing behavior — shifts in
    /// one pass. Pin as regression guard.
    #[test]
    fn dup_with_matching_stated_ref_shifts_in_one_pass() {
        let out = normalize_pass(provider_with_test_transcript(), "NM_TEST.1:c.10dupC");
        assert_eq!(out, "NM_TEST.1:c.14dup");
    }
}

// =============================================================================
// SECTION 2 — inside cis / trans compounds (the original surfacing context)
// =============================================================================
//
// The bug was originally surfaced inside a mixed-accession compound
// (`[NM_TEST.1:c.10dupA;NM_OTHER.1:c.5A>G]`). Pin both the cis and
// trans wrappers across both same-coord-different-accession and
// cross-coord-system partners.

mod inside_compound {
    use super::*;

    #[test]
    fn cis_same_coord_different_accession_idempotent() {
        let normalizer = Normalizer::new(provider_with_test_transcript());
        let parsed = parse_hgvs("[NM_TEST.1:c.10dupA;NM_OTHER.1:c.5A>G]").expect("parse");
        let p1 = normalizer.normalize(&parsed).expect("normalize pass 1");
        let p2 = normalizer.normalize(&p1).expect("normalize pass 2");
        assert_eq!(format!("{}", p1), format!("{}", p2));
        // Also pin the canonical output shape.
        assert_eq!(format!("{}", p1), "[NM_TEST.1:c.14dup;NM_OTHER.1:c.5A>G]");
    }

    #[test]
    fn cis_cross_coord_system_idempotent() {
        let normalizer = Normalizer::new(provider_with_test_transcript());
        let parsed = parse_hgvs("[NM_TEST.1:c.10dupA;NC_000001.11:g.100A>G]").expect("parse");
        let p1 = normalizer.normalize(&parsed).expect("normalize pass 1");
        let p2 = normalizer.normalize(&p1).expect("normalize pass 2");
        assert_eq!(format!("{}", p1), format!("{}", p2));
        assert_eq!(
            format!("{}", p1),
            "[NM_TEST.1:c.14dup;NC_000001.11:g.100A>G]"
        );
    }

    #[test]
    fn trans_same_coord_different_accession_idempotent() {
        let normalizer = Normalizer::new(provider_with_test_transcript());
        let parsed = parse_hgvs("[NM_TEST.1:c.10dupA];[NM_OTHER.1:c.5A>G]").expect("parse");
        let p1 = normalizer.normalize(&parsed).expect("normalize pass 1");
        let p2 = normalizer.normalize(&p1).expect("normalize pass 2");
        assert_eq!(format!("{}", p1), format!("{}", p2));
    }
}

// =============================================================================
// SECTION 3 — symmetric `del<base>` mismatched-stated-ref pin
// =============================================================================
//
// The same root cause (shuffle uses stated bases when present) could
// theoretically affect deletions with mismatched stated-ref. The
// existing `validate_reference` flags the mismatch with a warning,
// and the shuffle on Deletion in fact does not read alt bases (alt
// is empty for del). So del should always be idempotent regardless
// of stated-ref. Pin as a guard.

mod del_stated_ref {
    use super::*;

    #[test]
    fn del_with_mismatched_stated_ref_is_idempotent() {
        let normalizer = Normalizer::new(provider_with_test_transcript());
        // pos 10 is C; stating A is wrong. del is del regardless.
        let parsed = parse_hgvs("NM_TEST.1:c.10delA").expect("parse");
        let p1 = normalizer.normalize(&parsed).expect("normalize pass 1");
        let p2 = normalizer.normalize(&p1).expect("normalize pass 2");
        assert_eq!(format!("{}", p1), format!("{}", p2));
    }
}
