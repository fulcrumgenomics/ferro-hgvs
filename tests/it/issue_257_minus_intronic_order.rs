//! Audit for #81 § J3 — minus-strand intronic position ordering at the
//! parser/Display layer.
//!
//! The normalization-side rule (restore coding order after
//! genomic-space normalization, the `swapped` restore in
//! `normalize_intronic_{cds,tx,boundary_spanning_cds}`) is pinned by
//! `tests/coverage_gap_tests.rs` § minus-strand modules. This audit
//! covers the orthogonal parse + Display surface: position pairs in
//! ascending coding (5'→3') order must round-trip unchanged.

use ferro_hgvs::parse_hgvs;

const CDS: &str = "NM_000088.3";
const TX: &str = "NR_037639.1";

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
}

// =============================================================================
// SECTION 1 — Intronic insertion preserves ascending coding order
// =============================================================================

mod intronic_insertion_order {
    use super::*;

    /// Plus-offset pair: `30+2` < `30+3` on the coding axis.
    #[test]
    fn plus_offset_insertion_round_trips() {
        assert_round_trips(&format!("{CDS}:c.30+2_30+3insA"));
    }

    #[test]
    fn plus_offset_long_insertion_round_trips() {
        assert_round_trips(&format!("{CDS}:c.30+3_30+4insGGG"));
    }

    /// Minus-offset pair: `31-3` < `31-2` on the coding axis.
    #[test]
    fn minus_offset_insertion_round_trips() {
        assert_round_trips(&format!("{CDS}:c.31-3_31-2insT"));
    }

    /// n. analog of the same rule.
    #[test]
    fn nx_plus_offset_insertion_round_trips() {
        assert_round_trips(&format!("{TX}:n.30+2_30+3insA"));
    }

    #[test]
    fn nx_minus_offset_insertion_round_trips() {
        assert_round_trips(&format!("{TX}:n.31-3_31-2insT"));
    }
}

// =============================================================================
// SECTION 2 — Other edit types preserve coding order
// =============================================================================

mod intronic_other_edits {
    use super::*;

    #[test]
    fn plus_offset_del_round_trips() {
        assert_round_trips(&format!("{CDS}:c.30+1_30+5del"));
    }

    #[test]
    fn plus_offset_dup_round_trips() {
        assert_round_trips(&format!("{CDS}:c.30+1_30+5dup"));
    }

    #[test]
    fn plus_offset_inv_round_trips() {
        assert_round_trips(&format!("{CDS}:c.30+1_30+5inv"));
    }

    #[test]
    fn plus_offset_delins_round_trips() {
        assert_round_trips(&format!("{CDS}:c.30+1_30+5delinsACG"));
    }

    #[test]
    fn minus_offset_del_round_trips() {
        assert_round_trips(&format!("{CDS}:c.31-5_31-1del"));
    }

    #[test]
    fn minus_offset_dup_round_trips() {
        assert_round_trips(&format!("{CDS}:c.31-5_31-1dup"));
    }

    #[test]
    fn minus_offset_inv_round_trips() {
        assert_round_trips(&format!("{CDS}:c.31-5_31-1inv"));
    }

    /// n. analog spans.
    #[test]
    fn nx_plus_offset_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.30+1_30+5del"));
    }

    #[test]
    fn nx_minus_offset_dup_round_trips() {
        assert_round_trips(&format!("{TX}:n.31-5_31-1dup"));
    }
}

// =============================================================================
// SECTION 3 — Cross-intron range
// =============================================================================
//
// A range from `c.30+5` to `c.31-5` straddles the intron's interior in
// ascending coding order. Pin the round-trip.

mod cross_intron_range {
    use super::*;

    #[test]
    fn cross_intron_del_round_trips() {
        assert_round_trips(&format!("{CDS}:c.30+5_31-5del"));
    }
}

// =============================================================================
// SECTION 4 — Single intronic substitution (no pair to swap)
// =============================================================================
//
// Pin as a regression guard: single-position intronic subs never had
// the ordering bug, but they should keep round-tripping.

mod single_intronic_sub {
    use super::*;

    #[test]
    fn plus_offset_sub_round_trips() {
        assert_round_trips(&format!("{CDS}:c.30+1A>G"));
    }

    #[test]
    fn minus_offset_sub_round_trips() {
        assert_round_trips(&format!("{CDS}:c.31-1A>G"));
    }
}
