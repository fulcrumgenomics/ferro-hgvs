//! Audit for #81 § J2 — `n.` upstream (`-N`) and downstream (`*N`)
//! markers, canonical output.
//!
//! Per HGVS v21, the `n.` coordinate system uses:
//! - `-N` upstream of transcript start
//! - `N` inside the transcript
//! - `*N` downstream of transcript end
//! - `N+M` / `N-M` intronic offsets relative to position `N`
//!
//! All edit types should round-trip with each marker. No
//! canonicalization or marker-stripping is expected.

use ferro_hgvs::parse_hgvs;

const TX: &str = "NR_037639.1";

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
}

// =============================================================================
// SECTION 1 — Upstream `-N` marker
// =============================================================================

mod upstream {
    use super::*;

    #[test]
    fn neg_single_sub_round_trips() {
        assert_round_trips(&format!("{TX}:n.-1A>G"));
    }

    #[test]
    fn neg_deeper_sub_round_trips() {
        assert_round_trips(&format!("{TX}:n.-5A>G"));
    }

    #[test]
    fn neg_range_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.-100_-50del"));
    }

    #[test]
    fn neg_short_range_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.-3_-1del"));
    }

    #[test]
    fn neg_range_delins_round_trips() {
        assert_round_trips(&format!("{TX}:n.-3_-1delinsACG"));
    }

    #[test]
    fn neg_range_dup_round_trips() {
        assert_round_trips(&format!("{TX}:n.-3_-1dup"));
    }

    #[test]
    fn neg_range_inv_round_trips() {
        assert_round_trips(&format!("{TX}:n.-3_-1inv"));
    }
}

// =============================================================================
// SECTION 2 — Downstream `*N` marker
// =============================================================================

mod downstream {
    use super::*;

    #[test]
    fn star_single_sub_round_trips() {
        assert_round_trips(&format!("{TX}:n.*1A>G"));
    }

    #[test]
    fn star_deeper_sub_round_trips() {
        assert_round_trips(&format!("{TX}:n.*100A>G"));
    }

    #[test]
    fn star_range_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.*1_*5del"));
    }

    #[test]
    fn star_range_delins_round_trips() {
        assert_round_trips(&format!("{TX}:n.*3_*5delinsACG"));
    }

    #[test]
    fn star_range_dup_round_trips() {
        assert_round_trips(&format!("{TX}:n.*3_*5dup"));
    }

    #[test]
    fn star_range_inv_round_trips() {
        assert_round_trips(&format!("{TX}:n.*3_*5inv"));
    }
}

// =============================================================================
// SECTION 3 — Interior positions (positive controls)
// =============================================================================

mod interior {
    use super::*;

    #[test]
    fn interior_sub_round_trips() {
        assert_round_trips(&format!("{TX}:n.100A>G"));
    }

    #[test]
    fn interior_range_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.1_100del"));
    }
}

// =============================================================================
// SECTION 4 — Intronic offsets `N+M` / `N-M`
// =============================================================================

mod intronic_offsets {
    use super::*;

    #[test]
    fn plus_offset_sub_round_trips() {
        assert_round_trips(&format!("{TX}:n.100+1A>G"));
    }

    #[test]
    fn minus_offset_sub_round_trips() {
        assert_round_trips(&format!("{TX}:n.100-1A>G"));
    }

    #[test]
    fn plus_offset_range_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.100+5_100+10del"));
    }

    #[test]
    fn minus_offset_range_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.100-5_100-1del"));
    }
}

// =============================================================================
// SECTION 5 — Insertions and boundary-crossing ranges
// =============================================================================

mod insertions_and_crossings {
    use super::*;

    #[test]
    fn upstream_to_interior_ins_round_trips() {
        assert_round_trips(&format!("{TX}:n.-1_1insAAA"));
    }

    #[test]
    fn interior_to_downstream_ins_round_trips() {
        assert_round_trips(&format!("{TX}:n.100_*1insAAA"));
    }

    /// Boundary-spanning ranges are also covered by issue #253 (J1);
    /// pinned here to keep the J2 marker matrix self-contained.
    #[test]
    fn upstream_to_interior_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.-3_5del"));
    }

    #[test]
    fn interior_to_downstream_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.40_*3del"));
    }
}
