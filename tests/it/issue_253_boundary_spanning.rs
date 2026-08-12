//! Audit for #81 § J1 — user-supplied ranges spanning the 5'UTR↔CDS
//! or CDS↔3'UTR boundary (and the `n.` upstream↔transcript and
//! transcript↔downstream analogs).
//!
//! Policy: accept on input, preserve verbatim on Display, no warning.
//! The HGVS v21 spec does not explicitly forbid these forms. The
//! merge-barrier behavior from #103 / #89 is orthogonal — it only
//! prevents adjacent same-coord-system edits in opposite regions from
//! being folded into a single delins, never rejects user-supplied
//! boundary-spanning ranges.

use ferro_hgvs::parse_hgvs;

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
}

// =============================================================================
// SECTION 1 — c. 5'UTR↔CDS boundary
// =============================================================================

mod cds_5utr_boundary {
    use super::*;

    #[test]
    fn neg_to_pos_del_round_trips() {
        assert_round_trips("NM_000088.3:c.-3_5del");
    }

    #[test]
    fn neg_to_pos_del_with_ref_round_trips() {
        assert_round_trips("NM_000088.3:c.-3_5delAGCCC");
    }

    #[test]
    fn neg_to_pos_delins_round_trips() {
        assert_round_trips("NM_000088.3:c.-3_5delinsAGCCC");
    }

    #[test]
    fn neg_to_pos_dup_round_trips() {
        assert_round_trips("NM_000088.3:c.-3_5dup");
    }

    #[test]
    fn neg_to_pos_inv_round_trips() {
        assert_round_trips("NM_000088.3:c.-3_5inv");
    }

    /// `c.-1_1insAAA` — insertion across the boundary (no nucleotides
    /// removed, just an insertion at the UTR↔CDS junction).
    #[test]
    fn boundary_insertion_round_trips() {
        assert_round_trips("NM_000088.3:c.-1_1insAAA");
    }
}

// =============================================================================
// SECTION 2 — c. CDS↔3'UTR boundary
// =============================================================================

mod cds_3utr_boundary {
    use super::*;

    #[test]
    fn pos_to_star_del_round_trips() {
        assert_round_trips("NM_000088.3:c.4548_*3del");
    }

    #[test]
    fn pos_to_star_delins_round_trips() {
        assert_round_trips("NM_000088.3:c.4548_*3delinsAGC");
    }

    #[test]
    fn star_to_star_insertion_round_trips() {
        assert_round_trips("NM_000088.3:c.*1_*2insAAA");
    }
}

// =============================================================================
// SECTION 3 — c. ranges spanning both boundaries (entire CDS)
// =============================================================================

mod entire_cds_spanning {
    use super::*;

    /// `c.-3_*3del` — spans 5'UTR → CDS → 3'UTR. Accepted by the spec
    /// for completeness (no rule forbids the form).
    #[test]
    fn full_span_del_round_trips() {
        assert_round_trips("NM_000088.3:c.-3_*3del");
    }
}

// =============================================================================
// SECTION 4 — n. upstream↔transcript and transcript↔downstream
// =============================================================================

mod nx_boundaries {
    use super::*;

    #[test]
    fn n_upstream_to_body_del_round_trips() {
        assert_round_trips("NR_037639.1:n.-3_5del");
    }

    /// The 3' half of this boundary was pinned here as a round-trip
    /// (`NR_037639.1:n.40_*3del`) until #1748 refused `n.*N` at parse in every
    /// mode. `background/numbering.md:52` numbers this axis from the first to
    /// the last nucleotide of the reference sequence and `:53` grants intronic
    /// offsets as its only other zone, so there is no `*` region for a range to
    /// span *into* — the boundary this row was auditing does not exist on `n.`.
    ///
    /// Kept as a refusal rather than deleted, because that is the audit finding.
    /// The `c.` equivalent (`c.-3_*3del`, §3) is unaffected: the coding `*`
    /// zone is anchored to the stop codon and is still inside the transcript.
    /// The 5' half above is also unaffected — `n.-N` is refused in strict mode
    /// only, and this entry applies no mode (#1632).
    #[test]
    fn n_body_to_downstream_del_is_refused() {
        let err = parse_hgvs("NR_037639.1:n.40_*3del")
            .map(|v| v.to_string())
            .expect_err("`n.*N` is refused at parse in every mode (#1748)")
            .to_string();
        assert!(err.contains("E1003"), "expected the position code: {err}");
        assert!(
            err.contains("numbering.md:52"),
            "expected the clause: {err}"
        );
    }
}

// =============================================================================
// SECTION 5 — Same-region positive controls
// =============================================================================
//
// Pin that the same-region ranges keep round-tripping; the
// boundary-spanning audit cases above must not regress them.

mod same_region_positive_controls {
    use super::*;

    #[test]
    fn five_utr_only_del_round_trips() {
        assert_round_trips("NM_000088.3:c.-5_-3del");
    }

    #[test]
    fn cds_only_del_round_trips() {
        assert_round_trips("NM_000088.3:c.3_5del");
    }

    #[test]
    fn three_utr_only_del_round_trips() {
        assert_round_trips("NM_000088.3:c.*1_*5del");
    }
}
