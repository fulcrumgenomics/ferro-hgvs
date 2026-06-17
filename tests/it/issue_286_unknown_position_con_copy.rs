//! Audit for #286 — `?con<src>` / `?copy<N>` at unknown position.
//!
//! PR #240 (closes #239) accepted `?<edit>` on g./c./r./n. for
//! substitutions, identity, and the keyword-keyed edits (`del`,
//! `dup`, `ins`, `inv`). Conversion (`con<src>`) and copy-number
//! (`copy<N>`) edits paired with an unknown position were explicitly
//! deferred because the spec is silent on whether they are permitted
//! at `?`.
//!
//! Policy decision (issue #286): **accept**, on the precedent that
//! ferro already accepts every other `?<edit>` shape. An unknown
//! position with a defined edit semantics is consistent with the
//! "edit at unknown position" reading already in place for
//! substitution / identity / del / dup / ins / inv.
//!
//! This audit pins:
//! - `?con<src>` on c., g., n., r. — round-trip parity.
//! - `?copy<N>` on c., g., n., r. — round-trip parity.
//! - Negative shapes that must still error (e.g. `?con` with no
//!   source — a malformed conversion regardless of position).

use ferro_hgvs::parse_hgvs;

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
    // Idempotent: re-parse the rendered form and confirm Display is
    // stable across a second round.
    let v2 = parse_hgvs(&out).unwrap_or_else(|e| panic!("re-parse {out:?} failed: {e}"));
    let out2 = format!("{}", v2);
    assert_eq!(out2, out, "second round-trip mismatch for {s:?}");
}

#[track_caller]
fn assert_rejects(s: &str) {
    let result = parse_hgvs(s);
    assert!(
        result.is_err(),
        "expected parse failure for {s:?}, got {:?}",
        result.ok().map(|v| format!("{v}"))
    );
}

const CDS: &str = "NM_000088.3";
const GENOME: &str = "NC_000017.11";
const TX: &str = "NR_037639.1";
const MITO: &str = "NC_012920.1";
const CIRCULAR: &str = "NC_001802.1";

// =============================================================================
// SECTION 1 — `?con<src>` at unknown position
// =============================================================================
//
// Conversion takes a source HGVS expression (accession + colon +
// coord-system-specific interval). The intent of `?con<src>` is "a
// conversion happened from `<src>` into this transcript / genome /
// RNA, but the destination position is unknown". This is consistent
// with PR #240's framing for `?del` / `?dup` / `?A>G`.

mod conversion_at_unknown_position {
    use super::*;

    #[test]
    fn cds_unknown_position_conversion_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?conNM_000089.1:c.789_1011"));
    }

    #[test]
    fn genome_unknown_position_conversion_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.?conNC_000022.10:g.17179029_17179091"));
    }

    #[test]
    fn rna_unknown_position_conversion_round_trips() {
        assert_round_trips(&format!("{CDS}:r.?conNM_000089.1:r.789_1011"));
    }

    #[test]
    fn tx_unknown_position_conversion_round_trips() {
        assert_round_trips(&format!("{TX}:n.?conNR_000002.1:n.30_40"));
    }

    // m. and o. flow through `parse_uncertain_genome_pos` and already
    // accepted `?con<src>` / `?copy<N>` prior to #286. The tests below
    // pin that this behavior remains consistent with the c./g./n./r.
    // routes added by #286 — i.e. the four uncertain-position-bearing
    // coordinate systems all expose the same shape.

    #[test]
    fn mito_unknown_position_conversion_round_trips() {
        assert_round_trips(&format!("{MITO}:m.?conNC_012920.1:m.100_200"));
    }

    #[test]
    fn circular_unknown_position_conversion_round_trips() {
        assert_round_trips(&format!("{CIRCULAR}:o.?conNC_001802.1:o.100_200"));
    }
}

// =============================================================================
// SECTION 2 — `?copy<N>` at unknown position
// =============================================================================
//
// Copy-number takes an integer count (e.g. `copy2`, `copy4`). The
// intent of `?copy<N>` is "a copy-number variant exists but the
// affected interval is unknown".

mod copy_number_at_unknown_position {
    use super::*;

    #[test]
    fn cds_unknown_position_copy2_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?copy2"));
    }

    #[test]
    fn cds_unknown_position_copy_large_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?copy42"));
    }

    #[test]
    fn genome_unknown_position_copy3_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.?copy3"));
    }

    #[test]
    fn rna_unknown_position_copy4_round_trips() {
        assert_round_trips(&format!("{CDS}:r.?copy4"));
    }

    #[test]
    fn tx_unknown_position_copy2_round_trips() {
        assert_round_trips(&format!("{TX}:n.?copy2"));
    }

    #[test]
    fn mito_unknown_position_copy2_round_trips() {
        assert_round_trips(&format!("{MITO}:m.?copy2"));
    }

    #[test]
    fn circular_unknown_position_copy2_round_trips() {
        assert_round_trips(&format!("{CIRCULAR}:o.?copy2"));
    }
}

// =============================================================================
// SECTION 3 — Negative shapes
// =============================================================================
//
// `?con` with no source remains a malformed conversion. `?copy` with
// no integer remains a malformed copy-number edit. These are edit-level
// errors that exist independent of the position, and routing `?` to
// the position-based parser must not paper over them.

mod negative_shapes {
    use super::*;

    /// `?con` with no source after the keyword — malformed.
    #[test]
    fn cds_unknown_position_conversion_no_source_rejects() {
        assert_rejects(&format!("{CDS}:c.?con"));
    }

    /// `?copy` with no integer after the keyword — malformed.
    #[test]
    fn cds_unknown_position_copy_no_count_rejects() {
        assert_rejects(&format!("{CDS}:c.?copy"));
    }

    /// `?copy` with non-numeric trailing — malformed.
    #[test]
    fn genome_unknown_position_copy_non_numeric_rejects() {
        assert_rejects(&format!("{GENOME}:g.?copyABC"));
    }

    // The c.-only negatives above leave g./r./n. uncovered. The next
    // block pins that `?con` with no source and `?copy` with no count
    // are rejected at every coordinate system that #286 added, not just
    // the c. route.

    /// `?con` with no source on g. — malformed.
    #[test]
    fn genome_unknown_position_conversion_no_source_rejects() {
        assert_rejects(&format!("{GENOME}:g.?con"));
    }

    /// `?copy` with no count on g. — malformed.
    #[test]
    fn genome_unknown_position_copy_no_count_rejects() {
        assert_rejects(&format!("{GENOME}:g.?copy"));
    }

    /// `?con` with no source on r. — malformed.
    #[test]
    fn rna_unknown_position_conversion_no_source_rejects() {
        assert_rejects(&format!("{CDS}:r.?con"));
    }

    /// `?copy` with no count on r. — malformed.
    #[test]
    fn rna_unknown_position_copy_no_count_rejects() {
        assert_rejects(&format!("{CDS}:r.?copy"));
    }

    /// `?con` with no source on n. — malformed.
    #[test]
    fn tx_unknown_position_conversion_no_source_rejects() {
        assert_rejects(&format!("{TX}:n.?con"));
    }

    /// `?copy` with no count on n. — malformed.
    #[test]
    fn tx_unknown_position_copy_no_count_rejects() {
        assert_rejects(&format!("{TX}:n.?copy"));
    }
}

// =============================================================================
// SECTION 4 — Backwards-compatibility regression guards
// =============================================================================
//
// The known-position forms must continue to parse and round-trip
// exactly as before. Issue #286's change only adds new routes for
// `?` + (con | copy); existing position-bearing forms are untouched.

mod regression_guards {
    use super::*;

    #[test]
    fn cds_known_position_conversion_still_round_trips() {
        assert_round_trips(&format!("{CDS}:c.123_456conNM_000089.1:c.789_1011"));
    }

    #[test]
    fn cds_known_position_copy_still_round_trips() {
        assert_round_trips(&format!("{CDS}:c.123_456copy2"));
    }

    #[test]
    fn genome_known_position_conversion_still_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.12345_12400conNC_000022.10:g.100_155"));
    }

    #[test]
    fn genome_known_position_copy_still_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.12345_12400copy4"));
    }

    /// Whole-entity unknown (`?` alone) must still parse to the
    /// whole-entity-unknown edit, not be re-routed to the position
    /// parser by the new con/copy logic.
    #[test]
    fn cds_whole_entity_unknown_still_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?"));
    }

    #[test]
    fn genome_whole_entity_unknown_still_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.?"));
    }
}
