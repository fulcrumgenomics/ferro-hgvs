//! Audit for #245 — whole-entity predicted forms across NA coord
//! systems. Follow-up to #241 / PR #242 and #243 / PR #244.
//!
//! Per HGVS v21 `recommendations/uncertain.md`, the predicted-form
//! wrapper extends to whole-entity edits:
//!
//! - `c.(=)`, `g.(=)`, `r.(=)` — predicted whole-entity identity.
//! - `c.(?)`, `g.(?)`, `r.(?)` — predicted whole-entity unknown.
//! - `r.(0)` — predicted RNA no-product (RNA-only; no `c.(0)` form).
//!
//! Protein `p.(=)` and `p.(?)` already work; `p.(0)` is a separate
//! protein-level bug (tracked elsewhere).

use ferro_hgvs::parse_hgvs;

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
}

const CDS: &str = "NM_000088.3";
const GENOME: &str = "NC_000017.11";

// =============================================================================
// SECTION 1 — Predicted whole-entity identity `(=)`
// =============================================================================

mod predicted_identity {
    use super::*;

    #[test]
    fn cds_predicted_identity_round_trips() {
        assert_round_trips(&format!("{CDS}:c.(=)"));
    }

    #[test]
    fn genome_predicted_identity_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.(=)"));
    }

    #[test]
    fn rna_predicted_identity_round_trips() {
        assert_round_trips(&format!("{CDS}:r.(=)"));
    }
}

// =============================================================================
// SECTION 2 — Predicted whole-entity unknown `(?)`
// =============================================================================

mod predicted_unknown {
    use super::*;

    #[test]
    fn cds_predicted_unknown_round_trips() {
        assert_round_trips(&format!("{CDS}:c.(?)"));
    }

    #[test]
    fn genome_predicted_unknown_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.(?)"));
    }

    #[test]
    fn rna_predicted_unknown_round_trips() {
        assert_round_trips(&format!("{CDS}:r.(?)"));
    }
}

// =============================================================================
// SECTION 3 — Predicted RNA no-product `r.(0)`
// =============================================================================

mod predicted_no_product {
    use super::*;

    /// `r.0` exists (no transcript produced); `r.(0)` is its predicted form.
    #[test]
    fn rna_predicted_no_product_round_trips() {
        assert_round_trips(&format!("{CDS}:r.(0)"));
    }

    /// `c.0` is not a spec form — `0` is RNA-specific. The predicted
    /// `c.(0)` should likewise be rejected. Pinned as regression guard.
    #[test]
    fn cds_no_product_rejected() {
        assert!(parse_hgvs(&format!("{CDS}:c.0")).is_err());
        assert!(parse_hgvs(&format!("{CDS}:c.(0)")).is_err());
    }
}

// =============================================================================
// SECTION 4 — Certain forms remain working (regression guards)
// =============================================================================

mod regression_guards {
    use super::*;

    #[test]
    fn certain_identity_forms_still_round_trip() {
        assert_round_trips(&format!("{CDS}:c.="));
        assert_round_trips(&format!("{GENOME}:g.="));
        assert_round_trips(&format!("{CDS}:r.="));
    }

    #[test]
    fn certain_unknown_forms_still_round_trip() {
        assert_round_trips(&format!("{CDS}:c.?"));
        assert_round_trips(&format!("{GENOME}:g.?"));
        assert_round_trips(&format!("{CDS}:r.?"));
    }

    #[test]
    fn certain_rna_no_product_still_round_trips() {
        assert_round_trips(&format!("{CDS}:r.0"));
    }
}
