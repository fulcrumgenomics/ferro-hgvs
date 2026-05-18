//! Audit for #288 — whole-entity predicted forms across the
//! mitochondrial (`m.`), circular DNA (`o.`), and non-coding transcript
//! (`n.`) coordinate systems. Follow-up to #245 / PR #246, which wired
//! the same shapes for `c.` / `g.` / `r.`.
//!
//! Per HGVS v21 `recommendations/uncertain.md`, the predicted-form
//! wrapper extends to whole-entity edits on every coord system:
//!
//! - `m.(=)`, `o.(=)`, `n.(=)` — predicted whole-entity identity.
//! - `m.(?)`, `o.(?)`, `n.(?)` — predicted whole-entity unknown.
//! - `n.(0)` — predicted "no product" for the non-coding transcript.
//!   `m.0` / `o.0` are not spec forms (DNA isn't "produced"), so
//!   `m.(0)` / `o.(0)` stay rejected, matching the `c.(0)` guard
//!   pinned in #245.
//!
//! The bare forms (`m.=`, `m.?`, `o.=`, `o.?`, `n.=`, `n.?`, `n.0`)
//! also had no parser dispatch before this issue and are wired in the
//! same change. The certain/uncertain pair gets pinned together to
//! make Display round-trips meaningful.
//!
//! Accessions:
//!   - `NC_012920.1` — Homo sapiens mitochondrion (m.)
//!   - `NC_001416.1` — Enterobacteria phage lambda (o., circular DNA)
//!   - `NR_046018.2` — Human DDX11L1 non-coding RNA (n.)

use ferro_hgvs::parse_hgvs;

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
}

#[track_caller]
fn assert_rejected(s: &str) {
    assert!(
        parse_hgvs(s).is_err(),
        "{s:?} unexpectedly parsed — should be rejected"
    );
}

const MT: &str = "NC_012920.1";
const CIRC: &str = "NC_001416.1";
const TX: &str = "NR_046018.2";

// =============================================================================
// SECTION 1 — Predicted whole-entity identity `(=)`
// =============================================================================

mod predicted_identity {
    use super::*;

    #[test]
    fn mt_predicted_identity_round_trips() {
        assert_round_trips(&format!("{MT}:m.(=)"));
    }

    #[test]
    fn circular_predicted_identity_round_trips() {
        assert_round_trips(&format!("{CIRC}:o.(=)"));
    }

    #[test]
    fn tx_predicted_identity_round_trips() {
        assert_round_trips(&format!("{TX}:n.(=)"));
    }
}

// =============================================================================
// SECTION 2 — Predicted whole-entity unknown `(?)`
// =============================================================================

mod predicted_unknown {
    use super::*;

    #[test]
    fn mt_predicted_unknown_round_trips() {
        assert_round_trips(&format!("{MT}:m.(?)"));
    }

    #[test]
    fn circular_predicted_unknown_round_trips() {
        assert_round_trips(&format!("{CIRC}:o.(?)"));
    }

    #[test]
    fn tx_predicted_unknown_round_trips() {
        assert_round_trips(&format!("{TX}:n.(?)"));
    }
}

// =============================================================================
// SECTION 3 — Predicted "no product" `(0)` — n. only
// =============================================================================

mod predicted_no_product {
    use super::*;

    /// `n.0` (no transcript produced) exists for non-coding transcripts;
    /// `n.(0)` is its predicted form. Mirrors `r.(0)` from #245.
    #[test]
    fn tx_predicted_no_product_round_trips() {
        assert_round_trips(&format!("{TX}:n.(0)"));
    }

    /// `m.0` and `o.0` aren't spec forms — DNA isn't "produced" — so
    /// their predicted wrappers stay rejected. Regression guard
    /// mirroring `c.(0)` rejection from #245.
    #[test]
    fn mt_no_product_rejected() {
        assert_rejected(&format!("{MT}:m.0"));
        assert_rejected(&format!("{MT}:m.(0)"));
    }

    #[test]
    fn circular_no_product_rejected() {
        assert_rejected(&format!("{CIRC}:o.0"));
        assert_rejected(&format!("{CIRC}:o.(0)"));
    }
}

// =============================================================================
// SECTION 4 — Bare (certain) forms — pinned alongside the predicted
// pair, since the bare forms had no parser dispatch before this issue
// either.
// =============================================================================

mod certain_forms {
    use super::*;

    #[test]
    fn mt_certain_identity_round_trips() {
        assert_round_trips(&format!("{MT}:m.="));
    }

    #[test]
    fn mt_certain_unknown_round_trips() {
        assert_round_trips(&format!("{MT}:m.?"));
    }

    #[test]
    fn circular_certain_identity_round_trips() {
        assert_round_trips(&format!("{CIRC}:o.="));
    }

    #[test]
    fn circular_certain_unknown_round_trips() {
        assert_round_trips(&format!("{CIRC}:o.?"));
    }

    #[test]
    fn tx_certain_identity_round_trips() {
        assert_round_trips(&format!("{TX}:n.="));
    }

    #[test]
    fn tx_certain_unknown_round_trips() {
        assert_round_trips(&format!("{TX}:n.?"));
    }

    #[test]
    fn tx_certain_no_product_round_trips() {
        assert_round_trips(&format!("{TX}:n.0"));
    }
}

// =============================================================================
// SECTION 5 — Boundary regression guards: c./g./r. forms from PR #246
// still work, positional forms on m./o./n. still work.
// =============================================================================

mod regression_guards {
    use super::*;

    /// Pinned in #245 — must keep working on top of #288.
    #[test]
    fn cds_genome_rna_predicted_forms_still_round_trip() {
        assert_round_trips("NM_000088.3:c.(=)");
        assert_round_trips("NM_000088.3:c.(?)");
        assert_round_trips("NC_000017.11:g.(=)");
        assert_round_trips("NC_000017.11:g.(?)");
        assert_round_trips("NM_000088.3:r.(=)");
        assert_round_trips("NM_000088.3:r.(?)");
        assert_round_trips("NM_000088.3:r.(0)");
    }

    /// Sanity: positional m./o./n. forms unaffected.
    #[test]
    fn positional_forms_still_round_trip() {
        assert_round_trips(&format!("{MT}:m.100A>G"));
        assert_round_trips(&format!("{CIRC}:o.100A>G"));
        assert_round_trips(&format!("{TX}:n.100A>G"));
    }
}
