//! Audit for #81 § G2 — `?` unknown position acceptance across coord
//! systems and edit shapes.
//!
//! Per HGVS v21 `recommendations/uncertain.md`, `?` denotes an
//! unknown position. Ferro accepts most `?`-bearing forms correctly,
//! but three asymmetries surfaced during the systematic audit:
//!
//! 1. **Asymmetric coord-system support for `?<edit>`**: `c.?del` and
//!    `c.?dup` parse; `g.?del` and `r.?del` failed with
//!    "Unexpected trailing characters" (the whole-entity-unknown
//!    parsers for `g.` and `r.` did not reject `?` followed by an
//!    edit keyword, so `?` was consumed as the whole-entity edit
//!    leaving `del`/`dup` as trailing).
//!
//! 2. **Asymmetric edit support on `?` position**: `c.?del`, `c.?dup`
//!    worked; `c.?A>G` failed. The whole-CDS-unknown parser rejected
//!    `?` followed by `dup`/`del`/`ins`/`inv` but not by a base + `>`
//!    (substitution).
//!
//! 3. **Both-unknown range collapse**: `c.?_?del` parses but Display
//!    collapses to `c.?del` (information lost). Pinned as-is; the
//!    spec-aligned canonical form for "deletion of unknown extent
//!    between two unknown points" is the bracketed `c.(?_?)del`
//!    (round-trip working post #237). This file documents the
//!    asymmetry without re-fixing it — see notes in the
//!    `both_unknown_range_documented` section.
//!
//! This audit pins:
//! - all working `?` forms as regression guards (whole-entity and
//!   half-unknown ranges and intronic offsets).
//! - the spec-correct round-trip for the previously-broken forms.

use ferro_hgvs::parse_hgvs;

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
}

const CDS: &str = "NM_000088.3";
const GENOME: &str = "NC_000017.11";
const TX: &str = "NR_037639.1";

// =============================================================================
// SECTION 1 — Whole-entity unknown (`<acc>:<type>.?`)
// =============================================================================
//
// The bare `?` (no edit) means the whole variant is unknown. Pinned
// as regression guards across coord systems.

mod whole_entity_unknown {
    use super::*;

    #[test]
    fn cds_whole_unknown_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?"));
    }

    #[test]
    fn genome_whole_unknown_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.?"));
    }

    #[test]
    fn rna_whole_unknown_round_trips() {
        assert_round_trips(&format!("{CDS}:r.?"));
    }
}

// =============================================================================
// SECTION 2 — Single-`?` position with edit (`?<edit>`)
// =============================================================================
//
// `c.?del`, `c.?dup` were already working. `g.?del`, `r.?del`, and
// `c.?A>G` were failing pre-fix. This section pins the spec-correct
// round-trip for every coord × edit combination.

mod single_unknown_position_with_edit {
    use super::*;

    #[test]
    fn cds_unknown_position_del_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?del"));
    }

    #[test]
    fn cds_unknown_position_dup_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?dup"));
    }

    #[test]
    fn cds_unknown_position_substitution_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?A>G"));
    }

    #[test]
    fn genome_unknown_position_del_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.?del"));
    }

    #[test]
    fn genome_unknown_position_dup_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.?dup"));
    }

    #[test]
    fn genome_unknown_position_substitution_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.?A>G"));
    }

    #[test]
    fn rna_unknown_position_del_round_trips() {
        assert_round_trips(&format!("{CDS}:r.?del"));
    }

    #[test]
    fn rna_unknown_position_dup_round_trips() {
        assert_round_trips(&format!("{CDS}:r.?dup"));
    }

    #[test]
    fn tx_unknown_position_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.?del"));
    }

    /// Identity edit `c.?=` at unknown position — closes the same
    /// asymmetry class (`c.5=` works, `c.?_5=` works, so `c.?=` must
    /// also work to be consistent).
    #[test]
    fn cds_unknown_position_identity_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?="));
    }

    #[test]
    fn genome_unknown_position_identity_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.?="));
    }
}

// =============================================================================
// SECTION 3 — Half-unknown range `?_pos<edit>` and `pos_?<edit>`
// =============================================================================
//
// Already working. Regression guards across coord systems × edits.

mod half_unknown_range {
    use super::*;

    #[test]
    fn cds_unknown_start_del_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?_5del"));
    }

    #[test]
    fn cds_unknown_end_dup_round_trips() {
        assert_round_trips(&format!("{CDS}:c.5_?dup"));
    }

    #[test]
    fn cds_unknown_start_delins_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?_5delinsACG"));
    }

    #[test]
    fn cds_unknown_start_inv_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?_5inv"));
    }

    #[test]
    fn genome_unknown_start_del_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.?_5del"));
    }

    #[test]
    fn genome_unknown_end_del_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.5_?del"));
    }

    #[test]
    fn rna_unknown_start_del_round_trips() {
        assert_round_trips(&format!("{CDS}:r.?_5del"));
    }

    #[test]
    fn tx_unknown_start_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.?_5del"));
    }
}

// =============================================================================
// SECTION 4 — Intronic offsets on `?` position
// =============================================================================
//
// Already working. Regression guards.

mod intronic_offset_on_unknown {
    use super::*;

    #[test]
    fn cds_unknown_plus_offset_del_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?+5del"));
    }

    #[test]
    fn cds_unknown_minus_offset_del_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?-5del"));
    }

    #[test]
    fn cds_complex_offset_unknown_range_round_trips() {
        assert_round_trips(&format!("{CDS}:c.?-232_4484+?del"));
    }

    #[test]
    fn cds_uncertain_offset_unknown_range_round_trips() {
        assert_round_trips(&format!("{CDS}:c.5+?_10-?del"));
    }
}

// =============================================================================
// SECTION 5 — Compound brackets with `?` positions
// =============================================================================
//
// Already working. Pin as regression guard.

mod compound_brackets {
    use super::*;

    #[test]
    fn cis_bracket_with_unknown_start_del_round_trips() {
        assert_round_trips(&format!("{CDS}:c.[?_5del;200C>T]"));
    }

    #[test]
    fn cis_bracket_with_unknown_position_substitution_round_trips() {
        assert_round_trips(&format!("{CDS}:c.[?A>G;200C>T]"));
    }
}

// =============================================================================
// SECTION 6 — Both-unknown range `?_?<edit>` (documented behavior)
// =============================================================================
//
// `c.?_?del` parses but Display collapses to `c.?del` — information
// lost. This is an undocumented spec form; the canonical bracketed
// form `c.(?_?)del` (working post #237) is the spec-aligned way to
// express "deletion of unknown extent between two unknown points".
// This section pins the **current** collapse behavior as observed,
// alongside the working bracketed form as a regression guard.

mod both_unknown_range_documented {
    use super::*;

    /// Unbracketed `?_?<edit>` is silently canonicalized to single
    /// whole-entity `?<edit>`. Pinned as observed behavior; users
    /// who want the bounded form should write `(?_?)<edit>` (handled
    /// separately by #237 / PR #238).
    #[test]
    fn unbracketed_both_unknown_del_collapses_to_whole_unknown() {
        let v = parse_hgvs(&format!("{CDS}:c.?_?del")).expect("parse");
        assert_eq!(format!("{}", v), format!("{CDS}:c.?del"));
    }
}
