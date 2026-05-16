//! Audit for #81 § G1 — uncertain position ranges.
//!
//! Per HGVS v21.0 § `recommendations/uncertain.md` the form
//! `c.(a_b)<edit>` describes a single variant at one position
//! somewhere in `[a, b]` — *not* a multi-position range edit.
//! The spec example is verbatim:
//!
//! > `NC_000023.10:g.(33038277_33038278)C>T` (`LRG_199t1:c.(71_72)G>A`)
//! >  describes the variant `p.Trp24Ter` in the _DMD_ gene,
//! >  reported on protein level only.
//!
//! Before the fix, ferro emitted `c.(a)_(b)<edit>` (semantic loss:
//! reads as a multi-base range edit, not a single uncertain
//! position), and `g.(?_b)del` / `c.(a+1_b-1)del` rendered with a
//! duplicated boundary (`(?_b)_(?_b)del`). This file pins the
//! spec-correct round-trip for every coord system × edit shape
//! after the fix.
//!
//! The semantic distinction is:
//!
//! - `(a_b)<edit>`   → *single* uncertain position bounded by `[a, b]`.
//! - `(a)_(b)<edit>` → range from "approximately a" to "approximately b".
//!
//! Both forms remain valid; they are not interconvertible.

use ferro_hgvs::parse_hgvs;

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
}

// =============================================================================
// SECTION 1 — Substitution at a single uncertain position
// =============================================================================
//
// Spec: `LRG_199t1:c.(71_72)G>A` and
// `NC_000023.10:g.(33038277_33038278)C>T` from
// `recommendations/uncertain.md`.

mod substitution {
    use super::*;

    #[test]
    fn cds_single_uncertain_position_round_trips() {
        assert_round_trips("NM_000088.3:c.(123_127)A>G");
    }

    #[test]
    fn genome_single_uncertain_position_round_trips() {
        assert_round_trips("NC_000017.11:g.(123_127)A>G");
    }

    #[test]
    fn rna_single_uncertain_position_round_trips() {
        assert_round_trips("NM_000088.3:r.(123_127)a>g");
    }

    #[test]
    fn ncbi_tx_single_uncertain_position_round_trips() {
        assert_round_trips("NR_037639.1:n.(123_127)A>G");
    }

    /// Verbatim spec example.
    #[test]
    fn spec_example_dmd_protein_level_round_trips() {
        assert_round_trips("LRG_199t1:c.(71_72)G>A");
    }
}

// =============================================================================
// SECTION 1b — Protein (p.) single uncertain residue bounded by a range
// =============================================================================
//
// Protein uses three-letter amino acid codes but the surrounding `(a_b)`
// uncertain-position semantics are identical: `p.(Ala123_Pro131)Ter`
// means "a Ter substitution at one residue somewhere in [Ala123, Pro131]".

mod protein {
    use super::*;

    #[test]
    fn protein_ter_at_single_uncertain_residue_round_trips() {
        assert_round_trips("NP_003997.1:p.(Ala123_Pro131)Ter");
    }

    #[test]
    fn protein_fs_at_single_uncertain_residue_round_trips() {
        assert_round_trips("NP_003997.1:p.(Ala123_Pro131)fs");
    }

    #[test]
    fn protein_del_at_single_uncertain_residue_round_trips() {
        assert_round_trips("NP_003997.1:p.(Ala123_Pro131)del");
    }
}

// =============================================================================
// SECTION 2 — Deletion at a single uncertain position
// =============================================================================

mod deletion {
    use super::*;

    #[test]
    fn cds_single_uncertain_position_del_round_trips() {
        assert_round_trips("NM_000088.3:c.(123_127)del");
    }

    #[test]
    fn genome_single_uncertain_position_del_round_trips() {
        assert_round_trips("NC_000017.11:g.(123_127)del");
    }
}

// =============================================================================
// SECTION 3 — Duplication, inversion, delins at a single uncertain position
// =============================================================================

mod other_edits {
    use super::*;

    #[test]
    fn cds_single_uncertain_position_dup_round_trips() {
        assert_round_trips("NM_000088.3:c.(123_127)dup");
    }

    #[test]
    fn cds_single_uncertain_position_inv_round_trips() {
        assert_round_trips("NM_000088.3:c.(123_127)inv");
    }

    #[test]
    fn cds_single_uncertain_position_delins_round_trips() {
        assert_round_trips("NM_000088.3:c.(123_127)delinsACG");
    }
}

// =============================================================================
// SECTION 4 — Half-unknown bounds (?_pos), (pos_?)
// =============================================================================
//
// The lower or upper bound is unknown ("somewhere ≤ b" or "somewhere ≥ a").
// Spec uses these in array / FISH / MLPA descriptions where flanking probes
// were not tested.

mod half_unknown {
    use super::*;

    #[test]
    fn cds_lower_unknown_del_round_trips() {
        assert_round_trips("NM_000088.3:c.(?_127)del");
    }

    #[test]
    fn cds_upper_unknown_dup_round_trips() {
        assert_round_trips("NM_000088.3:c.(123_?)dup");
    }

    #[test]
    fn genome_lower_unknown_del_round_trips() {
        assert_round_trips("NC_000017.11:g.(?_127)del");
    }

    #[test]
    fn genome_upper_unknown_del_round_trips() {
        assert_round_trips("NC_000017.11:g.(123_?)del");
    }

    #[test]
    fn cds_both_unknown_del_round_trips() {
        assert_round_trips("NM_000088.3:c.(?_?)del");
    }
}

// =============================================================================
// SECTION 5 — Intronic-offset uncertain bounds
// =============================================================================

mod intronic_offset {
    use super::*;

    #[test]
    fn cds_intronic_offset_uncertain_del_round_trips() {
        assert_round_trips("NM_000088.3:c.(123+1_127-1)del");
    }

    #[test]
    fn cds_intronic_uncertain_offset_uncertain_bound_round_trips() {
        assert_round_trips("NM_000088.3:c.(123+?_127-?)del");
    }
}

// =============================================================================
// SECTION 6 — Spec-form preservation across compound brackets
// =============================================================================
//
// The spec form must survive inside a compound `c.[...;...]` allele
// bracket as well.

mod compound_brackets {
    use super::*;

    #[test]
    fn uncertain_range_inside_cis_bracket_round_trips() {
        assert_round_trips("NM_000088.3:c.[(123_127)A>G;200C>T]");
    }
}

// =============================================================================
// SECTION 7 — Nested forms (regression guards — already worked before #237)
// =============================================================================
//
// These are NOT the bug; they are the working "range_range" form. Pin so
// the fix does not regress them.

mod nested_range_regression_guards {
    use super::*;

    #[test]
    fn nested_two_ranges_round_trips() {
        assert_round_trips("NM_000088.3:c.(123_127)_(200_210)del");
    }

    #[test]
    fn nested_half_unknown_round_trips() {
        assert_round_trips("NM_000088.3:c.(?_127)_(200_?)del");
    }

    #[test]
    fn nested_range_then_certain_round_trips() {
        assert_round_trips("NM_000088.3:c.(123_127)_200del");
    }

    #[test]
    fn nested_certain_then_range_round_trips() {
        assert_round_trips("NM_000088.3:c.123_(200_210)del");
    }

    /// RNA two-paren range must round-trip — the `(a_b)` arm must not
    /// greedily consume the first paren when `_(c_d)` follows.
    #[test]
    fn rna_nested_two_ranges_round_trips() {
        assert_round_trips("NM_000088.3:r.(100_200)_(300_400)del");
    }

    /// RNA half-unknown two-paren range.
    #[test]
    fn rna_nested_half_unknown_round_trips() {
        assert_round_trips("NM_000088.3:r.(?_127)_(200_?)del");
    }

    /// RNA range, then certain, two-paren form.
    #[test]
    fn rna_nested_range_then_certain_round_trips() {
        assert_round_trips("NM_000088.3:r.(123_127)_200del");
    }

    /// RNA certain, then range, two-paren form.
    #[test]
    fn rna_nested_certain_then_range_round_trips() {
        assert_round_trips("NM_000088.3:r.123_(200_210)del");
    }
}

// =============================================================================
// SECTION 8 — Semantic distinction (c.(a)_(b) NOT confused with c.(a_b))
// =============================================================================
//
// The two-paren form `(a)_(b)<edit>` means "range from approximately
// a to approximately b" — a multi-base edit with both endpoints
// uncertain. This is distinct from `(a_b)<edit>` (single uncertain
// position). After the fix, both forms must round-trip cleanly and
// not be conflated.

mod two_paren_range_distinct {
    use super::*;

    /// The two-paren form `c.(a)_(b)del` must round-trip — it is a
    /// legitimate spec form, distinct from `c.(a_b)del`.
    #[test]
    fn two_paren_range_del_round_trips() {
        assert_round_trips("NC_000001.11:g.(100)_(200)del");
    }

    /// Half-paren mix: certain start, uncertain end.
    #[test]
    fn certain_start_uncertain_end_round_trips() {
        assert_round_trips("NC_000001.11:g.100_(200)del");
    }

    /// The two-paren form is NOT canonicalized to the single-paren
    /// form (they have different meanings).
    #[test]
    fn two_paren_form_is_not_collapsed_to_single_paren_form() {
        let a = parse_hgvs("NC_000017.11:g.(123_127)A>G").expect("parse a");
        let b = parse_hgvs("NC_000017.11:g.(123)_(127)A>G").expect("parse b");
        assert_ne!(
            format!("{}", a),
            format!("{}", b),
            "single-paren and two-paren forms must Display distinctly"
        );
        assert_eq!(format!("{}", a), "NC_000017.11:g.(123_127)A>G");
        assert_eq!(format!("{}", b), "NC_000017.11:g.(123)_(127)A>G");
    }
}
