//! Audit for #81 § G3 — uncertain / predicted-edit wrapper
//! `<acc>:<type>.(<pos><edit>)`.
//!
//! Per HGVS v21 `recommendations/uncertain.md`, parentheses around
//! the entire variant denote a predicted change with no experimental
//! proof. The protein form (`p.(Arg248Gln)`) is canonical and was
//! already working before this fix; the nucleic acid coord systems
//! (c./g./r./n./m.) had asymmetric / broken parser + Display paths.
//!
//! After the fix, every coord × edit × wrapper combination should
//! round-trip cleanly, and CDS-style uncertain-edit cross-system
//! parity holds with protein.

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
const MT: &str = "NC_012920.1";
const PROT: &str = "NP_003997.1";

// =============================================================================
// SECTION 1 — Substitution under the uncertain wrapper
// =============================================================================

mod substitution {
    use super::*;

    #[test]
    fn cds_uncertain_sub_round_trips() {
        assert_round_trips(&format!("{CDS}:c.(123A>G)"));
    }

    #[test]
    fn genome_uncertain_sub_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.(123A>G)"));
    }

    #[test]
    fn rna_uncertain_sub_round_trips() {
        assert_round_trips(&format!("{CDS}:r.(123a>g)"));
    }

    #[test]
    fn tx_uncertain_sub_round_trips() {
        assert_round_trips(&format!("{TX}:n.(123A>G)"));
    }

    #[test]
    fn mt_uncertain_sub_round_trips() {
        assert_round_trips(&format!("{MT}:m.(123A>G)"));
    }

    #[test]
    fn circular_uncertain_sub_round_trips() {
        assert_round_trips(&format!("{MT}:o.(123A>G)"));
    }

    /// Protein is the spec-canonical predicted form; already worked
    /// before #237 / G3 fix — pin as a regression guard.
    #[test]
    fn protein_uncertain_sub_round_trips() {
        assert_round_trips(&format!("{PROT}:p.(Arg248Gln)"));
    }
}

// =============================================================================
// SECTION 2 — Non-substitution edits under the uncertain wrapper
// =============================================================================

mod non_sub_edits {
    use super::*;

    #[test]
    fn cds_uncertain_del_round_trips() {
        assert_round_trips(&format!("{CDS}:c.(123_125del)"));
    }

    #[test]
    fn cds_uncertain_dup_round_trips() {
        assert_round_trips(&format!("{CDS}:c.(123_125dup)"));
    }

    #[test]
    fn cds_uncertain_inv_round_trips() {
        assert_round_trips(&format!("{CDS}:c.(123_125inv)"));
    }

    #[test]
    fn cds_uncertain_delins_round_trips() {
        assert_round_trips(&format!("{CDS}:c.(123_125delinsACG)"));
    }

    #[test]
    fn cds_uncertain_ins_round_trips() {
        assert_round_trips(&format!("{CDS}:c.(123_124insACG)"));
    }

    #[test]
    fn genome_uncertain_del_round_trips() {
        assert_round_trips(&format!("{GENOME}:g.(123_125del)"));
    }

    #[test]
    fn rna_uncertain_del_round_trips() {
        assert_round_trips(&format!("{CDS}:r.(123_125del)"));
    }

    #[test]
    fn rna_uncertain_dup_round_trips() {
        assert_round_trips(&format!("{CDS}:r.(123_125dup)"));
    }
}

// =============================================================================
// SECTION 3 — RNA uppercase canonicalization inside the uncertain wrapper
// =============================================================================
//
// Uppercase input inside `r.(...)` must lowercase to `r.(...)` per E1.

mod rna_lowercase_in_wrapper {
    use super::*;

    #[track_caller]
    fn assert_canonicalizes_to(input: &str, expected: &str) {
        let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
        let out = format!("{}", v);
        assert_eq!(out, expected, "input={input:?}");
    }

    #[test]
    fn uppercase_sub_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(&format!("{CDS}:r.(123A>G)"), &format!("{CDS}:r.(123a>g)"));
    }

    #[test]
    fn uppercase_del_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(
            &format!("{CDS}:r.(123_125delAUG)"),
            &format!("{CDS}:r.(123_125delaug)"),
        );
    }
}

// =============================================================================
// SECTION 4 — Canonical re-emission across the uncertain wrapper
// =============================================================================
//
// The wrapper-internal edit obeys the same canonicalization rules as
// outside the wrapper. E.g., a sub stays a sub; the position rendering
// is unchanged.

mod canonical_re_emission {
    use super::*;

    /// The certain and uncertain forms must Display distinctly —
    /// the wrapper is semantic, not cosmetic.
    #[test]
    fn certain_and_uncertain_subs_are_distinct() {
        let certain = parse_hgvs(&format!("{CDS}:c.123A>G")).expect("certain parse");
        let uncertain = parse_hgvs(&format!("{CDS}:c.(123A>G)")).expect("uncertain parse");
        assert_eq!(format!("{}", certain), format!("{CDS}:c.123A>G"));
        assert_eq!(format!("{}", uncertain), format!("{CDS}:c.(123A>G)"));
        assert_ne!(format!("{}", certain), format!("{}", uncertain));
    }

    /// G1 form `c.(123_127)A>G` (uncertain *position*, certain edit) and
    /// G3 form `c.(123A>G)` (certain position, uncertain *edit*) are
    /// structurally distinct and must not be conflated by the parser
    /// — both start with `(`; the parser distinguishes by whether the
    /// edit is *inside* or *outside* the parens. The exact `c.(123_127)`
    /// canonical shape is pinned by the G1 audit (#237 / PR #238); this
    /// test only asserts the two forms parse to different `Display`s and
    /// that the G3 form preserves its wrapper.
    #[test]
    fn g1_uncertain_position_distinct_from_g3_uncertain_edit() {
        let g1 = parse_hgvs(&format!("{CDS}:c.(123_127)A>G")).expect("G1 parse");
        let g3 = parse_hgvs(&format!("{CDS}:c.(123A>G)")).expect("G3 parse");
        assert_eq!(format!("{}", g3), format!("{CDS}:c.(123A>G)"));
        assert_ne!(
            format!("{}", g1),
            format!("{}", g3),
            "G1 and G3 forms must Display distinctly"
        );
    }
}

// =============================================================================
// SECTION 5 — Compound brackets with uncertain edits (deferred)
// =============================================================================
//
// `c.[(123A>G);125C>T]` puts a predicted member inside a cis bracket.
// The bracket-member parser does not yet support the predicted-edit
// wrapper at the bracket level. Pinned as currently rejected — a
// follow-up will extend the bracket parser to honor the wrapper.

mod compound_brackets {
    use super::*;

    #[test]
    fn cis_bracket_with_uncertain_first_member_is_currently_rejected() {
        assert!(parse_hgvs(&format!("{CDS}:c.[(123A>G);125C>T]")).is_err());
    }

    #[test]
    fn cis_bracket_with_uncertain_both_members_is_currently_rejected() {
        assert!(parse_hgvs(&format!("{CDS}:c.[(123A>G);(125C>T)]")).is_err());
    }
}
