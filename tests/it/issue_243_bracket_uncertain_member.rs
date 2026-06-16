//! Audit for #243 — compound brackets with uncertain-edit members,
//! follow-up to #241 (PR #242).
//!
//! At the variant level the predicted-edit wrapper `c.(<pos><edit>)`
//! is now supported across every NA coord system (#241). The bracket
//! parsers (`parse_cds_allele_shorthand`, `parse_genome_compound_allele`,
//! …) did not yet recognise the wrapper on per-member parses, so
//! inputs like `c.[(123A>G);125C>T]` failed.
//!
//! After this fix:
//!
//! - bracketed cis alleles round-trip with one or more predicted members,
//! - trans alleles round-trip with one or more predicted members,
//! - mixed forms (cis + trans + unknown-phase) preserve member-level
//!   wrappers,
//! - the working non-predicted forms are pinned as regression guards.

use ferro_hgvs::parse_hgvs;

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
}

const CDS: &str = "NM_000088.3";
const GENOME: &str = "NC_000017.11";
const MITO: &str = "NC_012920.1";
const CIRCULAR: &str = "AC_X.1";

// =============================================================================
// SECTION 1 — Cis bracket with one predicted member
// =============================================================================

mod cis_one_predicted {
    use super::*;

    #[test]
    fn cds_first_member_predicted() {
        assert_round_trips(&format!("{CDS}:c.[(123A>G);125C>T]"));
    }

    #[test]
    fn cds_second_member_predicted() {
        assert_round_trips(&format!("{CDS}:c.[123A>G;(125C>T)]"));
    }

    #[test]
    fn genome_first_member_predicted() {
        assert_round_trips(&format!("{GENOME}:g.[(123A>G);125C>T]"));
    }
}

// =============================================================================
// SECTION 2 — Cis bracket with both members predicted
// =============================================================================

mod cis_both_predicted {
    use super::*;

    #[test]
    fn cds_both_members_predicted() {
        assert_round_trips(&format!("{CDS}:c.[(123A>G);(125C>T)]"));
    }

    #[test]
    fn genome_both_members_predicted() {
        assert_round_trips(&format!("{GENOME}:g.[(123A>G);(125C>T)]"));
    }
}

// =============================================================================
// SECTION 3 — Trans bracket with predicted member(s)
// =============================================================================

mod trans_with_predicted {
    use super::*;

    #[test]
    fn cds_first_allele_predicted() {
        assert_round_trips(&format!("{CDS}:c.[(123A>G)];[125C>T]"));
    }

    #[test]
    fn cds_both_alleles_predicted() {
        assert_round_trips(&format!("{CDS}:c.[(123A>G)];[(125C>T)]"));
    }

    // g./m./o. route `];[` inputs through their own *_trans_allele_shorthand
    // parsers (not the cis compound-allele loop), so the predicted-wrapper
    // recognition has to be wired in there too. These pin that path.

    #[test]
    fn genome_first_allele_predicted() {
        assert_round_trips(&format!("{GENOME}:g.[(123A>G)];[125C>T]"));
    }

    #[test]
    fn genome_both_alleles_predicted() {
        assert_round_trips(&format!("{GENOME}:g.[(123A>G)];[(125C>T)]"));
    }

    #[test]
    fn mito_first_allele_predicted() {
        assert_round_trips(&format!("{MITO}:m.[(100A>G)];[200T>C]"));
    }

    #[test]
    fn mito_both_alleles_predicted() {
        assert_round_trips(&format!("{MITO}:m.[(100A>G)];[(200T>C)]"));
    }

    #[test]
    fn circular_first_allele_predicted() {
        assert_round_trips(&format!("{CIRCULAR}:o.[(100A>G)];[200T>C]"));
    }

    #[test]
    fn circular_both_alleles_predicted() {
        assert_round_trips(&format!("{CIRCULAR}:o.[(100A>G)];[(200T>C)]"));
    }
}

// =============================================================================
// SECTION 4 — Non-sub edits as predicted bracket members
// =============================================================================

mod non_sub_predicted_in_bracket {
    use super::*;

    #[test]
    fn cds_predicted_del_in_cis_bracket() {
        assert_round_trips(&format!("{CDS}:c.[(123_125del);200C>T]"));
    }

    #[test]
    fn cds_predicted_dup_in_cis_bracket() {
        assert_round_trips(&format!("{CDS}:c.[(123_125dup);200C>T]"));
    }

    #[test]
    fn cds_predicted_delins_in_cis_bracket() {
        assert_round_trips(&format!("{CDS}:c.[(123_125delinsACG);200C>T]"));
    }
}

// =============================================================================
// SECTION 5 — Regression guards (non-predicted bracket forms)
// =============================================================================
//
// Pin the working forms so adding predicted-member support doesn't
// accidentally regress them.

mod regression_guards {
    use super::*;

    #[test]
    fn cis_two_certain_subs_still_round_trips() {
        assert_round_trips(&format!("{CDS}:c.[123A>G;125C>T]"));
    }

    #[test]
    fn trans_two_certain_subs_still_round_trips() {
        assert_round_trips(&format!("{CDS}:c.[123A>G];[125C>T]"));
    }
}
