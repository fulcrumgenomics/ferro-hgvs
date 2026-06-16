//! Audit for #287 — RNA / TX trans-allele predicted-edit wrapper inside
//! compound brackets. Follow-up to PR #244 (which closed #243).
//!
//! PR #244 wired the predicted-edit wrapper `(<pos><edit>)` into the
//! `c.` shorthand bracket parsers (`parse_cds_allele_shorthand` and
//! `parse_cds_trans_allele_shorthand`). The RNA and TX counterparts —
//! `parse_rna_trans_allele_shorthand`, `parse_rna_allele_shorthand`,
//! `parse_tx_trans_allele_shorthand`, `parse_tx_compound_allele` —
//! were left on the bare `parse_<coord>_interval` + `parse_na_edit`
//! path, so inputs like `r.[(123a>g);125c>u]` and
//! `[NM_…:n.(123A>G)];[NM_…:n.125C>T]` failed.
//!
//! After this fix:
//!
//! - r. cis brackets round-trip with one or both members predicted,
//! - r. trans brackets round-trip with one or both members predicted,
//! - n. cis brackets round-trip with one or both members predicted,
//! - n. trans brackets round-trip with one or both members predicted,
//! - non-sub predicted bracket members (del / dup / delins) round-trip,
//! - the previously-working c. forms from PR #244 still round-trip
//!   (regression guards).
//!
//! Stacked on PR #244 — not mergeable until #244 lands.

use ferro_hgvs::parse_hgvs;

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
}

const TX: &str = "NM_000088.3";

// =============================================================================
// SECTION 1 — RNA cis bracket with predicted member(s)
// =============================================================================

mod rna_cis {
    use super::*;

    #[test]
    fn rna_first_member_predicted() {
        assert_round_trips(&format!("{TX}:r.[(123a>g);125c>u]"));
    }

    #[test]
    fn rna_second_member_predicted() {
        assert_round_trips(&format!("{TX}:r.[123a>g;(125c>u)]"));
    }

    #[test]
    fn rna_both_members_predicted() {
        assert_round_trips(&format!("{TX}:r.[(123a>g);(125c>u)]"));
    }
}

// =============================================================================
// SECTION 2 — RNA trans bracket with predicted member(s)
// =============================================================================

mod rna_trans {
    use super::*;

    #[test]
    fn rna_first_allele_predicted() {
        assert_round_trips(&format!("{TX}:r.[(123a>g)];[125c>u]"));
    }

    #[test]
    fn rna_second_allele_predicted() {
        assert_round_trips(&format!("{TX}:r.[123a>g];[(125c>u)]"));
    }

    #[test]
    fn rna_both_alleles_predicted() {
        assert_round_trips(&format!("{TX}:r.[(123a>g)];[(125c>u)]"));
    }
}

// =============================================================================
// SECTION 3 — TX (n.) cis bracket with predicted member(s)
// =============================================================================

mod tx_cis {
    use super::*;

    #[test]
    fn tx_first_member_predicted() {
        assert_round_trips(&format!("{TX}:n.[(123A>G);125C>T]"));
    }

    #[test]
    fn tx_second_member_predicted() {
        assert_round_trips(&format!("{TX}:n.[123A>G;(125C>T)]"));
    }

    #[test]
    fn tx_both_members_predicted() {
        assert_round_trips(&format!("{TX}:n.[(123A>G);(125C>T)]"));
    }
}

// =============================================================================
// SECTION 4 — TX (n.) trans bracket with predicted member(s)
// =============================================================================

mod tx_trans {
    use super::*;

    #[test]
    fn tx_first_allele_predicted() {
        assert_round_trips(&format!("{TX}:n.[(123A>G)];[125C>T]"));
    }

    #[test]
    fn tx_second_allele_predicted() {
        assert_round_trips(&format!("{TX}:n.[123A>G];[(125C>T)]"));
    }

    #[test]
    fn tx_both_alleles_predicted() {
        assert_round_trips(&format!("{TX}:n.[(123A>G)];[(125C>T)]"));
    }
}

// =============================================================================
// SECTION 5 — Non-sub edits as predicted bracket members
// =============================================================================

mod non_sub_predicted_in_bracket {
    use super::*;

    #[test]
    fn rna_predicted_del_in_cis_bracket() {
        assert_round_trips(&format!("{TX}:r.[(123_125del);200c>u]"));
    }

    #[test]
    fn rna_predicted_dup_in_cis_bracket() {
        assert_round_trips(&format!("{TX}:r.[(123_125dup);200c>u]"));
    }

    #[test]
    fn tx_predicted_delins_in_cis_bracket() {
        assert_round_trips(&format!("{TX}:n.[(123_125delinsACG);200C>T]"));
    }

    #[test]
    fn tx_predicted_del_in_trans_bracket() {
        assert_round_trips(&format!("{TX}:n.[(123_125del)];[200C>T]"));
    }
}

// =============================================================================
// SECTION 6 — Regression guards
// =============================================================================
//
// Pin the c. shorthand forms wired up by PR #244 plus the previously-working
// non-predicted RNA / TX forms so this fix doesn't accidentally regress them.

mod regression_guards {
    use super::*;

    #[test]
    fn cds_first_member_predicted_still_works() {
        assert_round_trips(&format!("{TX}:c.[(123A>G);125C>T]"));
    }

    #[test]
    fn cds_both_members_predicted_still_works() {
        assert_round_trips(&format!("{TX}:c.[(123A>G);(125C>T)]"));
    }

    #[test]
    fn cds_first_allele_predicted_still_works() {
        assert_round_trips(&format!("{TX}:c.[(123A>G)];[125C>T]"));
    }

    #[test]
    fn cds_both_alleles_predicted_still_works() {
        assert_round_trips(&format!("{TX}:c.[(123A>G)];[(125C>T)]"));
    }

    #[test]
    fn rna_cis_two_certain_subs_still_round_trips() {
        assert_round_trips(&format!("{TX}:r.[123a>g;125c>u]"));
    }

    #[test]
    fn rna_trans_two_certain_subs_still_round_trips() {
        assert_round_trips(&format!("{TX}:r.[123a>g];[125c>u]"));
    }

    #[test]
    fn tx_cis_two_certain_subs_still_round_trips() {
        assert_round_trips(&format!("{TX}:n.[123A>G;125C>T]"));
    }

    #[test]
    fn tx_trans_two_certain_subs_still_round_trips() {
        assert_round_trips(&format!("{TX}:n.[123A>G];[125C>T]"));
    }
}
