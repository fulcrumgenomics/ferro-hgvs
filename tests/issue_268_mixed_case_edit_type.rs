//! Tests for #268 — W1004 MixedCaseEditType soft-validation warning.
//!
//! Per HGVS general.md, edit-type tokens (`del`, `ins`, `dup`, `inv`,
//! `delins`, `con`) are spelled in lowercase. Before this PR ferro
//! produced a generic parse error in lenient mode for mixed-case
//! variants like `c.100Del`. This PR adds the full-string
//! `correct_edit_type_case_full` corrector and wires it as preprocessor
//! Phase 6b1 (before single-letter AA expansion).

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::{parse_hgvs_lenient, parse_hgvs_with_config};

#[track_caller]
fn assert_rewrites_with_w1004(input: &str, expected: &str) {
    let r = parse_hgvs_lenient(input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
    assert_eq!(
        format!("{}", r.result),
        expected,
        "rewrite mismatch for {input:?}"
    );
    assert!(
        r.warnings.iter().any(|w| w.error_type.code() == "W1004"),
        "expected W1004 for {input:?}; got {:?}",
        r.warnings
            .iter()
            .map(|w| w.error_type.code())
            .collect::<Vec<_>>()
    );
}

#[track_caller]
fn assert_no_w1004(input: &str) {
    let r = parse_hgvs_lenient(input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
    assert!(
        !r.warnings.iter().any(|w| w.error_type.code() == "W1004"),
        "unexpected W1004 for {input:?}"
    );
}

// =============================================================================
// SECTION 1 — Mixed-case edit tokens on nucleotide variants
// =============================================================================

mod nucleotide {
    use super::*;

    #[test]
    fn cds_titlecase_del() {
        assert_rewrites_with_w1004("NM_000088.3:c.100Del", "NM_000088.3:c.100del");
    }

    #[test]
    fn cds_uppercase_del() {
        assert_rewrites_with_w1004("NM_000088.3:c.100DEL", "NM_000088.3:c.100del");
    }

    #[test]
    fn cds_uppercase_ins() {
        assert_rewrites_with_w1004("NM_000088.3:c.100_101INSATG", "NM_000088.3:c.100_101insATG");
    }

    #[test]
    fn cds_uppercase_dup() {
        assert_rewrites_with_w1004("NM_000088.3:c.100_102DUPATG", "NM_000088.3:c.100_102dupATG");
    }

    #[test]
    fn cds_uppercase_inv() {
        assert_rewrites_with_w1004("NM_000088.3:c.100_102INV", "NM_000088.3:c.100_102inv");
    }

    #[test]
    fn cds_uppercase_delins() {
        assert_rewrites_with_w1004(
            "NM_000088.3:c.100_102DELINSAcg",
            "NM_000088.3:c.100_102delinsACG",
        );
    }

    /// `CON` chains with W3015 (deprecated `con` → `delins`).
    #[test]
    fn cds_uppercase_con_chains_to_delins() {
        let r = parse_hgvs_lenient("NM_000088.3:c.100_200CONNM_001.1:c.5_105").unwrap();
        let codes: Vec<_> = r
            .warnings
            .iter()
            .map(|w| w.error_type.code().to_string())
            .collect();
        assert!(codes.iter().any(|c| c == "W1004"));
        assert!(codes.iter().any(|c| c == "W3015"));
        assert_eq!(
            format!("{}", r.result),
            "NM_000088.3:c.100_200delins[NM_001.1:c.5_105]"
        );
    }

    /// Genomic accession also lowercases.
    #[test]
    fn genomic_uppercase_del() {
        assert_rewrites_with_w1004(
            "NC_000001.11:g.100_102DELAcg",
            "NC_000001.11:g.100_102delACG",
        );
    }
}

// =============================================================================
// SECTION 2 — Protein edit tokens
// =============================================================================
//
// Pin that `DEL` after a protein position (e.g. `p.Arg8_Lys10DEL`) is
// lowercased to `del`, *not* mis-expanded to `AspGluLeu` by the
// single-letter AA pass. This is the key reason Phase 6b1 runs before
// Phase 6c.

mod protein {
    use super::*;

    #[test]
    fn protein_uppercase_del() {
        assert_rewrites_with_w1004("NP_003997.1:p.Arg8_Lys10DEL", "NP_003997.1:p.Arg8_Lys10del");
    }

    #[test]
    fn protein_lowercase_del_unchanged() {
        assert_no_w1004("NP_003997.1:p.Arg8_Lys10del");
    }

    /// Asn/Glu/Ile don't get false-triggered.
    #[test]
    fn protein_three_letter_aa_unchanged() {
        assert_no_w1004("NP_003997.1:p.Asn10Ile");
    }
}

// =============================================================================
// SECTION 3 — Already-lowercase forms (no warning)
// =============================================================================

mod already_lowercase_no_warning {
    use super::*;

    #[test]
    fn cds_del_lowercase() {
        assert_no_w1004("NM_000088.3:c.100del");
    }

    #[test]
    fn cds_ins_lowercase() {
        assert_no_w1004("NM_000088.3:c.100_101insATG");
    }

    #[test]
    fn cds_dup_lowercase() {
        assert_no_w1004("NM_000088.3:c.100_102dup");
    }
}

// =============================================================================
// SECTION 4 — Strict mode rejects
// =============================================================================

mod strict_rejects {
    use super::*;

    #[test]
    fn strict_config_uppercase_del_rejected() {
        let result = parse_hgvs_with_config("NM_000088.3:c.100DEL", ErrorConfig::strict());
        assert!(result.is_err());
    }

    #[test]
    fn strict_config_uppercase_inv_rejected() {
        let result = parse_hgvs_with_config("NM_000088.3:c.100_102INV", ErrorConfig::strict());
        assert!(result.is_err());
    }
}
