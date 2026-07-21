//! MUST-level rejection of amino acids listed *after* a translation stop in
//! an inserted peptide — `p.Cys5_Ser6delinsTerGluAsp` (#1079).
//!
//! # Spec
//!
//! `recommendations/protein/delins.md:45`:
//!
//! > **NOTE**: the deletion-insertion is not described as
//! > `delinsSerSerTerAlaAsp`, amino acids after the translation termination
//! > codon are **not** listed.
//!
//! `recommendations/protein/insertion.md:43` states the same for insertions:
//!
//! > **NOTE**: the insertion is not described as `insSerSerTerAlaPro`; amino
//! > acids after the translation termination codon are not listed.
//!
//! Translation stops at the first `Ter`, so residues written after it are
//! not part of any protein product — listing them describes a sequence that
//! cannot exist. Under the spec's RFC 2119 reading
//! (`recommendations/style.md:9`) this is MUST-level.
//!
//! # Why reject rather than truncate
//!
//! Truncating at the first `Ter` is mechanical, and the diagnostic names the
//! truncated form. It is not applied automatically because it is not always
//! the spec's own correct description: when the insert *begins* with `Ter`,
//! `protein/substitution.md:20` says the variant is a nonsense substitution
//! at the preceding residue (`p.Cys5_Ser6delinsTerGluAsp` → `p.Tyr4Ter`),
//! which needs the reference sequence. Truncating would swap one
//! non-canonical description for another.

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::parse_hgvs;

const RESIDUES_AFTER_STOP: &[&str] = &[
    "NP_003997.1:p.Cys5_Ser6delinsTerGluAsp",
    "NP_004371.2:p.Asn47delinsSerSerTerAlaAsp",
    "NP_004371.2:p.Pro46_Asn47insSerSerTerAlaPro",
    "NP_003997.1:p.(Cys5_Ser6delinsTerGluAsp)",
];

#[test]
fn default_parse_rejects_residues_after_stop() {
    for input in RESIDUES_AFTER_STOP {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "delins.md:45 / insertion.md:43 forbid {input:?}; got {:?}",
            result.map(|v| v.to_string())
        );
    }
}

#[test]
fn every_mode_rejects_residues_after_stop() {
    for config in [
        ErrorConfig::strict(),
        ErrorConfig::lenient(),
        ErrorConfig::silent(),
    ] {
        assert!(
            parse_hgvs_with_config("NP_003997.1:p.Cys5_Ser6delinsTerGluAsp", config).is_err(),
            "no mode may accept residues listed after a stop"
        );
    }
}

#[test]
fn rejection_names_the_truncated_form() {
    let msg = parse_hgvs("NP_004371.2:p.Asn47delinsSerSerTerAlaAsp")
        .unwrap_err()
        .to_string();
    assert!(
        msg.contains("SerSerTer"),
        "the diagnostic should name the truncated peptide; got: {msg}"
    );
}

// =====================================================================
// Inserts that stop at their own end — every spec example must parse
// =====================================================================

#[test]
fn peptides_ending_at_the_stop_still_parse() {
    for input in [
        "NP_004371.2:p.(Asn47delinsSerSerTer)",
        "NP_003997.1:p.(Pro578_Lys579delinsLeuTer)",
        "NP_003997.1:p.(Met3_His4insGlyTer)",
        "NP_004371.2:p.(Pro46_Asn47insSerSerTer)",
        "NP_003997.1:p.Cys28delinsTrpVal",
        "NP_003997.1:p.Cys28_Lys29delinsTrp",
    ] {
        assert!(
            parse_hgvs(input).is_ok(),
            "{input:?} lists nothing after the stop and must parse"
        );
    }
}

/// `insTer<n>` / `ins*<n>` gives the *position* of the stop inside the
/// inserted sequence rather than spelling residues out, so nothing can
/// follow the stop and the rule does not apply.
#[test]
fn stop_position_insertion_form_is_untouched() {
    for input in [
        "NP_003997.1:p.Lys2_Leu3insTer12",
        "NP_060250.2:p.Gln746_Lys747ins*63",
    ] {
        assert!(parse_hgvs(input).is_ok(), "{input:?} must still parse");
    }
}
