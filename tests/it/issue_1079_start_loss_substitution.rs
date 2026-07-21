//! MUST-level rejection of a translation-initiation-codon variant described
//! as a plain amino-acid substitution — `p.Met1Thr`, `p.(Met1Val)` (#1079).
//!
//! # Spec
//!
//! `recommendations/protein/substitution.md:49`:
//!
//! > **no protein: `LRG_199p1:p.0`** … `LRG_199p1:p.0?` can be used when you
//! > predict that no protein is produced. Do not use descriptions like
//! > `p.Met1Thr`, this is for sure **not** the consequence of the effect on
//! > protein translation.
//!
//! `recommendations/checklist.md:65`:
//!
//! > the description `p.(Met1Val)` is not allowed (see Protein).
//!
//! "Do not use" / "is not allowed" is MUST-level under the spec's RFC 2119
//! reading (`recommendations/style.md:9`).
//!
//! `recommendations/protein/extension.md:28` rules out the sibling
//! `p.Met1Valext-4` form for the same reason:
//!
//! > this variant is **not** described as an extension (`p.Met1Valext-4`)
//! > since `Met1`, part of the normal amino acid sequence, is changed.
//!
//! # Why reject rather than canonicalise
//!
//! The spec names three correct forms — `p.0` (no protein produced), `p.0?`
//! (predicted no protein) and `p.(Met1?)` (consequence unknown) — and which
//! one applies depends on evidence the description does not carry. Picking
//! one for the author would assert a finding they never made, so every mode
//! rejects and the diagnostic lists all three.

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::parse_hgvs;

const START_LOSS_SUBSTITUTIONS: &[&str] = &[
    "NP_003997.1:p.Met1Thr",
    "NP_003997.1:p.(Met1Val)",
    "NP_003997.1:p.Met1Val",
    "NP_003997.1:p.M1T",
];

#[test]
fn default_parse_rejects_start_loss_substitution() {
    for input in START_LOSS_SUBSTITUTIONS {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "substitution.md:49 forbids {input:?}; got {:?}",
            result.map(|v| v.to_string())
        );
    }
}

#[test]
fn every_mode_rejects_start_loss_substitution() {
    for config in [
        ErrorConfig::strict(),
        ErrorConfig::lenient(),
        ErrorConfig::silent(),
    ] {
        assert!(
            parse_hgvs_with_config("NP_003997.1:p.Met1Thr", config).is_err(),
            "no mode may accept a start-loss substitution"
        );
    }
}

#[test]
fn rejection_lists_the_three_spec_forms() {
    let msg = parse_hgvs("NP_003997.1:p.Met1Thr").unwrap_err().to_string();
    for expected in ["p.0", "p.0?", "p.(Met1?)"] {
        assert!(
            msg.contains(expected),
            "the diagnostic should offer {expected}; got: {msg}"
        );
    }
}

// =====================================================================
// The spec-sanctioned start-codon descriptions still parse
// =====================================================================

#[test]
fn spec_sanctioned_start_codon_forms_parse() {
    for input in [
        "NP_003997.1:p.0",
        "NP_003997.1:p.0?",
        "NP_003997.1:p.(Met1?)",
        "NP_003997.1:p.Met1?",
        "NP_003997.1:p.Met1=",
        "NP_003997.1:p.Met1del",
        "NP_003997.1:p.Leu2_Met124del",
        "NP_003997.1:p.Met1_Leu2insArgSerThrVal",
    ] {
        assert!(
            parse_hgvs(input).is_ok(),
            "{input:?} is a spec example and must parse"
        );
    }
}

/// An N-terminal extension that *renames* `Met1` is forbidden for the same
/// reason (`protein/extension.md:28`: "this variant is **not** described as
/// an extension … since `Met1`, part of the normal amino acid sequence, is
/// changed"). Its replacement is the insertion form, which needs the
/// inserted residues, so there is nothing to canonicalise to.
#[test]
fn renaming_extension_at_the_initiator_is_rejected() {
    assert!(parse_hgvs("NP_003997.1:p.Met1Valext-4").is_err());
}

/// A plain N-terminal extension leaves `Met1` alone and is the spec's own
/// example, so it must keep parsing.
#[test]
fn plain_n_terminal_extension_is_untouched() {
    for input in [
        "NP_003997.1:p.Met1ext-5",
        "NP_003997.2:p.Met1ext-5",
        "NP_003997.1:p.(Met1ext-8)",
    ] {
        assert!(parse_hgvs(input).is_ok(), "{input:?} must still parse");
    }
}

/// Substitutions at any other position, and at position 1 of a reference
/// whose first residue is not the initiator methionine, are unaffected.
#[test]
fn ordinary_substitutions_are_untouched() {
    for input in [
        "NP_003997.1:p.Met2Thr",
        "NP_003997.1:p.Trp24Cys",
        "NP_003997.1:p.(Arg727Ser)",
    ] {
        assert!(parse_hgvs(input).is_ok(), "{input:?} must still parse");
    }
}
