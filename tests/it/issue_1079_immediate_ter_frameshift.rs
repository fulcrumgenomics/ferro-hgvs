//! MUST-level rejection of a frameshift that terminates immediately —
//! `p.Tyr4TerfsTer1` (#1079).
//!
//! # Spec
//!
//! `recommendations/protein/frameshift.md:22`:
//!
//! > **NOTE**: the shortest frameshift variant possible contains `fsTer2`;
//! > variants which introduce an **immediate** translation termination (stop)
//! > codon are described as nonsense variant, e.g., `p.Tyr4Ter` (or
//! > `p.Tyr4*`) not `p.Tyr4TerfsTer1`.
//!
//! `recommendations/protein/substitution.md:20` says the same from the
//! substitution side ("**NOTE**: not `p.Tyr4TerfsTer1`, but `p.Tyr4Ter`").
//!
//! # `fsTer2` is legal
//!
//! The spec names `fsTer2` as the *shortest possible* frameshift, so only
//! `fsTer1` — and the equivalent short form whose first new residue is
//! already `Ter` — is wrong. Everything from `fsTer2` up must keep parsing.
//!
//! # Why reject rather than canonicalise
//!
//! The repair is mechanical (`p.Tyr4Ter`) and the diagnostic names it
//! verbatim, but applying it silently would change the *edit kind* the
//! author wrote — a frameshift becomes a substitution — so it is surfaced
//! rather than performed.

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::parse_hgvs;

const IMMEDIATE_TER_FRAMESHIFTS: &[&str] = &[
    "NP_003997.1:p.Tyr4TerfsTer1",
    "NP_003997.1:p.(Tyr4TerfsTer1)",
    "NP_003997.1:p.Arg97ProfsTer1",
    "NP_003997.1:p.Tyr4Terfs",
    "NP_003997.1:p.Tyr4TerfsTer10",
];

#[test]
fn default_parse_rejects_immediate_ter_frameshift() {
    for input in IMMEDIATE_TER_FRAMESHIFTS {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "frameshift.md:22 forbids {input:?}; got {:?}",
            result.map(|v| v.to_string())
        );
    }
}

#[test]
fn every_mode_rejects_immediate_ter_frameshift() {
    for config in [
        ErrorConfig::strict(),
        ErrorConfig::lenient(),
        ErrorConfig::silent(),
    ] {
        assert!(
            parse_hgvs_with_config("NP_003997.1:p.Tyr4TerfsTer1", config).is_err(),
            "no mode may accept an immediately-terminating frameshift"
        );
    }
}

#[test]
fn rejection_names_the_nonsense_form() {
    let msg = parse_hgvs("NP_003997.1:p.Tyr4TerfsTer1")
        .unwrap_err()
        .to_string();
    assert!(
        msg.contains("p.Tyr4Ter"),
        "the diagnostic should name the nonsense form; got: {msg}"
    );
}

// =====================================================================
// `fsTer2` and beyond are legal and must not be touched
// =====================================================================

#[test]
fn shortest_legal_frameshift_still_parses() {
    for input in [
        "NP_003997.1:p.Tyr4ValfsTer2",
        "NP_003997.1:p.Arg456GlyfsTer17",
        "NP_003997.1:p.Arg97ProfsTer23",
        "NP_003997.1:p.Glu5ValfsTer5",
        "NP_003997.1:p.Ile327ArgfsTer?",
        "NP_003997.1:p.Arg97fs",
        "NP_003997.1:p.His150HisfsTer10",
    ] {
        assert!(
            parse_hgvs(input).is_ok(),
            "{input:?} is a legal frameshift and must parse"
        );
    }
}

/// The spec's own nonsense form — the thing `p.Tyr4TerfsTer1` should have
/// been written as — must of course parse.
#[test]
fn nonsense_substitution_still_parses() {
    for input in ["NP_003997.1:p.Tyr4Ter", "NP_003997.1:p.Tyr4*"] {
        assert!(parse_hgvs(input).is_ok(), "{input:?} must parse");
    }
}
