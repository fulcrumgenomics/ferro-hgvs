//! `*` (asterisk) is a spec-sanctioned stop-codon spelling and must be
//! accepted in strict mode — not rejected as a deprecated form (#1114).
//!
//! # Spec
//!
//! `recommendations/checklist.md:63`:
//!
//! > **`Ter` or `*` should be used** to indicate a translation stop codon;
//! > the `X` should not be used.
//!
//! `recommendations/general.md:54` ("the `*` can be used to indicate the
//! translation stop codon in both one- and three-letter amino acid code
//! descriptions") and `protein/frameshift.md` (which lists `p.Arg97Profs*23`
//! and `p.Gln151Thrfs*9` as valid examples) confirm `*` is co-equal with
//! `Ter`. Only `X` is deprecated.
//!
//! # Behavior
//!
//! ferro previously classified `*` as `DeprecatedStopCodonStar` (W3007) /
//! `DeprecatedFrameshiftStar` (W3009) and rejected it in strict mode (the
//! deprecated-form preprocessor pass turned `success = false`). This
//! reclassifies `*` as a valid `Ter` spelling: it parses to the `Ter` AST in
//! every mode (the core parser already accepts it) and Display canonicalizes
//! to the spec-preferred three-letter `Ter`. It is therefore *not* a
//! correction — no W3007/W3009 warning is emitted. The genuinely-deprecated
//! `X` spelling stays rejected.

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;

/// (input, expected canonical display) — spec-valid `*` stop spellings across
/// the nonsense-substitution (W3007) and frameshift-termination (W3009) slots.
const STAR_STOP_FORMS: &[(&str, &str)] = &[
    ("NP_003997.1:p.Trp26*", "NP_003997.1:p.Trp26Ter"),
    (
        "NP_003997.1:p.Arg97Profs*23",
        "NP_003997.1:p.Arg97ProfsTer23",
    ),
    (
        "NP_003997.1:p.Gln151Thrfs*9",
        "NP_003997.1:p.Gln151ThrfsTer9",
    ),
    // Extension: the `*` new-stop glyph is accepted and canonicalizes to `Ter`.
    (
        "NP_003997.1:p.Ter327Glnext*?",
        "NP_003997.1:p.Ter327GlnextTer?",
    ),
];

/// Strict mode accepts every `*` stop spelling and canonicalizes it to `Ter`.
#[test]
fn strict_accepts_star_stop_codon_as_ter() {
    for (input, expected) in STAR_STOP_FORMS {
        let parsed = parse_hgvs_with_config(input, ErrorConfig::strict())
            .unwrap_or_else(|e| panic!("strict must accept spec-valid {input:?}: {e}"));
        assert_eq!(parsed.result.to_string(), *expected, "input {input:?}");
    }
}

/// `*` is valid, not a deprecated form, so strict-accepting it emits no
/// W3007/W3009 deprecation warning.
#[test]
fn star_stop_codon_is_not_a_correction() {
    for (input, _) in STAR_STOP_FORMS {
        let parsed = parse_hgvs_with_config(input, ErrorConfig::strict())
            .unwrap_or_else(|e| panic!("strict must accept {input:?}: {e}"));
        assert!(
            parsed.warnings.is_empty(),
            "{input:?} is a valid `*` spelling; expected no warnings, got {:?}",
            parsed.warnings,
        );
    }
}

/// The genuinely-deprecated `X` stop spelling (checklist.md:63 "the `X`
/// should not be used") stays rejected in strict — the reclassification is
/// `*`-only, leaving `X` (W3008 / W3010) unchanged.
#[test]
fn strict_still_rejects_deprecated_x_stop() {
    // Strict rejects the deprecated `X` stop (W3008) and frameshift (W3010)
    // forms. Assert the rejection is specifically the deprecated-`X` diagnostic
    // — not an unrelated parse failure — so this can't silently pass on the
    // wrong error. (Strict surfaces the diagnostic message; the W-code itself is
    // carried on the warning in lenient/silent, exercised elsewhere.)
    for input in ["NP_003997.1:p.Trp26X", "NP_003997.1:p.Arg97ProfsX23"] {
        let err = parse_hgvs_with_config(input, ErrorConfig::strict())
            .expect_err(&format!(
                "strict must still reject deprecated `X` form {input:?}"
            ))
            .to_string();
        assert!(
            err.contains("Deprecated protein notation 'X"),
            "strict rejection of {input:?} must cite the deprecated `X` notation; got: {err}",
        );
    }
}

/// The 3'UTR `*` (`c.*10`) is a *coordinate* marker, a wholly separate meaning
/// from the protein stop-codon `*`. Accepting the stop-codon `*` must not
/// conflate the two: a UTR `*` still parses (in strict) and is preserved
/// verbatim, not rewritten to `Ter`.
#[test]
fn utr_star_coordinate_is_not_conflated_with_stop_codon() {
    let parsed = parse_hgvs_with_config("NM_004006.2:c.*10A>G", ErrorConfig::strict())
        .expect("strict must accept a 3'UTR `*` coordinate");
    assert_eq!(parsed.result.to_string(), "NM_004006.2:c.*10A>G");
}
