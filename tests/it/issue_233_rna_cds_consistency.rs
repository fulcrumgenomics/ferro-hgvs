//! Audit for #81 § E3 — r./c. canonicalization parity.
//!
//! HGVS v21.0 § 3.1 places `r.` and `c.` on the same coordinate axis:
//! same transcript, same position numbering, same structural shape
//! (intervals, sentinels, brackets, predicted wrappers). Only the
//! nucleotide alphabet differs (`A/C/G/T` ↔ `a/c/g/u`) and the prefix
//! (`c.` ↔ `r.`).
//!
//! This file pins, for matched pairs of `c.` and `r.` inputs on the
//! same transcript and the same numerical positions, that the only
//! difference between their `Display` outputs is:
//!
//! 1. the prefix `c.` vs `r.`, and
//! 2. the alphabet (`A/C/G/T` ↔ `a/c/g/u`).
//!
//! The two paths share `LocEdit` for position rendering, so any
//! divergence on numbering or punctuation between them is a bug.
//! The edit rendering routes through `NaEdit::Display` (uppercase)
//! and `NaEdit::to_rna_string` (lowercase) — pinning parity guards
//! against the two formatters drifting apart.

use ferro_hgvs::parse_hgvs;

const ACC: &str = "NM_000088.3";

/// Mechanically translate a `c.`-form `Display` string into the
/// expected `r.`-form by swapping the prefix and lowercasing the
/// nucleotide alphabet. Used as the ground truth for parity checks
/// — the actual `r.` `Display` output must equal this translation.
fn translate_cds_to_rna(s: &str) -> String {
    s.replace(":c.", ":r.")
        .replace('A', "a")
        .replace('C', "c")
        .replace('G', "g")
        .replace('T', "u")
}

/// Parse a `c.` input and its matched `r.` input, then assert that
/// the `r.` `Display` output equals the mechanical c→r translation
/// of the `c.` `Display` output.
#[track_caller]
fn assert_parity(c_input: &str, r_input: &str) {
    let c = parse_hgvs(c_input).unwrap_or_else(|e| panic!("parse {c_input:?}: {e}"));
    let r = parse_hgvs(r_input).unwrap_or_else(|e| panic!("parse {r_input:?}: {e}"));
    let c_out = format!("{}", c);
    let r_out = format!("{}", r);
    let translated = translate_cds_to_rna(&c_out);
    assert_eq!(
        r_out, translated,
        "r./c. parity broken:\n  c input    = {c_input}\n  c output   = {c_out}\n  c→r expect = {translated}\n  r input    = {r_input}\n  r output   = {r_out}"
    );
}

// =============================================================================
// SECTION 1 — Substitution
// =============================================================================

mod substitution {
    use super::*;

    #[test]
    fn exonic_sub_has_parity() {
        assert_parity(&format!("{ACC}:c.123A>G"), &format!("{ACC}:r.123a>g"));
    }

    #[test]
    fn intronic_offset_sub_has_parity() {
        assert_parity(&format!("{ACC}:c.123+5A>G"), &format!("{ACC}:r.123+5a>g"));
        assert_parity(&format!("{ACC}:c.123-5A>G"), &format!("{ACC}:r.123-5a>g"));
    }
}

// =============================================================================
// SECTION 2 — Deletion
// =============================================================================

mod deletion {
    use super::*;

    #[test]
    fn pos_only_has_parity() {
        assert_parity(
            &format!("{ACC}:c.123_125del"),
            &format!("{ACC}:r.123_125del"),
        );
    }

    #[test]
    fn stated_ref_has_parity() {
        assert_parity(
            &format!("{ACC}:c.123_125delACG"),
            &format!("{ACC}:r.123_125delacg"),
        );
    }
}

// =============================================================================
// SECTION 3 — Insertion
// =============================================================================

mod insertion {
    use super::*;

    #[test]
    fn insertion_has_parity() {
        assert_parity(
            &format!("{ACC}:c.123_124insACG"),
            &format!("{ACC}:r.123_124insacg"),
        );
    }
}

// =============================================================================
// SECTION 4 — Duplication
// =============================================================================

mod duplication {
    use super::*;

    #[test]
    fn pos_only_has_parity() {
        assert_parity(
            &format!("{ACC}:c.123_125dup"),
            &format!("{ACC}:r.123_125dup"),
        );
    }

    #[test]
    fn stated_ref_has_parity() {
        assert_parity(
            &format!("{ACC}:c.123_125dupACG"),
            &format!("{ACC}:r.123_125dupacg"),
        );
    }
}

// =============================================================================
// SECTION 5 — Delins
// =============================================================================

mod delins {
    use super::*;

    #[test]
    fn delins_with_stated_ref_has_parity() {
        assert_parity(
            &format!("{ACC}:c.123_125delACGinsTTT"),
            &format!("{ACC}:r.123_125delacginsuuu"),
        );
    }

    #[test]
    fn delins_without_stated_ref_has_parity() {
        assert_parity(
            &format!("{ACC}:c.123_125delinsTTT"),
            &format!("{ACC}:r.123_125delinsuuu"),
        );
    }
}

// =============================================================================
// SECTION 6 — Inversion
// =============================================================================

mod inversion {
    use super::*;

    #[test]
    fn pos_only_has_parity() {
        assert_parity(
            &format!("{ACC}:c.123_125inv"),
            &format!("{ACC}:r.123_125inv"),
        );
    }

    #[test]
    fn stated_ref_has_parity() {
        assert_parity(
            &format!("{ACC}:c.123_125invACG"),
            &format!("{ACC}:r.123_125invacg"),
        );
    }
}

// =============================================================================
// SECTION 7 — Repeat
// =============================================================================

mod repeat {
    use super::*;

    #[test]
    fn repeat_unit_has_parity() {
        assert_parity(&format!("{ACC}:c.123ACG[5]"), &format!("{ACC}:r.123acg[5]"));
    }
}

// =============================================================================
// SECTION 8 — Whole-entity sentinels
// =============================================================================
//
// `=` (identity) and `?` (unknown) are alphabet-free and structurally
// identical between `c.` and `r.`. `0` (no product) is intentionally
// **RNA-specific**: `r.0` ("no RNA produced") has no `c.` counterpart
// because `c.` does not have a "no DNA" concept. Pin both the parity
// for `=`/`?` and the asymmetry for `0`.

mod sentinels {
    use super::*;

    #[test]
    fn identity_has_parity() {
        assert_parity(&format!("{ACC}:c.="), &format!("{ACC}:r.="));
    }

    #[test]
    fn unknown_has_parity() {
        assert_parity(&format!("{ACC}:c.?"), &format!("{ACC}:r.?"));
    }

    /// `r.0` parses; `c.0` does not — this asymmetry is intentional
    /// per HGVS spec semantics ("no RNA produced" is RNA-specific).
    #[test]
    fn no_product_is_rna_only() {
        assert!(parse_hgvs(&format!("{ACC}:r.0")).is_ok(), "r.0 must parse");
        assert!(
            parse_hgvs(&format!("{ACC}:c.0")).is_err(),
            "c.0 must NOT parse (no `c.0` form in HGVS)"
        );
    }
}

// =============================================================================
// SECTION 9 — UTR markers
// =============================================================================
//
// `c.-N` (5'UTR) and `c.*N` (3'UTR) numbering rolls through to `r.`
// unchanged on the position side. Pin parity for both.

mod utr_markers {
    use super::*;

    #[test]
    fn five_prime_utr_has_parity() {
        assert_parity(&format!("{ACC}:c.-12A>G"), &format!("{ACC}:r.-12a>g"));
    }

    #[test]
    fn three_prime_utr_has_parity() {
        assert_parity(&format!("{ACC}:c.*5A>G"), &format!("{ACC}:r.*5a>g"));
    }
}

// =============================================================================
// SECTION 10 — Compound cis allele brackets
// =============================================================================

mod compound_cis {
    use super::*;

    #[test]
    fn two_sub_has_parity() {
        assert_parity(
            &format!("{ACC}:c.[123A>G;125C>T]"),
            &format!("{ACC}:r.[123a>g;125c>u]"),
        );
    }

    #[test]
    fn mixed_edit_shapes_have_parity() {
        assert_parity(
            &format!("{ACC}:c.[123A>G;200_202delACG;300_301insACG]"),
            &format!("{ACC}:r.[123a>g;200_202delacg;300_301insacg]"),
        );
    }
}
