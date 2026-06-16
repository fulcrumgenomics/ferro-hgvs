//! Pin the `r.` thymine policy (#81 § E follow-up to #232, closes #282).
//!
//! HGVS v21.0 RNA nomenclature uses the alphabet `a/c/g/u`; `t` is not
//! a valid RNA base. Prior to this issue, ferro accepted `t` leniently
//! and rendered it as lowercase `t` on `Display`. PR #232 pinned that
//! observed behavior but left the policy decision unresolved. PR #293
//! (issue #276) made `Display` emit `u` regardless of stored input,
//! closing half of the loop on the output side.
//!
//! This file pins the **input-side** policy: a `t` byte appearing as
//! an RNA base inside an `r.` description is **canonicalized to `u`**
//! at preprocess time, and a soft-validation warning **W3020
//! (RnaThymineCanonicalized)** is emitted in lenient mode. Strict
//! mode rejects with the same code visible in the error diagnostic.
//!
//! Boundary expectations:
//!
//! - Coverage spans every RNA edit shape where a base byte appears:
//!   substitution (ref/alt), deletion (stated ref), insertion (alt),
//!   duplication (stated ref), inversion (stated ref), delins
//!   (stated ref and alt), and repeat unit.
//! - Canonical `r.` input with `u` is unaffected (no W3020).
//! - `t` in `c.` / `g.` / `n.` / `m.` contexts is canonical DNA and
//!   must remain unaffected — W3020 fires only for `r.` descriptions.
//! - Strict mode propagates the policy: the diagnostic must surface
//!   W3020 (the parsed value never reaches the caller in strict mode).

use ferro_hgvs::error_handling::{ErrorConfig, ErrorMode, ErrorType};
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::parse_hgvs;

const ACC: &str = "NM_000088.3";

/// Parse with lenient mode; assert success, W3020 fired, and `Display`
/// of the parsed variant equals `expected`.
#[track_caller]
fn assert_canonicalizes_with_w3019(input: &str, expected: &str) {
    let parsed = parse_hgvs_with_config(input, ErrorConfig::lenient())
        .unwrap_or_else(|e| panic!("lenient parse {input:?} failed: {e}"));
    let warnings = &parsed.warnings;
    assert!(
        warnings.iter().any(|w| w.error_type.code() == "W3020"),
        "expected W3020 on {input:?}, got warnings: {:?}",
        warnings
            .iter()
            .map(|w| w.error_type.code())
            .collect::<Vec<_>>()
    );
    let out = format!("{}", parsed.result);
    assert_eq!(out, expected, "input={input:?}");
}

/// Parse with lenient mode; assert success and **no** W3020 fired.
#[track_caller]
fn assert_no_w3019(input: &str) {
    let parsed = parse_hgvs_with_config(input, ErrorConfig::lenient())
        .unwrap_or_else(|e| panic!("lenient parse {input:?} failed: {e}"));
    let warnings = &parsed.warnings;
    assert!(
        !warnings.iter().any(|w| w.error_type.code() == "W3020"),
        "unexpected W3020 on {input:?}: {:?}",
        warnings
            .iter()
            .map(|w| w.error_type.code())
            .collect::<Vec<_>>()
    );
}

/// Strict mode must reject; the error's display must mention W3020.
#[track_caller]
fn assert_strict_rejects_with_w3019(input: &str) {
    let result = parse_hgvs_with_config(input, ErrorConfig::strict());
    let err = match result {
        Ok(_) => panic!("strict parse {input:?} succeeded, expected rejection"),
        Err(e) => e,
    };
    let msg = format!("{err}");
    assert!(
        msg.contains("W3020") || msg.to_lowercase().contains("thymine"),
        "strict rejection for {input:?} did not surface W3020 / thymine: {msg}"
    );
}

/// Parse with lenient mode; assert success, the given warning code fired, and
/// the rewritten `Display` output equals `expected` (pins the seq-strip rewrite).
#[track_caller]
fn assert_lenient_strips_to(input: &str, code: &str, expected: &str) {
    let parsed = parse_hgvs_with_config(input, ErrorConfig::lenient())
        .unwrap_or_else(|e| panic!("lenient parse {input:?} failed: {e}"));
    let warnings = &parsed.warnings;
    assert!(
        warnings.iter().any(|w| w.error_type.code() == code),
        "expected {code} on {input:?}, got warnings: {:?}",
        warnings
            .iter()
            .map(|w| w.error_type.code())
            .collect::<Vec<_>>()
    );
    assert_eq!(
        format!("{}", parsed.result),
        expected,
        "rewrite mismatch for {input:?}"
    );
}

/// Strict mode must reject; the error display must mention the given code or keyword.
#[track_caller]
fn assert_strict_rejects_with_code(input: &str, code: &str, keyword: &str) {
    let result = parse_hgvs_with_config(input, ErrorConfig::strict());
    let err = match result {
        Ok(_) => panic!("strict parse {input:?} succeeded, expected rejection"),
        Err(e) => e,
    };
    let msg = format!("{err}");
    assert!(
        msg.contains(code) || msg.to_lowercase().contains(keyword),
        "strict rejection for {input:?} did not surface {code} / {keyword}: {msg}"
    );
}

// =============================================================================
// SECTION 1 — Substitution (sub)
// =============================================================================
//
// `r.<pos><ref>><alt>` — `t` in either ref or alt position triggers
// canonicalization to `u` with W3020.

mod substitution {
    use super::*;

    #[test]
    fn alt_t_canonicalizes_to_u_with_w3019() {
        assert_canonicalizes_with_w3019(&format!("{ACC}:r.123a>t"), &format!("{ACC}:r.123a>u"));
    }

    #[test]
    fn ref_t_canonicalizes_to_u_with_w3019() {
        assert_canonicalizes_with_w3019(&format!("{ACC}:r.123t>g"), &format!("{ACC}:r.123u>g"));
    }

    #[test]
    fn uppercase_alt_t_canonicalizes_to_u_with_w3019() {
        assert_canonicalizes_with_w3019(&format!("{ACC}:r.123A>T"), &format!("{ACC}:r.123a>u"));
    }

    #[test]
    fn canonical_u_input_does_not_fire_w3019() {
        assert_no_w3019(&format!("{ACC}:r.123a>u"));
        assert_no_w3019(&format!("{ACC}:r.123c>u"));
    }

    #[test]
    fn strict_mode_rejects_t_alt() {
        assert_strict_rejects_with_w3019(&format!("{ACC}:r.123a>t"));
    }
}

// =============================================================================
// SECTION 2 — Deletion (del) with stated reference bases
// =============================================================================

mod deletion {
    use super::*;

    #[test]
    fn stated_ref_with_t_canonicalizes_to_u_with_w3019() {
        // W3025 (DelExplicitSeq) fires before W3020 (RnaThymineCanonicalized):
        // the explicit sequence is stripped in Phase 13d before the thymine
        // check in Phase 13e can see the `t` bytes. W3025 supersedes W3020.
        // Pin the rewritten output too: the explicit seq is dropped entirely.
        assert_lenient_strips_to(
            &format!("{ACC}:r.123_125delaut"),
            "W3025",
            &format!("{ACC}:r.123_125del"),
        );
    }

    #[test]
    fn uppercase_stated_ref_with_t_canonicalizes() {
        // Same ordering as above: W3025 strips the seq (after case-folding),
        // so W3020 (thymine) never fires. The surfacing warning is W3025.
        assert_lenient_strips_to(
            &format!("{ACC}:r.123_125delAUT"),
            "W3025",
            &format!("{ACC}:r.123_125del"),
        );
    }

    #[test]
    // The explicit-sequence form fires W3025 (DelExplicitSeq), but this test
    // specifically validates that canonical RNA bases (`u`, not `t`) do not
    // trigger the thymine warning W3020.
    fn lowercase_canonical_aug_no_thymine_warning() {
        assert_no_w3019(&format!("{ACC}:r.123_125delaug"));
    }

    #[test]
    fn strict_mode_rejects_t_in_del_ref() {
        // Strict mode rejects with W3025 (DelExplicitSeq) because the explicit
        // sequence check in Phase 13d fires before the thymine check (Phase 13e).
        assert_strict_rejects_with_code(
            &format!("{ACC}:r.123_125delaut"),
            "W3025",
            "explicit sequence",
        );
    }
}

// =============================================================================
// SECTION 3 — Insertion (ins)
// =============================================================================

mod insertion {
    use super::*;

    #[test]
    fn inserted_seq_with_t_canonicalizes_to_u_with_w3019() {
        assert_canonicalizes_with_w3019(
            &format!("{ACC}:r.123_124insaut"),
            &format!("{ACC}:r.123_124insauu"),
        );
    }

    #[test]
    fn uppercase_inserted_t_canonicalizes() {
        assert_canonicalizes_with_w3019(
            &format!("{ACC}:r.123_124insAUT"),
            &format!("{ACC}:r.123_124insauu"),
        );
    }

    #[test]
    fn canonical_ins_no_warning() {
        assert_no_w3019(&format!("{ACC}:r.123_124insauc"));
    }
}

// =============================================================================
// SECTION 4 — Duplication (dup) with stated reference bases
// =============================================================================

mod duplication {
    use super::*;

    #[test]
    fn stated_ref_with_t_canonicalizes_to_u_with_w3019() {
        // W3024 (DupExplicitSeq) fires before W3020 (RnaThymineCanonicalized):
        // the explicit sequence is stripped in Phase 13c before the thymine
        // check can see the `t` bytes. W3024 supersedes W3020.
        // Pin the rewritten output too: the explicit seq is dropped entirely.
        assert_lenient_strips_to(
            &format!("{ACC}:r.123_125dupaut"),
            "W3024",
            &format!("{ACC}:r.123_125dup"),
        );
    }

    #[test]
    fn uppercase_dup_t_canonicalizes() {
        // Same ordering: W3024 strips the seq (after case-folding),
        // so W3020 (thymine) never fires. The surfacing warning is W3024.
        assert_lenient_strips_to(
            &format!("{ACC}:r.123_125dupAUT"),
            "W3024",
            &format!("{ACC}:r.123_125dup"),
        );
    }
}

// =============================================================================
// SECTION 5 — Inversion (inv) with stated reference bases
// =============================================================================

mod inversion {
    use super::*;

    #[test]
    fn stated_ref_with_t_canonicalizes_to_u_with_w3019() {
        assert_canonicalizes_with_w3019(
            &format!("{ACC}:r.123_125invaut"),
            &format!("{ACC}:r.123_125invauu"),
        );
    }
}

// =============================================================================
// SECTION 6 — Delins
// =============================================================================

mod delins {
    use super::*;

    #[test]
    fn alt_with_t_canonicalizes_to_u_with_w3019() {
        assert_canonicalizes_with_w3019(
            &format!("{ACC}:r.123_125delinsaut"),
            &format!("{ACC}:r.123_125delinsauu"),
        );
    }

    #[test]
    fn stated_ref_with_t_and_canonical_alt_canonicalizes() {
        // Only the lone `t` is rewritten — surrounding canonical RNA
        // bases (`a`, `u`) are not perturbed.
        assert_canonicalizes_with_w3019(
            &format!("{ACC}:r.123_125delautinsacg"),
            &format!("{ACC}:r.123_125delauuinsacg"),
        );
    }

    #[test]
    fn uppercase_delins_t_canonicalizes() {
        assert_canonicalizes_with_w3019(
            &format!("{ACC}:r.123_125delinsAUT"),
            &format!("{ACC}:r.123_125delinsauu"),
        );
    }
}

// =============================================================================
// SECTION 7 — Repeat unit
// =============================================================================
//
// `r.<pos><unit>[N]` — `t` inside the repeat unit triggers W3020.

mod repeat {
    use super::*;

    #[test]
    fn unit_with_t_canonicalizes_to_u_with_w3019() {
        assert_canonicalizes_with_w3019(
            &format!("{ACC}:r.123aut[5]"),
            &format!("{ACC}:r.123auu[5]"),
        );
    }

    #[test]
    fn uppercase_unit_t_canonicalizes() {
        assert_canonicalizes_with_w3019(
            &format!("{ACC}:r.123AUT[5]"),
            &format!("{ACC}:r.123auu[5]"),
        );
    }

    #[test]
    fn canonical_repeat_unit_no_warning() {
        assert_no_w3019(&format!("{ACC}:r.123aug[5]"));
    }
}

// =============================================================================
// SECTION 8 — Boundary: non-RNA contexts must NOT fire W3020
// =============================================================================
//
// `t` is canonical DNA in c./g./n./m. contexts. The W3020 detector
// must be scoped to `r.` descriptions and leave other coordinate
// types alone.

mod non_rna_unaffected {
    use super::*;

    #[test]
    fn cds_t_does_not_fire_w3019() {
        assert_no_w3019(&format!("{ACC}:c.123A>T"));
        assert_no_w3019(&format!("{ACC}:c.123_125delAUT"));
    }

    #[test]
    fn genome_t_does_not_fire_w3019() {
        assert_no_w3019("NC_000001.11:g.12345A>T");
    }

    #[test]
    fn noncoding_t_does_not_fire_w3019() {
        assert_no_w3019("NR_046018.2:n.100A>T");
    }

    #[test]
    fn mito_t_does_not_fire_w3019() {
        assert_no_w3019("NC_012920.1:m.3243A>T");
    }
}

// =============================================================================
// SECTION 9 — Compound allele (cis brackets)
// =============================================================================
//
// Bracketed alleles must propagate the canonicalization to each
// member edit; W3020 fires for every offending base.

mod compound_cis {
    use super::*;

    #[test]
    fn bracketed_member_t_canonicalizes() {
        assert_canonicalizes_with_w3019(
            &format!("{ACC}:r.[123a>t;125c>u]"),
            &format!("{ACC}:r.[123a>u;125c>u]"),
        );
    }

    /// The default `parse_hgvs` path does not run the preprocessor, so
    /// W3020 does not fire — but the parser still succeeds. The W3020
    /// canonicalization is observable through the lenient / silent
    /// paths (`parse_hgvs_with_config`). This test pins that the default
    /// `parse_hgvs` path remains tolerant of `t` (no rejection), so
    /// existing callers that bypass the preprocessor continue to work.
    /// Display-side T→u canonicalization is tracked separately in #276.
    #[test]
    fn parse_hgvs_default_path_still_accepts_t() {
        let v = parse_hgvs(&format!("{ACC}:r.123a>t"));
        assert!(v.is_ok(), "default parse_hgvs should still accept r.123a>t");
    }
}

// =============================================================================
// SECTION 10 — Mode behavior pinning
// =============================================================================

mod modes {
    use super::*;

    #[test]
    fn silent_mode_canonicalizes_without_warning() {
        let parsed =
            parse_hgvs_with_config(&format!("{ACC}:r.123a>t"), ErrorConfig::silent()).unwrap();
        assert!(!parsed
            .warnings
            .iter()
            .any(|w| w.error_type.code() == "W3020"));
        assert_eq!(parsed.preprocessed_input, format!("{ACC}:r.123a>u"));
    }

    #[test]
    fn strict_mode_default() {
        let config = ErrorConfig::new(ErrorMode::Strict);
        let result = parse_hgvs_with_config(&format!("{ACC}:r.123a>t"), config);
        assert!(result.is_err(), "strict should reject r.123a>t");
    }

    #[test]
    fn error_type_code_is_w3019() {
        assert_eq!(ErrorType::RnaThymineCanonicalized.code(), "W3020");
    }
}
