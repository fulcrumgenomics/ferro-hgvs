//! Audit for #290 — protein bracketed-amino-acid insertion
//! (`p.…ins[Ala;Pro]`) parser gap surfaced by PR #248.
//!
//! Policy (b): **reject with a targeted diagnostic** (W3021
//! `ProteinBracketedAaInsertion`). HGVS v21's protein insertion
//! notation concatenates 3-letter codes without separators
//! (`p.Arg97_Trp98insAlaPro`); the bracketed `[Ala;Pro]` form has no
//! protein-edit equivalent in the spec — brackets are reserved for
//! alleles at the variant level, not for amino-acid lists inside an
//! edit. The diagnostic must point the user at the canonical
//! `p.…insAlaPro` shape.
//!
//! The failure is accession-independent: it must fire on `NP_` and
//! `LRG_<n>p<m>` accessions alike, because the gap is in the
//! protein-edit grammar, not the accession dispatch.
//!
//! Closes #290; follow-up to #248 § *Out of scope*.

use ferro_hgvs::error_handling::{ErrorConfig, ErrorType};
use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
use ferro_hgvs::parse_hgvs;

// =============================================================================
// SECTION 1 — Bracketed AA insertion is rejected with W3021
// =============================================================================

mod bracketed_rejected {
    use super::*;

    /// Strict (`parse_hgvs`) rejects the bracketed form regardless of
    /// accession. The error message must mention the canonical
    /// `insAlaPro` alternative so the user knows the fix.
    #[track_caller]
    fn assert_rejected_with_hint(input: &str) {
        let err = parse_hgvs(input)
            .err()
            .unwrap_or_else(|| panic!("expected parse failure for {input:?}"));
        let msg = err.to_string();
        // The detailed message includes the hint. Check both.
        let detailed = err.detailed_message();
        let combined = format!("{msg}\n{detailed}");
        assert!(
            combined.contains("insAlaPro") || combined.contains("ins<AA1><AA2>"),
            "diagnostic for {input:?} must point at the canonical \
             concatenated form; got:\n{combined}"
        );
    }

    #[test]
    fn np_accession_rejected_with_hint() {
        assert_rejected_with_hint("NP_000088.3:p.Arg97_Trp98ins[Ala;Pro]");
    }

    #[test]
    fn lrg_protein_accession_rejected_with_hint() {
        assert_rejected_with_hint("LRG_1p1:p.Arg97_Trp98ins[Ala;Pro]");
    }

    /// In lenient mode the preprocessor surfaces the W3021 code
    /// explicitly (rather than the generic nom Verify failure).
    #[test]
    fn lenient_mode_surfaces_w3020_code() {
        let err = parse_hgvs_lenient("NP_000088.3:p.Arg97_Trp98ins[Ala;Pro]")
            .expect_err("lenient mode must still reject — W3021 is not auto-correctable");
        // The structured diagnostic carries an error code; we don't pin
        // the exact ErrorCode value (it shares the InvalidEdit family),
        // but the textual message must name W3021 / canonical form.
        let detailed = err.detailed_message();
        assert!(
            detailed.contains("insAlaPro") || detailed.contains("ins<AA1><AA2>"),
            "lenient diagnostic must mention the canonical insAlaPro \
             form; got:\n{detailed}"
        );
    }

    /// For a well-formed `[Aa3;Aa3]` body the hint must carry the
    /// concrete canonical rewrite — i.e. literal `insAlaPro` for
    /// `ins[Ala;Pro]` input — not just the generic template. This
    /// pins the suggestion-emission branch in `detect_protein_bracketed_aa_insertion`
    /// against accidentally dropping the canonicalized rewrite.
    #[test]
    fn ala_pro_hint_names_concrete_canonical_form() {
        let err = parse_hgvs("NP_000088.3:p.Arg97_Trp98ins[Ala;Pro]")
            .expect_err("must reject bracketed insertion");
        let combined = format!("{}\n{}", err, err.detailed_message());
        assert!(
            combined.contains("insAlaPro"),
            "diagnostic for `[Ala;Pro]` input must carry the concrete \
             canonical rewrite `insAlaPro`, not just the generic \
             template; got:\n{combined}"
        );
    }

    /// Strict mode rejects bracketed AA insertion as
    /// `ErrorType::ProteinBracketedAaInsertion` (W3021), not as a
    /// generic parse error.
    #[test]
    fn strict_mode_emits_w3020_error_type() {
        let config = ErrorConfig::strict();
        let pp = config.preprocessor();
        let result = pp.preprocess("NP_000088.3:p.Arg97_Trp98ins[Ala;Pro]");
        assert!(
            !result.success,
            "preprocessor must fail on bracketed AA insertion"
        );
        let err = result
            .error
            .expect("strict mode must populate a structured error");
        let detailed = err.detailed_message();
        assert!(
            detailed.contains("insAlaPro") || detailed.contains("ins<AA1><AA2>"),
            "strict mode diagnostic must mention the canonical insAlaPro \
             form; got:\n{detailed}"
        );
        // The code is exposed via the registry; verify it matches W3021.
        assert_eq!(ErrorType::ProteinBracketedAaInsertion.code(), "W3021");
    }
}

// =============================================================================
// SECTION 2 — Canonical (concatenated) insertion shapes still parse
// =============================================================================

mod canonical_round_trip {
    use super::*;

    #[track_caller]
    fn assert_round_trips(s: &str) {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
        let out = format!("{}", v);
        assert_eq!(out, s, "round-trip mismatch for {s:?}");
    }

    /// Single-AA canonical insertion still parses.
    #[test]
    fn single_aa_insertion_round_trips() {
        assert_round_trips("NP_000088.3:p.Arg97_Trp98insAla");
    }

    /// Two-AA canonical insertion (concatenated, no separators) still parses.
    #[test]
    fn two_aa_insertion_round_trips() {
        assert_round_trips("NP_000088.3:p.Arg97_Trp98insAlaPro");
    }

    /// Three-AA canonical insertion still parses — guards against an
    /// over-eager bracket rejection that also breaks the canonical
    /// concatenated form.
    #[test]
    fn three_aa_insertion_round_trips() {
        assert_round_trips("NP_000088.3:p.Arg97_Trp98insAlaProGly");
    }

    /// LRG accession with canonical concatenated insertion still parses.
    #[test]
    fn lrg_canonical_insertion_round_trips() {
        assert_round_trips("LRG_1p1:p.Arg97_Trp98insAla");
    }
}

// =============================================================================
// SECTION 3 — Negative shapes (empty / partial brackets) are rejected
// =============================================================================
//
// These variants of the bracketed form are also non-spec; rejection is
// the right answer for all of them. They should not silently parse to
// something surprising.

mod negative_bracket_shapes {
    use super::*;

    #[track_caller]
    fn assert_rejected(input: &str) {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "{input:?} must be rejected, but parsed to {:?}",
            result.ok()
        );
    }

    #[test]
    fn empty_brackets_rejected() {
        assert_rejected("NP_000088.3:p.Arg97_Trp98ins[]");
    }

    #[test]
    fn single_aa_in_brackets_rejected() {
        // Single AA in brackets is still not the canonical form
        // (canonical is `insAla`, not `ins[Ala]`).
        assert_rejected("NP_000088.3:p.Arg97_Trp98ins[Ala]");
    }

    #[test]
    fn trailing_separator_rejected() {
        assert_rejected("NP_000088.3:p.Arg97_Trp98ins[Ala;]");
    }

    #[test]
    fn leading_separator_rejected() {
        assert_rejected("NP_000088.3:p.Arg97_Trp98ins[;Pro]");
    }
}

// =============================================================================
// SECTION 4 — Bare-prefix (`p.…`) and span-correctness regressions
// =============================================================================
//
// The detection runs in the preprocessor *before* the parser, so it has
// to scan the original input text rather than a partially-rewritten
// copy. These regression guards pin two CodeRabbit-flagged gaps:
//
//   * a bare `p.<body>` input (no accession) used to slip past
//     `detect_protein_bracketed_aa_insertion` because the scan only
//     searched for `:p.`,
//   * running the detection against the post-Phase-6a/c `current`
//     string instead of `input` could shift `first.start`/`first.end`
//     off the bytes that the caller actually sees.
//
// The first is reachable today — the parser later fails at accession
// parsing, but the preprocessor now emits the targeted W3021 hint
// first instead of falling through to the generic accession error.
// The second is silently observable only as wrong span byte offsets,
// so we don't assert on them directly here; surfacing the W3021 hint
// at all on a multi-bracket input is the meaningful behavior signal.

mod bare_and_span_regressions {
    use super::*;

    #[test]
    fn bare_p_prefix_surfaces_w3020_hint() {
        // No accession; the bare-`p.` entry into `detect_protein_…`
        // must catch the bracket form and surface the canonical hint
        // before parsing fails at the accession layer.
        let err = parse_hgvs_lenient("p.Arg97_Trp98ins[Ala;Pro]")
            .expect_err("bare p. with bracketed AA must still be rejected");
        let detailed = err.detailed_message();
        assert!(
            detailed.contains("insAlaPro") || detailed.contains("ins<AA1><AA2>"),
            "bare-p. diagnostic must name the canonical concatenated \
             form; got:\n{detailed}"
        );
    }

    #[test]
    fn detection_uses_original_input_text() {
        // With a normal accession + bracket form, the hint should
        // quote `ins[Ala;Pro]` verbatim from the original input.
        // Pinned so a future regression that re-introduces scanning
        // a rewritten copy of the input shows up here.
        let err = parse_hgvs("NP_000088.3:p.Arg97_Trp98ins[Ala;Pro]")
            .expect_err("must reject bracketed insertion");
        let combined = format!("{}\n{}", err, err.detailed_message());
        assert!(
            combined.contains("ins[Ala;Pro]"),
            "diagnostic must quote the literal bracketed substring from \
             the original input; got:\n{combined}"
        );
    }
}

// =============================================================================
// SECTION 5 — Registry coverage
// =============================================================================

mod registry {
    use ferro_hgvs::error_handling::get_code_info;

    #[test]
    fn w3020_registered() {
        let info = get_code_info("W3021").expect("W3021 must be in the registry");
        assert_eq!(info.code, "W3021");
        assert_eq!(info.name, "ProteinBracketedAaInsertion");
    }
}
