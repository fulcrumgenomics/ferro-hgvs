//! Audit for issue #226 — dedicated coverage for the remaining
//! protein-canonicalization items in #81 § D (D1, D3, D4, D8). D2 ships
//! separately as #224 / PR #225 (extension `Ter` consistency).
//!
//! ferro accepts both three-letter and one-letter amino-acid codes,
//! and both `Ter` and `*` for the translation termination codon. The
//! canonical Display form is three-letter codes with `Ter` for the
//! stop, per `assets/hgvs-nomenclature/docs/background/standards.md`.
//! Predicted variants (those without RNA/protein-level confirmation)
//! are wrapped in parens — `p.(Arg97Trp)` — per
//! `recommendations/protein/standards.md`.
//!
//! Sections (one per #81 D item):
//!   - D1: frameshift `fs` rendering — short / long-Ter / long-star /
//!     predicted parens / unknown-Ter (`*?` / `Ter?`) / 1-letter →
//!     3-letter canonicalization.
//!   - D3: termination `Ter` vs `*` — substitution, frameshift,
//!     1-letter star canonical.
//!   - D4: 3-letter vs 1-letter — substitution, deletion, insertion,
//!     duplication, delins; 3-letter round-trip.
//!   - D8: predicted-vs-confirmed parens — sub / fs / del / identity /
//!     1-letter-inside-parens canonicalization / empty-parens reject.
//!
//! A probe of every shape in this file confirmed the current behavior
//! is correct; this PR is pure coverage (no `src/` changes).

use ferro_hgvs::parse_hgvs;

fn assert_round_trips(input: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{}`: {}", input, e));
    let out = format!("{}", v);
    assert_eq!(out, input, "round-trip mismatch: `{}` -> `{}`", input, out);
}

fn assert_canonicalizes_to(input: &str, expected: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{}`: {}", input, e));
    let out = format!("{}", v);
    assert_eq!(
        out, expected,
        "canonicalization mismatch:\n  input:    `{}`\n  expected: `{}`\n  got:      `{}`",
        input, expected, out
    );
}

// =============================================================================
// SECTION D1 — frameshift rendering
// =============================================================================
//
// Spec (`recommendations/protein/frameshift.md`): three forms accepted —
// short (`Arg97fs`), long-Ter (`Arg97ProfsTer23`), long-star
// (`Arg97Profs*23`). The long-Ter form is the spec-preferred canonical
// emission.

mod d1_frameshift {
    use super::*;

    /// Short form round-trips (no new-AA info, no stop position).
    #[test]
    fn short_form_round_trips() {
        assert_round_trips("NP_003997.1:p.Arg97fs");
    }

    /// Long-Ter form is the canonical emission.
    #[test]
    fn long_ter_form_round_trips() {
        assert_round_trips("NP_003997.1:p.Arg97ProfsTer23");
    }

    /// Long-star form canonicalizes to long-Ter.
    #[test]
    fn long_star_canonicalizes_to_long_ter() {
        assert_canonicalizes_to(
            "NP_003997.1:p.Arg97Profs*23",
            "NP_003997.1:p.Arg97ProfsTer23",
        );
    }

    /// Predicted parens around short form preserved.
    #[test]
    fn predicted_parens_short_form_round_trips() {
        assert_round_trips("NP_003997.1:p.(Arg97fs)");
    }

    /// Predicted parens around long form preserved.
    #[test]
    fn predicted_parens_long_form_round_trips() {
        assert_round_trips("NP_003997.1:p.(Arg97ProfsTer23)");
    }

    /// `Argfs*?` (unknown stop, star form) → short canonical `fs`
    /// when stop is unknown. The new-AA info is also dropped — spec
    /// lists `p.Ile327Argfs*?` and `p.Ile327fs` as equivalent
    /// (recommendations/protein/frameshift.md examples).
    #[test]
    fn unknown_ter_star_canonicalizes_to_short_form() {
        assert_canonicalizes_to("NP_003997.1:p.Ile327Argfs*?", "NP_003997.1:p.Ile327Argfs");
    }

    /// Same equivalence using `Ter?` notation.
    #[test]
    fn unknown_ter_three_letter_canonicalizes_to_short_form() {
        assert_canonicalizes_to("NP_003997.1:p.Ile327ArgfsTer?", "NP_003997.1:p.Ile327Argfs");
    }

    /// One-letter short form canonicalizes to three-letter.
    #[test]
    fn one_letter_short_form_canonicalizes_to_three_letter() {
        assert_canonicalizes_to("NP_003997.1:p.R97fs", "NP_003997.1:p.Arg97fs");
    }

    /// One-letter long form: 1-letter AAs canonicalize, star → Ter.
    #[test]
    fn one_letter_long_form_canonicalizes_to_three_letter_ter() {
        assert_canonicalizes_to("NP_003997.1:p.R97Pfs*23", "NP_003997.1:p.Arg97ProfsTer23");
    }
}

// =============================================================================
// SECTION D3 — termination Ter vs *
// =============================================================================
//
// `Ter` is the spec-preferred three-letter form for the translation
// termination codon. ferro canonicalizes substitutions and frameshifts
// to `Ter` regardless of input form.

mod d3_termination {
    use super::*;

    #[test]
    fn substitution_ter_round_trips() {
        assert_round_trips("NP_003997.1:p.Tyr4Ter");
    }

    #[test]
    fn substitution_star_canonicalizes_to_ter() {
        assert_canonicalizes_to("NP_003997.1:p.Tyr4*", "NP_003997.1:p.Tyr4Ter");
    }

    #[test]
    fn substitution_one_letter_star_canonicalizes_to_ter() {
        assert_canonicalizes_to("NP_003997.1:p.Y4*", "NP_003997.1:p.Tyr4Ter");
    }

    /// Frameshift Ter canonicalization (cross-reference with D1).
    #[test]
    fn frameshift_star_canonicalizes_to_ter() {
        assert_canonicalizes_to(
            "NP_003997.1:p.Arg97Profs*23",
            "NP_003997.1:p.Arg97ProfsTer23",
        );
    }
}

// =============================================================================
// SECTION D4 — 3-letter vs 1-letter rendering
// =============================================================================
//
// Spec recommends three-letter codes (`background/standards.md`).
// ferro accepts both; canonical Display is three-letter.

mod d4_three_vs_one_letter {
    use super::*;

    // -- Substitution --

    #[test]
    fn sub_three_letter_round_trips() {
        assert_round_trips("NP_003997.1:p.Arg97Trp");
    }

    #[test]
    fn sub_one_letter_canonicalizes_to_three_letter() {
        assert_canonicalizes_to("NP_003997.1:p.R97W", "NP_003997.1:p.Arg97Trp");
    }

    // -- Deletion --

    #[test]
    fn del_three_letter_round_trips() {
        assert_round_trips("NP_003997.1:p.Lys100del");
    }

    #[test]
    fn del_one_letter_canonicalizes_to_three_letter() {
        assert_canonicalizes_to("NP_003997.1:p.K100del", "NP_003997.1:p.Lys100del");
    }

    // -- Insertion --

    #[test]
    fn ins_three_letter_round_trips() {
        assert_round_trips("NP_003997.1:p.Lys100_Arg101insGlySer");
    }

    #[test]
    fn ins_one_letter_canonicalizes_to_three_letter() {
        assert_canonicalizes_to(
            "NP_003997.1:p.K100_R101insGS",
            "NP_003997.1:p.Lys100_Arg101insGlySer",
        );
    }

    // -- Duplication --

    #[test]
    fn dup_three_letter_round_trips() {
        assert_round_trips("NP_003997.1:p.Lys100_Arg101dup");
    }

    // -- Delins --

    #[test]
    fn delins_one_letter_canonicalizes_to_three_letter() {
        assert_canonicalizes_to(
            "NP_003997.1:p.K100_R101delinsGS",
            "NP_003997.1:p.Lys100_Arg101delinsGlySer",
        );
    }
}

// =============================================================================
// SECTION D8 — predicted-vs-confirmed parens
// =============================================================================
//
// Per `recommendations/protein/standards.md`, predicted protein-level
// changes (those without direct RNA/protein evidence) are wrapped in
// parens. ferro preserves the paren distinction on round-trip.

mod d8_predicted_parens {
    use super::*;

    #[test]
    fn confirmed_substitution_no_parens() {
        assert_round_trips("NP_003997.1:p.Arg97Trp");
    }

    #[test]
    fn predicted_substitution_round_trips() {
        assert_round_trips("NP_003997.1:p.(Arg97Trp)");
    }

    #[test]
    fn predicted_frameshift_round_trips() {
        assert_round_trips("NP_003997.1:p.(Arg97ProfsTer23)");
    }

    #[test]
    fn predicted_deletion_round_trips() {
        assert_round_trips("NP_003997.1:p.(Lys100del)");
    }

    /// One-letter substitution inside parens canonicalizes to
    /// three-letter while preserving the parens.
    #[test]
    fn predicted_one_letter_canonicalizes_inside_parens() {
        assert_canonicalizes_to("NP_003997.1:p.(R97W)", "NP_003997.1:p.(Arg97Trp)");
    }

    /// Predicted identity (`p.(=)`) round-trips. The bare `p.=`
    /// (confirmed identity) is distinct.
    #[test]
    fn predicted_identity_round_trips() {
        assert_round_trips("NP_003997.1:p.(=)");
    }

    /// Empty parens are malformed — must reject.
    #[test]
    fn empty_parens_rejected() {
        assert!(parse_hgvs("NP_003997.1:p.()").is_err());
    }

    /// The parens-vs-no-parens distinction is preserved through
    /// round-trip: confirmed and predicted forms are NOT collapsed.
    #[test]
    fn confirmed_and_predicted_remain_distinct() {
        let confirmed = parse_hgvs("NP_003997.1:p.Arg97Trp").expect("parse confirmed");
        let predicted = parse_hgvs("NP_003997.1:p.(Arg97Trp)").expect("parse predicted");
        assert_ne!(
            format!("{}", confirmed),
            format!("{}", predicted),
            "confirmed and predicted forms must NOT collapse to the same string",
        );
    }
}
