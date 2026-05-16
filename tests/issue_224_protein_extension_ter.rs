//! Audit for issue #224 — protein extension Display emits `extTerN`
//! (not `ext*N`) for C-terminal extensions, mirroring the `Ter`
//! canonicalization that Frameshift already uses.
//!
//! Spec: `assets/hgvs-nomenclature/docs/background/standards.md` —
//! three-letter codes (including `Ter`) are preferred; `*` is an
//! accepted alternative. ferro's canonical Display form is the
//! three-letter form for substitutions and frameshifts; this audit
//! pins the same canonicalization for extensions.
//!
//! Closes the D2 portion of #81 § D (also reaches into D3 since the
//! `Ter` canonicalization is a shared concern).

use ferro_hgvs::parse_hgvs;

fn assert_canonicalizes_to(input: &str, expected: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{}`: {}", input, e));
    let out = format!("{}", v);
    assert_eq!(
        out, expected,
        "canonicalization mismatch:\n  input:    `{}`\n  expected: `{}`\n  got:      `{}`",
        input, expected, out
    );
}

fn assert_round_trips(input: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{}`: {}", input, e));
    let out = format!("{}", v);
    assert_eq!(out, input, "round-trip mismatch: `{}` -> `{}`", input, out);
}

// =============================================================================
// SECTION 1 — C-terminal extension Ter canonicalization
// =============================================================================

mod c_terminal_extension {
    use super::*;

    /// Three-letter `Ter` input round-trips (canonical form).
    #[test]
    fn ter_ter_round_trips() {
        assert_round_trips("NP_003997.1:p.Ter315TyrextTer10");
    }

    /// Three-letter `Ter` start position + `*` extension count must
    /// canonicalize to `Ter` extension count.
    #[test]
    fn ter_star_canonicalizes_to_ter_ter() {
        assert_canonicalizes_to(
            "NP_003997.1:p.Ter315Tyrext*10",
            "NP_003997.1:p.Ter315TyrextTer10",
        );
    }

    /// `*` start position + `*` extension count: both canonicalize.
    #[test]
    fn star_star_canonicalizes_to_ter_ter() {
        assert_canonicalizes_to(
            "NP_003997.1:p.*315Tyrext*10",
            "NP_003997.1:p.Ter315TyrextTer10",
        );
    }

    /// Unknown extension count `*?` → `Ter?` after canonicalization.
    #[test]
    fn unknown_count_canonicalizes_star_to_ter() {
        assert_canonicalizes_to(
            "NP_003997.1:p.Ter315Tyrext*?",
            "NP_003997.1:p.Ter315TyrextTer?",
        );
    }

    /// `Ter?` already canonical.
    #[test]
    fn unknown_count_ter_form_round_trips() {
        assert_round_trips("NP_003997.1:p.Ter315TyrextTer?");
    }
}

// =============================================================================
// SECTION 2 — N-terminal extension is unaffected
// =============================================================================
//
// N-terminal extensions use a negative count (`ext-N`), not `*`/`Ter`,
// so the canonicalization change must NOT touch them. Pin to guard
// against an over-eager fix.

mod n_terminal_extension {
    use super::*;

    #[test]
    fn met_loss_short_round_trips() {
        assert_round_trips("NP_003997.1:p.Met1ext-5");
    }

    #[test]
    fn met_loss_with_new_aa_round_trips() {
        assert_round_trips("NP_003997.1:p.Met1Valext-12");
    }
}

// =============================================================================
// SECTION 3 — predicted-vs-confirmed parens preserved over extension
// =============================================================================

mod predicted_parens {
    use super::*;

    #[test]
    fn predicted_c_terminal_extension_canonicalizes() {
        assert_canonicalizes_to(
            "NP_003997.1:p.(Ter315Tyrext*10)",
            "NP_003997.1:p.(Ter315TyrextTer10)",
        );
    }

    #[test]
    fn predicted_c_terminal_extension_ter_round_trips() {
        assert_round_trips("NP_003997.1:p.(Ter315TyrextTer10)");
    }
}

// =============================================================================
// SECTION 4 — consistency with Frameshift Ter canonicalization
// =============================================================================
//
// Spec-relevant invariant: Display for protein edits emits `Ter` (never
// `*`) as the canonical stop-codon notation, regardless of input form
// and regardless of which edit kind contains it. Pin a representative
// frameshift case alongside an extension case to make the
// cross-edit-kind invariant explicit.

mod cross_edit_kind {
    use super::*;

    #[test]
    fn frameshift_star_canonicalizes_to_ter() {
        // Pre-existing behavior; pinned here so the invariant lives
        // in one place.
        assert_canonicalizes_to(
            "NP_003997.1:p.Arg97Profs*23",
            "NP_003997.1:p.Arg97ProfsTer23",
        );
    }

    #[test]
    fn substitution_star_canonicalizes_to_ter() {
        assert_canonicalizes_to("NP_003997.1:p.Tyr4*", "NP_003997.1:p.Tyr4Ter");
    }
}
