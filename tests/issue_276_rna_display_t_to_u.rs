//! Issue #276 — RNA `Display` must canonicalize DNA `T` to RNA `u`.
//!
//! Per HGVS v21.0, the RNA nucleotide alphabet is `a/c/g/u`. The
//! `Display` path on `RnaVariant` (and any `NaEdit::to_rna_string`
//! callsite) must emit `u` for any thymine byte, regardless of whether
//! that thymine entered the AST via:
//!
//! 1. lenient parsing of an `r.` input that contained `T`/`t` (the
//!    spec prefers `u`, but the parser is currently lenient about
//!    thymine input — see #282 for the parser-side policy decision),
//! 2. cross-coordinate translation (e.g. converting / decomposing
//!    from a `c.` representation that holds DNA bases), or
//! 3. shuffle/repeat-canonicalization that surfaces a previously
//!    unseen `T` in the emitted unit.
//!
//! Surfaced by PR #106 and PR #107's 3'-shift coverage matrices for
//! `del` and `dup`, where the RNA module's finer-periodicity sub-case
//! had to work around the bug by using an `AC` tandem rather than the
//! genomic/cds/noncoding modules' `AT` tandem. See
//! `tests/del_shift_matrix.rs` for the comment trail.
//!
//! Probe shapes covered here:
//! - `r.<pos><unit-with-T>[N]` (repeat-unit emission, the original
//!   regression),
//! - `r.<a>_<b>insT...` / `r.<pos>delT` / `r.<pos>dupT` (stated-base
//!   forms across edit shapes — all share the same `to_lowercase_char`
//!   / `to_lowercase_string` helper),
//! - `r.<pos>A>T` (substitution alt),
//! - `r.<a>_<b>delT...insT...` (delins both legs),
//! - `r.<a>_<b>invT...` (inversion stated ref),
//! - compound cis bracketed member with `T`.

use ferro_hgvs::parse_hgvs;

const ACC: &str = "NM_000088.3";

/// Parse `input`, format via `Display`, assert it equals `expected`.
#[track_caller]
fn assert_canonicalizes_to(input: &str, expected: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, expected, "input={input:?}");
}

// ----------------------------------------------------------------------
// Repeat units — the original regression from PR #106 / PR #107.
// ----------------------------------------------------------------------

mod repeat_unit {
    use super::*;

    /// Lowercase `at[5]` input — the unit is already lowercase but
    /// contains a DNA `T`. `Display` must emit `u`, not `t`.
    #[test]
    fn lowercase_unit_with_t_canonicalizes_to_u() {
        assert_canonicalizes_to(&format!("{ACC}:r.123at[5]"), &format!("{ACC}:r.123au[5]"));
    }

    /// Uppercase `AT[5]` input — mixed-case lenient input. `Display`
    /// lowercases A→a and translates T→u.
    #[test]
    fn uppercase_unit_with_t_canonicalizes_to_au() {
        assert_canonicalizes_to(&format!("{ACC}:r.123AT[5]"), &format!("{ACC}:r.123au[5]"));
    }

    /// Mixed-case `At[5]` input — A→a, T→u.
    #[test]
    fn mixed_case_unit_with_t_canonicalizes_to_au() {
        assert_canonicalizes_to(&format!("{ACC}:r.123At[5]"), &format!("{ACC}:r.123au[5]"));
    }

    /// Multi-byte unit `gat[3]` — all-lowercase but contains `t`.
    /// Pins that the translation is per-byte, not just leading.
    #[test]
    fn multi_byte_unit_with_t_canonicalizes_to_u() {
        assert_canonicalizes_to(&format!("{ACC}:r.123gat[3]"), &format!("{ACC}:r.123gau[3]"));
    }

    /// Multi-`T` unit `att[4]` — pins that every `t` byte in the
    /// emitted unit becomes `u`, not just the first.
    #[test]
    fn multi_t_unit_canonicalizes_each_t_to_u() {
        assert_canonicalizes_to(&format!("{ACC}:r.123att[4]"), &format!("{ACC}:r.123auu[4]"));
    }
}

// ----------------------------------------------------------------------
// Substitution alt with `T`.
// ----------------------------------------------------------------------

mod substitution {
    use super::*;

    #[test]
    fn lowercase_t_alt_canonicalizes_to_u() {
        assert_canonicalizes_to(&format!("{ACC}:r.123a>t"), &format!("{ACC}:r.123a>u"));
    }

    #[test]
    fn uppercase_t_alt_canonicalizes_to_u() {
        assert_canonicalizes_to(&format!("{ACC}:r.123A>T"), &format!("{ACC}:r.123a>u"));
    }

    /// `t>g` — ref with `T`. The ref position is the leg that
    /// previously rendered as lowercase `t` on r. Display.
    #[test]
    fn lowercase_t_ref_canonicalizes_to_u() {
        assert_canonicalizes_to(&format!("{ACC}:r.123t>g"), &format!("{ACC}:r.123u>g"));
    }
}

// ----------------------------------------------------------------------
// Deletion / Duplication / Inversion stated-ref forms.
// ----------------------------------------------------------------------

mod deletion {
    use super::*;

    #[test]
    fn stated_ref_with_t_canonicalizes_to_u() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125delaut"),
            &format!("{ACC}:r.123_125delauu"),
        );
    }

    #[test]
    fn uppercase_stated_ref_with_t_canonicalizes() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125delAUT"),
            &format!("{ACC}:r.123_125delauu"),
        );
    }

    #[test]
    fn single_base_del_t_canonicalizes_to_u() {
        assert_canonicalizes_to(&format!("{ACC}:r.123delT"), &format!("{ACC}:r.123delu"));
    }
}

mod duplication {
    use super::*;

    #[test]
    fn stated_ref_with_t_canonicalizes_to_u() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125dupaut"),
            &format!("{ACC}:r.123_125dupauu"),
        );
    }

    #[test]
    fn uppercase_stated_ref_with_t_canonicalizes() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125dupAUT"),
            &format!("{ACC}:r.123_125dupauu"),
        );
    }

    /// Single-base `dupT` — mirrors the single-base `delT` case in the
    /// `deletion` module. Pins that the dup stated-base path also
    /// translates a lone `T` byte to `u`.
    #[test]
    fn single_base_dup_t_canonicalizes_to_u() {
        assert_canonicalizes_to(&format!("{ACC}:r.123dupT"), &format!("{ACC}:r.123dupu"));
    }
}

mod inversion {
    use super::*;

    #[test]
    fn stated_ref_with_t_canonicalizes_to_u() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125invaut"),
            &format!("{ACC}:r.123_125invauu"),
        );
    }

    #[test]
    fn uppercase_stated_ref_with_t_canonicalizes() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125invAUT"),
            &format!("{ACC}:r.123_125invauu"),
        );
    }
}

// ----------------------------------------------------------------------
// Insertion / Delins — inserted alt with `T`.
// ----------------------------------------------------------------------

mod insertion {
    use super::*;

    #[test]
    fn inserted_alt_with_t_canonicalizes_to_u() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_124insAT"),
            &format!("{ACC}:r.123_124insau"),
        );
    }

    #[test]
    fn lowercase_inserted_alt_with_t_canonicalizes() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_124insat"),
            &format!("{ACC}:r.123_124insau"),
        );
    }
}

mod delins {
    use super::*;

    /// Both legs contain a `T`. Both must translate to `u`.
    #[test]
    fn t_in_both_ref_and_alt_canonicalize_to_u() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125delAUTinsT"),
            &format!("{ACC}:r.123_125delauuinsu"),
        );
    }

    #[test]
    fn t_only_in_alt_canonicalizes_to_u() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125delinsTTT"),
            &format!("{ACC}:r.123_125delinsuuu"),
        );
    }
}

// ----------------------------------------------------------------------
// Identity stated-ref — `r.<pos>...=` form. Pins canonicalization on
// the `=` (no-change) edit branch.
// ----------------------------------------------------------------------

mod identity {
    use super::*;

    /// Stated-ref `auT=` — Identity with a `T` byte in the ref-leg.
    /// `Display` must emit `auu=`, mirroring the `del`/`dup`/`inv`
    /// stated-ref branches.
    #[test]
    fn stated_ref_with_t_canonicalizes_to_u() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125auT="),
            &format!("{ACC}:r.123_125auu="),
        );
    }

    /// Uppercase `AUT=` — mixed-case lenient input through the
    /// Identity branch.
    #[test]
    fn uppercase_stated_ref_with_t_canonicalizes() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125AUT="),
            &format!("{ACC}:r.123_125auu="),
        );
    }
}

// ----------------------------------------------------------------------
// Multi-unit repeat — `<unit1>[N]<unit2>[M]...` form. Pins that every
// unit's bytes are canonicalized, not just the first.
// ----------------------------------------------------------------------

mod multi_repeat {
    use super::*;

    /// First unit `aT[3]` contains a `T`; second unit `ac[2]` is pure
    /// RNA. `Display` must translate the `T` in the first unit to `u`
    /// while leaving the second unit untouched.
    #[test]
    fn t_in_first_unit_canonicalizes_to_u() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123aT[3]ac[2]"),
            &format!("{ACC}:r.123au[3]ac[2]"),
        );
    }

    /// `T` in the second unit — pins that the per-unit translation
    /// runs across every unit, not just the first.
    #[test]
    fn t_in_second_unit_canonicalizes_to_u() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123ac[2]aT[3]"),
            &format!("{ACC}:r.123ac[2]au[3]"),
        );
    }
}

// ----------------------------------------------------------------------
// Compound cis allele bracket — pins canonicalization inside `[...;...]`.
// ----------------------------------------------------------------------

mod compound_cis {
    use super::*;

    #[test]
    fn member_with_t_canonicalizes_inside_brackets() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.[123A>T;200_202delinsTTT]"),
            &format!("{ACC}:r.[123a>u;200_202delinsuuu]"),
        );
    }
}
