//! Audit for issue #281 — mosaic/chimeric nesting policy and the
//! `[a/b]` bracketed mosaic form (follow-up to #217 / #216).
//!
//! Background.
//! PR #217 landed bracket-aware splitting for chimeric and non-compact
//! mosaic forms but explicitly deferred two related cases as
//! out-of-scope:
//!
//!   1. **Nesting** of `/` and `//` at the same level (e.g.
//!      `m.[3243A>G/T]//[3243A>C/G]` — chimeric whose first cell line
//!      is itself mosaic, or the mirror image mosaic-of-chimeric).
//!      HGVS v21 does not define this nesting; PR #217 rejects it
//!      with a generic "mixing mosaic and chimeric markers" message.
//!
//!   2. **`[a/b]`** — bracketed mosaic group. The HGVS spec does not
//!      define this form. The current parser routes through
//!      `parse_cis_allele` (which only knows about `;`) and emits a
//!      generic nom error.
//!
//! Policy decision pinned here.
//! For BOTH cases the parser continues to reject (the conservative,
//! spec-faithful choice — inventing a nesting policy that the spec
//! does not define would lock ferro into a semantics that might
//! conflict with a future spec revision). The rejection is upgraded
//! to a targeted diagnostic carrying the **W3019 NonSpecMosaicForm**
//! warning code and a hint pointing the user at the spec-supported
//! alternatives:
//!
//!   * compound brackets `[a;b]`,
//!   * dual fully-qualified slash `acc:c.X/acc:c.Y` (mosaic) or
//!     `acc:c.X//acc:c.Y` (chimeric),
//!   * compact mosaic short-hand `acc:c.<pos>=/<edit>`,
//!   * compact chimeric short-hand `acc:c.<pos>=//<edit>`.
//!
//! The same single code covers both sub-cases (one diagnostic family;
//! two detectors). Tests below pin both detectors plus boundary cases
//! that must continue to parse cleanly.

use ferro_hgvs::parse_hgvs;

/// Assert that `input` rejects AND the error message contains every
/// substring in `expected_substrings`. The substrings are tested
/// individually so a partial-match failure pinpoints the missing
/// fragment rather than dumping the whole message.
fn assert_rejects_with_all(input: &str, expected_substrings: &[&str]) {
    let err = match parse_hgvs(input) {
        Ok(v) => panic!("expected `{}` to reject, got {:?}", input, v),
        Err(e) => e,
    };
    let msg = err.to_string();
    for needle in expected_substrings {
        assert!(
            msg.contains(needle),
            "expected error for `{}` to contain `{}`, got `{}`",
            input,
            needle,
            msg,
        );
    }
}

// =============================================================================
// SECTION 1 — Nested mixed `/` + `//` at the same level
// =============================================================================
//
// HGVS v21 does not define chimeric-of-mosaic or mosaic-of-chimeric.
// The parser rejects with W3019 and points the user at the spec
// alternatives.

mod nested_mixed_slash {
    use super::*;

    /// Chimeric-of-mosaic: outer `//`, inner `/` inside one of the
    /// bracketed cell lines. Spec doesn't cover this nesting.
    #[test]
    fn chimeric_of_mosaic_rejected_with_w3019() {
        // Compact `m.` form from the issue body — Mt context, no
        // explicit accession, slash inside brackets on both sides.
        // Pin the EXACT rendering `[W3019 NonSpecMosaicForm]` so any
        // future drift in the structured-code prefix format (dropped
        // brackets, missing space, lower-case, etc.) trips this test
        // before reaching downstream consumers that key off the
        // bracketed form.
        let input = "m.[3243A>G/T]//[3243A>C/G]";
        let err = parse_hgvs(input).expect_err("expected rejection");
        let msg = err.to_string();
        assert!(
            msg.contains("[W3019 NonSpecMosaicForm]"),
            "expected error for `{}` to contain exact prefix `[W3019 NonSpecMosaicForm]`, got `{}`",
            input,
            msg,
        );
    }

    /// Mosaic-of-chimeric: outer `/`, inner `//` inside a bracketed
    /// group. The mirror image of the case above.
    #[test]
    fn mosaic_of_chimeric_rejected_with_w3019() {
        assert_rejects_with_all(
            "m.[3243A>G//T]/[3243A>C//G]",
            &["W3019", "NonSpecMosaicForm"],
        );
    }

    /// Non-bracketed mixed at the same level (the case PR #217 was
    /// already rejecting) MUST also surface the new W3019 code so a
    /// downstream tool can key off a single diagnostic family.
    #[test]
    fn flat_mosaic_then_chimeric_rejected_with_w3019() {
        assert_rejects_with_all(
            "NM_000088.3:c.1A>G/NM_000088.3:c.2A>G//NM_000088.3:c.3A>G",
            &["W3019", "NonSpecMosaicForm"],
        );
    }

    #[test]
    fn flat_chimeric_then_mosaic_rejected_with_w3019() {
        assert_rejects_with_all(
            "NM_000088.3:c.1A>G//NM_000088.3:c.2A>G/NM_000088.3:c.3A>G",
            &["W3019", "NonSpecMosaicForm"],
        );
    }

    /// The diagnostic message must mention the spec-supported
    /// alternatives so the user can rewrite the input.
    #[test]
    fn nested_mixed_diagnostic_mentions_spec_alternatives() {
        assert_rejects_with_all(
            "m.[3243A>G/T]//[3243A>C/G]",
            &[
                // hint must point the user somewhere actionable
                "=/",  // compact mosaic short-hand
                "=//", // compact chimeric short-hand
                "[",   // compound-bracket alternative reference
            ],
        );
    }

    /// Cross-coord coverage — the rejection is a parser-level
    /// structural rule, so it must fire regardless of which coord
    /// system the bracketed entries use. Pin `c.`, `g.`, and `r.` in
    /// addition to the `m.` flagship case above.
    #[test]
    fn cross_coord_c_chimeric_of_mosaic_rejected() {
        assert_rejects_with_all(
            "NM_000088.3:c.[100A>G/T]//NM_000088.3:c.[200A>C/G]",
            &["W3019", "NonSpecMosaicForm"],
        );
    }

    #[test]
    fn cross_coord_g_chimeric_of_mosaic_rejected() {
        assert_rejects_with_all(
            "NG_012232.1:g.[10A>G/T]//NG_012232.1:g.[20A>C/G]",
            &["W3019", "NonSpecMosaicForm"],
        );
    }

    #[test]
    fn cross_coord_r_chimeric_of_mosaic_rejected() {
        assert_rejects_with_all(
            "NM_000088.3:r.[5a>g/u]//NM_000088.3:r.[10c>g/a]",
            &["W3019", "NonSpecMosaicForm"],
        );
    }

    /// Trans-allele bracket with a flat `/` inside one arm. The outer
    /// separator is `;` (trans-allele) — not `/` or `//` — so the
    /// top-level mosaic/chimeric split never runs. The inner-`/`
    /// detector ("detector B") must still trip on the slash buried
    /// inside the first bracketed arm, confirming the detector
    /// generalises beyond the top-level mosaic/chimeric path.
    #[test]
    fn trans_allele_with_inner_slash_rejected_with_w3019() {
        assert_rejects_with_all(
            "NM_000088.3:c.[1A>G/2A>G];[3A>G]",
            &["W3019", "NonSpecMosaicForm"],
        );
    }
}

// =============================================================================
// SECTION 2 — `[a/b]` bracketed mosaic group
// =============================================================================
//
// HGVS doesn't define `[a/b]`. Today the parser routes through
// `parse_cis_allele` which knows only `;`, and the slash surfaces as a
// generic nom error. With W3019 the rejection becomes targeted.

mod bracketed_mosaic_form {
    use super::*;

    /// Bracketed `c.[100A>G/200T>C]` — the canonical shape called out
    /// in the issue body. Must reject with W3019, not a generic
    /// "Failed to parse variant" nom error.
    #[test]
    fn bracketed_mosaic_c_rejected_with_w3019() {
        assert_rejects_with_all(
            "NM_000088.3:c.[100A>G/200T>C]",
            &["W3019", "NonSpecMosaicForm"],
        );
    }

    /// Bracketed `g.` — same policy across coord systems.
    #[test]
    fn bracketed_mosaic_g_rejected_with_w3019() {
        assert_rejects_with_all(
            "NG_012232.1:g.[10A>G/20T>C]",
            &["W3019", "NonSpecMosaicForm"],
        );
    }

    /// Bracketed `n.` for the non-coding side.
    #[test]
    fn bracketed_mosaic_n_rejected_with_w3019() {
        assert_rejects_with_all(
            "NR_001234.1:n.[5A>G/10T>C]",
            &["W3019", "NonSpecMosaicForm"],
        );
    }

    /// Bracketed `r.` — mixed-case enforcement happens elsewhere; we
    /// just need the bracket+slash detector to fire here.
    #[test]
    fn bracketed_mosaic_r_rejected_with_w3019() {
        assert_rejects_with_all(
            "NM_000088.3:r.[5a>g/10u>c]",
            &["W3019", "NonSpecMosaicForm"],
        );
    }

    /// Bracketed mosaic across `m.` coord. The Mt case is the
    /// motivating example from the issue title; pin it explicitly.
    #[test]
    fn bracketed_mosaic_m_rejected_with_w3019() {
        assert_rejects_with_all("m.[3243A>G/3244T>C]", &["W3019", "NonSpecMosaicForm"]);
    }

    /// Bracketed chimeric `[a//b]` is the same family — `//` inside
    /// `[]` is also undefined. Pin the same diagnostic code.
    #[test]
    fn bracketed_chimeric_c_rejected_with_w3019() {
        assert_rejects_with_all(
            "NM_000088.3:c.[100A>G//200T>C]",
            &["W3019", "NonSpecMosaicForm"],
        );
    }

    /// The diagnostic for the bracketed form must also mention the
    /// spec alternatives, so users get a clear migration path
    /// regardless of which sub-case they tripped.
    #[test]
    fn bracketed_diagnostic_mentions_spec_alternatives() {
        assert_rejects_with_all(
            "NM_000088.3:c.[100A>G/200T>C]",
            &[
                "=/",  // compact mosaic
                "=//", // compact chimeric
                "[",   // compound-bracket alternative reference
            ],
        );
    }
}

// =============================================================================
// SECTION 3 — Boundary: legitimate forms must still parse
// =============================================================================
//
// Wiring W3019 must not regress the existing spec-supported shapes.

mod boundary_legit_forms_still_parse {
    use super::*;

    /// Fully-qualified mosaic: two complete variants joined by `/`.
    #[test]
    fn flat_mosaic_two_variants_parses() {
        let input = "NM_000088.3:c.1A>G/NM_000088.3:c.2A>G";
        let v = parse_hgvs(input).expect("legit mosaic must parse");
        assert_eq!(format!("{}", v), input);
    }

    /// Fully-qualified chimeric: two complete variants joined by `//`.
    #[test]
    fn flat_chimeric_two_variants_parses() {
        let input = "NM_000088.3:c.1A>G//NM_000088.3:c.2A>G";
        let v = parse_hgvs(input).expect("legit chimeric must parse");
        assert_eq!(format!("{}", v), input);
    }

    /// Compact mosaic short-hand `<pos>=/<edit>`.
    #[test]
    fn compact_mosaic_parses() {
        let input = "NM_000088.3:c.123=/A>G";
        let v = parse_hgvs(input).expect("compact mosaic must parse");
        assert_eq!(format!("{}", v), input);
    }

    /// Compact chimeric short-hand `<pos>=//<edit>`.
    #[test]
    fn compact_chimeric_parses() {
        let input = "NM_000088.3:c.123=//A>G";
        let v = parse_hgvs(input).expect("compact chimeric must parse");
        assert_eq!(format!("{}", v), input);
    }

    /// Bracketed inner alleles joined by `//` — spec-extension form
    /// landed in PR #217 (`[a;b]//[c;d]`). Slash is BETWEEN bracket
    /// groups (depth 0), not INSIDE one — so it must continue to
    /// parse cleanly. Canonicalizes to the compact `c.[…;…]/c.[…;…]`
    /// form, mirroring the PR #217 round-trip lattice.
    #[test]
    fn bracketed_groups_chimeric_parses() {
        let input =
            "[NM_000088.3:c.1A>G;NM_000088.3:c.2A>G]//[NM_000088.3:c.3A>G;NM_000088.3:c.4A>G]";
        let v = parse_hgvs(input).expect("bracketed chimeric must parse");
        let out = format!("{}", v);
        assert_eq!(
            out, "NM_000088.3:c.[1A>G;2A>G]//NM_000088.3:c.[3A>G;4A>G]",
            "expected canonicalization to compact form, got `{}`",
            out
        );
    }

    /// Bracketed inner alleles joined by `/` — mosaic equivalent.
    #[test]
    fn bracketed_groups_mosaic_parses() {
        let input =
            "[NM_000088.3:c.1A>G;NM_000088.3:c.2A>G]/[NM_000088.3:c.3A>G;NM_000088.3:c.4A>G]";
        let v = parse_hgvs(input).expect("bracketed mosaic must parse");
        let out = format!("{}", v);
        assert_eq!(
            out, "NM_000088.3:c.[1A>G;2A>G]/NM_000088.3:c.[3A>G;4A>G]",
            "expected canonicalization to compact form, got `{}`",
            out
        );
    }

    /// Cis allele `[a;b]` (no slash) must still parse. The W3019
    /// detector keys off `/` inside `[]` — a bare `;` cis allele
    /// must not trip it.
    #[test]
    fn cis_allele_no_slash_parses() {
        let input = "NM_000088.3:c.[100A>G;200T>C]";
        let v = parse_hgvs(input).expect("cis allele must parse");
        assert_eq!(format!("{}", v), input);
    }
}
