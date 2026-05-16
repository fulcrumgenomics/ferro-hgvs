//! Audit for issue #216 — chimeric (`//`) parse + round-trip and
//! non-compact mosaic (`/`) forms.
//!
//! Spec references:
//!   - `assets/hgvs-nomenclature/docs/general.md` (slash semantics).
//!   - `assets/hgvs-nomenclature/docs/recommendations/DNA/substitution.md`
//!     (compact mosaic + chimeric examples).
//!   - `assets/hgvs-nomenclature/docs/recommendations/DNA/alleles.md`
//!     (brackets + phase markers).
//!
//! Prior work landed:
//!   - #146 — trans-phase `[a];[b]` canonicalization.
//!   - #148 — unknown phase `[a(;)b]` (closes #123).
//!   - #153 — spec compact mosaic / chimeric `<pos>=/<edit>` and
//!     `<pos>=//<edit>` (closes #133).
//!
//! Remaining C2 coverage pinned here:
//!   - 2+ variant chimeric (`var//var`, `var//var//var`,
//!     `var//var//var//var`).
//!   - Cross-accession chimeric (`c.X//g.Y`).
//!   - Spec compact + `var//=` short-hand round-trip.
//!   - **Bracketed inner alleles**: `[a;b]//[c;d]` and `[a;b]/[c;d]` —
//!     chimerism / mosaicism between cis-allele groups. Defensible by
//!     compositional extension over the spec's cis production; the
//!     spec's formal grammar does not include `/` or `//` operators
//!     and the narrative only shows the compact `<pos>=/<edit>` form,
//!     so this is a parser extension. Previously rejected because the
//!     slash split was not bracket-aware.
//!   - **Standalone singleton bracket** `c.[a]` is REJECTED — the
//!     spec's cis production (`syntax.yaml:134-135`) requires `;` and
//!     ≥2 entries, and the committee explicitly classifies `c.[76A>C]`
//!     as invalid (`alleles.md:99-101`). Singleton-bracket phase
//!     operands `[a]//[b]` / `[a]/[b]` are also rejected — they add no
//!     expressiveness over `a//b` / `a/b` and have no spec referent.
//!   - Mosaic equivalents for the same lattice.
//!   - Mt + Circular coord systems chimeric / mosaic.
//!   - Triple-slash `c.X///c.Y` is malformed and must surface a parse
//!     error rather than silently round-trip.

use ferro_hgvs::hgvs::variant::{AllelePhase, HgvsVariant};
use ferro_hgvs::parse_hgvs;

fn assert_parse_roundtrips(input: &str) -> HgvsVariant {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{}`: {}", input, e));
    let out = format!("{}", v);
    assert_eq!(
        out, input,
        "round-trip mismatch:\n  input:  `{}`\n  output: `{}`",
        input, out
    );
    v
}

fn assert_parse_canonicalizes_to(input: &str, expected_output: &str) -> HgvsVariant {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{}`: {}", input, e));
    let out = format!("{}", v);
    assert_eq!(
        out, expected_output,
        "expected canonicalization mismatch:\n  input:    `{}`\n  expected: `{}`\n  got:      `{}`",
        input, expected_output, out
    );
    v
}

fn assert_rejects_with(input: &str, expected_substring: &str) {
    let err = match parse_hgvs(input) {
        Ok(v) => panic!("expected `{}` to reject, got {:?}", input, v),
        Err(e) => e,
    };
    let msg = err.to_string();
    assert!(
        msg.contains(expected_substring),
        "expected error for `{}` to contain `{}`, got `{}`",
        input,
        expected_substring,
        msg,
    );
}

fn expect_phase(variant: &HgvsVariant, phase: AllelePhase, count: usize) {
    let allele = match variant {
        HgvsVariant::Allele(a) => a,
        _ => panic!("expected AlleleVariant, got {:?}", variant),
    };
    assert_eq!(allele.phase, phase, "wrong AllelePhase");
    assert_eq!(
        allele.variants.len(),
        count,
        "wrong sub-variant count, got {}",
        allele.variants.len()
    );
}

// =============================================================================
// SECTION 1 — chimeric (`//`) full-form round-trip lattice
// =============================================================================

mod chimeric_long_form {
    use super::*;

    /// Two fully-qualified variants joined by `//`. Already exercised by
    /// existing parser tests; pinned here as the floor of the lattice.
    #[test]
    fn two_variants_same_accession_round_trips() {
        let v = assert_parse_roundtrips("NM_000088.3:c.100A>G//NM_000088.3:c.200C>T");
        expect_phase(&v, AllelePhase::Chimeric, 2);
    }

    /// Three-way chimeric (three distinct cell populations).
    #[test]
    fn three_variants_round_trips() {
        let v =
            assert_parse_roundtrips("NM_000088.3:c.1A>G//NM_000088.3:c.2A>G//NM_000088.3:c.3A>G");
        expect_phase(&v, AllelePhase::Chimeric, 3);
    }

    /// Four-way chimeric (real-world chimeras can carry more than three
    /// populations).
    #[test]
    fn four_variants_round_trips() {
        let v = assert_parse_roundtrips(
            "NM_000088.3:c.1A>G//NM_000088.3:c.2A>G//NM_000088.3:c.3A>G//NM_000088.3:c.4A>G",
        );
        expect_phase(&v, AllelePhase::Chimeric, 4);
    }

    /// Cross-accession chimeric (e.g. one variant on transcript, one on
    /// genomic context). Spec doesn't forbid this — it's analogous to
    /// the unknown-phase cross-reference case already pinned for
    /// `(;)`.
    #[test]
    fn cross_accession_round_trips() {
        let v = assert_parse_roundtrips("NM_000088.3:c.1A>G//NC_000017.11:g.2A>G");
        expect_phase(&v, AllelePhase::Chimeric, 2);
    }

    /// `var//=` short-hand: second cell population carries the
    /// reference. Existing parser test covers the basic shape — pinned
    /// here in audit context.
    #[test]
    fn eq_shorthand_rhs_round_trips() {
        let v = assert_parse_roundtrips("NM_000088.3:c.456C>T//=");
        expect_phase(&v, AllelePhase::Chimeric, 2);
    }

    /// Spec compact form: `<acc>:<type>.<pos>=//<edit>`.
    #[test]
    fn spec_compact_form_round_trips() {
        let v = assert_parse_roundtrips("NC_012920.1:m.3243=//A>G");
        expect_phase(&v, AllelePhase::Chimeric, 2);
    }

    /// Mt-context chimeric (heteroplasmy adjacent — mitochondrial cell
    /// populations carrying different alleles).
    #[test]
    fn mt_round_trips() {
        let v = assert_parse_roundtrips("NC_012920.1:m.3243A>G//NC_012920.1:m.3243A>T");
        expect_phase(&v, AllelePhase::Chimeric, 2);
    }

    /// Circular (o.) coord system chimeric.
    #[test]
    fn circular_round_trips() {
        let v = assert_parse_roundtrips("NC_000017.11:o.100A>G//NC_000017.11:o.200C>T");
        expect_phase(&v, AllelePhase::Chimeric, 2);
    }

    /// Inherited-accession RHS canonicalizes to fully-qualified form.
    /// (Already covered by existing tests; pinned here so a regression
    /// doesn't accidentally drop the canonicalization.)
    #[test]
    fn inherited_accession_canonicalizes() {
        let v = assert_parse_canonicalizes_to(
            "NM_000088.3:c.1A>G//c.2A>G",
            "NM_000088.3:c.1A>G//NM_000088.3:c.2A>G",
        );
        expect_phase(&v, AllelePhase::Chimeric, 2);
    }

    /// Inherited-accession RHS for the circular `o.` coord system —
    /// `parse_variant_with_inherited_accession` must accept bare `o.X`
    /// when the LHS already established the accession.
    #[test]
    fn inherited_accession_circular_canonicalizes() {
        let v = assert_parse_canonicalizes_to(
            "NC_000017.11:o.100A>G//o.200C>T",
            "NC_000017.11:o.100A>G//NC_000017.11:o.200C>T",
        );
        expect_phase(&v, AllelePhase::Chimeric, 2);
    }
}

// =============================================================================
// SECTION 2 — chimeric with BRACKETED inner alleles (`[a;b]//[c;d]`)
// =============================================================================
//
// The principled spec reading: each "allele" in a chimeric expression
// can itself be a multi-variant cis group `[a;b]`. Previously rejected
// because the slash split was not bracket-aware.

mod chimeric_bracketed_inner {
    use super::*;

    /// Two cis-allele groups joined by `//`: `[a;b]//[c;d]`.
    #[test]
    fn two_cis_groups_round_trips() {
        let input = "NM_000088.3:c.[1A>G;2A>G]//NM_000088.3:c.[3A>G;4A>G]";
        let v = assert_parse_roundtrips(input);
        expect_phase(&v, AllelePhase::Chimeric, 2);
        // Each sub-variant should itself be a Cis allele.
        let HgvsVariant::Allele(a) = &v else {
            unreachable!()
        };
        for sub in &a.variants {
            match sub {
                HgvsVariant::Allele(inner) => {
                    assert_eq!(inner.phase, AllelePhase::Cis);
                    assert_eq!(inner.variants.len(), 2);
                }
                _ => panic!("expected nested Cis AlleleVariant, got {:?}", sub),
            }
        }
    }

    /// Long-form expansion of the bracketed inner is accepted with the
    /// outer accession spelled per inner variant.
    #[test]
    fn long_form_expansion_round_trips() {
        let input =
            "[NM_000088.3:c.1A>G;NM_000088.3:c.2A>G]//[NM_000088.3:c.3A>G;NM_000088.3:c.4A>G]";
        // The canonical output uses the compact `c.[a;b]` form for each
        // inner cis group.
        assert_parse_canonicalizes_to(
            input,
            "NM_000088.3:c.[1A>G;2A>G]//NM_000088.3:c.[3A>G;4A>G]",
        );
    }
}

// =============================================================================
// SECTION 3 — mosaic (`/`) non-compact long-form lattice
// =============================================================================

mod mosaic_long_form {
    use super::*;

    #[test]
    fn two_variants_round_trips() {
        let v = assert_parse_roundtrips("NM_000088.3:c.100A>G/NM_000088.3:c.200C>T");
        expect_phase(&v, AllelePhase::Mosaic, 2);
    }

    #[test]
    fn three_variants_round_trips() {
        let v = assert_parse_roundtrips("NM_000088.3:c.1A>G/NM_000088.3:c.2A>G/NM_000088.3:c.3A>G");
        expect_phase(&v, AllelePhase::Mosaic, 3);
    }

    #[test]
    fn cross_accession_round_trips() {
        let v = assert_parse_roundtrips("NM_000088.3:c.1A>G/NC_000017.11:g.2A>G");
        expect_phase(&v, AllelePhase::Mosaic, 2);
    }

    #[test]
    fn eq_shorthand_rhs_round_trips() {
        let v = assert_parse_roundtrips("NM_000088.3:c.123A>G/=");
        expect_phase(&v, AllelePhase::Mosaic, 2);
    }

    #[test]
    fn spec_compact_form_round_trips() {
        let v = assert_parse_roundtrips("NM_000088.3:c.85=/T>C");
        expect_phase(&v, AllelePhase::Mosaic, 2);
    }

    #[test]
    fn mt_round_trips() {
        let v = assert_parse_roundtrips("NC_012920.1:m.3243A>G/NC_012920.1:m.3243A>T");
        expect_phase(&v, AllelePhase::Mosaic, 2);
    }

    #[test]
    fn circular_round_trips() {
        let v = assert_parse_roundtrips("NC_000017.11:o.100A>G/NC_000017.11:o.200C>T");
        expect_phase(&v, AllelePhase::Mosaic, 2);
    }

    #[test]
    fn inherited_accession_canonicalizes() {
        let v = assert_parse_canonicalizes_to(
            "NM_000088.3:c.1A>G/c.2A>G",
            "NM_000088.3:c.1A>G/NM_000088.3:c.2A>G",
        );
        expect_phase(&v, AllelePhase::Mosaic, 2);
    }

    /// Inherited-accession RHS for the circular `o.` coord system.
    #[test]
    fn inherited_accession_circular_canonicalizes() {
        let v = assert_parse_canonicalizes_to(
            "NC_000017.11:o.100A>G/o.200C>T",
            "NC_000017.11:o.100A>G/NC_000017.11:o.200C>T",
        );
        expect_phase(&v, AllelePhase::Mosaic, 2);
    }
}

// =============================================================================
// SECTION 4 — mosaic with BRACKETED inner alleles (`[a;b]/[c;d]`)
// =============================================================================

mod mosaic_bracketed_inner {
    use super::*;

    #[test]
    fn two_cis_groups_round_trips() {
        let input = "NM_000088.3:c.[1A>G;2A>G]/NM_000088.3:c.[3A>G;4A>G]";
        let v = assert_parse_roundtrips(input);
        expect_phase(&v, AllelePhase::Mosaic, 2);
        // Each sub-variant should itself be a Cis allele (mirrors the
        // chimeric version's structural assertion).
        let HgvsVariant::Allele(a) = &v else {
            unreachable!()
        };
        for sub in &a.variants {
            match sub {
                HgvsVariant::Allele(inner) => {
                    assert_eq!(inner.phase, AllelePhase::Cis);
                    assert_eq!(inner.variants.len(), 2);
                }
                _ => panic!("expected nested Cis AlleleVariant, got {:?}", sub),
            }
        }
    }

    #[test]
    fn long_form_expansion_round_trips() {
        let input =
            "[NM_000088.3:c.1A>G;NM_000088.3:c.2A>G]/[NM_000088.3:c.3A>G;NM_000088.3:c.4A>G]";
        assert_parse_canonicalizes_to(input, "NM_000088.3:c.[1A>G;2A>G]/NM_000088.3:c.[3A>G;4A>G]");
    }
}

// =============================================================================
// SECTION 5 — Malformed inputs (clean rejection)
// =============================================================================

mod malformed {
    use super::*;

    /// Triple slash `///` is not in the spec. The parser must reject
    /// rather than silently round-trip — three slashes is unambiguously
    /// a typo for either `/` or `//` and we should not pick for the
    /// user.
    #[test]
    fn triple_slash_rejected() {
        assert_rejects_with(
            "NM_000088.3:c.1A>G///NM_000088.3:c.2A>G",
            "triple-slash or stray separator",
        );
    }

    /// Empty chunks (`var//`) reject (existing behavior, pinned).
    #[test]
    fn empty_rhs_rejected() {
        assert_rejects_with("NM_000088.3:c.1A>G//", "Empty variant in chimeric allele");
        assert_rejects_with("NM_000088.3:c.1A>G/", "Empty variant in mosaic allele");
    }

    /// `=` cannot be the first variant in mosaic / chimeric (the
    /// reference must always come first per spec). Pinned.
    #[test]
    fn eq_first_rejected() {
        assert_rejects_with(
            "=//NM_000088.3:c.123A>G",
            "Cannot use '=' as first variant in chimeric notation",
        );
        assert_rejects_with(
            "=/NM_000088.3:c.123A>G",
            "Cannot use '=' as first variant in mosaic notation",
        );
    }

    /// Mixing mosaic and chimeric phase markers at the same nesting
    /// level is ambiguous — the spec doesn't define chimeric-of-mosaic
    /// or mosaic-of-chimeric semantics, so the parser rejects rather
    /// than silently picking an interpretation (the old code parsed it
    /// as `Chimeric([Mosaic([a, b]), c])`, which is not in the spec).
    #[test]
    fn mixed_slash_chimeric_with_inner_mosaic_rejected() {
        assert_rejects_with(
            "NM_000088.3:c.1A>G/NM_000088.3:c.2A>G//NM_000088.3:c.3A>G",
            "mixing mosaic and chimeric markers at the same nesting level",
        );
    }

    #[test]
    fn mixed_slash_chimeric_trailing_mosaic_rejected() {
        assert_rejects_with(
            "NM_000088.3:c.1A>G//NM_000088.3:c.2A>G/NM_000088.3:c.3A>G",
            "mixing mosaic and chimeric markers at the same nesting level",
        );
    }
}

// =============================================================================
// SECTION 6 — Stand-alone singleton bracket `[a]` is rejected
// =============================================================================
//
// The HGVS DNA cis grammar production is
// `"[" position_edit ";" position_edit "]"` (syntax.yaml lines 134-135)
// — it requires the `;` separator and ≥2 position_edits. The committee
// has also explicitly addressed `c.[76A>C]` standalone (alleles.md
// lines 99-101) and classified it as invalid: "the recommended
// description is `LRG_199t1:c.[76A>C];[76=]`". The cis detector must
// therefore require `;` inside the brackets — a singleton bracket
// expression carries no information beyond the bare inner variant
// and the spec does not generate it.

mod stand_alone_singleton_bracket_rejected {
    use super::*;

    #[test]
    fn singleton_bracket_rejected() {
        let result = parse_hgvs("[NM_000088.3:c.1A>G]");
        assert!(
            result.is_err(),
            "standalone singleton bracket `[a]` must reject per spec \
             (DNA alleles cis grammar requires `;`; alleles.md:99-101 \
             explicitly rejects `c.[76A>C]`), got {:?}",
            result
        );
    }

    /// Singleton-bracket chimeric / mosaic operands `[a]//[b]` and
    /// `[a]/[b]` are also rejected. They add no expressiveness over
    /// `a//b` / `a/b` (a single bracketed entry carries no cis
    /// grouping information), and the HGVS spec does not generate
    /// them — `/` and `//` have no formal productions in syntax.yaml,
    /// and the narrative only shows the compact `<pos>=/<edit>` form.
    #[test]
    fn singleton_bracket_chimeric_operand_rejected() {
        let result = parse_hgvs("[NM_000088.3:c.1A>G]//[NM_000088.3:c.2A>G]");
        assert!(
            result.is_err(),
            "singleton-bracket chimeric operand `[a]//[b]` must reject, got {:?}",
            result
        );
    }

    #[test]
    fn singleton_bracket_mosaic_operand_rejected() {
        let result = parse_hgvs("[NM_000088.3:c.1A>G]/[NM_000088.3:c.2A>G]");
        assert!(
            result.is_err(),
            "singleton-bracket mosaic operand `[a]/[b]` must reject, got {:?}",
            result
        );
    }
}
