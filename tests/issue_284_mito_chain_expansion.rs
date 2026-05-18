//! Tests for issue #284: mito chain `m.3243=/A>G/=/A>T` must not
//! synthesize an `m.1=` identity for the second `=` arm.
//!
//! Follow-up to #133 work-item-3 (PR #153), which cleaned up the
//! single-slash `var/=` shorthand so the synthetic whole-entity
//! identity is rendered as bare `=` rather than expanded to
//! `<acc>:<type>.1=`. The multi-slash chain form
//! `pos=/edit/=/edit2` was tracked separately under #284 because
//! the `=` second arm needs to inherit the *position* from the
//! leading position-identity arm (so the chain remains semantically
//! coherent), not be expanded to a synthetic `m.1=`.
//!
//! Closes #284.
//!
//! Reference shapes:
//!   - `NC_012920.1:m.3243=/A>G/=/A>T` (mosaic, single-slash chain)
//!   - `NC_012920.1:m.3243=//A>G//=//A>T` (chimeric, double-slash chain)

use ferro_hgvs::hgvs::AllelePhase;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

// =============================================================================
// SECTION 1 — Core pin
// =============================================================================
//
// The bug being fixed: a bare `=` chunk after a leading position-bound
// identity arm must inherit position from the LHS, not synthesize an
// `m.1=` (or `c.1=`, etc.) whole-entity identity.

mod core_pin {
    use super::*;

    /// Core bug pin: compact + slash-chained compact must NOT render a
    /// synthetic `m.1=` identity for the second `=` arm. The second `=`
    /// should inherit position from the preceding position-identity arm.
    #[test]
    fn compact_chain_with_identity_arm_inherits_position_mt() {
        let input = "NC_012920.1:m.3243=/A>G/=/A>T";
        let parsed = parse_hgvs(input).expect("chain must parse");
        let HgvsVariant::Allele(allele) = &parsed else {
            panic!("expected Allele, got {parsed:?}");
        };
        assert_eq!(allele.phase, AllelePhase::Mosaic);
        assert_eq!(allele.variants.len(), 4);

        let display = format!("{parsed}");
        let canonical =
            "NC_012920.1:m.3243=/NC_012920.1:m.3243A>G/NC_012920.1:m.3243=/NC_012920.1:m.3243A>T";
        assert_eq!(
            display, canonical,
            "Display must not synthesize an `m.1=` identity for the second `=` arm"
        );
        assert!(
            !display.contains("m.1="),
            "Display must not contain synthetic `m.1=`; got: {display}"
        );
    }

    /// Same as the mosaic core pin but for the chimeric (double-slash)
    /// chain shape.
    #[test]
    fn compact_chain_with_identity_arm_inherits_position_mt_chimeric() {
        let input = "NC_012920.1:m.3243=//A>G//=//A>T";
        let parsed = parse_hgvs(input).expect("chimeric chain must parse");
        let HgvsVariant::Allele(allele) = &parsed else {
            panic!("expected Allele, got {parsed:?}");
        };
        assert_eq!(allele.phase, AllelePhase::Chimeric);
        assert_eq!(allele.variants.len(), 4);

        let display = format!("{parsed}");
        let canonical =
            "NC_012920.1:m.3243=//NC_012920.1:m.3243A>G//NC_012920.1:m.3243=//NC_012920.1:m.3243A>T";
        assert_eq!(display, canonical);
        assert!(!display.contains("m.1="));
    }

    /// Generality across coord systems: the fix should not be mito-
    /// specific. The same chain shape on `c.`, `g.`, `r.`, and `n.`
    /// should also inherit position rather than synthesize a `<type>.1=`
    /// or whole-entity `<type>.=` arm. (Mirrors the cross-coord-system
    /// stance taken by #133 / PR #153.)
    #[test]
    fn compact_chain_with_identity_arm_inherits_position_cds_and_genome() {
        // c. chain.
        let cds_input = "NM_000088.3:c.85=/T>C/=/T>A";
        let cds_parsed = parse_hgvs(cds_input).expect("c. chain must parse");
        let cds_display = format!("{cds_parsed}");
        assert!(
            !cds_display.contains("c.1="),
            "no synthetic c.1= in c. chain; got {cds_display}"
        );
        assert!(
            !cds_display.contains("c.=/"),
            "no synthetic whole-entity c.= in c. chain; got {cds_display}"
        );
        assert!(
            cds_display.contains("c.85="),
            "c. chain must preserve position-bound c.85=; got {cds_display}"
        );

        // g. chain — positively assert the canonical position-preserving
        // form, not just the absence of the synthetic `g.1=`.
        let g_input = "NC_000023.11:g.33344590=/A>G/=/A>T";
        let g_parsed = parse_hgvs(g_input).expect("g. chain must parse");
        let g_display = format!("{g_parsed}");
        assert!(
            !g_display.contains("g.1="),
            "no synthetic g.1= in g. chain; got {g_display}"
        );
        assert!(
            g_display.contains("g.33344590="),
            "g. chain must preserve position-bound g.33344590=; got {g_display}"
        );

        // r. chain. The bare `=` arm must inherit `r.85` from the LHS.
        let r_input = "NM_000088.3:r.85=/u>c/=/u>a";
        let r_parsed = parse_hgvs(r_input).expect("r. chain must parse");
        let r_display = format!("{r_parsed}");
        assert!(
            !r_display.contains("r.1="),
            "no synthetic r.1= in r. chain; got {r_display}"
        );
        assert!(
            !r_display.contains("r.=/"),
            "no synthetic whole-entity r.= in r. chain; got {r_display}"
        );
        assert!(
            r_display.contains("r.85="),
            "r. chain must preserve position-bound r.85=; got {r_display}"
        );

        // n. chain.
        let n_input = "NR_046018.2:n.85=/T>C/=/T>A";
        let n_parsed = parse_hgvs(n_input).expect("n. chain must parse");
        let n_display = format!("{n_parsed}");
        assert!(
            !n_display.contains("n.1="),
            "no synthetic n.1= in n. chain; got {n_display}"
        );
        assert!(
            !n_display.contains("n.=/"),
            "no synthetic whole-entity n.= in n. chain; got {n_display}"
        );
        assert!(
            n_display.contains("n.85="),
            "n. chain must preserve position-bound n.85=; got {n_display}"
        );
    }

    /// Three bare `=` chunks in a row: `m.3243=/A>G/=/=/=` must expand
    /// every `=` arm to position 3243, not just the first one. The
    /// inheritance logic must not "use up" the LHS after the first
    /// inherited arm.
    #[test]
    fn compact_chain_three_identity_arms_all_inherit_position() {
        let input = "NC_012920.1:m.3243=/A>G/=/=/=";
        let parsed = parse_hgvs(input).expect("three-= chain must parse");
        let HgvsVariant::Allele(allele) = &parsed else {
            panic!("expected Allele, got {parsed:?}");
        };
        assert_eq!(allele.variants.len(), 5);

        let display = format!("{parsed}");
        let canonical = "NC_012920.1:m.3243=/NC_012920.1:m.3243A>G/NC_012920.1:m.3243=\
                         /NC_012920.1:m.3243=/NC_012920.1:m.3243=";
        assert_eq!(
            display, canonical,
            "all three `=` arms after the leading arm must inherit position 3243"
        );
        assert!(
            !display.contains("m.1="),
            "no synthetic m.1= in three-= chain; got {display}"
        );
    }

    /// Error path: a chain led by a bare `=` (no LHS to inherit from)
    /// must be rejected. The bare `=` cannot be the first arm in a
    /// mosaic/chimeric chain because there is nothing for the `=` to
    /// reference.
    #[test]
    fn compact_chain_bare_identity_as_first_arm_errors() {
        let inputs = [
            "=/NC_012920.1:m.3243A>G",
            "=//NC_012920.1:m.3243A>G",
            "=/NM_000088.3:c.85T>C",
        ];
        for input in inputs {
            let result = parse_hgvs(input);
            assert!(
                result.is_err(),
                "bare `=` as first arm must error for `{input}`; got {result:?}"
            );
        }
    }
}

// =============================================================================
// SECTION 2 — Round-trip stability
// =============================================================================
//
// Once parsed, the canonical Display form must reparse to itself, and
// the compact input form must converge to the canonical Display under
// a second pass.

mod round_trip {
    use super::*;

    /// Round-trip stability: parse the canonical form, Display, parse,
    /// Display — second pass must be idempotent.
    #[test]
    fn compact_chain_canonical_form_is_idempotent() {
        let canonical =
            "NC_012920.1:m.3243=/NC_012920.1:m.3243A>G/NC_012920.1:m.3243=/NC_012920.1:m.3243A>T";
        let parsed = parse_hgvs(canonical).expect("canonical form must parse");
        let first = format!("{parsed}");
        let reparsed = parse_hgvs(&first).expect("reparse of canonical Display must succeed");
        let second = format!("{reparsed}");
        assert_eq!(
            first, second,
            "canonical Display must be idempotent under reparse"
        );
        assert_eq!(
            second, canonical,
            "canonical Display must equal the canonical input verbatim"
        );
    }

    /// Compact-input round-trip: the compact input `m.3243=/A>G/=/A>T`
    /// parses, and the displayed canonical form re-parses identically
    /// across the second pass.
    #[test]
    fn compact_chain_compact_input_reparses_stably() {
        let compact = "NC_012920.1:m.3243=/A>G/=/A>T";
        let parsed1 = parse_hgvs(compact).expect("compact input must parse");
        let displayed = format!("{parsed1}");
        let parsed2 = parse_hgvs(&displayed).expect("first Display must reparse");
        assert_eq!(
            format!("{parsed2}"),
            displayed,
            "second-pass Display must match first-pass Display"
        );
    }
}

// =============================================================================
// SECTION 3 — Regression guards
// =============================================================================
//
// Pre-existing shapes that must NOT change as a side effect of the
// #284 fix: single-slash `var/=` shorthand, plain two-arm chains with
// no leading identity, and mixed chains where the second arm is fully
// qualified rather than bare `=`.

mod regression_guard {
    use super::*;

    /// Boundary: single-slash `var/=` shorthand (PR #153 / #133 work-item-3)
    /// still emits bare `=` on the RHS. The fix for #284 must not regress
    /// this case — in particular, the LHS here is a non-identity variant,
    /// so the `=` arm must keep producing the whole-entity identity that
    /// Display compresses back to `var/=`.
    #[test]
    fn single_slash_var_eq_shorthand_still_works() {
        let inputs = [
            "NC_012920.1:m.3243A>G/=",
            "NC_012920.1:m.3243A>G//=",
            "NM_000088.3:c.85T>C/=",
        ];
        for input in inputs {
            let parsed = parse_hgvs(input).expect("var/= must still parse");
            assert_eq!(
                format!("{parsed}"),
                input,
                "var/= shorthand round-trip for `{input}`"
            );
        }
    }

    /// Negative control: a legitimate two-arm chain where neither arm is
    /// a leading position-identity (`m.3243A>G/m.16569A>T`) is unchanged
    /// by the fix — neither rendering nor parsing should regress.
    #[test]
    fn negative_no_leading_identity_unchanged() {
        let input = "NC_012920.1:m.3243A>G/NC_012920.1:m.16569A>T";
        let parsed = parse_hgvs(input).expect("two-arm chain must parse");
        assert_eq!(format!("{parsed}"), input);
    }

    /// Boundary: a chain whose second arm is fully-qualified (not `=`)
    /// also remains unchanged — only the bare `=` chunk should resolve
    /// through the new path.
    #[test]
    fn chain_with_fully_qualified_arms_unchanged() {
        let input = "NC_012920.1:m.3243=/A>G/NC_012920.1:m.16569A>T";
        let parsed = parse_hgvs(input).expect("mixed chain must parse");
        let display = format!("{parsed}");
        assert!(
            display.contains("NC_012920.1:m.16569A>T"),
            "second fully-qualified arm preserved; got {display}"
        );
        assert!(
            !display.contains("m.1="),
            "no synthetic m.1= in mixed chain; got {display}"
        );
    }
}
