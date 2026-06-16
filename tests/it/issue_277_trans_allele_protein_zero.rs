//! Issue #277 — compound trans-allele Display→reparse drift for protein `p.0`.
//!
//! When a trans-allele compound contains `ProteinEdit::NoProtein` on one arm
//! (i.e. `[ACC:p.X];[ACC:p.0]`), the compact-form Display emits
//! `ACC:p.[X];[0]`. The bare-bracket `[0]` previously re-parsed as
//! `HgvsVariant::NullAllele`, NOT as the original `ProteinEdit::NoProtein` —
//! semantic drift between two structurally distinct concepts:
//!
//! - `NullAllele` (`[0]` in a non-coordinate-system context, e.g.
//!   `[ACC:p.X];[0]` with no shared `p.` prefix): the second X-chromosome is
//!   *absent* (alleles.md line 78).
//! - `ProteinEdit::NoProtein` (`p.0` whole-protein no-product): the second
//!   allele *exists* but produces no protein due to e.g. a start-codon
//!   variant (substitution.md line 42).
//!
//! Inside the protein-coordinate compact trans-allele form (`ACC:p.[X];[0]`),
//! the coordinate context is already pinned to `p.` by the shared compact
//! prefix; therefore `[0]` must resolve to `ProteinEdit::NoProtein`, not
//! `NullAllele`. This file pins the post-fix behavior.
//!
//! Cross-ref: PR #130, issue #81 § D6.

use ferro_hgvs::hgvs::edit::ProteinEdit;
use ferro_hgvs::hgvs::uncertainty::Mu;
use ferro_hgvs::hgvs::variant::{AllelePhase, ProteinVariant};
use ferro_hgvs::{parse_hgvs, HgvsVariant};

/// Extract the inner `ProteinEdit` from a `HgvsVariant::Protein`, asserting the
/// outer `Mu` wrapper is `Mu::Certain` (the canonical shape for `p.0` /
/// `p.0?`; see `protein_no_protein_roundtrip::p0_outer_mu_is_certain`).
fn protein_edit(variant: &HgvsVariant) -> &ProteinEdit {
    let pv: &ProteinVariant = match variant {
        HgvsVariant::Protein(pv) => pv,
        other => panic!("expected HgvsVariant::Protein, got {:?}", other),
    };
    match &pv.loc_edit.edit {
        Mu::Certain(edit) => edit,
        other => panic!(
            "expected outer Mu::Certain on protein no-product arm, got {:?}",
            other
        ),
    }
}

/// Compact-form behavior: bracketed `[0]` / `[0?]` inside a protein-coordinate
/// compact trans-allele routes to `ProteinEdit::NoProtein`, and the full
/// Display→reparse round-trip preserves the original structure.
mod compact_form {
    use super::*;

    /// Compact-form `[0]` inside a protein trans-allele must route to
    /// `ProteinEdit::NoProtein`, not `HgvsVariant::NullAllele`. This is the
    /// core bug.
    #[test]
    fn compact_trans_allele_zero_routes_to_no_protein() {
        let parsed = parse_hgvs("NP_003997.2:p.[Arg97Trp];[0]")
            .expect("parse compact trans-allele with bare-bracket zero");

        let allele = match &parsed {
            HgvsVariant::Allele(a) => a,
            other => panic!("expected Allele variant, got {:?}", other),
        };
        assert_eq!(allele.phase, AllelePhase::Trans, "must be trans phase");
        assert_eq!(allele.variants.len(), 2, "must have exactly two arms");

        // First arm: ACC:p.Arg97Trp — concrete protein substitution.
        assert!(
            matches!(&allele.variants[0], HgvsVariant::Protein(_)),
            "first arm must be Protein, got {:?}",
            allele.variants[0]
        );

        // Second arm: must be Protein(NoProtein), NOT NullAllele.
        assert!(
            matches!(&allele.variants[1], HgvsVariant::Protein(_)),
            "second arm must be Protein (carrying NoProtein); got {:?}",
            allele.variants[1]
        );
        match protein_edit(&allele.variants[1]) {
            ProteinEdit::NoProtein { predicted: false } => {}
            other => panic!(
                "second arm must carry ProteinEdit::NoProtein {{ predicted: false }}, got {:?}",
                other
            ),
        }
    }

    /// `[0]` as the *first* arm of a compact protein trans-allele is also
    /// `NoProtein`: coordinate context is shared, position-of-arm is
    /// irrelevant.
    #[test]
    fn compact_trans_allele_zero_first_arm_routes_to_no_protein() {
        let parsed = parse_hgvs("NP_003997.2:p.[0];[Arg97Trp]")
            .expect("parse compact trans-allele with bare-bracket zero first");

        let allele = match &parsed {
            HgvsVariant::Allele(a) => a,
            other => panic!("expected Allele variant, got {:?}", other),
        };
        assert_eq!(allele.phase, AllelePhase::Trans);
        assert_eq!(allele.variants.len(), 2);

        // First arm: Protein(NoProtein).
        assert!(
            matches!(&allele.variants[0], HgvsVariant::Protein(_)),
            "first arm must be Protein (carrying NoProtein); got {:?}",
            allele.variants[0]
        );
        match protein_edit(&allele.variants[0]) {
            ProteinEdit::NoProtein { predicted: false } => {}
            other => panic!(
                "first arm must carry ProteinEdit::NoProtein {{ predicted: false }}, got {:?}",
                other
            ),
        }

        // Second arm: ACC:p.Arg97Trp.
        assert!(matches!(&allele.variants[1], HgvsVariant::Protein(_)));
    }

    /// Three-member compact trans-allele with `[0]` in the *middle* position
    /// must route the middle arm to `ProteinEdit::NoProtein` (the shared
    /// `p.` prefix pins coordinate context for every bracketed member,
    /// regardless of arm count or arm index). Round-trip: Display→reparse
    /// preserves shape.
    #[test]
    fn compact_trans_allele_three_member_zero_middle_routes_to_no_protein() {
        let input = "NP_003997.2:p.[Arg97Trp];[0];[Met1Val]";
        let parsed = parse_hgvs(input).expect("parse 3-member compact trans-allele with [0] mid");

        let allele = match &parsed {
            HgvsVariant::Allele(a) => a,
            other => panic!("expected Allele, got {:?}", other),
        };
        assert_eq!(allele.phase, AllelePhase::Trans);
        assert_eq!(allele.variants.len(), 3, "must have exactly three arms");

        // All three arms must be Protein(_).
        for (i, arm) in allele.variants.iter().enumerate() {
            assert!(
                matches!(arm, HgvsVariant::Protein(_)),
                "arm {} must be Protein, got {:?}",
                i,
                arm
            );
        }

        // Middle arm must carry NoProtein.
        match protein_edit(&allele.variants[1]) {
            ProteinEdit::NoProtein { predicted: false } => {}
            other => panic!(
                "middle arm must carry ProteinEdit::NoProtein {{ predicted: false }}, got {:?}",
                other
            ),
        }

        // Round-trip Display→reparse must yield the structurally identical
        // variant.
        let displayed = format!("{}", parsed);
        let reparsed = parse_hgvs(&displayed).expect("reparse 3-member compact form");
        assert_eq!(
            parsed, reparsed,
            "3-member compact trans-allele round-trip drift; displayed = {:?}",
            displayed
        );
    }

    /// Round-trip: parsing the expanded form `[ACC:p.X];[ACC:p.0]` followed
    /// by Display→reparse must produce a structurally identical variant. The
    /// Display path emits the compact form `ACC:p.[X];[0]`; before the fix,
    /// reparse produced a different shape (`NullAllele` instead of
    /// `NoProtein`), silently breaking round-trip semantics.
    #[test]
    fn expanded_form_display_reparse_roundtrip_preserves_no_protein() {
        let input = "[NP_003997.2:p.Arg97Trp];[NP_003997.2:p.0]";
        let first = parse_hgvs(input).expect("parse expanded form");

        // Display should emit the compact form ACC:p.[X];[0].
        let displayed = format!("{}", first);
        assert_eq!(
            displayed, "NP_003997.2:p.[Arg97Trp];[0]",
            "compact-form Display drift"
        );

        // Reparsing the Display output must yield the same structure.
        let second = parse_hgvs(&displayed).expect("reparse compact form");
        assert_eq!(
            first, second,
            "Display->reparse round-trip must be structurally identical; this is the #277 drift"
        );

        // And the second arm must still be Protein(NoProtein), not
        // NullAllele.
        let allele = match &second {
            HgvsVariant::Allele(a) => a,
            other => panic!("expected Allele on reparse, got {:?}", other),
        };
        match protein_edit(&allele.variants[1]) {
            ProteinEdit::NoProtein { predicted: false } => {}
            other => panic!("expected NoProtein on reparsed second arm, got {:?}", other),
        }
    }

    /// Fully-qualified expanded form with `p.0` as the *first* arm —
    /// `[ACC:p.0];[ACC:p.X]` — must Display→reparse cleanly with the first
    /// arm carrying `ProteinEdit::NoProtein`. Arm order must not affect the
    /// round-trip.
    #[test]
    fn expanded_form_first_arm_zero_roundtrips_to_no_protein() {
        let input = "[NP_003997.2:p.0];[NP_003997.2:p.Arg97Trp]";
        let first = parse_hgvs(input).expect("parse expanded form with p.0 first");

        let displayed = format!("{}", first);
        let second = parse_hgvs(&displayed).expect("reparse Display output");
        assert_eq!(
            first, second,
            "first-arm p.0 expanded->compact round-trip drift; displayed = {:?}",
            displayed
        );

        let allele = match &second {
            HgvsVariant::Allele(a) => a,
            other => panic!("expected Allele on reparse, got {:?}", other),
        };
        assert_eq!(allele.phase, AllelePhase::Trans);
        assert_eq!(allele.variants.len(), 2);
        assert!(
            matches!(&allele.variants[0], HgvsVariant::Protein(_)),
            "first arm must be Protein (carrying NoProtein); got {:?}",
            allele.variants[0]
        );
        match protein_edit(&allele.variants[0]) {
            ProteinEdit::NoProtein { predicted: false } => {}
            other => panic!(
                "first arm must carry ProteinEdit::NoProtein {{ predicted: false }}, got {:?}",
                other
            ),
        }
    }

    /// Predicted variant: `p.0?` must round-trip through the compact form as
    /// well. `[ACC:p.X];[ACC:p.0?]` Display-emits `ACC:p.[X];[0?]`. Inside
    /// the protein compact trans-allele, `[0?]` must resolve to
    /// `ProteinEdit::NoProtein { predicted: true }`.
    #[test]
    fn compact_trans_allele_zero_predicted_roundtrips() {
        let input = "[NP_003997.2:p.Arg97Trp];[NP_003997.2:p.0?]";
        let first = parse_hgvs(input).expect("parse expanded form with p.0?");
        let displayed = format!("{}", first);
        assert_eq!(displayed, "NP_003997.2:p.[Arg97Trp];[0?]");

        let second = parse_hgvs(&displayed).expect("reparse compact form with [0?]");
        assert_eq!(first, second, "p.0? round-trip drift");

        let allele = match &second {
            HgvsVariant::Allele(a) => a,
            other => panic!("expected Allele, got {:?}", other),
        };
        match protein_edit(&allele.variants[1]) {
            ProteinEdit::NoProtein { predicted: true } => {}
            other => panic!("expected NoProtein predicted, got {:?}", other),
        }
    }
}

/// Boundary behavior: forms that must NOT route to `NoProtein` — cis
/// rejection, bare cross-coordinate `[0]`, and non-protein coordinate
/// systems. Guards against the fix over-reaching.
mod boundary {
    use super::*;

    /// Cis form `p.[X;0]` is explicitly forbidden by the spec — `[0]` is
    /// only allowed in trans context (`alleles.md` line 48). The cis-bracket
    /// parser must reject it. This test pins the rejection.
    #[test]
    fn cis_bracket_zero_is_rejected() {
        let bad: &[&str] = &[
            "NP_003997.2:p.[0]",           // singleton cis-bracket
            "NP_003997.2:p.[Arg97Trp;0]",  // 0 second
            "NP_003997.2:p.[0;Arg97Trp]",  // 0 first
            "NP_003997.2:p.[Arg97Trp;0?]", // predicted 0 second
            "NP_003997.2:p.[0?;Arg97Trp]", // predicted 0 first
        ];
        for input in bad {
            let result = parse_hgvs(input);
            assert!(
                result.is_err(),
                "expected parse_hgvs({:?}) to reject cis-form `[0]` as malformed, but got: {:?}",
                input,
                result.ok()
            );
        }
    }

    /// The expanded form `[ACC:p.X];[0]` (note: bare `[0]` with NO `ACC:p.`
    /// prefix) is the cross-coordinate "absent allele" marker. The outer
    /// `parse_trans_allele` path must continue to route this to
    /// `HgvsVariant::NullAllele`, untouched by this fix.
    ///
    /// This matches the existing pin
    /// `protein_no_protein_roundtrip::bare_bracketed_zero_is_null_allele_not_no_protein`
    /// and is duplicated here to guard against accidentally regressing the
    /// boundary while changing the compact-form code path.
    #[test]
    fn expanded_form_bare_zero_remains_null_allele() {
        let parsed =
            parse_hgvs("[NP_003997.2:p.Arg97Trp];[0]").expect("parse expanded form with bare [0]");
        let allele = match &parsed {
            HgvsVariant::Allele(a) => a,
            other => panic!("expected Allele, got {:?}", other),
        };
        assert_eq!(allele.phase, AllelePhase::Trans);
        assert_eq!(allele.variants.len(), 2);
        assert!(
            matches!(&allele.variants[1], HgvsVariant::NullAllele),
            "bare expanded-form [0] must remain NullAllele, got {:?}",
            allele.variants[1]
        );
    }

    /// Other coordinate systems (`c.`, `g.`, `n.`, `r.`, `m.`, `o.`) do not
    /// have a `0` whole-entity no-product edit. Their compact trans-allele
    /// `[0]` (if accepted at all) must continue to mean `NullAllele`, never
    /// silently route to a protein-only `NoProtein` concept. This guards
    /// against the fix over-reaching outside the protein dispatch.
    ///
    /// For DNA/RNA/MT/circular, `[0]` in compact trans-allele context
    /// (whether in the parser or in the round-trip pin) means the
    /// absent-allele marker; `NoProtein` is protein-only.
    #[test]
    fn non_protein_compact_trans_allele_zero_stays_null_allele() {
        // Each case is a compact-trans-allele input on a non-protein
        // coordinate system. The bare-bracket second arm is `[0]`; the
        // result must contain `HgvsVariant::NullAllele` as the second
        // sub-variant.
        let cases: &[&str] = &[
            "NM_000088.3:c.[459A>G];[0]",
            "NC_000001.11:g.[100A>G];[0]",
            "NR_002196.2:n.[100A>G];[0]",
            "NR_002196.2:r.[100a>g];[0]",
            "NC_012920.1:m.[100A>G];[0]",
            "NC_000001.11:o.[100A>G];[0]",
        ];
        for input in cases {
            let parsed = parse_hgvs(input)
                .unwrap_or_else(|e| panic!("parse_hgvs({:?}) failed: {}", input, e));
            let allele = match &parsed {
                HgvsVariant::Allele(a) => a,
                other => panic!("{:?}: expected Allele, got {:?}", input, other),
            };
            assert_eq!(
                allele.phase,
                AllelePhase::Trans,
                "{:?}: must be Trans",
                input
            );
            assert_eq!(allele.variants.len(), 2, "{:?}: must have 2 arms", input);
            assert!(
                matches!(&allele.variants[1], HgvsVariant::NullAllele),
                "{:?}: bare [0] on non-protein coord must remain NullAllele, got {:?}",
                input,
                allele.variants[1]
            );
        }
    }
}
