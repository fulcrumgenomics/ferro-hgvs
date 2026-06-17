//! Predicted / uncertain compound cis allele `[(a;b)]` (#545).
//!
//! `(...)` wraps a `;`-separated cis group inside the allele bracket, marking
//! the whole cis allele as predicted/uncertain — e.g.
//! `NP_003997.1:p.[(Ser68Arg;Asn594del)]`. Spec: `recommendations/uncertain.md`,
//! `recommendations/<axis>/alleles.md`.

use ferro_hgvs::hgvs::AllelePhase;
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, Normalizer};

#[test]
fn predicted_cis_allele_round_trips() {
    for s in [
        "NP_003997.1:p.[(Ser68Arg;Asn594del)]",
        "NP_003997.1:p.[(Ser73Arg;Asn103del)]",
        "LRG_199t1:r.[(578c>u;1339a>g;1680del)]",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        let HgvsVariant::Allele(a) = &v else {
            panic!("expected Allele for `{s}`, got {v:?}");
        };
        assert_eq!(a.phase, AllelePhase::Cis, "phase for `{s}`");
        assert!(a.uncertain, "uncertain flag for `{s}`");
        assert!(a.variants.len() >= 2, "member count for `{s}`");
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}

/// `p.[(Lys79*;Lys79Asn)]` parses as a predicted cis allele too; only the
/// `*` stop renders in its canonical `Ter` form (a pre-existing protein-stop
/// canonicalization, independent of the `(...)` wrapper). It moves off
/// `parse-error` to `diverges`, not `preserved`.
#[test]
fn predicted_cis_allele_with_stop_parses_ter_canonical() {
    let v = parse_hgvs("NP_003997.1:p.[(Lys79*;Lys79Asn)]").expect("must parse");
    let HgvsVariant::Allele(a) = &v else {
        panic!("expected Allele, got {v:?}");
    };
    assert_eq!(a.phase, AllelePhase::Cis);
    assert!(a.uncertain);
    assert_eq!(format!("{v}"), "NP_003997.1:p.[(Lys79Ter;Lys79Asn)]");
}

/// Normalize must be idempotent on predicted cis alleles: two consecutive
/// passes must produce identical output, and the `uncertain` flag and `Cis`
/// phase must survive both passes.  Mirrors `test_trans_normalize_is_idempotent`
/// in `allele_trans_phase.rs`.
///
/// This also pins the normalizer's singleton-merge branch: if normalization
/// were to collapse the two protein variants to one and then drop the
/// `uncertain` wrapper via the cis-unwrap path, `pass1.uncertain` would be
/// `false` and the Display would lose the `[(…)]` notation.
#[test]
fn predicted_cis_allele_normalize_is_idempotent() {
    let normalizer = Normalizer::new(MockProvider::new());
    let parsed = parse_hgvs("NP_003997.1:p.[(Ser68Arg;Asn594del)]").expect("must parse");
    let pass1 = normalizer
        .normalize(&parsed)
        .expect("normalize pass1 failed");
    let pass2 = normalizer
        .normalize(&pass1)
        .expect("normalize pass2 failed");
    assert_eq!(
        format!("{}", pass1),
        format!("{}", pass2),
        "normalize must be idempotent on predicted cis alleles",
    );
    let HgvsVariant::Allele(a) = pass2 else {
        panic!("expected Allele after second normalize");
    };
    assert_eq!(
        a.phase,
        AllelePhase::Cis,
        "phase must be Cis after normalize"
    );
    assert!(a.uncertain, "uncertain flag must survive normalize");
}

/// A single-member uncertain compound allele whose only member is an edit
/// containing a *nested* `;` (e.g. an insertion with a multi-segment inserted
/// allele `ins[A;T]`) must NOT be misclassified as a predicted cis allele.
///
/// `find_predicted_cis_bracket` requires a *top-level* `;` separator (depth 0
/// w.r.t. nested edit brackets / parens), not just any `;` anywhere in the
/// payload. The `;` in `ins[A;T]` is part of the inserted allele, at bracket
/// depth 1, so the input is not a `;`-separated cis group. Misclassifying it
/// reconstructs a single-member cis that collapses to a bare insertion,
/// dropping the `[(…)]` wrapper entirely (a round-trip failure).
///
/// Here the input is therefore a single-member compound allele that round-trips
/// to the bare insertion form, exactly as the non-predicted route already does
/// for `p.[(Ser68Arg)]` → `p.(Ser68Arg)`.
#[test]
fn nested_semicolon_in_edit_is_not_predicted_cis() {
    let s = "NM_004006.2:c.[(100_101ins[A;T])]";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    let HgvsVariant::Allele(a) = &v else {
        panic!("expected Allele for `{s}`, got {v:?}");
    };
    // Not a predicted cis allele: a single member whose nested `;` is part of
    // the inserted allele, not an allele-level cis separator.
    assert!(
        !(a.phase == AllelePhase::Cis && a.uncertain && a.variants.len() >= 2),
        "`{s}` must not be classified as a predicted cis allele, got {a:?}",
    );
    // Legit two-member predicted cis alleles still round-trip.
    for ok in [
        "NP_003997.1:p.[(Ser68Arg;Asn594del)]",
        "NM_004006.2:c.[(100A>G;200T>C)]",
    ] {
        let v = parse_hgvs(ok).unwrap_or_else(|e| panic!("must parse `{ok}`: {e}"));
        let HgvsVariant::Allele(a) = &v else {
            panic!("expected Allele for `{ok}`, got {v:?}");
        };
        assert_eq!(a.phase, AllelePhase::Cis, "phase for `{ok}`");
        assert!(a.uncertain, "uncertain flag for `{ok}`");
        assert_eq!(format!("{v}"), ok, "round-trip for `{ok}`");
    }
}

/// Predicted cis alleles round-trip correctly on the DNA coordinate axes
/// (`g.`, `c.`, `n.`).  All three axes route through the same
/// `find_predicted_cis_bracket` helper and must preserve the `[(…)]` wrapper.
#[test]
fn predicted_cis_allele_dna_axes_round_trip() {
    for s in [
        "NC_000001.11:g.[(100A>G;200T>C)]",
        "NM_000088.3:c.[(100A>G;200T>C)]",
        "NR_024540.1:n.[(100A>G;200T>C)]",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        let HgvsVariant::Allele(a) = &v else {
            panic!("expected Allele for `{s}`, got {v:?}");
        };
        assert_eq!(a.phase, AllelePhase::Cis, "phase for `{s}`");
        assert!(a.uncertain, "uncertain flag for `{s}`");
        assert!(a.variants.len() >= 2, "member count for `{s}`");
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}
