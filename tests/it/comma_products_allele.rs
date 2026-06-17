//! Comma allele `[a,b]` — different transcripts/proteins derived from one
//! allele (#545). Spec: `general.md`, `RNA/splicing.md`, `RNA/alleles.md`
//! ("two different transcripts ... derive from one variant on the DNA level").
//!
//! Distinct from cis (`;`), trans (`];[`), and-or (`^`), mosaic/chimeric
//! (`/`/`//`): the comma members are downstream products of a single DNA
//! event, modeled as `AllelePhase::Products`.

use ferro_hgvs::hgvs::AllelePhase;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

#[test]
fn comma_products_allele_round_trips() {
    // Accession-bearing forms. The bare-prefix examples (`r.[123a>u,...]`)
    // depend on bare-RNA-without-accession member parsing, which ferro does
    // not support — a limitation orthogonal to the comma relationship.
    for s in [
        "LRG_199t1:r.[897u>g,832_960del]",
        "NM_004006.3:r.[897u>g,832_960del]",
        "NP_003997.1:p.[Lys31Asn,Val25_Lys31del]",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        let HgvsVariant::Allele(a) = &v else {
            panic!("expected Allele for `{s}`, got {v:?}");
        };
        assert_eq!(a.phase, AllelePhase::Products, "phase for `{s}`");
        assert!(!a.uncertain, "must not be uncertain for `{s}`");
        assert_eq!(a.variants.len(), 2, "member count for `{s}`");
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}

/// Predicted comma allele `[(a,b)]`: the `(...)` wraps all members inside the
/// bracket (whole-group uncertainty) — `[(` ... `)]`.
#[test]
fn predicted_comma_products_allele_round_trips() {
    let s = "LRG_199t1:r.[(2603_2622del,2622g>c)]";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    let HgvsVariant::Allele(a) = &v else {
        panic!("expected Allele for `{s}`, got {v:?}");
    };
    assert_eq!(a.phase, AllelePhase::Products);
    assert!(a.uncertain, "predicted form must set uncertain");
    assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
}

/// A bracket may not mix the comma (product) separator with `;` (cis) —
/// they answer orthogonal questions, and the grammar admits only one
/// separator kind per bracket.
#[test]
fn mixed_comma_and_semicolon_separators_reject() {
    assert!(
        parse_hgvs("LRG_199t1:r.[897u>g,832_960del;900a>c]").is_err(),
        "mixing `,` and `;` in one bracket must reject"
    );
}

/// The comma (products) relationship is r./p.-only — different
/// transcripts/proteins derived from one allele. It is never valid on a
/// g./c. DNA-source description, so `c.[a,b]` / `g.[a,b]` must reject. This
/// guards the `find_products_bracket` axis restriction that keeps real
/// ClinVar/CMRG `c.[a,b]` inputs out of the products path. HGVS general.md,
/// RNA/alleles.md.
#[test]
fn comma_on_dna_axes_reject() {
    assert!(
        parse_hgvs("NM_014762.3:c.[881A>C,918G>C]").is_err(),
        "comma allele on the c. (CDS) axis must reject"
    );
    assert!(
        parse_hgvs("NC_000001.11:g.[100A>G,200C>T]").is_err(),
        "comma allele on the g. (genomic) axis must reject"
    );
}

/// A products allele requires an outer `r.`/`p.` coordinate stem; the
/// fully-expanded cross-accession bracket form `[ACC1:r.a,ACC2:r.b]` (no
/// shared outer stem) is out of spec scope and must reject rather than parse
/// into a shape that would not round-trip.
#[test]
fn bare_expanded_cross_accession_products_reject() {
    assert!(
        parse_hgvs("[NM_004006.3:r.897u>g,NM_004006.3:r.832_960del]").is_err(),
        "fully-expanded bare-bracket products form must reject"
    );
}
