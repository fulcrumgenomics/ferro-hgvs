//! `sup` supernumerary marker (#546). HGVS duplication.md:117 (whole-chromosome
//! supernumerary) and complex.md:161 (supernumerary ring chromosome).
use ferro_hgvs::parse_hgvs;

#[test]
fn supernumerary_whole_chromosome_round_trips() {
    let s = "NC_000023.11:g.pter_qtersup";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert_eq!(format!("{v}"), s);
}

#[test]
fn supernumerary_ring_round_trips() {
    let s = "NC_000022.11:g.[pter_(12200001_14700000)del::(37600001_410000000)_qterdel]sup";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert_eq!(format!("{v}"), s);
}

#[test]
#[allow(clippy::single_element_loop)] // kept as a table for future sup forms
fn supernumerary_with_gene_symbol_accession_round_trips() {
    for s in ["NC_000023.11(SOMEGENE):g.pter_qtersup"] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        assert_eq!(format!("{v}"), s);
    }
}

#[test]
fn sup_on_an_editful_bare_variant_is_rejected() {
    // `sup` applies only to the editless whole-chromosome form, not to an edit.
    assert!(parse_hgvs("NC_000023.11:g.12345A>Gsup").is_err());
    assert!(parse_hgvs("NC_000023.11:g.pter_qterdelsup").is_err());
}

#[test]
fn sup_on_non_genome_variant_is_rejected() {
    assert!(parse_hgvs("NM_000088.3:c.5A>Gsup").is_err());
}
