//! Issue #487 (one group): a genomic (`g.`) variant on the human
//! mitochondrial reference must normalize to mitochondrial (`m.`) notation.
//!
//! `NC_012920.1` (GRCh38 rCRS) and the deprecated `NC_001807` are mtDNA
//! references; HGVS requires `m.` for them. ferro parses `NC_012920.1:g.…` as
//! a genomic variant (preserved — see `tests/issue_261_*`), but `normalize()`
//! coerces it to the spec-required `m.` form.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

fn norm(input: &str) -> String {
    let normalizer = Normalizer::new(MockProvider::new());
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    format!(
        "{}",
        normalizer
            .normalize(&variant)
            .unwrap_or_else(|e| panic!("normalize {input}: {e}"))
    )
}

#[test]
fn genomic_substitution_on_mito_accession_coerces_to_m() {
    assert_eq!(norm("NC_012920.1:g.3243A>G"), "NC_012920.1:m.3243A>G");
}

#[test]
fn deprecated_mito_accession_also_coerces() {
    assert_eq!(norm("NC_001807.4:g.8993T>G"), "NC_001807.4:m.8993T>G");
}

#[test]
fn non_mito_genomic_accession_stays_g() {
    // A normal chromosomal accession must NOT be coerced.
    assert_eq!(norm("NC_000001.11:g.1000A>G"), "NC_000001.11:g.1000A>G");
}

#[test]
fn mito_input_already_m_is_unchanged() {
    assert_eq!(norm("NC_012920.1:m.3243A>G"), "NC_012920.1:m.3243A>G");
}
