//! Compound-reference accession pairing validation (#963).
//!
//! A compound reference `outer(inner)` names the annotated transcript/feature
//! (`inner`) in the context of a genomic/gene reference (`outer`), e.g.
//! `NC_000013.11(NM_004119.3):c.…` or `NG_012232.1(NM_004006.2):c.…`. The
//! parser previously accepted *any* pairing, including backwards forms where
//! the outer is a transcript or protein. These tests pin the end-to-end
//! behavior at the public `parse_hgvs` API:
//!
//! - backwards pairings (transcript/protein outer) are rejected, for both
//!   RefSeq and Ensembl (symmetric — the Ensembl compound branch was a faithful
//!   mirror of the permissive RefSeq one, PR #938);
//! - genuine genomic/gene-outer compounds still parse;
//! - the *inner* is intentionally left lenient — non-standard but real forms
//!   such as mutalyzer's `NG_x(NP_y)` protein wrapper must still parse so ferro
//!   can re-express them as the spec-correct bare form (see the
//!   `bare-np-protein` conformance cluster).

use ferro_hgvs::parse_hgvs;

/// Backwards compound references — the outer is a transcript or protein — are
/// rejected at the public parse API, for both RefSeq and Ensembl.
#[test]
fn backwards_compound_reference_is_rejected() {
    let backwards = [
        // RefSeq: transcript / non-coding / protein outer.
        "NM_004119.3(NG_012232.1):c.100A>G",
        "NR_046018.2(NM_000088.3):n.100A>G",
        "NP_004110.2(NM_000088.3):p.Val600Glu",
        // Ensembl: transcript / protein outer.
        "ENST00000375549.8(ENSG00000204370.13):c.100A>G",
        "ENSP00000012345.1(ENST00000375549.8):c.100A>G",
    ];
    for input in backwards {
        assert!(
            parse_hgvs(input).is_err(),
            "expected backwards compound reference to be rejected: {input}"
        );
    }
}

/// Genuine genomic/gene-outer compounds are unaffected and still parse.
#[test]
fn genomic_gene_outer_compound_reference_still_parses() {
    let valid = [
        "NC_000013.11(NM_004119.3):c.100A>G",
        "NG_007400.1(NM_000088.3):c.459A>G",
        "NC_000001.11(NR_046018.2):n.100A>G",
        "ENSG00000204370.13(ENST00000375549.8):c.100A>G",
        "LRG_1(NM_000088.3):c.459A>G",
    ];
    for input in valid {
        assert!(
            parse_hgvs(input).is_ok(),
            "expected genomic/gene-outer compound to parse: {input}"
        );
    }
}

/// The inner is intentionally lenient: a genomic/gene outer with a protein or
/// bare-genomic inner still parses (mutalyzer emits `NG_x(NP_y)` wrappers and
/// ferro must parse them to normalize to the spec-correct bare form).
#[test]
fn compound_reference_inner_stays_lenient() {
    let lenient = [
        "NC_000013.11(NP_004110.2):p.Val600Glu",
        "NG_007485.1(NP_478102.2):p.Asp68Glu",
        "NC_000023.11(LRG_199):g.1A>G",
    ];
    for input in lenient {
        assert!(
            parse_hgvs(input).is_ok(),
            "expected lenient inner compound to parse: {input}"
        );
    }
}
