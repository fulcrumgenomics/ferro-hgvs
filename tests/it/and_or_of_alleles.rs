//! And/or `^` joining bracketed allele operands (#544).
//!
//! `p.[(Asn158Asp)(;)(Asn158Ile)]^[(Asn158Val)]` — an and/or (`^`) between
//! two alleles: an unknown-phase predicted-protein group and a single-member
//! predicted-protein bracket (which inherits the accession). HGVS general.md.
//! Combines `^` (#547), predicted-protein members (#552), unknown-phase, and
//! single-member bracket operands.

use ferro_hgvs::hgvs::AllelePhase;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

#[test]
fn and_or_of_bracketed_alleles_round_trips() {
    let s = "NP_003997.2:p.[(Asn158Asp)(;)(Asn158Ile)]^[(Asn158Val)]";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    let HgvsVariant::Allele(a) = &v else {
        panic!("expected Allele for `{s}`, got {v:?}");
    };
    assert_eq!(a.phase, AllelePhase::AndOr, "phase");
    assert_eq!(a.variants.len(), 2, "operand count");
    assert_eq!(format!("{v}"), s, "round-trip");
}

/// A `c.` and/or form round-trips: a two-member cis bracket followed by a
/// bare single-member bracket that inherits the accession.
#[test]
fn and_or_of_cds_bracketed_alleles_round_trips() {
    let s = "NM_000088.3:c.[100A>G;200T>C]^[300del]";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    let HgvsVariant::Allele(a) = &v else {
        panic!("expected Allele for `{s}`, got {v:?}");
    };
    assert_eq!(a.phase, AllelePhase::AndOr, "phase");
    assert_eq!(a.variants.len(), 2, "operand count");
    assert_eq!(format!("{v}"), s, "round-trip");
}

/// A bare bracketed operand whose inner edit carries its *own* inner accession
/// (a cross-reference `delins[ACC:g.X_Y]`) must still inherit the bracket's
/// shared prefix. The chunk contains a `:` only via that inner accession, not a
/// leading accession stem, so it must route through the inherited-prefix path
/// rather than being misrouted to the fully-qualified parser.
#[test]
fn and_or_bare_operand_with_inner_accession_round_trips() {
    let s = "NC_000022.10:g.[5_10del]^[10_15delins[NC_000022.10:g.20_25]]";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    let HgvsVariant::Allele(a) = &v else {
        panic!("expected Allele for `{s}`, got {v:?}");
    };
    assert_eq!(a.phase, AllelePhase::AndOr, "phase");
    assert_eq!(a.variants.len(), 2, "operand count");
    assert_eq!(format!("{v}"), s, "round-trip");
}

/// A bare bracketed operand appearing *first* in the `^` chain (no accession
/// on the left-hand side) must produce a parse error — the guard at
/// `parse_and_or_allele` requires at least one accession-bearing operand
/// before any bare operand can inherit from it.
#[test]
fn and_or_bare_operand_first_is_rejected() {
    assert!(
        parse_hgvs("[100A>G]^[200T>C]").is_err(),
        "bare first operand (no accession) must be rejected"
    );
}
