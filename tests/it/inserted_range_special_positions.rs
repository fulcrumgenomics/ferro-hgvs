//! pter/qter endpoints in an inserted position-range with `inv` (#546).
//! HGVS DNA/complex.md:94-95 — homologous-chromosome translocation as delins.
use ferro_hgvs::parse_hgvs;

#[test]
fn delins_inserted_range_with_special_endpoints_round_trips() {
    for s in [
        "NC_000009.12:g.pter_26393001delins102425452_qterinv",
        "NC_000009.12:g.102425452_qterdelinspter_26393001inv",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}

#[test]
fn special_positions_in_inserted_range_are_rejected_outside_genome_axis() {
    // pter/qter are genome-only; a c. inserted range must not accept them.
    assert!(parse_hgvs("LRG_24t1:c.pter_qterdelinspter_qter").is_err());
}

#[test]
fn all_numeric_inserted_range_is_unchanged() {
    // Regression: numeric inserted ranges must keep their existing behavior.
    let s = "NC_000009.12:g.100_200delins850_900inv";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert_eq!(format!("{v}"), s);
}

#[test]
fn cen_endpoint_in_inserted_range_round_trips() {
    // The parser claims `cen` (centromere) support as a special inserted endpoint;
    // exercise both the inverted and the bare (`inverted: false`) Display branches.
    for s in [
        "NC_000009.12:g.cen_26393001delins102425452_qterinv",
        "NC_000009.12:g.100_200delins26393001_cen",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}

#[test]
fn special_inserted_range_round_trips_on_mito_and_circular_axes() {
    // The gate is the genome *family* (g./m./o.), not just `g.` — special-position
    // inserted ranges must parse on the mitochondrial (m.) and circular (o.) axes too.
    for s in [
        "NC_012920.1:m.100_200delins102425452_qterinv",
        "NC_012920.1:m.pter_200delins102425452_qter",
        "NC_012920.1:o.100_200delins102425452_qterinv",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}

#[test]
fn special_positions_in_inserted_range_are_rejected_on_rna_axis() {
    // Mirror the c. rejection: pter/qter inserted ranges are genome-only and must
    // be refused on the r. (RNA) axis as well.
    assert!(parse_hgvs("NM_004006.2:r.100_200delinspter_qter").is_err());
    assert!(parse_hgvs("NM_004006.2:r.pter_200delins102425452_qter").is_err());
}

#[test]
fn non_inverted_special_inserted_range_round_trips() {
    // Exercise the `inverted: false` Display branch on the genome axis (no trailing `inv`).
    let s = "NC_000009.12:g.102425452_qterdelinspter_26393001";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert_eq!(format!("{v}"), s);
}

#[test]
fn ins_inserted_range_with_special_endpoints_round_trips_on_genome_family() {
    // `allow_special_positions` is threaded into `parse_insertion` as well as
    // `parse_delins`, so the bare `ins<range>` form must accept special-position
    // inserted ranges (inverted and non-inverted) on the whole genome family
    // (g./m./o.), exactly like `delins`.
    for s in [
        "NC_000009.12:g.100_101inspter_qterinv",
        "NC_000009.12:g.100_101inspter_qter",
        "NC_000009.12:g.100_101ins102425452_qterinv",
        "NC_012920.1:m.100_101inspter_qterinv",
        "NC_012920.1:o.100_101inspter_qter",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}

#[test]
fn ins_special_positions_in_inserted_range_are_rejected_outside_genome_axis() {
    // Mirror the delins rejection: pter/qter inserted ranges are genome-only and
    // must be refused for the `ins` form on c./r. axes too.
    assert!(parse_hgvs("LRG_24t1:c.100_101inspter_qter").is_err());
    assert!(parse_hgvs("NM_004006.2:r.100_101inspter_qter").is_err());
}

#[test]
fn dupins_form_is_rejected_even_with_special_positions() {
    // `parse_dupins` threads `allow_special_positions`, but the `dupins<seq>` form
    // itself is not used in HGVS nomenclature (DNA/duplication.md:92) and is
    // rejected regardless of the inserted payload — special positions must not
    // provide a backdoor that lets `dupins` through.
    assert!(parse_hgvs("NC_000009.12:g.100_200dupinspter_qterinv").is_err());
    assert!(parse_hgvs("NC_000009.12:g.100_200dupins102425452_qter").is_err());
}

#[test]
fn inherited_accession_slash_form_accepts_special_inserted_range_on_genome_family() {
    // The slash mosaic form whose RHS omits the accession reconstructs the RHS as
    // a genome-family variant; that fallback path must accept special-position
    // inserted ranges (g./m./o.) consistently with the primary parse paths, and
    // still reject them on c./r.
    for s in [
        "NC_000009.12:g.50A>T/g.100_101inspter_qter",
        "NC_000009.12:g.50A>T/g.pter_26393001delins102425452_qterinv",
        "NC_012920.1:m.50A>T/m.100_101inspter_qter",
        "NC_012920.1:o.50A>T/o.100_101inspter_qter",
    ] {
        assert!(
            parse_hgvs(s).is_ok(),
            "genome-family inherited-accession must accept `{s}`"
        );
    }
    assert!(parse_hgvs("NM_004006.2:c.50A>T/c.100_101inspter_qter").is_err());
}
