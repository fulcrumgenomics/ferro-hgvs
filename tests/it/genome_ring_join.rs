//! Same-chromosome ring `::` deletion-join (#546). HGVS DNA/complex.md:127 —
//! the surviving ISCN2020 meaning of `::` (ring-chromosome break junction).
//! A ring requires at least two segments (DNA/complex.md:130); a single segment
//! or trailing empty segment is rejected by the grammar.
use ferro_hgvs::hgvs::variant::AlleleVariant;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

#[test]
fn ring_deletion_join_round_trips() {
    let s = "NC_000022.11:g.pter_(12200001_14700000)del::(37600001_410000000)_qterdel";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
}

#[test]
fn cross_chromosome_double_colon_is_rejected() {
    // Cross-chromosome `::` (inner accession) is ISCN2016, removed in ISCN2020
    // and spec-invalid — must NOT parse as a ring.
    assert!(parse_hgvs("NC_000002.12:g.pter_8247756::NC_000011.10:g.15825273_cen_qter").is_err());
}

#[test]
fn three_segment_ring_round_trips() {
    // exercises the join loop past the first separator
    let s = "NC_000022.11:g.pter_1000del::2000_3000del::4000_qterdel";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert_eq!(format!("{v}"), s);
}

#[test]
fn trailing_empty_segment_is_rejected() {
    assert!(parse_hgvs("NC_000022.11:g.pter_1000del::").is_err());
}

#[test]
fn plain_genome_variant_without_double_colon_is_unaffected() {
    let s = "NC_000022.11:g.1000_2000del";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert_eq!(format!("{v}"), s);
}

#[test]
fn genome_ring_variant_type_is_g_ring() {
    // GenomeRing must report "g::g" (not "g") so it is distinct from a plain
    // GenomeVariant and excluded from the compact-allele path.
    let s = "NC_000022.11:g.pter_1000del::2000_qterdel";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert_eq!(
        v.variant_type(),
        "g::g",
        "GenomeRing must report variant_type \"g::g\""
    );
}

#[test]
fn cis_allele_of_genome_rings_uses_expanded_form() {
    // Two GenomeRing sub-variants that share accession must NOT collapse into
    // compact allele form `ACC:g.[seg1;seg2]` — the segments already contain
    // `::` separators, so compact form would produce a malformed string.
    // The expanded form `[ACC:g.seg1::seg2;ACC:g.seg3::seg4]` must be used.
    let r1 = parse_hgvs("NC_000022.11:g.pter_1000del::2000_qterdel")
        .unwrap_or_else(|e| panic!("must parse ring 1: {e}"));
    let r2 = parse_hgvs("NC_000022.11:g.pter_3000del::4000_qterdel")
        .unwrap_or_else(|e| panic!("must parse ring 2: {e}"));

    // Confirm sub-variant types
    assert_eq!(r1.variant_type(), "g::g");
    assert_eq!(r2.variant_type(), "g::g");

    // Build a programmatic cis allele and confirm it uses the expanded form.
    let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![r1, r2]));
    let rendered = format!("{allele}");

    // The expanded form emits each sub-variant with its full `ACC:g.` prefix
    // inside the outer brackets.  A compact rendering would have only one
    // prefix before the `[` — check that both accessions appear.
    assert!(
        rendered.starts_with('['),
        "cis allele of rings must use expanded form (starts with `[`); got: {rendered}"
    );
    assert_eq!(
        rendered.matches("NC_000022.11:g.").count(),
        2,
        "expanded form must repeat the accession prefix for each sub-variant; got: {rendered}"
    );
}
