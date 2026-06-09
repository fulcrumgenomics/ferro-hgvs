//! Tests for trans-phase allele canonicalization.
//!
//! Tracking: issue #81 item C1 — "Trans-phase canonicalization — `[a];[b]` shape
//! vs malformed inputs. Round-trip fidelity unverified."
//!
//! Per the HGVS DNA alleles recommendation
//! (`assets/hgvs-nomenclature/docs/recommendations/DNA/alleles.md`), trans phase
//! is written as `[variant1];[variant2]` (each variant in its own bracket pair,
//! separated by `;`). Cis phase, by contrast, is `[variant1;variant2]` (a
//! single bracket pair containing semicolon-separated variants). The two shapes
//! describe biologically distinct genotypes and must not be confused by parser
//! or display.
//!
//! These tests pin:
//! 1. Round-trip fidelity for `[a];[b]` (both expanded and compact-prefix forms).
//! 2. That `[a;b]` (cis) does not collapse onto, nor parse as, `[a];[b]` (trans).
//! 3. That the `Cis`-only consecutive-edit merge pass (issue #72 / PR #80) is
//!    correctly excluded for trans alleles, even when the sub-variant positions
//!    are strictly adjacent (the case that would merge under cis).
//! 4. Mixed-accession trans round-trip (related to C4: mixed-accession alleles).
//!
//! See also: C2 (mosaic / chimeric), C3 (unknown-phase `[a(;)b]`),
//! C4 (mixed-accession alleles), C5 (3+ variants per allele).

use ferro_hgvs::hgvs::AllelePhase;
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, Normalizer};

fn parse_to_string(input: &str) -> String {
    let variant = parse_hgvs(input).expect("parse failed");
    format!("{}", variant)
}

fn normalize_to_string(input: &str) -> String {
    let normalizer = Normalizer::new(MockProvider::new());
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

fn parse_phase(input: &str) -> AllelePhase {
    match parse_hgvs(input).expect("parse failed") {
        HgvsVariant::Allele(a) => a.phase,
        other => panic!("expected Allele variant, got {:?}", other),
    }
}

// ---------------------------------------------------------------------------
// Round-trip fidelity for `[a];[b]`
// ---------------------------------------------------------------------------

/// `[ACC:c.a];[ACC:c.b]` (expanded form) round-trips through parse + display.
///
/// Same-accession trans alleles render in the spec's compact-prefix form
/// `ACC:c.[a];[b]`, so the second-pass round-trip is bit-identical.
#[test]
fn test_trans_expanded_form_round_trips_via_compact() {
    // First pass: expanded → compact (same accession, so compact-prefix form).
    assert_eq!(
        parse_to_string("[NM_000088.3:c.100A>G];[NM_000088.3:c.200T>C]"),
        "NM_000088.3:c.[100A>G];[200T>C]",
    );
    // Second pass: compact → compact (bit-identical, idempotent).
    assert_eq!(
        parse_to_string("NM_000088.3:c.[100A>G];[200T>C]"),
        "NM_000088.3:c.[100A>G];[200T>C]",
    );
    // Phase is Trans for both.
    assert_eq!(
        parse_phase("[NM_000088.3:c.100A>G];[NM_000088.3:c.200T>C]"),
        AllelePhase::Trans,
    );
    assert_eq!(
        parse_phase("NM_000088.3:c.[100A>G];[200T>C]"),
        AllelePhase::Trans,
    );
}

/// Trans round-trip for the spec's `LRG_199t1:c.[2376G>C];[3103del]` example.
///
/// Source: `assets/hgvs-nomenclature/docs/recommendations/DNA/alleles.md`,
/// "Variants in trans" section.
#[test]
fn test_trans_spec_example_round_trip() {
    assert_eq!(
        parse_to_string("NM_004006.2:c.[2376G>C];[3103del]"),
        "NM_004006.2:c.[2376G>C];[3103del]",
    );
    assert_eq!(
        parse_phase("NM_004006.2:c.[2376G>C];[3103del]"),
        AllelePhase::Trans,
    );
}

/// Trans alleles whose members are themselves multi-variant CIS groups:
/// `[a;b;c];[d;e]` (HGVS DNA alleles.md, compound heterozygote).
#[test]
fn test_trans_of_cis_groups_round_trips() {
    let s = "NM_004006.2:c.[296T>G;476T>C;1083A>C];[296T>G;1083A>C]";
    assert_eq!(parse_to_string(s), s, "round-trip for `{s}`");
    assert_eq!(parse_phase(s), AllelePhase::Trans, "phase for `{s}`");
}

/// The spec (DNA/alleles.md) marks the redundant cross-spelled `=` form as
/// invalid: "do not use `c.[2376G>C;3103=];[2376=;3103del]`". A multi-variant
/// cis group containing a position-identity (`=`) member must stay rejected
/// rather than round-trip the discouraged form.
#[test]
fn test_trans_of_cis_with_identity_member_is_rejected() {
    for s in [
        "LRG_199t1:c.[2376G>C;3103=];[2376=;3103del]",
        "NM_004006.2:c.[2376G>C;3103=];[2376=;3103del]",
    ] {
        assert!(
            parse_hgvs(s).is_err(),
            "spec-invalid cross-spelled `=` trans form must reject: `{s}`"
        );
    }
}

/// Inner members of a trans-of-cis allele are nested cis sub-alleles.
#[test]
fn test_trans_of_cis_nesting_shape() {
    let v = parse_hgvs("NM_004006.2:c.[296T>G;476T>C;1083A>C];[296T>G;1083A>C]").unwrap();
    let HgvsVariant::Allele(trans) = &v else {
        panic!("expected Allele, got {v:?}");
    };
    assert_eq!(trans.phase, AllelePhase::Trans);
    assert_eq!(trans.variants.len(), 2);
    let HgvsVariant::Allele(cis0) = &trans.variants[0] else {
        panic!(
            "expected nested cis allele member, got {:?}",
            trans.variants[0]
        );
    };
    assert_eq!(cis0.phase, AllelePhase::Cis);
    assert_eq!(cis0.variants.len(), 3);
}

/// Genomic trans round-trip — exercises the `g.` parser path.
#[test]
fn test_trans_genomic_round_trip() {
    assert_eq!(
        parse_to_string("[NC_000001.11:g.1000G>A];[NC_000001.11:g.2000T>C]"),
        "NC_000001.11:g.[1000G>A];[2000T>C]",
    );
    assert_eq!(
        parse_to_string("NC_000001.11:g.[1000G>A];[2000T>C]"),
        "NC_000001.11:g.[1000G>A];[2000T>C]",
    );
}

// ---------------------------------------------------------------------------
// Cis vs trans: distinct shapes, distinct phases
// ---------------------------------------------------------------------------

/// `[a;b]` (cis) and `[a];[b]` (trans) must not be confused.
///
/// They describe biologically distinct genotypes — cis = both variants on the
/// same chromosome, trans = one variant per chromosome (compound heterozygote).
/// The parser must classify them by shape, and the round-trip must preserve
/// the shape.
#[test]
fn test_cis_and_trans_are_not_confused() {
    let cis_input = "[NM_000088.3:c.100A>G;NM_000088.3:c.200T>C]";
    let trans_input = "[NM_000088.3:c.100A>G];[NM_000088.3:c.200T>C]";

    // Phase is correctly distinguished.
    assert_eq!(parse_phase(cis_input), AllelePhase::Cis);
    assert_eq!(parse_phase(trans_input), AllelePhase::Trans);

    // The two parse to non-equal HgvsVariant values.
    let cis_parsed = parse_hgvs(cis_input).expect("parse cis");
    let trans_parsed = parse_hgvs(trans_input).expect("parse trans");
    assert_ne!(
        cis_parsed, trans_parsed,
        "cis [a;b] and trans [a];[b] must not compare equal",
    );

    // Display preserves the shape: cis collapses inside one bracket pair,
    // trans keeps the bracket pair around each sub-variant.
    assert_eq!(
        format!("{}", cis_parsed),
        "NM_000088.3:c.[100A>G;200T>C]",
        "cis must render as [a;b]",
    );
    assert_eq!(
        format!("{}", trans_parsed),
        "NM_000088.3:c.[100A>G];[200T>C]",
        "trans must render as [a];[b]",
    );
}

// ---------------------------------------------------------------------------
// Trans is excluded from the Cis-only consecutive-edit merge (PR #80)
// ---------------------------------------------------------------------------

/// PR #80 / issue #72 introduced consecutive-edit merging for cis alleles
/// (e.g. `g.[1000G>A;1001A>C]` → `g.1000_1001delinsAC`). The merge is gated on
/// `AllelePhase::Cis` in `src/normalize/merge.rs`. This test pins that gating:
/// strictly adjacent variants in trans must NOT merge — they describe two
/// independent single-nucleotide events on opposite chromosomes, not one
/// dinucleotide event.
#[test]
fn test_trans_does_not_trigger_consecutive_edit_merge() {
    // Adjacent g. SNVs in trans — under cis these would merge to delins.
    let normalized = normalize_to_string("[NC_000001.11:g.1000G>A];[NC_000001.11:g.1001A>C]");
    assert!(
        !normalized.contains("delins"),
        "trans [1000G>A];[1001A>C] must not merge to delins, got: {}",
        normalized,
    );
    assert!(
        normalized.contains("1000G>A"),
        "first sub-variant must survive, got: {}",
        normalized,
    );
    assert!(
        normalized.contains("1001A>C"),
        "second sub-variant must survive, got: {}",
        normalized,
    );

    // Phase remains Trans after normalization (no unwrap to a single Genome).
    let normalizer = Normalizer::new(MockProvider::new());
    let parsed = parse_hgvs("[NC_000001.11:g.1000G>A];[NC_000001.11:g.1001A>C]").unwrap();
    let normalized_variant = normalizer.normalize(&parsed).unwrap();
    match normalized_variant {
        HgvsVariant::Allele(a) => assert_eq!(
            a.phase,
            AllelePhase::Trans,
            "phase must remain Trans after normalize",
        ),
        other => panic!("expected Allele after normalize, got {:?}", other),
    }
}

/// Sanity-paired with the test above: the same adjacent positions in cis DO
/// merge to delins under PR #80. If this assertion ever fails, the cross-check
/// in `test_trans_does_not_trigger_consecutive_edit_merge` becomes vacuous and
/// must be re-examined.
#[test]
fn test_cis_counterpart_does_merge_to_delins() {
    let normalized = normalize_to_string("NC_000001.11:g.[1000G>A;1001A>C]");
    assert_eq!(
        normalized, "NC_000001.11:g.1000_1001delinsAC",
        "cis adjacent SNVs must merge under PR #80; if not, the trans \
         exclusion test is vacuous",
    );
}

// ---------------------------------------------------------------------------
// Mixed-accession trans (related to C4)
// ---------------------------------------------------------------------------

/// A trans allele whose sub-variants share an accession renders in compact-
/// prefix form `ACC:c.[a];[b]`. Synthetic short accession exercises that path.
/// (See also C4 — mixed-accession alleles — for the cross-accession case.)
#[test]
fn test_mixed_accession_trans_same_accession_compact_round_trip() {
    // Both sub-variants share the synthetic accession `NM_X.1`, so the parser
    // accepts the expanded form and the display emits the compact-prefix form.
    let expanded = "[NM_X.1:c.1A>G];[NM_X.1:c.2T>C]";
    let compact = "NM_X.1:c.[1A>G];[2T>C]";
    assert_eq!(parse_to_string(expanded), compact);
    assert_eq!(parse_to_string(compact), compact);
    assert_eq!(parse_phase(expanded), AllelePhase::Trans);
}

/// When sub-variants come from different accessions or different coordinate
/// systems, the compact-prefix form does not apply: the expanded form
/// `[ACC1:...];[ACC2:...]` must round-trip bit-identically.
#[test]
fn test_mixed_accession_trans_distinct_accessions_round_trip() {
    let input = "[NC_000001.11:g.100A>G];[NM_000088.3:c.200T>C]";
    assert_eq!(parse_to_string(input), input);
    assert_eq!(parse_phase(input), AllelePhase::Trans);
}

// ---------------------------------------------------------------------------
// Malformed inputs — pin current rejection behavior
// ---------------------------------------------------------------------------

/// Whitespace inside the trans separator (`[a]; [b]`, `[a] ;[b]`,
/// `[a] ; [b]`) is currently REJECTED by the parser. The HGVS spec is silent
/// on whether human-readable whitespace should be tolerated, but pragmatic
/// implementations typically trim. This test pins ferro's current strict
/// behavior so any future relaxation is an explicit, reviewed change rather
/// than an accidental regression in either direction.
///
/// TODO(#81 C1 follow-up): consider tolerating cosmetic whitespace around the
/// `];` and `;[` separators in trans-allele inputs (the spec answer per the
/// issue description is "trim whitespace, accept"). Filing a follow-up issue
/// requires user approval for the title/body, so this is left as a TODO here.
#[test]
fn test_trans_with_whitespace_around_separator_is_rejected() {
    let cases = [
        "[NM_000088.3:c.100A>G]; [NM_000088.3:c.200T>C]",
        "[NM_000088.3:c.100A>G] ;[NM_000088.3:c.200T>C]",
        "[NM_000088.3:c.100A>G] ; [NM_000088.3:c.200T>C]",
        "NM_000088.3:c.[100A>G]; [200T>C]",
        "NM_000088.3:c.[100A>G] ;[200T>C]",
        "NM_000088.3:c.[100A>G] ; [200T>C]",
    ];
    for input in cases {
        assert!(
            parse_hgvs(input).is_err(),
            "expected parse error for whitespace-padded trans separator, \
             but parse succeeded: {:?}",
            input,
        );
    }
}

/// Unbalanced brackets and obviously-malformed shapes must be rejected.
#[test]
fn test_trans_obvious_malformed_shapes_rejected() {
    let cases = [
        // Missing closing bracket on second variant.
        "[NM_000088.3:c.100A>G];[NM_000088.3:c.200T>C",
        // Missing opening bracket on second variant.
        "[NM_000088.3:c.100A>G];NM_000088.3:c.200T>C]",
        // Empty trans bracket pair.
        "[];[NM_000088.3:c.100A>G]",
        // Trailing garbage after the trans allele.
        "[NM_000088.3:c.100A>G];[NM_000088.3:c.200T>C]xyz",
    ];
    for input in cases {
        assert!(
            parse_hgvs(input).is_err(),
            "expected parse error for malformed trans input, but parse \
             succeeded: {:?}",
            input,
        );
    }
}

// ---------------------------------------------------------------------------
// n. (non-coding transcript) — compact-prefix trans round-trip
// ---------------------------------------------------------------------------

/// Compact-prefix `n.[a];[b]` round-trips through parse + display.
///
/// Asymmetry fix paralleling PR #146's genomic case: the expanded form
/// `[ACC:n.a];[ACC:n.b]` already parsed via the top-level allele path, but
/// `Display` emits the compact form `ACC:n.[a];[b]` for same-accession
/// alleles, so without a compact-form parser path the `Display` → parse
/// round-trip was broken for `n.`.
#[test]
fn test_trans_n_compact_round_trip() {
    // First pass: expanded → compact (Display canonicalises same-acc trans).
    assert_eq!(
        parse_to_string("[NR_X.1:n.100A>G];[NR_X.1:n.200T>C]"),
        "NR_X.1:n.[100A>G];[200T>C]",
    );
    // Second pass: compact → compact (idempotent; was a parse error pre-fix).
    assert_eq!(
        parse_to_string("NR_X.1:n.[100A>G];[200T>C]"),
        "NR_X.1:n.[100A>G];[200T>C]",
    );
    assert_eq!(
        parse_phase("NR_X.1:n.[100A>G];[200T>C]"),
        AllelePhase::Trans,
    );
}

// ---------------------------------------------------------------------------
// m. (mitochondrial) — compact-prefix trans round-trip
// ---------------------------------------------------------------------------

/// Compact-prefix `m.[a];[b]` round-trips through parse + display.
///
/// Mitochondrial heteroplasmy is a canonical use case for trans (different
/// alleles in different mtDNA copies). Same Display → parse asymmetry as
/// `g.`/`n.`: fixed by paralleling `parse_genome_trans_allele_shorthand`
/// for `m.`.
#[test]
fn test_trans_m_compact_round_trip() {
    assert_eq!(
        parse_to_string("[NC_012920.1:m.100A>G];[NC_012920.1:m.200T>C]"),
        "NC_012920.1:m.[100A>G];[200T>C]",
    );
    assert_eq!(
        parse_to_string("NC_012920.1:m.[100A>G];[200T>C]"),
        "NC_012920.1:m.[100A>G];[200T>C]",
    );
    assert_eq!(
        parse_phase("NC_012920.1:m.[100A>G];[200T>C]"),
        AllelePhase::Trans,
    );
}

// ---------------------------------------------------------------------------
// o. (circular DNA, SVD-WG006) — compact-prefix trans round-trip
// ---------------------------------------------------------------------------

/// Compact-prefix `o.[a];[b]` round-trips through parse + display.
///
/// Circular DNA inherits the DNA allele grammar per SVD-WG006. Same
/// Display → parse asymmetry as `g.`/`n.`/`m.` fixed here.
#[test]
fn test_trans_o_compact_round_trip() {
    assert_eq!(
        parse_to_string("[AC_X.1:o.100A>G];[AC_X.1:o.200T>C]"),
        "AC_X.1:o.[100A>G];[200T>C]",
    );
    assert_eq!(
        parse_to_string("AC_X.1:o.[100A>G];[200T>C]"),
        "AC_X.1:o.[100A>G];[200T>C]",
    );
    assert_eq!(
        parse_phase("AC_X.1:o.[100A>G];[200T>C]"),
        AllelePhase::Trans,
    );
}

// ---------------------------------------------------------------------------
// p. (protein) — compact-prefix trans round-trip
// ---------------------------------------------------------------------------

/// Compact-prefix `p.[a];[b]` round-trips through parse + display.
///
/// Protein trans alleles use the same shape per
/// `assets/hgvs-nomenclature/docs/recommendations/protein/alleles.md`
/// (compound heterozygote). Pre-fix the parser only handled the cis form
/// `p.[a;b]`; the corresponding `parse_protein_trans_allele_shorthand`
/// added here parallels the DNA helpers.
#[test]
fn test_trans_p_compact_round_trip() {
    assert_eq!(
        parse_to_string("NP_003997.1:p.[Ser68Arg];[Asn594del]"),
        "NP_003997.1:p.[Ser68Arg];[Asn594del]",
    );
    assert_eq!(
        parse_to_string("[NP_003997.1:p.Ser68Arg];[NP_003997.1:p.Asn594del]"),
        "NP_003997.1:p.[Ser68Arg];[Asn594del]",
    );
    assert_eq!(
        parse_phase("NP_003997.1:p.[Ser68Arg];[Asn594del]"),
        AllelePhase::Trans,
    );
}

/// Spec line 38 (DNA alleles, applied to protein per the analogous protein
/// alleles section): homozygous trans (same variant on both alleles); both
/// chromosomes carry the same change. Phase still distinguishes from cis.
#[test]
fn test_trans_p_homozygous_round_trip() {
    assert_eq!(
        parse_to_string("NP_003997.1:p.[Ser68Arg];[Ser68Arg]"),
        "NP_003997.1:p.[Ser68Arg];[Ser68Arg]",
    );
    assert_eq!(
        parse_phase("NP_003997.1:p.[Ser68Arg];[Ser68Arg]"),
        AllelePhase::Trans,
    );
}

// ---------------------------------------------------------------------------
// Trans with special second-allele tokens: [?] (unknown), [0] (null),
// [pos=] (wild-type / no-change). All are spec-mandated forms.
// ---------------------------------------------------------------------------

/// `c.[a];[?]`: one allele has a variant, the other allele is unknown.
/// Spec: `assets/hgvs-nomenclature/docs/recommendations/DNA/alleles.md`,
/// "Variants in trans" section. Pre-existing functionality; pinned to catch
/// regressions and to confirm `[?]` lowers to the dedicated `UnknownAllele`
/// sentinel rather than to a coord-system variant (which would lose the
/// "second allele unknown" semantics).
#[test]
fn test_trans_c_with_unknown_allele_round_trip() {
    let input = "[NM_000088.3:c.100A>G];[?]";
    assert_eq!(parse_to_string(input), input);
    assert_eq!(parse_phase(input), AllelePhase::Trans);
    let parsed = parse_hgvs(input).expect("parse failed");
    if let HgvsVariant::Allele(a) = parsed {
        assert_eq!(a.variants.len(), 2);
        assert!(matches!(a.variants[1], HgvsVariant::UnknownAllele));
    } else {
        panic!("expected Allele");
    }
}

/// `c.[a];[0]`: one allele has a variant, the other allele is absent
/// (canonically used for hemizygous males on X-chromosome variants per
/// the "X-chromosome" Discussion note in
/// `recommendations/DNA/alleles.md`). Lowers to `NullAllele`.
#[test]
fn test_trans_c_with_null_allele_round_trip() {
    let input = "[NM_000088.3:c.100A>G];[0]";
    assert_eq!(parse_to_string(input), input);
    assert_eq!(parse_phase(input), AllelePhase::Trans);
    let parsed = parse_hgvs(input).expect("parse failed");
    if let HgvsVariant::Allele(a) = parsed {
        assert_eq!(a.variants.len(), 2);
        assert!(matches!(a.variants[1], HgvsVariant::NullAllele));
    } else {
        panic!("expected Allele");
    }
}

/// `c.[a];[pos=]`: one allele has a variant, the other is wild-type at the
/// SAME position. Spec line 41: `LRG_199t1:c.[2376G>C];[2376=]`. The note
/// on line 43 calls this out as the canonical pattern across edit types
/// (substitution, deletion, duplication, insertion).
#[test]
fn test_trans_c_with_wildtype_second_allele_round_trip() {
    // Same-accession, so Display emits compact form.
    assert_eq!(
        parse_to_string("NM_004006.2:c.[2376G>C];[2376=]"),
        "NM_004006.2:c.[2376G>C];[2376=]",
    );
    assert_eq!(
        parse_phase("NM_004006.2:c.[2376G>C];[2376=]"),
        AllelePhase::Trans,
    );
    // Spec note line 43 — same pattern for deletion.
    assert_eq!(
        parse_to_string("NM_004006.2:c.[2376del];[2376=]"),
        "NM_004006.2:c.[2376del];[2376=]",
    );
}

/// `g.[a];[?]` and `g.[a];[0]`: the same special tokens for the genomic
/// path. PR #146 added `parse_genome_trans_allele_shorthand`; this pin
/// asserts both special-token branches in that helper.
#[test]
fn test_trans_g_with_special_tokens_round_trip() {
    let unknown = "[NC_000001.11:g.100G>A];[?]";
    assert_eq!(parse_to_string(unknown), unknown);
    assert_eq!(parse_phase(unknown), AllelePhase::Trans);

    let null = "[NC_000001.11:g.100G>A];[0]";
    assert_eq!(parse_to_string(null), null);
    assert_eq!(parse_phase(null), AllelePhase::Trans);
}

/// `n.[a];[?]`: special tokens reach the new n.-trans helper added above.
/// The compact form parses; Display emits the expanded form for the [?]
/// branch (the Allele Display logic detects mixed-variant-kinds and falls
/// back). The Display→parse round-trip closes the loop.
#[test]
fn test_trans_n_with_unknown_allele_round_trip() {
    let parsed = parse_hgvs("NR_X.1:n.[100A>G];[?]").expect("parse failed");
    let displayed = format!("{}", parsed);
    assert_eq!(parse_to_string(&displayed), displayed);
    if let HgvsVariant::Allele(a) = parsed {
        assert_eq!(a.variants.len(), 2);
        assert!(matches!(a.variants[1], HgvsVariant::UnknownAllele));
        assert_eq!(a.phase, AllelePhase::Trans);
    } else {
        panic!("expected Allele");
    }
}

// ---------------------------------------------------------------------------
// Homozygous trans: same variant on both alleles. Distinct from cis-with-
// duplicate-sub-variants because the bracket shape encodes phase.
// ---------------------------------------------------------------------------

/// Spec line 38: `c.[2376G>C];[2376G>C]` — homozygous compound heterozygote.
/// Pre-existing in the c. trans-shorthand path; pinned for regression
/// coverage so a future parser change can't silently collapse a homozygous
/// trans into a cis duplicate.
#[test]
fn test_trans_c_homozygous_round_trip() {
    assert_eq!(
        parse_to_string("NM_004006.2:c.[2376G>C];[2376G>C]"),
        "NM_004006.2:c.[2376G>C];[2376G>C]",
    );
    assert_eq!(
        parse_phase("NM_004006.2:c.[2376G>C];[2376G>C]"),
        AllelePhase::Trans,
    );
}

/// Genomic homozygous trans — the `g.` counterpart of the c. case above,
/// exercising the `parse_genome_trans_allele_shorthand` path added in
/// PR #146.
#[test]
fn test_trans_g_homozygous_round_trip() {
    assert_eq!(
        parse_to_string("NC_000001.11:g.[100G>A];[100G>A]"),
        "NC_000001.11:g.[100G>A];[100G>A]",
    );
    assert_eq!(
        parse_phase("NC_000001.11:g.[100G>A];[100G>A]"),
        AllelePhase::Trans,
    );
}

// ---------------------------------------------------------------------------
// Normalize idempotency: normalize(normalize(x)) == normalize(x).
// ---------------------------------------------------------------------------

/// Trans alleles must be idempotent under repeated `normalize`. PR #146
/// asserted only that one-pass normalize doesn't break the trans
/// merge-barrier; this guards against rules that accidentally re-shuffle
/// trans on a second pass (e.g. by re-parsing through a path that drops
/// the phase wrapper).
#[test]
fn test_trans_normalize_is_idempotent() {
    let normalizer = Normalizer::new(MockProvider::new());
    let parsed =
        parse_hgvs("[NC_000001.11:g.1000G>A];[NC_000001.11:g.2000T>C]").expect("parse failed");
    let pass1 = normalizer.normalize(&parsed).expect("normalize pass1");
    let pass2 = normalizer.normalize(&pass1).expect("normalize pass2");
    assert_eq!(
        format!("{}", pass1),
        format!("{}", pass2),
        "normalize must be idempotent on trans alleles",
    );
    if let HgvsVariant::Allele(a) = pass2 {
        assert_eq!(a.phase, AllelePhase::Trans);
    } else {
        panic!("expected Allele after second normalize");
    }
}
