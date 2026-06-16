//! Property-based tests for HGVS variant parsing and normalization
//!
//! This module provides comprehensive property-based testing for the HGVS parser
//! using proptest. Tests cover all variant types and edit operations.

use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, Normalizer};
use proptest::prelude::*;
use proptest::test_runner::Config as ProptestConfig;

// =============================================================================
// Base strategies
// =============================================================================

/// Generate valid RefSeq genomic accession prefixes
#[allow(dead_code)]
fn genomic_accession_prefix() -> impl Strategy<Value = &'static str> {
    prop_oneof![Just("NC"), Just("NG"),]
}

/// Generate valid RefSeq transcript accession prefixes
#[allow(dead_code)]
fn transcript_accession_prefix() -> impl Strategy<Value = &'static str> {
    prop_oneof![Just("NM"), Just("NR"),]
}

/// Generate valid RefSeq protein accession prefixes
#[allow(dead_code)]
fn protein_accession_prefix() -> impl Strategy<Value = &'static str> {
    Just("NP")
}

/// Generate valid Ensembl accession
fn ensembl_transcript_accession() -> impl Strategy<Value = String> {
    "[0-9]{11}".prop_map(|n| format!("ENST{}", n))
}

/// Generate valid accession numbers (6 digits)
fn accession_number() -> impl Strategy<Value = String> {
    "[0-9]{6}".prop_map(|s| s)
}

/// Generate valid version numbers
fn version() -> impl Strategy<Value = u32> {
    1..20u32
}

/// Generate valid nucleotide bases
fn nucleotide() -> impl Strategy<Value = char> {
    prop_oneof![Just('A'), Just('C'), Just('G'), Just('T'),]
}

/// Generate valid RNA nucleotide bases
fn rna_nucleotide() -> impl Strategy<Value = char> {
    prop_oneof![Just('a'), Just('c'), Just('g'), Just('u'),]
}

/// Generate small positive position numbers
fn position() -> impl Strategy<Value = u64> {
    1..10000u64
}

/// Generate nucleotide sequence (1-10 bases)
#[allow(dead_code)]
fn nucleotide_sequence() -> impl Strategy<Value = String> {
    prop::collection::vec(nucleotide(), 1..=10).prop_map(|v| v.into_iter().collect())
}

/// Generate short nucleotide sequence (1-5 bases)
fn short_sequence() -> impl Strategy<Value = String> {
    prop::collection::vec(nucleotide(), 1..=5).prop_map(|v| v.into_iter().collect())
}

/// Generate 3-letter amino acid code
fn amino_acid() -> impl Strategy<Value = &'static str> {
    prop_oneof![
        Just("Ala"),
        Just("Arg"),
        Just("Asn"),
        Just("Asp"),
        Just("Cys"),
        Just("Gln"),
        Just("Glu"),
        Just("Gly"),
        Just("His"),
        Just("Ile"),
        Just("Leu"),
        Just("Lys"),
        Just("Met"),
        Just("Phe"),
        Just("Pro"),
        Just("Ser"),
        Just("Thr"),
        Just("Trp"),
        Just("Tyr"),
        Just("Val"),
        Just("Ter"),
    ]
}

/// Generate single-letter amino acid code
fn amino_acid_single() -> impl Strategy<Value = char> {
    prop_oneof![
        Just('A'),
        Just('R'),
        Just('N'),
        Just('D'),
        Just('C'),
        Just('Q'),
        Just('E'),
        Just('G'),
        Just('H'),
        Just('I'),
        Just('L'),
        Just('K'),
        Just('M'),
        Just('F'),
        Just('P'),
        Just('S'),
        Just('T'),
        Just('W'),
        Just('Y'),
        Just('V'),
        Just('*'), // Ter
    ]
}

/// Generate amino acid sequence (1-5 residues)
fn amino_acid_sequence() -> impl Strategy<Value = String> {
    prop::collection::vec(amino_acid(), 1..=5).prop_map(|v| v.join(""))
}

// =============================================================================
// Genomic variant strategies
// =============================================================================

/// Generate a genomic substitution variant string
fn genomic_substitution() -> impl Strategy<Value = String> {
    (
        accession_number(),
        version(),
        position(),
        nucleotide(),
        nucleotide(),
    )
        .prop_filter("ref != alt", |(_, _, _, r, a)| r != a)
        .prop_map(|(num, ver, pos, ref_base, alt_base)| {
            format!("NC_{}.{}:g.{}{}>{}", num, ver, pos, ref_base, alt_base)
        })
}

/// Generate a genomic deletion variant string
fn genomic_deletion() -> impl Strategy<Value = String> {
    (accession_number(), version(), position())
        .prop_map(|(num, ver, pos)| format!("NC_{}.{}:g.{}del", num, ver, pos))
}

/// Generate a range deletion variant string
fn genomic_range_deletion() -> impl Strategy<Value = String> {
    (accession_number(), version(), position(), 1..100u64).prop_map(|(num, ver, start, len)| {
        format!("NC_{}.{}:g.{}_{}del", num, ver, start, start + len)
    })
}

/// Generate a genomic duplication variant string
fn genomic_duplication() -> impl Strategy<Value = String> {
    (accession_number(), version(), position())
        .prop_map(|(num, ver, pos)| format!("NC_{}.{}:g.{}dup", num, ver, pos))
}

/// Generate a genomic range duplication variant string
fn genomic_range_duplication() -> impl Strategy<Value = String> {
    (accession_number(), version(), position(), 1..50u64).prop_map(|(num, ver, start, len)| {
        format!("NC_{}.{}:g.{}_{}dup", num, ver, start, start + len)
    })
}

/// Generate a genomic insertion variant string
fn genomic_insertion() -> impl Strategy<Value = String> {
    (accession_number(), version(), position(), short_sequence()).prop_map(
        |(num, ver, pos, seq)| format!("NC_{}.{}:g.{}_{}ins{}", num, ver, pos, pos + 1, seq),
    )
}

/// Generate a genomic delins variant string
fn genomic_delins() -> impl Strategy<Value = String> {
    (
        accession_number(),
        version(),
        position(),
        1..10u64,
        short_sequence(),
    )
        .prop_map(|(num, ver, start, len, seq)| {
            format!(
                "NC_{}.{}:g.{}_{}delins{}",
                num,
                ver,
                start,
                start + len,
                seq
            )
        })
}

/// Generate a genomic inversion variant string
fn genomic_inversion() -> impl Strategy<Value = String> {
    (accession_number(), version(), position(), 10..100u64).prop_map(|(num, ver, start, len)| {
        format!("NC_{}.{}:g.{}_{}inv", num, ver, start, start + len)
    })
}

/// Generate a genomic repeat variant string
fn genomic_repeat() -> impl Strategy<Value = String> {
    (accession_number(), version(), position(), 5..50u64)
        .prop_map(|(num, ver, pos, count)| format!("NC_{}.{}:g.{}[{}]", num, ver, pos, count))
}

/// Generate any genomic variant
fn any_genomic_variant() -> impl Strategy<Value = String> {
    prop_oneof![
        3 => genomic_substitution(),
        2 => genomic_deletion(),
        2 => genomic_range_deletion(),
        1 => genomic_duplication(),
        1 => genomic_range_duplication(),
        1 => genomic_insertion(),
        1 => genomic_delins(),
        1 => genomic_inversion(),
        1 => genomic_repeat(),
    ]
}

// =============================================================================
// CDS (coding) variant strategies
// =============================================================================

/// Generate a CDS substitution variant string (using known transcript)
fn cds_substitution() -> impl Strategy<Value = String> {
    (1..50i64, nucleotide(), nucleotide())
        .prop_filter("ref != alt", |(_, r, a)| r != a)
        .prop_map(|(pos, ref_base, alt_base)| {
            format!("NM_000088.3:c.{}{}>{}", pos, ref_base, alt_base)
        })
}

/// Generate a CDS deletion variant string
fn cds_deletion() -> impl Strategy<Value = String> {
    (accession_number(), version(), 1..1000i64)
        .prop_map(|(num, ver, pos)| format!("NM_{}.{}:c.{}del", num, ver, pos))
}

/// Generate a CDS range deletion variant string
fn cds_range_deletion() -> impl Strategy<Value = String> {
    (accession_number(), version(), 1..1000i64, 1..50i64).prop_map(|(num, ver, start, len)| {
        format!("NM_{}.{}:c.{}_{}del", num, ver, start, start + len)
    })
}

/// Generate a CDS duplication variant string
fn cds_duplication() -> impl Strategy<Value = String> {
    (accession_number(), version(), 1..1000i64)
        .prop_map(|(num, ver, pos)| format!("NM_{}.{}:c.{}dup", num, ver, pos))
}

/// Generate a CDS insertion variant string
fn cds_insertion() -> impl Strategy<Value = String> {
    (accession_number(), version(), 1..1000i64, short_sequence()).prop_map(
        |(num, ver, pos, seq)| format!("NM_{}.{}:c.{}_{}ins{}", num, ver, pos, pos + 1, seq),
    )
}

/// Generate a CDS delins variant string
fn cds_delins() -> impl Strategy<Value = String> {
    (
        accession_number(),
        version(),
        1..1000i64,
        1..10i64,
        short_sequence(),
    )
        .prop_map(|(num, ver, start, len, seq)| {
            format!(
                "NM_{}.{}:c.{}_{}delins{}",
                num,
                ver,
                start,
                start + len,
                seq
            )
        })
}

/// Generate a CDS variant with intronic offset
fn cds_intronic() -> impl Strategy<Value = String> {
    (
        accession_number(),
        version(),
        1..100i64,
        1..50i64,
        nucleotide(),
        nucleotide(),
    )
        .prop_filter("ref != alt", |(_, _, _, _, r, a)| r != a)
        .prop_map(|(num, ver, exon_pos, offset, ref_base, alt_base)| {
            format!(
                "NM_{}.{}:c.{}+{}{}>{}",
                num, ver, exon_pos, offset, ref_base, alt_base
            )
        })
}

/// Generate a CDS variant with negative intronic offset
fn cds_intronic_negative() -> impl Strategy<Value = String> {
    (
        accession_number(),
        version(),
        2..100i64,
        1..50i64,
        nucleotide(),
        nucleotide(),
    )
        .prop_filter("ref != alt", |(_, _, _, _, r, a)| r != a)
        .prop_map(|(num, ver, exon_pos, offset, ref_base, alt_base)| {
            format!(
                "NM_{}.{}:c.{}-{}{}>{}",
                num, ver, exon_pos, offset, ref_base, alt_base
            )
        })
}

/// Generate any CDS variant
fn any_cds_variant() -> impl Strategy<Value = String> {
    prop_oneof![
        3 => cds_substitution(),
        2 => cds_deletion(),
        2 => cds_range_deletion(),
        1 => cds_duplication(),
        1 => cds_insertion(),
        1 => cds_delins(),
        1 => cds_intronic(),
        1 => cds_intronic_negative(),
    ]
}

// =============================================================================
// Non-coding (n.) variant strategies
// =============================================================================

/// Generate a non-coding substitution variant string
fn noncoding_substitution() -> impl Strategy<Value = String> {
    (
        accession_number(),
        version(),
        position(),
        nucleotide(),
        nucleotide(),
    )
        .prop_filter("ref != alt", |(_, _, _, r, a)| r != a)
        .prop_map(|(num, ver, pos, ref_base, alt_base)| {
            format!("NR_{}.{}:n.{}{}>{}", num, ver, pos, ref_base, alt_base)
        })
}

/// Generate a non-coding deletion variant string
fn noncoding_deletion() -> impl Strategy<Value = String> {
    (accession_number(), version(), position())
        .prop_map(|(num, ver, pos)| format!("NR_{}.{}:n.{}del", num, ver, pos))
}

/// Generate any non-coding variant
fn any_noncoding_variant() -> impl Strategy<Value = String> {
    prop_oneof![
        2 => noncoding_substitution(),
        1 => noncoding_deletion(),
    ]
}

// =============================================================================
// RNA (r.) variant strategies
// =============================================================================

/// Generate an RNA substitution variant string
fn rna_substitution() -> impl Strategy<Value = String> {
    (
        accession_number(),
        version(),
        position(),
        rna_nucleotide(),
        rna_nucleotide(),
    )
        .prop_filter("ref != alt", |(_, _, _, r, a)| r != a)
        .prop_map(|(num, ver, pos, ref_base, alt_base)| {
            format!("NM_{}.{}:r.{}{}>{}", num, ver, pos, ref_base, alt_base)
        })
}

/// Generate an RNA deletion variant string
fn rna_deletion() -> impl Strategy<Value = String> {
    (accession_number(), version(), position())
        .prop_map(|(num, ver, pos)| format!("NM_{}.{}:r.{}del", num, ver, pos))
}

/// Generate any RNA variant
fn any_rna_variant() -> impl Strategy<Value = String> {
    prop_oneof![
        2 => rna_substitution(),
        1 => rna_deletion(),
    ]
}

// =============================================================================
// Protein (p.) variant strategies
// =============================================================================

/// Generate a protein substitution variant string
fn protein_substitution() -> impl Strategy<Value = String> {
    (
        accession_number(),
        version(),
        amino_acid(),
        1..1000u64,
        amino_acid(),
    )
        .prop_filter("ref != alt", |(_, _, r, _, a)| r != a)
        .prop_map(|(num, ver, ref_aa, pos, alt_aa)| {
            format!("NP_{}.{}:p.{}{}{}", num, ver, ref_aa, pos, alt_aa)
        })
}

/// Generate a protein substitution with single-letter codes
fn protein_substitution_single() -> impl Strategy<Value = String> {
    (
        accession_number(),
        version(),
        amino_acid_single(),
        1..1000u64,
        amino_acid_single(),
    )
        .prop_filter("ref != alt", |(_, _, r, _, a)| r != a)
        .prop_map(|(num, ver, ref_aa, pos, alt_aa)| {
            format!("NP_{}.{}:p.{}{}{}", num, ver, ref_aa, pos, alt_aa)
        })
}

/// Generate a protein deletion variant string
fn protein_deletion() -> impl Strategy<Value = String> {
    (accession_number(), version(), amino_acid(), 1..1000u64)
        .prop_map(|(num, ver, aa, pos)| format!("NP_{}.{}:p.{}{}del", num, ver, aa, pos))
}

/// Generate a protein range deletion variant string
fn protein_range_deletion() -> impl Strategy<Value = String> {
    (
        accession_number(),
        version(),
        amino_acid(),
        1..500u64,
        amino_acid(),
        1..100u64,
    )
        .prop_map(|(num, ver, aa1, pos1, aa2, len)| {
            format!(
                "NP_{}.{}:p.{}{}_{}{}del",
                num,
                ver,
                aa1,
                pos1,
                aa2,
                pos1 + len
            )
        })
}

/// Generate a protein duplication variant string
fn protein_duplication() -> impl Strategy<Value = String> {
    (accession_number(), version(), amino_acid(), 1..1000u64)
        .prop_map(|(num, ver, aa, pos)| format!("NP_{}.{}:p.{}{}dup", num, ver, aa, pos))
}

/// Generate a protein insertion variant string
fn protein_insertion() -> impl Strategy<Value = String> {
    (
        accession_number(),
        version(),
        amino_acid(),
        1..500u64,
        amino_acid(),
        amino_acid_sequence(),
    )
        .prop_map(|(num, ver, aa1, pos1, aa2, seq)| {
            format!(
                "NP_{}.{}:p.{}{}_{}{}ins{}",
                num,
                ver,
                aa1,
                pos1,
                aa2,
                pos1 + 1,
                seq
            )
        })
}

/// Generate a protein delins variant string
fn protein_delins() -> impl Strategy<Value = String> {
    (
        accession_number(),
        version(),
        amino_acid(),
        1..500u64,
        amino_acid(),
        1..10u64,
        amino_acid_sequence(),
    )
        .prop_map(|(num, ver, aa1, pos1, aa2, len, seq)| {
            format!(
                "NP_{}.{}:p.{}{}_{}{}delins{}",
                num,
                ver,
                aa1,
                pos1,
                aa2,
                pos1 + len,
                seq
            )
        })
}

/// Generate a protein frameshift variant string
fn protein_frameshift() -> impl Strategy<Value = String> {
    (
        accession_number(),
        version(),
        amino_acid(),
        1..1000u64,
        proptest::option::of(2..100u64),
    )
        .prop_map(|(num, ver, aa, pos, ter_pos)| match ter_pos {
            Some(t) => format!("NP_{}.{}:p.{}{}fsTer{}", num, ver, aa, pos, t),
            None => format!("NP_{}.{}:p.{}{}fs", num, ver, aa, pos),
        })
}

/// Generate a protein extension variant string
fn protein_extension() -> impl Strategy<Value = String> {
    (accession_number(), version(), 1..100i64)
        .prop_map(|(num, ver, count)| format!("NP_{}.{}:p.Met1ext-{}", num, ver, count))
}

/// Generate protein identity variant (p.=)
fn protein_identity() -> impl Strategy<Value = String> {
    (accession_number(), version()).prop_map(|(num, ver)| format!("NP_{}.{}:p.=", num, ver))
}

/// Generate predicted protein identity variant (p.(=))
fn protein_identity_predicted() -> impl Strategy<Value = String> {
    (accession_number(), version()).prop_map(|(num, ver)| format!("NP_{}.{}:p.(=)", num, ver))
}

/// Generate position-specific protein identity (p.Val600=)
fn protein_position_identity() -> impl Strategy<Value = String> {
    (accession_number(), version(), amino_acid(), 1..1000u64)
        .prop_map(|(num, ver, aa, pos)| format!("NP_{}.{}:p.{}{}=", num, ver, aa, pos))
}

/// Generate no protein variant (p.0)
fn protein_no_protein() -> impl Strategy<Value = String> {
    (accession_number(), version()).prop_map(|(num, ver)| format!("NP_{}.{}:p.0", num, ver))
}

/// Generate predicted no protein variant (p.0?)
fn protein_no_protein_predicted() -> impl Strategy<Value = String> {
    (accession_number(), version()).prop_map(|(num, ver)| format!("NP_{}.{}:p.0?", num, ver))
}

/// Generate any protein variant
fn any_protein_variant() -> impl Strategy<Value = String> {
    prop_oneof![
        3 => protein_substitution(),
        1 => protein_substitution_single(),
        2 => protein_deletion(),
        1 => protein_range_deletion(),
        1 => protein_duplication(),
        1 => protein_insertion(),
        1 => protein_delins(),
        1 => protein_frameshift(),
        1 => protein_extension(),
        1 => protein_identity(),
        1 => protein_identity_predicted(),
        1 => protein_position_identity(),
        1 => protein_no_protein(),
        1 => protein_no_protein_predicted(),
    ]
}

// =============================================================================
// Ensembl variant strategies
// =============================================================================

/// Generate Ensembl transcript variants
fn ensembl_cds_variant() -> impl Strategy<Value = String> {
    (
        ensembl_transcript_accession(),
        1..1000i64,
        nucleotide(),
        nucleotide(),
    )
        .prop_filter("ref != alt", |(_, _, r, a)| r != a)
        .prop_map(|(acc, pos, ref_base, alt_base)| {
            format!("{}:c.{}{}>{}", acc, pos, ref_base, alt_base)
        })
}

// =============================================================================
// All variant type composite strategy
// =============================================================================

/// Generate any parseable variant
fn any_variant() -> impl Strategy<Value = String> {
    prop_oneof![
        4 => any_genomic_variant(),
        4 => any_cds_variant(),
        2 => any_noncoding_variant(),
        2 => any_rna_variant(),
        4 => any_protein_variant(),
        1 => ensembl_cds_variant(),
    ]
}

// =============================================================================
// Property tests
// =============================================================================

proptest! {
    #![proptest_config(ProptestConfig::with_cases(500))]

    // -------------------------------------------------------------------------
    // Genomic variant roundtrip tests
    // -------------------------------------------------------------------------

    /// Test that genomic substitutions roundtrip through parse -> display -> parse
    #[test]
    fn test_genomic_substitution_roundtrip(variant in genomic_substitution()) {
        let parsed = parse_hgvs(&variant);
        prop_assert!(parsed.is_ok(), "Failed to parse: {}", variant);

        let displayed = format!("{}", parsed.unwrap());
        let reparsed = parse_hgvs(&displayed);
        prop_assert!(reparsed.is_ok(), "Failed to reparse: {}", displayed);
    }

    /// Test that genomic deletions roundtrip through parse -> display -> parse
    #[test]
    fn test_genomic_deletion_roundtrip(variant in genomic_deletion()) {
        let parsed = parse_hgvs(&variant);
        prop_assert!(parsed.is_ok(), "Failed to parse: {}", variant);

        let displayed = format!("{}", parsed.unwrap());
        let reparsed = parse_hgvs(&displayed);
        prop_assert!(reparsed.is_ok(), "Failed to reparse: {}", displayed);
    }

    /// Test that range deletions roundtrip through parse -> display -> parse
    #[test]
    fn test_range_deletion_roundtrip(variant in genomic_range_deletion()) {
        let parsed = parse_hgvs(&variant);
        prop_assert!(parsed.is_ok(), "Failed to parse: {}", variant);

        let displayed = format!("{}", parsed.unwrap());
        let reparsed = parse_hgvs(&displayed);
        prop_assert!(reparsed.is_ok(), "Failed to reparse: {}", displayed);
    }

    /// Test all genomic variants roundtrip
    #[test]
    fn test_any_genomic_roundtrip(variant in any_genomic_variant()) {
        let parsed = parse_hgvs(&variant);
        prop_assert!(parsed.is_ok(), "Failed to parse genomic variant: {}", variant);

        let displayed = format!("{}", parsed.unwrap());
        let reparsed = parse_hgvs(&displayed);
        prop_assert!(reparsed.is_ok(), "Failed to reparse genomic variant: {}", displayed);
    }

    // -------------------------------------------------------------------------
    // CDS variant roundtrip tests
    // -------------------------------------------------------------------------

    /// Test that CDS variants roundtrip
    #[test]
    fn test_cds_roundtrip(variant in any_cds_variant()) {
        let parsed = parse_hgvs(&variant);
        prop_assert!(parsed.is_ok(), "Failed to parse CDS variant: {}", variant);

        let displayed = format!("{}", parsed.unwrap());
        let reparsed = parse_hgvs(&displayed);
        prop_assert!(reparsed.is_ok(), "Failed to reparse CDS variant: {}", displayed);
    }

    // -------------------------------------------------------------------------
    // Non-coding variant roundtrip tests
    // -------------------------------------------------------------------------

    /// Test that non-coding variants roundtrip
    #[test]
    fn test_noncoding_roundtrip(variant in any_noncoding_variant()) {
        let parsed = parse_hgvs(&variant);
        prop_assert!(parsed.is_ok(), "Failed to parse non-coding variant: {}", variant);

        let displayed = format!("{}", parsed.unwrap());
        let reparsed = parse_hgvs(&displayed);
        prop_assert!(reparsed.is_ok(), "Failed to reparse non-coding variant: {}", displayed);
    }

    // -------------------------------------------------------------------------
    // RNA variant roundtrip tests
    // -------------------------------------------------------------------------

    /// Test that RNA variants roundtrip
    #[test]
    fn test_rna_roundtrip(variant in any_rna_variant()) {
        let parsed = parse_hgvs(&variant);
        prop_assert!(parsed.is_ok(), "Failed to parse RNA variant: {}", variant);

        let displayed = format!("{}", parsed.unwrap());
        let reparsed = parse_hgvs(&displayed);
        prop_assert!(reparsed.is_ok(), "Failed to reparse RNA variant: {}", displayed);
    }

    // -------------------------------------------------------------------------
    // Protein variant roundtrip tests
    // -------------------------------------------------------------------------

    /// Test that protein variants roundtrip
    #[test]
    fn test_protein_roundtrip(variant in any_protein_variant()) {
        let parsed = parse_hgvs(&variant);
        prop_assert!(parsed.is_ok(), "Failed to parse protein variant: {}", variant);

        let displayed = format!("{}", parsed.unwrap());
        let reparsed = parse_hgvs(&displayed);
        prop_assert!(reparsed.is_ok(), "Failed to reparse protein variant: {}", displayed);
    }

    // -------------------------------------------------------------------------
    // Ensembl variant roundtrip tests
    // -------------------------------------------------------------------------

    /// Test that Ensembl variants roundtrip
    #[test]
    fn test_ensembl_roundtrip(variant in ensembl_cds_variant()) {
        let parsed = parse_hgvs(&variant);
        prop_assert!(parsed.is_ok(), "Failed to parse Ensembl variant: {}", variant);

        let displayed = format!("{}", parsed.unwrap());
        let reparsed = parse_hgvs(&displayed);
        prop_assert!(reparsed.is_ok(), "Failed to reparse Ensembl variant: {}", displayed);
    }

    // -------------------------------------------------------------------------
    // Any variant roundtrip test
    // -------------------------------------------------------------------------

    /// Test that any variant can roundtrip
    #[test]
    fn test_any_variant_roundtrip(variant in any_variant()) {
        let parsed = parse_hgvs(&variant);
        prop_assert!(parsed.is_ok(), "Failed to parse variant: {}", variant);

        let displayed = format!("{}", parsed.unwrap());
        let reparsed = parse_hgvs(&displayed);
        prop_assert!(reparsed.is_ok(), "Failed to reparse variant: {}", displayed);
    }

    // -------------------------------------------------------------------------
    // Normalization tests
    // -------------------------------------------------------------------------

    /// Test that normalization is idempotent (normalize twice = normalize once)
    #[test]
    fn test_normalization_idempotent(variant in cds_substitution()) {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let parsed = parse_hgvs(&variant);
        prop_assume!(parsed.is_ok());
        let parsed = parsed.unwrap();

        // Normalize once
        let normalized = normalizer.normalize(&parsed);
        prop_assume!(normalized.is_ok());
        let normalized = normalized.unwrap();

        // Create new normalizer (to get fresh provider reference)
        let provider2 = MockProvider::with_test_data();
        let normalizer2 = Normalizer::new(provider2);

        // Normalize twice
        let normalized_twice = normalizer2.normalize(&normalized);
        prop_assert!(normalized_twice.is_ok());

        // Should produce same result
        let twice_str = format!("{}", normalized_twice.unwrap());
        let once_str = format!("{}", normalized);
        prop_assert_eq!(once_str, twice_str, "Normalization not idempotent");
    }

    // -------------------------------------------------------------------------
    // Invariant tests
    // -------------------------------------------------------------------------

    /// Test that parsing never panics on any input (resilience test)
    #[test]
    fn test_parser_does_not_panic(input in ".*") {
        // This should never panic, regardless of input
        let _ = parse_hgvs(&input);
    }

    /// Test that variant type is preserved through display
    #[test]
    fn test_variant_type_preserved(variant in any_variant()) {
        let parsed = parse_hgvs(&variant);
        prop_assume!(parsed.is_ok());
        let parsed = parsed.unwrap();

        let displayed = format!("{}", parsed);
        let reparsed = parse_hgvs(&displayed);
        prop_assume!(reparsed.is_ok());
        let reparsed = reparsed.unwrap();

        // Variant types should match
        let type_matches = matches!(
            (&parsed, &reparsed),
            (HgvsVariant::Genome(_), HgvsVariant::Genome(_))
            | (HgvsVariant::Cds(_), HgvsVariant::Cds(_))
            | (HgvsVariant::Tx(_), HgvsVariant::Tx(_))
            | (HgvsVariant::Rna(_), HgvsVariant::Rna(_))
            | (HgvsVariant::Protein(_), HgvsVariant::Protein(_))
            | (HgvsVariant::Mt(_), HgvsVariant::Mt(_))
            | (HgvsVariant::Allele(_), HgvsVariant::Allele(_))
        );
        prop_assert!(type_matches, "Variant type changed through display: {:?} vs {:?}", parsed, reparsed);
    }

    /// Test that display strings are non-empty
    #[test]
    fn test_display_nonempty(variant in any_variant()) {
        let parsed = parse_hgvs(&variant);
        prop_assume!(parsed.is_ok());
        let displayed = format!("{}", parsed.unwrap());
        prop_assert!(!displayed.is_empty(), "Display produced empty string");
    }

    /// Test that accession is preserved through roundtrip
    #[test]
    fn test_accession_preserved(variant in any_genomic_variant()) {
        let parsed = parse_hgvs(&variant);
        prop_assume!(parsed.is_ok());
        let parsed = parsed.unwrap();

        let displayed = format!("{}", parsed);
        let reparsed = parse_hgvs(&displayed);
        prop_assume!(reparsed.is_ok());
        let _reparsed = reparsed.unwrap(); // Keep for verification that reparse works

        // Extract accession from both (they should be in the displayed form)
        // For genomic variants, the accession is at the start
        let orig_acc: String = variant.chars().take_while(|c| *c != ':').collect();
        let repr_acc: String = displayed.chars().take_while(|c| *c != ':').collect();
        prop_assert_eq!(orig_acc, repr_acc, "Accession changed through roundtrip");
    }
}

/// Test that protein variants maintain their type through normalization
#[test]
fn test_protein_variant_type_preserved() {
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NP_000079.2:p.Val600Glu").unwrap();
    let normalized = normalizer.normalize(&variant).unwrap();

    assert!(matches!(normalized, HgvsVariant::Protein(_)));
}

/// Test that genomic variants maintain their type through normalization
#[test]
fn test_genomic_variant_type_preserved() {
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
    let normalized = normalizer.normalize(&variant).unwrap();

    assert!(matches!(normalized, HgvsVariant::Genome(_)));
}

/// Test that CDS variants maintain their type through normalization
#[test]
fn test_cds_variant_type_preserved() {
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
    let normalized = normalizer.normalize(&variant).unwrap();

    assert!(matches!(normalized, HgvsVariant::Cds(_)));
}
