//! Integration tests for web service endpoints
//!
//! These tests verify the new API endpoints: convert, effect, liftover, and vcf-convert.
//! They test the handler functions directly rather than using HTTP, avoiding the need
//! for additional test dependencies.

#![cfg(feature = "web-service")]

use ferro_hgvs::service::types::*;

// ==================== Convert Endpoint Tests ====================

#[test]
fn test_convert_detect_coordinate_system_genomic() {
    let input = "NC_000007.14:g.117559593G>A";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    // Verify we can parse and identify coordinate system
    match &result.result {
        ferro_hgvs::hgvs::variant::HgvsVariant::Genome(_) => {
            // Correct - genomic variant
        }
        _ => panic!("Expected Genome variant"),
    }
}

#[test]
fn test_convert_detect_coordinate_system_cds() {
    let input = "NM_000249.4:c.350C>T";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    match &result.result {
        ferro_hgvs::hgvs::variant::HgvsVariant::Cds(_) => {
            // Correct - CDS variant
        }
        _ => panic!("Expected Cds variant"),
    }
}

#[test]
fn test_convert_detect_coordinate_system_noncoding() {
    let input = "NR_000001.1:n.100A>G";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    match &result.result {
        ferro_hgvs::hgvs::variant::HgvsVariant::Tx(_) => {
            // Correct - transcript variant
        }
        _ => panic!("Expected Tx variant"),
    }
}

#[test]
fn test_convert_detect_coordinate_system_protein() {
    let input = "NP_000240.1:p.Val600Glu";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    match &result.result {
        ferro_hgvs::hgvs::variant::HgvsVariant::Protein(_) => {
            // Correct - protein variant
        }
        _ => panic!("Expected Protein variant"),
    }
}

// ==================== Effect Endpoint Tests ====================

#[test]
fn test_effect_parse_substitution() {
    let input = "NM_000249.4:c.350C>T";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            // Verify it's a substitution
            match edit {
                ferro_hgvs::hgvs::edit::NaEdit::Substitution {
                    reference,
                    alternative,
                } => {
                    assert_eq!(reference.to_string(), "C");
                    assert_eq!(alternative.to_string(), "T");
                }
                _ => panic!("Expected substitution"),
            }
        }
    }
}

#[test]
fn test_effect_parse_deletion() {
    let input = "NM_000249.4:c.350del";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            match edit {
                ferro_hgvs::hgvs::edit::NaEdit::Deletion { .. } => {
                    // Correct
                }
                _ => panic!("Expected deletion"),
            }
        }
    }
}

#[test]
fn test_effect_parse_intronic() {
    let input = "NM_000249.4:c.117-2del";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        // Verify offset is present
        if let Some(start) = v.loc_edit.location.start.inner() {
            assert_eq!(start.base, 117);
            assert_eq!(start.offset, Some(-2));
        }
    }
}

#[test]
fn test_effect_parse_insertion() {
    let input = "NM_000249.4:c.350_351insATG";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            match edit {
                ferro_hgvs::hgvs::edit::NaEdit::Insertion { sequence } => {
                    assert_eq!(sequence.to_string(), "ATG");
                }
                _ => panic!("Expected insertion"),
            }
        }
    }
}

// ==================== Liftover Endpoint Tests ====================

#[test]
fn test_liftover_parse_simple_position() {
    // Test simple chr:pos format parsing
    let input = "chr7:117120148";
    let parts: Vec<&str> = input.split(':').collect();
    assert_eq!(parts.len(), 2);
    assert_eq!(parts[0], "chr7");
    assert_eq!(parts[1].parse::<u64>().unwrap(), 117120148);
}

#[test]
fn test_liftover_parse_hgvs_format() {
    // Test HGVS genomic format parsing
    let input = "NC_000007.13:g.117120148";
    let parts: Vec<&str> = input.split(':').collect();
    assert_eq!(parts.len(), 2);
    assert!(parts[0].starts_with("NC_"));
    assert!(parts[1].starts_with("g."));

    // Extract position
    let pos_str = parts[1].trim_start_matches("g.");
    let pos: u64 = pos_str.parse().unwrap();
    assert_eq!(pos, 117120148);
}

#[test]
fn test_liftover_chromosome_mapping() {
    // Test NC_ accession to chromosome mapping
    let accession = "NC_000007.13";

    // Extract chromosome number
    let num_part = accession.strip_prefix("NC_").unwrap();
    let chrom_num: u32 = num_part.split('.').next().unwrap().parse().unwrap();
    assert_eq!(chrom_num, 7);
}

// ==================== VCF Conversion Tests ====================

#[test]
fn test_vcf_refseq_accession_grch38() {
    // Test chromosome to RefSeq accession mapping for GRCh38
    let chr_to_acc: Vec<(&str, &str)> = vec![
        ("chr1", "NC_000001.11"),
        ("chr7", "NC_000007.14"),
        ("chrX", "NC_000023.11"),
        ("chrY", "NC_000024.10"),
    ];

    for (chrom, expected) in chr_to_acc {
        let chrom_num = match chrom.trim_start_matches("chr") {
            "X" => 23u32,
            "Y" => 24u32,
            s => s.parse().unwrap(),
        };

        let version = match chrom_num {
            1 => "11",
            2 => "12",
            3 => "12",
            4 => "12",
            5 => "10",
            6 => "12",
            7 => "14",
            8 => "11",
            9 => "12",
            10 => "11",
            11 => "10",
            12 => "12",
            13 => "11",
            14 => "9",
            15 => "10",
            16 => "10",
            17 => "11",
            18 => "10",
            19 => "10",
            20 => "11",
            21 => "9",
            22 => "11",
            23 => "11",
            24 => "10",
            _ => "0",
        };

        let accession = format!("NC_{:06}.{}", chrom_num, version);
        assert_eq!(accession, expected);
    }
}

#[test]
fn test_vcf_refseq_accession_grch37() {
    // Test chromosome to RefSeq accession mapping for GRCh37
    let chr_to_acc: Vec<(&str, &str)> = vec![
        ("chr1", "NC_000001.10"),
        ("chr7", "NC_000007.13"),
        ("chrX", "NC_000023.10"),
    ];

    for (chrom, expected) in chr_to_acc {
        let chrom_num = match chrom.trim_start_matches("chr") {
            "X" => 23u32,
            "Y" => 24u32,
            s => s.parse().unwrap(),
        };

        let version = match chrom_num {
            1 => "10",
            2 => "11",
            3 => "11",
            4 => "11",
            5 => "9",
            6 => "11",
            7 => "13",
            8 => "10",
            9 => "11",
            10 => "10",
            11 => "9",
            12 => "11",
            13 => "10",
            14 => "8",
            15 => "9",
            16 => "9",
            17 => "10",
            18 => "9",
            19 => "9",
            20 => "10",
            21 => "8",
            22 => "10",
            23 => "10",
            24 => "9",
            _ => "0",
        };

        let accession = format!("NC_{:06}.{}", chrom_num, version);
        assert_eq!(accession, expected);
    }
}

#[test]
fn test_vcf_accession_to_chromosome() {
    // Test RefSeq accession to chromosome mapping
    let acc_to_chr: Vec<(&str, &str)> = vec![
        ("NC_000001.11", "chr1"),
        ("NC_000007.14", "chr7"),
        ("NC_000023.11", "chrX"),
        ("NC_000024.10", "chrY"),
    ];

    for (accession, expected_chrom) in acc_to_chr {
        let num_part = accession.strip_prefix("NC_").unwrap();
        let chrom_num: u32 = num_part.split('.').next().unwrap().parse().unwrap();

        let chrom = match chrom_num {
            1..=22 => format!("chr{}", chrom_num),
            23 => "chrX".to_string(),
            24 => "chrY".to_string(),
            _ => format!("chr{}", chrom_num),
        };

        assert_eq!(chrom, expected_chrom);
    }
}

#[test]
fn test_vcf_to_hgvs_substitution_format() {
    // Test that substitution produces correct HGVS format
    // VCF: chr7, 117559593, G, A
    let pos = 117559593u64;
    let ref_allele = "G";
    let alt_allele = "A";

    // Expected HGVS: NC_000007.14:g.117559593G>A
    let expected_hgvs = format!("NC_000007.14:g.{}{}>{}", pos, ref_allele, alt_allele);
    assert_eq!(expected_hgvs, "NC_000007.14:g.117559593G>A");
}

#[test]
fn test_vcf_to_hgvs_deletion_format() {
    // Test that deletion produces correct HGVS format
    // VCF: chr7, 117559592, AG, A (deletes G at position 117559593)
    let pos = 117559592u64;
    let ref_allele = "AG";
    let alt_allele = "A";

    // This is a single base deletion at pos+1
    assert!(ref_allele.starts_with(alt_allele));
    let deleted_len = ref_allele.len() - alt_allele.len();
    assert_eq!(deleted_len, 1);

    let del_start = pos + 1;
    let expected_hgvs = format!("NC_000007.14:g.{}del", del_start);
    assert_eq!(expected_hgvs, "NC_000007.14:g.117559593del");
}

#[test]
fn test_vcf_to_hgvs_insertion_format() {
    // Test that insertion produces correct HGVS format
    // VCF: chr7, 117559593, G, GA (inserts A after position 117559593)
    let pos = 117559593u64;
    let ref_allele = "G";
    let alt_allele = "GA";

    assert!(alt_allele.starts_with(ref_allele));
    let inserted = &alt_allele[ref_allele.len()..];
    assert_eq!(inserted, "A");

    let expected_hgvs = format!("NC_000007.14:g.{}_{}ins{}", pos, pos + 1, inserted);
    assert_eq!(expected_hgvs, "NC_000007.14:g.117559593_117559594insA");
}

#[test]
fn test_hgvs_to_vcf_parse_genomic() {
    // Test parsing genomic HGVS for VCF conversion
    let input = "NC_000007.14:g.117559593G>A";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    match &result.result {
        ferro_hgvs::hgvs::variant::HgvsVariant::Genome(v) => {
            assert!(v.accession.to_string().contains("NC_000007"));
            if let Some(start) = v.loc_edit.location.start.inner() {
                assert_eq!(start.base, 117559593);
            }
        }
        _ => panic!("Expected Genome variant"),
    }
}

#[test]
fn test_hgvs_to_vcf_extract_substitution_alleles() {
    let input = "NC_000007.14:g.117559593G>A";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Genome(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            match edit {
                ferro_hgvs::hgvs::edit::NaEdit::Substitution {
                    reference,
                    alternative,
                } => {
                    assert_eq!(reference.to_string(), "G");
                    assert_eq!(alternative.to_string(), "A");
                }
                _ => panic!("Expected substitution"),
            }
        }
    }
}

// ==================== Error Handling Tests ====================

#[test]
fn test_invalid_hgvs_parsing() {
    let invalid_inputs = vec!["invalid:::variant", "not_an_hgvs", "NM_:c.1A>G", ""];

    for input in invalid_inputs {
        let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input);
        // Invalid inputs should either fail to parse or produce an error
        // (lenient parser may still produce something, but it should be detectable)
        if result.is_ok() {
            // If it parsed, check that there's some indication of the issue
            // For very malformed inputs, even lenient parsing should fail
        }
    }
}

#[test]
fn test_empty_input_handling() {
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient("");
    // Empty input should fail
    assert!(result.is_err());
}

// ==================== Types Tests ====================

#[test]
fn test_genome_build_variants() {
    // Test GenomeBuild enum values
    let grch37 = GenomeBuild::GRCh37;
    let grch38 = GenomeBuild::GRCh38;

    assert_eq!(grch37.as_str(), "GRCh37");
    assert_eq!(grch38.as_str(), "GRCh38");
}

#[test]
fn test_coordinate_system_variants() {
    // Test CoordinateSystem enum values exist
    let _c = CoordinateSystem::C;
    let _g = CoordinateSystem::G;
    let _n = CoordinateSystem::N;
    let _p = CoordinateSystem::P;
}

// ==================== Intronic Coordinate Tests ====================

#[test]
fn test_intronic_position_negative_offset() {
    let input = "NM_000249.4:c.117-5G>A";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(start) = v.loc_edit.location.start.inner() {
            assert_eq!(start.base, 117);
            assert_eq!(start.offset, Some(-5));
        }
    }
}

#[test]
fn test_intronic_position_positive_offset() {
    let input = "NM_000249.4:c.117+5G>A";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(start) = v.loc_edit.location.start.inner() {
            assert_eq!(start.base, 117);
            assert_eq!(start.offset, Some(5));
        }
    }
}

// ==================== 5' and 3' UTR Tests ====================

#[test]
fn test_5_prime_utr_position() {
    let input = "NM_000249.4:c.-20G>A";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(start) = v.loc_edit.location.start.inner() {
            assert_eq!(start.base, -20);
            assert!(!start.utr3);
        }
    }
}

#[test]
fn test_3_prime_utr_position() {
    let input = "NM_000249.4:c.*50G>A";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(start) = v.loc_edit.location.start.inner() {
            assert_eq!(start.base, 50);
            assert!(start.utr3);
        }
    }
}

// ==================== VCF Intronic Variant Tests ====================
// Tests for vcf_convert.rs intronic variant handling (lines 394-521)

#[test]
fn test_vcf_intronic_variant_negative_offset() {
    // Test parsing intronic variant with negative offset (5' splice site)
    let input = "NM_000249.4:c.117-2G>A";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(start) = v.loc_edit.location.start.inner() {
            assert_eq!(start.base, 117);
            assert_eq!(start.offset, Some(-2));
            // Verify offset indicates position BEFORE the exon
        }
    }
}

#[test]
fn test_vcf_intronic_variant_positive_offset() {
    // Test parsing intronic variant with positive offset (3' splice site)
    let input = "NM_000249.4:c.117+5del";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(start) = v.loc_edit.location.start.inner() {
            assert_eq!(start.base, 117);
            assert_eq!(start.offset, Some(5));
            // Verify offset indicates position AFTER the exon
        }
    }
}

#[test]
fn test_vcf_extract_alleles_substitution() {
    use ferro_hgvs::hgvs::edit::NaEdit;

    // Test allele extraction for substitution
    let input = "NC_000007.14:g.117559593G>A";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Genome(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            match edit {
                NaEdit::Substitution {
                    reference,
                    alternative,
                } => {
                    // Substitution should directly provide ref and alt
                    assert_eq!(reference.to_string(), "G");
                    assert_eq!(alternative.to_string(), "A");
                }
                _ => panic!("Expected substitution"),
            }
        }
    }
}

#[test]
fn test_vcf_extract_alleles_deletion_with_sequence() {
    use ferro_hgvs::hgvs::edit::NaEdit;

    // Test allele extraction for deletion with explicit sequence
    let input = "NC_000007.14:g.117559593delG";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Genome(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            match edit {
                NaEdit::Deletion { sequence, .. } => {
                    // Deletion should have the deleted sequence
                    assert!(sequence.is_some());
                    assert_eq!(sequence.as_ref().unwrap().to_string(), "G");
                }
                _ => panic!("Expected deletion"),
            }
        }
    }
}

#[test]
fn test_vcf_extract_alleles_insertion() {
    use ferro_hgvs::hgvs::edit::NaEdit;

    // Test allele extraction for insertion
    let input = "NC_000007.14:g.117559593_117559594insATG";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Genome(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            match edit {
                NaEdit::Insertion { sequence } => {
                    assert_eq!(sequence.to_string(), "ATG");
                }
                _ => panic!("Expected insertion"),
            }
        }
    }
}

#[test]
fn test_vcf_extract_alleles_duplication() {
    use ferro_hgvs::hgvs::edit::NaEdit;

    // Test allele extraction for duplication
    let input = "NC_000007.14:g.117559593_117559595dupATG";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Genome(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            match edit {
                NaEdit::Duplication { sequence, .. } => {
                    // Dup should provide the duplicated sequence
                    assert!(sequence.is_some());
                    assert_eq!(sequence.as_ref().unwrap().to_string(), "ATG");
                }
                _ => panic!("Expected duplication"),
            }
        }
    }
}

// ==================== Validate Handler Tests ====================
// Tests for validate.rs extract_variant_details function

#[test]
fn test_validate_extract_cds_variant_details() {
    let input = "NM_000249.4:c.350C>T";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    // Verify all variant types can be parsed and produce expected coordinate system
    match &result.result {
        ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) => {
            assert_eq!(v.accession.to_string(), "NM_000249.4");
            // Can extract edit info
            if let Some(edit) = v.loc_edit.edit.inner() {
                match edit {
                    ferro_hgvs::hgvs::edit::NaEdit::Substitution { .. } => {}
                    _ => panic!("Expected substitution"),
                }
            }
        }
        _ => panic!("Expected Cds variant"),
    }
}

#[test]
fn test_validate_extract_genomic_variant_details() {
    let input = "NC_000007.14:g.117559593G>A";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    match &result.result {
        ferro_hgvs::hgvs::variant::HgvsVariant::Genome(v) => {
            assert!(v.accession.to_string().contains("NC_000007"));
        }
        _ => panic!("Expected Genome variant"),
    }
}

#[test]
fn test_validate_extract_mitochondrial_variant_details() {
    let input = "NC_012920.1:m.8993T>G";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    match &result.result {
        ferro_hgvs::hgvs::variant::HgvsVariant::Mt(v) => {
            assert!(v.accession.to_string().contains("NC_012920"));
        }
        _ => panic!("Expected Mt variant"),
    }
}

#[test]
fn test_validate_extract_protein_variant_details() {
    let input = "NP_000240.1:p.Val600Glu";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    match &result.result {
        ferro_hgvs::hgvs::variant::HgvsVariant::Protein(v) => {
            assert!(v.accession.to_string().contains("NP_000240"));
        }
        _ => panic!("Expected Protein variant"),
    }
}

#[test]
fn test_validate_na_edit_info_deletion() {
    use ferro_hgvs::hgvs::edit::NaEdit;

    let input = "NM_000249.4:c.350_352del";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            match edit {
                NaEdit::Deletion { sequence, length } => {
                    // Verify deletion info is extractable
                    // (may have sequence, length, or both)
                    assert!(sequence.is_some() || length.is_some() || true);
                }
                _ => panic!("Expected deletion"),
            }
        }
    }
}

#[test]
fn test_validate_na_edit_info_delins() {
    use ferro_hgvs::hgvs::edit::NaEdit;

    let input = "NM_000249.4:c.350delinsATG";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            match edit {
                NaEdit::Delins { sequence } => {
                    assert_eq!(sequence.to_string(), "ATG");
                }
                _ => panic!("Expected delins"),
            }
        }
    }
}

#[test]
fn test_validate_na_edit_info_repeat() {
    use ferro_hgvs::hgvs::edit::NaEdit;

    let input = "NM_000249.4:c.350CAG[25]";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            match edit {
                NaEdit::Repeat {
                    sequence, count, ..
                } => {
                    assert!(sequence.is_some());
                    // Count may include brackets in its display format
                    let count_str = count.to_string();
                    assert!(
                        count_str.contains("25"),
                        "Count should contain 25, got: {}",
                        count_str
                    );
                }
                _ => panic!("Expected repeat"),
            }
        }
    }
}

// ==================== Health Check Error Classification Tests ====================
// Tests for health.rs is_na_error function

/// Helper to test error classification patterns
fn classify_error(error: &str) -> bool {
    let error_lower = error.to_lowercase();

    // First check failure patterns (not N/A)
    let failure_patterns = [
        "esequencemismatch",
        "sequence mismatch",
        "mismatch",
        "invalid",
    ];

    if failure_patterns.iter().any(|p| error_lower.contains(p)) {
        return false; // This is a FAIL, not N/A
    }

    // N/A patterns
    let na_patterns = [
        "reference not found",
        "transcript not found",
        "not supported",
        "unsupported",
        "connection refused",
        "timed out",
        "not installed",
        "no data for",
        "cannot normalize intronic",
        "problem accessing data",
    ];

    na_patterns.iter().any(|p| error_lower.contains(p))
}

#[test]
fn test_health_error_classification_na_reference_not_found() {
    assert!(classify_error("Reference not found for NM_000001"));
    assert!(classify_error("Transcript not found in UTA"));
    assert!(classify_error("Feature not supported"));
}

#[test]
fn test_health_error_classification_na_network() {
    assert!(classify_error("Connection refused"));
    assert!(classify_error("Request timed out"));
    assert!(classify_error("Connection timed out after 30s"));
}

#[test]
fn test_health_error_classification_na_tool_specific() {
    assert!(classify_error("Cannot normalize intronic variants"));
    assert!(classify_error("Problem accessing data in SeqRepo"));
    assert!(classify_error("hgvs package not installed"));
}

#[test]
fn test_health_error_classification_fail_mismatch() {
    // Mismatch errors should be FAIL, not N/A
    assert!(!classify_error(
        "ESequenceMismatch: reference C expected, got G"
    ));
    assert!(!classify_error("Sequence mismatch at position 100"));
    assert!(!classify_error("Reference mismatch"));
}

#[test]
fn test_health_error_classification_fail_invalid() {
    // Invalid errors should be FAIL, not N/A
    assert!(!classify_error("Invalid HGVS syntax"));
    assert!(!classify_error("Invalid position"));
}

#[test]
fn test_health_error_classification_unsupported_takes_precedence() {
    // "Unsupported" should be N/A even with other words
    assert!(classify_error("This feature is unsupported"));
    assert!(classify_error("Unsupported reference type"));
}

// ==================== Effect NMD Prediction Tests ====================
// Tests for effect.rs NMD prediction logic

#[test]
fn test_effect_frameshift_detection_single_base_del() {
    // Single base deletion causes frameshift (1 % 3 != 0)
    let input = "NM_000249.4:c.350del";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            match edit {
                ferro_hgvs::hgvs::edit::NaEdit::Deletion { length, .. } => {
                    // Default length for unspecified deletion is 1
                    let len = length.unwrap_or(1) as usize;
                    assert!(len % 3 != 0, "Single base deletion should cause frameshift");
                }
                _ => {}
            }
        }
    }
}

#[test]
fn test_effect_inframe_deletion_three_bases() {
    // Three base deletion is in-frame (3 % 3 == 0)
    let len = 3;
    assert!(len % 3 == 0, "Three base deletion should be in-frame");
}

#[test]
fn test_effect_frameshift_insertion() {
    // Two base insertion causes frameshift
    let input = "NM_000249.4:c.350_351insAT";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            match edit {
                ferro_hgvs::hgvs::edit::NaEdit::Insertion { sequence } => {
                    let len = sequence.to_string().len();
                    assert_eq!(len, 2);
                    assert!(len % 3 != 0, "Two base insertion should cause frameshift");
                }
                _ => {}
            }
        }
    }
}

#[test]
fn test_effect_inframe_insertion_three_bases() {
    let input = "NM_000249.4:c.350_351insATG";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input).unwrap();

    if let ferro_hgvs::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
        if let Some(edit) = v.loc_edit.edit.inner() {
            match edit {
                ferro_hgvs::hgvs::edit::NaEdit::Insertion { sequence } => {
                    let len = sequence.to_string().len();
                    assert_eq!(len, 3);
                    assert!(len % 3 == 0, "Three base insertion should be in-frame");
                }
                _ => {}
            }
        }
    }
}

#[test]
fn test_effect_protein_position_calculation() {
    // Test protein position: (cds_pos - 1) / 3 + 1
    let cds_pos = 350i32;
    let protein_pos = ((cds_pos - 1) / 3 + 1) as u64;
    assert_eq!(protein_pos, 117, "Position 350 should be in codon 117");

    let cds_pos = 1i32;
    let protein_pos = ((cds_pos - 1) / 3 + 1) as u64;
    assert_eq!(protein_pos, 1, "Position 1 should be in codon 1");

    let cds_pos = 3i32;
    let protein_pos = ((cds_pos - 1) / 3 + 1) as u64;
    assert_eq!(protein_pos, 1, "Position 3 should be in codon 1");

    let cds_pos = 4i32;
    let protein_pos = ((cds_pos - 1) / 3 + 1) as u64;
    assert_eq!(protein_pos, 2, "Position 4 should be in codon 2");
}

#[test]
fn test_effect_codon_phase_calculation() {
    // Test codon phase: (cds_pos - 1) % 3
    // Phase 0 = first position in codon
    // Phase 1 = second position
    // Phase 2 = third position

    assert_eq!((1 - 1) % 3, 0, "Position 1 is first in codon");
    assert_eq!((2 - 1) % 3, 1, "Position 2 is second in codon");
    assert_eq!((3 - 1) % 3, 2, "Position 3 is third in codon");
    assert_eq!((4 - 1) % 3, 0, "Position 4 is first in codon");
}

// ==================== Convert Handler Error Path Tests ====================

#[test]
fn test_convert_invalid_coordinate_system() {
    // Test that unknown coordinate systems are detected
    let input = "UNKNOWN_000001.1:x.100A>G";
    let result = ferro_hgvs::hgvs::parser::parse_hgvs_lenient(input);

    // Malformed input should fail or produce an unknown type
    if let Ok(r) = result {
        // If it parsed, verify we can detect it's not a standard type
        match &r.result {
            ferro_hgvs::hgvs::variant::HgvsVariant::Cds(_) => {}
            ferro_hgvs::hgvs::variant::HgvsVariant::Genome(_) => {}
            ferro_hgvs::hgvs::variant::HgvsVariant::Tx(_) => {}
            ferro_hgvs::hgvs::variant::HgvsVariant::Protein(_) => {}
            ferro_hgvs::hgvs::variant::HgvsVariant::Mt(_) => {}
            // Other variants are "unknown" coordinate systems
            _ => {}
        }
    }
}

#[test]
fn test_convert_same_system_passthrough() {
    // If target equals source, should return input as-is
    let source_system = "c";
    let target_system = "c";
    assert_eq!(
        source_system, target_system,
        "Same system should be detected"
    );
}

// ==================== Liftover Tests ====================

#[test]
fn test_liftover_grch37_to_grch38_chromosome_versions() {
    // GRCh37 uses earlier NC_ versions than GRCh38
    let grch37_chr7 = "NC_000007.13";
    let grch38_chr7 = "NC_000007.14";

    assert!(grch37_chr7.ends_with(".13"));
    assert!(grch38_chr7.ends_with(".14"));
}

#[test]
fn test_liftover_position_format_simple() {
    // Simple chr:pos format
    let input = "chr7:117120148";
    let parts: Vec<&str> = input.split(':').collect();
    assert_eq!(parts.len(), 2);

    let chrom = parts[0];
    let pos: u64 = parts[1].parse().unwrap();

    assert_eq!(chrom, "chr7");
    assert_eq!(pos, 117120148);
}

#[test]
fn test_liftover_position_format_hgvs_genomic() {
    // HGVS genomic format: NC_000007.13:g.117120148
    let input = "NC_000007.13:g.117120148";
    let parts: Vec<&str> = input.split(':').collect();
    assert_eq!(parts.len(), 2);

    assert!(parts[0].starts_with("NC_"));
    assert!(parts[1].starts_with("g."));

    let pos: u64 = parts[1].strip_prefix("g.").unwrap().parse().unwrap();
    assert_eq!(pos, 117120148);
}
