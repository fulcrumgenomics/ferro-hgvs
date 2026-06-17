//! HGVS Standard Compliance Tests
//!
//! These tests verify compliance with the HGVS nomenclature standard
//! as defined at https://varnomen.hgvs.org/
//!
//! Tests are organized by HGVS variant type and edit type.

use ferro_hgvs::{parse_hgvs, HgvsVariant};

/// Test genomic (g.) variant parsing
mod genomic_variants {
    use super::*;

    #[test]
    fn test_genomic_substitution() {
        let variants = vec![
            "NC_000001.11:g.12345A>G",
            "NC_000001.11:g.1A>C",
            "NC_000023.10:g.33038255C>A",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
            assert!(matches!(result.unwrap(), HgvsVariant::Genome(_)));
        }
    }

    #[test]
    fn test_genomic_deletion() {
        let variants = vec![
            "NC_000001.11:g.12345del",
            "NC_000001.11:g.12345_12350del",
            "NG_012232.1:g.19del",
            "NG_012232.1:g.19_21del",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_genomic_duplication() {
        let variants = vec![
            "NC_000001.11:g.12345dup",
            "NC_000001.11:g.12345_12350dup",
            "NC_000023.10:g.33229407_33229410dup",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_genomic_insertion() {
        let variants = vec![
            "NC_000001.11:g.12345_12346insA",
            "NC_000001.11:g.12345_12346insATG",
            "NC_000023.10:g.32867861_32867862insT",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_genomic_delins() {
        let variants = vec![
            "NC_000001.11:g.12345delinsT",
            "NC_000001.11:g.12345_12350delinsATG",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_genomic_inversion() {
        let variants = vec![
            "NC_000001.11:g.12345_12350inv",
            "NC_000023.10:g.32361330_32361333inv",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }
}

/// Test coding DNA (c.) variant parsing
mod coding_variants {
    use super::*;

    #[test]
    fn test_cds_substitution() {
        let variants = vec![
            "NM_000088.3:c.10A>G",
            "NM_004006.1:c.145C>T",
            "NM_000088.3:c.589G>T",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
            assert!(matches!(result.unwrap(), HgvsVariant::Cds(_)));
        }
    }

    #[test]
    fn test_cds_intronic() {
        let variants = vec![
            "NM_000088.3:c.100+5G>A",
            "NM_000088.3:c.100-10del",
            "NM_004006.1:c.93+1G>T",
            "NM_000088.3:c.589-5G>T",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_cds_utr() {
        let variants = vec![
            "NM_000088.3:c.-20G>A", // 5' UTR
            "NM_000088.3:c.*50A>G", // 3' UTR
            "NM_004006.1:c.-20G>A",
            "NM_004006.1:c.*50A>G",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_cds_deletion() {
        let variants = vec![
            "NM_000088.3:c.100del",
            "NM_000088.3:c.100_102del",
            "NM_004006.1:c.720_991del",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_cds_duplication() {
        let variants = vec![
            "NM_000088.3:c.100dup",
            "NM_000088.3:c.100_102dup",
            "NM_004006.2:c.20dup",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_cds_insertion() {
        let variants = vec![
            "NM_000088.3:c.100_101insA",
            "NM_000088.3:c.100_101insATG",
            "NM_004006.2:c.169_170insA",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_cds_delins() {
        let variants = vec![
            "NM_000088.3:c.100delinsT",
            "NM_000088.3:c.100_102delinsATG",
            "NM_004006.1:c.79_80delinsTT",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_cds_inversion() {
        let variants = vec!["NM_000088.3:c.100_105inv", "NM_004006.2:c.5657_5660inv"];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }
}

/// Test transcript (n.) variant parsing
mod transcript_variants {
    use super::*;

    #[test]
    fn test_transcript_substitution() {
        let variants = vec!["NR_000001.1:n.100A>G", "NR_000001.1:n.1A>C"];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
            assert!(matches!(result.unwrap(), HgvsVariant::Tx(_)));
        }
    }

    #[test]
    fn test_transcript_deletion() {
        let variants = vec!["NR_000001.1:n.100del", "NR_000001.1:n.100_105del"];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }
}

/// Test protein (p.) variant parsing
mod protein_variants {
    use super::*;

    #[test]
    fn test_protein_substitution() {
        let variants = vec![
            "NP_000079.2:p.Val600Glu",
            "NP_000001.1:p.Arg100Gly",
            "NP_000079.2:p.Ala1Cys",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
            assert!(matches!(result.unwrap(), HgvsVariant::Protein(_)));
        }
    }

    #[test]
    fn test_protein_deletion() {
        let variants = vec!["NP_000079.2:p.Ala50del", "NP_000001.1:p.Arg100_Gly105del"];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_protein_frameshift() {
        let variants = vec!["NP_000079.2:p.Arg100fs", "NP_000079.2:p.Gly25fs"];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_protein_single_letter_codes() {
        // Single-letter amino acid codes (common in clinical reporting)
        let variants = vec![
            "NP_000079.2:p.V600E", // BRAF V600E
            "NP_000545.1:p.R248Q", // TP53 R248Q
            "NP_000138.2:p.G12D",  // KRAS G12D
            "NP_000138.2:p.G12V",  // KRAS G12V
            "NP_004314.2:p.L858R", // EGFR L858R
            "NP_000079.2:p.M1I",   // Start codon
            "NP_000079.2:p.R100*", // Nonsense with stop
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(
                result.is_ok(),
                "Failed to parse single-letter variant: {}",
                v
            );
            assert!(matches!(result.unwrap(), HgvsVariant::Protein(_)));
        }
    }

    #[test]
    fn test_protein_single_letter_display() {
        // Single-letter codes should display as three-letter codes
        let variant = parse_hgvs("NP_000079.2:p.V600E").unwrap();
        assert_eq!(variant.to_string(), "NP_000079.2:p.Val600Glu");

        let variant = parse_hgvs("NP_000138.2:p.G12D").unwrap();
        assert_eq!(variant.to_string(), "NP_000138.2:p.Gly12Asp");
    }
}

/// Test mitochondrial (m.) variant parsing
mod mitochondrial_variants {
    use super::*;

    #[test]
    fn test_mt_substitution() {
        let variants = vec!["NC_012920.1:m.3243A>G", "NC_012920.1:m.8344A>G"];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
            assert!(matches!(result.unwrap(), HgvsVariant::Mt(_)));
        }
    }

    #[test]
    fn test_mt_deletion() {
        let variants = vec!["NC_012920.1:m.8344del", "NC_012920.1:m.100_105del"];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }
}

/// Test accession formats
mod accession_formats {
    use super::*;

    #[test]
    fn test_ncbi_accessions() {
        let variants = vec![
            "NC_000001.11:g.100A>G",   // RefSeq genomic
            "NM_000088.3:c.100A>G",    // RefSeq mRNA
            "NR_000001.1:n.100A>G",    // RefSeq non-coding RNA
            "NP_000079.2:p.Val100Glu", // RefSeq protein
            "NG_012232.1:g.100A>G",    // RefSeq genomic region
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_lrg_accessions() {
        let variants = vec!["LRG_199:g.100A>G"];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_ensembl_accessions() {
        let variants = vec![
            "ENST00000380152.7:c.100A>G",
            "ENSG00000139618.15:g.100A>G",
            "ENSP00000369497.3:p.Val600Glu",
            "ENST00000012345:c.50del", // Without version
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_ensembl_roundtrip() {
        // Test that Ensembl accessions roundtrip correctly (no underscore)
        let variant = parse_hgvs("ENST00000380152.7:c.100A>G").unwrap();
        assert_eq!(variant.to_string(), "ENST00000380152.7:c.100A>G");

        let variant = parse_hgvs("ENSG00000139618.15:g.100A>G").unwrap();
        assert_eq!(variant.to_string(), "ENSG00000139618.15:g.100A>G");
    }
}

/// Test edge cases and special syntax
mod edge_cases {
    use super::*;

    #[test]
    fn test_gene_symbol_in_accession() {
        let variants = vec!["NG_012232.1(NM_004006.1):c.93+1G>T"];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse: {}", v);
        }
    }

    #[test]
    fn test_uncertain_positions() {
        // Note: Not all uncertainty notations are supported
        let variants = vec![
            "NC_000001.11:g.(100_200)del", // Uncertain range
        ];

        for v in variants {
            let result = parse_hgvs(v);
            // Just verify parsing doesn't panic
            let _ = result;
        }
    }
}

/// Tests based on HGVS nomenclature recommendations from hgvs-nomenclature.org
/// https://hgvs-nomenclature.org/stable/recommendations/
mod hgvs_recommendations {
    use super::*;

    /// Test that we follow the 3' rule for display (most 3' position preferred)
    /// Source: https://hgvs-nomenclature.org/stable/recommendations/general/
    #[test]
    fn test_3prime_rule_parsing() {
        // These variants should parse correctly - 3' shifting is a normalization concern
        let variants = vec![
            "NM_000088.3:c.10del",
            "NM_000088.3:c.10dup",
            "NM_000088.3:c.10_11insA",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse 3' rule example: {}", v);
        }
    }

    /// Test delins format compliance
    /// Source: https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/
    #[test]
    fn test_delins_format() {
        let variants = vec![
            // Single position delins
            "NC_000001.11:g.123delinsAC",
            // Range delins
            "NC_000001.11:g.123_129delinsAC",
            // CDS delins
            "NM_004006.2:c.6775_6777delinsC",
            // Codon replacement
            "LRG_199t1:c.145_147delinsTGG",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse delins: {}", v);
        }
    }

    /// Test protein extension variants
    /// Source: https://hgvs-nomenclature.org/stable/recommendations/protein/extension/
    #[test]
    fn test_protein_extension_format() {
        let variants = vec![
            // N-terminal extension
            "NP_003997.2:p.Met1ext-5",
            // C-terminal extension with Ter
            "NP_003997.2:p.Ter110extTer17",
            // C-terminal extension with asterisk
            "NP_003997.2:p.*110ext*17",
            // Unknown extension length
            "NP_003997.2:p.Ter327ext*?",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse extension: {}", v);
        }
    }

    /// Test frameshift variants
    /// Source: https://hgvs-nomenclature.org/stable/recommendations/protein/frameshift/
    #[test]
    fn test_frameshift_format() {
        let variants = vec![
            // Short form
            "NP_000079.2:p.Arg97fs",
            // Long form with termination position
            "NP_000079.2:p.Arg97fsTer23",
            // With asterisk
            "NP_000079.2:p.Arg97fs*23",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse frameshift: {}", v);
        }

        // Note: Full HGVS spec also supports "p.Arg97ProfsTer23" format
        // where Pro is the new amino acid, but we currently only store position+fs+ter_pos
    }

    /// Test RNA variant format (lowercase)
    /// Source: https://hgvs-nomenclature.org/stable/recommendations/RNA/substitution/
    #[test]
    fn test_rna_lowercase() {
        let variants = vec!["NM_004006.3:r.76a>c", "NM_004006.3:r.100u>a"];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse RNA: {}", v);
        }

        // Note: HGVS spec supports negative positions (r.-14a>c) for 5' UTR
        // and asterisk positions (r.*41u>a) for 3' UTR, but our RNA parser
        // currently uses simple positive positions only
    }

    /// Test that reference sequence versions are parsed correctly
    /// Source: https://hgvs-nomenclature.org/stable/recommendations/checklist/
    #[test]
    fn test_reference_versions() {
        let variants = vec![
            "NM_004006.2:c.123A>G",    // Version 2
            "NM_004006.3:c.123A>G",    // Version 3
            "NC_000001.11:g.12345A>G", // Chromosome version 11
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse versioned ref: {}", v);
        }
    }

    /// Test underscore range notation (not hyphen)
    /// Source: https://hgvs-nomenclature.org/stable/recommendations/checklist/
    #[test]
    fn test_underscore_ranges() {
        let variants = vec![
            "NM_000088.3:c.12_14del",
            "NC_000001.11:g.12345_12678del",
            "NP_000079.2:p.Lys23_Leu24del",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse range: {}", v);
        }
    }

    /// Test insertion format (between two positions)
    /// Source: https://hgvs-nomenclature.org/stable/recommendations/checklist/
    #[test]
    fn test_insertion_format() {
        // Insertions should be c.51_52insT, not c.52insT
        let variants = vec![
            "NM_000088.3:c.51_52insT",
            "NM_000088.3:c.100_101insATG",
            "NC_000001.11:g.12345_12346insGGG",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse insertion: {}", v);
        }
    }

    /// Test protein identity/synonymous format
    /// Source: https://hgvs-nomenclature.org/stable/recommendations/checklist/
    #[test]
    fn test_synonymous_format() {
        let variants = vec![
            "NP_000079.2:p.Leu54=",   // Synonymous (silent)
            "NP_000079.2:p.(Leu54=)", // Predicted synonymous
            "NP_000079.2:p.=",        // Whole protein identity
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse synonymous: {}", v);
        }
    }

    /// Test uncertain breakpoint notation for structural variants
    /// Source: https://hgvs-nomenclature.org/recommendations/uncertain/
    ///
    /// The (pos1_pos2)_(pos3_pos4) syntax describes variants with uncertain breakpoints:
    /// - pos1_pos4: maximal extent of the variant
    /// - pos2_pos3: minimal extent of the variant
    ///
    /// This notation is used when exact breakpoints haven't been sequenced (e.g., MLPA, PCR).
    /// Found in ClinGen Evidence Repository data.
    #[test]
    fn test_uncertain_breakpoint_ranges() {
        let variants = vec![
            // Uncertain intronic breakpoints - duplication spanning multiple exons
            // From ClinGen: NM_000527.5(LDLR):c.(313+1_314-1)_(1186+1_1187-1)dup
            // Breakpoints somewhere in introns flanking exons 4-11
            "NM_000527.5:c.(313+1_314-1)_(1186+1_1187-1)dup",
            // Uncertain deletion with unknown 5' breakpoint
            // The ? indicates position is completely unknown
            "NM_000088.3:c.(?_-187)_(190+1_191-1)del",
            // Genomic uncertain breakpoints (from HGVS spec example)
            "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup",
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(
                result.is_ok(),
                "Failed to parse uncertain breakpoint range: {} - {:?}",
                v,
                result.err()
            );
        }
    }
}

/// Summary test that counts coverage
#[test]
fn test_hgvs_standard_coverage() {
    println!("\nHGVS Standard Coverage:");
    println!("  Genomic (g.): substitution, deletion, duplication, insertion, delins, inversion");
    println!("  Coding (c.): substitution, intronic, UTR, deletion, duplication, insertion, delins, inversion");
    println!("  Transcript (n.): substitution, deletion");
    println!("  Protein (p.): substitution, deletion, frameshift, extension, identity");
    println!("  RNA (r.): substitution (lowercase)");
    println!("  Mitochondrial (m.): substitution, deletion");
    println!("  Accessions: NCBI (NC/NM/NR/NP/NG), LRG, Ensembl (ENST/ENSG)");
    println!("  HGVS Recommendations: 3' rule, delins, extensions, frameshifts");
}
