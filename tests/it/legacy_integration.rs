//! Integration tests for legacy transcript handling.
//!
//! Tests the preparation pipeline for legacy transcripts, including:
//! - Detection of legacy versions from patterns
//! - GenBank format parsing
//! - CDS metadata extraction
//! - Legacy transcript caching in ferro reference

#[cfg(feature = "benchmark")]
mod legacy_tests {
    use ferro_hgvs::benchmark::{LegacyMetadata, LegacyTranscript};
    use std::io::Write;
    use tempfile::TempDir;

    #[test]
    fn test_legacy_metadata_serialization() {
        // Verify metadata JSON structure matches expected format
        let metadata = LegacyMetadata {
            generated_at: "2024-01-01T00:00:00Z".to_string(),
            transcripts: [(
                "NM_000051.3".to_string(),
                LegacyTranscript {
                    id: "NM_000051.3".to_string(),
                    gene_symbol: Some("ATM".to_string()),
                    cds_start: Some(250),
                    cds_end: Some(9500),
                    sequence_length: 10000,
                },
            )]
            .into_iter()
            .collect(),
        };

        let json = serde_json::to_string_pretty(&metadata).unwrap();
        assert!(json.contains("NM_000051.3"));
        assert!(json.contains("ATM"));
        assert!(json.contains("cds_start"));
        assert!(json.contains("cds_end"));
    }

    #[test]
    fn test_legacy_metadata_deserialization() {
        let json = r#"{
            "generated_at": "2024-01-01T00:00:00Z",
            "transcripts": {
                "NM_000051.3": {
                    "id": "NM_000051.3",
                    "gene_symbol": "ATM",
                    "cds_start": 250,
                    "cds_end": 9500,
                    "sequence_length": 10000
                }
            }
        }"#;

        let metadata: LegacyMetadata = serde_json::from_str(json).unwrap();
        assert_eq!(metadata.transcripts.len(), 1);
        assert!(metadata.transcripts.contains_key("NM_000051.3"));

        let transcript = &metadata.transcripts["NM_000051.3"];
        assert_eq!(transcript.gene_symbol, Some("ATM".to_string()));
        assert_eq!(transcript.cds_start, Some(250));
        assert_eq!(transcript.cds_end, Some(9500));
        assert_eq!(transcript.sequence_length, 10000);
    }

    #[test]
    fn test_legacy_transcript_optional_fields() {
        // Test that optional fields can be None
        let transcript = LegacyTranscript {
            id: "NM_000001.1".to_string(),
            gene_symbol: None,
            cds_start: None,
            cds_end: None,
            sequence_length: 5000,
        };

        let json = serde_json::to_string(&transcript).unwrap();
        let parsed: LegacyTranscript = serde_json::from_str(&json).unwrap();

        assert_eq!(parsed.gene_symbol, None);
        assert_eq!(parsed.cds_start, None);
        assert_eq!(parsed.cds_end, None);
    }

    #[test]
    fn test_accession_sources_from_fai() {
        use ferro_hgvs::benchmark::AccessionSources;

        let dir = TempDir::new().unwrap();

        // Create subdirectories
        std::fs::create_dir_all(dir.path().join("transcripts")).unwrap();
        std::fs::create_dir_all(dir.path().join("supplemental")).unwrap();
        std::fs::create_dir_all(dir.path().join("genome")).unwrap();

        // Create mock FAI files
        let transcript_fai = dir.path().join("transcripts/refseq.fna.fai");
        let mut f = std::fs::File::create(&transcript_fai).unwrap();
        writeln!(f, "NM_000001.1\t1000\t0\t80\t81").unwrap();
        writeln!(f, "NM_000002.2\t2000\t1000\t80\t81").unwrap();

        let supplemental_fai = dir.path().join("supplemental/legacy.fna.fai");
        let mut f = std::fs::File::create(&supplemental_fai).unwrap();
        writeln!(f, "U31929.1\t500\t0\t80\t81").unwrap();

        let sources = AccessionSources::from_ferro_reference(dir.path()).unwrap();

        assert!(sources.contains("NM_000001.1"));
        assert!(sources.contains("NM_000002.2"));
        assert!(sources.contains("U31929.1"));
        assert!(!sources.contains("NM_999999.1"));

        assert_eq!(sources.source_for("NM_000001.1"), Some("transcripts"));
        assert_eq!(sources.source_for("U31929.1"), Some("supplemental"));
    }

    #[test]
    fn test_accession_sources_find_missing() {
        use ferro_hgvs::benchmark::AccessionSources;

        let mut sources = AccessionSources::default();
        sources.transcripts.insert("NM_000001.1".to_string());
        sources.transcripts.insert("NM_000002.2".to_string());

        let to_check = vec![
            "NM_000001.1".to_string(),
            "NM_000002.2".to_string(),
            "NM_000003.3".to_string(),
        ];
        let missing = sources.find_missing(&to_check);

        assert_eq!(missing.len(), 1);
        assert_eq!(missing[0], "NM_000003.3");
    }

    #[test]
    fn test_accession_sources_total_count_deduplicates() {
        use ferro_hgvs::benchmark::AccessionSources;

        let mut sources = AccessionSources::default();
        sources.transcripts.insert("NM_000001.1".to_string());
        sources.supplemental.insert("NM_000001.1".to_string()); // Same as transcripts
        sources.supplemental.insert("U31929.1".to_string());

        // Should deduplicate NM_000001.1
        assert_eq!(sources.total_count(), 2);
    }

    #[test]
    fn test_translate_cds_to_protein() {
        use ferro_hgvs::benchmark::translate::translate_cds_to_protein;

        // ATG = Met, TGG = Trp, TAA = Stop
        assert_eq!(
            translate_cds_to_protein("ATGTGGTAA"),
            Some("MW".to_string())
        );

        // No stop codon - translates all codons
        assert_eq!(translate_cds_to_protein("ATGTGG"), Some("MW".to_string()));

        // Incomplete codon at end is ignored
        assert_eq!(translate_cds_to_protein("ATGTGGA"), Some("MW".to_string()));
    }

    #[test]
    fn test_derive_protein_from_transcript() {
        use ferro_hgvs::benchmark::translate::derive_protein_from_transcript;

        // 5'UTR (3bp) + CDS (9bp) + 3'UTR (3bp)
        let transcript = "AAAATGTGGTAAGGG";
        assert_eq!(
            derive_protein_from_transcript(transcript, 3, 12),
            Some("MW".to_string())
        );

        // Invalid coordinates
        assert_eq!(derive_protein_from_transcript(transcript, 12, 3), None);
        assert_eq!(derive_protein_from_transcript(transcript, 0, 100), None);
    }

    #[test]
    fn test_cds_coords_conversion() {
        use ferro_hgvs::benchmark::translate::cds_coords_to_indices;

        // 1-based inclusive to 0-based half-open
        assert_eq!(cds_coords_to_indices(1, 100), (0, 100));
        assert_eq!(cds_coords_to_indices(251, 9500), (250, 9500));
    }
}
