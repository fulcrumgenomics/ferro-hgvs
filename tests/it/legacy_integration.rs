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
        // The derived length is no longer persisted (#2064).
        assert!(!json.contains("sequence_length"));
    }

    #[test]
    fn test_legacy_metadata_deserialization() {
        // A pre-#2064 on-disk file still carries `sequence_length`. It must still
        // load: serde ignores the now-unknown field rather than failing, so an
        // older reference is readable by the field-free type.
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

        let record = &metadata.transcripts["NM_000051.3"];
        assert_eq!(record.gene_symbol, Some("ATM".to_string()));
        assert_eq!(record.cds_start, Some(250));
        assert_eq!(record.cds_end, Some(9500));
    }

    #[test]
    fn test_legacy_transcript_optional_fields() {
        // Test that optional fields can be None on the on-disk record.
        let record = LegacyTranscript {
            id: "NM_000001.1".to_string(),
            gene_symbol: None,
            cds_start: None,
            cds_end: None,
        };

        let json = serde_json::to_string(&record).unwrap();
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

/// The derived length is no longer persisted for legacy-transcript metadata
/// (#2064, completing #2085 move 3).
///
/// These pin that the on-disk row [`LegacyTranscript`] carries only sourced
/// fields (accession, gene symbol, CDS bounds) and can no longer represent the
/// derived `sequence_length` — the field is *unrepresentable*, not merely
/// unread. Removing it changes the serialized bytes of the content-stamped
/// legacy metadata artifacts, and so deliberately BUMPS the blessed reference
/// identity from `aa8b3246d83055cc` to `f329e62bd95d8020` (#905/#933). This is
/// the opposite property to #2087, whose on-disk bytes were identity-neutral.
#[cfg(feature = "benchmark")]
mod legacy_metadata_no_derived_length {
    use ferro_hgvs::benchmark::{LegacyMetadata, LegacyTranscript};
    use ferro_hgvs::prepare::identity::fnv1a_hex;

    /// The exact post-#2064 serialization of a single-record metadata file, so
    /// the byte/stamp assertions below compare against a fixed golden rather than
    /// against the code under test. A single transcript keeps the map order
    /// deterministic; `cds_end: null` pins that `None` still serializes, and the
    /// absence of any `sequence_length` line is the change this whole PR is about.
    const GOLDEN_PRETTY: &str = "{\n  \"generated_at\": \"2024-01-01T00:00:00Z\",\n  \"transcripts\": {\n    \"NM_000051.3\": {\n      \"id\": \"NM_000051.3\",\n      \"gene_symbol\": \"ATM\",\n      \"cds_start\": 250,\n      \"cds_end\": null\n    }\n  }\n}";

    fn golden_metadata() -> LegacyMetadata {
        LegacyMetadata {
            generated_at: "2024-01-01T00:00:00Z".to_string(),
            transcripts: [(
                "NM_000051.3".to_string(),
                LegacyTranscript {
                    id: "NM_000051.3".to_string(),
                    gene_symbol: Some("ATM".to_string()),
                    cds_start: Some(250),
                    cds_end: None,
                },
            )]
            .into_iter()
            .collect(),
        }
    }

    /// The on-disk row cannot carry a derived length: it has no such field, so it
    /// never serializes one. This is the structural core of #2064 — a persisted
    /// length is no longer representable, not merely unread.
    #[test]
    fn legacy_metadata_carries_no_derived_length() {
        let record = LegacyTranscript {
            id: "NM_000051.3".to_string(),
            gene_symbol: Some("ATM".to_string()),
            cds_start: Some(250),
            cds_end: Some(9500),
        };
        let json = serde_json::to_string(&record).unwrap();
        assert!(
            !json.contains("sequence_length"),
            "LegacyTranscript must not serialize a derived length: {json}"
        );
        assert!(json.contains("\"id\":\"NM_000051.3\""));
        assert!(json.contains("\"cds_start\":250"));
    }

    /// The exact on-disk bytes after the field removal, pinned as a golden. This
    /// is the byte change that moves the content stamp and hence the identity.
    #[test]
    fn on_disk_bytes_drop_the_length_field() {
        let serialized = serde_json::to_string_pretty(&golden_metadata()).unwrap();
        assert_eq!(serialized, GOLDEN_PRETTY);
        assert!(!serialized.contains("sequence_length"));
    }

    /// The FNV-1a content stamp `src/prepare/identity.rs` injects into the
    /// reference identity, pinned as a literal digest.
    ///
    /// The literal is what makes this falsifiable. Comparing
    /// `fnv1a_hex(serialized)` against `fnv1a_hex(GOLDEN_PRETTY)` would be
    /// `f(x) == f(x)` once the sibling above has asserted `serialized ==
    /// GOLDEN_PRETTY`, so it would hold for *any* hash function and any golden.
    /// The digest below was derived independently of ferro (FNV-1a 64-bit over
    /// the golden's bytes, offset basis `cbf29ce484222325`, prime
    /// `100000001b3`), so it pins the stamp rather than restating the assertion
    /// above it. It differs from the pre-#2064 digest (`3135b86e16384d04`)
    /// because the removed `sequence_length` line changed the bytes.
    #[test]
    fn content_stamp_is_the_pinned_digest() {
        let serialized = serde_json::to_string_pretty(&golden_metadata()).unwrap();
        assert_eq!(fnv1a_hex(serialized.as_bytes()), "bdd80d2adb01c7c1");
    }
}
