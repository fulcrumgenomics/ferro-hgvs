//! Integration tests for legacy transcript handling.
//!
//! Tests the preparation pipeline for legacy transcripts, including:
//! - Detection of legacy versions from patterns
//! - GenBank format parsing
//! - CDS metadata extraction
//! - Legacy transcript caching in ferro reference

#[cfg(feature = "benchmark")]
mod legacy_tests {
    use ferro_hgvs::benchmark::{LegacyMetadata, LegacyTranscriptRecord};
    use std::io::Write;
    use tempfile::TempDir;

    #[test]
    fn test_legacy_metadata_serialization() {
        // Verify metadata JSON structure matches expected format
        let metadata = LegacyMetadata {
            generated_at: "2024-01-01T00:00:00Z".to_string(),
            transcripts: [(
                "NM_000051.3".to_string(),
                LegacyTranscriptRecord {
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

        let record = &metadata.transcripts["NM_000051.3"];
        assert_eq!(record.gene_symbol, Some("ATM".to_string()));
        assert_eq!(record.cds_start, Some(250));
        assert_eq!(record.cds_end, Some(9500));
        assert_eq!(record.sequence_length, 10000);
    }

    #[test]
    fn test_legacy_transcript_optional_fields() {
        // Test that optional fields can be None on the on-disk record.
        let record = LegacyTranscriptRecord {
            id: "NM_000001.1".to_string(),
            gene_symbol: None,
            cds_start: None,
            cds_end: None,
            sequence_length: 5000,
        };

        let json = serde_json::to_string(&record).unwrap();
        let parsed: LegacyTranscriptRecord = serde_json::from_str(&json).unwrap();

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

/// The sourced/derived split for legacy-transcript metadata (#2085).
///
/// These pin the structural invariant introduced by the first increment of
/// #2085: the *authority* type [`LegacyTranscript`] carries only sourced fields
/// and structurally cannot hold the derived `sequence_length`, while the on-disk
/// [`LegacyTranscriptRecord`] remains the sole carrier of the derived length so
/// the serialized bytes — and therefore the content stamps and the blessed
/// reference identity (`aa8b3246d83055cc`) — do not move.
#[cfg(feature = "benchmark")]
mod sourced_derived_split {
    use ferro_hgvs::benchmark::{LegacyMetadata, LegacyTranscript, LegacyTranscriptRecord};
    use ferro_hgvs::prepare::identity::fnv1a_hex;

    /// The exact pre-split serialization of a single-record metadata file, so the
    /// identity-neutrality assertions below compare against a fixed golden rather
    /// than against the code under test. A single transcript keeps the map order
    /// deterministic; `cds_end: null` pins that `None` still serializes as before.
    const GOLDEN_PRETTY: &str = "{\n  \"generated_at\": \"2024-01-01T00:00:00Z\",\n  \"transcripts\": {\n    \"NM_000051.3\": {\n      \"id\": \"NM_000051.3\",\n      \"gene_symbol\": \"ATM\",\n      \"cds_start\": 250,\n      \"cds_end\": null,\n      \"sequence_length\": 10000\n    }\n  }\n}";

    fn golden_metadata() -> LegacyMetadata {
        LegacyMetadata {
            generated_at: "2024-01-01T00:00:00Z".to_string(),
            transcripts: [(
                "NM_000051.3".to_string(),
                LegacyTranscriptRecord {
                    id: "NM_000051.3".to_string(),
                    gene_symbol: Some("ATM".to_string()),
                    cds_start: Some(250),
                    cds_end: None,
                    sequence_length: 10000,
                },
            )]
            .into_iter()
            .collect(),
        }
    }

    /// The sourced authority type cannot carry a derived length: it has no such
    /// field, so it never serializes one. This is the structural core of #2085 —
    /// a persisted length can no longer masquerade as a sourced fact.
    #[test]
    fn sourced_authority_type_has_no_length_field() {
        let sourced = LegacyTranscript {
            id: "NM_000051.3".to_string(),
            gene_symbol: Some("ATM".to_string()),
            cds_start: Some(250),
            cds_end: Some(9500),
        };
        let json = serde_json::to_string(&sourced).unwrap();
        assert!(
            !json.contains("sequence_length"),
            "sourced LegacyTranscript must not serialize a derived length: {json}"
        );
        assert!(json.contains("\"id\":\"NM_000051.3\""));
        assert!(json.contains("\"cds_start\":250"));
    }

    /// The on-disk bytes are byte-for-byte unchanged from the pre-split format,
    /// so the derived-artifact content stamps and the blessed reference identity
    /// do not move. This is what makes the increment identity-neutral (#905/#933).
    #[test]
    fn on_disk_bytes_are_identity_neutral() {
        let serialized = serde_json::to_string_pretty(&golden_metadata()).unwrap();
        assert_eq!(serialized, GOLDEN_PRETTY);
    }

    /// The FNV-1a content stamp `src/prepare/identity.rs` injects into the
    /// reference identity, pinned as a literal digest.
    ///
    /// The literal is what makes this falsifiable. Comparing
    /// `fnv1a_hex(serialized)` against `fnv1a_hex(GOLDEN_PRETTY)` — which this
    /// test used to do — is `f(x) == f(x)` once the sibling above has asserted
    /// `serialized == GOLDEN_PRETTY`, so it holds for *any* hash function and
    /// any golden. The digest below was derived independently of ferro (FNV-1a
    /// 64-bit over the golden's 228 bytes, offset basis `cbf29ce484222325`,
    /// prime `100000001b3`), so it pins the stamp rather than restating the
    /// assertion above it.
    #[test]
    fn content_stamp_is_the_pinned_digest() {
        let serialized = serde_json::to_string_pretty(&golden_metadata()).unwrap();
        assert_eq!(fnv1a_hex(serialized.as_bytes()), "3135b86e16384d04");
    }

    /// The persisted `sequence_length` is exactly the FASTA-derived length
    /// (`seq.len()`), i.e. fully redundant with its authority — the very property
    /// that makes it derived. `into_record` pairs the sourced authority with that
    /// derived length; `sourced()` projects back, dropping it.
    #[test]
    fn persisted_length_is_the_fasta_derived_length() {
        let seq = "ACGTACGTAC"; // stands in for the FASTA bases (authority)
        let sourced = LegacyTranscript {
            id: "NM_000001.1".to_string(),
            gene_symbol: None,
            cds_start: Some(1),
            cds_end: Some(9),
        };
        let record = sourced.clone().into_record(seq.len());

        // The persisted copy equals the authority's length, byte-for-byte.
        assert_eq!(record.sequence_length, seq.len());
        // Projecting back to the authority recovers exactly the sourced fields,
        // with no length surviving the round trip.
        assert_eq!(record.sourced(), sourced);
    }

    /// `sourced_transcripts()` is the length-free view consumers should read, so
    /// a persisted length can never be mistaken for the source of truth (#2085).
    #[test]
    fn sourced_transcripts_view_drops_the_derived_length() {
        let metadata = golden_metadata();
        let sourced = metadata.sourced_transcripts();
        assert_eq!(sourced.len(), 1);
        let view = &sourced["NM_000051.3"];
        assert_eq!(view.id, "NM_000051.3");
        assert_eq!(view.gene_symbol, Some("ATM".to_string()));
        assert_eq!(view.cds_start, Some(250));
        assert_eq!(view.cds_end, None);
        // A serialized sourced view carries no length key.
        assert!(!serde_json::to_string(view)
            .unwrap()
            .contains("sequence_length"));
    }
}
