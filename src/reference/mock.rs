//! Mock reference provider for testing

use crate::error::FerroError;
use crate::reference::provider::ReferenceProvider;
use crate::reference::transcript::Transcript;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::OnceLock;

/// Mock reference provider that loads transcripts from JSON
#[derive(Clone)]
pub struct MockProvider {
    transcripts: HashMap<String, Transcript>,
    proteins: HashMap<String, String>,
    /// Genomic sequences keyed by contig name
    genomic_sequences: HashMap<String, String>,
    /// Transcript ids this provider should report as *not* version-exact from
    /// [`ReferenceProvider::has_transcript_version_exact`]. Lets tests exercise
    /// the protein path's cross-version gate (#505) without a full FASTA/cdot
    /// manifest. Empty by default, so an unconfigured mock reports every
    /// transcript as version-exact (matching the trait default).
    non_version_exact: HashSet<String>,
}

impl MockProvider {
    /// Create an empty mock provider
    pub fn new() -> Self {
        Self {
            transcripts: HashMap::new(),
            proteins: HashMap::new(),
            genomic_sequences: HashMap::new(),
            non_version_exact: HashSet::new(),
        }
    }

    /// Mark `id` as not available at its exact requested version, so
    /// [`ReferenceProvider::has_transcript_version_exact`] returns `false` for
    /// it. Simulates a reference that only carries a sibling version's bases.
    pub fn mark_non_version_exact(&mut self, id: impl Into<String>) {
        self.non_version_exact.insert(id.into());
    }

    /// Load reference data from a JSON file.
    ///
    /// Accepts either a bare array of `Transcript` records or an object
    /// of the form `{ transcripts, proteins, genomic_sequences }`.
    pub fn from_json(path: &Path) -> Result<Self, FerroError> {
        // `deny_unknown_fields` so a typo'd key (e.g. `transripts`) produces
        // a clear error rather than silently defaulting to an empty provider.
        #[derive(serde::Deserialize)]
        #[serde(deny_unknown_fields)]
        struct ObjectForm {
            #[serde(default)]
            transcripts: Vec<Transcript>,
            #[serde(default)]
            proteins: HashMap<String, String>,
            #[serde(default)]
            genomic_sequences: HashMap<String, String>,
            /// Container metadata emitted by `ferro convert-gff`; accepted but ignored.
            #[serde(default, rename = "version")]
            _version: Option<String>,
            /// Container metadata emitted by `ferro convert-gff`; per-transcript
            /// `genome_build` on each `Transcript` is authoritative.
            #[serde(default, rename = "genome_build")]
            _genome_build: Option<String>,
        }

        let content = std::fs::read_to_string(path)?;
        let value: serde_json::Value = serde_json::from_str(&content)?;

        let (transcripts, proteins, genomic_sequences) = match value {
            serde_json::Value::Array(_) => {
                let transcripts: Vec<Transcript> = serde_json::from_value(value)?;
                (transcripts, HashMap::new(), HashMap::new())
            }
            serde_json::Value::Object(_) => {
                let obj: ObjectForm = serde_json::from_value(value)?;
                (obj.transcripts, obj.proteins, obj.genomic_sequences)
            }
            _ => {
                return Err(FerroError::Json {
                    msg: "MockProvider JSON root must be an array or object".to_string(),
                })
            }
        };

        let map: HashMap<String, Transcript> = transcripts
            .into_iter()
            .map(|tx| (tx.id.clone(), tx))
            .collect();

        Ok(Self {
            transcripts: map,
            proteins,
            genomic_sequences,
            non_version_exact: HashSet::new(),
        })
    }

    /// Add a transcript to the provider
    pub fn add_transcript(&mut self, transcript: Transcript) {
        self.transcripts.insert(transcript.id.clone(), transcript);
    }

    /// Add a protein sequence to the provider
    pub fn add_protein(&mut self, accession: impl Into<String>, sequence: impl Into<String>) {
        self.proteins.insert(accession.into(), sequence.into());
    }

    /// Add a genomic sequence for a contig/chromosome
    pub fn add_genomic_sequence(&mut self, contig: impl Into<String>, sequence: impl Into<String>) {
        self.genomic_sequences
            .insert(contig.into(), sequence.into());
    }

    /// Create a provider with some test transcripts
    pub fn with_test_data() -> Self {
        use crate::reference::transcript::{Exon, ManeStatus, Strand};

        let mut provider = Self::new();

        // Add a typical coding transcript (MANE Select)
        provider.add_transcript(Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some(
                "ATGCCCAAGGTGCTGCCCCAGATGCTGCCAGTGCTGCTGCTGCTGCTGCTGCTGCTGCTG".to_string(),
            ),
            cds_start: Some(1),
            cds_end: Some(60),
            exons: vec![
                Exon::new(1, 1, 20),
                Exon::new(2, 21, 40),
                Exon::new(3, 41, 60),
            ],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        // Add a transcript with UTRs
        provider.add_transcript(Transcript {
            id: "NM_001234.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("AAAAATGCCCAAGGGGGGGGGGGGGGGGGGGGGGGGGTAAAAAA".to_string()),
            cds_start: Some(5),
            cds_end: Some(38),
            exons: vec![
                Exon::new(1, 1, 15),
                Exon::new(2, 16, 30),
                Exon::new(3, 31, 44),
            ],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        // Add a minus strand transcript
        provider.add_transcript(Transcript {
            id: "NM_999999.1".to_string(),
            gene_symbol: Some("MINUS".to_string()),
            strand: Strand::Minus,
            sequence: Some("ATGCATGCATGCATGCATGCATGCATGCATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(30),
            exons: vec![Exon::new(1, 1, 32)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        // Add a transcript for testing 3' normalization of duplications
        // Sequence has "GAA" at c.8-c.9 to test that c.8dup -> c.9dup (2 A's)
        // and "GTTT" at c.21-c.24 to test that c.22dup -> c.24dup (3 T's)
        // Positions: 1234567890123456789012345678901234567890
        // Sequence:  ATGCCCGAAGCCCCCCCCCGTTTGCATGCATGCATGCAT
        //                  ^^           ^^^
        //                  c.8-9 (AA)   c.22-24 (TTT)
        provider.add_transcript(Transcript {
            id: "NM_888888.1".to_string(),
            gene_symbol: Some("DUPTEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGCCCGAAGCCCCCCCCCGTTTGCATGCATGCATGCAT".to_string()),
            cds_start: Some(1),
            cds_end: Some(39),
            exons: vec![Exon::new(1, 1, 39)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        // Add test protein sequences for validation testing
        // NP_000079.2 is used for BRAF V600E-style testing
        // We need Val at position 600, Arg at 97 for frameshift tests
        // Building a sequence with known positions:
        // - Position 1: M (Met)
        // - Position 97: R (Arg) for frameshift test
        // - Position 600: V (Val) for V600E test
        // - Position 23-25: K, A, E for range tests
        let mut seq = String::new();
        // Positions 1-96 (96 chars): just padding with Ala
        seq.push('M'); // Position 1
        for _ in 2..23 {
            seq.push('A'); // Positions 2-22
        }
        seq.push('K'); // Position 23 - for Lys23
        seq.push('A'); // Position 24 - for middle of range
        seq.push('E'); // Position 25 - for Glu25
        for _ in 26..97 {
            seq.push('A'); // Positions 26-96
        }
        seq.push('R'); // Position 97 - for Arg97 frameshift
        for _ in 98..600 {
            seq.push('A'); // Positions 98-599
        }
        seq.push('V'); // Position 600 - for Val600
        for _ in 601..=700 {
            seq.push('A'); // Positions 601-700
        }
        provider.add_protein("NP_000079.2", seq);

        // A simpler test protein for easier position testing
        // Position 1=M, 2=V, 3=L, 4=S, 5=P, etc.
        provider.add_protein("NP_TEST.1", "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFL");

        provider
    }

    /// Get the number of transcripts
    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    /// Check if provider is empty
    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }

    /// Get all transcript IDs
    pub fn transcript_ids(&self) -> Vec<&str> {
        self.transcripts.keys().map(|s| s.as_str()).collect()
    }

    /// Iterate over all transcripts registered in the provider.
    pub fn all_transcripts(&self) -> impl Iterator<Item = &Transcript> {
        self.transcripts.values()
    }

    /// True when `id` is declared as a chromosome/contig name — either an
    /// entry in `genomic_sequences` or the `chromosome` field on any
    /// transcript. Used by `get_sequence` to disambiguate genomic lookups
    /// from transcript-id lookups when the two name spaces collide
    /// (e.g. a transcript whose id is `chr1-gene.1`; see #311).
    fn is_known_contig(&self, id: &str) -> bool {
        if self.genomic_sequences.contains_key(id) {
            return true;
        }
        self.transcripts
            .values()
            .any(|tx| tx.chromosome.as_deref() == Some(id))
    }
}

impl Default for MockProvider {
    fn default() -> Self {
        Self::new()
    }
}

impl ReferenceProvider for MockProvider {
    fn get_transcript(&self, id: &str) -> Result<Transcript, FerroError> {
        // Handle versioned and unversioned lookups
        if let Some(tx) = self.transcripts.get(id) {
            return Ok(tx.clone());
        }

        // Try without version: only match a stored key whose base id (the
        // segment before the version dot) equals the requested id. A bare
        // `starts_with` would also match unrelated keys that merely share
        // a prefix (e.g. `chr1-gene.1` for a `chr1` lookup), causing a
        // genomic accession to be resolved as a transcript — see #311.
        //
        // Gated on an unversioned query: a versioned miss (`NM_000088.5`
        // against a stored `NM_000088.3`) must NOT silently return a
        // different version. The fallback only bridges `NM_000088` →
        // `NM_000088.3`.
        if !id.contains('.') {
            for (key, tx) in &self.transcripts {
                if key.split('.').next().unwrap_or(key) == id {
                    return Ok(tx.clone());
                }
            }
        }

        Err(FerroError::ReferenceNotFound { id: id.to_string() })
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        // If the requested id is registered as a contig — either as a key
        // in `genomic_sequences` or as the `chromosome` field on any
        // transcript — resolve it strictly against the genomic path.
        // Without this dispatch, a transcript whose own id (or, before
        // #311, whose id merely shared a prefix with the contig name)
        // would short-circuit the lookup and silently interpret genomic
        // coordinates as transcript-relative indices.
        if self.is_known_contig(id) {
            return self.get_genomic_sequence(id, start, end);
        }

        // Otherwise treat the id as a transcript accession (matches
        // FastaProvider's transcript-first ordering for transcript ids).
        if let Ok(transcript) = self.get_transcript(id) {
            return transcript
                .get_sequence(start, end)
                .map(|s| s.to_string())
                .ok_or_else(|| FerroError::InvalidCoordinates {
                    msg: format!("Position {}-{} out of range for {}", start, end, id),
                });
        }

        // Fall through to contig/chromosome lookup so genomic accessions
        // resolve against `genomic_sequences` instead of returning
        // ReferenceNotFound.
        self.get_genomic_sequence(id, start, end)
    }

    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        // Handle versioned and unversioned lookups
        let protein_seq = if let Some(seq) = self.proteins.get(accession) {
            seq.clone()
        } else if !accession.contains('.') {
            // Try without version: match a stored key whose base accession
            // (the segment before the version dot) equals the requested
            // accession. Strict base-id equality avoids the prefix-collision
            // failure mode that `get_transcript` had pre-#311.
            //
            // Gated on an unversioned query — a versioned miss must not
            // cross-version fall back (mirrors `get_transcript` above).
            let mut found = None;
            for (key, seq) in &self.proteins {
                if key.split('.').next().unwrap_or(key) == accession {
                    found = Some(seq.clone());
                    break;
                }
            }
            found.ok_or_else(|| FerroError::ProteinReferenceNotAvailable {
                accession: accession.to_string(),
                start,
                end,
            })?
        } else {
            return Err(FerroError::ProteinReferenceNotAvailable {
                accession: accession.to_string(),
                start,
                end,
            });
        };

        let start = start as usize;
        let end = end as usize;

        if start >= protein_seq.len() || end > protein_seq.len() || start > end {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "Position {}-{} out of range for {} (length {})",
                    start,
                    end,
                    accession,
                    protein_seq.len()
                ),
            });
        }

        Ok(protein_seq[start..end].to_string())
    }

    fn get_protein_length(&self, accession: &str) -> Result<u64, FerroError> {
        // Mirror the accession resolution in `get_protein_sequence`:
        // exact (versioned or stored) match first, then an unversioned
        // base-accession fallback, returning the stored length directly.
        if let Some(seq) = self.proteins.get(accession) {
            return Ok(seq.len() as u64);
        }
        if !accession.contains('.') {
            for (key, seq) in &self.proteins {
                if key.split('.').next().unwrap_or(key) == accession {
                    return Ok(seq.len() as u64);
                }
            }
        }
        // An accession this provider cannot resolve maps to length `0`,
        // matching the trait contract (the default probe leaves `lo = 0`
        // for an unresolvable protein). Returning an error here would
        // diverge from that contract.
        Ok(0)
    }

    fn has_protein_data(&self) -> bool {
        !self.proteins.is_empty()
    }

    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        let genomic_seq = self.genomic_sequences.get(contig).ok_or_else(|| {
            FerroError::GenomicReferenceNotAvailable {
                contig: contig.to_string(),
                start,
                end,
            }
        })?;

        let start = start as usize;
        let end = end as usize;

        if start >= genomic_seq.len() || end > genomic_seq.len() || start > end {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "Position {}-{} out of range for {} (length {})",
                    start,
                    end,
                    contig,
                    genomic_seq.len()
                ),
            });
        }

        Ok(genomic_seq[start..end].to_string())
    }

    fn has_genomic_data(&self) -> bool {
        !self.genomic_sequences.is_empty()
    }

    fn get_sequence_length(&self, id: &str) -> Result<u64, FerroError> {
        self.genomic_sequences
            .get(id)
            .map(|s| s.len() as u64)
            .ok_or_else(|| FerroError::ReferenceNotFound { id: id.to_string() })
    }

    fn has_transcript_version_exact(&self, id: &str) -> bool {
        !self.non_version_exact.contains(id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mock_provider_with_test_data() {
        let provider = MockProvider::with_test_data();
        assert!(!provider.is_empty());
        assert!(provider.len() >= 2);
    }

    #[test]
    fn test_get_transcript() {
        let provider = MockProvider::with_test_data();
        let tx = provider.get_transcript("NM_000088.3").unwrap();
        assert_eq!(tx.gene_symbol, Some("COL1A1".to_string()));
    }

    #[test]
    fn test_get_transcript_not_found() {
        let provider = MockProvider::with_test_data();
        let result = provider.get_transcript("NM_NONEXISTENT.1");
        assert!(result.is_err());
    }

    #[test]
    fn test_get_sequence() {
        let provider = MockProvider::with_test_data();
        let seq = provider.get_sequence("NM_000088.3", 0, 3).unwrap();
        assert_eq!(seq, "ATG");
    }

    #[test]
    fn test_get_sequence_falls_through_to_contig() {
        // Regression: get_sequence should fall through to contig lookup
        // when the id is not a transcript, matching FastaProvider behavior.
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("chr1", "ACGTACGT");
        let seq = provider.get_sequence("chr1", 0, 4).unwrap();
        assert_eq!(seq, "ACGT");
    }

    #[test]
    fn test_has_transcript() {
        let provider = MockProvider::with_test_data();
        assert!(provider.has_transcript("NM_000088.3"));
        assert!(!provider.has_transcript("NONEXISTENT"));
    }

    #[test]
    fn boxed_provider_forwards_has_transcript_version_exact() {
        // The protein gate (#505) is consulted through the provider type the
        // production CLI/benchmark paths use: `Box<dyn ReferenceProvider>`. If
        // the blanket boxed impl does not forward `has_transcript_version_exact`,
        // it silently falls through to the trait default `true` and the gate is
        // inert in production. Pin the forwarding here.
        let mut inner = MockProvider::new();
        inner.mark_non_version_exact("NM_TEST.1");
        let boxed: Box<dyn ReferenceProvider> = Box::new(inner);
        assert!(
            !boxed.has_transcript_version_exact("NM_TEST.1"),
            "boxed provider must forward to the inner provider's verdict (false), \
             not the trait default (true)"
        );
        assert!(boxed.has_transcript_version_exact("NM_OTHER.1"));
    }

    #[test]
    fn test_from_json_object_form_with_proteins_and_genomic() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"{
      "transcripts": [{
        "id": "NM_TEST.1",
        "gene_symbol": "TEST",
        "strand": "+",
        "sequence": "ATGCATGCAT",
        "cds_start": 1,
        "cds_end": 10,
        "exons": [{"number": 1, "start": 1, "end": 10}]
      }],
      "proteins": {
        "NP_TEST.1": "MAPLE"
      },
      "genomic_sequences": {
        "chr1": "ACGTACGTACGT"
      }
    }"#;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider = MockProvider::from_json(file.path()).unwrap();

        assert!(provider.has_transcript("NM_TEST.1"));
        assert!(provider.has_protein_data());
        assert_eq!(
            provider.get_protein_sequence("NP_TEST.1", 0, 5).unwrap(),
            "MAPLE"
        );
        assert!(provider.has_genomic_data());
        assert_eq!(provider.get_genomic_sequence("chr1", 0, 4).unwrap(), "ACGT");
    }

    #[test]
    fn test_get_protein_length_resolves_and_falls_back() {
        let mut provider = MockProvider::new();
        provider.add_protein("NP_TEST.1", "MAPLE");

        // Exact (versioned) match returns the stored length.
        assert_eq!(provider.get_protein_length("NP_TEST.1").unwrap(), 5);
        // Unversioned base accession falls back to the versioned entry.
        assert_eq!(provider.get_protein_length("NP_TEST").unwrap(), 5);
    }

    #[test]
    fn test_get_protein_length_unresolvable_is_zero() {
        // Contract: an accession the provider cannot resolve maps to a
        // length of `0` (matching the trait's default probe semantics),
        // NOT an error — otherwise `normalize()` is not behavior-preserving
        // when a provider switches from the probe loop to this API.
        let provider = MockProvider::with_test_data();
        assert_eq!(provider.get_protein_length("NP_NONEXISTENT.1").unwrap(), 0);
    }

    #[test]
    fn test_from_json_bare_array_form_still_works() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"[{
      "id": "NM_TEST.1",
      "gene_symbol": "TEST",
      "strand": "+",
      "sequence": "ATGCATGCAT",
      "cds_start": 1,
      "cds_end": 10,
      "exons": [{"number": 1, "start": 1, "end": 10}]
    }]"#;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider = MockProvider::from_json(file.path()).unwrap();

        assert!(provider.has_transcript("NM_TEST.1"));
        assert!(!provider.has_protein_data());
        assert!(!provider.has_genomic_data());
    }

    #[test]
    fn test_from_json_empty_object_form() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(b"{}").unwrap();

        let provider = MockProvider::from_json(file.path()).unwrap();

        assert!(provider.is_empty());
        assert!(!provider.has_protein_data());
        assert!(!provider.has_genomic_data());
    }

    #[test]
    fn test_from_json_rejects_scalar_root() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(b"42").unwrap();

        match MockProvider::from_json(file.path()) {
            Err(FerroError::Json { msg }) => {
                assert!(
                    msg.contains("array or object"),
                    "expected error message to mention 'array or object', got: {msg}",
                );
            }
            Err(e) => panic!("expected FerroError::Json, got {e}"),
            Ok(_) => panic!("expected scalar JSON root to be rejected"),
        }
    }

    #[test]
    fn test_from_json_rejects_unknown_field() {
        // A typo'd top-level key must error rather than silently parse as an
        // empty provider.
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(br#"{"transripts": []}"#).unwrap();

        match MockProvider::from_json(file.path()) {
            Err(e) => assert!(
                format!("{e}").contains("transripts"),
                "expected error to mention the unknown field, got {e}",
            ),
            Ok(_) => panic!("expected unknown field to be rejected"),
        }
    }

    #[test]
    fn from_json_accepts_convert_gff_output_with_metadata_keys() {
        // The container shape produced by `ferro convert-gff`.
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"{
            "version": "1.0",
            "genome_build": "GRCh38",
            "transcripts": []
        }"#;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider = MockProvider::from_json(file.path())
            .expect("convert-gff JSON with version/genome_build metadata should load");
        assert!(provider.is_empty());
    }

    #[test]
    fn test_mock_provider_returns_sequence_length_for_added_sequence() {
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_012920.1", "A".repeat(16569));
        assert_eq!(provider.get_sequence_length("NC_012920.1").unwrap(), 16569);
    }

    #[test]
    fn test_mock_provider_sequence_length_errors_for_unknown_id() {
        let provider = MockProvider::new();
        let err = provider.get_sequence_length("missing").unwrap_err();
        assert!(matches!(err, FerroError::ReferenceNotFound { .. }));
    }

    #[test]
    fn from_json_loads_convert_gff_output_with_transcript() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let json = r#"{
            "version": "1.0",
            "genome_build": "GRCh38",
            "transcripts": [
                {
                    "id": "NM_000001.1",
                    "strand": "+",
                    "sequence": "ATGC",
                    "exons": [{"number": 1, "start": 1, "end": 4}]
                }
            ]
        }"#;

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(json.as_bytes()).unwrap();

        let provider = MockProvider::from_json(file.path())
            .expect("convert-gff JSON with transcript should load");
        let tx = provider
            .get_transcript("NM_000001.1")
            .expect("tx not found");
        assert_eq!(tx.sequence.as_deref(), Some("ATGC"));
    }
}
