//! Reference manifest types and I/O operations.
//!
//! This module handles the definition, loading, and display of the `ReferenceManifest`
//! that tracks all prepared reference data files.

use crate::FerroError;
use std::fs::File;
use std::path::{Path, PathBuf};

/// Manifest of prepared reference data.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ReferenceManifest {
    /// When the data was prepared
    pub prepared_at: String,
    /// Transcript FASTA files
    pub transcript_fastas: Vec<PathBuf>,
    /// GRCh38 genome FASTA file (if downloaded)
    pub genome_fasta: Option<PathBuf>,
    /// GRCh37 genome FASTA file (if downloaded)
    #[serde(default)]
    pub genome_grch37_fasta: Option<PathBuf>,
    /// RefSeqGene FASTA files (NG_* accessions)
    #[serde(default)]
    pub refseqgene_fastas: Vec<PathBuf>,
    /// LRG FASTA files (LRG_* accessions)
    #[serde(default)]
    pub lrg_fastas: Vec<PathBuf>,
    /// LRG XML files with full annotation structure
    #[serde(default)]
    pub lrg_xmls: Vec<PathBuf>,
    /// LRG to RefSeq transcript mapping file
    #[serde(default)]
    pub lrg_refseq_mapping: Option<PathBuf>,
    /// cdot transcript metadata JSON for GRCh38 (if downloaded)
    pub cdot_json: Option<PathBuf>,
    /// cdot transcript metadata JSON for GRCh37 (if downloaded)
    #[serde(default)]
    pub cdot_grch37_json: Option<PathBuf>,
    /// Supplemental FASTA file (missing ClinVar transcripts fetched from NCBI)
    #[serde(default)]
    pub supplemental_fasta: Option<PathBuf>,
    /// Legacy transcript versions FASTA (older versions not in current RefSeq)
    #[serde(default)]
    pub legacy_transcripts_fasta: Option<PathBuf>,
    /// Legacy transcript metadata JSON (CDS coordinates, gene names)
    #[serde(default)]
    pub legacy_transcripts_metadata: Option<PathBuf>,
    /// Legacy GenBank sequences FASTA (non-RefSeq sequences like U31929.1)
    #[serde(default)]
    pub legacy_genbank_fasta: Option<PathBuf>,
    /// Legacy GenBank metadata JSON (CDS coordinates, gene names)
    #[serde(default)]
    pub legacy_genbank_metadata: Option<PathBuf>,
    /// Total number of transcripts
    pub transcript_count: usize,
    /// List of available accession prefixes
    pub available_prefixes: Vec<String>,
    /// Directory containing this manifest (runtime property, not serialized)
    #[serde(skip)]
    pub reference_dir: PathBuf,
}

impl Default for ReferenceManifest {
    fn default() -> Self {
        Self {
            prepared_at: chrono::Utc::now().to_rfc3339(),
            transcript_fastas: Vec::new(),
            genome_fasta: None,
            genome_grch37_fasta: None,
            refseqgene_fastas: Vec::new(),
            lrg_fastas: Vec::new(),
            lrg_xmls: Vec::new(),
            lrg_refseq_mapping: None,
            cdot_json: None,
            cdot_grch37_json: None,
            supplemental_fasta: None,
            legacy_transcripts_fasta: None,
            legacy_transcripts_metadata: None,
            legacy_genbank_fasta: None,
            legacy_genbank_metadata: None,
            transcript_count: 0,
            available_prefixes: Vec::new(),
            reference_dir: PathBuf::new(),
        }
    }
}

impl ReferenceManifest {
    /// Load manifest from directory, or create a fresh one if it doesn't exist.
    pub fn load_or_default(reference_dir: &Path) -> Result<Self, FerroError> {
        let manifest_path = reference_dir.join("manifest.json");

        let mut manifest = if manifest_path.exists() {
            let file = File::open(&manifest_path).map_err(|e| FerroError::Io {
                msg: format!("Failed to open manifest: {}", e),
            })?;

            serde_json::from_reader(file).map_err(|e| FerroError::Io {
                msg: format!("Failed to parse manifest: {}", e),
            })?
        } else {
            Self::default()
        };

        manifest.reference_dir = reference_dir.to_path_buf();
        manifest.make_paths_absolute();
        Ok(manifest)
    }

    /// Save manifest to its reference directory.
    ///
    /// Automatically deduplicates paths and converts them to relative (for portability)
    /// before serializing to JSON.
    pub fn save(&self) -> Result<(), FerroError> {
        let mut manifest = self.clone();
        manifest.deduplicate_paths();
        manifest.make_paths_relative();

        let manifest_path = self.reference_dir.join("manifest.json");
        let file = File::create(&manifest_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create manifest: {}", e),
        })?;
        serde_json::to_writer_pretty(file, &manifest).map_err(|e| FerroError::Io {
            msg: format!("Failed to write manifest: {}", e),
        })
    }

    /// Apply closures to all Vec<PathBuf> and Option<PathBuf> fields.
    fn for_each_path_mut(
        &mut self,
        mut vec_fn: impl FnMut(&mut Vec<PathBuf>),
        mut opt_fn: impl FnMut(&mut Option<PathBuf>),
    ) {
        // Vec<PathBuf> fields
        vec_fn(&mut self.transcript_fastas);
        vec_fn(&mut self.refseqgene_fastas);
        vec_fn(&mut self.lrg_fastas);
        vec_fn(&mut self.lrg_xmls);

        // Option<PathBuf> fields
        opt_fn(&mut self.genome_fasta);
        opt_fn(&mut self.genome_grch37_fasta);
        opt_fn(&mut self.lrg_refseq_mapping);
        opt_fn(&mut self.cdot_json);
        opt_fn(&mut self.cdot_grch37_json);
        opt_fn(&mut self.supplemental_fasta);
        opt_fn(&mut self.legacy_transcripts_fasta);
        opt_fn(&mut self.legacy_transcripts_metadata);
        opt_fn(&mut self.legacy_genbank_fasta);
        opt_fn(&mut self.legacy_genbank_metadata);
    }

    /// Convert all paths in the manifest to be relative to the reference directory.
    ///
    /// This ensures the manifest is portable - paths work when running from the
    /// directory containing the manifest, regardless of where `prepare` was run from.
    pub fn make_paths_relative(&mut self) {
        let base = self.reference_dir.clone();
        self.for_each_path_mut(
            |vec| {
                for p in vec {
                    if let Ok(stripped) = p.strip_prefix(&base) {
                        *p = stripped.to_path_buf();
                    }
                }
            },
            |opt| {
                if let Some(p) = opt {
                    if let Ok(stripped) = p.strip_prefix(&base) {
                        *p = stripped.to_path_buf();
                    }
                }
            },
        );
    }

    /// Convert all relative paths to absolute, resolved against the manifest's reference directory.
    ///
    /// Called when loading a manifest to ensure all paths are absolute for use in the program.
    pub fn make_paths_absolute(&mut self) {
        let base = self.reference_dir.clone();
        self.for_each_path_mut(
            |vec| {
                for p in vec {
                    if !p.is_absolute() {
                        *p = base.join(p.as_path());
                    }
                }
            },
            |opt| {
                if let Some(p) = opt {
                    if !p.is_absolute() {
                        *p = base.join(p.as_path());
                    }
                }
            },
        );
    }

    /// Deduplicate paths in all path lists.
    pub fn deduplicate_paths(&mut self) {
        self.for_each_path_mut(
            |vec| {
                vec.sort();
                vec.dedup();
            },
            |_opt| {
                // no-op: Option<PathBuf> is a single value, so dedup is irrelevant
            },
        );
    }
}

/// Check what reference data is available.
pub fn check_references(reference_dir: &Path) -> Result<ReferenceManifest, FerroError> {
    let manifest_path = reference_dir.join("manifest.json");

    if !manifest_path.exists() {
        return Err(FerroError::Io {
            msg: format!(
                "No reference data found at {}. Run 'ferro prepare' first.",
                reference_dir.display()
            ),
        });
    }

    ReferenceManifest::load_or_default(reference_dir)
}

/// Print a summary of reference data.
pub fn print_reference_summary(manifest: &ReferenceManifest) {
    eprintln!("=== Reference Data Summary ===");
    eprintln!("  Directory: {}", manifest.reference_dir.display());
    eprintln!("  Prepared at: {}", manifest.prepared_at);
    eprintln!("  Transcripts: {}", manifest.transcript_count);
    eprintln!(
        "  Available prefixes: {}",
        manifest.available_prefixes.join(", ")
    );

    if let Some(ref genome) = manifest.genome_fasta {
        eprintln!("  GRCh38 genome: {}", genome.display());
    }
    if let Some(ref genome) = manifest.genome_grch37_fasta {
        eprintln!("  GRCh37 genome: {}", genome.display());
    }
    if !manifest.refseqgene_fastas.is_empty() {
        eprintln!("  RefSeqGene files: {}", manifest.refseqgene_fastas.len());
    }
    if !manifest.lrg_fastas.is_empty() {
        eprintln!("  LRG files: {}", manifest.lrg_fastas.len());
    }
    if let Some(ref cdot) = manifest.cdot_json {
        eprintln!("  cdot metadata (GRCh38): {}", cdot.display());
    }
    if let Some(ref cdot) = manifest.cdot_grch37_json {
        eprintln!("  cdot metadata (GRCh37): {}", cdot.display());
    }
    if let Some(ref supp) = manifest.supplemental_fasta {
        eprintln!("  Supplemental transcripts: {}", supp.display());
    }
    if let Some(ref legacy) = manifest.legacy_transcripts_fasta {
        eprintln!("  Legacy transcripts: {}", legacy.display());
    }
    if let Some(ref genbank) = manifest.legacy_genbank_fasta {
        eprintln!("  Legacy GenBank: {}", genbank.display());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_paths_relative() {
        let mut manifest = ReferenceManifest::default();
        manifest.reference_dir = PathBuf::from("/ref/data");
        manifest.transcript_fastas = vec![PathBuf::from("/ref/data/transcripts.fa")];
        manifest.cdot_json = Some(PathBuf::from("/ref/data/cdot.json"));

        manifest.make_paths_relative();

        assert_eq!(
            manifest.transcript_fastas[0],
            PathBuf::from("transcripts.fa")
        );
        assert_eq!(manifest.cdot_json, Some(PathBuf::from("cdot.json")));
    }

    #[test]
    fn test_make_paths_absolute() {
        let mut manifest = ReferenceManifest::default();
        manifest.reference_dir = PathBuf::from("/ref/data");
        manifest.transcript_fastas = vec![PathBuf::from("transcripts.fa")];
        manifest.cdot_json = Some(PathBuf::from("cdot.json"));

        manifest.make_paths_absolute();

        assert_eq!(
            manifest.transcript_fastas[0],
            PathBuf::from("/ref/data/transcripts.fa")
        );
        assert_eq!(
            manifest.cdot_json,
            Some(PathBuf::from("/ref/data/cdot.json"))
        );
    }

    #[test]
    fn test_deduplicate_paths() {
        let mut manifest = ReferenceManifest::default();
        manifest.transcript_fastas = vec![
            PathBuf::from("b.fa"),
            PathBuf::from("a.fa"),
            PathBuf::from("b.fa"),
        ];
        manifest.lrg_fastas = vec![PathBuf::from("lrg.fa"), PathBuf::from("lrg.fa")];

        manifest.deduplicate_paths();

        assert_eq!(
            manifest.transcript_fastas,
            vec![PathBuf::from("a.fa"), PathBuf::from("b.fa")]
        );
        assert_eq!(manifest.lrg_fastas, vec![PathBuf::from("lrg.fa")]);
    }

    #[test]
    fn test_roundtrip_save_load_with_relative_paths() {
        use std::io::Read;
        use tempfile::TempDir;

        let dir = TempDir::new().unwrap();
        let ref_dir = dir.path();

        // Create manifest with absolute paths
        let manifest = ReferenceManifest {
            prepared_at: "2024-01-01T00:00:00Z".to_string(),
            transcript_fastas: vec![ref_dir.join("transcripts.fa")],
            genome_fasta: Some(ref_dir.join("genome.fa")),
            genome_grch37_fasta: Some(ref_dir.join("genome37.fa")),
            refseqgene_fastas: vec![ref_dir.join("ng.fa")],
            lrg_fastas: vec![ref_dir.join("lrg.fa")],
            lrg_xmls: vec![ref_dir.join("lrg.xml")],
            lrg_refseq_mapping: Some(ref_dir.join("lrg_mapping.txt")),
            cdot_json: Some(ref_dir.join("cdot.json")),
            cdot_grch37_json: Some(ref_dir.join("cdot37.json")),
            supplemental_fasta: Some(ref_dir.join("supplemental.fa")),
            legacy_transcripts_fasta: Some(ref_dir.join("legacy.fa")),
            legacy_transcripts_metadata: Some(ref_dir.join("legacy.json")),
            legacy_genbank_fasta: Some(ref_dir.join("genbank.fa")),
            legacy_genbank_metadata: Some(ref_dir.join("genbank.json")),
            transcript_count: 100,
            available_prefixes: vec!["NM".to_string()],
            reference_dir: ref_dir.to_path_buf(),
        };

        // Save the manifest (which should make paths relative)
        manifest.save().unwrap();

        // Verify on disk: paths are relative and reference_dir is not serialized
        let manifest_file = ref_dir.join("manifest.json");
        let mut contents = String::new();
        File::open(&manifest_file)
            .unwrap()
            .read_to_string(&mut contents)
            .unwrap();

        let json: serde_json::Value = serde_json::from_str(&contents).unwrap();

        // Check that reference_dir is not in the serialized JSON
        assert!(
            json.get("reference_dir").is_none(),
            "reference_dir should not be serialized"
        );

        // Check that paths are relative (not absolute)
        assert_eq!(
            json["transcript_fastas"][0],
            "transcripts.fa",
            "transcript_fastas should be stored as relative path"
        );
        assert_eq!(
            json["genome_fasta"],
            "genome.fa",
            "genome_fasta should be stored as relative path"
        );
        assert_eq!(
            json["cdot_json"],
            "cdot.json",
            "cdot_json should be stored as relative path"
        );

        // Load the manifest back
        let loaded = ReferenceManifest::load_or_default(ref_dir).unwrap();

        // Verify loaded manifest has reference_dir set
        assert_eq!(loaded.reference_dir, ref_dir, "reference_dir should be set after load");

        // Verify paths were converted back to absolute
        assert_eq!(
            loaded.transcript_fastas[0],
            ref_dir.join("transcripts.fa"),
            "transcript_fastas should be absolute after load"
        );
        assert_eq!(
            loaded.genome_fasta,
            Some(ref_dir.join("genome.fa")),
            "genome_fasta should be absolute after load"
        );
        assert_eq!(
            loaded.cdot_json,
            Some(ref_dir.join("cdot.json")),
            "cdot_json should be absolute after load"
        );

        // Verify all other fields are preserved
        assert_eq!(loaded.transcript_count, 100);
        assert_eq!(loaded.available_prefixes, vec!["NM"]);
    }
}
