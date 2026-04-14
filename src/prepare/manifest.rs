//! Reference manifest types and I/O operations.
//!
//! This module handles the definition, loading, and display of the `ReferenceManifest`
//! that tracks all prepared reference data files.

use crate::FerroError;
use std::fs::File;
use std::path::{Component, Path, PathBuf};

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
            // Empty string by default; `save()` populates this with the current time.
            prepared_at: String::new(),
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

    /// Validate the reference-root invariant before saving.
    ///
    /// Ensures that:
    /// 1. `reference_dir` is set (not empty PathBuf from default)
    /// 2. All tracked paths can be made relative to `reference_dir`
    /// 3. No absolute or out-of-root paths are written to manifest
    /// 4. No relative paths contain `..` components that could escape `reference_dir`
    ///
    /// Returns an Io error with a clear message if validation fails.
    fn validate_reference_root_invariant(&self) -> Result<(), FerroError> {
        // Check that reference_dir is set (not empty)
        if self.reference_dir.as_os_str().is_empty() {
            return Err(FerroError::Io {
                msg: "Invariant violation: reference_dir must be set before saving manifest. \
                       Manifest may have been created with default() and not properly initialized. \
                       Call load_or_default(reference_dir) instead."
                    .to_string(),
            });
        }

        let base = self.reference_dir.as_path();
        let mut out_of_root_paths = Vec::new();
        self.for_each_path(|p| {
            // Reject `..` components anywhere; otherwise after `make_paths_relative`
            // they would be persisted unchanged and resolve outside `reference_dir`
            // when loaded.
            if p.components().any(|c| matches!(c, Component::ParentDir)) {
                out_of_root_paths.push(format!(
                    "path '{}' contains a '..' component that could escape reference_dir",
                    p.display()
                ));
            } else if p.is_absolute() && !p.starts_with(base) {
                out_of_root_paths.push(format!(
                    "path '{}' is outside reference_dir '{}'",
                    p.display(),
                    base.display()
                ));
            }
        });

        if !out_of_root_paths.is_empty() {
            return Err(FerroError::Io {
                msg: format!(
                    "Invariant violation: {} path(s) are outside reference_dir. \
                     Manifest can only contain paths relative to reference_dir or within it:\n  {}",
                    out_of_root_paths.len(),
                    out_of_root_paths.join("\n  ")
                ),
            });
        }

        Ok(())
    }

    /// Save manifest to its reference directory.
    ///
    /// Automatically deduplicates paths in-place and refreshes `prepared_at` so that
    /// re-runs of `prepare` reflect when the manifest was last persisted. The on-disk
    /// JSON stores paths relative to `reference_dir` for portability, while the
    /// in-memory manifest retains its absolute paths.
    ///
    /// Validates the reference-root invariant to ensure all paths are within or can be
    /// made relative to the reference directory. Returns an Io error if validation fails.
    pub fn save(&mut self) -> Result<(), FerroError> {
        // Validate invariant before any modifications
        self.validate_reference_root_invariant()?;

        self.prepared_at = chrono::Utc::now().to_rfc3339();
        self.deduplicate_paths();

        // Serialize a relative-path view without mutating the in-memory absolute paths.
        let mut on_disk = self.clone();
        on_disk.make_paths_relative();

        let manifest_path = self.reference_dir.join("manifest.json");
        let file = File::create(&manifest_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create manifest: {}", e),
        })?;
        serde_json::to_writer_pretty(file, &on_disk).map_err(|e| FerroError::Io {
            msg: format!("Failed to write manifest: {}", e),
        })
    }

    /// Apply a closure to every tracked path, read-only.
    fn for_each_path(&self, mut f: impl FnMut(&Path)) {
        for p in self
            .transcript_fastas
            .iter()
            .chain(self.refseqgene_fastas.iter())
            .chain(self.lrg_fastas.iter())
            .chain(self.lrg_xmls.iter())
        {
            f(p);
        }
        for p in [
            &self.genome_fasta,
            &self.genome_grch37_fasta,
            &self.lrg_refseq_mapping,
            &self.cdot_json,
            &self.cdot_grch37_json,
            &self.supplemental_fasta,
            &self.legacy_transcripts_fasta,
            &self.legacy_transcripts_metadata,
            &self.legacy_genbank_fasta,
            &self.legacy_genbank_metadata,
        ]
        .into_iter()
        .flatten()
        {
            f(p);
        }
    }

    /// Apply a closure to every tracked path, mutably.
    fn for_each_path_mut(&mut self, mut f: impl FnMut(&mut PathBuf)) {
        for v in [
            &mut self.transcript_fastas,
            &mut self.refseqgene_fastas,
            &mut self.lrg_fastas,
            &mut self.lrg_xmls,
        ] {
            for p in v.iter_mut() {
                f(p);
            }
        }
        for o in [
            &mut self.genome_fasta,
            &mut self.genome_grch37_fasta,
            &mut self.lrg_refseq_mapping,
            &mut self.cdot_json,
            &mut self.cdot_grch37_json,
            &mut self.supplemental_fasta,
            &mut self.legacy_transcripts_fasta,
            &mut self.legacy_transcripts_metadata,
            &mut self.legacy_genbank_fasta,
            &mut self.legacy_genbank_metadata,
        ] {
            if let Some(p) = o.as_mut() {
                f(p);
            }
        }
    }

    /// Convert all paths in the manifest to be relative to the reference directory.
    ///
    /// This ensures the manifest is portable - paths work when running from the
    /// directory containing the manifest, regardless of where `prepare` was run from.
    /// Paths that cannot be stripped (e.g., outside `reference_dir`) are left
    /// unchanged; callers must run `validate_reference_root_invariant` first to
    /// guarantee a fully relative result.
    fn make_paths_relative(&mut self) {
        let base = self.reference_dir.clone();
        self.for_each_path_mut(|p| {
            if let Ok(stripped) = p.strip_prefix(&base) {
                *p = stripped.to_path_buf();
            }
        });
    }

    /// Convert all relative paths to absolute, resolved against the manifest's reference directory.
    ///
    /// Called when loading a manifest to ensure all paths are absolute for use in the program.
    /// Paths that are already absolute are left unchanged.
    fn make_paths_absolute(&mut self) {
        let base = self.reference_dir.clone();
        self.for_each_path_mut(|p| {
            if !p.is_absolute() {
                *p = base.join(p.as_path());
            }
        });
    }

    /// Deduplicate paths in all path lists.
    fn deduplicate_paths(&mut self) {
        for v in [
            &mut self.transcript_fastas,
            &mut self.refseqgene_fastas,
            &mut self.lrg_fastas,
            &mut self.lrg_xmls,
        ] {
            v.sort();
            v.dedup();
        }
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
        let mut manifest = ReferenceManifest {
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
            json["transcript_fastas"][0], "transcripts.fa",
            "transcript_fastas should be stored as relative path"
        );
        assert_eq!(
            json["genome_fasta"], "genome.fa",
            "genome_fasta should be stored as relative path"
        );
        assert_eq!(
            json["cdot_json"], "cdot.json",
            "cdot_json should be stored as relative path"
        );

        // Load the manifest back
        let loaded = ReferenceManifest::load_or_default(ref_dir).unwrap();

        // Verify loaded manifest has reference_dir set
        assert_eq!(
            loaded.reference_dir, ref_dir,
            "reference_dir should be set after load"
        );

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

    #[test]
    fn test_validate_reference_root_invariant_rejects_out_of_root_paths() {
        use tempfile::TempDir;

        let dir = TempDir::new().unwrap();
        let other = TempDir::new().unwrap();

        // In-root paths (absolute, under reference_dir) should pass.
        let mut manifest = ReferenceManifest {
            reference_dir: dir.path().to_path_buf(),
            transcript_fastas: vec![dir.path().join("transcripts.fa")],
            cdot_json: Some(dir.path().join("cdot.json")),
            ..ReferenceManifest::default()
        };
        manifest
            .validate_reference_root_invariant()
            .expect("in-root absolute paths should validate");

        // An absolute path outside reference_dir in a Vec field should fail.
        manifest
            .transcript_fastas
            .push(other.path().join("rogue.fa"));
        let err = manifest
            .validate_reference_root_invariant()
            .expect_err("out-of-root vec path must be rejected");
        assert!(format!("{}", err).contains("rogue.fa"));

        // An absolute path outside reference_dir in an Option field should also fail.
        manifest.transcript_fastas.pop();
        manifest.legacy_transcripts_fasta = Some(other.path().join("legacy.fa"));
        let err = manifest
            .validate_reference_root_invariant()
            .expect_err("out-of-root option path must be rejected");
        assert!(format!("{}", err).contains("legacy.fa"));

        // Empty reference_dir should fail with the dedicated message.
        let bad = ReferenceManifest::default();
        let err = bad
            .validate_reference_root_invariant()
            .expect_err("default manifest must fail validation");
        assert!(format!("{}", err).contains("reference_dir must be set"));
    }

    #[test]
    fn test_validate_reference_root_invariant_rejects_parent_dir_components() {
        use tempfile::TempDir;

        let dir = TempDir::new().unwrap();

        // A relative path with `..` would resolve outside reference_dir on load.
        let mut manifest = ReferenceManifest {
            reference_dir: dir.path().to_path_buf(),
            transcript_fastas: vec![PathBuf::from("../escape.fa")],
            ..ReferenceManifest::default()
        };
        let err = manifest
            .validate_reference_root_invariant()
            .expect_err("relative '..' path must be rejected");
        assert!(format!("{}", err).contains("escape.fa"));
        assert!(format!("{}", err).contains(".."));

        // Same check on Option fields.
        manifest.transcript_fastas.clear();
        manifest.cdot_json = Some(PathBuf::from("subdir/../../oops.json"));
        let err = manifest
            .validate_reference_root_invariant()
            .expect_err("Option path with '..' must be rejected");
        assert!(format!("{}", err).contains("oops.json"));
    }

    #[test]
    fn test_save_refreshes_prepared_at() {
        use std::io::Read;
        use tempfile::TempDir;

        let dir = TempDir::new().unwrap();
        let ref_dir = dir.path();

        let stale_timestamp = "2020-01-01T00:00:00+00:00".to_string();
        let mut manifest = ReferenceManifest {
            prepared_at: stale_timestamp.clone(),
            reference_dir: ref_dir.to_path_buf(),
            ..ReferenceManifest::default()
        };

        manifest.save().unwrap();

        let mut contents = String::new();
        File::open(ref_dir.join("manifest.json"))
            .unwrap()
            .read_to_string(&mut contents)
            .unwrap();
        let json: serde_json::Value = serde_json::from_str(&contents).unwrap();
        let saved_prepared_at = json["prepared_at"].as_str().unwrap();

        assert_ne!(
            saved_prepared_at, stale_timestamp,
            "save() should refresh prepared_at to a new timestamp on each call"
        );
        // Sanity check: the new timestamp should parse as RFC 3339
        chrono::DateTime::parse_from_rfc3339(saved_prepared_at)
            .expect("prepared_at should be a valid RFC 3339 timestamp");

        // save() should also update the in-memory `prepared_at` so it stays in
        // sync with the persisted value.
        assert_eq!(
            manifest.prepared_at, saved_prepared_at,
            "save() should refresh in-memory prepared_at to match what was persisted"
        );
    }
}
