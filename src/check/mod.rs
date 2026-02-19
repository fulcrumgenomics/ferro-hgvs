//! Reference data validation for ferro-hgvs.
//!
//! This module provides functionality to verify that reference data
//! is properly configured and available for normalization.

use crate::prepare::ReferenceManifest;
use crate::FerroError;
use std::fs::File;
use std::path::Path;

/// Result of checking reference data.
#[derive(Debug, Clone)]
pub struct CheckResult {
    /// Whether the reference data is valid
    pub valid: bool,
    /// The manifest if successfully loaded
    pub manifest: Option<ReferenceManifest>,
    /// Error messages if any
    pub errors: Vec<String>,
    /// Warning messages if any
    pub warnings: Vec<String>,
}

impl CheckResult {
    /// Create a successful check result.
    pub fn success(manifest: ReferenceManifest) -> Self {
        Self {
            valid: true,
            manifest: Some(manifest),
            errors: Vec::new(),
            warnings: Vec::new(),
        }
    }

    /// Create a failed check result.
    pub fn failure(error: String) -> Self {
        Self {
            valid: false,
            manifest: None,
            errors: vec![error],
            warnings: Vec::new(),
        }
    }

    /// Add a warning to the result.
    pub fn with_warning(mut self, warning: String) -> Self {
        self.warnings.push(warning);
        self
    }
}

/// Check reference data and return detailed result.
pub fn check_reference(reference_dir: &Path) -> CheckResult {
    let manifest_path = reference_dir.join("manifest.json");

    // Check manifest exists
    if !manifest_path.exists() {
        return CheckResult::failure(format!(
            "No reference data found at {}. Run 'ferro prepare' first.",
            reference_dir.display()
        ));
    }

    // Try to load manifest
    let manifest = match load_manifest(&manifest_path) {
        Ok(m) => m,
        Err(e) => return CheckResult::failure(format!("Failed to load manifest: {}", e)),
    };

    let mut result = CheckResult::success(manifest.clone());

    // Validate transcript files exist
    for fasta in &manifest.transcript_fastas {
        let full_path = reference_dir.join(fasta);
        let fna_path = full_path.with_extension("").with_extension("fna");
        if !fna_path.exists() {
            result.warnings.push(format!(
                "Transcript FASTA not found: {}",
                fna_path.display()
            ));
        }
    }

    // Validate genome file if specified
    if let Some(ref genome) = manifest.genome_fasta {
        let full_path = reference_dir.join(genome);
        if !full_path.exists() {
            result
                .warnings
                .push(format!("Genome FASTA not found: {}", full_path.display()));
        }
    }

    // Validate cdot file if specified
    if let Some(ref cdot) = manifest.cdot_json {
        let full_path = reference_dir.join(cdot);
        if !full_path.exists() {
            result
                .warnings
                .push(format!("cdot JSON not found: {}", full_path.display()));
        }
    }

    result
}

/// Load manifest from file.
fn load_manifest(manifest_path: &Path) -> Result<ReferenceManifest, FerroError> {
    let file = File::open(manifest_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open manifest: {}", e),
    })?;

    serde_json::from_reader(file).map_err(|e| FerroError::Io {
        msg: format!("Failed to parse manifest: {}", e),
    })
}

/// Print a detailed summary of reference data.
pub fn print_check_summary(result: &CheckResult, reference_dir: &Path) {
    if !result.valid {
        eprintln!("Reference check FAILED:");
        for error in &result.errors {
            eprintln!("  ERROR: {}", error);
        }
        return;
    }

    let manifest = result.manifest.as_ref().unwrap();

    eprintln!("=== Reference Data Check ===");
    eprintln!("  Directory: {}", reference_dir.display());
    eprintln!("  Status: OK");
    eprintln!();
    eprintln!("=== Contents ===");
    eprintln!("  Prepared at: {}", manifest.prepared_at);
    eprintln!("  Transcripts: {}", manifest.transcript_count);
    eprintln!("  Prefixes: {}", manifest.available_prefixes.join(", "));

    if let Some(ref genome) = manifest.genome_fasta {
        eprintln!("  GRCh38 genome: {}", genome.display());
    } else {
        eprintln!("  GRCh38 genome: not available");
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
        eprintln!("  cdot metadata: {}", cdot.display());
    } else {
        eprintln!("  cdot metadata: not available");
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

    // Print warnings
    if !result.warnings.is_empty() {
        eprintln!();
        eprintln!("=== Warnings ===");
        for warning in &result.warnings {
            eprintln!("  WARNING: {}", warning);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    #[test]
    fn test_check_missing_manifest() {
        let dir = TempDir::new().unwrap();
        let result = check_reference(dir.path());
        assert!(!result.valid);
        assert!(!result.errors.is_empty());
    }

    #[test]
    fn test_check_valid_manifest() {
        let dir = TempDir::new().unwrap();
        let manifest_path = dir.path().join("manifest.json");

        let manifest = ReferenceManifest {
            prepared_at: "2024-01-01T00:00:00Z".to_string(),
            transcript_fastas: Vec::new(),
            genome_fasta: None,
            genome_grch37_fasta: None,
            refseqgene_fastas: Vec::new(),
            lrg_fastas: Vec::new(),
            lrg_xmls: Vec::new(),
            lrg_refseq_mapping: None,
            cdot_json: None,
            supplemental_fasta: None,
            legacy_transcripts_fasta: None,
            legacy_transcripts_metadata: None,
            legacy_genbank_fasta: None,
            legacy_genbank_metadata: None,
            transcript_count: 0,
            available_prefixes: Vec::new(),
        };

        let mut file = File::create(&manifest_path).unwrap();
        write!(file, "{}", serde_json::to_string_pretty(&manifest).unwrap()).unwrap();

        let result = check_reference(dir.path());
        assert!(result.valid);
        assert!(result.manifest.is_some());
    }
}
