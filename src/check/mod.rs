//! Reference data validation for ferro-hgvs.
//!
//! This module provides functionality to verify that reference data
//! is properly configured and available for normalization.

use crate::prepare::ReferenceManifest;
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

    // Try to load manifest (paths are automatically made absolute)
    let manifest = match ReferenceManifest::load_or_default(reference_dir) {
        Ok(m) => m,
        Err(e) => return CheckResult::failure(format!("Failed to load manifest: {}", e)),
    };

    let mut result = CheckResult::success(manifest.clone());

    // Validate transcript files exist.
    //
    // `transcript_fastas` is populated by `prepare_references` with `.fna.gz`
    // RefSeq RNA paths; map those to the decompressed `.fna` companion to check.
    // Any other extension is checked as-is rather than silently rewritten.
    for fasta in &manifest.transcript_fastas {
        let check_path = fasta
            .file_name()
            .and_then(|n| n.to_str())
            .filter(|name| name.ends_with(".fna.gz"))
            .map(|name| fasta.with_file_name(&name[..name.len() - ".gz".len()]))
            .unwrap_or_else(|| fasta.clone());
        if !check_path.exists() {
            result.warnings.push(format!(
                "Transcript FASTA not found: {}",
                check_path.display()
            ));
        }
    }

    // Validate genome file if specified
    if let Some(ref genome) = manifest.genome_fasta {
        if !genome.exists() {
            result
                .warnings
                .push(format!("Genome FASTA not found: {}", genome.display()));
        }
    }

    // Validate cdot file if specified
    if let Some(ref cdot) = manifest.cdot_json {
        if !cdot.exists() {
            result
                .warnings
                .push(format!("cdot JSON not found: {}", cdot.display()));
        }
    }

    if let Some(ref cdot) = manifest.cdot_grch37_json {
        if !cdot.exists() {
            result
                .warnings
                .push(format!("cdot GRCh37 JSON not found: {}", cdot.display()));
        }
    }

    result
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
        eprintln!("  cdot metadata (GRCh38): {}", cdot.display());
    } else {
        eprintln!("  cdot metadata (GRCh38): not available");
    }

    if let Some(ref cdot) = manifest.cdot_grch37_json {
        eprintln!("  cdot metadata (GRCh37): {}", cdot.display());
    } else {
        eprintln!("  cdot metadata (GRCh37): not available");
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
    use std::fs::File;
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

        // Create a real transcript FASTA file with a relative path (as .fna, the expected format)
        let transcript_fasta = dir.path().join("example.fna");
        File::create(&transcript_fasta).unwrap();

        let mut manifest = ReferenceManifest {
            prepared_at: "2024-01-01T00:00:00Z".to_string(),
            transcript_fastas: vec![transcript_fasta.clone()],
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
            transcript_count: 1,
            available_prefixes: vec!["NM".to_string()],
            reference_dir: dir.path().to_path_buf(),
        };

        manifest.save().unwrap();

        let result = check_reference(dir.path());
        assert!(result.valid);
        assert!(result.manifest.is_some());
        assert!(
            result.warnings.is_empty(),
            "Expected no warnings for valid relative paths"
        );
    }
}
