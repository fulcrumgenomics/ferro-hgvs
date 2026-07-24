//! Reference data validation for ferro-hgvs.
//!
//! This module provides functionality to verify that reference data
//! is properly configured and available for normalization.

use crate::error::FerroError;
use crate::prepare::ReferenceManifest;
use crate::reference::mock::MockProvider;
use crate::reference::provider::ReferenceProvider;
use std::path::{Path, PathBuf};

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

/// Summary of a standalone `transcripts.json` reference (the
/// `convert-gff`/`build-transcript` → `MockProvider` path), as opposed to a
/// prepared manifest directory.
#[derive(Debug, Clone)]
pub struct TranscriptsJsonSummary {
    /// Path to the checked JSON file.
    pub path: PathBuf,
    /// Number of transcripts loaded.
    pub transcript_count: usize,
    /// Whether the reference carries **usable** genomic sequence (at least one
    /// non-empty contig), so the genome-dependent normalization rules can run — see
    /// the capability boundary in `docs/transcripts_json_schema.md`. An empty
    /// contig string does not count.
    pub genome_capable: bool,
    /// Whether the reference carries protein sequences.
    pub has_protein_data: bool,
}

/// Check a standalone `transcripts.json` file (rather than a prepared reference
/// directory). Loads it through [`MockProvider::from_json`], so it enforces the
/// same schema and version validation the runtime load does, and reports what
/// capability the reference has (#1012 comment 2).
///
/// Fails if the reference resolves nothing (no transcripts, no genomic bases, no
/// proteins) rather than green-lighting a dead reference.
pub fn check_transcripts_json(path: &Path) -> Result<TranscriptsJsonSummary, FerroError> {
    let provider = MockProvider::from_json(path)?;
    let transcript_count = provider.len();
    let genome_capable = provider.total_genomic_bases() > 0;
    let has_protein_data = provider.has_protein_data();

    if transcript_count == 0 && !genome_capable && !has_protein_data {
        return Err(FerroError::Json {
            msg: "reference has no usable data: no transcripts, genomic sequence, or proteins"
                .to_string(),
        });
    }

    Ok(TranscriptsJsonSummary {
        path: path.to_path_buf(),
        transcript_count,
        genome_capable,
        has_protein_data,
    })
}

/// Print a summary of a successful standalone `transcripts.json` check.
pub fn print_transcripts_json_summary(summary: &TranscriptsJsonSummary) {
    eprintln!("=== transcripts.json Check ===");
    eprintln!("  File: {}", summary.path.display());
    eprintln!("  Status: OK");
    eprintln!();
    eprintln!("=== Contents ===");
    eprintln!("  Transcripts: {}", summary.transcript_count);
    eprintln!(
        "  Genome-capable: {}",
        if summary.genome_capable {
            "yes"
        } else {
            "no (transcript-level normalization only)"
        }
    );
    eprintln!(
        "  Protein data: {}",
        if summary.has_protein_data {
            "yes"
        } else {
            "no"
        }
    );
}

/// Print a clean failure block for a standalone `transcripts.json` check, mirroring
/// the prepared-directory path's presentation instead of a raw error dump.
pub fn print_transcripts_json_failure(path: &Path, error: &FerroError) {
    eprintln!("=== transcripts.json Check ===");
    eprintln!("  File: {}", path.display());
    eprintln!("  Status: FAILED");
    eprintln!("  ERROR: {}", error);
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

    // Validate Ensembl cDNA files exist.
    //
    // `ensembl_transcript_fastas` holds the downloaded Ensembl cDNA `.fa.gz`
    // paths; `prepare_references` decompresses and `index_fasta`s the companion
    // `.fa`, which is what `MultiFastaProvider` actually resolves from. Mirror
    // the RefSeq loop above: map a trailing `.fa.gz`/`.fna.gz` to its
    // decompressed sidecar and warn if that file is missing, so a missing or
    // un-decompressed Ensembl FASTA surfaces here rather than as a silent
    // runtime resolution miss. Any other extension is checked as-is.
    for fasta in &manifest.ensembl_transcript_fastas {
        let check_path = fasta
            .file_name()
            .and_then(|n| n.to_str())
            .filter(|name| name.ends_with(".fa.gz") || name.ends_with(".fna.gz"))
            .map(|name| fasta.with_file_name(&name[..name.len() - ".gz".len()]))
            .unwrap_or_else(|| fasta.clone());
        if !check_path.exists() {
            result.warnings.push(format!(
                "Ensembl cDNA FASTA not found: {}",
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

    if let Some(ref cdot) = manifest.ensembl_cdot_json {
        if !cdot.exists() {
            result
                .warnings
                .push(format!("Ensembl cdot JSON not found: {}", cdot.display()));
        }
    }

    if let Some(ref cdot) = manifest.ensembl_cdot_grch37_json {
        if !cdot.exists() {
            result.warnings.push(format!(
                "Ensembl cdot GRCh37 JSON not found: {}",
                cdot.display()
            ));
        }
    }

    if let Some(ref overrides) = manifest.canonical_overrides {
        if !overrides.exists() {
            result.warnings.push(format!(
                "Canonical overrides JSON not found: {}",
                overrides.display()
            ));
        }
    }

    // A missing RefSeqGene alignments GFF3 silently disables NG_/LRG_
    // re-anchoring at load, so surface it like the other optional inputs above.
    if let Some(ref alignments) = manifest.refseqgene_alignments {
        if !alignments.exists() {
            result.warnings.push(format!(
                "RefSeqGene alignments GFF3 not found: {}",
                alignments.display()
            ));
        }
    }
    if let Some(ref alignments) = manifest.refseqgene_alignments_grch37 {
        if !alignments.exists() {
            result.warnings.push(format!(
                "GRCh37 RefSeqGene alignments GFF3 not found: {}",
                alignments.display()
            ));
        }
    }

    // A missing assembly report silently disables the data-driven build map
    // at load (#716), so surface it like the other optional inputs above.
    if let Some(ref report) = manifest.assembly_report {
        if !report.exists() {
            result
                .warnings
                .push(format!("Assembly report not found: {}", report.display()));
        }
    }
    if let Some(ref report) = manifest.assembly_report_grch37 {
        if !report.exists() {
            result.warnings.push(format!(
                "GRCh37 assembly report not found: {}",
                report.display()
            ));
        }
    }

    // A missing derived-placements JSON silently disables version-gap NG_/LRG_
    // projection at load, so surface it like the other optional inputs above.
    if let Some(ref placements) = manifest.derived_refseqgene_placements {
        if !placements.exists() {
            result.warnings.push(format!(
                "Derived RefSeqGene placements not found: {}",
                placements.display()
            ));
        }
    }

    // A missing RefSeqGene summary (LRG_RefSeqGene) silently disables legacy
    // gene-model selector resolution at load, so surface it like the other
    // optional inputs above.
    if let Some(ref summary) = manifest.refseqgene_summary {
        if !summary.exists() {
            result.warnings.push(format!(
                "RefSeqGene summary (LRG_RefSeqGene) not found: {}",
                summary.display()
            ));
        }
    }

    // A missing NG_ hosted-transcripts map silently drops parent-relative legacy
    // selector resolution at load (#792) — every NG_-parent legacy selector then
    // falls back to the global reference-standard map — so surface it like the
    // other optional inputs above.
    if let Some(ref ng_hosted) = manifest.ng_hosted_transcripts {
        if !ng_hosted.exists() {
            result.warnings.push(format!(
                "NG_ hosted-transcripts map not found: {}",
                ng_hosted.display()
            ));
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

    if !manifest.ensembl_transcript_fastas.is_empty() {
        eprintln!(
            "  Ensembl cDNA files: {}",
            manifest.ensembl_transcript_fastas.len()
        );
    }

    if let Some(ref cdot) = manifest.ensembl_cdot_json {
        eprintln!("  Ensembl cdot metadata (GRCh38): {}", cdot.display());
    }

    if let Some(ref cdot) = manifest.ensembl_cdot_grch37_json {
        eprintln!("  Ensembl cdot metadata (GRCh37): {}", cdot.display());
    }

    if !manifest.protein_fastas.is_empty() {
        eprintln!("  Protein FASTA files: {}", manifest.protein_fastas.len());
    }

    if let Some(ref overrides) = manifest.canonical_overrides {
        eprintln!("  Canonical overrides: {}", overrides.display());
    }

    if let Some(ref supp) = manifest.supplemental_fasta {
        eprintln!("  Supplemental transcripts: {}", supp.display());
    }

    if let Some(ref backfill) = manifest.backfill_transcripts_fasta {
        eprintln!("  Backfilled transcripts: {}", backfill.display());
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
    fn check_transcripts_json_reports_transcript_only_reference() {
        use std::io::Write;
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("transcripts.json");
        let json = r#"{
            "version": "1.0",
            "genome_build": "GRCh38",
            "transcripts": [
                {"id": "NM_1.1", "strand": "+", "sequence": "ACGTACGT",
                 "exons": [{"number": 1, "start": 1, "end": 8}]}
            ]
        }"#;
        std::fs::File::create(&path)
            .unwrap()
            .write_all(json.as_bytes())
            .unwrap();

        let summary = check_transcripts_json(&path).expect("loads");
        assert_eq!(summary.transcript_count, 1);
        assert!(
            !summary.genome_capable,
            "no genomic_sequences → not genome-capable"
        );
        assert!(!summary.has_protein_data);
    }

    #[test]
    fn check_transcripts_json_reports_genome_capable_reference() {
        use std::io::Write;
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("transcripts.json");
        let json = r#"{
            "version": "1.0",
            "transcripts": [
                {"id": "NM_1.1", "strand": "+", "sequence": "ACGT",
                 "chromosome": "chr1", "genomic_start": 1, "genomic_end": 4,
                 "exons": [{"number": 1, "start": 1, "end": 4,
                            "genomic_start": 1, "genomic_end": 4}]}
            ],
            "genomic_sequences": {"chr1": "ACGTACGT"}
        }"#;
        std::fs::File::create(&path)
            .unwrap()
            .write_all(json.as_bytes())
            .unwrap();

        let summary = check_transcripts_json(&path).expect("loads");
        assert_eq!(summary.transcript_count, 1);
        assert!(
            summary.genome_capable,
            "genomic_sequences present → genome-capable"
        );
    }

    #[test]
    fn check_transcripts_json_rejects_empty_reference() {
        use std::io::Write;
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("transcripts.json");
        std::fs::File::create(&path)
            .unwrap()
            .write_all(br#"{"transcripts": []}"#)
            .unwrap();
        assert!(
            check_transcripts_json(&path).is_err(),
            "a reference with no usable data must fail the check, not be green-lit"
        );
    }

    #[test]
    fn check_transcripts_json_empty_contig_is_not_genome_capable() {
        use std::io::Write;
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("transcripts.json");
        // An empty contig string keeps a transcript present but carries no bases.
        let json = r#"{
            "transcripts": [
                {"id": "NM_1.1", "strand": "+", "sequence": "ACGT",
                 "exons": [{"number": 1, "start": 1, "end": 4}]}
            ],
            "genomic_sequences": {"chr1": ""}
        }"#;
        std::fs::File::create(&path)
            .unwrap()
            .write_all(json.as_bytes())
            .unwrap();
        let summary = check_transcripts_json(&path).expect("loads (has a transcript)");
        assert!(
            !summary.genome_capable,
            "an empty genomic contig must not be reported as genome-capable"
        );
    }

    #[test]
    fn check_transcripts_json_propagates_incompatible_version_error() {
        use std::io::Write;
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("transcripts.json");
        std::fs::File::create(&path)
            .unwrap()
            .write_all(br#"{"version": "2.0", "transcripts": []}"#)
            .unwrap();

        assert!(
            check_transcripts_json(&path).is_err(),
            "an incompatible schema version must surface as a check failure"
        );
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
            protein_fastas: Vec::new(),
            genome_fasta: None,
            genome_grch37_fasta: None,
            refseqgene_fastas: Vec::new(),
            refseqgene_alignments: None,
            refseqgene_alignments_grch37: None,
            assembly_report: None,
            assembly_report_grch37: None,
            derived_refseqgene_placements: None,
            ng_hosted_transcripts: None,
            derived_transcript_placements: None,
            refseqgene_summary: None,
            lrg_fastas: Vec::new(),
            lrg_xmls: Vec::new(),
            lrg_refseq_mapping: None,
            cdot_json: None,
            cdot_grch37_json: None,
            ensembl_cdot_json: None,
            ensembl_cdot_grch37_json: None,
            ensembl_transcript_fastas: Vec::new(),
            supplemental_fasta: None,
            legacy_transcripts_fasta: None,
            legacy_transcripts_metadata: None,
            legacy_genbank_fasta: None,
            legacy_genbank_metadata: None,
            canonical_overrides: None,
            backfill_transcripts_fasta: None,
            transcript_count: 1,
            available_prefixes: vec!["NM".to_string()],
            reference_identity: None,
            manifest_schema_version: None,
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

    #[test]
    fn test_check_warns_on_missing_refseqgene_alignments() {
        let dir = TempDir::new().unwrap();

        // A real transcript FASTA so the manifest is otherwise valid.
        let transcript_fasta = dir.path().join("example.fna");
        File::create(&transcript_fasta).unwrap();

        // Point refseqgene_alignments at a file that does not exist.
        let missing_alignments = dir.path().join("missing_refseqgene_alignments.gff3");

        let mut manifest = ReferenceManifest {
            prepared_at: "2024-01-01T00:00:00Z".to_string(),
            transcript_fastas: vec![transcript_fasta],
            protein_fastas: Vec::new(),
            genome_fasta: None,
            genome_grch37_fasta: None,
            refseqgene_fastas: Vec::new(),
            refseqgene_alignments: Some(missing_alignments.clone()),
            refseqgene_alignments_grch37: None,
            assembly_report: None,
            assembly_report_grch37: None,
            derived_refseqgene_placements: None,
            ng_hosted_transcripts: None,
            derived_transcript_placements: None,
            refseqgene_summary: None,
            lrg_fastas: Vec::new(),
            lrg_xmls: Vec::new(),
            lrg_refseq_mapping: None,
            cdot_json: None,
            cdot_grch37_json: None,
            ensembl_cdot_json: None,
            ensembl_cdot_grch37_json: None,
            ensembl_transcript_fastas: Vec::new(),
            supplemental_fasta: None,
            legacy_transcripts_fasta: None,
            legacy_transcripts_metadata: None,
            legacy_genbank_fasta: None,
            legacy_genbank_metadata: None,
            canonical_overrides: None,
            backfill_transcripts_fasta: None,
            transcript_count: 1,
            available_prefixes: vec!["NM".to_string()],
            reference_identity: None,
            manifest_schema_version: None,
            reference_dir: dir.path().to_path_buf(),
        };

        // Create the artifact so `save()` (which now rejects a wired-but-
        // unreadable stamped artifact) succeeds, then remove it to reproduce the
        // "artifact went missing after prepare" state that `check` must warn on.
        File::create(&missing_alignments).unwrap();
        manifest.save().unwrap();
        std::fs::remove_file(&missing_alignments).unwrap();

        let result = check_reference(dir.path());
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.contains("RefSeqGene alignments GFF3 not found")),
            "Expected a warning for the missing alignments file, got {:?}",
            result.warnings
        );
    }

    #[test]
    fn test_check_warns_on_missing_refseqgene_alignments_grch37() {
        let dir = TempDir::new().unwrap();

        // A real transcript FASTA so the manifest is otherwise valid.
        let transcript_fasta = dir.path().join("example.fna");
        File::create(&transcript_fasta).unwrap();

        // Point refseqgene_alignments_grch37 at a file that does not exist.
        let missing_alignments = dir.path().join("missing_refseqgene_alignments_grch37.gff3");

        let mut manifest = ReferenceManifest {
            prepared_at: "2024-01-01T00:00:00Z".to_string(),
            transcript_fastas: vec![transcript_fasta],
            protein_fastas: Vec::new(),
            genome_fasta: None,
            genome_grch37_fasta: None,
            refseqgene_fastas: Vec::new(),
            refseqgene_alignments: None,
            refseqgene_alignments_grch37: Some(missing_alignments.clone()),
            assembly_report: None,
            assembly_report_grch37: None,
            derived_refseqgene_placements: None,
            ng_hosted_transcripts: None,
            derived_transcript_placements: None,
            refseqgene_summary: None,
            lrg_fastas: Vec::new(),
            lrg_xmls: Vec::new(),
            lrg_refseq_mapping: None,
            cdot_json: None,
            cdot_grch37_json: None,
            ensembl_cdot_json: None,
            ensembl_cdot_grch37_json: None,
            ensembl_transcript_fastas: Vec::new(),
            supplemental_fasta: None,
            legacy_transcripts_fasta: None,
            legacy_transcripts_metadata: None,
            legacy_genbank_fasta: None,
            legacy_genbank_metadata: None,
            canonical_overrides: None,
            backfill_transcripts_fasta: None,
            transcript_count: 1,
            available_prefixes: vec!["NM".to_string()],
            reference_identity: None,
            manifest_schema_version: None,
            reference_dir: dir.path().to_path_buf(),
        };

        // Create the artifact so `save()` (which now rejects a wired-but-
        // unreadable stamped artifact) succeeds, then remove it to reproduce the
        // "artifact went missing after prepare" state that `check` must warn on.
        File::create(&missing_alignments).unwrap();
        manifest.save().unwrap();
        std::fs::remove_file(&missing_alignments).unwrap();

        let result = check_reference(dir.path());
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.contains("GRCh37 RefSeqGene alignments GFF3 not found")),
            "Expected a warning for the missing GRCh37 alignments file, got {:?}",
            result.warnings
        );
    }

    #[test]
    fn test_check_warns_on_missing_assembly_report() {
        let dir = TempDir::new().unwrap();

        // A real transcript FASTA so the manifest is otherwise valid.
        let transcript_fasta = dir.path().join("example.fna");
        File::create(&transcript_fasta).unwrap();

        // Point assembly_report at a file that does not exist.
        let missing_report = dir.path().join("missing_assembly_report.txt");

        let mut manifest = ReferenceManifest {
            prepared_at: "2024-01-01T00:00:00Z".to_string(),
            transcript_fastas: vec![transcript_fasta],
            protein_fastas: Vec::new(),
            genome_fasta: None,
            genome_grch37_fasta: None,
            refseqgene_fastas: Vec::new(),
            refseqgene_alignments: None,
            refseqgene_alignments_grch37: None,
            assembly_report: Some(missing_report.clone()),
            assembly_report_grch37: None,
            derived_refseqgene_placements: None,
            ng_hosted_transcripts: None,
            derived_transcript_placements: None,
            refseqgene_summary: None,
            lrg_fastas: Vec::new(),
            lrg_xmls: Vec::new(),
            lrg_refseq_mapping: None,
            cdot_json: None,
            cdot_grch37_json: None,
            ensembl_cdot_json: None,
            ensembl_cdot_grch37_json: None,
            ensembl_transcript_fastas: Vec::new(),
            supplemental_fasta: None,
            legacy_transcripts_fasta: None,
            legacy_transcripts_metadata: None,
            legacy_genbank_fasta: None,
            legacy_genbank_metadata: None,
            canonical_overrides: None,
            backfill_transcripts_fasta: None,
            transcript_count: 1,
            available_prefixes: vec!["NM".to_string()],
            reference_identity: None,
            manifest_schema_version: None,
            reference_dir: dir.path().to_path_buf(),
        };

        // Create the artifact so `save()` (which now rejects a wired-but-
        // unreadable stamped artifact) succeeds, then remove it to reproduce the
        // "artifact went missing after prepare" state that `check` must warn on.
        File::create(&missing_report).unwrap();
        manifest.save().unwrap();
        std::fs::remove_file(&missing_report).unwrap();

        let result = check_reference(dir.path());
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.contains("Assembly report not found")),
            "Expected a warning for the missing assembly report file, got {:?}",
            result.warnings
        );
    }

    #[test]
    fn test_check_warns_on_missing_assembly_report_grch37() {
        let dir = TempDir::new().unwrap();

        // A real transcript FASTA so the manifest is otherwise valid.
        let transcript_fasta = dir.path().join("example.fna");
        File::create(&transcript_fasta).unwrap();

        // Point assembly_report_grch37 at a file that does not exist.
        let missing_report = dir.path().join("missing_assembly_report_grch37.txt");

        let mut manifest = ReferenceManifest {
            prepared_at: "2024-01-01T00:00:00Z".to_string(),
            transcript_fastas: vec![transcript_fasta],
            protein_fastas: Vec::new(),
            genome_fasta: None,
            genome_grch37_fasta: None,
            refseqgene_fastas: Vec::new(),
            refseqgene_alignments: None,
            refseqgene_alignments_grch37: None,
            assembly_report: None,
            assembly_report_grch37: Some(missing_report.clone()),
            derived_refseqgene_placements: None,
            ng_hosted_transcripts: None,
            derived_transcript_placements: None,
            refseqgene_summary: None,
            lrg_fastas: Vec::new(),
            lrg_xmls: Vec::new(),
            lrg_refseq_mapping: None,
            cdot_json: None,
            cdot_grch37_json: None,
            ensembl_cdot_json: None,
            ensembl_cdot_grch37_json: None,
            ensembl_transcript_fastas: Vec::new(),
            supplemental_fasta: None,
            legacy_transcripts_fasta: None,
            legacy_transcripts_metadata: None,
            legacy_genbank_fasta: None,
            legacy_genbank_metadata: None,
            canonical_overrides: None,
            backfill_transcripts_fasta: None,
            transcript_count: 1,
            available_prefixes: vec!["NM".to_string()],
            reference_identity: None,
            manifest_schema_version: None,
            reference_dir: dir.path().to_path_buf(),
        };

        // Create the artifact so `save()` (which now rejects a wired-but-
        // unreadable stamped artifact) succeeds, then remove it to reproduce the
        // "artifact went missing after prepare" state that `check` must warn on.
        File::create(&missing_report).unwrap();
        manifest.save().unwrap();
        std::fs::remove_file(&missing_report).unwrap();

        let result = check_reference(dir.path());
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.contains("GRCh37 assembly report not found")),
            "Expected a warning for the missing GRCh37 assembly report file, got {:?}",
            result.warnings
        );
    }

    #[test]
    fn test_check_warns_on_missing_refseqgene_summary() {
        let dir = TempDir::new().unwrap();

        // A real transcript FASTA so the manifest is otherwise valid.
        let transcript_fasta = dir.path().join("example.fna");
        File::create(&transcript_fasta).unwrap();

        // Point refseqgene_summary at a file that does not exist.
        let missing_summary = dir.path().join("missing_LRG_RefSeqGene");

        let mut manifest = ReferenceManifest {
            prepared_at: "2024-01-01T00:00:00Z".to_string(),
            transcript_fastas: vec![transcript_fasta],
            protein_fastas: Vec::new(),
            genome_fasta: None,
            genome_grch37_fasta: None,
            refseqgene_fastas: Vec::new(),
            refseqgene_alignments: None,
            refseqgene_alignments_grch37: None,
            assembly_report: None,
            assembly_report_grch37: None,
            derived_refseqgene_placements: None,
            ng_hosted_transcripts: None,
            derived_transcript_placements: None,
            refseqgene_summary: Some(missing_summary.clone()),
            lrg_fastas: Vec::new(),
            lrg_xmls: Vec::new(),
            lrg_refseq_mapping: None,
            cdot_json: None,
            cdot_grch37_json: None,
            ensembl_cdot_json: None,
            ensembl_cdot_grch37_json: None,
            ensembl_transcript_fastas: Vec::new(),
            supplemental_fasta: None,
            legacy_transcripts_fasta: None,
            legacy_transcripts_metadata: None,
            legacy_genbank_fasta: None,
            legacy_genbank_metadata: None,
            canonical_overrides: None,
            backfill_transcripts_fasta: None,
            transcript_count: 1,
            available_prefixes: vec!["NM".to_string()],
            reference_identity: None,
            manifest_schema_version: None,
            reference_dir: dir.path().to_path_buf(),
        };

        // Create the artifact so `save()` (which now rejects a wired-but-
        // unreadable stamped artifact) succeeds, then remove it to reproduce the
        // "artifact went missing after prepare" state that `check` must warn on.
        File::create(&missing_summary).unwrap();
        manifest.save().unwrap();
        std::fs::remove_file(&missing_summary).unwrap();

        let result = check_reference(dir.path());
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.contains("RefSeqGene summary (LRG_RefSeqGene) not found")),
            "Expected a warning for the missing summary file, got {:?}",
            result.warnings
        );
    }

    #[test]
    fn test_check_warns_on_missing_derived_refseqgene_placements() {
        let dir = TempDir::new().unwrap();

        // A real transcript FASTA so the manifest is otherwise valid.
        let transcript_fasta = dir.path().join("example.fna");
        File::create(&transcript_fasta).unwrap();

        // Point derived_refseqgene_placements at a file that does not exist.
        let missing_placements = dir.path().join("missing_derived_placements.json");

        let mut manifest = ReferenceManifest {
            prepared_at: "2024-01-01T00:00:00Z".to_string(),
            transcript_fastas: vec![transcript_fasta],
            protein_fastas: Vec::new(),
            genome_fasta: None,
            genome_grch37_fasta: None,
            refseqgene_fastas: Vec::new(),
            refseqgene_alignments: None,
            refseqgene_alignments_grch37: None,
            assembly_report: None,
            assembly_report_grch37: None,
            derived_refseqgene_placements: Some(missing_placements.clone()),
            ng_hosted_transcripts: None,
            derived_transcript_placements: None,
            refseqgene_summary: None,
            lrg_fastas: Vec::new(),
            lrg_xmls: Vec::new(),
            lrg_refseq_mapping: None,
            cdot_json: None,
            cdot_grch37_json: None,
            ensembl_cdot_json: None,
            ensembl_cdot_grch37_json: None,
            ensembl_transcript_fastas: Vec::new(),
            supplemental_fasta: None,
            legacy_transcripts_fasta: None,
            legacy_transcripts_metadata: None,
            legacy_genbank_fasta: None,
            legacy_genbank_metadata: None,
            canonical_overrides: None,
            backfill_transcripts_fasta: None,
            transcript_count: 1,
            available_prefixes: vec!["NM".to_string()],
            reference_identity: None,
            manifest_schema_version: None,
            reference_dir: dir.path().to_path_buf(),
        };

        // Create the artifact so `save()` (which now rejects a wired-but-
        // unreadable stamped artifact) succeeds, then remove it to reproduce the
        // "artifact went missing after prepare" state that `check` must warn on.
        File::create(&missing_placements).unwrap();
        manifest.save().unwrap();
        std::fs::remove_file(&missing_placements).unwrap();

        let result = check_reference(dir.path());
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.contains("Derived RefSeqGene placements not found")),
            "Expected a warning for the missing derived placements file, got {:?}",
            result.warnings
        );
    }

    #[test]
    fn test_check_warns_on_missing_ng_hosted_transcripts() {
        let dir = TempDir::new().unwrap();

        // A real transcript FASTA so the manifest is otherwise valid.
        let transcript_fasta = dir.path().join("example.fna");
        File::create(&transcript_fasta).unwrap();

        // Point ng_hosted_transcripts at a file that does not exist.
        let missing_ng_hosted = dir.path().join("missing_ng_hosted_transcripts.json");

        let mut manifest = ReferenceManifest {
            prepared_at: "2024-01-01T00:00:00Z".to_string(),
            transcript_fastas: vec![transcript_fasta],
            protein_fastas: Vec::new(),
            genome_fasta: None,
            genome_grch37_fasta: None,
            refseqgene_fastas: Vec::new(),
            refseqgene_alignments: None,
            refseqgene_alignments_grch37: None,
            assembly_report: None,
            assembly_report_grch37: None,
            derived_refseqgene_placements: None,
            ng_hosted_transcripts: Some(missing_ng_hosted.clone()),
            derived_transcript_placements: None,
            refseqgene_summary: None,
            lrg_fastas: Vec::new(),
            lrg_xmls: Vec::new(),
            lrg_refseq_mapping: None,
            cdot_json: None,
            cdot_grch37_json: None,
            ensembl_cdot_json: None,
            ensembl_cdot_grch37_json: None,
            ensembl_transcript_fastas: Vec::new(),
            supplemental_fasta: None,
            legacy_transcripts_fasta: None,
            legacy_transcripts_metadata: None,
            legacy_genbank_fasta: None,
            legacy_genbank_metadata: None,
            canonical_overrides: None,
            backfill_transcripts_fasta: None,
            transcript_count: 1,
            available_prefixes: vec!["NM".to_string()],
            reference_identity: None,
            manifest_schema_version: None,
            reference_dir: dir.path().to_path_buf(),
        };

        // Create the artifact so `save()` (which now rejects a wired-but-
        // unreadable stamped artifact) succeeds, then remove it to reproduce the
        // "artifact went missing after prepare" state that `check` must warn on.
        File::create(&missing_ng_hosted).unwrap();
        manifest.save().unwrap();
        std::fs::remove_file(&missing_ng_hosted).unwrap();

        let result = check_reference(dir.path());
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.contains("NG_ hosted-transcripts map not found")),
            "Expected a warning for the missing NG_ hosted-transcripts map, got {:?}",
            result.warnings
        );
    }
}
