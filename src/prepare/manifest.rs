//! Reference manifest types and I/O operations.
//!
//! This module handles the definition, loading, and display of the `ReferenceManifest`
//! that tracks all prepared reference data files.

use crate::FerroError;
use std::fs::File;
use std::path::{Component, Path, PathBuf};

/// Current reference-manifest schema version, written into every manifest by
/// [`ReferenceManifest::save`]. Bump this on any breaking change to the manifest
/// or derived-artifact format so an older build refuses a newer, incompatible
/// reference at load rather than silently misreading it. A manifest whose
/// `manifest_schema_version` exceeds this value is rejected (see
/// [`check_schema_version`]); an absent value means the reference predates schema
/// versioning and is accepted (its structure is still validated on load).
///
/// Invariant: `ReferenceManifest` uses `#[serde(deny_unknown_fields)]`, so any
/// *new* manifest field is an unknown field to older binaries. Every field
/// addition MUST bump this constant, so an older binary reading a newer
/// reference reports the actionable "prepared by a newer ferro — upgrade"
/// version message (the version gate runs before the strict deserialize) rather
/// than an opaque serde "unknown field" error. A future `MIN_SUPPORTED` floor
/// (currently unneeded — there is no older incompatible schema yet) would refuse
/// references older than this build can safely read.
pub const CURRENT_MANIFEST_SCHEMA_VERSION: u32 = 3;

/// Manifest of prepared reference data.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceManifest {
    /// When the data was prepared
    pub prepared_at: String,
    /// Transcript FASTA files
    pub transcript_fastas: Vec<PathBuf>,
    /// Protein FASTA files (`*.protein.faa`, NP_/XP_ accessions) for the
    /// translated-CDS-vs-canonical-protein check (issue #520).
    #[serde(default)]
    pub protein_fastas: Vec<PathBuf>,
    /// GRCh38 genome FASTA file (if downloaded)
    pub genome_fasta: Option<PathBuf>,
    /// GRCh37 genome FASTA file (if downloaded)
    #[serde(default)]
    pub genome_grch37_fasta: Option<PathBuf>,
    /// RefSeqGene FASTA files (NG_* accessions)
    #[serde(default)]
    pub refseqgene_fastas: Vec<PathBuf>,
    /// RefSeqGene→genome alignment GFF3 (NCBI `GCF_*_refseqgene_alignments.gff3`),
    /// giving each NG_ record's **GRCh38** chromosomal placement. Consumed to
    /// re-express transcript coordinates in an NG_ parent's own frame (#480).
    #[serde(default)]
    pub refseqgene_alignments: Option<PathBuf>,
    /// GRCh37 RefSeqGene→genome alignment GFF3 (the archived `GCF_000001405.25_105.*`
    /// snapshot). Merged with `refseqgene_alignments` so an NG_ resolves to its
    /// build-appropriate placement (#653/#713). Paired with `refseqgene_alignments`
    /// the same way `genome_grch37_fasta`/`cdot_grch37_json` pair their GRCh38 fields.
    #[serde(default)]
    pub refseqgene_alignments_grch37: Option<PathBuf>,
    /// Path to the GRCh38 NCBI `*_assembly_report.txt` (#716). Authoritative
    /// RefSeq-accession → assembly map for build inference; absent on manifests
    /// prepared before #716 (build inference then uses the hardcoded fallback).
    #[serde(default)]
    pub assembly_report: Option<PathBuf>,
    /// Path to the GRCh37 NCBI `*_assembly_report.txt` (#716). Present only when
    /// the GRCh37 genome was downloaded (`download_genome_grch37`).
    #[serde(default)]
    pub assembly_report_grch37: Option<PathBuf>,
    /// Derived `NG_`/`LRG_` placements (`derived_placement::DerivedPlacements`
    /// JSON) for versions absent from the alignment snapshots above — built at
    /// prepare time from each NG_'s GenBank record + cdot and validated by full
    /// sequence comparison (#728). Merged into the placement map *after* the
    /// authoritative GFF3 snapshots (first-source-wins, so a GFF3 entry is never
    /// overridden); fills only genuine version gaps.
    #[serde(default)]
    pub derived_refseqgene_placements: Option<PathBuf>,
    /// Per-`NG_`-version hosted-transcript map (`ng_hosted_transcripts.json`,
    /// #792), for NG_-parent-aware legacy `GENE_v001` selector resolution.
    #[serde(default)]
    pub ng_hosted_transcripts: Option<PathBuf>,
    /// Derived exon→genome structures (JSON) for old transcript versions absent
    /// from cdot — built at prepare time by sibling-anchoring + genome-readback
    /// validation (#790). Injected into the cdot mapper at load; fills only
    /// genuine version gaps (declines rather than mis-anchor).
    #[serde(default)]
    pub derived_transcript_placements: Option<PathBuf>,
    /// NCBI `LRG_RefSeqGene` association table, mapping each gene to its
    /// reference-standard transcript. Consumed to resolve legacy gene-model
    /// selectors (`NG_(GENE_v001):c.…`) to the transcript accession (#500/#637).
    #[serde(default)]
    pub refseqgene_summary: Option<PathBuf>,
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
    /// Ensembl cdot transcript metadata JSON for GRCh38 (ENST/ENSG/ENSP), if
    /// downloaded via `ferro prepare --ensembl`. Same schema as `cdot_json`.
    #[serde(default)]
    pub ensembl_cdot_json: Option<PathBuf>,
    /// Ensembl cdot transcript metadata JSON for GRCh37, if downloaded.
    #[serde(default)]
    pub ensembl_cdot_grch37_json: Option<PathBuf>,
    /// Ensembl cDNA FASTA files (ENST_* transcript sequences), if downloaded
    /// via `ferro prepare --ensembl`.
    #[serde(default)]
    pub ensembl_transcript_fastas: Vec<PathBuf>,
    /// Supplemental FASTA file (missing ClinVar transcripts fetched from NCBI)
    #[serde(default)]
    pub supplemental_fasta: Option<PathBuf>,
    /// Version-aware backfill FASTA (#842): deposited sequences for
    /// `accession.version`s that are present in cdot but absent from the bulk
    /// RefSeq RNA feed, fetched by `ferro prepare --backfill-transcripts`. The
    /// provider ingests it as a transcript FASTA so the primary path serves
    /// those versions directly instead of synthesizing lossy bases from the
    /// genome via the cdot exon CIGAR.
    #[serde(default)]
    pub backfill_transcripts_fasta: Option<PathBuf>,
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
    /// Canonical-overrides JSON: authoritative `(cds_start, cds_end, protein_id,
    /// tx_length)` per exact accession version, fetched by
    /// `ferro prepare --validate-canonical`. Consumed by the canonical-record
    /// validation / correction path (issue #520).
    #[serde(default)]
    pub canonical_overrides: Option<PathBuf>,
    /// Total number of transcripts
    pub transcript_count: usize,
    /// List of available accession prefixes
    pub available_prefixes: Vec<String>,
    /// The pinned cdot data release the reference's cdot artifacts came from
    /// (e.g. `"data_v0.2.32"`), recorded by `ferro prepare` (#1001). Provenance
    /// only — not verified at load and deliberately NOT part of the reference
    /// identity (the cdot basename, which is hashed, already encodes it).
    #[serde(default)]
    pub cdot_data_version: Option<String>,
    /// Content fingerprint (#1001), written by [`ReferenceManifest::save`].
    /// Absent on references prepared before content-hashing; verified fail-soft
    /// at load (warn by default, hard-fail under `--strict-reference`).
    #[serde(default)]
    pub reference_identity: Option<String>,
    /// Schema version of this manifest, written by [`ReferenceManifest::save`].
    /// Absent on references prepared before schema versioning; a value greater
    /// than [`CURRENT_MANIFEST_SCHEMA_VERSION`] means the reference was prepared
    /// by a newer, forward-incompatible `ferro` and is refused at load.
    #[serde(default)]
    pub manifest_schema_version: Option<u32>,
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
            transcript_count: 0,
            available_prefixes: Vec::new(),
            cdot_data_version: None,
            reference_identity: None,
            manifest_schema_version: None,
            reference_dir: PathBuf::new(),
        }
    }
}

impl ReferenceManifest {
    /// Load manifest from directory, or create a fresh one if it doesn't exist.
    pub fn load_or_default(reference_dir: &Path) -> Result<Self, FerroError> {
        let manifest_path = reference_dir.join("manifest.json");

        let mut manifest = if manifest_path.exists() {
            let bytes = std::fs::read(&manifest_path).map_err(|e| FerroError::Io {
                msg: format!("Failed to open manifest: {}", e),
            })?;
            // Parse untyped first so the load gate can check the schema VERSION
            // before the strict typed deserialize. A newer reference then reports
            // the actionable "upgrade ferro" message instead of an opaque serde
            // "unknown field" error (which `deny_unknown_fields` would raise first
            // if we deserialized straight into the struct). This mirrors the
            // runtime provider path (`MultiFastaProvider::from_manifest`), which
            // also validates the value before reading it.
            let value: serde_json::Value =
                serde_json::from_slice(&bytes).map_err(|e| FerroError::Io {
                    msg: format!("Failed to parse manifest: {}", e),
                })?;
            validate_loaded_manifest(&value)?;
            // Do the authoritative typed deserialize from the ORIGINAL BYTES, not
            // from `value`: a `serde_json::Value` round-trip silently collapses
            // duplicate JSON keys (last-wins), so `from_value` would lose serde's
            // duplicate-field rejection. Reading the bytes keeps a duplicate-key
            // manifest a hard error, consistent with this feature's "fail loud on
            // a misread reference" intent (#1001).
            serde_json::from_slice(&bytes).map_err(|e| FerroError::Io {
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
        self.manifest_schema_version = Some(CURRENT_MANIFEST_SCHEMA_VERSION);
        self.deduplicate_paths();

        // Serialize a relative-path view without mutating the in-memory absolute paths.
        let mut on_disk = self.clone();
        on_disk.make_paths_relative();

        // Stamp the content identity (#1001) onto the on-disk (relative-path)
        // view. Computed from a throwaway Value with the derived-artifact
        // content-stamps injected (reading the small artifacts from
        // reference_dir). The signature excludes `reference_identity` and
        // `prepared_at`, so this is self-consistent and idempotent on a re-bless
        // of byte-identical data. `derived_artifact_stamps` is NOT persisted —
        // it lives only in the throwaway value the hash consumes.
        let on_disk_value = serde_json::to_value(&on_disk).map_err(|e| FerroError::Io {
            msg: format!("Failed to serialize manifest for identity stamp: {e}"),
        })?;
        // Reject a wired-but-unreadable stamped artifact BEFORE hashing, so the
        // writer is symmetric with the load-time verifier: `inject_content_stamps`
        // silently drops an unreadable artifact, so without this pre-check `save()`
        // could mint a stamp over an artifact it couldn't read that the loader then
        // hard-rejects (an asymmetry the "re-run ferro prepare" remedy can't fix).
        crate::prepare::identity::ensure_stamped_artifacts_readable(
            &self.reference_dir,
            &on_disk_value,
        )?;
        let id = crate::prepare::identity::reference_identity(&self.reference_dir, &on_disk_value);
        on_disk.reference_identity = Some(id.clone());
        self.reference_identity = Some(id);

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
            .chain(self.protein_fastas.iter())
            .chain(self.refseqgene_fastas.iter())
            .chain(self.lrg_fastas.iter())
            .chain(self.lrg_xmls.iter())
            .chain(self.ensembl_transcript_fastas.iter())
        {
            f(p);
        }
        for p in [
            &self.genome_fasta,
            &self.genome_grch37_fasta,
            &self.refseqgene_alignments,
            &self.refseqgene_alignments_grch37,
            &self.assembly_report,
            &self.assembly_report_grch37,
            &self.derived_refseqgene_placements,
            &self.ng_hosted_transcripts,
            &self.derived_transcript_placements,
            &self.refseqgene_summary,
            &self.lrg_refseq_mapping,
            &self.cdot_json,
            &self.cdot_grch37_json,
            &self.ensembl_cdot_json,
            &self.ensembl_cdot_grch37_json,
            &self.supplemental_fasta,
            &self.legacy_transcripts_fasta,
            &self.legacy_transcripts_metadata,
            &self.legacy_genbank_fasta,
            &self.legacy_genbank_metadata,
            &self.canonical_overrides,
            &self.backfill_transcripts_fasta,
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
            &mut self.protein_fastas,
            &mut self.refseqgene_fastas,
            &mut self.lrg_fastas,
            &mut self.lrg_xmls,
            &mut self.ensembl_transcript_fastas,
        ] {
            for p in v.iter_mut() {
                f(p);
            }
        }
        for o in [
            &mut self.genome_fasta,
            &mut self.genome_grch37_fasta,
            &mut self.refseqgene_alignments,
            &mut self.refseqgene_alignments_grch37,
            &mut self.assembly_report,
            &mut self.assembly_report_grch37,
            &mut self.derived_refseqgene_placements,
            &mut self.ng_hosted_transcripts,
            &mut self.derived_transcript_placements,
            &mut self.refseqgene_summary,
            &mut self.lrg_refseq_mapping,
            &mut self.cdot_json,
            &mut self.cdot_grch37_json,
            &mut self.ensembl_cdot_json,
            &mut self.ensembl_cdot_grch37_json,
            &mut self.supplemental_fasta,
            &mut self.legacy_transcripts_fasta,
            &mut self.legacy_transcripts_metadata,
            &mut self.legacy_genbank_fasta,
            &mut self.legacy_genbank_metadata,
            &mut self.canonical_overrides,
            &mut self.backfill_transcripts_fasta,
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
            &mut self.protein_fastas,
            &mut self.refseqgene_fastas,
            &mut self.lrg_fastas,
            &mut self.lrg_xmls,
            &mut self.ensembl_transcript_fastas,
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

/// Reject a manifest prepared by a newer, forward-incompatible `ferro`.
///
/// An absent version (`None`) predates schema versioning and is accepted; a
/// version at or below [`CURRENT_MANIFEST_SCHEMA_VERSION`] is accepted; a higher
/// version is refused, because an older build cannot be trusted to read a newer
/// format correctly (it would otherwise silently misread it).
pub fn check_schema_version(version: Option<u32>) -> Result<(), FerroError> {
    if let Some(v) = version {
        if v > CURRENT_MANIFEST_SCHEMA_VERSION {
            return Err(FerroError::Io {
                msg: format!(
                    "reference manifest schema version {v} is newer than this build of ferro \
                     supports (maximum {CURRENT_MANIFEST_SCHEMA_VERSION}). Upgrade ferro, or \
                     re-run `ferro prepare` to regenerate the reference."
                ),
            });
        }
    }
    Ok(())
}

/// Validate a raw manifest JSON value before the permissive field reads in the
/// runtime loader ([`crate::reference::multi_fasta::MultiFastaProvider::from_manifest`]).
///
/// The runtime loader reads the manifest as an untyped `serde_json::Value` and
/// accesses each field with `.get(...).and_then(...)`, so a renamed, removed, or
/// wrong-typed field would otherwise degrade to `None` and be silently treated as
/// absent — producing wrong output across a run instead of failing. This gate:
///
/// 1. refuses a manifest from a newer, forward-incompatible `ferro`
///    ([`check_schema_version`]); and
/// 2. validates the value against the [`ReferenceManifest`] schema, catching a
///    missing required field, a wrong-typed field, or an unknown/typo'd field
///    (`#[serde(deny_unknown_fields)]`, #1001) — an unmodeled optional key would
///    otherwise silently deserialize as `None` and the run would proceed on
///    missing data.
pub fn validate_loaded_manifest(value: &serde_json::Value) -> Result<(), FerroError> {
    let version = value
        .get("manifest_schema_version")
        .and_then(|v| v.as_u64())
        .map(|v| v as u32);
    check_schema_version(version)?;

    serde_json::from_value::<ReferenceManifest>(value.clone()).map_err(|e| FerroError::Io {
        msg: format!(
            "reference manifest does not match the expected schema: {e}. The reference may have \
             been prepared by an incompatible version of ferro — re-run `ferro prepare` to \
             regenerate it."
        ),
    })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_paths_relative() {
        use tempfile::TempDir;

        let temp_dir = TempDir::new().unwrap();
        let ref_path = temp_dir.path().to_path_buf();

        let mut manifest = ReferenceManifest {
            reference_dir: ref_path.clone(),
            transcript_fastas: vec![ref_path.join("transcripts.fa")],
            protein_fastas: vec![ref_path.join("proteins.faa")],
            cdot_json: Some(ref_path.join("cdot.json")),
            canonical_overrides: Some(ref_path.join("canonical_overrides.json")),
            ..Default::default()
        };

        manifest.make_paths_relative();

        assert_eq!(
            manifest.transcript_fastas[0],
            PathBuf::from("transcripts.fa")
        );
        assert_eq!(
            manifest.protein_fastas[0],
            PathBuf::from("proteins.faa"),
            "protein_fastas must be normalized to a relative path"
        );
        assert_eq!(manifest.cdot_json, Some(PathBuf::from("cdot.json")));
        assert_eq!(
            manifest.canonical_overrides,
            Some(PathBuf::from("canonical_overrides.json")),
            "canonical_overrides must be normalized to a relative path"
        );
    }

    /// A manifest written before #713 (no `refseqgene_alignments_grch37` key)
    /// must still deserialize, with the new field defaulting to `None`
    /// (`#[serde(default)]`) — backward compatibility for prepared references.
    #[test]
    fn test_legacy_manifest_without_grch37_alignments_deserializes() {
        let json = r#"{
            "prepared_at": "2024-01-01T00:00:00Z",
            "transcript_fastas": ["transcripts.fa"],
            "genome_fasta": "genome.fa",
            "refseqgene_alignments": "refseqgene/aln.gff3",
            "transcript_count": 0,
            "available_prefixes": []
        }"#;
        let manifest: ReferenceManifest =
            serde_json::from_str(json).expect("legacy manifest must deserialize");
        assert_eq!(
            manifest.refseqgene_alignments,
            Some(PathBuf::from("refseqgene/aln.gff3")),
            "the GRCh38 alignment field is preserved"
        );
        assert_eq!(
            manifest.refseqgene_alignments_grch37, None,
            "absent GRCh37 alignment field defaults to None"
        );
    }

    #[test]
    fn test_make_paths_absolute() {
        use tempfile::TempDir;

        let temp_dir = TempDir::new().unwrap();
        let ref_path = temp_dir.path().to_path_buf();

        let mut manifest = ReferenceManifest {
            reference_dir: ref_path.clone(),
            transcript_fastas: vec![PathBuf::from("transcripts.fa")],
            protein_fastas: vec![PathBuf::from("proteins.faa")],
            cdot_json: Some(PathBuf::from("cdot.json")),
            canonical_overrides: Some(PathBuf::from("canonical_overrides.json")),
            ..Default::default()
        };

        manifest.make_paths_absolute();

        assert_eq!(
            manifest.transcript_fastas[0],
            ref_path.join("transcripts.fa")
        );
        assert_eq!(
            manifest.protein_fastas[0],
            ref_path.join("proteins.faa"),
            "protein_fastas must be resolved to an absolute path"
        );
        assert_eq!(manifest.cdot_json, Some(ref_path.join("cdot.json")));
        assert_eq!(
            manifest.canonical_overrides,
            Some(ref_path.join("canonical_overrides.json")),
            "canonical_overrides must be resolved to an absolute path"
        );
    }

    #[test]
    fn test_deduplicate_paths() {
        let mut manifest = ReferenceManifest {
            transcript_fastas: vec![
                PathBuf::from("b.fa"),
                PathBuf::from("a.fa"),
                PathBuf::from("b.fa"),
            ],
            protein_fastas: vec![
                PathBuf::from("q.faa"),
                PathBuf::from("p.faa"),
                PathBuf::from("q.faa"),
            ],
            lrg_fastas: vec![PathBuf::from("lrg.fa"), PathBuf::from("lrg.fa")],
            ..Default::default()
        };

        manifest.deduplicate_paths();

        assert_eq!(
            manifest.transcript_fastas,
            vec![PathBuf::from("a.fa"), PathBuf::from("b.fa")]
        );
        assert_eq!(
            manifest.protein_fastas,
            vec![PathBuf::from("p.faa"), PathBuf::from("q.faa")]
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
            protein_fastas: Vec::new(),
            genome_fasta: Some(ref_dir.join("genome.fa")),
            genome_grch37_fasta: Some(ref_dir.join("genome37.fa")),
            refseqgene_fastas: vec![ref_dir.join("ng.fa")],
            refseqgene_alignments: Some(ref_dir.join("refseqgene_alignments.gff3")),
            refseqgene_alignments_grch37: Some(ref_dir.join("refseqgene_alignments_grch37.gff3")),
            assembly_report: None,
            assembly_report_grch37: None,
            derived_refseqgene_placements: None,
            ng_hosted_transcripts: None,
            derived_transcript_placements: None,
            refseqgene_summary: Some(ref_dir.join("LRG_RefSeqGene")),
            lrg_fastas: vec![ref_dir.join("lrg.fa")],
            lrg_xmls: vec![ref_dir.join("lrg.xml")],
            lrg_refseq_mapping: Some(ref_dir.join("lrg_mapping.txt")),
            cdot_json: Some(ref_dir.join("cdot.json")),
            cdot_grch37_json: Some(ref_dir.join("cdot37.json")),
            ensembl_cdot_json: Some(ref_dir.join("ensembl_cdot.json")),
            ensembl_cdot_grch37_json: Some(ref_dir.join("ensembl_cdot37.json")),
            ensembl_transcript_fastas: vec![ref_dir.join("ensembl_cdna.fa")],
            supplemental_fasta: Some(ref_dir.join("supplemental.fa")),
            legacy_transcripts_fasta: Some(ref_dir.join("legacy.fa")),
            legacy_transcripts_metadata: Some(ref_dir.join("legacy.json")),
            legacy_genbank_fasta: Some(ref_dir.join("genbank.fa")),
            legacy_genbank_metadata: Some(ref_dir.join("genbank.json")),
            canonical_overrides: None,
            backfill_transcripts_fasta: Some(ref_dir.join("backfill/backfill_transcripts.fna")),
            transcript_count: 100,
            available_prefixes: vec!["NM".to_string()],
            cdot_data_version: None,
            reference_identity: None,
            manifest_schema_version: None,
            reference_dir: ref_dir.to_path_buf(),
        };

        // The content-stamped artifacts must exist on disk: `save()` now rejects a
        // wired-but-unreadable stamped artifact before hashing. (The bulk artifacts
        // — cdot/genome/transcript FASTAs — are not stamped, so they need not
        // exist for this path-roundtrip test.)
        for stamped in [
            "refseqgene_alignments.gff3",
            "refseqgene_alignments_grch37.gff3",
            "LRG_RefSeqGene",
            "lrg_mapping.txt",
            "legacy.fa",
            "legacy.json",
            "genbank.fa",
            "genbank.json",
        ] {
            File::create(ref_dir.join(stamped)).unwrap();
        }
        std::fs::create_dir_all(ref_dir.join("backfill")).unwrap();
        File::create(ref_dir.join("backfill/backfill_transcripts.fna")).unwrap();

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

        // Ensembl fields round-trip: relative on disk, absolute after load.
        assert_eq!(
            json["ensembl_cdot_json"], "ensembl_cdot.json",
            "ensembl_cdot_json should be stored as relative path"
        );
        assert_eq!(
            json["ensembl_transcript_fastas"][0], "ensembl_cdna.fa",
            "ensembl_transcript_fastas should be stored as relative path"
        );
        assert_eq!(
            loaded.ensembl_cdot_json,
            Some(ref_dir.join("ensembl_cdot.json")),
            "ensembl_cdot_json should be absolute after load"
        );
        assert_eq!(
            loaded.ensembl_cdot_grch37_json,
            Some(ref_dir.join("ensembl_cdot37.json")),
            "ensembl_cdot_grch37_json should be absolute after load"
        );
        assert_eq!(
            loaded.ensembl_transcript_fastas,
            vec![ref_dir.join("ensembl_cdna.fa")],
            "ensembl_transcript_fastas should be absolute after load"
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
    fn save_preserves_derived_refseqgene_placements_across_reload() {
        let dir = tempfile::tempdir().unwrap();
        let placements = dir.path().join("derived_refseqgene_placements.json");
        std::fs::write(&placements, "{}").unwrap();

        // First run wires the field (as a hand-edit or a --derive-ng-placements run would).
        let mut m = ReferenceManifest::load_or_default(dir.path()).unwrap();
        m.derived_refseqgene_placements = Some(placements.clone());
        m.save().unwrap();

        // A subsequent prepare loads the existing manifest — the field must survive.
        let reloaded = ReferenceManifest::load_or_default(dir.path()).unwrap();
        assert_eq!(
            reloaded
                .derived_refseqgene_placements
                .map(|p| p.file_name().unwrap().to_owned()),
            Some(std::ffi::OsString::from(
                "derived_refseqgene_placements.json"
            )),
            "derived_refseqgene_placements must be preserved across load_or_default/save"
        );
    }

    #[test]
    fn save_preserves_derived_transcript_placements_across_reload() {
        let dir = tempfile::tempdir().unwrap();
        let placements = dir.path().join("derived_transcript_placements.json");
        std::fs::write(&placements, "{}").unwrap();

        // First run wires the field (as a hand-edit or a --derive-transcript-placements run would).
        let mut m = ReferenceManifest::load_or_default(dir.path()).unwrap();
        m.derived_transcript_placements = Some(placements.clone());
        m.save().unwrap();

        // A subsequent prepare loads the existing manifest — the field must survive.
        let reloaded = ReferenceManifest::load_or_default(dir.path()).unwrap();
        assert_eq!(
            reloaded
                .derived_transcript_placements
                .map(|p| p.file_name().unwrap().to_owned()),
            Some(std::ffi::OsString::from(
                "derived_transcript_placements.json"
            )),
            "derived_transcript_placements must be preserved across load_or_default/save"
        );
    }

    #[test]
    fn assembly_report_paths_roundtrip_and_relativize() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_dir = tmp.path().to_path_buf();

        let mut manifest = ReferenceManifest {
            reference_dir: ref_dir.clone(),
            assembly_report: Some(ref_dir.join("genome/GRCh38.assembly_report.txt")),
            assembly_report_grch37: Some(ref_dir.join("genome/GRCh37.assembly_report.txt")),
            ..Default::default()
        };
        // The stamped artifacts must exist on disk: `save()` now rejects a
        // wired-but-unreadable stamped artifact before hashing.
        std::fs::create_dir_all(ref_dir.join("genome")).unwrap();
        std::fs::File::create(ref_dir.join("genome/GRCh38.assembly_report.txt")).unwrap();
        std::fs::File::create(ref_dir.join("genome/GRCh37.assembly_report.txt")).unwrap();
        manifest.save().unwrap();

        // On-disk JSON must store the report paths relative to reference_dir.
        let raw = std::fs::read_to_string(ref_dir.join("manifest.json")).unwrap();
        assert!(raw.contains("genome/GRCh38.assembly_report.txt"));
        assert!(!raw.contains(ref_dir.to_str().unwrap()));

        // Reloading absolutizes them back under reference_dir.
        let loaded = ReferenceManifest::load_or_default(&ref_dir).unwrap();
        assert_eq!(
            loaded.assembly_report,
            Some(ref_dir.join("genome/GRCh38.assembly_report.txt"))
        );
        assert_eq!(
            loaded.assembly_report_grch37,
            Some(ref_dir.join("genome/GRCh37.assembly_report.txt"))
        );
    }

    #[test]
    fn ng_hosted_transcripts_path_roundtrips() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_dir = tmp.path().to_path_buf();
        let mut m = ReferenceManifest {
            reference_dir: ref_dir.clone(),
            ng_hosted_transcripts: Some(ref_dir.join("ng_hosted_transcripts.json")),
            ..Default::default()
        };
        // The stamped artifact must exist on disk: `save()` now rejects a
        // wired-but-unreadable stamped artifact before hashing.
        std::fs::File::create(ref_dir.join("ng_hosted_transcripts.json")).unwrap();
        m.save().unwrap();
        let raw = std::fs::read_to_string(ref_dir.join("manifest.json")).unwrap();
        assert!(raw.contains("ng_hosted_transcripts.json"));
        assert!(!raw.contains(ref_dir.to_str().unwrap()));
        let loaded = ReferenceManifest::load_or_default(&ref_dir).unwrap();
        assert_eq!(
            loaded.ng_hosted_transcripts,
            Some(ref_dir.join("ng_hosted_transcripts.json"))
        );
    }

    #[test]
    fn backfill_transcripts_fasta_path_roundtrips() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_dir = tmp.path().to_path_buf();
        let mut m = ReferenceManifest {
            reference_dir: ref_dir.clone(),
            backfill_transcripts_fasta: Some(ref_dir.join("backfill/backfill_transcripts.fna")),
            ..Default::default()
        };
        // The stamped artifact must exist on disk: `save()` now rejects a
        // wired-but-unreadable stamped artifact before hashing.
        std::fs::create_dir_all(ref_dir.join("backfill")).unwrap();
        std::fs::File::create(ref_dir.join("backfill/backfill_transcripts.fna")).unwrap();
        m.save().unwrap();
        // On-disk JSON stores the path relative to reference_dir.
        let raw = std::fs::read_to_string(ref_dir.join("manifest.json")).unwrap();
        assert!(raw.contains("backfill/backfill_transcripts.fna"));
        assert!(!raw.contains(ref_dir.to_str().unwrap()));
        // Reloading absolutizes it back under reference_dir.
        let loaded = ReferenceManifest::load_or_default(&ref_dir).unwrap();
        assert_eq!(
            loaded.backfill_transcripts_fasta,
            Some(ref_dir.join("backfill/backfill_transcripts.fna"))
        );
    }

    #[test]
    fn pre_842_manifest_without_backfill_field_loads() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_dir = tmp.path();
        // A manifest JSON with NO backfill_transcripts_fasta key (pre-#842).
        std::fs::write(
            ref_dir.join("manifest.json"),
            r#"{"prepared_at":"","transcript_fastas":[],"genome_fasta":null,"cdot_json":null,"transcript_count":0,"available_prefixes":[]}"#,
        )
        .unwrap();
        let m = ReferenceManifest::load_or_default(ref_dir).unwrap();
        assert_eq!(m.backfill_transcripts_fasta, None);
    }

    #[test]
    fn pre_716_manifest_without_report_fields_loads() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_dir = tmp.path();
        // A manifest JSON with NO assembly_report keys (pre-#716).
        std::fs::write(
            ref_dir.join("manifest.json"),
            r#"{"prepared_at":"","transcript_fastas":[],"genome_fasta":null,"cdot_json":null,"transcript_count":0,"available_prefixes":[]}"#,
        )
        .unwrap();
        let m = ReferenceManifest::load_or_default(ref_dir).unwrap();
        assert_eq!(m.assembly_report, None);
        assert_eq!(m.assembly_report_grch37, None);
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

    /// A minimal but schema-valid manifest as a JSON value, for load-gate tests.
    fn minimal_manifest_json() -> serde_json::Value {
        serde_json::json!({
            "prepared_at": "2024-01-01T00:00:00Z",
            "transcript_fastas": ["transcripts.fa"],
            "genome_fasta": null,
            "cdot_json": null,
            "transcript_count": 0,
            "available_prefixes": [],
            "manifest_schema_version": CURRENT_MANIFEST_SCHEMA_VERSION
        })
    }

    #[test]
    fn save_stamps_current_schema_version() {
        let dir = tempfile::tempdir().unwrap();
        let mut m = ReferenceManifest {
            reference_dir: dir.path().to_path_buf(),
            ..Default::default()
        };
        m.save().unwrap();

        let raw = std::fs::read_to_string(dir.path().join("manifest.json")).unwrap();
        let json: serde_json::Value = serde_json::from_str(&raw).unwrap();
        assert_eq!(
            json["manifest_schema_version"].as_u64(),
            Some(CURRENT_MANIFEST_SCHEMA_VERSION as u64),
            "save() must stamp the current schema version into the manifest"
        );
        assert_eq!(
            m.manifest_schema_version,
            Some(CURRENT_MANIFEST_SCHEMA_VERSION),
            "save() must also update the in-memory schema version"
        );
    }

    #[test]
    fn check_schema_version_accepts_absent_and_current_but_rejects_newer() {
        assert!(
            check_schema_version(None).is_ok(),
            "a pre-versioning reference (absent version) is accepted"
        );
        assert!(check_schema_version(Some(CURRENT_MANIFEST_SCHEMA_VERSION)).is_ok());
        let err = check_schema_version(Some(CURRENT_MANIFEST_SCHEMA_VERSION + 1))
            .expect_err("a newer schema version must be refused");
        let msg = format!("{err}");
        assert!(
            msg.contains("newer than this build"),
            "actionable message: {msg}"
        );
        assert!(
            msg.contains("ferro prepare"),
            "message points to the remedy: {msg}"
        );
    }

    #[test]
    fn validate_loaded_manifest_accepts_clean_manifest() {
        validate_loaded_manifest(&minimal_manifest_json())
            .expect("a schema-valid manifest must load");
    }

    #[test]
    fn validate_loaded_manifest_rejects_unknown_field() {
        // A manifest carrying a key the struct does not model — e.g. a typo'd
        // optional key (`cdot_jsonn`) or drift-era metadata — must now be rejected,
        // not silently tolerated (#1001: an unmodeled optional key otherwise
        // deserializes as `None` and the run proceeds on missing data).
        let mut v = minimal_manifest_json();
        v.as_object_mut()
            .unwrap()
            .insert("cdot_jsonn".to_string(), serde_json::json!("typo.json"));
        let err = validate_loaded_manifest(&v)
            .expect_err("a manifest with an unknown field must be rejected");
        let msg = format!("{err}");
        assert!(
            msg.contains("does not match the expected schema"),
            "actionable schema message: {msg}"
        );
        assert!(
            msg.contains("ferro prepare"),
            "message points to the remedy: {msg}"
        );
    }

    #[test]
    fn validate_loaded_manifest_accepts_pre_versioning_manifest() {
        let mut v = minimal_manifest_json();
        v.as_object_mut().unwrap().remove("manifest_schema_version");
        validate_loaded_manifest(&v).expect("a manifest without a schema version must load");
    }

    #[test]
    fn validate_loaded_manifest_rejects_newer_schema_version() {
        let mut v = minimal_manifest_json();
        v["manifest_schema_version"] =
            serde_json::json!(CURRENT_MANIFEST_SCHEMA_VERSION as u64 + 1);
        let err = validate_loaded_manifest(&v).expect_err("a newer schema must be refused");
        assert!(format!("{err}").contains("newer than this build"));
    }

    #[test]
    fn validate_loaded_manifest_rejects_missing_required_field() {
        let mut v = minimal_manifest_json();
        v.as_object_mut().unwrap().remove("transcript_count");
        let err = validate_loaded_manifest(&v)
            .expect_err("a missing required field must be a hard error, not a silent None");
        assert!(format!("{err}").contains("does not match the expected schema"));
    }

    #[test]
    fn validate_loaded_manifest_rejects_wrong_typed_field() {
        let mut v = minimal_manifest_json();
        // transcript_fastas as a string instead of an array — the exact class the
        // untyped runtime loader would otherwise silently read as `None`.
        v["transcript_fastas"] = serde_json::json!("not-an-array");
        let err = validate_loaded_manifest(&v).expect_err("a wrong-typed field must be rejected");
        assert!(format!("{err}").contains("does not match the expected schema"));
    }

    #[test]
    fn load_or_default_reports_version_message_for_newer_schema_with_new_field() {
        // A reference prepared by a *newer* ferro carries both a higher schema
        // version and fields this build doesn't model. The version gate must win,
        // so the user sees an actionable "upgrade ferro" message rather than an
        // opaque serde "unknown field" error.
        let dir = tempfile::tempdir().unwrap();
        let mut v = minimal_manifest_json();
        let obj = v.as_object_mut().unwrap();
        obj.insert(
            "manifest_schema_version".to_string(),
            serde_json::json!(CURRENT_MANIFEST_SCHEMA_VERSION as u64 + 1),
        );
        obj.insert("a_future_field".to_string(), serde_json::json!("x"));
        std::fs::write(
            dir.path().join("manifest.json"),
            serde_json::to_string_pretty(&v).unwrap(),
        )
        .unwrap();

        let err = ReferenceManifest::load_or_default(dir.path())
            .expect_err("a newer-schema manifest must be refused");
        let msg = format!("{err}");
        assert!(
            msg.contains("newer than this build"),
            "version gate must win over the unknown-field error: {msg}"
        );
    }

    #[test]
    fn load_or_default_rejects_duplicate_json_keys() {
        // A literal duplicate JSON key is a malformed manifest. `load_or_default`
        // runs its authoritative typed deserialize over the ORIGINAL BYTES, so
        // serde's duplicate-field rejection fires. (A `serde_json::Value`
        // round-trip would instead silently collapse the duplicate, last-wins —
        // the regression this guards against.) Fail loud on a misread reference
        // (#1001).
        let dir = tempfile::tempdir().unwrap();
        // Minimal schema-valid manifest, but with `transcript_count` given twice.
        let raw = format!(
            r#"{{"prepared_at":"2024-01-01T00:00:00Z","transcript_fastas":[],"genome_fasta":null,"cdot_json":null,"transcript_count":0,"transcript_count":1,"available_prefixes":[],"manifest_schema_version":{}}}"#,
            CURRENT_MANIFEST_SCHEMA_VERSION
        );
        std::fs::write(dir.path().join("manifest.json"), raw).unwrap();

        let err = ReferenceManifest::load_or_default(dir.path())
            .expect_err("a manifest with a duplicate JSON key must be rejected");
        assert!(
            format!("{err}").contains("Failed to parse manifest"),
            "a duplicate key must surface as a parse error: {err}"
        );
    }

    #[test]
    fn reference_identity_field_round_trips() {
        let m = ReferenceManifest {
            reference_identity: Some("aa8b3246d83055cc".to_string()),
            ..Default::default()
        };
        let json = serde_json::to_value(&m).unwrap();
        assert_eq!(json["reference_identity"], "aa8b3246d83055cc");
        let back: ReferenceManifest = serde_json::from_value(json).unwrap();
        assert_eq!(back.reference_identity.as_deref(), Some("aa8b3246d83055cc"));
    }

    #[test]
    fn current_schema_version_is_three() {
        assert_eq!(CURRENT_MANIFEST_SCHEMA_VERSION, 3);
    }

    #[test]
    fn cdot_data_version_field_round_trips() {
        let m = ReferenceManifest {
            cdot_data_version: Some("data_v0.2.32".to_string()),
            ..Default::default()
        };
        let json = serde_json::to_value(&m).unwrap();
        assert_eq!(json["cdot_data_version"], "data_v0.2.32");
        let back: ReferenceManifest = serde_json::from_value(json).unwrap();
        assert_eq!(back.cdot_data_version.as_deref(), Some("data_v0.2.32"));
    }

    #[test]
    fn cdot_data_version_is_excluded_from_the_identity() {
        let base = serde_json::json!({
            "transcript_count": 3,
            "cdot_json": "cdot-0.2.32.refseq.GRCh38.json",
            "available_prefixes": ["NM_"],
        });
        let id = crate::prepare::identity::reference_identity_from_manifest(&base);
        let mut with_version = base.clone();
        with_version["cdot_data_version"] = serde_json::json!("data_v0.2.32");
        with_version["manifest_schema_version"] = serde_json::json!(3);
        assert_eq!(
            crate::prepare::identity::reference_identity_from_manifest(&with_version),
            id,
            "cdot_data_version and schema version must not enter the identity hash"
        );
    }

    #[test]
    fn save_stamps_a_reference_identity() {
        use tempfile::TempDir;
        let dir = TempDir::new().unwrap();
        let mut m = ReferenceManifest {
            reference_dir: dir.path().to_path_buf(),
            transcript_count: 1,
            available_prefixes: vec!["NM_".to_string()],
            ..Default::default()
        };
        m.save().unwrap();
        let saved: serde_json::Value =
            serde_json::from_slice(&std::fs::read(dir.path().join("manifest.json")).unwrap())
                .unwrap();
        let id = saved["reference_identity"]
            .as_str()
            .expect("identity stamped");
        assert_eq!(id.len(), 16); // 16-hex FNV-1a digest
                                  // Recompute matches what was written (excludes reference_identity itself).
        assert_eq!(
            crate::prepare::identity::reference_identity(dir.path(), &saved),
            id
        );
    }

    #[test]
    fn save_rejects_a_wired_but_unreadable_stamped_artifact() {
        use tempfile::TempDir;
        let dir = TempDir::new().unwrap();
        // Wire a stamped artifact that does not exist on disk. Previously `save()`
        // would silently drop its stamp (via `inject_content_stamps`) and mint an
        // identity the loader then hard-rejects; now it must fail at stamp time.
        let mut m = ReferenceManifest {
            reference_dir: dir.path().to_path_buf(),
            transcript_count: 1,
            derived_transcript_placements: Some(dir.path().join("missing.json")),
            ..Default::default()
        };
        let err = m.save().unwrap_err();
        assert!(
            matches!(err, FerroError::Io { .. }),
            "save() must reject a wired-but-unreadable artifact with an Io error, got {err:?}"
        );
        // And it must not have written a manifest with an inconsistent stamp.
        assert!(!dir.path().join("manifest.json").exists());
    }
}
