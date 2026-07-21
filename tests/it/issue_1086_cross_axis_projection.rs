//! Issue #1086: cross-axis projection defects.
//!
//! **D5** — the `r.` axis silently discarded intronic offsets, mapping an
//! intronic position onto a *different, existing* exonic RNA base
//! (`c.3675-45` → `r.3676`). Per the HGVS spec, an RNA reference sequence
//! contains no intron sequence and therefore cannot be used to describe a
//! variant affecting intronic bases (`background/numbering.md`, "RNA
//! reference sequences … do **not contain** intron sequences and can
//! therefore **not be used** to describe variants affecting these
//! sequences"). The correct behaviour is to decline the `r.` axis rather
//! than fabricate an exonic coordinate.
//!
//! Reference-backed: skips unless `FERRO_MANIFEST` points at a prepared
//! reference (the same manifest gate the other axis tests use).

use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::{MultiFastaProvider, ReferenceProvider, VariantProjection, VariantProjector};
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

// ---------------------------------------------------------------------------
// Manifest-gated provider — mirrors issue_498_whole_exon_deletion.rs
// ---------------------------------------------------------------------------

fn manifest_path() -> Option<PathBuf> {
    // FERRO_MANIFEST, when set, is authoritative — no fallback to well-known
    // paths, so CI can disable the runner with `FERRO_MANIFEST=/nonexistent`.
    if let Ok(path) = std::env::var("FERRO_MANIFEST") {
        let p = PathBuf::from(path);
        return if p.exists() { Some(p) } else { None };
    }
    let p = Path::new("benchmark-output/manifest.json");
    if p.exists() {
        return Some(p.to_path_buf());
    }
    None
}

fn provider() -> Option<Arc<MultiFastaProvider>> {
    static PROVIDER: OnceLock<Option<Arc<MultiFastaProvider>>> = OnceLock::new();
    PROVIDER
        .get_or_init(|| {
            let path = manifest_path()?;
            Some(Arc::new(
                MultiFastaProvider::from_manifest(&path)
                    .unwrap_or_else(|e| panic!("from_manifest({}) failed: {e}", path.display())),
            ))
        })
        .clone()
}

/// `Arc<MultiFastaProvider>` newtype implementing `ReferenceProvider + Clone`
/// so it can be plugged into [`VariantProjector`], which requires `Clone`.
#[derive(Clone)]
struct ArcProvider(Arc<MultiFastaProvider>);

impl ReferenceProvider for ArcProvider {
    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, ferro_hgvs::FerroError> {
        self.0.get_transcript(id)
    }
    fn get_sequence(
        &self,
        id: &str,
        start: u64,
        end: u64,
    ) -> Result<String, ferro_hgvs::FerroError> {
        self.0.get_sequence(id, start, end)
    }
    fn has_transcript(&self, id: &str) -> bool {
        self.0.has_transcript(id)
    }
    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, ferro_hgvs::FerroError> {
        self.0.get_genomic_sequence(contig, start, end)
    }
    fn has_genomic_data(&self) -> bool {
        self.0.has_genomic_data()
    }
    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, ferro_hgvs::FerroError> {
        self.0.get_protein_sequence(accession, start, end)
    }
    fn has_protein_data(&self) -> bool {
        self.0.has_protein_data()
    }
    fn genomic_placement(
        &self,
        parent: &ferro_hgvs::hgvs::variant::Accession,
    ) -> Option<ferro_hgvs::reference::GenomicPlacement> {
        self.0.genomic_placement(parent)
    }
    fn get_transcript_for_variant(
        &self,
        variant: &ferro_hgvs::hgvs::variant::HgvsVariant,
    ) -> Result<Arc<Transcript>, ferro_hgvs::FerroError> {
        self.0.get_transcript_for_variant(variant)
    }
    fn get_transcript_for_accession(
        &self,
        accession: &ferro_hgvs::hgvs::variant::Accession,
    ) -> Result<Arc<Transcript>, ferro_hgvs::FerroError> {
        self.0.get_transcript_for_accession(accession)
    }
    fn has_transcript_version_exact(&self, id: &str) -> bool {
        self.0.has_transcript_version_exact(id)
    }
    fn genomic_placement_on_build(
        &self,
        parent: &ferro_hgvs::hgvs::variant::Accession,
        build: Option<&str>,
    ) -> Option<ferro_hgvs::reference::GenomicPlacement> {
        self.0.genomic_placement_on_build(parent, build)
    }
    fn resolve_legacy_gene_selector(
        &self,
        selector: &str,
        ng_parent: Option<&ferro_hgvs::hgvs::variant::Accession>,
    ) -> Option<String> {
        self.0.resolve_legacy_gene_selector(selector, ng_parent)
    }
    fn sole_hosted_transcript(
        &self,
        ng_parent: &ferro_hgvs::hgvs::variant::Accession,
    ) -> Option<String> {
        self.0.sole_hosted_transcript(ng_parent)
    }
    fn get_protein_length(&self, accession: &str) -> Result<u64, ferro_hgvs::FerroError> {
        self.0.get_protein_length(accession)
    }
    fn get_sequence_length(&self, id: &str) -> Result<u64, ferro_hgvs::FerroError> {
        self.0.get_sequence_length(id)
    }
}

fn manifest_projector() -> Option<VariantProjector<ArcProvider>> {
    let p = provider()?;
    let cdot = p.cdot_mapper()?.clone();
    let projector = Projector::new(cdot);
    Some(VariantProjector::new(projector, ArcProvider(p)))
}

/// Project `input` onto `tx` and hand back the whole projection.
fn project(vp: &VariantProjector<ArcProvider>, input: &str, tx: &str) -> VariantProjection {
    let v = ferro_hgvs::parse_hgvs(input).expect("parse");
    vp.project_variant(&v, tx).expect("project")
}

// ---------------------------------------------------------------------------
// D5 — the r. axis must never silently map an intronic offset onto a plain
//      exonic r. position.
// ---------------------------------------------------------------------------

/// `NM_000552.3:c.3675-45_3692del` starts 45 nt inside an intron. That
/// nucleotide is absent from the mature RNA, so no plain exonic `r.` position
/// can name it: the projector must decline the `r.` axis. Before the fix it
/// emitted `NM_000552.3:r.(3676_3694del)` — a well-formed description of a
/// *different*, existing base.
#[test]
fn intronic_start_offset_declines_rna_axis() {
    let Some(vp) = manifest_projector() else {
        eprintln!("issue_1086: skipping — no manifest at FERRO_MANIFEST");
        return;
    };
    let proj = project(&vp, "NM_000552.3:c.3675-45_3692del", "NM_000552.3");
    assert!(
        proj.rna.is_none(),
        "intronic c.3675-45 must not project to a plain exonic r. position, got {:?}",
        proj.rna.map(|r| r.to_string()),
    );
    // The c. axis is correct here and must stay correct.
    assert_eq!(
        proj.coding.expect("coding axis").to_string(),
        "NM_000552.3:c.3675-44_3693del",
    );
}

/// `NM_004006.2:c.123-65_-50` likewise carries an intronic start offset; the
/// pre-fix output was the nonsensical `NM_004006.2:r.(123_-50)`.
#[test]
fn intronic_offset_utr_span_declines_rna_axis() {
    let Some(vp) = manifest_projector() else {
        eprintln!("issue_1086: skipping — no manifest at FERRO_MANIFEST");
        return;
    };
    let proj = project(&vp, "NM_004006.2:c.123-65_-50", "NM_004006.2");
    assert!(
        proj.rna.is_none(),
        "intronic c.123-65 must not project to a plain exonic r. position, got {:?}",
        proj.rna.map(|r| r.to_string()),
    );
    assert_eq!(
        proj.coding.expect("coding axis").to_string(),
        "NM_004006.2:c.123-65_-50",
    );
}

/// Guard against over-declining: a wholly exonic edit on the same transcript
/// still projects to `r.`, with the 3'-shift behaviour untouched.
#[test]
fn exonic_edit_still_projects_to_rna_axis() {
    let Some(vp) = manifest_projector() else {
        eprintln!("issue_1086: skipping — no manifest at FERRO_MANIFEST");
        return;
    };
    let proj = project(&vp, "NM_004006.2:c.897T>G", "NM_004006.2");
    assert_eq!(
        proj.rna.expect("exonic edits keep an r. axis").to_string(),
        "NM_004006.2:r.(897u>g)",
    );
}
