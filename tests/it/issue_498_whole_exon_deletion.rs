//! Real-world whole-exon-deletion protein prediction (#498), gated on
//! FERRO_MANIFEST (skips in manifest-less CI, like the other axis tests).
//!
//! Three cases:
//! - DMD (NM_004006.2) exon 10: in-frame whole-exon deletion (c.961_1149,
//!   189 nt = 63 codons) via the bare-`NM_` direct path. Oracle: HGVS spec
//!   protein/deletion.md "one or more exons" `p.(His321_Glu383del)`.
//! - DMD (NM_004006.2) exons 10–11: frameshift whole-exon deletion
//!   (c.961_1331, 371 nt) via the direct path. Spec exons-10–11 case
//!   (`p.(His321Leufs*3)` on NP_003997.2); ferro emits the version-`.1` value
//!   with matching anchor + stop distance.
//! - VSIR (NM_022153.2): a reverse-strand whole-exon deletion that currently
//!   DECLINES, blocked by the reverse-strand normalization bug #762 — a
//!   guard that ferro does not emit a *wrong* protein from mis-normalized
//!   endpoints. Predicts the `cases.json` oracle once #762 lands.

use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::{MultiFastaProvider, ReferenceProvider, VariantProjection, VariantProjector};
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

// ---------------------------------------------------------------------------
// Manifest-gated provider — mirrors the construction in mutalyzer_normalize_tests.rs
// ---------------------------------------------------------------------------

fn manifest_path() -> Option<PathBuf> {
    // FERRO_MANIFEST, when set, is authoritative — no fallback to well-known
    // paths. This lets CI explicitly disable the runner via
    // `FERRO_MANIFEST=/nonexistent` even on a host that happens to have one of
    // the well-known paths mounted.
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

/// `Arc<MultiFastaProvider>` newtype that implements `ReferenceProvider + Clone`
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
    fn resolve_legacy_gene_selector(&self, selector: &str) -> Option<String> {
        self.0.resolve_legacy_gene_selector(selector)
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

// ---------------------------------------------------------------------------
// Helper
// ---------------------------------------------------------------------------

fn protein_of(vp: &VariantProjector<ArcProvider>, input: &str, tx: &str) -> Option<String> {
    let v = ferro_hgvs::parse_hgvs(input).expect("parse");
    let proj: VariantProjection = vp.project_variant(&v, tx).expect("project");
    proj.protein.map(|p| p.to_string())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// VSIR: internal whole-exon deletion on a reverse-strand gene — currently
/// DECLINES, blocked by the reverse-strand coding-normalization bug (#762).
///
/// `NG_008835.1(NM_022153.2):c.677-21_704+62del` removes coding exon 7
/// (c.677_704, 28 nt → frameshift); mutalyzer predicts `p.(Arg226ProfsTer102)`
/// (the `cases.json` oracle). But ferro's reverse-strand normalization
/// mis-renders the interval as the crossed form `c.704-18_677+65del`, so the
/// protein gate is handed non-canonical endpoints and `whole_exon_deletion_span`
/// declines — rather than predict a *wrong* protein from the bad positions
/// (it would otherwise emit `p.(Ala227GlnfsTer6)`). This is the deliberate safe
/// behavior; it stays `None` until #762 corrects the normalization, at which
/// point the canonical endpoints `(677, -18)/(704, +65)` predict
/// `p.(Arg226ProfsTer102)`. This test guards against silently emitting a wrong
/// protein here.
#[test]
fn vsir_whole_exon_deletion_declines_pending_762() {
    let Some(vp) = manifest_projector() else {
        eprintln!("vsir_whole_exon_deletion_declines_pending_762: skipping — no manifest");
        return;
    };
    let got = protein_of(
        &vp,
        "NG_008835.1(NM_022153.2):c.677-21_704+62del",
        "NM_022153.2",
    );
    assert_eq!(
        got, None,
        "VSIR must decline (reverse-strand #762 mis-normalization) rather than \
         emit a wrong protein; unblocks once #762 lands"
    );
}

/// DMD exon 10 whole-exon in-frame deletion.
///
/// `NM_004006.2:c.961-1_1149+3del` removes exon 10 (c.961_1149, 189 nt =
/// 63 codons) → in-frame protein deletion.
///
/// Expected protein: HGVS spec protein/deletion.md oracle
/// `p.(His321_Glu383del)` on the DMD protein accession.
#[test]
fn dmd_exon10_whole_exon_deletion_inframe() {
    let Some(vp) = manifest_projector() else {
        eprintln!("dmd_exon10_whole_exon_deletion_inframe: skipping — no manifest");
        return;
    };
    let got = protein_of(&vp, "NM_004006.2:c.961-1_1149+3del", "NM_004006.2");
    // The exact NP_ accession and version are captured from ferro's live output.
    // The consequence `p.(His321_Glu383del)` is pinned to the HGVS spec oracle.
    assert_eq!(
        got.as_deref(),
        Some("NP_003997.1:p.(His321_Glu383del)"),
        "DMD exon 10 whole-exon in-frame deletion must match spec oracle"
    );
}

/// DMD exons 10–11 whole-exon deletion (forward-path frameshift).
///
/// `NM_004006.2:c.961-1_1331+1del` removes exons 10–11 (CDS c.961_1331, 371 nt,
/// 371 mod 3 ≠ 0) via the clean bare-`NM_` direct path → a frameshift starting
/// at codon 321. This is the spec's "deletion of exons 10 to 11" case
/// (protein/deletion.md gives `p.(His321Leufs*3)` on `NP_003997.2`). ferro
/// translates the version present in the manifest (`NP_003997.1`) and yields
/// `p.(His321PhefsTer3)`: same anchor (His321) and same stop distance (`fsTer3`)
/// as the spec; the substituted residue differs (Phe vs Leu) because the new
/// codon 321 is read from the `.1` CDS rather than the spec's `.2`. The value is
/// captured from ferro's live output; the anchor + frameshift structure match
/// the spec oracle. Complements the in-frame case by exercising the frameshift
/// branch of the predictor end-to-end on a real transcript.
#[test]
fn dmd_exons10_11_whole_exon_deletion_frameshift() {
    let Some(vp) = manifest_projector() else {
        eprintln!("dmd_exons10_11_whole_exon_deletion_frameshift: skipping — no manifest");
        return;
    };
    let got = protein_of(&vp, "NM_004006.2:c.961-1_1331+1del", "NM_004006.2");
    assert_eq!(
        got.as_deref(),
        Some("NP_003997.1:p.(His321PhefsTer3)"),
        "DMD exons 10–11 whole-exon frameshift (His321 anchor, fsTer3) per spec"
    );
}
