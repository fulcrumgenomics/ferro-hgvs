//! A `delins` that normalization decomposes into an in-cis allele of sub-edits
//! must still project to **one** protein consequence over the affected residue
//! range — never a bracketed protein allele.
//!
//! Spec basis (vendored `assets/hgvs-nomenclature/docs/recommendations/`):
//!
//! - `DNA/delins.md` — the coalescing rule exists to stop tools "making
//!   conflicting and incorrect predictions of two different substitutions at one
//!   position (e.g., `c.235_237delinsTAT` (`p.Lys79Tyr`) versus
//!   `c.[235A>T;237G>T]` (`p.[Lys79*;Lys79Asn]`)".
//! - `protein/substitution.md` — "the description `p.Arg76_Cys77delinsSerTrp` is
//!   correct, the description `p.[Arg76Ser;Cys77Trp]` is not correct".
//! - `DNA/alleles.md` — an allele must not mark a position both changed and
//!   unchanged (`c.[2376G>C;3103=]` is "not correct"), which rules out a
//!   bracketed protein allele carrying an `(=)` member.
//!
//! Describing each decomposed member independently is not merely mis-rendered:
//! each member is predicted against the *reference* codon, so the reported amino
//! acids can be wrong outright (the first case below yielded `Leu` where the
//! combined codon encodes `Ile`).
//!
//! Ground truth is Mutalyzer 3.1.1 on the same transcript. (Mutalyzer prints the
//! bare `p.(...)`; ferro prefixes the protein accession and, per #310, never
//! emits a `(GENE)` selector on a `p.` variant.)
//!
//! Reference-backed: skips unless `FERRO_MANIFEST` points at a prepared
//! reference (same gate as the other real-reference projection tests).

use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::multi_fasta::MultiFastaProvider;
use ferro_hgvs::VariantProjector;

/// WT1 — the transcript the four reported reproductions are on.
const WT1_TX: &str = "NM_024426.4";

fn manifest_path() -> Option<PathBuf> {
    // FERRO_MANIFEST, when set, is authoritative — no fallback to the
    // repo-relative default, so a deliberate `FERRO_MANIFEST=/nonexistent`
    // skips rather than silently picking up a stray local reference.
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

/// A projector over the manifest-loaded reference. `Arc<P>` has a blanket
/// `ReferenceProvider` impl and is `Clone`, so it plugs straight into
/// `VariantProjector` (same wiring as the `ferro project` CLI).
fn manifest_projector() -> Option<VariantProjector<Arc<MultiFastaProvider>>> {
    let provider = provider()?;
    let cdot = provider.cdot_mapper()?.clone();
    Some(VariantProjector::new(Projector::new(cdot), provider))
}

/// The projected protein consequence of `hgvs`, rendered.
fn project_protein(
    projector: &VariantProjector<Arc<MultiFastaProvider>>,
    hgvs: &str,
) -> Option<String> {
    let parsed = parse_hgvs(hgvs).unwrap_or_else(|e| panic!("parse {hgvs}: {e}"));
    let projection = projector
        .project_variant(&parsed, WT1_TX)
        .unwrap_or_else(|e| panic!("project {hgvs}: {e}"));
    projection.protein.map(|p| p.to_string())
}

/// Every reported `delins` must project to a single protein consequence matching
/// Mutalyzer 3.1.1 — and, whatever the rendering, never to a bracketed allele.
#[test]
fn decomposed_delins_projects_to_one_combined_protein_consequence() {
    let Some(projector) = manifest_projector() else {
        eprintln!("protein_cis_delins_decomposition: skipping — no manifest at FERRO_MANIFEST");
        return;
    };

    // (input, expected protein consequence per Mutalyzer 3.1.1)
    let cases = [
        // Two members at one residue (Phe139) plus an adjacent one: the exact
        // shape delins.md forbids. Per-member prediction also reported `Leu`
        // where the combined codon encodes `Ile`.
        (
            "NM_024426.4:c.415_419delinsATATG",
            "NP_077744.3:p.(Phe139_Ile140delinsIleCys)",
        ),
        // One member is silent: the unchanged residue is context, not a bracket
        // member (alleles.md), so only Lys141 is described — as a substitution.
        (
            "NM_024426.4:c.418_422delinsATTAC",
            "NP_077744.3:p.(Lys141Thr)",
        ),
        // Overlapping per-member ranges (both name Glu143) must collapse into one
        // delins across the whole affected range.
        (
            "NM_024426.4:c.426_430delinsCCATG",
            "NP_077744.3:p.(Gln142_Pro144delinsHisHisAla)",
        ),
        // A whole-protein `(=)` member must not be bracketed alongside a changed
        // one (alleles.md: no "changed and unchanged" allele).
        (
            "NM_024426.4:c.429_434delinsTACCTC",
            "NP_077744.3:p.(Glu143_Pro144delinsAspThr)",
        ),
    ];

    for (input, expected) in cases {
        let actual = project_protein(&projector, input)
            .unwrap_or_else(|| panic!("{input}: no protein consequence projected"));
        assert!(
            !actual.contains('['),
            "{input}: a decomposed delins must not project to a bracketed protein \
             allele (substitution.md, delins.md), got {actual}"
        );
        assert_eq!(actual, expected, "{input}: protein consequence");
    }
}
