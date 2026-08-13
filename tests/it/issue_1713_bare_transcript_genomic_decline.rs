//! Issue #1713 — `ferro project --axis g` declines every bare `NM_`-rooted
//! input, and the reason it gave said so in a way that pointed at the wrong
//! thing.
//!
//! **The filed claim does not reproduce, and that is the finding.** #1713
//! reports the CLI declining an axis "the library demonstrably can produce for
//! the same inputs on the same reference". Measured against the prepared
//! GRCh38 reference, the library declines identically:
//!
//! ```text
//! NM_004006.2:c.100del
//!   CLI  project_axis(.., Genomic, ..)  -> Unavailable
//!   LIB  project_variant(..).genomic    -> None
//!   LIB  project_to_genomic(..)         -> Err(UnsupportedProjection:
//!            "... has no parent reference (genomic_context) ... (see #327)")
//! ```
//!
//! There is no CLI/library asymmetry: `select_axis` reads
//! `VariantProjection::genomic`, so the CLI reports exactly the value the
//! library computed. The decline is specified — #327 requires
//! `UnsupportedProjection` for a transcript-coordinate input with no resolvable
//! parent, and `project_to_genomic_nc` deliberately refuses to synthesize a
//! parent from cdot.
//!
//! What *was* defective is the reporting, which is the issue's own third
//! option. The projector articulated the reason and dropped it, so the CLI fell
//! back to a string synthesized from the axis code alone — "no g.
//! representation for this variant" — asserting a property of the **variant**.
//! That is false here in a way anyone can check: the alignment is present, and
//! naming the genomic reference in the description renders the axis
//! immediately. The wording is what made a specified decline look like a
//! capability gap, and it is what this module pins.
//!
//! Reference-backed: skips unless `FERRO_MANIFEST` points at a prepared
//! reference.

use ferro_hgvs::cli::project::{project_axis, Axis, AxisOutcome};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::{MultiFastaProvider, ReferenceProvider, VariantProjector};
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

fn manifest_path() -> Option<PathBuf> {
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

/// Every bare-`NM_` input whose genomic axis the CLI declines must be declined
/// by the library too, for the same variant, on the same projector.
///
/// This is the invariant #1713 alleged was broken. It is asserted in both
/// directions — a CLI decline demands a library `None`, and a CLI render demands
/// a library `Some` with the identical string — so a future divergence fails
/// here rather than being re-filed as a capability gap.
///
/// `NM_004006.2:c.1704del` is deliberately in the corpus and deliberately not
/// expected to decline: #1708 re-parents a junction-crossing output onto the
/// genomic reference the crossing already resolved, so this row *renders*
/// despite arriving bare. Pinning agreement rather than a verdict is what lets
/// one corpus carry both classes.
///
/// **Agreement alone is too weak for that row**, though, so each input carries
/// its own expected verdict as well. A generic "either both render or both
/// decline" match accepts `(Unavailable, None)` for every row — so a regression
/// of #1708's re-parenting would drop `c.1704del` into the declining arm and
/// pass, while this doc comment went on claiming the row pins the render. The
/// `Verdict` column is what makes the two halves of the claim inseparable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Verdict {
    /// The genomic axis must be derivable for this input.
    Renders,
    /// The genomic axis must decline for this input.
    Declines,
}

#[test]
fn cli_and_library_agree_on_the_genomic_axis_for_bare_transcript_inputs() {
    let Some(vp) = manifest_projector() else {
        eprintln!(
            "issue_1713: skipping — no manifest at FERRO_MANIFEST. This is NOT the gate: \
             the hermetic half runs in CI as \
             src/cli/project.rs::project_axis_bare_nm_genomic_decline_names_cause_and_remedy"
        );
        return;
    };

    let mut checked = 0usize;
    for (input, verdict) in [
        ("NM_004006.2:c.100del", Verdict::Declines),
        ("NM_000492.4:c.1521_1523del", Verdict::Declines),
        // #1708: junction-crossing, re-parented onto the reference the crossing
        // already resolved, so this one renders despite arriving bare.
        ("NM_004006.2:c.1704del", Verdict::Renders),
        ("NC_000023.11(NM_004006.2):c.100del", Verdict::Renders),
    ] {
        let v = ferro_hgvs::parse_hgvs(input).expect("parse");
        let tx = v
            .accession()
            .map(|a| a.transcript_accession())
            .expect("transcript accession");

        let cli = project_axis(&vp, &v, Axis::Genomic, None)
            .unwrap_or_else(|e| panic!("{input}: CLI hard error {e}"));
        let lib = vp
            .project_variant(&v, &tx)
            .unwrap_or_else(|e| panic!("{input}: library hard error {e}"))
            .genomic
            .map(|g| g.to_string());

        match (&cli, &lib) {
            (AxisOutcome::Rendered { output, .. }, Some(g)) => {
                assert_eq!(
                    verdict,
                    Verdict::Renders,
                    "{input}: expected the genomic axis to decline, but it rendered {output}"
                );
                assert_eq!(
                    output, g,
                    "{input}: CLI and library rendered different genomic axes"
                );
            }
            (AxisOutcome::Unavailable { reason, .. }, None) => assert_eq!(
                verdict,
                Verdict::Declines,
                "{input}: expected the genomic axis to render, but it declined: {reason}"
            ),
            _ => panic!("{input}: CLI says {cli:?} but the library says {lib:?}"),
        }
        checked += 1;
    }
    // A corpus that silently shrank to nothing would pass every assertion above.
    assert_eq!(checked, 4, "the corpus must not shrink silently");
}

/// The decline a bare-`NM_` input receives must name its own cause and the
/// remedy, on the real reference.
///
/// The hermetic half of this lives in `src/cli/project.rs`
/// (`project_axis_bare_nm_genomic_decline_names_cause_and_remedy`) and runs in
/// CI. This half is what makes the wording checkable against *real* cdot: the
/// worked example it prints has to name the chromosome this transcript is
/// actually aligned to, which a mock cannot demonstrate.
#[test]
fn the_bare_transcript_genomic_decline_names_its_cause_and_the_remedy() {
    let Some(vp) = manifest_projector() else {
        eprintln!(
            "issue_1713: skipping — no manifest at FERRO_MANIFEST. This is NOT the gate: \
             the hermetic half runs in CI as \
             src/cli/project.rs::project_axis_bare_nm_genomic_decline_names_cause_and_remedy"
        );
        return;
    };

    // `NM_004006.2` (DMD) is on chrX; `NM_000492.4` (CFTR) on chr7. Two
    // different chromosomes, so a hardcoded accession in the message cannot
    // satisfy both.
    for (input, expect_parent) in [
        ("NM_004006.2:c.100del", "NC_000023.11"),
        ("NM_000492.4:c.1521_1523del", "NC_000007.14"),
    ] {
        let v = ferro_hgvs::parse_hgvs(input).expect("parse");
        let tx = v
            .accession()
            .map(|a| a.transcript_accession())
            .expect("transcript accession");

        let outcome = project_axis(&vp, &v, Axis::Genomic, None).expect("no hard error");
        let AxisOutcome::Unavailable { reason, .. } = outcome else {
            panic!("{input}: expected a decline, got {outcome:?}");
        };
        assert!(
            !reason.contains("no g. representation for this variant"),
            "{input}: the synthesized axis-code string states a property of the \
             variant, which is false here: {reason}"
        );
        assert!(
            reason.contains(&tx),
            "{input}: the decline must name the accession: {reason}"
        );
        assert!(
            reason.contains(&format!("{expect_parent}({tx})")),
            "{input}: the decline must show the remedy on this transcript's own \
             aligned genomic reference: {reason}"
        );

        // And the remedy has to actually work — otherwise the message is a
        // well-formed lie. Rebuild the description the decline recommends and
        // project it; the axis must render on the named reference.
        let suffix = input.split_once(':').expect("accession:axis").1;
        let anchored = format!("{expect_parent}({tx}):{suffix}");
        let av = ferro_hgvs::parse_hgvs(&anchored)
            .unwrap_or_else(|e| panic!("the recommended form must parse: {anchored}: {e}"));
        let outcome = project_axis(&vp, &av, Axis::Genomic, None).expect("no hard error");
        match outcome {
            AxisOutcome::Rendered { output, .. } => assert!(
                output.starts_with(&format!("{expect_parent}:g.")),
                "{anchored}: expected a g. axis on {expect_parent}, got {output}"
            ),
            other => panic!("{anchored}: the recommended form must render, got {other:?}"),
        }
    }
}
