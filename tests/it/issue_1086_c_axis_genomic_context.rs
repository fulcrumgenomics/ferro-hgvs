//! Issue #1086 (defect D3) — `ferro project --axis c` dropped the genomic
//! context (`NC_`/`NG_` parent) the input supplied, emitting a bare
//! `NM_004006.2:c.357+1G>A`. That string is the HGVS spec's own literal
//! **"not correct"** example.
//!
//! `assets/hgvs-nomenclature/docs/background/refseq.md:47-49`:
//!
//! > intronic sequences are considered to be **within the boundaries** of a
//! > transcript reference sequence and may only be used to describe a variant
//! > when a genomic reference sequence identifier is provided, e.g.,
//! > `NC_000023.11(NM_004006.2):c.93+1G>T` … **not correct:**
//! > `NM_004006.2:c.357+1G>A` … **correct:**
//! > `NC_000023.10(NM_004006.2):c.357+1G>A`,
//! > `NG_012232.1(NM_004006.2):c.357+1G>A`, `LRG_199t1:c.357+1G>A`.
//!
//! Restated for coding transcripts at `refseq.md:136` and for non-coding
//! transcripts at `refseq.md:157`; the parenthesised rendering itself is
//! specified at `refseq.md:138-140`:
//!
//! > when, based on a genomic reference sequence, variants are reported using a
//! > `c.` prefix, the transcript variant used should be indicated. … for `NC_`
//! > or `NG_` reference sequences, the annotated transcript used is given in
//! > parentheses directly following the accession.version number.
//!
//! ## The rule pinned here
//!
//! **The `c.`/`n.` axis retains the genomic context whenever the *input*
//! supplied one, and stays bare whenever it did not** — i.e. the projector
//! preserves the caller's reference framing rather than classifying the
//! position.
//!
//! Why this rule and not "wrap only intronic/flanking positions":
//!
//! 1. It is spec-sufficient: it is a superset of the mandatory case, so the
//!    "not correct" bare-intronic string can never be emitted from a parented
//!    input, and the only way to reach a bare intronic `c.` is a bare input,
//!    which was already spec-invalid before ferro touched it.
//! 2. It is spec-endorsed for exonic positions too: `refseq.md:138-140` frames
//!    the parenthesised form as what to do "when, based on a genomic reference
//!    sequence, variants are reported using a `c.` prefix" — no intronic
//!    qualifier. A `c.` description derived from `NC_000023.11(NM_004006.2)`
//!    *is* based on a genomic reference sequence whether or not the position
//!    happens to fall in an exon.
//! 3. It matches what ferro already does everywhere else: `ferro normalize`
//!    keeps the context, and `project --axis r` keeps it for exonic positions
//!    too (`NC_000023.11(NM_004006.2):r.(79a>g)`). Only `--axis c` differed.
//! 4. It is positional-classification-free, so it cannot drift out of sync with
//!    the intronic/UTR flag plumbing.
//!
//! `exonic_input_without_context_stays_bare` and
//! `exonic_input_with_context_keeps_context` pin both halves of the rule so a
//! future change to either direction is a deliberate one.

use ferro_hgvs::cli::project::{project_axis, Axis, AxisOutcome};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::{
    parse_hgvs, HgvsVariant, MultiFastaProvider, ReferenceProvider, VariantProjector,
};
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

/// `Arc<MultiFastaProvider>` newtype implementing `ReferenceProvider + Clone` so
/// it can be plugged into `VariantProjector` (mirrors
/// `tests/it/issue_843_allele_rna_build_scoping.rs`).
#[derive(Clone)]
struct ArcProvider(Arc<MultiFastaProvider>);

impl ReferenceProvider for ArcProvider {
    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, ferro_hgvs::FerroError> {
        self.0.get_transcript(id)
    }
    fn get_sequence(&self, id: &str, s: u64, e: u64) -> Result<String, ferro_hgvs::FerroError> {
        self.0.get_sequence(id, s, e)
    }
    fn has_transcript(&self, id: &str) -> bool {
        self.0.has_transcript(id)
    }
    fn get_genomic_sequence(
        &self,
        contig: &str,
        s: u64,
        e: u64,
    ) -> Result<String, ferro_hgvs::FerroError> {
        self.0.get_genomic_sequence(contig, s, e)
    }
    fn has_genomic_data(&self) -> bool {
        self.0.has_genomic_data()
    }
    fn get_protein_sequence(
        &self,
        accession: &str,
        s: u64,
        e: u64,
    ) -> Result<String, ferro_hgvs::FerroError> {
        self.0.get_protein_sequence(accession, s, e)
    }
    fn has_protein_data(&self) -> bool {
        self.0.has_protein_data()
    }
    fn get_transcript_for_variant(
        &self,
        variant: &HgvsVariant,
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
    fn infer_genome_build(
        &self,
        accession: &ferro_hgvs::hgvs::variant::Accession,
    ) -> Option<&'static str> {
        self.0.infer_genome_build(accession)
    }
    fn genomic_placement(
        &self,
        parent: &ferro_hgvs::hgvs::variant::Accession,
    ) -> Option<ferro_hgvs::reference::GenomicPlacement> {
        self.0.genomic_placement(parent)
    }
    fn genomic_placement_on_build(
        &self,
        parent: &ferro_hgvs::hgvs::variant::Accession,
        build: Option<&str>,
    ) -> Option<ferro_hgvs::reference::GenomicPlacement> {
        self.0.genomic_placement_on_build(parent, build)
    }
}

fn manifest_projector() -> Option<VariantProjector<ArcProvider>> {
    let p = provider()?;
    let cdot = p.cdot_mapper()?.clone();
    Some(VariantProjector::new(Projector::new(cdot), ArcProvider(p)))
}

/// Project `input` onto `axis` exactly as `ferro project --axis <axis>` does.
fn project(vp: &VariantProjector<ArcProvider>, input: &str, axis: Axis) -> AxisOutcome {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse({input}) failed: {e}"));
    project_axis(vp, &v, axis, None)
        .unwrap_or_else(|e| panic!("project_axis({input}, {}) failed: {e}", axis.code()))
}

/// Assert the `c.` axis renders exactly `expected`.
fn assert_coding(vp: &VariantProjector<ArcProvider>, input: &str, expected: &str) {
    match project(vp, input, Axis::Coding) {
        AxisOutcome::Rendered { output, .. } => assert_eq!(
            output, expected,
            "c. axis for {input}: expected {expected}, got {output}"
        ),
        AxisOutcome::Unavailable { reason, .. } => {
            panic!("c. axis for {input}: unexpectedly unavailable ({reason})")
        }
    }
}

macro_rules! skip_without_manifest {
    ($vp:ident) => {
        let Some($vp) = manifest_projector() else {
            eprintln!("issue_1086: skipping — no FERRO_MANIFEST / benchmark-output manifest");
            return;
        };
    };
}

/// The issue's headline reproduction: a two-endpoint deep-intronic deletion.
/// Both endpoints carry an intronic offset, so a bare `NM_` rendering is
/// spec-invalid (refseq.md:47-49).
#[test]
fn d3_intronic_deletion_keeps_genomic_context_on_c_axis() {
    skip_without_manifest!(vp);
    assert_coding(
        &vp,
        "NC_000023.11(NM_004006.2):c.4072-1234_5155-246del",
        "NC_000023.11(NM_004006.2):c.4072-1234_5155-246del",
    );
}

/// The spec's own `NC_`-parented intronic example (refseq.md:47).
#[test]
fn spec_nc_intronic_example_keeps_genomic_context_on_c_axis() {
    skip_without_manifest!(vp);
    assert_coding(
        &vp,
        "NC_000023.11(NM_004006.2):c.93+1G>T",
        "NC_000023.11(NM_004006.2):c.93+1G>T",
    );
}

/// The spec's own `NG_`-parented intronic example, listed verbatim under
/// **correct** at refseq.md:49.
#[test]
fn spec_ng_intronic_example_keeps_genomic_context_on_c_axis() {
    skip_without_manifest!(vp);
    assert_coding(
        &vp,
        "NG_012232.1(NM_004006.2):c.357+1G>A",
        "NG_012232.1(NM_004006.2):c.357+1G>A",
    );
}

/// The bare `NM_` intronic string the spec calls **not correct** must never be
/// produced from a parented input, on any of the reproductions.
#[test]
fn spec_not_correct_bare_intronic_string_is_never_emitted() {
    skip_without_manifest!(vp);
    for input in [
        "NC_000023.11(NM_004006.2):c.4072-1234_5155-246del",
        "NC_000023.11(NM_004006.2):c.93+1G>T",
        "NG_012232.1(NM_004006.2):c.357+1G>A",
    ] {
        let AxisOutcome::Rendered { output, .. } = project(&vp, input, Axis::Coding) else {
            panic!("c. axis for {input}: unexpectedly unavailable");
        };
        assert!(
            output.starts_with("NC_") || output.starts_with("NG_"),
            "c. axis for {input} emitted a bare transcript description ({output}); \
             refseq.md:47-49 lists that form under \"not correct\""
        );
    }
}

/// Rule half 1: an input that supplied **no** genomic context keeps rendering
/// bare. A purely exonic bare `NM_…:c.` description is spec-valid, and ferro
/// must not fabricate a parent the caller did not give (#121).
#[test]
fn exonic_input_without_context_stays_bare() {
    skip_without_manifest!(vp);
    assert_coding(&vp, "NM_004006.2:c.79A>G", "NM_004006.2:c.79A>G");
}

/// Rule half 2: an exonic input that **did** supply a genomic context keeps it,
/// matching `--axis r` (`NC_000023.11(NM_004006.2):r.(79a>g)`) and
/// `ferro normalize`. Pins the "preserve caller framing" rule against a future
/// narrowing to intronic-only.
#[test]
fn exonic_input_with_context_keeps_context() {
    skip_without_manifest!(vp);
    assert_coding(
        &vp,
        "NC_000023.11(NM_004006.2):c.79A>G",
        "NC_000023.11(NM_004006.2):c.79A>G",
    );
}

/// The `r.` axis was already correct — guard it against collateral damage.
/// (The `c.`/`r.` 3'-shift difference on this input is *required*, not a bug:
/// `RNA/deletion.md:20` excludes the exon/exon-junction exception to the 3'
/// rule for an RNA reference sequence. Deliberately not asserted here beyond
/// the framing.)
#[test]
fn rna_axis_framing_is_unchanged() {
    skip_without_manifest!(vp);
    let AxisOutcome::Rendered { output, .. } = project(
        &vp,
        "NC_000023.11(NM_004006.2):c.4072-1234_5155-246del",
        Axis::Rna,
    ) else {
        panic!("r. axis unexpectedly unavailable");
    };
    assert_eq!(output, "NC_000023.11(NM_004006.2):r.(4073_5156del)");
}
