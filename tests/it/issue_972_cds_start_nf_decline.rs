//! Issue #972 (Task 4): decline projection *to* the coding (`c.`) / protein
//! (`p.`) axis on a transcript whose 5' CDS is annotated incomplete (Ensembl
//! `cds_start_NF`) — no confirmed `ATG` initiation codon means `c.1`/`p.1` are
//! undefined, so HGVS does not recommend describing such a transcript on
//! either axis. The `n.` (and genomic) axes are unaffected: they don't depend
//! on knowing where translation starts.
//!
//! ## The `r.` (RNA) axis (review follow-up)
//!
//! `r.` is ALSO gated, not left valid. `src/project/rna.rs`'s module doc is
//! explicit: "Numbering is CDS-relative (identical to `c.`)" — for a coding
//! transcript, `predict_rna` derives the `r.` position straight from the
//! (internally-computed) `c.` value, so it inherits the exact same "no
//! confirmed ATG" problem `c.`/`p.` have. This is unlike `n.`, whose numbering
//! is transcript-native (position 1 = the first transcript base, no CDS/ATG
//! involved) and stays valid. All four "build coding" sites gate `rna`
//! alongside `coding`/`protein` now.
//!
//! ## Four guard sites, not one
//!
//! `HgvsVariant::Cds(`/`Tx(` output construction for the coding/protein/rna
//! axes happens in four independent places in `projector.rs`:
//!
//! 1. `project_single_inner` — the genome-pivot path (`g.` input).
//! 2. `project_coding_direct` — a bare `c.` input (no `genomic_context`).
//! 3. `project_noncoding_direct` — a bare `n.`/`r.` input.
//! 4. `polya_multiaxis_projection` / `terminus_multiaxis_projection` — two
//!    early-return short-circuits reached *from inside*
//!    `project_single_inner`, for a `c.*` poly-A input (#797) and a whole-arm
//!    `c.pter`/`c.qter` terminus input (#887) respectively. Both build a
//!    `coding: Some(...)` echo of their own, bypassing the outer function's
//!    guard entirely if not separately gated — this was the review-caught
//!    leak: a `c.*Ndel` (or `c.pterdel`) on a `cds_start_NF` transcript still
//!    emitted a naive `c.` before this fix.
//!
//! Four test groups below cover each site:
//!
//! - A deterministic `MockProvider` fixture (`fixture()`) exercises the
//!   primary guard site (site 1 above) with a **genomic** input — projecting
//!   a `g.` variant onto a synthetic `cds_start_NF` transcript and checking
//!   all five axes (`c.`/`p.`/`r.` decline, `n.`/`g.` still render). This runs
//!   unconditionally (no reference needed) and is the guardrail: a
//!   `cds_start_incomplete: false` sibling fixture proves the same inputs *do*
//!   predict c./p./r. normally, so the decline is specific to the flag, not a
//!   broken projector.
//!
//! - A `FERRO_MANIFEST`-gated test against the real Ensembl transcript
//!   `ENST00000011700` (VPS13D; genuinely `cds_start_NF` — chr1,
//!   NC_000001.11:12277121-…). This exercises the two *direct* (bare-accession)
//!   guard sites (sites 2/3 above) via a bare `c.`/`n.` input, and the CLI
//!   `project_axis`/`select_axis` layer.
//!
//!   NOTE: a **bare genomic** (`NC_000001.11:g.…`) input onto this transcript
//!   cannot currently be exercised end-to-end: `CdotMapper::
//!   transcript_genome_span[_on_build]` — consulted by the genome-pivot path
//!   before it ever reaches the #972 guard — returns `None` for every Ensembl
//!   (`ENST`) transcript on this prepared reference, because Ensembl cdot data
//!   is merged into a separate lazily-loaded sub-mapper that the genome-span
//!   side-table (built once from the primary/RefSeq map) never sees. That gap
//!   is pre-existing, unrelated to #972, and reproduces on `main` before this
//!   change too (verified: `project_variant` on a genomic input against
//!   `ENST00000011700` fails with `ReferenceNotFound` regardless of the
//!   `cds_start_incomplete` guard). The `MockProvider` fixture above covers the
//!   genomic-input pivot path that this real-transcript gap currently blocks.
//!
//! - A poly-A (`c.*`) `MockProvider` fixture (site 4a) reproduces the review's
//!   leak directly: a `c.*` endpoint whose genome pivot lands outside the
//!   cdot exon span short-circuits into `polya_multiaxis_projection` before
//!   the outer function's guard ever runs.
//!
//! - A whole-arm terminus (`c.pter`/`c.qter`) `MockProvider` fixture (site 4b)
//!   reproduces the same defect class in `terminus_multiaxis_projection`.

use ferro_hgvs::cli::project::{project_axis, Axis, AxisOutcome};
use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::{ReferenceProvider, Strand};
use ferro_hgvs::{parse_hgvs, MultiFastaProvider, VariantProjection, VariantProjector};
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

// =============================================================================
// MockProvider fixture — deterministic, no reference required.
// =============================================================================

/// Build a (Projector, MockProvider) pair for `NM_TEST.1`: CDS "ATGCGCTAA"
/// (Met-Arg-Stop) on chr1 plus-strand at genome 1000-1008 (1-based
/// inclusive), mirroring `tests/it/projection.rs::fixture` exactly except for
/// `cds_start_incomplete`, which this helper sets per `incomplete`.
fn fixture(incomplete: bool) -> (Projector, MockProvider) {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TEST.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: incomplete,
            gene_name: Some("TESTGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            // [genome_start(0-based), genome_end(0-based excl), tx_start, tx_end]
            exons: vec![[1000, 1009, 0, 9]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_TEST.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TESTGENE".to_string()),
        TxStrand::Plus,
        "ATGCGCTAA".to_string(),
        Some(1),
        Some(9),
        vec![Exon::new(1, 1, 9)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1008),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    let prefix = "N".repeat(999);
    let suffix = "N".repeat(100);
    provider.add_genomic_sequence("chr1", format!("{}ATGCGCTAA{}", prefix, suffix));
    (projector, provider)
}

/// `cds_start_incomplete: true` (the #972 case): projecting a genomic
/// substitution inside the CDS must decline `c.`/`p.` while `n.`/`g.` still
/// render.
#[test]
fn genomic_input_declines_coding_protein_but_not_noncoding_genomic() {
    let (projector, provider) = fixture(true);
    let vp = VariantProjector::new(projector, provider);
    let result: VariantProjection = vp
        .project("NC_000001.11:g.1003C>A", "NM_TEST.1")
        .expect("projection itself must succeed — only the coding/protein axes decline");

    assert!(
        result.coding.is_none(),
        "cds_start_NF transcript: c. must decline (no confirmed ATG), got {:?}",
        result.coding.map(|v| v.to_string())
    );
    assert!(
        result.protein.is_none(),
        "cds_start_NF transcript: p. must decline (no confirmed ATG), got {:?}",
        result.protein.map(|v| v.to_string())
    );
    assert!(
        result.rna.is_none(),
        "cds_start_NF transcript: r. must decline too — r. numbering is \
         CDS-relative (identical to c.) for a coding transcript, per \
         `project::rna`'s module doc, so it inherits the same missing-ATG \
         problem, got {:?}",
        result.rna.map(|v| v.to_string())
    );
    assert!(
        result.noncoding.is_some(),
        "n. does not depend on the initiation codon and must still render"
    );
    assert!(
        result.genomic.is_some(),
        "g. does not depend on the initiation codon and must still render"
    );
}

/// Guardrail contrast: the identical fixture with `cds_start_incomplete:
/// false` predicts c./p./r. normally (matches `tests/it/projection.rs::
/// end_to_end_missense`). Proves the decline above is specific to the flag,
/// not a broken/overly-aggressive projector.
#[test]
fn genomic_input_predicts_coding_protein_when_cds_start_is_complete() {
    let (projector, provider) = fixture(false);
    let vp = VariantProjector::new(projector, provider);
    let result: VariantProjection = vp
        .project("NC_000001.11:g.1003C>A", "NM_TEST.1")
        .expect("projection should succeed");

    let c = result.coding.as_ref().unwrap().to_string();
    assert!(c.contains(":c.4C>A"), "got c. = {c}");
    let p = result.protein.as_ref().unwrap().to_string();
    assert_eq!(p, "NP_TEST.1:p.(Arg2Ser)");
    assert!(
        result.rna.is_some(),
        "r. must still predict when cds_start is complete (guardrail — the \
         decline above is specific to the flag, not a broken r. predictor)"
    );
}

/// Same fixture, driven through the CLI `project_axis`/`select_axis` layer
/// (`ferro project --axis c|p|r|n|g`) rather than the raw library call, to lock
/// the user-facing contract: c./p./r. surface as `AxisOutcome::Unavailable`
/// (exit 0, not a hard error), n./g. as `AxisOutcome::Rendered`.
#[test]
fn cli_axis_c_and_p_are_unavailable_n_and_g_are_rendered() {
    let (projector, provider) = fixture(true);
    let vp = VariantProjector::new(projector, provider);
    let variant = parse_hgvs("NC_000001.11:g.1003C>A").expect("parse");

    for axis in [Axis::Coding, Axis::Protein, Axis::Rna] {
        let outcome =
            project_axis(&vp, &variant, axis, Some("NM_TEST.1")).expect("not a hard error");
        match outcome {
            AxisOutcome::Unavailable { transcript_id, .. } => {
                assert_eq!(transcript_id.as_deref(), Some("NM_TEST.1"));
            }
            other => panic!("axis {axis:?}: expected Unavailable, got {other:?}"),
        }
    }
    for axis in [Axis::Noncoding, Axis::Genomic] {
        let outcome =
            project_axis(&vp, &variant, axis, Some("NM_TEST.1")).expect("not a hard error");
        assert!(
            matches!(outcome, AxisOutcome::Rendered { .. }),
            "axis {axis:?}: expected Rendered, got {outcome:?}"
        );
    }
}

// =============================================================================
// Real-reference test — ENST00000011700 (VPS13D), genuinely `cds_start_NF`.
// =============================================================================

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

/// `Arc<MultiFastaProvider>` newtype implementing `ReferenceProvider + Clone`
/// so it can be plugged into `VariantProjector` (mirrors
/// `tests/it/issue_498_whole_exon_deletion.rs`).
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

const VPS13D_TX: &str = "ENST00000011700.10";

/// Build a bare `c.` substitution `{VPS13D_TX}:c.50<ref>><alt>` at an
/// arbitrary in-CDS position (cds_start=0, i.e. no 5'UTR at all — the
/// annotation could not place a true ATG, hallmark of cds_start_NF). The ref
/// base is read from the transcript sequence so the substitution is valid
/// regardless of reference-data drift.
fn vps13d_coding_variant(
    raw_provider: &MultiFastaProvider,
) -> ferro_hgvs::hgvs::variant::HgvsVariant {
    let base = raw_provider
        .get_sequence(VPS13D_TX, 49, 50)
        .expect("read VPS13D base at tx pos 50");
    let alt = if base == "A" { "T" } else { "A" };
    let hgvs = format!("{VPS13D_TX}:c.50{base}>{alt}");
    parse_hgvs(&hgvs).unwrap_or_else(|e| panic!("parse {hgvs}: {e}"))
}

/// A bare `c.` input on the real, genuinely-`cds_start_NF` VPS13D transcript
/// must decline `c.` (trivially — it *is* the input), `p.`, and `r.` (derived
/// from `c.` and thus CDS-relative too), while its derived `n.` form still
/// renders. Exercises `project_coding_direct`.
#[test]
fn vps13d_bare_coding_input_declines_protein_but_derives_noncoding() {
    let Some(raw_provider) = provider() else {
        eprintln!(
            "vps13d_bare_coding_input_declines_protein_but_derives_noncoding: skipping — no manifest"
        );
        return;
    };
    let Some(vp) = manifest_projector() else {
        eprintln!(
            "vps13d_bare_coding_input_declines_protein_but_derives_noncoding: skipping — no manifest"
        );
        return;
    };
    let variant = vps13d_coding_variant(&raw_provider);

    let projection = vp
        .project_variant(&variant, VPS13D_TX)
        .unwrap_or_else(|e| panic!("project {variant}: {e}"));

    assert!(
        projection.coding.is_none(),
        "VPS13D ({VPS13D_TX}) is cds_start_NF: c. must decline, got {:?}",
        projection.coding.map(|v| v.to_string())
    );
    assert!(
        projection.protein.is_none(),
        "VPS13D ({VPS13D_TX}) is cds_start_NF: p. must decline, got {:?}",
        projection.protein.map(|v| v.to_string())
    );
    assert!(
        projection.rna.is_none(),
        "VPS13D ({VPS13D_TX}) is cds_start_NF: r. must decline too (CDS-relative \
         numbering, same as c.), got {:?}",
        projection.rna.map(|v| v.to_string())
    );
    let n = projection
        .noncoding
        .as_ref()
        .unwrap_or_else(|| panic!("n. must still render for a cds_start_NF transcript"))
        .to_string();
    assert!(n.contains(":n.50"), "got n. = {n}");
}

/// The CLI axis layer on the same real transcript: `--axis c`/`--axis p`/
/// `--axis r` decline as `AxisOutcome::Unavailable` (exit 0), `--axis n`
/// renders.
#[test]
fn vps13d_cli_axis_c_and_p_unavailable_n_rendered() {
    let Some(raw_provider) = provider() else {
        eprintln!("vps13d_cli_axis_c_and_p_unavailable_n_rendered: skipping — no manifest");
        return;
    };
    let Some(vp) = manifest_projector() else {
        eprintln!("vps13d_cli_axis_c_and_p_unavailable_n_rendered: skipping — no manifest");
        return;
    };
    let variant = vps13d_coding_variant(&raw_provider);

    for axis in [Axis::Coding, Axis::Protein, Axis::Rna] {
        let outcome = project_axis(&vp, &variant, axis, None)
            .unwrap_or_else(|e| panic!("axis {axis:?} must not be a hard error: {e}"));
        match outcome {
            AxisOutcome::Unavailable {
                transcript_id,
                reason,
            } => {
                assert_eq!(transcript_id.as_deref(), Some(VPS13D_TX));
                assert!(!reason.is_empty());
            }
            other => panic!("axis {axis:?}: expected Unavailable, got {other:?}"),
        }
    }
    let outcome = project_axis(&vp, &variant, Axis::Noncoding, None)
        .unwrap_or_else(|e| panic!("axis n must not be a hard error: {e}"));
    assert!(
        matches!(outcome, AxisOutcome::Rendered { .. }),
        "axis n: expected Rendered, got {outcome:?}"
    );
}

// =============================================================================
// Poly-A multi-axis leak (review finding #1) — `polya_multiaxis_projection`'s
// early-return short-circuit inside `project_single_inner`.
// =============================================================================

/// Build a deterministic sequence by repeating `pat` up to `len` bases
/// (mirrors `tests/it/issue_868_genomic_nc_units.rs::cycle`).
fn cycle(pat: &str, len: usize) -> String {
    pat.chars().cycle().take(len).collect()
}

/// `cds_start_incomplete: incomplete` variant of `tests/it/
/// issue_868_genomic_nc_units.rs::polya_fixture` — identical CDS/exon/genome
/// layout (a coding transcript whose exon→genome map strips the
/// post-transcriptional poly-A tail), except this fixture is projected via
/// the **multi-axis** `VariantProjector::project`/`project_variant`, not the
/// single-axis `project_to_genomic` the #868 test uses. A `c.*4` endpoint
/// (tx 12, the exon's genome_end) lands outside the cdot exon span, so
/// `project_single_inner` short-circuits into `polya_multiaxis_projection`
/// (#797) — the site whose early return built `coding: Some(coding)`
/// unconditionally, bypassing the `cds_start_incomplete` guard (the review's
/// leak: a `c.*Ndel` on a `cds_start_NF` transcript still emitted a naive
/// `c.` before this fix).
fn polya_fixture(incomplete: bool) -> VariantProjector<MockProvider> {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_POLYA972.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: incomplete,
            gene_name: Some("POLYAGENE972".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            // Exon covers only the genomic core (poly-A tail stripped): tx 0..12.
            exons: vec![[1000, 1012, 0, 12]],
            cds_start: Some(3),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_POLYA972.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    // FASTA carries the full 15-base transcript incl. the 3-base poly-A tail;
    // the exon map above only spans tx 0..12, so c.*4..c.*6 (tx 12..15) is the
    // stripped poly-A region the fallback walks into.
    provider.add_transcript(Transcript::new(
        "NM_POLYA972.1".to_string(),
        Some("POLYAGENE972".to_string()),
        TxStrand::Plus,
        "TTTATGCGCGTAAAA".to_string(), // 15 bases, trailing "AAA" tail
        Some(4),                       // 1-based inclusive CDS start (tx index 3)
        Some(9),
        vec![Exon::with_genomic(1, 1, 12, 1000, 1011)],
        Some("NC_000001.11".to_string()),
        Some(1000),
        Some(1011),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    // Genomic backbone: 999 non-repeating leading bases so HGVS g.1000 == 0-based
    // index 999, then the 12-base exon core, then a non-repeating downstream tail.
    let genomic = format!(
        "{}{}{}",
        cycle("ACGT", 999),
        "TTTATGCGCGTA",
        cycle("CAGT", 100)
    );
    provider.add_genomic_sequence("NC_000001.11", genomic);
    VariantProjector::new(projector, provider)
}

/// RED→GREEN for the review's poly-A leak: before the fix,
/// `polya_multiaxis_projection`'s early return built `coding: Some(coding)`
/// unconditionally. `protein`/`rna`/`noncoding` are already `None` on this
/// path regardless (a poly-A position has no protein consequence and no
/// genomic exon context to derive `n.`/`r.` from), so only `coding` is
/// asserted here — it is the field the leak actually affected.
#[test]
fn polya_input_declines_coding_on_cds_start_nf_transcript() {
    let vp = polya_fixture(true);
    let projection = vp
        .project("NC_000001.11(NM_POLYA972.1):c.*4del", "NM_POLYA972.1")
        .expect("projection itself must succeed — only the coding axis declines");

    assert!(
        projection.coding.is_none(),
        "cds_start_NF transcript: c. must decline even via the poly-A \
         multi-axis short-circuit, got {:?}",
        projection.coding.map(|v| v.to_string())
    );
    assert!(
        projection.genomic.is_some(),
        "g. does not depend on the initiation codon and must still render"
    );
}

/// Guardrail: the identical fixture with `cds_start_incomplete: false` still
/// reports the `c.` echo on the poly-A short-circuit path — proves the
/// decline above is specific to the flag, not a broken poly-A fallback.
#[test]
fn polya_input_reports_coding_when_cds_start_is_complete() {
    let vp = polya_fixture(false);
    let projection = vp
        .project("NC_000001.11(NM_POLYA972.1):c.*4del", "NM_POLYA972.1")
        .expect("projection should succeed");

    let c = projection
        .coding
        .as_ref()
        .expect("c. must render when cds_start is complete")
        .to_string();
    assert!(c.contains(":c.*4del"), "got c. = {c}");
}

// =============================================================================
// Whole-arm terminus multi-axis leak (same defect class) —
// `terminus_multiaxis_projection`'s early-return short-circuit.
// =============================================================================

/// A projector with an NG_ parent (length `len`) and an `NM_TEST.1` cdot
/// record whose only field that matters here is `cds_start_incomplete`.
/// `project_cds_terminus_to_parent` only needs the parent's length (no
/// transcript record), so the `MockProvider` need not carry an `NM_TEST.1`
/// transcript — `terminus_multiaxis_projection`'s `cds_start_incomplete`
/// lookup resolves straight from cdot.
fn ng_parent_projector_with_transcript(
    len: usize,
    cds_start_incomplete: bool,
) -> VariantProjector<MockProvider> {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TEST.1".to_string(),
        CdotTranscript {
            cds_start_incomplete,
            gene_name: Some("TESTGENE".to_string()),
            contig: "NG_TEST.1".to_string(),
            strand: Strand::Plus,
            exons: vec![[0, len as u64, 0, len as u64]],
            cds_start: Some(0),
            cds_end: Some(len as u64),
            gene_id: None,
            protein: Some("NP_TEST.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NG_TEST.1", "A".repeat(len));
    VariantProjector::new(projector, provider)
}

/// RED→GREEN for the same defect class in `terminus_multiaxis_projection`:
/// before the fix, its early return built `coding: Some(...)` unconditionally
/// for a whole-arm `c.pter`/`c.qter` input, with no `cds_start_incomplete`
/// check at all.
#[test]
fn terminus_input_declines_coding_on_cds_start_nf_transcript() {
    let vp = ng_parent_projector_with_transcript(50, true);
    let projection = vp
        .project("NG_TEST.1(NM_TEST.1):c.pterdel", "NM_TEST.1")
        .expect("projection itself must succeed — only the coding axis declines");

    assert!(
        projection.coding.is_none(),
        "cds_start_NF transcript: c. must decline even via the whole-arm \
         terminus multi-axis short-circuit, got {:?}",
        projection.coding.map(|v| v.to_string())
    );
    assert!(
        projection.genomic.is_some(),
        "g. does not depend on the initiation codon and must still render"
    );
}

/// Guardrail: the identical fixture with `cds_start_incomplete: false` still
/// reports the `c.` echo on the terminus short-circuit path — proves the
/// decline above is specific to the flag, not a broken terminus resolver.
#[test]
fn terminus_input_reports_coding_when_cds_start_is_complete() {
    let vp = ng_parent_projector_with_transcript(50, false);
    let projection = vp
        .project("NG_TEST.1(NM_TEST.1):c.pterdel", "NM_TEST.1")
        .expect("projection should succeed");

    let c = projection
        .coding
        .as_ref()
        .expect("c. must render when cds_start is complete")
        .to_string();
    assert!(c.contains(":c.pter"), "got c. = {c}");
}

// =============================================================================
// Task 5 — INPUT-side c./p./r. over a `cds_start_incomplete` transcript,
// mode-gated `W5004` (`IncompleteCdsStartReference`).
//
// Everything above this section is Task 4: OUTPUT-side, i.e. *projecting*
// onto c./p./r. is declined regardless of which axis the input arrived on.
// This section is the complementary INPUT side: the user's OWN `c.`/`p.`/`r.`
// variant is described directly against a `cds_start_incomplete` transcript
// (e.g. `ferro normalize 'ENST00000011700.10:c.50A>G'`) — no projection
// involved at all. `Normalizer::normalize`/`normalize_with_diagnostics` must
// decline to re-number the coordinate (it cannot be trusted without a
// confirmed ATG) and gate acceptance by the active error mode:
//
// - strict:  reject with a `FerroError` carrying W5004.
// - lenient: warn (W5004) and pass the coordinate through UNCHANGED.
// - silent:  pass the coordinate through UNCHANGED, no warning.
//
// `n.` over the identical transcript is untouched in every mode — `n.`
// numbering is transcript-native and does not depend on the CDS start.
// =============================================================================

const CDS_NF_TX: &str = "NM_CDSNF972.1";

/// A minimal coding transcript for the input-side gate: 9-base CDS
/// `"ATGCGCTAA"` (Met-Arg-Stop, cds_start=1, cds_end=9, one exon spanning the
/// whole transcript), `cds_start_incomplete` set per `incomplete`. No genomic
/// sequence is registered — the gate fires (or, for the `incomplete: false`
/// guardrail, a plain substitution needs no reference window) before any
/// genomic bases would be consulted.
fn cds_start_nf_provider(incomplete: bool) -> MockProvider {
    let mut provider = MockProvider::new();
    let mut transcript = Transcript::new(
        CDS_NF_TX.to_string(),
        Some("CDSNFGENE972".to_string()),
        TxStrand::Plus,
        "ATGCGCTAA".to_string(),
        Some(1),
        Some(9),
        vec![Exon::new(1, 1, 9)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1008),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    );
    transcript.cds_start_incomplete = incomplete;
    provider.add_transcript(transcript);
    provider
}

/// RED→GREEN, mode 1/3: strict mode rejects a `c.` input against a
/// `cds_start_NF` transcript with a hard error carrying W5004. Before the
/// fix this normalized/renumbered naively instead.
#[test]
fn strict_mode_rejects_c_input_over_cds_start_nf_transcript() {
    use ferro_hgvs::{FerroError, NormalizeConfig, Normalizer};

    let normalizer =
        Normalizer::with_config(cds_start_nf_provider(true), NormalizeConfig::strict());
    let variant = parse_hgvs(&format!("{CDS_NF_TX}:c.5G>A")).expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("strict mode must reject c. input over a cds_start_NF transcript");
    match err {
        FerroError::InvalidCoordinates { msg } => {
            assert!(
                msg.contains("W5004"),
                "expected the W5004 code in the rejection message, got {msg:?}"
            );
        }
        other => panic!("expected FerroError::InvalidCoordinates (W5004), got {other:?}"),
    }
}

/// RED→GREEN, mode 2/3: lenient mode warns (W5004) and returns the `c.`
/// coordinate UNCHANGED — no re-numbering, since the transcript cannot
/// support it without a confirmed ATG.
#[test]
fn lenient_mode_warns_and_preserves_c_input_over_cds_start_nf_transcript() {
    use ferro_hgvs::{NormalizeConfig, Normalizer};

    let normalizer =
        Normalizer::with_config(cds_start_nf_provider(true), NormalizeConfig::lenient());
    let variant = parse_hgvs(&format!("{CDS_NF_TX}:c.5G>A")).expect("parse");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("lenient mode must accept (only warn), not reject");
    assert_eq!(
        format!("{}", result.result),
        format!("{CDS_NF_TX}:c.5G>A"),
        "lenient mode must pass the coordinate through UNCHANGED — no re-numbering"
    );
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "INCOMPLETE_CDS_START_REFERENCE"),
        "expected a W5004 (IncompleteCdsStartReference) warning in lenient mode, got {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>()
    );
}

/// RED→GREEN, mode 3/3: silent mode passes the `c.` coordinate through
/// UNCHANGED with no warning surfaced at all.
#[test]
fn silent_mode_accepts_c_input_over_cds_start_nf_transcript_without_warning() {
    use ferro_hgvs::{NormalizeConfig, Normalizer};

    let normalizer =
        Normalizer::with_config(cds_start_nf_provider(true), NormalizeConfig::silent());
    let variant = parse_hgvs(&format!("{CDS_NF_TX}:c.5G>A")).expect("parse");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("silent mode must accept");
    assert_eq!(
        format!("{}", result.result),
        format!("{CDS_NF_TX}:c.5G>A"),
        "silent mode must pass the coordinate through UNCHANGED"
    );
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| w.code() == "INCOMPLETE_CDS_START_REFERENCE"),
        "silent mode must not surface a W5004 warning, got {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>()
    );
}

/// `n.` over the identical `cds_start_NF` transcript is unaffected in every
/// mode: `n.` numbering is transcript-native and does not depend on a
/// confirmed CDS start, unlike `c.`/`p.`/`r.`.
#[test]
fn n_input_over_cds_start_nf_transcript_is_unaffected_in_all_modes() {
    use ferro_hgvs::{NormalizeConfig, Normalizer};

    for config in [
        NormalizeConfig::strict(),
        NormalizeConfig::lenient(),
        NormalizeConfig::silent(),
    ] {
        let normalizer = Normalizer::with_config(cds_start_nf_provider(true), config);
        let variant = parse_hgvs(&format!("{CDS_NF_TX}:n.5G>A")).expect("parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("n. must not be gated by cds_start_incomplete in any mode");
        assert_eq!(format!("{}", normalized), format!("{CDS_NF_TX}:n.5G>A"));
    }
}

/// `r.` shares `c.`'s CDS-relative numbering, so the normalize surface must
/// gate it the same way. Strict rejects with W5004; lenient warns and passes
/// the coordinate through unchanged (mirrors the `c.` gate above, exercising
/// the `normalize_rna` gate directly rather than only via projection).
#[test]
fn strict_and_lenient_gate_r_input_over_cds_start_nf_transcript() {
    use ferro_hgvs::{FerroError, NormalizeConfig, Normalizer};

    // Strict: hard error carrying W5004.
    let strict = Normalizer::with_config(cds_start_nf_provider(true), NormalizeConfig::strict());
    let variant = parse_hgvs(&format!("{CDS_NF_TX}:r.5g>a")).expect("parse");
    match strict
        .normalize(&variant)
        .expect_err("strict mode must reject r. input over a cds_start_NF transcript")
    {
        FerroError::InvalidCoordinates { msg } => assert!(
            msg.contains("W5004"),
            "expected W5004 in the rejection message, got {msg:?}"
        ),
        other => panic!("expected FerroError::InvalidCoordinates (W5004), got {other:?}"),
    }

    // Lenient: warn (W5004) and pass through unchanged.
    let lenient = Normalizer::with_config(cds_start_nf_provider(true), NormalizeConfig::lenient());
    let result = lenient
        .normalize_with_diagnostics(&variant)
        .expect("lenient mode must accept (only warn)");
    assert_eq!(
        format!("{}", result.result),
        format!("{CDS_NF_TX}:r.5g>a"),
        "lenient mode must pass the r. coordinate through UNCHANGED"
    );
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "INCOMPLETE_CDS_START_REFERENCE"),
        "expected a W5004 warning in lenient mode, got {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>()
    );
}

/// Guardrail: the identical fixture with `cds_start_incomplete: false`
/// normalizes `c.` normally even in strict mode — proves the rejection
/// above is specific to the flag, not a broken/overly-aggressive normalizer.
#[test]
fn strict_mode_accepts_c_input_when_cds_start_is_complete() {
    use ferro_hgvs::{NormalizeConfig, Normalizer};

    let normalizer =
        Normalizer::with_config(cds_start_nf_provider(false), NormalizeConfig::strict());
    let variant = parse_hgvs(&format!("{CDS_NF_TX}:c.5G>A")).expect("parse");
    let normalized = normalizer
        .normalize(&variant)
        .expect("must not reject when cds_start is complete");
    assert_eq!(format!("{}", normalized), format!("{CDS_NF_TX}:c.5G>A"));
}

// =============================================================================
// Real-reference test — ENST00000011700.10 (VPS13D), genuinely `cds_start_NF`
// (mirrors the `provider()`/`vps13d_coding_variant` helpers above, which back
// the Task 4 projection-side real-transcript tests).
// =============================================================================

/// A bare `c.` input against the real, genuinely-`cds_start_NF` VPS13D
/// transcript must be rejected in strict mode, carrying W5004.
#[test]
fn vps13d_strict_mode_rejects_bare_coding_input() {
    use ferro_hgvs::{FerroError, NormalizeConfig, Normalizer};

    let Some(raw_provider) = provider() else {
        eprintln!("vps13d_strict_mode_rejects_bare_coding_input: skipping — no manifest");
        return;
    };
    let variant = vps13d_coding_variant(&raw_provider);
    let normalizer = Normalizer::with_config(ArcProvider(raw_provider), NormalizeConfig::strict());
    let err = normalizer
        .normalize(&variant)
        .expect_err("strict mode must reject c. input over VPS13D (cds_start_NF)");
    match err {
        FerroError::InvalidCoordinates { msg } => {
            assert!(
                msg.contains("W5004"),
                "expected the W5004 code in the rejection message, got {msg:?}"
            );
        }
        other => panic!("expected FerroError::InvalidCoordinates (W5004), got {other:?}"),
    }
}
