//! #843 — the predicted-`r.` allele transcript fetch must be build-scoped by the
//! input's `NC_*` build, not resolved build-agnostically against the primary
//! (GRCh38) transcript. Real-reference, end-to-end routing guard, gated on
//! FERRO_MANIFEST (skips in manifest-less CI like the other manifest tests).
//!
//! The deterministic red→green for the bug lives in the hermetic unit tests
//! `multi_build_projection_tests::allele_predicted_rna_is_build_scoped` /
//! `…_does_not_pollute_cache_with_parentless_key` (they model a transcript whose
//! GRCh37/GRCh38 bases differ at a deletion 3'-shift boundary — a divergence only
//! genome-reconstructed transcripts exhibit, and only at a shift boundary, which
//! is vanishingly rare on real data). This integration test instead pins the
//! issue's literal acceptance — a GRCh37 `g.` allele → `r.` — end-to-end on the
//! real dual-build reference: it confirms the allele path ROUTES GRCh37
//! (resolves the GRCh37 transcript and emits a GRCh37-framed `r.`), and that the
//! `g.`-allele prediction agrees with the equivalent `NC_(NM_)`-parented c.
//! allele on the same build (the parented form was already build-scoped).
//!
//! Transcript: `XM_005254401.2` (on `NC_000015.9` GRCh37 / `NC_000015.10`
//! GRCh38) — a genome-reconstructed transcript whose bases differ at a deletion
//! 3'-shift boundary, so a real `g.`-allele `c.1191del` 3'-shifts to a
//! **build-distinct** predicted `r.`: `r.1193del` on GRCh37, `r.1191del` on
//! GRCh38. Pre-#843 the bare `g.` allele read the primary (GRCh38) transcript
//! for `predict_rna`, so the GRCh37 allele did NOT emit its own `r.1193del`
//! (it degraded to no/primary-build prediction) — a genuine, reproducible
//! red→green on the real reference, not just a hermetic mock.

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

/// `Arc<MultiFastaProvider>` newtype that forwards every build-aware method to
/// the inner provider, so build routing (`get_transcript_for_accession`,
/// `infer_genome_build`) reaches the real `MultiFastaProvider` logic.
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

const TX: &str = "XM_005254401.2";

/// Project a single-deletion `g.` allele `NC_*(.ver):g.[<pos>del]` and return the
/// predicted `r.` string (or `None` if no prediction).
fn allele_rna(vp: &VariantProjector<ArcProvider>, nc: &str, pos: u64) -> Option<String> {
    let input = format!("{nc}:g.[{pos}del]");
    let v = parse_hgvs(&input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    vp.project_variant(&v, TX).ok()?.rna.map(|r| r.to_string())
}

/// The equivalent NC-parented c.-allele prediction on the same build. The
/// parented form ALREADY carries the build via `genomic_context`, so it is the
/// build-correct oracle the bare `g.` allele must match once #843 build-scopes
/// the predict_rna fetch.
fn parented_c_allele_rna(
    vp: &VariantProjector<ArcProvider>,
    nc: &str,
    cpos: u64,
) -> Option<String> {
    let input = format!("{nc}({TX}):c.[{cpos}del]");
    let v = parse_hgvs(&input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    vp.project_variant(&v, TX).ok()?.rna.map(|r| r.to_string())
}

/// Normalize a predicted-`r.` string to just its edit body (`(1193del)` /
/// `1193del`), dropping the accession prefix and any predicted `( )` wrapper, so
/// a bare `g.` allele's `r.<edit>` and a parented c. allele's `r.(<edit>)`
/// compare on coordinate/edit alone (the build-sensitive part).
fn rna_edit_body(rna: &str) -> String {
    let after = rna.split(":r.").nth(1).unwrap_or(rna);
    after
        .trim_start_matches('(')
        .trim_end_matches(')')
        .to_string()
}

/// The build-divergence oracle for `XM_005254401.2:c.1191del`: GRCh37 3'-shifts
/// to `1193del`, GRCh38 stays `1191del`. Carries each build's genomic position
/// for `c.1191` and its parented-c predicted-`r.` edit body; the helper returns
/// `None` (skip) only if the manifest's placement drifted so the parented-c
/// oracle no longer predicts.
struct Oracle {
    g37_pos: u64,
    g38_pos: u64,
    c37_body: String,
    c38_body: String,
}

fn oracle(vp: &VariantProjector<ArcProvider>) -> Option<Oracle> {
    // Genomic coordinates of c.1191 on each build (from the parented c. allele's
    // `.genomic`), so the bare `g.` allele uses each build's own coordinate.
    let g_of = |nc: &str| -> Option<u64> {
        let v = parse_hgvs(&format!("{nc}({TX}):c.1191del")).ok()?;
        let g = vp.project_variant(&v, TX).ok()?.genomic?;
        format!("{g}")
            .split(":g.")
            .nth(1)?
            .split(|c: char| !c.is_ascii_digit())
            .next()?
            .parse::<u64>()
            .ok()
    };
    let c37_body = rna_edit_body(&parented_c_allele_rna(vp, "NC_000015.9", 1191)?);
    let c38_body = rna_edit_body(&parented_c_allele_rna(vp, "NC_000015.10", 1191)?);
    Some(Oracle {
        g37_pos: g_of("NC_000015.9")?,
        g38_pos: g_of("NC_000015.10")?,
        c37_body,
        c38_body,
    })
}

/// A bare GRCh37 `g.` allele must predict the GRCh37-shifted `r.` — its own
/// `r.1193del`, NOT the GRCh38 primary's `r.1191del`. Pre-#843 the allele's
/// `predict_rna` read the primary (GRCh38) transcript, so the GRCh37 allele did
/// not produce its build-correct `r.` (it degraded to none/primary). Post-fix it
/// matches the build-scoped parented-c oracle. (The GRCh38 companion is checked
/// in the same test to lock the builds apart.)
#[test]
fn g_allele_predict_rna_is_build_scoped_on_real_reference() {
    let Some(vp) = manifest_projector() else {
        eprintln!("issue_843: skipping — no manifest");
        return;
    };
    let Some(o) = oracle(&vp) else {
        eprintln!("issue_843: skipping — parented-c oracle unavailable (data drift)");
        return;
    };

    // Pin the divergence premise: the two builds' parented oracles must differ,
    // else this transcript no longer exercises the bug (guards against silent
    // data drift that would make the test a tautology).
    assert_ne!(
        o.c37_body, o.c38_body,
        "test premise broken: XM_005254401.2 c.1191del no longer 3'-shifts \
         differently across builds (37={}, 38={}); pick a new divergent locus",
        o.c37_body, o.c38_body
    );

    // Bare g. alleles, each on its own build's genomic coordinate.
    let g37 = allele_rna(&vp, "NC_000015.9", o.g37_pos)
        .map(|r| rna_edit_body(&r))
        .expect(
            "GRCh37 bare g. allele MUST predict an r.; a missing prediction means \
             predict_rna read the wrong (primary-build) transcript (#843 pre-fix)",
        );
    let g38 = allele_rna(&vp, "NC_000015.10", o.g38_pos)
        .map(|r| rna_edit_body(&r))
        .expect("GRCh38 bare g. allele MUST predict an r.");

    // Each bare g. allele must equal its OWN build's oracle, not the primary's.
    assert_eq!(
        g37, o.c37_body,
        "GRCh37 g. allele r. ({g37}) must match the GRCh37 oracle ({}); \
         a mismatch (esp. equal to the GRCh38 value {}) means predict_rna read \
         the primary (GRCh38) transcript",
        o.c37_body, o.c38_body
    );
    assert_eq!(
        g38, o.c38_body,
        "GRCh38 g. allele r. ({g38}) must match the GRCh38 oracle ({})",
        o.c38_body
    );
    // And the two builds must remain distinct end-to-end.
    assert_ne!(
        g37, g38,
        "GRCh37 and GRCh38 g. alleles must predict build-distinct r."
    );
}
