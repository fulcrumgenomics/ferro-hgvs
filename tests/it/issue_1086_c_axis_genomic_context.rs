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

use ferro_hgvs::{parse_hgvs, MultiFastaProvider, VariantProjector};
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

// `Arc<MultiFastaProvider>` is used directly as the projector's provider: the
// library ships a complete blanket `impl ReferenceProvider for Arc<T>`
// (`reference/provider.rs:434`) and `Arc` is `Clone`, so it satisfies
// `VariantProjector<P: ReferenceProvider + Clone>` as-is. A hand-written
// newtype would have to re-delegate every method, and any it missed would
// silently fall back to the trait's stub default — e.g. `get_sequence_length`,
// whose default is `Err(ReferenceNotFound)`, which the poly-A and terminus
// paths below both depend on.
type ArcProvider = Arc<MultiFastaProvider>;

fn manifest_projector() -> Option<VariantProjector<ArcProvider>> {
    let p = provider()?;
    let cdot = p.cdot_mapper()?.clone();
    Some(VariantProjector::new(Projector::new(cdot), p))
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

/// The `r.` axis retains the input's genomic context too — guard it against
/// collateral damage from the `c.`/`n.` framing change. An exonic input renders
/// on `r.` (lowercase RNA alphabet) with the genomic parent echoed, exactly as
/// the module doc notes (`NC_000023.11(NM_004006.2):r.(79a>g)`).
///
/// (An earlier revision pinned the deep-intronic `c.4072-1234_5155-246del`
/// here; after #1089 the `r.` axis declines partial-intronic spans — it accepts
/// only whole-exon-skip shapes — so that input is now `Unavailable`. This test
/// uses an exonic input, which still renders, to keep asserting the framing.)
#[test]
fn rna_axis_framing_is_unchanged() {
    skip_without_manifest!(vp);
    let AxisOutcome::Rendered { output, .. } =
        project(&vp, "NC_000023.11(NM_004006.2):c.79A>G", Axis::Rna)
    else {
        panic!("r. axis unexpectedly unavailable");
    };
    assert_eq!(output, "NC_000023.11(NM_004006.2):r.(79a>g)");
}

/// The non-coding (`n.`) half of the "c./n." guarantee: an input that supplied
/// a genomic context keeps it on the `n.` axis too. `NM_004006.2:c.79` is
/// `n.323` (the 5'UTR is 244 bp), framed under the same `NC_` parent.
#[test]
fn noncoding_axis_retains_context() {
    skip_without_manifest!(vp);
    let AxisOutcome::Rendered { output, .. } =
        project(&vp, "NC_000023.11(NM_004006.2):c.79A>G", Axis::Noncoding)
    else {
        panic!("n. axis unexpectedly unavailable");
    };
    assert_eq!(output, "NC_000023.11(NM_004006.2):n.323A>G");
}

// ---------------------------------------------------------------------------
// The two multi-axis echo paths (`polya_multiaxis_projection`,
// `terminus_multiaxis_projection`), deliberately left out of scope by the
// commit above and folded in here as the same defect.
//
// Both short-circuit the normal genome→CDS re-derivation and instead *echo the
// input's own `c.` form* as the coding axis — and both then explicitly did
// `c.accession.genomic_context = None`. Each justified that stripping in its
// own doc comment as "matching how `coding` is reported on the normal path".
// That justification was true when written and is no longer: the normal path
// now preserves the caller's framing, so the two echoes are the last places
// that still drop it — and an echo is the case where dropping it is least
// defensible, since the parent being discarded is one the caller typed on this
// very description.
//
// Both echo an axis in the *same* reference frame the genomic axis is anchored
// to (the input's `NG_`/`NC_` parent), so neither has a "different reference
// frame" reason to render bare; see the per-test notes.
// ---------------------------------------------------------------------------

/// `polya_multiaxis_projection` (#797): a `c.*` position in the transcript's
/// post-transcriptional poly-A tail has no genomic alignment, so the coding
/// axis is the input echoed rather than re-derived. The echo dropped the `NG_`
/// parent; the genomic axis it is reported alongside is anchored in exactly
/// that parent's frame (`NG_012337.1:g.13954del`), so the bare form made the
/// two axes of one projection disagree on their reference.
#[test]
fn polya_echo_keeps_genomic_context_on_c_axis() {
    skip_without_manifest!(vp);
    assert_coding(
        &vp,
        "NG_012337.1(NM_003002.2):c.*830del",
        "NG_012337.1(NM_003002.2):c.*830del",
    );
}

/// `terminus_multiaxis_projection` (#887): a whole-arm terminus marker
/// (`pter`/`qter`) does not round-trip through cdot, so the coding axis is
/// again the input echoed. A terminus marker is inherently a *genomic*-frame
/// landmark — it is only resolvable at all because the input named a parent
/// (`project_cds_terminus_to_parent` declines outright without one) — so the
/// bare `NM_…:c.pterdel` names a position its own reference cannot locate.
#[test]
fn terminus_echo_keeps_genomic_context_on_c_axis() {
    skip_without_manifest!(vp);
    for marker in ["pter", "qter", "pter_qter"] {
        let input = format!("NG_012337.1(NM_003002.2):c.{marker}del");
        assert_coding(&vp, &input, &input);
    }
}

/// Guard the other half of the rule on both echo paths: an input that supplied
/// no genomic context must not have one invented (#121). Neither bare input
/// may acquire a parenthesised parent.
#[test]
fn echo_paths_never_invent_a_genomic_context() {
    skip_without_manifest!(vp);
    for input in [
        // Poly-A tail position, bare transcript reference.
        "NM_003002.2:c.*830del",
        // Whole-arm terminus on a bare LRG transcript.
        "LRG_9t1:c.qterdel",
        "LRG_9t1:c.pter_qterdel",
    ] {
        let AxisOutcome::Rendered { output, .. } = project(&vp, input, Axis::Coding) else {
            panic!("c. axis for {input}: unexpectedly unavailable");
        };
        assert!(
            !output.contains('('),
            "c. axis for {input} fabricated a genomic context ({output}); a bare \
             input must stay bare (#121)"
        );
    }
}

/// The genomic axis of both echo paths is untouched by the framing change —
/// it stays the parent-frame coordinate each path already reported.
#[test]
fn echo_paths_genomic_axis_is_unchanged() {
    skip_without_manifest!(vp);
    for (input, expected) in [
        // Poly-A: the re-anchored walked coordinate.
        (
            "NG_012337.1(NM_003002.2):c.*830del",
            "NG_012337.1:g.13954del",
        ),
        // Terminus: pter → the parent's 5' end, qter → its 3' end.
        ("NG_012337.1(NM_003002.2):c.pterdel", "NG_012337.1:g.3del"),
        (
            "NG_012337.1(NM_003002.2):c.qterdel",
            "NG_012337.1:g.15948del",
        ),
        (
            "NG_012337.1(NM_003002.2):c.pter_qterdel",
            "NG_012337.1:g.1_15948del",
        ),
    ] {
        let AxisOutcome::Rendered { output, .. } = project(&vp, input, Axis::Genomic) else {
            panic!("g. axis for {input}: unexpectedly unavailable");
        };
        assert_eq!(
            output, expected,
            "g. axis for {input}: expected {expected}, got {output}"
        );
    }
}
