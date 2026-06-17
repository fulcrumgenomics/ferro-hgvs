//! Issue #537 — emit the genomic (g.) axis for `c.pter`/`c.qter` inputs on an
//! NG/LRG/NC-parent reference.
//!
//! PR #534 made ferro spec-correctly *refuse* to number a transcript-flank
//! `c.pter`/`c.qter` marker on the `c.` axis. The genomic axis, however, has a
//! concrete spec-valid answer: the marker denotes the parent reference's own
//! terminus, so `pter` → the 5'-most genomic coordinate (`g.1`) and `qter` →
//! the 3'-most (`g.<length>`); a `pter_qter` range spans the whole reference.
//! These map straight to the parent reference's termini with no transcript
//! mapping — `project_to_genomic` previously rejected them as unresolvable
//! special sentinels.
//!
//! The emitted `g.1`/`g.<length>` is the pre-normalization projection; the
//! normalizer's 3' rule rolls it further (e.g. `g.1del` through a leading
//! homopolymer run, which is how mutalyzer's `NG_012337.1:g.3del` arises). That
//! normalization step is exercised by the conformance harness, not here.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/537>.

use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, FerroError, VariantProjector};

/// A projector whose provider knows only the NG parent's *length* — pter/qter
/// resolution needs nothing else (no cdot transcript, no chromosomal
/// placement). The parent reference is `len` bases long.
fn ng_parent_projector(len: usize) -> VariantProjector<MockProvider> {
    let projector = Projector::new(CdotMapper::new());
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NG_TEST.1", "A".repeat(len));
    VariantProjector::new(projector, provider)
}

/// A projector whose provider knows only the genomic LRG parent's *length*.
/// A bare `LRG_<n>t<m>` transcript has no `genomic_context`, so the terminus
/// path derives its parent structurally (`LRG_<n>`, #480); the parent reference
/// is `len` bases long.
fn lrg_parent_projector(len: usize) -> VariantProjector<MockProvider> {
    let projector = Projector::new(CdotMapper::new());
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("LRG_24", "A".repeat(len));
    VariantProjector::new(projector, provider)
}

#[test]
fn pter_resolves_to_parent_five_prime_terminus() {
    let vp = ng_parent_projector(50);
    let v = parse_hgvs("NG_TEST.1(NM_TEST.1):c.pterdel").expect("parse");
    let g = vp
        .project_to_genomic(&v)
        .expect("pter must project to the parent terminus");
    assert_eq!(g.to_string(), "NG_TEST.1:g.1del");
}

#[test]
fn qter_resolves_to_parent_three_prime_terminus() {
    let vp = ng_parent_projector(50);
    let v = parse_hgvs("NG_TEST.1(NM_TEST.1):c.qterdel").expect("parse");
    let g = vp
        .project_to_genomic(&v)
        .expect("qter must project to the parent terminus");
    // qter → the parent reference's last base (g.<length>).
    assert_eq!(g.to_string(), "NG_TEST.1:g.50del");
}

#[test]
fn pter_qter_range_spans_the_whole_parent_reference() {
    let vp = ng_parent_projector(50);
    let v = parse_hgvs("NG_TEST.1(NM_TEST.1):c.pter_qterdel").expect("parse");
    let g = vp
        .project_to_genomic(&v)
        .expect("pter_qter must span the whole parent reference");
    assert_eq!(g.to_string(), "NG_TEST.1:g.1_50del");
}

#[test]
fn qter_uses_the_actual_parent_length() {
    // A different reference length must change the qter coordinate — confirming
    // it is read from the provider, not a constant.
    let vp = ng_parent_projector(12345);
    let v = parse_hgvs("NG_TEST.1(NM_TEST.1):c.qterdel").expect("parse");
    let g = vp.project_to_genomic(&v).expect("qter projects");
    assert_eq!(g.to_string(), "NG_TEST.1:g.12345del");
}

#[test]
fn lrg_bare_transcript_pter_resolves_to_its_structural_lrg_parent() {
    // A bare `LRG_<n>t<m>` transcript carries no `genomic_context`, so the
    // terminus path derives the genomic parent structurally (`LRG_24`, #480)
    // and resolves `pter` to its 5'-most genomic coordinate (g.1) — exercising
    // the `lrg_genomic_parent` arm, not the `genomic_context` arm.
    let vp = lrg_parent_projector(50);
    let v = parse_hgvs("LRG_24t1:c.pterdel").expect("parse");
    let g = vp
        .project_to_genomic(&v)
        .expect("LRG pter must project to the structural LRG parent terminus");
    assert_eq!(g.to_string(), "LRG_24:g.1del");
}

#[test]
fn lrg_bare_transcript_qter_resolves_to_its_structural_lrg_parent() {
    // The qter mirror of the pter case: the 3'-most genomic coordinate is the
    // structural LRG parent's last base (g.<length>).
    let vp = lrg_parent_projector(50);
    let v = parse_hgvs("LRG_24t1:c.qterdel").expect("parse");
    let g = vp
        .project_to_genomic(&v)
        .expect("LRG qter must project to the structural LRG parent terminus");
    assert_eq!(g.to_string(), "LRG_24:g.50del");
}

#[test]
fn unknown_parent_length_is_an_honest_error() {
    // The provider does not know the parent's length (the pinned parent version
    // is absent — a reference-coverage gap, #645/#672, not a projection bug).
    // `project_to_genomic` surfaces the provider error rather than fabricating a
    // coordinate.
    let projector = Projector::new(CdotMapper::new());
    let vp = VariantProjector::new(projector, MockProvider::new());
    let v = parse_hgvs("NG_TEST.1(NM_TEST.1):c.pterdel").expect("parse");
    let err = vp
        .project_to_genomic(&v)
        .expect_err("absent parent reference must error, not fabricate a coordinate");
    assert!(
        matches!(err, FerroError::ReferenceNotFound { .. }),
        "expected ReferenceNotFound for the absent parent, got: {err:?}"
    );
}

#[test]
fn bare_transcript_pter_falls_through_to_the_no_parent_path() {
    // No `genomic_context` parent and not an LRG transcript: the terminus path
    // declines (returns to the normal projection), which raises the canonical
    // "no parent reference" error rather than the terminus path inventing one.
    let vp = ng_parent_projector(50);
    let v = parse_hgvs("NM_TEST.1:c.pterdel").expect("parse");
    let err = vp
        .project_to_genomic(&v)
        .expect_err("a bare-transcript pter has no parent reference to anchor on");
    assert!(
        matches!(err, FerroError::UnsupportedProjection { .. }),
        "expected UnsupportedProjection (no parent), got: {err:?}"
    );
}
