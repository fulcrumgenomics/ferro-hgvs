//! Issue #952 — c.→p. consequence prediction must *decline observably* (never
//! silently) for the edit kinds that have no protein predictor yet.
//!
//! `predict_protein_consequence` returns `Ok(None)` for Repeat / MultiRepeat /
//! DupIns / Conversion / NPaddedDeletion / SubstitutionNoRef /
//! BreakpointInsertion / Identity. Before #952 that decline was a bare
//! `_ => {}` with no comment and no log, indistinguishable from a deliberate
//! "no consequence". #952 keeps the `None` (implementing these consequences is
//! future work) but makes the gap explicit + logged. These tests lock the
//! contract: such a coding input projects *without error or panic* and simply
//! carries no protein axis — the projection itself still succeeds.

use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::VariantProjector;

/// Mock reference with a clean `ATG`-start CDS: `ATG CGC GCG TAA`
/// (Met-Arg-Ala-Stop), coding bases 1..12 on a single exon.
fn make_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let tx = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("MYGENE".to_string()),
        TxStrand::Plus,
        "ATGCGCGCGTAA".to_string(),
        Some(1),
        Some(12),
        vec![Exon::with_genomic(1, 1, 12, 1000, 1011)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1011),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    )
    .with_protein_id(Some("NP_TEST.1".to_string()));
    provider.add_transcript(tx);
    let prefix = "N".repeat(999);
    let suffix = "N".repeat(100);
    provider.add_genomic_sequence("chr1", format!("{}ATGCGCGCGTAA{}", prefix, suffix));
    provider
}

fn projector() -> VariantProjector<MockProvider> {
    let provider = make_provider();
    let cdot = CdotMapper::from_transcripts(provider.all_transcripts());
    let proj = Projector::new(cdot);
    VariantProjector::new(proj, provider)
}

/// Each of these coding inputs carries an edit kind with no c.→p. predictor
/// (Repeat / MultiRepeat / NPaddedDeletion / SubstitutionNoRef / Identity), so
/// `predict_protein_consequence` reaches the #952 decline arm. The projection
/// must still succeed — no error, no panic — and simply carry no protein axis.
///
/// (`Conversion` is deliberately absent: it is rewritten to a `Delins` with the
/// fetched source sequence *before* this dispatch, so it does get a protein
/// consequence and never reaches the decline arm. `DupIns` is absent because the
/// `dupins` notation is parse-rejected — DNA/duplication.md:92 — so no string can
/// reach that arm.)
#[test]
fn unpredictable_edit_kinds_decline_gracefully_without_protein() {
    let vp = projector();
    for input in [
        "NM_TEST.1:c.4_6[3]",        // Repeat
        "NM_TEST.1:c.4CGC[2]GCG[2]", // MultiRepeat
        "NM_TEST.1:c.4_6delN[3]",    // NPaddedDeletion
        "NM_TEST.1:c.4>A",           // SubstitutionNoRef
        "NM_TEST.1:c.4=",            // Identity (position-specific)
    ] {
        let variant =
            ferro_hgvs::parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
        let projection = vp
            .project_variant(&variant, "NM_TEST.1")
            .unwrap_or_else(|e| panic!("project {input} must not error: {e}"));
        assert!(
            projection.protein.is_none(),
            "{input}: no c.→p. predictor exists for this edit kind, so the projection must \
             decline (protein=None), not fabricate a consequence — got {:?}",
            projection.protein.as_ref().map(ToString::to_string),
        );
    }
}

/// Guardrail contrast: a plain substitution on the same transcript *does* get a
/// protein consequence. This proves the decline above is specific to the
/// unpredictable edit kinds and not a broken projector that drops every protein.
#[test]
fn substitution_still_predicts_protein() {
    let vp = projector();
    let variant = ferro_hgvs::parse_hgvs("NM_TEST.1:c.4C>A").expect("parse");
    let projection = vp
        .project_variant(&variant, "NM_TEST.1")
        .expect("substitution projects");
    assert_eq!(
        projection
            .protein
            .as_ref()
            .map(ToString::to_string)
            .as_deref(),
        Some("NP_TEST.1:p.(Arg2Ser)"),
        "a plain substitution must still predict its protein consequence",
    );
}
