//! The five-axis projection pass behind the HGVS spec test enumeration's
//! `project-{g,c,n,r,p}` dimensions.
//!
//! Three callers must walk the spec corpus identically, so the pass lives here
//! once rather than three times:
//!
//! - `examples/extract_spec_enumeration_windows.rs` runs it against a real
//!   prepared reference through a recording provider, to capture the committed
//!   hermetic reference slice;
//! - `examples/generate_spec_enumeration.rs` replays it against that committed
//!   slice to produce the enumeration rows;
//! - `tests/it/spec_enumeration_tests.rs` replays it again to verify each row
//!   still reproduces its recorded behaviour.
//!
//! Any divergence between them would mean the fixture is missing bases the
//! replay needs, or that a row's `observed` is not reproducible.
//!
//! The pass is deliberately thin over the public `ferro project` core
//! ([`crate::cli::project`]) — the same [`project_axis`] entry point the
//! `ferro project --axis` subcommand uses. Going through it (rather than
//! reading [`VariantProjection`](crate::project::VariantProjection) fields
//! directly) means the rows also pin the shipped decline-vs-hard-failure
//! classification, which is part of the CLI's contract.

use std::collections::BTreeMap;
use std::path::Path;

use crate::cli::project::{project_axis, select_axis, Axis, AxisOutcome};
use crate::conformance::reference_snapshot::{parse_fasta, render_fasta};
use crate::conformance::reference_window::WindowFixture;
use crate::hgvs::variant::HgvsVariant;
use crate::project::VariantProjector;
use crate::reference::ReferenceProvider;
use crate::FerroError;

/// The five projection axes, in the order the enumeration reports them.
pub const AXES: [(char, Axis); 5] = [
    ('g', Axis::Genomic),
    ('c', Axis::Coding),
    ('n', Axis::Noncoding),
    ('r', Axis::Rna),
    ('p', Axis::Protein),
];

/// Sidecar FASTA carrying the committed slice's **transcript** bases.
pub const SLICE_TRANSCRIPTS_FASTA: &str = "spec_enumeration_transcripts.fna";
/// Sidecar FASTA carrying the committed slice's **genomic window** bases.
pub const SLICE_WINDOWS_FASTA: &str = "spec_enumeration_windows.fna";

/// Move every base string out of `fixture` into two FASTA sidecars, returning
/// `(transcripts, windows)` keyed by transcript id and `"<contig>:<start>"`.
///
/// The reference slice is ~90% raw bases. Left inline, the JSON is a
/// half-megabyte blob that is neither reviewable nor diffable, and trips the
/// repo's large-file guard. Split out, the JSON stays small enough to read as
/// structure (which transcripts, which windows, which placements) and the bases
/// live in the standard format for bases — greppable, and each sidecar
/// independently reviewable. Same division of labour as
/// [`reference_snapshot`](super::reference_snapshot), which pairs
/// `transcripts.metadata.json` with `transcripts.fna`.
pub fn split_slice_sequences(
    fixture: &mut WindowFixture,
) -> (BTreeMap<String, String>, BTreeMap<String, String>) {
    let mut transcripts = BTreeMap::new();
    for tx in &mut fixture.transcripts {
        if let Some(seq) = tx.sequence.take() {
            transcripts.insert(tx.id.clone(), seq);
        }
    }
    let mut windows = BTreeMap::new();
    for w in &mut fixture.genomic {
        windows.insert(window_key(&w.contig, w.start), std::mem::take(&mut w.bases));
    }
    (transcripts, windows)
}

/// Put the sidecar bases back, the inverse of [`split_slice_sequences`].
///
/// A window whose bases are missing is an error, not an empty window: serving
/// zero bases would silently turn a real projection into a spurious decline.
/// A transcript whose bases are missing keeps `sequence: None`, which is a
/// legitimate state for a transcript record.
pub fn attach_slice_sequences(
    fixture: &mut WindowFixture,
    transcripts: &BTreeMap<String, String>,
    windows: &BTreeMap<String, String>,
) -> Result<(), FerroError> {
    for tx in &mut fixture.transcripts {
        if let Some(seq) = transcripts.get(&tx.id) {
            tx.sequence = Some(seq.clone());
        }
    }
    for w in &mut fixture.genomic {
        let key = window_key(&w.contig, w.start);
        let bases = windows
            .get(&key)
            .ok_or(FerroError::ReferenceNotFound { id: key })?;
        w.bases = bases.clone();
    }
    Ok(())
}

/// The FASTA record name for a genomic window.
fn window_key(contig: &str, start: u64) -> String {
    format!("{contig}:{start}")
}

/// Render a split fixture as its three committed files' contents:
/// `(metadata_json, transcripts_fasta, windows_fasta)`.
pub fn render_slice(fixture: &WindowFixture) -> Result<(String, String, String), FerroError> {
    let mut fixture = fixture.clone();
    let (transcripts, windows) = split_slice_sequences(&mut fixture);
    Ok((
        fixture.to_json()?,
        render_fasta(&transcripts),
        render_fasta(&windows),
    ))
}

/// Load a committed slice: the metadata JSON plus its two sidecar FASTAs, which
/// live beside it.
pub fn load_slice<P: AsRef<Path>>(metadata_path: P) -> Result<WindowFixture, FerroError> {
    let metadata_path = metadata_path.as_ref();
    let dir = metadata_path.parent().unwrap_or(Path::new("."));
    let mut fixture = WindowFixture::from_json_path(metadata_path)?;
    let transcripts = parse_fasta(&std::fs::read_to_string(dir.join(SLICE_TRANSCRIPTS_FASTA))?);
    let windows = parse_fasta(&std::fs::read_to_string(dir.join(SLICE_WINDOWS_FASTA))?);
    attach_slice_sequences(&mut fixture, &transcripts, &windows)?;
    Ok(fixture)
}

/// What one axis produced for one input.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AxisResult {
    /// The axis rendered an HGVS string — the informative case.
    Rendered(String),
    /// The engine (or the projection) declined this axis, with a reason. Still
    /// a real behavioural assertion: the `ferro project` contract is that a
    /// decline is an exit-0 outcome carrying a stable reason.
    Unavailable(String),
    /// A hard engine failure (reference not found, ambiguous input, ...).
    Error(String),
    /// Structurally inapplicable: the input names no transcript to project on
    /// (a bare `p.`/`m.`/`o.` input, an allele list, ...). Not emitted as a row.
    NotApplicable(String),
}

impl AxisResult {
    /// A single-line rendering used as the row's `observed` value.
    pub fn as_observed(&self) -> String {
        match self {
            AxisResult::Rendered(s) => s.clone(),
            AxisResult::Unavailable(r) => format!("unavailable: {r}"),
            AxisResult::Error(e) => format!("error: {e}"),
            AxisResult::NotApplicable(r) => format!("not applicable: {r}"),
        }
    }
}

/// The transcript a variant projects on, plus the per-axis results.
#[derive(Debug, Clone)]
pub struct PassResult {
    /// The transcript the projection ran on, when one could be chosen.
    pub transcript: Option<String>,
    /// One entry per axis letter in [`AXES`] order.
    pub axes: BTreeMap<char, AxisResult>,
}

impl PassResult {
    /// Every axis reports the same structural non-applicability.
    fn all_not_applicable(reason: &str) -> Self {
        PassResult {
            transcript: None,
            axes: AXES
                .iter()
                .map(|(c, _)| (*c, AxisResult::NotApplicable(reason.to_string())))
                .collect(),
        }
    }
}

/// Choose the transcript to project `variant` onto.
///
/// Only a transcript-coordinate input (`c.`/`n.`/`r.`) has an answer: its own
/// accession *is* the transcript. Everything else returns `None`.
///
/// A **bare genomic** input is deliberately included in that `None`, and the
/// reason is not performance. `ferro project` refuses a bare `g.` input without
/// an explicit `--transcript`, because which transcript to read a chromosomal
/// position against is a caller's clinical choice — and the HGVS spec states no
/// such choice for its examples. Silently picking the engine's first overlapping
/// transcript would manufacture an expectation the spec never wrote, which is
/// exactly what these rows must not do. (The `g.` axis is still exercised, in
/// the direction that carries information: `c.`/`n.`/`r.` → `g.`.)
///
/// Protein, mitochondrial, circular and allele-list inputs carry no transcript
/// frame at all; `ferro project` cannot project them either.
pub fn transcript_for<P: ReferenceProvider + Clone>(
    _projector: &VariantProjector<P>,
    variant: &HgvsVariant,
) -> Option<String> {
    match variant {
        HgvsVariant::Cds(_) | HgvsVariant::Tx(_) | HgvsVariant::Rna(_) => {
            variant.accession().map(|a| a.transcript_accession())
        }
        _ => None,
    }
}

/// Run the five-axis pass for one parsed variant.
pub fn project_all_axes<P: ReferenceProvider + Clone>(
    projector: &VariantProjector<P>,
    variant: &HgvsVariant,
) -> PassResult {
    let Some(transcript) = transcript_for(projector, variant) else {
        return PassResult::all_not_applicable(
            "input carries no transcript frame to project from (bare g./p./m./o. \
             or an allele list); `ferro project` declines these too",
        );
    };

    // One projection serves all five axes — they describe one variant, so
    // re-projecting per axis would be five times the work for the same answer.
    let axes = match projector.project_variant(variant, &transcript) {
        Ok(projection) => AXES
            .iter()
            .map(|(code, axis)| {
                let result = match select_axis(&projection, *axis) {
                    AxisOutcome::Rendered { output, .. } => AxisResult::Rendered(output),
                    AxisOutcome::Unavailable { reason, .. } => AxisResult::Unavailable(reason),
                };
                (*code, result)
            })
            .collect(),
        // A whole-projection failure applies identically to all five axes, but
        // decline (exit 0) and hard failure (exit 1) are different contracts.
        // Rather than restate the shipped classification here — where it could
        // silently drift — ask `project_axis` once and reuse its verdict. The
        // probe is genomic, so its non-rendered *message* is genomic-flavoured;
        // re-frame it axis-neutrally before fanning it out so a c./n./r./p. row
        // does not read as though the failure were specific to the g. axis.
        Err(_) => {
            let verdict = match project_axis(projector, variant, Axis::Genomic, Some(&transcript)) {
                Ok(AxisOutcome::Rendered { output, .. }) => AxisResult::Rendered(output),
                Ok(AxisOutcome::Unavailable { reason, .. }) => AxisResult::Unavailable(format!(
                    "projection unavailable for all axes: {reason}"
                )),
                Err(e) => AxisResult::Error(format!("projection failed for all axes: {e}")),
            };
            AXES.iter()
                .map(|(code, _)| (*code, verdict.clone()))
                .collect()
        }
    };
    PassResult {
        transcript: Some(transcript),
        axes,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::cdot::{CdotMapper, CdotTranscript};
    use crate::data::projection::Projector;
    use crate::reference::mock::MockProvider;
    use crate::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
    use crate::reference::Strand;

    /// One plus-strand coding transcript `NM_TEST.1` on `NC_000001.11`
    /// [1000..1008], CDS = the whole 9-base exon `ATGCGCTAA`.
    fn fixture() -> VariantProjector<MockProvider> {
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_TEST.1".to_string(),
            CdotTranscript {
                cds_start_incomplete: false,
                gene_name: Some("TESTGENE".to_string()),
                contig: "NC_000001.11".to_string(),
                strand: Strand::Plus,
                exons: vec![[1000, 1009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: Some("NP_TEST.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript::new(
            "NM_TEST.1".to_string(),
            Some("TESTGENE".to_string()),
            TxStrand::Plus,
            "ATGCGCTAA".to_string(),
            Some(1),
            Some(9),
            vec![Exon::new(1, 1, 9)],
            Some("NC_000001.11".to_string()),
            Some(1000),
            Some(1008),
            Default::default(),
            ManeStatus::default(),
            None,
            None,
        ));
        provider.add_genomic_sequence(
            "NC_000001.11",
            format!("{}ATGCGCTAA{}", "N".repeat(1000), "N".repeat(100)),
        );
        VariantProjector::new(Projector::new(cdot), provider)
    }

    #[test]
    fn axes_cover_the_five_projection_targets_exactly_once() {
        let mut codes: Vec<char> = AXES.iter().map(|(c, _)| *c).collect();
        codes.sort_unstable();
        assert_eq!(codes, vec!['c', 'g', 'n', 'p', 'r']);
    }

    #[test]
    fn coding_input_projects_every_axis_from_one_pass() {
        let vp = fixture();
        let v = crate::parse_hgvs("NM_TEST.1:c.4C>A").expect("parses");
        let pass = project_all_axes(&vp, &v);

        assert_eq!(pass.transcript.as_deref(), Some("NM_TEST.1"));
        // Every axis is reported — the pass never silently omits one.
        assert_eq!(pass.axes.len(), AXES.len());
        // The informative axes render real HGVS strings on the right accessions.
        assert_eq!(
            pass.axes[&'c'],
            AxisResult::Rendered("NM_TEST.1:c.4C>A".to_string())
        );
        // This mock transcript carries no genome alignment, so the `g.` axis
        // declines — and the decline is reported as a first-class outcome with
        // its reason, not dropped.
        assert_eq!(
            pass.axes[&'g'],
            AxisResult::Unavailable("no g. representation for this variant".to_string())
        );
        assert_eq!(
            pass.axes[&'p'],
            AxisResult::Rendered("NP_TEST.1:p.(Arg2Ser)".to_string())
        );
    }

    #[test]
    fn bare_genomic_input_is_not_applicable_not_an_arbitrary_transcript() {
        // A bare `g.` input overlaps NM_TEST.1, but the spec states no
        // transcript choice for it, so the pass must decline rather than pick
        // one and manufacture an expectation.
        let vp = fixture();
        let v = crate::parse_hgvs("NC_000001.11:g.1004C>A").expect("parses");
        let pass = project_all_axes(&vp, &v);
        assert_eq!(pass.transcript, None);
        for (code, _) in AXES {
            assert!(
                matches!(pass.axes[&code], AxisResult::NotApplicable(_)),
                "axis {code} should be not-applicable for a bare genomic input"
            );
        }
    }

    #[test]
    fn protein_input_has_no_transcript_frame() {
        let vp = fixture();
        let v = crate::parse_hgvs("NP_TEST.1:p.Arg2Ser").expect("parses");
        assert_eq!(transcript_for(&vp, &v), None);
    }

    #[test]
    fn unknown_transcript_is_a_pinned_error_not_a_silent_drop() {
        // An accession absent from the hermetic slice must surface as a real
        // outcome on every axis, never as a missing row.
        let vp = fixture();
        let v = crate::parse_hgvs("NM_ABSENT.1:c.4C>A").expect("parses");
        let pass = project_all_axes(&vp, &v);
        assert_eq!(pass.transcript.as_deref(), Some("NM_ABSENT.1"));
        for (code, _) in AXES {
            assert!(
                matches!(pass.axes[&code], AxisResult::Error(_)),
                "axis {code} must surface the unknown transcript as a hard error, \
                 not drop it or soften it to unavailable: got {:?}",
                pass.axes[&code]
            );
        }
    }

    #[test]
    fn observed_rendering_distinguishes_the_outcome_kinds() {
        assert_eq!(
            AxisResult::Rendered("NM_TEST.1:c.4C>A".to_string()).as_observed(),
            "NM_TEST.1:c.4C>A"
        );
        assert_eq!(
            AxisResult::Unavailable("no p. representation".to_string()).as_observed(),
            "unavailable: no p. representation"
        );
        assert_eq!(
            AxisResult::Error("boom".to_string()).as_observed(),
            "error: boom"
        );
    }
}
