//! CLI-agnostic core for the `ferro project` subcommand: select one output
//! axis from a [`VariantProjection`] and classify engine errors into
//! "unavailable" (exit 0) vs "hard failure" (exit 1).

use crate::error::FerroError;
use crate::hgvs::variant::{CoordinateAxis, HgvsVariant};
use crate::normalize::NormalizationWarning;
use crate::project::{VariantProjection, VariantProjector};
use crate::reference::ReferenceProvider;

/// The output axis selected by `--axis`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    Genomic,
    Coding,
    Noncoding,
    Protein,
    Rna,
}

impl Axis {
    /// Parse the clap-validated `--axis` value (`g`/`c`/`n`/`p`/`r`).
    pub fn parse(s: &str) -> Option<Axis> {
        match s {
            "g" => Some(Axis::Genomic),
            "c" => Some(Axis::Coding),
            "n" => Some(Axis::Noncoding),
            "p" => Some(Axis::Protein),
            "r" => Some(Axis::Rna),
            _ => None,
        }
    }

    /// The one-letter axis code.
    pub fn code(self) -> &'static str {
        match self {
            Axis::Genomic => "g",
            Axis::Coding => "c",
            Axis::Noncoding => "n",
            Axis::Protein => "p",
            Axis::Rna => "r",
        }
    }
}

/// A normalization warning raised while projecting, flattened for reporting.
///
/// [`NormalizationWarning`] is `#[non_exhaustive]` and implements neither
/// `PartialEq` nor `Eq`, so it cannot be embedded in [`AxisOutcome`] (which
/// derives both, and whose equality several tests rely on). Flattening to the
/// code plus the rendered message keeps those derives and is all the CLI needs
/// — the code for machine consumers, the message for humans. #1182.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProjectionWarning {
    /// The warning's stable code (e.g. `POSITION_PAST_END`).
    pub code: String,
    /// The rendered, human-readable message.
    pub message: String,
}

impl ProjectionWarning {
    /// Flatten a normalizer warning for reporting.
    pub fn from_normalization(warning: &NormalizationWarning) -> Self {
        Self {
            code: warning.code().to_string(),
            message: warning.to_string(),
        }
    }

    /// Flatten a whole set, in the order the normalizer raised them.
    pub fn from_normalization_set(warnings: &[NormalizationWarning]) -> Vec<Self> {
        warnings.iter().map(Self::from_normalization).collect()
    }
}

/// The result of selecting one axis from a projection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AxisOutcome {
    /// The axis rendered to an HGVS string.
    Rendered {
        transcript_id: String,
        output: String,
        /// Warnings raised while normalizing the input for this projection.
        /// Empty when it normalized cleanly. Previously discarded outright, so
        /// diagnostics such as W4004 (position past the end of the reference)
        /// could not reach the user from `ferro project` in any configuration.
        warnings: Vec<ProjectionWarning>,
    },
    /// The axis is legitimately not available (exit 0); carries a reason and,
    /// when known, the transcript the projection was attempted on.
    Unavailable {
        transcript_id: Option<String>,
        reason: String,
        /// Warnings raised before the axis turned out to be unavailable. An
        /// unavailable axis can still have produced diagnostics worth showing.
        warnings: Vec<ProjectionWarning>,
    },
}

/// A hard failure for `ferro project` — maps to a nonzero exit (parse error,
/// reference/transcript not found, `--transcript` mismatch, ambiguous bare-g,
/// IO). Distinct from [`AxisOutcome::Unavailable`], which is exit 0.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProjectCliError(pub String);

impl std::fmt::Display for ProjectCliError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

/// Whether an engine `FerroError` means "this representation is not available"
/// (a clean decline → exit 0) rather than a hard failure. The #662 NG_/LRG_
/// re-anchor decline and the p./m./o./allele / missing-`genomic_context` / `?`
/// cases all surface as `UnsupportedProjection`.
fn engine_error_is_unavailable(e: &FerroError) -> bool {
    matches!(
        e,
        FerroError::UnsupportedProjection { .. }
            | FerroError::GenomicReferenceNotAvailable { .. }
            | FerroError::ProteinReferenceNotAvailable { .. }
            | FerroError::ProteinSequenceUnavailable { .. }
            | FerroError::IntronicVariant { .. }
            | FerroError::TranscriptVersionNotExact { .. }
            | FerroError::TranscriptSequenceUnreconstructable { .. }
    )
}

/// Read the chosen axis off a completed projection, rendering its HGVS string
/// or reporting it unavailable.
pub fn select_axis(projection: &VariantProjection, axis: Axis) -> AxisOutcome {
    let declines = &projection.axis_decline_reasons;
    // Read the axis and its recorded decline reason together, so the reason can
    // never be taken from a *different* axis than the one being reported.
    let (field, decline) = match axis {
        Axis::Genomic => (&projection.genomic, &declines.genomic),
        Axis::Coding => (&projection.coding, &declines.coding),
        Axis::Noncoding => (&projection.noncoding, &declines.noncoding),
        Axis::Protein => (&projection.protein, &declines.protein),
        Axis::Rna => (&projection.rna, &declines.rna),
    };
    let tx = projection.transcript_id.clone();
    // #1182: the projection's normalization warnings were read off nowhere and
    // dropped here, which is what made W4004 (and every other normalizer
    // diagnostic) structurally unreachable from `ferro project`. Carry them onto
    // both outcomes.
    let warnings = ProjectionWarning::from_normalization_set(&projection.normalization_warnings);
    match field {
        Some(v) => AxisOutcome::Rendered {
            transcript_id: tx,
            output: v.to_string(),
            warnings,
        },
        None => AxisOutcome::Unavailable {
            transcript_id: Some(tx),
            reason: unavailable_reason(axis, decline.as_deref()),
            warnings,
        },
    }
}

/// Phrase why `axis` is unavailable.
///
/// The engine's own explanation *replaces* the synthesized string when one was
/// recorded (#1198). The synthesized string is derived from nothing but the axis
/// code, so it restates what both output formats already say — text writes
/// `[<axis>]: unavailable:` and JSON carries `"axis"` plus
/// `"status": "unavailable"` — while the real explanation was computed, logged
/// at `trace`, and dropped.
///
/// With no recorded explanation the wording is unchanged. An axis can be absent
/// because it was never applicable, or because its decline is not one the
/// projector records (see [`crate::project::AxisDeclineReasons`]) — and
/// inventing a cause for
/// either would trade one generic answer for a wrong one.
///
/// This does **not** widen what exits 0. Which errors are a clean decline is
/// decided by [`engine_error_is_unavailable`] for a failed *projection*, and by
/// the projector's best-effort axis contract for a failed *axis* — an
/// underivable axis has always been reported as unavailable at exit 0, whatever
/// its error's class, so that the other four axes survive. This function only
/// changes what that record says, never whether it is emitted.
fn unavailable_reason(axis: Axis, decline: Option<&str>) -> String {
    match decline {
        Some(explanation) => explanation.to_string(),
        None => format!("no {}. representation for this variant", axis.code()),
    }
}

/// The accession base (everything before the version dot), for version-flex
/// matching of a user-supplied `--transcript` against enumerated results.
fn accession_base(acc: &str) -> &str {
    acc.split('.').next().unwrap_or(acc)
}

/// Project `variant` to the requested `axis`.
///
/// - Transcript-coordinate input (`c.`/`n.`/`r.`): the transcript is taken from
///   the input accession; a supplied `--transcript` is cross-checked and a
///   mismatch is a hard error.
/// - Bare genomic input (`g.`): `--transcript` is required; the engine
///   enumerates overlapping transcripts and the chosen one is selected.
pub fn project_axis<P: ReferenceProvider + Clone>(
    projector: &VariantProjector<P>,
    variant: &HgvsVariant,
    axis: Axis,
    transcript: Option<&str>,
) -> Result<AxisOutcome, ProjectCliError> {
    // Route every *genomic-axis* input through the genomic dispatch, not just a
    // lone `HgvsVariant::Genome`. A genomic cis allele
    // (`NC_…:g.[261del;263_264insA]`) parses to `HgvsVariant::Allele`, whose
    // accession is a *chromosome*, not a transcript — so the transcript-
    // coordinate branch below would compare `--transcript` against that
    // chromosome accession and reject a valid projection (#1562). Keying on the
    // coordinate axis rather than the variant kind catches the allele (and any
    // other genomic-axis shape) the same way a lone `g.` variant is caught.
    if variant.coordinate_axis() == Some(CoordinateAxis::Genomic) {
        return project_axis_genomic(projector, variant, axis, transcript);
    }

    // Transcript-coordinate input: the transcript is named in the accession.
    let input_tx = variant
        .accession()
        .map(|a| a.transcript_accession())
        .ok_or_else(|| ProjectCliError(format!("input {variant} has no transcript to project")))?;

    if let Some(requested) = transcript {
        if accession_base(requested) != accession_base(&input_tx) {
            return Err(ProjectCliError(format!(
                "--transcript {requested} does not match the input's transcript {input_tx}"
            )));
        }
    }

    match projector.project_variant(variant, &input_tx) {
        Ok(projection) => Ok(select_axis(&projection, axis)),
        Err(e) if engine_error_is_unavailable(&e) => Ok(AxisOutcome::Unavailable {
            transcript_id: Some(input_tx),
            reason: e.to_string(),
            // The projection never completed, so there is no warning set to carry.
            warnings: Vec::new(),
        }),
        Err(e) => Err(ProjectCliError(e.to_string())),
    }
}

/// Bare-genomic dispatch: enumerate overlapping transcripts and pick one.
fn project_axis_genomic<P: ReferenceProvider + Clone>(
    projector: &VariantProjector<P>,
    variant: &HgvsVariant,
    axis: Axis,
    transcript: Option<&str>,
) -> Result<AxisOutcome, ProjectCliError> {
    let Some(requested) = transcript else {
        // Ambiguous: list the overlapping transcripts so the user can choose.
        let ids = enumerate_transcript_ids(projector, variant)?;
        let listed = if ids.is_empty() {
            "none overlap this position".to_string()
        } else {
            ids.join(", ")
        };
        return Err(ProjectCliError(format!(
            "a genomic (g.) input requires --transcript; overlapping transcripts \
             (longest-CDS, then alphabetical): {listed}"
        )));
    };

    // Project directly onto the requested transcript rather than enumerating and
    // filtering: `project_variant_all` silently drops per-transcript errors, so a
    // transcript that overlaps but whose projection *declines* would otherwise be
    // misreported as "does not overlap". Going direct surfaces a genuine decline
    // as Unavailable (exit 0) and a non-overlap / hard error as a hard failure.
    match projector.project_variant(variant, requested) {
        Ok(projection) => Ok(select_axis(&projection, axis)),
        Err(e) if engine_error_is_unavailable(&e) => Ok(AxisOutcome::Unavailable {
            transcript_id: Some(requested.to_string()),
            reason: e.to_string(),
            // The projection never completed, so there is no warning set to carry.
            warnings: Vec::new(),
        }),
        Err(e) => Err(ProjectCliError(e.to_string())),
    }
}

/// Enumerate the transcript ids overlapping a genomic input (for the
/// "needs --transcript" error message), in the engine's clinical-priority order
/// (MANE/canonical sets are empty in CLI context, so this is longest-CDS then
/// alphabetical).
fn enumerate_transcript_ids<P: ReferenceProvider + Clone>(
    projector: &VariantProjector<P>,
    variant: &HgvsVariant,
) -> Result<Vec<String>, ProjectCliError> {
    match projector.project_variant_all(variant) {
        Ok(ps) => Ok(ps.into_iter().map(|p| p.transcript_id).collect()),
        Err(e) if engine_error_is_unavailable(&e) => Ok(Vec::new()),
        Err(e) => Err(ProjectCliError(e.to_string())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::cdot::{CdotMapper, CdotTranscript};
    use crate::data::projection::Projector;
    use crate::project::{AxisDeclineReasons, VariantProjection};
    use crate::reference::mock::MockProvider;
    use crate::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
    use crate::reference::Strand;

    #[test]
    fn axis_parse_roundtrips_codes() {
        for code in ["g", "c", "n", "p", "r"] {
            assert_eq!(Axis::parse(code).unwrap().code(), code);
        }
        assert_eq!(Axis::parse("m"), None);
        assert_eq!(Axis::parse(""), None);
    }

    #[test]
    fn unsupported_projection_is_unavailable_not_hard() {
        let e = FerroError::UnsupportedProjection {
            reason: "no chromosomal placement is known for NG_TEST.1".to_string(),
        };
        assert!(engine_error_is_unavailable(&e));
    }

    #[test]
    fn transcript_sequence_unreconstructable_is_unavailable_not_hard() {
        let e = FerroError::TranscriptSequenceUnreconstructable {
            id: "NM_TEST.1".to_string(),
            insertions: 1,
        };
        assert!(engine_error_is_unavailable(&e));
    }

    #[test]
    fn parse_and_not_found_are_hard() {
        assert!(!engine_error_is_unavailable(
            &FerroError::ReferenceNotFound {
                id: "NM_X".to_string()
            }
        ));
        assert!(!engine_error_is_unavailable(&FerroError::Parse {
            msg: "bad".to_string(),
            pos: 0,
            diagnostic: None,
        }));
    }

    fn projection_with(genomic: Option<&str>, protein: Option<&str>) -> VariantProjection {
        VariantProjection {
            normalization_warnings: Vec::new(),
            genomic: genomic.map(|s| crate::parse_hgvs(s).unwrap()),
            coding: None,
            noncoding: None,
            protein: protein.map(|s| crate::parse_hgvs(s).unwrap()),
            rna: None,
            transcript_id: "NM_000088.3".to_string(),
            gene_symbol: None,
            is_frameshift: false,
            is_intronic: false,
            is_utr: false,
            affects_init: false,
            axis_decline_reasons: AxisDeclineReasons::default(),
        }
    }

    /// #1182: `select_axis` must carry the projection's normalization warnings
    /// onto the outcome. They were previously read by nobody, which is what made
    /// W4004 structurally unreachable from `ferro project`.
    #[test]
    fn select_axis_carries_warnings_onto_a_rendered_outcome() {
        let mut proj = projection_with(Some("NC_000017.11:g.50198003C>A"), None);
        proj.normalization_warnings = vec![NormalizationWarning::PositionPastEnd {
            accession: "NM_000088.3".to_string(),
            coordinate_system: "c".to_string(),
            position: "c.-238".to_string(),
            bound_kind: "5utr-start".to_string(),
            bound_value: 237,
        }];
        match select_axis(&proj, Axis::Genomic) {
            AxisOutcome::Rendered { warnings, .. } => {
                assert_eq!(warnings.len(), 1, "the warning must survive selection");
                assert_eq!(warnings[0].code, "POSITION_PAST_END");
            }
            other => panic!("expected Rendered, got {other:?}"),
        }
    }

    /// The unavailable arm matters as much: an axis can be unavailable *because*
    /// of what the warning describes, so dropping it there would discard the
    /// explanation exactly when it is most useful.
    #[test]
    fn select_axis_carries_warnings_onto_an_unavailable_outcome() {
        let mut proj = projection_with(None, None);
        proj.normalization_warnings = vec![NormalizationWarning::PositionPastEnd {
            accession: "NM_000088.3".to_string(),
            coordinate_system: "c".to_string(),
            position: "c.-238".to_string(),
            bound_kind: "5utr-start".to_string(),
            bound_value: 237,
        }];
        match select_axis(&proj, Axis::Genomic) {
            AxisOutcome::Unavailable { warnings, .. } => {
                assert_eq!(warnings.len(), 1);
                assert_eq!(warnings[0].code, "POSITION_PAST_END");
            }
            other => panic!("expected Unavailable, got {other:?}"),
        }
    }

    /// A clean projection must not grow a spurious diagnostic.
    #[test]
    fn select_axis_carries_no_warnings_for_a_clean_projection() {
        let proj = projection_with(Some("NC_000017.11:g.50198003C>A"), None);
        match select_axis(&proj, Axis::Genomic) {
            AxisOutcome::Rendered { warnings, .. } => assert!(warnings.is_empty()),
            other => panic!("expected Rendered, got {other:?}"),
        }
    }

    /// #1198: when the engine recorded *why* an axis declined, that explanation
    /// must reach the outcome — verbatim, and *instead of* the string
    /// synthesized from the axis code alone. Decorating it with the synthesized
    /// string would put the axis and its unavailability into the reason twice
    /// over, since both output formats already state each of them.
    #[test]
    fn select_axis_surfaces_a_recorded_decline_reason_verbatim() {
        let engine_explanation =
            "transcript position n.0 (start) is outside NM_000088.3: `n.` numbering starts at 1";
        let mut proj = projection_with(None, None);
        proj.axis_decline_reasons.noncoding = Some(engine_explanation.to_string());
        match select_axis(&proj, Axis::Noncoding) {
            AxisOutcome::Unavailable { reason, .. } => {
                assert!(
                    !reason.contains("representation for this variant"),
                    "the synthesized string must be replaced, not decorated; got {reason:?}"
                );
                assert_eq!(reason, engine_explanation);
            }
            other => panic!("expected Unavailable, got {other:?}"),
        }
    }

    /// A reason recorded for one axis must not be reported against another —
    /// the whole point is that the explanation is *this* axis's, so borrowing
    /// it would be a new way of lying about the cause.
    #[test]
    fn select_axis_does_not_borrow_another_axis_decline_reason() {
        let mut proj = projection_with(None, None);
        proj.axis_decline_reasons.noncoding = Some("the n. axis's own problem".to_string());
        match select_axis(&proj, Axis::Protein) {
            AxisOutcome::Unavailable { reason, .. } => {
                assert_eq!(reason, "no p. representation for this variant");
            }
            other => panic!("expected Unavailable, got {other:?}"),
        }
    }

    /// With nothing recorded the wording is unchanged: most absent axes were
    /// never applicable to the input, and inventing a cause would be worse than
    /// saying nothing.
    #[test]
    fn select_axis_keeps_the_generic_reason_when_nothing_was_recorded() {
        let proj = projection_with(None, None);
        match select_axis(&proj, Axis::Noncoding) {
            AxisOutcome::Unavailable { reason, .. } => {
                assert_eq!(reason, "no n. representation for this variant");
            }
            other => panic!("expected Unavailable, got {other:?}"),
        }
    }

    /// The axis value, not the presence of a reason, decides which outcome you
    /// get. Checking the reason first — plausible, since it is the new input to
    /// this function — would report a rendered axis as unavailable.
    #[test]
    fn a_decline_reason_does_not_override_a_rendered_axis() {
        let mut proj = projection_with(Some("NC_000017.11:g.50198003C>A"), None);
        proj.axis_decline_reasons.genomic = Some("stale".to_string());
        assert_eq!(
            select_axis(&proj, Axis::Genomic),
            AxisOutcome::Rendered {
                transcript_id: "NM_000088.3".to_string(),
                output: "NC_000017.11:g.50198003C>A".to_string(),
                warnings: Vec::new(),
            }
        );
    }

    #[test]
    fn select_axis_renders_present_field() {
        let proj = projection_with(Some("NC_000017.11:g.50198003C>A"), None);
        let outcome = select_axis(&proj, Axis::Genomic);
        assert_eq!(
            outcome,
            AxisOutcome::Rendered {
                transcript_id: "NM_000088.3".to_string(),
                output: "NC_000017.11:g.50198003C>A".to_string(),
                warnings: Vec::new(),
            }
        );
    }

    #[test]
    fn select_axis_absent_field_is_unavailable() {
        let proj = projection_with(None, None);
        let outcome = select_axis(&proj, Axis::Protein);
        match outcome {
            AxisOutcome::Unavailable {
                transcript_id,
                reason,
                ..
            } => {
                assert_eq!(transcript_id.as_deref(), Some("NM_000088.3"));
                assert!(reason.contains("protein") || reason.contains("p."));
            }
            other => panic!("expected Unavailable, got {other:?}"),
        }
    }

    /// A single plus-strand coding transcript NM_TEST.1 on NC_000001.11
    /// [1000..1008], CDS = the whole 9-base exon "ATGCGCTAA". The contig is the
    /// chromosome accession (NOT "chr1") because a bare `NC_000001.11:g.` input
    /// resolves its contig from the accession, so the cdot/transcript must be
    /// indexed under the same `NC_000001.11` name for overlap enumeration.
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

        VariantProjector::new(projector, provider)
    }

    #[test]
    fn project_axis_coding_input_renders_protein() {
        let vp = fixture();
        let variant = crate::parse_hgvs("NM_TEST.1:c.4C>A").unwrap();
        let outcome = project_axis(&vp, &variant, Axis::Protein, None).unwrap();
        match outcome {
            AxisOutcome::Rendered {
                transcript_id,
                output,
                ..
            } => {
                assert_eq!(transcript_id, "NM_TEST.1");
                assert!(output.starts_with("NP_TEST.1:p."), "got {output}");
            }
            other => panic!("expected Rendered, got {other:?}"),
        }
    }

    #[test]
    fn project_axis_transcript_mismatch_is_hard_error() {
        let vp = fixture();
        let variant = crate::parse_hgvs("NM_TEST.1:c.4C>A").unwrap();
        let err = project_axis(&vp, &variant, Axis::Protein, Some("NM_OTHER.1")).unwrap_err();
        assert!(
            err.0.contains("NM_OTHER.1") && err.0.to_lowercase().contains("match"),
            "got {}",
            err.0
        );
    }

    #[test]
    fn project_axis_bare_nm_genomic_axis_unavailable() {
        // A bare `NM_` coding description names no genomic parent, so there is
        // no reference for the `g.` axis to be written against → `g.` is None.
        //
        // This comment previously read "has no genome alignment", which is not
        // what the fixture models and not why the axis declines — see
        // `project_axis_bare_nm_genomic_decline_names_cause_and_remedy` below,
        // which pins that this very fixture aligns the transcript to
        // `NC_000001.11` and projects against it happily.
        let vp = fixture();
        let variant = crate::parse_hgvs("NM_TEST.1:c.4C>A").unwrap();
        let outcome = project_axis(&vp, &variant, Axis::Genomic, None).unwrap();
        assert!(
            matches!(outcome, AxisOutcome::Unavailable { .. }),
            "got {outcome:?}"
        );
    }

    /// #1713: the `g.` decline on a bare transcript input must state its own
    /// cause and the remedy, rather than the string synthesized from the axis
    /// code alone.
    ///
    /// The synthesized wording — "no g. representation for this variant" —
    /// asserts a property of the **variant**. That is false here and measurably
    /// misleading: the reference data required to derive the axis is present
    /// (this fixture aligns `NM_TEST.1` to `NC_000001.11`, and the second half
    /// of this test renders that very axis), and what is missing is a genomic
    /// reference in the **description**. #1713 was filed against ferro on the
    /// strength of that wording, reading it as evidence that the CLI could not
    /// reach a capability the library had.
    ///
    /// `project_to_genomic` declines a parentless accession by design — #327
    /// specifies it, and `projector.rs` refuses to synthesize a parent from
    /// cdot — so the decline itself is correct and is not what this pins. What
    /// is owed is that the reason reaches the user, which is a decline the
    /// projector articulates and then discards (`AxisDeclineReasons`).
    #[test]
    fn project_axis_bare_nm_genomic_decline_names_cause_and_remedy() {
        let vp = fixture();

        let bare = crate::parse_hgvs("NM_TEST.1:c.4C>A").unwrap();
        let outcome = project_axis(&vp, &bare, Axis::Genomic, None).unwrap();
        let AxisOutcome::Unavailable { reason, .. } = outcome else {
            panic!("expected Unavailable, got {outcome:?}");
        };
        assert!(
            !reason.contains("no g. representation for this variant"),
            "the decline must not fall back to the axis-code string, which \
             states a property of the variant that is false here: {reason}"
        );
        assert!(
            reason.contains("NM_TEST.1"),
            "the decline must name the accession it is about: {reason}"
        );
        assert!(
            reason.contains("genomic reference"),
            "the decline must name the missing thing: {reason}"
        );
        assert!(
            reason.contains("NC_000001.11(NM_TEST.1)"),
            "the decline must show the remedy on this transcript's own genomic \
             reference, which the projector can see: {reason}"
        );

        // The other half of the claim, in the same test so the two cannot drift:
        // name the genomic reference and the identical variant renders. Nothing
        // about the reference data changed between these two calls, which is
        // what makes "no g. representation for this variant" the wrong story.
        let anchored = crate::parse_hgvs("NC_000001.11(NM_TEST.1):c.4C>A").unwrap();
        let outcome = project_axis(&vp, &anchored, Axis::Genomic, None).unwrap();
        match outcome {
            AxisOutcome::Rendered { output, .. } => {
                assert!(output.starts_with("NC_000001.11:g."), "got {output}");
            }
            other => panic!("expected Rendered once the parent is named, got {other:?}"),
        }
    }

    /// The remedy is spelled on the axis the user wrote, not always `c.`.
    ///
    /// `project_coding_direct` passes a literal `"c."`; `project_noncoding_direct`
    /// threads its own `axis_label` (`"n."` / `"r."`). That parameter is the only
    /// behavioural difference between the two paths' declines, and the `c.`-only
    /// tests beside this one cannot see it — a `project_noncoding_direct` that
    /// hardcoded `"c."` would pass every other assertion in the tree. So the two
    /// non-coding axes are pinned here, and pinned as *not* `c.`, since a
    /// substring check for the accession alone would survive the regression.
    #[test]
    fn project_axis_bare_nm_genomic_decline_spells_the_remedy_on_the_inputs_own_axis() {
        let vp = fixture();

        for (input, axis_label) in [("NM_TEST.1:n.4C>A", "n."), ("NM_TEST.1:r.4c>a", "r.")] {
            let variant = crate::parse_hgvs(input).unwrap();
            let outcome = project_axis(&vp, &variant, Axis::Genomic, None).unwrap();
            let AxisOutcome::Unavailable { reason, .. } = outcome else {
                panic!("{input}: expected Unavailable, got {outcome:?}");
            };
            assert!(
                reason.contains(&format!("NC_000001.11(NM_TEST.1):{axis_label}")),
                "{input}: the remedy must be spelled on the input's own axis \
                 ({axis_label}): {reason}"
            );
            assert!(
                !reason.contains("NC_000001.11(NM_TEST.1):c."),
                "{input}: the remedy must not fall back to the coding path's \
                 hardcoded `c.`: {reason}"
            );
        }
    }

    #[test]
    fn project_axis_coding_input_n_axis_renders() {
        let vp = fixture();
        let variant = crate::parse_hgvs("NM_TEST.1:c.4C>A").unwrap();
        let outcome = project_axis(&vp, &variant, Axis::Noncoding, None).unwrap();
        assert!(
            matches!(outcome, AxisOutcome::Rendered { .. }),
            "got {outcome:?}"
        );
    }

    #[test]
    fn project_axis_bare_genomic_without_transcript_is_hard_error_listing_overlaps() {
        let vp = fixture();
        let variant = crate::parse_hgvs("NC_000001.11:g.1003C>A").unwrap();
        let err = project_axis(&vp, &variant, Axis::Coding, None).unwrap_err();
        assert!(err.0.contains("requires --transcript"), "got {}", err.0);
        assert!(
            err.0.contains("NM_TEST.1"),
            "should list the overlapping transcript: {}",
            err.0
        );
    }

    #[test]
    fn project_axis_bare_genomic_with_transcript_renders() {
        let vp = fixture();
        let variant = crate::parse_hgvs("NC_000001.11:g.1003C>A").unwrap();
        let outcome = project_axis(&vp, &variant, Axis::Coding, Some("NM_TEST.1")).unwrap();
        match outcome {
            AxisOutcome::Rendered {
                transcript_id,
                output,
                ..
            } => {
                assert_eq!(transcript_id, "NM_TEST.1");
                assert!(output.contains(":c."), "got {output}");
            }
            other => panic!("expected Rendered, got {other:?}"),
        }
    }

    /// #1562: a genomic (`g.`) *cis allele* must project onto the named
    /// transcript exactly as a lone genomic variant does. Before the fix the CLI
    /// only routed `HgvsVariant::Genome(_)` through the genomic dispatch; a
    /// genomic allele parses to `HgvsVariant::Allele`, so it fell through to the
    /// transcript-coordinate branch, which read the input's *chromosome*
    /// accession as its "transcript" and rejected `--transcript` for not
    /// matching it — refusing a valid projection outright.
    ///
    /// The assertion is on the projected *member list*, not merely that the call
    /// succeeds: the whole point of projecting a cis allele is that every member
    /// is carried through onto the target axis. Two substitutions six bases apart
    /// stay individual (`general.md:34` — separated by more than one unchanged
    /// nucleotide), so a projection that dropped or merged a member would fail
    /// this. The exact string is measured, not inferred.
    #[test]
    fn project_axis_genomic_cis_allele_projects_all_members() {
        let vp = fixture();
        let variant = crate::parse_hgvs("NC_000001.11:g.[1002T>A;1008A>C]").unwrap();
        let outcome = project_axis(&vp, &variant, Axis::Coding, Some("NM_TEST.1")).unwrap();
        match outcome {
            AxisOutcome::Rendered {
                transcript_id,
                output,
                ..
            } => {
                assert_eq!(transcript_id, "NM_TEST.1");
                assert_eq!(output, "NC_000001.11(NM_TEST.1):c.[2T>A;9A>C]");
            }
            other => panic!("expected Rendered, got {other:?}"),
        }
    }

    /// The discriminating half of #1562: routing genomic alleles through the
    /// genomic dispatch must NOT weaken the mismatch guard for a genuine
    /// *transcript-coordinate* cis allele. A `c.` allele names its transcript in
    /// the accession, so a `--transcript` that disagrees with it is still a hard
    /// error — exactly as for a lone `c.` input.
    #[test]
    fn project_axis_coding_cis_allele_transcript_mismatch_is_hard_error() {
        let vp = fixture();
        let variant = crate::parse_hgvs("NM_TEST.1:c.[3dup;6del]").unwrap();
        let err = project_axis(&vp, &variant, Axis::Protein, Some("NM_OTHER.1")).unwrap_err();
        assert!(
            err.0.contains("NM_OTHER.1") && err.0.to_lowercase().contains("match"),
            "got {}",
            err.0
        );
    }
}
