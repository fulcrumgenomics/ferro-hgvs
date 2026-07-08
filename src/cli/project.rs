//! CLI-agnostic core for the `ferro project` subcommand: select one output
//! axis from a [`VariantProjection`] and classify engine errors into
//! "unavailable" (exit 0) vs "hard failure" (exit 1).

use crate::error::FerroError;
use crate::hgvs::variant::HgvsVariant;
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

/// The result of selecting one axis from a projection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AxisOutcome {
    /// The axis rendered to an HGVS string.
    Rendered {
        transcript_id: String,
        output: String,
    },
    /// The axis is legitimately not available (exit 0); carries a reason and,
    /// when known, the transcript the projection was attempted on.
    Unavailable {
        transcript_id: Option<String>,
        reason: String,
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
    let field = match axis {
        Axis::Genomic => &projection.genomic,
        Axis::Coding => &projection.coding,
        Axis::Noncoding => &projection.noncoding,
        Axis::Protein => &projection.protein,
        Axis::Rna => &projection.rna,
    };
    let tx = projection.transcript_id.clone();
    match field {
        Some(v) => AxisOutcome::Rendered {
            transcript_id: tx,
            output: v.to_string(),
        },
        None => AxisOutcome::Unavailable {
            transcript_id: Some(tx),
            reason: format!("no {}. representation for this variant", axis.code()),
        },
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
    if matches!(variant, HgvsVariant::Genome(_)) {
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
    use crate::project::VariantProjection;
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
        }
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
        // A bare NM_ coding input has no genome alignment → g. is None.
        let vp = fixture();
        let variant = crate::parse_hgvs("NM_TEST.1:c.4C>A").unwrap();
        let outcome = project_axis(&vp, &variant, Axis::Genomic, None).unwrap();
        assert!(
            matches!(outcome, AxisOutcome::Unavailable { .. }),
            "got {outcome:?}"
        );
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
            } => {
                assert_eq!(transcript_id, "NM_TEST.1");
                assert!(output.contains(":c."), "got {output}");
            }
            other => panic!("expected Rendered, got {other:?}"),
        }
    }
}
