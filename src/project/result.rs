//! `VariantProjection` result type.

use crate::hgvs::variant::HgvsVariant;
use crate::normalize::NormalizationWarning;

/// A variant resolved across coordinate systems.
///
/// Each coordinate axis is independently optional: a projection carries
/// whatever representations are derivable from the input and the available
/// reference data, and `None` means "this representation is not available"
/// (not an error). In particular `genomic` is `None` for a bare transcript
/// (e.g. bare-`NM_`) coding input that carries no genome alignment — protein
/// is still predicted directly from the transcript's CDS (#498). For a
/// genome-anchored input `genomic` is the canonical g. representation.
#[derive(Debug, Clone)]
pub struct VariantProjection {
    pub genomic: Option<HgvsVariant>,
    pub coding: Option<HgvsVariant>,
    /// The `n.` (transcript-relative) representation on this transcript, when
    /// derivable. Populated for both coding transcripts (derived genome-free
    /// from the `c.` form via CDS-offset arithmetic) and non-coding transcripts
    /// (the same `Tx` form also carried in `coding`). `None` when no transcript
    /// coordinate is available (e.g. an empty allele).
    ///
    /// Note: for a non-coding (`NR_`) transcript `coding` *also* holds this
    /// `Tx` form — a pre-existing quirk of the `coding` field's name. This axis
    /// gives callers a single field that always means "the `n.` form".
    pub noncoding: Option<HgvsVariant>,
    pub protein: Option<HgvsVariant>,
    /// The predicted RNA consequence `r.(…)` on this transcript, when derivable.
    /// CDS-relative numbering (matches `c.`); `None` when not representable
    /// (no transcript sequence, non-c./n. input, or an unresolved payload).
    pub rna: Option<HgvsVariant>,
    pub transcript_id: String,
    pub gene_symbol: Option<String>,
    pub is_frameshift: bool,
    pub is_intronic: bool,
    pub is_utr: bool,
    /// Whether this projection's edit reaches the translation initiation codon
    /// (CDS 1–3, exonic and non-offset). Set on the single-variant coding path
    /// via `edit_reaches_initiation_codon`; the allele-protein path reads it
    /// to collapse a cis allele's whole-protein consequence to the
    /// initiation-codon-unknown form (`p.(Met1?)` / `p.?`) when a member
    /// disrupts the start codon (#771 follow-on). On an aggregate (allele)
    /// projection this is the OR of its members — `true` when any member's edit
    /// reaches the initiation codon — and that value propagates outward so an
    /// enclosing nested allele's cis-collapse check can read it. `false` for
    /// non-coding and empty-allele projections, and for any aggregate whose
    /// members all leave the start codon untouched.
    pub affects_init: bool,
    /// Warnings emitted while normalizing the input for this projection (e.g.
    /// an auto-corrected reference-sequence mismatch, or — once the reduced-
    /// capability warning lands on the normalize path — a no-genome degrade).
    /// Empty when the input normalized cleanly. Populated at the projector's
    /// normalize entry points; deep projection-building helpers default it to
    /// empty and the entry point attaches the captured set.
    pub normalization_warnings: Vec<NormalizationWarning>,
}
