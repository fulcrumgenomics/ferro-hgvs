//! `VariantProjection` result type.

use crate::hgvs::variant::HgvsVariant;

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
    pub transcript_id: String,
    pub gene_symbol: Option<String>,
    pub is_frameshift: bool,
    pub is_intronic: bool,
    pub is_utr: bool,
}
