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
    pub protein: Option<HgvsVariant>,
    pub transcript_id: String,
    pub gene_symbol: Option<String>,
    pub is_frameshift: bool,
    pub is_intronic: bool,
    pub is_utr: bool,
}
