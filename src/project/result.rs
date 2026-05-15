//! `VariantProjection` result type.

use crate::hgvs::variant::HgvsVariant;

/// The result of projecting a g. variant onto a transcript.
#[derive(Debug, Clone)]
pub struct VariantProjection {
    pub genomic: HgvsVariant,
    pub coding: Option<HgvsVariant>,
    pub protein: Option<HgvsVariant>,
    pub transcript_id: String,
    pub gene_symbol: Option<String>,
    pub is_frameshift: bool,
    pub is_intronic: bool,
    pub is_utr: bool,
}
