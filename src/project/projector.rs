//! `VariantProjector` orchestrator.

use crate::data::projection::Projector;
use crate::error::FerroError;
use crate::hgvs::variant::HgvsVariant;
use crate::normalize::NormalizeConfig;
use crate::project::result::VariantProjection;
use crate::reference::ReferenceProvider;

#[allow(dead_code)] // TODO(issue-200): fields used in Task 7
pub struct VariantProjector<P: ReferenceProvider + Clone> {
    projector: Projector,
    provider: P,
    normalize_config: NormalizeConfig,
}

impl<P: ReferenceProvider + Clone> VariantProjector<P> {
    pub fn new(projector: Projector, provider: P) -> Self {
        Self {
            projector,
            provider,
            normalize_config: NormalizeConfig::default(),
        }
    }

    pub fn with_normalize_config(mut self, config: NormalizeConfig) -> Self {
        self.normalize_config = config;
        self
    }

    /// Project an HGVS variant onto a transcript. (Stub — implemented in later task.)
    pub fn project_variant(
        &self,
        _variant: &HgvsVariant,
        _transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        todo!("project_variant not yet implemented")
    }

    /// Parse, normalize, and project an HGVS string onto a transcript.
    pub fn project(
        &self,
        _hgvs_string: &str,
        _transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        todo!("project not yet implemented")
    }
}
