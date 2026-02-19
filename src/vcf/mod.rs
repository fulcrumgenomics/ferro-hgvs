//! VCF (Variant Call Format) support
//!
//! This module provides types and utilities for working with VCF files,
//! including parsing, writing, and converting between VCF and HGVS formats.

mod annotate;
mod annotator;
mod batch;
mod from_hgvs;
mod parser;
mod record;
mod to_hgvs;

pub use annotate::{
    determine_consequence, generate_info_header_lines, Consequence, VcfAnnotator, INFO_ANN,
    INFO_CONSEQUENCE, INFO_GENE, INFO_HGVS_C, INFO_HGVS_P, INFO_TRANSCRIPT,
};
pub use annotator::{
    AnnotationConfig, BatchAnnotationResult, MultiIsoformAnnotation, MultiIsoformAnnotator,
};
pub use batch::{BatchConfig, BatchProcessor, BatchStats};
pub use from_hgvs::{genomic_hgvs_to_vcf, HgvsToVcfConverter, VcfConversionResult};
pub use parser::{open_vcf, parse_vcf_string, split_multiallelic, VcfHeader, VcfReader};
pub use record::{InfoValue, VcfRecord};
pub use to_hgvs::{vcf_to_genomic_hgvs, HgvsAnnotation, VcfToHgvsConverter};
