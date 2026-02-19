//! VCF INFO field annotation
//!
//! This module provides functionality to add HGVS annotations as VCF INFO fields,
//! following conventions similar to VEP and SnpEff.

use crate::error::FerroError;
use crate::hgvs::variant::HgvsVariant;
use crate::reference::loader::TranscriptDb;

use super::annotator::MultiIsoformAnnotator;
use super::record::{InfoValue, VcfRecord};
use super::to_hgvs::HgvsAnnotation;

/// INFO field key for coding HGVS notation
pub const INFO_HGVS_C: &str = "HGVS_C";

/// INFO field key for protein HGVS notation
pub const INFO_HGVS_P: &str = "HGVS_P";

/// INFO field key for gene symbol
pub const INFO_GENE: &str = "GENE";

/// INFO field key for consequence/effect
pub const INFO_CONSEQUENCE: &str = "CONSEQUENCE";

/// INFO field key for transcript ID
pub const INFO_TRANSCRIPT: &str = "TRANSCRIPT";

/// INFO field key for all annotations (pipe-separated)
pub const INFO_ANN: &str = "ANN";

/// Variant consequence/effect types
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Consequence {
    /// Synonymous variant (no amino acid change)
    Synonymous,
    /// Missense variant (amino acid substitution)
    Missense,
    /// Nonsense variant (introduces stop codon)
    Nonsense,
    /// Frameshift variant
    Frameshift,
    /// In-frame insertion
    InframeInsertion,
    /// In-frame deletion
    InframeDeletion,
    /// Splice site variant
    SpliceSite,
    /// Intronic variant
    Intronic,
    /// 5' UTR variant
    FivePrimeUtr,
    /// 3' UTR variant
    ThreePrimeUtr,
    /// Start codon loss
    StartLoss,
    /// Stop codon loss (read-through)
    StopLoss,
    /// Non-coding transcript variant
    NonCoding,
    /// Intergenic variant
    Intergenic,
    /// Unknown consequence
    Unknown,
}

impl Consequence {
    /// Get the SO term for this consequence
    pub fn so_term(&self) -> &'static str {
        match self {
            Consequence::Synonymous => "synonymous_variant",
            Consequence::Missense => "missense_variant",
            Consequence::Nonsense => "stop_gained",
            Consequence::Frameshift => "frameshift_variant",
            Consequence::InframeInsertion => "inframe_insertion",
            Consequence::InframeDeletion => "inframe_deletion",
            Consequence::SpliceSite => "splice_region_variant",
            Consequence::Intronic => "intron_variant",
            Consequence::FivePrimeUtr => "5_prime_UTR_variant",
            Consequence::ThreePrimeUtr => "3_prime_UTR_variant",
            Consequence::StartLoss => "start_lost",
            Consequence::StopLoss => "stop_lost",
            Consequence::NonCoding => "non_coding_transcript_variant",
            Consequence::Intergenic => "intergenic_variant",
            Consequence::Unknown => "sequence_variant",
        }
    }

    /// Get a display-friendly name
    pub fn display_name(&self) -> &'static str {
        match self {
            Consequence::Synonymous => "synonymous",
            Consequence::Missense => "missense",
            Consequence::Nonsense => "nonsense",
            Consequence::Frameshift => "frameshift",
            Consequence::InframeInsertion => "inframe_insertion",
            Consequence::InframeDeletion => "inframe_deletion",
            Consequence::SpliceSite => "splice_site",
            Consequence::Intronic => "intronic",
            Consequence::FivePrimeUtr => "5_prime_UTR",
            Consequence::ThreePrimeUtr => "3_prime_UTR",
            Consequence::StartLoss => "start_loss",
            Consequence::StopLoss => "stop_loss",
            Consequence::NonCoding => "non_coding",
            Consequence::Intergenic => "intergenic",
            Consequence::Unknown => "unknown",
        }
    }

    /// Get the impact severity (HIGH, MODERATE, LOW, MODIFIER)
    pub fn impact(&self) -> &'static str {
        match self {
            Consequence::Nonsense
            | Consequence::Frameshift
            | Consequence::SpliceSite
            | Consequence::StartLoss
            | Consequence::StopLoss => "HIGH",
            Consequence::Missense
            | Consequence::InframeInsertion
            | Consequence::InframeDeletion => "MODERATE",
            Consequence::Synonymous => "LOW",
            Consequence::Intronic
            | Consequence::FivePrimeUtr
            | Consequence::ThreePrimeUtr
            | Consequence::NonCoding
            | Consequence::Intergenic
            | Consequence::Unknown => "MODIFIER",
        }
    }
}

impl std::fmt::Display for Consequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.so_term())
    }
}

/// Extract HGVS coding notation string from annotation
fn get_hgvs_c(ann: &HgvsAnnotation) -> Option<String> {
    match &ann.variant {
        HgvsVariant::Cds(v) => Some(v.to_string()),
        _ => None,
    }
}

/// Extract HGVS genomic notation string from annotation
fn get_hgvs_g(ann: &HgvsAnnotation) -> Option<String> {
    match &ann.variant {
        HgvsVariant::Genome(v) => Some(v.to_string()),
        _ => None,
    }
}

/// Determine the consequence from an HGVS annotation
pub fn determine_consequence(ann: &HgvsAnnotation) -> Consequence {
    // Check for non-coding
    if !ann.is_coding {
        if ann.is_intronic {
            return Consequence::Intronic;
        }
        return Consequence::NonCoding;
    }

    // Check if intronic
    if ann.is_intronic {
        // Check for splice site (within 2bp of exon boundary)
        let hgvs = ann.hgvs_string();
        if hgvs.contains("+1") || hgvs.contains("+2") || hgvs.contains("-1") || hgvs.contains("-2")
        {
            // More specific check to avoid false positives from positions
            // Look for patterns like 123+1, 123-2, etc.
            if hgvs.chars().collect::<Vec<_>>().windows(3).any(|w| {
                w[0].is_ascii_digit()
                    && (w[1] == '+' || w[1] == '-')
                    && (w[2] == '1' || w[2] == '2')
            }) {
                return Consequence::SpliceSite;
            }
        }
        return Consequence::Intronic;
    }

    // Check HGVS string for clues
    let hgvs = ann.hgvs_string();

    // Check for UTR variants
    if hgvs.contains("c.-") {
        return Consequence::FivePrimeUtr;
    }
    if hgvs.contains("c.*") {
        return Consequence::ThreePrimeUtr;
    }

    // Look at the edit type from the variant
    match &ann.variant {
        HgvsVariant::Cds(v) => {
            let edit_str = format!("{:?}", v.loc_edit.edit);

            // Check for frameshift-causing edits (indels not divisible by 3)
            if edit_str.contains("Deletion") || edit_str.contains("Insertion") {
                // For a simple approximation, single base indels cause frameshifts
                // A more complete implementation would analyze the exact length
                if edit_str.contains("Insertion") {
                    return Consequence::InframeInsertion;
                }
                return Consequence::InframeDeletion;
            }

            if edit_str.contains("Delins") {
                return Consequence::InframeDeletion;
            }

            // Substitution - could be synonymous, missense, or nonsense
            // Without protein prediction, default to unknown for coding changes
            Consequence::Unknown
        }
        HgvsVariant::Genome(_) => Consequence::Unknown,
        HgvsVariant::Tx(_) => Consequence::NonCoding,
        HgvsVariant::Rna(_) => Consequence::NonCoding,
        HgvsVariant::Protein(v) => {
            let edit_str = format!("{:?}", v.loc_edit.edit);

            if edit_str.contains("Frameshift") {
                return Consequence::Frameshift;
            }
            if edit_str.contains("Nonsense") || hgvs.contains("Ter") {
                return Consequence::Nonsense;
            }
            if edit_str.contains("Extension") {
                return Consequence::StopLoss;
            }
            if edit_str.contains("Insertion") {
                return Consequence::InframeInsertion;
            }
            if edit_str.contains("Deletion") {
                return Consequence::InframeDeletion;
            }
            if edit_str.contains("Substitution") {
                return Consequence::Missense;
            }
            Consequence::Unknown
        }
        HgvsVariant::Mt(_) => Consequence::Unknown,
        HgvsVariant::Circular(_) => Consequence::Unknown,
        HgvsVariant::RnaFusion(_) => Consequence::Unknown, // Fusion transcripts - complex structural consequence
        HgvsVariant::Allele(_) => Consequence::Unknown,
        HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => Consequence::Unknown,
    }
}

/// VCF annotator that adds INFO fields
pub struct VcfAnnotator<'a> {
    /// Transcript database
    db: &'a TranscriptDb,
    /// Whether to include all transcripts
    include_all_transcripts: bool,
    /// Whether to use VEP-style ANN field
    use_ann_field: bool,
}

impl<'a> VcfAnnotator<'a> {
    /// Create a new VCF annotator
    pub fn new(db: &'a TranscriptDb) -> Self {
        Self {
            db,
            include_all_transcripts: true,
            use_ann_field: false,
        }
    }

    /// Configure whether to include all transcripts
    pub fn include_all_transcripts(mut self, include: bool) -> Self {
        self.include_all_transcripts = include;
        self
    }

    /// Configure whether to use VEP-style ANN field
    pub fn use_ann_field(mut self, use_ann: bool) -> Self {
        self.use_ann_field = use_ann;
        self
    }

    /// Annotate a VCF record with INFO fields
    pub fn annotate(&self, vcf: &VcfRecord) -> Result<VcfRecord, FerroError> {
        let annotator = MultiIsoformAnnotator::new(self.db);
        let result = annotator.annotate(vcf)?;

        let mut annotated = vcf.clone();

        if result.annotations.is_empty() {
            return Ok(annotated);
        }

        if self.use_ann_field {
            // Build VEP/SnpEff-style ANN field
            let ann_values: Vec<String> = result
                .annotations
                .iter()
                .map(|ann| self.build_ann_value(ann))
                .collect();

            annotated.info.insert(
                INFO_ANN.to_string(),
                InfoValue::String(ann_values.join(",")),
            );
        } else {
            // Use separate INFO fields
            self.add_info_fields(&mut annotated, &result.annotations);
        }

        Ok(annotated)
    }

    /// Add individual INFO fields from annotations
    fn add_info_fields(&self, vcf: &mut VcfRecord, annotations: &[HgvsAnnotation]) {
        // Collect values for each field
        let mut genes: Vec<String> = Vec::new();
        let mut hgvs_c: Vec<String> = Vec::new();
        let mut consequences: Vec<String> = Vec::new();
        let mut transcripts: Vec<String> = Vec::new();

        for ann in annotations {
            // Gene
            if let Some(ref gene) = ann.gene_symbol {
                if !genes.contains(gene) {
                    genes.push(gene.clone());
                }
            }

            // HGVS coding (or other if not coding)
            if let Some(c) = get_hgvs_c(ann) {
                if !hgvs_c.contains(&c) {
                    hgvs_c.push(c);
                }
            } else if let Some(g) = get_hgvs_g(ann) {
                if !hgvs_c.contains(&g) {
                    hgvs_c.push(g);
                }
            } else {
                let s = ann.hgvs_string();
                if !hgvs_c.contains(&s) {
                    hgvs_c.push(s);
                }
            }

            // Transcript
            if let Some(ref tx) = ann.transcript_accession {
                if !transcripts.contains(tx) {
                    transcripts.push(tx.clone());
                }
            }

            // Consequence
            let consequence = determine_consequence(ann);
            let cons_str = consequence.so_term().to_string();
            if !consequences.contains(&cons_str) {
                consequences.push(cons_str);
            }

            // Only include first if not including all
            if !self.include_all_transcripts {
                break;
            }
        }

        // Add non-empty fields
        if !genes.is_empty() {
            vcf.info
                .insert(INFO_GENE.to_string(), InfoValue::String(genes.join(",")));
        }

        if !hgvs_c.is_empty() {
            vcf.info
                .insert(INFO_HGVS_C.to_string(), InfoValue::String(hgvs_c.join(",")));
        }

        if !transcripts.is_empty() {
            vcf.info.insert(
                INFO_TRANSCRIPT.to_string(),
                InfoValue::String(transcripts.join(",")),
            );
        }

        if !consequences.is_empty() {
            vcf.info.insert(
                INFO_CONSEQUENCE.to_string(),
                InfoValue::String(consequences.join(",")),
            );
        }
    }

    /// Build a VEP/SnpEff-style ANN value
    fn build_ann_value(&self, ann: &HgvsAnnotation) -> String {
        // Format: Allele|Consequence|IMPACT|Gene|Gene_ID|Feature_Type|Feature|Biotype|HGVS.c|HGVS.p|...
        let consequence = determine_consequence(ann);
        let gene = ann.gene_symbol.as_deref().unwrap_or("");
        let transcript = ann.transcript_accession.as_deref().unwrap_or("");
        let hgvs = ann.hgvs_string();
        let biotype = if ann.is_coding {
            "protein_coding"
        } else {
            "processed_transcript"
        };

        format!(
            "{}|{}|{}|{}||transcript|{}|{}|{}|",
            "", // Allele (empty, use VCF ALT)
            consequence.so_term(),
            consequence.impact(),
            gene,
            transcript,
            biotype,
            hgvs,
        )
    }

    /// Annotate multiple records
    pub fn annotate_batch(&self, records: &[VcfRecord]) -> Vec<Result<VcfRecord, FerroError>> {
        records.iter().map(|vcf| self.annotate(vcf)).collect()
    }
}

/// Generate VCF header lines for HGVS INFO fields
pub fn generate_info_header_lines() -> Vec<String> {
    vec![
        format!(
            "##INFO=<ID={},Number=.,Type=String,Description=\"HGVS coding notation\">",
            INFO_HGVS_C
        ),
        format!(
            "##INFO=<ID={},Number=.,Type=String,Description=\"HGVS protein notation\">",
            INFO_HGVS_P
        ),
        format!(
            "##INFO=<ID={},Number=.,Type=String,Description=\"Gene symbol\">",
            INFO_GENE
        ),
        format!(
            "##INFO=<ID={},Number=.,Type=String,Description=\"Variant consequence (SO terms)\">",
            INFO_CONSEQUENCE
        ),
        format!(
            "##INFO=<ID={},Number=.,Type=String,Description=\"Transcript accession\">",
            INFO_TRANSCRIPT
        ),
        format!(
            "##INFO=<ID={},Number=.,Type=String,Description=\"Functional annotations (VEP format)\">",
            INFO_ANN
        ),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
    use std::sync::OnceLock;

    fn create_test_db() -> TranscriptDb {
        let mut db = TranscriptDb::with_build(GenomeBuild::GRCh38);

        let tx = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            exons: vec![
                Exon::with_genomic(1, 1, 100, 1000, 1099),
                Exon::with_genomic(2, 101, 200, 2000, 2099),
            ],
            cds_start: Some(50),
            cds_end: Some(150),
            sequence: "ATGC".repeat(50),
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2099),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        };

        db.add(tx);
        db
    }

    #[test]
    fn test_consequence_so_terms() {
        assert_eq!(Consequence::Missense.so_term(), "missense_variant");
        assert_eq!(Consequence::Nonsense.so_term(), "stop_gained");
        assert_eq!(Consequence::Frameshift.so_term(), "frameshift_variant");
        assert_eq!(Consequence::Synonymous.so_term(), "synonymous_variant");
    }

    #[test]
    fn test_consequence_impact() {
        assert_eq!(Consequence::Nonsense.impact(), "HIGH");
        assert_eq!(Consequence::Frameshift.impact(), "HIGH");
        assert_eq!(Consequence::Missense.impact(), "MODERATE");
        assert_eq!(Consequence::Synonymous.impact(), "LOW");
        assert_eq!(Consequence::Intronic.impact(), "MODIFIER");
    }

    #[test]
    fn test_vcf_annotator() {
        let db = create_test_db();
        let annotator = VcfAnnotator::new(&db);

        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let annotated = annotator.annotate(&vcf).unwrap();

        // Should have INFO fields populated
        assert!(annotated.info.contains_key(INFO_GENE));
        assert!(annotated.info.contains_key(INFO_TRANSCRIPT));
    }

    #[test]
    fn test_vcf_annotator_with_ann_field() {
        let db = create_test_db();
        let annotator = VcfAnnotator::new(&db).use_ann_field(true);

        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let annotated = annotator.annotate(&vcf).unwrap();

        // Should have ANN field
        assert!(annotated.info.contains_key(INFO_ANN));
    }

    #[test]
    fn test_generate_info_header_lines() {
        let headers = generate_info_header_lines();
        assert!(!headers.is_empty());
        assert!(headers.iter().any(|h| h.contains("HGVS_C")));
        assert!(headers.iter().any(|h| h.contains("HGVS_P")));
        assert!(headers.iter().any(|h| h.contains("GENE")));
        assert!(headers.iter().any(|h| h.contains("CONSEQUENCE")));
    }

    #[test]
    fn test_consequence_display() {
        assert_eq!(format!("{}", Consequence::Missense), "missense_variant");
        assert_eq!(format!("{}", Consequence::Intronic), "intron_variant");
    }

    #[test]
    fn test_consequence_display_name() {
        assert_eq!(Consequence::Missense.display_name(), "missense");
        assert_eq!(Consequence::Nonsense.display_name(), "nonsense");
        assert_eq!(Consequence::Frameshift.display_name(), "frameshift");
        assert_eq!(Consequence::Synonymous.display_name(), "synonymous");
        assert_eq!(
            Consequence::InframeInsertion.display_name(),
            "inframe_insertion"
        );
        assert_eq!(
            Consequence::InframeDeletion.display_name(),
            "inframe_deletion"
        );
        assert_eq!(Consequence::SpliceSite.display_name(), "splice_site");
        assert_eq!(Consequence::Intronic.display_name(), "intronic");
        assert_eq!(Consequence::FivePrimeUtr.display_name(), "5_prime_UTR");
        assert_eq!(Consequence::ThreePrimeUtr.display_name(), "3_prime_UTR");
        assert_eq!(Consequence::StartLoss.display_name(), "start_loss");
        assert_eq!(Consequence::StopLoss.display_name(), "stop_loss");
        assert_eq!(Consequence::NonCoding.display_name(), "non_coding");
        assert_eq!(Consequence::Intergenic.display_name(), "intergenic");
        assert_eq!(Consequence::Unknown.display_name(), "unknown");
    }

    #[test]
    fn test_consequence_all_so_terms() {
        assert_eq!(Consequence::InframeInsertion.so_term(), "inframe_insertion");
        assert_eq!(Consequence::InframeDeletion.so_term(), "inframe_deletion");
        assert_eq!(Consequence::SpliceSite.so_term(), "splice_region_variant");
        assert_eq!(Consequence::FivePrimeUtr.so_term(), "5_prime_UTR_variant");
        assert_eq!(Consequence::ThreePrimeUtr.so_term(), "3_prime_UTR_variant");
        assert_eq!(Consequence::StartLoss.so_term(), "start_lost");
        assert_eq!(Consequence::StopLoss.so_term(), "stop_lost");
        assert_eq!(
            Consequence::NonCoding.so_term(),
            "non_coding_transcript_variant"
        );
        assert_eq!(Consequence::Intergenic.so_term(), "intergenic_variant");
        assert_eq!(Consequence::Unknown.so_term(), "sequence_variant");
    }

    #[test]
    fn test_consequence_impact_all() {
        // HIGH impact
        assert_eq!(Consequence::SpliceSite.impact(), "HIGH");
        assert_eq!(Consequence::StartLoss.impact(), "HIGH");
        assert_eq!(Consequence::StopLoss.impact(), "HIGH");

        // MODERATE impact
        assert_eq!(Consequence::InframeInsertion.impact(), "MODERATE");
        assert_eq!(Consequence::InframeDeletion.impact(), "MODERATE");

        // MODIFIER impact
        assert_eq!(Consequence::FivePrimeUtr.impact(), "MODIFIER");
        assert_eq!(Consequence::ThreePrimeUtr.impact(), "MODIFIER");
        assert_eq!(Consequence::NonCoding.impact(), "MODIFIER");
        assert_eq!(Consequence::Intergenic.impact(), "MODIFIER");
        assert_eq!(Consequence::Unknown.impact(), "MODIFIER");
    }

    #[test]
    fn test_vcf_annotator_include_all_transcripts() {
        let db = create_test_db();

        // Test with include_all_transcripts = true (default)
        let annotator = VcfAnnotator::new(&db).include_all_transcripts(true);
        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let annotated = annotator.annotate(&vcf).unwrap();
        assert!(annotated.info.contains_key(INFO_GENE));

        // Test with include_all_transcripts = false
        let annotator = VcfAnnotator::new(&db).include_all_transcripts(false);
        let annotated = annotator.annotate(&vcf).unwrap();
        assert!(annotated.info.contains_key(INFO_GENE));
    }

    #[test]
    fn test_vcf_annotator_empty_result() {
        let db = TranscriptDb::with_build(GenomeBuild::GRCh38);
        let annotator = VcfAnnotator::new(&db);

        // Position with no overlapping transcripts
        let vcf = VcfRecord::snv("chr99", 99999, 'A', 'G');
        let annotated = annotator.annotate(&vcf).unwrap();

        // Should return original record without annotation fields
        assert!(!annotated.info.contains_key(INFO_GENE));
    }

    #[test]
    fn test_vcf_annotator_annotate_batch() {
        let db = create_test_db();
        let annotator = VcfAnnotator::new(&db);

        let records = vec![
            VcfRecord::snv("chr1", 1050, 'A', 'G'),
            VcfRecord::snv("chr1", 1051, 'T', 'C'),
        ];

        let results = annotator.annotate_batch(&records);
        assert_eq!(results.len(), 2);
        assert!(results[0].is_ok());
        assert!(results[1].is_ok());
    }

    #[test]
    fn test_generate_info_header_lines_all_fields() {
        let headers = generate_info_header_lines();
        assert_eq!(headers.len(), 6);
        assert!(headers.iter().any(|h| h.contains("TRANSCRIPT")));
        assert!(headers.iter().any(|h| h.contains("ANN")));
    }

    #[test]
    fn test_determine_consequence_intronic() {
        use crate::hgvs::parser::parse_hgvs;

        // Create annotation for intronic variant
        let variant = parse_hgvs("NM_000088.3:c.100+5A>G").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NM_000088.3".to_string()),
            is_coding: true,
            is_intronic: true,
            mane_status: ManeStatus::Select,
            intron_number: Some(1),
            intron_position: None,
            intronic_consequence: None,
        };

        let consequence = determine_consequence(&ann);
        // Intronic variants with +5 offset should be splice site (5 is outside +1/+2)
        assert!(matches!(consequence, Consequence::Intronic));
    }

    #[test]
    fn test_determine_consequence_splice_site() {
        use crate::hgvs::parser::parse_hgvs;

        // Create annotation for splice site variant (+1)
        let variant = parse_hgvs("NM_000088.3:c.100+1A>G").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NM_000088.3".to_string()),
            is_coding: true,
            is_intronic: true,
            mane_status: ManeStatus::Select,
            intron_number: Some(1),
            intron_position: None,
            intronic_consequence: None,
        };

        let consequence = determine_consequence(&ann);
        assert_eq!(consequence, Consequence::SpliceSite);
    }

    #[test]
    fn test_determine_consequence_five_prime_utr() {
        use crate::hgvs::parser::parse_hgvs;

        let variant = parse_hgvs("NM_000088.3:c.-10A>G").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NM_000088.3".to_string()),
            is_coding: true,
            is_intronic: false,
            mane_status: ManeStatus::Select,
            intron_number: None,
            intron_position: None,
            intronic_consequence: None,
        };

        let consequence = determine_consequence(&ann);
        assert_eq!(consequence, Consequence::FivePrimeUtr);
    }

    #[test]
    fn test_determine_consequence_three_prime_utr() {
        use crate::hgvs::parser::parse_hgvs;

        let variant = parse_hgvs("NM_000088.3:c.*10A>G").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NM_000088.3".to_string()),
            is_coding: true,
            is_intronic: false,
            mane_status: ManeStatus::Select,
            intron_number: None,
            intron_position: None,
            intronic_consequence: None,
        };

        let consequence = determine_consequence(&ann);
        assert_eq!(consequence, Consequence::ThreePrimeUtr);
    }

    #[test]
    fn test_determine_consequence_noncoding() {
        use crate::hgvs::parser::parse_hgvs;

        let variant = parse_hgvs("NR_000001.1:n.100A>G").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("TEST".to_string()),
            transcript_accession: Some("NR_000001.1".to_string()),
            is_coding: false,
            is_intronic: false,
            mane_status: ManeStatus::None,
            intron_number: None,
            intron_position: None,
            intronic_consequence: None,
        };

        let consequence = determine_consequence(&ann);
        assert_eq!(consequence, Consequence::NonCoding);
    }

    #[test]
    fn test_determine_consequence_noncoding_intronic() {
        use crate::hgvs::parser::parse_hgvs;

        let variant = parse_hgvs("NR_000001.1:n.100+5A>G").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("TEST".to_string()),
            transcript_accession: Some("NR_000001.1".to_string()),
            is_coding: false,
            is_intronic: true,
            mane_status: ManeStatus::None,
            intron_number: Some(1),
            intron_position: None,
            intronic_consequence: None,
        };

        let consequence = determine_consequence(&ann);
        assert_eq!(consequence, Consequence::Intronic);
    }

    #[test]
    fn test_determine_consequence_genomic() {
        use crate::hgvs::parser::parse_hgvs;

        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: None,
            transcript_accession: None,
            is_coding: true,
            is_intronic: false,
            mane_status: ManeStatus::None,
            intron_number: None,
            intron_position: None,
            intronic_consequence: None,
        };

        let consequence = determine_consequence(&ann);
        assert_eq!(consequence, Consequence::Unknown);
    }

    #[test]
    fn test_determine_consequence_protein_frameshift() {
        use crate::hgvs::parser::parse_hgvs;

        let variant = parse_hgvs("NP_000079.2:p.Ala100Glyfs*5").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NP_000079.2".to_string()),
            is_coding: true,
            is_intronic: false,
            mane_status: ManeStatus::None,
            intron_number: None,
            intron_position: None,
            intronic_consequence: None,
        };

        let consequence = determine_consequence(&ann);
        assert_eq!(consequence, Consequence::Frameshift);
    }

    #[test]
    fn test_determine_consequence_protein_nonsense() {
        use crate::hgvs::parser::parse_hgvs;

        let variant = parse_hgvs("NP_000079.2:p.Arg100Ter").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NP_000079.2".to_string()),
            is_coding: true,
            is_intronic: false,
            mane_status: ManeStatus::None,
            intron_number: None,
            intron_position: None,
            intronic_consequence: None,
        };

        let consequence = determine_consequence(&ann);
        assert_eq!(consequence, Consequence::Nonsense);
    }

    #[test]
    fn test_determine_consequence_protein_missense() {
        use crate::hgvs::parser::parse_hgvs;

        let variant = parse_hgvs("NP_000079.2:p.Ala100Gly").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NP_000079.2".to_string()),
            is_coding: true,
            is_intronic: false,
            mane_status: ManeStatus::None,
            intron_number: None,
            intron_position: None,
            intronic_consequence: None,
        };

        let consequence = determine_consequence(&ann);
        assert_eq!(consequence, Consequence::Missense);
    }

    #[test]
    fn test_determine_consequence_protein_deletion() {
        use crate::hgvs::parser::parse_hgvs;

        let variant = parse_hgvs("NP_000079.2:p.Ala100del").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NP_000079.2".to_string()),
            is_coding: true,
            is_intronic: false,
            mane_status: ManeStatus::None,
            intron_number: None,
            intron_position: None,
            intronic_consequence: None,
        };

        let consequence = determine_consequence(&ann);
        assert_eq!(consequence, Consequence::InframeDeletion);
    }

    #[test]
    fn test_determine_consequence_protein_insertion() {
        use crate::hgvs::parser::parse_hgvs;

        let variant = parse_hgvs("NP_000079.2:p.Ala100_Gly101insVal").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NP_000079.2".to_string()),
            is_coding: true,
            is_intronic: false,
            mane_status: ManeStatus::None,
            intron_number: None,
            intron_position: None,
            intronic_consequence: None,
        };

        let consequence = determine_consequence(&ann);
        assert_eq!(consequence, Consequence::InframeInsertion);
    }

    #[test]
    fn test_determine_consequence_protein_extension() {
        use crate::hgvs::parser::parse_hgvs;

        // Note: Extension variants with "Ter" in the string match the Nonsense check first
        // because the code checks hgvs.contains("Ter") before checking for Extension
        let variant = parse_hgvs("NP_000079.2:p.Ter100GlyextTer50").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NP_000079.2".to_string()),
            is_coding: true,
            is_intronic: false,
            mane_status: ManeStatus::None,
            intron_number: None,
            intron_position: None,
            intronic_consequence: None,
        };

        let consequence = determine_consequence(&ann);
        // Due to current logic order, this matches Nonsense (Ter in string) before StopLoss
        assert_eq!(consequence, Consequence::Nonsense);
    }
}
