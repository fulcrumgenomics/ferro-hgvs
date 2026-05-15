//! Typed feature representation. See spec §5, §6 Stage 3.

use smol_str::SmolStr;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FeatureType {
    Gene,
    PseudoGene,
    Mrna,
    Transcript,
    NcRna,
    LncRna,
    MiRna,
    TRna,
    RRna,
    PrimaryTranscript,
    Exon,
    Cds,
    StartCodon,
    StopCodon,
    FivePrimeUtr,
    ThreePrimeUtr,
    Utr,
    Other(SmolStr),
}

impl FeatureType {
    pub fn from_so_term(s: &str) -> Self {
        match s {
            "gene" => Self::Gene,
            "pseudogene" => Self::PseudoGene,
            "mRNA" | "messenger_RNA" | "protein_coding_transcript" => Self::Mrna,
            "transcript" => Self::Transcript,
            "ncRNA" | "non_coding_transcript" => Self::NcRna,
            "lncRNA" | "lnc_RNA" | "lincRNA" | "long_noncoding_RNA" => Self::LncRna,
            "miRNA" => Self::MiRna,
            "tRNA" => Self::TRna,
            "rRNA" => Self::RRna,
            "primary_transcript" => Self::PrimaryTranscript,
            "exon" => Self::Exon,
            "CDS" => Self::Cds,
            "start_codon" => Self::StartCodon,
            "stop_codon" => Self::StopCodon,
            "five_prime_UTR" | "5'UTR" | "5UTR" => Self::FivePrimeUtr,
            "three_prime_UTR" | "3'UTR" | "3UTR" => Self::ThreePrimeUtr,
            "UTR" => Self::Utr,
            other => Self::Other(SmolStr::new(other)),
        }
    }

    pub fn is_transcript_like(&self) -> bool {
        matches!(
            self,
            Self::Mrna
                | Self::Transcript
                | Self::NcRna
                | Self::LncRna
                | Self::MiRna
                | Self::TRna
                | Self::RRna
                | Self::PrimaryTranscript
        )
    }

    pub fn is_gene_like(&self) -> bool {
        matches!(self, Self::Gene | Self::PseudoGene)
    }
}

pub type AttributeMap = Vec<(SmolStr, SmolStr)>;

pub fn attr_get<'a>(attrs: &'a AttributeMap, key: &str) -> Option<&'a str> {
    attrs
        .iter()
        .find(|(k, _)| k.eq_ignore_ascii_case(key))
        .map(|(_, v)| v.as_str())
}

#[derive(Debug, Clone)]
pub struct Feature {
    pub id: Option<SmolStr>,
    pub parent_ids: Vec<SmolStr>,
    pub feature_type: FeatureType,
    pub seqid: SmolStr,
    pub range: (u64, u64),
    pub strand: crate::reference::transcript::Strand,
    /// CDS phase (0/1/2); populated by the parser but consumed in Phase 2+.
    #[allow(dead_code)]
    pub phase: Option<u8>,
    /// Raw key-value attributes; populated from the parsed record.
    pub attrs: AttributeMap,
    pub source_line: u64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn so_term_maps_to_known_variants() {
        assert_eq!(FeatureType::from_so_term("mRNA"), FeatureType::Mrna);
        assert_eq!(
            FeatureType::from_so_term("messenger_RNA"),
            FeatureType::Mrna
        );
        assert_eq!(
            FeatureType::from_so_term("transcript"),
            FeatureType::Transcript
        );
        assert_eq!(FeatureType::from_so_term("CDS"), FeatureType::Cds);
        assert_eq!(FeatureType::from_so_term("exon"), FeatureType::Exon);
        assert_eq!(
            FeatureType::from_so_term("start_codon"),
            FeatureType::StartCodon
        );
        assert_eq!(
            FeatureType::from_so_term("stop_codon"),
            FeatureType::StopCodon
        );
        assert_eq!(
            FeatureType::from_so_term("five_prime_UTR"),
            FeatureType::FivePrimeUtr
        );
        assert_eq!(
            FeatureType::from_so_term("3'UTR"),
            FeatureType::ThreePrimeUtr
        );
        assert_eq!(FeatureType::from_so_term("gene"), FeatureType::Gene);
        assert_eq!(FeatureType::from_so_term("lncRNA"), FeatureType::LncRna);
    }

    #[test]
    fn unknown_so_term_preserves_string() {
        match FeatureType::from_so_term("some_weird_thing") {
            FeatureType::Other(s) => assert_eq!(s, "some_weird_thing"),
            _ => panic!("expected Other"),
        }
    }

    #[test]
    fn transcript_like_predicates() {
        assert!(FeatureType::Mrna.is_transcript_like());
        assert!(FeatureType::Transcript.is_transcript_like());
        assert!(FeatureType::NcRna.is_transcript_like());
        assert!(FeatureType::LncRna.is_transcript_like());
        assert!(!FeatureType::Gene.is_transcript_like());
        assert!(!FeatureType::Exon.is_transcript_like());
        assert!(!FeatureType::Cds.is_transcript_like());
    }

    #[test]
    fn attr_get_is_case_insensitive() {
        let attrs: AttributeMap = vec![(
            smol_str::SmolStr::new("Gene"),
            smol_str::SmolStr::new("BRCA1"),
        )];
        assert_eq!(attr_get(&attrs, "gene"), Some("BRCA1"));
        assert_eq!(attr_get(&attrs, "GENE"), Some("BRCA1"));
        assert_eq!(attr_get(&attrs, "missing"), None);
    }
}
