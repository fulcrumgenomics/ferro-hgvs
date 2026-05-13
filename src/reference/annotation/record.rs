//! AnnotationRecord trait + GFF3 and GTF wrappers. See spec §6 Stage 2.

use crate::reference::transcript::Strand;
use smol_str::SmolStr;

use super::feature::{AttributeMap, FeatureType};

pub trait AnnotationRecord: Sized {
    /// Parse one line of the underlying format. Returns `Ok(None)` for
    /// comments / blank lines so the dispatch loop can skip them without
    /// allocating a diagnostic.
    fn parse(line: &str, source_line: u64) -> Result<Option<Self>, RecordParseError>;
    fn seqid(&self) -> &str;
    fn feature_type(&self) -> &FeatureType;
    fn start(&self) -> u64;
    fn end(&self) -> u64;
    fn strand(&self) -> Strand;
    fn phase(&self) -> Option<u8>;
    fn id(&self) -> Option<&str>;
    fn parents(&self) -> &[String];
    /// Look up a single attribute by key; used in Phase 2.5+.
    #[allow(dead_code)]
    fn attribute(&self, key: &str) -> Option<&str>;
    fn source_line(&self) -> u64;
    /// Consume the record and return its attribute map; used in Phase 2.5+.
    #[allow(dead_code)]
    fn into_attrs(self) -> AttributeMap;
}

#[derive(Debug, thiserror::Error)]
pub enum RecordParseError {
    #[error("too few tab-separated columns: {0}")]
    TooFewColumns(usize),
    #[error("invalid coordinate '{0}' in column {1}")]
    BadCoordinate(String, usize),
    #[error("invalid phase '{0}'")]
    BadPhase(String),
    #[error("invalid strand '{0}'")]
    BadStrand(String),
}

#[derive(Debug, Clone)]
pub struct Gff3Record {
    seqid: String,
    feature_type: FeatureType,
    start: u64,
    end: u64,
    strand: Strand,
    phase: Option<u8>,
    id: Option<String>,
    parents: Vec<String>,
    /// Raw attributes; consumed by `into_attrs` in Phase 2.5+.
    #[allow(dead_code)]
    attrs: AttributeMap,
    source_line: u64,
}

impl Gff3Record {
    /// Parse one GFF3 line. Returns Ok(None) for comments / blank lines.
    pub fn parse(line: &str, source_line: u64) -> Result<Option<Self>, RecordParseError> {
        let trimmed = line.trim_end_matches('\n').trim_end_matches('\r');
        if trimmed.is_empty() || trimmed.starts_with('#') {
            return Ok(None);
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 9 {
            return Err(RecordParseError::TooFewColumns(fields.len()));
        }
        let start = fields[3]
            .parse::<u64>()
            .map_err(|_| RecordParseError::BadCoordinate(fields[3].into(), 4))?;
        let end = fields[4]
            .parse::<u64>()
            .map_err(|_| RecordParseError::BadCoordinate(fields[4].into(), 5))?;
        let strand = match fields[6] {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            "." | "?" => Strand::Unknown,
            other => return Err(RecordParseError::BadStrand(other.into())),
        };
        let phase = match fields[7] {
            "." => None,
            s => Some(
                s.parse::<u8>()
                    .map_err(|_| RecordParseError::BadPhase(s.into()))?,
            ),
        };
        let (id, parents, attrs) = parse_gff3_attrs(fields[8]);
        Ok(Some(Self {
            seqid: fields[0].to_string(),
            feature_type: FeatureType::from_so_term(fields[2]),
            start,
            end,
            strand,
            phase,
            id,
            parents,
            attrs,
            source_line,
        }))
    }
}

fn parse_gff3_attrs(s: &str) -> (Option<String>, Vec<String>, AttributeMap) {
    let mut id = None;
    let mut parents = Vec::new();
    let mut attrs: AttributeMap = Vec::new();
    for part in s.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let (k, v) = match part.split_once('=') {
            Some(kv) => kv,
            None => continue,
        };
        let v = url_decode(v);
        match k {
            "ID" => id = Some(v),
            "Parent" => parents = v.split(',').map(|p| p.to_string()).collect(),
            _ => attrs.push((SmolStr::new(k), SmolStr::new(v))),
        }
    }
    (id, parents, attrs)
}

fn url_decode(s: &str) -> String {
    let bytes = s.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(s.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' && i + 2 < bytes.len() {
            if let Ok(b) =
                u8::from_str_radix(std::str::from_utf8(&bytes[i + 1..i + 3]).unwrap_or(""), 16)
            {
                out.push(b);
                i += 3;
                continue;
            }
        }
        out.push(bytes[i]);
        i += 1;
    }
    String::from_utf8(out).unwrap_or_else(|_| s.to_string())
}

impl AnnotationRecord for Gff3Record {
    fn parse(line: &str, source_line: u64) -> Result<Option<Self>, RecordParseError> {
        Self::parse(line, source_line)
    }
    fn seqid(&self) -> &str {
        &self.seqid
    }
    fn feature_type(&self) -> &FeatureType {
        &self.feature_type
    }
    fn start(&self) -> u64 {
        self.start
    }
    fn end(&self) -> u64 {
        self.end
    }
    fn strand(&self) -> Strand {
        self.strand
    }
    fn phase(&self) -> Option<u8> {
        self.phase
    }
    fn id(&self) -> Option<&str> {
        self.id.as_deref()
    }
    fn parents(&self) -> &[String] {
        &self.parents
    }
    fn attribute(&self, key: &str) -> Option<&str> {
        super::feature::attr_get(&self.attrs, key)
    }
    fn source_line(&self) -> u64 {
        self.source_line
    }
    fn into_attrs(self) -> AttributeMap {
        self.attrs
    }
}

#[derive(Debug, Clone)]
pub struct GtfRecord {
    seqid: String,
    feature_type: FeatureType,
    start: u64,
    end: u64,
    strand: Strand,
    phase: Option<u8>,
    id: Option<String>,
    parents: Vec<String>,
    /// Raw attributes; consumed by `into_attrs` in Phase 2.5+.
    #[allow(dead_code)]
    attrs: AttributeMap,
    source_line: u64,
}

impl GtfRecord {
    /// Parse one GTF line. Returns Ok(None) for comments / blank lines.
    pub fn parse(line: &str, source_line: u64) -> Result<Option<Self>, RecordParseError> {
        let trimmed = line.trim_end_matches('\n').trim_end_matches('\r');
        if trimmed.is_empty() || trimmed.starts_with('#') {
            return Ok(None);
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 9 {
            return Err(RecordParseError::TooFewColumns(fields.len()));
        }
        let start = fields[3]
            .parse::<u64>()
            .map_err(|_| RecordParseError::BadCoordinate(fields[3].into(), 4))?;
        let end = fields[4]
            .parse::<u64>()
            .map_err(|_| RecordParseError::BadCoordinate(fields[4].into(), 5))?;
        let strand = match fields[6] {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            "." | "?" => Strand::Unknown,
            other => return Err(RecordParseError::BadStrand(other.into())),
        };
        let phase = match fields[7] {
            "." => None,
            s => Some(
                s.parse::<u8>()
                    .map_err(|_| RecordParseError::BadPhase(s.into()))?,
            ),
        };
        let attrs = parse_gtf_attrs(fields[8]);
        let feature_type = FeatureType::from_so_term(fields[2]);

        let transcript_id = super::feature::attr_get(&attrs, "transcript_id").map(String::from);
        let gene_id = super::feature::attr_get(&attrs, "gene_id").map(String::from);

        let (id, parents) = match &feature_type {
            FeatureType::Gene | FeatureType::PseudoGene => (gene_id.clone(), Vec::new()),
            ft if ft.is_transcript_like() => (transcript_id.clone(), gene_id.into_iter().collect()),
            _ => (None, transcript_id.into_iter().collect()),
        };
        Ok(Some(Self {
            seqid: fields[0].to_string(),
            feature_type,
            start,
            end,
            strand,
            phase,
            id,
            parents,
            attrs,
            source_line,
        }))
    }
}

fn parse_gtf_attrs(s: &str) -> AttributeMap {
    let mut out = AttributeMap::new();
    for part in s.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let (k, v) = match part.split_once(' ') {
            Some(kv) => kv,
            None => continue,
        };
        let v = v.trim().trim_matches('"');
        out.push((SmolStr::new(k), SmolStr::new(v)));
    }
    out
}

impl AnnotationRecord for GtfRecord {
    fn parse(line: &str, source_line: u64) -> Result<Option<Self>, RecordParseError> {
        Self::parse(line, source_line)
    }
    fn seqid(&self) -> &str {
        &self.seqid
    }
    fn feature_type(&self) -> &FeatureType {
        &self.feature_type
    }
    fn start(&self) -> u64 {
        self.start
    }
    fn end(&self) -> u64 {
        self.end
    }
    fn strand(&self) -> Strand {
        self.strand
    }
    fn phase(&self) -> Option<u8> {
        self.phase
    }
    fn id(&self) -> Option<&str> {
        self.id.as_deref()
    }
    fn parents(&self) -> &[String] {
        &self.parents
    }
    fn attribute(&self, key: &str) -> Option<&str> {
        super::feature::attr_get(&self.attrs, key)
    }
    fn source_line(&self) -> u64 {
        self.source_line
    }
    fn into_attrs(self) -> AttributeMap {
        self.attrs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::Strand;

    #[test]
    fn gff3_record_parses_basic_mrna_line() {
        let line = "chr1\t.\tmRNA\t100\t200\t.\t+\t.\tID=tx1;Parent=gene01;gene=GENE1";
        let rec = Gff3Record::parse(line, 5)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.seqid(), "chr1");
        assert_eq!(rec.feature_type(), &FeatureType::Mrna);
        assert_eq!(rec.start(), 100);
        assert_eq!(rec.end(), 200);
        assert_eq!(rec.strand(), Strand::Plus);
        assert_eq!(rec.id(), Some("tx1"));
        assert_eq!(rec.parents(), &["gene01"]);
        assert_eq!(rec.attribute("gene"), Some("GENE1"));
        assert_eq!(rec.source_line(), 5);
    }

    #[test]
    fn gff3_record_handles_multi_parent() {
        let line = "chr1\t.\texon\t100\t200\t.\t+\t.\tID=ex1;Parent=tx1,tx2";
        let rec = Gff3Record::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.parents(), &["tx1", "tx2"]);
    }

    #[test]
    fn gff3_record_handles_unknown_strand() {
        let line = "chr1\t.\texon\t100\t200\t.\t.\t.\tID=ex1";
        let rec = Gff3Record::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.strand(), Strand::Unknown);
    }

    #[test]
    fn gff3_record_url_decodes_attributes() {
        let line = "chr1\t.\tgene\t100\t200\t.\t+\t.\tID=g1;Name=My%20Gene";
        let rec = Gff3Record::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.attribute("Name"), Some("My Gene"));
    }

    #[test]
    fn gff3_record_url_decodes_multibyte_utf8() {
        // "%C3%A9" is the percent-encoding of the two UTF-8 bytes for "é".
        let line = "chr1\t.\tgene\t100\t200\t.\t+\t.\tID=g1;Name=caf%C3%A9";
        let rec = Gff3Record::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.attribute("Name"), Some("café"));
    }

    #[test]
    fn gff3_record_skips_comments() {
        assert!(Gff3Record::parse("# header", 1).unwrap().is_none());
        assert!(Gff3Record::parse("", 1).unwrap().is_none());
    }

    #[test]
    fn gff3_record_rejects_malformed_coordinate() {
        let line = "chr1\t.\texon\tabc\t200\t.\t+\t.\tID=ex1";
        assert!(Gff3Record::parse(line, 1).is_err());
    }

    #[test]
    fn gtf_record_parses_exon_with_transcript_id() {
        let line = "chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_id \"ENSG1\"; transcript_id \"ENST1\"; exon_number \"1\";";
        let rec = GtfRecord::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.seqid(), "chr1");
        assert_eq!(rec.feature_type(), &FeatureType::Exon);
        assert_eq!(rec.start(), 100);
        assert_eq!(rec.end(), 200);
        assert_eq!(rec.parents(), &["ENST1"]);
        assert_eq!(rec.attribute("gene_id"), Some("ENSG1"));
        assert_eq!(rec.attribute("exon_number"), Some("1"));
    }

    #[test]
    fn gtf_record_transcript_row_uses_self_id() {
        let line = "chr1\tHAVANA\ttranscript\t100\t500\t.\t+\t.\tgene_id \"ENSG1\"; transcript_id \"ENST1\";";
        let rec = GtfRecord::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.id(), Some("ENST1"));
        assert_eq!(rec.parents(), &["ENSG1"]);
    }

    #[test]
    fn gtf_record_gene_row_uses_gene_id() {
        let line = "chr1\tHAVANA\tgene\t100\t500\t.\t+\t.\tgene_id \"ENSG1\";";
        let rec = GtfRecord::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.id(), Some("ENSG1"));
        assert!(rec.parents().is_empty());
    }

    #[test]
    fn gtf_record_parses_minus_strand() {
        let line = "chr1\tHAVANA\texon\t100\t200\t.\t-\t.\tgene_id \"g1\"; transcript_id \"tx1\";";
        let rec = GtfRecord::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.strand(), Strand::Minus);
    }

    #[test]
    fn gtf_record_parses_phase_for_cds() {
        let line = "chr1\tHAVANA\tCDS\t100\t200\t.\t+\t2\tgene_id \"g1\"; transcript_id \"tx1\";";
        let rec = GtfRecord::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.phase(), Some(2));
    }
}
