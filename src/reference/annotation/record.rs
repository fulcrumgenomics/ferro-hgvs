//! AnnotationRecord trait + GFF3 and GTF wrappers. See spec §6 Stage 2.

use std::io::Cursor;

use crate::reference::transcript::Strand;
use noodles_gff::feature::record::{Phase, Strand as NoodlesStrand};
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
    /// Look up a single attribute by key.
    #[allow(dead_code)]
    fn attribute(&self, key: &str) -> Option<&str>;
    fn source_line(&self) -> u64;
    /// Consume the record and return its attribute map.
    fn into_attrs(self) -> AttributeMap;
}

/// Map a noodles strand parse result to our `Strand`, falling back to a
/// `BadStrand` error that carries the *actual* column-7 token. Both GFF3 and
/// GTF use the same `NoodlesStrand` type (re-exported from `noodles-gff`), so
/// the helper is generic over the noodles error type.
///
/// Note: in GFF3, `'?'` is a valid strand (folded into `NoodlesStrand::Unknown`),
/// so reaching the `Err` arm means column 7 held something else (e.g. `x`,
/// `++`); reporting `"?"` would be misleading. Recovering the raw token gives
/// the actionable diagnostic.
fn parse_strand<E>(
    nstrand: Result<NoodlesStrand, E>,
    trimmed: &str,
) -> Result<Strand, RecordParseError> {
    match nstrand {
        Ok(NoodlesStrand::Forward) => Ok(Strand::Plus),
        Ok(NoodlesStrand::Reverse) => Ok(Strand::Minus),
        // GFF3 spec: '.' means unstranded, '?' means relevant-but-unknown.
        // noodles-gtf rejects '?' outright, so for GTF we only see None here.
        Ok(NoodlesStrand::None | NoodlesStrand::Unknown) => Ok(Strand::Unknown),
        Err(_) => {
            let raw = trimmed.split('\t').nth(6).unwrap_or("").to_string();
            Err(RecordParseError::BadStrand(raw))
        }
    }
}

/// Map a noodles phase parse result to our `Option<u8>`, applying the lenient
/// fallback contract: out-of-spec phases (e.g. `'5'`) survive parsing so the
/// builder can emit W-LOAD-111 and treat the phase as unavailable. Noodles'
/// `Phase` enum strictly rejects anything outside `{'.', '0', '1', '2'}`, so on
/// `Err` we re-read column 8 of the raw line and accept any `u8`; otherwise
/// `None`.
fn parse_phase<E>(nphase: Option<Result<Phase, E>>, trimmed: &str) -> Option<u8> {
    match nphase {
        None => None,
        Some(Ok(Phase::Zero)) => Some(0u8),
        Some(Ok(Phase::One)) => Some(1u8),
        Some(Ok(Phase::Two)) => Some(2u8),
        Some(Err(_)) => trimmed
            .split('\t')
            .nth(7)
            .and_then(|s| s.parse::<u8>().ok()),
    }
}

/// Decode a GFF3 attribute value as a list of IDs. noodles-gff parses
/// comma-separated values as `Value::Array` and single values as
/// `Value::String`; both shapes flow into the same `Vec<String>` here.
/// Used for both `Parent` and `Derives_from` parent-edge attributes.
fn parse_gff_id_list(value: &noodles_gff::record::attributes::field::Value) -> Vec<String> {
    use noodles_gff::record::attributes::field::Value as GffValue;
    match value {
        GffValue::String(s) => vec![s.to_string()],
        GffValue::Array(arr) => arr.iter().map(|cow| cow.to_string()).collect(),
    }
}

#[derive(Debug, thiserror::Error)]
pub enum RecordParseError {
    #[error("invalid coordinate '{0}' in column {1}")]
    BadCoordinate(String, usize),
    #[error("invalid strand '{0}'")]
    BadStrand(String),
    /// Wraps any error returned by the underlying noodles parser
    /// (`Reader::read_line`, `Line::as_record`, or `Attributes::iter`).
    /// The wrapped string preserves the noodles error message so the
    /// `E-LOAD-001 MalformedRecord` diagnostic carries actionable detail
    /// rather than a generic "parse failed".
    #[error("malformed record: {0}")]
    Malformed(String),
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

        // Parse using noodles-gff. The Reader strips any trailing newline from
        // the buffer, so wrapping `trimmed` (already newline-stripped above) in
        // a Cursor is safe.
        let cursor = Cursor::new(trimmed.as_bytes());
        let mut reader = noodles_gff::io::Reader::new(cursor);
        let mut nline = noodles_gff::Line::default();
        reader
            .read_line(&mut nline)
            .map_err(|e| RecordParseError::Malformed(e.to_string()))?;

        let nrecord = match nline.as_record() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Err(RecordParseError::Malformed(e.to_string())),
            // noodles treated the line as a comment/directive — shouldn't happen
            // because we already guarded against '#' above.
            None => return Ok(None),
        };

        let seqid = nrecord.reference_sequence_name().to_string();
        let feature_type = FeatureType::from_so_term(&nrecord.ty().to_string());

        let start = nrecord
            .start()
            .map(|p| usize::from(p) as u64)
            .map_err(|_| RecordParseError::BadCoordinate("start".into(), 4))?;

        let end = nrecord
            .end()
            .map(|p| usize::from(p) as u64)
            .map_err(|_| RecordParseError::BadCoordinate("end".into(), 5))?;

        let strand = parse_strand(nrecord.strand(), trimmed)?;
        let phase = parse_phase(nrecord.phase(), trimmed);

        // Extract ID, Parent (+ FlyBase `Derives_from`), and remaining
        // attributes. `Derives_from` indicates a derives-from edge rather
        // than the standard part-of relationship; for ferro's feature
        // graph both establish a parent edge, so the two are merged into
        // a single parent list (dedup'd, Parent values first, Derives_from
        // values appended). #397 item 1.
        let mut id: Option<String> = None;
        let mut parents: Vec<String> = Vec::new();
        let mut derives_from: Vec<String> = Vec::new();
        let mut attrs: AttributeMap = Vec::new();

        for result in nrecord.attributes().iter() {
            let (tag, value) = result.map_err(|e| RecordParseError::Malformed(e.to_string()))?;
            let key_str = tag.to_string();
            match key_str.as_str() {
                "ID" => {
                    // ID is always a single string in GFF3.
                    id = Some(value.as_ref().to_string());
                }
                "Parent" => {
                    parents = parse_gff_id_list(&value);
                }
                "Derives_from" => {
                    derives_from = parse_gff_id_list(&value);
                }
                _ => {
                    // value.as_ref() gives &BStr; convert to str for SmolStr.
                    attrs.push((
                        SmolStr::new(&key_str),
                        SmolStr::new(value.as_ref().to_string()),
                    ));
                }
            }
        }

        // Merge `Derives_from` IDs into `parents`, preserving order and
        // de-duplicating against the existing parents.
        for d in derives_from {
            if !parents.contains(&d) {
                parents.push(d);
            }
        }

        Ok(Some(Self {
            seqid,
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

        // Parse using noodles-gtf.
        let cursor = Cursor::new(trimmed.as_bytes());
        let mut reader = noodles_gtf::io::Reader::new(cursor);
        let mut nline = noodles_gtf::Line::default();
        reader
            .read_line(&mut nline)
            .map_err(|e| RecordParseError::Malformed(e.to_string()))?;

        let nrecord = match nline.as_record() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Err(RecordParseError::Malformed(e.to_string())),
            None => return Ok(None),
        };

        let seqid = nrecord.reference_sequence_name().to_string();
        let feature_type = FeatureType::from_so_term(&nrecord.ty().to_string());

        let start = nrecord
            .start()
            .map(|p| usize::from(p) as u64)
            .map_err(|_| RecordParseError::BadCoordinate("start".into(), 4))?;

        let end = nrecord
            .end()
            .map(|p| usize::from(p) as u64)
            .map_err(|_| RecordParseError::BadCoordinate("end".into(), 5))?;

        // GTF strand: noodles maps '.' to Strand::None, '+' to Forward, '-' to Reverse.
        // GTF does not define '?'; noodles-gtf returns an error for anything
        // other than '.', '+', '-'.
        let strand = parse_strand(nrecord.strand(), trimmed)?;
        let phase = parse_phase(nrecord.phase(), trimmed);

        // Build the attribute map using noodles-gtf's parser.
        let nattrs = nrecord
            .attributes()
            .map_err(|e| RecordParseError::Malformed(e.to_string()))?;

        let mut attrs: AttributeMap = Vec::new();
        for result in nattrs.iter() {
            let (key, value) = result.map_err(|e| RecordParseError::Malformed(e.to_string()))?;
            use noodles_gtf::record::attributes::field::Value as GtfValue;
            // Flatten array values to their first element (GTF duplicates are rare).
            let val_str = match value {
                GtfValue::String(s) => s.to_string(),
                GtfValue::Array(parts) => parts
                    .iter()
                    .next()
                    .map(|c| c.to_string())
                    .unwrap_or_default(),
            };
            attrs.push((SmolStr::new(key.to_string()), SmolStr::new(val_str)));
        }

        let transcript_id = super::feature::attr_get(&attrs, "transcript_id").map(String::from);
        let gene_id = super::feature::attr_get(&attrs, "gene_id").map(String::from);

        let (id, parents) = match &feature_type {
            FeatureType::Gene | FeatureType::PseudoGene => (gene_id.clone(), Vec::new()),
            ft if ft.is_transcript_like() => (transcript_id.clone(), gene_id.into_iter().collect()),
            _ => (None, transcript_id.into_iter().collect()),
        };

        Ok(Some(Self {
            seqid,
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
    fn gff3_record_honors_derives_from_as_parent_edge() {
        // FlyBase emits `Derives_from` for transcripts derived from
        // pseudogenes / alleles. Treat it as an additional parent edge
        // so the feature graph captures the part-of relationship. #397
        // item 1.
        let line = "2L\tFlyBase\tmRNA\t100\t200\t.\t+\t.\tID=FBtr0000001;Derives_from=FBgn0000001";
        let rec = Gff3Record::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.parents(), &["FBgn0000001"]);
    }

    #[test]
    fn gff3_record_merges_parent_and_derives_from() {
        // Both attributes present — parent set is the union (Parent
        // values first, Derives_from values appended). De-duplicate
        // identical IDs.
        let line = "2L\tFlyBase\tmRNA\t100\t200\t.\t+\t.\tID=FBtr0000001;Parent=FBgn0000001;Derives_from=FBgn0000002";
        let rec = Gff3Record::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.parents(), &["FBgn0000001", "FBgn0000002"]);
    }

    #[test]
    fn gff3_record_handles_multi_value_derives_from() {
        // `Derives_from=A,B` (comma-separated multi-value) — both
        // values flow into the parent set.
        let line = "2L\tFlyBase\tmRNA\t100\t200\t.\t+\t.\tID=FBtr0000001;Derives_from=FBgn0000001,FBgn0000002";
        let rec = Gff3Record::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.parents(), &["FBgn0000001", "FBgn0000002"]);
    }

    #[test]
    fn gff3_record_dedups_identical_parent_and_derives_from() {
        // If a record lists the same ID under both `Parent` and
        // `Derives_from`, only one parent edge should land in the
        // graph.
        let line = "2L\tFlyBase\tmRNA\t100\t200\t.\t+\t.\tID=FBtr0000001;Parent=FBgn0000001;Derives_from=FBgn0000001";
        let rec = Gff3Record::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.parents(), &["FBgn0000001"]);
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
    fn gff3_record_bad_strand_reports_actual_value() {
        // Invalid strand 'x' — error should report 'x', not the misleading '?'
        // (since '?' is itself a valid GFF3 strand).
        let line = "chr1\t.\texon\t100\t200\t.\tx\t.\tID=ex1";
        match Gff3Record::parse(line, 1) {
            Err(RecordParseError::BadStrand(v)) => assert_eq!(v, "x"),
            other => panic!("expected BadStrand(\"x\"), got {:?}", other),
        }
    }

    #[test]
    fn gff3_record_lenient_phase_fallback() {
        // CDS row with out-of-spec phase '5' must reach the builder. Noodles'
        // strict Phase enum rejects '5', so the lenient fallback re-reads
        // column 8 and yields Some(5).
        let line = "chr1\t.\tCDS\t100\t200\t.\t+\t5\tID=cds1";
        let rec = Gff3Record::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.phase(), Some(5));
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

    #[test]
    fn gtf_record_bad_strand_reports_actual_value() {
        // Invalid strand 'x' — error should report 'x', not '?'.
        let line = "chr1\tHAVANA\texon\t100\t200\t.\tx\t.\tgene_id \"g1\"; transcript_id \"tx1\";";
        match GtfRecord::parse(line, 1) {
            Err(RecordParseError::BadStrand(v)) => assert_eq!(v, "x"),
            other => panic!("expected BadStrand(\"x\"), got {:?}", other),
        }
    }

    #[test]
    fn gtf_record_lenient_phase_fallback() {
        // Out-of-spec phase '5' must reach the builder via the lenient fallback.
        let line = "chr1\tHAVANA\tCDS\t100\t200\t.\t+\t5\tgene_id \"g1\"; transcript_id \"tx1\";";
        let rec = GtfRecord::parse(line, 1)
            .expect("parse")
            .expect("not comment");
        assert_eq!(rec.phase(), Some(5));
    }
}
