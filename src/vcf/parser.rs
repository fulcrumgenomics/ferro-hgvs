//! VCF file parsing using noodles-vcf
//!
//! This module provides utilities for parsing VCF files into VcfRecord types.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use noodles_vcf as nvcf;
use nvcf::variant::record::{AlternateBases, Filters, Ids};

use crate::error::FerroError;
use crate::reference::transcript::GenomeBuild;

use super::record::{InfoValue, VcfRecord};

/// A parsed VCF file header
#[derive(Debug, Clone)]
pub struct VcfHeader {
    /// Contigs defined in the header (##contig lines)
    pub contigs: Vec<String>,
    /// Sample names from the header line
    pub samples: Vec<String>,
    /// Raw header for reference
    inner: nvcf::Header,
}

impl VcfHeader {
    /// Get the number of samples in the VCF
    pub fn sample_count(&self) -> usize {
        self.samples.len()
    }

    /// Check if a contig is defined in the header
    pub fn has_contig(&self, name: &str) -> bool {
        self.contigs.iter().any(|c| c == name)
    }
}

/// VCF file reader that yields VcfRecord instances
pub struct VcfReader<R> {
    inner: nvcf::io::Reader<R>,
    header: VcfHeader,
    genome_build: GenomeBuild,
}

impl<R: BufRead> VcfReader<R> {
    /// Create a new VCF reader from a buffered reader
    pub fn new(reader: R) -> Result<Self, FerroError> {
        let mut inner = nvcf::io::Reader::new(reader);
        let noodles_header = inner.read_header().map_err(|e| FerroError::Io {
            msg: format!("Failed to parse VCF header: {}", e),
        })?;

        // Extract contig names
        let contigs: Vec<String> = noodles_header
            .contigs()
            .keys()
            .map(|k| k.to_string())
            .collect();

        // Extract sample names
        let samples: Vec<String> = noodles_header
            .sample_names()
            .iter()
            .map(|s| s.to_string())
            .collect();

        let header = VcfHeader {
            contigs,
            samples,
            inner: noodles_header,
        };

        Ok(Self {
            inner,
            header,
            genome_build: GenomeBuild::default(),
        })
    }

    /// Get a reference to the parsed header
    pub fn header(&self) -> &VcfHeader {
        &self.header
    }

    /// Set the genome build for parsed records
    pub fn with_genome_build(mut self, build: GenomeBuild) -> Self {
        self.genome_build = build;
        self
    }

    /// Read the next VCF record
    pub fn read_record(&mut self) -> Result<Option<VcfRecord>, FerroError> {
        let mut record = nvcf::variant::RecordBuf::default();

        match self.inner.read_record_buf(&self.header.inner, &mut record) {
            Ok(0) => Ok(None), // EOF
            Ok(_) => {
                let vcf_record = convert_record(&record, self.genome_build, &self.header)?;
                Ok(Some(vcf_record))
            }
            Err(e) => Err(FerroError::Io {
                msg: format!("Failed to parse VCF record: {}", e),
            }),
        }
    }

    /// Iterate over all records in the VCF file
    pub fn records(self) -> VcfRecordIterator<R> {
        VcfRecordIterator {
            reader: self,
            done: false,
        }
    }
}

/// Open a VCF file from a path
pub fn open_vcf<P: AsRef<Path>>(path: P) -> Result<VcfReader<BufReader<File>>, FerroError> {
    let file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
        msg: format!("Failed to open VCF file: {}", e),
    })?;
    VcfReader::new(BufReader::new(file))
}

/// Parse VCF from a string
pub fn parse_vcf_string(vcf_content: &str) -> Result<VcfReader<BufReader<&[u8]>>, FerroError> {
    let reader = BufReader::new(vcf_content.as_bytes());
    VcfReader::new(reader)
}

/// Iterator over VCF records
pub struct VcfRecordIterator<R> {
    reader: VcfReader<R>,
    done: bool,
}

impl<R: BufRead> Iterator for VcfRecordIterator<R> {
    type Item = Result<VcfRecord, FerroError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        match self.reader.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(e) => {
                self.done = true;
                Some(Err(e))
            }
        }
    }
}

/// Convert a noodles VCF record to our VcfRecord type
fn convert_record(
    record: &nvcf::variant::RecordBuf,
    genome_build: GenomeBuild,
    header: &VcfHeader,
) -> Result<VcfRecord, FerroError> {
    // Chromosome
    let chrom = record.reference_sequence_name().to_string();

    // Position (noodles uses 1-based positions)
    let pos = record
        .variant_start()
        .map(|p| p.get() as u64)
        .ok_or_else(|| FerroError::Io {
            msg: "Missing position in VCF record".to_string(),
        })?;

    // ID field
    let id = {
        let ids = record.ids();
        if ids.is_empty() {
            None
        } else {
            Some(
                ids.iter()
                    .map(|id| id.to_string())
                    .collect::<Vec<_>>()
                    .join(";"),
            )
        }
    };

    // Reference allele
    let reference = record.reference_bases().to_string();

    // Alternate alleles
    let alternate: Vec<String> = record
        .alternate_bases()
        .iter()
        .map(|a| match a {
            Ok(allele) => allele.to_string(),
            Err(_) => ".".to_string(),
        })
        .collect();

    // Quality score
    let quality = record.quality_score();

    // Filter
    let filter = {
        let filters = record.filters();
        if filters.is_pass() {
            None // PASS
        } else {
            let filter_strs: Vec<_> = filters
                .iter(&header.inner)
                .filter_map(|f| f.ok())
                .map(|filter| filter.to_string())
                .collect();
            if filter_strs.is_empty() {
                None
            } else {
                Some(filter_strs.join(";"))
            }
        }
    };

    // INFO fields - use the RecordBuf's typed info field directly
    let mut info = std::collections::HashMap::new();
    for (key, value) in record.info().as_ref() {
        let key_str = key.to_string();
        if let Some(v) = value {
            let info_value = convert_info_value(v);
            info.insert(key_str, info_value);
        } else {
            info.insert(key_str, InfoValue::Flag);
        }
    }

    // FORMAT keys (e.g. "GT:DP:GQ") and per-sample values.
    let samples_buf = record.samples();
    let format_keys: Vec<String> = samples_buf
        .keys()
        .as_ref()
        .iter()
        .map(|k| k.to_string())
        .collect();

    let format = if format_keys.is_empty() {
        None
    } else {
        Some(format_keys.join(":"))
    };

    // One map per sample column, in sample-column order (matching the header's
    // `sample_names()`), keyed by FORMAT field name. A field whose value is
    // missing (`.`) or absent from a short sample column is omitted from the
    // map; `VcfRecord`'s `Display` renders any absent FORMAT key as `.`, so the
    // omission round-trips correctly and keeps the maps minimal.
    let samples: Vec<std::collections::HashMap<String, String>> = samples_buf
        .values()
        .map(|sample| {
            sample
                .keys()
                .as_ref()
                .iter()
                .zip(sample.values())
                .filter_map(|(key, value)| {
                    value
                        .as_ref()
                        .map(|v| (key.to_string(), format_sample_value(v)))
                })
                .collect()
        })
        .collect();

    Ok(VcfRecord {
        chrom,
        pos,
        id,
        reference,
        alternate,
        quality,
        filter,
        info,
        format,
        samples,
        genome_build,
    })
}

/// Convert a noodles INFO value to our InfoValue type
fn convert_info_value(value: &nvcf::variant::record_buf::info::field::Value) -> InfoValue {
    use nvcf::variant::record_buf::info::field::Value;

    match value {
        Value::Integer(v) => InfoValue::Integer(*v as i64),
        Value::Float(v) => InfoValue::Float(*v as f64),
        Value::Flag => InfoValue::Flag,
        Value::Character(v) => InfoValue::Character(*v),
        Value::String(v) => InfoValue::String(v.clone()),
        Value::Array(arr) => {
            use nvcf::variant::record_buf::info::field::value::Array;
            match arr {
                Array::Integer(vals) => {
                    let ints: Vec<i64> = vals.iter().filter_map(|v| v.map(|i| i as i64)).collect();
                    InfoValue::IntegerArray(ints)
                }
                Array::Float(vals) => {
                    let floats: Vec<f64> =
                        vals.iter().filter_map(|v| v.map(|f| f as f64)).collect();
                    InfoValue::FloatArray(floats)
                }
                Array::Character(vals) => {
                    let chars: Vec<String> = vals
                        .iter()
                        .filter_map(|v| v.map(|c| c.to_string()))
                        .collect();
                    InfoValue::StringArray(chars)
                }
                Array::String(vals) => {
                    let strs: Vec<String> = vals.iter().filter_map(|v| v.clone()).collect();
                    InfoValue::StringArray(strs)
                }
            }
        }
    }
}

/// Format a per-sample FORMAT field value as canonical VCF text.
///
/// Mirrors noodles' own sample-value writer: scalars use their natural string
/// form, arrays are comma-joined with missing elements rendered as `.`, and
/// genotypes are rendered VCF 4.0/4.3-style (see [`format_genotype`]).
fn format_sample_value(value: &nvcf::variant::record_buf::samples::sample::Value) -> String {
    use nvcf::variant::record_buf::samples::sample::value::Array;
    use nvcf::variant::record_buf::samples::sample::Value;

    /// Render an array of optional scalars, joining with `,` and rendering
    /// missing elements as `.`.
    fn join_array<T: std::fmt::Display>(values: &[Option<T>]) -> String {
        values
            .iter()
            .map(|v| {
                v.as_ref()
                    .map_or_else(|| ".".to_string(), |x| x.to_string())
            })
            .collect::<Vec<_>>()
            .join(",")
    }

    match value {
        Value::Integer(v) => v.to_string(),
        Value::Float(v) => v.to_string(),
        Value::Character(v) => v.to_string(),
        Value::String(v) => v.clone(),
        Value::Genotype(genotype) => format_genotype(genotype),
        Value::Array(array) => match array {
            Array::Integer(values) => join_array(values),
            Array::Float(values) => join_array(values),
            Array::Character(values) => join_array(values),
            Array::String(values) => join_array(values),
        },
    }
}

/// Render a genotype value as VCF 4.0/4.3-style text (e.g. `0/1`, `0|1`,
/// `./.`, `0|1|2`).
///
/// The first allele contributes only its position (a missing position is `.`);
/// each subsequent allele contributes its phasing separator (`|` phased, `/`
/// unphased) followed by its position. This matches noodles' `write_genotype`
/// for file formats below VCF 4.4 (the project emits VCF 4.3).
fn format_genotype(
    genotype: &nvcf::variant::record_buf::samples::sample::value::Genotype,
) -> String {
    use nvcf::variant::record::samples::series::value::genotype::Phasing;

    let mut rendered = String::new();
    for (i, allele) in genotype.as_ref().iter().enumerate() {
        if i > 0 {
            rendered.push(match allele.phasing() {
                Phasing::Phased => '|',
                Phasing::Unphased => '/',
            });
        }
        match allele.position() {
            Some(position) => rendered.push_str(&position.to_string()),
            None => rendered.push('.'),
        }
    }
    rendered
}

/// Split a multi-allelic VCF record into multiple bi-allelic records
pub fn split_multiallelic(record: &VcfRecord) -> Vec<VcfRecord> {
    if !record.is_multiallelic() {
        return vec![record.clone()];
    }

    record
        .alternate
        .iter()
        .map(|alt| {
            let mut new_record = record.clone();
            new_record.alternate = vec![alt.clone()];
            new_record
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINIMAL_VCF: &str = r#"##fileformat=VCFv4.3
##contig=<ID=chr1,length=249250621>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	12345	rs123	A	G	30	PASS	DP=100
chr1	12346	.	AT	A	20	.	.
chr1	12347	.	A	G,T	40	PASS	DP=50
"#;

    #[test]
    fn test_parse_vcf_string() {
        let reader = parse_vcf_string(MINIMAL_VCF).unwrap();
        assert_eq!(reader.header().contigs, vec!["chr1"]);
        assert!(reader.header().samples.is_empty());
    }

    #[test]
    fn test_read_records() {
        let mut reader = parse_vcf_string(MINIMAL_VCF).unwrap();

        // First record: SNV
        let record1 = reader.read_record().unwrap().unwrap();
        assert_eq!(record1.chrom, "chr1");
        assert_eq!(record1.pos, 12345);
        assert_eq!(record1.id, Some("rs123".to_string()));
        assert_eq!(record1.reference, "A");
        assert_eq!(record1.alternate, vec!["G"]);
        assert!(record1.is_snv());
        assert!(record1.passes_filters());

        // Second record: Deletion
        let record2 = reader.read_record().unwrap().unwrap();
        assert_eq!(record2.pos, 12346);
        assert_eq!(record2.reference, "AT");
        assert_eq!(record2.alternate, vec!["A"]);
        assert!(record2.is_deletion());
        assert!(record2.id.is_none());

        // Third record: Multi-allelic
        let record3 = reader.read_record().unwrap().unwrap();
        assert!(record3.is_multiallelic());
        assert_eq!(record3.alternate, vec!["G", "T"]);

        // EOF
        assert!(reader.read_record().unwrap().is_none());
    }

    #[test]
    fn test_records_iterator() {
        let reader = parse_vcf_string(MINIMAL_VCF).unwrap();
        let records: Vec<_> = reader.records().collect();
        assert_eq!(records.len(), 3);
        assert!(records.iter().all(|r| r.is_ok()));
    }

    #[test]
    fn test_split_multiallelic() {
        let record = VcfRecord::new(
            "chr1".to_string(),
            100,
            "A".to_string(),
            vec!["G".to_string(), "T".to_string()],
        );

        let split = split_multiallelic(&record);
        assert_eq!(split.len(), 2);
        assert_eq!(split[0].alternate, vec!["G"]);
        assert_eq!(split[1].alternate, vec!["T"]);
    }

    #[test]
    fn test_split_biallelic() {
        let record = VcfRecord::snv("chr1", 100, 'A', 'G');
        let split = split_multiallelic(&record);
        assert_eq!(split.len(), 1);
    }

    #[test]
    fn test_genome_build_propagation() {
        let reader = parse_vcf_string(MINIMAL_VCF)
            .unwrap()
            .with_genome_build(GenomeBuild::GRCh37);

        let records: Vec<_> = reader.records().filter_map(|r| r.ok()).collect();
        assert!(records
            .iter()
            .all(|r| r.genome_build == GenomeBuild::GRCh37));
    }

    #[test]
    fn test_info_parsing() {
        let reader = parse_vcf_string(MINIMAL_VCF).unwrap();
        let records: Vec<_> = reader.records().filter_map(|r| r.ok()).collect();

        // First record should have DP=100
        assert_eq!(records[0].info.get("DP"), Some(&InfoValue::Integer(100)));

        // Second record should have empty INFO
        assert!(records[1].info.is_empty() || !records[1].info.contains_key("DP"));
    }

    /// VCF with FORMAT and per-sample genotype columns for sample-value tests.
    const SAMPLE_VCF: &str = r#"##fileformat=VCFv4.3
##contig=<ID=chr1,length=249250621>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	30	PASS	.	GT:DP:GQ	0/1:35:99	1|1:40:50
chr1	200	.	C	T	30	PASS	.	GT:DP:GQ	./.:.:.	0|1:12:88
chr1	300	.	G	A	30	PASS	.	GT:DP:GQ	1/1	0/0:20:40
"#;

    #[test]
    fn test_parse_single_sample_format_values() {
        // A single record's first sample carries each FORMAT field, keyed by name.
        let mut reader = parse_vcf_string(SAMPLE_VCF).unwrap();
        let record = reader.read_record().unwrap().unwrap();

        assert_eq!(record.format.as_deref(), Some("GT:DP:GQ"));
        assert_eq!(record.samples.len(), 2);

        let s1 = &record.samples[0];
        assert_eq!(s1.get("GT").map(String::as_str), Some("0/1"));
        assert_eq!(s1.get("DP").map(String::as_str), Some("35"));
        assert_eq!(s1.get("GQ").map(String::as_str), Some("99"));
    }

    #[test]
    fn test_parse_multi_sample_values() {
        // Per-sample maps are populated in sample-column order (S1 then S2),
        // and phased genotypes (`|`) are preserved verbatim.
        let mut reader = parse_vcf_string(SAMPLE_VCF).unwrap();
        let record = reader.read_record().unwrap().unwrap();

        let s2 = &record.samples[1];
        assert_eq!(s2.get("GT").map(String::as_str), Some("1|1"));
        assert_eq!(s2.get("DP").map(String::as_str), Some("40"));
        assert_eq!(s2.get("GQ").map(String::as_str), Some("50"));

        // The two columns must not be conflated.
        assert_ne!(record.samples[0].get("GT"), record.samples[1].get("GT"));
    }

    #[test]
    fn test_parse_missing_sample_values() {
        // Missing field values (`.`) are omitted from the map; a missing-allele
        // genotype (`./.`) is a present GT value and is rendered as such.
        let mut reader = parse_vcf_string(SAMPLE_VCF).unwrap();
        reader.read_record().unwrap().unwrap(); // skip record 1
        let record = reader.read_record().unwrap().unwrap();

        let s1 = &record.samples[0];
        assert_eq!(s1.get("GT").map(String::as_str), Some("./."));
        assert!(!s1.contains_key("DP"), "missing DP should be omitted");
        assert!(!s1.contains_key("GQ"), "missing GQ should be omitted");

        let s2 = &record.samples[1];
        assert_eq!(s2.get("GT").map(String::as_str), Some("0|1"));
        assert_eq!(s2.get("DP").map(String::as_str), Some("12"));
    }

    #[test]
    fn test_parse_ragged_sample_columns() {
        // A sample column with fewer fields than FORMAT declares: trailing keys
        // are simply absent from the map (and Display renders them as `.`).
        let mut reader = parse_vcf_string(SAMPLE_VCF).unwrap();
        for _ in 0..2 {
            reader.read_record().unwrap().unwrap();
        }
        let record = reader.read_record().unwrap().unwrap();

        // S1 has only "1/1" (GT), no DP/GQ columns.
        let s1 = &record.samples[0];
        assert_eq!(s1.get("GT").map(String::as_str), Some("1/1"));
        assert!(!s1.contains_key("DP"));
        assert!(!s1.contains_key("GQ"));

        // S2 is fully populated.
        let s2 = &record.samples[1];
        assert_eq!(s2.get("GT").map(String::as_str), Some("0/0"));
        assert_eq!(s2.get("DP").map(String::as_str), Some("20"));
    }

    #[test]
    fn test_sites_only_still_parses() {
        // A VCF with no FORMAT/sample columns parses with no format and no
        // sample maps.
        let mut reader = parse_vcf_string(MINIMAL_VCF).unwrap();
        let record = reader.read_record().unwrap().unwrap();
        assert!(record.format.is_none());
        assert!(record.samples.is_empty());
    }

    #[test]
    fn test_sample_roundtrip_display() {
        // FORMAT and sample columns survive a parse -> Display round-trip. Uses
        // only integer/GT fields to avoid float text-formatting ambiguity.
        let mut reader = parse_vcf_string(SAMPLE_VCF).unwrap();
        let record = reader.read_record().unwrap().unwrap();

        let rendered = format!("{}", record);
        let columns: Vec<&str> = rendered.split('\t').collect();
        // CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1 S2
        assert_eq!(columns[8], "GT:DP:GQ");
        assert_eq!(columns[9], "0/1:35:99");
        assert_eq!(columns[10], "1|1:40:50");
    }
}
