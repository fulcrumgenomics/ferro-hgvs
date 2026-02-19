//! Batch processing for VCF files
//!
//! This module provides efficient batch processing of VCF records,
//! supporting streaming, parallel processing, and memory-efficient operations.

use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use crate::error::FerroError;
use crate::reference::loader::TranscriptDb;

use super::annotate::VcfAnnotator;
use super::record::VcfRecord;

/// Statistics from batch processing
#[derive(Debug, Clone, Default)]
pub struct BatchStats {
    /// Total records processed
    pub total: usize,
    /// Successfully annotated records
    pub annotated: usize,
    /// Records that failed annotation
    pub failed: usize,
    /// Records with no overlapping transcripts
    pub no_transcripts: usize,
    /// Error messages from failed records (when skip_failed is true)
    /// Contains up to 100 error messages to prevent memory issues with large files
    pub errors: Vec<String>,
}

/// Maximum number of error messages to store to prevent memory issues
const MAX_STORED_ERRORS: usize = 100;

impl BatchStats {
    /// Get the success rate
    pub fn success_rate(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            self.annotated as f64 / self.total as f64
        }
    }

    /// Record an error, storing up to MAX_STORED_ERRORS messages
    pub fn record_error(&mut self, error: &crate::error::FerroError) {
        self.failed += 1;
        if self.errors.len() < MAX_STORED_ERRORS {
            self.errors.push(error.to_string());
        }
    }
}

/// Progress callback type
pub type ProgressCallback = Box<dyn Fn(usize, usize) + Send + Sync>;

/// Configuration for batch processing
#[derive(Default)]
pub struct BatchConfig {
    /// Number of records to process in each batch
    pub batch_size: usize,
    /// Whether to skip records that fail annotation
    pub skip_failed: bool,
    /// Progress callback (called with current count and total estimate)
    pub progress: Option<ProgressCallback>,
}

impl BatchConfig {
    /// Create a new batch configuration
    pub fn new() -> Self {
        Self {
            batch_size: 1000,
            skip_failed: true,
            progress: None,
        }
    }

    /// Set the batch size
    pub fn with_batch_size(mut self, size: usize) -> Self {
        self.batch_size = size;
        self
    }

    /// Set whether to skip failed records
    pub fn skip_failed(mut self, skip: bool) -> Self {
        self.skip_failed = skip;
        self
    }

    /// Set the progress callback
    pub fn with_progress<F>(mut self, callback: F) -> Self
    where
        F: Fn(usize, usize) + Send + Sync + 'static,
    {
        self.progress = Some(Box::new(callback));
        self
    }
}

/// Batch processor for VCF records
pub struct BatchProcessor<'a> {
    /// Transcript database
    db: &'a TranscriptDb,
    /// Configuration
    config: BatchConfig,
    /// Statistics
    stats: BatchStats,
}

impl<'a> BatchProcessor<'a> {
    /// Create a new batch processor
    pub fn new(db: &'a TranscriptDb) -> Self {
        Self {
            db,
            config: BatchConfig::new(),
            stats: BatchStats::default(),
        }
    }

    /// Set configuration
    pub fn with_config(mut self, config: BatchConfig) -> Self {
        self.config = config;
        self
    }

    /// Get current statistics
    pub fn stats(&self) -> &BatchStats {
        &self.stats
    }

    /// Process a single record
    pub fn process_record(&self, record: &VcfRecord) -> Result<VcfRecord, FerroError> {
        let annotator = VcfAnnotator::new(self.db);
        annotator.annotate(record)
    }

    /// Process a batch of records
    pub fn process_batch(&mut self, records: &[VcfRecord]) -> Vec<Result<VcfRecord, FerroError>> {
        let annotator = VcfAnnotator::new(self.db);

        let results: Vec<_> = records.iter().map(|r| annotator.annotate(r)).collect();

        // Update stats
        for result in &results {
            self.stats.total += 1;
            match result {
                Ok(_) => self.stats.annotated += 1,
                Err(_) => self.stats.failed += 1,
            }
        }

        results
    }

    /// Stream process from reader to writer
    pub fn process_stream<R: Read, W: Write>(
        &mut self,
        reader: R,
        writer: W,
    ) -> Result<BatchStats, FerroError> {
        let buf_reader = BufReader::new(reader);
        let mut buf_writer = BufWriter::new(writer);

        let mut batch = Vec::with_capacity(self.config.batch_size);

        for line in buf_reader.lines() {
            let line = line?;

            if line.starts_with('#') {
                // Pass through header lines
                writeln!(buf_writer, "{}", line)?;
                continue;
            }

            // Empty lines are valid in VCF files and are silently skipped.
            // This matches the behavior of most VCF parsers (bcftools, htslib).
            if line.is_empty() {
                continue;
            }

            // Parse record
            let record = match parse_vcf_line(&line) {
                Ok(r) => r,
                Err(e) => {
                    if !self.config.skip_failed {
                        return Err(e);
                    }
                    // Store error for reporting even when skipping
                    self.stats.record_error(&e);
                    continue;
                }
            };

            batch.push(record);

            // Process batch when full
            if batch.len() >= self.config.batch_size {
                self.write_batch(&batch, &mut buf_writer)?;

                // Report progress
                if let Some(ref callback) = self.config.progress {
                    callback(self.stats.total, 0); // 0 = unknown total
                }

                batch.clear();
            }
        }

        // Process remaining records
        if !batch.is_empty() {
            self.write_batch(&batch, &mut buf_writer)?;
        }

        buf_writer.flush()?;
        Ok(self.stats.clone())
    }

    /// Write a batch of records to the writer
    fn write_batch<W: Write>(
        &mut self,
        records: &[VcfRecord],
        writer: &mut W,
    ) -> Result<(), FerroError> {
        let results = self.process_batch(records);

        for (result, original) in results.iter().zip(records.iter()) {
            match result {
                Ok(annotated) => {
                    writeln!(writer, "{}", annotated)?;
                }
                Err(_) => {
                    if !self.config.skip_failed {
                        // Write original if not skipping failed
                        writeln!(writer, "{}", original)?;
                    }
                }
            }
        }

        Ok(())
    }

    /// Process a file
    pub fn process_file<P: AsRef<Path>>(
        &mut self,
        input: P,
        output: P,
    ) -> Result<BatchStats, FerroError> {
        let reader = std::fs::File::open(input.as_ref())?;
        let writer = std::fs::File::create(output.as_ref())?;
        self.process_stream(reader, writer)
    }
}

/// Parallel batch processor (requires rayon feature)
#[cfg(feature = "parallel")]
#[allow(dead_code)]
pub struct ParallelBatchProcessor<'a> {
    /// Transcript database
    db: &'a TranscriptDb,
    /// Configuration
    config: BatchConfig,
}

#[cfg(feature = "parallel")]
#[allow(dead_code)]
impl<'a> ParallelBatchProcessor<'a> {
    /// Create a new parallel batch processor
    pub fn new(db: &'a TranscriptDb) -> Self {
        Self {
            db,
            config: BatchConfig::new(),
        }
    }

    /// Process records in parallel
    pub fn process_parallel(&self, records: &[VcfRecord]) -> Vec<Result<VcfRecord, FerroError>> {
        use rayon::prelude::*;

        let annotator = VcfAnnotator::new(self.db);

        records.par_iter().map(|r| annotator.annotate(r)).collect()
    }
}

/// Simple iterator-based streaming processor
#[allow(dead_code)]
pub struct StreamingProcessor<'a, R: BufRead> {
    /// Reader
    reader: R,
    /// Transcript database
    db: &'a TranscriptDb,
    /// Header lines
    header: Vec<String>,
    /// Whether header has been consumed
    header_consumed: bool,
}

#[allow(dead_code)]
impl<'a, R: BufRead> StreamingProcessor<'a, R> {
    /// Create a new streaming processor
    pub fn new(reader: R, db: &'a TranscriptDb) -> Self {
        Self {
            reader,
            db,
            header: Vec::new(),
            header_consumed: false,
        }
    }

    /// Get the VCF header lines
    pub fn header(&mut self) -> Result<&[String], FerroError> {
        if !self.header_consumed {
            self.consume_header()?;
        }
        Ok(&self.header)
    }

    /// Consume header lines from the reader
    fn consume_header(&mut self) -> Result<(), FerroError> {
        let mut line = String::new();
        while self.reader.read_line(&mut line)? > 0 {
            if line.starts_with('#') {
                self.header.push(line.trim_end().to_string());
                line.clear();
            } else {
                // Put back the first data line by keeping it in a buffer
                // Note: This is a simplified implementation
                break;
            }
        }
        self.header_consumed = true;
        Ok(())
    }
}

/// Parse a single VCF line into a record
fn parse_vcf_line(line: &str) -> Result<VcfRecord, FerroError> {
    use crate::reference::transcript::GenomeBuild;

    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 5 {
        return Err(FerroError::Parse {
            msg: "Invalid VCF line: not enough fields".to_string(),
            pos: 0,
            diagnostic: None,
        });
    }

    let chrom = fields[0].to_string();
    let pos: u64 = fields[1].parse().map_err(|_| FerroError::Parse {
        msg: "Invalid position".to_string(),
        pos: 0,
        diagnostic: None,
    })?;
    let id = if fields[2] == "." {
        None
    } else {
        Some(fields[2].to_string())
    };
    let reference = fields[3].to_string();
    let alternate: Vec<String> = fields[4].split(',').map(|s| s.to_string()).collect();

    Ok(VcfRecord {
        chrom,
        pos,
        id,
        reference,
        alternate,
        quality: None,
        filter: None,
        info: Default::default(),
        format: None,
        samples: Vec::new(),
        genome_build: GenomeBuild::GRCh38,
    })
}

/// Count-based progress reporter
#[allow(dead_code)]
pub fn count_progress(count: &Arc<AtomicUsize>) -> impl Fn(usize, usize) + Send + Sync + Clone {
    let count = Arc::clone(count);
    move |current, _total| {
        count.store(current, Ordering::SeqCst);
    }
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
    fn test_batch_stats() {
        let mut stats = BatchStats::default();
        assert_eq!(stats.success_rate(), 0.0);

        stats.total = 100;
        stats.annotated = 90;
        stats.failed = 10;
        assert!((stats.success_rate() - 0.9).abs() < 0.001);
    }

    #[test]
    fn test_batch_config() {
        let config = BatchConfig::new().with_batch_size(500).skip_failed(false);

        assert_eq!(config.batch_size, 500);
        assert!(!config.skip_failed);
    }

    #[test]
    fn test_batch_processor() {
        let db = create_test_db();
        let processor = BatchProcessor::new(&db);

        let record = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let result = processor.process_record(&record);
        assert!(result.is_ok());
    }

    #[test]
    fn test_parse_vcf_line() {
        let line = "chr1\t12345\t.\tA\tG\t.\t.\t.";
        let record = parse_vcf_line(line).unwrap();
        assert_eq!(record.chrom, "chr1");
        assert_eq!(record.pos, 12345);
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["G"]);
    }

    #[test]
    fn test_parse_vcf_line_multiallelic() {
        let line = "chr1\t12345\trs123\tA\tG,C\t.\t.\t.";
        let record = parse_vcf_line(line).unwrap();
        assert_eq!(record.id, Some("rs123".to_string()));
        assert_eq!(record.alternate, vec!["G", "C"]);
    }

    // ===== BatchStats Extended Tests =====

    #[test]
    fn test_batch_stats_success_rate_zero_total() {
        let stats = BatchStats::default();
        assert_eq!(stats.success_rate(), 0.0);
    }

    #[test]
    fn test_batch_stats_success_rate_all_successful() {
        let stats = BatchStats {
            total: 100,
            annotated: 100,
            failed: 0,
            no_transcripts: 0,
            errors: Vec::new(),
        };
        assert!((stats.success_rate() - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_batch_stats_success_rate_partial() {
        let stats = BatchStats {
            total: 100,
            annotated: 75,
            failed: 25,
            no_transcripts: 10,
            errors: Vec::new(),
        };
        assert!((stats.success_rate() - 0.75).abs() < 0.001);
    }

    // ===== BatchConfig Extended Tests =====

    #[test]
    fn test_batch_config_default() {
        let config = BatchConfig::default();
        assert_eq!(config.batch_size, 0);
        assert!(!config.skip_failed);
        assert!(config.progress.is_none());
    }

    #[test]
    fn test_batch_config_new() {
        let config = BatchConfig::new();
        assert_eq!(config.batch_size, 1000);
        assert!(config.skip_failed);
        assert!(config.progress.is_none());
    }

    #[test]
    fn test_batch_config_builder_chain() {
        let config = BatchConfig::new().with_batch_size(500).skip_failed(false);

        assert_eq!(config.batch_size, 500);
        assert!(!config.skip_failed);
    }

    #[test]
    fn test_batch_config_with_progress() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        use std::sync::Arc;

        let counter = Arc::new(AtomicUsize::new(0));
        let counter_clone = counter.clone();

        let config = BatchConfig::new().with_progress(move |current, _total| {
            counter_clone.store(current, Ordering::SeqCst);
        });

        assert!(config.progress.is_some());

        // Call the callback
        if let Some(ref callback) = config.progress {
            callback(42, 100);
        }

        assert_eq!(counter.load(Ordering::SeqCst), 42);
    }

    // ===== BatchProcessor Extended Tests =====

    #[test]
    fn test_batch_processor_with_config() {
        let db = create_test_db();
        let config = BatchConfig::new().with_batch_size(100);
        let processor = BatchProcessor::new(&db).with_config(config);

        let stats = processor.stats();
        assert_eq!(stats.total, 0);
    }

    #[test]
    fn test_batch_processor_stats() {
        let db = create_test_db();
        let processor = BatchProcessor::new(&db);

        let stats = processor.stats();
        assert_eq!(stats.total, 0);
        assert_eq!(stats.annotated, 0);
        assert_eq!(stats.failed, 0);
    }

    #[test]
    fn test_batch_processor_process_batch() {
        let db = create_test_db();
        let mut processor = BatchProcessor::new(&db);

        let records = vec![
            VcfRecord::snv("chr1", 1050, 'A', 'G'),
            VcfRecord::snv("chr1", 1060, 'C', 'T'),
        ];

        let results = processor.process_batch(&records);
        assert_eq!(results.len(), 2);

        let stats = processor.stats();
        assert_eq!(stats.total, 2);
    }

    #[test]
    fn test_batch_processor_process_record_outside_transcript() {
        let db = create_test_db();
        let processor = BatchProcessor::new(&db);

        // Position outside any transcript
        let record = VcfRecord::snv("chr1", 5000, 'A', 'G');
        let result = processor.process_record(&record);

        // Should still succeed but without annotation
        assert!(result.is_ok());
    }

    // ===== parse_vcf_line Extended Tests =====

    #[test]
    fn test_parse_vcf_line_too_few_fields() {
        let line = "chr1\t12345\t.\tA";
        let result = parse_vcf_line(line);
        assert!(result.is_err());

        if let Err(FerroError::Parse { msg, .. }) = result {
            assert!(msg.contains("not enough fields"));
        }
    }

    #[test]
    fn test_parse_vcf_line_invalid_position() {
        let line = "chr1\tnotanumber\t.\tA\tG";
        let result = parse_vcf_line(line);
        assert!(result.is_err());

        if let Err(FerroError::Parse { msg, .. }) = result {
            assert!(msg.contains("Invalid position"));
        }
    }

    #[test]
    fn test_parse_vcf_line_with_id() {
        let line = "chr1\t12345\trs12345\tA\tG";
        let record = parse_vcf_line(line).unwrap();
        assert_eq!(record.id, Some("rs12345".to_string()));
    }

    #[test]
    fn test_parse_vcf_line_deletion() {
        let line = "chr1\t12345\t.\tATG\tA";
        let record = parse_vcf_line(line).unwrap();
        assert_eq!(record.reference, "ATG");
        assert_eq!(record.alternate, vec!["A"]);
    }

    #[test]
    fn test_parse_vcf_line_insertion() {
        let line = "chr1\t12345\t.\tA\tATGC";
        let record = parse_vcf_line(line).unwrap();
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["ATGC"]);
    }

    #[test]
    fn test_parse_vcf_line_symbolic_alt() {
        let line = "chr1\t12345\t.\tA\t<DEL>";
        let record = parse_vcf_line(line).unwrap();
        assert_eq!(record.alternate, vec!["<DEL>"]);
    }

    // ===== Stream Processing Tests =====

    #[test]
    fn test_batch_processor_process_stream() {
        let db = create_test_db();
        let mut processor = BatchProcessor::new(&db)
            .with_config(BatchConfig::new().with_batch_size(2).skip_failed(true));

        let input = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1050\t.\tA\tG\t.\t.\t.\nchr1\t1060\t.\tC\tT\t.\t.\t.\n";
        let output = Vec::new();

        let result = processor.process_stream(&input[..], output);
        assert!(result.is_ok());

        let stats = result.unwrap();
        assert!(stats.total > 0 || stats.failed > 0);
    }

    #[test]
    fn test_batch_processor_process_stream_empty() {
        let db = create_test_db();
        let mut processor = BatchProcessor::new(&db);

        let input = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n";
        let output = Vec::new();

        let result = processor.process_stream(&input[..], output);
        assert!(result.is_ok());
    }

    #[test]
    fn test_batch_processor_process_stream_with_invalid_lines() {
        let db = create_test_db();
        let mut processor =
            BatchProcessor::new(&db).with_config(BatchConfig::new().skip_failed(true));

        // Include an invalid line
        let input = b"##header\nchr1\tinvalidpos\t.\tA\tG\t.\t.\t.\nchr1\t1050\t.\tA\tG\t.\t.\t.\n";
        let output = Vec::new();

        let result = processor.process_stream(&input[..], output);
        assert!(result.is_ok());

        let stats = result.unwrap();
        assert!(stats.failed > 0);
    }

    #[test]
    fn test_batch_processor_process_stream_skip_empty_lines() {
        let db = create_test_db();
        let mut processor = BatchProcessor::new(&db);

        let input = b"##header\n\n\nchr1\t1050\t.\tA\tG\t.\t.\t.\n\n";
        let output = Vec::new();

        let result = processor.process_stream(&input[..], output);
        assert!(result.is_ok());
    }

    // ===== count_progress Tests =====

    #[test]
    fn test_count_progress() {
        let count = Arc::new(AtomicUsize::new(0));
        let progress = count_progress(&count);

        progress(100, 1000);
        assert_eq!(count.load(Ordering::SeqCst), 100);

        progress(200, 1000);
        assert_eq!(count.load(Ordering::SeqCst), 200);
    }

    // ===== StreamingProcessor Tests =====

    #[test]
    fn test_streaming_processor_new() {
        let db = create_test_db();
        let input = b"##header\n#CHROM\tPOS\tID\tREF\tALT\nchr1\t1050\t.\tA\tG";
        let reader = BufReader::new(&input[..]);

        let processor = StreamingProcessor::new(reader, &db);
        assert!(!processor.header_consumed);
        assert!(processor.header.is_empty());
    }

    #[test]
    fn test_streaming_processor_header() {
        let db = create_test_db();
        let input = b"##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1>\n#CHROM\tPOS\tID\tREF\tALT\nchr1\t1050\t.\tA\tG";
        let reader = BufReader::new(&input[..]);

        let mut processor = StreamingProcessor::new(reader, &db);
        let header = processor.header().unwrap();

        assert_eq!(header.len(), 3);
        assert!(header[0].starts_with("##fileformat"));
        assert!(header[2].starts_with("#CHROM"));
    }
}
