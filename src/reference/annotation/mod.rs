//! Annotation file loading (GFF3 / GTF).
//!
//! The primary entry point is [`load_annotations`], which reads a GFF3 or GTF
//! file and produces a `TranscriptDb` plus a [`LoaderReport`] aggregating
//! per-record diagnostics. Format is auto-detected from the path extension
//! and content unless overridden via [`LoaderConfig::with_format`].
//!
//! # Example
//!
//! ```no_run
//! use ferro_hgvs::reference::annotation::{load_annotations, LoaderConfig};
//!
//! let cfg = LoaderConfig::new();
//! let (db, report) = load_annotations("annotations.gff3", &cfg)?;
//! println!("{}", report.summary_line());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! See `docs/superpowers/specs/2026-05-13-gff-gtf-loader-rewrite-design.md`
//! for the full design.

pub(crate) mod builder;
pub mod diagnostics;
pub(crate) mod feature;
pub mod format_detect;
pub(crate) mod graph;
pub(crate) mod record;

pub use diagnostics::{
    DiagnosticPayload, LoaderDiagnostic, LoaderReport, Severity, SourceLocation,
};
pub use format_detect::AnnotationFormat;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use crate::error::FerroError;
use crate::reference::loader::TranscriptDb;
use crate::reference::transcript::GenomeBuild;

use builder::build_transcripts;
use graph::FeatureGraph;
use record::{AnnotationRecord, Gff3Record, GtfRecord};

/// Controls how malformed records are handled during loading.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorMode {
    /// Any `Severity::Error` diagnostic causes [`load_annotations`] to return `Err`.
    Strict,
    /// Bad records are dropped and their diagnostics recorded; loading continues (default).
    Lenient,
    /// Bad records are silently dropped; no diagnostics are recorded.
    Silent,
}

/// Configuration for [`load_annotations`].
///
/// Construct with [`LoaderConfig::new`] and chain builder methods as needed.
#[non_exhaustive]
pub struct LoaderConfig {
    pub format: Option<AnnotationFormat>,
    pub error_mode: ErrorMode,
    pub genome_build: Option<GenomeBuild>,
}

impl Default for LoaderConfig {
    fn default() -> Self {
        Self::new()
    }
}

impl LoaderConfig {
    /// Create a new `LoaderConfig` with default settings:
    /// format auto-detected, [`ErrorMode::Lenient`], genome build GRCh38.
    pub fn new() -> Self {
        Self {
            format: None,
            error_mode: ErrorMode::Lenient,
            genome_build: None,
        }
    }

    /// Enable strict error mode: any parse error causes [`load_annotations`] to return `Err`.
    pub fn with_strict(mut self) -> Self {
        self.error_mode = ErrorMode::Strict;
        self
    }

    /// Enable lenient error mode (default): bad records are dropped, diagnostics recorded.
    pub fn with_lenient(mut self) -> Self {
        self.error_mode = ErrorMode::Lenient;
        self
    }

    /// Enable silent error mode: bad records are dropped without recording diagnostics.
    pub fn with_silent(mut self) -> Self {
        self.error_mode = ErrorMode::Silent;
        self
    }

    /// Override auto-detected annotation format with an explicit value.
    pub fn with_format(mut self, fmt: AnnotationFormat) -> Self {
        self.format = Some(fmt);
        self
    }

    /// Set the genome build / assembly version used when constructing transcripts.
    pub fn with_genome_build(mut self, b: GenomeBuild) -> Self {
        self.genome_build = Some(b);
        self
    }
}

/// Parse one line through the format-specific `AnnotationRecord` impl and
/// either ingest the resulting record into `graph` or record an
/// `E-LOAD-001 MalformedRecord` diagnostic. Hoisting this out of the
/// per-format match arms keeps the diagnostic-construction logic in one
/// place — the two arms were byte-for-byte identical except for the
/// parser type.
fn parse_and_ingest<R: AnnotationRecord>(
    line: &str,
    line_no: u64,
    source_path: &Path,
    graph: &mut FeatureGraph,
    report: &mut LoaderReport,
) {
    match R::parse(line, line_no) {
        Ok(Some(rec)) => graph.ingest(&rec),
        Ok(None) => {}
        Err(e) => report.record(LoaderDiagnostic::error(
            "E-LOAD-001",
            format!("MalformedRecord: {}", e),
            SourceLocation {
                path: source_path.to_path_buf(),
                line: line_no,
            },
            None,
            DiagnosticPayload::MalformedRecord {
                field: "row",
                value: line.to_string(),
            },
        )),
    }
}

/// Load transcripts from a GFF3 or GTF file.
///
/// Returns a populated `TranscriptDb` and a [`LoaderReport`] containing
/// any diagnostics produced during parsing. In [`ErrorMode::Strict`], the
/// presence of any `Severity::Error` diagnostic causes this function to
/// return `Err`; in [`ErrorMode::Lenient`] (the default) bad records are
/// dropped and their diagnostics recorded in the report.
///
/// The annotation format is auto-detected from the file extension and content
/// unless overridden via [`LoaderConfig::with_format`].
pub fn load_annotations<P: AsRef<Path>>(
    path: P,
    config: &LoaderConfig,
) -> Result<(TranscriptDb, LoaderReport), FerroError> {
    let path = path.as_ref();
    let format = match config.format {
        Some(f) => f,
        None => format_detect::detect_format(path)?,
    };
    let genome_build = config.genome_build.unwrap_or(GenomeBuild::GRCh38);

    let file = File::open(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", path.display(), e),
    })?;
    let reader = BufReader::new(file);

    let mut graph = FeatureGraph::new();
    // Silent mode is enforced at the report layer so every recording site
    // (`Gff3Record::parse`, `GtfRecord::parse`, orphan resolution, the
    // transcript builder) observes the same gate without each call site
    // needing its own conditional.
    let mut report = LoaderReport {
        silent: config.error_mode == ErrorMode::Silent,
        ..Default::default()
    };
    let source_path: PathBuf = path.into();

    // Read raw bytes per line and decode as UTF-8 so a single bad byte
    // becomes an `E-LOAD-001` diagnostic rather than aborting the whole
    // load. Lossy decoding (`String::from_utf8_lossy`) preserves the rest
    // of the line for the parser; the diagnostic still fires on the
    // original undecoded bytes so strict mode rejects the file.
    let mut reader = reader;
    let mut buf: Vec<u8> = Vec::new();
    let mut line_no: u64 = 0;
    loop {
        buf.clear();
        let n = reader
            .read_until(b'\n', &mut buf)
            .map_err(|e| FerroError::Io {
                msg: format!("read: {}", e),
            })?;
        if n == 0 {
            break;
        }
        line_no += 1;
        // Strip trailing \n and \r from the chunk before decoding so the
        // diagnostic's `value` matches what a UTF-8-clean line would carry.
        let trimmed = match buf.last() {
            Some(b'\n') => {
                let len = buf.len() - 1;
                let len = if len > 0 && buf[len - 1] == b'\r' {
                    len - 1
                } else {
                    len
                };
                &buf[..len]
            }
            _ => &buf[..],
        };
        let (line, malformed_utf8) = match std::str::from_utf8(trimmed) {
            Ok(s) => (s.to_string(), false),
            Err(_) => (String::from_utf8_lossy(trimmed).into_owned(), true),
        };
        if malformed_utf8 {
            report.record(LoaderDiagnostic::error(
                "E-LOAD-001",
                format!("MalformedRecord: invalid UTF-8 on line {}", line_no),
                SourceLocation {
                    path: source_path.clone(),
                    line: line_no,
                },
                None,
                DiagnosticPayload::MalformedRecord {
                    field: "row",
                    value: line,
                },
            ));
            continue;
        }
        match format {
            AnnotationFormat::Gff3 => {
                parse_and_ingest::<Gff3Record>(
                    &line,
                    line_no,
                    &source_path,
                    &mut graph,
                    &mut report,
                );
            }
            AnnotationFormat::Gtf => {
                parse_and_ingest::<GtfRecord>(
                    &line,
                    line_no,
                    &source_path,
                    &mut graph,
                    &mut report,
                );
            }
        }
    }

    graph.resolve();
    for (line_no, parent_id, child_id) in &graph.orphan_diagnostics {
        report.record(LoaderDiagnostic::warning(
            "W-LOAD-010",
            format!("OrphanFeature: parent '{}' not found", parent_id),
            SourceLocation {
                path: source_path.clone(),
                line: *line_no,
            },
            child_id.clone(),
            DiagnosticPayload::OrphanFeature {
                parent_id: parent_id.clone(),
                child_id: child_id.clone(),
            },
        ));
    }
    report.records_dropped += graph.orphan_count;

    let transcripts = build_transcripts(&graph, format, genome_build, source_path, &mut report);

    let mut db = TranscriptDb::with_build(genome_build);
    for tx in transcripts {
        db.add(tx);
    }

    if config.error_mode == ErrorMode::Strict && report.has_error {
        return Err(FerroError::Io {
            msg: format!("Strict loader rejected {} due to errors", path.display()),
        });
    }

    Ok((db, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_gff3(content: &str) -> NamedTempFile {
        let mut tf = tempfile::Builder::new().suffix(".gff3").tempfile().unwrap();
        tf.write_all(content.as_bytes()).unwrap();
        tf.flush().unwrap();
        tf
    }

    fn write_gtf(content: &str) -> NamedTempFile {
        let mut tf = tempfile::Builder::new().suffix(".gtf").tempfile().unwrap();
        tf.write_all(content.as_bytes()).unwrap();
        tf.flush().unwrap();
        tf
    }

    #[test]
    fn load_gff3_issue_183_single_exon_no_exon_line() {
        let f = write_gff3(
            "##gff-version 3\n\
             seq1\t.\tgene\t100\t1200\t.\t+\t.\tID=gene01;Name=gene01\n\
             seq1\t.\tmRNA\t100\t1200\t.\t+\t.\tID=gene01.1;Parent=gene01\n\
             seq1\t.\tCDS\t100\t1200\t.\t+\t0\tParent=gene01.1\n",
        );
        let cfg = LoaderConfig::new();
        let (db, report) = load_annotations(f.path(), &cfg).unwrap();
        assert_eq!(db.len(), 1, "expected 1 transcript from no-exon GFF3");
        assert_eq!(report.transcripts_loaded, 1);
        let tx = db.get("gene01.1").unwrap();
        assert_eq!(tx.exons.len(), 1);
        assert_eq!(tx.exons[0].genomic_start, Some(100));
    }

    #[test]
    fn load_gtf_multi_exon_transcript_end_to_end() {
        let f = write_gtf(
            "chr1\tHAVANA\tgene\t100\t500\t.\t+\t.\tgene_id \"ENSG1\"; gene_name \"GENE1\";\n\
             chr1\tHAVANA\ttranscript\t100\t500\t.\t+\t.\tgene_id \"ENSG1\"; transcript_id \"ENST1\";\n\
             chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_id \"ENSG1\"; transcript_id \"ENST1\";\n\
             chr1\tHAVANA\texon\t300\t500\t.\t+\t.\tgene_id \"ENSG1\"; transcript_id \"ENST1\";\n\
             chr1\tHAVANA\tCDS\t150\t450\t.\t+\t0\tgene_id \"ENSG1\"; transcript_id \"ENST1\";\n",
        );
        let cfg = LoaderConfig::new();
        let (db, _report) = load_annotations(f.path(), &cfg).unwrap();
        assert_eq!(db.len(), 1);
        let tx = db.get("ENST1").unwrap();
        assert_eq!(tx.exons.len(), 2);
        assert_eq!(tx.exons[0].genomic_start, Some(100));
        assert_eq!(tx.exons[1].genomic_end, Some(500));
        assert_eq!(tx.chromosome.as_deref(), Some("chr1"));
    }

    #[test]
    fn auto_detect_gff3_from_content_when_extension_is_generic() {
        let mut tf = tempfile::Builder::new().suffix(".txt").tempfile().unwrap();
        writeln!(tf, "##gff-version 3").unwrap();
        writeln!(tf, "chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1").unwrap();
        writeln!(tf, "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=g1").unwrap();
        writeln!(tf, "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1").unwrap();
        tf.flush().unwrap();
        let cfg = LoaderConfig::new();
        let (db, _report) = load_annotations(tf.path(), &cfg).unwrap();
        assert_eq!(db.len(), 1);
    }

    #[test]
    fn strict_mode_returns_err_on_malformed_record() {
        let f = write_gff3(
            "##gff-version 3\n\
             chr1\t.\tgene\tnot_a_number\t500\t.\t+\t.\tID=g1\n",
        );
        let cfg = LoaderConfig::new().with_strict();
        let result = load_annotations(f.path(), &cfg);
        assert!(result.is_err(), "strict mode must reject malformed record");
    }

    #[test]
    fn strict_mode_rejects_when_errors_come_after_warning_sample_fills() {
        // `sample_diagnostics` is capped at 100. If 100+ warning-level events fill the
        // buffer first (here: orphans recorded by `graph.resolve()`), a later
        // error-level diagnostic (here: an E-LOAD-103 from `build_transcripts`) is
        // counted in `diagnostics_by_code` but not present in the sample. Strict
        // mode must still reject the file in that case.
        let mut content = String::from("##gff-version 3\n");
        // 150 ID-less features referencing a missing parent → 150 orphan warnings.
        for _ in 0..150 {
            content.push_str("chr1\t.\texon\t100\t500\t.\t+\t.\tParent=missing\n");
        }
        // One transcript with unknown strand → builder records E-LOAD-103 error.
        content.push_str("chr1\t.\tmRNA\t100\t500\t.\t.\t.\tID=tx_unknown\n");
        content.push_str("chr1\t.\texon\t100\t500\t.\t.\t.\tParent=tx_unknown\n");
        let f = write_gff3(&content);
        let cfg = LoaderConfig::new().with_strict();
        let result = load_annotations(f.path(), &cfg);
        assert!(
            result.is_err(),
            "strict mode must reject when error-level diagnostic exists beyond the \
             sample cap (got Ok)"
        );
    }

    #[test]
    fn lenient_mode_records_diagnostic_without_failing() {
        let f = write_gff3(
            "##gff-version 3\n\
             chr1\t.\tgene\tnot_a_number\t500\t.\t+\t.\tID=g1\n\
             chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g2\n\
             chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx2;Parent=g2\n\
             chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx2\n",
        );
        let cfg = LoaderConfig::new(); // Lenient by default
        let (db, report) = load_annotations(f.path(), &cfg).unwrap();
        assert_eq!(db.len(), 1, "good records should still load");
        assert!(
            report.diagnostics_by_code.contains_key("E-LOAD-001"),
            "malformed record should be recorded"
        );
    }

    #[test]
    fn minus_strand_multi_exon_stored_in_transcript_order() {
        let f = write_gff3(
            "##gff-version 3\n\
             chr1\t.\tmRNA\t100\t500\t.\t-\t.\tID=tx1\n\
             chr1\t.\texon\t100\t200\t.\t-\t.\tParent=tx1\n\
             chr1\t.\texon\t300\t500\t.\t-\t.\tParent=tx1\n",
        );
        let cfg = LoaderConfig::new();
        let (db, _report) = load_annotations(f.path(), &cfg).unwrap();
        let tx = db.get("tx1").unwrap();
        assert_eq!(tx.strand, crate::reference::transcript::Strand::Minus);
        assert_eq!(tx.exons.len(), 2);
        // Exons are stored in transcript order (5'→3'), matching the
        // convention used by cdot and `multi_fasta.rs`. On minus strand the
        // 5'-most exon has the *higher* genomic coordinates.
        assert!(tx.exons[0].genomic_start.unwrap() > tx.exons[1].genomic_start.unwrap());
    }

    #[test]
    fn non_utf8_line_in_lenient_records_diagnostic_and_continues() {
        // Lenient mode must survive a single non-UTF-8 line: previously
        // `reader.lines()` surfaced decode failures as a hard
        // `FerroError::Io`, aborting the entire load. Now the bad line is
        // recorded as `E-LOAD-001` and the surrounding records still load.
        use std::io::Write;
        let mut tf = tempfile::Builder::new().suffix(".gff3").tempfile().unwrap();
        writeln!(tf, "##gff-version 3").unwrap();
        // First good record:
        writeln!(tf, "chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1;Name=g1").unwrap();
        // Inject a non-UTF-8 byte (lone 0xFF) on its own line:
        tf.write_all(b"chr1\t.\texon\t100\t200\t.\t+\t.\tName=bad_\xff\n")
            .unwrap();
        // Second good record:
        writeln!(tf, "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx2;Parent=g1").unwrap();
        writeln!(tf, "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx2").unwrap();
        tf.flush().unwrap();

        let cfg = LoaderConfig::new(); // Lenient by default
        let (db, report) = load_annotations(tf.path(), &cfg).unwrap();
        assert_eq!(
            db.len(),
            1,
            "good records around a bad UTF-8 line must still load"
        );
        assert!(
            report.diagnostics_by_code.contains_key("E-LOAD-001"),
            "non-UTF-8 line must be recorded as E-LOAD-001, got {:?}",
            report.diagnostics_by_code
        );
    }

    #[test]
    fn non_utf8_line_in_strict_returns_err() {
        // Strict mode promotes the non-UTF-8 diagnostic to a hard error,
        // matching the behavior for any other E-LOAD-001 malformed record.
        use std::io::Write;
        let mut tf = tempfile::Builder::new().suffix(".gff3").tempfile().unwrap();
        writeln!(tf, "##gff-version 3").unwrap();
        tf.write_all(b"chr1\t.\texon\t100\t200\t.\t+\t.\tName=bad_\xff\n")
            .unwrap();
        tf.flush().unwrap();
        let cfg = LoaderConfig::new().with_strict();
        assert!(
            load_annotations(tf.path(), &cfg).is_err(),
            "strict mode must reject a file with a non-UTF-8 line"
        );
    }

    #[test]
    fn silent_mode_drops_records_without_recording_diagnostics() {
        // Input mixes three diagnostic sources: a malformed record
        // (E-LOAD-001), an orphan reference (W-LOAD-010), and an
        // unknown-strand transcript (E-LOAD-103 from the builder).
        // Silent mode must record none of them and must not flip
        // `has_error`.
        let f = write_gff3(
            "##gff-version 3\n\
             chr1\t.\tgene\tnot_a_number\t500\t.\t+\t.\tID=g1\n\
             chr1\t.\texon\t100\t500\t.\t+\t.\tParent=missing_parent\n\
             chr1\t.\tmRNA\t100\t500\t.\t.\t.\tID=tx_no_strand\n\
             chr1\t.\texon\t100\t500\t.\t.\t.\tParent=tx_no_strand\n",
        );
        let cfg = LoaderConfig::new().with_silent();
        let (_db, report) = load_annotations(f.path(), &cfg).unwrap();
        assert!(
            report.diagnostics_by_code.is_empty(),
            "silent mode must not aggregate diagnostics by code, got {:?}",
            report.diagnostics_by_code
        );
        assert!(
            report.sample_diagnostics.is_empty(),
            "silent mode must not populate the diagnostic sample"
        );
        assert!(
            !report.has_error,
            "silent mode must not flip `has_error`; strict+silent is a configuration error"
        );
    }

    #[test]
    fn empty_file_loads_to_empty_db() {
        let f = write_gff3("##gff-version 3\n# comment\n");
        let cfg = LoaderConfig::new();
        let (db, report) = load_annotations(f.path(), &cfg).unwrap();
        assert_eq!(db.len(), 0);
        assert_eq!(report.transcripts_loaded, 0);
    }
}
