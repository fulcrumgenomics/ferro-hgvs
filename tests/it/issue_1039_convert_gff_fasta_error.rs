//! #1039: `convert_gff` must not silently emit a wrong transcript sequence when
//! a supplied FASTA cannot cover an exon's genomic coordinates.
//!
//! Two failure shapes exist, and the tolerant default masked both:
//!
//! 1. **Hard extraction error** — the contig is absent from the FASTA, or the
//!    exon starts past the contig end. `FastaProvider::get_sequence` returns
//!    `Err`, and the old code replaced it with `"N".repeat(len)`.
//! 2. **Silent truncation** — the exon overhangs the contig end. The default
//!    `FastaProvider` clamps the read and returns *fewer* real bases than the
//!    exon span (an `Ok` shorter than requested), which the old code pushed
//!    verbatim.
//!
//! Either way the emitted `transcripts.json` carried a plausible-but-wrong
//! reference (a wrong FASTA, a sub-region FASTA paired with whole-genome
//! coordinates, or a mismatched assembly) with no error and no diagnostic. The
//! default path must now fail fast for BOTH shapes (mirroring the `#1026`
//! `--emit-genomic-sequences` coverage check); the tolerant behavior must be
//! reachable only when the user opts out of FASTA validation via
//! `no_validate_fasta`.

use std::io::Write;

use ferro_hgvs::reference::annotation::{convert_gff, ConvertGffConfig};
use ferro_hgvs::reference::transcript::GenomeBuild;
use tempfile::NamedTempFile;

fn write_gff3(content: &str) -> NamedTempFile {
    let mut tf = tempfile::Builder::new().suffix(".gff3").tempfile().unwrap();
    tf.write_all(content.as_bytes()).unwrap();
    tf.flush().unwrap();
    tf
}

/// Write a single-contig FASTA of `len` deterministic bases to a temp file.
fn write_fasta(name: &str, len: usize) -> NamedTempFile {
    let bases = b"ACGT";
    let mut tf = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
    let seq: String = (0..len).map(|i| bases[i % 4] as char).collect();
    writeln!(tf, ">{}", name).unwrap();
    writeln!(tf, "{}", seq).unwrap();
    tf.flush().unwrap();
    tf
}

/// A single-exon transcript on `chr1` spanning genomic 100..500. Paired below
/// with a FASTA that lacks `chr1`, so extracting the exon sequence fails — the
/// stand-in for a wrong FASTA / mismatched assembly.
fn single_exon_chr1_gff3() -> NamedTempFile {
    write_gff3(
        "##gff-version 3\n\
         chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1;Name=GENE1\n\
         chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=g1;gene=GENE1\n\
         chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1\n\
         chr1\t.\tCDS\t150\t450\t.\t+\t0\tParent=tx1\n\
         chr1\t.\tstop_codon\t448\t450\t.\t+\t0\tParent=tx1\n",
    )
}

/// A single-exon transcript on `chr1` spanning genomic 100..900, with its CDS
/// (150..450) safely inside a 520-base contig. Loader CDS/start-codon validation
/// therefore passes, but the exon overhangs the contig end (900 > 520), so the
/// FASTA read is clamped to fewer bases than the exon span — the stand-in for a
/// sub-region FASTA paired with whole-genome coordinates.
fn exon_overhangs_contig_gff3() -> NamedTempFile {
    write_gff3(
        "##gff-version 3\n\
         chr1\t.\tgene\t100\t900\t.\t+\t.\tID=g1;Name=GENE1\n\
         chr1\t.\tmRNA\t100\t900\t.\t+\t.\tID=tx1;Parent=g1;gene=GENE1\n\
         chr1\t.\texon\t100\t900\t.\t+\t.\tParent=tx1\n\
         chr1\t.\tCDS\t150\t450\t.\t+\t0\tParent=tx1\n\
         chr1\t.\tstop_codon\t448\t450\t.\t+\t0\tParent=tx1\n",
    )
}

#[test]
fn extraction_error_fails_fast_by_default() {
    let gff = single_exon_chr1_gff3();
    // FASTA has chr2 only — the exon's chr1 contig is absent (wrong FASTA).
    let fasta = write_fasta("chr2", 520);
    let config = ConvertGffConfig {
        genome_build: GenomeBuild::GRCh38,
        ..ConvertGffConfig::new()
    };

    let err = convert_gff(gff.path(), Some(fasta.path()), &config)
        .expect_err("a FASTA missing the exon's contig must be a hard error, not an N-run");
    let msg = err.to_string();
    assert!(
        msg.contains("tx1"),
        "error should name the offending transcript, got: {msg}"
    );
    assert!(
        msg.contains("does not cover the annotation"),
        "error should explain the coverage mismatch (as the #1026 emit path does), got: {msg}"
    );
}

#[test]
fn extraction_error_tolerated_with_no_validate_fasta() {
    // The user has explicitly opted out of FASTA-aware validation, so the tolerant
    // N-fill is preserved: convert-gff succeeds and the transcript sequence is the
    // all-`N` placeholder covering the exon span (500 - 100 + 1 = 401 bases).
    let gff = single_exon_chr1_gff3();
    let fasta = write_fasta("chr2", 520);
    let config = ConvertGffConfig {
        genome_build: GenomeBuild::GRCh38,
        no_validate_fasta: true,
        ..ConvertGffConfig::new()
    };

    let outcome = convert_gff(gff.path(), Some(fasta.path()), &config)
        .expect("no_validate_fasta opts into the tolerant N-fill");
    let sequence = outcome.json["transcripts"][0]["sequence"]
        .as_str()
        .expect("transcript sequence emitted");
    assert_eq!(sequence.len(), 401, "N-fill covers the full exon span");
    assert!(
        sequence.chars().all(|c| c == 'N'),
        "the tolerated placeholder is an all-N run, got: {sequence}"
    );
}

#[test]
fn partial_coverage_fails_fast_by_default() {
    // The exon spans chr1:100..900 but the contig is only 520 bases, so the FASTA
    // read is clamped to 421 real bases — a silent truncation, not an N-run. It
    // must still fail: the FASTA does not cover the annotation's coordinates.
    let gff = exon_overhangs_contig_gff3();
    let fasta = write_fasta("chr1", 520);
    let config = ConvertGffConfig {
        genome_build: GenomeBuild::GRCh38,
        ..ConvertGffConfig::new()
    };

    let err = convert_gff(gff.path(), Some(fasta.path()), &config)
        .expect_err("a FASTA too short to cover the exon must be a hard error, not a truncation");
    let msg = err.to_string();
    assert!(
        msg.contains("tx1"),
        "error should name the offending transcript, got: {msg}"
    );
    assert!(
        msg.contains("does not cover the annotation"),
        "error should explain the coverage mismatch (as the #1026 emit path does), got: {msg}"
    );
}

#[test]
fn partial_coverage_tolerated_with_no_validate_fasta() {
    // Opted out of FASTA validation, the historical best-effort behavior is
    // preserved: the transcript keeps the real bases the FASTA *does* cover
    // (chr1:100..520 = 421 bases), rather than failing or being blanked to Ns.
    let gff = exon_overhangs_contig_gff3();
    let fasta = write_fasta("chr1", 520);
    let config = ConvertGffConfig {
        genome_build: GenomeBuild::GRCh38,
        no_validate_fasta: true,
        ..ConvertGffConfig::new()
    };

    let outcome = convert_gff(gff.path(), Some(fasta.path()), &config)
        .expect("no_validate_fasta tolerates partial coverage");
    let sequence = outcome.json["transcripts"][0]["sequence"]
        .as_str()
        .expect("transcript sequence emitted");
    // The FASTA fills base `ACGT[i % 4]` at each absolute 0-based index `i`. The
    // clamped read covers indices 99..=519 (chr1:100..520): index 99 → `T`, then
    // indices 100..=519 are 420 bases = `ACGT` × 105. Pin the exact bases so a
    // wrong-offset extraction (not just an N-run) would fail the test too.
    let expected: String = format!("T{}", "ACGT".repeat(105));
    assert_eq!(
        sequence, expected,
        "the covered bases (chr1:100..520) are kept as-is, in order"
    );
}
