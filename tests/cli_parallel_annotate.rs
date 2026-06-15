//! Integration test for parallel batch `ferro annotate-vcf` (`-j/--workers`).
//!
//! Invariance property: the annotated output for any worker count must be
//! byte-identical to the serial (`-j1`) output, including record ordering across
//! chunk boundaries. Uses the committed micro chr21 GFF fixture (one transcript,
//! `NR_027228.1`, two exons) so the test needs no external data.

use std::io::Write;
use std::process::Command;
use tempfile::NamedTempFile;

// Exon ranges of NR_027228.1 in the micro GFF.
const EXONS: &[(u64, u64)] = &[(5_016_902, 5_017_178), (5_011_799, 5_012_256)];
const BASES: [&str; 4] = ["A", "C", "G", "T"];

/// Write a VCF with `n` distinct records placed inside the transcript's exons
/// (so every record annotates). `n` exceeds the driver's 4096-record chunk size
/// so ordering is exercised across multiple chunks. Each record has a unique ID
/// so any reordering changes the output byte stream.
fn write_vcf(n: usize) -> (NamedTempFile, std::path::PathBuf) {
    let mut tf = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
    writeln!(tf, "##fileformat=VCFv4.2").unwrap();
    writeln!(tf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    // Walk positions across both exons, cycling as needed to reach `n`, but give
    // every record a unique ID so duplicate positions stay distinguishable.
    let mut written = 0usize;
    let mut id = 0usize;
    while written < n {
        for &(start, end) in EXONS {
            for pos in start..=end {
                if written >= n {
                    break;
                }
                let r = BASES[(pos % 4) as usize];
                let a = BASES[((pos + 1) % 4) as usize];
                writeln!(tf, "NC_000021.9\t{pos}\tv{id}\t{r}\t{a}\t.\tPASS\t.").unwrap();
                written += 1;
                id += 1;
            }
        }
    }
    tf.flush().unwrap();
    let path = tf.path().to_path_buf();
    (tf, path)
}

fn annotate(vcf: &std::path::Path, workers: usize) -> Vec<u8> {
    let bin = env!("CARGO_BIN_EXE_ferro");
    let out = Command::new(bin)
        .args([
            "annotate-vcf",
            "-i",
            vcf.to_str().unwrap(),
            "--gff",
            "tests/fixtures/annotation/refseq_chr21_micro.gff3",
            "--build",
            "GRCh38",
            "-j",
            &workers.to_string(),
            "-o",
            "-",
        ])
        .output()
        .expect("run ferro annotate-vcf");
    assert!(
        out.status.success(),
        "annotate-vcf -j{workers} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

#[test]
fn parallel_annotate_output_is_byte_identical_to_serial() {
    let (_keep, vcf) = write_vcf(10_000); // > 2 chunks
    let serial = annotate(&vcf, 1);
    // Sanity: every record actually got annotated (one `HGVS_C=` per record),
    // not just one record passing through.
    assert_eq!(
        String::from_utf8_lossy(&serial).matches("HGVS_C=").count(),
        10_000,
        "expected all 10,000 records to be annotated"
    );
    for workers in [2usize, 4, 8] {
        assert_eq!(
            serial,
            annotate(&vcf, workers),
            "annotate-vcf -j{workers} output differs from serial (-j1)"
        );
    }
}

#[test]
fn parallel_annotate_is_deterministic() {
    let (_keep, vcf) = write_vcf(8_000);
    assert_eq!(annotate(&vcf, 8), annotate(&vcf, 8));
}
