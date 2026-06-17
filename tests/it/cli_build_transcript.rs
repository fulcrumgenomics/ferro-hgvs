//! Integration tests for the `ferro build-transcript` subcommand (#184).

use std::io::Write;
use std::process::Command;
use tempfile::NamedTempFile;

/// Write a single-contig FASTA to a temporary file.
///
/// Returns the `NamedTempFile` (to keep it alive) and its path.
fn write_fasta(name: &str, seq: &str) -> (NamedTempFile, std::path::PathBuf) {
    let mut tf = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
    writeln!(tf, ">{}", name).unwrap();
    writeln!(tf, "{}", seq).unwrap();
    tf.flush().unwrap();
    let path = tf.path().to_path_buf();
    (tf, path)
}

#[test]
fn build_transcript_emits_valid_json_for_single_exon() {
    // 60 bp synthetic construct, CDS 1..60.
    let (_keep, fasta_path) = write_fasta(
        "construct1",
        "ATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCTAA",
    );
    let out = tempfile::Builder::new().suffix(".json").tempfile().unwrap();
    let out_path = out.path().to_path_buf();
    drop(out); // close handle so build-transcript can overwrite

    let bin = env!("CARGO_BIN_EXE_ferro");
    let status = Command::new(bin)
        .args([
            "build-transcript",
            "--fasta",
            fasta_path.to_str().unwrap(),
            "--cds-start",
            "1",
            "--cds-end",
            "60",
            "--output",
            out_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        status.status.success(),
        "build-transcript failed: stderr={}",
        String::from_utf8_lossy(&status.stderr)
    );

    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&out_path).unwrap()).unwrap();
    assert_eq!(json["version"], "1.0");
    assert_eq!(json["genome_build"], "GRCh38");
    let tx = &json["transcripts"][0];
    assert_eq!(tx["id"], "construct1");
    assert_eq!(tx["chromosome"], "construct1");
    assert_eq!(tx["strand"], "+");
    assert_eq!(tx["cds_start"], 1);
    assert_eq!(tx["cds_end"], 60);
    assert_eq!(
        tx["sequence"],
        "ATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCTAA"
    );
    assert_eq!(tx["genomic_start"], 1);
    assert_eq!(tx["genomic_end"], 60);
    assert_eq!(tx["exons"][0]["genomic_start"], 1);
    assert_eq!(tx["exons"][0]["genomic_end"], 60);
    assert_eq!(tx["exons"][0]["start"], 1);
    assert_eq!(tx["exons"][0]["end"], 60);
    assert_eq!(tx["exons"][0]["number"], 1);
}

#[test]
fn build_transcript_rejects_cds_beyond_contig() {
    let (_keep, fasta_path) = write_fasta("small", "ATGAAA");
    let out = tempfile::Builder::new().suffix(".json").tempfile().unwrap();
    let out_path = out.path().to_path_buf();
    drop(out);

    let bin = env!("CARGO_BIN_EXE_ferro");
    let status = Command::new(bin)
        .args([
            "build-transcript",
            "--fasta",
            fasta_path.to_str().unwrap(),
            "--cds-start",
            "1",
            "--cds-end",
            "100",
            "--output",
            out_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        !status.status.success(),
        "expected failure for CDS beyond contig"
    );
}

#[test]
fn build_transcript_supports_minus_strand_and_custom_id() {
    let (_keep, fasta_path) = write_fasta("construct2", "ATGAAACCCGGGTTTAAACCCGGGTTTTAA");
    let out = tempfile::Builder::new().suffix(".json").tempfile().unwrap();
    let out_path = out.path().to_path_buf();
    drop(out);

    let bin = env!("CARGO_BIN_EXE_ferro");
    let status = Command::new(bin)
        .args([
            "build-transcript",
            "--fasta",
            fasta_path.to_str().unwrap(),
            "--cds-start",
            "1",
            "--cds-end",
            "30",
            "--strand",
            "-",
            "--id",
            "tx_custom",
            "--gene",
            "MYGENE",
            "--output",
            out_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        status.status.success(),
        "build-transcript failed: stderr={}",
        String::from_utf8_lossy(&status.stderr)
    );

    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&out_path).unwrap()).unwrap();
    let tx = &json["transcripts"][0];
    assert_eq!(tx["id"], "tx_custom");
    assert_eq!(tx["gene_symbol"], "MYGENE");
    assert_eq!(tx["strand"], "-");
    // Sequence should be reverse-complemented
    // Original: ATGAAACCCGGGTTTAAACCCGGGTTTTAA (30 bp)
    // RC: TTAAAACCCGGGTTTAAACCCGGGTTTCAT
    assert_eq!(tx["sequence"], "TTAAAACCCGGGTTTAAACCCGGGTTTCAT");
}
