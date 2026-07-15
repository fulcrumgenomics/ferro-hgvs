//! The extracted library `convert_gff` (`ferro_hgvs::reference::annotation::convert_gff`)
//! must produce output byte-identical to the `ferro convert-gff` CLI for the
//! same inputs and flags. The CLI is now a thin wrapper over this function, and
//! the Python binding wraps the same function — so this test pins the shared
//! serializer that guarantees CLI ↔ library ↔ Python parity (#1035).

use std::io::Write;
use std::path::Path;
use std::process::Command;

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

/// A single-exon transcript on chr1 with a CDS.
fn single_exon_gff3() -> NamedTempFile {
    write_gff3(
        "##gff-version 3\n\
         chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1;Name=GENE1\n\
         chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=g1;gene=GENE1\n\
         chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1\n\
         chr1\t.\tCDS\t150\t450\t.\t+\t0\tParent=tx1\n\
         chr1\t.\tstop_codon\t448\t450\t.\t+\t0\tParent=tx1\n",
    )
}

/// Run `ferro convert-gff` and return the exact bytes it writes to `--output`.
fn cli_output_bytes(args: &[&str]) -> Vec<u8> {
    let out = tempfile::Builder::new()
        .suffix(".json")
        .tempfile()
        .unwrap()
        .into_temp_path();
    let bin = env!("CARGO_BIN_EXE_ferro");
    let mut full: Vec<&str> = vec!["convert-gff"];
    full.extend_from_slice(args);
    full.extend_from_slice(&["--output", out.to_str().unwrap()]);
    let result = Command::new(bin).args(&full).output().unwrap();
    assert!(
        result.status.success(),
        "convert-gff failed: stderr={}",
        String::from_utf8_lossy(&result.stderr)
    );
    std::fs::read(&out).unwrap()
}

/// Serialize a library `convert_gff` outcome exactly as the CLI writes it:
/// pretty JSON followed by a trailing newline (the CLI's `writeln!`).
fn library_output_bytes(gff: &Path, fasta: Option<&Path>, config: &ConvertGffConfig) -> Vec<u8> {
    let outcome = convert_gff(gff, fasta, config).expect("library convert_gff succeeds");
    let mut bytes = serde_json::to_string_pretty(&outcome.json)
        .unwrap()
        .into_bytes();
    bytes.push(b'\n');
    bytes
}

#[test]
fn library_matches_cli_transcript_only() {
    let gff = single_exon_gff3();
    let config = ConvertGffConfig {
        genome_build: GenomeBuild::GRCh38,
        ..ConvertGffConfig::new()
    };

    let cli = cli_output_bytes(&["--gff", gff.path().to_str().unwrap()]);
    let lib = library_output_bytes(gff.path(), None, &config);

    assert_eq!(
        lib,
        cli,
        "library output must be byte-identical to the CLI\nlib:\n{}\ncli:\n{}",
        String::from_utf8_lossy(&lib),
        String::from_utf8_lossy(&cli),
    );
}

#[test]
fn library_matches_cli_with_emit_genomic_sequences() {
    let gff = single_exon_gff3();
    let fasta = write_fasta("chr1", 520);
    let config = ConvertGffConfig {
        genome_build: GenomeBuild::GRCh38,
        no_validate_fasta: true,
        emit_genomic_sequences: true,
        ..ConvertGffConfig::new()
    };

    let cli = cli_output_bytes(&[
        "--gff",
        gff.path().to_str().unwrap(),
        "--fasta",
        fasta.path().to_str().unwrap(),
        "--no-validate-fasta",
        "--emit-genomic-sequences",
    ]);
    let lib = library_output_bytes(gff.path(), Some(fasta.path()), &config);

    assert_eq!(
        lib,
        cli,
        "library output with --emit-genomic-sequences must be byte-identical to the CLI\nlib:\n{}\ncli:\n{}",
        String::from_utf8_lossy(&lib),
        String::from_utf8_lossy(&cli),
    );
}

#[test]
fn library_matches_cli_with_gene_filter() {
    let gff = single_exon_gff3();
    let config = ConvertGffConfig {
        genome_build: GenomeBuild::GRCh38,
        genes: Some(vec!["GENE1".to_string()]),
        ..ConvertGffConfig::new()
    };

    let cli = cli_output_bytes(&["--gff", gff.path().to_str().unwrap(), "--genes", "GENE1"]);
    let lib = library_output_bytes(gff.path(), None, &config);

    assert_eq!(
        lib,
        cli,
        "library output with a gene filter must be byte-identical to the CLI\nlib:\n{}\ncli:\n{}",
        String::from_utf8_lossy(&lib),
        String::from_utf8_lossy(&cli),
    );

    // The filter is real: a non-matching gene yields zero transcripts.
    let empty = ConvertGffConfig {
        genome_build: GenomeBuild::GRCh38,
        genes: Some(vec!["NOSUCHGENE".to_string()]),
        ..ConvertGffConfig::new()
    };
    let outcome = convert_gff(gff.path(), None, &empty).unwrap();
    assert_eq!(outcome.json["transcripts"].as_array().unwrap().len(), 0);
}
