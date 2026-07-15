//! The extracted library `build_transcript`
//! (`ferro_hgvs::reference::annotation::build_transcript`) must produce output
//! byte-identical to the `ferro build-transcript` CLI for the same inputs and
//! flags. The CLI is now a thin wrapper over this function, and the Python
//! binding wraps the same function — so this test pins the shared serializer
//! that guarantees CLI ↔ library ↔ Python parity (#1035).

use std::io::Write;
use std::path::Path;
use std::process::Command;

use ferro_hgvs::reference::annotation::{build_transcript, BuildTranscriptConfig};
use tempfile::NamedTempFile;

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

/// Run `ferro build-transcript` and return the exact bytes it writes to `-o`.
fn cli_output_bytes(fasta: &Path, extra: &[&str]) -> Vec<u8> {
    let out = tempfile::Builder::new()
        .suffix(".json")
        .tempfile()
        .unwrap()
        .into_temp_path();
    let bin = env!("CARGO_BIN_EXE_ferro");
    let mut args: Vec<&str> = vec!["build-transcript", "--fasta", fasta.to_str().unwrap()];
    args.extend_from_slice(extra);
    args.extend_from_slice(&["-o", out.to_str().unwrap()]);
    let result = Command::new(bin).args(&args).output().unwrap();
    assert!(
        result.status.success(),
        "build-transcript failed: stderr={}",
        String::from_utf8_lossy(&result.stderr)
    );
    std::fs::read(&out).unwrap()
}

/// Serialize a library `build_transcript` outcome exactly as the CLI writes it:
/// pretty JSON with NO trailing newline (the CLI's `std::fs::write`).
fn library_output_bytes(fasta: &Path, config: &BuildTranscriptConfig) -> Vec<u8> {
    let outcome = build_transcript(fasta, config).expect("library build_transcript succeeds");
    serde_json::to_string_pretty(&outcome.json)
        .unwrap()
        .into_bytes()
}

#[test]
fn library_matches_cli_plus_strand() {
    let fasta = write_fasta("construct1", 60);
    let config = BuildTranscriptConfig::new(1, 60);

    let cli = cli_output_bytes(fasta.path(), &["--cds-start", "1", "--cds-end", "60"]);
    let lib = library_output_bytes(fasta.path(), &config);

    assert_eq!(
        String::from_utf8_lossy(&lib),
        String::from_utf8_lossy(&cli),
        "library output must be byte-identical to the CLI"
    );
}

#[test]
fn library_matches_cli_minus_strand_custom_id_gene() {
    let fasta = write_fasta("chrX", 90);
    let config = BuildTranscriptConfig {
        id: Some("MYTX.1".to_string()),
        strand: "-".to_string(),
        gene: Some("GENEX".to_string()),
        ..BuildTranscriptConfig::new(4, 87)
    };

    let cli = cli_output_bytes(
        fasta.path(),
        &[
            "--cds-start",
            "4",
            "--cds-end",
            "87",
            "--strand",
            "-",
            "--id",
            "MYTX.1",
            "--gene",
            "GENEX",
        ],
    );
    let lib = library_output_bytes(fasta.path(), &config);

    assert_eq!(
        String::from_utf8_lossy(&lib),
        String::from_utf8_lossy(&cli),
        "library output (minus strand, custom id/gene) must be byte-identical to the CLI"
    );
}

#[test]
fn library_matches_cli_with_emit_genomic_sequences() {
    let fasta = write_fasta("construct1", 120);
    let config = BuildTranscriptConfig {
        emit_genomic_sequences: true,
        ..BuildTranscriptConfig::new(1, 120)
    };

    let cli = cli_output_bytes(
        fasta.path(),
        &[
            "--cds-start",
            "1",
            "--cds-end",
            "120",
            "--emit-genomic-sequences",
        ],
    );
    let lib = library_output_bytes(fasta.path(), &config);

    assert_eq!(
        String::from_utf8_lossy(&lib),
        String::from_utf8_lossy(&cli),
        "library output with --emit-genomic-sequences must be byte-identical to the CLI"
    );
}
