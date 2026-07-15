//! Integration test for `ferro convert-gff` flags added in Phase 3.

use std::io::Write;
use std::process::Command;
use tempfile::NamedTempFile;

fn write_gff3(content: &str) -> NamedTempFile {
    let mut tf = tempfile::Builder::new().suffix(".gff3").tempfile().unwrap();
    tf.write_all(content.as_bytes()).unwrap();
    tf.flush().unwrap();
    tf
}

/// Write a single-contig FASTA of `len` deterministic bases to a temp file.
fn write_fasta(name: &str, len: usize) -> NamedTempFile {
    write_multi_fasta(&[(name, len)])
}

/// Write a multi-contig FASTA (each `(name, len)`) of deterministic bases.
fn write_multi_fasta(contigs: &[(&str, usize)]) -> NamedTempFile {
    let bases = b"ACGT";
    let mut tf = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
    for (name, len) in contigs {
        let seq: String = (0..*len).map(|i| bases[i % 4] as char).collect();
        writeln!(tf, ">{}", name).unwrap();
        writeln!(tf, "{}", seq).unwrap();
    }
    tf.flush().unwrap();
    tf
}

#[test]
fn strict_mode_fails_on_malformed_coordinate() {
    let f = write_gff3(
        "##gff-version 3\n\
         chr1\t.\tgene\tnot_a_number\t500\t.\t+\t.\tID=g1\n",
    );
    let bin = env!("CARGO_BIN_EXE_ferro");
    let out = Command::new(bin)
        .args([
            "convert-gff",
            "--gff",
            f.path().to_str().unwrap(),
            "--strict",
        ])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "strict mode should fail on malformed record: stdout={}, stderr={}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn diagnostics_json_writes_file() {
    let f = write_gff3(
        "##gff-version 3\n\
         chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1\n",
    );
    let json_path = tempfile::Builder::new()
        .suffix(".json")
        .tempfile()
        .unwrap()
        .into_temp_path();
    let bin = env!("CARGO_BIN_EXE_ferro");
    let out = Command::new(bin)
        .args([
            "convert-gff",
            "--gff",
            f.path().to_str().unwrap(),
            "--diagnostics-json",
            json_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "convert-gff failed: stderr={}",
        String::from_utf8_lossy(&out.stderr)
    );
    let body = std::fs::read_to_string(&json_path).expect("diagnostics file");
    let diags: Vec<serde_json::Value> =
        serde_json::from_str(&body).expect("diagnostics JSON parses as array");
    assert!(
        diags
            .iter()
            .any(|d| d.get("code").and_then(|v| v.as_str()) == Some("W-LOAD-100")),
        "expected diagnostic with code=W-LOAD-100, got: {}",
        body
    );
}

/// A single-exon transcript on chr1, plus a chr1 FASTA long enough to cover it.
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

#[test]
fn emit_genomic_sequences_makes_output_genome_capable() {
    let gff = single_exon_gff3();
    let fasta = write_fasta("chr1", 520);
    let out = tempfile::Builder::new()
        .suffix(".json")
        .tempfile()
        .unwrap()
        .into_temp_path();

    let bin = env!("CARGO_BIN_EXE_ferro");
    let result = Command::new(bin)
        .args([
            "convert-gff",
            "--gff",
            gff.path().to_str().unwrap(),
            "--fasta",
            fasta.path().to_str().unwrap(),
            "--no-validate-fasta",
            "--emit-genomic-sequences",
            "--output",
            out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "convert-gff failed: stderr={}",
        String::from_utf8_lossy(&result.stderr)
    );

    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&out).unwrap()).unwrap();
    let genomic = &json["genomic_sequences"]["chr1"];
    let genomic = genomic.as_str().expect("chr1 genomic sequence emitted");
    // Full contig bytes are emitted, matching the FASTA contig exactly (not just
    // its length — a wrong-offset or N-filled slice would share the length).
    let expected: String = (0..520).map(|i| b"ACGT"[i % 4] as char).collect();
    assert_eq!(genomic, expected);
    // The exon spans genomic 100..500, so the transcript placement is fully backed.
    let provider =
        ferro_hgvs::reference::mock::MockProvider::from_json(&out).expect("loads and validates");
    use ferro_hgvs::reference::provider::ReferenceProvider;
    assert!(provider.has_genomic_data());
}

#[test]
fn without_emit_flag_no_genomic_sequences() {
    let gff = single_exon_gff3();
    let fasta = write_fasta("chr1", 520);
    let out = tempfile::Builder::new()
        .suffix(".json")
        .tempfile()
        .unwrap()
        .into_temp_path();

    let bin = env!("CARGO_BIN_EXE_ferro");
    let result = Command::new(bin)
        .args([
            "convert-gff",
            "--gff",
            gff.path().to_str().unwrap(),
            "--fasta",
            fasta.path().to_str().unwrap(),
            "--no-validate-fasta",
            "--output",
            out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(result.status.success());

    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&out).unwrap()).unwrap();
    assert!(
        json.get("genomic_sequences").is_none(),
        "genomic_sequences must be absent without --emit-genomic-sequences"
    );
    // A transcripts-only reference stays transcript-capable only.
    let provider = ferro_hgvs::reference::mock::MockProvider::from_json(&out).unwrap();
    use ferro_hgvs::reference::provider::ReferenceProvider;
    assert!(!provider.has_genomic_data());
}

#[test]
fn emit_genomic_sequences_warns_when_nothing_is_placed() {
    // A transcript filter that matches nothing leaves no emitted transcript with a
    // chromosome, so --emit-genomic-sequences has nothing to back: convert-gff must
    // succeed, warn, and NOT write an empty (misleading) genomic_sequences map.
    let gff = single_exon_gff3();
    let fasta = write_fasta("chr1", 520);
    let out = tempfile::Builder::new()
        .suffix(".json")
        .tempfile()
        .unwrap()
        .into_temp_path();

    let bin = env!("CARGO_BIN_EXE_ferro");
    let result = Command::new(bin)
        .args([
            "convert-gff",
            "--gff",
            gff.path().to_str().unwrap(),
            "--fasta",
            fasta.path().to_str().unwrap(),
            "--no-validate-fasta",
            "--transcripts",
            "NO_SUCH_TRANSCRIPT",
            "--emit-genomic-sequences",
            "--output",
            out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "convert-gff failed: stderr={}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(
        String::from_utf8_lossy(&result.stderr)
            .contains("no emitted transcript has a genomic placement"),
        "expected a warning about nothing to place, stderr={}",
        String::from_utf8_lossy(&result.stderr)
    );

    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&out).unwrap()).unwrap();
    assert!(
        json.get("genomic_sequences").is_none(),
        "an empty genomic_sequences map must not be written"
    );
}

#[test]
fn emit_genomic_sequences_fails_fast_when_fasta_too_short() {
    // The GFF places the exon at chr1:100..500, but the FASTA contig is only 200
    // bases (a sub-region FASTA with whole-genome coordinates). convert-gff must
    // FAIL AT EMIT — not produce a file the loader would later reject (#1026).
    let gff = single_exon_gff3(); // exon 100..500 on chr1
    let fasta = write_fasta("chr1", 200);
    let out = tempfile::Builder::new()
        .suffix(".json")
        .tempfile()
        .unwrap()
        .into_temp_path();

    let bin = env!("CARGO_BIN_EXE_ferro");
    let result = Command::new(bin)
        .args([
            "convert-gff",
            "--gff",
            gff.path().to_str().unwrap(),
            "--fasta",
            fasta.path().to_str().unwrap(),
            "--no-validate-fasta",
            "--emit-genomic-sequences",
            "--output",
            out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        !result.status.success(),
        "convert-gff must fail when the FASTA cannot cover the placement"
    );
    assert!(
        String::from_utf8_lossy(&result.stderr).contains("does not cover the annotation"),
        "expected a coordinate-coverage error, got: {}",
        String::from_utf8_lossy(&result.stderr)
    );
}

#[test]
fn emit_genomic_sequences_handles_multiple_contigs() {
    // Two transcripts on two contigs; both contigs must be emitted and the result
    // must load as genome-capable.
    let gff = write_gff3(
        "##gff-version 3\n\
         chr1\t.\tgene\t100\t200\t.\t+\t.\tID=g1\n\
         chr1\t.\tmRNA\t100\t200\t.\t+\t.\tID=tx1;Parent=g1\n\
         chr1\t.\texon\t100\t200\t.\t+\t.\tParent=tx1\n\
         chr2\t.\tgene\t50\t150\t.\t+\t.\tID=g2\n\
         chr2\t.\tmRNA\t50\t150\t.\t+\t.\tID=tx2;Parent=g2\n\
         chr2\t.\texon\t50\t150\t.\t+\t.\tParent=tx2\n",
    );
    let fasta = write_multi_fasta(&[("chr1", 260), ("chr2", 200)]);
    let out = tempfile::Builder::new()
        .suffix(".json")
        .tempfile()
        .unwrap()
        .into_temp_path();

    let bin = env!("CARGO_BIN_EXE_ferro");
    let result = Command::new(bin)
        .args([
            "convert-gff",
            "--gff",
            gff.path().to_str().unwrap(),
            "--fasta",
            fasta.path().to_str().unwrap(),
            "--no-validate-fasta",
            "--emit-genomic-sequences",
            "--output",
            out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "convert-gff failed: stderr={}",
        String::from_utf8_lossy(&result.stderr)
    );

    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&out).unwrap()).unwrap();
    // Assert exact per-contig bytes (using the same deterministic generator as
    // write_multi_fasta), so a contig-swap or wrong-offset read can't pass.
    let expect = |len: usize| -> String { (0..len).map(|i| b"ACGT"[i % 4] as char).collect() };
    assert_eq!(
        json["genomic_sequences"]["chr1"].as_str().unwrap(),
        expect(260)
    );
    assert_eq!(
        json["genomic_sequences"]["chr2"].as_str().unwrap(),
        expect(200)
    );

    let provider = ferro_hgvs::reference::mock::MockProvider::from_json(&out).expect("loads");
    use ferro_hgvs::reference::provider::ReferenceProvider;
    assert!(provider.has_genomic_data());
}

#[test]
fn emit_genomic_sequences_requires_fasta() {
    let gff = single_exon_gff3();
    let bin = env!("CARGO_BIN_EXE_ferro");
    let result = Command::new(bin)
        .args([
            "convert-gff",
            "--gff",
            gff.path().to_str().unwrap(),
            "--emit-genomic-sequences",
        ])
        .output()
        .unwrap();
    assert!(
        !result.status.success(),
        "--emit-genomic-sequences without --fasta should be rejected by the CLI"
    );
}
