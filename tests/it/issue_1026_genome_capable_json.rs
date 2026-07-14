//! Issue #1026: a `transcripts.json` that carries a top-level `genomic_sequences`
//! map (as emitted by `ferro convert-gff --emit-genomic-sequences` /
//! `ferro build-transcript --emit-genomic-sequences`) is genome-capable: a
//! `Normalizer` built from it via `MockProvider::from_json` runs the genome-aware
//! normalization rules, whereas the same reference without `genomic_sequences`
//! does not.
//!
//! `Normalizer(reference_json=…)` (Python) and
//! `Normalizer::new(MockProvider::from_json(…))` (Rust) go through the same
//! normalize pipeline, so proving the Rust path covers the Python binding too.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::provider::ReferenceProvider;
use ferro_hgvs::{parse_hgvs, FerroError, Normalizer};
use std::io::Write;
use tempfile::NamedTempFile;

/// Build a two-exon, plus-strand non-coding transcript with a small intron and
/// (optionally) the backing genomic contig, then serialize it to a temp JSON file
/// in the same object form the CLI emits.
///
/// Geometry (plus strand):
///   Exon 1: tx 1..20,  genomic 1001..1020
///   Exon 2: tx 21..40, genomic 1101..1120   (intron genomic 1021..1100)
fn write_reference_json(with_genomic_sequences: bool) -> NamedTempFile {
    // 40-base transcript sequence.
    let tx_seq = "ACGT".repeat(10);

    let transcript = serde_json::json!({
        "id": "NM_1026TEST.1",
        "gene_symbol": "GENE1026",
        "strand": "+",
        "sequence": tx_seq,
        "chromosome": "chrT",
        "genomic_start": 1001,
        "genomic_end": 1120,
        "exons": [
            {"number": 1, "start": 1,  "end": 20, "genomic_start": 1001, "genomic_end": 1020},
            {"number": 2, "start": 21, "end": 40, "genomic_start": 1101, "genomic_end": 1120}
        ]
    });

    let mut obj = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [transcript],
    });

    if with_genomic_sequences {
        // 1200-base contig, long enough to back the placement (max coord 1120).
        let contig = "ACGT".repeat(300);
        let mut genomic = serde_json::Map::new();
        genomic.insert("chrT".to_string(), serde_json::json!(contig));
        obj["genomic_sequences"] = serde_json::Value::Object(genomic);
    }

    let mut file = NamedTempFile::new().unwrap();
    file.write_all(serde_json::to_string_pretty(&obj).unwrap().as_bytes())
        .unwrap();
    file.flush().unwrap();
    file
}

/// Without `genomic_sequences`, the reference is transcript-only: `has_genomic_data`
/// is false and an intronic normalization cannot run the genome-aware path — it
/// reports the intronic variant as unresolvable (the pre-#1015 main behavior).
#[test]
fn transcripts_only_reference_cannot_normalize_intronic() {
    let file = write_reference_json(false);
    let provider = MockProvider::from_json(file.path()).unwrap();
    assert!(
        !provider.has_genomic_data(),
        "a reference without genomic_sequences must not report genomic data"
    );

    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NM_1026TEST.1:n.20+3del").expect("parse");
    let result = normalizer.normalize(&variant);
    assert!(
        matches!(result, Err(FerroError::IntronicVariant { .. })),
        "transcript-only reference should report the intronic variant as unresolvable, got: {result:?}"
    );
}

/// With `genomic_sequences` present, the same reference is genome-capable: the
/// intronic normalization reaches the genome-aware path and no longer fails with
/// `IntronicVariant`. This is exactly the capability `--emit-genomic-sequences`
/// unlocks (#1026).
#[test]
fn genome_capable_reference_normalizes_intronic() {
    let file = write_reference_json(true);
    let provider = MockProvider::from_json(file.path()).unwrap();
    assert!(
        provider.has_genomic_data(),
        "a reference with genomic_sequences must report genomic data"
    );

    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NM_1026TEST.1:n.20+3del").expect("parse");
    // The genome-aware intronic path runs and SUCCEEDS (vs the transcript-only
    // reference above, which errors with IntronicVariant). On this synthetic
    // reference the deletion is already in canonical position, so the normalized
    // value is the input unchanged — but critically it is an Ok result from the
    // genome-aware path, not an error.
    let normalized = normalizer
        .normalize(&variant)
        .expect("genome-capable reference must reach the genome-aware path and succeed");
    assert_eq!(normalized.to_string(), "NM_1026TEST.1:n.20+3del");
}
