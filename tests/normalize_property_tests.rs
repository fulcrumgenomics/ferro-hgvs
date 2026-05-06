//! Transcript-driven normalize() property tests.
//!
//! Two fundamental properties:
//!   1. Normalized output is itself valid HGVS (re-parses cleanly).
//!   2. Normalization is idempotent: normalize(normalize(x)) == normalize(x).
//!
//! These tests use real transcript metadata from
//! tests/fixtures/sequences/normalization_transcripts.json - they are
//! independent of the v21.0 spec fixture (which uses MockProvider::new()).
//!
//! Originally lived in tests/hgvs_spec_compliance_tests.rs (deleted in #84).

use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::Normalizer;
use serde::Deserialize;
use std::collections::HashMap;
use std::sync::OnceLock;

#[derive(Debug, Deserialize, Clone)]
struct TranscriptData {
    id: String,
    gene_symbol: String,
    strand: String,
    sequence: String,
    cds_start: u64,
    cds_end: u64,
    exons: Vec<ExonData>,
}

#[derive(Debug, Deserialize, Clone)]
struct ExonData {
    number: u32,
    start: u64,
    end: u64,
}

static TRANSCRIPTS: OnceLock<HashMap<String, TranscriptData>> = OnceLock::new();

fn get_transcripts() -> &'static HashMap<String, TranscriptData> {
    TRANSCRIPTS.get_or_init(|| {
        let path = format!(
            "{}/tests/fixtures/sequences/normalization_transcripts.json",
            env!("CARGO_MANIFEST_DIR")
        );
        let json =
            std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
        let txs: Vec<TranscriptData> =
            serde_json::from_str(&json).unwrap_or_else(|e| panic!("failed to parse {path}: {e}"));
        txs.into_iter().map(|t| (t.id.clone(), t)).collect()
    })
}

fn provider_for(accession: &str) -> Option<MockProvider> {
    let d = get_transcripts().get(accession)?;
    let mut provider = MockProvider::new();
    let strand = if d.strand == "+" {
        Strand::Plus
    } else {
        Strand::Minus
    };
    let exons: Vec<Exon> = d
        .exons
        .iter()
        .map(|e| Exon::new(e.number, e.start, e.end))
        .collect();
    let tx = Transcript::new(
        d.id.clone(),
        Some(d.gene_symbol.clone()),
        strand,
        d.sequence.clone(),
        Some(d.cds_start),
        Some(d.cds_end),
        exons,
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(tx);
    Some(provider)
}

const TRANSCRIPT_INPUTS: &[&str] = &[
    "NM_000492.4:c.3846G>A",
    "NM_007294.4:c.5266A>G",
    "NM_000546.6:c.215C>G",
    "NM_000546.6:c.743G>A",
    "NM_000492.4:c.1521_1523del",
    "NM_007294.4:c.68_69del",
    "NM_000059.4:c.5946del",
    "NM_000546.6:c.532del",
    "NM_007294.4:c.5266dup",
    "NM_001089.2:c.3057_3058insT",
    "NR_153405.1:n.3650G>A",
    "NR_153405.1:n.3799C>T",
];

#[test]
fn normalized_output_is_valid_hgvs() {
    let mut failures = Vec::new();
    let mut tested = 0usize;
    let mut skipped = 0usize;
    for input in TRANSCRIPT_INPUTS {
        let acc = match input.split(':').next() {
            Some(a) => a,
            None => {
                skipped += 1;
                continue;
            }
        };
        let provider = match provider_for(acc) {
            Some(p) => p,
            None => {
                skipped += 1;
                continue;
            }
        };
        let v = match parse_hgvs(input) {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!("parse failed for {input}: {e}"));
                continue;
            }
        };
        let n = match Normalizer::new(provider).normalize(&v) {
            Ok(n) => n,
            Err(e) => {
                failures.push(format!("normalize failed for {input}: {e}"));
                continue;
            }
        };
        let f = format!("{n}");
        if let Err(e) = parse_hgvs(&f) {
            failures.push(format!("{input} -> {f}: {e}"));
        }
        tested += 1;
    }
    eprintln!("normalized_output_is_valid_hgvs: tested {tested}, skipped {skipped}");
    assert!(
        tested > 0,
        "normalized_output_is_valid_hgvs exercised no cases (all skipped)"
    );
    assert!(
        failures.is_empty(),
        "Normalized output failed to re-parse:\n{}",
        failures.join("\n")
    );
}

#[test]
fn normalization_is_idempotent() {
    const INPUTS: &[&str] = &[
        "NM_000492.4:c.1521_1523del",
        "NM_000492.4:c.3846G>A",
        "NM_007294.4:c.68_69del",
        "NM_007294.4:c.5266dup",
        "NM_000059.4:c.5946del",
        "NM_000546.6:c.215C>G",
        "NM_000546.6:c.532del",
        "NR_153405.1:n.3650G>A",
    ];
    let mut failures = Vec::new();
    let mut tested = 0usize;
    let mut skipped = 0usize;
    for input in INPUTS {
        let acc = input.split(':').next().unwrap();
        let provider = match provider_for(acc) {
            Some(p) => p,
            None => {
                skipped += 1;
                continue;
            }
        };
        let v = match parse_hgvs(input) {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!("parse failed for {input}: {e}"));
                continue;
            }
        };
        let n1 = match Normalizer::new(provider.clone()).normalize(&v) {
            Ok(n) => n,
            Err(e) => {
                failures.push(format!("normalize failed for {input}: {e}"));
                continue;
            }
        };
        let f1 = format!("{n1}");
        let v2 = match parse_hgvs(&f1) {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!("re-parse {input} -> {f1}: {e}"));
                continue;
            }
        };
        let n2 = match Normalizer::new(provider).normalize(&v2) {
            Ok(n) => n,
            Err(e) => {
                failures.push(format!("re-normalize {input} -> {f1}: {e}"));
                continue;
            }
        };
        let f2 = format!("{n2}");
        if f1 != f2 {
            failures.push(format!("not idempotent: {input} -> {f1} -> {f2}"));
        }
        tested += 1;
    }
    eprintln!("normalization_is_idempotent: tested {tested}, skipped {skipped}");
    assert!(
        tested > 0,
        "normalization_is_idempotent exercised no cases (all skipped)"
    );
    assert!(
        failures.is_empty(),
        "Idempotency failed:\n{}",
        failures.join("\n")
    );
}
