//! Tests against hand-curated real-format slice fixtures. These fixtures
//! exercise format-specific quirks (RefSeq attribute conventions, GENCODE
//! tag fields, FlyBase Derives_from) that purely programmatic fixtures miss.

use ferro_hgvs::reference::annotation::{load_annotations, LoaderConfig};
use std::path::Path;

#[test]
fn refseq_micro_loads_minus_strand_transcript_with_mane_select() {
    let (db, report) = load_annotations(
        Path::new("tests/fixtures/annotation/refseq_chr21_micro.gff3"),
        &LoaderConfig::new(),
    )
    .unwrap();
    assert_eq!(db.len(), 1, "expected one transcript from RefSeq fixture");
    let tx = db.get("NR_027228.1").expect("RefSeq transcript loaded");
    assert_eq!(tx.exons.len(), 2);
    assert!(matches!(
        tx.strand,
        ferro_hgvs::reference::transcript::Strand::Minus
    ));
    assert!(matches!(
        tx.mane_status,
        ferro_hgvs::reference::transcript::ManeStatus::Select
    ));
    assert_eq!(report.transcripts_loaded, 1);
}

#[test]
fn gencode_micro_loads_minus_strand_transcript() {
    let (db, _report) = load_annotations(
        Path::new("tests/fixtures/annotation/gencode_chr21_micro.gtf"),
        &LoaderConfig::new(),
    )
    .unwrap();
    assert_eq!(db.len(), 1);
    let tx = db
        .get("ENST00000620911")
        .expect("GENCODE transcript loaded");
    assert_eq!(tx.exons.len(), 2);
    assert!(matches!(
        tx.strand,
        ferro_hgvs::reference::transcript::Strand::Minus
    ));
    assert_eq!(tx.gene_symbol.as_deref(), Some("MIR99AHG"));
    assert!(matches!(
        tx.mane_status,
        ferro_hgvs::reference::transcript::ManeStatus::Select
    ));
}

#[test]
fn flybase_micro_loads_two_exon_protein_coding() {
    let (db, _report) = load_annotations(
        Path::new("tests/fixtures/annotation/flybase_micro.gff3"),
        &LoaderConfig::new(),
    )
    .unwrap();
    assert_eq!(db.len(), 1);
    let tx = db.get("FBtr0300689").expect("FlyBase transcript loaded");
    assert_eq!(tx.exons.len(), 2);
    // CDS spans both exons; cds_start should fall within exon 1 and cds_end within exon 2
    // (both in transcript coordinates).
    let cds_start = tx.cds_start.expect("cds_start present");
    let cds_end = tx.cds_end.expect("cds_end present");
    assert!(
        cds_start >= tx.exons[0].start && cds_start <= tx.exons[0].end,
        "cds_start {} not in exon[0] [{}, {}]",
        cds_start,
        tx.exons[0].start,
        tx.exons[0].end
    );
    assert!(
        cds_end >= tx.exons[1].start && cds_end <= tx.exons[1].end,
        "cds_end {} not in exon[1] [{}, {}]",
        cds_end,
        tx.exons[1].start,
        tx.exons[1].end
    );
}
