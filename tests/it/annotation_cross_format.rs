//! Cross-format consistency invariant: GFF3 and GTF representations of the
//! same biological annotation should produce equivalent Transcript objects.

use ferro_hgvs::reference::annotation::{load_annotations, LoaderConfig};
use ferro_hgvs::reference::transcript::Transcript;
use std::io::Write;
use tempfile::NamedTempFile;

fn assert_transcripts_equivalent(a: &Transcript, b: &Transcript) {
    assert_eq!(a.id, b.id, "id mismatch");
    assert_eq!(a.strand, b.strand, "strand mismatch");
    assert_eq!(a.chromosome, b.chromosome, "chromosome mismatch");
    assert_eq!(a.genomic_start, b.genomic_start, "genomic_start mismatch");
    assert_eq!(a.genomic_end, b.genomic_end, "genomic_end mismatch");
    assert_eq!(a.cds_start, b.cds_start, "cds_start mismatch");
    assert_eq!(a.cds_end, b.cds_end, "cds_end mismatch");
    assert_eq!(a.exons.len(), b.exons.len(), "exon count mismatch");
    for (i, (ea, eb)) in a.exons.iter().zip(&b.exons).enumerate() {
        assert_eq!(
            ea.genomic_start, eb.genomic_start,
            "exon {} genomic_start",
            i
        );
        assert_eq!(ea.genomic_end, eb.genomic_end, "exon {} genomic_end", i);
        assert_eq!(ea.start, eb.start, "exon {} tx start", i);
        assert_eq!(ea.end, eb.end, "exon {} tx end", i);
    }
}

fn write_with_suffix(suffix: &str, content: &str) -> NamedTempFile {
    let mut tf = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
    tf.write_all(content.as_bytes()).unwrap();
    tf.flush().unwrap();
    tf
}

#[test]
fn single_exon_transcript_equivalent_across_formats() {
    let gff = write_with_suffix(
        ".gff3",
        "##gff-version 3\n\
         chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1;Name=GENE1\n\
         chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=g1;gene=GENE1\n\
         chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1\n\
         chr1\t.\tCDS\t150\t450\t.\t+\t0\tParent=tx1\n\
         chr1\t.\tstop_codon\t448\t450\t.\t+\t0\tParent=tx1\n",
    );

    // GTF: same gene, stop included in CDS by convention (so CDS extends to 450).
    let gtf = write_with_suffix(
        ".gtf",
        "chr1\tHAVANA\tgene\t100\t500\t.\t+\t.\tgene_id \"g1\"; gene_name \"GENE1\";\n\
         chr1\tHAVANA\ttranscript\t100\t500\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\"; gene_name \"GENE1\";\n\
         chr1\tHAVANA\texon\t100\t500\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n\
         chr1\tHAVANA\tCDS\t150\t450\t.\t+\t0\tgene_id \"g1\"; transcript_id \"tx1\";\n",
    );

    let (db_gff, _) = load_annotations(gff.path(), &LoaderConfig::new()).unwrap();
    let (db_gtf, _) = load_annotations(gtf.path(), &LoaderConfig::new()).unwrap();

    let tx_gff = db_gff.get("tx1").expect("GFF3 transcript loaded");
    let tx_gtf = db_gtf.get("tx1").expect("GTF transcript loaded");
    assert_transcripts_equivalent(tx_gff, tx_gtf);
}

#[test]
fn multi_exon_transcript_equivalent_across_formats() {
    // Two-exon transcript with an intron; CDS spans both exons.
    let gff = write_with_suffix(
        ".gff3",
        "##gff-version 3\n\
         chr1\t.\tgene\t100\t800\t.\t+\t.\tID=g1\n\
         chr1\t.\tmRNA\t100\t800\t.\t+\t.\tID=tx1;Parent=g1\n\
         chr1\t.\texon\t100\t300\t.\t+\t.\tParent=tx1\n\
         chr1\t.\texon\t500\t800\t.\t+\t.\tParent=tx1\n\
         chr1\t.\tCDS\t150\t300\t.\t+\t0\tParent=tx1\n\
         chr1\t.\tCDS\t500\t650\t.\t+\t0\tParent=tx1\n\
         chr1\t.\tstop_codon\t648\t650\t.\t+\t0\tParent=tx1\n",
    );

    let gtf = write_with_suffix(
        ".gtf",
        "chr1\tHAVANA\ttranscript\t100\t800\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n\
         chr1\tHAVANA\texon\t100\t300\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n\
         chr1\tHAVANA\texon\t500\t800\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n\
         chr1\tHAVANA\tCDS\t150\t300\t.\t+\t0\tgene_id \"g1\"; transcript_id \"tx1\";\n\
         chr1\tHAVANA\tCDS\t500\t650\t.\t+\t0\tgene_id \"g1\"; transcript_id \"tx1\";\n",
    );

    let (db_gff, _) = load_annotations(gff.path(), &LoaderConfig::new()).unwrap();
    let (db_gtf, _) = load_annotations(gtf.path(), &LoaderConfig::new()).unwrap();

    let tx_gff = db_gff.get("tx1").expect("GFF3 transcript loaded");
    let tx_gtf = db_gtf.get("tx1").expect("GTF transcript loaded");
    assert_transcripts_equivalent(tx_gff, tx_gtf);
}
