//! Issue #1139 — a `gene_symbol` on a **non-RefSeq** transcript was rendered as
//! a second, stacked parenthesised specification on the `c.`/`n.` axes:
//! `Template(Template-gene.1)(GENE1):c.5A>G`.
//!
//! Two things are wrong with that string. `background/refseq.md:40-41` permits a
//! **single** specification and restricts its value to a transcript/protein
//! accession ("Gene symbols should not be used as specification"), so both the
//! doubled `(...)(...)` group and the gene symbol in a specification slot are
//! non-conformant. Worse, ferro's own parser rejects it (`Expected ':' after
//! accession`), so `project()` emits a coding name that cannot be fed back to
//! `Normalizer`/`parse`.
//!
//! Root cause: `Display` suppressed the gene-symbol selector only when the
//! accession's prefix was a *recognised* transcript prefix
//! (`Accession::is_transcript_reference` — `NM_`/`NR_`/`ENST`/… ), which a
//! custom id like `Template-gene.1` is not. Once #1086 started stamping a
//! genomic input's own accession as the `genomic_context`, the accession already
//! rendered its one permitted specification (`Template(Template-gene.1)`) and
//! the selector became a second group.
//!
//! PR-CI-runnable: `JsonProvider` fixture over a synthetic ORF, no reference
//! manifest.

use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::{parse_hgvs, JsonProvider, VariantProjector};
use std::io::Write;

/// Synthetic 81 bp ORF (no real locus). Base 5 is `A`, and `c.5A>G` rewrites
/// codon 2 `GAT` (Asp) to `GGT` (Gly).
const ORF: &str =
    "ATGGATTATACTTGTACTTGCAAACAGCCATTTACTAGTGCATGGTTTTTGCTGCAACATGCTCAGAACACTCATGGTTAA";

/// A single custom (non-RefSeq) transcript `Template-gene.1` spanning the whole
/// of contig `Template`, optionally carrying `gene_symbol`.
///
/// `gene_symbol` is what `convert-gff` populates from a GFF3 `gene_name`/`gene`/
/// `Name` attribute on the transcript feature, so this is the shape any
/// custom-reference pipeline produces — not just hand-authored JSON.
fn projector_for(gene_symbol: Option<&str>) -> VariantProjector<JsonProvider> {
    let mut transcript = serde_json::json!({
        "id": "Template-gene.1",
        "chromosome": "Template",
        "strand": "+",
        "sequence": ORF,
        "cds_start": 1,
        "cds_end": 81,
        "genomic_start": 1,
        "genomic_end": 81,
        "protein_id": "Template-gene.1",
        "exons": [{"number": 1, "start": 1, "end": 81, "genomic_start": 1, "genomic_end": 81}],
    });
    if let Some(symbol) = gene_symbol {
        transcript["gene_symbol"] = serde_json::json!(symbol);
    }
    let document = serde_json::json!({
        "version": "1.0",
        "transcripts": [transcript],
        "genomic_sequences": { "Template": ORF },
    });

    let mut file = tempfile::NamedTempFile::new().expect("temp reference file");
    file.write_all(document.to_string().as_bytes())
        .expect("write reference JSON");
    let provider = JsonProvider::from_json(file.path()).expect("load reference JSON");
    let cdot = CdotMapper::from_transcripts(provider.all_transcripts());
    VariantProjector::new(Projector::new(cdot), provider)
}

/// The `c.` axis of `Template:g.5A>G` projected onto `Template-gene.1`.
fn coding_name(gene_symbol: Option<&str>) -> String {
    projector_for(gene_symbol)
        .project("Template:g.5A>G", "Template-gene.1")
        .expect("project g. onto the custom transcript")
        .coding
        .expect("coding axis")
        .to_string()
}

/// The `n.` axis of the same projection.
fn noncoding_name(gene_symbol: Option<&str>) -> String {
    projector_for(gene_symbol)
        .project("Template:g.5A>G", "Template-gene.1")
        .expect("project g. onto the custom transcript")
        .noncoding
        .expect("noncoding axis")
        .to_string()
}

#[test]
fn gene_symbol_is_not_stacked_as_a_second_specification_on_the_coding_axis() {
    // The genomic context IS the one permitted specification; the gene symbol
    // must not be appended as a second group.
    assert_eq!(
        coding_name(Some("GENE1")),
        "Template(Template-gene.1):c.5A>G"
    );
}

#[test]
fn gene_symbol_is_not_stacked_as_a_second_specification_on_the_noncoding_axis() {
    assert_eq!(
        noncoding_name(Some("GENE1")),
        "Template(Template-gene.1):n.5A>G"
    );
}

#[test]
fn projected_coding_name_round_trips_through_the_parser() {
    // The defect surfaced away from its cause: `project()` returned silently and
    // the malformed name only blew up when something re-parsed it.
    let c = coding_name(Some("GENE1"));
    parse_hgvs(&c)
        .unwrap_or_else(|e| panic!("projector emitted an unparseable c. name {c:?}: {e}"));
}

#[test]
fn projected_noncoding_name_round_trips_through_the_parser() {
    let n = noncoding_name(Some("GENE1"));
    parse_hgvs(&n)
        .unwrap_or_else(|e| panic!("projector emitted an unparseable n. name {n:?}: {e}"));
}

#[test]
fn a_reference_without_a_gene_symbol_is_unchanged() {
    // The gene-symbol-free reference always rendered correctly; pin that the fix
    // does not disturb it.
    assert_eq!(coding_name(None), "Template(Template-gene.1):c.5A>G");
    assert_eq!(noncoding_name(None), "Template(Template-gene.1):n.5A>G");
}

#[test]
fn the_protein_axis_stays_gene_symbol_free() {
    // #310 already dropped the symbol from the `p.` axis; the `c.`/`n.` axes now
    // agree with it for the same reference.
    let p = projector_for(Some("GENE1"))
        .project("Template:g.5A>G", "Template-gene.1")
        .expect("project g. onto the custom transcript")
        .protein
        .expect("protein axis")
        .to_string();
    assert_eq!(p, "Template-gene.1:p.(Asp2Gly)");
}
