//! End-to-end projector coverage matrix.
//!
//! Closes the gap identified after c14:
//!
//! 1. Pre-c16 the projector had exactly ONE minus-strand end-to-end test
//!    (a single-base substitution on a single-exon transcript). The
//!    component-level tests (`src/project/edit.rs`, `src/data/cdot.rs`)
//!    exercise pieces of the pipeline but don't compose them through
//!    `VariantProjector::project`. The cases here run the full pipeline
//!    for both strands and every NaEdit category.
//!
//! 2. Pre-c17 the indel-protein tests covered each edit *type* but not
//!    each edit *position* (start of CDS, end of CDS, exon junction,
//!    UTR-adjacent, etc.). The cases here pin those structural edges.
//!
//! Fixture conventions: minimal hand-built `CdotTranscript` +
//! `MockProvider` pairs so the test reads the same way the existing
//! `tests/projection.rs` fixtures do. No external manifest required.

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::sequence::reverse_complement;
use ferro_hgvs::VariantProjector;

// ─── shared fixture builders ─────────────────────────────────────────────────

/// Plus-strand single-exon transcript on chr1:
///   genome 1000..1009 (1-based 1000-1008 inclusive)
///   CDS: full 9 bases — Met-Arg-Stop (ATG CGC TAA)
fn plus_single_exon() -> (Projector, MockProvider) {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_PLUS1.1".to_string(),
        CdotTranscript {
            gene_name: Some("PLUSGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            exons: vec![[1000, 1009, 0, 9]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_PLUS1.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_PLUS1.1".to_string(),
        Some("PLUSGENE".to_string()),
        TxStrand::Plus,
        "ATGCGCTAA".to_string(),
        Some(1),
        Some(9),
        vec![Exon::new(1, 1, 9)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1008),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    let prefix = "N".repeat(999);
    let suffix = "N".repeat(100);
    provider.add_genomic_sequence("chr1", format!("{}ATGCGCTAA{}", prefix, suffix));
    (projector, provider)
}

/// Minus-strand single-exon transcript on chr1:
///   genome 1000..1009 same coordinates, but the transcript's 5'→3'
///   direction is the reverse-complement.
///   Transcript reads "ATGCGCTAA" (same Met-Arg-Stop) — that means the
///   genomic + strand at this range carries the reverse-complement of
///   the transcript: "TTAGCGCAT".
fn minus_single_exon() -> (Projector, MockProvider) {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_MINUS1.1".to_string(),
        CdotTranscript {
            gene_name: Some("MINUSGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Minus,
            exons: vec![[1000, 1009, 0, 9]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_MINUS1.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_MINUS1.1".to_string(),
        Some("MINUSGENE".to_string()),
        TxStrand::Minus,
        "ATGCGCTAA".to_string(),
        Some(1),
        Some(9),
        vec![Exon::new(1, 1, 9)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1008),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    // Genomic + strand carries revcomp("ATGCGCTAA") = "TTAGCGCAT" at
    // 1-based positions 1000-1008.
    let prefix = "N".repeat(999);
    let suffix = "N".repeat(100);
    let genomic = format!("{}{}{}", prefix, reverse_complement("ATGCGCTAA"), suffix);
    provider.add_genomic_sequence("chr1", genomic);
    (projector, provider)
}

/// Plus-strand two-exon transcript on chr2:
///   Exon 1: genome 2000..2006, tx 0..6   ("ATGCGC")  — Met + first codon Arg
///   Intron : genome 2006..2100
///   Exon 2: genome 2100..2103, tx 6..9   ("TAA")     — Stop
///   CDS spans the whole transcript.
fn plus_two_exon() -> (Projector, MockProvider) {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_PLUS2.1".to_string(),
        CdotTranscript {
            gene_name: Some("PLUS2GENE".to_string()),
            contig: "chr2".to_string(),
            strand: Strand::Plus,
            exons: vec![[2000, 2006, 0, 6], [2100, 2103, 6, 9]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_PLUS2.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_PLUS2.1".to_string(),
        Some("PLUS2GENE".to_string()),
        TxStrand::Plus,
        "ATGCGCTAA".to_string(),
        Some(1),
        Some(9),
        vec![
            Exon::with_genomic(1, 1, 6, 2000, 2005),
            Exon::with_genomic(2, 7, 9, 2100, 2102),
        ],
        Some("chr2".to_string()),
        Some(2000),
        Some(2102),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    // Build genome: exon1 + intron + exon2.
    let mut genomic = "N".repeat(1999);
    genomic.push_str("ATGCGC"); // exon 1
    genomic.push_str(&"N".repeat(94)); // intron 2006..2100
    genomic.push_str("TAA"); // exon 2
    genomic.push_str(&"N".repeat(100));
    provider.add_genomic_sequence("chr2", genomic);
    (projector, provider)
}

// ─── minus-strand × full edit matrix ─────────────────────────────────────────

/// Minus-strand substitution: c.4C>A.
///
/// Transcript reads "ATG CGC TAA". c.4 is the first base of codon 2 (Arg
/// codon CGC), so c.4C>A changes CGC → AGC = Ser. On the genome (+ strand)
/// this is at the position carrying the revcomp of c.4 = revcomp(C) = G,
/// changing to revcomp(A) = T.
///
///   genome 1004 carries reverse-complement of c.4 (1004 = 1008-4):
///   1008 ← c.1 / 1007 ← c.2 / 1006 ← c.3 / 1005 ← c.4 / …
///   Wait — for a 9-base exon on `-` strand the mapping is
///     c.i = (genome_end - 1) - g_pos_0based = 1008 - g_pos
///   so c.4 ↔ g.1005 (genome 1-based).
#[test]
fn minus_strand_substitution() {
    let (projector, provider) = minus_single_exon();
    let vp = VariantProjector::new(projector, provider);
    let res = vp.project("chr1:g.1005G>T", "NM_MINUS1.1").unwrap();
    let c = res.coding.unwrap().to_string();
    assert!(c.contains(":c.4C>A"), "got c. = {c}");
    let p = res.protein.unwrap().to_string();
    assert!(p.contains("Arg2Ser"), "got p. = {p}");
    assert!(!res.is_frameshift && !res.is_intronic && !res.is_utr);
}

/// Minus-strand single-base deletion at c.4 — frameshift.
#[test]
fn minus_strand_deletion_frameshift() {
    let (projector, provider) = minus_single_exon();
    let vp = VariantProjector::new(projector, provider);
    let res = vp.project("chr1:g.1005del", "NM_MINUS1.1").unwrap();
    let c = res.coding.unwrap().to_string();
    assert!(c.contains("del"), "got c. = {c}");
    assert!(res.is_frameshift);
}

/// Minus-strand codon-aligned deletion (whole codon CGC at c.4_6).
#[test]
fn minus_strand_codon_aligned_deletion() {
    let (projector, provider) = minus_single_exon();
    let vp = VariantProjector::new(projector, provider);
    // Delete c.4_6 = CGC. On the genome (+ strand) that's bases carrying
    // revcomp(CGC) = GCG, which are at g.1003_1005.
    let res = vp.project("chr1:g.1003_1005del", "NM_MINUS1.1").unwrap();
    let c = res.coding.unwrap().to_string();
    assert!(c.contains(":c.4_6del"), "got c. = {c}");
    let p = res.protein.unwrap().to_string();
    assert!(p.contains("Arg2del"), "got p. = {p}");
    assert!(!res.is_frameshift);
}

/// Minus-strand insertion of three bases between c.3 and c.4 (codon-aligned).
/// Inserted "GGG" on the transcript = "CCC" on the genome.
///
/// The HGVS canonical form may render this as a `repeat` (`c.<pos><base>[<n>]`)
/// rather than `ins` when the inserted bases adjoin a homopolymer run — both
/// representations are valid. We assert structural properties (success +
/// in-frame + non-empty c.) rather than the exact lexical form.
#[test]
fn minus_strand_codon_aligned_insertion() {
    let (projector, provider) = minus_single_exon();
    let vp = VariantProjector::new(projector, provider);
    let res = vp.project("chr1:g.1005_1006insCCC", "NM_MINUS1.1").unwrap();
    let c = res.coding.unwrap().to_string();
    assert!(c.starts_with("NM_MINUS1.1"), "got c. = {c}");
    assert!(c.contains(":c."), "got c. = {c}");
    assert!(!res.is_frameshift);
}

/// Minus-strand single-base insertion — frameshift.
#[test]
fn minus_strand_insertion_frameshift() {
    let (projector, provider) = minus_single_exon();
    let vp = VariantProjector::new(projector, provider);
    let res = vp.project("chr1:g.1005_1006insT", "NM_MINUS1.1").unwrap();
    let c = res.coding.unwrap().to_string();
    assert!(c.contains("ins"), "got c. = {c}");
    assert!(res.is_frameshift);
}

/// Minus-strand duplication of c.4_6 (whole Arg codon).
#[test]
fn minus_strand_codon_aligned_duplication() {
    let (projector, provider) = minus_single_exon();
    let vp = VariantProjector::new(projector, provider);
    let res = vp.project("chr1:g.1003_1005dup", "NM_MINUS1.1").unwrap();
    let c = res.coding.unwrap().to_string();
    assert!(c.contains("dup"), "got c. = {c}");
    let p = res.protein.unwrap().to_string();
    // Duplicating Arg → inserts an extra Arg after position 2.
    assert!(p.contains("Arg") && p.contains("dup"), "got p. = {p}");
}

/// Minus-strand delins replacing the Arg codon with two codons.
#[test]
fn minus_strand_delins() {
    let (projector, provider) = minus_single_exon();
    let vp = VariantProjector::new(projector, provider);
    // c.4_6delinsAAATTT (replace CGC with AAATTT = Lys+Phe, in-frame +3 bases).
    // On the genome (+ strand) the replacement bases are revcomp(AAATTT) =
    // "AAATTT" → "AAATTT" reversed = "TTTAAA"? Actually
    // revcomp("AAATTT") = "AAATTT" (palindrome-like). Compute:
    //   AAATTT → complement TTTAAA → reverse AAATTT. Yes palindromic.
    let res = vp
        .project("chr1:g.1003_1005delinsAAATTT", "NM_MINUS1.1")
        .unwrap();
    let c = res.coding.unwrap().to_string();
    assert!(c.contains("delins"), "got c. = {c}");
    assert!(!res.is_frameshift);
}

/// Minus-strand inversion of the Arg codon CGC → GCG (in-frame, same length).
#[test]
fn minus_strand_inversion() {
    let (projector, provider) = minus_single_exon();
    let vp = VariantProjector::new(projector, provider);
    let res = vp.project("chr1:g.1003_1005inv", "NM_MINUS1.1").unwrap();
    let c = res.coding.unwrap().to_string();
    assert!(c.contains("inv"), "got c. = {c}");
    assert!(!res.is_frameshift);
}

// ─── multi-exon spans (plus strand) ──────────────────────────────────────────

/// Substitution at the last base of exon 1 on a two-exon plus-strand
/// transcript. This exercises the per-exon mapping logic at an exon
/// boundary.
#[test]
fn plus_two_exon_substitution_in_exon1() {
    let (projector, provider) = plus_two_exon();
    let vp = VariantProjector::new(projector, provider);
    // g.2005 is the last base of exon 1 → c.6 (the 'C' of codon 2 CGC).
    let res = vp.project("chr2:g.2005C>A", "NM_PLUS2.1").unwrap();
    let c = res.coding.unwrap().to_string();
    assert!(c.contains(":c.6C>A"), "got c. = {c}");
}

/// Substitution at the first base of exon 2.
#[test]
fn plus_two_exon_substitution_in_exon2() {
    let (projector, provider) = plus_two_exon();
    let vp = VariantProjector::new(projector, provider);
    // g.2100 is the first base of exon 2 → c.7 (the 'T' of the stop codon).
    let res = vp.project("chr2:g.2100T>A", "NM_PLUS2.1").unwrap();
    let c = res.coding.unwrap().to_string();
    assert!(c.contains(":c.7T>A"), "got c. = {c}");
}

/// Intronic substitution should project to c.N+offset notation and flag
/// intronic.
#[test]
fn plus_two_exon_intronic_substitution() {
    let (projector, provider) = plus_two_exon();
    let vp = VariantProjector::new(projector, provider);
    // g.2050 is well inside the intron between exon 1 (ends at 2005) and
    // exon 2 (starts at 2100).
    let res = vp.project("chr2:g.2050N>A", "NM_PLUS2.1").unwrap();
    let c = res.coding.unwrap().to_string();
    // Intronic offset notation: c.6+44 (44 bp past the last base of exon 1).
    assert!(
        c.contains('+') || c.contains('-'),
        "expected intronic notation, got {c}"
    );
    assert!(res.is_intronic);
}

// ─── indel structural edges (plus-strand, simple exon) ───────────────────────

/// Edit at the very first base of CDS (c.1): substitution of the start
/// codon's first base.
#[test]
fn substitution_at_c1_first_cds_base() {
    let (projector, provider) = plus_single_exon();
    let vp = VariantProjector::new(projector, provider);
    let res = vp.project("chr1:g.1000A>G", "NM_PLUS1.1").unwrap();
    let c = res.coding.unwrap().to_string();
    assert!(c.contains(":c.1A>G"), "got c. = {c}");
    // c.1A>G changes ATG (Met) → GTG (Val). HGVS p. recommendation: report
    // as p.(Met1?) since we can't tell if translation still starts there.
    // We accept either the strict ? or the explicit Met1Val depending on
    // implementation; just assert SOMETHING with Met1 was emitted.
    let p = res.protein.unwrap().to_string();
    assert!(p.contains("Met1") || p.contains("M1"), "got p. = {p}");
}

/// Edit at the last base before the stop codon — codon 2 last base
/// substitution.
#[test]
fn substitution_at_last_cds_base_before_stop() {
    let (projector, provider) = plus_single_exon();
    let vp = VariantProjector::new(projector, provider);
    // c.6 is the last base of Arg codon (CGC). Change C→A → CGA (still Arg).
    let res = vp.project("chr1:g.1005C>A", "NM_PLUS1.1").unwrap();
    let c = res.coding.unwrap().to_string();
    assert!(c.contains(":c.6C>A"), "got c. = {c}");
    let p = res.protein.unwrap().to_string();
    // Synonymous substitution: HGVS uses p.(Arg2=) or similar identity form.
    assert!(
        p.contains('='),
        "expected synonymous identity, got p. = {p}"
    );
}

/// Indel that spans an exon-intron junction — the c./n. position lands at
/// the splice site. We just assert the projection succeeds and the result
/// reflects the junction (intronic-flagged or splice-distance reported).
#[test]
fn indel_spanning_exon_intron_junction() {
    let (projector, provider) = plus_two_exon();
    let vp = VariantProjector::new(projector, provider);
    // Delete g.2005_2006 — last base of exon 1 + first base of intron.
    let res = vp.project("chr2:g.2005_2006del", "NM_PLUS2.1");
    // Some projectors will treat this as intronic with `del` notation;
    // others as a splice-affecting deletion. Accept either path as long
    // as it didn't error.
    assert!(
        res.is_ok(),
        "exon/intron-junction del must project, got {:?}",
        res.err()
    );
    let r = res.unwrap();
    let c = r.coding.unwrap().to_string();
    assert!(c.contains("del"), "got c. = {c}");
}
