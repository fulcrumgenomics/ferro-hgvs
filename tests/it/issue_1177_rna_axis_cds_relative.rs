//! Issue #1177 — an `r.` input's position must be read as **CDS-relative**.
//!
//! On a coding transcript the `r.` axis is CDS-relative — the same axis as
//! `c.`, so `r.N` describes the same base as `c.N` (HGVS
//! `background/numbering.md` L58/L61; pinned for the normalizer by #469 and
//! `issue_291_rna_axis_convention.rs`). `VariantProjection::rna` already
//! documents its own output that way ("CDS-relative numbering (matches `c.`)").
//!
//! The projector's *input* handling disagreed: it copied `RnaPos.base` straight
//! into a transcript-relative position, on the premise that the two share a
//! `base`/`offset`/`utr3` struct shape. They share a shape but not an *origin*,
//! so every `r.` input on a coding transcript was silently mis-projected by
//! `cds_start` — no error, no warning. On the genomic axis of a real transcript
//! with large introns that put the result hundreds of kb from the truth.
//!
//! The invariant these tests pin is **cross-axis agreement**: three spellings of
//! the same physical base — `c.N`, `n.(cds_start + N - 1)` and `r.N` — must
//! project to identical results on every axis. That invariant was entirely
//! absent, which is why the defect survived; it is not expressible as a
//! single-axis expectation, and the `FERRO_ASSERT_IDEMPOTENT` oracle cannot see
//! it either (it checks `norm(norm(x)) == norm(x)`, not cross-axis agreement).

use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::{parse_hgvs, VariantProjection, VariantProjector};

/// 5'UTR length, and therefore the 1-based transcript position of `c.1`.
///
/// Deliberately non-trivial: with `cds_start = 1` the CDS-relative and
/// transcript-relative axes coincide and the bug is invisible.
const CDS_START: i64 = 21;

/// 20-base 5'UTR, then a 30-base CDS (`ATG` + 8 codons + `TAA`), then a 10-base
/// 3'UTR — 60 bases total.
///
/// | tx    | region | notes                          |
/// |-------|--------|--------------------------------|
/// | 1..20 | 5'UTR  | `tx 11` = `G` = `c.-10` = `r.-10` |
/// | 21..50| CDS    | `tx 21` = `A` of `ATG` = `c.1`; `tx 25` = `G` = `c.5` = `r.5` |
/// | 51..60| 3'UTR  |                                |
const TX_SEQUENCE: &str = concat!(
    "CTGACTGACTGACTGACTGA",     // 1..20   5'UTR
    "ATG",                      // 21..23  c.1..c.3
    "CGCGCGCGCGCGCGCGCGCGCGCG", // 24..47  c.4..c.27
    "TAA",                      // 48..50  c.28..c.30 (stop)
    "GTCAGTCAGT",               // 51..60  3'UTR
);

fn projector() -> VariantProjector<MockProvider> {
    let mut provider = MockProvider::new();
    let tx = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("MYGENE".to_string()),
        TxStrand::Plus,
        TX_SEQUENCE.to_string(),
        Some(CDS_START as u64),
        Some(50),
        vec![Exon::with_genomic(1, 1, 60, 1000, 1059)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1059),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    )
    .with_protein_id(Some("NP_TEST.1".to_string()));
    provider.add_transcript(tx);
    let prefix = "N".repeat(999);
    let suffix = "N".repeat(100);
    provider.add_genomic_sequence("chr1", format!("{prefix}{TX_SEQUENCE}{suffix}"));
    let cdot = CdotMapper::from_transcripts(provider.all_transcripts());
    VariantProjector::new(Projector::new(cdot), provider)
}

/// Every axis of a projection, rendered for comparison.
fn axes(p: &VariantProjection) -> Vec<(&'static str, Option<String>)> {
    vec![
        ("g", p.genomic.as_ref().map(ToString::to_string)),
        ("c", p.coding.as_ref().map(ToString::to_string)),
        ("n", p.noncoding.as_ref().map(ToString::to_string)),
        ("r", p.rna.as_ref().map(ToString::to_string)),
        ("p", p.protein.as_ref().map(ToString::to_string)),
    ]
}

fn project(input: &str) -> VariantProjection {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    projector()
        .project_variant(&variant, "NM_TEST.1")
        .unwrap_or_else(|e| panic!("project {input} must not error: {e}"))
}

/// Assert that three spellings of one base agree on every axis.
///
/// `c.` is the reference spelling: it is the axis the projector has always read
/// correctly, so a disagreement localises the fault to the other input.
fn assert_axes_agree(coding: &str, noncoding: &str, rna: &str) {
    let expected = axes(&project(coding));
    for other in [noncoding, rna] {
        assert_eq!(
            axes(&project(other)),
            expected,
            "`{other}` and `{coding}` describe the same base, so every projected axis must \
             match; a mismatch means the input's coordinate axis was misread",
        );
    }
}

/// A positive, in-CDS `r.` coordinate. `r.5` is `c.5` is `n.25`
/// (`CDS_START + 5 - 1`).
///
/// Pre-fix this reported `c.-16` for the `r.` input (`5 - 21`), i.e. an
/// ordinary in-CDS RNA description was relocated into the 5'UTR.
#[test]
fn in_cds_r_position_agrees_with_c_and_n() {
    assert_axes_agree("NM_TEST.1:c.5G>A", "NM_TEST.1:n.25G>A", "NM_TEST.1:r.5g>a");
}

/// A negative 5'UTR `r.` coordinate. `r.-10` is `c.-10` is `n.11`
/// (`CDS_START - 10`) — note HGVS `c.`/`r.` numbering has no zero, so the
/// mapping is `cds_start + N` for negative `N`, not `cds_start + N - 1`.
#[test]
fn five_prime_utr_r_position_agrees_with_c_and_n() {
    assert_axes_agree(
        "NM_TEST.1:c.-10G>A",
        "NM_TEST.1:n.11G>A",
        "NM_TEST.1:r.-10g>a",
    );
}

/// The same invariant on the **genomic** axis, which a bare `NM_` input cannot
/// reach (no genome alignment, so `genomic` is `None`). Supplying an explicit
/// genomic context routes the projection through the genome pivot
/// (`map_cnr_position_to_genome`) — a second, independent place where the input
/// axis is interpreted.
///
/// This matters most in practice: on a real transcript with large introns a
/// `cds_start`-sized error in transcript space lands the genomic coordinate
/// hundreds of kb away from the truth.
#[test]
fn genomic_axis_of_r_position_agrees_with_c_and_n() {
    assert_axes_agree(
        "chr1(NM_TEST.1):c.5G>A",
        "chr1(NM_TEST.1):n.25G>A",
        "chr1(NM_TEST.1):r.5g>a",
    );
}

/// The `r.` axis of an `r.` input must be a fixed point: re-projecting an
/// already-`r.` description must not move it.
///
/// This is the property whose failure originally exposed the bug — re-feeding
/// the `r.` output walked by exactly `-cds_start` per pass, forever.
#[test]
fn r_axis_projection_is_a_fixed_point() {
    let first = project("NM_TEST.1:r.5g>a")
        .rna
        .expect("an r. input must have an r. axis")
        .to_string();
    // Pin the position, not just the stability. A fixed-point assertion alone is
    // satisfied by any bug that shifts once and then stays put, which is exactly
    // the shape of a constant `cds_start` offset applied at a single conversion
    // site. The position must come back as `5` — the r. axis of an r. input is
    // the same base — while the `( )` is expected: the r. axis is *predicted*
    // from the c. form (`predict_rna`), and HGVS parenthesises a predicted RNA
    // consequence.
    assert_eq!(
        first, "NM_TEST.1:r.(5g>a)",
        "the r. axis of an r. input must reproduce the input's position, as a predicted form",
    );
    let second = project(&first)
        .rna
        .expect("an r. input must have an r. axis")
        .to_string();
    assert_eq!(
        second, first,
        "projecting an r. description onto the r. axis must be idempotent; drift means the \
         input's coordinate axis is re-interpreted on each pass",
    );
}

/// Guardrail: on a **non-coding** transcript `r.` and `n.` genuinely are the
/// same axis (there is no CDS to be relative to), so the verbatim copy is
/// correct there and must stay untouched by the fix.
#[test]
fn noncoding_transcript_r_position_is_transcript_relative() {
    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NR_TEST.1".to_string(),
        Some("MYGENE".to_string()),
        TxStrand::Plus,
        TX_SEQUENCE.to_string(),
        None,
        None,
        vec![Exon::with_genomic(1, 1, 60, 1000, 1059)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1059),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    let prefix = "N".repeat(999);
    provider.add_genomic_sequence("chr1", format!("{prefix}{TX_SEQUENCE}{}", "N".repeat(100)));
    let cdot = CdotMapper::from_transcripts(provider.all_transcripts());
    let vp = VariantProjector::new(Projector::new(cdot), provider);

    let rna = parse_hgvs("NR_TEST.1:r.25g>a").expect("parse");
    let tx = parse_hgvs("NR_TEST.1:n.25G>A").expect("parse");
    let from_rna = vp.project_variant(&rna, "NR_TEST.1").expect("project r.");
    let from_tx = vp.project_variant(&tx, "NR_TEST.1").expect("project n.");
    assert_eq!(
        axes(&from_rna),
        axes(&from_tx),
        "with no CDS, r.25 and n.25 are the same position and must project identically",
    );
}
