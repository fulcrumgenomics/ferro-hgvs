//! #879 — a bare `NG_`/`LRG_` genomic input projected onto a *single* requested
//! transcript must de-anchor into the `NC_` chromosome frame before the overlap
//! guard runs, or it declines with `TranscriptNotOverlapping` (the parent-frame
//! `g.` coordinate is compared against `NC_`-frame cdot exon spans).
//!
//! Split from #857. The fan-out (`project_variant_all`) already de-anchors +
//! re-frames; this locks in that the single-transcript paths (`project_variant`,
//! `project_normalized`) do too — while leaving `NG_(NM_):c.` transcript inputs
//! on their existing pivot path (the Genome-gate, review M1).
//!
//! CI-runnable (MockProvider, no manifest). The `GenomicPlacement` is
//! **non-identity** (`parent_start != nc_start`) so a frame bug cannot hide
//! behind an identity transform (cf. the `issue_480_*` pattern), and the
//! transcript has a terminal 3'UTR (`cds_end < exon_end`) so a boundary deletion
//! spans `c.<n>_*1` and exercises the #857 C-terminal-extension protein shape.

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::provider::GenomicPlacement;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{parse_hgvs, FerroError, VariantProjection, VariantProjector};

/// CDS "ATGAAATAA" (Met-Lys-Ter, stop at c.7_9) followed by the 3'UTR "TCTAA",
/// so the full transcript is "ATGAAATAATCTAA" (14 bp). The CDS ends at `c.9`
/// (`cds_end = 9`) while the exon runs to tx position 14, so `c.9_*1` straddles
/// the CDS→3'UTR boundary. Placed on chr1 plus strand at HGVS `[1000, 1014)`.
const FULL_TX: &str = "ATGAAATAATCTAA";

/// Build a projector + provider for `NM_TEST.1`/`NP_TEST.1` on chr1 (plus
/// strand). Placements are added per-test so each can pick a strand/offset.
fn base_fixture() -> (Projector, MockProvider) {
    let len = FULL_TX.len() as u64;
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TEST.1".to_string(),
        CdotTranscript {
            gene_name: Some("TESTGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            // [g_start, g_end, tx_start, tx_end]; exon spans HGVS [1000, 1014).
            exons: vec![[1000, 1000 + len, 0, len]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_TEST.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TESTGENE".to_string()),
        TxStrand::Plus,
        FULL_TX.to_string(),
        Some(1),
        Some(9),
        vec![Exon::new(1, 1, len)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1000 + len - 1),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    // 999-N pad so HGVS exon coord 1000 lands at 0-based index 999; the exon's
    // genome bases then match the transcript and the #644 sequence-aware
    // projection sees no phantom indel (see the note in issue_480's base_fixture).
    let prefix = "N".repeat(999);
    let suffix = "N".repeat(100);
    provider.add_genomic_sequence("chr1", format!("{prefix}{FULL_TX}{suffix}"));
    (projector, provider)
}

fn parse_acc(s: &str) -> ferro_hgvs::hgvs::variant::Accession {
    match parse_hgvs(&format!("{s}:g.1=")).expect("parse accession") {
        ferro_hgvs::HgvsVariant::Genome(g) => g.accession,
        other => panic!("expected genome variant, got {other:?}"),
    }
}

/// Plus-strand, non-identity placement: `NG_TEST.1` base 5001 == chr 1000, so
/// the CDS→3'UTR boundary deletion `c.9_*1del` (chr 1008_1009) is written in the
/// parent frame as `NG_TEST.1:g.5009_5010del`.
///
/// The expected `c./p./g.` strings are a deliberate regression pin from the
/// fixture's known placement + CDS layout, cross-validated against the real
/// manifest reproducer in #879 (`NG_012337.1:g.13124_13125del` →
/// `c.480_*1del` / `p.(Ter160CysextTer29)`), which exercises the same code path.
#[test]
fn bare_ng_genomic_projects_onto_single_transcript_plus_strand() {
    let (projector, mut provider) = base_fixture();
    provider.add_genomic_placement(
        "NG_TEST.1",
        GenomicPlacement {
            nc: parse_acc("NC_000001.11"),
            parent_start: 5001,
            nc_start: 1000,
            nc_end: 1013,
            strand: Strand::Plus,
        },
    );
    let vp = VariantProjector::new(projector, provider);

    // Parent-frame input straddling the CDS→3'UTR boundary.
    let variant = parse_hgvs("NG_TEST.1:g.5009_5010del").expect("parse");
    let r = vp
        .project_variant(&variant, "NM_TEST.1")
        .expect("bare NG_ genomic input must de-anchor and project, not decline");

    assert_eq!(
        r.coding.as_ref().map(ToString::to_string).as_deref(),
        Some("NG_TEST.1(NM_TEST.1):c.9_*1del"),
        "coding must be framed under the NG_ parent"
    );
    // The protein axis (the #857 C-terminal-extension shape) must survive the
    // re-frame and be relabeled under the NG_ parent (review N1).
    assert_eq!(
        r.protein.as_ref().map(ToString::to_string).as_deref(),
        Some("NG_TEST.1(NP_TEST.1):p.(Ter3TyrextTer1)"),
        "protein axis must survive + be reframed under the NG_ parent"
    );
    // The genomic axis is re-anchored back into the NG_ parent frame.
    assert_eq!(
        r.genomic.as_ref().map(ToString::to_string).as_deref(),
        Some("NG_TEST.1:g.5009_5010del"),
        "genomic axis must be re-anchored into the NG_ parent frame"
    );
}

/// Minus-strand placement (revcomp de-anchor path): `NG_REV.1` runs antiparallel
/// to the chromosome, so a parent-frame substitution reverse-complements into
/// the `NC_` frame before projection, then re-frames back.
#[test]
fn bare_ng_genomic_projects_onto_single_transcript_minus_strand() {
    let (projector, mut provider) = base_fixture();
    provider.add_genomic_placement(
        "NG_REV.1",
        GenomicPlacement {
            nc: parse_acc("NC_000001.11"),
            parent_start: 1,
            nc_start: 1000,
            nc_end: 1013,
            strand: Strand::Minus,
        },
    );
    let vp = VariantProjector::new(projector, provider);

    // chr 1003 (= c.4, base A) maps to parent frame 11 (1 + 1013 - 1003); the
    // parent reads the complement (T), so `T>G` here revcomps to `A>C` on the
    // chromosome → c.4A>C → Lys2Gln.
    let variant = parse_hgvs("NG_REV.1:g.11T>G").expect("parse");
    let r = vp
        .project_variant(&variant, "NM_TEST.1")
        .expect("minus-strand bare NG_ genomic input must de-anchor and project");

    assert_eq!(
        r.coding.as_ref().map(ToString::to_string).as_deref(),
        Some("NG_REV.1(NM_TEST.1):c.4A>C"),
        "coding must be framed under the NG_ parent (minus strand)"
    );
    assert_eq!(
        r.protein.as_ref().map(ToString::to_string).as_deref(),
        Some("NG_REV.1(NP_TEST.1):p.(Lys2Gln)"),
        "protein axis must survive + be reframed under the NG_ parent (minus strand)"
    );
    assert_eq!(
        r.genomic.as_ref().map(ToString::to_string).as_deref(),
        Some("NG_REV.1:g.11T>G"),
        "genomic axis must be re-anchored into the NG_ parent frame (minus strand)"
    );
}

/// `project_normalized` on a (pre-normalized) bare-`NG_` genomic input must
/// yield the same result as `project_variant` — the parent branch re-normalizes
/// in the `NC_` frame (review M2), a documented exception to its skip-normalize
/// contract, since the caller cannot pre-normalize in a frame it doesn't know.
#[test]
fn project_normalized_de_anchors_bare_ng_genomic() {
    let (projector, mut provider) = base_fixture();
    provider.add_genomic_placement(
        "NG_TEST.1",
        GenomicPlacement {
            nc: parse_acc("NC_000001.11"),
            parent_start: 5001,
            nc_start: 1000,
            nc_end: 1013,
            strand: Strand::Plus,
        },
    );
    let vp = VariantProjector::new(projector, provider);

    let variant = parse_hgvs("NG_TEST.1:g.5009_5010del").expect("parse");
    let via_normalized = vp
        .project_normalized(&variant, "NM_TEST.1")
        .expect("project_normalized must de-anchor a bare NG_ genomic input too");
    // Cross-path check (not a duplicated string pin): `project_normalized` (whose
    // parent branch re-normalizes) must agree with `project_variant` on every axis
    // for the same bare-`NG_` genomic input — the two are distinct entry points.
    let via_variant = vp
        .project_variant(&variant, "NM_TEST.1")
        .expect("project_variant must de-anchor a bare NG_ genomic input");
    let axes = |r: &VariantProjection| {
        (
            r.coding.as_ref().map(ToString::to_string),
            r.protein.as_ref().map(ToString::to_string),
            r.genomic.as_ref().map(ToString::to_string),
            r.noncoding.as_ref().map(ToString::to_string),
            r.rna.as_ref().map(ToString::to_string),
        )
    };
    assert_eq!(
        axes(&via_normalized),
        axes(&via_variant),
        "project_normalized must match project_variant on a bare NG_ genomic input",
    );
    // Pin the shared value so a change to BOTH paths still surfaces.
    assert_eq!(
        via_normalized
            .coding
            .as_ref()
            .map(ToString::to_string)
            .as_deref(),
        Some("NG_TEST.1(NM_TEST.1):c.9_*1del"),
    );
}

/// Lock-in (review M1): a `c.` input carrying an `NG_` `genomic_context` must be
/// **unchanged** by this PR — the Genome-gate keeps it on its existing pivot
/// path (it must NOT be routed through the de-anchor helper, whose cnr branch
/// would otherwise fire). Pin all three axes so a regression is caught.
#[test]
fn ng_context_coding_input_keeps_pivot_path_unchanged() {
    let (projector, mut provider) = base_fixture();
    provider.add_genomic_placement(
        "NG_TEST.1",
        GenomicPlacement {
            nc: parse_acc("NC_000001.11"),
            parent_start: 5001,
            nc_start: 1000,
            nc_end: 1013,
            strand: Strand::Plus,
        },
    );
    let vp = VariantProjector::new(projector, provider);

    // c. input with NG_ context — NOT a bare Genome input.
    let variant = parse_hgvs("NG_TEST.1(NM_TEST.1):c.4A>C").expect("parse");
    let r = vp
        .project_variant(&variant, "NM_TEST.1")
        .expect("c. input with NG_ context projects via the pivot path");

    // Today's pivot-path output, pinned verbatim: the single-transcript c. path
    // renders the coding axis with its gene-symbol selector (not reframed under
    // the NG_ parent), a bare protein, and an NG_-framed genomic axis. The
    // Genome-gate must leave every one of these unchanged.
    assert_eq!(
        r.coding.as_ref().map(ToString::to_string).as_deref(),
        Some("NM_TEST.1(TESTGENE):c.4A>C"),
    );
    assert_eq!(
        r.protein.as_ref().map(ToString::to_string).as_deref(),
        Some("NP_TEST.1:p.(Lys2Gln)"),
    );
    assert_eq!(
        r.genomic.as_ref().map(ToString::to_string).as_deref(),
        Some("NG_TEST.1:g.5004A>C"),
    );
}

/// Regression: a genuinely-non-overlapping `NC_` genomic input still declines —
/// de-anchoring happens BEFORE the overlap guard, it does not weaken the guard.
#[test]
fn genuine_non_overlap_still_declines() {
    let (projector, provider) = base_fixture();
    let vp = VariantProjector::new(projector, provider);
    let variant = parse_hgvs("NC_000001.11:g.5000A>G").expect("parse");
    let err = vp
        .project_variant(&variant, "NM_TEST.1")
        .expect_err("a non-overlapping NC_ input must still decline");
    assert!(
        matches!(err, FerroError::TranscriptNotOverlapping { .. }),
        "expected TranscriptNotOverlapping, got: {err:?}"
    );
}
