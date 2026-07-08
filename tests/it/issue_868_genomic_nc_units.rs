//! Issue #868 — characterization tests for the two units extracted from
//! `project_to_genomic_nc`:
//!   1. the #785 silent-version-substitution gate
//!      (`enforce_exact_transcript_version`), and
//!   2. the step-5 coordinate map (`map_cnr_position_to_genome`), including the
//!      #797 poly-A `c.*` fallback and the non-coding `n.` arm.
//!
//! These exercise the outbound `c./n./r. → g.` path via the public
//! `VariantProjector::project_to_genomic` with a hermetic `MockProvider` + cdot
//! fixture (no manifest). They pin today's behavior so the #868 extraction is
//! provably byte-identical.
//!
//! Scope: these are **refactor-safety characterization tests**, not
//! independent-oracle conformance checks — they assert ferro's own output for
//! inputs that traverse the two extracted units, so they fail loudly if the
//! extraction changes behavior. Independent-oracle conformance (vs mutalyzer)
//! lives in the manifest-gated `mutalyzer_normalize_tests`.
//! <https://github.com/fulcrumgenomics/ferro-hgvs/issues/868>

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{parse_hgvs, FerroError, VariantProjector};

/// Build a deterministic sequence by repeating `pat` up to `len` bases (shared
/// by the fixtures that need a long genomic flank).
fn cycle(pat: &str, len: usize) -> String {
    pat.chars().cycle().take(len).collect()
}

/// (a) #785 gate — the input names an explicit transcript version the reference
/// does NOT carry exactly (`NM_000123.1`), but a *sibling* version resolves in
/// cdot (`NM_000123.2`). The outbound `c. → g.` projection must decline with
/// `TranscriptVersionNotExact` rather than silently projecting onto the
/// sibling's frame.
fn version_gate_fixture() -> VariantProjector<MockProvider> {
    let mut cdot = CdotMapper::new();
    // Only the SIBLING version (.2) lives in cdot; a request for .1 resolves to
    // it via cdot's base-accession version fallback (`base_to_versioned`).
    cdot.add_transcript(
        "NM_000123.2".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("GENE123".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            exons: vec![[1000, 1010, 0, 10]],
            cds_start: Some(0),
            cds_end: Some(10),
            gene_id: None,
            protein: Some("NP_000123.2".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    // The reference does not carry the exact requested version .1 — only a
    // sibling's bases. This is the precondition the #785 gate keys off.
    provider.mark_non_version_exact("NM_000123.1");
    VariantProjector::new(projector, provider)
}

#[test]
fn version_gate_declines_cross_version_substitution_outbound() {
    let vp = version_gate_fixture();
    // NG_ parent carries no build tag, so `build_hint` is None and the gate uses
    // the build-agnostic cdot lookup — which resolves .1 → .2 (sibling present),
    // firing the decline.
    let v = parse_hgvs("NG_000123.1(NM_000123.1):c.1A>G").expect("parse");
    let err = vp.project_to_genomic(&v).expect_err(
        "an explicit non-exact transcript version with a resolving sibling must decline",
    );
    match err {
        FerroError::TranscriptVersionNotExact { requested } => {
            assert_eq!(
                requested, "NM_000123.1",
                "the requested version must be surfaced verbatim"
            );
        }
        other => panic!("expected TranscriptVersionNotExact, got: {other:?}"),
    }
}

#[test]
fn version_gate_permits_exact_version_outbound() {
    // Control: the SAME fixture, but requesting the exact sibling version .2
    // (which the reference does carry exactly) must NOT trip the gate — it
    // projects the c. position onto the genome. This proves the gate is a
    // cross-version guard, not a blanket rejection. An NC_ parent is already in
    // its own (chromosome) frame, so the pivot passes through without needing a
    // placement to re-anchor.
    let vp = version_gate_fixture();
    let v = parse_hgvs("NC_000001.11(NM_000123.2):c.1A>G").expect("parse");
    let out = vp
        .project_to_genomic(&v)
        .expect("an exact-version request must project, not decline");
    // NC_ pivot returned in chromosome frame; c.1 (cds_start=0) → genome 1000.
    let s = match out {
        HgvsVariant::Genome(ref g) => g.to_string(),
        other => panic!("expected Genome variant, got: {other:?}"),
    };
    assert!(
        s.contains(":g.1000"),
        "expected g.1000 for c.1 on the exact version, got: {s}"
    );
}

/// (b) #797 poly-A `c.*` → g. — a coding transcript whose exon→genome map strips
/// the post-transcriptional poly-A tail. A `c.*4` endpoint lands in the poly-A
/// region, so the CDS-aware map (`cds_to_genome`) declines; the fallback
/// (`try_extend_polya_to_genome`) resolves it to the contiguous downstream
/// genome coordinate. Mirrors the in-module `make_polya_test_provider_and_projector`.
fn polya_fixture() -> VariantProjector<MockProvider> {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_POLYA.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("POLYAGENE".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            // Exon covers only the genomic core (poly-A tail stripped): tx 0..12.
            exons: vec![[1000, 1012, 0, 12]],
            cds_start: Some(3),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_POLYA.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    // FASTA carries the full 15-base transcript incl. the 3-base poly-A tail;
    // the exon map above only spans tx 0..12, so c.*4..c.*6 (tx 12..15) is the
    // stripped poly-A region the fallback walks into.
    provider.add_transcript(Transcript::new(
        "NM_POLYA.1".to_string(),
        Some("POLYAGENE".to_string()),
        TxStrand::Plus,
        "TTTATGCGCGTAAAA".to_string(), // 15 bases, trailing "AAA" tail
        Some(4),                       // 1-based inclusive CDS start (tx index 3)
        Some(9),
        vec![Exon::with_genomic(1, 1, 12, 1000, 1011)],
        Some("NC_000001.11".to_string()),
        Some(1000),
        Some(1011),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    // Genomic backbone: 999 non-repeating leading bases so HGVS g.1000 == 0-based
    // index 999, then the 12-base exon core, then a non-repeating downstream tail
    // (so a single-base del at the walked coordinate stays put under genomic-frame
    // renormalization). Registered under the NC_ output accession too so the
    // public projector's renormalization can fetch bases.
    let genomic = format!(
        "{}{}{}",
        cycle("ACGT", 999),
        "TTTATGCGCGTA",
        cycle("CAGT", 100)
    );
    provider.add_genomic_sequence("NC_000001.11", genomic);
    VariantProjector::new(projector, provider)
}

#[test]
fn polya_core_3utr_projects_normally() {
    // Control: a `c.*` in the genomic CORE of the 3'UTR (tx 9, still inside the
    // exon) projects via the normal CDS-aware map, NOT the fallback. c.*1 → g.1009.
    let vp = polya_fixture();
    let v = parse_hgvs("NC_000001.11(NM_POLYA.1):c.*1del").expect("parse");
    let out = vp
        .project_to_genomic(&v)
        .expect("genomic-core 3'UTR must project");
    let s = match out {
        HgvsVariant::Genome(ref g) => g.to_string(),
        other => panic!("expected Genome variant, got: {other:?}"),
    };
    assert!(
        s.contains(":g.1009del"),
        "expected g.1009del for c.*1del, got: {s}"
    );
}

#[test]
fn polya_region_walks_into_downstream_genome() {
    // c.*4 = tx 12 = the exon's genome_end (exclusive) → the first stripped
    // poly-A base → genome 1012, resolved only by the #797 fallback walk.
    let vp = polya_fixture();
    let v = parse_hgvs("NC_000001.11(NM_POLYA.1):c.*4del").expect("parse");
    let out = vp
        .project_to_genomic(&v)
        .expect("c.*4del in the poly-A region must project via the contiguous walk");
    let s = match out {
        HgvsVariant::Genome(ref g) => g.to_string(),
        other => panic!("expected Genome variant, got: {other:?}"),
    };
    assert!(
        s.contains(":g.1012del"),
        "expected g.1012del for c.*4del (poly-A walk), got: {s}"
    );
}

/// (c) non-coding `n.` → g. outbound — exercises the non-coding arm of the
/// coordinate map (`tx_to_genome`), distinct from the coding CDS-aware path.
#[test]
fn noncoding_n_position_projects_via_tx_to_genome() {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NR_000200.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("NCGENE".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            exons: vec![[1000, 1010, 0, 10]],
            cds_start: None,
            cds_end: None,
            gene_id: None,
            protein: None,
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NR_000200.1".to_string(),
        Some("NCGENE".to_string()),
        TxStrand::Plus,
        "ACGTACGTAC".to_string(), // 10 bases
        None,
        None,
        vec![Exon::with_genomic(1, 1, 10, 1000, 1009)],
        Some("NC_000001.11".to_string()),
        Some(1000),
        Some(1009),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    let genomic = format!(
        "{}{}{}",
        cycle("ACGT", 999),
        "ACGTACGTAC",
        cycle("CAGT", 100)
    );
    provider.add_genomic_sequence("NC_000001.11", genomic);
    let vp = VariantProjector::new(projector, provider);

    // n.5 (1-based tx 5 → 0-based tx 4) → genome 1000 + 4 = 1004 on the plus strand.
    let v = parse_hgvs("NC_000001.11(NR_000200.1):n.5del").expect("parse");
    let out = vp
        .project_to_genomic(&v)
        .expect("non-coding n. must project");
    let s = match out {
        HgvsVariant::Genome(ref g) => g.to_string(),
        other => panic!("expected Genome variant, got: {other:?}"),
    };
    assert!(
        s.contains(":g.1004del"),
        "expected g.1004del for n.5del, got: {s}"
    );
}
