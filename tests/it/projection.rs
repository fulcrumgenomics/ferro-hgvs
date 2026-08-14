//! End-to-end tests for the variant-projection module.

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::{ReferenceProvider, Strand};
use ferro_hgvs::{FerroError, VariantProjection, VariantProjector};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

/// Build a test (Projector, MockProvider) pair for NM_TEST.1.
///
/// The transcript encodes "ATGCGCTAA" (Met-Arg-Stop, 3 codons) on chr1 plus
/// strand at genomic positions 1000-1008 (1-based inclusive).
///
/// Cdot coordinate layout (0-based):
///   exon: genome [1000, 1009) — tx [0, 9)
///   cds_start = 0 (no 5' UTR), cds_end = 9
///
/// Coordinate mapping (cdot 0-based genome → HGVS c. 1-based):
///   g.1000 → tx_pos 0 → c.1  (A)
///   g.1001 → tx_pos 1 → c.2  (T)
///   g.1002 → tx_pos 2 → c.3  (G)
///   g.1003 → tx_pos 3 → c.4  (C ← ref base for the substitution test)
///   ...
fn fixture() -> (Projector, MockProvider) {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TEST.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("TESTGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            // [genome_start(0-based), genome_end(0-based excl), tx_start, tx_end]
            exons: vec![[1000, 1009, 0, 9]],
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
    // Genomic sequence: 999 N's + "ATGCGCTAA" + 100 N's.
    // 0-based index 999 = 'A', 1000 = 'T', 1001 = 'G', 1002 = 'C', ...
    let prefix = "N".repeat(999);
    let suffix = "N".repeat(100);
    provider.add_genomic_sequence("chr1", format!("{}ATGCGCTAA{}", prefix, suffix));
    (projector, provider)
}

#[test]
fn end_to_end_missense() {
    let (projector, provider) = fixture();
    let vp = VariantProjector::new(projector, provider);
    // g.1003C>A: the 4th base (c.4) of codon 2 "CGC" (Arg) → "AGC" (Ser).
    let result: VariantProjection = vp
        .project("NC_000001.11:g.1003C>A", "NM_TEST.1")
        .expect("projection should succeed");
    assert_eq!(result.transcript_id, "NM_TEST.1");
    let c = result.coding.as_ref().unwrap().to_string();
    assert!(c.contains(":c.4C>A"), "got c. = {}", c);
    let p = result.protein.as_ref().unwrap().to_string();
    assert_eq!(p, "NP_TEST.1:p.(Arg2Ser)");
    assert!(!result.is_frameshift);
    assert!(!result.is_intronic);
    assert!(!result.is_utr);
}

#[test]
fn end_to_end_no_overlap() {
    let (projector, provider) = fixture();
    let vp = VariantProjector::new(projector, provider);
    // chr1:5000 is far outside the transcript at genome [1000, 1009).
    let err = vp.project("NC_000001.11:g.5000A>G", "NM_TEST.1");
    assert!(err.is_err(), "expected error for off-transcript position");
}

// =============================================================================
// Issue #332: `VariantProjector` must also route transcript lookups through
// the variant-aware `ReferenceProvider::get_transcript_for_variant` so the
// projector benefits from the same NG/NC-parent build resolution as
// `Normalizer::normalize`. Confirms the projector resolves an intronic c.
// input that carries an NG_* parent without erroring on the codon-fetch step.
// =============================================================================

/// Intronic-aware fixture: a 3-exon transcript on chr1 plus strand, intron
/// homopolymers, so intronic 3' shifting is well-defined. The projector path
/// only requires that the substitution-codon-fetch and indel-codon-fetch
/// sites use the variant-aware lookup; this fixture exercises an exonic g.
/// substitution against an NG-parented HGVS form to confirm the projector's
/// `provider.get_transcript_for_variant(variant)` plumbing is in place.
fn intronic_fixture() -> (Projector, MockProvider) {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TEST.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("TESTGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            // 3 exons, 30bp total, plus 2 introns of 50bp each.
            // Exon 1: g.[1001,1010] tx [0,10)
            // Exon 2: g.[1061,1070] tx [10,20)
            // Exon 3: g.[1121,1130] tx [20,30)
            exons: vec![
                [1000, 1010, 0, 10],
                [1060, 1070, 10, 20],
                [1120, 1130, 20, 30],
            ],
            cds_start: Some(0),
            cds_end: Some(30),
            gene_id: None,
            protein: Some("NP_TEST.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    let tx_seq = "ATGCGCTAAATGCGCTAAATGCGCTAACGT"; // 30 bases
    provider.add_transcript(Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TESTGENE".to_string()),
        TxStrand::Plus,
        tx_seq.to_string(),
        Some(1),
        Some(30),
        vec![
            Exon::with_genomic(1, 1, 10, 1001, 1010),
            Exon::with_genomic(2, 11, 20, 1061, 1070),
            Exon::with_genomic(3, 21, 30, 1121, 1130),
        ],
        Some("chr1".to_string()),
        Some(1001),
        Some(1130),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    // Plant a deterministic genomic sequence on chr1.
    let mut bytes = vec![b'N'; 2000];
    let plant = |bytes: &mut Vec<u8>, start_1: usize, seq: &str| {
        for (i, b) in seq.bytes().enumerate() {
            bytes[start_1 - 1 + i] = b;
        }
    };
    plant(&mut bytes, 1001, "ATGCGCTAAA");
    plant(&mut bytes, 1011, &"A".repeat(50));
    plant(&mut bytes, 1061, "TGCGCTAAAT");
    plant(&mut bytes, 1071, &"T".repeat(50));
    plant(&mut bytes, 1121, "GCGCTAACGT");
    provider.add_genomic_sequence("chr1", String::from_utf8(bytes).unwrap());
    (projector, provider)
}

#[test]
fn project_intronic_c_dot_succeeds() {
    // Pins that intronic g.→c. projection runs through `Normalizer::normalize`
    // (intronic path, which uses the variant-aware lookup after #332) and the
    // projector's per-codon transcript fetches without erroring.
    //
    // The input here is `g.`-form (the projector's primary supported entry
    // point); the per-codon transcript lookups inside the projector consult
    // `get_transcript_for_variant` on internally-constructed `CdsVariant`s
    // (which do NOT carry `genomic_context` today — that's a separate gap
    // tracked by #328). Direct NG-parented projection coverage will land
    // once #328 ships.
    let (projector, provider) = intronic_fixture();
    let vp = VariantProjector::new(projector, provider);
    let result = vp
        .project("NC_000001.11:g.1012A>T", "NM_TEST.1")
        .expect("projection must succeed for intronic input");
    assert_eq!(result.transcript_id, "NM_TEST.1");
    assert!(result.is_intronic);
    let c = result.coding.as_ref().unwrap().to_string();
    // Pin the position number — a regression that mis-identifies the
    // intronic offset would still match a loose `c.10+` substring.
    assert!(
        c.contains("c.10+3"),
        "expected intronic offset c.10+3; got {}",
        c
    );
    // Intronic substitutions skip protein prediction.
    assert!(result.protein.is_none(), "no p. for intronic substitutions");
}

// ─── transcript-cache regression test ────────────────────────────────────────
//
// A counting wrapper around any `ReferenceProvider` that tallies how many times
// each downstream method is invoked. Used below to assert that the projector's
// transcript cache collapses N identical lookups into one.

#[derive(Clone)]
struct CountingProvider<P: ReferenceProvider + Clone> {
    inner: P,
    transcript_calls: Arc<AtomicUsize>,
}

impl<P: ReferenceProvider + Clone> CountingProvider<P> {
    fn new(inner: P) -> Self {
        Self {
            inner,
            transcript_calls: Arc::new(AtomicUsize::new(0)),
        }
    }
}

impl<P: ReferenceProvider + Clone> ReferenceProvider for CountingProvider<P> {
    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError> {
        self.transcript_calls.fetch_add(1, Ordering::SeqCst);
        self.inner.get_transcript(id)
    }
    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        self.inner.get_sequence(id, start, end)
    }
    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        self.inner.get_genomic_sequence(contig, start, end)
    }
    fn has_genomic_data(&self) -> bool {
        self.inner.has_genomic_data()
    }
    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        self.inner.get_protein_sequence(accession, start, end)
    }
    fn has_protein_data(&self) -> bool {
        self.inner.has_protein_data()
    }
}

/// Project many substitutions on a single transcript and assert that the
/// projector's own fetches of that transcript collapse.
///
/// Without the cache, every substitution variant pulls the transcript from the
/// provider during protein prediction. With the cache, a `VariantProjector`
/// fetches each unique `transcript_id` at most once *through
/// `cached_get_transcript_for_variant`* for the lifetime of the projector.
///
/// # The pin is 10, not 1, and that is a REGRESSION recorded rather than fixed
///
/// The count is not the cache's hit rate; it is the total over the provider,
/// and the projector's cache does not sit in front of every path that reaches
/// it. The **normalizer** the projector holds fetches from the provider
/// directly, so any axis the projector re-derives by normalizing costs one
/// transcript fetch per variant however warm the projector's own cache is.
///
/// #1712 made the `n.` axis one of those axes
/// (`rulings[noncoding-axis-is-re-derived-on-its-own-reference]`, operator
/// ruling 2026-08-13), which took this count from **1** to **10** across the
/// same three variants. The pin is re-blessed rather than the reach gap closed,
/// deliberately: the gap is a performance defect in where the cache sits, it
/// predates #1712 on the genomic axis (`reanchor_and_normalize_genomic` and
/// `project_to_genomic_normalized` have normalized since #737/#867), and fixing
/// it inside an adjudication PR would mix a performance change into a ruling.
/// It is filed as **#1860**.
///
/// **This is not a weakened guard.** It was an exact equality and still is, so
/// it fails in both directions. More calls is a new regression. **Fewer means
/// #1860 was closed** — which is a good outcome and still a failure here, so
/// that closing cannot land without re-pinning this number and saying so.
/// Ten is not a round number anybody would reach for; it is what was measured,
/// and reading it as a target rather than as a record is the mistake to avoid.
#[test]
fn transcript_cache_collapses_repeat_lookups() {
    let (projector, provider) = fixture();
    let counting = CountingProvider::new(provider);
    let counter = counting.transcript_calls.clone();
    let vp = VariantProjector::new(projector, counting);

    // Three substitutions at three CDS positions on the same transcript. Each
    // hits the protein-prediction path, which is the user of `get_transcript`.
    for s in [
        "NC_000001.11:g.1003C>A",
        "NC_000001.11:g.1004G>A",
        "NC_000001.11:g.1005C>A",
    ] {
        vp.project(s, "NM_TEST.1")
            .expect("projection should succeed");
    }

    let calls = counter.load(Ordering::SeqCst);
    assert_eq!(
        calls, 10,
        "expected exactly ten get_transcript calls across 3 variants (got \
         {calls}). The projector's own fetches still collapse to one; the rest \
         are the normalizer's, which do not route through that cache — see \
         #1860 and this test's doc comment. FEWER than ten means #1860 was \
         closed: re-pin this number in that change rather than here."
    );
}

// =============================================================================
// Issue #1712: every axis the projector renders must be a fixed point of bare
// normalization — asserted MANIFEST-FREE.
// =============================================================================

/// A transcript with a 5'UTR and two repeat tracts, on `MockProvider`.
///
/// The 5'UTR is what makes this fixture say anything: with `cds_start = 1` the
/// `n.` and `c.` coordinates coincide and a reframing bug is invisible. Here
/// `c.1` sits at `n.3`, so the two axes carry different numbers.
///
/// ```text
/// n.  1  2 | 3  4  5 | 6  7  8 | 9 10 11 |12 13 14 |15 16 17 |18 19 20
/// c.       | 1  2  3 | 4  5  6 | 7  8  9 |10 11 12 |13 14 15 |16 17 18
///     G  G | A  T  G | T  T  T | T  G  C | G  C  G | C  G  C | T  A  A
/// ```
///
/// `c.4_7` is a four-base `T` homopolymer and `c.8_15` a `GC` dinucleotide
/// tract — the two shapes whose re-derivation on the `n.` reference the ruling
/// `rulings[noncoding-axis-is-re-derived-on-its-own-reference]` moves.
fn repeat_tract_fixture() -> (Projector, MockProvider) {
    const TX: &str = "GGATGTTTTGCGCGCTAA";
    let len = TX.len() as u64; // 18

    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_RPT.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("RPTGENE".to_string()),
            contig: "chr9".to_string(),
            strand: Strand::Plus,
            exons: vec![[2000, 2000 + len, 0, len]],
            cds_start: Some(2), // 0-based tx offset of `c.1` → n.3
            cds_end: Some(len),
            gene_id: None,
            protein: Some("NP_RPT.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_RPT.1".to_string(),
        Some("RPTGENE".to_string()),
        TxStrand::Plus,
        TX.to_string(),
        Some(3), // 1-based `c.1` position on the transcript
        Some(len),
        vec![Exon::new(1, 1, len)],
        Some("chr9".to_string()),
        Some(2001),
        Some(2000 + len),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    provider.add_genomic_sequence(
        "chr9",
        format!("{}{}{}", "N".repeat(2000), TX, "N".repeat(100)),
    );
    (projector_from(cdot), provider)
}

fn projector_from(cdot: CdotMapper) -> Projector {
    Projector::new(cdot)
}

/// **Manifest-free counterpart to `mutalyzer_normalize_tests::axis_noncoding_idempotent`.**
///
/// That test measures the real residue over 681 projected axes and is the one
/// that proved the ruling closes it — but it resolves its reference from
/// `FERRO_MANIFEST`, which **PR CI does not provision**, so it takes a skip path
/// `nextest` reports as PASS. It is nightly-only coverage. This test asserts the
/// same contract on a `MockProvider` fixture built programmatically, so the
/// invariant is guarded pre-merge.
///
/// The contract: whatever the projector hands a caller must already be what bare
/// normalization would produce. A caller that normalizes a projected axis — the
/// obvious defensive thing to do — must not see the string change under them.
///
/// Reverting `noncoding_axis_with_reason`'s `normalize_or_fallback` call fails
/// this test, which is what makes it a guard rather than a restatement.
#[test]
fn every_projected_noncoding_axis_is_a_fixed_point() {
    use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

    let (projector, provider) = repeat_tract_fixture();
    let normalizer = Normalizer::with_config(provider.clone(), NormalizeConfig::default());
    let vp = VariantProjector::new(projector, provider);

    // Deletions and insertions inside each tract, plus a plain substitution as
    // the negative control — a shape with nothing to maximize, which must be a
    // fixed point both before and after the ruling.
    let inputs = [
        "NM_RPT.1:c.5del",
        "NM_RPT.1:c.6del",
        "NM_RPT.1:c.4_5del",
        "NM_RPT.1:c.9_10del",
        "NM_RPT.1:c.8_9insGC",
        "NM_RPT.1:c.10_11insCG",
        "NM_RPT.1:c.2T>A",
    ];

    // The fixture must actually be able to PRODUCE the shape, or "nothing
    // drifted" below is a fact about the fixture rather than about the
    // projector. `c.4_5del` is the load-bearing row, and the two axes below are
    // the whole point of the ruling side by side:
    //
    // ```text
    // c.  c.6_7del      3'-shifted within the T tract, still a plain `del`
    // n.  n.6_9T[2]     re-derived on the `n.` reference, tract MAXIMIZED
    // ```
    //
    // Both denote the same two bases. The `c.` axis shifts but does not
    // maximize; the `n.` axis does, because it is normalized on its own
    // reference rather than reframed from the `c.` form. Reverting the fix makes
    // the `n.` row read `n.8_9del` — the plain CDS-offset image of the `c.`
    // answer — which is the repeat-form class: 16 of the 20 rows the ruling
    // moved on the real reference, and the majority shape the reference-backed
    // test measures.
    let anchor = vp
        .project_variant(
            &parse_hgvs("NM_RPT.1:c.4_5del").expect("anchor must parse"),
            "NM_RPT.1",
        )
        .expect("anchor must project");
    assert_eq!(
        anchor
            .noncoding
            .as_ref()
            .map(ToString::to_string)
            .as_deref(),
        Some("NM_RPT.1:n.6_9T[2]"),
        "the `n.` axis of `c.4_5del` must be the repeat form re-derived on the \
         transcript's own reference; if this fixture stopped producing that \
         shape the fixed-point assertion below would be quantifying over \
         nothing that ever moved"
    );
    // …and the `c.` axis keeps a plain `del`, which is what makes the pair a
    // divergence rather than a change of normalizer behaviour on both axes.
    // Pinned as the exact string, so a future change that started maximizing the
    // tract on `c.` too would fail here rather than quietly making the `n.`
    // assertion above trivial.
    assert_eq!(
        anchor.coding.as_ref().map(ToString::to_string).as_deref(),
        Some("NM_RPT.1:c.6_7del"),
        "the coding axis must be untouched by this ruling"
    );

    let mut checked = 0usize;
    let mut drifted = Vec::new();
    for input in inputs {
        let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));
        let projection = vp
            .project_variant(&variant, "NM_RPT.1")
            .unwrap_or_else(|e| panic!("{input} must project: {e}"));
        let Some(axis) = projection.noncoding.as_ref() else {
            panic!("{input}: the `n.` axis must render for this fixture to test anything");
        };
        let rendered = axis.to_string();
        checked += 1;
        let renormalized = normalizer
            .normalize(axis)
            .unwrap_or_else(|e| panic!("{rendered} must re-normalize: {e}"))
            .to_string();
        if renormalized != rendered {
            drifted.push(format!("{input} : {rendered} -> {renormalized}"));
        }
    }

    // Honest-zero discipline, the same as the reference-backed test's: a fixture
    // that stopped rendering `n.` axes would otherwise report no drift and read
    // as a pass. The `let ... else` above already panics per input; this asserts
    // the loop ran at all.
    assert_eq!(
        checked,
        inputs.len(),
        "every input must contribute an `n.` axis"
    );
    assert!(
        drifted.is_empty(),
        "the projector emitted `n.` axes that bare normalization rewrites (#1712). \
         The projector re-derives this axis on its own reference, so there is no \
         exemption list to add to:\n{}",
        drifted.join("\n")
    );
}
