//! Audit (item B4 remaining of tracking issue #81 / closes #209):
//! 3' rule for repeat unit boundary placement across the rest of the
//! matrix.
//!
//! The cyclic-rotation `ins` branch landed in #155 (closes #132). For
//! `del`/`dup` the shuffle phase-alignment lemma (`shuffle.rs` /
//! `rules.rs::deletion_to_repeat` doc) ensures that a post-shift slice
//! lands on the natural unit boundary, so the existing del/dup shift
//! matrices (#106, #107) cover plus/minus strand × g./c./n./r. ×
//! cyclic-rotation phase. The remaining gaps this file pins:
//!
//! 1. **Intronic `c.` positions** — `normalize_intronic_cds` previously
//!    passed `is_coding=true` into `normalize_na_edit`, which fired the
//!    codon-frame gate (`unit_len % 3 != 0`) on repeat normalization
//!    *inside introns*. Per `assets/hgvs-nomenclature/docs/
//!    recommendations/DNA/repeated.md`:
//!
//!       > This restriction only applies to the coding sequence,
//!       > which does not include the introns or the UTR sequence.
//!
//!    So an intronic A-homopolymer dup must emit `A[N+k]`, not the
//!    gated `ins<literal>` fallback. The companion intronic `n.`
//!    parser already passed `is_coding=false` (correct).
//!
//! 2. **UTR `c.` positions** — `normalize_cds` (exonic path) previously
//!    passed `is_coding=true` unconditionally, including for 5' UTR
//!    (`c.-N`) and 3' UTR (`c.*N`) positions which are also exempt from
//!    the codon-frame restriction per the same spec clause. Spec
//!    example: `NM_024312.4:c.-6_-3G[6]` is valid even though
//!    `unit_len=1`.
//!
//! 3. **Boundary-spanning `c.` positions** — `normalize_boundary_spanning_cds`
//!    previously passed `is_coding=true`. Variants that cross an
//!    exon/intron boundary do not sit purely in coding sequence, so
//!    they follow the same intronic exemption.
//!
//! (Mitochondrial `m.` normalization is a separate stub-level gap
//! discovered while writing this audit — `normalize_mt` does not run
//! window-based shuffling or repeat-notation canonicalization. That
//! work is tracked under #210 and intentionally out of scope here so
//! this PR stays focused on the codon-frame gate fix.)

mod common;

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// 256 bases of padding on each side — matches the
/// `common::synthetic::PAD_OFFSET` contract so the normalizer's 100bp
/// shuffle window stays in range. Pure `ACGT` repeats so no accidental
/// tract overlap with the test cores.
const PAD: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
                   ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
                   ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
                   ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
const PAD_OFFSET: u64 = 256;

/// Round-trip helper: normalize, then re-parse the result and re-normalize
/// to confirm the canonical form is a fixed point.
fn assert_canonical_round_trip(provider: MockProvider, input: &str, expected: &str) {
    let normalizer = Normalizer::new(provider);
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for {:?}: {:?}", input, e));
    let n = normalizer
        .normalize(&v)
        .unwrap_or_else(|e| panic!("normalize failed for {:?}: {:?}", input, e));
    let n_str = format!("{}", n);
    assert_eq!(n_str, expected, "first-pass normalization mismatch");

    let v2 =
        parse_hgvs(&n_str).unwrap_or_else(|e| panic!("re-parse failed for {:?}: {:?}", n_str, e));
    let n2 = normalizer
        .normalize(&v2)
        .unwrap_or_else(|e| panic!("re-normalize failed for {:?}: {:?}", n_str, e));
    assert_eq!(
        format!("{}", n2),
        expected,
        "canonical form is not idempotent"
    );
}

/// Reverse-complement a DNA string (for minus-strand provider construction).
fn rc(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            other => other,
        })
        .collect()
}

// =============================================================================
// Section 1: Intronic c. positions — codon-frame restriction does NOT apply.
//
// Build a small 2-exon plus-strand transcript whose intron carries an
// `A[5]` homopolymer (unit_len=1, never codon-aligned) and a `TG[4]`
// dinucleotide tract (unit_len=2, also blocked by the gate). Both
// would currently route to gated `ins<literal>` because
// `normalize_intronic_cds` passes `is_coding=true`. Per `repeated.md`
// the restriction applies only to the coding sequence — introns are
// exempt.
// =============================================================================

mod intronic_plus_strand {
    use super::*;

    /// 2-exon plus-strand transcript layout (`p = PAD_OFFSET + 1`; genomic
    /// coordinates are 1-based, matching real cdot data — `Exon::with_genomic`
    /// genomic positions are 1-based HGVS positions, so the first core base
    /// (0-based slice index `PAD_OFFSET`) is at 1-based genomic position
    /// `p = PAD_OFFSET + 1`; exon 1 spans 1-based `[p, p+29]`):
    /// ```text
    ///   Exon 1: tx 1-30   genomic p    .. p+29
    ///   Intron 1 (30 bp): genomic p+30 .. p+59
    ///     content (1-based intron offset):
    ///       offset  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 ...
    ///       base    C  C  A  A  A  A  A  C  C  T  G  T  G  T  G  T  G  C..C
    ///                     `--- A[5] ---`     `------ TG[4] ------`
    ///   Exon 2: tx 31-60  genomic p+60 .. p+89
    /// ```
    fn make_provider() -> MockProvider {
        let intron_content = "CCAAAAACCTGTGTGTGCCCCCCCCCCCCC";
        assert_eq!(intron_content.len(), 30, "intron must be exactly 30 bp");
        let exon1 = "GAATTCAAAAGGCCTTCCGGAACCGGTAAA";
        let exon2 = "TTGGAACCGGAATTCCAAGGCCAATTGGTT";
        let core = format!("{}{}{}", exon1, intron_content, exon2);
        let genomic = format!("{}{}{}", PAD, core, PAD);

        let p = PAD_OFFSET + 1;
        let tx_seq = format!("{}{}", exon1, exon2);
        let transcript = Transcript::new(
            "NM_INTRON.1".to_string(),
            Some("INTRONGENE".to_string()),
            Strand::Plus,
            tx_seq,
            Some(1),
            Some(60),
            vec![
                Exon::with_genomic(1, 1, 30, p, p + 29),
                Exon::with_genomic(2, 31, 60, p + 60, p + 89),
            ],
            Some("chr_intron".to_string()),
            Some(p),
            Some(p + 89),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );

        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("chr_intron", genomic);
        provider.add_transcript(transcript);
        provider
    }

    #[test]
    fn intronic_homopolymer_dup_emits_repeat_notation() {
        // Intronic A[5] tract at c.30+3..c.30+7. Dup of "AA" at
        // c.30+3_30+4 should emit `A[7]` over the reference-tract
        // extent. Pre-fix: codon-frame gate fires (`is_coding=true`
        // for intronic CDS), routing to `c.30+7_30+8insAA`. Post-fix:
        // intronic positions pass `is_coding=false`, allowing repeat
        // notation.
        let p = make_provider();
        assert_canonical_round_trip(
            p,
            "NM_INTRON.1:c.30+3_30+4dup",
            "NM_INTRON.1:c.30+3_30+7A[7]",
        );
    }

    #[test]
    fn intronic_homopolymer_del_emits_repeat_notation() {
        // Intronic A[5] at c.30+3..c.30+7. Del of "AA" at c.30+3_30+4
        // should emit `A[3]`. Pre-fix: `deletion_to_repeat`
        // short-circuits at the codon-frame gate and the dup arm
        // falls through to plain `del` at the post-shift position.
        let p = make_provider();
        assert_canonical_round_trip(
            p,
            "NM_INTRON.1:c.30+3_30+4del",
            "NM_INTRON.1:c.30+3_30+7A[3]",
        );
    }

    // Note: a single-copy multi-base intronic dup (e.g. `c.30+10_30+11dup`
    // of one TG unit adjacent to TG[4]) currently emits the 3'-shifted
    // dup at a position one base earlier than the expected 3'-most
    // unit-aligned boundary. The discrepancy is in the intronic
    // position-conversion pipeline rather than the codon-frame gate
    // (single-copy paths don't enter the gate), so it is out of scope
    // for this PR and tracked separately. We do NOT pin it here.

    #[test]
    fn intronic_dinucleotide_dup_multi_copy_emits_repeat_notation() {
        // Intronic TG[4] at intron offsets 10..17 (8 bases). Dup of
        // "TGTG" at c.30+10_30+13 (2 copies, phase-aligned). Should
        // emit `TG[6]` over the reference-tract extent. unit_len=2 —
        // codon-frame gate fires pre-fix (`is_coding=true` rejects
        // unit_len % 3 != 0); post-fix allows repeat notation.
        //
        // The tract's end position (intron offset 17 of 30) sits past
        // the midpoint, so ferro renders it acceptor-relative as
        // `c.31-14` (per `find_intron_at_genomic`'s
        // `dist_to_5prime <= dist_to_3prime` tie-breaker).
        let p = make_provider();
        assert_canonical_round_trip(
            p,
            "NM_INTRON.1:c.30+10_30+13dup",
            "NM_INTRON.1:c.30+10_31-14TG[6]",
        );
    }

    #[test]
    fn intronic_dinucleotide_del_multi_copy_emits_repeat_notation() {
        // Intronic TG[4] at intron offsets 10..17. Del of "TGTG" at
        // c.30+10_30+13. Should emit `TG[2]`. Same codon-frame gate
        // issue as the dup case; same acceptor-relative rendering of
        // the past-midpoint endpoint.
        let p = make_provider();
        assert_canonical_round_trip(
            p,
            "NM_INTRON.1:c.30+10_30+13del",
            "NM_INTRON.1:c.30+10_31-14TG[2]",
        );
    }

    #[test]
    fn intronic_homopolymer_multi_copy_ins_emits_repeat_notation() {
        // Intronic A[5] tract at c.30+3..c.30+7. Inserting two more As
        // adjacent to the tract (`c.30+7_30+8insAA`) should emit
        // `A[7]` over the reference-tract extent via
        // `insertion_to_repeat`. The `is_coding` flag flows into that
        // helper too (line 195 of rules.rs); the intronic exemption
        // must cover the `ins` arm as well as `dup`/`del`.
        let p = make_provider();
        assert_canonical_round_trip(
            p,
            "NM_INTRON.1:c.30+7_30+8insAA",
            "NM_INTRON.1:c.30+3_30+7A[7]",
        );
    }
}

// =============================================================================
// Section 2: Intronic minus-strand variants — orientation + codon-frame
// exemption together.
//
// Pins the spec-correct form for an intronic dup on a minus-strand
// transcript. The orientation flip in `normalize_intronic_cds` already
// lands the bytes in transcript view (#98); the only remaining bug is
// the codon-frame gate firing in intronic context.
// =============================================================================

mod intronic_minus_strand {
    use super::*;

    /// 2-exon minus-strand transcript with the intron carrying an
    /// `A[5]` homopolymer at c.30+3..c.30+7 *in transcript view*.
    ///
    /// Genomic layout (left-to-right, with `p = PAD_OFFSET`):
    /// ```text
    ///   Exon 2 RC: genomic p+1   .. p+30   (tx 31..60 maps here)
    ///   Intron:    genomic p+31  .. p+60   (tx-view RC of the desired content)
    ///   Exon 1 RC: genomic p+61  .. p+90   (tx 1..30 maps here)
    /// ```
    fn make_provider() -> MockProvider {
        // Transcript-view intron content (5'→3' on the transcript):
        // 2 bp filler + A[5] homopolymer + 23 bp filler = 30 bp.
        let intron_tx_view = "CCAAAAACCCCCCCCCCCCCCCCCCCCCCC";
        assert_eq!(intron_tx_view.len(), 30);
        // Genomic plus-strand intron sequence is the RC of the
        // transcript view.
        let intron_genomic = rc(intron_tx_view);

        let exon1_tx_view = "GAATTCAAAAGGCCTTCCGGAACCGGTAAA";
        let exon1_genomic = rc(exon1_tx_view);

        let exon2_tx_view = "TTGGAACCGGAATTCCAAGGCCAATTGGTT";
        let exon2_genomic = rc(exon2_tx_view);

        let core = format!("{}{}{}", exon2_genomic, intron_genomic, exon1_genomic);
        let genomic = format!("{}{}{}", PAD, core, PAD);

        let p = PAD_OFFSET + 1;
        let tx_seq = format!("{}{}", exon1_tx_view, exon2_tx_view);
        let transcript = Transcript::new(
            "NM_MINTRON.1".to_string(),
            Some("MINTRONGENE".to_string()),
            Strand::Minus,
            tx_seq,
            Some(1),
            Some(60),
            vec![
                // Exon 1 (coding start) maps to the rightmost genomic block.
                // 1-based genomic ranges (matching real cdot data).
                Exon::with_genomic(1, 1, 30, p + 60, p + 89),
                Exon::with_genomic(2, 31, 60, p, p + 29),
            ],
            Some("chr_mintron".to_string()),
            Some(p),
            Some(p + 89),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );

        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("chr_mintron", genomic);
        provider.add_transcript(transcript);
        provider
    }

    #[test]
    fn minus_strand_intronic_homopolymer_dup_emits_repeat_notation() {
        // Transcript-view A[5] in intron at c.30+3..c.30+7. Dup of
        // "AA" at c.30+3_30+4 → `A[7]`. Mirrors the plus-strand test
        // in section 1; the only difference is the orientation. The
        // existing `tests/coverage_gap_tests.rs::
        // test_minus_strand_intronic_multi_base_dup` now pins the
        // spec-correct repeat form (`c.30+2_30+4A[6]`), matching this
        // gate-exemption behavior for intronic `c.` positions.
        let p = make_provider();
        assert_canonical_round_trip(
            p,
            "NM_MINTRON.1:c.30+3_30+4dup",
            "NM_MINTRON.1:c.30+3_30+7A[7]",
        );
    }

    #[test]
    fn minus_strand_intronic_homopolymer_del_emits_repeat_notation() {
        // Del of "AA" at c.30+3_30+4 of intronic A[5] → `A[3]`.
        let p = make_provider();
        assert_canonical_round_trip(
            p,
            "NM_MINTRON.1:c.30+3_30+4del",
            "NM_MINTRON.1:c.30+3_30+7A[3]",
        );
    }
}

// =============================================================================
// Section 3: 5' UTR c. positions — codon-frame restriction does NOT apply.
//
// Per `repeated.md`: `NM_024312.4:c.-6_-3G[6]` is valid even though
// unit_len=1 (homopolymer is never codon-aligned), because the variant
// is in 5' UTR (before `c.1`). We pin the equivalent rule on the
// dup/del side.
// =============================================================================

mod utr_5prime {
    use super::*;

    /// Single-exon plus-strand transcript:
    /// - tx 1..30 is 5' UTR (`c.-30..c.-1`); contains `AAAAA` at tx
    ///   pos 5..9 = c.-26..c.-22
    /// - tx 31..60 is CDS (`c.1..c.30`)
    fn make_provider() -> MockProvider {
        let utr5 = "GCCGAAAAACCGGTACCGGTACCGGTACCG"; // 30 bp
        let cds = "ATGCATGCATGCATGCATGCATGCATGCAT"; // 30 bp
        let core = format!("{}{}", utr5, cds);
        let genomic = format!("{}{}{}", PAD, core, PAD);

        let p = PAD_OFFSET + 1;
        let tx_seq = format!("{}{}", utr5, cds);
        let transcript = Transcript::new(
            "NM_UTR5.1".to_string(),
            Some("UTR5GENE".to_string()),
            Strand::Plus,
            tx_seq,
            Some(31), // CDS starts at tx 31 → c.1
            Some(60), // CDS ends at tx 60 → c.30
            // 1-based genomic range: 60 bases at [p, p+60).
            vec![Exon::with_genomic(1, 1, 60, p, p + 59)],
            Some("chr_utr5".to_string()),
            Some(p),
            Some(p + 59),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );

        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("chr_utr5", genomic);
        provider.add_transcript(transcript);
        provider
    }

    #[test]
    fn utr5_homopolymer_dup_emits_repeat_notation() {
        // 5' UTR A[5] tract at c.-26..c.-22 (tx pos 5..9). Dup of
        // "AA" at c.-26_-25 → `A[7]` over the reference-tract extent.
        // Pre-fix: `is_coding=true` is passed unconditionally for
        // exonic CDS variants, including UTR — codon-frame gate fires
        // and routes to `c.-22_-21insAA`. Post-fix: UTR positions get
        // `is_coding=false`.
        let p = make_provider();
        assert_canonical_round_trip(p, "NM_UTR5.1:c.-26_-25dup", "NM_UTR5.1:c.-26_-22A[7]");
    }

    #[test]
    fn utr5_homopolymer_del_emits_repeat_notation() {
        // Del of "AA" at c.-26_-25 of A[5] tract → `A[3]`.
        let p = make_provider();
        assert_canonical_round_trip(p, "NM_UTR5.1:c.-26_-25del", "NM_UTR5.1:c.-26_-22A[3]");
    }
}

// =============================================================================
// Section 4: 3' UTR c. positions — codon-frame restriction does NOT apply.
// =============================================================================

mod utr_3prime {
    use super::*;

    /// Single-exon plus-strand transcript:
    /// - tx 1..30 is CDS (`c.1..c.30`)
    /// - tx 31..60 is 3' UTR (`c.*1..c.*30`); contains `AAAAA` at
    ///   tx pos 35..39 = c.*5..c.*9
    fn make_provider() -> MockProvider {
        let cds = "ATGCATGCATGCATGCATGCATGCATGCAT"; // 30 bp
        let utr3 = "CGCGAAAAACGGTACCGGTACCGGTACCGG"; // 30 bp
        let core = format!("{}{}", cds, utr3);
        let genomic = format!("{}{}{}", PAD, core, PAD);

        let p = PAD_OFFSET + 1;
        let tx_seq = format!("{}{}", cds, utr3);
        let transcript = Transcript::new(
            "NM_UTR3.1".to_string(),
            Some("UTR3GENE".to_string()),
            Strand::Plus,
            tx_seq,
            Some(1),
            Some(30),
            // 1-based genomic range: 60 bases at [p, p+60).
            vec![Exon::with_genomic(1, 1, 60, p, p + 59)],
            Some("chr_utr3".to_string()),
            Some(p),
            Some(p + 59),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );

        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("chr_utr3", genomic);
        provider.add_transcript(transcript);
        provider
    }

    #[test]
    fn utr3_homopolymer_dup_emits_repeat_notation() {
        // 3' UTR A[5] tract at c.*5..c.*9. Dup of "AA" at c.*5_*6 →
        // `A[7]`. Same codon-frame-gate-shouldn't-apply rationale.
        let p = make_provider();
        assert_canonical_round_trip(p, "NM_UTR3.1:c.*5_*6dup", "NM_UTR3.1:c.*5_*9A[7]");
    }

    #[test]
    fn utr3_homopolymer_del_emits_repeat_notation() {
        // Del of "AA" at c.*5_*6 of A[5] tract → `A[3]`.
        let p = make_provider();
        assert_canonical_round_trip(p, "NM_UTR3.1:c.*5_*6del", "NM_UTR3.1:c.*5_*9A[3]");
    }
}

// =============================================================================
// Section 5: Sanity — CDS-proper context still gates the codon-frame
// restriction (regression guard).
// =============================================================================

mod cds_proper_still_gates {
    use super::*;

    /// Single-exon plus-strand transcript whose entire 30 bp is CDS,
    /// with `AAAAA` (A[5]) at tx pos 4..8 = `c.4..c.8`.
    fn make_provider() -> MockProvider {
        let tx = "GCCAAAAACCGGTACCGGTACCGGTACCGT"; // 30 bp; A[5] at tx 4..8
        let genomic = format!("{}{}{}", PAD, tx, PAD);
        let p = PAD_OFFSET + 1;
        let transcript = Transcript::new(
            "NM_CDSP.1".to_string(),
            Some("CDSPGENE".to_string()),
            Strand::Plus,
            tx.to_string(),
            Some(1),
            Some(30),
            // 1-based genomic range: 30 bases at [p, p+30).
            vec![Exon::with_genomic(1, 1, 30, p, p + 29)],
            Some("chr_cdsp".to_string()),
            Some(p),
            Some(p + 29),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );

        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("chr_cdsp", genomic);
        provider.add_transcript(transcript);
        provider
    }

    #[test]
    fn cds_proper_homopolymer_dup_still_gated() {
        // CDS-proper A[5] at c.4..c.8. Dup of "AA" at c.4_5 must
        // emit `c.8_9insAA` (codon-frame gate fires because
        // unit_len=1 in a coding context). Pins existing behavior
        // identical to `tests/dup_shift_matrix.rs::cds_minus::
        // multi_base_dup_of_multiple_units::case_1` (modulo strand)
        // — locks in that the B4-remaining fix does NOT regress
        // CDS-proper gating.
        let p = make_provider();
        assert_canonical_round_trip(p, "NM_CDSP.1:c.4_5dup", "NM_CDSP.1:c.8_9insAA");
    }

    #[test]
    fn cds_proper_homopolymer_del_still_gated() {
        // CDS-proper A[5] at c.4..c.8. Del of "AA" at c.4_5 — gate
        // fires, `deletion_to_repeat` returns None, falls back to
        // plain `del` at the post-shift position. The shift lands
        // on the 3'-most 2-base window inside the tract.
        let p = make_provider();
        assert_canonical_round_trip(p, "NM_CDSP.1:c.4_5del", "NM_CDSP.1:c.7_8del");
    }
}
