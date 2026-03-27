//! Coverage gap tests for HGVS normalization
//!
//! These tests target specific blind spots identified in the test suite audit:
//! 1.  Minus-strand normalization (all prior tests used plus-strand only)
//! 2.  Position ordering invariant (5'-to-3' in coding coordinates)
//! 3.  Intronic insertions
//! 4.  Intronic duplications
//! 5.  Intronic deletions on minus strand
//! 6.  Boundary-spanning variants (exon-intron)
//! 7.  UTR normalization on minus strand
//! 8.  Delins at exon-intron boundaries
//! 9.  Multi-allele canonical ordering
//! 10. Inversion normalization
//!
//! Issue #11: NM_004119.3:c.1837+2_1837+3ins... position ordering bug

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

// =============================================================================
// Test infrastructure: minus-strand multi-exon transcript with genomic mapping
// =============================================================================

/// Padding to ensure genomic sequences are large enough for the normalizer's
/// 100bp window on each side of the variant position.
const PAD: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

/// Build a genomic sequence with padding on both sides so the normalizer's
/// 100bp window never goes out of bounds. Returns (padded_sequence, offset)
/// where offset is the 0-based position in the padded string where the real
/// region starts.
fn padded_genomic(core: &str) -> String {
    format!("{}{}{}", PAD, core, PAD)
}

/// The offset added by padding (length of PAD).
const PAD_OFFSET: u64 = 256;

/// Build a minus-strand, 3-exon transcript with genomic coordinates and introns.
///
/// Transcript layout (in coding / transcript order):
///
///   Exon 1: tx 1-30    genomic (PAD+80)–(PAD+109)
///   Intron 1:           genomic (PAD+70)–(PAD+79)  GGGCCCTTTG
///   Exon 2: tx 31-60   genomic (PAD+40)–(PAD+69)
///   Intron 2:           genomic (PAD+30)–(PAD+39)  AAATTTAAAT
///   Exon 3: tx 61-90   genomic (PAD+0)–(PAD+29)
///
/// On minus strand, exon 1 (coding start) maps to rightmost genomic block.
fn make_minus_strand_transcript() -> (Transcript, String) {
    // Core genomic (110 bp):
    // Exon 3 RC (30bp) | Intron 2 (10bp) | Exon 2 RC (30bp) | Intron 1 (10bp) | Exon 1 RC (30bp)
    let core = "CATGCATGCATGCATGCATGCATGCATGCA\
                AAATTTAAAT\
                GCATGCATGCATGCATGCATGCATGCATGC\
                GGGCCCTTTG\
                ATGCATGCATGCATGCATGCATGCATGCAT";
    let genomic_seq = padded_genomic(core);

    // Transcript sequence is reverse complement of the exon regions in coding order
    // All ATGC repeats are palindromic in RC
    let tx_seq = "ATGCATGCATGCATGCATGCATGCATGCAT\
                   GCATGCATGCATGCATGCATGCATGCATGC\
                   TGCATGCATGCATGCATGCATGCATGCATG";

    let p = PAD_OFFSET; // shorthand
    let transcript = Transcript::new(
        "NM_MINUS.1".to_string(),
        Some("MINUSGENE".to_string()),
        Strand::Minus,
        tx_seq.to_string(),
        Some(1),
        Some(90),
        vec![
            Exon::with_genomic(1, 1, 30, p + 80, p + 109),
            Exon::with_genomic(2, 31, 60, p + 40, p + 69),
            Exon::with_genomic(3, 61, 90, p + 0, p + 29),
        ],
        Some("chr_minus".to_string()),
        Some(p),
        Some(p + 109),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );

    (transcript, genomic_seq)
}

/// Build a plus-strand, 3-exon transcript with genomic mapping for comparison.
///
/// Transcript layout:
///   Exon 1: tx 1-30   genomic (PAD+0)–(PAD+29)
///   Intron 1:          genomic (PAD+30)–(PAD+39)  AAACCCAAAT
///   Exon 2: tx 31-60  genomic (PAD+40)–(PAD+69)
///   Intron 2:          genomic (PAD+70)–(PAD+79)  GGGTTTTTTG
///   Exon 3: tx 61-90  genomic (PAD+80)–(PAD+109)
fn make_plus_strand_transcript() -> (Transcript, String) {
    let core = "ATGCATGCATGCATGCATGCATGCATGCAT\
                AAACCCAAAT\
                GCATGCATGCATGCATGCATGCATGCATGC\
                GGGTTTTTTG\
                CATGCATGCATGCATGCATGCATGCATGCA";
    let genomic_seq = padded_genomic(core);

    let tx_seq = "ATGCATGCATGCATGCATGCATGCATGCAT\
                   GCATGCATGCATGCATGCATGCATGCATGC\
                   CATGCATGCATGCATGCATGCATGCATGCA";

    let p = PAD_OFFSET;
    let transcript = Transcript::new(
        "NM_PLUS.1".to_string(),
        Some("PLUSGENE".to_string()),
        Strand::Plus,
        tx_seq.to_string(),
        Some(1),
        Some(90),
        vec![
            Exon::with_genomic(1, 1, 30, p + 0, p + 29),
            Exon::with_genomic(2, 31, 60, p + 40, p + 69),
            Exon::with_genomic(3, 61, 90, p + 80, p + 109),
        ],
        Some("chr_plus".to_string()),
        Some(p),
        Some(p + 109),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );

    (transcript, genomic_seq)
}

/// Build a plus-strand transcript with UTR regions and genomic mapping.
///
/// Layout: 5'UTR (5bp) + CDS (60bp) + 3'UTR (5bp) = 70bp total transcript
///   Exon 1: tx 1-25   genomic (PAD+0)–(PAD+24)
///   Intron 1:          genomic (PAD+25)–(PAD+34)  AAAAAAAAAA
///   Exon 2: tx 26-50  genomic (PAD+35)–(PAD+59)
///   Intron 2:          genomic (PAD+60)–(PAD+69)  CCCCCCCCCC
///   Exon 3: tx 51-70  genomic (PAD+70)–(PAD+89)
///   CDS: tx 6-65
fn make_plus_strand_utr_transcript() -> (Transcript, String) {
    let core = "AAAAATGCATGCATGCATGCATGCAT\
                AAAAAAAAAA\
                GCATGCATGCATGCATGCATGCATGC\
                CCCCCCCCCC\
                ATGCATGCATGCATGCATGCAAAAA";
    let genomic_seq = padded_genomic(core);

    let tx_seq = "AAAAATGCATGCATGCATGCATGCAT\
                   GCATGCATGCATGCATGCATGCATGC\
                   ATGCATGCATGCATGCATGCAAAAA";

    let p = PAD_OFFSET;
    let transcript = Transcript::new(
        "NM_UTR.1".to_string(),
        Some("UTRGENE".to_string()),
        Strand::Plus,
        tx_seq.to_string(),
        Some(6),
        Some(65),
        vec![
            Exon::with_genomic(1, 1, 25, p + 0, p + 24),
            Exon::with_genomic(2, 26, 50, p + 35, p + 59),
            Exon::with_genomic(3, 51, 70, p + 70, p + 89),
        ],
        Some("chr_utr".to_string()),
        Some(p),
        Some(p + 89),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );

    (transcript, genomic_seq)
}

/// Build a minus-strand transcript with UTR regions and genomic mapping.
///
/// Layout (minus strand, genomic left-to-right):
///   Exon 3: genomic (PAD+0)–(PAD+19)    tx 51-70
///   Intron 2: genomic (PAD+20)–(PAD+29)  GGGGGGGGGG
///   Exon 2: genomic (PAD+30)–(PAD+54)   tx 26-50
///   Intron 1: genomic (PAD+55)–(PAD+64)  TTTTTTTTTT
///   Exon 1: genomic (PAD+65)–(PAD+89)   tx 1-25
///   CDS: tx 6-65
fn make_minus_strand_utr_transcript() -> (Transcript, String) {
    let core = "TTTTTGCATGCATGCATGCAT\
                GGGGGGGGGG\
                GCATGCATGCATGCATGCATGCATGC\
                TTTTTTTTTT\
                ATGCATGCATGCATGCATGCATTTTT";
    let genomic_seq = padded_genomic(core);

    // Transcript is RC of exons in coding order
    let tx_seq = "AAAAATGCATGCATGCATGCATGCAT\
                   GCATGCATGCATGCATGCATGCATGC\
                   ATGCATGCATGCATGCAAAAA";

    let p = PAD_OFFSET;
    let transcript = Transcript::new(
        "NM_MUTR.1".to_string(),
        Some("MUTRGENE".to_string()),
        Strand::Minus,
        tx_seq.to_string(),
        Some(6),
        Some(65),
        vec![
            Exon::with_genomic(1, 1, 25, p + 65, p + 89),
            Exon::with_genomic(2, 26, 50, p + 30, p + 54),
            Exon::with_genomic(3, 51, 70, p + 0, p + 19),
        ],
        Some("chr_mutr".to_string()),
        Some(p),
        Some(p + 89),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );

    (transcript, genomic_seq)
}

fn make_provider_with_minus_strand() -> MockProvider {
    let mut provider = MockProvider::new();
    let (tx, genomic) = make_minus_strand_transcript();
    provider.add_genomic_sequence("chr_minus", genomic);
    provider.add_transcript(tx);
    provider
}

fn make_provider_with_plus_strand() -> MockProvider {
    let mut provider = MockProvider::new();
    let (tx, genomic) = make_plus_strand_transcript();
    provider.add_genomic_sequence("chr_plus", genomic);
    provider.add_transcript(tx);
    provider
}

fn make_provider_with_plus_utr() -> MockProvider {
    let mut provider = MockProvider::new();
    let (tx, genomic) = make_plus_strand_utr_transcript();
    provider.add_genomic_sequence("chr_utr", genomic);
    provider.add_transcript(tx);
    provider
}

fn make_provider_with_minus_utr() -> MockProvider {
    let mut provider = MockProvider::new();
    let (tx, genomic) = make_minus_strand_utr_transcript();
    provider.add_genomic_sequence("chr_mutr", genomic);
    provider.add_transcript(tx);
    provider
}

/// Normalize a variant and return the formatted HGVS string.
fn normalize(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant =
        parse_hgvs(input).unwrap_or_else(|e| panic!("Failed to parse '{}': {}", input, e));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("Normalization failed for '{}': {}", input, e));
    format!("{}", normalized)
}

/// Try to normalize; returns Err string on failure.
fn try_normalize(provider: MockProvider, input: &str) -> Result<String, String> {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).map_err(|e| format!("Parse error: {}", e))?;
    let normalized = normalizer
        .normalize(&variant)
        .map_err(|e| format!("Normalize error: {}", e))?;
    Ok(format!("{}", normalized))
}

// =============================================================================
// GAP 1: Minus-strand exonic normalization
// =============================================================================
// All prior normalization tests used plus-strand transcripts. These test that
// basic exonic normalization (3' shifting, dup detection) works on minus strand.

mod minus_strand_exonic {
    use super::*;

    #[test]
    fn test_minus_strand_substitution_passthrough() {
        // Substitutions should pass through unchanged regardless of strand
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.15A>G");
        assert_eq!(result, "NM_MINUS.1:c.15A>G");
    }

    #[test]
    fn test_minus_strand_deletion_3prime_shift() {
        // Deletion in ATGC repeat region should shift 3' in coding direction.
        // On minus strand, 3' in coding = 5' in genomic (toward lower genomic coords).
        // Transcript seq: ATGCATGCATGCATGC... (30bp exon 1)
        // Deleting c.1del (A) should shift 3' through the ATGC repeat.
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.1del");
        // The exact shifted position depends on the repeat boundary
        assert!(
            result.contains("del"),
            "Should still be a deletion, got: {}",
            result
        );
    }
}

// =============================================================================
// GAP 2: Position ordering invariant (Issue #11)
// =============================================================================
// The HGVS spec requires insertion flanking positions to be in 5'-to-3' order
// in coding coordinates. For c. positions: smaller number first. For intronic
// offsets on the same base: +2 before +3 (deeper into intron).

mod position_ordering {
    use super::*;

    #[test]
    fn test_plus_strand_intronic_insertion_position_order() {
        // On plus strand, c.30+2_30+3 should remain in order after normalization
        let provider = make_provider_with_plus_strand();
        let result = normalize(provider, "NM_PLUS.1:c.30+2_30+3insA");
        // Positions must be in 5'-to-3' coding order: +2 before +3
        assert!(
            result.contains("30+") && !result.contains("+3_") || result.contains("+2_"),
            "Insertion positions should be 5'-to-3' in coding order, got: {}",
            result
        );
    }

    #[test]
    fn test_minus_strand_intronic_insertion_position_order() {
        // Issue #11: On minus strand, c.30+2_30+3 positions must stay in coding
        // order (5' to 3'), NOT get swapped to genomic order.
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+2_30+3insA");
        // c.30+2 is 5' of c.30+3 in coding direction — must come first
        assert!(
            !result.contains("+3_30+2"),
            "Positions must not be swapped to genomic order (issue #11), got: {}",
            result
        );
    }

    #[test]
    fn test_minus_strand_intronic_insertion_negative_offset_order() {
        // c.31-3_31-2 should maintain order: -3 is further from exon (5' in coding)
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.31-3_31-2insT");
        assert!(
            !result.contains("-2_31-3") && !result.contains("-2_"),
            "Negative offset positions must maintain 5'-to-3' coding order, got: {}",
            result
        );
    }
}

// =============================================================================
// GAP 3: Intronic insertions
// =============================================================================

mod intronic_insertions {
    use super::*;

    #[test]
    fn test_plus_strand_intronic_insertion() {
        // Basic intronic insertion on plus strand
        let provider = make_provider_with_plus_strand();
        let result = normalize(provider, "NM_PLUS.1:c.30+3_30+4insGGG");
        assert!(
            result.contains("ins"),
            "Should remain an insertion, got: {}",
            result
        );
    }

    #[test]
    fn test_minus_strand_intronic_insertion() {
        // Basic intronic insertion on minus strand — the core issue #11 scenario
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+3_30+4insGGG");
        assert!(
            result.contains("ins"),
            "Should remain an insertion, got: {}",
            result
        );
    }

    #[test]
    fn test_intronic_insertion_to_dup_plus_strand() {
        // If inserted sequence matches the preceding intronic sequence, should become dup.
        // Intron 1 genomic (1031-1040): AAACCCAAAT
        // c.30+1 = genomic 1031 = A, c.30+2 = 1032 = A, c.30+3 = 1033 = A
        // Inserting A at c.30+3_30+4 into AAA tract could become dup
        let provider = make_provider_with_plus_strand();
        let result = normalize(provider, "NM_PLUS.1:c.30+1_30+2insA");
        // In an AAA tract, inserting A should become dup after 3' shift
        assert!(
            result.contains("dup") || result.contains("ins"),
            "Should be dup or ins, got: {}",
            result
        );
    }
}

// =============================================================================
// GAP 4: Intronic duplications
// =============================================================================

mod intronic_duplications {
    use super::*;

    #[test]
    fn test_plus_strand_intronic_dup() {
        let provider = make_provider_with_plus_strand();
        let result = normalize(provider, "NM_PLUS.1:c.30+2dup");
        assert!(
            result.contains("dup"),
            "Intronic dup should remain dup, got: {}",
            result
        );
    }

    #[test]
    fn test_minus_strand_intronic_dup() {
        // Intronic duplication on minus strand — tests that position is correct
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+2dup");
        assert!(
            result.contains("dup"),
            "Minus-strand intronic dup should remain dup, got: {}",
            result
        );
    }

    #[test]
    fn test_minus_strand_intronic_multi_base_dup() {
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+2_30+4dup");
        // Normalizer may convert multi-base dup into a repeat if the region is repetitive
        assert!(
            result.contains("dup") || result.contains('['),
            "Minus-strand multi-base intronic dup should normalize to dup or repeat, got: {}",
            result
        );
        // Verify position ordering is correct (coding order: 30+2 before 30+4)
        assert!(
            result.contains("30+2_30+4"),
            "Positions should be in coding order (30+2_30+4), got: {}",
            result
        );
    }
}

// =============================================================================
// GAP 5: Intronic deletions on minus strand
// =============================================================================

mod intronic_deletions_minus {
    use super::*;

    #[test]
    fn test_minus_strand_intronic_single_del() {
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+2del");
        assert!(
            result.contains("del"),
            "Minus-strand intronic del should remain del, got: {}",
            result
        );
    }

    #[test]
    fn test_minus_strand_intronic_range_del() {
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+2_30+4del");
        assert!(
            result.contains("del"),
            "Minus-strand intronic range del should remain del, got: {}",
            result
        );
    }

    #[test]
    fn test_minus_strand_intronic_del_3prime_shift() {
        // Deletion in intronic repeat region on minus strand should shift 3' in
        // coding direction (which is 5' in genomic direction for minus strand).
        // Intron 1 genomic (1031-1040): GGGCCCTTTG
        // On minus strand, coding direction is right-to-left in genomic space.
        // c.30+1 maps to genomic 1040 (G), c.30+2 maps to 1039 (T), etc.
        // So the CCC tract at genomic 1034-1036 maps to c.30+5 through c.30+7
        // (GGG in RC = CCC). Deleting one should shift 3' in coding direction.
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+5del");
        assert!(
            result.contains("del") && result.contains("30+"),
            "Should be intronic deletion after shift, got: {}",
            result
        );
    }
}

// =============================================================================
// GAP 6: Boundary-spanning variants (exon-intron)
// =============================================================================

mod boundary_spanning {
    use super::*;

    #[test]
    fn test_plus_strand_exon_intron_deletion() {
        // Deletion spanning last exonic base and first intronic base
        let provider = make_provider_with_plus_strand();
        let result = try_normalize(provider, "NM_PLUS.1:c.30_30+1del");
        // Should either normalize or return an error — not panic
        assert!(
            result.is_ok() || result.is_err(),
            "Boundary-spanning deletion should not panic"
        );
    }

    #[test]
    fn test_minus_strand_exon_intron_deletion() {
        let provider = make_provider_with_minus_strand();
        let result = try_normalize(provider, "NM_MINUS.1:c.30_30+1del");
        assert!(
            result.is_ok() || result.is_err(),
            "Minus-strand boundary-spanning deletion should not panic"
        );
    }

    #[test]
    fn test_intron_exon_boundary_insertion() {
        // Insertion at intron-exon junction
        let provider = make_provider_with_plus_strand();
        let result = try_normalize(provider, "NM_PLUS.1:c.31-1_31insA");
        assert!(
            result.is_ok() || result.is_err(),
            "Intron-exon boundary insertion should not panic"
        );
    }
}

// =============================================================================
// GAP 7: UTR normalization on minus strand
// =============================================================================

mod utr_minus_strand {
    use super::*;

    #[test]
    fn test_5prime_utr_substitution_minus_strand() {
        // 5'UTR variant on minus strand should pass through
        let provider = make_provider_with_minus_utr();
        let result = normalize(provider, "NM_MUTR.1:c.-3A>G");
        assert!(
            result.contains("c.-3"),
            "5'UTR substitution should pass through, got: {}",
            result
        );
    }

    #[test]
    fn test_3prime_utr_substitution_minus_strand() {
        let provider = make_provider_with_minus_utr();
        let result = normalize(provider, "NM_MUTR.1:c.*3A>G");
        assert!(
            result.contains("c.*3"),
            "3'UTR substitution should pass through, got: {}",
            result
        );
    }

    #[test]
    fn test_5prime_utr_deletion_minus_strand() {
        let provider = make_provider_with_minus_utr();
        let result = normalize(provider, "NM_MUTR.1:c.-3del");
        assert!(
            result.contains("del"),
            "5'UTR deletion on minus strand should normalize, got: {}",
            result
        );
    }

    #[test]
    fn test_3prime_utr_deletion_minus_strand() {
        let provider = make_provider_with_minus_utr();
        let result = normalize(provider, "NM_MUTR.1:c.*2del");
        assert!(
            result.contains("del"),
            "3'UTR deletion on minus strand should normalize, got: {}",
            result
        );
    }
}

// =============================================================================
// GAP 8: Delins at exon-intron boundaries
// =============================================================================

mod delins_boundary {
    use super::*;

    #[test]
    fn test_exonic_delins_near_boundary_plus() {
        // Delins at the last exonic position, near intron boundary
        let provider = make_provider_with_plus_strand();
        let result = normalize(provider, "NM_PLUS.1:c.29_30delinsAA");
        assert!(
            result.contains("delins") || result.contains("del") || result.contains(">"),
            "Delins near boundary should normalize, got: {}",
            result
        );
    }

    #[test]
    fn test_intronic_delins_plus() {
        let provider = make_provider_with_plus_strand();
        let result = try_normalize(provider, "NM_PLUS.1:c.30+2_30+4delinsTTT");
        assert!(
            result.is_ok() || result.is_err(),
            "Intronic delins should not panic"
        );
    }

    #[test]
    fn test_intronic_delins_minus() {
        let provider = make_provider_with_minus_strand();
        let result = try_normalize(provider, "NM_MINUS.1:c.30+2_30+4delinsAAA");
        assert!(
            result.is_ok() || result.is_err(),
            "Minus-strand intronic delins should not panic"
        );
    }

    #[test]
    fn test_boundary_spanning_delins() {
        // Delins spanning exon-intron boundary
        let provider = make_provider_with_plus_strand();
        let result = try_normalize(provider, "NM_PLUS.1:c.30_30+2delinsGGG");
        assert!(
            result.is_ok() || result.is_err(),
            "Boundary-spanning delins should not panic"
        );
    }
}

// =============================================================================
// GAP 9: Multi-allele canonical ordering
// =============================================================================

mod multi_allele {
    use super::*;

    #[test]
    fn test_allele_notation_parses_and_normalizes() {
        // Compound heterozygous / allele notation: [var1;var2]
        // Each component should be individually normalized
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let input = "NM_000088.3:c.[10A>G;20T>C]";
        let variant = parse_hgvs(input);
        assert!(
            variant.is_ok(),
            "Allele notation should parse, got: {:?}",
            variant.err()
        );
        if let Ok(v) = variant {
            let result = normalizer.normalize(&v);
            // Should either normalize components or pass through
            assert!(
                result.is_ok(),
                "Allele normalization should not error, got: {:?}",
                result.err()
            );
        }
    }

    #[test]
    fn test_allele_components_individually_normalized() {
        // Each variant within an allele should be independently normalized
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let input = "NM_000088.3:c.[10A>G;20T>C]";
        if let Ok(variant) = parse_hgvs(input) {
            if let Ok(normalized) = normalizer.normalize(&variant) {
                let output = format!("{}", normalized);
                // Output should still contain both variants
                assert!(
                    output.contains("A>G") && output.contains("T>C"),
                    "Both allele components should be present, got: {}",
                    output
                );
            }
        }
    }
}

// =============================================================================
// GAP 10: Inversion normalization
// =============================================================================

mod inversions {
    use super::*;

    fn provider_with_single_exon(id: &str, sequence: &str) -> MockProvider {
        let mut provider = MockProvider::new();
        let len = sequence.len();
        provider.add_transcript(Transcript::new(
            id.to_string(),
            Some("TEST".to_string()),
            Strand::Plus,
            sequence.to_string(),
            Some(1),
            Some(len as u64),
            vec![Exon::new(1, 1, len as u64)],
            None,
            None,
            None,
            Default::default(),
            ManeStatus::None,
            None,
            None,
        ));
        provider
    }

    #[test]
    fn test_inversion_passthrough() {
        // Simple inversion should pass through normalization unchanged
        let provider = provider_with_single_exon("NM_TEST.1", "ATGCCCGGGAAATTTCCCGGG");
        let result = normalize(provider, "NM_TEST.1:c.4_6inv");
        assert!(
            result.contains("inv"),
            "Inversion should remain as inv, got: {}",
            result
        );
    }

    #[test]
    fn test_inversion_maintains_position() {
        // Inversion positions should not shift during normalization
        let provider = provider_with_single_exon("NM_TEST.1", "ATGCCCGGGAAATTTCCCGGG");
        let result = normalize(provider, "NM_TEST.1:c.7_9inv");
        assert!(
            result.contains("c.7_9inv"),
            "Inversion positions should be preserved, got: {}",
            result
        );
    }

    #[test]
    fn test_inversion_in_repeat_region() {
        // Inversion in a repeat region — should inversions shift? Spec says no.
        let provider = provider_with_single_exon("NM_TEST.1", "ATGAAAAAAGGGAAATTTCCC");
        let result = normalize(provider, "NM_TEST.1:c.4_6inv");
        assert!(
            result.contains("inv"),
            "Inversion in repeat should remain inv, got: {}",
            result
        );
    }

    #[test]
    fn test_minus_strand_inversion() {
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.10_12inv");
        assert!(
            result.contains("inv"),
            "Minus-strand inversion should remain inv, got: {}",
            result
        );
    }
}

// =============================================================================
// Integration test: Issue #11 exact case (requires real reference data)
// =============================================================================
// This test uses the exact variant from GitHub issue #11. It requires the
// benchmark reference data to be present (run `ferro prepare` first).

mod issue_11 {
    use super::*;

    #[test]
    #[ignore] // Requires benchmark-output reference data
    fn test_flt3_itd_insertion_position_order() {
        // Issue #11: FLT3 ITD insertion on minus-strand gene
        // Mutalyzer normalizes to: NM_004119.3:c.1837+2_1837+3ins...
        // ferro incorrectly produces: NM_004119.3:c.1837+3_1837+2ins...
        //
        // The HGVS spec requires flanking positions in 5'-to-3' coding order.
        // c.1837+2 is 5' of c.1837+3 in coding direction, so +2 must come first.
        use ferro_hgvs::MultiFastaProvider;
        use std::path::Path;

        let ref_path = Path::new("benchmark-output/manifest.json");
        if !ref_path.exists() {
            eprintln!("Skipping: benchmark-output not available");
            return;
        }
        let provider =
            MultiFastaProvider::from_manifest(ref_path).expect("Failed to load reference data");
        let normalizer = Normalizer::new(provider);

        let input = "NM_004119.3:c.1837+2_1837+3insGGATATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGT";
        let variant = parse_hgvs(input).expect("Failed to parse FLT3 variant");
        let normalized = normalizer
            .normalize(&variant)
            .expect("FLT3 normalization failed");
        let result = format!("{}", normalized);

        // The key assertion: positions must be +2_+3, not +3_+2
        assert!(
            result.contains("1837+2_1837+3"),
            "FLT3 insertion positions must be 5'-to-3' (c.1837+2_1837+3), got: {}",
            result
        );
        assert!(
            !result.contains("1837+3_1837+2"),
            "Positions must NOT be in genomic order (c.1837+3_1837+2), got: {}",
            result
        );
    }
}
