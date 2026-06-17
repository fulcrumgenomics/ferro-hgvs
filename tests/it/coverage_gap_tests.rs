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
//!
//! Issue #94: every gap test that previously asserted only `contains(...)`
//! or `is_ok() || is_err()` is locked with `assert_eq!` against the
//! exact normalizer output. Two root causes have since been fixed:
//! #98 (minus-strand intronic ref-base orientation) and #97 (5'UTR
//! CDS↔tx off-by-one that collapsed UTR del positions to `c.?`); the
//! affected tests now lock spec-correct outputs (full-tract position
//! with the transcript-view repeat unit, per the HGVS DNA repeat
//! spec; UTR positions resolved per the HGVS no-c.0 numbering rule).
//! Boundary-spanning panic-canary tests are intentionally left as
//! `is_ok() || is_err()` per #94 scope.

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
/// Genomic coordinates are 1-based (matching real cdot data): the core region's
/// first base is at 1-based position `p = PAD_OFFSET + 1`. Transcript layout (in
/// coding / transcript order):
///
///   Exon 1: tx 1-30    genomic (p+80)–(p+109)
///   Intron 1:           genomic (p+70)–(p+79)  GGGCCCTTTG
///   Exon 2: tx 31-60   genomic (p+40)–(p+69)
///   Intron 2:           genomic (p+30)–(p+39)  AAATTTAAAT
///   Exon 3: tx 61-90   genomic (p+0)–(p+29)
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

    let p = PAD_OFFSET + 1; // shorthand
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
            Exon::with_genomic(3, 61, 90, p, p + 29),
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
/// Genomic coordinates are 1-based; the core's first base is at `p = PAD_OFFSET + 1`.
/// Transcript layout:
///   Exon 1: tx 1-30   genomic (p+0)–(p+29)
///   Intron 1:          genomic (p+30)–(p+39)  AAACCCAAAT
///   Exon 2: tx 31-60  genomic (p+40)–(p+69)
///   Intron 2:          genomic (p+70)–(p+79)  GGGTTTTTTG
///   Exon 3: tx 61-90  genomic (p+80)–(p+109)
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

    let p = PAD_OFFSET + 1;
    let transcript = Transcript::new(
        "NM_PLUS.1".to_string(),
        Some("PLUSGENE".to_string()),
        Strand::Plus,
        tx_seq.to_string(),
        Some(1),
        Some(90),
        vec![
            Exon::with_genomic(1, 1, 30, p, p + 29),
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

/// Build a plus-strand, 3-exon transcript whose exon1 3' tail and intron1 5'
/// form a single homopolymer run that straddles the junction — the fixture for
/// #704 sub-problem A (purely-exonic indel that, under the 3' rule, must shift
/// out of the exon and into the following intron).
///
/// `cds_start == 1` so `c.N == n.N`; the same `n.30del` / `c.30del` exercise
/// both the `n.` and `c.` post-checks against one genomic layout.
///
/// Genomic coordinates are 1-based; the core's first base is at `p = PAD_OFFSET + 1`.
/// Layout:
///   Exon 1: tx 1-30   genomic (p+0)–(p+29)   ends `...AAAA` (tx 27-30)
///   Intron 1:          genomic (p+30)–(p+39)  `AAAGCGCGCG` (starts `AAA`)
///   Exon 2: tx 31-60  genomic (p+40)–(p+69)
///   Intron 2:          genomic (p+70)–(p+79)  `CTCTCTCTCT`
///   Exon 3: tx 61-90  genomic (p+80)–(p+109)
///
/// The junction run is genomic (p+26)–(p+32) = `AAAAAAA` (4 exonic + 3 intronic
/// A's). A single-base deletion/duplication at the exon's 3' edge (`n.30`) is
/// part of that run, so the 3' rule shifts it to the run's 3'-most base at
/// genomic p+32 = `n.30+3` (intron position 3).
fn make_edge_homopolymer_transcript() -> (Transcript, String) {
    // Exon 1 (30) ends in AAAA; intron 1 (10) starts AAA → 7-A junction run.
    let core = "GCGCGCGCGCGCGCGCGCGCGCGCGCAAAA\
                AAAGCGCGCG\
                TGTGTGTGTGTGTGTGTGTGTGTGTGTGTG\
                CTCTCTCTCT\
                ACACACACACACACACACACACACACACAC";
    let genomic_seq = padded_genomic(core);

    // Plus strand: transcript sequence is the spliced exons in order.
    let tx_seq = "GCGCGCGCGCGCGCGCGCGCGCGCGCAAAA\
                   TGTGTGTGTGTGTGTGTGTGTGTGTGTGTG\
                   ACACACACACACACACACACACACACACAC";

    let p = PAD_OFFSET + 1;
    let transcript = Transcript::new(
        "NM_EDGE.1".to_string(),
        Some("EDGEGENE".to_string()),
        Strand::Plus,
        tx_seq.to_string(),
        Some(1),
        Some(90),
        vec![
            Exon::with_genomic(1, 1, 30, p, p + 29),
            Exon::with_genomic(2, 31, 60, p + 40, p + 69),
            Exon::with_genomic(3, 61, 90, p + 80, p + 109),
        ],
        Some("chr_edge".to_string()),
        Some(p),
        Some(p + 109),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );

    (transcript, genomic_seq)
}

fn make_provider_with_edge() -> MockProvider {
    let mut provider = MockProvider::new();
    let (tx, genomic) = make_edge_homopolymer_transcript();
    provider.add_genomic_sequence("chr_edge", genomic);
    provider.add_transcript(tx);
    provider
}

/// Build a minus-strand transcript with UTR regions and genomic mapping.
///
/// Genomic coordinates are 1-based; the core's first base is at `p = PAD_OFFSET + 1`.
/// Layout (minus strand, genomic left-to-right):
///   Exon 3: genomic (p+0)–(p+19)    tx 51-70
///   Intron 2: genomic (p+20)–(p+29)  GGGGGGGGGG
///   Exon 2: genomic (p+30)–(p+54)   tx 26-50
///   Intron 1: genomic (p+55)–(p+64)  TTTTTTTTTT
///   Exon 1: genomic (p+65)–(p+89)   tx 1-25
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

    let p = PAD_OFFSET + 1;
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
            Exon::with_genomic(3, 51, 70, p, p + 19),
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
    fn test_minus_strand_deletion_no_shift() {
        // Single-A deletion in transcript-view ATGC tandem repeat.
        // The deleted A is not part of an A homopolymer (it sits between
        // C and T), so no equivalent representation exists; canonical form
        // is c.1del.
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.1del");
        assert_eq!(result, "NM_MINUS.1:c.1del");
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
        // Plus-strand intron 1 has an AAA tract at c.30+1..c.30+3.
        // ins A canonicalizes to dup at the 3'-most A in the run (#81 A1).
        let provider = make_provider_with_plus_strand();
        let result = normalize(provider, "NM_PLUS.1:c.30+2_30+3insA");
        assert_eq!(result, "NM_PLUS.1:c.30+3dup");
    }

    #[test]
    fn test_minus_strand_intronic_insertion_position_order() {
        // Issue #11: positions must stay in 5'-to-3' coding order (no
        // swap). Transcript-view bases at c.30+1..c.30+4 are C, A, A, A
        // — analogous to the plus-strand AAA tract — so ins A
        // canonicalizes to dup at the 3'-most A in the run (c.30+4) per
        // #81 A1, with the issue-#11 no-swap invariant preserved.
        // Issue #98 closed the minus-strand intronic ref-base
        // orientation gap that had prevented this canonicalization
        // from firing.
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+2_30+3insA");
        assert_eq!(result, "NM_MINUS.1:c.30+4dup");
    }

    #[test]
    fn test_minus_strand_intronic_insertion_negative_offset_order() {
        // Transcript-view bases at c.31-3..c.31-1 are C, C, C; inserted T
        // does not match the flank, so no equivalent shift exists. Locks
        // both the no-swap invariant and the no-shift behavior.
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.31-3_31-2insT");
        assert_eq!(result, "NM_MINUS.1:c.31-3_31-2insT");
    }
}

// =============================================================================
// GAP 3: Intronic insertions
// =============================================================================
// Restored from PR #93 / issue #94: the ins-shift matrix in
// `tests/ins_shift_matrix.rs` has zero intronic cells (its synthetic
// builder places the CDS over the whole transcript, no intron offsets),
// so these gap tests are still the only intronic-ins coverage we have.

mod intronic_insertions {
    use super::*;

    #[test]
    fn test_plus_strand_intronic_insertion() {
        // Plus-strand intron 1 genomic sequence is AAACCCAAAT at
        // c.30+1..c.30+10. Inserting GGG at c.30+3_30+4 sits between
        // c.30+3=A and c.30+4=C; ins[0]=G matches neither flank, so no
        // 3'-shift and no dup canonicalization — canonical form is the
        // unchanged input.
        let provider = make_provider_with_plus_strand();
        let result = normalize(provider, "NM_PLUS.1:c.30+3_30+4insGGG");
        assert_eq!(result, "NM_PLUS.1:c.30+3_30+4insGGG");
    }

    #[test]
    fn test_minus_strand_intronic_insertion() {
        // Issue #11 core scenario. Minus-strand intron 1 reads as
        // CAAAGGGCCC in transcript orientation at c.30+1..c.30+10 (RC
        // of genomic GGGCCCTTTG). Inserting GGG at c.30+3_30+4 sits
        // between c.30+3=A and c.30+4=A; the would-be alternative
        // representation at c.30+4_30+5 yields a different sequence
        // (AAA|GGG|GGG vs AA|GGG|AGGG), so no equivalent shift exists.
        // Positions also remain in coding (5'→3') order per #11.
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+3_30+4insGGG");
        assert_eq!(result, "NM_MINUS.1:c.30+3_30+4insGGG");
    }

    #[test]
    fn test_intronic_insertion_to_dup_plus_strand() {
        // Plus-strand intron 1: AAACCCAAAT. Inserting A at c.30+1_30+2
        // falls inside the leading AAA tract; ins A in a homopolymer
        // canonicalizes to dup at the 3'-most A in the run per #81 A1
        // (3'-shift walks 30+1_30+2 → 30+2_30+3 → 30+3_30+4, then
        // collapses to dup of the preceding base, c.30+3dup).
        let provider = make_provider_with_plus_strand();
        let result = normalize(provider, "NM_PLUS.1:c.30+1_30+2insA");
        assert_eq!(result, "NM_PLUS.1:c.30+3dup");
    }
}

// =============================================================================
// GAP 4: Intronic duplications
// =============================================================================

mod intronic_duplications {
    use super::*;

    #[test]
    fn test_plus_strand_intronic_dup() {
        // Plus-strand intron AAA tract at c.30+1..c.30+3; dup of A at +2
        // shifts 3' to the rightmost A in the run.
        let provider = make_provider_with_plus_strand();
        let result = normalize(provider, "NM_PLUS.1:c.30+2dup");
        assert_eq!(result, "NM_PLUS.1:c.30+3dup");
    }

    #[test]
    fn test_minus_strand_intronic_dup() {
        // Transcript-view AAA tract at c.30+2..c.30+4; dup of A at +2
        // shifts to c.30+4dup per #81 A6. Issue #98 closed the
        // minus-strand intronic ref-base orientation gap that had
        // prevented this shift from firing.
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+2dup");
        assert_eq!(result, "NM_MINUS.1:c.30+4dup");
    }

    #[test]
    fn test_minus_strand_intronic_multi_base_dup() {
        // Transcript-view AAA tract at c.30+2..c.30+4 duplicated to
        // AAAAAA. Per HGVS spec (repeated.md):
        //
        //   > This restriction only applies to the coding sequence,
        //   > which does not include the introns or the UTR sequence.
        //
        // The codon-frame `unit_len % 3 == 0` restriction does NOT apply
        // to intronic positions, so `normalize_intronic_cds` passes
        // `is_coding=false` and the canonical form is `A[6]` over the
        // reference-tract extent, not the gated `ins<literal>` fallback.
        // Issue #98 ensured the minus-strand orientation produces
        // transcript-view `A`s rather than genomic-view `T`s.
        // (B4-remaining, issue #209.)
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+2_30+4dup");
        assert_eq!(result, "NM_MINUS.1:c.30+2_30+4A[6]");
    }
}

// =============================================================================
// GAP 5: Intronic deletions on minus strand
// =============================================================================

mod intronic_deletions_minus {
    use super::*;

    #[test]
    fn test_minus_strand_intronic_single_del() {
        // Transcript-view AAA tract at c.30+2..c.30+4; single-A del at
        // +2 shifts to the 3'-most A (c.30+4) per #81 A5. Issue #98
        // closed the minus-strand intronic ref-base orientation gap
        // that had prevented this shift from firing.
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+2del");
        assert_eq!(result, "NM_MINUS.1:c.30+4del");
    }

    #[test]
    fn test_minus_strand_intronic_range_del() {
        // Transcript-view bases c.30+1..c.30+5 are C, A, A, A, G. Deleting
        // the entire AAA run (c.30+2..c.30+4) leaves C-G with no equivalent
        // representation; canonical form is unchanged.
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+2_30+4del");
        assert_eq!(result, "NM_MINUS.1:c.30+2_30+4del");
    }

    #[test]
    fn test_minus_strand_intronic_del_3prime_shift() {
        // Transcript-view GGG tract at c.30+5..c.30+7 (RC of genomic
        // CCC at PAD+73..PAD+75); single-G del at +5 shifts to the
        // 3'-most G in the run per #81 A5. For this 10-base intron
        // ferro renders the resulting position with acceptor-relative
        // notation (positions past the midpoint use `c.31-N`); c.30+7
        // and c.31-4 designate the same intronic base. Issue #98
        // closed the orientation gap that had prevented this shift
        // from firing.
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+5del");
        assert_eq!(result, "NM_MINUS.1:c.31-4del");
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

    #[test]
    fn test_minus_strand_boundary_spanning_delins_identity() {
        // Minus-strand boundary-spanning delins discriminator:
        // c.30 (last exonic base of exon 1, tx view) = T;
        // c.30+1 (first intronic base in tx-3' direction, tx view) = C
        // (RC of genomic 'G' at the high end of intron 1).
        // Inserting "TC" therefore matches the tx-view reference and PR #78
        // collapses it to identity.
        //
        // If `normalize_boundary_spanning_cds` failed to flip the genomic
        // window into transcript view before calling the delins-identity
        // rule, the rule would compare "TC" against the genomic-strand
        // bytes (G,A) and the insertion would not collapse — the result
        // would stay as `c.30_30+1delinsTC`. Asserting `c.30_30+1=` locks
        // the strand-aware path.
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30_30+1delinsTC");
        assert_eq!(result, "NM_MINUS.1:c.30_30+1=");
    }
}

// =============================================================================
// GAP 7: UTR normalization on minus strand
// =============================================================================

mod utr_minus_strand {
    use super::*;

    #[test]
    fn test_5prime_utr_substitution_minus_strand() {
        // Substitutions never shift; round-trip through normalize.
        let provider = make_provider_with_minus_utr();
        let result = normalize(provider, "NM_MUTR.1:c.-3A>G");
        assert_eq!(result, "NM_MUTR.1:c.-3A>G");
    }

    #[test]
    fn test_3prime_utr_substitution_minus_strand() {
        let provider = make_provider_with_minus_utr();
        let result = normalize(provider, "NM_MUTR.1:c.*3A>G");
        assert_eq!(result, "NM_MUTR.1:c.*3A>G");
    }

    #[test]
    fn test_5prime_utr_deletion_minus_strand() {
        // Single-base 5'UTR del on minus strand. The 5'UTR is the
        // homopolymer `AAAAA` at c.-5..c.-1 (transcript view); a
        // single-A del at c.-3 3'-shifts to the rightmost A in the
        // run, c.-1. Issue #97 closed the off-by-one in the CDS↔tx
        // 5'UTR coordinate mapping that had previously caused the
        // result to collapse to `c.?del`.
        let provider = make_provider_with_minus_utr();
        let result = normalize(provider, "NM_MUTR.1:c.-3del");
        assert_eq!(result, "NM_MUTR.1:c.-1del");
    }

    #[test]
    fn test_3prime_utr_deletion_minus_strand() {
        // Single-base 3'UTR del on minus strand. Locks current behavior;
        // whether further shift is expected depends on the UTR tract
        // shape — the fixture's 3'UTR is short and not a clean
        // homopolymer at this offset.
        let provider = make_provider_with_minus_utr();
        let result = normalize(provider, "NM_MUTR.1:c.*2del");
        assert_eq!(result, "NM_MUTR.1:c.*2del");
    }
}

// =============================================================================
// GAP 8: Delins at exon-intron boundaries
// =============================================================================

mod delins_boundary {
    use super::*;

    #[test]
    fn test_exonic_delins_near_boundary_plus_collapses_to_substitution() {
        // Delins at the last two exonic bases of exon 1 (c.29=A, c.30=T).
        // delinsAA shares the leading `A` with the reference, so under the
        // HGVS minimal-form rules (sub > del > inv > dup > ins) the canonical
        // form is the single-base substitution c.30T>A. The shared-affix
        // trimming added in #81 A2 makes this collapse fire on plus-strand
        // exonic boundaries the same way it fires elsewhere.
        let provider = make_provider_with_plus_strand();
        let result = normalize(provider, "NM_PLUS.1:c.29_30delinsAA");
        assert_eq!(result, "NM_PLUS.1:c.30T>A");
    }

    #[test]
    fn test_intronic_delins_plus() {
        // Plus-strand intronic bases at c.30+2..c.30+4 are A, A, C.
        // delinsTTT shares no leading/trailing match with AAC, so the
        // canonical form is unchanged.
        let provider = make_provider_with_plus_strand();
        let result = normalize(provider, "NM_PLUS.1:c.30+2_30+4delinsTTT");
        assert_eq!(result, "NM_PLUS.1:c.30+2_30+4delinsTTT");
    }

    #[test]
    fn test_intronic_delins_minus() {
        // Transcript-view bases at c.30+2..c.30+4 are A, A, A; insertion
        // is also AAA — this is an identity delins and PR #78 rewrites
        // it to identity (`c.30+2_30+4=`). Issue #98 closed the
        // minus-strand intronic ref-base orientation gap that had
        // prevented PR #78's identity check from succeeding.
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.30+2_30+4delinsAAA");
        assert_eq!(result, "NM_MINUS.1:c.30+2_30+4=");
    }

    #[test]
    fn test_boundary_spanning_delins() {
        // Delins spanning exon-intron boundary (c.30 is last exonic base,
        // c.30+1, +2 are intronic). True boundary-spanning case where
        // either canonical form or an explicit error is acceptable —
        // panic-canary only, per #94 scope.
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
        // Compound allele notation [var1;var2]. Both components are
        // simple substitutions at separated positions, which never shift,
        // so the round-trip is identity.
        let provider = MockProvider::with_test_data();
        let result = normalize(provider, "NM_000088.3:c.[10A>G;20T>C]");
        assert_eq!(result, "NM_000088.3:c.[10A>G;20T>C]");
    }

    #[test]
    fn test_allele_components_individually_normalized() {
        // Each variant within an allele is independently normalized.
        // The deletion c.34_36del (CTG) sits at the 5' end of a long CTG
        // tandem repeat and 3'-shifts to c.38_40del; the substitution
        // c.10A>G never shifts. Locking the mixed-shift outcome verifies
        // that components are normalized in isolation rather than the
        // bracket being treated as opaque.
        let provider = MockProvider::with_test_data();
        let result = normalize(provider, "NM_000088.3:c.[34_36del;10A>G]");
        assert_eq!(result, "NM_000088.3:c.[38_40del;10A>G]");
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
        // Inversions are not subject to 3'-shifting per HGVS spec; round-trip.
        let provider = provider_with_single_exon("NM_TEST.1", "ATGCCCGGGAAATTTCCCGGG");
        let result = normalize(provider, "NM_TEST.1:c.4_6inv");
        assert_eq!(result, "NM_TEST.1:c.4_6inv");
    }

    #[test]
    fn test_inversion_maintains_position() {
        let provider = provider_with_single_exon("NM_TEST.1", "ATGCCCGGGAAATTTCCCGGG");
        let result = normalize(provider, "NM_TEST.1:c.7_9inv");
        assert_eq!(result, "NM_TEST.1:c.7_9inv");
    }

    #[test]
    fn test_inversion_in_repeat_region() {
        // Inversion in a repeat region — spec says inversions don't shift,
        // even when the surrounding sequence is repetitive.
        let provider = provider_with_single_exon("NM_TEST.1", "ATGAAAAAAGGGAAATTTCCC");
        let result = normalize(provider, "NM_TEST.1:c.4_6inv");
        assert_eq!(result, "NM_TEST.1:c.4_6inv");
    }

    #[test]
    fn test_minus_strand_inversion() {
        let provider = make_provider_with_minus_strand();
        let result = normalize(provider, "NM_MINUS.1:c.10_12inv");
        assert_eq!(result, "NM_MINUS.1:c.10_12inv");
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

// =============================================================================
// #704: exon/intron 3' rule parity for n. (and r.) — boundary-spanning n.
//
// `NM_PLUS.1` has cds_start=1, so c.N == n.N: a boundary-spanning n. variant
// (one exonic + one intronic endpoint) must normalize through genomic space
// exactly as the c. analog does (#670), rather than erroring "not supported
// for n. coords". The c. analogs (c.30_30+1del, c.31-1_31insA) already
// normalize to themselves; the n. forms must too (not Err).
// =============================================================================

mod issue_704_nr_boundary {
    use super::*;

    #[test]
    fn n_boundary_spanning_deletion_normalizes() {
        // Exonic last base (n.30) + first intronic base (n.30+1): no equivalent
        // shifted position, so it normalizes to itself — but it must NOT error.
        let result = try_normalize(make_provider_with_plus_strand(), "NM_PLUS.1:n.30_30+1del");
        assert_eq!(result, Ok("NM_PLUS.1:n.30_30+1del".to_string()));
    }

    #[test]
    fn n_intron_exon_junction_insertion_normalizes() {
        let result = try_normalize(make_provider_with_plus_strand(), "NM_PLUS.1:n.31-1_31insA");
        assert_eq!(result, Ok("NM_PLUS.1:n.31-1_31insA".to_string()));
    }

    #[test]
    fn n_both_intronic_same_intron_still_shifts() {
        // Regression guard: the already-working both-intronic same-intron path
        // (routed to normalize_intronic_tx) must keep shifting (#670 sub-problem
        // B must not break it). c.30+1del -> c.30+3del; n. mirrors it.
        let result = normalize(make_provider_with_plus_strand(), "NM_PLUS.1:n.30+1del");
        assert_eq!(result, "NM_PLUS.1:n.30+3del");
    }

    #[test]
    fn n_multi_intron_span_routes_to_boundary_spanning() {
        // Both endpoints intronic but in DIFFERENT introns (intron1 -> intron2,
        // spanning exon2). The old `n.` path sent both-intronic variants to the
        // single-intron `normalize_intronic_tx` (no `same_intron` split); after
        // #704 the `same_intron_tx` split routes this to boundary-spanning. No
        // equivalent shift exists, so it normalizes to itself — matching the c.
        // analog — but must be handled, not mis-bounded to intron1.
        let result = try_normalize(make_provider_with_plus_strand(), "NM_PLUS.1:n.30+2_60+2del");
        assert_eq!(result, Ok("NM_PLUS.1:n.30+2_60+2del".to_string()));
    }
}

// =============================================================================
// #704 sub-problem A: purely-exonic exon→intron edge post-check for n. (and c.)
//
// On `NM_EDGE.1`, exon1's 3' tail (`...AAAA`) and intron1's 5' (`AAA...`) form a
// 7-A run spanning the junction. A purely-exonic single-base del/dup at the
// exon's 3' edge (`n.30` / `c.30`) is part of that run, so the 3' rule must
// shift it OUT of the exon and into the intron — landing at the run's 3'-most
// base, intron position 3 (`n.30+3` / `c.30+3`). The exon-confined shuffle only
// sees spliced bases and stops at the exon edge; the genomic-space post-check
// (mirror of #670's `normalize_cds` block) is what carries it across.
// =============================================================================

mod issue_704_nr_exon_to_intron {
    use super::*;

    #[test]
    fn c_exonic_edge_deletion_shifts_into_intron() {
        // Validates the genomic fixture against the already-shipped `c.`
        // post-check (#670): the c. analog must already cross into the intron.
        let result = normalize(make_provider_with_edge(), "NM_EDGE.1:c.30del");
        assert_eq!(result, "NM_EDGE.1:c.30+3del");
    }

    #[test]
    fn n_exonic_edge_deletion_shifts_into_intron() {
        // The `n.` parity case (#704 sub-problem A). Before the post-check is
        // ported to `normalize_tx` this stays exon-confined as `n.30del`.
        let result = normalize(make_provider_with_edge(), "NM_EDGE.1:n.30del");
        assert_eq!(result, "NM_EDGE.1:n.30+3del");
    }

    #[test]
    fn c_exonic_edge_duplication_shifts_into_intron() {
        let result = normalize(make_provider_with_edge(), "NM_EDGE.1:c.30dup");
        assert_eq!(result, "NM_EDGE.1:c.30+3dup");
    }

    #[test]
    fn n_exonic_edge_duplication_shifts_into_intron() {
        let result = normalize(make_provider_with_edge(), "NM_EDGE.1:n.30dup");
        assert_eq!(result, "NM_EDGE.1:n.30+3dup");
    }
}

// =============================================================================
// #704 increment r.: route intronic / boundary-spanning `r.` through the tx
// machinery, plus the sub-problem-A exon→intron edge post-check.
//
// `NM_PLUS.1` has cds_start=1 so `r.N == c.N == n.N`; each `r.` case mirrors the
// already-passing `n.` case. Pre-#704 every intronic/boundary `r.` returned
// `Err(IntronicVariant)`. The display carries the `r.` prefix and U/T rendering.
// =============================================================================

mod issue_704_r_intronic {
    use super::*;

    #[test]
    fn r_boundary_spanning_deletion_normalizes() {
        // One exonic + one intronic endpoint: no equivalent shift, normalizes to
        // itself — but must NOT error (was Err pre-#704). Mirror of the `n.` case.
        let result = try_normalize(make_provider_with_plus_strand(), "NM_PLUS.1:r.30_30+1del");
        assert_eq!(result, Ok("NM_PLUS.1:r.30_30+1del".to_string()));
    }

    #[test]
    fn r_intron_exon_junction_insertion_normalizes() {
        let result = try_normalize(make_provider_with_plus_strand(), "NM_PLUS.1:r.31-1_31insa");
        assert_eq!(result, Ok("NM_PLUS.1:r.31-1_31insa".to_string()));
    }

    #[test]
    fn r_both_intronic_same_intron_shifts() {
        // Both intronic in the same intron → normalize_intronic_tx; the homopolymer
        // shifts. Mirror of `n.30+1del -> n.30+3del`.
        let result = normalize(make_provider_with_plus_strand(), "NM_PLUS.1:r.30+1del");
        assert_eq!(result, "NM_PLUS.1:r.30+3del");
    }

    #[test]
    fn r_multi_intron_span_routes_to_boundary_spanning() {
        // Both intronic but in DIFFERENT introns → boundary-spanning (the
        // same_intron_tx split). No equivalent shift; normalizes to itself.
        let result = try_normalize(make_provider_with_plus_strand(), "NM_PLUS.1:r.30+2_60+2del");
        assert_eq!(result, Ok("NM_PLUS.1:r.30+2_60+2del".to_string()));
    }

    #[test]
    fn r_exonic_edge_deletion_shifts_into_intron() {
        // Sub-problem A for `r.`: a purely-exonic edge del must cross into the
        // intron. Mirror of `c./n.30del -> 30+3del` on the homopolymer fixture.
        let result = normalize(make_provider_with_edge(), "NM_EDGE.1:r.30del");
        assert_eq!(result, "NM_EDGE.1:r.30+3del");
    }

    #[test]
    fn r_exonic_edge_duplication_shifts_into_intron() {
        let result = normalize(make_provider_with_edge(), "NM_EDGE.1:r.30dup");
        assert_eq!(result, "NM_EDGE.1:r.30+3dup");
    }

    #[test]
    fn r_minus_strand_cds_relative_intronic_shifts() {
        // The shipped cases above all use NM_PLUS.1/NM_EDGE.1 with cds_start=1,
        // where r.N == n.N — so they never exercise the CDS-relative arm of
        // `rna_pos_to_txpos`/`txpos_to_rnapos`. NM_MUTR.1 is minus-strand with
        // cds_start=6, so r.N != n.N (r.20 == n.25, r.21 == n.26) and the anchor
        // base must be remapped through the CDS offset on the way to tx and back.
        //
        // Each `r.` result is the trusted part-1 `n.` result (same genomic
        // position) re-expressed in CDS-relative numbering with lowercase U/T
        // rendering: n.26-3 -> r.21-3, n.24_26-3A[8] -> r.19_21-3a[8]. This also
        // covers boundary-spanning (`r.20_20+1`), the genomic-space repeat
        // notation, and the minus-strand intronic orientation in one place.
        let p = make_provider_with_minus_utr;

        // same-intron shuffle (both intronic) on a minus strand, cds-relative
        assert_eq!(
            try_normalize(p(), "NM_MUTR.1:r.20+1del"),
            Ok("NM_MUTR.1:r.21-3del".to_string())
        );
        // boundary-spanning (one exonic + one intronic endpoint)
        assert_eq!(
            try_normalize(p(), "NM_MUTR.1:r.20_20+1del"),
            Ok("NM_MUTR.1:r.19_21-3a[8]".to_string())
        );
        // multi-base intronic span in the same intron
        assert_eq!(
            try_normalize(p(), "NM_MUTR.1:r.20+1_20+2del"),
            Ok("NM_MUTR.1:r.19_21-3a[8]".to_string())
        );

        // Cross-check the r.<->n. correspondence the assertions above encode:
        // r.20 == n.25 and r.21 == n.26, so the r. outputs must be the n. outputs
        // shifted by the same cds offset (5) — proving the CDS-relative round trip.
        assert_eq!(
            try_normalize(p(), "NM_MUTR.1:n.25+1del"),
            Ok("NM_MUTR.1:n.26-3del".to_string())
        );
        assert_eq!(
            try_normalize(p(), "NM_MUTR.1:n.25_25+1del"),
            Ok("NM_MUTR.1:n.24_26-3A[8]".to_string())
        );
    }

    #[test]
    fn r_intronic_edit_with_uracil_base() {
        // Cross-cut of #704 (intronic `r.` routing) and #736 (the `u`->`T`
        // normalization): an intronic `r.` delins whose inserted bases are
        // written with the RNA base `u`. The genomic-space engine must see DNA,
        // so without the #736 mapping the reverse-complement inversion would be
        // missed and this would stay `delinsuu`. `r.30+1_30+2` is the intron-1
        // `aa`; `uu` is its reverse complement, so the canonical form is `inv`,
        // matching the `n.`-with-`T` analog.
        assert_eq!(
            normalize(
                make_provider_with_plus_strand(),
                "NM_PLUS.1:r.30+1_30+2delinsuu"
            ),
            "NM_PLUS.1:r.30+1_30+2inv",
        );
        assert_eq!(
            normalize(
                make_provider_with_plus_strand(),
                "NM_PLUS.1:n.30+1_30+2delinsTT"
            ),
            "NM_PLUS.1:n.30+1_30+2inv",
        );
        // A plain intronic `r.` insertion carrying `u` must also normalize
        // (not error) and match its `n.`-with-`T` analog.
        assert_eq!(
            try_normalize(
                make_provider_with_plus_strand(),
                "NM_PLUS.1:r.30+1_30+2insu"
            ),
            Ok("NM_PLUS.1:r.30+1_30+2insu".to_string()),
        );
    }
}
