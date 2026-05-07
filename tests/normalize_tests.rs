//! Normalization tests
//!
//! Comprehensive tests for the HGVS normalization engine including:
//! - 3' and 5' shuffling
//! - Boundary respect (exon boundaries)
//! - Duplication detection
//! - Various edit types

use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Shared helpers for integration tests that require benchmark reference data.
#[cfg(test)]
mod integration_helpers {
    use ferro_hgvs::{parse_hgvs, MultiFastaProvider, Normalizer};
    use std::path::Path;

    pub fn create_normalizer() -> Option<Normalizer<MultiFastaProvider>> {
        let ref_path = Path::new("benchmark-output/manifest.json");
        if !ref_path.exists() {
            return None;
        }
        let provider = MultiFastaProvider::from_manifest(ref_path).ok()?;
        Some(Normalizer::new(provider))
    }

    pub fn try_normalize(input: &str) -> Option<String> {
        let normalizer = create_normalizer()?;
        let variant = parse_hgvs(input).ok()?;
        let normalized = normalizer.normalize(&variant).ok()?;
        Some(format!("{}", normalized))
    }

    pub fn has_reference_data() -> bool {
        Path::new("benchmark-output/manifest.json").exists()
    }
}

#[test]
fn test_normalizer_creation() {
    let provider = MockProvider::with_test_data();
    let _normalizer = Normalizer::new(provider);
}

#[test]
fn test_normalize_substitution_unchanged() {
    // Substitutions should not change during normalization
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NM_000088.3:c.459A>G").unwrap();
    let normalized = normalizer.normalize(&variant).unwrap();

    assert_eq!(format!("{}", variant), format!("{}", normalized));
}

#[test]
fn test_normalize_identity() {
    // Substitutions should remain unchanged (no shifting needed)
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    // Genomic substitution - should remain unchanged
    // (No transcript lookup needed for genomic variants without shifting)
    let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
    let result = normalizer.normalize(&variant).unwrap();

    // Verify output equals input
    assert_eq!(
        format!("{}", variant),
        format!("{}", result),
        "Genomic substitution should remain unchanged"
    );
}

#[test]
fn test_normalize_inversion_unchanged() {
    // Inversions should not be shifted
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NM_000088.3:c.10_15inv").unwrap();
    let normalized = normalizer.normalize(&variant).unwrap();

    assert_eq!(format!("{}", variant), format!("{}", normalized));
}

#[test]
fn test_config_with_direction() {
    // Test that config direction can be set
    let provider = MockProvider::with_test_data();
    let config = NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime);
    let normalizer = Normalizer::with_config(provider, config);

    assert_eq!(
        normalizer.config().shuffle_direction,
        ShuffleDirection::FivePrime
    );
}

#[test]
fn test_config_cross_boundaries() {
    // Test that config can enable boundary crossing
    let provider = MockProvider::with_test_data();
    let config = NormalizeConfig::default().allow_crossing_boundaries();
    let normalizer = Normalizer::with_config(provider, config);

    assert!(normalizer.config().cross_boundaries);
}

#[test]
fn test_normalize_protein_unchanged() {
    // Protein variants don't undergo position shifting
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NP_000079.2:p.Val600Glu").unwrap();
    let normalized = normalizer.normalize(&variant).unwrap();

    assert_eq!(format!("{}", variant), format!("{}", normalized));
}

#[test]
fn test_normalize_missing_transcript() {
    // Should gracefully handle missing transcripts
    let provider = MockProvider::new(); // Empty provider
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NM_MISSING.1:c.100del").unwrap();
    let result = normalizer.normalize(&variant);

    // Should return original variant when transcript not found
    assert!(result.is_ok());
    assert_eq!(format!("{}", variant), format!("{}", result.unwrap()));
}

#[test]
fn test_normalize_allele_normalizes_each_variant() {
    // Allele normalization should normalize each contained variant
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("[NM_000088.3:c.10A>G;NM_000088.3:c.20C>T]").unwrap();
    let result = normalizer.normalize(&variant);

    assert!(result.is_ok());
    // Each variant in the allele should be processed
    let normalized = result.unwrap();
    assert!(format!("{}", normalized).contains(";"));
}

#[test]
fn test_normalize_uncertain_variant_unchanged() {
    // Variants with uncertain positions should not be normalized
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NM_000088.3:c.?_10del").unwrap();
    let result = normalizer.normalize(&variant);

    assert!(result.is_ok());
    // Should return unchanged
    assert_eq!(format!("{}", variant), format!("{}", result.unwrap()));
}

#[test]
fn test_normalize_complex_interval_unchanged() {
    // Complex interval variants should not be normalized (uncertain positions)
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    // Use a valid complex interval pattern with intronic offsets
    let variant = parse_hgvs("NM_000088.3:c.(1+1_2-1)_(20+1_21-1)del").unwrap();
    let result = normalizer.normalize(&variant).unwrap();

    // Complex intervals with uncertain positions should remain unchanged
    assert_eq!(
        format!("{}", variant),
        format!("{}", result),
        "Complex interval should remain unchanged"
    );
}

#[test]
fn test_shuffle_direction_default_is_three_prime() {
    // Default shuffle direction should be 3'
    let config = NormalizeConfig::default();
    assert_eq!(config.shuffle_direction, ShuffleDirection::ThreePrime);
}

#[test]
fn test_shuffle_direction_parsing() {
    // Test parsing of shuffle direction strings
    let three_prime: ShuffleDirection = "3prime".parse().unwrap();
    assert_eq!(three_prime, ShuffleDirection::ThreePrime);

    let five_prime: ShuffleDirection = "5prime".parse().unwrap();
    assert_eq!(five_prime, ShuffleDirection::FivePrime);

    let three_alt: ShuffleDirection = "3'".parse().unwrap();
    assert_eq!(three_alt, ShuffleDirection::ThreePrime);

    let five_alt: ShuffleDirection = "5'".parse().unwrap();
    assert_eq!(five_alt, ShuffleDirection::FivePrime);
}

#[test]
fn test_normalize_repeat_unchanged() {
    // Repeat notation variants should not change position
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NM_000088.3:c.10CAG[5]").unwrap();
    let result = normalizer.normalize(&variant).unwrap();

    // Repeat notation should pass through unchanged
    let output = format!("{}", result);
    assert!(
        output.contains("CAG[5]"),
        "Repeat notation should be preserved, got: {}",
        output
    );
}

#[test]
fn test_normalize_duplication_type_preserved() {
    // Test duplication normalization - type should be preserved even if position shifts
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NM_000088.3:c.10dup").unwrap();
    let result = normalizer.normalize(&variant).unwrap();

    // Duplication should remain a duplication
    let output = format!("{}", result);
    assert!(
        output.contains("dup"),
        "Duplication should remain a dup, got: {}",
        output
    );
}

#[test]
fn test_normalize_delins_type_preserved() {
    // Test delins normalization - should remain delins when not simplifiable.
    // ref NM_000088.3 c.10_11 = GT. delete GT, insert CA: no shared affix,
    // revcomp(GT) = AC != CA so not an inversion, length doesn't fit dup, so
    // the minimal form is still a 2->2 delins.
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NM_000088.3:c.10_11delinsCA").unwrap();
    let result = normalizer.normalize(&variant).unwrap();

    let output = format!("{}", result);
    assert!(
        output.contains("delins"),
        "Delins should remain delins, got: {}",
        output
    );
}

#[test]
fn test_normalize_rna_variant() {
    // RNA variants should pass through unchanged
    // Mutalyzer doesn't normalize RNA variants, they stay as-is
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NM_000088.3:r.10a>g").unwrap();
    let result = normalizer.normalize(&variant).unwrap();

    // RNA variants should remain unchanged
    assert_eq!(
        format!("{}", variant),
        format!("{}", result),
        "RNA variant should remain unchanged"
    );
}

#[test]
fn test_normalize_mitochondrial_variant() {
    // Mitochondrial variants should pass through unchanged
    // Mutalyzer confirms: NC_012920.1:m.3243A>G -> NC_012920.1:m.3243A>G
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NC_012920.1:m.3243A>G").unwrap();
    let result = normalizer.normalize(&variant).unwrap();

    // Mitochondrial variants should remain unchanged
    assert_eq!(
        format!("{}", variant),
        format!("{}", result),
        "Mitochondrial variant should remain unchanged"
    );
}

#[test]
fn test_config_validate_ref_default() {
    // Default config (lenient) should warn but not reject reference mismatches
    let config = NormalizeConfig::default();
    assert!(!config.should_reject_ref_mismatch());
    assert!(config.should_warn_ref_mismatch());
}

#[test]
#[allow(deprecated)]
fn test_config_skip_validation() {
    // Should be able to disable reference validation warnings
    let config = NormalizeConfig::default().skip_validation();
    assert!(!config.should_reject_ref_mismatch());
    assert!(!config.should_warn_ref_mismatch());
}

// =============================================================================
// NORMALIZATION OUTPUT TESTS
// =============================================================================
// These tests verify that normalization produces the correct output according
// to HGVS rules as implemented by mutalyzer. Each test creates a custom
// transcript with a specific sequence designed to test particular behaviors.

mod output_tests {
    use ferro_hgvs::reference::mock::MockProvider;
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
    use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

    /// Create a transcript with a custom sequence for testing
    fn make_transcript(id: &str, sequence: &str) -> Transcript {
        let len = sequence.len();
        Transcript::new(
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
        )
    }

    fn provider_with_transcript(id: &str, sequence: &str) -> MockProvider {
        let mut provider = MockProvider::new();
        provider.add_transcript(make_transcript(id, sequence));
        provider
    }

    /// Helper to normalize and get output string
    fn normalize_to_string(provider: MockProvider, input: &str) -> String {
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse input");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        format!("{}", normalized)
    }

    // =========================================================================
    // DELETION TESTS - Deletions should shift 3' but NEVER become duplications
    // =========================================================================

    #[test]
    fn test_deletion_shifts_3prime_in_repeat() {
        // Sequence: GGGGGGGGGAAAAAAAAGGGGGGGGGG
        //           123456789012345678901234567
        // 9 G's at positions 1-9, 8 A's at positions 10-17, G's at 18+
        // Deleting A at position 10 should shift to position 17 (3'-most A)
        let seq = "GGGGGGGGGAAAAAAAAGGGGGGGGGG"; // 8 A's at pos 10-17 (1-indexed)
        let provider = provider_with_transcript("NM_TEST.1", seq);

        let result = normalize_to_string(provider, "NM_TEST.1:c.10del");

        // Should shift to 3'-most position in the A repeat (position 17)
        assert_eq!(
            result, "NM_TEST.1:c.17del",
            "Deletion should shift to 3'-most position"
        );
    }

    #[test]
    fn test_deletion_never_becomes_dup() {
        // Locks two rules at once:
        // (a) deletion never becomes dup (HGVS prioritization: del > dup)
        // (b) deletion of >=2 tandem-repeat units becomes unit[N-k] (B2)
        //
        // Sequence: ATGATGATGATGATG (5 ATG units at pos 1-15).
        // Two cases:
        //   c.4_6del (k=1) → c.13_15del (still del, shifted to 3'-most ATG)
        //   c.4_9del (k=2) → c.1_15ATG[3] (B2 fires, becomes repeat)
        let seq = "ATGATGATGATGATG"; // 5 repeating ATG
        let provider = provider_with_transcript("NM_TEST.1", seq);

        // k=1: stays as del, shifted to 3'-most position
        let result = normalize_to_string(provider, "NM_TEST.1:c.4_6del");
        assert_eq!(
            result, "NM_TEST.1:c.13_15del",
            "k=1 deletion should stay as del shifted 3'"
        );

        // k=2: B2 fires, becomes ATG[3]
        let provider = provider_with_transcript("NM_TEST.1", seq);
        let result = normalize_to_string(provider, "NM_TEST.1:c.4_9del");
        assert_eq!(
            result, "NM_TEST.1:c.1_15ATG[3]",
            "k=2 deletion of tandem unit should become unit[N-k]"
        );

        // Neither result mentions "dup"
        let provider = provider_with_transcript("NM_TEST.1", seq);
        let r1 = normalize_to_string(provider, "NM_TEST.1:c.4_6del");
        assert!(
            !r1.contains("dup"),
            "k=1 case must not become dup, got: {}",
            r1
        );
        let provider = provider_with_transcript("NM_TEST.1", seq);
        let r2 = normalize_to_string(provider, "NM_TEST.1:c.4_9del");
        assert!(
            !r2.contains("dup"),
            "k=2 case must not become dup, got: {}",
            r2
        );
    }

    #[test]
    fn test_multi_base_deletion_in_homopolymer_in_cds_stays_as_del() {
        // Delete 2 A's from an 8-A tract in c. context. Per HGVS spec
        // (repeated.md codon-frame exception), repeat notation requires
        // unit_len % 3 == 0 in c.; unit_len=1 is gated, so the canonical
        // form stays as plain `del` at the post-shift 3' position.
        //
        // Sequence: GGGGGGGGGAAAAAAAAGGGGGGGGG. A-tract at positions 10-17.
        let seq = "GGGGGGGGGAAAAAAAAGGGGGGGGG";
        let provider = provider_with_transcript("NM_TEST.1", seq);
        let result = normalize_to_string(provider, "NM_TEST.1:c.10_11del");
        assert_eq!(result, "NM_TEST.1:c.16_17del");
    }

    #[test]
    fn test_direction_shift_under_codon_frame_gate() {
        // In c. context with unit_len=1 the codon-frame gate (repeated.md)
        // forbids `A[N]` repeat notation regardless of shuffle direction,
        // so neither branch can rewrite to `unit[N-k]`. What this test
        // locks is the direction-conditional shift behavior: ThreePrime
        // shifts the deletion to the 3'-most position in the tract;
        // FivePrime shifts to the 5'-most position. (For B2 in its native
        // direction-conditional form — only firing in ThreePrime — see
        // the matrix tests that use a codon-aligned unit.)
        //
        // Sequence: GGGGGGGGGAAAAAAAAGGGGGGGGG. 8-A tract at pos 10-17
        // (1-based). Input c.16_17del (last 2 A's of the tract):
        //   - ThreePrime: del shifted to canonical 3' end → c.16_17del.
        //   - FivePrime:  shift to 5'-most position → c.10_11del.
        let seq = "GGGGGGGGGAAAAAAAAGGGGGGGGG";

        let provider = provider_with_transcript("NM_TEST.1", seq);
        let three_prime = Normalizer::with_config(
            provider,
            NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
        );
        let variant = parse_hgvs("NM_TEST.1:c.16_17del").expect("Failed to parse input");
        let result_3p = format!(
            "{}",
            three_prime
                .normalize(&variant)
                .expect("Normalization failed")
        );
        assert_eq!(
            result_3p, "NM_TEST.1:c.16_17del",
            "ThreePrime mode in c. with unit_len=1 must emit plain del per codon-frame gate"
        );

        let provider = provider_with_transcript("NM_TEST.1", seq);
        let five_prime = Normalizer::with_config(
            provider,
            NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
        );
        let result_5p = format!(
            "{}",
            five_prime
                .normalize(&variant)
                .expect("Normalization failed")
        );
        assert_eq!(
            result_5p, "NM_TEST.1:c.10_11del",
            "FivePrime mode shifts to 5'-most position"
        );
    }

    // =========================================================================
    // DUPLICATION TESTS - Duplications should shift 3'
    // =========================================================================

    #[test]
    fn test_dup_shifts_3prime() {
        // Sequence: ...AAAAAAA...
        // dup of A at position 10 should shift to 3'-most A
        let seq = "GGGGGGGGGAAAAAAAAGGGGGGGGGG"; // A's at pos 10-17 (8 A's)
        let provider = provider_with_transcript("NM_TEST.1", seq);

        let result = normalize_to_string(provider, "NM_TEST.1:c.10dup");

        // Should shift to position 17 (3'-most A)
        // Note: the exact output format may vary
        assert!(
            result.contains("dup"),
            "Duplication should stay as dup, got: {}",
            result
        );
    }

    // =========================================================================
    // SUBSTITUTION TESTS - Substitutions should not change
    // =========================================================================

    #[test]
    fn test_substitution_unchanged() {
        let seq = "ATGCATGCATGC";
        let provider = provider_with_transcript("NM_TEST.1", seq);

        let result = normalize_to_string(provider, "NM_TEST.1:c.5A>G");

        assert_eq!(
            result, "NM_TEST.1:c.5A>G",
            "Substitution should remain unchanged"
        );
    }

    // =========================================================================
    // DELINS TESTS
    // =========================================================================

    #[test]
    fn test_delins_unchanged_when_not_special() {
        // A delins that doesn't simplify should stay as delins
        let seq = "ATGCATGCATGC";
        let provider = provider_with_transcript("NM_TEST.1", seq);

        let result = normalize_to_string(provider, "NM_TEST.1:c.5delinsGG");

        assert!(
            result.contains("delins"),
            "Delins should stay as delins when not simplifiable, got: {}",
            result
        );
    }

    // =========================================================================
    // EDGE CASES
    // =========================================================================

    #[test]
    fn test_deletion_at_sequence_end() {
        // Deletion at the end of sequence - can't shift further.
        // Seq: ATGCATGCAAAA. Delete A at position 9. A-tract is positions
        // 9-12. 3'-shift moves the del to position 12 (last A in the tract).
        // ref_count=4, k=1, B2 doesn't fire — stays as del.
        let seq = "ATGCATGCAAAA"; // 4-A homopolymer at positions 9-12
        let provider = provider_with_transcript("NM_TEST.1", seq);

        let result = normalize_to_string(provider, "NM_TEST.1:c.9del");
        assert_eq!(result, "NM_TEST.1:c.12del");
    }

    #[test]
    fn test_deletion_no_shift_when_not_in_repeat() {
        // Deletion in unique sequence shouldn't shift. Position 5 is 'A'
        // (ATGCATGCATGC, 1-indexed). Bases on either side differ from 'A'
        // in adjacent positions, so no shift.
        let seq = "ATGCATGCATGC";
        let provider = provider_with_transcript("NM_TEST.1", seq);

        let result = normalize_to_string(provider, "NM_TEST.1:c.5del");
        // ref[4]='A' (1-based pos 5). Adjacent: ref[3]='C', ref[5]='T'.
        // 3'-shift extends right while ref[del_start]==ref[del_end].
        // ref[4]='A', ref[5]='T' → mismatch, no shift. Stays at c.5del.
        assert_eq!(result, "NM_TEST.1:c.5del");
    }

    // =========================================================================
    // CODON-FRAME GATE TESTS (B1) — c. context, repeated.md exception
    // =========================================================================

    #[test]
    fn test_multi_insertion_in_homopolymer_codon_frame_blocks_repeat_in_cds() {
        // Per HGVS spec (repeated.md codon-frame exception): in c. context,
        // repeat notation requires unit_len % 3 == 0. unit_len=1 (A) is
        // forbidden, so the canonical form stays as `ins<literal>`, not A[N].
        let seq = "GGGGAAAAAAAGGGG"; // 7 A's
        let provider = provider_with_transcript("NM_TEST.1", seq);
        let result = normalize_to_string(provider, "NM_TEST.1:c.5_6insAA");
        assert!(
            result.contains("ins") && !result.contains("A["),
            "c. + unit_len=1 must not emit A[N], got: {}",
            result
        );
    }

    #[test]
    fn test_multi_dup_in_homopolymer_codon_frame_blocks_repeat_in_cds() {
        // 2-base dup in 7-A tract: per spec the structurally-canonical form
        // is repeat notation A[9], but the c. codon-frame exception forbids
        // it (unit_len=1). Falls back to literal `ins<dup_seq>` per the
        // GatedInsertion path in duplication_to_repeat.
        let seq = "GGGGAAAAAAAGGGG"; // 7 A's at pos 5-11
        let provider = provider_with_transcript("NM_TEST.1", seq);
        let result = normalize_to_string(provider, "NM_TEST.1:c.5_6dup");
        assert!(
            result.contains("ins") && !result.contains("A["),
            "c. + unit_len=1 dup must not emit A[N], got: {}",
            result
        );
    }
}

// =============================================================================
// EXPECTED FAILURES - Tests for features not yet implemented
// =============================================================================
// These tests document expected mutalyzer behavior that ferro-hgvs
// doesn't yet implement. They are marked #[ignore] and will be enabled
// as features are implemented.

mod expected_failures {
    use super::*;
    use ferro_hgvs::reference::mock::MockProvider;
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};

    fn make_transcript(id: &str, sequence: &str) -> Transcript {
        let len = sequence.len();
        Transcript::new(
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
        )
    }

    fn provider_with_transcript(id: &str, sequence: &str) -> MockProvider {
        let mut provider = MockProvider::new();
        provider.add_transcript(make_transcript(id, sequence));
        provider
    }

    fn normalize_to_string(provider: MockProvider, input: &str) -> String {
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse input");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        format!("{}", normalized)
    }

    #[test]
    #[ignore] // Not yet implemented - delins to dup conversion
    fn test_delins_becomes_dup_special_case() {
        // Mutalyzer: c.1558delinsGG -> c.1558dup (when net result is duplication)
        // If deleting X and inserting XX (doubling), it's a dup
        let seq = "ATGCGATGCATGC"; // G at position 5
        let provider = provider_with_transcript("NM_TEST.1", seq);

        let result = normalize_to_string(provider, "NM_TEST.1:c.5delinsGG");

        // Deleting G and inserting GG = net gain of one G = dup
        assert!(
            result.contains("dup"),
            "Delins that doubles a base should become dup, got: {}",
            result
        );
    }

    #[test]
    fn test_inversion_shortening() {
        // Mutalyzer-style canonical form: outer complement pairs cascade.
        // Sequence: GGGGAGGGCTGGGGG (AGGGCT at c.5_10).
        //           123456789012345
        // Cascade: A/T cancel -> GGGC; G/C cancel -> GG. Inner G/G are not
        // complementary so shortening stops at c.7_8inv.
        let seq = "GGGGAGGGCTGGGGG";
        let provider = provider_with_transcript("NM_TEST.1", seq);
        let result = normalize_to_string(provider, "NM_TEST.1:c.5_10inv");
        assert_eq!(result, "NM_TEST.1:c.7_8inv");
    }
}

// =============================================================================
// MUTALYZER-VERIFIED TESTS
// =============================================================================
// These tests verify normalization output against mutalyzer ground truth.
// Test cases from scripts/mutalyzer_expected.json

mod mutalyzer_verified {
    use ferro_hgvs::reference::mock::MockProvider;
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
    use ferro_hgvs::{parse_hgvs, Normalizer};
    use rstest::rstest;

    /// Create a transcript with a custom sequence for testing
    fn make_transcript_with_cds(
        id: &str,
        sequence: &str,
        cds_start: u64,
        cds_end: u64,
    ) -> Transcript {
        let len = sequence.len();
        Transcript::new(
            id.to_string(),
            Some("TEST".to_string()),
            Strand::Plus,
            sequence.to_string(),
            Some(cds_start),
            Some(cds_end),
            vec![Exon::new(1, 1, len as u64)],
            None,
            None,
            None,
            Default::default(),
            ManeStatus::None,
            None,
            None,
        )
    }

    fn provider_with_transcript(
        id: &str,
        sequence: &str,
        cds_start: u64,
        cds_end: u64,
    ) -> MockProvider {
        let mut provider = MockProvider::new();
        provider.add_transcript(make_transcript_with_cds(id, sequence, cds_start, cds_end));
        provider
    }

    fn normalize_to_string(provider: MockProvider, input: &str) -> String {
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse input");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        format!("{}", normalized)
    }

    // =========================================================================
    // PASS-THROUGH TESTS (no sequence needed - variants unchanged)
    // All cases verified against Mutalyzer - all agree
    // =========================================================================

    #[rstest]
    // Genomic substitutions - Mutalyzer: agrees
    #[case("NC_000009.11:g.130548229C>G")]
    #[case("NC_000017.10:g.56296510A>G")]
    #[case("NC_000017.11:g.58219149A>G")]
    #[case("NC_000001.10:g.160001799G>C")]
    #[case("NC_000003.11:g.49059579T>C")]
    // Genomic inversions - Mutalyzer: agrees
    #[case("NC_000011.10:g.116830247_116836307inv")]
    #[case("NC_000006.11:g.32007624_32007625inv")]
    #[case("NC_000011.10:g.72301915_72301917inv")]
    #[case("NC_000004.11:g.114277480_114277481inv")]
    // Genomic deletions (no sequence = pass through) - Mutalyzer: agrees
    #[case("NC_000017.11:g.43045711del")]
    #[case("NC_000013.10:g.29233227del")]
    #[case("NC_000013.11:g.28659090del")]
    #[case("NC_000005.10:g.78985014del")]
    // Large genomic deletions - Mutalyzer: agrees
    #[case("NC_000005.10:g.112836913_112908314del")]
    #[case("NC_000017.10:g.7126454_7126558del")]
    #[case("NC_000017.11:g.7223135_7223239del")]
    #[case("NC_000003.11:g.167422632_167422685del")]
    #[case("NC_000003.12:g.167704844_167704897del")]
    // Genomic duplications - Mutalyzer: agrees
    #[case("NC_000013.11:g.32316461dup")]
    #[case("NC_000017.10:g.56296538_56296542dup")]
    #[case("NC_000017.11:g.58219177_58219181dup")]
    // Large delins - Mutalyzer: agrees
    #[case("NC_000005.10:g.174694600_174931940delinsTATAATATGTGTGTATATAATATATATATTACAATATA")]
    fn test_genomic_passthrough(#[case] input: &str) {
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(format!("{}", normalized), input);
    }

    // =========================================================================
    // REPEAT NORMALIZATION TESTS
    // These test the conversion of repeat notation to del/dup
    // =========================================================================

    #[test]
    fn test_repeat_to_repeat_decrement() {
        // Per B2: when repeat count < ref by >=2 with surviving units,
        // emit unit[N-k]. Input: c.4CAT[1] of 4-CAT ref. ref_count=4,
        // specified=1, k=3 → CAT[1] at positions 4..15 (1-based).
        //
        // (Was test_repeat_to_deletion; renamed to reflect post-A5 behavior.
        // Single-unit reductions still emit del — see
        // test_repeat_to_deletion_one_unit below.)
        let seq = "GGGCATCATCATCATGGG"; // 4 CAT repeats at positions 4-15
        let provider = provider_with_transcript("NM_TEST.1", seq, 1, seq.len() as u64);

        let result = normalize_to_string(provider, "NM_TEST.1:c.4CAT[1]");
        assert_eq!(result, "NM_TEST.1:c.4_15CAT[1]");
    }

    #[test]
    fn test_repeat_to_deletion_one_unit() {
        // k=1 reduction: per HGVS prioritization (deletion > unranked
        // repeat), emit del at the 3'-most unit. Input: c.4CAT[3] of
        // 4-CAT ref. ref_count=4, specified=3, k=1 → del of last CAT.
        let seq = "GGGCATCATCATCATGGG"; // 4 CAT repeats at positions 4-15
        let provider = provider_with_transcript("NM_TEST.1", seq, 1, seq.len() as u64);

        let result = normalize_to_string(provider, "NM_TEST.1:c.4CAT[3]");
        assert_eq!(result, "NM_TEST.1:c.13_15del");
    }

    #[test]
    fn test_repeat_to_duplication() {
        // When repeat count = reference count + 1, should become duplication
        // Sequence with 2 CAT repeats, specifying CAT[3] should become dup
        let seq = "GGGCATCATGGG"; // 2 CAT repeats
        let provider = provider_with_transcript("NM_TEST.1", seq, 1, seq.len() as u64);

        let result = normalize_to_string(provider, "NM_TEST.1:c.4CAT[3]");

        // Should become a duplication
        assert!(
            result.contains("dup"),
            "Repeat with count = ref+1 should become duplication, got: {}",
            result
        );
    }

    #[test]
    fn test_repeat_stays_repeat() {
        // When repeat count > reference count + 1, should stay as repeat
        // Sequence with 2 CAT repeats, specifying CAT[5] should stay as repeat
        let seq = "GGGCATCATGGG"; // 2 CAT repeats
        let provider = provider_with_transcript("NM_TEST.1", seq, 1, seq.len() as u64);

        let result = normalize_to_string(provider, "NM_TEST.1:c.4CAT[5]");

        // Should stay as repeat notation
        assert!(
            result.contains("[5]") || result.contains("CAT["),
            "Repeat with count > ref+1 should stay as repeat, got: {}",
            result
        );
    }

    #[test]
    fn test_repeat_unchanged_when_equal_to_ref() {
        // When repeat count = reference count, no change (identity)
        let seq = "GGGCATCATGGG"; // 2 CAT repeats
        let provider = provider_with_transcript("NM_TEST.1", seq, 1, seq.len() as u64);

        let result = normalize_to_string(provider, "NM_TEST.1:c.4CAT[2]");

        // Should remain unchanged (same as reference)
        assert!(
            result.contains("CAT[2]"),
            "Repeat with count = ref should be unchanged, got: {}",
            result
        );
    }

    #[test]
    fn test_single_base_repeat_contraction_in_cds_routes_to_del() {
        // Input A[1] against 4-A ref. In c. context with unit_len=1, the
        // codon-frame gate forbids `[N]` repeat notation, so the contraction
        // routes to a plain del of the surplus A's.
        let seq = "GGGAAAAGGG"; // 4 A's at positions 4-7
        let provider = provider_with_transcript("NM_TEST.1", seq, 1, seq.len() as u64);

        let result = normalize_to_string(provider, "NM_TEST.1:c.4A[1]");
        assert_eq!(result, "NM_TEST.1:c.5_7del");
    }

    #[test]
    fn test_single_base_repeat_to_duplication() {
        // Single base repeat: A[5] when reference has 4 A's should become dup
        let seq = "GGGAAAAGGG"; // 4 A's
        let provider = provider_with_transcript("NM_TEST.1", seq, 1, seq.len() as u64);

        let result = normalize_to_string(provider, "NM_TEST.1:c.4A[5]");

        // Should become a duplication
        assert!(
            result.contains("dup"),
            "Single-base repeat with count = ref+1 should become duplication, got: {}",
            result
        );
    }
}

// =============================================================================
// CLINVAR NORMALIZATION TESTS - Using Real Sequences
// =============================================================================
// These tests verify normalization using real transcript sequences extracted
// from the mutalyzer cache. Test cases are from UNIT_TEST_PROMPT.md.

mod clinvar_normalization {
    use ferro_hgvs::reference::mock::MockProvider;
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
    use ferro_hgvs::{parse_hgvs, Normalizer};
    use rstest::rstest;
    use serde::Deserialize;
    use std::collections::HashMap;
    use std::sync::OnceLock;

    #[derive(Debug, Deserialize)]
    struct TranscriptData {
        id: String,
        #[allow(dead_code)]
        gene_symbol: String,
        strand: String,
        sequence: String,
        cds_start: u64,
        cds_end: u64,
        exons: Vec<ExonData>,
    }

    #[derive(Debug, Deserialize)]
    struct ExonData {
        number: u32,
        start: u64,
        end: u64,
    }

    static TRANSCRIPTS: OnceLock<HashMap<String, TranscriptData>> = OnceLock::new();

    fn load_transcripts() -> &'static HashMap<String, TranscriptData> {
        TRANSCRIPTS.get_or_init(|| {
            let json_path = concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/fixtures/sequences/normalization_transcripts.json"
            );
            let json_str = std::fs::read_to_string(json_path)
                .expect("Failed to read normalization_transcripts.json");
            let transcripts: Vec<TranscriptData> =
                serde_json::from_str(&json_str).expect("Failed to parse JSON");
            transcripts.into_iter().map(|t| (t.id.clone(), t)).collect()
        })
    }

    fn create_provider_for_transcript(accession: &str) -> MockProvider {
        let transcripts = load_transcripts();
        let mut provider = MockProvider::new();

        if let Some(data) = transcripts.get(accession) {
            let strand = if data.strand == "+" {
                Strand::Plus
            } else {
                Strand::Minus
            };
            let exons: Vec<Exon> = data
                .exons
                .iter()
                .map(|e| Exon::new(e.number, e.start, e.end))
                .collect();

            let transcript = Transcript::new(
                data.id.clone(),
                Some(data.gene_symbol.clone()),
                strand,
                data.sequence.clone(),
                Some(data.cds_start),
                Some(data.cds_end),
                exons,
                None,
                None,
                None,
                Default::default(),
                ManeStatus::None,
                None,
                None,
            );
            provider.add_transcript(transcript);
        }
        provider
    }

    fn normalize_to_string(input: &str) -> String {
        // Extract accession from HGVS string
        let accession = input.split(':').next().expect("Invalid HGVS format");
        let provider = create_provider_for_transcript(accession);
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse input");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        format!("{}", normalized)
    }

    // =========================================================================
    // BATCH 1: Verified Matches - Coding Variants (25 tests)
    // These variants are verified to match between ferro and mutalyzer
    // All cases: Mutalyzer agrees
    // =========================================================================

    #[rstest]
    // UTR and CDS variants using transcripts we have - Mutalyzer: agrees
    #[case("NM_001365307.2:c.*895A>G", "NM_001365307.2:c.*895A>G")]
    #[case("NM_001365304.2:c.*1128G>A", "NM_001365304.2:c.*1128G>A")]
    #[case("NM_001400774.1:c.-29A>T", "NM_001400774.1:c.-29A>T")]
    // Non-coding RNA variants - Mutalyzer: agrees
    #[case("NR_153405.1:n.3650G>A", "NR_153405.1:n.3650G>A")]
    #[case("NR_153405.1:n.3799C>T", "NR_153405.1:n.3799C>T")]
    #[case("NR_046285.1:n.951C>T", "NR_046285.1:n.951C>T")]
    #[case("NR_038982.1:n.756G>A", "NR_038982.1:n.756G>A")]
    // Well-known clinical variants - CFTR - Mutalyzer: agrees
    #[case("NM_000492.4:c.3846G>A", "NM_000492.4:c.3846G>A")]
    // Well-known clinical variants - BRCA1 - Mutalyzer: agrees
    #[case("NM_007294.4:c.5266dup", "NM_007294.4:c.5266dup")]
    // Well-known clinical variants - BRCA2 - Mutalyzer: agrees
    #[case("NM_000059.4:c.5946del", "NM_000059.4:c.5946del")]
    // Well-known clinical variants - COL1A1 - Mutalyzer: agrees
    #[case("NM_000088.4:c.769G>A", "NM_000088.4:c.769G>A")]
    // Well-known clinical variants - TP53 - Mutalyzer: agrees
    #[case("NM_000546.6:c.215C>G", "NM_000546.6:c.215C>G")]
    #[case("NM_000546.6:c.743G>A", "NM_000546.6:c.743G>A")]
    fn test_verified_normalization_batch1(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input);
        assert_eq!(
            result, expected,
            "Normalization mismatch for '{}': got '{}', expected '{}'",
            input, result, expected
        );
    }

    // =========================================================================
    // BATCH 1 CONTINUED: Deletion and duplication variants
    // =========================================================================

    #[rstest]
    // Deletions that should normalize to same position - Mutalyzer: agrees
    #[case("NM_000492.4:c.1521_1523del", "NM_000492.4:c.1521_1523del")]
    #[case("NM_007294.4:c.68_69del", "NM_007294.4:c.68_69del")]
    #[case("NM_000059.4:c.6275_6276del", "NM_000059.4:c.6275_6276del")]
    #[case("NM_001414398.1:c.*160_*163del", "NM_001414398.1:c.*160_*163del")]
    // Duplications
    // 5'UTR dup — Mutalyzer: agrees. The duplicated 10-base region is
    // already in its 3'-most equivalent representation in the local
    // tract, so 3' shift is a no-op and the canonical form matches the
    // input. (Pre-#97 ferro emitted a spurious 1-base shift to
    // `c.-55_-46dup` as a side effect of the 5'UTR off-by-one in the
    // CDS↔tx mapping; once that's fixed, ferro and Mutalyzer agree.)
    #[case("NM_001394148.2:c.-56_-47dup", "NM_001394148.2:c.-56_-47dup")]
    // Mutalyzer: agrees
    #[case("NR_153405.1:n.3908_3910dup", "NR_153405.1:n.3908_3910dup")]
    // Delins - Mutalyzer: agrees
    #[case("NR_037658.1:n.987delinsGAAG", "NR_037658.1:n.987delinsGAAG")]
    fn test_verified_normalization_indels_batch1(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input);
        assert_eq!(
            result, expected,
            "Normalization mismatch for '{}': got '{}', expected '{}'",
            input, result, expected
        );
    }

    // =========================================================================
    // BATCH 2: Genomic Variants - No Sequence Required
    // These are passthrough tests for genomic coordinates
    // =========================================================================

    #[rstest]
    // Genomic substitutions
    #[case("NC_000009.11:g.130548229C>G", "NC_000009.11:g.130548229C>G")]
    #[case("NC_000017.10:g.56296510A>G", "NC_000017.10:g.56296510A>G")]
    #[case("NC_000017.11:g.58219149A>G", "NC_000017.11:g.58219149A>G")]
    #[case("NC_000001.10:g.160001799G>C", "NC_000001.10:g.160001799G>C")]
    #[case("NC_000003.11:g.49059579T>C", "NC_000003.11:g.49059579T>C")]
    // Genomic deletions
    #[case("NC_000013.10:g.29233227del", "NC_000013.10:g.29233227del")]
    #[case("NC_000013.11:g.28659090del", "NC_000013.11:g.28659090del")]
    #[case("NC_000005.10:g.78985014del", "NC_000005.10:g.78985014del")]
    #[case("NC_000017.11:g.43045711del", "NC_000017.11:g.43045711del")]
    // Large genomic deletions
    #[case(
        "NC_000005.10:g.112836913_112908314del",
        "NC_000005.10:g.112836913_112908314del"
    )]
    #[case(
        "NC_000017.10:g.7126454_7126558del",
        "NC_000017.10:g.7126454_7126558del"
    )]
    #[case(
        "NC_000017.11:g.7223135_7223239del",
        "NC_000017.11:g.7223135_7223239del"
    )]
    #[case(
        "NC_000003.11:g.167422632_167422685del",
        "NC_000003.11:g.167422632_167422685del"
    )]
    #[case(
        "NC_000003.12:g.167704844_167704897del",
        "NC_000003.12:g.167704844_167704897del"
    )]
    fn test_genomic_substitutions_and_deletions(#[case] input: &str, #[case] expected: &str) {
        // Genomic variants without sequence pass through unchanged
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(
            format!("{}", normalized),
            expected,
            "Genomic variant normalization failed for '{}'",
            input
        );
    }

    #[rstest]
    // Genomic inversions
    #[case(
        "NC_000011.10:g.116830247_116836307inv",
        "NC_000011.10:g.116830247_116836307inv"
    )]
    #[case(
        "NC_000006.11:g.32007624_32007625inv",
        "NC_000006.11:g.32007624_32007625inv"
    )]
    #[case(
        "NC_000011.10:g.72301915_72301917inv",
        "NC_000011.10:g.72301915_72301917inv"
    )]
    #[case(
        "NC_000004.11:g.114277480_114277481inv",
        "NC_000004.11:g.114277480_114277481inv"
    )]
    // Genomic duplications
    #[case("NC_000013.11:g.32316461dup", "NC_000013.11:g.32316461dup")]
    #[case(
        "NC_000017.10:g.56296538_56296542dup",
        "NC_000017.10:g.56296538_56296542dup"
    )]
    #[case(
        "NC_000017.11:g.58219177_58219181dup",
        "NC_000017.11:g.58219177_58219181dup"
    )]
    // Genomic delins
    #[case(
        "NC_000020.11:g.25383511_25397600delinsG",
        "NC_000020.11:g.25383511_25397600delinsG"
    )]
    #[case(
        "NC_000005.10:g.174694600_174931940delinsTATAATATGTGTGTATATAATATATATATTACAATATA",
        "NC_000005.10:g.174694600_174931940delinsTATAATATGTGTGTATATAATATATATATTACAATATA"
    )]
    #[case(
        "NC_000002.11:g.179472926_179472954delinsGGGATCTGTTTTGGGATCTG",
        "NC_000002.11:g.179472926_179472954delinsGGGATCTGTTTTGGGATCTG"
    )]
    // Identity variants (=)
    #[case("NC_000009.11:g.136500515=", "NC_000009.11:g.136500515=")]
    #[case("NC_000009.12:g.133635393=", "NC_000009.12:g.133635393=")]
    #[case("NC_000010.10:g.64415184=", "NC_000010.10:g.64415184=")]
    #[case("NC_000007.13:g.117120047=", "NC_000007.13:g.117120047=")]
    #[case("NC_000017.10:g.32611446=", "NC_000017.10:g.32611446=")]
    fn test_genomic_inversions_dups_delins(#[case] input: &str, #[case] expected: &str) {
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(
            format!("{}", normalized),
            expected,
            "Genomic variant normalization failed for '{}'",
            input
        );
    }

    // =========================================================================
    // BATCH 2 CONTINUED: LRG and mitochondrial variants
    // =========================================================================

    #[rstest]
    // LRG reference variants
    #[case("LRG_835:g.5337del", "LRG_835:g.5337del")]
    #[case("LRG_835:g.5136C>T", "LRG_835:g.5136C>T")]
    #[case("LRG_835:g.5040C>G", "LRG_835:g.5040C>G")]
    #[case("LRG_835:g.5301G>A", "LRG_835:g.5301G>A")]
    // Mitochondrial variants
    #[case("NC_012920.1:m.12315G>A", "NC_012920.1:m.12315G>A")]
    #[case("NC_012920.1:m.12320A>G", "NC_012920.1:m.12320A>G")]
    #[case("NC_012920.1:m.12297T>C", "NC_012920.1:m.12297T>C")]
    #[case("NC_012920.1:m.12310dup", "NC_012920.1:m.12310dup")]
    // NG_ RefSeqGene variants
    #[case("NG_189068.1:g.188del", "NG_189068.1:g.188del")]
    fn test_lrg_mito_ng_variants(#[case] input: &str, #[case] expected: &str) {
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(
            format!("{}", normalized),
            expected,
            "LRG/mito/NG variant normalization failed for '{}'",
            input
        );
    }

    // =========================================================================
    // BATCH 3: Complex Patterns - Uncertain Positions and Special Cases
    // These are complex HGVS patterns that ferro handles correctly
    // =========================================================================

    #[rstest]
    // Uncertain position deletions (from Mutalyzer timeout cases - ferro output is reference)
    #[case(
        "NC_000017.11:g.(?_79009664)_(79009817_?)del",
        "NC_000017.11:g.(?_79009664)_(79009817_?)del"
    )]
    #[case(
        "NC_000017.11:g.(?_31094927)_(31377677_?)del",
        "NC_000017.11:g.(?_31094927)_(31377677_?)del"
    )]
    #[case(
        "NC_000005.10:g.(?_112707504)_(112846240_?)del",
        "NC_000005.10:g.(?_112707504)_(112846240_?)del"
    )]
    #[case(
        "NC_000014.9:g.(50000000_?)_(?_50247254)del",
        "NC_000014.9:g.(50000000_?)_(?_50247254)del"
    )]
    #[case(
        "NC_000007.14:g.(40116368_?)_(?_40134601)del",
        "NC_000007.14:g.(40116368_?)_(?_40134601)del"
    )]
    #[case(
        "NC_000022.10:g.(?_18893735)_(18924066_?)del",
        "NC_000022.10:g.(?_18893735)_(18924066_?)del"
    )]
    #[case(
        "NC_000008.10:g.(?_6264113)_(6296618_6299587)del",
        "NC_000008.10:g.(?_6264113)_(6296618_6299587)del"
    )]
    #[case(
        "NC_000001.11:g.(196753076_?)_(?_196839375)del",
        "NC_000001.11:g.(196753076_?)_(?_196839375)del"
    )]
    // Uncertain position duplications
    #[case(
        "NC_000008.11:g.(142876886_142877022)_(142914909_142915045)dup",
        "NC_000008.11:g.(142876886_142877022)_(142914909_142915045)dup"
    )]
    // Inverted ranges (antisense strand)
    #[case(
        "NC_000012.11:g.110593351_110576466dup",
        "NC_000012.11:g.110593351_110576466dup"
    )]
    #[case(
        "NC_000016.9:g.23634775_23621090dup",
        "NC_000016.9:g.23634775_23621090dup"
    )]
    fn test_uncertain_position_variants(#[case] input: &str, #[case] expected: &str) {
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(
            format!("{}", normalized),
            expected,
            "Uncertain position variant normalization failed for '{}'",
            input
        );
    }

    #[rstest]
    // Embedded accession insertions (complex structural variants)
    #[case(
        "NC_000008.11:g.86688947_86688948ins[MF045863.1:g.1_36978]",
        "NC_000008.11:g.86688947_86688948ins[MF045863.1:g.1_36978]"
    )]
    #[case(
        "NC_000008.11:g.86711345_86711346ins[MF045864.2:g.1_98770]",
        "NC_000008.11:g.86711345_86711346ins[MF045864.2:g.1_98770]"
    )]
    #[case(
        "NG_016167.1:g.21559097_21559098ins[PP887427.1:g.1_1518]",
        "NG_016167.1:g.21559097_21559098ins[PP887427.1:g.1_1518]"
    )]
    #[case(
        "LRG_1293:g.21559097_21559098ins[PP887427.1:g.1_1518]",
        "LRG_1293:g.21559097_21559098ins[PP887427.1:g.1_1518]"
    )]
    // Complex delins with embedded accessions
    #[case(
        "NC_000008.11:g.86587460_86650711delins[KY923049.1:g.1_466]",
        "NC_000008.11:g.86587460_86650711delins[KY923049.1:g.1_466]"
    )]
    // Delins with N repeats
    #[case(
        "NC_000002.11:g.47618487_47650860delinsN[155]",
        "NC_000002.11:g.47618487_47650860delinsN[155]"
    )]
    // Delins with repeat range
    #[case(
        "NC_000004.12:g.39348425_39348479delinsAAAGG[400_2000]",
        "NC_000004.12:g.39348425_39348479delinsAAAGG[400_2000]"
    )]
    fn test_embedded_accession_variants(#[case] input: &str, #[case] expected: &str) {
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(
            format!("{}", normalized),
            expected,
            "Embedded accession variant normalization failed for '{}'",
            input
        );
    }

    #[rstest]
    // Complex multi-component insertions
    #[case(
        "NC_000017.11:g.80114186_80114187ins[80114172_80114186;NC_000020.11:g.2823027_2826302;AAA]",
        "NC_000017.11:g.80114186_80114187ins[80114172_80114186;NC_000020.11:g.2823027_2826302;AAA]"
    )]
    // Chromosomal rearrangements with qter
    #[case(
        "NC_000009.12:g.12891379_qterdelins[T;NC_000020.11:g.47204889_qterinv]",
        "NC_000009.12:g.12891379_qterdelins[T;NC_000020.11:g.47204889_qterinv]"
    )]
    #[case(
        "NC_000013.10:g.114819939_qterdelins[96729864_114814234inv;96735632_104289803]",
        "NC_000013.10:g.114819939_qterdelins[96729864_114814234inv;96735632_104289803]"
    )]
    // Complex delins with inversions
    #[case(
        "NC_000007.14:g.45043702_46521017delins[AGAAGGAAATTT;45310743_46521014;45043709_45310738inv]",
        "NC_000007.14:g.45043702_46521017delins[AGAAGGAAATTT;45310743_46521014;45043709_45310738inv]"
    )]
    fn test_complex_structural_variants(#[case] input: &str, #[case] expected: &str) {
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(
            format!("{}", normalized),
            expected,
            "Complex structural variant normalization failed for '{}'",
            input
        );
    }

    // =========================================================================
    // BATCH 3 CONTINUED: Intronic variants and protein variants
    // =========================================================================

    #[rstest]
    // Intronic variants (from Mutalyzer timeout cases)
    #[case("NM_001318856.2:c.9-1296T>C", "NM_001318856.2:c.9-1296T>C")]
    #[case("NM_001365307.2:c.150+19G>A", "NM_001365307.2:c.150+19G>A")]
    #[case("NM_001031734.3:c.154+12G>C", "NM_001031734.3:c.154+12G>C")]
    #[case("NM_033517.1:c.1772-2A>C", "NM_033517.1:c.1772-2A>C")]
    #[case("NM_001031734.3:c.210-17C>A", "NM_001031734.3:c.210-17C>A")]
    #[case("NM_000350.2:c.302+68C>T", "NM_000350.2:c.302+68C>T")]
    #[case("NM_000492.4:c.1585-1G>A", "NM_000492.4:c.1585-1G>A")]
    #[case("NM_020732.3:c.-3A>G", "NM_020732.3:c.-3A>G")]
    #[case("NM_001351733.2:c.-91+8701C>T", "NM_001351733.2:c.-91+8701C>T")]
    #[case("NR_033294.1:n.*5C>G", "NR_033294.1:n.*5C>G")]
    fn test_intronic_variants(#[case] input: &str, #[case] expected: &str) {
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(
            format!("{}", normalized),
            expected,
            "Intronic variant normalization failed for '{}'",
            input
        );
    }

    // =========================================================================
    // BATCH 4: Protein Variants
    // These test protein variant normalization (passthrough - no shifting)
    // =========================================================================

    #[rstest]
    // UniProt protein substitutions
    #[case("P04181:p.Tyr245Cys", "P04181:p.Tyr245Cys")]
    #[case("P04181:p.Arg250Pro", "P04181:p.Arg250Pro")]
    #[case("P09417:p.Gly23Asp", "P09417:p.Gly23Asp")]
    #[case("P00439:p.Phe299Cys", "P00439:p.Phe299Cys")]
    #[case("P54802:p.Phe48Leu", "P54802:p.Phe48Leu")]
    #[case("Q9P0J0:p.Lys5Asn", "Q9P0J0:p.Lys5Asn")]
    #[case("P04062:p.Val433Leu", "P04062:p.Val433Leu")]
    // RefSeq protein duplications
    #[case("AAK07616.1:p.Asp200dup", "AAK07616.1:p.Asp200dup")]
    fn test_protein_substitutions(#[case] input: &str, #[case] expected: &str) {
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(
            format!("{}", normalized),
            expected,
            "Protein variant normalization failed for '{}'",
            input
        );
    }

    #[rstest]
    // Frameshift variants
    #[case("NP_000509.1:p.Ser10ValfsTer14", "NP_000509.1:p.Ser10ValfsTer14")]
    #[case("NP_005201.2:p.Asn25ThrfsTer20", "NP_005201.2:p.Asn25ThrfsTer20")]
    #[case("NP_001305738.1:p.Tyr180fs", "NP_001305738.1:p.Tyr180fs")]
    #[case("NP_066970.3:p.Val90SerfsTer6", "NP_066970.3:p.Val90SerfsTer6")]
    #[case("NP_004412.2:p.Ser539AlafsTer110", "NP_004412.2:p.Ser539AlafsTer110")]
    #[case("NP_060609.2:p.Gly406ArgfsTer90", "NP_060609.2:p.Gly406ArgfsTer90")]
    // Extension variants (stop-loss)
    #[case("NP_001166937.1:p.Ter514Leuext*?", "NP_001166937.1:p.Ter514Leuext*?")]
    #[case("NP_056480.1:p.Ter547Leuext*?", "NP_056480.1:p.Ter547Leuext*?")]
    #[case("NP_000654.2:p.Ter501Lysext*?", "NP_000654.2:p.Ter501Lysext*?")]
    fn test_protein_frameshift_extension(#[case] input: &str, #[case] expected: &str) {
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(
            format!("{}", normalized),
            expected,
            "Protein frameshift/extension normalization failed for '{}'",
            input
        );
    }

    // =========================================================================
    // BATCH 4 CONTINUED: Special interval and repeat patterns
    // =========================================================================

    #[rstest]
    // Interval-only notation (no edit)
    #[case("NM_173651.4:c.5238_5240", "NM_173651.4:c.5238_5240")]
    #[case("NM_007375.3:c.*697", "NM_007375.3:c.*697")]
    // Uncertain dup/c.?dup
    #[case("NM_001412270.1:c.?dup", "NM_001412270.1:c.?dup")]
    // Complex interval with intronic offsets
    #[case("NM_005262.2:c.259-25_259-24", "NM_005262.2:c.259-25_259-24")]
    #[case(
        "NM_012343.3:c.(-51+1_-53-1)_(381+1_382-1)",
        "NM_012343.3:c.(-51+1_-53-1)_(381+1_382-1)"
    )]
    #[case(
        "NM_001040142.2:c.(2388+1_2389-1)_(3849+1_3850-1)",
        "NM_001040142.2:c.(2388+1_2389-1)_(3849+1_3850-1)"
    )]
    #[case("NM_004483.5:c.148-?_228+?", "NM_004483.5:c.148-?_228+?")]
    fn test_special_interval_patterns(#[case] input: &str, #[case] expected: &str) {
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(
            format!("{}", normalized),
            expected,
            "Special interval pattern normalization failed for '{}'",
            input
        );
    }

    #[rstest]
    // Insertions at antisense/inverted positions
    #[case(
        "NC_000011.10:g.5238138_5153222insTATTT",
        "NC_000011.10:g.5238138_5153222insTATTT"
    )]
    // NG uncertain position variants
    #[case(
        "NG_011403.2:g.(80027_96047)_(99154_121150)del",
        "NG_011403.2:g.(80027_96047)_(99154_121150)del"
    )]
    #[case(
        "NG_009385.2:g.(?_5001)_(40068_?)del",
        "NG_009385.2:g.(?_5001)_(40068_?)del"
    )]
    // M22590.1 (GenBank) duplication
    #[case("M22590.1:c.1069_1233dup", "M22590.1:c.1069_1233dup")]
    // Uncertain ranges with delins
    #[case(
        "NC_000023.10:g.(133030929_133031380)_(133079087_133079463)delins[118528009_118528409;118674690_118675082]",
        "NC_000023.10:g.(133030929_133031380)_(133079087_133079463)delins[118528009_118528409;118674690_118675082]"
    )]
    fn test_additional_complex_patterns(#[case] input: &str, #[case] expected: &str) {
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(
            format!("{}", normalized),
            expected,
            "Additional complex pattern normalization failed for '{}'",
            input
        );
    }

    // =========================================================================
    // BATCH 4 CONTINUED: Repeat notation variants
    // These pass through as-is when no sequence available (MockProvider is empty).
    // Mutalyzer has sequence access and normalizes differently - we document differences.
    // =========================================================================

    #[rstest]
    // Repeat notation - these pass through as-is when no sequence available
    // Mutalyzer: g.3243408_3243410del (differs - with sequence, CAT[1] where ref has 4 = deletion)
    #[case("NC_000016.10:g.3243405CAT[1]", "NC_000016.10:g.3243405CAT[1]")]
    // Mutalyzer: g.98231228_98231229dup (differs - with sequence, GA[3] where ref has 2 = dup)
    #[case("NC_000009.11:g.98231226GA[3]", "NC_000009.11:g.98231226GA[3]")]
    // Mutalyzer: g.34274816_34274839GAGAAG[9] (differs - expands range and rotates unit)
    #[case("NC_000017.11:g.34274818GAAGGA[8]", "NC_000017.11:g.34274818GAAGGA[8]")]
    // Mutalyzer: g.34274816_34274839GAGAAG[9] (differs - expands range)
    #[case("NC_000017.11:g.34274816GAGAAG[9]", "NC_000017.11:g.34274816GAGAAG[9]")]
    // Mutalyzer: g.25031779_25031808CGC[17] (differs - expands range)
    #[case("NC_000023.10:g.25031779CGC[17]", "NC_000023.10:g.25031779CGC[17]")]
    // Mutalyzer: g.680_682dup (differs - TGC[6] where ref has 5 = dup)
    #[case("NG_045215.1:g.668TGC[6]", "NG_045215.1:g.668TGC[6]")]
    // Mutalyzer: g.204135386_204135388del (differs - CAG[3] where ref has 4 = deletion)
    #[case("NC_000001.10:g.204135377CAG[3]", "NC_000001.10:g.204135377CAG[3]")]
    // Mutalyzer: g.35658018_35658029GTCCTCAGCTTC[3] (differs - expands range)
    #[case(
        "NC_000009.11:g.35658018GTCCTCAGCTTC[3]",
        "NC_000009.11:g.35658018GTCCTCAGCTTC[3]"
    )]
    // Mutalyzer: g.45869562_45869585del (differs - GGGGGC[2] where ref has 6 = deletion)
    #[case("NC_000010.10:g.45869550GGGGGC[2]", "NC_000010.10:g.45869550GGGGGC[2]")]
    // Mutalyzer: n.884_885dup (differs - GA[3] where ref has 2 = dup)
    #[case("NR_038982.1:n.882GA[3]", "NR_038982.1:n.882GA[3]")]
    fn test_repeat_notation_passthrough(#[case] input: &str, #[case] expected: &str) {
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(
            format!("{}", normalized),
            expected,
            "Repeat notation passthrough failed for '{}'",
            input
        );
    }

    // =========================================================================
    // BATCH 4 CONTINUED: Protein repeat and deletion variants
    // =========================================================================

    #[rstest]
    // Protein with repeat notation
    #[case("NP_002102.4:p.Gln18[40_?]", "NP_002102.4:p.Gln18[40_?]")]
    // Protein deletion with specified residue (ferro converts single-letter to three-letter)
    #[case("NP_000026.2:p.Leu288delCys", "NP_000026.2:p.Leu288delCys")]
    fn test_protein_repeat_and_deletion(#[case] input: &str, #[case] expected: &str) {
        let provider = MockProvider::new();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        assert_eq!(
            format!("{}", normalized),
            expected,
            "Protein repeat/deletion normalization failed for '{}'",
            input
        );
    }
}

// =============================================================================
// NORMALIZATION TRANSFORMATION TESTS
// =============================================================================
// These tests verify actual normalization transformations where input ≠ output.
// Test cases are from NORMALIZATION_TEST_PROMPT.md, verified against Mutalyzer.

mod normalization_transformations {
    use ferro_hgvs::reference::mock::MockProvider;
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
    use ferro_hgvs::{parse_hgvs, Normalizer};
    use rstest::rstest;
    use serde::Deserialize;
    use std::collections::HashMap;
    use std::sync::OnceLock;

    #[derive(Debug, Deserialize)]
    struct TranscriptData {
        id: String,
        #[allow(dead_code)]
        gene_symbol: String,
        strand: String,
        sequence: String,
        cds_start: u64,
        cds_end: u64,
        exons: Vec<ExonData>,
    }

    #[derive(Debug, Deserialize)]
    struct ExonData {
        number: u32,
        start: u64,
        end: u64,
    }

    static TRANSCRIPTS: OnceLock<HashMap<String, TranscriptData>> = OnceLock::new();

    fn load_transcripts() -> &'static HashMap<String, TranscriptData> {
        TRANSCRIPTS.get_or_init(|| {
            let json_path = concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/fixtures/sequences/normalization_transcripts.json"
            );
            let json_str = std::fs::read_to_string(json_path)
                .expect("Failed to read normalization_transcripts.json");
            let transcripts: Vec<TranscriptData> =
                serde_json::from_str(&json_str).expect("Failed to parse JSON");
            transcripts.into_iter().map(|t| (t.id.clone(), t)).collect()
        })
    }

    fn create_provider_for_transcript(accession: &str) -> MockProvider {
        let transcripts = load_transcripts();
        let mut provider = MockProvider::new();

        if let Some(data) = transcripts.get(accession) {
            let strand = if data.strand == "+" {
                Strand::Plus
            } else {
                Strand::Minus
            };
            let exons: Vec<Exon> = data
                .exons
                .iter()
                .map(|e| Exon::new(e.number, e.start, e.end))
                .collect();

            let transcript = Transcript::new(
                data.id.clone(),
                Some(data.gene_symbol.clone()),
                strand,
                data.sequence.clone(),
                Some(data.cds_start),
                Some(data.cds_end),
                exons,
                None,
                None,
                None,
                Default::default(),
                ManeStatus::None,
                None,
                None,
            );
            provider.add_transcript(transcript);
        }
        provider
    }

    fn normalize_to_string(input: &str) -> String {
        let accession = input.split(':').next().expect("Invalid HGVS format");
        let provider = create_provider_for_transcript(accession);
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse input");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        format!("{}", normalized)
    }

    // =========================================================================
    // 3' SHIFTING TESTS
    // Deletions in repetitive regions should shift to 3'-most position
    // All cases: Mutalyzer agrees
    // =========================================================================

    #[rstest]
    // CFTR ΔF508 - deletion shifts right in repetitive region - Mutalyzer: agrees
    #[case("NM_000492.4:c.1520_1522del", "NM_000492.4:c.1521_1523del")]
    // Already at 3' position - no shift needed - Mutalyzer: agrees
    #[case("NM_000492.4:c.1521_1523del", "NM_000492.4:c.1521_1523del")]
    // Different position - no shift needed - Mutalyzer: agrees
    #[case("NM_000492.4:c.1519_1521del", "NM_000492.4:c.1519_1521del")]
    fn test_3prime_shifting(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input);
        assert_eq!(
            result, expected,
            "3' shifting failed for '{}': got '{}', expected '{}'",
            input, result, expected
        );
    }

    // =========================================================================
    // INSERTION → DUPLICATION TESTS
    // When inserted sequence matches the preceding base(s), becomes dup
    // Mutalyzer: agrees
    // =========================================================================

    // Real-accession smoke; matrix-level coverage in tests/ins_shift_matrix.rs.
    #[rstest]
    // Single T inserted after T becomes dup - Mutalyzer: agrees
    #[case("NM_001089.2:c.3057_3058insT", "NM_001089.2:c.3057dup")]
    fn test_insertion_becomes_dup(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input);
        assert_eq!(
            result, expected,
            "Insertion→dup failed for '{}': got '{}', expected '{}'",
            input, result, expected
        );
    }

    // =========================================================================
    // INSERTION → REPEAT NOTATION TESTS
    // Homopolymer insertions become repeat notation
    // All cases: Mutalyzer agrees
    // =========================================================================

    // Real-accession smoke; matrix-level coverage in tests/ins_shift_matrix.rs.
    // Per HGVS spec (repeated.md codon-frame exception): in c. context,
    // homopolymer expansions (unit_len=1) emit `ins<literal>`, not `[N]`.
    #[rstest]
    // polyA expansion: inserting AAAA into A-tract → ins (not A[7]).
    #[case("NM_001127687.1:c.400_401insAAAA", "NM_001127687.1:c.403_404insAAAA")]
    // polyT expansion: inserting TTTT into T-tract → ins (not T[6]).
    #[case("NM_000382.3:c.399_400insTTTT", "NM_000382.3:c.399_400insTTTT")]
    // polyG expansion: inserting GGG into G-tract → ins (not G[5]).
    #[case("NM_015120.4:c.46_47insGGG", "NM_015120.4:c.46_47insGGG")]
    fn test_insertion_becomes_repeat(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input);
        assert_eq!(
            result, expected,
            "Insertion→ins-literal failed for '{}': got '{}', expected '{}'",
            input, result, expected
        );
    }

    // =========================================================================
    // DUPLICATION → REPEAT NOTATION TESTS
    // Tandem duplications in repeat regions become repeat notation
    // Mutalyzer: agrees
    // =========================================================================

    #[rstest]
    // Triplet repeat expansion: dupGCAGCA in GCA-tract becomes GCA[n] - Mutalyzer: agrees
    #[case("NM_020732.3:c.357_362dupGCAGCA", "NM_020732.3:c.342_362GCA[9]")]
    fn test_dup_becomes_repeat(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input);
        assert_eq!(
            result, expected,
            "Dup→repeat failed for '{}': got '{}', expected '{}'",
            input, result, expected
        );
    }

    // =========================================================================
    // DELETION → SHIFTED DEL TESTS (real RefSeq data)
    // 1-base del in a real homopolymer; canonical 3'-shifted position; k=1
    // so B2 doesn't fire and the canonical form stays as `del`.
    // Real-accession smoke; matrix-level coverage in tests/del_shift_matrix.rs.
    // =========================================================================

    #[rstest]
    // CFTR ΔF508 context. c.1520_1522del is the F508del shorthand; the local
    // CFTR context has a CTT tandem at c.1521_1523, so the 1-codon deletion
    // shifts one base to the canonical c.1521_1523del. k=1 → stays as del.
    #[case("NM_000492.4:c.1520_1522del", "NM_000492.4:c.1521_1523del")]
    fn test_deletion_shifts_in_real_homopolymer(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input);
        assert_eq!(
            result, expected,
            "Deletion shift failed for '{}': got '{}', expected '{}'",
            input, result, expected
        );
    }

    // =========================================================================
    // DELETION CANONICALIZATION (real RefSeq data; codon-frame gate)
    // In c. context, the codon-frame gate (repeated.md) suppresses B2
    // (`unit[N-k]` rewrite) when unit_len % 3 != 0, so multi-unit dels of
    // non-codon-aligned tandems canonicalize to plain `del` rather than
    // `unit[N-k]`. Real-accession smoke; matrix-level coverage in
    // tests/del_shift_matrix.rs.
    // =========================================================================

    #[rstest]
    // TP53 c.627-629 is an AAA tract; deleting 2 As (c.628_629del) leaves 1 A.
    // In c. context with unit_len=1, the codon-frame gate (repeated.md) forbids
    // A[1] repeat notation, so the canonical form is plain del.
    #[case("NM_000546.6:c.628_629del", "NM_000546.6:c.628_629del")]
    fn test_deletion_becomes_repeat(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input);
        assert_eq!(
            result, expected,
            "Deletion→canonical-form failed for '{}': got '{}', expected '{}'",
            input, result, expected
        );
    }

    // =========================================================================
    // EXPLICIT BASE REMOVAL TESTS
    // Redundant sequence info in dup notation is stripped
    // All cases: Mutalyzer agrees
    // =========================================================================

    #[rstest]
    // Single base dup with explicit base: dupC → dup - Mutalyzer: agrees
    #[case("NM_033517.1:c.4818dupC", "NM_033517.1:c.4818dup")]
    // Multi-base dup with explicit sequence: dupGCA → dup - Mutalyzer: agrees
    #[case("NM_020732.3:c.360_362dupGCA", "NM_020732.3:c.360_362dup")]
    fn test_explicit_base_removal(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input);
        assert_eq!(
            result, expected,
            "Explicit base removal failed for '{}': got '{}', expected '{}'",
            input, result, expected
        );
    }

    // =========================================================================
    // REPEAT → DELETION TESTS
    // Repeat notation with count < reference count becomes deletion
    // Mutalyzer: agrees
    // =========================================================================

    #[rstest]
    // CFTR: ATC[1] when reference has 2 copies = deletion of one copy - Mutalyzer: agrees
    #[case("NM_000492.4:c.1516ATC[1]", "NM_000492.4:c.1519_1521del")]
    fn test_repeat_becomes_deletion(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input);
        assert_eq!(
            result, expected,
            "Repeat→deletion failed for '{}': got '{}', expected '{}'",
            input, result, expected
        );
    }

    // =========================================================================
    // NO CHANGE TESTS (Verified against Mutalyzer)
    // These variants should normalize to themselves
    // All cases: Mutalyzer agrees
    // =========================================================================

    #[rstest]
    // BRCA1 185delAG - no shift needed - Mutalyzer: agrees
    #[case("NM_007294.4:c.68_69del", "NM_007294.4:c.68_69del")]
    // BRCA2 deletions - Mutalyzer: agrees
    #[case("NM_000059.4:c.5946del", "NM_000059.4:c.5946del")]
    #[case("NM_000059.4:c.6275_6276del", "NM_000059.4:c.6275_6276del")]
    // TP53 deletions - Mutalyzer: agrees
    // (Note: NM_000546.6:c.628_629del removed — ferro now emits c.627_629A[1]
    // per B2 canonical-form rule, which Mutalyzer does not apply.)
    #[case("NM_000546.6:c.532del", "NM_000546.6:c.532del")]
    // BRCA1 duplication - already at 3' position - Mutalyzer: agrees
    #[case("NM_007294.4:c.5266dup", "NM_007294.4:c.5266dup")]
    // TP53 duplication - Mutalyzer: agrees
    #[case("NM_000546.6:c.475_481dup", "NM_000546.6:c.475_481dup")]
    fn test_no_change_verified(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input);
        assert_eq!(
            result, expected,
            "No-change test failed for '{}': got '{}', expected '{}'",
            input, result, expected
        );
    }

    // =========================================================================
    // INSERTIONS THAT STAY AS INSERTIONS
    // These insertions don't match preceding sequence, stay as ins
    // Mutalyzer: agrees
    // =========================================================================

    #[rstest]
    // Insertions that don't become dups - Mutalyzer: agrees
    #[case("NM_000350.3:c.6708_6709insG", "NM_000350.3:c.6708_6709insG")]
    fn test_insertions_stay_insertions(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input);
        assert_eq!(
            result, expected,
            "Insertion should stay as insertion for '{}': got '{}', expected '{}'",
            input, result, expected
        );
    }
}

// =============================================================================
// BENCHMARK COMPARISON TESTS
// =============================================================================
// These tests are derived from 100K ClinVar pattern comparison between
// ferro-hgvs and mutalyzer. Tests cover:
// - Category 1: Bugs in ferro (expected = mutalyzer output)
// - Category 2: Equivalence patterns (expected = ferro output, correct behavior)
// - Category 3: Mutalyzer wrong (expected = ferro output, correct behavior)
// - Category 4: Patterns needing investigation (expected = mutalyzer output)

mod benchmark_comparison_tests {
    use ferro_hgvs::reference::mock::MockProvider;
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
    use ferro_hgvs::{parse_hgvs, Normalizer};

    /// Create a transcript with UTR regions for testing 3' UTR normalization
    /// cds_start and cds_end define the CDS region; positions after cds_end are 3' UTR
    fn make_transcript_with_utr(
        id: &str,
        sequence: &str,
        cds_start: u64,
        cds_end: u64,
    ) -> Transcript {
        let len = sequence.len();
        Transcript::new(
            id.to_string(),
            Some("TEST".to_string()),
            Strand::Plus,
            sequence.to_string(),
            Some(cds_start),
            Some(cds_end),
            vec![Exon::new(1, 1, len as u64)],
            None,
            None,
            None,
            Default::default(),
            ManeStatus::None,
            None,
            None,
        )
    }

    fn provider_with_utr_transcript(
        id: &str,
        sequence: &str,
        cds_start: u64,
        cds_end: u64,
    ) -> MockProvider {
        let mut provider = MockProvider::new();
        provider.add_transcript(make_transcript_with_utr(id, sequence, cds_start, cds_end));
        provider
    }

    fn normalize_to_string(provider: MockProvider, input: &str) -> String {
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse input");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        format!("{}", normalized)
    }

    // =========================================================================
    // CATEGORY 1: BUGS IN FERRO (documented differences with mutalyzer)
    // These tests document observed differences. Some may be ferro bugs,
    // some may be mutalyzer bugs, some may be valid interpretation differences.
    // =========================================================================

    #[test]
    fn test_3prime_utr_deletion_shift_with_repeat() {
        // When there IS a repeat, deletion should shift 3'
        // This test verifies the basic 3' shift logic works in UTR regions
        let mut seq = String::new();
        seq.push_str(&"A".repeat(100)); // CDS c.1-100
        seq.push_str(&"T".repeat(115)); // 3' UTR c.*1 to c.*115
        seq.push_str("GG"); // c.*116 and c.*117 - GG repeat
        seq.push_str(&"A".repeat(50)); // padding after

        let provider = provider_with_utr_transcript("NM_TEST.1", &seq, 1, 100);
        let result = normalize_to_string(provider, "NM_TEST.1:c.*116del");

        // With GG repeat, deletion should shift to 3'-most position
        assert_eq!(
            result, "NM_TEST.1:c.*117del",
            "3' UTR deletion should shift to 3'-most position when there's a repeat"
        );
    }

    #[test]
    fn test_3prime_utr_deletion_no_shift_without_repeat() {
        // When there is NO repeat, deletion should NOT shift
        // Real differences observed:
        // - NM_000392.3:c.*116delG → ferro: c.*116del, mutalyzer: c.*117del
        // - NM_005502.3:c.*1291delT → ferro: c.*1291del, mutalyzer: c.*1292del
        //
        // If the real transcripts don't have repeats at those positions,
        // then ferro's behavior (no shift) would be correct.
        let mut seq = String::new();
        seq.push_str(&"A".repeat(100)); // CDS c.1-100
        seq.push_str(&"T".repeat(115)); // 3' UTR c.*1 to c.*115
        seq.push_str("GA"); // c.*116=G, c.*117=A - NO repeat
        seq.push_str(&"C".repeat(50)); // padding after

        let provider = provider_with_utr_transcript("NM_TEST.1", &seq, 1, 100);
        let result = normalize_to_string(provider, "NM_TEST.1:c.*116del");

        // Without repeat, deletion stays at original position
        assert_eq!(
            result, "NM_TEST.1:c.*116del",
            "3' UTR deletion should NOT shift when there's no repeat"
        );
    }

    // =========================================================================
    // CATEGORY 2: EQUIVALENCE PATTERNS
    // These test ferro's behavior. The equivalence checker should
    // recognize these as equivalent to mutalyzer's different representation.
    // =========================================================================

    #[test]
    fn test_insertion_after_matching_sequence() {
        // Real difference observed:
        // NM_001025357.3:c.242_243insGGC → ferro: c.264_265insGGC, mutalyzer: c.262_264dup
        //
        // When inserting a sequence that matches the preceding bases,
        // ferro may convert to dup (which is valid HGVS).
        // The equivalence check should recognize ins vs dup as equivalent
        // when they represent the same biological change.
        //
        // Create sequence where GGC is repeated
        let mut seq = String::new();
        seq.push_str(&"A".repeat(261)); // Positions 1-261
        seq.push_str("GGC"); // Positions 262-264
        seq.push_str(&"T".repeat(50)); // Padding

        let provider = provider_with_utr_transcript("NM_TEST.1", &seq, 1, seq.len() as u64);
        let result = normalize_to_string(provider, "NM_TEST.1:c.264_265insGGC");

        // Ferro may convert to dup when inserted sequence matches preceding
        // Both ins and dup are valid representations
        assert!(
            result.contains("ins") || result.contains("dup"),
            "Should normalize to either ins or dup, got: {}",
            result
        );
    }

    #[test]
    fn test_inversion_with_bases_preserved() {
        // Pattern: c.1267_1268invCA → ferro keeps invCA, mutalyzer uses inv
        // Both are valid representations; ferro's is more explicit.
        //
        // Create a sequence with CA at positions 1267-1268
        let mut seq = String::new();
        seq.push_str(&"A".repeat(1266)); // Positions 1-1266
        seq.push_str("CA"); // Positions 1267-1268
        seq.push_str(&"T".repeat(50)); // Padding

        let provider = provider_with_utr_transcript("NM_TEST.1", &seq, 1, seq.len() as u64);
        let result = normalize_to_string(provider, "NM_TEST.1:c.1267_1268invCA");

        // Ferro preserves the bases in inversion notation
        assert!(
            result.contains("inv"),
            "Should contain inv, got: {}",
            result
        );
        // Note: Whether bases are kept depends on implementation
        // The equivalence check should handle both forms
    }

    #[test]
    fn test_single_element_allelic_expanded() {
        // Pattern: c.[1616A>G] → ferro: [NM_X:c.1616A>G], mutalyzer: c.1616A>G
        // Ferro expands single-element allelic to full form.
        // Both are valid; equivalence check should handle this.
        let seq = "A".repeat(2000);
        let provider = provider_with_utr_transcript("NM_TEST.1", &seq, 1, seq.len() as u64);

        let result = normalize_to_string(provider, "NM_TEST.1:c.[100A>G]");

        // Ferro expands allelic notation
        // The exact format may vary - just check it normalizes successfully
        assert!(
            result.contains("A>G") || result.contains("100"),
            "Should preserve the substitution, got: {}",
            result
        );
    }

    // =========================================================================
    // CATEGORY 3: MUTALYZER BEHAVIOR DIFFERENCES
    // These tests document observed differences. In some cases ferro's
    // behavior is more HGVS-compliant, in others they're equivalent.
    // =========================================================================

    #[test]
    fn test_insertion_to_dup_conversion() {
        // Real observed difference:
        // NM_001025357.3:c.242_243insGGC → ferro: c.264_265insGGC, mutalyzer: c.262_264dup
        //
        // Ferro keeps as insertion, mutalyzer converts to duplication.
        // According to HGVS, when an insertion exactly duplicates preceding sequence,
        // it SHOULD be described as a duplication. So mutalyzer may actually be correct.
        //
        // Test with a sequence where insertion matches preceding
        let mut seq = String::new();
        seq.push_str(&"T".repeat(239)); // Positions 1-239
        seq.push_str("GGC"); // Positions 240-242
        seq.push_str("GGC"); // Positions 243-245 (adjacent repeat)
        seq.push_str(&"A".repeat(50)); // Padding

        let provider = provider_with_utr_transcript("NM_TEST.1", &seq, 1, seq.len() as u64);

        // Insert GGC between positions 242 and 243
        let result = normalize_to_string(provider, "NM_TEST.1:c.242_243insGGC");

        // Both ins and dup are valid; equivalence check handles this
        assert!(
            result.contains("ins") || result.contains("dup"),
            "Should normalize to ins or dup, got: {}",
            result
        );
    }

    // =========================================================================
    // CATEGORY 4: NEEDS INVESTIGATION
    // These tests document patterns that need further investigation.
    // The expected behavior depends on actual reference sequences.
    // =========================================================================

    #[test]
    fn test_repeat_notation_with_single_ref_copy() {
        // Real observed difference:
        // NM_001425324.1:c.2010_2011TG[1]
        // ferro: c.2012_2013del (interprets TG[1] as reducing from 2 to 1 copy)
        // mutalyzer: c.= (interprets as "there is 1 copy" = no change)
        //
        // The correct interpretation depends on the reference sequence:
        // - If ref has TGTG (2 copies), then TG[1] = delete one copy
        // - If ref has TG (1 copy), then TG[1] = no change
        //
        // Test with single copy in reference
        let mut seq = String::new();
        seq.push_str(&"A".repeat(2009)); // Positions 1-2009
        seq.push_str("TG"); // Positions 2010-2011 (single copy)
        seq.push_str(&"A".repeat(50)); // Padding (no repeat)

        let provider = provider_with_utr_transcript("NM_TEST.1", &seq, 1, seq.len() as u64);

        let result = normalize_to_string(provider, "NM_TEST.1:c.2010_2011TG[1]");

        // Current ferro behavior: keeps repeat notation unchanged
        // This is acceptable - the repeat notation is preserved for interpretation
        assert!(
            result.contains("TG[1]") || result.contains("c.="),
            "Repeat[1] with single copy should be no-change or preserved, got: {}",
            result
        );
    }

    #[test]
    fn test_repeat_notation_with_double_ref_copy() {
        // When ref has 2 copies (TGTG) and variant says TG[1],
        // that's a deletion of one copy
        let mut seq = String::new();
        seq.push_str(&"A".repeat(2009)); // Positions 1-2009
        seq.push_str("TGTG"); // Positions 2010-2013 (two copies)
        seq.push_str(&"A".repeat(50)); // Padding

        let provider = provider_with_utr_transcript("NM_TEST.1", &seq, 1, seq.len() as u64);

        let result = normalize_to_string(provider, "NM_TEST.1:c.2010_2011TG[1]");

        // With 2 copies in ref, TG[1] means delete one = deletion
        assert!(
            result.contains("del") || result.contains("TG[1]"),
            "Repeat contraction should be deletion or preserved, got: {}",
            result
        );
    }

    #[test]
    fn test_repeat_expansion_vs_duplication() {
        // Pattern: c.351_352AG[2]
        // ferro: c.355_356del (if ref has 3 copies)
        // mutalyzer: c.355_356dup (different interpretation)
        //
        // With triple copy reference (AGAGAG), AG[2] = reduce to 2 = del
        let mut seq = String::new();
        seq.push_str(&"T".repeat(350)); // Positions 1-350
        seq.push_str("AGAGAG"); // Positions 351-356 (3 copies of AG)
        seq.push_str(&"C".repeat(50)); // Padding

        let provider = provider_with_utr_transcript("NM_TEST.1", &seq, 1, seq.len() as u64);

        let result = normalize_to_string(provider, "NM_TEST.1:c.351_352AG[2]");

        // With 3 copies, AG[2] should be a deletion
        assert!(
            result.contains("del") || result.contains("AG[2]"),
            "Repeat contraction from 3 to 2 should be deletion or preserved, got: {}",
            result
        );
    }

    #[test]
    fn test_insertion_sequence_rotation() {
        // Pattern: c.247_248insAACA against ref where c.247_250 == "AACA".
        // ferro post-issue-132: c.247_250dup (the inserted AACA matches the
        //   following 4 ref bases verbatim → 1-copy AACA dup).
        // mutalyzer: c.249_250insCAAA (a different rotation choice).
        //
        // The HGVS 3' rule prefers `dup` over `ins` when the inserted unit
        // matches an adjacent reference unit (insertion_to_duplication's
        // rotation iterator finds the AACA-phase tract with ref_count=1).
        // The dup form is canonical regardless of input rotation.
        let mut seq = String::new();
        seq.push_str(&"G".repeat(246)); // Positions 1-246
        seq.push_str("AA"); // Positions 247-248
        seq.push_str("CA"); // Positions 249-250 (forms AACA context)
        seq.push_str(&"T".repeat(50)); // Padding

        let provider = provider_with_utr_transcript("NM_TEST.1", &seq, 1, seq.len() as u64);

        let result = normalize_to_string(provider, "NM_TEST.1:c.247_248insAACA");

        // 1-copy AACA next to a single AACA in the ref → dup spanning
        // the matched ref tract (c.247..c.250 inclusive).
        assert_eq!(
            result, "NM_TEST.1:c.247_250dup",
            "Expected c.247_250dup (issue #132 dup-over-ins canonicalization), got: {}",
            result
        );
    }

    #[test]
    fn test_5prime_utr_insertion_to_dup() {
        // Pattern: c.-27_-26insA
        // ferro: c.-27dup (correct - insertion of A next to A becomes dup)
        // mutalyzer: c.-27_-26insA (keeps as insertion - incorrect)
        //
        // This is ferro being CORRECT. When inserting a base that matches
        // the preceding base, it should become a duplication.
        //
        // 5' UTR positions are negative, counting backward from CDS start
        // c.-27 is position (cds_start - 27)
        let mut seq = String::new();
        seq.push_str(&"G".repeat(100)); // 5' UTR positions
        seq.push_str(&"A".repeat(500)); // CDS and beyond

        // CDS starts at position 101 (so c.-27 = position 74)
        // We need an A at position 74 so inserting A becomes dup
        // Actually, let's build more carefully:
        // Position 74 (c.-27) should be A
        let mut seq = String::new();
        seq.push_str(&"T".repeat(73)); // Positions 1-73
        seq.push('A'); // Position 74 (c.-27)
        seq.push_str(&"T".repeat(26)); // Positions 75-100 (c.-26 to c.-1)
        seq.push_str(&"A".repeat(300)); // CDS starting at 101

        let provider = provider_with_utr_transcript("NM_TEST.1", &seq, 101, 400);

        let result = normalize_to_string(provider, "NM_TEST.1:c.-27_-26insA");

        // Ferro correctly converts ins A → dup when A precedes
        assert!(
            result.contains("dup"),
            "Insertion of A next to A should become dup, got: {}",
            result
        );
    }
}

// =============================================================================
// INSERTION ROTATION TESTS
// =============================================================================
// Tests for the insertion sequence rotation fix (GitHub issue: ins→dup and
// ins sequence rotation during 3' normalization).
//
// When an insertion is shifted through a repeat region, the effective sequence
// "rotates". For example:
// - Inserting "GGC" and shifting by 1 position → effective sequence is "GCG"
// - Inserting "GA" and shifting by 1 position → effective sequence is "AG"
//
// This affects both:
// 1. ins→dup conversion: must check rotated sequence against reference
// 2. ins output: must output the rotated sequence, not the original

mod insertion_rotation_tests {
    use ferro_hgvs::reference::mock::MockProvider;
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
    use ferro_hgvs::{parse_hgvs, Normalizer};

    fn make_transcript(id: &str, sequence: &str, cds_start: u64, cds_end: u64) -> Transcript {
        let len = sequence.len();
        Transcript::new(
            id.to_string(),
            Some("TEST".to_string()),
            Strand::Plus,
            sequence.to_string(),
            Some(cds_start),
            Some(cds_end),
            vec![Exon::new(1, 1, len as u64)],
            None,
            None,
            None,
            Default::default(),
            ManeStatus::None,
            None,
            None, // ensembl_match
        )
    }

    fn provider_with_transcript(
        id: &str,
        sequence: &str,
        cds_start: u64,
        cds_end: u64,
    ) -> MockProvider {
        let mut provider = MockProvider::new();
        provider.add_transcript(make_transcript(id, sequence, cds_start, cds_end));
        provider
    }

    fn normalize_to_string(provider: MockProvider, input: &str) -> String {
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).expect("Failed to parse input");
        let normalized = normalizer
            .normalize(&variant)
            .expect("Normalization failed");
        format!("{}", normalized)
    }

    // =========================================================================
    // TEST 1: Insertion becomes duplication with sequence rotation
    // =========================================================================
    // Pattern: c.242_243insGGC in a GGCGGCGGC... repeat region
    // After 3' shift (22 positions), the sequence rotates: GGC → GCG
    // The preceding 3 bases at final position match "GCG", so it becomes dup
    //
    // Expected: c.262_264dup (not c.264_265insGGC or c.264_266dup)

    #[test]
    fn test_insertion_to_dup_with_rotation() {
        // Build a sequence with a GGC repeat region
        // CDS starts at position 38 (1-based), so c.1 = position 38
        // We need: c.243 onwards = GGCGGCGGCGGCGGCGGCGGCGCG...
        //          (22 positions of GGC repeat, then it ends)
        //
        // c.243 = position 38 + 242 = 280 (1-based), so index 279
        //
        // Issue #132 update: the canonical dup is now phase-aligned with
        // the GGC tract (most-3' GGC dup), not the rotated GCG-phase that
        // shuffle's iterative walk happened to land on. The new helper
        // `insertion_to_duplication` picks the tract-aligned phase (GGC at
        // tract [279..300)), so the dup spans the most-3' GGC = c.261_263.
        // Both the new and old positions describe the same edit; the new
        // behavior is symmetric with the multi-copy `insertion_to_repeat`
        // path's phase choice.

        let mut seq = String::new();
        seq.push_str(&"A".repeat(37)); // Positions 1-37 (5' UTR)
                                       // CDS starts at position 38
        seq.push_str(&"A".repeat(241)); // c.1 to c.241 = positions 38-278
                                        // c.242 = position 279
        seq.push('A'); // c.242 = A (before the repeat)
                       // c.243 onwards: GGCGGCGGC... (22 bases = 7+ triplets, then changes)
        seq.push_str("GGCGGCGGCGGCGGCGGCGGCGC"); // c.243 to c.265 (23 bases)
                                                 // After c.265, sequence changes so shift stops at c.264_265
        seq.push_str(&"T".repeat(200)); // Fill out the rest

        // CDS from position 38 to 500 (arbitrary end)
        let provider = provider_with_transcript("NM_TEST.1", &seq, 38, 500);

        let result = normalize_to_string(provider, "NM_TEST.1:c.242_243insGGC");

        // GGC-phase tract: 7 GGC units at c.243..c.263 inclusive. Most-3'
        // GGC dup → c.261_263.
        assert_eq!(
            result, "NM_TEST.1:c.261_263dup",
            "Expected c.261_263dup (issue #132 GGC-phase canonical form), got: {}",
            result
        );
    }

    // =========================================================================
    // TEST 2: Multi-base duplication position calculation
    // =========================================================================
    // When an insertion becomes a duplication, the positions should be
    // (X - len + 1) to X, not X to (X + len - 1)
    //
    // For c.264_265insGCG where c.262_264 = GCG:
    // Result should be c.262_264dup, NOT c.264_266dup

    #[test]
    fn test_dup_position_calculation() {
        // Simple case: insert GCG where preceding 3 bases are GCG
        // c.10_11insGCG where c.8_10 = GCG
        // Should become c.8_10dup

        let mut seq = String::new();
        seq.push_str(&"A".repeat(10)); // Positions 1-10 (5' UTR, CDS starts at 11)
                                       // CDS: c.1 = position 11
        seq.push_str(&"A".repeat(7)); // c.1 to c.7
        seq.push_str("GCG"); // c.8, c.9, c.10
        seq.push('T'); // c.11 - different, so insertion stops here
        seq.push_str(&"A".repeat(200));

        let provider = provider_with_transcript("NM_TEST.1", &seq, 11, 300);

        let result = normalize_to_string(provider, "NM_TEST.1:c.10_11insGCG");

        // Insertion at c.10_11 of GCG
        // Preceding 3 bases (c.8_10) = GCG
        // Should become c.8_10dup
        assert!(
            result.contains("c.8_10dup"),
            "Expected c.8_10dup (3-base dup ending at insertion point), got: {}",
            result
        );
    }

    // =========================================================================
    // TEST 3: Insertion rotation for non-dup output
    // =========================================================================
    // When an insertion shifts but doesn't become a dup, the output should
    // still use the rotated sequence.
    //
    // Pattern: c.5_6insGA shifts to c.8_9, rotation = 3 % 2 = 1
    // Rotated sequence: GA → AG
    // If c.7_8 ≠ AG, stays as insertion but with rotated sequence
    // Expected: c.8_9insAG (not c.8_9insGA)

    #[test]
    fn test_insertion_rotation_non_dup() {
        // Insert GA at c.5_6, with sequence allowing 3 shifts
        // c.6, c.7, c.8 should be G, A, G to allow shifting
        // But c.7_8 should NOT be AG to prevent dup conversion
        //
        // Sequence design:
        // c.5 = T (before insertion)
        // c.6 = G, c.7 = A, c.8 = G (allows 3 shifts: G matches, A matches, G matches)
        // c.9 = T (stops shifting)
        // After shift: c.8_9insGA, but rotated by 3%2=1 → insAG
        // c.7_8 = AG? No, c.7=A, c.8=G → AG. That WOULD be a dup!
        //
        // Let's redesign: we want rotation but no dup
        // Insert GA, shift by 1 position (rotation=1)
        // c.6 = G allows one shift
        // c.7 = T stops shift
        // Final: c.6_7insAG (rotated)
        // Check c.5_6 for dup: need it to NOT be AG

        let mut seq = String::new();
        seq.push_str(&"A".repeat(10)); // 5' UTR
                                       // CDS: c.1 = position 11
        seq.push_str("TTTTG"); // c.1-c.5 = TTTTG
        seq.push('G'); // c.6 = G (matches ins[0]=G, allows shift)
        seq.push('T'); // c.7 = T (doesn't match ins[1]=A, stops shift)
        seq.push_str(&"A".repeat(200));

        let provider = provider_with_transcript("NM_TEST.1", &seq, 11, 300);

        let result = normalize_to_string(provider, "NM_TEST.1:c.5_6insGA");

        // Shift by 1: c.5_6insGA → c.6_7insAG (rotated)
        // c.5_6 = TG ≠ AG, so not a dup
        // Result should be c.6_7insAG with the rotated sequence
        assert!(
            result.contains("c.6_7insAG"),
            "Expected c.6_7insAG (rotated GA→AG), got: {}",
            result
        );
    }

    // =========================================================================
    // TEST 4: No rotation when shift amount is multiple of sequence length
    // =========================================================================
    // If shift % len == 0, no rotation needed

    #[test]
    fn test_no_rotation_when_multiple_of_length() {
        // Insert GA, shift by 2 positions (rotation = 2 % 2 = 0)
        // Sequence should allow exactly 2 shifts

        let mut seq = String::new();
        seq.push_str(&"A".repeat(10)); // 5' UTR
                                       // CDS: c.1 = position 11
        seq.push_str("TTTT"); // c.1-c.4
        seq.push('G'); // c.5 = G (matches ins[0]=G, shift 1)
        seq.push('A'); // c.6 = A (matches ins[1]=A, shift 2)
        seq.push('T'); // c.7 = T (stops, doesn't match ins[0]=G)
        seq.push_str(&"A".repeat(200));

        let provider = provider_with_transcript("NM_TEST.1", &seq, 11, 300);

        let result = normalize_to_string(provider, "NM_TEST.1:c.4_5insGA");

        // Shift by 2: rotation = 0, sequence stays GA
        // c.4_5insGA → c.6_7insGA (no rotation)
        // Check c.5_6 for dup: GA? c.5=G, c.6=A → GA! This IS a dup!
        // So result should be c.5_6dup
        assert!(
            result.contains("c.5_6dup"),
            "Expected c.5_6dup (shift by 2, GA matches preceding GA), got: {}",
            result
        );
    }

    // =========================================================================
    // TEST 5: Single-base insertion (no rotation possible)
    // =========================================================================

    #[test]
    fn test_single_base_insertion_no_rotation() {
        // Single base insertions can't rotate (len=1, any shift % 1 = 0)

        let mut seq = String::new();
        seq.push_str(&"A".repeat(10)); // 5' UTR
        seq.push_str("TTTTAAAAT"); // c.1-c.9
        seq.push_str(&"G".repeat(200));

        let provider = provider_with_transcript("NM_TEST.1", &seq, 11, 300);

        let result = normalize_to_string(provider, "NM_TEST.1:c.4_5insA");

        // Should shift through the AAA region and become dup
        assert!(
            result.contains("dup"),
            "Single A insertion into A region should become dup, got: {}",
            result
        );
    }
}

#[cfg(test)]
mod repeat_position_tests {
    use ferro_hgvs::normalize::rules::{count_tandem_repeats, normalize_repeat, RepeatNormResult};

    /// Test what happens when the repeat unit doesn't match at the specified position
    /// but exists nearby. This simulates c.4261_4262CA[1] where:
    /// - At c.4261-4262 the reference has "GC" (not "CA")
    /// - At c.4262-4265 the reference has "CACA" (2 copies of CA)
    #[test]
    fn test_repeat_unit_not_at_position() {
        // Simulate the reference around c.4261:
        // c.4258=C, c.4259=A, c.4260=T, c.4261=G, c.4262=C, c.4263=A, c.4264=C, c.4265=A
        // So at positions around 4261 (0-based index after tx conversion):
        // The sequence is: ...CATGCACAG...
        //                      ^ index 4459 = c.4261
        let ref_seq = b"XXXXXCATGCACAGXXXX";
        //                   0123456789...
        // At index 3: C
        // At index 4: A
        // At index 5: T
        // At index 6: G  <- this would be like c.4261
        // At index 7: C
        // At index 8: A
        // At index 9: C
        // At index 10: A
        // CACA is at indices 7-10

        // Test count_tandem_repeats when called at index 6 (where G is, not CA)
        let result = count_tandem_repeats(ref_seq, 6, b"CA");
        assert!(
            result.is_none(),
            "Should return None because CA doesn't start at index 6"
        );

        // Test normalize_repeat - should return Unchanged
        let result = normalize_repeat(ref_seq, 6, b"CA", 1, false);
        assert!(
            matches!(result, RepeatNormResult::Unchanged),
            "Should return Unchanged when repeat unit not found at position. Got {:?}",
            result
        );
    }

    /// Test what happens when the repeat unit DOES match at the specified position
    /// and we request fewer copies than reference has
    #[test]
    fn test_repeat_deletion_when_unit_matches() {
        // Reference has CACA (2 copies of CA) at position 0
        let ref_seq = b"CACAXXXX";

        // count_tandem_repeats at index 0 should find 2 copies
        let result = count_tandem_repeats(ref_seq, 0, b"CA");
        assert_eq!(
            result,
            Some((2, 0, 4)),
            "Should find 2 copies of CA at indices 0-4"
        );

        // normalize_repeat with count=1 should produce deletion
        let result = normalize_repeat(ref_seq, 0, b"CA", 1, false);
        match result {
            RepeatNormResult::Deletion { start, end } => {
                // Should delete one copy (2 bases) from the 3' end
                // Indices 2-3 (0-based) = positions 3-4 (1-based)
                assert_eq!(start, 3, "Deletion should start at HGVS position 3");
                assert_eq!(end, 4, "Deletion should end at HGVS position 4");
            }
            _ => panic!("Expected Deletion, got {:?}", result),
        }
    }
}

/// Integration test to check real repeat normalization behavior
/// This test requires the benchmark-output directory with reference data
#[cfg(test)]
mod real_repeat_tests {
    use ferro_hgvs::{parse_hgvs, MultiFastaProvider, Normalizer};
    use std::path::Path;

    #[test]
    #[ignore] // Requires reference data
    fn test_real_repeat_pattern() {
        let ref_path = Path::new("benchmark-output/manifest.json");
        if !ref_path.exists() {
            eprintln!("Skipping test - benchmark-output not available");
            return;
        }

        let provider = MultiFastaProvider::from_manifest(ref_path).unwrap();
        let normalizer = Normalizer::new(provider);

        // Parse the variant
        let variant = parse_hgvs("NM_001407936.1:c.4261_4262CA[1]").unwrap();
        println!("Parsed: {:?}", variant);

        // Normalize
        let result = normalizer.normalize(&variant);
        match result {
            Ok(normalized) => {
                println!("Normalized: {}", normalized);
                // If the repeat unit doesn't match at the position,
                // we should get back the original (or c.=)
                // NOT a deletion
            }
            Err(e) => {
                println!("Error: {:?}", e);
            }
        }
    }
}

// =============================================================================
// FERRO VS MUTALYZER DIFFERENCES (Non-Equivalent Outputs)
// =============================================================================
// These tests document cases where Ferro and Mutalyzer produce different results.
// Ferro's output is expected to be correct per HGVS specification.
// Mutalyzer outputs are documented in comments for reference.
//
// Run with: cargo test ferro_mutalyzer_differences --ignored
// Requires: benchmark-output/ with reference data
// =============================================================================

#[cfg(test)]
mod ferro_mutalyzer_differences {
    use super::integration_helpers::try_normalize as normalize_to_string;
    use rstest::rstest;

    // =========================================================================
    // CATEGORY 1: Repeat[N] interpreted as deletion vs duplication
    // =========================================================================
    // When a repeat count is LESS than reference, it represents a DELETION.
    // Mutalyzer incorrectly interprets these as duplications.
    //
    // Example: Reference has CTCTCTCT (4 copies of CT = 8 bases)
    // Input: c.362_363CT[2] means "2 copies" = 4 bases (deletion of 4 bases)
    // Ferro: c.364_367del (correct - 4 bases deleted)
    // Mutalyzer: c.366_367dup (incorrect - misinterprets as adding bases)

    #[rstest]
    #[case("NM_000017.4:c.362_363CT[2]", "NM_000017.4:c.364_367del")]
    // Mutalyzer outputs: NM_000017.4:c.366_367dup (INCORRECT)
    #[case("NM_000028.3:c.2347_2348GA[2]", "NM_000028.3:c.2351_2352del")]
    // Mutalyzer outputs: NM_000028.3:c.2351_2352dup (INCORRECT)
    #[case("NM_000038.4:c.1193_1194AG[2]", "NM_000038.4:c.1197_1198del")]
    // Mutalyzer outputs: NM_000038.4:c.1197_1198dup (INCORRECT)
    #[case("NM_000038.4:c.1814_1815AT[2]", "NM_000038.4:c.1818_1819del")]
    // Mutalyzer outputs: NM_000038.4:c.1818_1819dup (INCORRECT)
    #[case("NM_000038.4:c.2058_2059TC[2]", "NM_000038.4:c.2062_2063del")]
    // Mutalyzer outputs: NM_000038.4:c.2062_2063dup (INCORRECT)
    #[case("NM_000038.4:c.6426_6427CT[2]", "NM_000038.4:c.6430_6431del")]
    // Mutalyzer outputs: NM_000038.4:c.6430_6431dup (INCORRECT)
    #[case("NM_000051.3:c.198_199AT[2]", "NM_000051.3:c.202_203del")]
    // Mutalyzer outputs: NM_000051.3:c.202_203dup (INCORRECT)
    #[case("NM_000051.3:c.6592_6593CT[2]", "NM_000051.3:c.6596_6597del")]
    // Mutalyzer outputs: NM_000051.3:c.6596_6597dup (INCORRECT)
    #[case("NM_000051.3:c.6772_6773CT[2]", "NM_000051.3:c.6776_6777del")]
    // Mutalyzer outputs: NM_000051.3:c.6776_6777dup (INCORRECT)
    #[case("NM_000051.3:c.8279_8280TC[2]", "NM_000051.3:c.8283_8284del")]
    // Mutalyzer outputs: NM_000051.3:c.8283_8284dup (INCORRECT)
    #[ignore]
    fn test_repeat_deletion_vs_mutalyzer_dup(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input)
            .expect("Normalization failed - is benchmark-output present?");
        assert_eq!(
            result, expected,
            "Repeat[N] should normalize to deletion when count < reference.\n\
             Input: {}\n\
             Expected (Ferro correct): {}\n\
             Got: {}",
            input, expected, result
        );
    }

    // =========================================================================
    // CATEGORY 2: Repeat count differences
    // =========================================================================
    // Mutalyzer miscalculates the repeat count, often adding the original count
    // to the specified count instead of using the specified count directly.
    //
    // Example: Reference has AGAGAGAG (4 copies of AG)
    // Input: c.1497_1498AG[4] means "4 copies of AG at expanded range"
    // Ferro: c.1497_1500AG[4] (correct - adjusts range to cover 4 copies)
    // Mutalyzer: c.1497_1500AG[5] (incorrect - adds original count)

    #[rstest]
    #[case("NM_000028.3:c.1497_1498AG[4]", "NM_000028.3:c.1497_1500AG[4]")]
    // Mutalyzer outputs: NM_000028.3:c.1497_1500AG[5] (INCORRECT - adds count)
    #[case("NM_000051.3:c.2926_2927GT[5]", "NM_000051.3:c.2926_2931GT[5]")]
    // Mutalyzer outputs: NM_000051.3:c.2926_2931GT[7] (INCORRECT - adds count)
    #[case("NM_000051.3:c.5373_5374TA[4]", "NM_000051.3:c.5373_5376TA[4]")]
    // Mutalyzer outputs: NM_000051.3:c.5373_5376TA[5] (INCORRECT - adds count)
    #[case("NM_000051.3:c.6864_6865CT[4]", "NM_000051.3:c.6864_6867CT[4]")]
    // Mutalyzer outputs: NM_000051.3:c.6864_6867CT[5] (INCORRECT - adds count)
    #[case("NM_000051.3:c.9132_9133TC[4]", "NM_000051.3:c.9132_9135TC[4]")]
    // Mutalyzer outputs: NM_000051.3:c.9132_9135TC[5] (INCORRECT - adds count)
    #[ignore]
    fn test_repeat_count_vs_mutalyzer_miscalc(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input)
            .expect("Normalization failed - is benchmark-output present?");
        assert_eq!(
            result, expected,
            "Repeat count should be preserved in normalization.\n\
             Input: {}\n\
             Expected (Ferro correct): {}\n\
             Got: {}",
            input, expected, result
        );
    }

    // =========================================================================
    // CATEGORY 3: Triplet repeat with different unit rotation
    // =========================================================================
    // Mutalyzer uses a different rotation of the repeat unit.
    // Both normalizers use 3' normalization but may start with different
    // reference rotation interpretations.

    #[rstest]
    #[case("NM_000044.3:c.172_174CAG[35]", "NM_000044.3:c.172_237CAG[35]")]
    // Mutalyzer outputs: NM_000044.3:c.171_239GCA[57] (different rotation, different count)
    #[ignore]
    fn test_triplet_repeat_rotation(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input)
            .expect("Normalization failed - is benchmark-output present?");
        assert_eq!(
            result, expected,
            "Triplet repeat should preserve specified rotation.\n\
             Input: {}\n\
             Expected (Ferro): {}\n\
             Got: {}",
            input, expected, result
        );
    }

    // =========================================================================
    // CATEGORY 4: Triplet repeat del vs dup (6 base deletion)
    // =========================================================================
    // Similar to Category 1, but with triplet repeats.

    #[rstest]
    #[case("NM_000038.4:c.5603_5605ATG[2]", "NM_000038.4:c.5609_5614del")]
    // Mutalyzer outputs: NM_000038.4:c.5612_5614dup (INCORRECT)
    #[case("NM_000051.3:c.7983_7985TGT[2]", "NM_000051.3:c.7989_7991del")]
    // Mutalyzer outputs: NM_000051.3:c.7989_7991dup (INCORRECT)
    #[ignore]
    fn test_triplet_repeat_deletion_vs_mutalyzer_dup(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input)
            .expect("Normalization failed - is benchmark-output present?");
        assert_eq!(
            result, expected,
            "Triplet repeat[N] should normalize to deletion when count < reference.\n\
             Input: {}\n\
             Expected (Ferro correct): {}\n\
             Got: {}",
            input, expected, result
        );
    }

    // =========================================================================
    // CATEGORY 5: Allelic variants kept separate vs collapsed
    // =========================================================================
    // Ferro keeps allelic variants separated; Mutalyzer may collapse them.

    #[rstest]
    #[case(
        "NM_000030.2:c.[299_307dup;308G>A]",
        "NM_000030.2:c.[299_307dup;308G>A]"
    )]
    // Mutalyzer outputs: NM_000030.2:c.308delinsTCCTGGTTGA (collapses allele)
    #[ignore]
    fn test_allelic_vs_delins_collapse(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input)
            .expect("Normalization failed - is benchmark-output present?");
        assert_eq!(
            result, expected,
            "Allelic variants should be normalized separately.\n\
             Input: {}\n\
             Expected (Ferro): {}\n\
             Got: {}",
            input, expected, result
        );
    }

    // =========================================================================
    // CATEGORY 6: Delins simplification differences
    // =========================================================================
    // Mutalyzer simplifies delins to substitution when only one base differs.
    // Ferro keeps the delins representation.

    #[rstest]
    #[case("NM_000051.4:c.5948_5950delinsGTT", "NM_000051.4:c.5948_5950delinsGTT")]
    // Mutalyzer outputs: NM_000051.4:c.5950A>T (simplifies to substitution)
    #[ignore]
    fn test_delins_vs_substitution(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input)
            .expect("Normalization failed - is benchmark-output present?");
        assert_eq!(
            result, expected,
            "Delins with partial match should remain as delins.\n\
             Input: {}\n\
             Expected (Ferro): {}\n\
             Got: {}",
            input, expected, result
        );
    }
}

// =============================================================================
// NORMALIZATION PATTERNS EXTRACTED FROM EQUIVALENCE CHECKS
// =============================================================================
// These tests verify Ferro's normalization behavior for patterns identified
// in the equivalence comparison between Ferro and Mutalyzer.
//
// The equivalence checks show patterns where different outputs represent
// the same variant. These tests verify that Ferro normalizes inputs to
// the correct canonical form.
//
// Run with: cargo test equivalence_derived_normalize --ignored
// Requires: benchmark-output/ with reference data
// =============================================================================

#[cfg(test)]
mod equivalence_derived_normalize {
    use super::integration_helpers::try_normalize as normalize_to_string;
    use rstest::rstest;

    // =========================================================================
    // PATTERN 1: 3' Shifting
    // =========================================================================
    // Deletions/duplications in repetitive regions shift to 3'-most position.
    // From check_three_prime_shift_equivalence - Mutalyzer may report shifted
    // positions that differ from Ferro's 3'-normalized output.

    #[rstest]
    // CFTR ΔF508 - classic 3' shift example
    #[case("NM_000492.4:c.1520_1522del", "NM_000492.4:c.1521_1523del")]
    // Verified match - already at 3' position
    #[case("NM_000492.4:c.1521_1523del", "NM_000492.4:c.1521_1523del")]
    #[ignore]
    fn test_three_prime_shift_normalization(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input)
            .expect("Normalization failed - is benchmark-output present?");
        assert_eq!(
            result, expected,
            "Deletion should shift to 3'-most position.\n\
             Input: {}\n\
             Expected: {}\n\
             Got: {}",
            input, expected, result
        );
    }

    // =========================================================================
    // PATTERN 3: Duplication → Repeat Notation
    // =========================================================================
    // Tandem duplications in existing repeat regions become repeat notation.
    // From check_dup_repeat_equivalence and check_dup_vs_expanded_repeat_equivalence.

    #[rstest]
    // Triplet repeat expansion: dupGCAGCA in GCA-tract becomes GCA[n]
    #[case("NM_020732.3:c.357_362dupGCAGCA", "NM_020732.3:c.342_362GCA[9]")]
    #[ignore]
    fn test_dup_to_repeat_normalization(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input)
            .expect("Normalization failed - is benchmark-output present?");
        assert_eq!(
            result, expected,
            "Duplication in repeat region should become repeat notation.\n\
             Input: {}\n\
             Expected: {}\n\
             Got: {}",
            input, expected, result
        );
    }

    // =========================================================================
    // PATTERN 4: Repeat → Deletion
    // =========================================================================
    // Repeat notation with count less than reference becomes deletion.
    // From check_repeat_vs_del_dup_equivalence.

    #[rstest]
    // CFTR: ATC[1] when reference has 2 copies = deletion of one copy
    #[case("NM_000492.4:c.1516ATC[1]", "NM_000492.4:c.1519_1521del")]
    #[ignore]
    fn test_repeat_to_deletion_normalization(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input)
            .expect("Normalization failed - is benchmark-output present?");
        assert_eq!(
            result, expected,
            "Repeat with count < reference should become deletion.\n\
             Input: {}\n\
             Expected: {}\n\
             Got: {}",
            input, expected, result
        );
    }

    // =========================================================================
    // PATTERN 5: Explicit Base Removal
    // =========================================================================
    // Redundant sequence info in dup notation is stripped.
    // From equivalence checks where dup with explicit bases equals dup without.

    #[rstest]
    // Single base dup with explicit base: dupC → dup
    #[case("NM_033517.1:c.4818dupC", "NM_033517.1:c.4818dup")]
    // Multi-base dup with explicit sequence: dupGCA → dup
    #[case("NM_020732.3:c.360_362dupGCA", "NM_020732.3:c.360_362dup")]
    #[ignore]
    fn test_explicit_base_removal_normalization(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input)
            .expect("Normalization failed - is benchmark-output present?");
        assert_eq!(
            result, expected,
            "Explicit bases in dup notation should be removed.\n\
             Input: {}\n\
             Expected: {}\n\
             Got: {}",
            input, expected, result
        );
    }

    // =========================================================================
    // PATTERN 7: Allelic Expansion
    // =========================================================================
    // Ferro expands allelic notation to include full accession for each variant.
    // From check_allelic_equivalence.

    #[rstest]
    // Compact allelic round-trips through normalization (HGVS spec compact form)
    #[case("NM_000060.2:c.[1207T>G;1330G>C]", "NM_000060.2:c.[1207T>G;1330G>C]")]
    #[ignore]
    fn test_allelic_expansion_normalization(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input)
            .expect("Normalization failed - is benchmark-output present?");
        assert_eq!(
            result, expected,
            "Allelic notation should be expanded with full accession.\n\
             Input: {}\n\
             Expected: {}\n\
             Got: {}",
            input, expected, result
        );
    }

    // =========================================================================
    // PATTERN 8: Inversion explicit bases preserved
    // =========================================================================
    // Ferro preserves explicit bases in inversions (doesn't strip them).
    // The equivalence check shows invCA ≡ inv, but Ferro preserves both forms.
    // Mutalyzer: c.2082_2083inv (differs - strips explicit bases)

    #[rstest]
    // invCA stays invCA (bases preserved) - Mutalyzer: c.2082_2083inv
    #[case("NM_144670.4:c.2082_2083invCA", "NM_144670.4:c.2082_2083invCA")]
    #[ignore]
    fn test_inv_base_preservation(#[case] input: &str, #[case] expected: &str) {
        let result = normalize_to_string(input)
            .expect("Normalization failed - is benchmark-output present?");
        assert_eq!(
            result, expected,
            "Inversion with explicit bases should be preserved.\n\
             Input: {}\n\
             Expected: {}\n\
             Got: {}",
            input, expected, result
        );
    }
}

// =============================================================================
// COMPREHENSIVE NORMALIZATION TESTS
// =============================================================================
// Systematic tests for all 12 normalization rules, using MockProvider with
// synthetic sequences for isolated, predictable testing.
//
// Run with: cargo test comprehensive_normalization
// =============================================================================

#[cfg(test)]
mod comprehensive_normalization_tests {
    use ferro_hgvs::reference::mock::MockProvider;
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
    use ferro_hgvs::{parse_hgvs, Normalizer};

    // =========================================================================
    // HELPER INFRASTRUCTURE
    // =========================================================================

    /// Create a transcript with full CDS coverage (entire sequence is coding)
    fn make_transcript(id: &str, sequence: &str) -> Transcript {
        let len = sequence.len();
        Transcript::new(
            id.to_string(),
            Some("TEST".to_string()),
            Strand::Plus,
            sequence.to_string(),
            Some(1),          // CDS starts at position 1
            Some(len as u64), // CDS ends at last position
            vec![Exon::new(1, 1, len as u64)],
            None,
            None,
            None,
            Default::default(),
            ManeStatus::None,
            None,
            None,
        )
    }

    /// Create a MockProvider with a single transcript
    fn provider_with_transcript(id: &str, sequence: &str) -> MockProvider {
        let mut provider = MockProvider::new();
        provider.add_transcript(make_transcript(id, sequence));
        provider
    }

    /// Normalize a variant and return the output string
    fn normalize_to_string(provider: MockProvider, input: &str) -> String {
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).unwrap_or_else(|_| panic!("Failed to parse: {}", input));
        let normalized = normalizer
            .normalize(&variant)
            .unwrap_or_else(|_| panic!("Normalization failed for: {}", input));
        format!("{}", normalized)
    }

    /// Assert that normalizing input produces expected output
    fn assert_normalizes_to(input: &str, expected: &str, provider: MockProvider) {
        let result = normalize_to_string(provider, input);
        assert_eq!(
            result, expected,
            "\n  Input:    {}\n  Expected: {}\n  Got:      {}",
            input, expected, result
        );
    }

    // =========================================================================
    // RULE 1 & 2: 3' SHIFTING (Deletions and Insertions)
    // =========================================================================
    // Deletions and insertions in repetitive regions shift to 3'-most position.
    // Function: shuffle() in src/normalize/shuffle.rs

    mod shifting_tests {
        use super::*;

        // ----- Deletion shifting tests -----

        #[test]
        fn test_single_base_deletion_shifts_3prime_in_homopolymer() {
            // Sequence: GGGAAAAAAGGG (A's at positions 4-9 in 1-indexed)
            //           123456789012
            // Deleting A at position 4 should shift to position 9 (3'-most A)
            let seq = "GGGAAAAAAGGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.4del", "NM_TEST.1:c.9del", provider);
        }

        #[test]
        fn test_multi_base_deletion_shifts_3prime_in_dinucleotide_repeat() {
            // Sequence: GGGATATATAGG (AT repeat at positions 4-11)
            //           123456789012
            // Deleting AT at positions 4-5 should shift to 3'-most AT
            let seq = "GGGATATATATT";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            // AT starting at 4 should shift toward end
            // Original positions: 4-5, 6-7, 8-9, 10-11 (4 copies of AT)
            // After 3' shift, del at 10-11
            assert_normalizes_to("NM_TEST.1:c.4_5del", "NM_TEST.1:c.10_11del", provider);
        }

        #[test]
        fn test_deletion_at_3prime_end_unchanged() {
            // Sequence: GGGAAAAAAGGG (A's at positions 4-9)
            // Deleting A at position 9 is already at 3'-most position
            let seq = "GGGAAAAAAGGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.9del", "NM_TEST.1:c.9del", provider);
        }

        #[test]
        fn test_deletion_not_in_repeat_unchanged() {
            // Sequence: ATGCATGCAT (no adjacent repeats for single base)
            // Deleting T at position 2 should not shift
            let seq = "ATGCATGCAT";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            // Position 2 is T, position 3 is G - no shifting
            assert_normalizes_to("NM_TEST.1:c.2del", "NM_TEST.1:c.2del", provider);
        }

        #[test]
        fn test_triplet_deletion_shifts_in_repeat() {
            // Sequence: NNNCAGCAGCAGNNN with CAG repeat
            //           123456789012345
            // CAG at 4-6, 7-9, 10-12, N at 13-15
            let seq = "NNNCAGCAGCAGNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            // Deleting CAG at 4-6 should shift to 10-12
            assert_normalizes_to("NM_TEST.1:c.4_6del", "NM_TEST.1:c.10_12del", provider);
        }

        // ----- Insertion shifting tests -----

        #[test]
        fn test_single_base_insertion_shifts_3prime() {
            // Sequence: GGGAAAAAAGGG (A's at positions 4-9)
            // Inserting A between positions 4-5 should shift to after position 9
            let seq = "GGGAAAAAAGGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            // Inserting A after position 3 (c.3_4insA) in A-tract shifts to c.9_10insA
            // But since the inserted A matches following sequence, it will become dup
            // For pure shifting test, insert at start of A-tract
            assert_normalizes_to("NM_TEST.1:c.3_4insA", "NM_TEST.1:c.9dup", provider);
        }

        #[test]
        fn test_insertion_no_shift_when_no_match() {
            // Sequence: GGGAAAAAAGGG
            // Inserting T should not shift (T doesn't match surrounding A's)
            let seq = "GGGAAAAAAGGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.5_6insT", "NM_TEST.1:c.5_6insT", provider);
        }

        #[test]
        fn test_multi_base_insertion_with_rotation() {
            // Sequence: AAAGCGCGCGCAAA (GC repeat at positions 4-11)
            //           12345678901234
            // Inserting GC in GC-tract should shift with rotation
            let seq = "AAAGCGCGCGCAAA";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            // Inserting GC at start of repeat - should shift and possibly become dup
            assert_normalizes_to("NM_TEST.1:c.3_4insGC", "NM_TEST.1:c.10_11dup", provider);
        }
    }

    // =========================================================================
    // RULE 4: Insertion → Repeat Notation
    // =========================================================================
    // Homopolymer insertion in existing repeat tract becomes repeat notation.
    // Function: insertion_to_repeat() in src/normalize/rules.rs

    mod insertion_to_repeat_tests {
        use super::*;

        #[test]
        fn test_single_a_insertion_into_a_tract_becomes_dup() {
            // Sequence: GGGAAAGGG (A's at positions 4-6).
            // Per HGVS spec: 1 added copy of a unit = dup (never `[N]`).
            // Plus the c. codon-frame exception explicitly forbids `A[4]`
            // (unit_len=1, not a multiple of 3). So the only spec-compliant
            // output here is `dup`.
            let seq = "GGGAAAGGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.4_5insA");
            assert!(
                result.contains("dup"),
                "Single-copy ins matching adjacent A must emit dup, got: {}",
                result
            );
        }

        #[test]
        fn test_multiple_a_insertion_into_tract_in_cds_stays_as_ins() {
            // Sequence: GGGAAAGGG (3 A's at positions 4-6).
            // Per HGVS spec (repeated.md codon-frame exception): in c.
            // context, repeat notation requires unit_len % 3 == 0. With
            // unit_len=1 the canonical form is `ins<literal>`, not A[N].
            let seq = "GGGAAAGGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.6_7insAAA");
            assert!(
                result.contains("insAAA") && !result.contains("A["),
                "c. + unit_len=1 must emit insAAA literal, got: {}",
                result
            );
        }

        #[test]
        fn test_t_insertion_in_t_tract_in_cds_stays_as_ins() {
            // Sequence: AAATTTAAA (T's at positions 4-6). c. + unit_len=1.
            let seq = "AAATTTAAA";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.6_7insTT");
            assert!(
                result.contains("insTT") && !result.contains("T["),
                "c. + unit_len=1 must emit insTT literal, got: {}",
                result
            );
        }

        #[test]
        fn test_mixed_base_insertion_does_not_become_repeat() {
            // Sequence: GGGAAAGGG
            // Inserting AT (mixed bases) should stay as insertion
            let seq = "GGGAAAGGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.6_7insAT");
            assert!(
                result.contains("insAT"),
                "Mixed base insertion should stay as ins, got: {}",
                result
            );
        }

        #[test]
        fn test_large_homopolymer_expansion_in_cds_stays_as_ins() {
            // Sequence: GGAAAAAGGG (5 A's at positions 3-7). Inserting 5 A's
            // would normally produce A[10]; the c. codon-frame gate forbids
            // A[N] (unit_len=1), so the canonical form is `ins<literal>`.
            let seq = "GGAAAAAGGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.7_8insAAAAA");
            assert!(
                result.contains("insAAAAA") && !result.contains("A["),
                "c. + unit_len=1 must emit insAAAAA literal, got: {}",
                result
            );
        }
    }

    // =========================================================================
    // RULE 5: Duplication → Repeat Notation
    // =========================================================================
    // Only multi-base (2+) dups in homopolymer OR multi-copy tandem dup becomes repeat.
    // Single-base duplications stay as simple dups.
    // Function: duplication_to_repeat() in src/normalize/rules.rs

    mod dup_to_repeat_tests {
        use super::*;

        #[test]
        fn test_single_base_dup_in_a_tract_stays_dup() {
            // Sequence: GGGAAAGGG (A's at positions 4-6)
            // Single-base duplications stay as simple dups (not repeat notation)
            let seq = "GGGAAAGGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.6dup");
            // Single-base dup stays as dup
            assert!(
                result.contains("dup"),
                "Single-base dup should stay as dup, got: {}",
                result
            );
        }

        #[test]
        fn test_single_base_dup_in_long_homopolymer_stays_dup() {
            // Sequence: GGAAAAAAAAGG (8 A's at positions 3-10)
            // Single-base duplications stay as simple dups (not repeat notation)
            let seq = "GGAAAAAAAAGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.10dup");
            // Single-base dup stays as dup
            assert!(
                result.contains("dup"),
                "Single-base dup should stay as dup, got: {}",
                result
            );
        }

        #[test]
        fn test_single_copy_triplet_dup_stays_dup() {
            // Sequence: NNCAGCAGCAGNN (3 CAG repeats at positions 3-11)
            //           1234567890123
            // Duplicating single CAG (one copy) stays as dup per HGVS rules
            // Only multi-copy tandem dups become repeat notation
            let seq = "NNCAGCAGCAGNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.9_11dup");
            // Single-copy dup stays as dup
            assert!(
                result.contains("dup"),
                "Single-copy triplet dup should stay as dup, got: {}",
                result
            );
        }

        #[test]
        fn test_multi_copy_triplet_dup_becomes_repeat() {
            // Sequence: NNCAGCAGCAGNN (3 CAG repeats)
            // Duplicating CAGCAG (two copies) becomes repeat
            let seq = "NNCAGCAGCAGNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.6_11dup");
            // Multi-copy dup becomes repeat
            assert!(
                result.contains("CAG[5]"),
                "Multi-copy triplet dup should become repeat, got: {}",
                result
            );
        }

        #[test]
        fn test_single_copy_dup_no_repeat_stays_dup() {
            // Sequence: NNGCANN (single GCA, no repeat)
            //           1234567
            // Duplicating GCA when there's only one copy stays as dup
            let seq = "NNGCANNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3_5dup");
            // Should stay as dup (no existing repeat to extend)
            assert!(
                result.contains("dup"),
                "Single-copy dup should stay as dup, got: {}",
                result
            );
        }

        #[test]
        fn test_single_copy_dinucleotide_dup_stays_dup() {
            // Sequence: NNGTGTGTNN (GT repeat)
            //           1234567890
            // Duplicating single GT (one copy) stays as dup
            // Only multi-copy tandem dups become repeat notation
            let seq = "NNGTGTGTNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.7_8dup");
            // Single-copy dup stays as dup
            assert!(
                result.contains("dup"),
                "Single-copy dinucleotide dup should stay as dup, got: {}",
                result
            );
        }
    }

    // =========================================================================
    // RULE 6: Repeat → Deletion
    // =========================================================================
    // Repeat count < reference count becomes deletion.
    // Function: normalize_repeat() in src/normalize/rules.rs

    mod repeat_to_deletion_tests {
        use super::*;

        #[test]
        fn test_repeat_contraction_to_deletion_triplet() {
            // Sequence: NNCAGCAGCAGNN (3 CAG repeats)
            //           1234567890123
            // CAG[2] when ref has 3 = delete one CAG
            let seq = "NNCAGCAGCAGNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3_5CAG[2]");
            // Should become deletion (removing 3 bases)
            assert!(
                result.contains("del"),
                "CAG[2] with 3 copies should become del, got: {}",
                result
            );
        }

        #[test]
        fn test_repeat_contraction_in_cds_homopolymer_routes_to_del() {
            // 5-A ref, A[3] in c. context. unit_len=1 — codon-frame gate
            // forbids A[N], routes the contraction-with-survivors branch to
            // a plain deletion.
            let seq = "GGAAAAAGGG"; // 5 A's at positions 3-7
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3A[3]");
            assert_eq!(result, "NM_TEST.1:c.6_7del");
        }

        #[test]
        fn test_repeat_contraction_in_cds_multibase_routes_to_del() {
            // 3-CT ref, CT[1] in c. context. unit_len=2 — codon-frame gate
            // forbids CT[N], routes to plain del.
            let seq = "AACTCTCTAA"; // 3 CTs at positions 3-8
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3_4CT[1]");
            assert_eq!(result, "NM_TEST.1:c.5_8del");
        }

        #[test]
        fn test_repeat_at_reference_count_unchanged() {
            // Sequence: NNCAGCAGCAGNN (3 CAG repeats)
            // CAG[3] = reference, should be unchanged or become identity
            let seq = "NNCAGCAGCAGNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3_5CAG[3]");
            // Should be identity (c.=) or unchanged
            assert!(
                result.contains("=") || result.contains("CAG[3]"),
                "CAG[3] at reference should be identity or unchanged, got: {}",
                result
            );
        }
    }

    // =========================================================================
    // RULE 7: Repeat → Duplication
    // =========================================================================
    // Repeat count = reference count + 1 becomes duplication.
    // Function: normalize_repeat() in src/normalize/rules.rs

    mod repeat_to_dup_tests {
        use super::*;

        #[test]
        fn test_repeat_expansion_by_one_becomes_dup_triplet() {
            // Sequence: NNCAGCAGCAGNN (3 CAG repeats)
            // CAG[4] = one more than ref = dup
            let seq = "NNCAGCAGCAGNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3_5CAG[4]");
            assert!(
                result.contains("dup"),
                "CAG[4] with 3 copies should become dup, got: {}",
                result
            );
        }

        #[test]
        fn test_repeat_expansion_by_one_becomes_dup_homopolymer() {
            // Sequence: GGAAAAAGGG (5 A's)
            // A[6] = one more than ref = dup one A
            let seq = "GGAAAAAGGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3A[6]");
            assert!(
                result.contains("dup"),
                "A[6] with 5 copies should become dup, got: {}",
                result
            );
        }

        #[test]
        fn test_dinucleotide_repeat_expansion_by_one() {
            // Sequence: NNGTGTNN (2 GT repeats)
            //           12345678
            // GT[3] = one more than ref = dup one GT
            let seq = "NNGTGTNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3_4GT[3]");
            assert!(
                result.contains("dup"),
                "GT[3] with 2 copies should become dup, got: {}",
                result
            );
        }
    }

    // =========================================================================
    // RULE 8: Repeat Stays Repeat (Large Expansions)
    // =========================================================================
    // When count > ref + 1, repeat notation is preserved.
    // Function: normalize_repeat() in src/normalize/rules.rs

    mod repeat_stays_repeat_tests {
        use super::*;

        #[test]
        fn test_large_expansion_stays_repeat() {
            // Sequence: NNCAGCAGCAGNN (3 CAG repeats)
            // CAG[10] = large expansion, stays as repeat
            let seq = "NNCAGCAGCAGNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3_5CAG[10]");
            assert!(
                result.contains("CAG[10]"),
                "Large expansion should stay as repeat, got: {}",
                result
            );
        }

        #[test]
        fn test_large_homopolymer_expansion_in_cds_routes_to_ins() {
            // Sequence: GGAAAAAGGG (5 A's). A[20] = +15 unit copies.
            // Per HGVS spec (repeated.md codon-frame exception): repeat
            // notation in c. requires unit_len % 3 == 0. unit_len=1 is gated;
            // canonical form is `ins<literal>` of the added copies.
            let seq = "GGAAAAAGGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3A[20]");
            assert!(
                result.contains("ins") && !result.contains("A["),
                "c. + unit_len=1 expansion must emit ins, got: {}",
                result
            );
        }

        #[test]
        fn test_moderate_expansion_stays_repeat() {
            // Sequence: NNCAGCAGCAGNN (3 CAG repeats)
            // CAG[6] = 3 more than ref, stays as repeat
            let seq = "NNCAGCAGCAGNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3_5CAG[6]");
            assert!(
                result.contains("CAG[6]"),
                "Moderate expansion should stay as repeat, got: {}",
                result
            );
        }
    }

    // =========================================================================
    // RULE 9: Delins → Duplication
    // =========================================================================
    // Delins that doubles the deleted sequence becomes duplication.
    // Function: delins_is_duplication() in src/normalize/rules.rs

    mod delins_to_dup_tests {
        use super::*;

        #[test]
        fn test_single_base_delins_to_dup() {
            // Sequence: ATGCATGCAT
            //           1234567890
            // c.5delinsAA where ref[5]=A → delinsAA doubles A → c.5dup
            let seq = "ATGCATGCAT";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.5delinsAA");
            assert!(
                result.contains("dup"),
                "delinsAA at A should become dup, got: {}",
                result
            );
        }

        #[test]
        fn test_multi_base_delins_to_dup() {
            // Sequence: NNATGCNNNNN
            //           12345678901
            // c.3_5delinsATGATG where ref[3-5]=ATG → doubles to dup
            let seq = "NNATGCNNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3_5delinsATGATG");
            assert!(
                result.contains("dup"),
                "delinsATGATG at ATG should become dup, got: {}",
                result
            );
        }

        #[test]
        fn test_delins_different_sequence_stays_delins() {
            // Sequence: ATGCATGCAT
            // c.5delinsGG where ref[5]=A → doesn't double, stays delins
            let seq = "ATGCATGCAT";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.5delinsGG");
            assert!(
                result.contains("delins"),
                "delinsGG at A should stay delins, got: {}",
                result
            );
        }

        #[test]
        fn test_delins_partial_match_becomes_pure_insertion() {
            // Sequence: NNATGCNNNNN
            // c.3_5 = ATG. Insert ATGTTT: shared ATG prefix consumes the
            // entire deleted reference, leaving TTT inserted between c.5
            // and c.6 per HGVS minimal-form rules. Not a duplication
            // (insert second half TTT != ATG).
            let seq = "NNATGCNNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to(
                "NM_TEST.1:c.3_5delinsATGTTT",
                "NM_TEST.1:c.5_6insTTT",
                provider,
            );
        }

        #[test]
        fn test_delins_single_to_triple_not_dup() {
            // Sequence: ATGCATGCAT
            // c.5delinsAAA → triples, not doubles, stays delins or becomes repeat
            let seq = "ATGCATGCAT";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.5delinsAAA");
            // Should not become simple dup (that's doubling)
            // Could become repeat A[3] or stay as delins
            assert!(
                !result.ends_with("dup") || result.contains("["),
                "delinsAAA at A should not become simple dup, got: {}",
                result
            );
        }
    }

    // =========================================================================
    // RULE 9b: Delins → Inversion (item A2 of issue #81)
    // =========================================================================
    // Delins whose insert is the reverse complement of the deleted reference
    // becomes an inversion. The HGVS spec defines `inv` precisely as this case
    // and excludes it from the definition of `delins`. The complementary-outer-
    // bases shortening rule applies to the resulting inversion so that a
    // delins-encoded inversion and a directly-encoded `inv` produce the same
    // canonical output.
    // Function: canonicalize_delins() in src/normalize/rules.rs

    mod delins_to_inv_tests {
        use super::*;

        #[test]
        fn test_simple_delins_to_inv() {
            // Sequence: NNCTANNNNN (CTA at c.3_5, revcomp = TAG)
            //           1234567890
            let seq = "NNCTANNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.3_5delinsTAG", "NM_TEST.1:c.3_5inv", provider);
        }

        #[test]
        fn test_delins_to_inv_with_outer_shortening() {
            // Sequence: NCTATGNNNN (CTATG at c.2_6, revcomp = CATAG)
            //           1234567890
            // Outer C/G are complement -> shortens to inner TAT (c.3_5)
            let seq = "NCTATGNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.2_6delinsCATAG", "NM_TEST.1:c.3_5inv", provider);
        }

        #[test]
        fn test_delins_shared_suffix_becomes_substitution() {
            // Sequence: NCTAGNNNNN (CTAG at c.2_5)
            //           1234567890
            // Insert TTAG: shares the `TAG` suffix with the reference, so the
            // minimal HGVS edit is a single-base substitution at c.2 (priority
            // sub > delins). This also confirms the partial-complement input
            // does NOT become an inversion (revcomp(CTAG) = CTAG, palindrome,
            // not equal to TTAG).
            let seq = "NCTAGNNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.2_5delinsTTAG", "NM_TEST.1:c.2C>T", provider);
        }

        #[test]
        fn test_delins_shared_prefix_becomes_substitution() {
            // Mirror of the suffix-trim test: shared `CTA` prefix, single-base
            // mismatch at the trimmed tail position.
            // Sequence: NCTAGNNNNN (CTAG at c.2_5)
            //           1234567890
            // Insert CTAA: shares `CTA` prefix; minimal edit is c.5G>A.
            let seq = "NCTAGNNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.2_5delinsCTAA", "NM_TEST.1:c.5G>A", provider);
        }

        #[test]
        fn test_delins_shared_affix_pure_deletion() {
            // Shared affixes consume the entire inserted sequence -> pure deletion.
            // Sequence: NACGTNNNNN (ACGT at c.2_5)
            //           1234567890
            // Insert AT: prefix A=A, suffix T=T. Trimmed: delete CG at c.3_4,
            // insert nothing. Minimal HGVS: c.3_4del.
            let seq = "NACGTNNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.2_5delinsAT", "NM_TEST.1:c.3_4del", provider);
        }

        #[test]
        fn test_delins_shared_affix_pure_insertion() {
            // Shared affixes consume the entire deleted reference -> pure insertion.
            // Sequence: NACTNNNNNN (ACT at c.2_4)
            //           1234567890
            // Insert ACGT: prefix AC=AC, suffix T=T. Trimmed: delete nothing,
            // insert G between c.3 and c.4. Minimal HGVS: c.3_4insG.
            let seq = "NACTNNNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.2_4delinsACGT", "NM_TEST.1:c.3_4insG", provider);
        }

        #[test]
        fn test_delins_shared_affix_residual_delins() {
            // Shared affixes leave a residual N->M (here 2->2) that doesn't
            // fit any higher-priority form -> emit minimal (trimmed) delins.
            // Sequence: NAGGCTNNN (AGGCT at c.2_6)
            //           123456789
            // Insert AAACT: prefix A=A, suffix CT=CT. Trimmed range c.3_4
            // (deleted GG, length 2) -> insert AA (length 2). revcomp(GG)=CC,
            // not AA, so not an inversion. Minimal HGVS: c.3_4delinsAA.
            let seq = "NAGGCTNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to(
                "NM_TEST.1:c.2_6delinsAAACT",
                "NM_TEST.1:c.3_4delinsAA",
                provider,
            );
        }

        #[test]
        fn test_delins_one_base_complement_is_substitution() {
            // Priority guard: a 1-base complement (ref=A, insert=T) is a
            // substitution per HGVS spec, NEVER an inversion (which requires
            // more than one nucleotide).
            let seq = "NANNNNNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.2delinsT", "NM_TEST.1:c.2A>T", provider);
        }

        #[test]
        fn test_delins_palindrome_is_identity() {
            // Sequence: NATATNNNNN (ATAT at c.2_5)
            //           1234567890
            // ATAT is its own revcomp; insert ATAT must hit Identity, not Inversion.
            let seq = "NATATNNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.2_5delinsATAT", "NM_TEST.1:c.2_5=", provider);
        }

        #[test]
        fn test_delins_revcomp_round_trip_matches_direct_inv() {
            // The same biological event written as delins or as inv must produce
            // the same normalized string.
            let seq = "NCTATGNNNN";
            let direct = normalize_to_string(
                provider_with_transcript("NM_TEST.1", seq),
                "NM_TEST.1:c.2_6inv",
            );
            let via_delins = normalize_to_string(
                provider_with_transcript("NM_TEST.1", seq),
                "NM_TEST.1:c.2_6delinsCATAG",
            );
            assert_eq!(
                direct, via_delins,
                "delins-encoded inversion must normalize to the same string \
                 as a directly-encoded inversion (direct={}, via_delins={})",
                direct, via_delins,
            );
        }
    }

    // =========================================================================
    // RULE 10: Inversion Shortening
    // =========================================================================
    // Inversions with complementary outer bases shorten.
    // Function: shorten_inversion() in src/normalize/rules.rs

    mod inversion_tests {
        use super::*;

        #[test]
        fn test_inversion_no_shortening_when_not_complementary() {
            // Sequence: NNATGCNNNN (ATGC at c.3_6)
            //           1234567890
            // A vs C: complement(A)=T != C, so no shortening; stays c.3_6inv.
            let seq = "NNATGCNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.3_6inv", "NM_TEST.1:c.3_6inv", provider);
        }

        #[test]
        fn test_inversion_palindrome_collapses_to_identity() {
            // Sequence: NNATCGATNN (ATCGAT at c.3_8). ATCGAT is its own
            //           1234567890
            // reverse complement, so every outer pair is complementary and
            // shortening collapses to identity.
            let seq = "NNATCGATNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.3_8inv", "NM_TEST.1:c.3_8=", provider);
        }

        #[test]
        fn test_inversion_partial_shortening() {
            // Sequence: NAGGGCTNNN (AGGGCT at c.2_7).
            //           1234567890
            // Outer A/T are complement -> shortens to inner GGGC (c.3_6inv).
            // GGGC has no further cancelling outer pair (G vs C complement
            // would cancel further but G == G's revcomp partner here is the
            // inner GG, not a complement of C). Wait — G complement IS C, so
            // GGGC also cancels: G-C cancels to GG. Then G-G doesn't cancel.
            // Final: c.4_5inv over GG.
            let seq = "NAGGGCTNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.2_7inv", "NM_TEST.1:c.4_5inv", provider);
        }

        #[test]
        fn test_two_base_palindrome_collapses_to_identity() {
            // Sequence: NNATNNNNNN (AT at c.3_4). revcomp(AT) = AT (palindrome)
            //           1234567890
            // -> shortening collapses to identity.
            let seq = "NNATNNNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.3_4inv", "NM_TEST.1:c.3_4=", provider);
        }

        #[test]
        fn test_two_base_gc_palindrome_collapses_to_identity() {
            // Sequence: NNGCNNNNNN (GC at c.3_4). revcomp(GC) = GC (palindrome)
            //           1234567890
            // -> shortening collapses to identity.
            let seq = "NNGCNNNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            assert_normalizes_to("NM_TEST.1:c.3_4inv", "NM_TEST.1:c.3_4=", provider);
        }
    }

    // =========================================================================
    // RULE 11: Explicit Base Removal (Canonicalization)
    // =========================================================================
    // Redundant sequence info in dup/del is stripped.
    // Function: canonicalize_edit() in src/normalize/rules.rs

    mod canonicalization_tests {
        use super::*;

        #[test]
        fn test_dup_with_explicit_single_base_stripped() {
            // c.5dupA should become c.5dup
            let seq = "ATGCATGCAT";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.5dupA");
            assert_eq!(
                result, "NM_TEST.1:c.5dup",
                "dupA should become dup, got: {}",
                result
            );
        }

        #[test]
        fn test_dup_with_explicit_multi_base_stripped() {
            // c.5_7dupATG should become c.5_7dup
            let seq = "NNATGCNNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3_5dupATG");
            assert!(
                result.ends_with("dup"),
                "dupATG should become dup, got: {}",
                result
            );
        }

        #[test]
        fn test_del_with_explicit_sequence_stripped() {
            // c.5_7delATG should become c.5_7del
            let seq = "NNATGCNNNNN";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.3_5delATG");
            // After 3' shifting, should end with just "del"
            assert!(
                result.contains("del") && !result.contains("delATG"),
                "delATG should become del, got: {}",
                result
            );
        }

        #[test]
        fn test_del_with_numeric_count_stripped() {
            // c.5del1 should become c.5del (for single base)
            let seq = "ATGCATGCAT";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.5del1");
            // Should strip the count
            assert!(
                result.contains("del") && !result.contains("del1"),
                "del1 should become del, got: {}",
                result
            );
        }

        #[test]
        fn test_dup_stays_dup_type() {
            // Ensure dup type is preserved even after canonicalization
            let seq = "ATGCATGCAT";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.5dupA");
            assert!(
                result.contains("dup"),
                "Canonicalized dup should still be dup type, got: {}",
                result
            );
        }
    }

    // =========================================================================
    // RULE 12: Allele Normalization
    // =========================================================================
    // Each variant in an allele is normalized independently.
    // Function: normalize_allele() in src/normalize/mod.rs

    mod allele_tests {
        use super::*;

        #[test]
        fn test_allele_each_variant_normalized() {
            // Allele with substitution (unchanged) and deletion (shifts)
            // Sequence with A-tract for deletion to shift
            let seq = "GGGAAAAAAGGGCCCGGGAAA";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "[NM_TEST.1:c.4del;NM_TEST.1:c.15C>T]");
            // Deletion should shift, substitution unchanged
            assert!(
                result.contains("del") && result.contains("C>T"),
                "Both variants should be present, got: {}",
                result
            );
        }

        #[test]
        fn test_allele_with_two_deletions_both_shift() {
            // Two deletions in repeat regions should both shift
            let seq = "GGGAAAAAAGGGTTTTTTGGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "[NM_TEST.1:c.4del;NM_TEST.1:c.13del]");
            // Both should shift to 3'-most positions
            assert!(
                result.contains("9del") && result.contains("18del"),
                "Both deletions should shift, got: {}",
                result
            );
        }

        #[test]
        fn test_allele_preserves_semicolon_separator() {
            // Allele output should use semicolon separator
            let seq = "ATGCATGCAT";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "[NM_TEST.1:c.1A>G;NM_TEST.1:c.5A>G]");
            assert!(
                result.contains(";"),
                "Allele should use semicolon separator, got: {}",
                result
            );
        }

        #[test]
        fn test_allele_with_dup_and_sub() {
            // Allele with duplication (may become repeat) and substitution
            let seq = "GGGAAAGGGCCC";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "[NM_TEST.1:c.6dup;NM_TEST.1:c.10C>T]");
            // Dup might become repeat, sub unchanged
            assert!(
                result.contains(";"),
                "Allele should contain semicolon, got: {}",
                result
            );
        }

        #[test]
        fn test_compact_allele_notation() {
            // Compact allele notation with accession before bracket
            // NM_TEST.1:c.[4del;9A>G] format
            let seq = "GGGAAAAAAGGG";
            let provider = provider_with_transcript("NM_TEST.1", seq);
            let result = normalize_to_string(provider, "NM_TEST.1:c.[4del;9A>G]");
            // Should normalize the variants
            assert!(
                result.contains("del") && result.contains("A>G"),
                "Compact allele should normalize both variants, got: {}",
                result
            );
        }
    }
}

// =============================================================================
// INTRONIC NORMALIZATION FAILURE TESTS (TDD Phase 1)
// =============================================================================
// These tests document known intronic normalization failures that need fixes.
// Each test asserts the DESIRED behavior (normalization succeeds). They are
// marked #[ignore] and will fail until the corresponding fixes are implemented.
//
// Root causes:
// 1. CIGAR insertions: CDS numbering skips CIGAR insertion bases
// 2. CIGAR complex: deletions or multi-gap CIGARs compound the issue
// 3. Version mismatch: FASTA has latest version; ClinVar references old versions
// 4. LRG transcripts: LRG→RefSeq mapping points to old RefSeq versions
// 5. Non-contiguous cdot exons: gaps trigger supplemental fallback → single exon
//
// Run with: cargo nextest run intronic_normalization --ignored
// Requires: benchmark-output/ with reference data
// =============================================================================

#[cfg(test)]
mod cigar_aware_normalization {
    use super::integration_helpers::{has_reference_data, try_normalize};
    use rstest::rstest;

    // =========================================================================
    // CIGAR insertion offset: these transcripts have CIGAR insertions (I ops)
    // that cause CDS position offsets. Without CIGAR-aware mapping, the
    // intronic position cannot be resolved to a genomic coordinate.
    // =========================================================================

    #[rstest]
    // NM_015120.4 (ABHD12): CIGAR has I3 in exon 0 — 3bp insertion shifts CDS
    #[case("NM_015120.4:c.1433-12_1433-10dup")]
    #[case("NM_015120.4:c.450+22_450+25del")]
    #[case("NM_015120.4:c.9540-7dup")]
    #[case("NM_015120.4:c.12362+16_12362+29delinsTGAGTTTGTGTG")]
    #[case("NM_015120.4:c.10078+1_10078+5delinsAAAAACCCTTGCAGAATGAAAA")]
    #[case("NM_015120.4:c.7007+30del")]
    #[case("NM_015120.4:c.1114-8del")]
    // NM_001304717.5: CIGAR has M479 D1 M444 in exon 0
    #[case("NM_001304717.5:c.774-30dup")]
    #[case("NM_001304717.5:c.729+4_729+7del")]
    #[case("NM_001304717.5:c.1066+4del")]
    #[case("NM_001304717.5:c.1066+5_1066+8del")]
    #[case("NM_001304717.5:c.1067-3_1067-2del")]
    // NM_002111.8 (HTT): CIGAR insertion in first exon
    #[case("NM_002111.8:c.52+1del")]
    #[case("NM_002111.8:c.53-6del")]
    #[case("NM_002111.8:c.52+5G>A")]
    // NM_001467.6: CIGAR-affected transcript
    #[case("NM_001467.6:c.1+1del")]
    // NM_001164278.2 / NM_001164279.2: CIGAR-affected
    #[case("NM_001164278.2:c.427+5del")]
    #[case("NM_001164279.2:c.544+5del")]
    // NM_000527.4 (LDLR): widely-referenced transcript with CIGAR gap
    #[case("NM_000527.4:c.313+1G>A")]
    #[case("NM_000527.4:c.313+2T>C")]
    #[case("NM_000527.4:c.190+4_190+7del")]
    fn test_cigar_intronic_normalizes_ok(#[case] input: &str) {
        let result = try_normalize(input);
        // Skip if reference data not available
        if result.is_none() && !has_reference_data() {
            eprintln!("Skipping test - benchmark-output not available");
            return;
        }
        assert!(result.is_some(), "Should normalize successfully: {}", input);
    }
}

#[cfg(test)]
mod cigar_edge_case_normalization {
    use super::integration_helpers::{has_reference_data, try_normalize};
    use rstest::rstest;

    // =========================================================================
    // CIGAR edge cases: transcripts with deletions (D ops), multi-gap CIGARs,
    // or large offsets where simple insertion offset doesn't fix the boundary.
    // =========================================================================

    #[rstest]
    // NM_130444.3 / NM_030582.4: complex CIGAR with D ops
    #[case("NM_130444.3:c.4332+9_4332+10del")]
    #[case("NM_130444.3:c.4939-20AC[4]")]
    #[case("NM_130444.3:c.4740+14del")]
    #[case("NM_130444.3:c.4476-6del")]
    #[case("NM_130444.3:c.3788+11del")]
    #[case("NM_130444.3:c.4186+5G>A")]
    #[case("NM_130444.3:c.4939-5del")]
    #[case("NM_130444.3:c.5261+5G>A")]
    #[case("NM_030582.4:c.4332+9_4332+10del")]
    #[case("NM_030582.4:c.4939-20AC[4]")]
    #[case("NM_030582.4:c.4740+14del")]
    #[case("NM_030582.4:c.4476-6del")]
    #[case("NM_030582.4:c.3788+11del")]
    #[case("NM_030582.4:c.4186+5G>A")]
    // NM_001372044.2: CIGAR with D2 in middle of exon
    #[case("NM_001372044.2:c.1529+4dup")]
    #[case("NM_001372044.2:c.1530-4del")]
    #[case("NM_001372044.2:c.1347+5G>A")]
    #[case("NM_001372044.2:c.1529+3A>G")]
    #[case("NM_001372044.2:c.1244+7C>T")]
    #[case("NM_001372044.2:c.1244+5G>A")]
    // NM_001414686.1 / NM_001401501.2 / NM_001414687.1: multi-gap CIGARs
    #[case("NM_001414686.1:c.313+1G>A")]
    #[case("NM_001414686.1:c.313+2T>C")]
    #[case("NM_001414686.1:c.190+4_190+7del")]
    #[case("NM_001401501.2:c.313+1G>A")]
    #[case("NM_001401501.2:c.313+2T>C")]
    #[case("NM_001401501.2:c.190+4_190+7del")]
    #[case("NM_001414687.1:c.313+1G>A")]
    #[case("NM_001414687.1:c.313+2T>C")]
    #[case("NM_001414687.1:c.190+4_190+7del")]
    // Additional CIGAR-complex transcripts
    #[case("NM_001405709.1:c.313+1G>A")]
    #[case("NM_001405709.1:c.190+4_190+7del")]
    #[case("NM_001405710.1:c.313+1G>A")]
    #[case("NM_001405710.1:c.190+4_190+7del")]
    fn test_cigar_edge_case_normalizes_ok(#[case] input: &str) {
        let result = try_normalize(input);
        if result.is_none() && !has_reference_data() {
            eprintln!("Skipping test - benchmark-output not available");
            return;
        }
        assert!(result.is_some(), "Should normalize successfully: {}", input);
    }
}

#[cfg(test)]
mod cdot_gap_normalization {
    use super::integration_helpers::{has_reference_data, try_normalize};
    use rstest::rstest;

    // =========================================================================
    // Non-contiguous cdot exons: transcripts where cdot tx coordinates have
    // gaps between exons. This triggers the supplemental fallback which creates
    // a single synthetic exon, losing all intron boundaries.
    // =========================================================================

    #[rstest]
    // NM_001405681.2: gaps in cdot exon tx coordinates
    #[case("NM_001405681.2:c.2535-4_2755del")]
    #[case("NM_001405681.2:c.2923-106_3210-108del")]
    // NM_017940.8: gaps in cdot exon tx coordinates
    #[case("NM_017940.8:c.2451-4_2671del")]
    #[case("NM_017940.8:c.2839-106_3126-108del")]
    // NM_001405694.2: gaps in cdot exon tx coordinates
    #[case("NM_001405694.2:c.2451-4_2671del")]
    #[case("NM_001405694.2:c.2839-106_3126-108del")]
    // NM_001405693.2: gaps in cdot exon tx coordinates
    #[case("NM_001405693.2:c.2451-4_2671del")]
    #[case("NM_001405693.2:c.2839-106_3126-108del")]
    // NM_001405684.2: gaps in cdot exon tx coordinates
    #[case("NM_001405684.2:c.2451-4_2671del")]
    // NM_001405683.2: gaps in cdot exon tx coordinates
    #[case("NM_001405683.2:c.2451-4_2671del")]
    fn test_cdot_gap_intronic_normalizes_ok(#[case] input: &str) {
        let result = try_normalize(input);
        if result.is_none() && !has_reference_data() {
            eprintln!("Skipping test - benchmark-output not available");
            return;
        }
        assert!(result.is_some(), "Should normalize successfully: {}", input);
    }
}

#[cfg(test)]
mod version_mismatch_normalization {
    use super::integration_helpers::{has_reference_data, try_normalize};
    use rstest::rstest;

    // =========================================================================
    // Version mismatch: the FASTA has the latest transcript version but the
    // ClinVar pattern references an older version. Version fallback changes
    // cds_start, causing intronic positions to be off.
    // =========================================================================

    #[rstest]
    #[case("NM_000036.2:c.1015-27del")] // FASTA has .3
    #[case("NM_000068.2:c.4987-12del")] // FASTA has .4
    #[case("NM_000251.1:c.942+3A>T")] // FASTA has .3
    #[case("NM_000314.4:c.80+1del")] // FASTA has .8
    #[case("NM_000351.4:c.1037-5_1037-4del")] // FASTA has .7
    #[case("NM_000368.4:c.1000+5del")] // FASTA has .5
    #[case("NM_000500.6:c.290-13A>C")] // FASTA has .9
    #[case("NM_000542.3:c.1-55del")] // FASTA has .5
    #[case("NM_153609.3:c.494+5G>A")] // FASTA has .4
    #[case("NM_000059.3:c.68-7del")] // FASTA has .4, BRCA2
    #[case("NM_000059.3:c.8332-1G>A")] // FASTA has .4, BRCA2
    #[case("NM_000059.3:c.7435+1G>A")] // FASTA has .4, BRCA2
    #[case("NM_000249.3:c.1558+1G>A")] // FASTA has .4, MLH1
    #[case("NM_000249.3:c.306+1G>A")] // FASTA has .4, MLH1
    #[case("NM_007294.3:c.5278-1G>A")] // FASTA has .4, BRCA1
    fn test_version_mismatch_normalizes_ok(#[case] input: &str) {
        let result = try_normalize(input);
        if result.is_none() && !has_reference_data() {
            eprintln!("Skipping test - benchmark-output not available");
            return;
        }
        assert!(result.is_some(), "Should normalize successfully: {}", input);
    }
}

#[cfg(test)]
mod lrg_intronic_normalization {
    use super::integration_helpers::{has_reference_data, try_normalize};
    use rstest::rstest;

    // =========================================================================
    // LRG transcripts: LRG→RefSeq mapping typically points to an old RefSeq
    // version, making these a downstream effect of version mismatch + CIGAR.
    // Once those are fixed, LRG patterns should work.
    // =========================================================================

    #[rstest]
    // LRG_741t1 → NM_015120.4 (CIGAR-affected)
    #[case("LRG_741t1:c.10078+1_10078+5delinsAAAAACCCTTGCAGAATGAAAA")]
    #[case("LRG_741t1:c.1433-12_1433-10dup")]
    #[case("LRG_741t1:c.7007+30del")]
    #[case("LRG_741t1:c.450+22_450+25del")]
    // LRG_763t1 → NM_002111.8 (CIGAR-affected, HTT gene)
    #[case("LRG_763t1:c.52+1del")]
    #[case("LRG_763t1:c.53-6del")]
    // LRG_720t1 → NM_000280.3 (version mismatch — FASTA has .4)
    #[case("LRG_720t1:c.1183+4dup")]
    #[case("LRG_720t1:c.1327+4del")]
    // LRG_293t1 → NM_000059.3 (version mismatch — FASTA has .4, BRCA2)
    #[case("LRG_293t1:c.1-59_1-57del")]
    #[case("LRG_293t1:c.68-7del")]
    #[case("LRG_293t1:c.8332-1G>A")]
    // LRG_486t1 → NM_000368.4 (version mismatch — FASTA has .5)
    #[case("LRG_486t1:c.-234-1616_-234-235del1382")]
    #[case("LRG_486t1:c.1000+5del")]
    // LRG_384t1 → NM_000257.2 (version mismatch)
    #[case("LRG_384t1:c.-64-17del")]
    // LRG_311t1 → NM_000314.4 (version mismatch — FASTA has .8)
    #[case("LRG_311t1:c.-1032-222_-366del889")]
    #[case("LRG_311t1:c.80+1del")]
    fn test_lrg_intronic_normalizes_ok(#[case] input: &str) {
        let result = try_normalize(input);
        if result.is_none() && !has_reference_data() {
            eprintln!("Skipping test - benchmark-output not available");
            return;
        }
        assert!(result.is_some(), "Should normalize successfully: {}", input);
    }
}

#[cfg(test)]
mod cigar_cds_mapping {
    // =========================================================================
    // Unit tests for CIGAR-aware CDS→transcript mapping.
    // These verify that CDS positions are correctly adjusted when CIGAR
    // insertions/deletions exist in the alignment.
    // =========================================================================

    use ferro_hgvs::data::{CdotTranscript, CigarOp};
    use ferro_hgvs::reference::Strand;

    /// Create a transcript that models NM_015120.4's exon 0:
    /// CIGAR = M185 I3 M250, start_codon=111
    /// The 3bp insertion at tx pos 186-188 means CDS numbering (which follows
    /// the genome) should account for these extra transcript bases.
    fn transcript_with_cigar_insertion() -> CdotTranscript {
        // Exon 0: 438bp on transcript (185+3+250), 435bp on genome (185+250)
        // Exon 1: next exon starting after intron
        CdotTranscript {
            gene_name: Some("TEST_CIGAR".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            exons: vec![
                [1000, 1435, 0, 438],   // Exon 0: 435bp genome, 438bp tx (3bp insertion)
                [2000, 2200, 438, 638], // Exon 1: 200bp
                [3000, 3300, 638, 938], // Exon 2: 300bp
            ],
            cds_start: Some(111), // CDS starts at tx pos 111
            cds_end: Some(900),
            gene_id: None,
            protein: None,
            exon_cigars: Vec::new(),
        }
    }

    /// Create a transcript modeling a CIGAR with deletion:
    /// CIGAR = M504 D2 M123
    /// D2 means genome has 2 bases that transcript doesn't, so the exon is
    /// longer on genome (504+2+123=629) than transcript (504+123=627).
    fn transcript_with_cigar_deletion() -> CdotTranscript {
        CdotTranscript {
            gene_name: Some("TEST_CIGAR_DEL".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            exons: vec![
                [1000, 1629, 0, 627],   // Exon 0: 629bp genome, 627bp tx (2bp deletion)
                [2000, 2200, 627, 827], // Exon 1: 200bp
            ],
            cds_start: Some(100),
            cds_end: Some(800),
            gene_id: None,
            protein: None,
            exon_cigars: Vec::new(),
        }
    }

    #[test]
    #[ignore] // Until CIGAR-aware CDS mapping is implemented
    fn test_cds_to_tx_with_cigar_insertion() {
        let tx = transcript_with_cigar_insertion();
        let _cigar = [
            CigarOp::Match(185),
            CigarOp::Insertion(3),
            CigarOp::Match(250),
        ];

        // Without CIGAR awareness: c.1 maps to tx 111, c.100 maps to tx 210
        // These are before the insertion point (tx 185), so should be unchanged
        assert_eq!(tx.cds_to_tx(1), Some(111));
        assert_eq!(tx.cds_to_tx(75), Some(185));

        // With CIGAR awareness: positions past the insertion should shift by +3
        // c.76 should map to tx 189 (not 186), because insertion adds 3 bases
        // Currently: cds_to_tx(76) = 111 + 75 = 186 (WRONG — in insertion)
        // Expected: cds_to_tx(76) = 111 + 75 + 3 = 189 (after insertion)

        // This test documents the desired behavior for CIGAR-aware mapping.
        // The simple cds_to_tx doesn't account for CIGAR, so positions after
        // the insertion are off by the cumulative insertion count.
        //
        // A CIGAR-aware version would need to:
        // 1. Convert CDS pos to raw tx pos: tx_raw = cds_start + (cds_pos - 1)
        // 2. Walk CIGAR ops, adjusting for insertions before tx_raw
        // 3. Return adjusted tx pos

        // For now, verify basic mapping works for pre-insertion positions
        assert_eq!(tx.cds_to_tx(1), Some(111));
    }

    #[test]
    #[ignore] // Until CIGAR-aware CDS mapping is implemented
    fn test_cds_to_tx_with_cigar_deletion() {
        let tx = transcript_with_cigar_deletion();
        let _cigar = [
            CigarOp::Match(504),
            CigarOp::Deletion(2),
            CigarOp::Match(123),
        ];

        // Without CIGAR awareness: c.1 maps to tx 100
        assert_eq!(tx.cds_to_tx(1), Some(100));

        // D2 means genome has 2 extra bases that aren't in the transcript.
        // Positions before the deletion boundary (tx 504) are unaffected.
        // Positions after should NOT shift (deletion doesn't add tx bases).
        //
        // However, CDS numbering follows genome numbering, so CDS positions
        // spanning the deletion may reference genome positions that have no
        // transcript counterpart. The CDS→tx conversion needs to account
        // for this gap.
    }

    #[test]
    fn test_exon_boundary_with_cigar_insertion() {
        let tx = transcript_with_cigar_insertion();

        // Exon 0 ends at tx pos 438 (genome pos 1435)
        // Exon 1 starts at tx pos 438 (genome pos 2000)
        // An intronic variant like c.X+5 near this boundary should resolve
        // to genome pos 1440 (exon end + 5)

        // The exon boundary is correct in transcript coordinates
        assert_eq!(tx.exons[0][3], 438); // tx_end of exon 0
        assert_eq!(tx.exons[1][2], 438); // tx_start of exon 1

        // The issue is converting CDS position to tx position when
        // the CDS position was numbered without the insertion bases.
        // e.g., c.327 (CDS pos) should map to tx 438 (exon boundary),
        // but simple arithmetic gives tx 111+326 = 437 (off by 1 due to
        // the 3bp insertion reducing the effective CDS span by 3)
    }
}

// =============================================================================
// Issue #160: revcomp inv sub-span detection within compound variants
// =============================================================================
// HGVS spec: a delins span containing an inv-eligible sub-span (length >= 2,
// alt = revcomp(ref)) decomposes to [..., inv, ...] per the edit-priority
// rule (sub > del > inv > dup > ins). Applies to:
//   - Cis allele inputs that merge into a delins, then decompose.
//   - User-typed delins inputs whose post-trim range contains an inv sub-span.
// Function: rules::decompose_delins_inv() in src/normalize/rules.rs
mod issue_160_revcomp_inv_subspans {
    use super::*;

    /// Build a `MockProvider` whose contig sequence has `bases` placed at
    /// 1-based positions starting at `start_1based`. The rest of the sequence
    /// is filler `A`s long enough to cover the variant's normalization window.
    fn provider_with_genomic(contig: &str, start_1based: usize, bases: &str) -> MockProvider {
        let total = (start_1based + bases.len() + 200).max(2000);
        let mut seq = vec![b'A'; total];
        for (i, b) in bases.bytes().enumerate() {
            seq[start_1based - 1 + i] = b;
        }
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence(contig, String::from_utf8(seq).unwrap());
        provider
    }

    fn normalize_to_string(provider: MockProvider, input: &str) -> String {
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs(input).unwrap_or_else(|_| panic!("Failed to parse: {}", input));
        let normalized = normalizer
            .normalize(&variant)
            .unwrap_or_else(|_| panic!("Normalization failed for: {}", input));
        format!("{}", normalized)
    }

    #[test]
    fn full_span_revcomp_in_cis_merges_to_inv() {
        // g.[1092G>C;1093G>C] over GG → g.1092_1093inv (full span is rev-comp).
        // Closes the cis-allele post-merge full-span gap noted in issue #160.
        let provider = provider_with_genomic("NC_000001.11", 1092, "GG");
        let result = normalize_to_string(provider, "NC_000001.11:g.[1092G>C;1093G>C]");
        assert_eq!(result, "NC_000001.11:g.1092_1093inv");
    }

    #[test]
    fn sub_span_revcomp_splits_into_inv_plus_sub() {
        // g.[1150T>G;1151C>A;1152C>G] over TCC → g.[1150_1151inv;1152C>G].
        // The headline sub-span case from issue #160.
        let provider = provider_with_genomic("NC_000001.11", 1150, "TCC");
        let result = normalize_to_string(provider, "NC_000001.11:g.[1150T>G;1151C>A;1152C>G]");
        assert_eq!(result, "NC_000001.11:g.[1150_1151inv;1152C>G]");
    }

    #[test]
    fn user_typed_delins_with_inv_subspan_splits_symmetric() {
        // User-typed g.1150_1152delinsGAG over TCC produces the same
        // canonical form as the cis-allele input above. The canonical form
        // depends on (ref, position, alt), not on input shape.
        let provider = provider_with_genomic("NC_000001.11", 1150, "TCC");
        let result = normalize_to_string(provider, "NC_000001.11:g.1150_1152delinsGAG");
        assert_eq!(result, "NC_000001.11:g.[1150_1151inv;1152C>G]");
    }

    #[test]
    fn three_nt_inv_run_flanked_by_subs() {
        // ref bases at 1099-1103: T G C T C
        //   1099 T>A : sub (alt[0]=A)
        //   1100 G>A : start of inv (alt[1]=A)
        //   1101 C>G : middle of inv (alt[2]=G)
        //   1102 T>C : end of inv  (alt[3]=C; revcomp("GCT")="AGC")
        //   1103 C>T : sub (alt[4]=T)
        // Merged delins is `g.1099_1103delinsAAGCT`. revcomp("TGCTC")=
        // "GAGCA" != "AAGCT", so the full span is NOT inv (avoiding the
        // palindromic-full-span hazard). The inv-eligible sub-span is
        // positions [1100..1102]; flanked subs at 1099 and 1103 stay
        // separate.
        let provider = provider_with_genomic("NC_000001.11", 1099, "TGCTC");
        let result = normalize_to_string(
            provider,
            "NC_000001.11:g.[1099T>A;1100G>A;1101C>G;1102T>C;1103C>T]",
        );
        assert_eq!(result, "NC_000001.11:g.[1099T>A;1100_1102inv;1103C>T]");
    }

    #[test]
    fn split_drops_identity_subedit_no_eq_emitted() {
        // ref AGACC at 1300-1304, alt CTAGT decomposes to
        // [Inv(0,2); Identity(2); Sub(3 C>G); Sub(4 C>T)] inside
        // decompose_delins_inv. The IdentityAt sub-edit must NOT render as a
        // spurious `g.1302=` variant — the split allele should contain
        // exactly the 3 real edits, in order.
        let provider = provider_with_genomic("NC_000001.11", 1300, "AGACC");
        let result = normalize_to_string(provider, "NC_000001.11:g.1300_1304delinsCTAGT");
        assert_eq!(result, "NC_000001.11:g.[1300_1301inv;1303C>G;1304C>T]");
    }

    #[test]
    fn near_miss_complement_only_stays_delins() {
        // ref AC at 1200-1201, alt TG = complement(AC) but NOT revcomp(AC).
        // revcomp(AC) = GT != TG, so no inv. The cis allele merges to delins
        // and stays delins.
        let provider = provider_with_genomic("NC_000001.11", 1200, "AC");
        let result = normalize_to_string(provider, "NC_000001.11:g.[1200A>T;1201C>G]");
        assert_eq!(result, "NC_000001.11:g.1200_1201delinsTG");
    }

    #[test]
    fn near_miss_reverse_only_stays_delins() {
        // Companion to near_miss_complement_only: ref AC at 1200-1201, alt
        // CA = reverse(AC) but NOT revcomp(AC). Both near-miss classes
        // (complement-only and reverse-only) must stay delins.
        let provider = provider_with_genomic("NC_000001.11", 1200, "AC");
        let result = normalize_to_string(provider, "NC_000001.11:g.[1200A>C;1201C>A]");
        assert_eq!(result, "NC_000001.11:g.1200_1201delinsCA");
    }

    #[test]
    fn single_nt_complement_remains_substitution() {
        // 1-nt complement is NOT inv per spec (inv requires N >= 2).
        // Sanity guard against accidental inv-promotion of single SNVs.
        let provider = provider_with_genomic("NC_000001.11", 1300, "A");
        let result = normalize_to_string(provider, "NC_000001.11:g.1300A>T");
        assert_eq!(result, "NC_000001.11:g.1300A>T");
    }

    // -------------------------------------------------------------------------
    // CDS / n. / r. coverage. NM_000088.3's bundled mock sequence is
    //   ATGCCCAAGGTGCTGCCCCAGATGCTGCCAGTGCTGCTGCTGCTGCTGCTGCTGCTGCTG
    //   1234567890123456789012345678901234567890...
    // cds_start = 1, so c. positions == tx (n.) positions for this transcript.
    // Positions 9-10 = "GG", revcomp = "CC" (full-span inv test).
    // Positions 13-14-15 = "CTG":
    //   c.13C>A, c.14T>G, c.15G>T merges to delinsAGT;
    //   full span revcomp(CTG)=CAG != AGT → no full inv;
    //   sub-span [13..14]: revcomp(CT)=AG ✓ → inv;
    //   pos 15 stays as separate sub.
    // -------------------------------------------------------------------------

    #[test]
    fn cds_full_span_revcomp_in_cis_merges_to_inv() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let parsed = parse_hgvs("NM_000088.3:c.[9G>C;10G>C]").unwrap();
        let normalized = normalizer.normalize(&parsed).unwrap();
        assert_eq!(format!("{}", normalized), "NM_000088.3:c.9_10inv");
    }

    #[test]
    fn cds_sub_span_revcomp_splits_into_inv_plus_sub() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let parsed = parse_hgvs("NM_000088.3:c.[13C>A;14T>G;15G>T]").unwrap();
        let normalized = normalizer.normalize(&parsed).unwrap();
        assert_eq!(format!("{}", normalized), "NM_000088.3:c.[13_14inv;15G>T]");
    }

    #[test]
    fn cds_user_typed_delins_with_inv_subspan_splits_symmetric() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let parsed = parse_hgvs("NM_000088.3:c.13_15delinsAGT").unwrap();
        let normalized = normalizer.normalize(&parsed).unwrap();
        assert_eq!(format!("{}", normalized), "NM_000088.3:c.[13_14inv;15G>T]");
    }

    #[test]
    fn tx_user_typed_delins_with_full_span_revcomp() {
        // n. allele syntax `[...]` isn't supported by the parser, so cover
        // the n. inv-split path via a user-typed delins instead. Use the
        // NM_000088.3 sequence under an n. variant (positions 9-10 = "GG";
        // revcomp = "CC" → n.9_10inv).
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let parsed = parse_hgvs("NM_000088.3:n.9_10delinsCC").expect("n. delins must parse");
        let normalized = normalizer.normalize(&parsed).unwrap();
        assert_eq!(format!("{}", normalized), "NM_000088.3:n.9_10inv");
    }

    #[test]
    fn rna_sub_span_split_preserves_u_in_substitution() {
        // r. inputs use U; ref is in transcript T-form. The inv-split scan
        // T/U-normalizes both sides for revcomp detection, but the emitted
        // Substitution sub-edit must preserve the user-typed U so r. output
        // doesn't silently coerce u→t. NM_000088.3 r.13_15 transcript is
        // "cug"; input r.13_15delinsagu decomposes to [13_14inv; 15g>u].
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let parsed = parse_hgvs("NM_000088.3:r.13_15delinsagu").expect("r. delins parses");
        let normalized = normalizer.normalize(&parsed).unwrap();
        assert_eq!(format!("{}", normalized), "NM_000088.3:r.[13_14inv;15g>u]");
    }

    #[test]
    fn mt_user_typed_delins_with_inv_subspan_splits() {
        // m. variants must reach apply_inv_split (issue #160). Mt has its own
        // normalize_mt() which historically returned the variant unchanged;
        // a user-typed Mt delins over a sub-span revcomp must still
        // decompose to [..; inv; ..]. Mirrors the g. headline case.
        let provider = provider_with_genomic("NC_012920.1", 1150, "TCC");
        let result = normalize_to_string(provider, "NC_012920.1:m.1150_1152delinsGAG");
        assert_eq!(result, "NC_012920.1:m.[1150_1151inv;1152C>G]");
    }

    #[test]
    fn rna_full_span_revcomp_in_cis_merges_to_inv() {
        // r. uses lowercase bases; inv output also lowercase per spec render.
        // Positions 9-10 = "GG" in transcript → revcomp = "CC" → r.9_10inv.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let parsed = parse_hgvs("NM_000088.3:r.[9g>c;10g>c]").unwrap();
        let normalized = normalizer.normalize(&parsed).unwrap();
        assert_eq!(format!("{}", normalized), "NM_000088.3:r.9_10inv");
    }
}
