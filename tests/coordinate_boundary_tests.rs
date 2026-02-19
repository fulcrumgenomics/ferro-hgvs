//! Comprehensive coordinate system boundary tests
//!
//! These tests are specifically designed to catch off-by-one errors
//! at coordinate system boundaries, covering scenarios from recent bugs:
//!
//! - cdot tx coordinates misinterpreted as 0-based (commit 944a4e9)
//! - Insertion shuffle boundary checking (commit d8391aa)
//! - 5' UTR coordinate conversion off-by-one
//! - CDS coordinate indexing bug

use ferro_hgvs::coords::{
    cdot_genomic_to_closed, cdot_tx_coords, hgvs_pos_to_index, hgvs_to_spdi_pos, index_to_hgvs_pos,
    spdi_to_hgvs_pos, OneBasedInterval, OneBasedPos, ZeroBasedInterval, ZeroBasedPos,
};

// ============================================================================
// Position conversion tests
// ============================================================================

mod position_conversions {
    use super::*;

    #[test]
    fn test_zero_to_one_based_first_element() {
        // The most critical case: first element
        // Index 0 (first element) -> position 1
        assert_eq!(ZeroBasedPos::new(0).to_one_based().value(), 1);
    }

    #[test]
    fn test_one_to_zero_based_first_element() {
        // Position 1 -> index 0
        assert_eq!(OneBasedPos::new(1).to_zero_based().value(), 0);
    }

    #[test]
    fn test_conversion_at_common_boundaries() {
        // Test at powers of 10 - common mistake points
        let test_cases = [(0u64, 1u64), (9, 10), (99, 100), (999, 1000), (9999, 10000)];

        for (zb_val, ob_val) in test_cases {
            let zb = ZeroBasedPos::new(zb_val);
            let ob = zb.to_one_based();
            assert_eq!(
                ob.value(),
                ob_val,
                "ZeroBasedPos({}) should convert to OneBasedPos({})",
                zb_val,
                ob_val
            );

            let ob = OneBasedPos::new(ob_val);
            let zb = ob.to_zero_based();
            assert_eq!(
                zb.value(),
                zb_val,
                "OneBasedPos({}) should convert to ZeroBasedPos({})",
                ob_val,
                zb_val
            );
        }
    }

    #[test]
    fn test_roundtrip_preserves_value() {
        for i in 0..1000 {
            let zb = ZeroBasedPos::new(i);
            let ob = zb.to_one_based();
            let zb_again = ob.to_zero_based();
            assert_eq!(
                zb.value(),
                zb_again.value(),
                "Roundtrip failed for value {}",
                i
            );
        }

        for i in 1..=1000 {
            let ob = OneBasedPos::new(i);
            let zb = ob.to_zero_based();
            let ob_again = zb.to_one_based();
            assert_eq!(
                ob.value(),
                ob_again.value(),
                "Roundtrip failed for value {}",
                i
            );
        }
    }
}

// ============================================================================
// One-based position validation tests
// ============================================================================

mod one_based_validation {
    use super::*;

    #[test]
    #[should_panic(expected = "1-based position cannot be 0")]
    fn test_one_based_rejects_zero() {
        let _ = OneBasedPos::new(0);
    }

    #[test]
    fn test_try_new_returns_none_for_zero() {
        assert!(OneBasedPos::try_new(0).is_none());
    }

    #[test]
    fn test_try_new_returns_some_for_valid() {
        assert!(OneBasedPos::try_new(1).is_some());
        assert!(OneBasedPos::try_new(100).is_some());
    }

    #[test]
    fn test_new_unchecked_allows_zero() {
        // This should not panic - it's unchecked
        let pos = OneBasedPos::new_unchecked(0);
        assert_eq!(pos.value(), 0);
        // Note: using this position for conversion would be incorrect
    }
}

// ============================================================================
// Interval conversion tests
// ============================================================================

mod interval_conversions {
    use super::*;

    #[test]
    fn test_half_open_to_closed_basic() {
        // [0, 10) in 0-based half-open = [1, 10] in 1-based closed
        let ho = ZeroBasedInterval {
            start: ZeroBasedPos::new(0),
            end: ZeroBasedPos::new(10),
        };
        let closed = ho.to_one_based_closed();

        assert_eq!(closed.start.value(), 1);
        assert_eq!(closed.end.value(), 10);
    }

    #[test]
    fn test_closed_to_half_open_basic() {
        // [1, 10] in 1-based closed = [0, 10) in 0-based half-open
        let closed = OneBasedInterval {
            start: OneBasedPos::new(1),
            end: OneBasedPos::new(10),
        };
        let ho = closed.to_zero_based_half_open();

        assert_eq!(ho.start.value(), 0);
        assert_eq!(ho.end.value(), 10);
    }

    #[test]
    fn test_interval_length_consistency() {
        // Both should represent the same length
        let ho = ZeroBasedInterval {
            start: ZeroBasedPos::new(5),
            end: ZeroBasedPos::new(15),
        };
        let closed = ho.to_one_based_closed();

        assert_eq!(ho.len(), 10);
        assert_eq!(closed.len(), 10);
        assert_eq!(ho.len(), closed.len());
    }

    #[test]
    fn test_single_element_interval() {
        // [5, 6) in half-open = [6, 6] in closed (single element)
        let ho = ZeroBasedInterval {
            start: ZeroBasedPos::new(5),
            end: ZeroBasedPos::new(6),
        };
        assert_eq!(ho.len(), 1);

        let closed = ho.to_one_based_closed();
        assert_eq!(closed.start.value(), 6);
        assert_eq!(closed.end.value(), 6);
        assert_eq!(closed.len(), 1);
    }
}

// ============================================================================
// cdot coordinate tests (covers bug 944a4e9)
// ============================================================================

mod cdot_coordinates {
    use super::*;

    #[test]
    fn test_cdot_genomic_to_closed() {
        // cdot genomic: [0, 100) 0-based half-open
        // Expected: [1, 100] 1-based closed
        let (start, end) = cdot_genomic_to_closed(0, 100);
        assert_eq!(start, 1);
        assert_eq!(end, 100);
    }

    #[test]
    fn test_cdot_tx_coords_unchanged() {
        // cdot tx coordinates are already 1-based
        // This is a documentation function - values should be unchanged
        // This tests the fix for bug 944a4e9

        let (tx_start, tx_end) = cdot_tx_coords(1, 100);
        assert_eq!(tx_start, 1, "cdot tx_start should remain 1-based");
        assert_eq!(tx_end, 100, "cdot tx_end should remain 1-based");

        // First exon should start at 1, not 0
        let (tx_start, _) = cdot_tx_coords(1, 50);
        assert_eq!(
            tx_start, 1,
            "First exon tx_start should be 1, not 0 (bug 944a4e9)"
        );
    }

    #[test]
    fn test_cdot_mixed_coordinate_awareness() {
        // cdot uses DIFFERENT coordinate systems for genomic vs transcript:
        // - Genomic: 0-based half-open
        // - Transcript: 1-based

        // Simulate cdot exon: [1000, 1100, 1, 1, 100]
        // genomic_start=1000, genomic_end=1100 (0-based half-open)
        // exon_num=1, tx_start=1, tx_end=100 (1-based)

        let g_start = 1000u64;
        let g_end = 1100u64;
        let tx_start = 1u64;
        let tx_end = 100u64;

        // Convert genomic to 1-based closed
        let (g_start_1b, g_end_1b) = cdot_genomic_to_closed(g_start, g_end);
        assert_eq!(g_start_1b, 1001);
        assert_eq!(g_end_1b, 1100);

        // Transcript coords should NOT be converted
        let (tx_start_out, tx_end_out) = cdot_tx_coords(tx_start, tx_end);
        assert_eq!(tx_start_out, 1); // NOT 2
        assert_eq!(tx_end_out, 100); // NOT 101

        // Lengths should match
        let genomic_len = g_end - g_start; // 100 (half-open)
        let tx_len = tx_end - tx_start + 1; // 100 (closed)
        assert_eq!(genomic_len, tx_len);
    }
}

// ============================================================================
// SPDI/HGVS conversion tests
// ============================================================================

mod spdi_hgvs_conversions {
    use super::*;

    #[test]
    fn test_spdi_to_hgvs_first_position() {
        // SPDI position 0 = HGVS position 1
        assert_eq!(spdi_to_hgvs_pos(0), 1);
    }

    #[test]
    fn test_hgvs_to_spdi_first_position() {
        // HGVS position 1 = SPDI position 0
        assert_eq!(hgvs_to_spdi_pos(1), 0);
    }

    #[test]
    fn test_spdi_hgvs_roundtrip() {
        for i in 0..1000 {
            let hgvs = spdi_to_hgvs_pos(i);
            let back = hgvs_to_spdi_pos(hgvs);
            assert_eq!(i, back, "SPDI->HGVS->SPDI roundtrip failed for {}", i);
        }

        for i in 1..=1000 {
            let spdi = hgvs_to_spdi_pos(i);
            let back = spdi_to_hgvs_pos(spdi);
            assert_eq!(i, back, "HGVS->SPDI->HGVS roundtrip failed for {}", i);
        }
    }

    #[test]
    fn test_spdi_hgvs_example_from_docs() {
        // From SPDI documentation:
        // For a variant at genomic position 12345 (1-based HGVS),
        // SPDI position is 12344 (0-based)
        assert_eq!(hgvs_to_spdi_pos(12345), 12344);
        assert_eq!(spdi_to_hgvs_pos(12344), 12345);
    }
}

// ============================================================================
// Array indexing tests
// ============================================================================

mod array_indexing {
    use super::*;

    #[test]
    fn test_hgvs_pos_to_index_first_element() {
        // HGVS position 1 (first base) -> index 0
        assert_eq!(hgvs_pos_to_index(1), 0);
    }

    #[test]
    fn test_index_to_hgvs_pos_first_element() {
        // Index 0 -> HGVS position 1
        assert_eq!(index_to_hgvs_pos(0), 1);
    }

    #[test]
    fn test_array_access_with_zero_based() {
        let seq = b"ACGTACGT";

        // Using ZeroBasedPos for array access
        let pos = ZeroBasedPos::new(0);
        assert_eq!(seq[pos.as_index()], b'A');

        let pos = ZeroBasedPos::new(3);
        assert_eq!(seq[pos.as_index()], b'T');

        let pos = ZeroBasedPos::new(7);
        assert_eq!(seq[pos.as_index()], b'T');
    }

    #[test]
    fn test_array_access_with_one_based_conversion() {
        let seq = b"ACGTACGT";

        // HGVS position 1 (first base) should access index 0
        let hgvs_pos = 1u64;
        let idx = hgvs_pos_to_index(hgvs_pos);
        assert_eq!(seq[idx], b'A');

        // HGVS position 4 should access index 3 (T)
        let hgvs_pos = 4u64;
        let idx = hgvs_pos_to_index(hgvs_pos);
        assert_eq!(seq[idx], b'T');
    }
}

// ============================================================================
// Insertion position tests (covers bug d8391aa)
// ============================================================================

mod insertion_positions {
    use super::*;

    #[test]
    fn test_insertion_flanking_positions() {
        // In HGVS, an insertion like g.100_101insATG means:
        // - Insert between positions 100 and 101
        // - The flanking bases are at positions 100 and 101

        let hgvs_pos_left = 100u64;
        let hgvs_pos_right = 101u64;

        // To check sequence at these positions, convert to 0-based indices
        let idx_left = hgvs_pos_to_index(hgvs_pos_left);
        let idx_right = hgvs_pos_to_index(hgvs_pos_right);

        // Index 99 is position 100, index 100 is position 101
        assert_eq!(idx_left, 99);
        assert_eq!(idx_right, 100);

        // These are adjacent
        assert_eq!(idx_right - idx_left, 1);
    }

    #[test]
    fn test_insertion_shuffle_boundary() {
        // Bug d8391aa: insertion normalization was checking ref_seq[end]
        // instead of ref_seq[end-1] for the second flanking position

        // For g.100_101insA, when checking if insert matches flanking:
        // - Left flank is ref_seq[99] (HGVS 100)
        // - Right flank is ref_seq[100] (HGVS 101)

        let seq = b"ACGTACGTACGT";
        let hgvs_start = 5u64; // g.5
        let hgvs_end = 6u64; // g.6

        // Correct: check ref_seq[4] and ref_seq[5]
        let idx_start = hgvs_pos_to_index(hgvs_start);
        let idx_end = hgvs_pos_to_index(hgvs_end);

        assert_eq!(idx_start, 4);
        assert_eq!(idx_end, 5);
        assert_eq!(seq[idx_start], b'A');
        assert_eq!(seq[idx_end], b'C');
    }
}

// ============================================================================
// Arithmetic safety tests
// ============================================================================

mod arithmetic_safety {
    use super::*;

    #[test]
    fn test_zero_based_checked_sub_prevents_underflow() {
        let pos = ZeroBasedPos::new(5);

        // Valid subtraction
        assert_eq!(pos.checked_sub(3), Some(ZeroBasedPos::new(2)));

        // Would underflow
        assert_eq!(pos.checked_sub(10), None);

        // Exact boundary
        assert_eq!(pos.checked_sub(5), Some(ZeroBasedPos::new(0)));
    }

    #[test]
    fn test_one_based_checked_sub_prevents_zero() {
        let pos = OneBasedPos::new(5);

        // Valid subtraction
        assert_eq!(pos.checked_sub(3), Some(OneBasedPos::new(2)));

        // Would become 0 (invalid for 1-based)
        assert_eq!(pos.checked_sub(5), None);

        // Would underflow
        assert_eq!(pos.checked_sub(10), None);
    }

    #[test]
    fn test_saturating_sub() {
        let pos = ZeroBasedPos::new(5);

        assert_eq!(pos.saturating_sub(3), ZeroBasedPos::new(2));
        assert_eq!(pos.saturating_sub(10), ZeroBasedPos::new(0)); // Saturates at 0
    }
}

// ============================================================================
// Edge case tests
// ============================================================================

mod edge_cases {
    use super::*;

    #[test]
    fn test_large_positions() {
        // Test with large values that might cause overflow issues
        let large = u64::MAX - 1;

        let zb = ZeroBasedPos::new(large);
        // This should overflow
        assert!(zb.checked_add(10).is_none());

        // But this should work
        let zb = ZeroBasedPos::new(1000);
        assert!(zb.checked_add(1000).is_some());
    }

    #[test]
    fn test_display_formatting() {
        let zb = ZeroBasedPos::new(42);
        assert_eq!(format!("{}", zb), "42");

        let ob = OneBasedPos::new(42);
        assert_eq!(format!("{}", ob), "42");
    }

    #[test]
    fn test_ordering() {
        let a = ZeroBasedPos::new(5);
        let b = ZeroBasedPos::new(10);
        assert!(a < b);
        assert!(b > a);
        assert_eq!(a, ZeroBasedPos::new(5));

        let a = OneBasedPos::new(5);
        let b = OneBasedPos::new(10);
        assert!(a < b);
        assert!(b > a);
        assert_eq!(a, OneBasedPos::new(5));
    }

    #[test]
    fn test_from_conversions() {
        let zb = ZeroBasedPos::new(42);
        let raw: u64 = zb.into();
        assert_eq!(raw, 42);

        let ob = OneBasedPos::new(42);
        let raw: u64 = ob.into();
        assert_eq!(raw, 42);
    }
}
