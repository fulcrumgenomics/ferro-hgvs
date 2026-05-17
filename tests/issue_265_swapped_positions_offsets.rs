//! Audit + regression tests for #265 — extend W4001 (SwappedPositions)
//! to offset-bearing positions (`c.100+5_99+3del`) and 3'UTR `*N`
//! markers (`c.*5_*1del`).
//!
//! #264 wired the integer-only corrector into the preprocessor; it
//! explicitly bailed on offset / `*N` forms to avoid silently dropping
//! offsets. #265 replaces the integer-only parser with a full endpoint
//! tokenizer that captures the raw substring and a comparable sort key.

use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;

#[track_caller]
fn assert_rewrites_with_w4001(input: &str, expected: &str) {
    let r = parse_hgvs_lenient(input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
    assert_eq!(
        format!("{}", r.result),
        expected,
        "rewrite mismatch for {input:?}"
    );
    assert!(
        r.warnings.iter().any(|w| w.error_type.code() == "W4001"),
        "expected W4001 for {input:?}; got {:?}",
        r.warnings
            .iter()
            .map(|w| w.error_type.code())
            .collect::<Vec<_>>()
    );
}

#[track_caller]
fn assert_no_rewrite_no_w4001(input: &str) {
    let r = parse_hgvs_lenient(input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
    assert_eq!(
        format!("{}", r.result),
        input,
        "input must round-trip unchanged for {input:?}"
    );
    assert!(
        !r.warnings.iter().any(|w| w.error_type.code() == "W4001"),
        "unexpected W4001 for {input:?}"
    );
}

// =============================================================================
// SECTION 1 — Intronic-offset swaps
// =============================================================================

mod intronic_offset {
    use super::*;

    /// Different main coords with `+offset`. Sort key is
    /// `(main_axis, offset)` so 100+5 > 99+3 swaps to 99+3 < 100+5.
    #[test]
    fn cross_main_plus_offset_swapped() {
        assert_rewrites_with_w4001("NM_000088.3:c.100+5_99+3del", "NM_000088.3:c.99+3_100+5del");
    }

    /// Same main coord, swapped `+offset` magnitude.
    #[test]
    fn same_main_swapped_plus_offsets() {
        assert_rewrites_with_w4001(
            "NM_000088.3:c.100+5_100+3del",
            "NM_000088.3:c.100+3_100+5del",
        );
    }

    /// `-offset` (acceptor side) — same main coord, swapped magnitudes.
    /// `100-2` < `100-5` because `-2 > -5` on the signed offset axis.
    #[test]
    fn same_main_swapped_minus_offsets() {
        assert_rewrites_with_w4001(
            "NM_000088.3:c.100-2_100-5del",
            "NM_000088.3:c.100-5_100-2del",
        );
    }

    /// Ascending intronic offsets stay unchanged.
    #[test]
    fn ascending_intronic_unchanged() {
        assert_no_rewrite_no_w4001("NM_000088.3:c.99+3_100+5del");
        assert_no_rewrite_no_w4001("NM_000088.3:c.100-5_100-2del");
    }

    /// Main coord without offset vs main coord with positive offset:
    /// `100` < `100+3` (zero offset less than positive offset).
    #[test]
    fn bare_main_versus_plus_offset_ascending_unchanged() {
        assert_no_rewrite_no_w4001("NM_000088.3:c.100_100+3del");
    }

    /// `100+3` > `100` so it would swap; pinned.
    #[test]
    fn plus_offset_versus_bare_main_swapped() {
        assert_rewrites_with_w4001("NM_000088.3:c.100+3_100del", "NM_000088.3:c.100_100+3del");
    }
}

// =============================================================================
// SECTION 2 — 3'UTR `*N` marker swaps
// =============================================================================

mod three_prime_utr_star {
    use super::*;

    /// Two `*N` endpoints, swapped.
    #[test]
    fn star_to_star_swapped() {
        assert_rewrites_with_w4001("NM_000088.3:c.*5_*1del", "NM_000088.3:c.*1_*5del");
    }

    /// Ascending `*N` pair unchanged.
    #[test]
    fn ascending_star_pair_unchanged() {
        assert_no_rewrite_no_w4001("NM_000088.3:c.*1_*5del");
    }

    /// `*N` always sorts after non-`*` per the spec — CDS-start to
    /// 3'UTR-end ascending stays unchanged.
    #[test]
    fn cds_to_star_ascending_unchanged() {
        assert_no_rewrite_no_w4001("NM_000088.3:c.100_*5del");
    }

    /// `*N` start with non-`*` end is a class swap.
    #[test]
    fn star_to_cds_swapped() {
        assert_rewrites_with_w4001("NM_000088.3:c.*1_100del", "NM_000088.3:c.100_*1del");
    }
}

// =============================================================================
// SECTION 3 — 5'UTR `-N` marker swaps
// =============================================================================

mod five_prime_utr_minus {
    use super::*;

    /// CDS start, 5'UTR end (`5_-3`) is a class swap; `-3 < 5` so swap.
    #[test]
    fn cds_to_five_prime_utr_swapped() {
        assert_rewrites_with_w4001("NM_000088.3:c.5_-3del", "NM_000088.3:c.-3_5del");
    }

    /// 5'UTR start, CDS end ascending unchanged.
    #[test]
    fn five_prime_utr_to_cds_ascending_unchanged() {
        assert_no_rewrite_no_w4001("NM_000088.3:c.-3_5del");
    }

    /// Two negative endpoints, swapped magnitudes. `-3 > -5` so
    /// `-3_-5` is descending → rewrites to `-5_-3`.
    #[test]
    fn five_prime_utr_two_negatives_swapped() {
        assert_rewrites_with_w4001("NM_000088.3:c.-3_-5del", "NM_000088.3:c.-5_-3del");
    }

    /// Two negative endpoints, ascending (`-5 < -3`) unchanged.
    #[test]
    fn five_prime_utr_two_negatives_ascending_unchanged() {
        assert_no_rewrite_no_w4001("NM_000088.3:c.-5_-3del");
    }
}

// =============================================================================
// SECTION 4 — Bare integer (regression guards from #264)
// =============================================================================

mod bare_integer_regression {
    use super::*;

    #[test]
    fn cds_bare_swapped() {
        assert_rewrites_with_w4001("NM_000088.3:c.200_100del", "NM_000088.3:c.100_200del");
    }

    #[test]
    fn genome_bare_swapped() {
        assert_rewrites_with_w4001("NC_000001.11:g.500_100dup", "NC_000001.11:g.100_500dup");
    }
}
