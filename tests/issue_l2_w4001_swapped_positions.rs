//! Audit + regression tests for #81 § L2 — wire W4001
//! (SwappedPositions) into the preprocessor.
//!
//! The corrector `correct_swapped_positions` existed in
//! `src/error_handling/corrections.rs` from #137's L1 audit but the
//! preprocessor never called it; the audit fixture classified W4001 as
//! `registered`. This PR adds Phase 15a to the preprocessor and pins
//! the resulting emission surface.

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::{parse_hgvs_lenient, parse_hgvs_with_config};
use ferro_hgvs::parse_hgvs;

fn warning_codes(input: &str) -> Vec<String> {
    let r = parse_hgvs_lenient(input).unwrap_or_else(|e| panic!("lenient parse {input:?}: {e}"));
    r.warnings
        .iter()
        .map(|w| w.error_type.code().to_string())
        .collect()
}

// =============================================================================
// SECTION 1 — Lenient mode emits W4001 and rewrites
// =============================================================================

mod lenient_emits_w4001 {
    use super::*;

    #[test]
    fn cds_swapped_emits_w4001_and_rewrites() {
        let r = parse_hgvs_lenient("NM_000088.3:c.200_100del").unwrap();
        assert_eq!(format!("{}", r.result), "NM_000088.3:c.100_200del");
        assert!(r.warnings.iter().any(|w| w.error_type.code() == "W4001"));
    }

    #[test]
    fn genome_swapped_dup_emits_w4001_and_rewrites() {
        let r = parse_hgvs_lenient("NC_000001.11:g.500_100dup").unwrap();
        assert_eq!(format!("{}", r.result), "NC_000001.11:g.100_500dup");
        assert!(r.warnings.iter().any(|w| w.error_type.code() == "W4001"));
    }

    #[test]
    fn ascending_input_emits_no_w4001() {
        let codes = warning_codes("NM_000088.3:c.100_200del");
        assert!(
            !codes.iter().any(|c| c == "W4001"),
            "ascending range must not emit W4001, got {codes:?}"
        );
    }
}

// =============================================================================
// SECTION 2 — Strict-mode rejection through the preprocessor
// =============================================================================
//
// `parse_hgvs_with_config(ErrorConfig::strict())` runs the preprocessor
// (which now includes Phase 15a) and rejects swapped ranges with a
// W4001 diagnostic. The bare `parse_hgvs` entry point bypasses the
// preprocessor and goes straight to the nom parser; the parser
// accepts swapped CDS ranges as a Verify-failure and accepts swapped
// genomic ranges silently. Pin both behaviors so a future "always run
// the preprocessor" change is surfaced.

mod strict_through_preprocessor {
    use super::*;

    #[test]
    fn strict_config_cds_swapped_rejected() {
        let result = parse_hgvs_with_config("NM_000088.3:c.200_100del", ErrorConfig::strict());
        assert!(result.is_err());
    }

    #[test]
    fn strict_config_genome_swapped_rejected() {
        let result = parse_hgvs_with_config("NC_000001.11:g.500_100dup", ErrorConfig::strict());
        assert!(result.is_err());
    }

    /// Bare `parse_hgvs` bypasses the preprocessor; the CDS parser
    /// rejects swapped ranges via the underlying nom Verify check.
    #[test]
    fn bare_parse_cds_swapped_rejected() {
        assert!(parse_hgvs("NM_000088.3:c.200_100del").is_err());
    }

    /// Bare `parse_hgvs` bypasses the preprocessor; the genomic parser
    /// accepts swapped ranges silently. Pinned to surface this gap if
    /// the parser is ever tightened.
    #[test]
    fn bare_parse_genome_swapped_accepted_by_parser() {
        assert!(parse_hgvs("NC_000001.11:g.500_100dup").is_ok());
    }
}

// =============================================================================
// SECTION 3 — Partial coverage: offset-bearing forms not handled
// =============================================================================
//
// Pin the documented partial-coverage limit: the corrector handles only
// simple integer pairs (`c.A_Bdel`). Offset-bearing inverted forms fall
// through to the parser. These tests pin the gap so a future enhancement
// surfaces them.

mod partial_offset_coverage {
    use super::*;

    /// `c.100+5_99+3` has inverted positions across an intron offset;
    /// the corrector explicitly bails on offset-bearing forms so the
    /// input is not silently rewritten (which would drop the offset).
    #[test]
    fn offset_swapped_not_rewritten() {
        let r = parse_hgvs_lenient("NM_000088.3:c.100+5_99+3del").unwrap();
        // No W4001 emitted; offset-bearing form is left unchanged.
        assert!(
            !r.warnings.iter().any(|w| w.error_type.code() == "W4001"),
            "offset-bearing swapped form must not emit W4001 (corrector would drop the offset)"
        );
        // Round-trip is unchanged.
        assert_eq!(format!("{}", r.result), "NM_000088.3:c.100+5_99+3del");
    }

    /// `c.*5_*1del` (3'UTR `*N`) — also bails because `*` is not a bare
    /// integer; the corrector is integer-only.
    #[test]
    fn three_prime_utr_swapped_not_rewritten() {
        let r = parse_hgvs_lenient("NM_000088.3:c.*5_*1del");
        if let Ok(r) = r {
            assert!(
                !r.warnings.iter().any(|w| w.error_type.code() == "W4001"),
                "three-prime-UTR swapped form unexpectedly emitted W4001"
            );
        }
    }
}
