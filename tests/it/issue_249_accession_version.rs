//! Audit for #81 § I2 — versioned vs unversioned accession policy.
//!
//! Policy (already implemented):
//! 1. Strict parse accepts both versioned (`NM_000088.3`) and
//!    unversioned (`NM_000088`).
//! 2. Display preserves verbatim — never synthesizes or strips.
//! 3. Trailing empty-version dot (`NM_000088.`) is rejected.
//! 4. `PartialEq` treats versioned and unversioned as distinct.
//! 5. Lenient mode emits soft-validation warning W3001
//!    (`MissingVersion`) per occurrence; never mutates the parsed
//!    accession.

use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
use ferro_hgvs::parse_hgvs;

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
}

// =============================================================================
// SECTION 1 — Versioned accessions round-trip verbatim
// =============================================================================

mod versioned {
    use super::*;

    #[test]
    fn refseq_cds_versioned_round_trips() {
        assert_round_trips("NM_000088.3:c.100A>G");
    }

    #[test]
    fn refseq_genome_versioned_round_trips() {
        assert_round_trips("NC_000017.11:g.100A>G");
    }

    #[test]
    fn refseq_protein_versioned_round_trips() {
        assert_round_trips("NP_003997.1:p.Arg8Gln");
    }

    #[test]
    fn refseq_noncoding_versioned_round_trips() {
        assert_round_trips("NR_037639.1:n.100A>G");
    }

    #[test]
    fn ensembl_versioned_round_trips() {
        assert_round_trips("ENST00000456328.2:c.100A>G");
    }

    /// Multi-digit version numbers round-trip too.
    #[test]
    fn refseq_large_version_round_trips() {
        assert_round_trips("NM_000088.30:c.100A>G");
    }
}

// =============================================================================
// SECTION 2 — Unversioned accessions round-trip verbatim (no synthesis)
// =============================================================================

mod unversioned_preserve {
    use super::*;

    #[test]
    fn refseq_cds_unversioned_round_trips() {
        assert_round_trips("NM_000088:c.100A>G");
    }

    #[test]
    fn refseq_genome_unversioned_round_trips() {
        assert_round_trips("NC_000017:g.100A>G");
    }

    #[test]
    fn refseq_protein_unversioned_round_trips() {
        assert_round_trips("NP_003997:p.Arg8Gln");
    }

    #[test]
    fn refseq_noncoding_unversioned_round_trips() {
        assert_round_trips("NR_037639:n.100A>G");
    }

    #[test]
    fn ensembl_unversioned_round_trips() {
        assert_round_trips("ENST00000456328:c.100A>G");
    }
}

// =============================================================================
// SECTION 3 — Pathological forms
// =============================================================================

mod pathological {
    use super::*;

    /// Trailing empty version dot must be rejected.
    #[test]
    fn trailing_empty_dot_rejected() {
        assert!(parse_hgvs("NM_000088.:c.100A>G").is_err());
    }

    /// LRG accessions carry no version segment — they're inherently
    /// versioned by definition. Pin that the parser accepts the bare
    /// form (not "missing version").
    #[test]
    fn lrg_has_no_version_segment() {
        assert_round_trips("LRG_199t1:c.100A>G");
        assert_round_trips("LRG_199:g.100A>G");
        assert_round_trips("LRG_199p1:p.Arg8Gln");
    }
}

// =============================================================================
// SECTION 4 — Equality semantics
// =============================================================================

mod equality {
    use super::*;

    /// Versioned and unversioned are distinct identifiers — `PartialEq`
    /// must not collapse them.
    #[test]
    fn versioned_and_unversioned_are_distinct() {
        let versioned = parse_hgvs("NM_000088.3:c.100A>G").expect("versioned parse");
        let unversioned = parse_hgvs("NM_000088:c.100A>G").expect("unversioned parse");
        assert_ne!(versioned, unversioned);
    }

    /// Different version numbers are likewise distinct.
    #[test]
    fn different_versions_are_distinct() {
        let v3 = parse_hgvs("NM_000088.3:c.100A>G").expect("v3 parse");
        let v4 = parse_hgvs("NM_000088.4:c.100A>G").expect("v4 parse");
        assert_ne!(v3, v4);
    }
}

// =============================================================================
// SECTION 5 — Lenient mode emits W3001 for unversioned accessions
// =============================================================================

mod lenient_warnings {
    use super::*;

    fn warning_codes(input: &str) -> Vec<String> {
        let r = parse_hgvs_lenient(input)
            .unwrap_or_else(|e| panic!("lenient parse {input:?} failed: {e}"));
        r.warnings
            .iter()
            .map(|w| w.error_type.code().to_string())
            .collect()
    }

    #[test]
    fn versioned_emits_no_w3001() {
        let codes = warning_codes("NM_000088.3:c.100A>G");
        assert!(
            !codes.iter().any(|c| c == "W3001"),
            "versioned should not emit W3001, got {codes:?}"
        );
    }

    #[test]
    fn unversioned_emits_w3001() {
        let codes = warning_codes("NM_000088:c.100A>G");
        assert!(
            codes.iter().any(|c| c == "W3001"),
            "unversioned should emit W3001, got {codes:?}"
        );
    }

    /// Lenient mode does not mutate the parsed accession — the warning
    /// is emitted, but the variant's accession remains unversioned.
    #[test]
    fn lenient_does_not_synthesize_version() {
        let r = parse_hgvs_lenient("NM_000088:c.100A>G").unwrap();
        assert_eq!(format!("{}", r.result), "NM_000088:c.100A>G");
    }

    /// LRG accessions have no missing-version concept — no W3001.
    #[test]
    fn lrg_emits_no_w3001() {
        let codes = warning_codes("LRG_199t1:c.100A>G");
        assert!(
            !codes.iter().any(|c| c == "W3001"),
            "LRG should not emit W3001, got {codes:?}"
        );
    }
}
