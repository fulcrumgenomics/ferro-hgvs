//! Audit for #81 § I4 — chromosome alias parse + canonical output.
//!
//! Policy (already implemented):
//! 1. Strict parse accepts the full alias surface: bare `chr1`, bare
//!    numeric `1`, case variants, special chromosomes (`X`, `Y`, `M`,
//!    `MT`), assembly-wrapped (`GRCh38(chr1)`), canonical `NC_…`.
//! 2. Display preserves verbatim — no canonicalization to `NC_` from a
//!    bare alias (would require assembly context the parser lacks).
//! 3. PartialEq treats different alias forms as distinct accessions,
//!    even when they refer to the same physical chromosome.

use ferro_hgvs::parse_hgvs;

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
}

// =============================================================================
// SECTION 1 — Bare `chr<N>` aliases
// =============================================================================

mod chr_prefixed {
    use super::*;

    #[test]
    fn chr1_round_trips() {
        assert_round_trips("chr1:g.100A>G");
    }

    #[test]
    fn chr17_round_trips() {
        assert_round_trips("chr17:g.100A>G");
    }

    #[test]
    fn chr_x_round_trips() {
        assert_round_trips("chrX:g.100A>G");
    }

    #[test]
    fn chr_y_round_trips() {
        assert_round_trips("chrY:g.100A>G");
    }

    #[test]
    fn chr_m_round_trips() {
        assert_round_trips("chrM:g.100A>G");
    }

    #[test]
    fn chr_mt_round_trips() {
        assert_round_trips("chrMT:g.100A>G");
    }
}

// =============================================================================
// SECTION 2 — Bare numeric/single-letter aliases
// =============================================================================

mod bare_aliases {
    use super::*;

    #[test]
    fn bare_1_round_trips() {
        assert_round_trips("1:g.100A>G");
    }

    #[test]
    fn bare_17_round_trips() {
        assert_round_trips("17:g.100A>G");
    }

    #[test]
    fn bare_x_round_trips() {
        assert_round_trips("X:g.100A>G");
    }
}

// =============================================================================
// SECTION 3 — Case variants
// =============================================================================

mod case_variants {
    use super::*;

    /// `CHR1` (all-uppercase) is accepted verbatim.
    #[test]
    fn upper_chr1_round_trips() {
        assert_round_trips("CHR1:g.100A>G");
    }

    /// `Chr1` (title-case) is accepted verbatim.
    #[test]
    fn title_chr1_round_trips() {
        assert_round_trips("Chr1:g.100A>G");
    }
}

// =============================================================================
// SECTION 4 — Assembly-wrapped aliases
// =============================================================================

mod assembly_wrapped {
    use super::*;

    #[test]
    fn grch37_chr1_round_trips() {
        assert_round_trips("GRCh37(chr1):g.100A>G");
    }

    #[test]
    fn grch38_chr1_round_trips() {
        assert_round_trips("GRCh38(chr1):g.100A>G");
    }

    #[test]
    fn grch38_chrx_round_trips() {
        assert_round_trips("GRCh38(chrX):g.100A>G");
    }
}

// =============================================================================
// SECTION 5 — Canonical `NC_…` accessions
// =============================================================================

mod canonical_refseq {
    use super::*;

    #[test]
    fn nc_chr1_grch38_round_trips() {
        assert_round_trips("NC_000001.11:g.100A>G");
    }

    #[test]
    fn nc_chr17_grch38_round_trips() {
        assert_round_trips("NC_000017.11:g.100A>G");
    }

    #[test]
    fn nc_chrx_grch38_round_trips() {
        assert_round_trips("NC_000023.11:g.100A>G");
    }

    #[test]
    fn nc_mitochondrion_round_trips() {
        assert_round_trips("NC_012920.1:m.100A>G");
    }
}

// =============================================================================
// SECTION 6 — Equality semantics
// =============================================================================
//
// Different alias forms referring to the same physical chromosome are
// distinct accessions — without an assembly context the parser cannot
// canonicalize them, and PartialEq operates on the parsed accession
// string, not on biological identity.

mod equality {
    use super::*;

    #[test]
    fn chr1_distinct_from_bare_1() {
        let a = parse_hgvs("chr1:g.100A>G").unwrap();
        let b = parse_hgvs("1:g.100A>G").unwrap();
        assert_ne!(a, b);
    }

    #[test]
    fn chr1_distinct_from_nc_000001() {
        let a = parse_hgvs("chr1:g.100A>G").unwrap();
        let b = parse_hgvs("NC_000001.11:g.100A>G").unwrap();
        assert_ne!(a, b);
    }

    #[test]
    fn case_variants_distinct() {
        let a = parse_hgvs("chr1:g.100A>G").unwrap();
        let b = parse_hgvs("CHR1:g.100A>G").unwrap();
        assert_ne!(a, b);
    }
}

// =============================================================================
// SECTION 7 — Display does not synthesize a canonical accession
// =============================================================================
//
// Bare aliases stay bare on Display. We pin this because the alternative
// would require assembly context the parser does not have and would
// silently fabricate biologically-meaningful identifiers.

mod no_synthesis {
    use super::*;

    fn display_of(s: &str) -> String {
        format!("{}", parse_hgvs(s).unwrap())
    }

    #[test]
    fn chr1_stays_chr1() {
        assert_eq!(display_of("chr1:g.100A>G"), "chr1:g.100A>G");
    }

    #[test]
    fn bare_1_stays_bare_1() {
        assert_eq!(display_of("1:g.100A>G"), "1:g.100A>G");
    }

    #[test]
    fn grch38_chr1_stays_assembly_wrapped() {
        assert_eq!(display_of("GRCh38(chr1):g.100A>G"), "GRCh38(chr1):g.100A>G");
    }
}
