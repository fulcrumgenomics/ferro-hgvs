//! Audit for #81 § K1 (remaining) — HGVS ↔ SPDI coverage matrix.
//!
//! Three categories of behavior are pinned:
//! 1. Syntactic round-trip: HGVS → SPDI → HGVS preserves the input verbatim.
//! 2. Semantic round-trip: SPDI canonicalizes the HGVS form (e.g. dup→ins,
//!    delins drops the deleted seq). Pinned as the expected canonical
//!    re-emission, not a bug. `inv` and `m.` are recovered by inspection on
//!    the reverse leg and round-trip verbatim (see syntactic section).
//! 3. Provider-required / unsupported: returns a typed error.
//!
//! Out of scope: length-mismatch input (`g.100_110delAAAATTTGCC` with
//! 10-base seq across 11 positions) — the parser silently accepts it and
//! the SPDI conversion clips to the seq length. Tracked separately.

use ferro_hgvs::parse_hgvs;
use ferro_hgvs::spdi::{hgvs_to_spdi_simple, spdi_to_hgvs, ConversionError};

#[track_caller]
fn assert_syntactic_round_trip(s: &str) {
    let hgvs = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?}: {e}"));
    let spdi = hgvs_to_spdi_simple(&hgvs).unwrap_or_else(|e| panic!("to_spdi {s:?}: {e}"));
    let back = spdi_to_hgvs(&spdi).unwrap_or_else(|e| panic!("from_spdi {spdi:?}: {e}"));
    let back_s = format!("{}", back);
    assert_eq!(back_s, s, "syntactic round-trip mismatch for {s:?}");
}

#[track_caller]
fn assert_semantic_round_trip(input: &str, canonical: &str) {
    let hgvs = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
    let spdi = hgvs_to_spdi_simple(&hgvs).unwrap_or_else(|e| panic!("to_spdi {input:?}: {e}"));
    let back = spdi_to_hgvs(&spdi).unwrap_or_else(|e| panic!("from_spdi {spdi:?}: {e}"));
    let back_s = format!("{}", back);
    assert_eq!(
        back_s, canonical,
        "semantic round-trip mismatch for {input:?}; expected canonical {canonical:?}"
    );
}

// =============================================================================
// SECTION 1 — Syntactic round-trip
// =============================================================================

mod syntactic_round_trip {
    use super::*;

    #[test]
    fn substitution_round_trips() {
        assert_syntactic_round_trip("NC_000001.11:g.12345A>G");
    }

    #[test]
    fn insertion_round_trips() {
        assert_syntactic_round_trip("NC_000001.11:g.100_101insATG");
    }

    #[test]
    fn single_base_del_with_seq_round_trips() {
        assert_syntactic_round_trip("NC_000001.11:g.100delA");
    }

    #[test]
    fn range_del_with_seq_round_trips() {
        assert_syntactic_round_trip("NC_000001.11:g.100_102delATG");
    }

    #[test]
    fn identity_round_trips() {
        assert_syntactic_round_trip("NC_000001.11:g.100A=");
    }

    #[test]
    fn lrg_substitution_round_trips() {
        assert_syntactic_round_trip("LRG_199:g.100A>G");
    }

    #[test]
    fn lrg_del_with_seq_round_trips() {
        assert_syntactic_round_trip("LRG_199:g.100_102delATG");
    }

    /// SPDI itself has no `inv` edit type — the forward leg encodes it as
    /// the reverse-complement substring. The reverse leg recognises the
    /// substring as a palindromic inversion and recovers `inv` verbatim.
    #[test]
    fn inv_with_seq_round_trips() {
        assert_syntactic_round_trip("NC_000001.11:g.100_104invTAGCA");
    }

    /// SPDI has no coord-system tag, but `NC_012920.1` is the canonical
    /// mitochondrial accession, so the reverse leg recovers `m.` verbatim.
    #[test]
    fn mt_substitution_round_trips() {
        assert_syntactic_round_trip("NC_012920.1:m.100A>G");
    }
}

// =============================================================================
// SECTION 2 — Semantic round-trip (canonicalization through SPDI)
// =============================================================================

mod semantic_round_trip {
    use super::*;

    /// SPDI has no `dup` concept; a duplication encodes as an insertion
    /// after the duplicated region. Reverse path emits the insertion form.
    #[test]
    fn dup_canonicalizes_to_insertion_after_region() {
        assert_semantic_round_trip(
            "NC_000001.11:g.100_102dupATG",
            "NC_000001.11:g.102_103insATG",
        );
    }

    /// `delATGinsTTCC` SPDI-encodes as `ATG:TTCC`; reverse drops the
    /// explicit deleted seq (still stored in SPDI but not re-emitted).
    #[test]
    fn delins_canonicalizes_to_implicit_deleted_seq() {
        assert_semantic_round_trip(
            "NC_000001.11:g.100_102delATGinsTTCC",
            "NC_000001.11:g.100_102delinsTTCC",
        );
    }
}

// =============================================================================
// SECTION 3 — Provider-required (errors with MissingReferenceData)
// =============================================================================

mod provider_required {
    use super::*;

    fn try_to_spdi(s: &str) -> Result<(), ConversionError> {
        let hgvs = parse_hgvs(s).unwrap();
        hgvs_to_spdi_simple(&hgvs).map(|_| ())
    }

    #[test]
    fn del_without_seq_requires_provider() {
        assert!(matches!(
            try_to_spdi("NC_000001.11:g.100_102del"),
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    #[test]
    fn dup_without_seq_requires_provider() {
        assert!(matches!(
            try_to_spdi("NC_000001.11:g.100_102dup"),
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    #[test]
    fn inv_without_seq_requires_provider() {
        assert!(matches!(
            try_to_spdi("NC_000001.11:g.100_105inv"),
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    #[test]
    fn delins_without_del_seq_requires_provider() {
        assert!(matches!(
            try_to_spdi("NC_000001.11:g.100_102delinsTTCC"),
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    #[test]
    fn single_base_del_without_seq_requires_provider() {
        assert!(matches!(
            try_to_spdi("NC_000001.11:g.100del"),
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    #[test]
    fn cds_variant_requires_provider() {
        assert!(matches!(
            try_to_spdi("NM_000088.3:c.100A>G"),
            Err(ConversionError::ProviderRequired { .. })
        ));
    }

    #[test]
    fn repeat_without_unit_requires_provider() {
        assert!(matches!(
            try_to_spdi("NC_000001.11:g.100_102[5]"),
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }
}

// =============================================================================
// SECTION 4 — Unsupported variant types (errors with UnsupportedVariantType)
// =============================================================================

mod unsupported_types {
    use super::*;

    fn try_to_spdi(s: &str) -> Result<(), ConversionError> {
        let hgvs = parse_hgvs(s).unwrap();
        hgvs_to_spdi_simple(&hgvs).map(|_| ())
    }

    #[test]
    fn protein_variant_unsupported() {
        assert!(matches!(
            try_to_spdi("NP_003997.1:p.Arg8Gln"),
            Err(ConversionError::UnsupportedVariantType { .. })
        ));
    }

    #[test]
    fn compound_allele_unsupported() {
        assert!(matches!(
            try_to_spdi("NC_000001.11:g.[100A>G;200C>T]"),
            Err(ConversionError::UnsupportedVariantType { .. })
        ));
    }
}

// =============================================================================
// SECTION 5 — Direct SPDI inputs (no HGVS source)
// =============================================================================
//
// Pin the SPDI → HGVS direction for the canonical SPDI shapes, locking
// the inverse of the syntactic round-trip section.

mod spdi_to_hgvs_direct {
    use ferro_hgvs::spdi::{spdi_to_hgvs, SpdiVariant};

    #[test]
    fn sub_emits_sub() {
        let s = SpdiVariant::new("NC_000001.11", 12344, "A", "G");
        assert_eq!(
            spdi_to_hgvs(&s).unwrap().to_string(),
            "NC_000001.11:g.12345A>G"
        );
    }

    #[test]
    fn range_deletion_emits_del_with_seq() {
        let s = SpdiVariant::deletion("NC_000001.11", 99, "ATG");
        assert_eq!(
            spdi_to_hgvs(&s).unwrap().to_string(),
            "NC_000001.11:g.100_102delATG"
        );
    }

    #[test]
    fn single_base_deletion_emits_del_with_seq() {
        let s = SpdiVariant::deletion("NC_000001.11", 99, "A");
        assert_eq!(
            spdi_to_hgvs(&s).unwrap().to_string(),
            "NC_000001.11:g.100delA"
        );
    }

    #[test]
    fn insertion_emits_ins() {
        // SPDI 100 is the 0-based interbase position AFTER 1-based base
        // 100 → HGVS `g.100_101ins…` (#390 corrects a prior off-by-one
        // that emitted `g.101_102ins…`).
        let s = SpdiVariant::insertion("NC_000001.11", 100, "ATG");
        assert_eq!(
            spdi_to_hgvs(&s).unwrap().to_string(),
            "NC_000001.11:g.100_101insATG"
        );
    }

    #[test]
    fn delins_emits_delins() {
        let s = SpdiVariant::delins("NC_000001.11", 99, "ATG", "TTCC");
        assert_eq!(
            spdi_to_hgvs(&s).unwrap().to_string(),
            "NC_000001.11:g.100_102delinsTTCC"
        );
    }

    #[test]
    fn identity_emits_eq() {
        let s = SpdiVariant::new("NC_000001.11", 99, "A", "A");
        assert_eq!(
            spdi_to_hgvs(&s).unwrap().to_string(),
            "NC_000001.11:g.100A="
        );
    }
}
