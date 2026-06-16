//! Tests for #270 — preserve `inv` and `m.` coord system through
//! HGVS ↔ SPDI round-trip.
//!
//! SPDI has no `inv` edit type and no coord-system tag, so the
//! forward leg lost these shapes. The reverse leg now recovers:
//!
//! - **inv recovery**: a delins where `ins == reverse_complement(del)`
//!   is canonically an inversion.
//! - **m. recovery**: NC_012920 / NC_001807 accessions emit `m.` instead
//!   of `g.`.
//!
//! `dup` recovery is unchanged (still requires reference data via
//! `spdi_to_hgvs_with_ref`).

use ferro_hgvs::parse_hgvs;
use ferro_hgvs::spdi::{hgvs_to_spdi_simple, spdi_to_hgvs, SpdiVariant};

#[track_caller]
fn assert_round_trips(input: &str) {
    let hgvs = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
    let spdi = hgvs_to_spdi_simple(&hgvs).unwrap_or_else(|e| panic!("to_spdi {input:?}: {e}"));
    let back = spdi_to_hgvs(&spdi).unwrap_or_else(|e| panic!("from_spdi {spdi:?}: {e}"));
    assert_eq!(
        format!("{}", back),
        input,
        "round-trip mismatch for {input:?}"
    );
}

// =============================================================================
// SECTION 1 — Inversion preservation
// =============================================================================

mod inversion_preservation {
    use super::*;

    /// Canonical 4-base inversion round-trips through SPDI.
    #[test]
    fn genomic_inv_round_trips() {
        assert_round_trips("NC_000001.11:g.100_103invTAGC");
    }

    /// Five-base inversion with non-palindromic seq.
    #[test]
    fn genomic_inv_5bp_round_trips() {
        assert_round_trips("NC_000001.11:g.100_104invATAGC");
    }

    /// Direct SPDI → HGVS where deleted = reverse_complement(inserted).
    #[test]
    fn direct_spdi_with_rc_emits_inversion() {
        let spdi = SpdiVariant::delins("NC_000001.11", 99, "TAGC", "GCTA");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(format!("{}", hgvs), "NC_000001.11:g.100_103invTAGC");
    }

    /// Length-1 inversion is just an SNV; identity SPDI doesn't
    /// canonicalize to inv (we require length >= 2).
    #[test]
    fn single_base_not_inversion() {
        let spdi = SpdiVariant::new("NC_000001.11", 99, "A", "T");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(format!("{}", hgvs), "NC_000001.11:g.100A>T");
    }

    /// Non-RC delins still emits delins.
    #[test]
    fn non_rc_delins_unchanged() {
        let spdi = SpdiVariant::delins("NC_000001.11", 99, "ATG", "TTCC");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(format!("{}", hgvs), "NC_000001.11:g.100_102delinsTTCC");
    }
}

// =============================================================================
// SECTION 2 — Mitochondrial coord-system preservation
// =============================================================================

mod mt_preservation {
    use super::*;

    /// Mitochondrial substitution preserves `m.`.
    #[test]
    fn mt_sub_round_trips() {
        assert_round_trips("NC_012920.1:m.100A>G");
    }

    /// Mitochondrial deletion preserves `m.`.
    #[test]
    fn mt_del_round_trips() {
        assert_round_trips("NC_012920.1:m.100_102delATG");
    }

    /// Mitochondrial inversion (combines both #270 features).
    #[test]
    fn mt_inv_round_trips() {
        assert_round_trips("NC_012920.1:m.100_103invTAGC");
    }

    /// Mitochondrial identity preserves `m.`.
    #[test]
    fn mt_identity_round_trips() {
        assert_round_trips("NC_012920.1:m.100A=");
    }

    /// Non-mitochondrial accessions still emit `g.`.
    #[test]
    fn nuclear_chrom_still_genomic() {
        assert_round_trips("NC_000001.11:g.12345A>G");
    }

    /// Older mitochondrial accession (NC_001807) also recognized.
    #[test]
    fn older_mt_accession_recognized() {
        let spdi = SpdiVariant::new("NC_001807.4", 99, "A", "G");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(format!("{}", hgvs), "NC_001807.4:m.100A>G");
    }
}

// =============================================================================
// SECTION 3 — Out-of-scope shapes unchanged
// =============================================================================

mod unchanged {
    use super::*;

    /// `dup` still canonicalizes to `ins` on the no-provider path.
    /// (Provider-required dup recovery is covered by `spdi_to_hgvs_with_ref`.)
    #[test]
    fn dup_canonicalizes_to_ins_without_provider() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_102dupATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(format!("{}", back), "NC_000001.11:g.102_103insATG");
    }

    /// `delins` with explicit del seq still drops the deleted seq on
    /// Display (canonical HGVS form).
    #[test]
    fn delins_with_explicit_del_drops_del_seq() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_102delATGinsTTCC").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(format!("{}", back), "NC_000001.11:g.100_102delinsTTCC");
    }
}
