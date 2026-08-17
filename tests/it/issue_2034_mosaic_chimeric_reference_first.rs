//! Mosaic (`/`) and chimeric (`//`) substitution alleles must render the
//! reference allele first (#2034).
//!
//! `DNA/substitution.md:49`: "irrespective of the frequency in which each
//! nucleotide was found, the reference is always described first." The clause is
//! non-normative prose (a house-choice the spec leaves ferro to answer), and the
//! operator adjudicated it **reference-first** — see the
//! `mosaic-chimeric-substitution-reference-allele-first` ruling in
//! `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`.
//!
//! So a variant-first spelling normalizes to reference-first:
//!   `g.5A>G/=`  → `g.5=/A>G`   (mosaic)
//!   `g.5A>G//=` → `g.5=//A>G`  (chimeric)
//! while a reference-first spelling is already a fixed point.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, NormalizeConfig, Normalizer};

use crate::common::hg38_window::{local_desc, provider};

/// Normalize a `NC_TESTWIN.1:g.<body>` description against the genomic window
/// provider and return the rendered string.
fn normalize_genomic(body: &str) -> String {
    let input = local_desc(body);
    let variant: HgvsVariant =
        parse_hgvs(&input).unwrap_or_else(|e| panic!("parse `{input}`: {e}"));
    Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
    .normalize(&variant)
    .unwrap_or_else(|e| panic!("normalize `{input}`: {e}"))
    .to_string()
}

// Local window base at local position 5 is `A` (window starts `CAAAAG…`), so
// `5A>G` states the reference base correctly and passes strict validation.

#[test]
fn genomic_mosaic_substitution_reorders_reference_first() {
    assert_eq!(normalize_genomic("5A>G/="), "NC_TESTWIN.1:g.5=/A>G");
}

#[test]
fn genomic_chimeric_substitution_reorders_reference_first() {
    assert_eq!(normalize_genomic("5A>G//="), "NC_TESTWIN.1:g.5=//A>G");
}

#[test]
fn genomic_reference_first_substitution_is_a_fixed_point() {
    assert_eq!(normalize_genomic("5=/A>G"), "NC_TESTWIN.1:g.5=/A>G");
    assert_eq!(normalize_genomic("5=//A>G"), "NC_TESTWIN.1:g.5=//A>G");
}

/// A mosaic of two *real* edits carries no reference allele, so there is
/// nothing to move to the front — it must be left in authored order.
#[test]
fn genomic_mosaic_of_two_real_edits_is_untouched() {
    // Bases 5 (`A`) and 6 (`G`) both state their reference correctly.
    let input = "NC_TESTWIN.1:g.5A>G/NC_TESTWIN.1:g.6G>C";
    let variant: HgvsVariant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse `{input}`: {e}"));
    let out = Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
    .normalize(&variant)
    .unwrap_or_else(|e| panic!("normalize `{input}`: {e}"))
    .to_string();
    assert_eq!(out, input);
}

/// Build a single-exon coding transcript whose CDS position `c.85` is a `T`,
/// so `c.85T>C` states its reference base correctly.
fn cds_provider() -> MockProvider {
    const UTR5_LEN: usize = 40;
    // 90-base CDS: `A`×84 then `T` at c.85 then `A`×5.
    let cds: String = format!("{}T{}", "A".repeat(84), "A".repeat(5));
    let sequence = format!("{}{cds}", "C".repeat(UTR5_LEN));
    let cds_start = (UTR5_LEN + 1) as u64;
    let cds_end = (UTR5_LEN + cds.len()) as u64;
    let tx_len = sequence.len() as u64;

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_TEST2034.1".to_string(),
        Some("TEST2034".to_string()),
        Strand::Plus,
        sequence,
        Some(cds_start),
        Some(cds_end),
        vec![Exon::new(1, 1, tx_len)],
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

fn normalize_cds(input: &str) -> String {
    let variant: HgvsVariant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse `{input}`: {e}"));
    Normalizer::with_config(
        cds_provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
    .normalize(&variant)
    .unwrap_or_else(|e| panic!("normalize `{input}`: {e}"))
    .to_string()
}

/// The issue's exact example, on a coding transcript.
#[test]
fn cds_mosaic_substitution_reorders_reference_first() {
    assert_eq!(
        normalize_cds("NM_TEST2034.1:c.85T>C/="),
        "NM_TEST2034.1:c.85=/T>C"
    );
}

#[test]
fn cds_chimeric_substitution_reorders_reference_first() {
    assert_eq!(
        normalize_cds("NM_TEST2034.1:c.85T>C//="),
        "NM_TEST2034.1:c.85=//T>C"
    );
}

#[test]
fn cds_reference_first_substitution_is_a_fixed_point() {
    assert_eq!(
        normalize_cds("NM_TEST2034.1:c.85=/T>C"),
        "NM_TEST2034.1:c.85=/T>C"
    );
}

/// Regression guard (#2107 review): the reorder must NOT fire when the second
/// member is a *ranged, non-co-located* identity (a screened region), because
/// reordering would rebuild the first member from the substitution's location
/// and **silently discard** that screened region. `c.85T>C/c.60_70=` denotes
/// two distinct facts — a substitution at c.85 and an identity assertion over
/// c.60_70 — and must round-trip verbatim (it did pre-#2034).
#[test]
fn cds_mosaic_ranged_identity_second_member_is_untouched() {
    let input = "NM_TEST2034.1:c.85T>C/NM_TEST2034.1:c.60_70=";
    assert_eq!(normalize_cds(input), input);
}

#[test]
fn cds_chimeric_ranged_identity_second_member_is_untouched() {
    let input = "NM_TEST2034.1:c.85T>C//NM_TEST2034.1:c.60_70=";
    assert_eq!(normalize_cds(input), input);
}
