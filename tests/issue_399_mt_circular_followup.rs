//! End-to-end coverage for issue #399 — the F1 follow-up to PR #380.
//!
//! These tests exercise wraparound-aware span math through the public
//! `_with_provider` APIs. The no-provider `get_indel_length` path
//! returns `None` for wraparound on `m.`/`o.` (pinned by
//! `tests/mito_circular_audit.rs`); these tests pin the spec-correct
//! values when a provider can supply contig length.
//!
//! The second section (vcf_conversion_*) verifies that `HgvsToVcfConverter`
//! rejects wraparound `m.`/`o.` variants with a clear error rather than
//! silently producing garbage from calling `get_reference_sequence` with
//! `start > end`.
//!
//! The third section (spdi_conversion_*) verifies that `hgvs_to_spdi` and
//! `hgvs_to_spdi_simple` reject wraparound `m.`/`o.` variants: SPDI is a
//! single-edit format with no native representation for circular-contig
//! wraparound.

use ferro_hgvs::parse_hgvs;
use ferro_hgvs::python_helpers::get_indel_length_with_provider;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::spdi::convert::{hgvs_to_spdi, hgvs_to_spdi_simple};
use ferro_hgvs::vcf::HgvsToVcfConverter;

fn mt_provider() -> MockProvider {
    let mut p = MockProvider::new();
    // NC_012920.1 is 16569 bp; for span math we only need the length,
    // not the bases. Use a placeholder string of the right length.
    p.add_genomic_sequence("NC_012920.1", "A".repeat(16569));
    p
}

#[test]
fn wraparound_del_indel_length_is_spec_correct_with_provider() {
    let v = parse_hgvs("NC_012920.1:m.16569_1del").unwrap();
    let p = mt_provider();
    // 2-nt deletion (positions 16569 and 1); del reports -span.
    assert_eq!(get_indel_length_with_provider(&v, &p), Some(-2));
}

#[test]
fn wraparound_dup_indel_length_is_spec_correct_with_provider() {
    let v = parse_hgvs("NC_012920.1:m.16560_5dup").unwrap();
    let p = mt_provider();
    // Wraparound span = (16569 − 16560 + 1) + 5 = 15.
    assert_eq!(get_indel_length_with_provider(&v, &p), Some(15));
}

#[test]
fn wraparound_delins_indel_length_is_spec_correct_with_provider() {
    let v = parse_hgvs("NC_012920.1:m.16569_1delinsT").unwrap();
    let p = mt_provider();
    // 2-nt window replaced with 1-nt insert; net = 1 − 2 = −1.
    assert_eq!(get_indel_length_with_provider(&v, &p), Some(-1));
}

// =============================================================================
// VCF conversion: wraparound m./o. variants must be rejected at the boundary
// =============================================================================

/// Minimal transcript stub for VCF converter tests.
///
/// `HgvsToVcfConverter` needs a `&Transcript` but for `Mt`/`Circular` variants
/// `convert_genome` derives the chromosome from the accession, not from the
/// transcript. The stub only needs to be structurally valid.
fn minimal_transcript() -> Transcript {
    Transcript::new(
        "NC_012920.1".to_string(),
        None,
        Strand::Plus,
        None,
        None,
        None,
        vec![],
        Some("chrM".to_string()),
        None,
        None,
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    )
}

#[test]
fn vcf_conversion_rejects_wraparound_mt_del() {
    let v = parse_hgvs("NC_012920.1:m.16569_1del").unwrap();
    let tx = minimal_transcript();
    let p = mt_provider();
    let converter = HgvsToVcfConverter::new(&tx, &p);
    let err = converter.convert(&v).unwrap_err();
    let msg = format!("{err}");
    assert!(
        msg.contains("wraparound"),
        "expected error mentioning 'wraparound', got: {msg}"
    );
}

#[test]
fn vcf_conversion_rejects_wraparound_circular_dup() {
    // SVD-WG006's o.-axis example: J01749.1:o.4344_197dup wraps the origin
    // on a 5386-bp plasmid. The Circular arm of the converter must reject
    // it with the same shape of error as the Mt arm.
    let v = parse_hgvs("J01749.1:o.4344_197dup").unwrap();
    let tx = minimal_transcript();
    let mut p = MockProvider::new();
    p.add_genomic_sequence("J01749.1", "A".repeat(5386));
    let converter = HgvsToVcfConverter::new(&tx, &p);
    let err = converter.convert(&v).unwrap_err();
    let msg = format!("{err}");
    assert!(
        msg.contains("wraparound") && msg.contains("o."),
        "expected error mentioning 'wraparound' and 'o.', got: {msg}"
    );
}

#[test]
fn vcf_conversion_rejects_wraparound_mt_delins() {
    let v = parse_hgvs("NC_012920.1:m.16569_1delinsT").unwrap();
    let tx = minimal_transcript();
    let p = mt_provider();
    let converter = HgvsToVcfConverter::new(&tx, &p);
    assert!(converter.convert(&v).is_err());
}

#[test]
fn vcf_conversion_rejects_wraparound_mt_dup() {
    let v = parse_hgvs("NC_012920.1:m.16560_5dup").unwrap();
    let tx = minimal_transcript();
    let p = mt_provider();
    let converter = HgvsToVcfConverter::new(&tx, &p);
    assert!(converter.convert(&v).is_err());
}

#[test]
fn vcf_conversion_still_accepts_linear_mt_sub() {
    let v = parse_hgvs("NC_012920.1:m.100A>G").unwrap();
    let tx = minimal_transcript();
    let p = mt_provider();
    let converter = HgvsToVcfConverter::new(&tx, &p);
    assert!(converter.convert(&v).is_ok());
}

// =============================================================================
// SPDI conversion: wraparound m./o. variants must be rejected at the boundary
// =============================================================================
// SPDI is a single-edit format defined as (sequence, position, deleted, inserted).
// It has no native representation for circular-contig wraparound, so the
// conversion must reject those inputs with a clear error rather than emitting
// a SPDI with start > end that downstream consumers cannot interpret.

#[test]
fn spdi_conversion_rejects_wraparound_mt_del() {
    let v = parse_hgvs("NC_012920.1:m.16569_1del").unwrap();
    let p = mt_provider();
    let err = hgvs_to_spdi(&v, &p).unwrap_err();
    let msg = format!("{err}");
    assert!(
        msg.contains("wraparound"),
        "expected error mentioning 'wraparound', got: {msg}"
    );
}

#[test]
fn spdi_conversion_rejects_wraparound_mt_dup() {
    let v = parse_hgvs("NC_012920.1:m.16560_5dup").unwrap();
    let p = mt_provider();
    assert!(hgvs_to_spdi(&v, &p).is_err());
}

#[test]
fn spdi_conversion_still_accepts_linear_mt() {
    let v = parse_hgvs("NC_012920.1:m.100A>G").unwrap();
    let p = mt_provider();
    assert!(hgvs_to_spdi(&v, &p).is_ok());
}

#[test]
fn spdi_conversion_simple_accepts_linear_circular() {
    // Non-wraparound o. variants now convert successfully via the new
    // circular_to_spdi_simple path added alongside the wraparound guard.
    // Pins the positive side of the new o. dispatch.
    let v = parse_hgvs("J01749.1:o.100A>G").unwrap();
    assert!(hgvs_to_spdi_simple(&v).is_ok());
}

#[test]
fn spdi_conversion_simple_rejects_wraparound_mt_del() {
    // hgvs_to_spdi_simple also converts m. variants and must apply the same guard.
    let v = parse_hgvs("NC_012920.1:m.16569_1del").unwrap();
    let err = hgvs_to_spdi_simple(&v).unwrap_err();
    let msg = format!("{err}");
    assert!(
        msg.contains("wraparound"),
        "expected error mentioning 'wraparound', got: {msg}"
    );
}

#[test]
fn spdi_conversion_rejects_wraparound_circular() {
    // SVD-WG006's o.-axis example: J01749.1:o.4344_197dup wraps the origin
    // on a 5386-bp plasmid.
    let v = parse_hgvs("J01749.1:o.4344_197dup").unwrap();
    let mut p = MockProvider::new();
    p.add_genomic_sequence("J01749.1", "A".repeat(5386));
    let err = hgvs_to_spdi(&v, &p).unwrap_err();
    let msg = format!("{err}");
    assert!(
        msg.contains("wraparound") && msg.contains("o."),
        "expected error mentioning 'wraparound' and 'o.', got: {msg}"
    );
}

// =============================================================================
// Validate handler: wraps_origin field on ValidateResponse
// =============================================================================

#[cfg(feature = "web-service")]
mod validate_response_tests {
    use ferro_hgvs::service::handlers::validate::validate_hgvs;

    #[test]
    fn validate_response_wraps_origin_true_for_wraparound_mt() {
        let response = validate_hgvs("NC_012920.1:m.16569_1del");
        assert!(
            response.wraps_origin,
            "expected wraps_origin=true on wraparound m. del"
        );
    }

    #[test]
    fn validate_response_wraps_origin_true_for_wraparound_circular() {
        let response = validate_hgvs("J01749.1:o.4344_197dup");
        assert!(
            response.wraps_origin,
            "expected wraps_origin=true on wraparound o. dup"
        );
    }

    #[test]
    fn validate_response_wraps_origin_false_for_linear_mt() {
        let response = validate_hgvs("NC_012920.1:m.100A>G");
        assert!(
            !response.wraps_origin,
            "expected wraps_origin=false on linear m. sub"
        );
    }
}
