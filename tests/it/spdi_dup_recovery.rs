//! Integration test for SPDI→HGVS dup-form recovery (issue #119).
//!
//! Verifies that the public crate API (`ferro_hgvs::spdi_to_hgvs_with_ref`)
//! recovers HGVS `dup` form from canonical SPDI insertions when given a
//! reference, and leaves true insertions as `ins`.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::spdi::SpdiVariant;
use ferro_hgvs::{hgvs_to_spdi_simple, parse_hgvs, spdi_to_hgvs, spdi_to_hgvs_with_ref};

fn provider_with_atg_at_100() -> MockProvider {
    let mut contig = "N".repeat(99);
    contig.push_str("ATG"); // 1-based 100..102 = "ATG"
    contig.push_str(&"N".repeat(50));
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_000001.11", contig);
    p
}

#[test]
fn public_api_recovers_dup_from_spdi_insertion() {
    let provider = provider_with_atg_at_100();
    // SPDI 102::ATG is the canonical interbase form of g.100_102dupATG
    // under the post-#390 convention: position 102 is the boundary
    // AFTER 1-based base 102, matching the equivalent g.102_103ins.
    let spdi = SpdiVariant::insertion("NC_000001.11", 102, "ATG");
    let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
    assert_eq!(hgvs.to_string(), "NC_000001.11:g.100_102dupATG");
}

#[test]
fn public_api_full_roundtrip_dup() {
    let provider = provider_with_atg_at_100();
    let original = "NC_000001.11:g.100_102dupATG";

    let hgvs = parse_hgvs(original).unwrap();
    let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
    let recovered = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
    assert_eq!(recovered.to_string(), original);
}

#[test]
fn public_api_no_ref_path_unchanged() {
    // Audit pin at the public-API level: spdi_to_hgvs (no reference)
    // still emits ins form for a dup-shaped SPDI, even after issue
    // #119. SPDI input updated to 102::ATG to match the post-#390
    // canonical form of g.100_102dupATG; the ins-form rendering
    // remains g.102_103ins (SPDI 102 = boundary AFTER 1-based 102).
    let spdi = SpdiVariant::insertion("NC_000001.11", 102, "ATG");
    let hgvs = spdi_to_hgvs(&spdi).unwrap();
    assert_eq!(hgvs.to_string(), "NC_000001.11:g.102_103insATG");
}
