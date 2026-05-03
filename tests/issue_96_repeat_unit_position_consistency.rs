//! Regression test for issue #96 (residual after #98).
//!
//! Per HGVS spec for repeated sequences, the position in
//! `prefix:position_or_range unit[N]` must locate the **first repeat
//! unit** — its length must match the unit length. ferro previously
//! emitted the location of the **entire reference tract** (ref count ×
//! unit length), which is a longer span than the unit.
//!
//! Examples (from `varnomen.hgvs.org/recommendations/DNA/variant/repeated/`):
//!
//! - `NC_000001.11:g.123456_123458GCA[12]` — 3-base position for a
//!   3-base unit.
//! - Single-base unit: position is a single base.
//!
//! These tests fail before the fix in `duplication_to_repeat` and
//! `insertion_to_repeat` (which both used `start..ref_end_of_full_tract`)
//! and pass after the fix sets `end = start + unit_len - 1`.
//!
//! Independent of the symptom-locked tests in
//! `tests/coverage_gap_tests.rs` and `tests/ins_shift_matrix.rs`.
//!
//! The fixture is a single-exon transcript whose sequence contains a
//! contiguous tandem run, so the position is computable purely from the
//! transcript-view bases.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

fn provider_with_single_exon(id: &str, sequence: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    let len = sequence.len() as u64;
    provider.add_transcript(Transcript::new(
        id.to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        sequence.to_string(),
        Some(1),
        Some(len),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        GenomeBuild::Unknown,
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

fn normalize(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).expect("parse should succeed");
    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize should succeed");
    format!("{}", normalized)
}

/// **Diagnostic 1 — single-base unit, dup → repeat.**
///
/// Tx `CAAACG` (c.1..6 = C,A,A,A,C,G); AAA tract at c.2..4. Dup
/// `c.2_4` adds 3 more As → 6 As total → repeat notation with unit
/// `A`. The position must locate ONE A — a single base — not the
/// original 3-base tract.
#[test]
fn dup_single_base_unit_position_is_single_base() {
    let provider = provider_with_single_exon("NM_TEST.1", "CAAACG");
    let result = normalize(provider, "NM_TEST.1:c.2_4dup");
    assert_eq!(result, "NM_TEST.1:c.2A[6]");
}

/// **Diagnostic 2 — multi-base unit, dup → repeat.**
///
/// Tx `CGCAGCAGCAT` (c.1..11). GCA tract at c.2..10 (3 units).
/// Duplicating two GCA units (`c.3_8dup` = `AGCAGC` shifted to one
/// of the GCA boundaries) triggers repeat notation with unit `GCA`.
/// The position must locate the FIRST GCA unit — 3 bases starting
/// at c.2 — i.e. `c.2_4`, not the full 9-base tract.
#[test]
fn dup_three_base_unit_position_is_three_bases() {
    let provider = provider_with_single_exon("NM_TEST.1", "CGCAGCAGCAT");
    let result = normalize(provider, "NM_TEST.1:c.3_8dup");
    assert_eq!(result, "NM_TEST.1:c.2_4GCA[5]");
}

/// **Diagnostic 3 — single-base unit, ins → repeat.**
///
/// Tx `CAAACG`, AAA tract at c.2..4. Inserting two more As
/// (`c.3_4insAA`) adds 2 units of `A` → 5 As total → repeat
/// notation. Position must be a single base.
#[test]
fn ins_single_base_unit_position_is_single_base() {
    let provider = provider_with_single_exon("NM_TEST.1", "CAAACG");
    let result = normalize(provider, "NM_TEST.1:c.3_4insAA");
    assert_eq!(result, "NM_TEST.1:c.2A[5]");
}

/// **Diagnostic 4 — multi-base unit, ins → repeat.**
///
/// Tx `CGCAGCAGCAT`, GCA tract at c.2..10 (3 units). Inserting two
/// more GCA units (`c.5_6insGCAGCA`) adds 2 units → 5 GCA total →
/// repeat. Position must be 3 bases (the first GCA at c.2_4), not
/// the full 9-base tract.
#[test]
fn ins_three_base_unit_position_is_three_bases() {
    let provider = provider_with_single_exon("NM_TEST.1", "CGCAGCAGCAT");
    let result = normalize(provider, "NM_TEST.1:c.5_6insGCAGCA");
    assert_eq!(result, "NM_TEST.1:c.2_4GCA[5]");
}
