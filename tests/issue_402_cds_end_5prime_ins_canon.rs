//! Regression for #402.
//!
//! Surfaced by `tests/fixtures/biocommons-normalize/baseline-failures/normalized.txt`:
//! `NM_000051.3:c.9171_*1insA` (5'+cross) should canonicalise to
//! `NM_000051.3:c.9171dup` per biocommons; ferro currently emits
//! `NM_000051.3:c.9171_*1insA` (the input unchanged).
//!
//! Geometry of the failing case (derived from the adjacent passing row
//! `NM_000051.3:c.9170_9171insAT` 3'+cross → `c.9171_*1dup`):
//!
//!   axis: ... c.9170 c.9171 | c.*1 ...
//!   base: ...  ?      A     |  T   ...
//!                            ^
//!                            cds_end boundary
//!
//! Inserting `A` between c.9171 (A) and c.*1 (T) duplicates the
//! immediately-5' CDS base — the canonical form is the CDS-resident
//! `c.9171dup`. ferro emits the input unchanged: the post-shift
//! ins → dup recogniser does not fire for an insertion straddling the
//! CDS-end boundary even though the inserted alt exactly matches the
//! immediately-5' base.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Synthetic transcript whose c.cds_end = A and c.*1 = T — same shape
/// as `NM_000051.3` around the failing input's boundary.
///
///   axis: c.-3 c.-2 c.-1 | c.1 c.2 c.3 c.4 c.5 | c.*1 c.*2 c.*3 ...
///   base:  G    G    G   | C   G   T   C   A   | T    G    G    ...
///
/// CDS spans c.1..c.5 (synthetic; the bug is independent of CDS
/// length). The c.-1 base is G so the synthetic does not accidentally
/// trip the spanning-dup canon at the CDS-start (irrelevant here).
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let mut seq = String::from("GGG"); // 5'UTR (c.-3..c.-1)
    seq.push_str("CGTCA"); // c.1..c.5 (CDS); c.5 = A is cds_end
    seq.push_str("TGG"); // c.*1..c.*3 (3'UTR); c.*1 = T
    while seq.len() < 203 {
        seq.push('G');
    }
    let len = seq.len() as u64;
    // c.1 byte index 3 → 1-based tx pos 4; c.5 byte index 7 → tx pos 8.
    let transcript = Transcript::new(
        "NM_TEST402.1".to_string(),
        Some("TEST402".to_string()),
        Strand::Plus,
        seq,
        Some(4),
        Some(8),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

fn normalize_with(direction: ShuffleDirection, input: &str) -> String {
    let normalizer = Normalizer::with_config(
        provider(),
        NormalizeConfig::default()
            .with_direction(direction)
            .allow_crossing_boundaries(),
    );
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{normalized}")
}

#[test]
fn five_prime_insertion_at_cds_end_boundary_canonicalises_to_cds_dup() {
    // c.5_*1insA on a transcript whose c.5 = A: biocommons emits
    // c.5dup (the inserted A duplicates the immediately-5' CDS base);
    // ferro currently leaves the input unchanged.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST402.1:c.5_*1insA"),
        "NM_TEST402.1:c.5dup",
    );
}

#[test]
fn three_prime_direction_unchanged_on_same_input() {
    // 3'-direction sanity. The same insertion under 3'-direction has
    // no further-3' tract to walk into (c.*1 = T ≠ A), so the canonical
    // form is also c.5dup — the immediately-5' tandem unit. Locks in
    // direction symmetry for this fix.
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, "NM_TEST402.1:c.5_*1insA"),
        "NM_TEST402.1:c.5dup",
    );
}
