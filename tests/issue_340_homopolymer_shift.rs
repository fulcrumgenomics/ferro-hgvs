//! Regression test for issue #340 subgroup (i): 5'-direction shuffle of
//! a single-base insertion in a homopolymer tract picks the 3'-most
//! position instead of the 5'-most. Biocommons cases (from #324):
//!
//!   NM_000051.3:c.14_15insT   bio: c.14dup   ferro: c.15dup   (5prime)
//!   NM_000051.3:c.-4_-3insA   bio: c.-4dup   ferro: c.-3dup   (5prime)
//!   NM_000051.3:c.9171_*1insA bio: c.9171dup ferro: c.9170dup (5prime)
//!   NM_000051.3:c.*4_*5insT   bio: c.*3dup   ferro: c.*4dup   (5prime)
//!
//! Self-contained synthetic transcript pins the pattern without needing
//! the real NM_000051.3 manifest. The (ii) delins → dup and (iii)
//! CDS-start clamp subgroups in the issue body are left for follow-up
//! PRs per the issue author's own scoping note.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Synthetic transcript:
///   tx positions 1-20  (CDS-only, no 5'UTR / 3'UTR — avoids overlap
///                       with the axis-clamp work in #337)
///   c.1..c.20 == "ATGAATTAAATAGAAAAATA"
///                       ↑   ↑↑      — c.5=A, c.6=T, c.7=T, c.8=A
///
/// `c.6_7insT` inserts another T into the 2-T tract at c.6-c.7,
/// growing it to a 3-T tract at c.6, c.7, c.8 (or 5, 6, 7 depending
/// on shift direction). Under 5prime, the canonical dup is c.6dup;
/// under 3prime, c.7dup. ferro currently picks 3prime regardless.
///
/// We use a homopolymer of length 2 specifically so the shift has
/// exactly one position to cover — making the off-by-one isolable.
fn provider_with_t_tract_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    //                       12345678901234567890
    let sequence = String::from("ATGAATTAAATAGAAAAATA");
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(1),
        Some(20),
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

fn normalize_with(direction: ShuffleDirection, cross: bool, input: &str) -> String {
    let mut config = NormalizeConfig::default().with_direction(direction);
    if cross {
        config = config.allow_crossing_boundaries();
    }
    let normalizer = Normalizer::with_config(provider_with_t_tract_transcript(), config);
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{}", normalized)
}

#[test]
fn five_prime_insertion_in_homopolymer_picks_leftmost_dup() {
    // c.6_7insT into the c.6-c.7 T-tract (c.5=A, c.6=T, c.7=T, c.8=A).
    // 5prime-most dup is c.6dup. ferro pre-fix emits c.7dup.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, false, "NM_TEST.1:c.6_7insT"),
        "NM_TEST.1:c.6dup",
    );
}

#[test]
fn three_prime_insertion_in_homopolymer_picks_rightmost_dup() {
    // Sanity check that the 3prime direction still picks c.7dup.
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, false, "NM_TEST.1:c.6_7insT"),
        "NM_TEST.1:c.7dup",
    );
}

#[test]
fn five_prime_cross_insertion_in_homopolymer_picks_leftmost_dup() {
    // cross_boundaries=true must NOT change the homopolymer dup position;
    // it only affects exon-intron crossing. Single-exon transcript here
    // has no introns, but the test pins behavior under both cross modes.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, true, "NM_TEST.1:c.6_7insT"),
        "NM_TEST.1:c.6dup",
    );
}
