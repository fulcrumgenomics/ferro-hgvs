//! Regression test for issue #351 (folded into PR #343): the
//! `Boundaries::new(1, sequence_length())` callers in `normalize_tx` and
//! `normalize_genome` previously passed a 1-based-input left bound to
//! `shuffle()`'s 0-based contract, making transcript position 1 the
//! leftmost reachable 0-based new_start (i.e. transcript position 1 was
//! unreachable for 5'-direction shuffles).
//!
//! After the audit, every `Boundaries::new(1, X)` was converted to
//! `Boundaries::new(0, X)`. This test exercises a 5'-direction shuffle
//! on the `n.` axis where the canonical position SHOULD reach `n.1` —
//! pre-#351 the shuffle would stop one position short.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Non-coding transcript whose first 3 bases are all `A`. A 5'-direction
/// shuffle of `n.3del` walks back through the `AAA` tract and must reach
/// `n.1del`.
fn provider_with_aaa_prefix_noncoding() -> MockProvider {
    let mut provider = MockProvider::new();
    //                          123456789012
    let sequence = "AAACGTACGTAC".to_string();
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NR_LEFTBOUND.1".to_string(),
        Some("LEFTBOUND".to_string()),
        Strand::Plus,
        sequence,
        None,
        None,
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

#[test]
fn five_prime_shuffle_on_n_axis_reaches_position_1() {
    let normalizer = Normalizer::with_config(
        provider_with_aaa_prefix_noncoding(),
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let variant = parse_hgvs("NR_LEFTBOUND.1:n.3del").expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    assert_eq!(
        format!("{}", normalized),
        "NR_LEFTBOUND.1:n.1del",
        "5'-direction shuffle through the n.1..n.3 AAA-tract must reach n.1, \
         not stop at n.2 (the pre-#351 off-by-one)",
    );
}
