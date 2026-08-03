//! Regression for #1349: a 5' normalization that was not a fixed point.
//!
//! `NM_TEST1349.1:c.2_3insTCC` normalized (5') to `c.1delinsCCTC`, and
//! re-normalizing *that* gave `c.-1_2dup` — a stable answer only on the third
//! pass. The 3' direction reached `c.-1_2dup` immediately, so the two directions
//! also disagreed on the first pass about an edit that is direction-independent.
//!
//! Root cause is in `rules::insertion_to_duplication`, not in the CDS-start
//! clamp. That function tries every cyclic rotation of the inserted unit, scores
//! each by the size of the reference tandem tract it abuts, and keeps the single
//! highest-scoring one — applying the #882 equivalence check
//! (`tandem_edit_preserves_insertion`) only afterwards. Here the alt `TCC`
//! abuts a one-copy tract (the `TCC` at c.-1_2, which is the correct answer),
//! while its rotation `CCT` abuts a *two-copy* tract further 3'. `CCT` wins on
//! count, fails the equivalence check because that tract is out of phase with
//! the insertion, and the function returns `None` — discarding the correct
//! `TCC` candidate, which was never reconsidered.
//!
//! With no `ins → dup` promotion the 5' path emitted a plain insertion resting
//! at the CDS boundary, which #383's clamp then rewrote as `c.1delinsCCTC`. The
//! clamp was doing its job on the input it was handed; the promotion upstream is
//! what was missing. Once the dup is recognised, #401's `spanning_dup_exception`
//! already lets it through the clamp unchanged.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// 5'UTR `…GTACGT` followed by a CDS opening `CCCCTCCTCCTCCTCCCGGGG`.
///
///   axis: c.-6 … c.-2 c.-1 | c.1 c.2 c.3 c.4 c.5 c.6 c.7 c.8 …
///   base:  G   A  C   G    T  | C   C   C   C   T   C   C   T  …
///
/// The two tracts that compete for `c.2_3insTCC`:
///   - `TCC` at c.-1_2 (one copy) — abuts the insertion point, in phase.
///   - `CCT` at c.3_5 and c.6_8 (two copies) — out of phase with the insertion.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let mut seq = String::from("GTACGT"); // 5'UTR: c.-6 … c.-1
    seq.push_str("CCCCTCCTCCTCCTCCCGGGG"); // c.1 … c.21
    while seq.len() < 6 + 30 {
        seq.push('A'); // pad the CDS out to 30 bases (10 codons)
    }
    let cds_end = seq.len() as u64;
    seq.push_str("AAAAAAAAAA"); // 3'UTR
    let len = seq.len() as u64;
    let transcript = Transcript::new(
        "NM_TEST1349.1".to_string(),
        Some("TEST1349".to_string()),
        Strand::Plus,
        seq,
        Some(7), // c.1 is the 7th base (1-based)
        Some(cds_end),
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

/// The inserted `TCC` duplicates the `TCC` at c.-1_2, so the canonical form is
/// that spanning dup — `dup` outranks both `ins` and `delins` in the spec's
/// prioritisation, and the dup-source genuinely equals the alt.
#[test]
fn five_prime_insertion_promotes_to_the_in_phase_dup() {
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST1349.1:c.2_3insTCC"),
        "NM_TEST1349.1:c.-1_2dup",
    );
}

/// The 3' direction already reached this answer (via the post-shuffle simple-dup
/// path rather than the rotation search), and must keep reaching it.
#[test]
fn three_prime_insertion_promotes_to_the_in_phase_dup() {
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, "NM_TEST1349.1:c.2_3insTCC"),
        "NM_TEST1349.1:c.-1_2dup",
    );
}

/// The defect's actual signature: `norm(norm(x)) != norm(x)`. Pinned separately
/// from the value assertions above so a future change that moves the canonical
/// answer still has to keep it a fixed point.
#[test]
fn five_prime_normalization_is_a_fixed_point() {
    let once = normalize_with(ShuffleDirection::FivePrime, "NM_TEST1349.1:c.2_3insTCC");
    let twice = normalize_with(ShuffleDirection::FivePrime, &once);
    assert_eq!(once, twice, "normalization must be idempotent");
}

/// The boundary `delins` that the clamp used to emit is still reachable by
/// writing it directly, and still normalizes to the dup — so the fix did not
/// simply disable the clamp.
#[test]
fn the_equivalent_boundary_delins_still_normalizes_to_the_same_dup() {
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST1349.1:c.1delinsCCTC"),
        "NM_TEST1349.1:c.-1_2dup",
    );
}
