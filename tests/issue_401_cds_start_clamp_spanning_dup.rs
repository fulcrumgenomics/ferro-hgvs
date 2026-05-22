//! Regression for #401.
//!
//! Surfaced by `tests/fixtures/biocommons-normalize/baseline-failures/normalized.txt`:
//! `NM_000051.3:c.1_2insCA` (5'+cross) should canonicalise to
//! `NM_000051.3:c.-1_1dup` per biocommons; ferro currently emits
//! `NM_000051.3:c.1delinsACA` (the PR #385 / #383 clamp output, applied
//! where the canonical form is actually a *spanning* dup — one endpoint
//! in 5'UTR (c.-1), one endpoint in CDS (c.1)).
//!
//! Root cause: PR #385's CDS-start canon clamp fires unconditionally on
//! `start_axis == Cds && new_tx_start < cds_start`. When the
//! canonicalisation produces a `Duplication` whose `tx_end >= cds_start`,
//! the dup-source spans the c.-1/c.1 boundary and IS the spec-canonical
//! form (HGVS DNA §general; edit-type priority `dup > ins`). The clamp
//! must skip these spanning duplications and keep the dup output.
//!
//! The PR #385 clamp continues to fire for:
//!   - Entirely-UTR rewrites (`new_tx_end < cds_start`) — e.g.
//!     `c.-2_-1dup` on `NM_212556.2` (GGG-ATG context) → `c.1delinsACA`.
//!   - Spanning Insertion / Delins outputs — e.g. `c.-1_1insCAT` on
//!     `NM_212556.2:c.2_3insCAT` → `c.1delinsATCA`. (Spanning *insertions*
//!     are not the canonical form; only spanning *duplications* are
//!     preserved.)
//!
//! Synthetic transcript here reproduces `NM_000051.3`'s c.-1..c.3 = "C |
//! A T G" context, the only difference vs the existing
//! `tests/issue_383_canon_cds_start_clamp.rs` (which uses "G | A T G"
//! UTR-flank) is `c.-1 = C`. That single base difference makes
//! `ref[c.-1] ++ ref[c.1] = "CA" == alt`, so the canon ins → dup
//! recogniser emits `c.-1_1dup` — which must survive the clamp.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Synthetic transcript whose c.-1 = C, c.1..c.3 = "ATG" — same shape as
/// `NM_000051.3` around the failing input.
///
///   axis: c.-3 c.-2 c.-1 | c.1 c.2 c.3 c.4 c.5 ...
///   base:  G    G    C   | A   T   G   C   ...
///
/// CDS spans c.1..c.200 (synthetic; CDS-end side is irrelevant here).
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let mut seq = String::from("GGC"); // 5'UTR (c.-3 = G, c.-2 = G, c.-1 = C)
    seq.push_str("ATGC"); // c.1..c.4 = A, T, G, C
    while seq.len() < 203 {
        seq.push('G');
    }
    let len = seq.len() as u64;
    // c.1 is at byte index 3 (0-based) = byte position 4 (1-based);
    // matches tests/issue_383_canon_cds_start_clamp.rs layout.
    let transcript = Transcript::new(
        "NM_TEST401.1".to_string(),
        Some("TEST401".to_string()),
        Strand::Plus,
        seq,
        Some(4),
        Some(203),
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
fn five_prime_insertion_emits_spanning_dup_when_alt_matches_boundary_bases() {
    // c.1_2insCA on a transcript whose c.-1 = C, c.1 = A: the inserted
    // "CA" equals ref[c.-1] ++ ref[c.1], so the canonical form is the
    // spanning dup `c.-1_1dup` (dup-source = the two boundary bases).
    // PR #385's clamp must NOT collapse this to `c.1delinsACA`.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST401.1:c.1_2insCA"),
        "NM_TEST401.1:c.-1_1dup",
    );
}

#[test]
fn three_prime_direction_also_emits_spanning_dup() {
    // 3'-direction sanity. On this transcript the only adjacent tandem
    // unit matching the inserted "CA" is the spanning c.-1_1 pair; the
    // 3'-direction tie-break still anchors there because there is no
    // 3'-side tract to prefer. The PR #385 clamp must not collapse this
    // either.
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, "NM_TEST401.1:c.1_2insCA"),
        "NM_TEST401.1:c.-1_1dup",
    );
}
