//! Regression tests for issue #418 subgroup (a): NM_212556.2 panic +
//! canonicalisation misses at the c.0/c.1 boundary on transcripts
//! without a 5'UTR (`cds_start == 1`).
//!
//! Setup mirrors issue #383's synthetic transcript geometry but with
//! `cds_start = 1` (no 5'UTR), reproducing the real NM_212556.2 layout
//! whose c.1..c.3 = "ATG…" surfaced post-PR #417 when the actual `.2`
//! bases reached `axis_normalized` (pre-#417, the version-strip fallback
//! resolved `.2` against `.4` bases, masking the bug).
//!
//! Three failure modes covered:
//!
//!   NM_212556.2:c.1delinsCA    panic       -> c.1delinsCA   (input verbatim)
//!   NM_212556.2:c.1_2insCA     c.1insCA    -> c.1delinsACA  (boundary absorb)
//!   NM_212556.2:c.2_3insCAT    c.1insCAT   -> c.1delinsATCA (boundary absorb)
//!
//! Root causes:
//!
//! * The panic surfaces in `hgvs_pos_to_index(0)` when
//!   `canonicalize_delins`'s shared-affix trim collapses `c.1delinsCA` to
//!   `Insertion("C")` at `after_index = 0` (i.e. an insertion before the
//!   first transcript byte). The recursive `normalize_na_edit` call then
//!   passes `tx_start = 0` — not a valid 1-based HGVS position — and the
//!   `pos - 1` conversion underflows. Fixed by suppressing the recursion
//!   when `after_index = 0` and restoring the input Delins form (the same
//!   disposition the #383 post-shift clamp already gives Delins inputs
//!   whose canonicalisation would cross the CDS-start boundary).
//!
//! * The canon mismatch for `c.1_2insCA` / `c.2_3insCAT` (5'+cross) is
//!   the same boundary issue from the Insertion-input side: the 5'-shift
//!   saturates at the transcript start (`result.start = 0` in 0-based
//!   coords), but the 1-based HGVS conversion clamps both `new_tx_start`
//!   and `new_tx_end` back to `cds_start = 1`, producing a degenerate
//!   `c.1ins<seq>` shape that bypasses the #383 clamp's
//!   `new_tx_start < cds_start` gate. Fixed by extending the gate to
//!   detect the degenerate-Insertion saturation case and apply the same
//!   clamp formula.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Synthetic transcript whose `cds_start = 1` (no 5'UTR) and whose
/// first three CDS bases are `c.1..c.3 = "ATG"` — the same start-codon
/// context as NM_212556.2.
///
///   axis:  c.1 c.2 c.3 c.4 c.5 c.6 ...
///   base:   A   T   G   C   G   G  ...
///
/// CDS spans c.1..c.1000 (synthetic); 3'UTR is empty.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let mut seq = String::from("ATGC"); // c.1..c.4
    while seq.len() < 1500 {
        seq.push('G');
    }
    let len = seq.len() as u64;
    let transcript = Transcript::new(
        "NM_TEST418.1".to_string(),
        Some("TEST418".to_string()),
        Strand::Plus,
        seq,
        Some(1),    // cds_start: c.1 at tx byte 1 (1-based) = index 0 (no 5'UTR)
        Some(1003), // cds_end (synthetic)
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
    format!("{}", normalized)
}

#[test]
fn delins_at_c1_no_utr5_three_prime_no_panic() {
    // Pre-fix this panicked with `attempt to subtract with overflow` at
    // `hgvs_pos_to_index(0)` because the shared-affix trim of
    // `c.1delinsCA` (delete c.1=A, insert "CA") collapses to an
    // `Insertion("C")` at `after_index = 0` and the recursive
    // `normalize_na_edit` then calls `hgvs_pos_to_index(0)` which
    // underflows. Spec-canonical: restore input verbatim.
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, "NM_TEST418.1:c.1delinsCA"),
        "NM_TEST418.1:c.1delinsCA",
    );
}

#[test]
fn delins_at_c1_no_utr5_five_prime_no_panic() {
    // Same shape, 5'-direction. Both directions must restore the input
    // since `canonicalize_delins` is direction-agnostic and produces the
    // same panicking `Insertion(after_index=0)` shape regardless of
    // direction.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST418.1:c.1delinsCA"),
        "NM_TEST418.1:c.1delinsCA",
    );
}

#[test]
fn insertion_at_c1_no_utr5_five_prime_clamps_to_delins_at_c1() {
    // `c.1_2insCA` 5'+cross on a no-5'UTR transcript whose c.1=A: the
    // 5'-shuffle walks left by one position (alt rotates "CA" → "AC"),
    // saturating at the transcript start (`result.start = 0`). The
    // coordinate conversion clamps both `new_tx_start` and `new_tx_end`
    // back to `cds_start = 1`, producing a degenerate `c.1insCA` shape
    // that the original `new_tx_start < cds_start` gate missed. Spec-
    // canonical clamp absorbs ref[c.1] = "A" into the alt: `c.1delinsACA`.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST418.1:c.1_2insCA"),
        "NM_TEST418.1:c.1delinsACA",
    );
}

#[test]
fn insertion_two_bases_in_no_utr5_five_prime_clamps_to_delins_at_c1() {
    // `c.2_3insCAT` 5'+cross: 5'-shift rotates "CAT" twice walking left
    // through c.2=T and c.1=A, saturating at the transcript start. Same
    // saturation as above. Spec-canonical clamp emits
    // `c.1delins<ref[c.1..c.3]>+<alt[..2]> = ATCA`.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST418.1:c.2_3insCAT"),
        "NM_TEST418.1:c.1delinsATCA",
    );
}

#[test]
fn insertion_at_c1_no_utr5_three_prime_unaffected() {
    // Negative: 3'-direction shuffle on the same input must NOT clamp.
    // ref[c.3] = G != alt rotation, so 3'-shift doesn't move and the
    // input is preserved. Pins the directional-specificity of the
    // saturation gate so a stray clamp on the 3'-direction would be
    // caught.
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, "NM_TEST418.1:c.1_2insCA"),
        "NM_TEST418.1:c.1_2insCA",
    );
}

#[test]
fn insertion_two_bases_in_no_utr5_three_prime_unaffected() {
    // Negative companion to the 5'-direction case for `c.2_3insCAT`.
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, "NM_TEST418.1:c.2_3insCAT"),
        "NM_TEST418.1:c.2_3insCAT",
    );
}

#[test]
fn cds_interior_insertion_no_utr5_unaffected() {
    // Negative: a CDS-interior insertion far from the c.1 boundary must
    // not be touched by the saturation gate. The synthetic transcript's
    // c.10..c.11 are both G; inserting "C" between them doesn't shift,
    // and `new_tx_start = 10 != cds_start = 1`, so the gate stays off.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST418.1:c.10_11insC"),
        "NM_TEST418.1:c.10_11insC",
    );
}

#[test]
fn delins_cds_interior_no_utr5_unaffected() {
    // Negative: a CDS-interior delins not at the c.1 boundary must not
    // trigger the post-shift Delins-arm clamp. Pin the result so the
    // suppression branch staying off on interior inputs is regressed
    // here, not just in `issue_383_canon_cds_start_clamp`.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST418.1:c.10delinsCA"),
        "NM_TEST418.1:c.10delinsCA",
    );
}
