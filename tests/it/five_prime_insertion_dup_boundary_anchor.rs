//! Regression tests for the 5'-direction insertion→duplication anchor at the
//! transcript-start boundary (idempotency-campaign follow-up to PR #1169/#1170).
//!
//! When a 5'-shuffled insertion (or a delins that reduces to one) comes to rest
//! at the transcript start (`result.start == 0`, i.e. inserting immediately
//! before c.1) and its rotated sequence duplicates the following reference
//! tract, the emitted `dup` must be anchored on that tract. The 3'-side dup
//! anchor previously used `new_start`, which the insertion coordinate
//! adjust-back (`saturating_sub(1)`) collapses to c.1 at the boundary — so the
//! dup was mislabelled one position 3', producing a DIFFERENT haplotype
//! (`c.2_3insGA` in `G` + A-run → `c.2_3dup` (unit `AA`) instead of the correct
//! `c.1_2dup` (unit `GA`)). 3' was already correct; the fix derives the anchor
//! from the matched-tract position `pos_idx` so both directions agree.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

fn norm(core: &str, input: &str, dir: ShuffleDirection) -> String {
    Normalizer::with_config(
        SyntheticBuilder::cds(core, 1, core.len() as u64, Strand::Plus).build(),
        NormalizeConfig::default().with_direction(dir),
    )
    .normalize(&parse_hgvs(input).expect("parse"))
    .expect("normalize")
    .to_string()
}

/// `G` then a 6-A run: `c.2_3insGA` is the same haplotype as duplicating the
/// leading `GA` (`c.1_2`), not the interior `AA` (`c.2_3`). Both directions must
/// agree on `c.1_2dup`.
#[test]
fn five_prime_boundary_insertion_dups_matched_tract_not_offset() {
    assert_eq!(
        norm(
            "GAAAAAAC",
            "NM_TEST.1:c.2_3insGA",
            ShuffleDirection::FivePrime
        ),
        "NM_TEST.1:c.1_2dup",
    );
}

#[test]
fn three_prime_agrees_on_same_boundary_dup() {
    assert_eq!(
        norm(
            "GAAAAAAC",
            "NM_TEST.1:c.2_3insGA",
            ShuffleDirection::ThreePrime
        ),
        "NM_TEST.1:c.1_2dup",
    );
}

/// The resulting `c.1_2dup` must be a 5' fixed point (it previously drifted to
/// the wrong-haplotype `c.2_3dup`).
#[test]
fn five_prime_boundary_dup_is_idempotent() {
    let once = norm(
        "GAAAAAAC",
        "NM_TEST.1:c.2_3insGA",
        ShuffleDirection::FivePrime,
    );
    assert_eq!(
        norm("GAAAAAAC", &once, ShuffleDirection::FivePrime),
        once,
        "norm(norm(x)) must equal norm(x)"
    );
}

/// A delins that canonicalizes to the same boundary insertion takes the same
/// path: `c.1delinsGCG` in `GCC…` is `c.1_2dup` (dup of leading `GC`), idempotent.
#[test]
fn five_prime_boundary_delins_reduces_to_matched_tract_dup() {
    let out = norm(
        "GCCTGTGTGTGC",
        "NM_TEST.1:c.1delinsGCG",
        ShuffleDirection::FivePrime,
    );
    assert_eq!(out, "NM_TEST.1:c.1_2dup");
    assert_eq!(
        norm("GCCTGTGTGTGC", &out, ShuffleDirection::FivePrime),
        out,
        "delins-derived boundary dup must be idempotent"
    );
}
