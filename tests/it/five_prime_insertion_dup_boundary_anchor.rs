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

fn norm_genomic(core: &str, input: &str, dir: ShuffleDirection) -> String {
    Normalizer::with_config(
        SyntheticBuilder::genomic(core).build(),
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

/// 5' single-copy dup into a dinucleotide tract whose run continues one base
/// further 5' in the *opposite* phase than the abutting-unit anchor. Inserting
/// `TG` at the 3' edge of an odd-length `GTGTGTGTG` run: `find_tandem_extent`
/// anchors on the abutting `TG` phase (one base short of the run's true 5' end),
/// so the naive `ref_start` slot emits `g.259_260dup` — but the run's 5'-most
/// phase is `GT` one base further 5'. Without the mirror `five_prime_align_tract`
/// (the sibling of the 3' branch's `three_prime_align_tract`), the emitted dup
/// under-shifts by one and re-shuffles on the second pass (`g.259_260dup` →
/// `g.258_259dup`), breaking idempotency. The fix pushes to the true 5'-most
/// slot: `g.258_259dup`, a fixed point. Regression for the ins→dup sibling of #8.
#[test]
fn five_prime_dinucleotide_dup_crosses_off_phase_extension_and_is_idempotent() {
    // core: A | GTGTGTGTG (odd 9-base GT run, g.258..266) | C
    let core = "AGTGTGTGTGC";
    let once = norm_genomic(
        core,
        "NC_TEST.1:g.266_267insTG",
        ShuffleDirection::FivePrime,
    );
    assert_eq!(
        once, "NC_TEST.1:g.258_259dup",
        "must reach the true 5'-most GT slot"
    );
    assert_eq!(
        norm_genomic(core, &once, ShuffleDirection::FivePrime),
        once,
        "5' ins→dup over a dinucleotide tract must be a fixed point"
    );
}

/// 3' is unaffected: it phases to the 3'-most slot and is idempotent.
#[test]
fn three_prime_dinucleotide_dup_unchanged_and_idempotent() {
    let core = "AGTGTGTGTGC";
    let once = norm_genomic(
        core,
        "NC_TEST.1:g.266_267insTG",
        ShuffleDirection::ThreePrime,
    );
    assert_eq!(once, "NC_TEST.1:g.265_266dup");
    assert_eq!(
        norm_genomic(core, &once, ShuffleDirection::ThreePrime),
        once,
        "3' must remain idempotent"
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
