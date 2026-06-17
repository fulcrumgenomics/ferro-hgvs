//! Issue #132: 3' rule misses cyclic-rotation single-copy insertions.
//!
//! When a literal `ins` payload is a non-zero cyclic rotation of an adjacent
//! reference repeat unit, the 3' rule should shift the variant to the most-3'
//! dup-aligned position. The 2+-copy path (`insertion_to_repeat`) already
//! handles this; this file pins the single-copy fix and the round-trip
//! equivalence between rotation-mismatched alts and the matched alt.

use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};

/// First synthetic core from the issue: `AAGTGTGTTTTTTTAA`.
/// GT[3] tract at core pos 3..8 (HGVS pos 259..264 with PAD_OFFSET=256).
/// Most-3' GT in the tract is at core 7..8 = HGVS 263..264.
const CORE_GT: &str = "AAGTGTGTTTTTTTAA";

/// Second synthetic core derived from the issue's CAG example. The issue
/// body uses `AAACAGCAGCAGCAACAGAA` (CAG[3] starting at core idx 3, two
/// bases past the insertion point at HGVS 258_259) — that layout has a
/// non-tandem 'AA' between the insertion point and the tract, so the
/// rotation helper correctly returns None and the variant stays as `ins`.
/// To exercise the cyclic-rotation path on a 3-base unit, we shift the
/// tract to start AT the insertion point: `AACAGCAGCAGCAACAGAA`.
/// CAG[3] tract at core 3..12 → HGVS 259..267 — abutting HGVS 258_259.
const CORE_CAG: &str = "AACAGCAGCAGCAACAGAA";

/// HGVS sequence-id used by `SyntheticBuilder::genomic`.
const SEQID: &str = "NC_TEST.1";

#[test]
fn ins_tg_rotates_to_most_3prime_dup() {
    // alt = TG, ref tract = GT (rotation r=1). Single copy added.
    // Expected: 3'-shift through GT[3] then dup at most-3' GT.
    let p = SyntheticBuilder::genomic(CORE_GT).build();
    let result = normalize_to_string(p, &format!("{}:g.258_259insTG", SEQID));
    assert_eq!(result, format!("{}:g.263_264dup", SEQID));
}

#[test]
fn ins_gt_matched_rotation_already_works() {
    // alt = GT, ref tract = GT (rotation r=0). Already shifts via shuffle.
    // This is the "already works" baseline from the issue table — keep it
    // as a regression lock so any change to the rotation logic does not
    // silently flip this.
    let p = SyntheticBuilder::genomic(CORE_GT).build();
    let result = normalize_to_string(p, &format!("{}:g.258_259insGT", SEQID));
    assert_eq!(result, format!("{}:g.263_264dup", SEQID));
}

#[test]
fn ins_tgtg_two_copies_already_works() {
    // alt = TGTG, ref tract = GT[3]. Multi-copy path
    // (`insertion_to_repeat`) already handles this; lock the behavior.
    let p = SyntheticBuilder::genomic(CORE_GT).build();
    let result = normalize_to_string(p, &format!("{}:g.258_259insTGTG", SEQID));
    assert_eq!(result, format!("{}:g.259_264GT[5]", SEQID));
}

#[test]
fn ins_idempotent_post_fix() {
    // Round-trip: normalize(normalize(insTG)) == normalize(insTG).
    let p = SyntheticBuilder::genomic(CORE_GT).build();
    let once = normalize_to_string(p, &format!("{}:g.258_259insTG", SEQID));
    let p2 = SyntheticBuilder::genomic(CORE_GT).build();
    let twice = normalize_to_string(p2, &once);
    assert_eq!(once, twice);
    assert_eq!(once, format!("{}:g.263_264dup", SEQID));
}

#[test]
fn ins_agc_cag_tract_rotates_to_dup() {
    // CAG[3] tract starting at HGVS 259 (core idx 2). Insert AGC at
    // HGVS 258_259 (just before tract's first C). alt "AGC" is the r=1
    // rotation of canonical "CAG".
    //
    // Rotation iteration on `smallest_repeat_unit("AGC") = "AGC"` →
    // "AGC" (r=0): no anchor matches the abutting window;
    // "GCA" (r=1): no anchor matches;
    // "CAG" (r=2): anchor at ins_point=258 matches; tract walks right
    //   through CAGCAGCAG, ref_count=3.
    // Only CAG-rotation matches. Most-3' CAG dup → tract end=267,
    // dup_start_idx=264, dup_end_idx=266 → HGVS [265..267].
    let p = SyntheticBuilder::genomic(CORE_CAG).build();
    let result = normalize_to_string(p, &format!("{}:g.258_259insAGC", SEQID));
    assert_eq!(result, format!("{}:g.265_267dup", SEQID));
}

#[test]
fn ins_gca_cag_tract_rotates_to_dup() {
    // Same layout as above, alt = GCA (r=2 of canonical "GCA"). Rotation
    // iteration on `smallest_repeat_unit("GCA") = "GCA"`: only the "CAG"
    // rotation finds an anchored tract abutting the ins point. Result is
    // identical to the AGC case — both spec-equivalent rotations of the
    // same tandem unit produce the same canonical dup.
    let p = SyntheticBuilder::genomic(CORE_CAG).build();
    let result = normalize_to_string(p, &format!("{}:g.258_259insGCA", SEQID));
    assert_eq!(result, format!("{}:g.265_267dup", SEQID));
}
