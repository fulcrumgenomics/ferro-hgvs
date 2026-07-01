//! Issue #132 / #882: cyclic-rotation single-copy insertions and the in-phase
//! gate on dup/repeat conversion.
//!
//! Corrected understanding (#882): a cyclic rotation of the reference repeat
//! unit only 3'-shifts to a `dup`/repeat when the rotation is IN PHASE with the
//! tract at the insertion cut — i.e. when the emitted tandem edit decodes to
//! exactly the same sequence as the input insertion. An IN-PHASE rotation
//! (e.g. `ins_gt_matched_rotation_already_works`) becomes a dup. An OUT-OF-PHASE
//! rotation would decode to a DIFFERENT sequence, so per HGVS duplication.md:19
//! ("no evidence the extra copy is in tandem … it should be described as an
//! insertion") the #882 phase gate rejects the conversion and the variant stays
//! a plain `ins`. The synthetic cores below (`CORE_GT`, `CORE_CAG`) are phased
//! so the tested rotations are OUT of phase, so they stay `ins` — matching live
//! mutalyzer 3.1.1.

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
fn ins_tg_out_of_phase_stays_ins() {
    // alt = TG, ref tract = GT (rotation r=1) — OUT OF PHASE at the 258_259
    // cut. A candidate dup at the most-3' GT would decode to a different
    // sequence than inserting "TG" here, so the #882 phase gate rejects it and
    // the variant stays a plain `ins` (duplication.md:19). (Mutalyzer keeps
    // out-of-phase insertions as ins; verified on real loci — this synthetic
    // NC_TEST.1 accession isn't resolvable by mutalyzer.)
    let p = SyntheticBuilder::genomic(CORE_GT).build();
    let result = normalize_to_string(p, &format!("{}:g.258_259insTG", SEQID));
    assert_eq!(result, format!("{}:g.258_259insTG", SEQID));
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
fn ins_tgtg_out_of_phase_stays_ins() {
    // alt = TGTG, ref tract = GT[3] — OUT OF PHASE at the 258_259 cut. The
    // candidate GT[5] repeat would decode to a different sequence than
    // inserting "TGTG" here, so the #882 phase gate rejects it and the variant
    // stays a plain `ins` (duplication.md:19).
    let p = SyntheticBuilder::genomic(CORE_GT).build();
    let result = normalize_to_string(p, &format!("{}:g.258_259insTGTG", SEQID));
    assert_eq!(result, format!("{}:g.258_259insTGTG", SEQID));
}

#[test]
fn ins_idempotent_post_fix() {
    // Round-trip: normalize(normalize(insTG)) == normalize(insTG). The
    // out-of-phase insTG stays a plain `ins` (#882) and is idempotent.
    let p = SyntheticBuilder::genomic(CORE_GT).build();
    let once = normalize_to_string(p, &format!("{}:g.258_259insTG", SEQID));
    let p2 = SyntheticBuilder::genomic(CORE_GT).build();
    let twice = normalize_to_string(p2, &once);
    assert_eq!(once, twice);
    assert_eq!(once, format!("{}:g.258_259insTG", SEQID));
}

#[test]
fn ins_agc_out_of_phase_stays_ins() {
    // CAG[3] tract starting at HGVS 259 (core idx 2). Insert AGC at
    // HGVS 258_259 (just before tract's first C). alt "AGC" is the r=1
    // rotation of canonical "CAG" — OUT OF PHASE at this cut. A candidate dup
    // anywhere in the tract would decode to a different sequence than inserting
    // "AGC" here, so the #882 phase gate rejects the conversion and the variant
    // stays a plain `ins` (duplication.md:19). (Mutalyzer keeps out-of-phase
    // insertions as ins; verified on real loci — NC_TEST.1 is synthetic.)
    let p = SyntheticBuilder::genomic(CORE_CAG).build();
    let result = normalize_to_string(p, &format!("{}:g.258_259insAGC", SEQID));
    assert_eq!(result, format!("{}:g.258_259insAGC", SEQID));
}

#[test]
fn ins_gca_out_of_phase_stays_ins() {
    // Same layout/flank as `ins_agc_out_of_phase_stays_ins`, alt = GCA (a
    // different out-of-phase rotation of the same tandem unit). It likewise
    // stays a plain `ins` at the input position per the #882 phase gate.
    let p = SyntheticBuilder::genomic(CORE_CAG).build();
    let result = normalize_to_string(p, &format!("{}:g.258_259insGCA", SEQID));
    assert_eq!(result, format!("{}:g.258_259insGCA", SEQID));
}
