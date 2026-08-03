//! Commuting is tested against the payload a member would carry **at** the
//! junction it lands on, not the one it carries now.
//!
//! `clamp_sibling_crossing_junctions` lets a member meet a sibling whose payload
//! commutes with its own, because sharing a junction is well-defined there and
//! the pair then merges (#1286, #1308). Landing rotates the payload into phase —
//! and a rotation is exactly what can destroy the commuting identity:
//!
//! ```text
//! reference  ("ACGT" x 64) + "TAAAACCA" + ("ACGT" x 64)
//! g.[260_261insAC;261_262insAC]
//!   `AC` and `AC` commute, so the first member was allowed onto junction 261
//!   at 261 its payload is `CA`, and `CA ++ AC != AC ++ CA`
//!
//!   intended  T A A A A C A A C C C A
//!   emitted   T A A A A A C C A C C A
//! ```
//!
//! Two members then share a junction with no defined order between them, and the
//! merge downstream concatenates them in *rendered* order — producing `ACCA`,
//! which is neither payload nor any rotation of the pair.
//!
//! Testing the landed payload keeps them apart: the bound falls back to one
//! position short, and the allele renders as its input does.

use crate::common::synthetic::assert_padded_preserving;

#[test]
fn a_rotation_that_breaks_commuting_keeps_the_pair_apart() {
    // #1312. `AC` at 260 becomes `CA` at 261, which does not commute with the
    // `AC` already there, so the member may not land on that junction.
    //
    // The clamp still refuses that landing — this file's property is unchanged —
    // but since #1235 the allele's output is derived from the sequence the pair
    // denotes rather than assembled from the clamped members, and that sequence
    // is one four-base insertion. It is the `#1312` row of
    // `cis_spelling_confluence_gap.rs`, now converged. The `ACCA` corruption the
    // module doc describes is still excluded: `assert_padded_preserving` applies
    // both spellings with an applier that is not the normalizer.
    let output = assert_padded_preserving("TAAAACCA", "NC_TEST.1:g.[260_261insAC;261_262insAC]");
    assert_eq!(output, "NC_TEST.1:g.262_263insAACC");
}

#[test]
fn a_rotation_that_preserves_commuting_still_lets_them_merge() {
    // The guard. In a poly-A tract every rotation of `A` is `A`, so the landed
    // payload still commutes and the pair merges into one member.
    let output = assert_padded_preserving("AAAAAA", "NC_TEST.1:g.[258_259insA;259_260insA]");
    assert_eq!(output, "NC_TEST.1:g.257_263A[9]");
}
