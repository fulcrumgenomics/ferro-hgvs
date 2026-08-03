//! Moving a junction-occupying member must not rewrite the bases it adds.
//!
//! A member that adds sequence at a junction denotes its payload *at that
//! junction*. Two of the sibling-repair passes move such a member, and a
//! `Duplication` reads its payload from the reference **under its own span** —
//! so sliding the span silently changes what the member copies:
//!
//! ```text
//! reference  ("ACGT" x 64) + "GCAAATAACATCCAGA" + ("ACGT" x 64)
//! g.[263_264insAC;265del]
//!   per-member ->  g.[265_266dup;265del]      dup payload "CA"
//!   clamped    ->  g.[263_264dup;265del]      payload silently became "AA"
//!   merged     ->  g.265delinsAA              the inserted C is gone
//!
//!   intended  T A A C A A T
//!   emitted   T A A A A A T
//! ```
//!
//! The repair is to rotate the payload into phase at the destination:
//!
//! ```text
//! insert P at junction j  ==  insert (ref[j] ++ P[..n-1]) at junction j - 1
//!                             only when P[n-1] == ref[j]
//! ```
//!
//! and to re-spell the member as a plain insertion carrying it whenever the
//! translated spelling would denote something else. #1292 is that defect at
//! `clamp_sibling_crossing_junctions`; #1280 is the same root cause reached
//! through `demote_repeats_spanning_siblings`, whose duplication is then pulled
//! back by the very same clamp.
//!
//! The rotation is what makes this safe in a tandem repeat. #1280 records an
//! attempted fix — demote the repeat to an insertion carrying its **literal**
//! payload — that was reverted precisely because it skipped the rotation:
//! invisible in a homopolymer, wrong in a `(CAG)` tract. The last test here is
//! that counterexample.

use crate::common::synthetic::assert_padded_preserving;

#[test]
fn a_clamped_duplication_keeps_the_bases_it_copies() {
    // #1292. The inserted `AC` cancels against the deleted `C`, leaving one
    // added `A`; the payload must survive being pulled back two positions.
    let output = assert_padded_preserving("GCAAATAACATCCAGA", "NC_TEST.1:g.[263_264insAC;265del]");
    assert_eq!(output, "NC_TEST.1:g.266dup");
}

#[test]
fn a_demoted_repeat_keeps_its_payload_when_the_clamp_pulls_it_back() {
    // #1280 (1). An insertion grows a nine-`T` tract; the repeat is demoted to
    // a duplication over the tract's 3'-most bases, which reaches across the
    // substituted base, and the junction clamp then pulls it back.
    let output = assert_padded_preserving(
        "CTTTTTTTTTAATATATTTTG",
        "NC_TEST.1:g.[259_260insTTTTT;262T>A]",
    );
    assert_eq!(output, "NC_TEST.1:g.263_267dup");
}

#[test]
fn the_pulled_back_payload_does_not_pick_up_a_flanking_base() {
    // #1280 (2). The same shape with a `GC` run 5' of the tract, which is what
    // made the corruption visible: the translated duplication read a `C` from
    // outside the tract and emitted `g.272delinsCTTTTA` — a base that appears
    // nowhere in the sequence the input denotes.
    let output = assert_padded_preserving(
        "CGCGCGCGCGCTTTTTTTTTAATATATTTTG",
        "NC_TEST.1:g.[269_270insTTTTT;272T>A]",
    );
    assert_eq!(output, "NC_TEST.1:g.273_277dup");
}

#[test]
fn a_tandem_repeat_insertion_beside_a_deletion_merges_to_one_insertion() {
    // #1280's counterexample, and the reason the payload is *rotated* rather
    // than carried literally. A duplication is phase-correct by construction
    // where it already sits, so a member that does not have to move keeps it.
    //
    // Since #1235 the pair does not survive to be spelled per member: the merged
    // form is derived from the sequence, and a `-1` deletion with a `+6`
    // insertion six bases away has a five-column single-gap explanation against
    // the input's seven, so the derivation takes it. The payload lands 5' of the
    // `CAG` tract with a `C` to its 3', so it cannot rotate onto the tract and is
    // not a tandem duplication of what precedes it — the `dup` mandate is not in
    // play on the merged form.
    let output = assert_padded_preserving(
        "CTTTTCAGCAGCAGCAGCAGCAGTTTTG",
        "NC_TEST.1:g.[268_269insAGCAGC;262del]",
    );
    assert_eq!(output, "NC_TEST.1:g.261_262insAGCAG");
}
