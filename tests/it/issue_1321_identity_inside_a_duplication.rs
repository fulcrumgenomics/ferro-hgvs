//! A duplication covers the positions it copies, so an identity member inside
//! its span is a residue to drop.
//!
//! `drop_identity_members_covered_by_siblings` removes a member that cancelled
//! to `=` while a sibling grew over the bases it names (#1297). It gates on
//! `claims_reference_bases`, which deliberately excludes `Duplication` — a dup
//! adds a copy at a junction and consumes nothing — so a duplication that grew
//! across the cancelled member's position did not count:
//!
//! ```text
//! reference  ("ACGT" x 64) + "TCCCAGAAAAT" + ("ACGT" x 64)
//!
//! g.[261_262insGA;262_263insA;263del]  ->  g.[262_263dup;263=]
//!   the dup spans 262-263 and the identity member names 263, inside it
//! ```
//!
//! `g.262_263dup` on its own already denotes the input — it copies `G` at 262
//! and `A` at 263, giving `CAGA` + `GA` + `AAAT`, which is exactly what the
//! three input members describe. The `263=` adds nothing and overlaps, so the
//! pair violates #1235's second criterion.
//!
//! This is #1297's shape one spelling over, and the gate widens from
//! `claims_reference_bases` to `blocks_sibling_shift` — the predicate that is
//! already "claims bases, plus a duplication", and which exists because a
//! duplication reads its payload from the reference under its own span. A
//! duplication therefore *names* the positions it copies even though it consumes
//! none of them, which is exactly the property that makes an identity member
//! inside it redundant (#1321).
//!
//! Insertions stay excluded, and must: an insertion's span is the gap it sits
//! in, so `g.264_265insA` nominally covers position 264 while contradicting
//! nothing there. `blocks_sibling_shift` excludes them for the same reason
//! `claims_reference_bases` did.

use crate::common::synthetic::assert_padded_preserving;

/// The core the proptest shrank to, and the reference the seed describes.
const CORE: &str = "TCCCAGAAAAT";

#[test]
fn an_identity_member_inside_a_duplication_is_dropped() {
    // #1321. The seed's own case. `262_263dup` alone denotes the input; the
    // `263=` it was left beside adds nothing and overlaps it.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[261_262insGA;262_263insA;263del]");
    assert_eq!(output, "NC_TEST.1:g.262_263dup");
}

#[test]
fn an_identity_member_beside_an_unrelated_deletion_survives() {
    // The guard #1297 installed and this must not weaken. An identity member is
    // real information when nothing covers it: `g.[1002=;1005del]` keeps both,
    // and dropping identities outright would lose what the `=` records.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[258=;260del]");
    assert!(
        output.contains("258="),
        "an uncovered identity member must survive: `{output}`"
    );
}

#[test]
fn an_identity_member_beside_an_insertion_survives() {
    // The reason the gate cannot simply be "a sibling's span covers it". An
    // insertion's span is the gap it occupies, so `259_260insG` nominally spans
    // 259-260 while changing nothing at 259 — it contradicts an identity there,
    // and `blocks_sibling_shift` excludes it exactly as `claims_reference_bases`
    // did.
    // The identity sits *inside* the insertion's nominal span, which is the
    // case the paragraph above is about — an identity outside it was never at
    // risk from this gate and so proves nothing.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[259=;259_260insG]");
    assert!(
        output.contains("259="),
        "an insertion must not drop an identity member inside its span: `{output}`"
    );
}
