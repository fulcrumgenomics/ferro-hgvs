//! The examined-`=` provenance channel: annotation-only re-attach (design §3 ③).
//!
//! After a strategy has produced its definite partition, the caller's asserted
//! examined regions ([`ExaminedRegion`]) are re-attached as identity members —
//! the `=` a sequence-first derivation drops, because an unchanged stretch
//! between two edits emits no member of its own. This is the bakeoff analog of
//! the production plan's `reattach_examined` (2026-08-19 provenance-channel-plan
//! §step 3): pure over `(partition, examined, reference, resulting)`, it NEVER
//! reads an input spelling (confluence-by-construction), and it is annotation-only
//! — an identity member changes no base, so round-trip is preserved and the
//! fragmentation metrics (`member_count` / `rendered_cost`, in the harness's
//! `bakeoff::metrics`) that filter out
//! `Identity` are unmoved.
//!
//! A region is honoured only if it is **validated against the sequence**: in
//! range, and touched by no edit member (an edit over it would contradict
//! "unchanged"). Because [`RefApplier`](crate::partition::metrics::RefApplier)
//! copies every untouched reference gap verbatim, an edit-free region is
//! unchanged in the output by construction; the explicit byte comparison against
//! `resulting` makes that independent of whether the partition round-trips, so an
//! examined-but-changed region is dropped rather than re-attached.
//!
//! ∅ (empty) examined ⟹ the pass is a no-op, so a corpus that supplies no
//! provenance is byte-identical with and without it.

// The re-attach helpers below are consumed only by the dev-gated arm layer / bakeoff
// analog; the production ruled core uses a subset. In a non-`dev` build the rest is
// dead, which the required `cargo clippy --release -- -D warnings` job (ci.yml) rejects
// — allow it only there. A genuinely-dead item is still caught by the `dev` clippy pass.
#![cfg_attr(not(feature = "dev"), allow(dead_code))]

use crate::partition::block_ctx::ExaminedRegion;
use crate::partition::output::{EditKind, Member, Partition};

/// The number of bases a member contributes to the resulting sequence. Mirrors
/// [`RefApplier`](crate::partition::metrics::RefApplier)'s match arms exactly:
/// `Inv`/`Identity` are length-preserving with an EMPTY `inserted` (the payload
/// is implied), so their output length is the span, not `inserted.len()` — the
/// trap a naive `inserted.len()` projection falls into.
fn member_output_len(m: &Member) -> usize {
    match m.kind {
        EditKind::Inv | EditKind::Identity => m.ref_end - m.ref_start,
        EditKind::Del => 0,
        EditKind::Ins | EditKind::Sub | EditKind::Delins | EditKind::Dup => m.inserted.len(),
    }
}

/// Whether an examined region is honourable: in range, non-empty, touched by no
/// member, and byte-for-byte unchanged in `resulting`. The resulting offset is
/// the reference offset shifted by the net length delta of every member lying
/// entirely 5' of the region (members are, by the no-overlap check, either
/// wholly before or wholly after — never within).
fn region_is_confirmed_unchanged(
    region: &ExaminedRegion,
    members: &[Member],
    reference: &[u8],
    resulting: &[u8],
) -> bool {
    if region.start >= region.end || region.end > reference.len() {
        return false;
    }
    // An edit (or an already-attached identity) over the region contradicts the
    // "unchanged" claim, or would double-annotate it.
    let overlaps = members
        .iter()
        .any(|m| m.ref_start < region.end && region.start < m.ref_end);
    if overlaps {
        return false;
    }
    // Project the reference offset onto the resulting axis and compare bytes.
    let delta: isize = members
        .iter()
        .filter(|m| m.ref_end <= region.start)
        .map(|m| member_output_len(m) as isize - (m.ref_end - m.ref_start) as isize)
        .sum();
    let res_start = region.start as isize + delta;
    let len = region.end - region.start;
    if res_start < 0 || res_start as usize + len > resulting.len() {
        return false;
    }
    let res_start = res_start as usize;
    reference[region.start..region.end] == resulting[res_start..res_start + len]
}

/// Re-attach the caller's validated examined-`=` regions as identity members
/// (design §3 ③). Annotation-only and pure over its inputs; ∅ examined is a
/// no-op. See the module docs.
pub(crate) fn reattach_examined(
    partition: Partition,
    examined: &[ExaminedRegion],
    reference: &[u8],
    resulting: &[u8],
) -> Partition {
    if examined.is_empty() {
        return partition; // ∅ ⟹ no-op: the default corpus is untouched.
    }
    let mut members = partition.members;
    for region in examined {
        if region_is_confirmed_unchanged(region, &members, reference, resulting) {
            members.push(Member {
                kind: EditKind::Identity,
                ref_start: region.start,
                ref_end: region.end,
                inserted: Vec::new(),
            });
        }
    }
    // Canonical order (RefApplier's own key): a zero-width edit sorts before a
    // span at a shared start. Deterministic regardless of the caller's region
    // order, so confluence-within-class is preserved.
    members.sort_by_key(|m| (m.ref_start, m.ref_end));
    Partition { members }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn member(kind: EditKind, ref_start: usize, ref_end: usize, ins: &str) -> Member {
        Member {
            kind,
            ref_start,
            ref_end,
            inserted: ins.as_bytes().to_vec(),
        }
    }

    fn identity(ref_start: usize, ref_end: usize) -> Member {
        member(EditKind::Identity, ref_start, ref_end, "")
    }

    #[test]
    fn empty_examined_is_a_no_op() {
        // reference AAAAA -> resulting AACAA (a Sub at offset 2)
        let p = Partition {
            members: vec![member(EditKind::Sub, 2, 3, "C")],
        };
        let got = reattach_examined(p.clone(), &[], b"AAAAA", b"AACAA");
        assert_eq!(
            got, p,
            "no examined regions ⟹ the partition is returned as-is"
        );
    }

    #[test]
    fn an_unchanged_examined_region_beside_an_edit_is_reattached() {
        // reference AAAAA -> resulting AACAA: Sub at [2,3), region [3,5) unchanged.
        let p = Partition {
            members: vec![member(EditKind::Sub, 2, 3, "C")],
        };
        let region = ExaminedRegion { start: 3, end: 5 };
        let got = reattach_examined(p, &[region], b"AAAAA", b"AACAA");
        assert_eq!(
            got.members,
            vec![member(EditKind::Sub, 2, 3, "C"), identity(3, 5)],
            "the edit-free, unchanged region is re-attached as an identity member"
        );
    }

    #[test]
    fn a_region_an_edit_touches_is_dropped() {
        // The Sub at [2,3) overlaps the claimed region [1,4): contradicts "unchanged".
        let p = Partition {
            members: vec![member(EditKind::Sub, 2, 3, "C")],
        };
        let region = ExaminedRegion { start: 1, end: 4 };
        let got = reattach_examined(p, &[region], b"AAAAA", b"AACAA");
        assert_eq!(
            got.members,
            vec![member(EditKind::Sub, 2, 3, "C")],
            "a region an edit touches is not re-attached"
        );
    }

    #[test]
    fn an_out_of_range_region_is_dropped() {
        let p = Partition {
            members: vec![member(EditKind::Sub, 2, 3, "C")],
        };
        let region = ExaminedRegion { start: 4, end: 9 }; // end past the reference
        let got = reattach_examined(p, &[region], b"AAAAA", b"AACAA");
        assert_eq!(got.members, vec![member(EditKind::Sub, 2, 3, "C")]);
    }

    #[test]
    fn the_resulting_offset_is_projected_across_an_indel() {
        // reference AAAAAA -> resulting AAAA (a Del of [1,3)). Region [3,6) is
        // unchanged; on the resulting axis it lands at [1,4), which the byte check
        // must find by shifting -2 for the deletion — not at reference offset 3.
        let p = Partition {
            members: vec![member(EditKind::Del, 1, 3, "")],
        };
        let region = ExaminedRegion { start: 3, end: 6 };
        let got = reattach_examined(p, &[region], b"AAAAAA", b"AAAA");
        assert_eq!(
            got.members,
            vec![member(EditKind::Del, 1, 3, ""), identity(3, 6)],
            "the region projects across the deletion and is confirmed unchanged"
        );
    }

    #[test]
    fn the_projection_survives_a_length_preserving_inversion() {
        // reference ACGTAA -> resulting ACGATA is NOT what an inv of [3,5) gives;
        // build the resulting from the members so the byte check is honest.
        // Inv of [1,4) (CGT -> revcomp ACG) then an unchanged tail [4,6).
        // reference: A C G T A A ; inv[1,4) => A + revcomp(CGT)=ACG + T A A
        //                                     = A A C G T A A ? recompute below.
        // Use a clean case: reference ACGTT -> inv[0,3) (ACG->CGT) => CGTTT tail T T.
        let reference = b"ACGTT";
        // inv of [0,3): revcomp(ACG) = CGT ; tail [3,5) = TT unchanged.
        let resulting = b"CGTTT";
        let p = Partition {
            members: vec![member(EditKind::Inv, 0, 3, "")],
        };
        let region = ExaminedRegion { start: 3, end: 5 }; // the unchanged TT tail
        let got = reattach_examined(p, &[region], reference, resulting);
        assert_eq!(
            got.members,
            vec![member(EditKind::Inv, 0, 3, ""), identity(3, 5)],
            "an Inv is length-preserving, so the tail projects with zero shift"
        );
    }

    #[test]
    fn regions_are_reattached_in_canonical_order_regardless_of_input_order() {
        // Two unchanged regions [0,1) and [4,5) around a Sub at [2,3), supplied
        // out of order — the result must be sorted by (ref_start, ref_end).
        let p = Partition {
            members: vec![member(EditKind::Sub, 2, 3, "C")],
        };
        let regions = [
            ExaminedRegion { start: 4, end: 5 },
            ExaminedRegion { start: 0, end: 1 },
        ];
        let got = reattach_examined(p, &regions, b"AAAAA", b"AACAA");
        assert_eq!(
            got.members,
            vec![
                identity(0, 1),
                member(EditKind::Sub, 2, 3, "C"),
                identity(4, 5),
            ],
            "members come back in canonical (ref_start, ref_end) order"
        );
    }
}
