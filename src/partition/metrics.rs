//! Soundness primitives: round-trip and renderability — the two hard gates
//! (design §6.3). Pure over `(reference, resulting)` / a partition, so they live
//! in `partition` and are shared by both the harness and any production caller.
//! The PROFILE / scoring metrics (agreement, member counts, idempotence) stay in
//! `bakeoff::metrics`, since they are a measurement concern.

use crate::partition::block_ctx::BlockCtx;
use crate::partition::output::{EditKind, Member, Partition};
use std::borrow::Cow;

/// Reconstruct the resulting sequence by applying a partition to a reference.
/// Injected so metrics are testable without ferro's production renderer.
pub trait SequenceApplier {
    fn apply(&self, reference: &[u8], p: &Partition) -> Option<Vec<u8>>;

    /// Whether applying `p` to `reference` reconstructs `resulting` exactly — the
    /// round-trip soundness question ([`round_trips`]), factored so an implementor
    /// can answer it without materialising the spliced window. The default routes
    /// through [`apply`](Self::apply) and is the oracle every specialised path is
    /// diffed against — [`RefApplier`]'s allocation-free override is pinned to it
    /// by `reconstructs_streaming_matches_apply`.
    fn reconstructs(&self, reference: &[u8], resulting: &[u8], p: &Partition) -> bool {
        self.apply(reference, p).as_deref() == Some(resulting)
    }
}

/// Complement over the full IUPAC alphabet, `None` on a byte it does not model.
/// Routed through [`crate::sequence::complement_base`] so an `Inv` round-trips
/// under the same definition `merge::is_inversion` used to type it: a span holding
/// an ambiguity code (`R` pairs with `Y`) is an inversion there, so it must apply
/// as one here, and a span holding an unmodelled byte is no inversion at all. Both
/// [`RefApplier::apply`] and its streaming [`RefApplier::reconstructs`] use it.
fn complement(b: u8) -> Option<u8> {
    crate::sequence::complement_base(b)
}

/// Reverse complement under [`complement`]; `None` when any byte is unmodelled.
fn revcomp(b: &[u8]) -> Option<Vec<u8>> {
    b.iter().rev().map(|&c| complement(c)).collect()
}

/// Members in ascending `(ref_start, ref_end)` order — at a shared start a
/// zero-width edit (Identity/Ins) must precede a span so it does not read as an
/// overlap against the cursor. A `Cut`'s members already arrive in that order
/// (both the ruled and canonical paths build them in reference order), which is
/// every call on the production round-trip gate, so the common case borrows them
/// and pays for neither the clone nor the sort; only an out-of-order caller owns
/// and sorts a copy. Both [`RefApplier::apply`] and [`RefApplier::reconstructs`]
/// splice through this, so the fast path has one definition (see
/// `apply_is_order_independent`).
fn in_reference_order(p: &Partition) -> Cow<'_, [Member]> {
    if p.members.is_sorted_by_key(|m| (m.ref_start, m.ref_end)) {
        Cow::Borrowed(&p.members)
    } else {
        let mut v = p.members.clone();
        v.sort_by_key(|m| (m.ref_start, m.ref_end));
        Cow::Owned(v)
    }
}

/// Splices members into the reference in reference order. Members must be
/// sorted-able, non-overlapping, and in-bounds; otherwise returns None.
pub struct RefApplier;

impl SequenceApplier for RefApplier {
    fn apply(&self, reference: &[u8], p: &Partition) -> Option<Vec<u8>> {
        let mut out = Vec::with_capacity(reference.len());
        let mut cursor = 0usize;
        for m in in_reference_order(p).iter() {
            if m.ref_start < cursor || m.ref_end > reference.len() || m.ref_start > m.ref_end {
                return None; // overlap or out of bounds
            }
            out.extend_from_slice(&reference[cursor..m.ref_start]);
            match m.kind {
                EditKind::Identity => out.extend_from_slice(&reference[m.ref_start..m.ref_end]),
                EditKind::Del => {}
                EditKind::Ins | EditKind::Delins | EditKind::Dup | EditKind::Sub => {
                    out.extend_from_slice(&m.inserted)
                }
                EditKind::Inv => {
                    out.extend_from_slice(&revcomp(&reference[m.ref_start..m.ref_end])?)
                }
            }
            cursor = m.ref_end;
        }
        out.extend_from_slice(&reference[cursor..]);
        Some(out)
    }

    /// Allocation-free round-trip: walk the members in reference order and compare
    /// each spliced segment directly against `resulting`, never materialising the
    /// window and short-circuiting on the first mismatch. Behaviour is identical
    /// to `self.apply(reference, p).as_deref() == Some(resulting)` — every branch
    /// mirrors [`apply`](Self::apply), the same structural guards reject an
    /// overlapping/out-of-range member (there returning `None`, here `false`), and
    /// the trailing `matched == resulting.len()` is the exact-length equality
    /// `apply`'s `==` enforces. The equivalence is pinned by
    /// `reconstructs_streaming_matches_apply`.
    fn reconstructs(&self, reference: &[u8], resulting: &[u8], p: &Partition) -> bool {
        let mut cursor = 0usize; // consumed prefix of `reference`
        let mut matched = 0usize; // consumed prefix of `resulting`

        // Compare `expected` against the next unmatched bytes of `resulting`,
        // advancing `matched`; `false` on overrun or mismatch (early exit).
        let take_segment = |matched: &mut usize, expected: &[u8]| -> bool {
            let end = *matched + expected.len();
            if end > resulting.len() || &resulting[*matched..end] != expected {
                return false;
            }
            *matched = end;
            true
        };

        // Same reference-order walk as `apply`, sharing its fast path.
        for m in in_reference_order(p).iter() {
            if m.ref_start < cursor || m.ref_end > reference.len() || m.ref_start > m.ref_end {
                return false; // overlap or out of bounds — `apply` returns None here
            }
            // The identity gap preceding the member.
            if !take_segment(&mut matched, &reference[cursor..m.ref_start]) {
                return false;
            }
            match m.kind {
                EditKind::Identity => {
                    if !take_segment(&mut matched, &reference[m.ref_start..m.ref_end]) {
                        return false;
                    }
                }
                EditKind::Del => {}
                EditKind::Ins | EditKind::Delins | EditKind::Dup | EditKind::Sub => {
                    if !take_segment(&mut matched, &m.inserted) {
                        return false;
                    }
                }
                EditKind::Inv => {
                    // revcomp(reference[m.ref_start..m.ref_end]) compared without
                    // allocating: the k-th output base is complement of the k-th
                    // base from the span's end.
                    let span = &reference[m.ref_start..m.ref_end];
                    let end = matched + span.len();
                    if end > resulting.len() {
                        return false;
                    }
                    for (k, &out_byte) in resulting[matched..end].iter().enumerate() {
                        if complement(span[span.len() - 1 - k]) != Some(out_byte) {
                            return false;
                        }
                    }
                    matched = end;
                }
            }
            cursor = m.ref_end;
        }
        // The trailing identity segment, then exact-length equality.
        if !take_segment(&mut matched, &reference[cursor..]) {
            return false;
        }
        matched == resulting.len()
    }
}

/// Round-trip soundness: does the partition, applied to the reference, rebuild
/// the block's resulting sequence? A hard gate — a failing arm is broken, not ugly.
///
/// Delegates to [`SequenceApplier::reconstructs`], so for [`RefApplier`] the
/// answer is computed by a streaming byte-compare that never materialises the
/// spliced window (the default impl still routes through [`apply`](SequenceApplier::apply)).
pub fn round_trips(applier: &dyn SequenceApplier, ctx: &BlockCtx, p: &Partition) -> bool {
    applier.reconstructs(ctx.reference, ctx.resulting, p)
}

/// Renderability (a hard gate, distinct from round-trip): can this partition be
/// written as valid HGVS at all? The one way it can fail that [`round_trips`]
/// does NOT catch is the #486 coincident-slot violation (the #2203 rederive bug):
/// two content-bearing insertions anchored at the *same* interbase. In this model
/// both `Ins` and a peeled `Dup` are zero-width members (`ref_start == ref_end`),
/// so two of them at one offset are two insertions into one slot — invalid HGVS,
/// which production's #486 guard rejects. [`RefApplier`] is BLIND to it: at a
/// shared offset `p` a zero-width member passes `p < cursor` (`p < p` is false)
/// and its payload is byte-concatenated with the other's, so the partition
/// "round-trips" while being unrenderable.
///
/// Zero-width `Identity` members (Step-5 examined re-attachments) are EXEMPT —
/// they render nothing at the slot (annotation, not content), so they never
/// collide with an insertion there. Span members (`ref_start < ref_end`) occupy a
/// range, not a slot; `RefApplier` already rejects a zero-width member landing
/// inside a span as an overlap, so only the coincident zero-width pair is blind.
pub fn is_renderable(p: &Partition) -> bool {
    // A member counts here iff it is content-bearing and zero-width (an `Ins` or a
    // peeled `Dup` at an interbase); the partition is unrenderable iff two of them
    // share one interbase. Member counts per block are tiny, so this pairwise scan
    // is cheaper than the `HashSet` it replaces — no allocation, and it returns on
    // the common "no collision" case after one pass. It is O(n²), but n is tiny;
    // order-independent, exactly like the set: it fires iff some later counted
    // member repeats an earlier one's offset.
    let counted = |m: &Member| m.ref_start == m.ref_end && m.kind != EditKind::Identity;
    let members = &p.members;
    for (i, a) in members.iter().enumerate() {
        if counted(a)
            && members[i + 1..]
                .iter()
                .any(|b| counted(b) && b.ref_start == a.ref_start)
        {
            return false; // a second insertion at an already-occupied interbase
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::block_ctx::{BlockCtx, FrameContext, Molecule, Provenance};
    use crate::partition::output::Member;

    fn m(kind: EditKind, ref_start: usize, ref_end: usize, inserted: &str) -> Member {
        Member {
            kind,
            ref_start,
            ref_end,
            inserted: inserted.as_bytes().to_vec(),
        }
    }

    /// An `Inv` applies under the IUPAC complement `merge::is_inversion` types it
    /// with: `R` pairs with `Y`, so `RT` inverts to `AY`. A span holding a byte no
    /// complement models is not an inversion, so the partition does not apply.
    #[test]
    fn inversion_applies_under_the_iupac_complement() {
        let inv = |reference: &[u8]| {
            RefApplier.apply(
                reference,
                &Partition {
                    members: vec![m(EditKind::Inv, 0, 2, "")],
                },
            )
        };
        assert_eq!(inv(b"RT"), Some(b"AY".to_vec()));
        assert_eq!(inv(b"ACGG"), Some(b"GTGG".to_vec()));
        // Soft-masked reference bases complement as `is_inversion` compares them.
        assert_eq!(inv(b"ac"), Some(b"GT".to_vec()));
        assert_eq!(inv(b"XT"), None, "an unmodelled byte has no complement");
    }

    /// The #2203 shape: a peeled zero-width `Dup` and a plain zero-width `Ins` at
    /// the SAME interbase. It round-trips (RefApplier byte-concatenates the two
    /// payloads, blind to the collision) yet is unrenderable (#486 coincident
    /// slot). The gate catches exactly what the applier cannot.
    #[test]
    fn coincident_insertions_round_trip_but_are_unrenderable() {
        let reference = b"AC";
        // Dup of "AC" at interbase 2, and an Ins of "G" at interbase 2.
        let p = Partition {
            members: vec![m(EditKind::Dup, 2, 2, "AC"), m(EditKind::Ins, 2, 2, "G")],
        };
        let frame = FrameContext::NonCoding;
        let provenance = Provenance::none();
        let ctx = BlockCtx {
            reference,
            resulting: b"ACACG",
            frame: &frame,
            molecule: Molecule::Dna,
            provenance: &provenance,
        };
        assert!(
            round_trips(&RefApplier, &ctx, &p),
            "RefApplier is blind: the two payloads concatenate to ACACG"
        );
        assert!(
            !is_renderable(&p),
            "two content-bearing insertions at interbase 2 violate #486"
        );
    }

    /// A zero-width `Identity` (a Step-5 examined re-attachment) is annotation, not
    /// content — it may share an interbase with an insertion (Flag A: exempt).
    #[test]
    fn a_zero_width_identity_does_not_collide_with_an_insertion() {
        let p = Partition {
            members: vec![m(EditKind::Ins, 2, 2, "G"), m(EditKind::Identity, 2, 2, "")],
        };
        assert!(is_renderable(&p), "Identity is exempt from the slot rule");
    }

    /// One insertion per interbase is renderable, including two insertions at
    /// DIFFERENT interbases.
    #[test]
    fn insertions_at_distinct_interbases_are_renderable() {
        let one = Partition {
            members: vec![m(EditKind::Ins, 2, 2, "G")],
        };
        let two_distinct = Partition {
            members: vec![m(EditKind::Ins, 1, 1, "G"), m(EditKind::Dup, 3, 3, "AC")],
        };
        assert!(is_renderable(&one));
        assert!(is_renderable(&two_distinct));
    }

    /// The allocation-free `is_renderable` must agree with the `HashSet`-based
    /// reference on every shape — collisions, spans, exempt identities, and
    /// reversed member order — so the pairwise scan is order-independent exactly
    /// like the set it replaced.
    #[test]
    fn is_renderable_matches_a_set_based_reference() {
        fn reference(p: &Partition) -> bool {
            let mut occupied = std::collections::HashSet::new();
            for m in &p.members {
                if m.ref_start == m.ref_end
                    && m.kind != EditKind::Identity
                    && !occupied.insert(m.ref_start)
                {
                    return false;
                }
            }
            true
        }
        let cases: Vec<Vec<Member>> = vec![
            vec![],
            vec![m(EditKind::Ins, 2, 2, "G")],
            vec![m(EditKind::Ins, 2, 2, "G"), m(EditKind::Dup, 2, 2, "AC")],
            vec![m(EditKind::Ins, 2, 2, "G"), m(EditKind::Identity, 2, 2, "")],
            vec![m(EditKind::Ins, 1, 1, "G"), m(EditKind::Dup, 3, 3, "AC")],
            vec![m(EditKind::Del, 1, 3, ""), m(EditKind::Ins, 3, 3, "A")],
            vec![
                m(EditKind::Dup, 5, 5, "AC"),
                m(EditKind::Ins, 2, 2, "G"),
                m(EditKind::Ins, 5, 5, "T"),
            ],
            vec![
                m(EditKind::Ins, 2, 2, "A"),
                m(EditKind::Ins, 2, 2, "B"),
                m(EditKind::Ins, 2, 2, "C"),
            ],
        ];
        for members in &cases {
            let p = Partition {
                members: members.clone(),
            };
            assert_eq!(is_renderable(&p), reference(&p), "mismatch for {members:?}");
        }
    }

    /// `RefApplier::apply` is order-independent: members handed to it out of
    /// reference order splice to the same sequence as the ascending list. This
    /// pins the reordering fallback the fast path relies on — a `Cut`'s members
    /// already arrive ascending (both the ruled and canonical paths build them in
    /// reference order), so `apply` skips the clone-and-sort, but an out-of-order
    /// caller must still get the reference-order splice.
    #[test]
    fn apply_is_order_independent() {
        let reference = b"ACGTACGT";
        let ascending = Partition {
            members: vec![m(EditKind::Del, 2, 4, ""), m(EditKind::Ins, 6, 6, "XX")],
        };
        let shuffled = Partition {
            members: vec![m(EditKind::Ins, 6, 6, "XX"), m(EditKind::Del, 2, 4, "")],
        };
        let expected = b"ACACXXGT";
        assert_eq!(
            RefApplier.apply(reference, &ascending).as_deref(),
            Some(&expected[..]),
            "ascending members splice in reference order",
        );
        assert_eq!(
            RefApplier.apply(reference, &shuffled).as_deref(),
            Some(&expected[..]),
            "shuffled members are reordered before splicing",
        );
    }

    /// The allocation-free `RefApplier::reconstructs` must agree with the
    /// `apply`-based oracle on every input — sound splices, unsound ones,
    /// out-of-range members, coincident insertions, and length mismatches alike.
    /// `apply` stays the reference implementation; this differential test pins the
    /// streaming fast path to it (the same keep-the-oracle-optimise-beside-it
    /// discipline the DAG aligner uses).
    #[test]
    fn reconstructs_streaming_matches_apply() {
        let cases: Vec<(&[u8], Partition)> = vec![
            (b"ACGT", part(&[(EditKind::Sub, 1, 2, "T")])),
            (b"ACGT", part(&[(EditKind::Del, 1, 3, "")])),
            (b"ACGT", part(&[(EditKind::Ins, 2, 2, "GG")])),
            (b"ACGTACGT", part(&[(EditKind::Delins, 2, 5, "TT")])),
            (b"ACGT", part(&[(EditKind::Dup, 2, 2, "AC")])),
            (b"ACGTAC", part(&[(EditKind::Inv, 1, 5, "")])),
            (
                b"ACGT",
                part(&[(EditKind::Identity, 0, 2, ""), (EditKind::Sub, 3, 4, "A")]),
            ),
            (
                b"ACGTACGT",
                part(&[
                    (EditKind::Del, 1, 2, ""),
                    (EditKind::Ins, 4, 4, "X"),
                    (EditKind::Sub, 6, 7, "A"),
                ]),
            ),
            // coincident zero-width insertions: round-trips via byte-concat
            (
                b"AC",
                part(&[(EditKind::Dup, 2, 2, "AC"), (EditKind::Ins, 2, 2, "G")]),
            ),
            // overlapping members -> apply returns None
            (
                b"ACGT",
                part(&[(EditKind::Sub, 1, 3, "XX"), (EditKind::Sub, 2, 4, "YY")]),
            ),
            // out of bounds -> apply returns None
            (b"ACGT", part(&[(EditKind::Sub, 3, 9, "Z")])),
            // descending input order (exercises the reordering fallback)
            (
                b"ACGTACGT",
                part(&[(EditKind::Sub, 6, 7, "A"), (EditKind::Sub, 1, 2, "T")]),
            ),
            // empty partition -> identity of the whole reference
            (b"ACGT", part(&[])),
            // a full inversion at the start
            (b"ACGT", part(&[(EditKind::Inv, 0, 4, "")])),
            // IUPAC and soft-masked inversion spans, and an unmodelled byte (apply
            // returns None): the streaming path must share `apply`'s complement.
            (b"ARTC", part(&[(EditKind::Inv, 1, 3, "")])),
            (b"acGT", part(&[(EditKind::Inv, 0, 2, "")])),
            (b"AXTC", part(&[(EditKind::Inv, 1, 3, "")])),
        ];

        for (reference, p) in &cases {
            let applied = RefApplier.apply(reference, p);

            // Candidate `resulting` strings: shared shapes plus, when the splice
            // succeeds, the true reconstruction and several perturbations of it so
            // both the true and false branches of the oracle are exercised.
            let mut candidates: Vec<Vec<u8>> = vec![
                b"".to_vec(),
                b"A".to_vec(),
                reference.to_vec(),
                b"ACGT".to_vec(),
                b"ACACG".to_vec(),
            ];
            if let Some(seq) = &applied {
                candidates.push(seq.clone());
                let mut longer = seq.clone();
                longer.push(b'N');
                candidates.push(longer);
                if !seq.is_empty() {
                    candidates.push(seq[..seq.len() - 1].to_vec());
                    let mut flipped = seq.clone();
                    flipped[0] = if flipped[0] == b'A' { b'C' } else { b'A' };
                    candidates.push(flipped);
                }
            }

            for cand in &candidates {
                let oracle = applied.as_deref() == Some(cand.as_slice());
                let streaming = RefApplier.reconstructs(reference, cand, p);
                assert_eq!(
                    streaming, oracle,
                    "reconstructs disagreed with the apply oracle\n  ref={:?}\n  members={:?}\n  \
                     cand={:?}\n  oracle={oracle} streaming={streaming}",
                    reference, p.members, cand,
                );
            }
        }
    }

    fn part(specs: &[(EditKind, usize, usize, &str)]) -> Partition {
        Partition {
            members: specs
                .iter()
                .map(|&(kind, s, e, ins)| m(kind, s, e, ins))
                .collect(),
        }
    }

    proptest::proptest! {
        /// The randomised half of the differential pin: `reconstructs` equals the
        /// `apply`-based oracle on random references, arbitrary member lists (which
        /// may overlap, be unsorted, or be malformed with `ref_start > ref_end`),
        /// and candidate `resulting` strings — matching and not. Both sides see the
        /// same members, so the equality must hold whatever the partition's shape.
        #[test]
        fn reconstructs_streaming_matches_apply_on_random_inputs(
            reference in proptest::collection::vec(
                proptest::sample::select(vec![b'A', b'C', b'G', b'T']), 0..24usize),
            raw_members in proptest::collection::vec(
                (
                    proptest::sample::select(vec![
                        EditKind::Identity, EditKind::Del, EditKind::Ins,
                        EditKind::Delins, EditKind::Dup, EditKind::Sub, EditKind::Inv,
                    ]),
                    0..12usize,
                    0..12usize,
                    proptest::collection::vec(
                        proptest::sample::select(vec![b'A', b'C', b'G', b'T']), 0..6usize),
                ),
                0..5usize),
            random_cand in proptest::collection::vec(
                proptest::sample::select(vec![b'A', b'C', b'G', b'T']), 0..12usize),
        ) {
            let p = Partition {
                members: raw_members
                    .into_iter()
                    .map(|(kind, ref_start, ref_end, inserted)| Member {
                        kind,
                        ref_start,
                        ref_end,
                        inserted,
                    })
                    .collect(),
            };
            let applied = RefApplier.apply(&reference, &p);
            let mut candidates: Vec<Vec<u8>> = vec![reference.clone(), random_cand, Vec::new()];
            if let Some(seq) = &applied {
                candidates.push(seq.clone());
                let mut longer = seq.clone();
                longer.push(b'N');
                candidates.push(longer);
                if !seq.is_empty() {
                    candidates.push(seq[..seq.len() - 1].to_vec());
                }
            }
            for cand in &candidates {
                let oracle = applied.as_deref() == Some(cand.as_slice());
                let streaming = RefApplier.reconstructs(&reference, cand, &p);
                proptest::prop_assert_eq!(
                    streaming,
                    oracle,
                    "ref={:?} members={:?} cand={:?}",
                    reference,
                    p.members,
                    cand
                );
            }
        }
    }
}
