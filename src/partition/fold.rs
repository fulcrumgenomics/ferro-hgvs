//! `fold_coincident_insertions` — the by-construction renderability invariant
//! (design §D6 Layer 1; the #2203 lesson).
//!
//! Two or more content-bearing insertions anchored at the SAME interbase are
//! never valid HGVS (production's #486 coincident-slot guard rejects them), yet
//! [`RefApplier`](crate::partition::metrics::RefApplier) byte-concatenates them and
//! reports a round-trip — so the shape is unrenderable-but-round-tripping (the
//! #2203 `342_343dup;343_344ins` regression). Rather than merely FLAG it (that is
//! [`is_renderable`](crate::partition::metrics::is_renderable)'s job, the tripwire),
//! the pipeline never PRODUCES it: coincident zero-width insertions denote the one
//! insertion they concatenate to, so they are folded into a single `Ins`.
//!
//! In this model both `Ins` and a peeled tandem `Dup` are zero-width
//! (`ref_start == ref_end`), so a `dup`+`ins` collision — the literal #2203 shape —
//! is two fragments of one insertion event and is folded, cited to
//! `duplication-must-ranks-the-label-not-the-partition` (the `dup` label ranks a
//! *change*, never a *partition* that exposes one; the folded member may re-earn
//! `dup` afterward). A reference-consuming member (`Sub`/`Del`/`Delins`/`Inv`,
//! `ref_start < ref_end`) can never occupy an interbase slot, so the peeled
//! `[dup; sub]` (#2175) is protected by geometry — the collision predicate never
//! matches it. Round-trip preservation is definitional: the payloads concatenate
//! in the applier's application order (stable sort by `(ref_start, ref_end)` keeps
//! coincident members in member-vector order), which is exactly the order this
//! fold concatenates them.

use crate::partition::output::{EditKind, Member, Partition};
use std::collections::HashMap;

/// Whether any interbase carries two or more content-bearing zero-width members —
/// the shape [`fold_coincident_insertions`] folds and [`is_renderable`] rejects.
/// The single collision predicate, so the fold and the `Cut::new` fast-path check
/// (design §10 step 2) cannot disagree on what "collides" means.
///
/// [`is_renderable`]: crate::partition::metrics::is_renderable
pub(crate) fn has_coincident_collision(members: &[Member]) -> bool {
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for m in members {
        if m.ref_start == m.ref_end && m.kind != EditKind::Identity {
            let n = counts.entry(m.ref_start).or_default();
            *n += 1;
            if *n >= 2 {
                return true;
            }
        }
    }
    false
}

/// Fold every interbase carrying two or more content-bearing zero-width members
/// into a single `Ins`. A no-op (returns the partition unchanged) when no
/// interbase collides — so the ∅ path, and every already-renderable partition,
/// is byte-identical.
pub(crate) fn fold_coincident_insertions(partition: Partition) -> Partition {
    let members = partition.members;

    if !has_coincident_collision(&members) {
        return Partition { members }; // no collision: exact passthrough
    }

    // Interbase offset -> member indices of the content-bearing zero-width members
    // anchored there, in ascending (member-vector) order = applier order.
    let mut slots: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, m) in members.iter().enumerate() {
        if m.ref_start == m.ref_end && m.kind != EditKind::Identity {
            slots.entry(m.ref_start).or_default().push(i);
        }
    }

    // For each colliding slot: the folded `Ins` replaces the group's FIRST member;
    // the rest are dropped. Payloads concatenate in applier order.
    let mut folded_at_first: HashMap<usize, Member> = HashMap::new();
    let mut dropped: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for (&p, idxs) in &slots {
        if idxs.len() < 2 {
            continue;
        }
        let mut payload = Vec::new();
        for &i in idxs {
            payload.extend_from_slice(&members[i].inserted);
        }
        folded_at_first.insert(
            idxs[0],
            Member {
                kind: EditKind::Ins,
                ref_start: p,
                ref_end: p,
                inserted: payload,
            },
        );
        for &i in &idxs[1..] {
            dropped.insert(i);
        }
    }

    let mut out = Vec::with_capacity(members.len());
    for (i, m) in members.into_iter().enumerate() {
        if let Some(folded) = folded_at_first.remove(&i) {
            out.push(folded);
        } else if !dropped.contains(&i) {
            out.push(m);
        }
    }
    Partition { members: out }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::metrics::{is_renderable, RefApplier, SequenceApplier};

    fn ins(p: usize, s: &str) -> Member {
        Member {
            kind: EditKind::Ins,
            ref_start: p,
            ref_end: p,
            inserted: s.as_bytes().to_vec(),
        }
    }
    fn dup(p: usize, s: &str) -> Member {
        Member {
            kind: EditKind::Dup,
            ref_start: p,
            ref_end: p,
            inserted: s.as_bytes().to_vec(),
        }
    }

    /// The literal #2203 shape: a peeled `dup` and a plain `ins` at one interbase.
    /// Folds to one renderable `Ins`; the applied bytes are unchanged.
    #[test]
    fn dup_ins_collision_folds_to_one_renderable_ins() {
        let reference = b"AC";
        let before = Partition {
            members: vec![dup(2, "AC"), ins(2, "G")],
        };
        assert_eq!(RefApplier.apply(reference, &before).unwrap(), b"ACACG");
        assert!(!is_renderable(&before));
        let after = fold_coincident_insertions(before);
        assert_eq!(after.members, vec![ins(2, "ACG")]);
        assert!(is_renderable(&after));
        assert_eq!(
            RefApplier.apply(reference, &after).unwrap(),
            b"ACACG",
            "round-trip preserved"
        );
    }

    /// `[dup; sub]` (#2175): the `dup` is zero-width, the `sub` consumes a
    /// reference span and can never share the interbase slot — protected by
    /// geometry (the collision predicate never matches it), no label carve-out.
    #[test]
    fn dup_beside_a_sub_is_untouched() {
        let p = Partition {
            members: vec![
                dup(14, "CA"),
                Member {
                    kind: EditKind::Sub,
                    ref_start: 14,
                    ref_end: 15,
                    inserted: b"C".to_vec(),
                },
            ],
        };
        let folded = fold_coincident_insertions(p.clone());
        assert_eq!(folded.members, p.members, "no collision: peel preserved");
    }

    /// Two plain insertions at one interbase fold.
    #[test]
    fn ins_ins_collision_folds() {
        let folded = fold_coincident_insertions(Partition {
            members: vec![ins(3, "T"), ins(3, "A")],
        });
        assert_eq!(folded.members, vec![ins(3, "TA")]);
    }

    /// Order sensitivity: payloads concatenate in applier order (member-vector
    /// order), so "AB" then "CD" folds to "ABCD" — never "CDAB". This is what
    /// makes round-trip preservation definitional.
    #[test]
    fn fold_concatenates_in_applier_order() {
        let folded = fold_coincident_insertions(Partition {
            members: vec![ins(5, "AB"), ins(5, "CD")],
        });
        assert_eq!(folded.members, vec![ins(5, "ABCD")]);
    }

    /// Three coincident members (mixed Dup/Ins) fold to one Ins.
    #[test]
    fn three_coincident_members_fold_to_one() {
        let folded = fold_coincident_insertions(Partition {
            members: vec![dup(7, "A"), ins(7, "T"), dup(7, "G")],
        });
        assert_eq!(folded.members, vec![ins(7, "ATG")]);
    }

    /// No collision => exact passthrough (the ∅ / already-renderable invariant).
    #[test]
    fn no_collision_is_identity() {
        let p = Partition {
            members: vec![
                ins(1, "G"),
                dup(3, "AC"),
                Member {
                    kind: EditKind::Del,
                    ref_start: 5,
                    ref_end: 8,
                    inserted: Vec::new(),
                },
            ],
        };
        assert_eq!(fold_coincident_insertions(p.clone()).members, p.members);
    }
}
