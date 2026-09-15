//! Soundness primitives: round-trip and renderability — the two hard gates
//! (design §6.3). Pure over `(reference, resulting)` / a partition, so they live
//! in `partition` and are shared by both the harness and any production caller.
//! The PROFILE / scoring metrics (agreement, member counts, idempotence) stay in
//! `bakeoff::metrics`, since they are a measurement concern.

use crate::partition::block_ctx::BlockCtx;
use crate::partition::output::{EditKind, Partition};

/// Reconstruct the resulting sequence by applying a partition to a reference.
/// Injected so metrics are testable without ferro's production renderer.
pub trait SequenceApplier {
    fn apply(&self, reference: &[u8], p: &Partition) -> Option<Vec<u8>>;
}

fn revcomp(b: &[u8]) -> Vec<u8> {
    b.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => *other,
        })
        .collect()
}

/// Splices members into the reference in reference order. Members must be
/// sorted-able, non-overlapping, and in-bounds; otherwise returns None.
pub struct RefApplier;

impl SequenceApplier for RefApplier {
    fn apply(&self, reference: &[u8], p: &Partition) -> Option<Vec<u8>> {
        let mut out = Vec::with_capacity(reference.len());
        let mut cursor = 0usize;
        let mut members = p.members.clone();
        // (ref_start, ref_end): at a shared start, a zero-width edit (Identity/Ins)
        // must precede a span so it does not read as an overlap against the cursor.
        members.sort_by_key(|m| (m.ref_start, m.ref_end));
        for m in &members {
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
                    out.extend_from_slice(&revcomp(&reference[m.ref_start..m.ref_end]))
                }
            }
            cursor = m.ref_end;
        }
        out.extend_from_slice(&reference[cursor..]);
        Some(out)
    }
}

/// Round-trip soundness: does the partition, applied to the reference, rebuild
/// the block's resulting sequence? A hard gate — a failing arm is broken, not ugly.
pub fn round_trips(applier: &dyn SequenceApplier, ctx: &BlockCtx, p: &Partition) -> bool {
    applier.apply(ctx.reference, p).as_deref() == Some(ctx.resulting)
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
    let mut occupied = std::collections::HashSet::new();
    for m in &p.members {
        let content_bearing_zero_width = m.ref_start == m.ref_end && m.kind != EditKind::Identity;
        if content_bearing_zero_width && !occupied.insert(m.ref_start) {
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
}
