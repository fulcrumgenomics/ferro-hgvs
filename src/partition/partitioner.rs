//! The `Partitioner` trait: a name, a fallible `(BlockCtx) -> Partition`, and a
//! declared canonical-coalesce eligibility — plus the round-trip soundness gate
//! every arm's output must clear before it can be returned as `Ok` (design §3/§6).
//!
//! This module is purely additive: nothing here is wired into production, and no
//! shipped behavior changes. It gives Task 3's legacy-rule wrappers and Task 4's
//! registry a common interface to implement against.

use crate::partition::block_ctx::BlockCtx;
use crate::partition::metrics::{is_renderable, round_trips, RefApplier};
use crate::partition::output::Partition;
#[cfg(feature = "dev")]
use crate::partition::strategy::Strategy;

/// Why a partitioner failed to produce a partition. `Unsound` and `Unrenderable`
/// are the two hard gates [`validate_sound`] enforces; `GridTooLarge` is reserved
/// for a grid-cell entrant declining a window beyond its complexity budget (not
/// produced by anything in this task).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PartitionError {
    /// The partition, applied back to the reference, does not reconstruct the
    /// block's resulting sequence (see [`round_trips`]).
    Unsound,
    /// The partition round-trips but cannot be rendered as valid HGVS — e.g. two
    /// content-bearing insertions sharing one interbase slot (see
    /// [`is_renderable`]).
    Unrenderable,
    /// The block exceeds a grid-cell arm's complexity budget.
    GridTooLarge { ref_len: usize, alt_len: usize },
}

/// A partitioning strategy: derives a [`Partition`] from a [`BlockCtx`], subject
/// to the round-trip soundness gate.
pub trait Partitioner {
    /// A stable, human-readable name for this arm.
    fn name(&self) -> &str;
    /// Partition `ctx`, or decline with a [`PartitionError`].
    fn partition(&self, ctx: &BlockCtx) -> Result<Partition, PartitionError>;
    /// Whether this arm takes the canonical coalesce passes (the `delins.md:44-47`
    /// payload-coincidence merge-back). No default: every arm must declare it —
    /// the compile-time half of the enum-exhaustiveness replacement (design §6).
    fn cuts_with_canonical(&self) -> bool;
}

/// Validate a candidate partition against the two hard gates: it must round-trip
/// to `ctx.resulting` when applied to `ctx.reference`, and it must be renderable
/// as valid HGVS. Returns the partition unchanged on success.
pub fn validate_sound(ctx: &BlockCtx, p: Partition) -> Result<Partition, PartitionError> {
    if !round_trips(&RefApplier, ctx, &p) {
        return Err(PartitionError::Unsound);
    }
    if !is_renderable(&p) {
        return Err(PartitionError::Unrenderable);
    }
    Ok(p)
}

#[cfg(feature = "dev")]
impl Partitioner for Strategy {
    fn name(&self) -> &str {
        Strategy::name(self)
    }

    fn partition(&self, ctx: &BlockCtx) -> Result<Partition, PartitionError> {
        validate_sound(ctx, self.run(ctx))
    }

    /// `Grid`: eligible whenever the cell's Dial-B config applies at least one of
    /// the canonical merge-back carve-outs (C1 codon gap-1, C2 net-deletion
    /// payload-coincidence, C3 unequal-length placed-gap) — the Grid arm's
    /// canonical-coalesce eligibility. The exact production semantics for bakeoff
    /// arms (which of these correspond to ferro's shipped `canonical` vs
    /// `canonical-coalesced` rules) are settled in Task 4.
    /// `Mono`: delegates to the wrapped [`Monolithic`](crate::partition::arm::Monolithic)
    /// arm's own declaration.
    fn cuts_with_canonical(&self) -> bool {
        match self {
            Strategy::Grid { dial_b, .. } => dial_b.c1 || dial_b.c2 || dial_b.c3,
            Strategy::Mono { arm, .. } => arm.cuts_with_canonical(),
        }
    }
}

#[cfg(all(test, feature = "dev"))]
mod tests {
    use super::*;
    use crate::partition::arms::r2_winner;
    use crate::partition::block_ctx::{BlockCtx, FrameContext, Molecule, Provenance};
    use crate::partition::output::{EditKind, Member, Partition};

    fn ctx<'a>(r: &'a [u8], s: &'a [u8], f: &'a FrameContext, p: &'a Provenance) -> BlockCtx<'a> {
        BlockCtx {
            reference: r,
            resulting: s,
            frame: f,
            molecule: Molecule::from_axis("g"),
            provenance: p,
        }
    }

    #[test]
    fn a_sound_arm_returns_ok() {
        let (f, p) = (FrameContext::NonCoding, Provenance::none());
        let arm = r2_winner::peel_comp(); // coalesce-structural not required; any sound arm
        let out = arm.partition(&ctx(b"ACGT", b"AGT", &f, &p));
        assert!(out.is_ok(), "a sound partition must be Ok");
    }

    #[test]
    fn an_unsound_partition_is_rejected() {
        // validate_sound is the unit under test; feed it a hand-built partition that
        // cannot reconstruct `resulting` (identity over a changed block).
        let (f, p) = (FrameContext::NonCoding, Provenance::none());
        let c = ctx(b"ACGT", b"AGGT", &f, &p);
        let bad = Partition {
            members: vec![Member {
                kind: EditKind::Identity,
                ref_start: 0,
                ref_end: 4,
                inserted: vec![],
            }],
        };
        assert!(matches!(
            validate_sound(&c, bad),
            Err(PartitionError::Unsound)
        ));
    }

    #[test]
    fn a_slot_collision_is_unrenderable() {
        // Deviates from the brief's literal `resulting` fixture (`b"ACGT"`): the
        // brief's verbatim listing does not round-trip against this partition, so
        // `validate_sound` returns `Err(Unsound)` before `is_renderable` is ever
        // consulted, and the test cannot observe `Unrenderable` at all — see the
        // Task 2 report. `RefApplier` byte-concatenates the two zero-width Ins
        // payloads blindly (that blindness is the whole point of this gate, per
        // `metrics.rs`'s own `coincident_insertions_round_trip_but_are_unrenderable`),
        // so `resulting` must be the reference with both payloads spliced in at
        // interbase 2 for `round_trips` to pass and hand off to `is_renderable`.
        let (f, p) = (FrameContext::NonCoding, Provenance::none());
        let c = ctx(b"ACGT", b"ACACGT", &f, &p);
        let coincident = Partition {
            members: vec![
                Member {
                    kind: EditKind::Ins,
                    ref_start: 2,
                    ref_end: 2,
                    inserted: b"A".to_vec(),
                },
                Member {
                    kind: EditKind::Ins,
                    ref_start: 2,
                    ref_end: 2,
                    inserted: b"C".to_vec(),
                },
            ],
        };
        assert!(matches!(
            validate_sound(&c, coincident),
            Err(PartitionError::Unrenderable)
        ));
    }
}
