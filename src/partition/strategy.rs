//! Strategy composition: a grid cell (L1 × L2 × DialB) or a monolithic entrant.

use crate::partition::arm::{DialBConfig, L1Segmenter, L2Typer, MergeContext, Monolithic};
use crate::partition::block_ctx::BlockCtx;
use crate::partition::examined::reattach_examined;
use crate::partition::fold::fold_coincident_insertions;
use crate::partition::output::Partition;

/// A named strategy: either a grid cell (L1 × L2 × DialB) or a monolithic entrant.
///
/// The trait objects are `Send + Sync` so `run_matrix` can share `&[Strategy]`
/// across the rayon pool AND so a `Strategy` (and any wrapping `Box`/`Arc`) is
/// itself `Send + Sync` — which the process-wide `partition::registry` needs to
/// hold one in a `static OnceLock`, and the typed-config handle needs to move one
/// across threads. Every arm this crate builds is a plain data struct, so `Send`
/// costs nothing.
pub enum Strategy {
    Grid {
        name: String,
        l1: Box<dyn L1Segmenter + Send + Sync>,
        l2: Box<dyn L2Typer + Send + Sync>,
        dial_b: DialBConfig,
    },
    Mono {
        name: String,
        arm: Box<dyn Monolithic + Send + Sync>,
    },
}

impl Strategy {
    pub fn name(&self) -> &str {
        match self {
            Strategy::Grid { name, .. } => name,
            Strategy::Mono { name, .. } => name,
        }
    }

    /// Run the strategy on one block context. For a grid cell: L1 segments, L2 types each
    /// segment, then Dial-B merge-back runs over the full member list. Whatever an
    /// arm produces, coincident zero-width insertions are folded into one (the
    /// by-construction renderability invariant, design §D6 Layer 1 — a no-op unless
    /// the arm exposed a #486 slot collision), then the examined-`=` provenance
    /// channel is re-attached last, as an annotation over the definite partition
    /// (design §3 ③) — both no-ops under the ∅ path.
    pub fn run(&self, ctx: &BlockCtx) -> Partition {
        let partition = match self {
            Strategy::Grid { l1, l2, dial_b, .. } => {
                let segs = l1.segment(ctx.reference, ctx.resulting, ctx.molecule, ctx.provenance);
                let mut members = Vec::new();
                for s in &segs {
                    members.extend(l2.type_segment(s, ctx.reference, ctx.resulting, ctx.frame));
                }
                let members = l2.merge_back(
                    members,
                    ctx.reference,
                    ctx.resulting,
                    &MergeContext::of(ctx),
                    dial_b,
                );
                Partition { members }
            }
            Strategy::Mono { arm, .. } => arm.partition(ctx.reference, ctx.resulting, ctx.frame),
        };
        let partition = fold_coincident_insertions(partition);
        reattach_examined(
            partition,
            &ctx.provenance.examined,
            ctx.reference,
            ctx.resulting,
        )
    }
}
