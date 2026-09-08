//! Neutral, arm-facing partitioner layer (design §7).
//!
//! `partition` is the home of the arm machinery relocated out of `bakeoff/`: the
//! output types ([`output`]), the composition layer ([`arm`], [`strategy`]), the
//! per-block context ([`block_ctx`]) every arm is a pure function of, the
//! soundness gates ([`metrics`]), and the concrete arm implementations
//! ([`arms`]). The composition layer takes a [`BlockCtx`] — never a
//! measurement-side `Case` — so this layer depends on nothing in `bakeoff`.
//!
//! # Dependency direction
//!
//! `bakeoff → partition` and `normalize → partition`; `partition` reaches back
//! into `normalize` only through the seams below. Each is intentional and
//! confined to the module named, and none of them pulls the high-level
//! normalize/HGVS logic (variant types, providers, rendering) into this layer:
//!
//! * `arms::arms` (dev-only) uses
//!   `crate::normalize::seqfirst` (`align` / `clcs`) — a self-contained, low-level
//!   alignment primitive that has no intra-crate `use crate::` imports of its own
//!   and already serves both production `merge.rs` and these arms. Relocating
//!   `seqfirst` to a neutral home is a deferred follow-on.
//! * `registry` (dev-only) wraps the four legacy `FERRO_PARTITION` rules, whose
//!   `Monolithic` bodies must reach `merge`-private `partition_block*`/`Piece`
//!   internals — so those wrappers are DEFINED in `merge` and only *registered*
//!   there (Task 3, design §8).
//! * [`bridge`] adapts between a [`ruled::Cut`] and `merge`'s `Piece`, and takes
//!   `normalize::ShuffleDirection` for placement.
//! * [`seed`] starts the ruled core from `merge`'s canonical partitioner
//!   primitives (flank trim, whole-span reverse complement).
//! * [`adapters`] wraps the existing `merge` passes and predicates as ruled rules,
//!   so the ruled core executes the shipping pass rather than a re-implementation.
//! * [`driver`] takes `normalize::ShuffleDirection` in its public signature.
//!
//! Every other `partition` module — `block_ctx`/`output`/`metrics`/`partitioner`/
//! `ruled`/`rules`/`fold`/`examined`/`coincidence`, and the dev-only
//! `arm`/`strategy` — stays free of `normalize`.

pub mod adapters;
pub mod block_ctx;
pub mod bridge;
pub mod coincidence;
pub mod driver;
pub mod examined;
pub mod fold;
pub mod metrics;
pub mod output;
pub mod partitioner;
pub mod ruled;
pub mod rules;
pub mod seed;

// The bake-off layer (design §7 composition + arm implementations). Dev-only:
// the step-9 flip production-wires the ruled core above, but the L1×L2×Dial-B
// grid, the arm registry, and the `Strategy` composition remain a measurement
// tool that ships in no release build.
#[cfg(feature = "dev")]
pub mod arm;
#[cfg(feature = "dev")]
pub mod arms;
#[cfg(feature = "dev")]
pub mod registry;
#[cfg(feature = "dev")]
pub mod strategy;

// The arm-facing surface (design §7 "Produces").
#[cfg(feature = "dev")]
pub use arm::{DialBConfig, L1Segmenter, L2Typer, MergeContext, Monolithic};
pub use block_ctx::{
    BlockCtx, ExaminedRegion, FrameContext, IndividuationPolicy, Molecule, Provenance,
};
pub use output::{EditKind, Member, Partition, Segment};
pub use partitioner::{validate_sound, PartitionError, Partitioner};
#[cfg(feature = "dev")]
pub use strategy::Strategy;
