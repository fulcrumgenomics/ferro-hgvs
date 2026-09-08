//! Concrete partitioner arm implementations (design §7).
//!
//! Each module is one L1 segmenter, L2 typer, Dial-B pass, or named entrant,
//! relocated unchanged from `bakeoff/`. [`arms`] additionally uses
//! `crate::normalize::seqfirst` (`align` / `clcs`) — the single intentional
//! `partition → normalize` edge, treated as a shared low-level alignment
//! primitive (see the crate-level note on [`crate::partition`]).

pub mod anchor_peel;
// The `arms` module keeps its name (the control-arm implementations) inside the
// `arms/` directory; the inner-module-inception lint is expected here.
#[allow(clippy::module_inception)]
pub mod arms;
pub mod dial_b;
pub mod l1_align;
pub mod l2_typers;
pub mod op_extract;
pub mod r2_winner;
