//! Sequence-first cis-allele canonicalization: alignment and partitioning.
//!
//! This module derives member boundaries from the **denoted sequence** rather
//! than from how the input allele was spelled. Given a reference block and the
//! alternate block it becomes, it builds the DAG of every minimal-cost
//! alignment and cuts the reference at the edges common to all of them.
//!
//! Two spellings of one variant produce the same `(ref_block, alt_block)` pair,
//! so they produce the same partition. Confluence is structural here, not a
//! property that has to be tested for and repaired.
//!
//! # The two rules
//!
//! 1. **Dominance.** A reference position counts as unchanged only if *every*
//!    minimal alignment matches it; a junction has a forced insertion only if
//!    *every* minimal alignment inserts there. A coincidental match is not in
//!    every minimal alignment, so coincidental splits are refused structurally
//!    — without gap penalties, a block-size cap, or a net-insertion
//!    restriction.
//! 2. **Separation.** Two runs of change merge when separated by fewer than
//!    [`MIN_SEPARATION`] *unchanged reference bases*.
//!
//! Both are load-bearing and neither suffices alone. See
//! `partition::partition_members` for the measured evidence.
//!
//! # Coordinates
//!
//! 0-based half-open throughout, matching [`crate::normalize::shuffle`].
//! Reference offsets are relative to the start of `ref_block`.

// The sequence-first path IS compiled into production: `normalize_allele` runs
// it as a shadow comparison. It is gated behind `FERRO_SEQFIRST_SHADOW=1` and
// never affects output, so the migration step that remains is promoting it from
// shadow to authoritative — not wiring it up at all, which the earlier version
// of this comment claimed.
//
// The allowance is therefore scoped to the items that genuinely have no
// non-test, non-benchmark caller yet, rather than blanketed over the module: a
// module-wide `allow` would also hide a helper that becomes unreachable when the
// shadow is promoted, which is exactly when the warning is worth having.
//
// The module (and `align`'s `AlignmentDag::build`/`edit_distance`, needed by
// the `seqfirst_align` criterion benchmark) is `pub` only under the `dev`
// feature so an external benchmark crate can reach it; it stays `pub(crate)`
// otherwise, since there is no public API here yet.

#[cfg(feature = "dev")]
pub mod align;
#[cfg(not(feature = "dev"))]
pub(crate) mod align;
pub(crate) mod partition;

/// Unchanged reference bases two runs of change must be separated by before the
/// split between them is believed.
///
/// `2`, not `1`. The HGVS spec's general rule (`general.md:34`) says variants
/// "separated by one or more nucleotides" are described individually, but its
/// own exception (`general.md:35`) already carves out "separated by one
/// nucleotide, together affecting one amino acid" as a delins — i.e. a single
/// unchanged base between two runs merges them. `DNA/delins.md:42` applies
/// exactly that exception, rejecting the individually-described
/// `c.[145C>T;147C>G]` in favor of `LRG_199t1:c.145_147delinsTGG`. The spec's
/// worked example at `delins.md:44` (`LRG_199t1:c.850_901`) has that same
/// one-unchanged-base gap and is only reproducible when it merges, so
/// `MIN_SEPARATION` is `2` — merge when fewer than 2 unchanged bases separate
/// two runs — not `1`.
pub(crate) const MIN_SEPARATION: u32 = 2;
