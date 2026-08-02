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

// This module is complete and tested but not yet wired into `normalize_allele`;
// that is the migration step. Remove this allow when the wiring lands.
#![allow(dead_code)]

pub(crate) mod align;
pub(crate) mod partition;

/// Unchanged reference bases two runs of change must be separated by before the
/// split between them is believed.
///
/// `2`, not `1`. The HGVS spec contradicts itself here: `general.md:34` says
/// variants "separated by one or more nucleotides" are described individually,
/// while its own NOTE says the rule is moving to "separated by less than two →
/// delins". The spec's worked example (`delins.md:44`, `LRG_199t1:c.850_901`)
/// is only reproducible under the second reading, so that is the one used.
pub(crate) const MIN_SEPARATION: u32 = 2;
