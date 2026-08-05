//! A block the partitioner collapses *by policy* must be reachable from the
//! split spelling too, or one variant has two normal forms.
//!
//! `partition_block` deliberately returns the whole block when every gap
//! between pieces is a single unchanged base and the split buys no
//! higher-priority type — `every_separation_is_a_single_base` plus
//! `split_buys_no_higher_priority_type`. The reasoning is that a lone matched
//! base inside a length-changing block is more likely coincidence than
//! structure, so the spanning `delins` is the better description
//! (`delins.md:44-47`).
//!
//! That is a policy, and it was unreachable from half its own inputs:
//!
//! ```text
//! reference  ("ACGT" x 64) + "ACTAC" + ("ACGT" x 64)
//!
//! g.257_261delinsGGGTGG              stable — the policy's own preferred form
//! g.[257_258delinsGGG;260_261delinsGG]  stable — and denotes the same GGGTGG
//! ```
//!
//! Both were fixed points, so one variant had two stable normal forms. The
//! spanning form is necessarily *heavier* than the split it replaces — merging
//! across an unchanged base always is — and `canonicalize_from_sequence`'s
//! weight bound refuses any derivation heavier than the input. So the policy
//! could never be applied to an input that had already arrived split.
//!
//! The bound guards against derivations that are *accidentally* worse: the
//! single-gap aligner parking its lone gap at one end and reading the offset as
//! a run of substitutions. A rule that widens on purpose is not that, so
//! `partition_block` now reports whether its collapse was `ByPolicy`, and only
//! that case is exempt.
//!
//! The narrowness is the point. #999's `g.[306dup;308C>A]` is also a two-member
//! input whose derived pieces collapse to one, and it must stay split — its
//! split buys a `dup`, a higher-priority type, so `split_buys_no_higher_priority_type`
//! is false, nothing is recorded as `ByPolicy`, and the bound still refuses.

use crate::common::synthetic::assert_padded_preserving;

/// Five bases with a single `T` at 259, which is the coincidental match.
const CORE: &str = "ACTAC";

#[test]
fn the_split_spelling_reaches_the_policys_own_spanning_form() {
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[257_258delinsGGG;260_261delinsGG]");
    assert_eq!(output, "NC_TEST.1:g.257_261delinsGGGTGG");
}

#[test]
fn the_spanning_spelling_is_still_a_fixed_point() {
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.257_261delinsGGGTGG");
    assert_eq!(output, "NC_TEST.1:g.257_261delinsGGGTGG");
}
