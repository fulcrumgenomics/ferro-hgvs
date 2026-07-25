//! Cross-document rules for changes that are *contiguous* on the reference.
//!
//! Companion to `cross_doc_compliance.rs`. Those probes need no reference
//! sequence; these two do, because the rule only becomes observable once
//! the reference says the changes abut.
//!
//! Neither case is described where an implementer would look for it.
//! `DNA/substitution.md` says nothing about two substitutions that happen
//! to be adjacent, and `DNA/deletion.md` says nothing about two deletions
//! that happen to be adjacent. The governing text is elsewhere:
//!
//!   - `DNA/inversion.md:5` — an inversion is "a sequence change where,
//!     compared to a reference sequence, **more than one nucleotide**
//!     replacing the original sequence is the reverse complement of the
//!     original sequence." Two adjacent substitutions can satisfy that
//!     definition, and when they do the pair *is* an inversion.
//!   - `general.md:41` — the 3'rule: "the most 3' position possible of the
//!     reference sequence is arbitrarily assigned to have been changed",
//!     which applies to the merged deletion, not to the input members.
//!
//! Deliberately **not** relied on: `general.md:33-39` discusses variants
//! *separated by* one or more nucleotides, and the broader "less than two
//! nucleotides should be described as a delins" line there is flagged in
//! the spec itself as an SVD-WG proposal, not the current recommendation.
//! Both cases below are separated by **zero** nucleotides — the changes are
//! contiguous, so they describe one edit under the definitions above and do
//! not depend on that pending change.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

/// Normalize `input` 3'-ward against a synthetic genomic reference whose
/// `core` begins at 1-based position 257 on contig `NC_TEST.1`.
fn norm3(core: &str, input: &str) -> String {
    let normalizer = Normalizer::with_config(
        SyntheticBuilder::genomic(core).build(),
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    );
    normalizer
        .normalize(&parse_hgvs(input).expect("parse"))
        .expect("normalize")
        .to_string()
}

/// Two adjacent substitutions whose combined replacement is the reverse
/// complement of the reference pair are an inversion (`DNA/inversion.md:5`).
///
/// Core `GCTCGTTAGCT` puts `ref[259..=260]` = `TC`. Substituting `259T>G`
/// and `260C>A` replaces `TC` with `GA`, and `revcomp(TC) == GA`, so the
/// pair satisfies the inversion definition and must render as `259_260inv`
/// rather than as two substitutions or a `delins`.
///
/// The third member is the control: `267T>A` sits seven nucleotides away,
/// so it is *not* contiguous and must survive as its own member. A merge
/// rule that swept up every member of the allele would fail here.
#[test]
fn adjacent_substitutions_forming_reverse_complement_become_inversion() {
    let out = norm3("GCTCGTTAGCT", "NC_TEST.1:g.[259T>G;260C>A;267T>A]");
    assert_eq!(out, "NC_TEST.1:g.[259_260inv;267T>A]");
}

/// The discriminator for the test above: adjacency alone does **not** make
/// an inversion. `DNA/inversion.md:5` requires the replacement to be the
/// reverse complement of the original, so a contiguous pair that fails that
/// test is a plain `delins`.
///
/// Same core, same positions, one base different: `259T>G` with `260C>G`
/// replaces `TC` with `GG`, and `revcomp(TC) == GA != GG`. The pair is still
/// contiguous — so it still merges into one edit — but it must merge to
/// `259_260delinsGG`, not to `inv`.
///
/// Without this case, the `inv` assertion above would also be satisfied by
/// an implementation that rendered *every* adjacent substitution pair as an
/// inversion.
#[test]
fn adjacent_substitutions_without_reverse_complement_become_delins() {
    let out = norm3("GCTCGTTAGCT", "NC_TEST.1:g.[259T>G;260C>G]");
    assert_eq!(out, "NC_TEST.1:g.259_260delinsGG");
}

/// Two adjacent single-base deletions are one two-base deletion, and the
/// merged deletion is then placed by the 3'rule (`general.md:41`) — not
/// left at the coordinates the input happened to name.
///
/// Core `TCACACAG` puts a `CA` dinucleotide tract at `ref[258..=263]`
/// (`CACACA`), bounded by `T` at 257 and `G` at 264. Deleting 258 and 259
/// individually removes one `CA` unit; merged that is `258_259del`, and
/// since deleting any one of the three `CA` units yields the same sequence
/// (`T` + `CACA` + `G`), the 3'rule assigns the most 3' unit: `262_263del`.
///
/// 264 is `G`, so the shift stops there and never reaches the padding —
/// the expected position is a property of the core, not of the fixture.
#[test]
fn adjacent_single_base_deletions_merge_then_shift_three_prime() {
    let out = norm3("TCACACAG", "NC_TEST.1:g.[258del;259del]");
    assert_eq!(out, "NC_TEST.1:g.262_263del");
}
