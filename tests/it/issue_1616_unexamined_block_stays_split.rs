//! A block `MAX_SPLIT_BLOCK` refused to examine may not be emitted as one
//! spanning member off the coding DNA axis.
//!
//! # The mechanism
//!
//! `partition_block` short-circuits on length:
//!
//! ```text
//! if past_the_split_cap(reference, result) { return whole(); }
//! ```
//!
//! What comes back is **not** a finding that the block holds one variant. It is
//! the absence of a finding — no rule about the derived pieces ever ran, because
//! no pieces were derived. Emitting it regardless asserts "there is no separation
//! here" about a block nothing looked at.
//!
//! Off `c.` that assertion has no licence. `general.md:34` — "two variants
//! separated by one or more nucleotides should be described individually and
//! **not** as a 'delins'" — governs unopposed there, and the one passage that
//! could override it for this shape, `DNA/delins.md:44-47`, reaches `c.` and
//! nothing else per the decided
//! `rulings[delins-payload-coincidence-carve-out-is-coding-dna-scoped]`. So
//! `canonicalize_from_sequence` declines and the per-member pipeline answers,
//! which is the individual form `:34` asks for.
//!
//! # Why it is not the deleted weight bound under another name
//!
//! `rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]` deleted a
//! refusal whose comparand was **the input's own member list**. The gate here has
//! no term for the input: it reads the block's two lengths against a constant and
//! the axis off `AxisFrame`. The `SPANNING` and `SPLIT` inputs below are the same
//! variant spelled two ways and reach the gate identically — which the deleted
//! bound could not do, and which is why both are pinned here.
//!
//! # Why the gate is the axis KIND
//!
//! `AxisFrame::reading_frame` is **true** for a coding `r.`, so a gate written in
//! terms of the frame would extend a `DNA/` clause onto the RNA axis — the trap
//! that ruling names by name, and the defect #1711 already was once. The `c.`
//! control below is the positive half: it must keep the spanning form, or the
//! gate has become a blanket refusal rather than an axis scope.
//!
//! # This test is not vacuous, and that was measured rather than assumed
//!
//! With the call-site gate removed (`if false`), the `g.` row comes back as the
//! single spanning `delins` — the answer of the **rejected** SVD-WG010 proposal,
//! which is what `spec_conformance_axis`'s negative guard counts. See the
//! mutation table in the PR description.

use crate::common::synthetic::{normalize_to_string, padded, SyntheticBuilder, PAD_OFFSET};
use ferro_hgvs::reference::transcript::Strand;

/// Core length, chosen so the trimmed block clears `merge::MAX_SPLIT_BLOCK`
/// (1 024) with room to spare — the short-circuit is the whole subject.
const SPAN: usize = 1_100;

/// Core offsets, 1-based, of the two members and the base between them.
///
/// Stated once here and mapped onto each axis below, so the `g.` case and its
/// `c.` control are provably the same block rather than two hand-written ones.
const DELINS_START: usize = 1;
/// Last reference base the wide `delins` member consumes. 1 032 bases, so the
/// trimmed block clears the 1 024 cap.
const DELINS_END: usize = 1_032;
/// The single unchanged base between the two members — `general.md:34`'s
/// "separated by one or more nucleotides", at its minimum.
const SEPARATOR: usize = DELINS_END + 1;
/// The lone deleted base 3' of the separator.
const DEL_AT: usize = SEPARATOR + 1;

/// What the wide member puts in place of the reference bases it removes.
///
/// Shorter than the span it replaces, so the block is **length-changing** and
/// therefore inside the cap's regime at all — an equal-length block is exempt
/// from the cap by construction and is pinned by
/// `residual_above_cap_confluence`.
const PAYLOAD: &str = "TTGCAACGTTGCAACGTTGCAACGTTGCAA";

/// A deterministic core, built rather than committed.
///
/// The unit's length is co-prime with the pad's period-4 `ACGT`, so no tract
/// runs out of the core into the padding and the 3'-shift has nothing to slide
/// along.
fn core() -> String {
    const UNIT: &[u8] = b"ACGTTGCAAGCT";
    let bases: Vec<u8> = (0..SPAN).map(|i| UNIT[i % UNIT.len()]).collect();
    String::from_utf8(bases).expect("ASCII bases")
}

/// The two members, spelled at `offset` added to each core position.
///
/// `offset` is [`PAD_OFFSET`] on the padded genomic contig and `0` on the `c.`
/// axis, where `cds_start = PAD_OFFSET + 1` already puts core base `p` at `c.p`.
fn split_members(offset: usize) -> String {
    format!(
        "[{}_{}delins{PAYLOAD};{}del]",
        DELINS_START + offset,
        DELINS_END + offset,
        DEL_AT + offset
    )
}

/// The variant spelled as its two individual members, on the genomic contig.
fn genomic_split() -> String {
    format!("NC_TEST.1:g.{}", split_members(PAD_OFFSET as usize))
}

/// The same variant spelled as the one spanning `delins` that covers both
/// members **and the unchanged base between them**.
///
/// Built from the core the split is applied to, so the two spellings cannot
/// drift apart.
fn genomic_spanning(core: &str) -> String {
    let separator = core.as_bytes()[SEPARATOR - 1] as char;
    format!(
        "NC_TEST.1:g.{}_{}delins{PAYLOAD}{separator}",
        DELINS_START + PAD_OFFSET as usize,
        DEL_AT + PAD_OFFSET as usize
    )
}

#[test]
fn a_genomic_block_past_the_cap_keeps_its_members_individual() {
    let core = core();
    let provider = SyntheticBuilder::genomic(&core).build();
    let output = normalize_to_string(provider, &genomic_split());
    assert!(
        output.contains(';'),
        "`general.md:34` governs a frameless axis, and the block was never \
         examined — the two members must stay individual, got {output}"
    );
    assert_eq!(
        output,
        genomic_split(),
        "nothing licenses a re-spelling here: the block is past \
         `MAX_SPLIT_BLOCK`, so no derivation ran at all"
    );
}

#[test]
fn the_spanning_spelling_of_the_same_variant_is_left_alone() {
    let core = core();
    let provider = SyntheticBuilder::genomic(&core).build();
    let output = normalize_to_string(provider, &genomic_spanning(&core));
    assert_eq!(
        output,
        genomic_spanning(&core),
        "the gate declines the derivation; it does not rewrite an input that \
         already arrived spanning. The two spellings reach it identically, which \
         is what makes this an axis gate rather than the input-relative weight \
         bound `rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]` \
         deleted. **This is the confluence the gate costs**, and it is disclosed \
         rather than hidden: one variant, two normal forms, because a block past \
         the cap cannot be split and `:34` forbids merging it off `c.`. Rule 2 \
         outranks rule 3, so the trade is the ruled one"
    );
}

#[test]
fn the_same_block_on_the_coding_dna_axis_still_spans() {
    // `padded(core)` as the transcript with `cds_start = PAD_OFFSET + 1` puts
    // core base `p` at `c.p`, so this is the same block as above with a real
    // 5'UTR ahead of it rather than the identity CDS `cds_start == 1` gives.
    let core = core();
    let provider = SyntheticBuilder::cds(
        &padded(&core),
        PAD_OFFSET + 1,
        PAD_OFFSET + SPAN as u64,
        Strand::Plus,
    )
    .build();
    let output = normalize_to_string(provider, &format!("NM_TEST.1:c.{}", split_members(0)));
    // Pinned as the exact spanning description rather than as "has no `;`".
    // Absence of a separator is satisfied by any single member, including one
    // with the wrong bounds or the wrong payload, so it would pass a build that
    // spanned the block incorrectly — and the whole subject here is which span a
    // block past the cap comes back as.
    let separator = core.as_bytes()[SEPARATOR - 1] as char;
    assert_eq!(
        output,
        format!("NM_TEST.1:c.{DELINS_START}_{DEL_AT}delins{PAYLOAD}{separator}"),
        "`c.` is the one axis `DNA/delins.md:44-47` reaches, so the spanning \
         form `:47` recommends stands — this is the positive half of the axis \
         scope, and without it the gate could be a blanket refusal"
    );
    assert!(
        !output.contains(';'),
        "and it is one member, not a bracket: {output}"
    );
}
