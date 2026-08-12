//! A delins must split at its unchanged bases regardless of how long it is.
//!
//! The canonicalizer originally refused to partition any block longer than
//! 32 nt, on the reading that `delins.md:44-47` keeps a long replacement as one
//! spanning delins. That conflated two different things. The spec's concern
//! there is *coincidental* alignment — "parts of the inserted sequence align" —
//! not length as such, and the separation rule (`delins.md:17`) carries no
//! length qualifier at all.
//!
//! The coincidence concern is handled by `separations_are_meaningful`, which
//! since #1271 covers **every** length-changing regime — net deletions as well
//! as net insertions. Equal-length blocks need no cover at all, having no
//! alignment choice to get wrong. So one bound (`MAX_SPLIT_BLOCK`) now applies
//! everywhere, and the coincidence question is answered by a rule about the
//! derived pieces rather than by length.
//!
//! **And that rule is now scoped by AXIS as well**, which is what the second
//! case below records. Disbelieving a separation is `delins.md:44-47`'s
//! payload-coincidence carve-out, and the operator ruling
//! `delins-payload-coincidence-carve-out-is-coding-dna-scoped` scopes that
//! passage to the coding DNA axis. Both cases here are on `NC_TEST.1:g.`, so
//! both now split — the second one having pinned the spanning form until the
//! scope was applied.
//!
//! It was not always so: net deletions were once the unguarded regime and were
//! held instead by a second, smaller bound of 32 nt
//! (`MAX_UNGUARDED_SPLIT_BLOCK`) that guarded by accident. #1271 retired it.
//! The cases below pin the outcomes, which are unchanged by that move.

use crate::common::cis_apply_oracle::normalized_preserving;
use crate::common::synthetic::hgvs;

/// A 40 nt delins — comfortably past the old 32 nt cap — whose interior leaves
/// one base untouched. It must be described as two members, exactly as the same
/// shape is at 5 nt.
#[test]
fn a_delins_longer_than_the_old_cap_still_splits() {
    // Core: 40 bases. The replacement changes everything except position 20.
    let core = "ACGT".repeat(10);
    let mut replacement: Vec<char> = core
        .chars()
        .map(|base| match base {
            'A' => 'C',
            'C' => 'A',
            'G' => 'T',
            _ => 'G',
        })
        .collect();
    // Keep index 19 (HGVS position 20 within the core) identical to the
    // reference, so the block has a genuine unchanged interior base.
    replacement[19] = core.as_bytes()[19] as char;
    let replacement: String = replacement.into_iter().collect();

    let input = hgvs(
        &format!("NC_TEST.1:g.{{0}}_{{1}}delins{replacement}"),
        &[1, core.len() as u64],
    );
    // Through the shared oracle, so the output is checked to denote the input's
    // own bases BEFORE the string below is compared. The string pin alone is a
    // re-bless waiting to happen: whoever next moves this expectation only has
    // to make the literal agree with whatever came out, and a partition that
    // dropped or duplicated a base would pass. `normalized_preserving` reaches
    // the bases through `hgvs_to_spdi` and a splice, never through the
    // normalizer, so it cannot agree merely because normalization produced it.
    let out = normalized_preserving(&core, &input);

    // Pin the whole result, not just "it split somewhere". `contains(';')` plus
    // "the spanning payload is gone" is satisfied by *any* partition, including
    // a coincidental one this periodic `ACGT` fixture could plausibly produce —
    // which would let the very defect this file guards against pass as a fix.
    //
    // The expected answer is fully determined: the block is equal-length, so
    // `best_alignment` short-circuits to a positional comparison with no gap to
    // place and no alignment choice to make. Every column differs except index
    // 19, giving exactly two pieces split at that base.
    let first: String = replacement.chars().take(19).collect();
    let second: String = replacement.chars().skip(20).collect();
    let expected = format!("NC_TEST.1:g.[1_19delins{first};21_40delins{second}]");
    assert_eq!(
        out, expected,
        "a 40 nt delins with an unchanged interior base must split at exactly that base"
    );
}

/// A long **net deletion** on a FRAMELESS axis splits at its coincidental
/// matches, because `delins.md:44-47` does not reach that axis.
///
/// # This test asserted the opposite until the carve-out was scoped
///
/// It pinned the spanning form, on the reading that `delins.md:44-47`'s worked
/// example — `LRG_199t1:c.850_901delinsTTCCTCGATGCCTG`, the very block this
/// test plants here — settles the shape wherever it occurs. The operator ruling
/// `delins-payload-coincidence-carve-out-is-coding-dna-scoped` says it does
/// not: `:47`'s stated reason is preventing "incorrect predictions for the
/// consequences on protein level", which has nothing to bite on where no
/// protein is coded, so the passage reaches `c.` and nothing else. The spec's
/// own example is on `c.`; this test is on `NC_TEST.1:g.`, and the axis is the
/// whole difference.
///
/// So the two halves are pinned on the two axes they belong to. The `c.` half —
/// the spec's own block, staying one member — is
/// `merge::tests::the_spec_delins_example_is_never_audited_as_a_split` and the
/// `LRG_199` row of the reproducer corpus, both of which pass the block through
/// `partition_block` with the carve-out in reach. This is the frameless half.
///
/// The payload here is drawn from the reference's own alphabet on a periodic
/// core, so coincidental matches are abundant. That makes it the strongest
/// available statement of the frameless answer rather than a corner case: if
/// `general.md:34` governs unopposed off `c.`, it governs here.
#[test]
fn a_long_net_deletion_splits_on_a_frameless_axis() {
    // 52 nt core, replaced by 14 — well past the 32 nt unguarded bound.
    let core = "ACGT".repeat(13);
    assert_eq!(core.len(), 52);
    let replacement = "TTCCTCGATGCCTG";
    assert_eq!(replacement.len(), 14);

    let input = hgvs(
        &format!("NC_TEST.1:g.{{0}}_{{1}}delins{replacement}"),
        &[1, core.len() as u64],
    );
    // Through the shared oracle, for the reason the sibling above gives: the
    // five-member pin below is a string comparison and nothing in it checks the
    // members still reassemble to the input's bases. That matters most here,
    // where the payload is drawn from the reference's own alphabet on a periodic
    // core — a wrong partition of an `ACGT` tract is exactly the shape that
    // looks plausible in a diff.
    let out = normalized_preserving(&core, &input);

    // Pin the whole result, for the reason the sibling test above gives: a
    // predicate over the string is satisfied by outputs whose span or trimming
    // regressed, so the guard would keep passing after it had stopped working.
    //
    // The expected answer is fully determined. The block reaches the aligner,
    // whose proposed split `separations_are_meaningful` used to refuse — every
    // gap between consecutive pieces is one base wide, against the
    // `RAISED_PIECE_SEPARATION` this block's net change required. That raise is
    // the payload-coincidence carve-out, so it now applies on the coding DNA
    // axis only and the derived pieces stand here. Five members, every
    // separation exactly one unchanged base, which is what `general.md:34` and
    // `delins.md:17` ask for.
    let expected = "NC_TEST.1:g.[1_39delinsT;41A>C;43del;45_49delinsCGATGC;51_52delinsTG]";
    assert_eq!(
        out, expected,
        "a 52 nt -> 14 nt delins on a frameless axis splits at its unchanged bases"
    );
}
