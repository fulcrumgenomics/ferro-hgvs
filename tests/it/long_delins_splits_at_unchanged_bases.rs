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
//! covers net insertions. It does not reach net *deletions*, so those keep the
//! original 32 nt bound (`MAX_UNGUARDED_SPLIT_BLOCK`) as their only guard, and
//! the raise applies to the two regimes that are safe: equal-length blocks,
//! which have no alignment choice to get wrong, and net insertions, which the
//! rule covers. Both halves of that split are pinned below.

use crate::common::synthetic::{hgvs, SyntheticBuilder};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

fn normalize(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
    format!(
        "{}",
        normalizer
            .normalize(&variant)
            .unwrap_or_else(|e| panic!("normalize failed for `{input}`: {e}"))
    )
}

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

    let provider = SyntheticBuilder::genomic(&core).build();
    let input = hgvs(
        &format!("NC_TEST.1:g.{{0}}_{{1}}delins{replacement}"),
        &[1, core.len() as u64],
    );
    let out = normalize(provider, &input);

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

/// A long **net deletion** must stay one spanning delins.
///
/// `separations_are_meaningful` only covers net insertions, so nothing else
/// stops `best_alignment` seizing on a coincidentally-surviving base here. The
/// block is 52 nt replaced by 14 — the shape of `delins.md:44-47`'s own worked
/// example, `LRG_199t1:c.850_901delinsTTCCTCGATGCCTG`, which raising the bound
/// across the board split into three members and turned one correct protein
/// consequence into three bogus ones.
///
/// The payload here is drawn from the reference's own alphabet on a periodic
/// core, so coincidental matches are abundant — exactly the regime the spec
/// warns about, and the case a length-blind cap gets wrong.
#[test]
fn a_long_net_deletion_stays_one_spanning_delins() {
    // 52 nt core, replaced by 14 — well past the 32 nt unguarded bound.
    let core = "ACGT".repeat(13);
    assert_eq!(core.len(), 52);
    let replacement = "TTCCTCGATGCCTG";
    assert_eq!(replacement.len(), 14);

    let provider = SyntheticBuilder::genomic(&core).build();
    let input = hgvs(
        &format!("NC_TEST.1:g.{{0}}_{{1}}delins{replacement}"),
        &[1, core.len() as u64],
    );
    let out = normalize(provider, &input);

    // Pin the whole result, for the reason the sibling test above gives:
    // `!contains(';')` plus `contains(replacement)` is satisfied by *any*
    // single-member output, including one whose span or trimming regressed —
    // `1_51`, `2_52`, or a payload with something spliced onto it would all
    // pass while the guard under test had stopped working.
    //
    // The expected answer is fully determined. The block is a net deletion, so
    // `partition_block` caps it at `MAX_UNGUARDED_SPLIT_BLOCK` (32) and returns
    // the whole block as one piece; and neither endpoint trims, because the
    // core begins `A`/ends `T` while the payload begins `T`/ends `G`. So the
    // span is the full core and the payload is untouched.
    let expected = format!("NC_TEST.1:g.1_{}delins{replacement}", core.len());
    assert_eq!(
        out, expected,
        "a 52 nt -> 14 nt delins must stay one spanning member with its payload intact"
    );
}
