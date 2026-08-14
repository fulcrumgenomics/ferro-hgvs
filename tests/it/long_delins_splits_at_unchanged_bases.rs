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
//! It was not always so: net deletions were once the unguarded regime and were
//! held instead by a second, smaller bound of 32 nt
//! (`MAX_UNGUARDED_SPLIT_BLOCK`) that guarded by accident. #1271 retired it.
//! The cases below pin the outcomes, which are unchanged by that move.

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
    // # MOVED BY #1835 TO 21 MEMBERS, AND RETURNED BY #1878
    //
    // The expectation here is two members split at the unchanged base, on the
    // ground that the block is equal-length, so there is no gap to place and no
    // alignment choice to make: every column differs except index 19.
    //
    // #1835 replaced that with a 21-member split, and its arithmetic was right as
    // far as it went. The core is `ACGT` x10 and the payload is that with A<->C
    // and G<->T swapped, i.e. `CATG` x10, with index 19 restored. Read
    // position-wise, 39 of 40 columns differ: cost 39. Delete the leading `A`
    // instead and the frames line up, so the same block is 1 deletion, 19
    // substitutions and 1 insertion: cost 21. The minimal alignment is the gapped
    // one, by a wide margin rather than by a tie.
    //
    // WHAT `rulings[equal-length-block-column-correspondence-is-unique]` (decided)
    // says is that cost is the wrong question here. Forty reference bases against
    // a forty-base payload have exactly ONE column correspondence, so the
    // changed-column set is a fact about the variant; the gapped alignment is
    // cheaper but it is a re-frame, and its "unchanged" bases are unchanged only
    // relative to a shift the variant does not have. `general.md:34` keys on
    // nucleotides separating two variants, and under the correspondence there is
    // exactly one such nucleotide — index 19. `delins.md:16` types the runs each
    // side of it, and `:17` keeps them individual.
    //
    // The sentence the flip retired — "the block is equal-length, so there is no
    // alignment choice to make" — is restored as a property of the SHAPE rather
    // than of `partition_block`. That is the whole of #1878.
    //
    // The file's thesis is untouched either way: a delins must split at its
    // unchanged bases regardless of length, and both pinnings assert it.
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
    // # RE-PINNED BY THE PARTITION DEFAULT FLIP (#1835) — AND THE TITLE OF THIS
    // TEST IS NOW WRONG FOR A REASON WORTH READING
    //
    // On the `live` arm the block reached the aligner, whose proposed split
    // `separations_are_meaningful` refused — every gap between consecutive pieces
    // is one base wide, against the `RAISED_PIECE_SEPARATION` this block's net
    // change requires — so `partition_block` returned the whole block as one
    // piece. That guard is `partition_block`'s and is not on the canonical path,
    // so under the new default the block splits into fourteen pure deletions.
    //
    // LICENSED TWICE OVER, and both legs are worth stating because either alone
    // would settle it.
    //
    // (1) THE AXIS. `DNA/delins.md:47` is what recommends the spanning form for a
    // payload-coincidence split, and `delins-payload-coincidence-carve-out-is-coding-dna-scoped`
    // (decided) scopes it to the coding DNA axis. This is a `g.` row, so `:47`
    // does not reach it and `general.md:34` governs unqualified. Note the doc
    // above says raising the bound "turned one correct protein consequence into
    // three bogus ones" — that is `:47`'s own stated ground, and the axis ruling's
    // reasoning is precisely that it has nothing to bite on off the coding axis:
    // a genomic reference predicts no protein consequence, correctly or
    // otherwise. The comment describes the `c.`-sited original
    // (`LRG_199t1:c.850_901delinsTTCCTCGATGCCTG`), which is a different row on a
    // different axis and is unaffected.
    //
    // (2) THE MECHANISM, which would settle it even on a `c.` axis. Every one of
    // the fourteen members below is a PURE DELETION — no member supplies a base.
    // `delins-recommendation-reach-when-the-input-arrives-split` (decided,
    // operator ruling 2026-08-12) rules that `:47` reaches a split only where an
    // INSERTED SEQUENCE re-aligned, i.e. where some derived member supplies bases
    // while consuming a different number of reference bases. "A split of pure
    // deletions inserts nothing, so nothing re-aligned, so `general.md:34`/`:17`
    // govern it unqualified." That is this row exactly.
    //
    // WHY THIS DIFFERS FROM `:44-47`'s OWN EXAMPLE, which merges. The record
    // records that `c.850_901delinsTTCCTCGATGCCTG` derives to
    // `c.[850_866del;870_880del;887_892del;895_896delinsC;899_901del]`, whose
    // `895_896delinsC` consumes two reference bases to place one — gap-bearing, so
    // it merges. This fixture's periodic `ACGT` core against `TTCCTCGATGCCTG`
    // happens to derive with no gap-bearing member at all, so the same rule sends
    // it the other way. The two outcomes are consistent; the sequences differ.
    //
    // The test name is left as it is so the history stays greppable, but read the
    // assertion rather than the name: this row now pins that a long net deletion
    // SPLITS on a genomic reference.
    let expected = "NC_TEST.1:g.[1_3del;5_7del;9del;11_13del;15del;17del;20del;22_23del;\
                    25_26del;28_29del;31_33del;35del;37_38del;40_52del]";
    assert_eq!(
        out, expected,
        "a 52 nt -> 14 nt net deletion on a genomic reference splits at its \
         unchanged bases: `delins.md:47` is coding-axis-scoped, and every member \
         here is a pure deletion so `:46`'s re-alignment mechanism never occurred"
    );
}
