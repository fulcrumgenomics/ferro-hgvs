//! Acceptance corpus for the cis-allele normalizer's sequence-first rewrite (v1, genomic
//! axis).
//!
//! The rewrite (#1235) replaces "normalize each member independently, then repair the
//! damage with seven passes" with "apply the allele, align the result to the reference,
//! partition at the alignment's dominators, render." This module pins the open defects
//! that rewrite is meant to fix, and tries to say honestly what it does and does not cover.
//! An earlier draft of this file claimed to be the one place every open defect the rewrite
//! is not meant to fix was recorded; an independent audit found six in-class defects it had
//! simply omitted, one BOUNDARY exclusion it had misclassified, and two more BOUNDARY
//! exclusions whose stated *reason* was wrong even though the exclusion itself was right.
//! Corrected below; this paragraph is left in as the record of that.
//!
//! ## PARTITION targets asserted here
//!
//! Every assertion below is sequence-level — never a hardcoded "correct" output string,
//! since these are open defects and nobody has stated what the fixed rendering should be.
//! Three forms:
//!
//! - *Confluence* (#1260, #1262): two spellings of one variant must reach the same
//!   normalized string. **Both now do**, so `the_confluence_targets_converge` is a
//!   positive pin rather than a pinned failure.
//!   `every_confluence_target_denotes_one_variant` is **not**
//!   itself a pinned failure — it is a sanity guard on this file's own rows (both
//!   spellings really are the same variant, checked by `cis_apply_oracle::apply`,
//!   independent of the normalizer under test) and it keeps passing after the rewrite
//!   lands; a malformed pin here would prove nothing about the normalizer, which is why the
//!   guard exists at all.
//! - *Denotation* (#1394): normalizing must not change what the description denotes. It
//!   does, today. `assert_denotation_currently_broken` asserts the normalized output
//!   denotes no sequence at all, because #1394's members overlap and the independent
//!   applier declines to splice them. **#1267 and #1325 were the other halves of this row
//!   and are fixed** — the 5' half of `clamp_sibling_crossing_junctions` bounds the
//!   junction, and a repeat outgrowing its tract now demotes to the insertion that growth
//!   stands for; both cases left the PARTITION list rather than being re-pinned, see the
//!   notes where their tests used to be.
//! - *Boundary-validity* (#1274): normalizing must not name a position past the contig's
//!   end. **Fixed** — by #1327's terminal re-spelling, which lands the repair's junction
//!   inside the sequence instead of one past it. Flipped from a pinned failure to a
//!   positive pin, `a_cancelling_edit_at_the_contig_end_stays_inside_the_contig`, per this
//!   file's own contract: when the assertion goes red, the case moves out of the target
//!   list rather than being silenced.
//!
//! Apart from the guard noted above, every test here is a **pinned failure**, not a
//! `#[ignore]`d one: it currently passes *because* it asserts the defect reproduces. The
//! moment the rewrite fixes a case, the assertion flips and the test goes red — that's the
//! signal to delete it (or move it to a "fixed" note), never something to silence.
//!
//! ## The other half of the acceptance evidence
//!
//! `tests/it/cis_spelling_confluence_gap.rs` holds the *other* harvested spelling pairs, on
//! the same `cis_apply_oracle` + `padded()` machinery and the same assert-then-flip-red
//! contract as this file. Read the two files together — this one is not the sole record of
//! confluence targets, only of the ones harvested separately from #1235's descendant chain.
//!
//! It started at eight diverging pairs and now has **none**: #1287, #1290, #1301, #1308,
//! #1312 and #1320 converged once the sequence-first pass stopped refusing a derivation
//! that collapses to a single pure insertion; #1304 followed once `main`'s removal of the
//! input-separator veto (#1345) let its three-member spelling merge as well; and #1296 —
//! the last — converged once the derivation was allowed to read a repeat member instead of
//! refusing the group that carried it.
//!
//! ## #1235 itself: what's actually running
//!
//! #1235 is the tracking/design issue for the whole rewrite and supplies no reproducer of
//! its own (a table of one-line sketches, not full HGVS strings). Its three acceptance
//! criteria are confluence, no overlapping/out-of-order members, and idempotence. An
//! earlier draft of this doc said these are "exactly what `cis_allele_confluence_proptest.rs`
//! ... already exercise" — true only for half of it:
//!
//! - The **substitution model** (`Haplotype`/`haplotype_strategy`) backs three properties
//!   that run on every test invocation: `all_encodings_of_one_haplotype_converge`,
//!   `the_converged_form_is_a_fixed_point`, and `members_are_disjoint_and_ascending`. These
//!   genuinely exercise confluence, idempotence, and ordering — for substitution-only,
//!   length-preserving haplotypes.
//! - The **indel model** (`IndelHaplotype`/`indel_haplotype_strategy`) is what reaches
//!   insertions, deletions, and the shifting machinery — and its confluence/no-overlap
//!   property, `an_indel_haplotype_normalizes_to_its_own_sequence`, is `#[ignore]`d and its
//!   own doc comment says plainly: "this currently fails, and the failure is real." It is
//!   the property that found the entire #1286 -> #1287 -> #1290 -> #1292 -> #1296 -> #1301
//!   -> #1304 -> #1297 -> #1308 -> #1312 -> #1316 -> #1320 -> #1321 -> #1323 -> #1325 ->
//!   #1394 chain (seventeen defects). It does not run by default; enabling it is
//!   `cargo test --features dev -- --ignored an_indel_haplotype`.
//!
//! So criterion 2 (no overlapping/out-of-order members) for indel haplotypes is currently
//! enforced by no test that runs in CI. #1394 below is the live failure that property last
//! found; there is nothing further to assert for #1235 itself beyond that, but "already
//! exercised" overstated what's actually wired in.
//!
//! ## PARTITION issues not asserted
//!
//! - **#1271** — closed. The issue's own worked example (the HGVS spec's `LRG_199` delins)
//!   was never a live defect here: the regression that used to split it four ways was
//!   already held by a regime-aware bound (`MAX_UNGUARDED_SPLIT_BLOCK`, per the issue's own
//!   account). What #1271 tracked was the *principled* follow-up — extending
//!   `separations_are_meaningful` to net deletions and retiring that bound — and that has
//!   since landed, so the block is now refused on its merits rather than by length. It is
//!   pinned at the unit level by `partition_block_refuses_a_coincidental_net_deletion_split`
//!   and end-to-end by `long_delins_splits_at_unchanged_bases`, and the spec enumeration
//!   records the move (`projection-splits-single-member` 10 -> 9). Still nothing to add to
//!   *this* corpus: the illustrative example is stated as a coding-axis (`c.`) description
//!   whose split coordinates are relative rather than full genomic positions, so there is no
//!   sequence content to build a synthetic genomic reproducer from without inventing it.
//! - **#1284** — three cis repair passes (`respell_colliding_duplications`,
//!   `coalesce_members_at_one_junction`, and the junction arm of
//!   `demote_repeats_spanning_siblings`) are gated `CisKind::Genome | Mt`, so a `c.`/`n.`/
//!   `r.` cis allele that needs one of them emits a description ferro's own parser rejects
//!   — and with `FERRO_ASSERT_REPARSE` armed, that panics and aborts the run instead of
//!   failing with a diff. Out of a genomic-only v1's *axis* scope, so not asserted here —
//!   recorded, not silently dropped, because a v1 rewrite that only wires up `g.` leaves
//!   this exact hole open on the other axes.
//! - **#1328** — `junction_rank` (`src/hgvs/variant.rs`) has no `HgvsVariant::Protein` arm
//!   (falls through its `match` to `_ => None`), so protein cis members sharing a span tie
//!   and fall back to an alphabetical tie-break that renders `dup` before `ins` — reversing
//!   apply order. Same shape as the closed #1301, on the protein axis. Out of a
//!   genomic-only v1's scope, recorded rather than omitted.
//!
//! ## Coverage gaps in the acceptance evidence itself (not defects in the normalizer)
//!
//! - **#1268** — no committed sweep generates a **three**-member cis allele; both committed
//!   sweeps (`cis_junction_crossing_shift.rs`'s
//!   `no_two_member_allele_normalizes_to_a_different_sequence` and
//!   `repeat_span_sibling_overlap.rs`'s
//!   `no_two_member_allele_normalizes_to_overlapping_members`) are explicitly two-member.
//!   Every recent defect in the #1286 chain, #1394 included, is three-member. This is a gap
//!   in the corpus's own evidence, flagged prominently rather than quietly worked around:
//!   nothing here asserts three-member coverage either, so the gap remains open after this
//!   file lands.
//! - **#1283** — the ~276k-case sibling-crossing sweep is read as covering sibling-crossing
//!   shifts broadly, but cannot reach five named shapes it does not generate (a second
//!   member that is `dup`/`ins` rather than a plain edit, more than two members, siblings
//!   in *trans* rather than in phase, seeds shorter than 20bp, and non-`g.` axes). Same
//!   status as #1268: a gap in the acceptance evidence, not asserted or closed here.
//!
//! `cis_allele_confluence_proptest.rs`'s own module doc already names both gaps ("#1268/
//! #1283's complaint about coverage that reads wider than it is") — this section is not the
//! first place they're written down, but they were missing from this corpus specifically.
//!
//! ## Marginal / render-path — recorded, not asserted
//!
//! - **#1318** — `canonicalize_delins` (`src/normalize/rules.rs`) does case-sensitive byte
//!   comparisons on the delins-to-inversion typing path. A soft-masked (lowercase)
//!   reference can therefore stop a genuine reverse complement from being typed `inv`, so
//!   two spellings of what should be one variant diverge by case alone. Single-member, so
//!   not a cis-allele partitioning defect strictly, but it is confluence-shaped and lives in
//!   the *render* step v1's rewrite inherits unchanged rather than replaces — worth carrying
//!   here rather than dropping for being adjacent to the class instead of squarely in it.
//! - **#1264** — the parser accepts an `ins` whose two anchors are not adjacent positions,
//!   but normalizing or re-rendering the same variant is rejected outbound — an inbound/
//!   outbound asymmetry. Currently recorded only as an oracle-exemption comment in
//!   `normalize_reparse_invariant.rs` (why that file's re-parse oracle passes over it), with
//!   no assertion anywhere that the asymmetry itself holds or has closed. Noted here so it
//!   isn't invisible to a reader of this corpus; not asserted, since fixing it is a parser
//!   grammar question orthogonal to member-partitioning.
//!
//! ## BOUNDARY — excluded on axis/coordinate-layer grounds, not because the rewrite can't
//! ## incidentally affect them
//!
//! None of the three below are asserted here. But two of the three reproducers currently
//! reach the pipeline **through a repair pass the rewrite deletes** — so "the rewrite has no
//! reason to touch this" is the wrong justification for excluding them, even though excluding
//! them from this corpus's *asserted* set is still the right call (no genuine partition
//! defect is being deferred; see each entry).
//!
//! - **#1282** — 5'-shifting a member to position 1 underflows `hgvs_pos_to_index`'s
//!   `(pos - 1) as usize` at the coordinate layer and panics in debug builds (silently wraps
//!   to a garbage index in release). This is a real coordinate-layer bug independent of any
//!   repair pass, and excluding it from the partition-defect set is correct. **Caveat:** its
//!   reproducer (`g.[3_4insT;1T>A]`) is two members the rewrite will very plausibly merge
//!   into one block before the 5'-shift ever reaches position 1 alone — so the panic may
//!   simply stop reproducing once the rewrite lands, with the underflow itself untouched.
//!   That would be **masking**, not fixing: a third outcome distinct from "fixed" and
//!   "still broken" that this corpus does not have machinery to tell apart from a real fix,
//!   because nothing here re-tests `hgvs_pos_to_index` directly. The underflow needs its own
//!   test against the coordinate function, independent of this reproducer, so a masked
//!   #1282 cannot be mistaken for a closed one.
//! - **#1307** — respelling a duplication that ends at a contig's last base as an insertion
//!   places it at a gap that does not exist (`24_25insC` on a 24-base contig). The out-of-
//!   bounds gap is produced by `respell_colliding_duplications` -> `respell_at_gap` — a
//!   repair pass that exists only to fix a collision independent per-member normalization
//!   manufactures. Under sequence-first, `24dup` and `24C>G` both touch position 24 and land
//!   in one changed block, so no collision is manufactured and `respell_at_gap` is never
//!   called: the rewrite plausibly retires this reproducer **incidentally**, not because
//!   anyone decided it was in scope. Not asserted here regardless, since there is no
//!   partition defect being deferred — the repair-pass bug and the rewrite's likely
//!   side-effect on it are both worth knowing, which is why they're written down instead of
//!   just "BOUNDARY, out of scope."
//! - **#1327** — the same `respell_at_gap` gap-placement defect as #1307, reached through
//!   the same two repair-pass callers (`respell_colliding_duplications` and
//!   `coalesce_members_at_one_junction`), but on the mitochondrial/circular (`m.`) axis: a
//!   junction at the contig's last base should wrap to `16569_1` instead of running off the
//!   end. Excluded on **axis scope alone** — circular wraparound is out of a genomic-only
//!   v1's scope regardless of which pass produces the bug or whether the rewrite happens to
//!   retire it too. No concrete reproducer exists in the issue (a code-level report only).
//!
//! ## OTHER
//!
//! - **#1270** — an unconfirmed coding-axis (`c.`) codon-frame-gate asymmetry between the
//!   deletion-to-repeat and insertion-to-repeat paths. The issue's own author states no
//!   case has been constructed that reaches it ("an asymmetry with a plausible consequence
//!   rather than a demonstrated defect"). Not genomic, not reproduced, nothing to assert.
//!
//! ## An unexpected pass, found while building this corpus
//!
//! #1267's issue body gives three lines over reference `ACAAAAAAAACGTACGTACG`. When this
//! corpus was built, only the insertion form (`g.[4A>G;9_10insA]`) still reproduced; the
//! two duplication-form lines — `g.[4A>G;9dup]` -> `g.[4A>G;5dup]` and
//! `g.[5A>G;10dup]` -> `g.[5A>G;6dup]` — were independently checked here (via this
//! module's own applier) and already normalized correctly, `blocks_sibling_shift` having
//! fixed them. They were reported rather than asserted, per the corpus's own rules.
//!
//! All three now pass, the insertion form included, so nothing from #1267 is asserted
//! here any more.

use crate::common::cis_apply_oracle::{apply, normalize, normalize_in};
use crate::common::synthetic::padded;
use ferro_hgvs::ShuffleDirection;

/// `(issue, core, split spelling, spanning/merged spelling)` — the two spellings that
/// currently reach two different fixed points instead of one.
const CONFLUENCE_TARGETS: &[(&str, &str, &str, &str)] = &[
    (
        "#1260",
        "AAAAAA",
        "TEMPLATE:g.[258_259insC;259_260insC]",
        "TEMPLATE:g.258_260delinsACACA",
    ),
    (
        "#1262",
        "AAAAAA",
        "TEMPLATE:g.[258A>C;260del]",
        "TEMPLATE:g.258_260delinsCA",
    ),
];

/// Both spellings in every row must denote the same sequence, or the row is a broken pin
/// rather than evidence about the normalizer (this has happened before in this campaign).
///
/// This is a sanity guard on this file's own rows, not a pinned failure: it passes because
/// the rows are well-formed, and it keeps passing after the rewrite lands.
#[test]
fn every_confluence_target_denotes_one_variant() {
    for (issue, core, a, b) in CONFLUENCE_TARGETS {
        let seq = padded(core);
        let from_a = apply(&seq, a).unwrap_or_else(|| panic!("{issue}: `{a}` does not apply"));
        let from_b = apply(&seq, b).unwrap_or_else(|| panic!("{issue}: `{b}` does not apply"));
        assert_eq!(
            from_a, from_b,
            "{issue}: `{a}` and `{b}` are not the same variant"
        );
    }
}

/// **Both confluence targets have converged.** Positive pins, not pinned failures.
///
/// #1260 converged first, on the split form PR #1285 named, which the two-gap
/// alignment made expressible. #1262 followed when the input-separator veto was
/// removed: the derivation had always reached the same piece set from both of its
/// spellings, and the veto was the only thing answering them differently.
#[test]
fn the_confluence_targets_converge() {
    let expected = [
        ("#1260", "TEMPLATE:g.[258_259insC;259_260insC]"),
        ("#1262", "TEMPLATE:g.258_259delinsC"),
    ];
    for (issue, core, a, b) in CONFLUENCE_TARGETS {
        let seq = padded(core);
        let (norm_a, norm_b) = (normalize(&seq, a), normalize(&seq, b));
        assert_eq!(norm_a, norm_b, "{issue}: `{a}` and `{b}` must converge");
        let (_, want) = expected
            .iter()
            .find(|(id, _)| id == issue)
            .unwrap_or_else(|| panic!("no expected form recorded for {issue}"));
        assert_eq!(norm_a, *want, "{issue}: converged on an unexpected form");
    }
}

/// Assert that normalizing `input` in `direction` currently changes what it denotes: the
/// output either splices to a different sequence than the input, or (when its members now
/// overlap) does not splice to a well-defined sequence at all. Flips to a failing
/// assertion, on purpose, the moment the rewrite fixes the case — that is the point of
/// this corpus.
fn assert_denotation_currently_broken(
    seq: &str,
    input: &str,
    direction: ShuffleDirection,
    issue: &str,
) {
    let actual = normalize_in(seq, input, direction);
    let from_input = apply(seq, input)
        .unwrap_or_else(|| panic!("{issue}: `{input}` has no well-defined resulting sequence"));
    let from_output = apply(seq, &actual);
    assert_ne!(
        from_output,
        Some(from_input),
        "{issue} appears fixed — `{input}` now normalizes to `{actual}` (under {direction:?} \
         shuffle), which denotes the same sequence as the input. Move this case out of the \
         PARTITION target list and delete this test."
    );
}

// #1267 was pinned here as a denotation target and is **fixed** — the 5' half of
// `clamp_sibling_crossing_junctions` bounds the junction, so `g.[4A>G;9_10insA]`
// now reaches `g.3_4insG`, which denotes the input's bases. Per this file's
// contract the case moves out of the PARTITION list rather than being re-pinned;
// the positive assertions live with the rest of their shape in
// `cis_junction_crossing_shift.rs`
// (`an_insertion_junction_does_not_cross_an_upstream_sibling`,
// `an_insertion_junction_does_not_cross_an_upstream_junction_sibling` and
// `a_duplication_junction_does_not_cross_an_upstream_junction_sibling`).

// #1325 was pinned here as a denotation target and is **fixed** — a repeat whose
// growth exceeds its tract is now re-spelled as the insertion that growth stands
// for, which restores a junction `clamp_sibling_crossing_junctions` can bound, so
// `g.[262_263insAA;263_264insAA;264_265insC]` reaches
// `g.[263_264insAAAA;264_265insC]` and denotes the input's bases. Per this file's
// contract the case moves out of the PARTITION list rather than being re-pinned;
// the positive assertions live with the rest of their shape in
// `issue_1325_repeat_growth_swallows_junction.rs`.

/// #1394: the same defect on the other branch of the same guard, and what #1325 was
/// masking. Two insertions collapse into a tract-wide repeat correctly, but a sibling
/// *deletion* 3'-shifts into that tract; the tract grew by more bases than it can express
/// as a duplication, so the demotion pass declines and the repeat is left spanning
/// (swallowing) the deletion's bases. The independent applier refuses to splice the
/// resulting members at all, since they overlap — there is no single well-defined
/// resulting sequence for them, which is itself the defect this pins (#1235's
/// criterion 2).
///
/// #1325's repair does not extend here: an insertion claims no bases and so blocks no
/// sibling's shift, and widening its route to the base-claiming branch costs
/// `issue_1296_repeat_claims_its_bases`' two tests (measured, not assumed).
#[test]
fn a_repeat_growth_exceeding_its_tract_still_swallows_a_sibling_deletion() {
    let seq = padded("CAGGCAAACAGTGAAG");
    assert_denotation_currently_broken(
        &seq,
        "TEMPLATE:g.[262del;263_264insAA;264_265insAA]",
        ShuffleDirection::ThreePrime,
        "#1394",
    );
}

/// The highest 1-based HGVS position named anywhere in `output`'s coordinate part (the
/// substring after the accession's `:`). Used only to check a position against a contig
/// length, not to parse HGVS in general — `TEMPLATE` (this module's only accession) has no
/// digits, so every digit run found belongs to a position, not to the accession.
fn max_hgvs_position(output: &str) -> u64 {
    let coord_part = output.rsplit(':').next().unwrap_or(output);
    coord_part
        .split(|c: char| !c.is_ascii_digit())
        .filter_map(|token| token.parse::<u64>().ok())
        .max()
        .unwrap_or(0)
}

/// #1274, **fixed**: an insertion and a deletion that cancel exactly at a contig's last
/// base used to normalize to a coordinate one past the end (`g.10_11=` on a 10-base
/// contig), because a repair pass derived a span for the (correctly empty) residue
/// independently of the contig's actual length.
///
/// #1327's terminal re-spelling closed it — the repair's landing junction is now expressed
/// through the boundary identity at the last base rather than at the junction past it — and
/// this file's contract says a case that goes green moves out of the target list. So the
/// assertion is inverted rather than deleted: the boundary-validity property is worth a
/// standing guard, and inverting keeps the reproducer that found the defect attached to it.
///
/// Two things are asserted, because the position bound alone is weak. The independent
/// applier never caught this one the way it catches #1394 — an `=` edit's SPDI triple does
/// not exercise the out-of-bounds check a deletion's does, so `apply()` happily returned the
/// (correct) reference — which is exactly why "stays in range" needs pairing with "still
/// denotes the reference", or a fix that clamped the span while corrupting the sequence
/// would pass.
#[test]
fn a_cancelling_edit_at_the_contig_end_stays_inside_the_contig() {
    let seq = "ACGTACGTAA";
    let input = "TEMPLATE:g.[8_9insA;10del]";
    let actual = normalize(seq, input);
    let max_pos = max_hgvs_position(&actual);
    assert!(
        max_pos <= seq.len() as u64,
        "`{input}` -> `{actual}` names position {max_pos}, past the {}-base contig (#1274)",
        seq.len()
    );
    // The pair cancels exactly, so the allele denotes the reference unchanged. Asserted
    // through the independent applier, not by re-normalizing.
    let want = apply(&padded(seq), input).expect("input applies");
    let got = apply(&padded(seq), &actual).expect("output applies");
    assert_eq!(
        got, want,
        "`{input}` -> `{actual}` changed the denoted sequence"
    );
}
