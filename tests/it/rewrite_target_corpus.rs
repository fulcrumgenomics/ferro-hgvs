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
//! - *Confluence* (#1260, #1262): two spellings of one variant used to be required
//!   to reach the same normalized string, and both did. **Under the partition
//!   model #1260 no longer does, and is no longer required to** — its two
//!   spellings assert different partitions of the reference (two adjacent
//!   insertions against one spanning delins), so they are two assertions rather
//!   than two spellings of one question. **#1262 no longer converges either** —
//!   it did until `partition_block_preserving` stopped declining, which is a
//!   later measurement than the one this paragraph originally carried; see
//!   `TARGET_FORMS`. `the_confluence_targets_converge` is therefore replaced by
//!   `each_confluence_target_reaches_one_form_per_asserted_partition`, which pins
//!   **both** strings of both rows exactly (strictly stronger than the
//!   convergence it replaces, which pinned one per row) plus the number of rows
//!   that still agree. See `TARGET_FORMS` for the argument, including what
//!   #1262's convergence used to rest on and why that support is gone.
//!   `every_confluence_target_denotes_one_variant` is **not**
//!   itself a pinned failure — it is a sanity guard on this file's own rows (both
//!   spellings really are the same variant, checked by `cis_apply_oracle::apply`,
//!   independent of the normalizer under test) and it keeps passing after the rewrite
//!   lands; a malformed pin here would prove nothing about the normalizer, which is why the
//!   guard exists at all.
//! - *Denotation* (#1267, #1325, #1394): normalizing must not change what the description
//!   denotes. **All three are now fixed**, so this row has no live target and asserts
//!   nothing here — the 5' half of `clamp_sibling_crossing_junctions` bounds the junction
//!   (#1267), and a repeat outgrowing its tract demotes to the insertion that growth
//!   stands for, on both the junction branch (#1325) and the trapped-deletion branch
//!   (#1394). Each left the PARTITION list rather than being re-pinned; see the notes
//!   where their tests used to be.
//! - *Boundary-validity* (#1274): normalizing must not name a position past the contig's
//!   end. **Fixed** — by #1327's terminal re-spelling, which lands the repair's junction
//!   inside the sequence instead of one past it. Flipped from a pinned failure to a
//!   positive pin, `a_cancelling_edit_at_the_contig_end_stays_inside_the_contig`, per this
//!   file's own contract: when the assertion goes red, the case moves out of the target
//!   list rather than being silenced.
//!
//! Tests here are **pinned failures**, not `#[ignore]`d ones: one passes *because* it
//! asserts the defect reproduces. The moment the rewrite fixes a case, the assertion flips
//! and the test goes red — that's the signal to delete it (or move it to a "fixed" note),
//! never something to silence. As of the #1394 fix every pinned failure has flipped, so
//! what remains is the confluence guard plus three positive pins.
//!
//! ## The other half of the acceptance evidence
//!
//! `tests/it/cis_spelling_confluence_gap.rs` holds the *other* harvested spelling pairs, on
//! the same `cis_apply_oracle` + `padded()` machinery and the same assert-then-flip-red
//! contract as this file. Read the two files together — this one is not the sole record of
//! confluence targets, only of the ones harvested separately from #1235's descendant chain.
//!
//! It started at eight diverging pairs and its `DIVERGENT` *defect* table is still
//! **empty** — but four of those convergences have since been given back on purpose.
//! #1290, #1296, #1308 and #1312 still converge; #1287, #1301, #1320 and #1304 do not,
//! and they sit in that file's `PARTITION_DIVERGENT` table rather than in `DIVERGENT`,
//! for the same reason #1260 and #1262 sit here: under
//! `partition-is-the-unit-of-normalization` their two spellings assert different
//! partitions of one sequence, so two canonical strings is the model working rather
//! than a defect awaiting a fix.
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
//!   #1394 chain (seventeen defects). With #1394 fixed it passes and **no longer carries
//!   `#[ignore]`**, so it now runs on every invocation and criterion 2 is enforced in CI
//!   for indel haplotypes — see its doc comment for what the enabling soak does and does
//!   not establish.
//!
//! So criterion 2 (no overlapping/out-of-order members) for indel haplotypes was, for the
//! whole life of this file, enforced by no test that ran in CI. That is no longer true: the
//! indel property is enabled as of the #1394 fix. The "already exercised" claim the earlier
//! draft made is now roughly accurate — it just was not when it was written.
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

use crate::common::cis_apply_oracle::{apply, normalize};
use crate::common::synthetic::padded;

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

/// `(issue, split spelling's output, spanning spelling's output)` — one pinned
/// string per spelling, because the two no longer agree.
///
/// # One of the two convergences survived the partition rule, and one did not
///
/// Both rows were pinned as *converged* until `FERRO_PARTITION` defaulted to
/// `preserve`:
///
/// ```text
/// #1260  both spellings -> g.[258_259insC;259_260insC]     LOST
/// #1262  both spellings -> g.258_259delinsC                HELD, same string
/// ```
///
/// **#1260 diverges now, and that is the model working rather than a defect.**
/// `g.[258_259insC;259_260insC]` asserts **two** changed blocks — insertions at
/// two adjacent gaps, one unchanged nucleotide apart, and the genomic floor
/// merges only members closer than that — while `g.258_260delinsACACA` asserts
/// **one**. The old partitioner reached one answer by discarding both assertions
/// and re-deriving from the resulting sequence, which converges any two
/// spellings of one sequence by construction. `preserve` keeps each, and the
/// spanning spelling's block is unequal-length (3 reference bases against a
/// 5-base payload) so the split move cannot reach inside it either.
///
/// So this is a **lost convergence on a row filed as a confluence defect**, and
/// a representation change for anyone who stored the merged form. It is recorded
/// as an adjudication rather than absorbed: what the partition model says is that
/// #1260's two spellings were never one question. The applier still agrees they
/// denote one sequence — [`every_confluence_target_denotes_one_variant`] pins
/// that independently of the normalizer — so `EquivalenceChecker` remains the
/// fallback for a consumer that needs to know they are the same variant.
///
/// **#1262 has now lost its convergence too, and the note this replaces
/// explains why it had kept it.** That note read: "the reason is a *decline*
/// rather than a preservation: `258A>C` and `260del` are one unchanged
/// nucleotide apart in the input, but the deletion 3'-shifts out of the trimmed
/// block, so `partition_block_preserving` cannot place every member inside the
/// window, returns `None`, and the caller falls back to `partition_block` — the
/// same re-derivation the default is moving away from." That fallback no longer
/// fires. `partition_block_preserving` returns **window** coordinates instead of
/// forcing each member into the trimmed block — the change
/// `partition-is-the-unit-of-normalization` records as driving the arm's decline
/// rate to 0.0% in both directions — so there is nothing left to refuse, and the
/// split spelling keeps its two members:
///
/// ```text
/// #1262  split     was TEMPLATE:g.258_259delinsC   now TEMPLATE:g.[258A>C;263del]
///        spanning     TEMPLATE:g.258_259delinsC       TEMPLATE:g.258_259delinsC
/// ```
///
/// So the old note's warning — "do not read this row as evidence that the
/// partition rule preserves two-member alleles here" — has been overtaken: it
/// now *is* that evidence, and the row it was warning about no longer exists.
/// The stale claim it left behind has now been corrected at its source rather
/// than annotated around. It was **not** in the ruling record — an earlier
/// revision of this comment said "the ruling record's own rationale still says
/// `#1262 still converges`", and `hgvs_spec_normalization_overrides.json`
/// contains no mention of #1262 at all. The sentence lived in
/// `src/normalize/merge/splitter_reproducer_corpus.rs`'s module doc, which said
/// "End to end, though, **#1262 converges**"; that paragraph now records the
/// non-convergence and points here for the pinned strings.
///
/// **#1260 is unmoved by this change**, both strings byte-identical to what they
/// were, which is what says the movement is #1262's fallback disappearing rather
/// than a partitioner-wide shift.
const TARGET_FORMS: &[(&str, &str, &str)] = &[
    (
        "#1260",
        "TEMPLATE:g.[258_259insC;259_260insC]",
        "TEMPLATE:g.259delinsCAC",
    ),
    (
        "#1262",
        "TEMPLATE:g.[258A>C;263del]",
        "TEMPLATE:g.258_259delinsC",
    ),
];

/// How many of [`TARGET_FORMS`] still reach one string from both spellings.
///
/// **Zero of two**, down from one of two and originally two of two. Pinned
/// separately from the strings so neither direction can happen quietly: a row
/// converging again means an assertion was overruled somewhere and needs saying,
/// and a row losing its convergence is a migration for whoever stored the merged
/// form.
///
/// Zero is a stronger claim than one, not a weaker one — it says every row in
/// this table is now a *deliberate* non-confluence licensed by
/// `partition-is-the-unit-of-normalization`, with no row left whose agreement
/// depends on a fallback path. Both rows' split spellings still denote their
/// inputs' bases ([`every_confluence_target_denotes_one_variant`] proves the two
/// spellings are one variant; the fixed-point loop below proves each answer is
/// stable), so this is criterion 1 of #1235 and nothing worse.
const CONVERGING_TARGETS: usize = 0;

/// Each confluence target reaches the pinned form from each of its two
/// spellings, and exactly [`CONVERGING_TARGETS`] of them still agree.
///
/// Replaces `the_confluence_targets_converge`, whose name asserted a property
/// that now holds for one row of the two — see [`TARGET_FORMS`] for the argument
/// and for what each row used to print.
///
/// This is a **stronger** guard than the one it replaces, not a relaxed one:
/// convergence pinned one string per row, this pins two, so a change that moves
/// either spelling anywhere fails here. The convergence count is pinned on top,
/// so a row converging or diverging is caught as well as a string moving.
#[test]
fn each_confluence_target_reaches_one_form_per_asserted_partition() {
    assert_eq!(
        TARGET_FORMS.len(),
        CONFLUENCE_TARGETS.len(),
        "every target row needs its two forms"
    );
    let mut converging = 0usize;
    for ((issue, core, a, b), (echo, want_a, want_b)) in CONFLUENCE_TARGETS.iter().zip(TARGET_FORMS)
    {
        assert_eq!(issue, echo, "the two tables fell out of order");
        let seq = padded(core);
        let (norm_a, norm_b) = (normalize(&seq, a), normalize(&seq, b));
        assert_eq!(
            norm_a, *want_a,
            "{issue}: the split spelling `{a}` moved. Both spellings are \
             pinned here; a move is a representation change whichever way it \
             goes."
        );
        assert_eq!(
            norm_b, *want_b,
            "{issue}: the spanning spelling `{b}` moved. Both spellings are \
             pinned here; a move is a representation change whichever way it \
             goes."
        );
        // Both pinned forms are fixed points. Preservation that is not a fixed
        // point means the second pass re-partitions the member the first pass kept,
        // and neither exact-string assertion above can see that: they only pin
        // where the *first* pass landed. #1260's spanning spelling is the site most
        // exposed to it — it now lands on the trimmed `TEMPLATE:g.259delinsCAC`
        // rather than the authored `g.258_260` span, so a further trim on a second
        // pass is the plausible failure.
        for (label, form) in [("split", &norm_a), ("spanning", &norm_b)] {
            let again = normalize(&seq, form);
            assert_eq!(
                &again, form,
                "{issue}: the {label} spelling's answer `{form}` is not a fixed point \
                 (re-normalized to `{again}`), so a second pass moves what the first \
                 pass settled"
            );
        }
        if norm_a == norm_b {
            converging += 1;
        }
    }
    assert_eq!(
        converging, CONVERGING_TARGETS,
        "the number of target rows whose two spellings agree moved. Converging \
         one is not a free win — the two spellings assert different partitions, \
         so agreement means one assertion was overruled and the PR must say \
         which. Diverging one is a migration for whoever stored the merged form."
    );
}

// The *denotation* row is now empty. Both of its targets are fixed, and per this file's
// contract they moved out of the PARTITION list rather than being re-pinned:
//
// - **#1325** — a repeat whose growth exceeds its tract is re-spelled as the insertion that
//   growth stands for, which restores a junction `clamp_sibling_crossing_junctions` can
//   bound, so `g.[262_263insAA;263_264insAA;264_265insC]` reaches
//   `g.[263_264insAAAA;264_265insC]`.
// - **#1394** — the base-claiming branch of that same guard, where the swallowed sibling is
//   a deletion the tract holds in both snapshots, so the clamp has nowhere to pull it back
//   to; `g.[262del;263_264insAA;264_265insAA]` reaches `g.262_264A[6]`.
//
// The positive assertions live with the rest of their shape, in
// `issue_1325_repeat_growth_swallows_junction.rs` and
// `issue_1394_repeat_growth_swallows_deletion.rs`. `assert_denotation_currently_broken`
// went with them: it had no remaining caller, and a helper kept alive for a row with no
// members is exactly the corpus rot this file's self-checks exist to catch.

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
