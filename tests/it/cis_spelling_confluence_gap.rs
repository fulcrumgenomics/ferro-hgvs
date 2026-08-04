//! Two spellings of one variant that normalize to two different strings.
//!
//! #1235's criterion 1 requires every encoding of a variant to reach one
//! normalized string. This pair does not, and its *split* spelling is an
//! expectation this repository blessed and shipped.
//!
//! Six of the original eight converged when the sequence-first pass stopped
//! refusing a derivation that collapses to a single pure insertion — #1287,
//! #1290, #1301, #1308, #1312 and #1320 have moved to [`CONVERGED`] below. The
//! refusal existed to protect two capabilities the derivation lacked (the
//! terminal-insertion clamp and `dup` typing); both now live in the piece
//! builder, so every one of those six merged forms is derived from the sequence
//! rather than assembled per member.
//!
//! #1304 makes seven. It needed the input-separator veto gone as well (#1345,
//! landed separately on `main`), because its split spelling is three members
//! whose derived form merges across a base the input left between two of them.
//!
//! #1296 was the eighth and last, and it needed something different from all of
//! them: not a capability the derivation lacked, but permission to *read* its
//! own output. The per-member pipeline promotes `274_275insAA` to `274A[3]`,
//! and `collect_canonical_edits` refused any group carrying a repeat — so the
//! derivation was locked out of that variant forever. With repeats lowered to
//! the `ins`/`del` they denote, both spellings reach `g.273C>A`. [`DIVERGENT`]
//! is now empty.
//!
//! The pairs were found by deriving each variant's minimal-alignment partition
//! from the resulting sequence, rendering the derived single-block form, and
//! normalizing both spellings. Every pair below was verified to denote the same
//! sequence by an applier independent of the normalizer.
//!
//! Each row asserts the two spellings still DIVERGE. This test failing means a
//! pair converged, and that pair's blessed expectation should be re-blessed to
//! the converged form.
//!
//! Both tables are measured in both directions. Every row above was blessed
//! against the 3' rule, but confluence is a property of the normalizer rather
//! than of one shuffle direction, and `--direction 5prime` is a supported public
//! option — so [`DIVERGENT_UNDER_FIVE_PRIME`] records the [`CONVERGED`] rows that
//! converge under 3' and diverge under 5', and
//! [`the_eight_spelling_pairs_still_diverge_under_five_prime`] pins that the
//! [`DIVERGENT`] eight diverge under 5' too. Both use the same assert-then-flip
//! idiom: a red assertion is a result, not a maintenance chore.

use crate::common::cis_apply_oracle::{apply, normalize, normalize_in};
use crate::common::synthetic::padded;
use ferro_hgvs::ShuffleDirection;

/// `(issue, core, split spelling, merged spelling)`.
///
/// Empty: every harvested pair has converged. Kept — rather than deleted along
/// with its last row — because it is where the next divergence found by
/// `cis_allele_confluence_proptest`'s indel model gets pinned, and because
/// [`every_pinned_pair_denotes_one_variant`] and
/// [`every_divergent_row_still_diverges`] are the contract a new row must
/// satisfy — the latter is dormant while this table is empty and starts
/// checking the moment a row is added, and
/// [`the_divergent_table_is_empty`] pins the count that says so.
const DIVERGENT: &[(&str, &str, &str, &str)] = &[];

/// `(issue, core, split spelling, merged spelling, the string they agree on)`
/// for pairs that converge today.
///
/// The fifth column exists because agreement alone is not a value pin: both
/// spellings regressing together onto the same *wrong* string would satisfy an
/// agreement check silently. Pinning what they agree on is what makes a row
/// evidence rather than a tautology. Note it is often neither of the two input
/// spellings — several settle on repeat notation.
const CONVERGED: &[(&str, &str, &str, &str, &str)] = &[
    // The first six converged when the lone-pure-insertion refusal was removed:
    // each merged form is now derived from the sequence rather than assembled
    // per member, so the split spelling reaches it too.
    (
        "#1287",
        "ATACAGAAAATCAGGGCATA",
        "TEMPLATE:g.[261_262insGA;263_264insAA]",
        "TEMPLATE:g.263_264insGAAA",
        "TEMPLATE:g.263_264insGAAA",
    ),
    (
        "#1290",
        "ATACAGAAAATCAGGGCATA",
        "TEMPLATE:g.[263_264insA;265_266insC]",
        "TEMPLATE:g.266_267insCA",
        "TEMPLATE:g.266_267insCA",
    ),
    (
        "#1301",
        "GCATGAAAAT",
        "TEMPLATE:g.[263_264insAC;264_265insAA]",
        "TEMPLATE:g.264_265insCAAA",
        "TEMPLATE:g.264_265insCAAA",
    ),
    (
        "#1308",
        "CAGAAGATGAATAA",
        "TEMPLATE:g.[263_264insTG;264_265insTG]",
        "TEMPLATE:g.265_266insTTGG",
        "TEMPLATE:g.265_266insTTGG",
    ),
    (
        "#1312",
        "TAAAACCA",
        "TEMPLATE:g.[260_261insAC;261_262insAC]",
        "TEMPLATE:g.262_263insAACC",
        "TEMPLATE:g.262_263insAACC",
    ),
    (
        "#1320",
        "AACAGTAAAATAT",
        "TEMPLATE:g.[263_264insAC;265_266insAA;266_267insAA]",
        "TEMPLATE:g.264_265insCAAAAA",
        "TEMPLATE:g.264_265insCAAAAA",
    ),
    // #1304 converged on top of those six, and needed one more thing they did
    // not: `origin/main`'s removal of the input-separator veto (#1345). Its
    // split spelling is three members whose derived form is a single piece
    // covering a base the input left between two of them, so the veto refused it
    // even with the lone-insertion refusal gone. With both removed the pair
    // reaches the `dup` the sequence states, and the sequence-first pass — not
    // the junction barrier — is what now decides this allele's output. See
    // `issue_1304_junction_barrier_snapshot.rs`, where the barrier keeps its own
    // coverage on a shape the derivation declines to merge.
    (
        "#1304",
        "GCATGAAAAT",
        "TEMPLATE:g.[260_261insGA;261_262insA;264del]",
        "TEMPLATE:g.262_263insGA",
        "TEMPLATE:g.261_262dup",
    ),
    // #1296 is the eighth and last, and the only one that converged by letting
    // the derivation *read* a form rather than by teaching it to build one. Its
    // agreed string is neither input spelling: two members denoting one changed
    // base come back as that substitution. Verified identical under
    // `--direction 5prime`.
    (
        "#1296",
        "AAAAAAATAATCGCAACAGAAG",
        "TEMPLATE:g.[272_273del;274_275insAA]",
        "TEMPLATE:g.273delinsA",
        "TEMPLATE:g.273C>A",
    ),
    // Converged before that change, and unmoved by it.
    (
        "#1286",
        "AAAAAA",
        "TEMPLATE:g.[258_259insA;259_260insA]",
        "TEMPLATE:g.263_264insAA",
        "TEMPLATE:g.257_263A[9]",
    ),
    (
        "#1297",
        "GCATGAAAAT",
        "TEMPLATE:g.[261_262insAA;263del;264_265insA]",
        "TEMPLATE:g.265_266insAA",
        "TEMPLATE:g.262_265A[6]",
    ),
    (
        "#1316",
        "CAGCCAGTCAGCGCATCAG",
        "TEMPLATE:g.[261_262insAA;262_263insAA]",
        "TEMPLATE:g.262_263insAAAA",
        "TEMPLATE:g.262A[5]",
    ),
    (
        "#1321",
        "TCCCAGAAAAT",
        "TEMPLATE:g.[261_262insGA;262_263insA;263del]",
        "TEMPLATE:g.263_264insGA",
        "TEMPLATE:g.262_263dup",
    ),
    (
        "#1323",
        "CAGGGATCAT",
        "TEMPLATE:g.[260del;261_262insGA;262_263insGA]",
        "TEMPLATE:g.262_263insAGA",
        "TEMPLATE:g.262_263insAGA",
    ),
];

/// Both spellings in every row must denote the same sequence, or the row is a
/// broken pin rather than evidence about the normalizer.
#[test]
fn every_pinned_pair_denotes_one_variant() {
    let check = |issue: &str, core: &str, split: &str, merged: &str| {
        let seq = padded(core);
        let a = apply(&seq, split)
            .unwrap_or_else(|| panic!("{issue}: split spelling `{split}` does not apply"));
        let b = apply(&seq, merged)
            .unwrap_or_else(|| panic!("{issue}: merged spelling `{merged}` does not apply"));
        assert_eq!(
            a, b,
            "{issue}: `{split}` and `{merged}` are not the same variant"
        );
    };
    for (issue, core, split, merged) in DIVERGENT {
        check(issue, core, split, merged);
    }
    for (issue, core, split, merged, _) in CONVERGED {
        check(issue, core, split, merged);
    }
}

#[test]
fn the_divergent_table_is_empty() {
    // Named for what actually executes. The predecessor was called
    // `no_pinned_pair_still_diverges`, which promised a divergence check its
    // body cannot perform while `DIVERGENT` is empty, so the name described
    // coverage that was not there. Renamed rather than deleted: the emptiness
    // claim is itself worth pinning. The divergence check it promised lives in
    // [`every_divergent_row_still_diverges`], deliberately apart so this
    // assertion cannot abort it.
    //
    // The count "none" is asserted in three places' prose — this test's name,
    // the module doc above, and `tests/it/rewrite_target_corpus.rs` — and in none
    // of them executably. Adding or removing a row would leave all three wrong
    // and silent. `splitter_reproducer_corpus.rs` guards its own table the same
    // way.
    assert!(
        DIVERGENT.is_empty(),
        "row count changed; update this test's name, the module doc, and \
         tests/it/rewrite_target_corpus.rs's reference to this table"
    );
}

/// Every [`DIVERGENT`] row must still diverge, or it belongs in [`CONVERGED`].
///
/// **Dormant while the table is empty — zero iterations, by construction.**
/// Kept, not deleted, because it is the contract the next row must satisfy: a
/// pair added to `DIVERGENT` must actually diverge when it is added, and must
/// redden here when it converges. Deleting it would drop that check exactly
/// when nobody is looking for it; leaving it labelled costs nothing and stops
/// it reading as live coverage.
///
/// It is a **separate test** from [`the_divergent_table_is_empty`] on purpose.
/// Held together, the emptiness assertion runs first and aborts the whole test
/// on the very change that makes this loop live: adding the first row fails the
/// count pin, and the divergence check below never executes. Split, adding a
/// row reddens the count pin *and* exercises the row.
#[test]
fn every_divergent_row_still_diverges() {
    for (issue, core, split, merged) in DIVERGENT {
        let seq = padded(core);
        let (a, b) = (normalize(&seq, split), normalize(&seq, merged));
        assert_ne!(
            a, b,
            "{issue} appears fixed — `{split}` and `{merged}` now both normalize to \
             `{a}`. Re-bless {issue}'s expectation to the converged form and delete \
             this row."
        );
    }
}

/// #1235's criterion 1 stated directly rather than as a table row: a **lone**
/// spelling of a variant and its multi-member spelling must reach one string.
///
/// This is `#1287`'s pair. It gets its own named test rather than only a row in
/// [`CONVERGED`] because it is the acceptance criterion itself and should not
/// have to be read out of a table: the two spellings must denote one variant
/// (proved below with an applier that is not the normalizer), reach one string,
/// and that string must be a fixed point.
///
/// **Which half moved.** The multi-member spelling is the one that moved onto
/// the lone one, by entering the sequence-first pass as a cis allele. The lone
/// `g.263_264insGAAA` is a fixed point either way, so this pair does *not*
/// exercise `is_splittable_single_member`'s widening from `delins`/`inv` to any
/// edit type — that gate is exercised by
/// `issue_1205_genome_contig_bounds_clamp::a_lone_insertion_written_past_the_contig_end_is_clamped`
/// and its `m.` twin, which are the shapes where a lone member's derived answer
/// differs from its per-member one at all.
#[test]
fn a_lone_insertion_and_its_multi_member_spelling_converge() {
    let seq = padded("ATACAGAAAATCAGGGCATA");
    let lone = "TEMPLATE:g.263_264insGAAA";
    let split = "TEMPLATE:g.[261_262insGA;263_264insAA]";

    // Not assumed: proved with an applier that is not the normalizer, so a
    // convergence onto a *wrong* shared string cannot pass this test.
    assert_eq!(
        apply(&seq, lone).expect("lone spelling applies"),
        apply(&seq, split).expect("split spelling applies"),
        "`{lone}` and `{split}` must denote one variant"
    );

    let from_lone = normalize(&seq, lone);
    let from_split = normalize(&seq, split);
    assert_eq!(
        from_lone, from_split,
        "one variant, two spellings, two normalized strings"
    );
    assert_eq!(from_lone, lone, "the shared answer must be a fixed point");
}

#[test]
fn converged_pairs_stay_converged() {
    for (issue, core, split, merged, expected) in CONVERGED {
        let seq = padded(core);
        let (a, b) = (normalize(&seq, split), normalize(&seq, merged));
        assert_eq!(
            a, b,
            "{issue} regressed — `{split}` -> `{a}` but `{merged}` -> `{b}`; these \
             two spellings of one variant must agree."
        );
        // Agreement alone is also satisfied by both spellings moving together
        // onto a wrong string, so the value they agree on is pinned too.
        assert_eq!(
            &a, expected,
            "{issue} moved — `{split}` and `{merged}` now agree on `{a}`, not `{expected}`"
        );
    }
}

/// Rows from [`CONVERGED`] that converge under the 3' rule but **not** under 5'.
///
/// Confluence is a property of the normalizer, not of one shuffle direction, and
/// `--direction 5prime` is a supported public option on both the CLI
/// (`src/bin/ferro.rs`) and the Python bindings. Every row above was blessed
/// against the 3' direction only, so these gaps were never measured.
///
/// **Both original rows have since left this set, for unrelated reasons.**
///
/// #1321's 5' divergence was the cancelled-member residue: `g.[262_263insA;263del]`
/// merges to `g.263delinsA`, which restates the reference and so renders as
/// `g.263=`, and the split spelling kept that `=` while the merged spelling never
/// grew one. Dropping the residue where the merge creates it converges the pair on
/// `g.261_262dup` — see `issue_1321_identity_inside_a_duplication.rs`.
///
/// #1323 closed here: lowering a repeat member to the `ins`/`del` it denotes let
/// `collect_canonical_edits` read its 5' spelling, which is the gate that had been
/// refusing it; both spellings now reach one string under 5' as they already did
/// under 3'.
///
/// The set is kept rather than deleted with its last row — it is where the next
/// direction-dependent gap gets pinned, and `the_five_prime_confluence_gap_is_unchanged`
/// asserts every unlisted row converges, so an empty set is a stronger claim than a
/// populated one, not a dormant check.
///
/// What a row here means, when there is one: both spellings still denote the
/// input's bases and both are stable fixed points, so it is criterion 1 of #1235
/// — non-confluence — and nothing else. That is why neither oracle sees such a
/// row: `FERRO_ASSERT_IDEMPOTENT` re-normalizes a single spelling and finds it
/// stable, and `FERRO_ASSERT_REPARSE` finds it well-formed. Only comparing two
/// spellings of one variant exposes it.
const DIVERGENT_UNDER_FIVE_PRIME: &[&str] = &[];

#[test]
fn the_five_prime_confluence_gap_is_unchanged() {
    // Assert-then-flip, the same contract the rows above carry: each listed row
    // is asserted to still DIVERGE, and each unlisted row to still CONVERGE. A
    // red assertion here is a *result*, not a maintenance chore — move the row
    // between the two sets rather than widening the list.
    //
    // The cause is not the shuffle direction reaching the derivation, which was
    // measured and refuted: threading `ShuffleDirection` into
    // `canonicalize_from_sequence` changes 298 of 39,600 swept outputs and all of
    // them are sequence-preserving both before and after, so it is a lateral
    // re-spelling with no correctness content. For the one row still listed the
    // derivation instead *refuses*: `canonicalize_from_sequence` returns `None`,
    // behind two stacked gates — the `collect_canonical_edits` catch-all
    // (reached via a redundant `=` member the 5' pipeline emits), then the
    // `changed_columns_of_pieces` weight bound. Each masks the next, which is
    // the pattern #1235 and #1345 both describe. A third gate stood behind those
    // two when this was written — `needs_unsupported_form`, since removed — so
    // the ladder is two deep today, not three.
    //
    // Every listed id must still name a `CONVERGED` row. The loop below reaches
    // an entry only through that table, so an id naming no row is inert: it
    // asserts nothing, and the set goes on claiming to pin a divergence that is
    // no longer measured. Deleting a row and leaving its id here is exactly how
    // that happens, and it is silent — verified by adding a bogus id, which left
    // all six tests in this file green. Same class of drift as the row count
    // pinned in `the_eight_spelling_pairs_still_diverge`.
    for issue in DIVERGENT_UNDER_FIVE_PRIME {
        assert!(
            CONVERGED.iter().any(|(row, ..)| row == issue),
            "{issue} is listed in DIVERGENT_UNDER_FIVE_PRIME but names no CONVERGED \
             row, so it pins nothing. Delete the entry, or restore the row it names."
        );
    }
    for (issue, core, split, merged, _three_prime_expected) in CONVERGED {
        let seq = padded(core);
        let (a, b) = (
            normalize_in(&seq, split, ShuffleDirection::FivePrime),
            normalize_in(&seq, merged, ShuffleDirection::FivePrime),
        );
        if DIVERGENT_UNDER_FIVE_PRIME.contains(issue) {
            assert_ne!(
                a, b,
                "{issue} now CONVERGES under 5' shuffle — both spellings reach \
                 `{a}`. Remove it from DIVERGENT_UNDER_FIVE_PRIME; the gap is \
                 closing."
            );
        } else {
            assert_eq!(
                a, b,
                "{issue} has started diverging under 5' shuffle — `{split}` -> \
                 `{a}` but `{merged}` -> `{b}`. Two spellings of one variant must \
                 agree in both directions; add it to DIVERGENT_UNDER_FIVE_PRIME \
                 only with a measured reason."
            );
        }
    }
}

#[test]
fn the_eight_spelling_pairs_still_diverge_under_five_prime() {
    // The 5' half of `the_eight_spelling_pairs_still_diverge`. Those eight rows
    // were blessed against the 3' rule only, so without this the DIVERGENT table
    // had no 5' coverage at all and the module doc's "both directions" claim
    // covered only CONVERGED.
    //
    // All eight diverge under 5' as well, but they reach *different* spellings
    // than under 3' (#1296's merged form settles at `g.273C>A` here against
    // `g.273delinsA` there), so this is measuring the 5' pipeline rather than
    // restating the 3' result. Assert-then-flip: a row that starts converging is
    // the gap closing, and belongs in the message below rather than deleted.
    for (issue, core, split, merged) in DIVERGENT {
        let seq = padded(core);
        let (a, b) = (
            normalize_in(&seq, split, ShuffleDirection::FivePrime),
            normalize_in(&seq, merged, ShuffleDirection::FivePrime),
        );
        assert_ne!(
            a, b,
            "{issue} now CONVERGES under 5' shuffle — `{split}` and `{merged}` both \
             reach `{a}`. The gap is closing; re-bless {issue} and move it to \
             CONVERGED if it converges under 3' too."
        );
    }
}

#[test]
fn the_five_prime_divergences_are_non_confluence_and_nothing_worse() {
    // The severity claim in DIVERGENT_UNDER_FIVE_PRIME's doc, made executable. If
    // one of these ever starts changing the denoted sequence or ceases to be a
    // fixed point, that is a different and much worse defect than the one
    // recorded here, and it must not hide behind a known non-confluence.
    //
    // Every output is now checked directly, and two changes were needed to get
    // there. The one output that could not be put to the applier was `#1321`'s,
    // which carried a cancelled-member `=` residue: #1362 gave an identity's SPDI
    // triple the span it claims, so it no longer collides with a sibling's
    // insertion junction and such an output is expressible at all; #1379 then
    // stopped the residue being emitted, so this corpus no longer produces one.
    // The identity-stripping rewrite and the tolerated-skip count this loop used
    // to carry are both gone, and an inexpressible output is now a failure rather
    // than a known exception.
    for (issue, core, split, merged, _three_prime_expected) in CONVERGED {
        if !DIVERGENT_UNDER_FIVE_PRIME.contains(issue) {
            continue;
        }
        let seq = padded(core);
        let want = apply(&seq, split).expect("input applies");
        for input in [split, merged] {
            let out = normalize_in(&seq, input, ShuffleDirection::FivePrime);
            assert_eq!(
                normalize_in(&seq, &out, ShuffleDirection::FivePrime),
                out,
                "{issue}: `{out}` is no longer a fixed point"
            );
            let got = apply(&seq, &out).unwrap_or_else(|| {
                panic!(
                    "{issue}: `{out}` is inexpressible to the applier. Since #1362 \
                     every output here is expressible, and since #1379 none carries \
                     a cancelled-member `=` residue at all, so this is a new blind \
                     spot rather than a known one — find what shape the applier \
                     declines before treating it as tolerable."
                )
            });
            assert_eq!(
                got, want,
                "{issue}: `{input}` -> `{out}` no longer denotes the input's bases"
            );
        }
    }
}
