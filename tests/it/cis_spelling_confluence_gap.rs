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
//! ## Four of those convergences have since been given back, deliberately
//!
//! `partition-is-the-unit-of-normalization` (DECIDED, 2026-08-08,
//! `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`) says the unit
//! of normalization is the partition the input asserts, not the pair of sequences
//! it denotes — and every convergence above was reached by re-deriving a
//! partition from the resulting sequence. So `#1287`, `#1301`, `#1320` and
//! `#1304` no longer converge under 3', and they have moved to
//! [`PARTITION_DIVERGENT`], **not** back to [`DIVERGENT`]: a row there is a
//! defect awaiting a fix, and these are the ruling working. Their split
//! spellings' members are one or more unchanged reference nucleotides apart,
//! which `general.md:34` says are "described individually and **not** as a
//! \"delins\"".
//!
//! Nine rows still converge, and four more ([`DIVERGENT_UNDER_FIVE_PRIME`])
//! converge under 3' but not under 5'. Both tables' row counts are pinned by
//! [`the_two_tables_have_the_shapes_the_docs_claim`], because these paragraphs
//! cannot fail.
//!
//! This is a representation change for anyone who stored a merged form, and it is
//! written down as such rather than absorbed. What every affected pair keeps is
//! `EquivalenceLevel::SequenceMatch`, which still answers "are these the same
//! variant"; only a consumer keying on the normalized string sees a difference.
//!
//! The pairs were found by deriving each variant's minimal-alignment partition
//! from the resulting sequence, rendering the derived single-block form, and
//! normalizing both spellings. Every pair below was verified to denote the same
//! sequence by an applier independent of the normalizer.
//!
//! Each [`DIVERGENT`] row asserts the two spellings still DIVERGE. That test
//! failing means a pair converged, and that pair's blessed expectation should be
//! re-blessed to the converged form. [`PARTITION_DIVERGENT`] is the other way
//! round — its rows are *meant* to diverge — which is why the two are separate
//! tables rather than one with a flag.
//!
//! All three tables are measured in both directions. Every row above was blessed
//! against the 3' rule, but confluence is a property of the normalizer rather
//! than of one shuffle direction, and `--direction 5prime` is a supported public
//! option — so [`DIVERGENT_UNDER_FIVE_PRIME`] records the [`CONVERGED`] rows that
//! converge under 3' and diverge under 5', and
//! [`every_divergent_row_still_diverges_under_five_prime`] holds the same
//! contract for [`DIVERGENT`] — **dormant, since that table is empty**, exactly
//! as its 3' sibling is — while
//! [`every_partition_divergent_row_reaches_its_pinned_forms`] pins all four of
//! its rows' strings in both directions, sixteen strings in total. All use the
//! same assert-then-flip idiom: a red assertion is a result, not a maintenance
//! chore.

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

/// One [`PARTITION_DIVERGENT`] row:
/// `(issue, core, split, merged, split@3', merged@3', split@5', merged@5')`.
///
/// Named because the tuple is eight wide and `clippy::type_complexity` is right
/// that the bare form is unreadable. The width is deliberate — see the table's
/// own doc for why every one of the four outputs is pinned separately.
type PartitionDivergentRow = (
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
);

/// `(issue, core, split, merged, split@3', merged@3', split@5', merged@5')` for
/// pairs the **partition ruling** deliberately keeps apart.
///
/// # Why this is a separate table from [`DIVERGENT`], and not a reopening of it
///
/// [`DIVERGENT`] is the *defect* table: a row there is a pair that ought to
/// converge and does not, and its contract
/// ([`every_divergent_row_still_diverges`]) is assert-then-flip — a row going
/// green is a fix. These four rows are the opposite. Under
/// `partition-is-the-unit-of-normalization` (DECIDED, 2026-08-08,
/// `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`) the two
/// spellings of each are **two different assertions about how many variants
/// there are**, not two spellings of one question, so their reaching two
/// canonical strings is the model working. The ruling says so directly: "Two
/// spellings that assert different partitions of the same bases reach different
/// canonical forms by design."
///
/// In every row the split spelling's members are separated by one or more
/// unchanged reference nucleotides, which `general.md:34` says are "described
/// individually and **not** as a \"delins\"", while the merged spelling asserts
/// one block. Each row's previously-blessed agreement was reached by discarding
/// both assertions and re-deriving the partition from the resulting sequence —
/// exactly the move the ruling removes.
///
/// # Both strings are pinned, in both directions
///
/// Four columns rather than one shared answer, because agreement is no longer
/// the property. This is strictly stronger than the [`CONVERGED`] contract these
/// rows came from, which pinned one string per row: a move on either side, in
/// either direction, reddens here.
///
/// The 3'/5' split is not decoration. **Three of the four still converge under
/// 5'**, because 5'-shifting walks the two payloads onto one junction, where
/// separation is zero and the ruling's own MERGE move applies. So the partition
/// is preserved in both directions; what differs is where canonical placement
/// puts the members before the floor is measured. That is the reading campaign
/// issue #65 records as still open (D2) — measured here, not adjudicated here.
///
/// # Every row is non-confluence and nothing worse
///
/// [`every_pinned_pair_denotes_one_variant`] proves the two spellings denote one
/// variant with `cis_apply_oracle::apply` (a `hgvs_to_spdi` walk the normalizer
/// does not consult), and
/// [`every_partition_divergent_row_reaches_its_pinned_forms`] proves each of the
/// sixteen outputs still denotes those bases and is a fixed point — four rows
/// times the four strings each row pins (`split`/`merged` in each of the two
/// shuffle directions), which is the count that test's own doc states. So a consumer
/// asking "are these the same variant" is still answered — by
/// `EquivalenceLevel::SequenceMatch` — and only a consumer keying on the string
/// sees a difference.
const PARTITION_DIVERGENT: &[PartitionDivergentRow] = &[
    // Separation two: the reference bases at 262 and 263 sit between the two
    // insertions. This row is #1235's acceptance criterion 1 stated as a table
    // entry; the named test
    // [`a_lone_insertion_and_its_multi_member_spelling_assert_different_partitions`]
    // states it directly.
    (
        "#1287",
        "ATACAGAAAATCAGGGCATA",
        "TEMPLATE:g.[261_262insGA;263_264insAA]",
        "TEMPLATE:g.263_264insGAAA",
        "TEMPLATE:g.[262_263dup;265_266dup]",
        "TEMPLATE:g.263_264insGAAA",
        "TEMPLATE:g.[261_262dup;263_266A[6]]",
        "TEMPLATE:g.262_263insAGAA",
    ),
    // Separation one: the reference base at 264. Diverges at 3'; at 5' both
    // payloads shift onto one junction and the MERGE move closes the gap.
    (
        "#1301",
        "GCATGAAAAT",
        "TEMPLATE:g.[263_264insAC;264_265insAA]",
        "TEMPLATE:g.264_265insCAAA",
        "TEMPLATE:g.[264_265insCA;264_265dup]",
        "TEMPLATE:g.264_265insCAAA",
        "TEMPLATE:g.261_262insAAAC",
        "TEMPLATE:g.261_262insAAAC",
    ),
    // Three authored members. The first stays its own block (separation two);
    // the second and third land on one junction at 3' and merge, so the split
    // spelling answers with two members rather than three.
    (
        "#1320",
        "AACAGTAAAATAT",
        "TEMPLATE:g.[263_264insAC;265_266insAA;266_267insAA]",
        "TEMPLATE:g.264_265insCAAAAA",
        "TEMPLATE:g.[264_265insCA;266_267insAAAA]",
        "TEMPLATE:g.264_265insCAAAAA",
        "TEMPLATE:g.262_263insAACAAA",
        "TEMPLATE:g.262_263insAACAAA",
    ),
    // **CHARACTERIZATION PIN, NOT AN ADJUDICATION — read this row differently
    // from the three above.** Its 3' split answer carries a redundant identity
    // member, `265=`: the `264del` cancels against part of the insertions and
    // the residue is emitted rather than dropped. This module's own
    // [`DIVERGENT_UNDER_FIVE_PRIME`] doc records that #1379 stopped such a
    // residue being emitted at all, so its reappearance is a suspected
    // regression of that fix, **not** a consequence of the partition ruling.
    //
    // It is pinned exactly rather than softened so that the fix reddens this row
    // instead of passing silently, and it is deliberately not blessed as
    // correct. The same shape is under separate adjudication as
    // `issue_1304_junction_barrier_snapshot::a_moved_sibling_shape_now_merges_from_the_sequence`;
    // when that lands, re-measure this row and move it to [`CONVERGED`] if the
    // residue's removal converges the pair (it would: the merged spelling
    // already answers `g.261_262dup`, and the two agree under 5').
    (
        "#1304",
        "GCATGAAAAT",
        "TEMPLATE:g.[260_261insGA;261_262insA;264del]",
        "TEMPLATE:g.262_263insGA",
        "TEMPLATE:g.[261_262dup;265=]",
        "TEMPLATE:g.261_262dup",
        "TEMPLATE:g.261_262dup",
        "TEMPLATE:g.261_262dup",
    ),
];

/// `(issue, core, split spelling, merged spelling, the string they agree on)`
/// for pairs that converge today.
///
/// The fifth column exists because agreement alone is not a value pin: both
/// spellings regressing together onto the same *wrong* string would satisfy an
/// agreement check silently. Pinning what they agree on is what makes a row
/// evidence rather than a tautology. Note it is often neither of the two input
/// spellings — several settle on repeat notation.
const CONVERGED: &[(&str, &str, &str, &str, &str)] = &[
    (
        "#1290",
        "ATACAGAAAATCAGGGCATA",
        "TEMPLATE:g.[263_264insA;265_266insC]",
        "TEMPLATE:g.266_267insCA",
        "TEMPLATE:g.266_267insCA",
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
    for (issue, core, split, merged, ..) in PARTITION_DIVERGENT {
        check(issue, core, split, merged);
    }
}

/// Every [`PARTITION_DIVERGENT`] row reaches its four pinned strings, each of
/// which denotes the input's bases and is a fixed point.
///
/// This is the guard the four rows carry in place of the convergence check they
/// left behind, and it is stronger than that check in three ways: it pins two
/// strings per row instead of one, it does so in both shuffle directions instead
/// of only 3', and it proves the bases of every output through an applier the
/// normalizer does not consult instead of inferring them from agreement.
///
/// It deliberately does **not** assert that the two spellings disagree. Three of
/// the four rows agree under 5', so a blanket `assert_ne` would be false; the
/// per-direction string pins say exactly which rows agree where, which is more
/// information than a boolean and cannot rot into a claim nobody re-measured.
#[test]
fn every_partition_divergent_row_reaches_its_pinned_forms() {
    for (issue, core, split, merged, split3, merged3, split5, merged5) in PARTITION_DIVERGENT {
        let seq = padded(core);
        let want = apply(&seq, split)
            .unwrap_or_else(|| panic!("{issue}: split spelling `{split}` does not apply"));
        for (direction, input, expected) in [
            (ShuffleDirection::ThreePrime, split, split3),
            (ShuffleDirection::ThreePrime, merged, merged3),
            (ShuffleDirection::FivePrime, split, split5),
            (ShuffleDirection::FivePrime, merged, merged5),
        ] {
            let got = normalize_in(&seq, input, direction);
            assert_eq!(
                &got, expected,
                "{issue}: `{input}` moved under {direction:?}. Both spellings are \
                 pinned in both directions; a move is a representation change \
                 whichever way it goes."
            );
            assert_eq!(
                normalize_in(&seq, &got, direction),
                got,
                "{issue}: `{got}` is not a fixed point under {direction:?}"
            );
            let denoted = apply(&seq, &got).unwrap_or_else(|| {
                panic!(
                    "{issue}: `{input}` -> `{got}` denotes no sequence under \
                     {direction:?} — its members overlap or claim one interbase \
                     twice, which is a defect and not a partition difference"
                )
            });
            assert_eq!(
                denoted, want,
                "{issue}: `{input}` -> `{got}` no longer denotes the input's bases \
                 under {direction:?}"
            );
        }
    }
}

/// The row counts of both tables, pinned so neither can move in silence.
///
/// Nine rows converge and four are kept apart by the partition ruling. The
/// module doc and [`DIVERGENT_UNDER_FIVE_PRIME`]'s doc both quote those numbers
/// in prose, and prose does not fail.
#[test]
fn the_two_tables_have_the_shapes_the_docs_claim() {
    assert_eq!(
        CONVERGED.len(),
        9,
        "row count changed; update the module doc and this test"
    );
    assert_eq!(
        PARTITION_DIVERGENT.len(),
        4,
        "row count changed; update the module doc and this test"
    );
    for (issue, ..) in PARTITION_DIVERGENT {
        assert!(
            !CONVERGED.iter().any(|(row, ..)| row == issue),
            "{issue} is in both tables; it can only be in one"
        );
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

/// #1235's criterion 1 stated directly rather than as a table row — and the one
/// place in this file where the partition ruling's cost is written down as a
/// **loss**, because that is what it is.
///
/// This is `#1287`'s pair. It gets its own named test rather than only a row in
/// [`PARTITION_DIVERGENT`] because it is the acceptance criterion itself and
/// should not have to be read out of a table.
///
/// RE-BLESSED under `partition-is-the-unit-of-normalization` (DECIDED,
/// 2026-08-08, `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`),
/// and renamed, because the old name asserted the convergence that is gone:
///
/// ```text
/// lone   TEMPLATE:g.263_264insGAAA              -> itself                       (unmoved)
/// split  TEMPLATE:g.[261_262insGA;263_264insAA] -> TEMPLATE:g.[262_263dup;265_266dup]
///                                       was      -> TEMPLATE:g.263_264insGAAA
/// ```
///
/// **State the cost plainly: #1235's criterion 1 does not hold for this pair any
/// more, and the ruling says it is not meant to.** Its own text: "This does not
/// make ferro confluent with respect to sequence equivalence, and it is not
/// meant to. Two spellings that assert different partitions of the same bases
/// reach different canonical forms by design." The split spelling's two members
/// sit at the junctions `261|262` and `263|264` with the reference bases at 262
/// and 263 unchanged between them — separation two, which `general.md:34`
/// describes individually — while the lone spelling asserts one block. The
/// previous agreement was reached by discarding both assertions and re-deriving
/// the partition from the resulting sequence.
///
/// **Which half moved: the multi-member one.** The lone `g.263_264insGAAA` is a
/// fixed point either way and is byte-identical to what it has always been, so a
/// consumer holding the lone spelling's answer is unaffected; a consumer holding
/// the split spelling's is looking at a migration.
///
/// What the pair still has, and what a consumer should use instead, is
/// `EquivalenceLevel::SequenceMatch` — asserted below on the two *outputs*, not
/// only on the inputs, so this test says the divergence is a representation
/// difference and nothing worse.
#[test]
fn a_lone_insertion_and_its_multi_member_spelling_assert_different_partitions() {
    use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};
    use ferro_hgvs::parse_hgvs;

    let seq = padded("ATACAGAAAATCAGGGCATA");
    let lone = "TEMPLATE:g.263_264insGAAA";
    let split = "TEMPLATE:g.[261_262insGA;263_264insAA]";

    // Not assumed: proved with an applier that is not the normalizer, so the
    // divergence below is a representation difference rather than the two
    // spellings having been different variants all along.
    let want = apply(&seq, lone).expect("lone spelling applies");
    assert_eq!(
        want,
        apply(&seq, split).expect("split spelling applies"),
        "`{lone}` and `{split}` must denote one variant"
    );

    let from_lone = normalize(&seq, lone);
    let from_split = normalize(&seq, split);

    // Both strings pinned exactly. Pinning only the divergence would be
    // satisfied by either side regressing anywhere at all.
    assert_eq!(
        from_lone, lone,
        "the lone spelling is unmoved and must stay a fixed point"
    );
    assert_eq!(
        from_split, "TEMPLATE:g.[262_263dup;265_266dup]",
        "the multi-member spelling keeps both members (separation 2, general.md:34)"
    );
    assert_ne!(
        from_lone, from_split,
        "the two spellings assert different partitions and must stay distinct; \
         re-converging them in either direction is a representation change that \
         needs adjudicating, not a silent improvement"
    );

    // Neither answer may drift off the bases or off its own fixed point — a
    // known non-confluence must not become cover for something worse.
    for form in [&from_lone, &from_split] {
        assert_eq!(
            normalize(&seq, form),
            *form,
            "`{form}` is not a fixed point"
        );
        assert_eq!(
            apply(&seq, form).unwrap_or_else(|| panic!("`{form}` denotes no sequence")),
            want,
            "`{form}` no longer denotes the input's bases"
        );
    }

    // …and the property #1235 is really about is still available to a consumer,
    // on the rung built for it. `EquivalenceChecker` reaches `SequenceMatch` only
    // after `Identical` and `NormalizedMatch` have both declined, so this
    // simultaneously re-states the divergence and shows it is harmless to a
    // consumer that asks the right question.
    let provider = crate::common::cis_apply_oracle::provider(&seq);
    let verdict = EquivalenceChecker::new(provider)
        .check(
            &parse_hgvs(&from_lone).expect("parse the normalized lone spelling"),
            &parse_hgvs(&from_split).expect("parse the normalized multi-member spelling"),
        )
        .expect("the two normalized forms must be comparable");
    assert_eq!(
        verdict.level,
        EquivalenceLevel::SequenceMatch,
        "the two canonical forms must still be recognisable as one variant"
    );
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
/// **Four rows entered this set with the partition ruling** — `#1290`, `#1308`,
/// `#1312` and `#1323`. All four converge under 3' and diverge under 5', for one
/// reason: under 3' the split spelling's payloads shift onto the same junction
/// as the merged form's single block, so the two agree, while under 5' they walk
/// the other way and the split spelling keeps the members `general.md:34` says
/// to describe individually. That is `partition-is-the-unit-of-normalization`
/// (DECIDED, 2026-08-08,
/// `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`) working, not
/// a defect — but it is recorded here rather than absorbed, because it is a
/// direction-dependent representation difference and
/// [`the_five_prime_divergences_are_non_confluence_and_nothing_worse`] is what
/// keeps it from becoming cover for something worse. Measured: all eight
/// outputs denote their inputs' bases and all eight are fixed points.
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
/// `the_five_prime_confluence_gap_is_unchanged` asserts every unlisted row
/// converges, so the set is a two-sided ledger rather than a suppression list:
/// adding a row is as loud as removing one.
///
/// What a row here means, when there is one: both spellings still denote the
/// input's bases and both are stable fixed points, so it is criterion 1 of #1235
/// — non-confluence — and nothing else. That is why neither oracle sees such a
/// row: `FERRO_ASSERT_IDEMPOTENT` re-normalizes a single spelling and finds it
/// stable, and `FERRO_ASSERT_REPARSE` finds it well-formed. Only comparing two
/// spellings of one variant exposes it.
const DIVERGENT_UNDER_FIVE_PRIME: &[&str] = &["#1290", "#1308", "#1312", "#1323"];

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
    // pinned in `the_divergent_table_is_empty`.
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

/// The 5' half of [`every_divergent_row_still_diverges`], and **dormant on the
/// same terms** — zero iterations while [`DIVERGENT`] is empty.
///
/// Renamed from `the_eight_spelling_pairs_still_diverge_under_five_prime`, whose
/// name and body both promised coverage of eight rows that are no longer there.
/// The 3' sibling above was labelled when the table emptied and this one was
/// missed, which is how a name outlived its table.
///
/// Kept for the same reason the sibling is: it is the contract the next row must
/// satisfy in the 5' direction. When the rows existed they diverged under 5' too,
/// but reached *different* spellings than under 3' (#1296's merged form settled
/// at `g.273C>A` here against `g.273delinsA` there), so this measures the 5'
/// pipeline rather than restating the 3' result — which is why it stays a
/// separate test rather than being folded into the sibling.
#[test]
fn every_divergent_row_still_diverges_under_five_prime() {
    // Assert-then-flip: a row that starts converging is the gap closing, and
    // belongs in the message below rather than deleted.
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
