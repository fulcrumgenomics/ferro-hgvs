//! Pins `rulings[unchanged-is-read-over-every-minimal-alignment]`'s own data
//! against the enumerator that implements it.
//!
//! The record decides that **a reference base is unchanged iff it is matched in
//! every minimal alignment of the block**, and argues that decision through
//! four worked cases. Until `common::minimal_alignment` existed, those cases
//! were arithmetic nobody could re-run: the rule quantifies over *all* minimal
//! alignments and nothing enumerated them. Every case the record states is
//! reproduced below, so a change to the instrument that stopped agreeing with
//! the ruling would go red rather than quietly redefining "unchanged".
//!
//! # This module changes no behaviour
//!
//! Nothing here normalizes, and nothing under `src/` reads the enumerator. The
//! record's rule is *applied* by `src/normalize/`; what was missing was a way
//! to **check** an application of it. These tests are that check and nothing
//! more.
//!
//! # The record does not say which cost model governs
//!
//! It is inferable from the record's own worked example — `CAG -> AGA`,
//! "position-wise: cost 3" over three substitutions, so a substitution costs 1
//! — but it is never written down, and the two plausible models disagree on a
//! real reported block. [`the_two_cost_models_disagree_on_1420_v4`] is the
//! record of that gap. It picks no winner, deliberately: which model governs is
//! an open adjudication, and pinning both answers is what keeps the question
//! answerable later instead of being settled by whichever model an instrument
//! happened to hardcode.

use crate::common::minimal_alignment::{
    minimal_alignments, minimal_alignments_capped, position_wise_matches, CostModel, Step,
    DEFAULT_ALIGNMENT_CAP,
};

/// Assert cost, alignment count and unchanged set at once, so a failure names
/// all three rather than stopping at the first.
#[track_caller]
fn assert_block(
    reference: &[u8],
    alternate: &[u8],
    model: CostModel,
    expected_cost: u32,
    expected_count: usize,
    expected_unchanged: &[u32],
) {
    let report = minimal_alignments(reference, alternate, model).unwrap_or_else(|exceeded| {
        panic!(
            "{} -> {} under {model:?}: {exceeded}",
            String::from_utf8_lossy(reference),
            String::from_utf8_lossy(alternate),
        )
    });
    assert_eq!(
        (report.cost(), report.count(), report.unchanged()),
        (expected_cost, expected_count, expected_unchanged),
        "{} -> {} under {model:?}: (cost, alignments, unchanged)",
        String::from_utf8_lossy(reference),
        String::from_utf8_lossy(alternate),
    );
}

/// **The record's headline case.** `CAG -> AGA` is equal length, and equal
/// length does not imply no indel.
///
/// The record states it directly:
///
/// ```text
/// position-wise:            C->A, A->G, G->A            cost 3
/// del C, A=A, G=G, ins A                                cost 2   <- minimal
/// ```
///
/// So reference offsets 1 and 2 are unchanged, the two changes *are* separated
/// by them, and `general.md:34` describes them individually. The position-wise
/// reading, which matches nothing at all here, is the one the record says has no
/// authority behind it.
///
/// Identical under both cost models, which is why the record could argue from
/// it without stating a model.
#[test]
fn cag_to_aga_is_a_del_plus_ins_around_two_unchanged_bases() {
    for model in CostModel::ALL {
        assert_block(b"CAG", b"AGA", model, 2, 1, &[1, 2]);
    }

    assert_eq!(
        position_wise_matches(b"CAG", b"AGA"),
        Vec::<u32>::new(),
        "the position-wise reading matches nothing, which is the whole point"
    );

    let report = minimal_alignments(b"CAG", b"AGA", CostModel::Levenshtein).unwrap();
    assert!(
        report.is_unique(),
        "one minimal alignment, so nothing to agree on"
    );
    assert_eq!(
        report.alignments()[0],
        vec![Step::Del, Step::Match, Step::Match, Step::Ins],
        "the record spells this alignment out: del C, A=A, G=G, ins A"
    );
}

/// **The case that forced the column/cell distinction.** `GACA -> AGAT`.
///
/// The record's argument: the minimal alignments "both match ref offset 1, but
/// to different alt offsets, so the cell-based notion sees nothing while the
/// base is unchanged either way." The intersection here is exactly one
/// reference offset, and no single match *edge* is common to all the
/// alignments — which is why `Dominators::matched_ref()` (cell-based) and the
/// record's clause (column-based) give different answers on it, and why the
/// record insists that is a notion mismatch rather than a bug.
///
/// # The record says "the two cost-3 alignments"; there are three
///
/// Measured, not inferred: under substitution cost 1 this block has **three**
/// distinct minimal alignments, whose matched-reference sets are `{1, 3}`,
/// `{0, 1}` and `{0, 1}`. Two of the three agree on their matched set, which is
/// the most likely reading behind the record's "two" — it is a count of
/// distinct outcomes, not of alignments.
///
/// The discrepancy is pinned rather than resolved because **the record's
/// conclusion is unaffected**: reference offset 1 is matched in all three, so
/// the intersection is `{1}` on either count. This test records the arithmetic
/// so nobody re-derives it a fourth time.
#[test]
fn gaca_to_agat_intersects_to_one_base_across_three_alignments() {
    assert_block(b"GACA", b"AGAT", CostModel::Levenshtein, 3, 3, &[1]);

    let report = minimal_alignments(b"GACA", b"AGAT", CostModel::Levenshtein).unwrap();
    let matched: Vec<Vec<u32>> = (0..report.count())
        .map(|i| report.matched_reference_offsets(i))
        .collect();
    assert_eq!(
        matched,
        vec![vec![1, 3], vec![0, 1], vec![0, 1]],
        "three alignments, two distinct matched-reference outcomes"
    );
    assert!(
        !report.is_unique(),
        "this is the case where the unchanged set is a genuine intersection"
    );

    // Under the other model the block costs 4 and agrees on nothing at all.
    assert_block(b"GACA", b"AGAT", CostModel::SubstitutionCostsTwo, 4, 9, &[]);
}

/// `1420-v2`: the `TEMPLATE` block at 38-41, `TTGC -> ATTG`.
///
/// The `g.[37dup;41del]` / `g.38_41delinsATTG` pair in
/// `reported_confluence_pairs::REPORTED_PAIRS`. Exactly one optimal alignment,
/// so the unchanged set is pinned rather than intersected, and the two cost
/// models agree — this block says nothing about the open model question.
///
/// Reference offsets 0-2 are unchanged, i.e. `g.38`, `g.39` and `g.40`; the
/// change is an insertion before the block and a deletion of `g.41`.
#[test]
fn template_38_41_is_one_alignment_over_three_unchanged_bases() {
    for model in CostModel::ALL {
        assert_block(b"TTGC", b"ATTG", model, 2, 1, &[0, 1, 2]);
    }

    let report = minimal_alignments(b"TTGC", b"ATTG", CostModel::Levenshtein).unwrap();
    assert_eq!(
        report.alignments()[0],
        vec![Step::Ins, Step::Match, Step::Match, Step::Match, Step::Del],
        "insert A, then TTG matches, then delete C"
    );
    assert_eq!(
        position_wise_matches(b"TTGC", b"ATTG"),
        vec![1],
        "the position-wise reading finds one agreeing base, while three of the \
         four are unchanged"
    );
}

/// `1420-v3`: the `TEMPLATE` block at 37-40, `ATTG -> CATT`.
///
/// The same shape as `1420-v2` one base upstream — an insertion before the
/// block and a deletion at its end — and like it, unique and model-independent.
/// Kept as its own case because it is a separate reported row, not because it
/// exercises anything `1420-v2` does not.
#[test]
fn template_37_40_is_one_alignment_over_three_unchanged_bases() {
    for model in CostModel::ALL {
        assert_block(b"ATTG", b"CATT", model, 2, 1, &[0, 1, 2]);
    }

    let report = minimal_alignments(b"ATTG", b"CATT", CostModel::Levenshtein).unwrap();
    assert_eq!(
        report.alignments()[0],
        vec![Step::Ins, Step::Match, Step::Match, Step::Match, Step::Del],
    );
}

/// **`1420-v4` is where the two cost models disagree, and the record does not
/// say which one governs.**
///
/// The `TEMPLATE` block at 21-24, `ATGC -> GCTG`
/// (`g.[21delinsGC;24del]` / `g.21_24delinsGCTG`).
///
/// | model | cost | alignments | unchanged (block) | unchanged (genomic) |
/// |---|---:|---:|---|---|
/// | substitution costs 1 | 3 | 2 | `{1, 2}` | `g.22`, `g.23` |
/// | substitution costs 2 | 4 | 6 | `{2}` | `g.23` |
///
/// Under one model `g.22` is an unchanged base and the changes either side of
/// it are "separated by one or more nucleotides", so `general.md:34` describes
/// them individually. Under the other it is not, and they are consecutive.
/// **That is a difference in what the spec clause says about this block, from a
/// parameter the record never fixes.**
///
/// # Nothing here picks a winner
///
/// The record's own worked example implies substitution cost 1, but implication
/// is not a ruling, and `rulings[adjudication-precedence-order]` is the register
/// for escalating a residue like this — not a test module. So both answers are
/// pinned, in full, and this doc comment is the record that the question is
/// open. If it is later decided, this test is where the decision lands and what
/// it has to be reconciled against.
#[test]
fn the_two_cost_models_disagree_on_1420_v4() {
    assert_block(b"ATGC", b"GCTG", CostModel::Levenshtein, 3, 2, &[1, 2]);
    assert_block(
        b"ATGC",
        b"GCTG",
        CostModel::SubstitutionCostsTwo,
        4,
        6,
        &[2],
    );

    let unit = minimal_alignments(b"ATGC", b"GCTG", CostModel::Levenshtein).unwrap();
    let paired = minimal_alignments(b"ATGC", b"GCTG", CostModel::SubstitutionCostsTwo).unwrap();
    assert_ne!(
        unit.unchanged(),
        paired.unchanged(),
        "if these ever agree the disagreement this test exists to record has \
         been resolved by a change to the instrument, which is the one way it \
         must not be resolved"
    );

    // And the position-wise reading agrees with neither, matching nothing at
    // all — so the block has three readings and no two of them coincide.
    assert_eq!(position_wise_matches(b"ATGC", b"GCTG"), Vec::<u32>::new());
}

/// Identical sequences: zero cost, one alignment, every base unchanged.
///
/// The degenerate end of the rule, and the one shape where the position-wise
/// reading and the minimal-alignment reading are *guaranteed* to coincide
/// rather than merely happening to. Elsewhere the two are incomparable — see
/// [`neither_notion_contains_the_other`].
#[test]
fn identical_sequences_leave_every_base_unchanged() {
    for model in CostModel::ALL {
        assert_block(
            b"ACGTACGT",
            b"ACGTACGT",
            model,
            0,
            1,
            &[0, 1, 2, 3, 4, 5, 6, 7],
        );
    }
    assert_eq!(
        position_wise_matches(b"ACGTACGT", b"ACGTACGT"),
        (0..8).collect::<Vec<u32>>()
    );
}

/// A block with no shared base at all: nothing is unchanged under either model.
///
/// The counts diverge sharply even though the sets agree. Under substitution
/// cost 1 the four substitutions are the only minimal alignment; under
/// substitution cost 2 each of them ties with a deletion/insertion pair, and
/// every interleaving of the result is optimal — 321 alignments for a four-base
/// block, the central Delannoy number `D(4, 4)`. That is the growth
/// [`DEFAULT_ALIGNMENT_CAP`] exists to bound.
#[test]
fn a_block_with_no_shared_base_leaves_nothing_unchanged() {
    assert_block(b"AAAA", b"TTTT", CostModel::Levenshtein, 4, 1, &[]);
    assert_block(
        b"AAAA",
        b"TTTT",
        CostModel::SubstitutionCostsTwo,
        8,
        321,
        &[],
    );
}

/// Empty blocks, a pure deletion and a pure insertion.
#[test]
fn empty_and_single_sided_blocks_are_handled() {
    for model in CostModel::ALL {
        assert_block(b"", b"", model, 0, 1, &[]);
        assert_block(b"A", b"", model, 1, 1, &[]);
        assert_block(b"", b"A", model, 1, 1, &[]);
        assert_block(b"ACGT", b"", model, 4, 1, &[]);
    }
}

/// **Neither notion contains the other**, so the record's rule cannot be
/// approximated by the position-wise reading in either direction.
///
/// Measured by brute force over every `{A, C}` pair up to length 5 under
/// substitution cost 1: of 1364 equal-length pairs, 384 disagree — 178 have an
/// unchanged base that is not a position-wise match, and 250 have a
/// position-wise match that is not unchanged. Both witnesses below come from
/// that sweep.
///
/// This is worth pinning because the natural guess is that one set bounds the
/// other, and a property test written on that guess passes vacuously or fails
/// mysteriously. It does neither: the two are simply incomparable.
#[test]
fn neither_notion_contains_the_other() {
    // Unchanged but not position-wise matched: offset 1 is matched in every
    // minimal alignment of `ACA -> CAC`, and `A != C` position-wise.
    let report = minimal_alignments(b"ACA", b"CAC", CostModel::Levenshtein).unwrap();
    assert_eq!(report.unchanged(), &[1]);
    assert_eq!(position_wise_matches(b"ACA", b"CAC"), Vec::<u32>::new());

    // Position-wise matched but not unchanged: `AAC -> ACA` agrees at offset 0,
    // yet no reference base survives every minimal alignment.
    let report = minimal_alignments(b"AAC", b"ACA", CostModel::Levenshtein).unwrap();
    assert_eq!(report.unchanged(), &[] as &[u32]);
    assert_eq!(position_wise_matches(b"AAC", b"ACA"), vec![0]);
}

/// **The two models are incomparable too** — neither's unchanged set contains
/// the other's.
///
/// `1420-v4` gives one direction (substitution cost 1 sees `{1, 2}`, cost 2 sees
/// `{2}`). `AAC -> CAA` gives the other, which is the direction that is easy to
/// assume away: the *stricter-looking* model reports **more** unchanged bases,
/// because pricing a substitution at 2 removes alignments from the optimal set
/// and an intersection over fewer alignments is larger.
#[test]
fn neither_cost_model_bounds_the_other() {
    let unit = minimal_alignments(b"AAC", b"CAA", CostModel::Levenshtein).unwrap();
    let paired = minimal_alignments(b"AAC", b"CAA", CostModel::SubstitutionCostsTwo).unwrap();
    assert_eq!(unit.unchanged(), &[1]);
    assert_eq!(paired.unchanged(), &[0, 1]);
}

/// **The four rows another decided record calls SPEC-CONFORMANT are all
/// position-wise readings, and not one of them is minimal.**
///
/// `rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]` justifies four
/// moved rows on the ground that "the four are EQUAL-LENGTH blocks, where the
/// column correspondence is unique and so the changed-column set is a fact
/// rather than a choice", and then applies `DNA/delins.md:15`/`:16`/`:17` to the
/// changed columns that ground yields. The record *this* module pins falsifies
/// it: equal length does not make the correspondence unique, and on these four
/// it does not even make position-wise minimal.
///
/// Two `decided` records resting on contradictory premises is the thing no guard
/// was watching, so this is that guard. It asserts the arithmetic rather than the
/// disposition — which of the two forms ferro should emit is an adjudication and
/// is not settled here.
///
/// | row's block | position-wise | minimal | unchanged |
/// |---|---:|---:|---|
/// | `TTTTTTTAAT -> ATTTTTTTAA` | cost 3 | **cost 2**, 1 alignment | `{0..=8}` |
/// | `CAG -> AGA` | cost 3 | **cost 2**, 1 alignment | `{1, 2}` |
/// | `AATA -> TAAT` | cost 3 | **cost 2**, 1 alignment | `{0, 1, 2}` |
/// | `GTAAAA -> TAAAAG` | cost 3 | **cost 2**, 1 alignment | `{1..=5}` |
///
/// All four are one shape: one base inserted at one end of the block and a
/// **different** base deleted at the other, everything between matching. That is
/// why every row is **model-independent** — with no substitution in the optimal
/// alignment there is no substitution to price, so the open question
/// [`the_two_cost_models_disagree_on_1420_v4`] records cannot arise here.
///
/// Not a *rotation*, which would re-insert the base it removed: only the
/// `GTAAAA -> TAAAAG` row below is one. An earlier revision of this comment
/// called all of them rotations, which is false for the four that matter and
/// true only for the row that is not one of them.
///
/// The fourth block is not one of the four. The record files it among "the
/// fifteen ... unequal-length blocks", and it is equal length 6/6 — pinned here
/// so that miscount cannot be re-derived as a fact.
#[test]
fn the_weight_bound_records_four_rows_are_position_wise_readings() {
    /// One block, and the two readings of it that the two records disagree
    /// about. Named fields rather than a tuple because `unchanged` and
    /// `position_wise` are the same type and transposing them would invert
    /// exactly the claim this test exists to make.
    struct Row {
        reference: &'static [u8],
        alternate: &'static [u8],
        /// Matched in **every** minimal alignment — the ruling's notion.
        unchanged: &'static [u32],
        /// Matched by pairing column *i* against column *i*.
        position_wise: &'static [u32],
    }

    let rows = [
        // `cis_junction_crossing_shift::the_five_prime_junction_barrier_does_not_over_clamp`
        // — `TEMPLATE:g.3_12`, from `TRACT`.
        Row {
            reference: b"TTTTTTTAAT",
            alternate: b"ATTTTTTTAA",
            unchanged: &[0, 1, 2, 3, 4, 5, 6, 7, 8],
            position_wise: &[1, 2, 3, 4, 5, 6, 8],
        },
        // `issue_1284_transcript_axis_collision::a_noncoding_del_dup_collision_is_repaired`
        // and `normalize_reparse_invariant::a_del_beside_a_dup_re_spells_instead_of_colliding`
        // — one block reached from two unrelated inputs, and the record's own
        // headline case above.
        Row {
            reference: b"CAG",
            alternate: b"AGA",
            unchanged: &[1, 2],
            position_wise: &[],
        },
        // `cis_confluence_adjudication::two_adjacent_members_that_both_consume_reference_are_one_delins`
        // — `NM_TEST.1:c.10_13`, from `CORE`. Note this row's position-wise
        // changed columns are {0, 2, 3}: TWO groups separated by one matching
        // base, not one run. So the position-wise reading does not yield the
        // spanning `c.10_13delinsTAAT` either — it yields
        // `c.[10A>T;12_13delinsAT]`. This row is therefore not the clean
        // position-wise-against-minimal contrast the other three are.
        Row {
            reference: b"AATA",
            alternate: b"TAAT",
            unchanged: &[0, 1, 2],
            position_wise: &[1],
        },
        // Filed by the record among the UNEQUAL-length fifteen, and equal 6/6.
        Row {
            reference: b"GTAAAA",
            alternate: b"TAAAAG",
            unchanged: &[1, 2, 3, 4, 5],
            position_wise: &[2, 3, 4],
        },
    ];

    for Row {
        reference,
        alternate,
        unchanged,
        position_wise,
    } in rows
    {
        assert_eq!(
            reference.len(),
            alternate.len(),
            "{} -> {}: the record's premise is about EQUAL-length blocks",
            String::from_utf8_lossy(reference),
            String::from_utf8_lossy(alternate),
        );

        for model in CostModel::ALL {
            // Cost 2 with equal lengths and no substitution: one insertion and
            // one deletion, so the answer cannot depend on the substitution
            // price.
            assert_block(reference, alternate, model, 2, 1, unchanged);
        }

        assert_eq!(
            position_wise_matches(reference, alternate),
            position_wise,
            "{} -> {}: position-wise matches",
            String::from_utf8_lossy(reference),
            String::from_utf8_lossy(alternate),
        );

        // The whole point: the position-wise reading is NOT one of the minimal
        // alignments, so the changed-column set
        // `rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]`
        // derived from it is a choice and not a fact. That is
        // `rulings[unchanged-is-read-over-every-minimal-alignment]` applied to
        // these four rows, not a reading taken here — read that record, which
        // is the only place the rule is stated. Position-wise costs one
        // substitution per mismatched column; minimal costs 2.
        let position_wise_cost = (reference.len() - position_wise.len()) as u32;
        assert!(
            position_wise_cost > 2,
            "{} -> {}: position-wise costs {position_wise_cost}, which must \
             exceed the minimal cost of 2 for this row to witness anything",
            String::from_utf8_lossy(reference),
            String::from_utf8_lossy(alternate),
        );

        // And it disagrees about which bases survive, in the direction that
        // matters: the minimal reading finds MORE unchanged bases, so the
        // record's derivation describes more change than the sequence requires.
        let minimal: std::collections::BTreeSet<u32> = unchanged.iter().copied().collect();
        let positional: std::collections::BTreeSet<u32> = position_wise.iter().copied().collect();
        assert!(
            positional.is_subset(&minimal) && positional != minimal,
            "{} -> {}: the position-wise matches must be a strict subset of the \
             unchanged set",
            String::from_utf8_lossy(reference),
            String::from_utf8_lossy(alternate),
        );
    }
}

/// Exceeding the cap is an error, never a truncated answer.
///
/// A silently truncated enumeration would make the intersection **too
/// permissive** — bases would be reported unchanged that some unenumerated
/// minimal alignment does not match — which is wrong in exactly the direction
/// the record was written to prevent. So the refusal is checked directly, and
/// checked to still report the cost, which the polynomial grid always knows.
#[test]
fn exceeding_the_cap_refuses_rather_than_truncating() {
    let exceeded = minimal_alignments_capped(
        b"AAAA",
        b"TTTT",
        CostModel::SubstitutionCostsTwo,
        320, // one below this block's 321
    )
    .expect_err("321 alignments must not fit under a cap of 320");
    assert_eq!(exceeded.cap, 320);
    assert_eq!(
        exceeded.cost, 8,
        "the cost is known even when the set is not"
    );
    assert_eq!(exceeded.model, CostModel::SubstitutionCostsTwo);

    // One more, and it fits.
    let report = minimal_alignments_capped(b"AAAA", b"TTTT", CostModel::SubstitutionCostsTwo, 321)
        .expect("321 alignments fit under a cap of 321");
    assert_eq!(report.count(), 321);

    // A cap of zero is a legitimate way to ask only for the cost.
    let cost_only = minimal_alignments_capped(b"CAG", b"AGA", CostModel::Levenshtein, 0)
        .expect_err("a cap of zero admits no alignment");
    assert_eq!(cost_only.cost, 2);
}

/// The default cap is large enough for every case this module pins, and the
/// margin is stated rather than assumed.
#[test]
fn every_pinned_case_fits_well_inside_the_default_cap() {
    let worst = [
        (&b"CAG"[..], &b"AGA"[..]),
        (b"GACA", b"AGAT"),
        (b"TTGC", b"ATTG"),
        (b"ATTG", b"CATT"),
        (b"ATGC", b"GCTG"),
        (b"AAAA", b"TTTT"),
        (b"TTTTTTTAAT", b"ATTTTTTTAA"),
        (b"AATA", b"TAAT"),
        (b"GTAAAA", b"TAAAAG"),
    ]
    .into_iter()
    .flat_map(|(r, a)| CostModel::ALL.map(move |m| (r, a, m)))
    .map(|(r, a, m)| {
        minimal_alignments(r, a, m)
            .expect("every pinned case fits the default cap")
            .count()
    })
    .max()
    .expect("the list is not empty");

    assert_eq!(worst, 321, "`AAAA -> TTTT` under substitution cost 2");
    assert!(
        worst * 4 < DEFAULT_ALIGNMENT_CAP,
        "the default cap ({DEFAULT_ALIGNMENT_CAP}) should leave headroom over \
         the worst pinned case ({worst})"
    );
}
