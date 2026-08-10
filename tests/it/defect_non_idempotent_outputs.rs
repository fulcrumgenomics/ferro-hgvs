//! The corpus's `non_idempotent_outputs` class, dissected: **4 outputs at 3' and
//! 4 at 5'** that are not their own fixed point.
//!
//! # HALF OF THIS CLASS IS FIXED — the three `scale-*` rows are gone
//!
//! It was **7 at 3'**. The branch that adds `general.md:35`'s second conjunct
//! ("together affecting one amino acid") to `coalesce_coding_frame_separation`
//! removed the whole of "Mechanism 2 of 2" below: the three
//! `scale-c3p-sep{120,128,136}-del-del` rows now settle in **one** pass, and the
//! coding-axis onset sweep that found the class on five of six `ACGT` cores finds
//! it on **none**.
//!
//! The mechanism is the one this file identified and is worth stating in its own
//! terms, because it is the confirmation: the second pass was reaching a shorter
//! partition than the first *by merging across codon boundaries the exception
//! never reached*. Remove those merges and the first pass lands where the second
//! used to. The three second-pass answers this file pinned are unchanged strings
//! — they are simply now the **first**-pass answers, which is why
//! [`pass_one_is_already_a_fixed_point_of_per_member_normalization`] can keep the
//! same table of three expected strings, moved one pass earlier.
//!
//! **The pins below are flipped, not deleted.** Every one of them now asserts the
//! fixed point, and says in its own doc what it used to assert and what changed.
//! `KNOWN_DIVERGENT_INPUTS`'s discipline applies here too: a deviation that has
//! been fixed must stay visible, so that it cannot silently come back. What is
//! **not** claimed is that the class is fixed — mechanism 1, the four `cds-end`
//! families, is untouched at both 3' and 5'.
//!
//! # What this file is for
//!
//! `spec_conformance_axis.rs` pins the class as two integers. Two integers
//! cannot say *why*, and the two hypotheses on offer when this was opened were
//! both wrong: that the three `scale-*-sep{120,128,136}` rows are a block-ladder
//! effect at 128 (they are not — see
//! [`the_scale_onset_is_not_the_corpus_scale_band`]), and that the class follows
//! from ferro normalizing cis members independently and concatenating (it does
//! not — see [`pass_one_is_already_a_fixed_point_of_per_member_normalization`]).
//! Both are pinned here as refutations, in this repository's "record what was
//! refuted, not only what was decided" tradition.
//!
//! # These pin DEFECTS, not answers
//!
//! Every expectation below is what ferro produces on
//! `feat/exhaustive-spec-corpus` at `a8c415ba`, and the idempotency failures are
//! **wrong**. Idempotency is not a spec clause — it is an invariant of any
//! normalizer, and `FERRO_ASSERT_IDEMPOTENT` exists to police it — so the
//! citations below bear on the *mechanism*, not on the invariant. When one is
//! fixed, the assertion moves to the correct value in the same commit, the
//! `spec_conformance_axis` census is re-blessed in that same commit, and the
//! comment says so.
//!
//! # Why this module is in `ORACLE_EXCLUDE`
//!
//! It measures over the spec corpus, and `FERRO_ASSERT_IDEMPOTENT` **panics** on
//! exactly the condition it measures. A panicking row contributes no output, so
//! an armed run would not redden this file — it would silently empty it. See
//! `tests/it/oracle_exclude_invariant.rs`, which fails if this module is not
//! named in `ci.yml`'s `ORACLE_EXCLUDE`.

use std::sync::OnceLock;

use ferro_hgvs::conformance::spec_corpus::{
    corpus, corpus_cores, denotation_of, CorpusBounds, Denotation, Frame, RefShape, Row, SpecCorpus,
};
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

// ---------------------------------------------------------------------------
// The class, by row id
// ---------------------------------------------------------------------------

/// The four `cds-end` families. These are the whole 5' set, and four of the
/// seven at 3'.
const CDS_END_FAMILIES: [&str; 4] = [
    "s01-c3m-cds-end-dup-dup-p1-sep0",
    "s01-c3m-cds-end-dup-dup-p1-sep1",
    "s01-c3p-cds-end-dup-dup-p1-sep0",
    "s01-c3p-cds-end-dup-dup-p1-sep1",
];

/// The three `scale-separation` families. **They were the 3'-only half of the
/// class and they are FIXED** — kept as a named list, because the tests that used
/// to pin their instability now pin that they are fixed points, and the list is
/// what stops a later generator edit from quietly un-generating the rows that
/// prove it.
const SCALE_FAMILIES: [&str; 3] = [
    "scale-c3p-sep120-del-del",
    "scale-c3p-sep128-del-del",
    "scale-c3p-sep136-del-del",
];

/// Every family holding a non-idempotent output, in id order.
///
/// **`SCALE_FAMILIES` is deliberately no longer chained in.** It was, and the
/// three rows it names are the ones the amino-acid precondition fixed; leaving
/// them here would make every caller iterate an empty spelling list, which is the
/// silent-vacuity failure `every_scale_family_is_now_a_fixed_point` exists to
/// make loud instead.
fn affected_families() -> Vec<&'static str> {
    let mut ids: Vec<&'static str> = CDS_END_FAMILIES.to_vec();
    ids.sort_unstable();
    ids
}

/// The one spelling of `row` whose output is not its own fixed point.
///
/// Fails with a message naming the flip rather than panicking on an out-of-bounds
/// index. A bare `[0]` on this list reported "the defect got fixed" as
/// `index out of bounds: the len is 0 but the index is 0`, which is the opposite
/// of the repo's rule that a fixed deviation must fail *legibly* — the same rule
/// `ferro_produces_the_form_the_spec_states` implements by failing when a listed
/// input starts matching the spec.
fn sole_non_idempotent_spelling(row: &Row, direction: ShuffleDirection) -> String {
    let spellings = non_idempotent_spellings(row, direction);
    assert_eq!(
        spellings.len(),
        1,
        "{} ({direction:?}): expected exactly one non-idempotent spelling. An EMPTY list means \
         the defect this test pins is FIXED — flip the pin to assert the fixed point and \
         re-bless `spec_conformance_axis`'s census in the same commit; more than one means the \
         class grew. Found {spellings:#?}",
        row.id
    );
    spellings[0].to_string()
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// The corpus, built once. Enumerating it is the dominant cost in this file, and
/// every test here reads the same rows, so sharing is worth more than the
/// per-test independence `spec_conformance_axis` prefers for its two directions.
fn built() -> &'static SpecCorpus {
    static BUILT: OnceLock<SpecCorpus> = OnceLock::new();
    BUILT.get_or_init(|| corpus(&CorpusBounds::default()))
}

fn row(id: &str) -> &'static Row {
    built()
        .rows
        .iter()
        .find(|row| row.id == id)
        .unwrap_or_else(|| panic!("the corpus no longer generates `{id}`"))
}

fn normalizer_for(frame: &Frame, direction: ShuffleDirection) -> Normalizer<MockProvider> {
    Normalizer::with_config(
        frame.provider().clone(),
        NormalizeConfig::default().with_direction(direction),
    )
}

/// One normalization pass.
fn normalize(frame: &Frame, input: &str, direction: ShuffleDirection) -> String {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));
    normalizer_for(frame, direction)
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("{input} must normalize: {e}"))
        .to_string()
}

/// `steps` successive passes, **without** stopping early on a repeat — the point
/// of this file is what the third pass does, and a chain that short-circuits at
/// the second cannot answer it.
fn iterate(frame: &Frame, input: &str, direction: ShuffleDirection, steps: usize) -> Vec<String> {
    let mut outputs = Vec::with_capacity(steps);
    let mut current = input.to_string();
    for _ in 0..steps {
        current = normalize(frame, &current, direction);
        outputs.push(current.clone());
    }
    outputs
}

/// Whether one pass moves `input`'s own output.
fn is_fixed_point(frame: &Frame, input: &str, direction: ShuffleDirection) -> bool {
    let once = normalize(frame, input, direction);
    normalize(frame, &once, direction) == once
}

/// The spellings of `row` whose output is not its own fixed point.
fn non_idempotent_spellings(row: &Row, direction: ShuffleDirection) -> Vec<&str> {
    let frame = row.frame();
    row.spellings
        .iter()
        .filter(|spelling| !is_fixed_point(&frame, spelling, direction))
        .map(String::as_str)
        .collect()
}

/// The `index`th core the corpus draws for `seeds` seeds, with a message rather
/// than a bare panic when the generator stops drawing that many.
///
/// The tests below are written on specific cores (`corpus_cores` emits an `AT`
/// core then an `ACGT` one per seed), so a generator that drew fewer would make
/// them panic inside `Vec::remove` with nothing naming the cause. That is the
/// #1460/#1478 drift class rather than a normalization change, and it should say
/// so.
fn core_at(seeds: u32, len: usize, index: usize) -> String {
    let mut cores = corpus_cores(seeds, len);
    assert!(
        index < cores.len(),
        "the corpus draws {} cores for {seeds} seeds at length {len}, so core {index} — the one \
         this test is written on — is no longer generated. That is generator drift, not a \
         normalization change.",
        cores.len()
    );
    cores.remove(index)
}

/// The 0-based index into `frame`'s served sequence whose HGVS label is `label`.
fn index_of(frame: &Frame, label: &str) -> usize {
    (0..frame.served().len())
        .find(|&i| frame.label(i) == label)
        .unwrap_or_else(|| panic!("no position labelled `{label}` on this frame"))
}

/// Every interbase of `frame` at which a two-base insertion comes back typed as
/// something other than an insertion, by the label of its 5' anchor.
fn interbases_that_retype(frame: &Frame, shape: RefShape) -> Vec<String> {
    let normalizer = normalizer_for(frame, ShuffleDirection::ThreePrime);
    (0..frame.served().len() - 1)
        .filter_map(|i| {
            let anchor = frame.label(i);
            let input = format!(
                "{}:{}.{anchor}_{}insCT",
                frame.accession(),
                shape.prefix(),
                frame.label(i + 1)
            );
            let parsed = parse_hgvs(&input).ok()?;
            let output = normalizer.normalize(&parsed).ok()?.to_string();
            output.contains("delins").then_some(anchor)
        })
        .collect()
}

/// A two-deletion design over `sep` unchanged bases, respelled as the single
/// spanning `delins` the corpus generates for it.
fn spanning_delins(frame: &Frame, shape: RefShape, start: usize, sep: usize) -> Option<String> {
    let served = frame.served();
    let end = start + 3 + sep;
    if end + 2 >= served.len() {
        return None;
    }
    Some(format!(
        "{}:{}.{}_{}delins{}",
        frame.accession(),
        shape.prefix(),
        frame.label(start),
        frame.label(end),
        &served[start + 2..start + 2 + sep]
    ))
}

/// The smallest separation at which the spanning `delins` stops being a fixed
/// point, anchored 100 bases into the drawn core so every shape addresses the
/// same bases.
fn onset(shape: RefShape, core: &str) -> Option<usize> {
    let frame = Frame::build(shape, core);
    let start = frame.core_offset() + 100;
    (1..=160).find(|&sep| {
        spanning_delins(&frame, shape, start, sep)
            .is_some_and(|input| !is_fixed_point(&frame, &input, ShuffleDirection::ThreePrime))
    })
}

/// `input`'s members, each normalized on its own and reassembled in place.
///
/// String surgery rather than reconstruction through the AST, so what is
/// compared is exactly what a reader sees in the assertion message.
fn each_member_normalized(frame: &Frame, input: &str, direction: ShuffleDirection) -> String {
    let (accession, body) = input.split_once(":c.").expect("a `c.` description");
    let members = body.trim_start_matches('[').trim_end_matches(']');
    let normalized: Vec<String> = members
        .split(';')
        .map(|member| {
            let alone = format!("{accession}:c.{member}");
            normalize(frame, &alone, direction)
                .split_once(":c.")
                .expect("a `c.` output")
                .1
                .to_string()
        })
        .collect();
    format!("{accession}:c.[{}]", normalized.join(";"))
}

// ---------------------------------------------------------------------------
// What the class IS
// ---------------------------------------------------------------------------

/// **Question.** Which spelling of an affected family is the non-idempotent one?
///
/// **The spanning `delins` respelling, in all 8 cases — never an authored
/// allele, and never a description with more than one member.** Each of the
/// four `cds-end` families holds exactly one non-idempotent spelling at 3' and
/// one at 5' as well, and every one of those 8 is a single-member `delins`
/// that the corpus generated as a *respelling* of a two-member design.
///
/// **It was 11 across seven families.** The three `scale-*` families no longer
/// hold one at all; that half is checked by
/// [`every_scale_family_is_now_a_fixed_point`], immediately below, so the fix
/// stays asserted rather than merely absent from this loop.
///
/// That is the first thing to know about this class, because it disposes of the
/// most natural guess. The adjacent, already-pinned defect is that ferro
/// normalizes cis members independently and concatenates
/// (`the_cds_end_flush_pair_is_its_two_members_normalized_separately` in
/// `spec_corpus_regressions.rs`). A defect of member *interaction* cannot be
/// reached from an input that has no members. Whatever moves these outputs, it
/// moves them before an allele exists.
///
/// The corpus's own census is the denominator: 4 and 4 are the only
/// non-idempotent outputs across all 58,552 spellings, pinned in
/// `spec_conformance_axis.rs`. This test does not re-scan the corpus — it names
/// the rows so a later generator edit cannot un-generate them silently, which is
/// the #1456/#1460/#1478 blindness class.
#[test]
fn the_non_idempotent_output_is_always_the_spanning_delins_respelling() {
    let mut three_prime = 0usize;
    let mut five_prime_families: Vec<&str> = Vec::new();

    for id in affected_families() {
        let row = row(id);

        let at_three = non_idempotent_spellings(row, ShuffleDirection::ThreePrime);
        assert_eq!(
            at_three.len(),
            1,
            "{id} must hold exactly one non-idempotent spelling at 3', found {at_three:#?}"
        );
        three_prime += at_three.len();

        let spelling = at_three[0];
        assert!(
            !spelling.contains('['),
            "{id}: the non-idempotent spelling must be single-member, got `{spelling}`"
        );
        assert!(
            spelling.contains("delins"),
            "{id}: the non-idempotent spelling must be the spanning delins, got `{spelling}`"
        );
        assert_ne!(
            spelling,
            row.authored_spelling(),
            "{id}: the non-idempotent spelling must be a RESPELLING, not the authored design"
        );

        let at_five = non_idempotent_spellings(row, ShuffleDirection::FivePrime);
        if !at_five.is_empty() {
            assert_eq!(
                at_five, at_three,
                "{id}: the two directions must agree on WHICH spelling moves"
            );
            five_prime_families.push(id);
        }
    }

    // Compared against the AXIS census, not against a literal. A literal `7` is
    // tautological here: the loop runs once per `affected_families()` entry and
    // each is asserted to hold exactly one spelling, so the count is the length
    // of a hardcoded list and equals itself whatever the corpus does. The axis
    // figure is derived independently — by scanning every row rather than the
    // seven named ones — so a family the corpus NEWLY makes non-idempotent moves
    // it while this list stays put, and the two disagree. That disagreement is
    // the coverage; it is the same cross-check
    // `corpus_prohibited_inputs::every_prohibition_violating_output_is_a_re_emitted_prohibited_input`
    // makes with `CENSUS_TOTAL`.
    assert_eq!(
        three_prime,
        crate::spec_conformance_axis::THREE_PRIME.non_idempotent_outputs,
        "PINNED DEFECT — the named families no longer account for the axis census's 3' \
         non-idempotent outputs. If the corpus made a new family non-idempotent, name it in \
         `CDS_END_FAMILIES` here and re-bless the axis in the same commit; if the count FELL, \
         a family was fixed — move it out of `affected_families` and flip its pins"
    );
    assert_eq!(
        five_prime_families,
        CDS_END_FAMILIES.to_vec(),
        "PINNED DEFECT — the 5' class is the four cds-end families and nothing else"
    );
    assert_eq!(
        five_prime_families.len(),
        crate::spec_conformance_axis::FIVE_PRIME.non_idempotent_outputs,
        "the 5' family list no longer accounts for the axis census's 5' non-idempotent outputs"
    );
}

/// **Question.** The three `scale-*` families were the 3'-only half of this
/// class. Are they still?
///
/// **No — every spelling of all three is now its own fixed point, in both
/// directions.** This is the FLIPPED pin: it used to be
/// [`the_non_idempotent_output_is_always_the_spanning_delins_respelling`]'s job
/// to find one moving spelling per `scale` row, and there is none to find. The
/// cause is `coalesce_coding_frame_separation` gaining `general.md:35`'s second
/// conjunct — "two variants … together affecting one amino acid" — so it stops
/// merging pairs whose merged span crosses a codon boundary. Those merges were
/// exactly what the second pass was doing that the first did not.
///
/// It exists as its own test rather than as an absence, because an absence is
/// invisible: dropping `SCALE_FAMILIES` from `affected_families` alone would have
/// left three rows with no assertion on them at all, and a fix that regressed
/// would then show up only as the axis census moving — a struct diff, in a
/// different file, with no mechanism attached. That is the
/// "a committed test that has never executed reads as coverage" trap pointed
/// backwards.
///
/// Asserted over **every** spelling of each row and in both directions, not just
/// over the one that used to move, since the merge the precondition now declines
/// is reachable from any respelling of the design.
#[test]
fn every_scale_family_is_now_a_fixed_point() {
    let mut checked = 0usize;
    for id in SCALE_FAMILIES {
        let row = row(id);
        let frame = row.frame();
        assert!(
            !row.spellings.is_empty(),
            "{id}: VACUOUS — the corpus generates no spellings for this row"
        );
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            for spelling in &row.spellings {
                let once = normalize(&frame, spelling, direction);
                assert_eq!(
                    normalize(&frame, &once, direction),
                    once,
                    "{id} ({direction:?}): FIXED, and it must stay fixed — `{spelling}` reached \
                     `{once}`, which must be its own fixed point. Before the amino-acid \
                     precondition the spanning-delins spelling of this row took two passes to \
                     settle. If this fires the defect is back: re-bless \
                     `spec_conformance_axis`'s `non_idempotent_outputs` and say so in the PR."
                );
                checked += 1;
            }
        }
    }
    assert!(
        checked >= SCALE_FAMILIES.len() * 2,
        "VACUOUS — only {checked} spellings were checked across {} rows and two directions",
        SCALE_FAMILIES.len()
    );
}

/// **Question.** Does the class reach a fixed point at all, or does it
/// oscillate? Nobody had checked, and an oscillation would be far more serious
/// than a two-step settle.
///
/// **It settles, at the second pass. `norm^3 == norm^2` for all 8.** The chain
/// is `x -> a -> b -> b`, with `a != b`: one extra pass, never more, and no
/// two-cycle. So a caller who normalizes twice gets a stable answer, and the
/// blast radius is "the first answer is wrong", not "there is no answer".
///
/// Asserted in both directions and on both strands, and asserted `a != b` as
/// well — without that the test would pass vacuously the moment the defect is
/// fixed, which is the wrong way round for a pin. When it *is* fixed, `a == b`
/// and this test must be deleted rather than relaxed.
///
/// **It was "all 11", and the `a != b` guard did not save it — the LOOP did the
/// vacating.** When the three `scale-*` rows became fixed points,
/// `non_idempotent_spellings` returned empty for them, so the guard sat inside a
/// loop body that never ran: three rows silently stopped being measured and this
/// test stayed green. The `chains` counter below is the fix, and it is the reason
/// a guard against vacuity has to be outside the iteration it is guarding.
#[test]
fn norm_cubed_equals_norm_squared_for_every_row_in_the_class() {
    let mut chains = 0usize;
    for id in affected_families() {
        let row = row(id);
        let frame = row.frame();
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            for spelling in non_idempotent_spellings(row, direction) {
                let chain = iterate(&frame, spelling, direction, 4);
                assert_ne!(
                    chain[0], chain[1],
                    "{id} ({direction:?}): PINNED DEFECT — `{spelling}` must still move on the \
                     second pass for this test to be measuring anything. If it no longer moves, \
                     the defect is FIXED: delete this file's pins and re-bless the census."
                );
                assert_eq!(
                    chain[1], chain[2],
                    "{id} ({direction:?}): the second pass is NOT a fixed point — this is an \
                     oscillation or a longer chain, which is a materially worse defect than the \
                     two-step settle recorded here. Chain: {chain:#?}"
                );
                assert_eq!(
                    chain[2], chain[3],
                    "{id} ({direction:?}): chain: {chain:#?}"
                );
                chains += 1;
            }
        }
    }
    // Outside the loop, deliberately: the `assert_ne!` above is the pin's own
    // anti-vacuity guard and it cannot fire when the list it iterates is empty.
    assert_eq!(
        chains,
        CDS_END_FAMILIES.len() * 2,
        "VACUOUS — {chains} chains were walked, but the four `cds-end` families each hold one \
         non-idempotent spelling in each of the two directions. A LOWER number means a family \
         was fixed and this test stopped measuring it silently; flip its pin and re-bless \
         `spec_conformance_axis` rather than letting the count drift"
    );
}

/// **Question.** Does the instability cost bases, or only spelling?
///
/// **Only spelling. Every step of every chain denotes the row's own sequence.**
/// That is the honest scoping of this class and it is worth stating plainly: it
/// is representation churn, not the sequence-preservation failure that
/// `a_deletion_flush_against_an_insertion_at_the_cds_end_changes_the_sequence`
/// pins. `background/basics.md:38` — "The recommendations for the description of
/// sequence variants are designed to be **stable**, **meaningful**,
/// **memorable**, and **unequivocal**" — puts stability first among the spec's
/// four stated values, so a normalizer whose answer depends on how many times it
/// has been run fails the first of them; but no base is lost.
///
/// The consequence for a fix: any of the three forms in a chain is a *legal*
/// description of the same variant, so choosing between them is the
/// `canonical-form-choice-when-both-legal` ruling, which is `decided` —
/// `general.md:157` governs, and ferro derives the description from the
/// RESULTING SEQUENCE rather than preserving the input's spelling or the
/// previously-shipped string. That names the direction a fix takes; it does not
/// make the move free. Moving these strings is still a declared representation
/// change under this repo's `Representation-Change:` trailer.
#[test]
fn the_chain_never_changes_the_sequence_it_denotes() {
    let mut steps_checked = 0usize;
    for id in affected_families() {
        let row = row(id);
        let frame = row.frame();
        let expected = row.denoted.as_deref().unwrap_or_else(|| {
            panic!("{id} must carry a denotation for this test to mean anything")
        });
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            for spelling in non_idempotent_spellings(row, direction) {
                for step in iterate(&frame, spelling, direction, 3) {
                    assert_eq!(
                        denotation_of(frame.provider(), frame.served(), &step),
                        Denotation::Sequence(expected.to_string()),
                        "{id} ({direction:?}): `{step}` denotes different bases from its row — \
                         that would move this row into the sequence-changing class, which is a \
                         rank-1 defect and must be reported as such rather than as churn"
                    );
                    steps_checked += 1;
                }
            }
        }
    }
    // The same anti-vacuity guard `norm_cubed_equals_norm_squared_for_every_row_in_the_class`
    // carries, and for the same reason: every assertion here lives inside a loop
    // over `non_idempotent_spellings`, which goes empty the moment a family is
    // fixed. Three chains stopped being walked once the `scale-*` rows settled
    // and nothing went red.
    assert_eq!(
        steps_checked,
        CDS_END_FAMILIES.len() * 2 * 3,
        "VACUOUS — {steps_checked} chain steps were checked; expected 3 per direction per \
         `cds-end` family. A lower number means a family was fixed and dropped out of the \
         iteration silently"
    );
}

// ---------------------------------------------------------------------------
// Mechanism 1 of 2: the CDS/3'UTR boundary re-types an insertion
// ---------------------------------------------------------------------------

/// **Question.** The four `cds-end` chains never contain an allele —
/// `c.*1delinsCTT -> c.72_*1insCT -> c.72delinsCCT`, three single-member forms
/// of one variant. So what actually moves between the first pass and the second?
///
/// **The *typing* of an insertion at one interbase.** Both passes place the
/// inserted `CT` at the same interbase, between `c.72` and `c.*1`; they disagree
/// only on whether to spell it as an insertion or as a `delins` that re-states
/// `c.72`. Sweeping a two-base insertion across **all 95 interbases** of the
/// transcript, **exactly one** re-types: the interbase immediately 3' of the last
/// CDS base. Everywhere else — including the 5'UTR/CDS boundary `c.-1_1`, mid-CDS
/// `c.71_72`, and inside the 3'UTR at `c.*1_*2` — an insertion stays an
/// insertion.
///
/// So this is neither a shift nor a partition: same interbase, different edit
/// type. And it is a **confluence failure surfacing as non-idempotency** — the
/// output typing depends on the input typing, so a `delins` input yields the
/// insertion form and the insertion form yields the `delins`.
///
/// `background/numbering.md:30` defines the zone the boundary separates —
/// "nucleotides downstream (3') of the translation termination codon (stop) are
/// marked with a `*` (asterisk) and numbered `c.*1`, `c.*2`, `c.*3`, etc., going
/// further downstream" — but it states **no rule for comparing positions across
/// zones**, and contains the word "allele" zero times. Ferro's cross-zone
/// arithmetic is therefore house policy, not compliance; see this repository's
/// `CLAUDE.md`, "Comparing `c.` positions across numbering zones is OUR policy".
/// That is the same locus the sequence-changing and denotes-no-sequence classes
/// localise to, which is why fixing one should be expected to move the others.
#[test]
fn an_insertion_immediately_after_the_last_cds_base_is_retyped_as_a_delins() {
    for id in CDS_END_FAMILIES {
        let row = row(id);
        let frame = row.frame();
        let shape = row.shape;

        // The sweep: of every interbase on the transcript, only one re-types.
        let retyped = interbases_that_retype(&frame, shape);
        assert_eq!(
            retyped,
            vec!["72".to_string()],
            "{id}: PINNED DEFECT — exactly one interbase re-types a two-base insertion, and it \
             is the one immediately 3' of the last CDS base (c.72). Correct behaviour: emit the \
             insertion, as all 94 other interbases do."
        );

        // …and c.72 really is the last CDS base, so the pin above is about the
        // zone boundary rather than about a magic number.
        let (_, cds_end) = frame.cds().expect("a coding frame has a CDS");
        assert_eq!(
            frame.label(cds_end - 1),
            "72",
            "{id}: c.72 must BE the last CDS base for the sweep to mean what it says"
        );

        // The named controls, spelled out so a partial fix is visible.
        for (region, anchor) in [("mid-CDS", "71"), ("3'UTR", "*1"), ("5'UTR/CDS", "-1")] {
            let at = index_of(&frame, anchor);
            let input = format!("NM_TEST.1:c.{anchor}_{}insCT", frame.label(at + 1));
            let output = normalize(&frame, &input, ShuffleDirection::ThreePrime);
            assert!(
                !output.contains("delins"),
                "{id}: the {region} control must NOT re-type — if it does, the defect is wider \
                 than the CDS end and this file's diagnosis is wrong. {input} -> {output}"
            );
        }

        // The exemplar chain, which is what the census counts.
        let spelling = sole_non_idempotent_spelling(row, ShuffleDirection::ThreePrime);
        let chain = iterate(&frame, &spelling, ShuffleDirection::ThreePrime, 2);
        assert_eq!(
            chain,
            vec!["NM_TEST.1:c.72_*1insCT", "NM_TEST.1:c.72delinsCCT"],
            "{id}: PINNED DEFECT — one variant, one interbase, two typings. Correct behaviour: \
             both passes emit the same one."
        );

        // Feeding the FIRST pass's answer back is what re-types it, so the two
        // passes disagree about a description ferro itself produced.
        assert_eq!(
            normalize(
                &frame,
                "NM_TEST.1:c.72_*1insCT",
                ShuffleDirection::ThreePrime
            ),
            "NM_TEST.1:c.72delinsCCT",
            "{id}: the insertion form is not a fixed point"
        );
        assert_eq!(
            normalize(
                &frame,
                "NM_TEST.1:c.72delinsCCT",
                ShuffleDirection::ThreePrime
            ),
            "NM_TEST.1:c.72delinsCCT",
            "{id}: the delins form IS a fixed point, which is why the chain settles at two"
        );
    }
}

/// **Question.** Is the re-typing about the *exon structure* — the boundary is
/// also near the end of the transcript — or about the CDS?
///
/// **The CDS, and specifically its last base.** On one drawn core, four
/// references built from the identical sequence:
///
/// ```text
/// c. three exons, CDS 13..84   -> re-types at c.72   (the last CDS base)
/// c. one exon,    CDS 1..95    -> re-types at c.95   (the last CDS base)
/// n. three exons, no CDS       -> re-types nowhere
/// g. no transcript at all      -> re-types nowhere
/// ```
///
/// The single-exon row is the load-bearing one: it has no exon junction anywhere,
/// a completely different CDS extent, and the re-typing tracks the CDS end
/// rather than a fixed coordinate. The `n.` row is the other half — same exons,
/// same bases, same coordinates, no CDS, no re-typing.
///
/// The two normalization stages a coding axis has and the other two do not are
/// `coalesce_coding_frame_separation` and `apply_coding_codon_exception`
/// (`src/normalize/merge.rs`), both keyed on the frame's reading frame, which
/// `general.md:35` is the clause for: "**exception**: two variants separated by
/// one nucleotide, together affecting one amino acid, should be described as a
/// 'delins'". Named as the place to look, not as a proven cause — this test
/// measures the discriminator, not the call site.
#[test]
fn the_retyping_needs_a_cds_and_tracks_the_cds_end_not_the_exon_structure() {
    let core = core_at(1, 96, 1);

    for (shape, expected) in [
        (RefShape::CodingMultiExon(Strand::Plus), vec!["72"]),
        (RefShape::CodingMultiExon(Strand::Minus), vec!["72"]),
        (RefShape::CodingSingleExon, vec!["95"]),
        (RefShape::NonCodingMultiExon(Strand::Plus), vec![]),
        (RefShape::Genomic, vec![]),
    ] {
        let frame = Frame::build(shape, &core);
        let retyped = interbases_that_retype(&frame, shape);
        assert_eq!(
            retyped,
            expected.iter().map(ToString::to_string).collect::<Vec<_>>(),
            "{}: PINNED DEFECT — an insertion re-types at the last CDS base and nowhere else. \
             A shape with no CDS re-types nowhere, which is what makes the CDS the variable.",
            shape.label()
        );
        if let (Some((_, cds_end)), Some(anchor)) = (frame.cds(), expected.first()) {
            assert_eq!(
                frame.label(cds_end - 1),
                *anchor,
                "{}: the re-typing interbase must be the CDS end on this shape too",
                shape.label()
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Mechanism 2 of 2: a long delins is partitioned twice, differently
// ---------------------------------------------------------------------------

/// **Question.** The three `scale-*` rows look nothing like the `cds-end` four:
/// their chains run `one delins -> an 9-to-13-member allele -> a 2-to-3-member
/// allele`. Same mechanism, or a second one?
///
/// **A second one, and it is axis-specific rather than length-specific.** The
/// spanning `delins` is *partitioned* on the first pass and *re-partitioned* on
/// the second, and the second answer is the shorter one. Holding the drawn core,
/// the anchor and the block geometry fixed and varying only the reference shape,
/// over the twelve cores `corpus_cores(6, 528)` draws — `corpus_cores` emits two
/// per seed, an `AT` core then an `ACGT` one, so the even indices below are the
/// two-letter alphabet:
///
/// ```text
/// core        c. (one exon)   n. (three exons)   g.
/// 0,2,4,…  AT   never             never          never
/// 1  ACGT       sep >= 84         never          never
/// 3  ACGT       sep >= 77         never          never
/// 5  ACGT       sep >= 71         never          never
/// 7  ACGT       sep >= 62         never          never
/// 9  ACGT       never             never          never
/// 11 ACGT       sep >= 69         never          never
/// ```
///
/// Twelve cores, three shapes, separations 1 through 160: **36 sweeps, and the
/// only ones that ever move are on the coding axis** — five of the six `ACGT`
/// cores, none of the six `AT` cores, and none of the twenty-four `n.`/`g.`
/// sweeps. The `n.` axis is not merely quieter — the block at the end of this
/// test pins one geometry where it splits the same block **further** than the
/// coding axis does and is a fixed point anyway, while the coding axis's
/// smaller partition is not. So the defect is not "long blocks are partitioned
/// badly"; it is that the coding axis's extra stages reach a different answer
/// the second time they see the same variant.
///
/// The `n.` shape is the load-bearing control rather than `g.`: it is built from
/// the same bytes with the same exon structure and the same numbering, so the
/// CDS is the only thing that differs. `g.` is included because it also has no
/// transcript at all.
///
/// # FLIPPED — the coding column is now `never` too, on all twelve cores
///
/// The diagnosis above was right and the table above is the evidence for it, so
/// both are kept verbatim as the record of what was measured. **The measurement
/// no longer reproduces**, and that is the point: the whole coding column is now
/// `never`, so the sweep is `0 of 36` rather than `5 of 36`.
///
/// The reason is the diagnosis itself, taken seriously. "The coding axis's extra
/// stages reach a different answer the second time they see the same variant" —
/// the extra stage was `coalesce_coding_frame_separation` merging pairs one
/// unchanged base apart *without* checking `general.md:35`'s second conjunct,
/// "together affecting one amino acid". A merge across a codon boundary is not an
/// answer the exception reaches, so the pass declining it removes the disagreement
/// between the two passes rather than papering over it.
///
/// **This test is now a negative guard and is worth more as one than as a pin.**
/// Every assertion below reads `None`/fixed-point, so it fails the moment any
/// axis — coding or not — reacquires a separation at which a spanning `delins`
/// stops being its own fixed point. The `AT`/`ACGT` and `c.`/`n.`/`g.` structure
/// is retained deliberately: it is what makes a re-acquisition attributable to an
/// axis, and re-deriving it from scratch cost the measurement that produced the
/// table above.
///
/// The `n.` half is untouched by the change and is still asserted the same way,
/// including the "splits it FURTHER and is stable anyway" contrast — which is now
/// a contrast between two stable partitions of different sizes rather than
/// between a stable and an unstable one.
#[test]
fn the_scale_class_is_a_coding_axis_repartition_not_a_long_block_effect() {
    let mut coding_onsets: Vec<Option<usize>> = Vec::new();
    for core in corpus_cores(6, 528) {
        coding_onsets.push(onset(RefShape::CodingSingleExon, &core));
        for shape in [
            RefShape::NonCodingMultiExon(Strand::Plus),
            RefShape::Genomic,
        ] {
            assert_eq!(
                onset(shape, &core),
                None,
                "{}: PINNED — no separation from 1 to 160 makes a spanning delins \
                 non-idempotent on an axis with no CDS. If this fires, the class is wider than \
                 the coding axis and the diagnosis in this file is wrong.",
                shape.label()
            );
        }
    }
    assert!(
        !coding_onsets.is_empty(),
        "VACUOUS — no cores were swept, so the coding column below is a claim about the \
         generator rather than about the normalizer"
    );
    assert_eq!(
        coding_onsets,
        // Was `[None, Some(84), None, Some(77), None, Some(71), None, Some(62),
        // None, None, None, Some(69)]` — the five ACGT onsets in the table above.
        vec![None; 12],
        "FIXED, and it must stay fixed — no separation from 1 to 160 makes a spanning delins \
         non-idempotent on the CODING axis either, now that \
         `coalesce_coding_frame_separation` tests `general.md:35`'s amino-acid conjunct. A \
         `Some(_)` here is the class coming back: it names the core index and the separation, so \
         re-measure that cell before touching anything else, and re-bless \
         `spec_conformance_axis`'s `non_idempotent_outputs` in the same commit."
    );

    // The `n.` axis is stable because its partition is a FIXED POINT, not
    // because it declines to partition. Without this the sweep above would be
    // equally consistent with "only the coding axis splits at all", which would
    // point a fix at the splitter rather than at the second pass. It is also what
    // keeps the flipped sweep above from being satisfiable by a normalizer that
    // simply stopped partitioning: both axes still split this block, into 8 and
    // 39 members respectively.
    //
    // One geometry rather than the sweep, because it takes a block long enough
    // that both axes split it: `scale-c3p-sep136`'s core, anchored inside the
    // second exon so nothing is masked by an exon junction, over 160 unchanged
    // bases. `c.` and `n.` here are the same bytes, the same three exons and the
    // same numbering; the CDS is the only difference between them.
    let core = &row("scale-c3p-sep136-del-del").core;
    let (start, sep) = (180usize, 160usize);
    let coding_shape = RefShape::CodingMultiExon(Strand::Plus);
    let coding = Frame::build(coding_shape, core);
    let noncoding_shape = RefShape::NonCodingMultiExon(Strand::Plus);
    let noncoding = Frame::build(noncoding_shape, core);

    let coding_input = spanning_delins(&coding, coding_shape, start, sep).expect("in range");
    let coding_output = normalize(&coding, &coding_input, ShuffleDirection::ThreePrime);
    let input = spanning_delins(&noncoding, noncoding_shape, start, sep).expect("in range");
    let output = normalize(&noncoding, &input, ShuffleDirection::ThreePrime);

    assert!(
        output.matches(';').count() > coding_output.matches(';').count(),
        "the `n.` axis must split this block FURTHER than the coding axis for the point to \
         hold — it is not stable because it is doing less work:\n  n.: {output}\n  c.: {coding_output}"
    );
    assert!(
        is_fixed_point(&noncoding, &input, ShuffleDirection::ThreePrime),
        "…and that larger partition must still be a fixed point: {input} -> {output}"
    );
    assert!(
        is_fixed_point(&coding, &coding_input, ShuffleDirection::ThreePrime),
        "FIXED — the coding axis's smaller partition is a fixed point too. This assertion was \
         `!is_fixed_point` and pinned the defect; it inverted when \
         `coalesce_coding_frame_separation` stopped merging across codon boundaries. If it fires \
         again the two passes disagree once more: {coding_input} -> {coding_output}"
    );
}

/// **Question.** The three affected `scale-*` rows sit at separations 120, 128
/// and 136, straddling 128 — the block-ladder boundary the corpus deliberately
/// varies (#1460). Is 128 the trigger?
///
/// **No, and the corpus's own rows say so twice.** `scale-c3p-sep120-del-del`
/// carries **no** scale band at all and is non-idempotent;
/// `scale-g-sep128-del-del` carries `canonical-pad-128` and is a fixed point. The
/// onsets measured above — 84, 77 and 71 on three different cores — are not a
/// constant, let alone 128, and the corpus's own coding scale rows simply have no
/// separation between 8 and 120 to have caught it earlier.
///
/// Recorded as a refutation because the reading is natural and would send a fix
/// at the wrong knob: a length threshold is a one-line change and would move
/// nothing. Per this repository's `CLAUDE.md`, "record what was refuted, not only
/// what was decided" — a measurement that kills a plausible belief is worth as
/// much as the ruling, because the belief recurs.
///
/// # FLIPPED — and the refutation is what the fix confirmed
///
/// The unbanded sep-120 row is now a fixed point too, so both halves of the
/// contrast read the same and the *contrast* is spent. The refutation is not: the
/// knob that fixed the class was `general.md:35`'s amino-acid conjunct, which is
/// a property of codon phase and has nothing to do with block length or with
/// `canonical-pad-128`. A length threshold would indeed have moved nothing. The
/// belief this test killed stayed dead, and it was killed for the right reason.
///
/// **What it still measures**, which is why it is flipped rather than deleted: the
/// two rows' `scale_bands` metadata, which is the corpus's own record of what it
/// is varying and is not otherwise asserted anywhere; and that neither row
/// reacquires the instability. If the class comes back, this says whether it came
/// back on a banded row, a bandless one, or both — the question the original
/// measurement had to be built to answer.
#[test]
fn the_scale_onset_is_not_the_corpus_scale_band() {
    // No band, and it no longer moves either.
    let unbanded = row("scale-c3p-sep120-del-del");
    assert!(
        unbanded.scale_bands.is_empty(),
        "the sep-120 row must carry no scale band for this contrast to hold, got {:?}",
        unbanded.scale_bands
    );
    assert_eq!(
        non_idempotent_spellings(unbanded, ShuffleDirection::ThreePrime),
        Vec::<&str>::new(),
        "FIXED — the bandless sep-120 row was non-idempotent and is not any more. This \
         assertion pinned `1`; the amino-acid conjunct took it to `0`. A non-empty list is the \
         class returning on a row that carries no scale band at all, which would again rule the \
         band out as the trigger"
    );

    // The band, and it never moved.
    let banded = row("scale-g-sep128-del-del");
    assert!(
        banded.scale_bands.contains(&"canonical-pad-128"),
        "the genomic sep-128 row must carry the 128 band, got {:?}",
        banded.scale_bands
    );
    assert!(
        non_idempotent_spellings(banded, ShuffleDirection::ThreePrime).is_empty(),
        "a row that DOES straddle 128 is a fixed point, so the band is not the trigger"
    );
}

/// **Question.** Is the re-partition the known cis-member-independence defect —
/// ferro normalizing each member on its own and concatenating — showing up a pass
/// later?
///
/// **No, and the measurement points the other way.** Take the first pass's
/// many-member output, normalize **each member on its own**, and reassemble:
/// the result is the first pass's output, unchanged, for all three rows. So the
/// first answer already *is* a fixed point of per-member normalization. What the
/// second pass does is something per-member normalization cannot: it re-derives
/// the partition **across** members, collapsing nine, eleven or thirteen of them
/// into two or three.
///
/// That inverts the usual reading. `the_cds_end_flush_pair_is_its_two_members_normalized_separately`
/// pins a defect where ferro does too *little* between members; this one is a
/// defect where it does something between members on the second pass that it did
/// not do on the first. A fix that only made member handling independent would
/// make this class **worse**, by removing the pass that reaches the shorter form.
///
/// This is also #1420's thesis in miniature — that the corrections "are
/// reconciled only by … apply the edit set to the reference, then re-derive the
/// minimal canonical partition from the resulting sequence". Ferro does exactly
/// that on pass two and not on pass one, and the gap between them is the class.
///
/// # FLIPPED — pass one now lands where pass two used to, and the table is the SAME
///
/// The three expected strings below have not changed. What changed is which pass
/// produces them: the spanning `delins` reaches each of them on the **first**
/// pass, and that answer is its own fixed point, so there is no second pass to
/// compare against. `each_member_normalized` is therefore now checked against the
/// single answer rather than against `^1` of a two-step chain, and the
/// `chain[1] != chain[0]` and "the second pass's partition is SHORTER" assertions
/// are gone, because there is no longer a `chain[1]` distinct from `chain[0]` to
/// make them about.
///
/// **The refutation this test exists for survives intact and is now stronger.**
/// The claim was: the class is not ferro normalizing cis members independently and
/// concatenating, because the first pass's output was already a fixed point of
/// per-member normalization. That is still asserted, and the second pass that used
/// to do "something per-member normalization cannot" is precisely the pass that
/// was merging across codon boundaries — a merge `general.md:35` never authorised.
/// So the mechanism that made the two passes disagree was not member interaction
/// at all, which is what the test said. The `matches_authored` column is the part
/// worth keeping most: it is unchanged too, so **sep-136 still disagrees with its
/// own authored design's output**, which is why that family remains non-confluent
/// in `spec_conformance_axis`'s census even though it is now idempotent. Idempotent
/// and confluent are different properties and this row is the demonstration.
#[test]
fn pass_one_is_already_a_fixed_point_of_per_member_normalization() {
    // The spanning `delins` respelling, named the way the corpus names it rather
    // than via `non_idempotent_spellings` — which now returns EMPTY for these
    // rows, that being the whole point of the flip.
    let spanning_delins_of = |row: &Row| -> String {
        let candidates: Vec<&String> = row
            .spellings
            .iter()
            .filter(|s| !s.contains('[') && s.contains("delins"))
            .collect();
        assert_eq!(
            candidates.len(),
            1,
            "{}: the corpus must still respell this design as exactly one single-member spanning \
             delins for this test to be about the same input, found {candidates:#?}",
            row.id
        );
        candidates[0].clone()
    };

    for id in SCALE_FAMILIES {
        let row = row(id);
        let frame = row.frame();
        let spelling = spanning_delins_of(row);
        let chain = iterate(&frame, &spelling, ShuffleDirection::ThreePrime, 2);

        assert!(
            !spelling.contains('['),
            "{id}: the input is a single delins, so nothing about members applies to it"
        );
        assert_eq!(
            each_member_normalized(&frame, &chain[0], ShuffleDirection::ThreePrime),
            chain[0],
            "{id}: PINNED — the answer is unchanged by normalizing each member independently, so \
             member independence was never what moved it. This is the refutation this test \
             exists for and it is unaffected by the fix"
        );
        assert_eq!(
            chain[1], chain[0],
            "{id}: FIXED — the second pass must now agree with the first. This assertion was \
             `assert_ne!` and pinned the re-partition; if it fires again the two passes disagree \
             once more and `spec_conformance_axis`'s `non_idempotent_outputs` needs re-blessing \
             UP, which is the one direction it forbids"
        );
    }

    // The exact answers, unchanged strings measured one pass earlier. Two of the
    // three land on the authored design's own output and the third does not, which
    // is why sep-136's family is still non-confluent.
    for (id, expected, matches_authored) in [
        (
            "scale-c3p-sep120-del-del",
            "NM_TEST.1:c.[107_108del;229_230del]",
            true,
        ),
        (
            "scale-c3p-sep128-del-del",
            "NM_TEST.1:c.[110_111del;240_241del]",
            true,
        ),
        (
            "scale-c3p-sep136-del-del",
            "NM_TEST.1:c.[114_115del;250_251del;252_253inv]",
            false,
        ),
    ] {
        let row = row(id);
        let frame = row.frame();
        let spelling = spanning_delins_of(row);
        let once = normalize(&frame, &spelling, ShuffleDirection::ThreePrime);
        assert_eq!(
            once, expected,
            "{id}: the answer moved. It is the same string this test pinned as the SECOND pass's \
             answer before the amino-acid conjunct; only the pass that produces it changed"
        );
        let authored = normalize(
            &frame,
            row.authored_spelling(),
            ShuffleDirection::ThreePrime,
        );
        assert_eq!(
            once == authored,
            matches_authored,
            "{id}: PINNED DEFECT — whether the answer agrees with the authored design's own \
             output is still unstable across the three rows, and that is a CONFLUENCE defect \
             surviving the idempotency fix. {once}, authored -> {authored}"
        );
    }
}

/// **Question.** Would sweeping this class over a multi-exon transcript find it?
///
/// **Only if the block stays inside one exon — otherwise it is masked.** A
/// spanning `delins` whose block crosses an exon junction is not partitioned at
/// all: it comes back as the single member it went in as, and is therefore
/// trivially a fixed point. Slide the same 100-base block along one core and the
/// defect switches on precisely when the block clears the junction:
///
/// ```text
/// block 1-based 121..220   crosses exon1/exon2 (176/177)   1 member   fixed point
/// block 1-based 171..270   crosses                          1 member   fixed point
/// block 1-based 181..280   inside exon 2                    6 members  NOT a fixed point
/// ```
///
/// Recorded because it is the same trap as
/// `the_five_prime_boundary_masks_the_same_per_member_defect`: a sweep that
/// happened to place its blocks across junctions would measure **zero** and read
/// as evidence of safety, when it is evidence that the corpus could not build the
/// thing (`CLAUDE.md`: "a corpus zero is a claim about the corpus, not about the
/// change"). It is also why the affected rows are all `mid-cds` — the only region
/// where a 124-base block fits inside one exon of this fixture.
///
/// # FLIPPED — the third row, and ONLY the third row
///
/// The two junction-crossing rows are byte-for-byte unchanged: still one member,
/// still trivially fixed points. The junction-clearing row moves from
/// **6 members, NOT a fixed point** to **2 members, a fixed point**, landing on
/// `c.[115_116del;213_214del]`.
///
/// **That is the same fix as everywhere else in this file, and its direction is
/// not a surprise once the mechanism is named.** 6 -> 2 is *fewer* members, i.e.
/// the shorter partition — which is precisely what the second pass used to
/// produce and what `pass_one_is_already_a_fixed_point_of_per_member_normalization`
/// now pins the first pass producing. Nothing is being merged here that was not
/// merged before; the block is being partitioned once, correctly, instead of twice.
/// Read as "something is now collapsed that used to stay whole", the row would
/// look like the opposite of the rest of the file — hence this paragraph.
///
/// **The corpus-design lesson this test was written for is untouched**, which is
/// why it keeps its name. A junction-crossing block is still left whole, so it
/// still cannot exhibit a re-partition of any kind, so a sweep confined to
/// junction-crossing blocks would still have measured zero and still would have
/// meant nothing. What has changed is only that there is no longer a defect on the
/// other side of the junction for it to have been blind to — and the guard now
/// fails if either half moves: if a crossing block starts being split, or if the
/// clearing block stops being a fixed point.
#[test]
fn a_block_crossing_an_exon_junction_is_not_split_so_the_defect_is_masked_there() {
    let core = core_at(3, 528, 1);
    let shape = RefShape::CodingMultiExon(Strand::Plus);
    let frame = Frame::build(shape, &core);
    let exons = frame.exons();
    assert!(
        !exons.is_empty(),
        "VACUOUS — a multi-exon fixture with no exons cannot have a junction for this test to \
         slide a block across"
    );
    let junction = exons[0].1;
    assert_eq!(
        junction, 176,
        "the fixture's first exon must end where this test assumes"
    );

    let sep = 96usize;
    // `(start, crosses the junction, members emitted, is a fixed point)`. The third
    // row was `(180, false, 6, false)` before the amino-acid conjunct.
    for (start, crosses, members, fixed) in [
        (120usize, true, 1, true),
        (170, true, 1, true),
        (180, false, 2, true),
    ] {
        let input = spanning_delins(&frame, shape, start, sep).expect("in range");
        // `Frame::exons` is 1-based inclusive; `start` and the block's last index
        // are 0-based, so both get converted before being compared to `junction`.
        let block_first = start + 1;
        let block_last = start + 3 + sep + 1;
        assert_eq!(
            block_first <= junction && block_last > junction,
            crosses,
            "the block from {start} must {} the junction for this row to mean what it says",
            if crosses { "cross" } else { "clear" }
        );
        let output = normalize(&frame, &input, ShuffleDirection::ThreePrime);
        assert_eq!(
            output.matches(';').count() + 1,
            members,
            "PINNED — a junction-crossing block is left whole, and a junction-clearing one is \
             partitioned once into its minimal form ({output})"
        );
        assert_eq!(
            is_fixed_point(&frame, &input, ShuffleDirection::ThreePrime),
            fixed,
            "a block that is never split cannot show a re-partition defect, so a sweep confined \
             to junction-crossing blocks would measure zero and mean nothing — that is the \
             corpus-design lesson, and it is unaffected by the junction-clearing block having \
             become a fixed point ({input} -> {output})"
        );
    }
}
