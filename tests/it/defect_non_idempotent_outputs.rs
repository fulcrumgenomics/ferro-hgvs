//! The corpus's `non_idempotent_outputs` class, dissected — and re-measured
//! against the partition stack, which **halved it without fixing either of its
//! two mechanisms**.
//!
//! # Verdict: MIXED. Read this before citing anything below.
//!
//! This file was written against the pre-partition corpus branch and
//! characterised **7 outputs at 3' and 4 at 5'**, in two mechanisms. On the
//! stack the census reads **4 and 4**. That drop is real and measured twice,
//! but it is not one story:
//!
//! | | mechanism A — CDS-end re-typing | mechanism B — coding-axis re-partition |
//! |---|---|---|
//! | before | 4 rows at 3', 4 at 5' | 3 rows at 3', 0 at 5' |
//! | after | **4 at 3', 4 at 5' — unchanged** | **0 — unreachable, not fixed** |
//!
//! **Mechanism A is not fixed.** Every measurement in it reproduces on the stack
//! to the byte: the same one interbase re-types, on the same shapes, with the
//! same chain. [`an_insertion_immediately_after_the_last_cds_base_is_retyped_as_a_delins`]
//! and [`the_retyping_needs_a_cds_and_tracks_the_cds_end_not_the_exon_structure`]
//! were written as pins against the old revision and pass unmodified.
//!
//! **Mechanism B is out of reach, not repaired.** The shipped partition rule is
//! now `PartitionRule::Preserve`: it keeps the members the *input* asserted
//! rather than re-deriving a partition from the resulting sequence. A
//! single-member spanning `delins` therefore comes back as the single member it
//! went in as, so there is no first partition for a second pass to disagree
//! with. The re-derivation partitioner that *did* disagree is untouched — see
//! [`the_scale_class_is_out_of_reach_of_the_shipped_rule_rather_than_fixed`],
//! which reaches it through `dev_partitioners::live` and shows it still cuts the
//! very same block into many members.
//!
//! # The revisions, and one correction
//!
//! - **Before:** `feat/exhaustive-spec-corpus` at `a8c415ba` — the corpus branch
//!   forked from `main` at `35de96c8`, before it was rebased onto the partition
//!   stack. Every "before" figure quoted below was re-measured there, not
//!   copied.
//! - **After:** the stack tip, `feat/exhaustive-spec-corpus` rebased onto #1571.
//! - **`main` at `35de96c8` itself cannot run any of this.** It has no
//!   `conformance::spec_corpus` module — the corpus arrives with this branch —
//!   so "the numbers on main" has always meant "on main plus the corpus". Said
//!   plainly because an earlier draft of this note cited `35de96c8` as if the
//!   sweeps had been run there.
//!
//! # The one-switch control, which is why "unreachable" is the honest word
//!
//! `FERRO_PARTITION=live` selects the retired re-derivation partitioner. On the
//! **same binary at the stack tip** it restores every figure this file
//! originally pinned:
//!
//! ```text
//! shipped (FERRO_PARTITION unset)   population 4 at 3'   onsets [None x12]
//! FERRO_PARTITION=live              population 7 at 3'   onsets [_,84,_,77,_,71,_,62,_,_,_,69]
//! ```
//!
//! Nothing about mechanism B was corrected. What changed is which partitioner is
//! handed the block. That distinction is the whole reason this file was not
//! deleted.
//!
//! # What actually makes it unreachable, and what to watch
//!
//! `Preserve` falls back to the re-derivation partitioner whenever it declines
//! (an overlapping member, a member outside the window, a surviving repeat, or a
//! failed rebuild check). That fallback is mechanism B's only remaining route
//! under the shipped rule — and it is measured at **zero**:
//! [`the_shipped_rule_never_hands_a_corpus_block_to_the_re_derivation_partitioner`]
//! pins 0 declines over the whole corpus. **If that number ever rises, this
//! class comes back**, so it is pinned as the reachability bound rather than
//! left as a remark.
//!
//! # Evidence that the disappearance is not a blind corpus
//!
//! Per this repository's `CLAUDE.md` — "a corpus zero is a claim about the
//! corpus, not about the change" — a zero here was checked three ways rather
//! than reported:
//!
//! 1. **The corpus's own census**, which does not read this file:
//!    `spec_conformance_axis` measures `non_idempotent_outputs: 4` at 3' and 4
//!    at 5' (against pins of 7 and 4).
//! 2. **The seam oracle**, which fires on *any* non-idempotent normalization
//!    anywhere rather than only on shapes this file builds. Armed with all three
//!    flags over CI's selection — 9,809 tests, corpus modules excluded —
//!    `FERRO_ASSERT_IDEMPOTENT`, `FERRO_ASSERT_REPARSE` and
//!    `FERRO_ASSERT_IN_BOUNDS` fired **0 times**. The instrument was confirmed
//!    live rather than assumed: re-running *this* module with
//!    `FERRO_ASSERT_IDEMPOTENT=1` reddens **5 of its 10 tests**, and every fire
//!    names the same surviving cds-end input (`c.*1delinsCTT ->
//!    c.72_*1insCT`) — mechanism A, and nothing else.
//! 3. **A construction attempt**, 33,720 inputs over six shape families the
//!    committed sweep cannot build (overlapping, nested and coincident members
//!    at scale; equal-length spanning `delins`; comb and reverse-complement
//!    payloads; repeat-spelled members), across five reference shapes, twelve
//!    cores and both directions: **0 hits**. That number is reported with its
//!    calibration, because on its own it would be worthless — the same harness
//!    under `FERRO_PARTITION=live` finds **412**, all of them the *original*
//!    spanning-`delins` shape. So the harness reaches the mechanism, and the six
//!    new shape families reach it under neither rule; they are blind and are
//!    named here as blind rather than counted as coverage.
//!
//! # These pin a DEFECT and a REACHABILITY BOUND, not answers
//!
//! Mechanism A's expectations are what ferro produces and they are **wrong** —
//! idempotency is not a spec clause, it is an invariant of any normalizer, so
//! the citations bear on the *mechanism* and not on the invariant. When it is
//! fixed, the assertion moves in the same commit, `spec_conformance_axis`'s
//! census is re-blessed in that same commit, and the comment says so.
//!
//! Mechanism B's expectations are the opposite kind: they pin that the class
//! **stays** out of reach, and they fail loudly if the shipped rule starts
//! handing blocks back to the re-derivation partitioner.
//!
//! # What the disappearance cost
//!
//! Worth stating so nobody reads the 7 -> 4 as free. The same rule change that
//! made the spanning `delins` a fixed point made it a fixed point *at its own
//! spelling* — so a family whose members used to converge on one output now
//! reports two. Over the corpus that is `converged` 9,140 -> 3,509 at 3'
//! (`split_two` 2,435 -> 7,943). The class did not go away because the
//! normalizer got better at it; it went away because the normalizer stopped
//! re-deriving, and confluence paid.
//!
//! # Why this module is in `ORACLE_EXCLUDE`
//!
//! It measures over the spec corpus, and `FERRO_ASSERT_IDEMPOTENT` **panics** on
//! exactly the condition it measures — verified above, 5 of 10 tests red under
//! the flag. A panicking row contributes no output, so an armed run would not
//! redden this file for the right reason: it would silently empty it. See
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

/// The four `cds-end` families — mechanism A. Before the stack these were the
/// whole 5' set and four of the seven at 3'; **after it they are the whole
/// class in both directions.**
const CDS_END_FAMILIES: [&str; 4] = [
    "s01-c3m-cds-end-dup-dup-p1-sep0",
    "s01-c3m-cds-end-dup-dup-p1-sep1",
    "s01-c3p-cds-end-dup-dup-p1-sep0",
    "s01-c3p-cds-end-dup-dup-p1-sep1",
];

/// The three `scale-separation` families — mechanism B. Non-idempotent at 3'
/// before the stack, fixed points under the shipped rule.
const SCALE_FAMILIES: [&str; 3] = [
    "scale-c3p-sep120-del-del",
    "scale-c3p-sep128-del-del",
    "scale-c3p-sep136-del-del",
];

/// Every family that held a non-idempotent output on `a8c415ba`, in id order.
///
/// Deliberately still the union of both mechanisms even though three of the
/// seven now hold none: a test that only looked at the survivors could not
/// notice mechanism B coming back.
fn affected_families() -> Vec<&'static str> {
    let mut ids: Vec<&'static str> = CDS_END_FAMILIES.into_iter().chain(SCALE_FAMILIES).collect();
    ids.sort_unstable();
    ids
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
///
/// Note the payload is four bases shorter than the block it replaces, so the
/// member is never equal-length — which is exactly why the shipped `Preserve`
/// rule leaves it whole (its interior column split is scoped to equal-length
/// members). This helper *is* the shape mechanism B needed, which is why it is
/// kept even though its sweep now measures a wall of `None`.
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
// What the class IS, now
// ---------------------------------------------------------------------------

/// **Question.** After the partition stack, which spellings are still not their
/// own fixed point, and which stopped being one?
///
/// **The four `cds-end` respellings survive, in both directions; the three
/// `scale-*` respellings do not.** So the class is 4 at 3' and 4 at 5', down
/// from 7 and 4, and the entire drop is mechanism B.
///
/// Every survivor is still the spanning `delins` respelling — never an authored
/// allele, never a description with more than one member — which is the first
/// thing to know about it, because it disposes of the most natural guess. The
/// adjacent, already-pinned defect is that ferro normalizes cis members
/// independently and concatenates
/// (`the_cds_end_flush_pair_is_its_two_members_normalized_separately` in
/// `spec_corpus_regressions.rs`). A defect of member *interaction* cannot be
/// reached from an input that has no members. Whatever moves these outputs, it
/// moves them before an allele exists.
///
/// The corpus's own census is the denominator: 4 and 4 are the only
/// non-idempotent outputs across all 58,552 spellings, and
/// `spec_conformance_axis.rs` measures them independently of this file. This
/// test does not re-scan the corpus — it names the rows so a later generator
/// edit cannot un-generate them silently, which is the #1456/#1460/#1478
/// blindness class.
///
/// The three `scale-*` rows are asserted at **zero rather than dropped**, for
/// the same reason: they are the guard that mechanism B has not come back.
#[test]
fn the_surviving_class_is_the_four_cds_end_respellings_in_both_directions() {
    let mut three_prime: Vec<&str> = Vec::new();
    let mut five_prime: Vec<&str> = Vec::new();

    for id in affected_families() {
        let row = row(id);
        let at_three = non_idempotent_spellings(row, ShuffleDirection::ThreePrime);
        let at_five = non_idempotent_spellings(row, ShuffleDirection::FivePrime);

        if SCALE_FAMILIES.contains(&id) {
            // MECHANISM B, out of reach under the shipped rule. Before the
            // stack each of these held exactly one at 3'.
            assert!(
                at_three.is_empty() && at_five.is_empty(),
                "{id}: PINNED — mechanism B is out of reach under the shipped partition rule. \
                 A non-idempotent spelling here means the re-derivation partitioner is being \
                 reached again; 3' {at_three:#?}, 5' {at_five:#?}"
            );
            continue;
        }

        // MECHANISM A, unchanged.
        assert_eq!(
            at_three.len(),
            1,
            "{id} must hold exactly one non-idempotent spelling at 3', found {at_three:#?}"
        );
        three_prime.push(id);

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

        assert_eq!(
            at_five, at_three,
            "{id}: the two directions must agree on WHICH spelling moves"
        );
        five_prime.push(id);
    }

    assert_eq!(
        three_prime,
        CDS_END_FAMILIES.to_vec(),
        "PINNED DEFECT — the 3' class is the four cds-end families and nothing else (it was \
         those four PLUS the three scale families on `a8c415ba`). If this moved, re-bless \
         `spec_conformance_axis`'s census in the same commit"
    );
    assert_eq!(
        five_prime,
        CDS_END_FAMILIES.to_vec(),
        "PINNED DEFECT — the 5' class is the four cds-end families and nothing else"
    );
}

/// **Question.** Does the class reach a fixed point at all, or does it
/// oscillate? Nobody had checked, and an oscillation would be far more serious
/// than a two-step settle.
///
/// **It settles, at the second pass. `norm^3 == norm^2` for all of them.** The
/// chain is `x -> a -> b -> b`, with `a != b`: one extra pass, never more, and
/// no two-cycle. So a caller who normalizes twice gets a stable answer, and the
/// blast radius is "the first answer is wrong", not "there is no answer".
///
/// Asserted in both directions and on both strands, and asserted `a != b` as
/// well — without that the test would pass vacuously the moment the defect is
/// fixed, which is the wrong way round for a pin. When it *is* fixed, `a == b`
/// and this test must be deleted rather than relaxed.
///
/// It also asserts its own denominator. On `a8c415ba` this covered 11 chains;
/// under the shipped rule it covers 8, and a version of this test that merely
/// looped over whatever it found would have gone green at 0 — which is the
/// failure this whole file exists to not repeat.
#[test]
fn norm_cubed_equals_norm_squared_for_every_row_in_the_class() {
    let mut chains = 0usize;
    for id in affected_families() {
        let row = row(id);
        let frame = row.frame();
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            for spelling in non_idempotent_spellings(row, direction) {
                chains += 1;
                let chain = iterate(&frame, spelling, direction, 4);
                assert_ne!(
                    chain[0], chain[1],
                    "{id} ({direction:?}): PINNED DEFECT — `{spelling}` must still move on the \
                     second pass for this test to be measuring anything. If it no longer moves, \
                     mechanism A is FIXED: delete this file's pins and re-bless the census."
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
            }
        }
    }
    assert_eq!(
        chains, 8,
        "PINNED — 8 chains (4 cds-end families x 2 directions). It was 11 on `a8c415ba`, the \
         extra 3 being the scale rows at 3'. A change in this number is a change in the class, \
         and a ZERO would make every assertion above vacuous"
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
/// `canonical-form-choice-when-both-legal` ruling (still `undecided`), and
/// whichever way it lands, moving these strings is a declared representation
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
                    steps_checked += 1;
                    assert_eq!(
                        denotation_of(frame.provider(), frame.served(), &step),
                        Denotation::Sequence(expected.to_string()),
                        "{id} ({direction:?}): `{step}` denotes different bases from its row — \
                         that would move this row into the sequence-changing class, which is a \
                         rank-1 defect and must be reported as such rather than as churn"
                    );
                }
            }
        }
    }
    assert_eq!(
        steps_checked, 24,
        "PINNED — 8 chains x 3 steps. A zero here would make this test a claim about nothing"
    );
}

// ---------------------------------------------------------------------------
// Mechanism 1 of 2: the CDS/3'UTR boundary re-types an insertion — SURVIVES
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
/// **Unchanged by the partition stack.** This test was written as a pin against
/// `a8c415ba` and passes on the stack without a single value moved: the same
/// interbase, the same three controls, the same two-step chain, the same fixed
/// point. Mechanism A is not what the partition rule touches — which is the
/// evidence that the 7 -> 4 drop is entirely mechanism B.
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
        let spelling = non_idempotent_spellings(row, ShuffleDirection::ThreePrime)[0].to_string();
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
/// **Unchanged by the partition stack**, to the same four rows and the same two
/// anchors. Re-measured on the stack rather than assumed.
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
    let core = corpus_cores(1, 96).remove(1);

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
// Mechanism 2 of 2: a long delins was partitioned twice, differently —
// now UNREACHABLE, and this section is about proving that word
// ---------------------------------------------------------------------------

/// **Question — the whole point of this file's second half.** The three
/// `scale-*` rows were non-idempotent on `a8c415ba` and are fixed points on the
/// stack. Was the defect **fixed**, or did the shapes that produced it stop
/// reaching the code that produced it?
///
/// **The second. The re-derivation partitioner still cuts the very same block
/// into the very same many members; the shipped rule simply never asks it to.**
///
/// The mechanism, as measured on `a8c415ba`: a spanning `delins` was
/// *partitioned* on the first pass and *re-partitioned* on the second, and the
/// second answer was the shorter one. Holding the drawn core, the anchor and the
/// block geometry fixed and varying only the reference shape, over the twelve
/// cores `corpus_cores(6, 528)` draws — `corpus_cores` emits two per seed, an
/// `AT` core then an `ACGT` one, so the even indices below are the two-letter
/// alphabet:
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
/// Under the shipped rule **all thirty-six of those sweeps are `None`**,
/// including the coding column. That is what this test pins.
///
/// And here is why "unreachable" rather than "fixed", in three measurements:
///
/// 1. **The same binary reproduces the old table.** `FERRO_PARTITION=live` at
///    the stack tip gives back `[None, 84, None, 77, None, 71, None, 62, None,
///    None, None, 69]`, i.e. every original onset. The knob is process-global
///    and read once into a `OnceLock`, and `merge.rs`'s own suite asserts it is
///    unset, so it cannot be flipped inside a test — it is quoted here as a
///    measurement and reproduced by hand with
///    `FERRO_PARTITION=live cargo run --release --features dev --example dump_partitions`.
/// 2. **The partitioner is reachable from this test and still disagrees.** The
///    block below is handed to `dev_partitioners::live`, which *is*
///    `partition_block`, the rule that was the default before the stack. It
///    still cuts the block into many members while the shipped normalization
///    returns one. No environment variable involved; this is an assertion.
/// 3. **What changed is upstream of the partitioner.** `PartitionRule::Preserve`
///    keeps the members the input asserted. A spanning `delins` asserts exactly
///    one member and its payload is four bases shorter than its block, so
///    `Preserve`'s interior column split — which is scoped to *equal-length*
///    members — does not apply, and the member survives whole. Measured: a block
///    that came back as **6 members and was not a fixed point** on `a8c415ba`
///    comes back as **1 member and is** one now.
#[cfg(feature = "dev")]
#[test]
fn the_scale_class_is_out_of_reach_of_the_shipped_rule_rather_than_fixed() {
    use ferro_hgvs::normalize::dev_partitioners;

    // (1) The sweep, now silent on every axis including the coding one.
    for core in corpus_cores(6, 528) {
        for shape in [
            RefShape::CodingSingleExon,
            RefShape::NonCodingMultiExon(Strand::Plus),
            RefShape::Genomic,
        ] {
            assert_eq!(
                onset(shape, &core),
                None,
                "{}: PINNED — under the shipped partition rule no separation from 1 to 160 \
                 makes a spanning delins non-idempotent, on ANY axis. The coding axis lost \
                 idempotency at 84/77/71/62/69 on five of these cores before the stack; if a \
                 number reappears here, the re-derivation partitioner is being reached again.",
                shape.label()
            );
        }
    }

    // (2) The partitioner itself, reached directly. `dev_partitioners::live` is
    // `partition_block` — the pre-stack default — and it still disagrees with
    // the shipped answer on the identical block.
    let row = row("scale-c3p-sep136-del-del");
    let core = &row.core;
    let shape = RefShape::CodingMultiExon(Strand::Plus);
    let frame = Frame::build(shape, core);
    let (start, sep) = (180usize, 160usize);
    let input = spanning_delins(&frame, shape, start, sep).expect("in range");
    let output = normalize(&frame, &input, ShuffleDirection::ThreePrime);
    assert_eq!(
        output.matches(';').count(),
        0,
        "PINNED — the shipped rule preserves the single member the input asserted: {input} -> \
         {output}"
    );
    assert!(
        is_fixed_point(&frame, &input, ShuffleDirection::ThreePrime),
        "…and a member that is never split cannot be re-split, so it is a fixed point: {input}"
    );

    // The same block, through the retired partitioner. `spanning_delins` builds
    // `label(start) .. label(start + 3 + sep)` inclusive with a payload drawn
    // from `served[start + 2 ..]`, so those are the two slices to compare.
    let served = frame.served().as_bytes();
    let reference = &served[start..=start + 3 + sep];
    let alt = &served[start + 2..start + 2 + sep];
    let pieces = dev_partitioners::live(reference, alt);
    assert!(
        pieces.len() > 1,
        "PINNED — `partition_block` is untouched and still re-derives a many-member partition \
         of the very block the shipped rule now leaves whole. If this ever returns one piece, \
         mechanism B really has been fixed and this test should say so instead. Got {pieces:#?}"
    );
}

/// **Question.** `Preserve` falls back to the re-derivation partitioner whenever
/// it declines. That fallback is mechanism B's only remaining route under the
/// shipped rule — so how often does it fire?
///
/// **Never, over the whole corpus.** Measured across every spelling in both
/// directions, `partition_block_preserving` is kept on every block it is offered
/// and declines **zero** times. That, and not any repair, is what currently
/// makes mechanism B unreachable.
///
/// This is pinned as a **reachability bound**, which is a different kind of
/// assertion from everything above it. The others say "this defect is still
/// here" or "this output is right". This one says: *the only door back to the
/// defect is shut, and here is the number that will tell you when it opens.* A
/// change that makes `Preserve` decline — a new member geometry, a repeat the
/// window cannot lower, a failed rebuild check — hands blocks back to
/// `partition_block`, and the class returns without anything in this file
/// otherwise noticing.
///
/// The census is process-global and never reset between tests, which is safe
/// here only because nextest runs each test in its own process; it is reset
/// explicitly anyway so corpus construction cannot contribute.
///
/// Why a counter and not a log: `dev_partitioners::preserve_census`'s own doc
/// records that no logger is installed anywhere in this crate, so all 19
/// `log::` call sites are no-ops and a grep-based measurement returns 0 whether
/// or not the branch ran. This arm was once reported as preserving 100% of
/// blocks when it preserved 66.6%.
#[cfg(feature = "dev")]
#[test]
fn the_shipped_rule_never_hands_a_corpus_block_to_the_re_derivation_partitioner() {
    use ferro_hgvs::normalize::dev_partitioners;

    let built = built();
    dev_partitioners::reset_preserve_census();

    let mut attempted = 0usize;
    let mut normalizations = 0usize;
    for row in &built.rows {
        let frame = row.frame();
        for spelling in &row.spellings {
            for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
                attempted += 1;
                // The count of what actually reached the normalizer is taken
                // *after* the parse, not before it. With one counter
                // incremented ahead of this `continue`, an unparseable spelling
                // still raised the denominator, so the pin below could stay
                // green while strictly fewer blocks were offered to
                // `partition_block_preserving` — which weakens the
                // `declined == 0` claim in exactly the way this test's doc says
                // it must not.
                //
                // A parse failure is NOT a corpus defect here, which is why this
                // is two counters rather than a panic: the corpus deliberately
                // carries prohibited inputs, and the parser rejects some of them
                // outright (`NM_TEST.1:c.[9_12del;11_14dup]` is refused as a
                // self-cancelling allele). Both totals are pinned so neither the
                // corpus nor the parser's acceptance can shift unobserved.
                let Ok(parsed) = parse_hgvs(spelling) else {
                    continue;
                };
                normalizations += 1;
                let _ = normalizer_for(&frame, direction).normalize(&parsed);
            }
        }
    }

    let (kept, declined) = dev_partitioners::preserve_census();
    assert_eq!(
        attempted, 117_104,
        "PINNED — the corpus. 58,552 spellings x 2 directions; if the corpus shrinks, the \
         zero below becomes a claim about the corpus"
    );
    assert_eq!(
        normalizations, 116_872,
        "PINNED — the denominator that matters: how many of the {attempted} attempts the \
         parser accepted and so actually reached `partition_block_preserving`. MEASURED: 232 \
         of the corpus total (116 spellings x 2 directions) never get here, because the \
         corpus carries prohibited inputs the parser refuses outright. Pinned apart from the \
         corpus total because a parser change that started rejecting rows would otherwise \
         shrink this silently while the total above stayed green — which is exactly what one \
         combined counter allowed"
    );
    assert!(
        kept > 100_000,
        "VACUOUS — the preserving arm was reached only {kept} times, so `declined == 0` would \
         say nothing about it"
    );
    assert_eq!(
        declined, 0,
        "PINNED REACHABILITY BOUND — `partition_block_preserving` declined {declined} times, \
         which hands that many blocks to `partition_block`. That is the retired re-derivation \
         partitioner, and it is still non-idempotent at scale (see \
         `the_scale_class_is_out_of_reach_of_the_shipped_rule_rather_than_fixed`). A non-zero \
         here means mechanism B has a live route again — re-run the onset sweep before \
         re-blessing this number. kept={kept}"
    );
}

/// **Question.** The three affected `scale-*` rows sat at separations 120, 128
/// and 136, straddling 128 — the block-ladder boundary the corpus deliberately
/// varies (#1460). Was 128 the trigger?
///
/// **No, and the corpus's own rows said so twice.** On `a8c415ba`,
/// `scale-c3p-sep120-del-del` carried **no** scale band at all and was
/// non-idempotent, while `scale-g-sep128-del-del` carried `canonical-pad-128`
/// and was a fixed point. The onsets measured there — 84, 77, 71, 62 and 69 on
/// five different cores — were not a constant, let alone 128.
///
/// Recorded as a refutation because the reading is natural and would send a fix
/// at the wrong knob: a length threshold is a one-line change and would have
/// moved nothing. Per this repository's `CLAUDE.md`, "record what was refuted,
/// not only what was decided" — a measurement that kills a plausible belief is
/// worth as much as the ruling, because the belief recurs.
///
/// **The refutation survives the stack from the other side.** Under the shipped
/// rule *both* rows are fixed points, so the band still fails to discriminate —
/// it just now fails by predicting instability that is not there rather than
/// stability that was not. Asserted that way round rather than deleted, because
/// the pairing (`no band + moves` against `band + does not`) is the whole
/// argument and it has to remain visible when mechanism B comes back.
#[test]
fn the_scale_onset_is_not_the_corpus_scale_band() {
    // No band. It moved on `a8c415ba`; it does not now.
    let unbanded = row("scale-c3p-sep120-del-del");
    assert!(
        unbanded.scale_bands.is_empty(),
        "the sep-120 row must carry no scale band for this contrast to hold, got {:?}",
        unbanded.scale_bands
    );

    // The band, and it did not move either — before the stack or after.
    let banded = row("scale-g-sep128-del-del");
    assert!(
        banded.scale_bands.contains(&"canonical-pad-128"),
        "the genomic sep-128 row must carry the 128 band, got {:?}",
        banded.scale_bands
    );

    // The refutation, stated as the comparison rather than as two counts: the
    // band never separated the two rows. Before the stack the unbanded row moved
    // and the banded one did not; now neither does. Either way the band is not
    // the variable.
    let unbanded_moves =
        !non_idempotent_spellings(unbanded, ShuffleDirection::ThreePrime).is_empty();
    let banded_moves = !non_idempotent_spellings(banded, ShuffleDirection::ThreePrime).is_empty();
    assert!(
        !banded_moves,
        "PINNED — a row that DOES straddle 128 has never been non-idempotent, so the band is \
         not the trigger"
    );
    assert!(
        !unbanded_moves,
        "PINNED — under the shipped rule the unbanded sep-120 row is a fixed point too. It was \
         NOT one on `a8c415ba`, which is where this refutation's other half came from; if it \
         moves again, mechanism B is back and the onset sweep is the thing to re-run"
    );
}

/// **Question.** Was the re-partition the known cis-member-independence defect —
/// ferro normalizing each member on its own and concatenating — showing up a
/// pass later?
///
/// **No, and the measurement pointed the other way.** On `a8c415ba`: take the
/// first pass's many-member output, normalize **each member on its own**, and
/// reassemble — the result was the first pass's output, unchanged, for all three
/// rows. So the first answer already *was* a fixed point of per-member
/// normalization. What the second pass did was something per-member
/// normalization cannot: it re-derived the partition **across** members,
/// collapsing nine, eleven or thirteen of them into two or three. For the
/// record, those first-pass answers were, at 3':
///
/// ```text
/// sep120  c.[107_108del;194_195del;196T>G;200_201delinsTG;203_205delinsTAC;
///           207_208delinsGT;213_221delinsCCAGGACTA;223_228delinsGGCCAA;230A>C]
///          -> c.[107_108del;229_230del]
/// sep128  11 members -> c.[110_111del;240_241del]
/// sep136  13 members -> c.[114_115del;250_251del;252_253inv]
/// ```
///
/// That inverted the usual reading.
/// `the_cds_end_flush_pair_is_its_two_members_normalized_separately` pins a
/// defect where ferro does too *little* between members; this one was a defect
/// where it did something between members on the second pass that it did not do
/// on the first. A fix that only made member handling independent would have
/// made this class **worse**, by removing the pass that reached the shorter form.
///
/// It was also #1420's thesis in miniature — that the corrections "are
/// reconciled only by … apply the edit set to the reference, then re-derive the
/// minimal canonical partition from the resulting sequence". Ferro did exactly
/// that on pass two and not on pass one, and the gap between them was the class.
///
/// **What the shipped rule does instead: neither pass.** `Preserve` asserts the
/// input's own partition, so the spanning `delins` is returned unchanged and the
/// question of which pass re-derives never arises. That is what this test now
/// pins, per row and by exact string — the strings are the tripwire, because
/// "is a fixed point" alone would stay true if the rule started splitting into
/// some other stable partition.
#[test]
fn the_shipped_rule_returns_the_scale_respelling_unchanged_instead_of_partitioning_it() {
    for id in SCALE_FAMILIES {
        let row = row(id);
        let frame = row.frame();

        // The respelling: the one single-member spelling in the family, and the
        // one that was non-idempotent before the stack.
        let respelling = row
            .spellings
            .iter()
            .find(|spelling| !spelling.contains('['))
            .unwrap_or_else(|| panic!("{id} must carry a single-member respelling"))
            .clone();
        assert!(
            respelling.contains("delins"),
            "{id}: the respelling must be the spanning delins, got `{respelling}`"
        );

        let chain = iterate(&frame, &respelling, ShuffleDirection::ThreePrime, 2);
        assert_eq!(
            chain[0], respelling,
            "{id}: PINNED — the shipped rule returns the spanning delins UNCHANGED. It used to \
             partition it into 9, 11 or 13 members on the first pass and re-derive 2 or 3 on \
             the second; both passes are now no-ops. A different first answer here means the \
             partition rule moved."
        );
        assert_eq!(chain[0], chain[1], "{id}: …and so it is a fixed point");

        // Per-member normalization is a no-op on it too — trivially, since it
        // has one member. Asserted so the helper stays exercised and so the
        // contrast with the pre-stack many-member output is explicit.
        assert_eq!(
            each_member_normalized(&frame, &chain[0], ShuffleDirection::ThreePrime),
            format!(
                "{}:c.[{}]",
                chain[0].split_once(":c.").expect("a `c.` output").0,
                chain[0].split_once(":c.").expect("a `c.` output").1
            ),
            "{id}: the single member normalizes to itself"
        );
    }

    // The authored design's own output is unchanged from `a8c415ba` — so what
    // moved is the respelling, not the design. This is also where the family's
    // confluence went: the two spellings used to reach one output and now reach
    // two, which is the trade the module docs quantify (converged 9,140 ->
    // 3,509 at 3').
    for (id, authored_output) in [
        (
            "scale-c3p-sep120-del-del",
            "NM_TEST.1:c.[107_108del;229_230del]",
        ),
        (
            "scale-c3p-sep128-del-del",
            "NM_TEST.1:c.[110_111del;240_241del]",
        ),
        (
            "scale-c3p-sep136-del-del",
            "NM_TEST.1:c.[114_115del;253_254del]",
        ),
    ] {
        let row = row(id);
        let frame = row.frame();
        let authored = normalize(
            &frame,
            row.authored_spelling(),
            ShuffleDirection::ThreePrime,
        );
        assert_eq!(
            authored, authored_output,
            "{id}: the authored design's output moved"
        );

        let respelling = row
            .spellings
            .iter()
            .find(|spelling| !spelling.contains('['))
            .expect("a single-member respelling");
        let respelled = normalize(&frame, respelling, ShuffleDirection::ThreePrime);
        assert_ne!(
            respelled, authored,
            "{id}: PINNED — the design and its respelling reach DIFFERENT outputs under the \
             shipped rule. That is the confluence this class was traded for, and it is stated \
             here so the 7 -> 4 improvement is never quoted without it"
        );
    }
}

/// **Question.** Would sweeping this class over a multi-exon transcript find it?
///
/// **On `a8c415ba`: only if the block stayed inside one exon.** A spanning
/// `delins` whose block crossed an exon junction was not partitioned at all — it
/// came back as the single member it went in as, and was therefore trivially a
/// fixed point. Sliding the same 100-base block along one core, the defect
/// switched on precisely when the block cleared the junction:
///
/// ```text
///                                          a8c415ba              shipped rule
/// block 1-based 121..220  crosses 176/177   1 member  fixed       1 member  fixed
/// block 1-based 171..270  crosses           1 member  fixed       1 member  fixed
/// block 1-based 181..280  inside exon 2     6 members NOT fixed   1 member  fixed
/// ```
///
/// **Under the shipped rule the junction is no longer the discriminator,
/// because nothing is split.** All three positions come back as one member. The
/// masking is therefore not gone — it has been superseded by a rule that leaves
/// every such block whole, which produces the same *measurement* for a
/// completely different reason. Pinned in that form, with the third row asserted
/// as the one that changed, so the distinction cannot be lost.
///
/// Recorded because the original finding is the same trap as
/// `the_five_prime_boundary_masks_the_same_per_member_defect`: a sweep that
/// happened to place its blocks across junctions would measure **zero** and read
/// as evidence of safety, when it is evidence that the corpus could not build the
/// thing (`CLAUDE.md`: "a corpus zero is a claim about the corpus, not about the
/// change"). It is also why the affected rows were all `mid-cds` — the only
/// region where a 124-base block fits inside one exon of this fixture.
#[test]
fn a_block_crossing_an_exon_junction_is_not_split_and_now_neither_is_one_inside_an_exon() {
    let core = corpus_cores(3, 528).remove(1);
    let shape = RefShape::CodingMultiExon(Strand::Plus);
    let frame = Frame::build(shape, &core);
    let junction = frame.exons()[0].1;
    assert_eq!(
        junction, 176,
        "the fixture's first exon must end where this test assumes"
    );

    let sep = 96usize;
    for (start, crosses, members_before, fixed_before) in [
        (120usize, true, 1, true),
        (170, true, 1, true),
        (180, false, 6, false),
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
            1,
            "PINNED — under the shipped rule EVERY one of these blocks is left whole, junction \
             or no junction (this one was {members_before} member(s) on `a8c415ba`): {output}"
        );
        assert!(
            is_fixed_point(&frame, &input, ShuffleDirection::ThreePrime),
            "PINNED — and so every one is a fixed point (this one was fixed={fixed_before} on \
             `a8c415ba`): {input}"
        );

        // The row that changed, named as such. Without this the block above
        // reads as three identical results and the reader cannot tell that one
        // of them used to be the only one that measured anything.
        if !crosses {
            assert_eq!(
                (members_before, fixed_before),
                (6, false),
                "the in-exon row is the one whose behaviour the partition rule moved; if this \
                 record is edited, the contrast above stops meaning anything"
            );
        }
    }
}
