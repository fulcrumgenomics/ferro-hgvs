//! The equivalence relation, made fit to serve as the confluence release gate.
//!
//! # Why the relation cannot be defined by `normalize`
//!
//! The gate the project wants to state is "`normalize` is constant on each
//! equivalence class". `normalize` maps description -> description, so if the
//! class were itself defined by `normalize` the statement would be a tautology.
//! The relation has to come from `apply`, which maps description -> **sequence**
//! — a different codomain, where equality is byte equality on bases.
//!
//! # Why apply-equality on ONE axis is not enough
//!
//! `LRG_199t1:c.3921dup` and `LRG_199t1:c.3922dup` denote the same *transcript*
//! sequence (both positions carry `T`), so single-axis apply-equality calls them
//! equivalent. They project 2,790 bp apart, into different exons, and the spec
//! says so in as many words — `DNA/duplication.md:148`: "one would end up at the
//! wrong nucleotide, in the wrong exon". The relation therefore has to hold on
//! **every determined axis**, which for a `c.`/`n.`/`r.` description is the
//! transcript *and* the genome its exon alignment carries it to.
//!
//! These tests pin that with a two-exon synthetic transcript, because
//! `SyntheticBuilder` builds single-exon transcripts only and so cannot express
//! the geometry at all.
//!
//! # Why the denotational tests do not go through `check` (#1800)
//!
//! `check` is the whole ladder, and it short-circuits at `NormalizedMatch` the
//! moment two spellings normalize to one string — so a denotational rung is
//! reached only for a pair that **diverges**. Driving these tests through
//! `check` therefore made every one of them depend on the normalizer staying
//! non-confluent over some pair, and reducing divergence is the project's own
//! goal (README rule 3). Twice now a confluence PR has destroyed the pair a
//! test here was standing on: #1649 converged #1419's `[19_23del;27_33del]` /
//! `19_33delinsCGG`, and #1616 converges #1420-v2's `c.[48dup;52del]` /
//! `c.49_52delinsATTG` and its genomic sibling.
//!
//! Swapping in a fresh diverging pair buys one PR and guarantees a third
//! occurrence, and it also requires proving the replacement's divergence is not
//! itself a defect tracked elsewhere — otherwise the guard is pinned on a bug
//! and retires the day it is fixed.
//!
//! So the denotational tests below call
//! [`EquivalenceChecker::compare_denotations`], which runs the same rung logic
//! over the pair **as written** and never asks whether the two normalize alike.
//! The precondition is then *built* rather than *found*: the two forms differ
//! because this file spells them differently, which no normalization change can
//! undo. `check`'s own ladder is still pinned, by
//! `check_stops_at_normalized_match_and_so_cannot_reach_a_denotational_rung` —
//! which asserts the short-circuit that made the old strategy fragile, and is
//! the one test here that *wants* a converging pair.
//!
//! **Re-pinning to `NormalizedMatch` was never an option.**
//! `rulings[confluence-gate-is-apply-equality-on-every-determined-axis]` holds
//! that equivalence is apply-equality on every determined axis and never
//! `NormalizedMatch` — "a relation over `apply`, whose codomain is bases, so it
//! cannot collapse into the normalizer it gates".

use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel, EquivalenceResult};
use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider};

// ---------------------------------------------------------------------------
// A two-exon fixture, laid out so a duplication straddles the junction
// ---------------------------------------------------------------------------

const CONTIG: &str = "chr_two_exon";
const TX: &str = "NM_TWOEXON.1";

/// Exon 1 holds transcript bases 1..=10 at genomic 1001..=1010;
/// exon 2 holds transcript bases 11..=79 at genomic 2001..=2069 (69 bases).
///
/// Exon 1 ends in `T` and exon 2 begins with `T`, so transcript bases 10 and 11
/// are a two-base run that the junction splits. Duplicating either one yields
/// the same transcript sequence and two different genomic sequences, which is
/// the whole point of the fixture.
///
/// The rest of exon 2 is the `reported_partition_verdicts` contig, so
/// #1420-v2's partition pair — `[37dup;41del]` against `38_41delinsATTG` — can
/// be written on this transcript. Transcript position = that contig's position
/// + 11.
///
/// **What that pair is used for here has changed (#1800).** It is a *partition*
/// of one edit: two spellings whose SPDI triple lists have different shapes (two
/// triples against one) and whose reconstructed windows agree. That is what
/// makes it a real exercise of the sequence rung and of the genomic
/// re-derivation above it — a pair spelled identically would compare a triple
/// list against itself and measure nothing. Whether the *normalizer* keeps the
/// two apart is no longer load-bearing: these tests never normalize.
///
/// Two earlier readings of this fixture are recorded so they are not
/// re-attempted. #1419's `[19_23del;27_33del]` / `19_33delinsCGG` stood here
/// until #1649 converged it, and #1420-v2 replaced it on the strength of
/// staying divergent — which #1616 then ended. Divergence is not a property to
/// build on; see the module header.
const EXON1_TX: &str = "ACGCACGCAT";
const EXON2_TX: &str = "TATGCACCAGTCACCAGTCTGATGCGGATCACGTGCAATTGCACGTGCAATTGGATCCGATCGTACGTA";

/// Where exon 2's junction base sits on the genome, and the offset from a
/// transcript position in exon 2 to its genomic position.
const EXON2_GENOMIC_START: u64 = 2001;
const EXON2_TX_START: u64 = 11;

fn two_exon_provider() -> MockProvider {
    two_exon_provider_with_contig(Some(CONTIG.to_string()))
}

/// The same fixture with the transcript's contig withheld, so
/// `transcript_triples_to_genomic` cannot name a genomic axis and the second
/// axis comes back `NotComputed`.
///
/// The transcript *sequence* is still served, so the own-axis comparison
/// succeeds — which is the whole point: this isolates "the second axis could
/// not be computed" from "nothing could be computed".
fn two_exon_provider_without_a_contig() -> MockProvider {
    two_exon_provider_with_contig(None)
}

fn two_exon_provider_with_contig(chromosome: Option<String>) -> MockProvider {
    let tx_sequence = format!("{EXON1_TX}{EXON2_TX}");
    assert_eq!(tx_sequence.len(), 79, "transcript length");
    assert_eq!(
        &tx_sequence[9..11],
        "TT",
        "the junction splits a two-base run"
    );

    // Filler, exon 1, intron, exon 2, filler. The intron is a `G` homopolymer
    // so neither exon's terminal `T` can extend into it.
    let mut genomic = String::with_capacity(2160);
    genomic.push_str(&"A".repeat(1000)); // g.1..=1000
    genomic.push_str(EXON1_TX); // g.1001..=1010
    genomic.push_str(&"G".repeat(990)); // g.1011..=2000
    genomic.push_str(EXON2_TX); // g.2001..=2069
    genomic.push_str(&"A".repeat(91)); // trailing filler
    assert_eq!(genomic.len(), 2160);

    let exon2_end_tx = EXON2_TX_START + EXON2_TX.len() as u64 - 1;
    let exon2_end_genomic = EXON2_GENOMIC_START + EXON2_TX.len() as u64 - 1;
    let transcript = Transcript::new(
        TX.to_string(),
        Some("TWOEXON".to_string()),
        Strand::Plus,
        tx_sequence,
        Some(1),
        Some(78),
        vec![
            Exon::with_genomic(1, 1, 10, 1001, 1010),
            Exon::with_genomic(
                2,
                EXON2_TX_START,
                exon2_end_tx,
                EXON2_GENOMIC_START,
                exon2_end_genomic,
            ),
        ],
        chromosome,
        Some(1001),
        Some(exon2_end_genomic),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );

    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(CONTIG, genomic);
    provider.add_transcript(transcript);
    provider
}

fn level(provider: &MockProvider, a: &str, b: &str) -> EquivalenceLevel {
    check(provider, a, b).level
}

fn check(provider: &MockProvider, a: &str, b: &str) -> EquivalenceResult {
    let checker = EquivalenceChecker::new(provider.clone());
    checker
        .check(&parse_hgvs(a).unwrap(), &parse_hgvs(b).unwrap())
        .unwrap()
}

/// Drive the denotational ladder over two **constructed** forms, and assert it
/// was reached over exactly those two.
///
/// This is the whole of #1800's restructure. `compare_denotations` runs the
/// sequence rung and the cross-axis strengthening above it over the pair as
/// written, so the precondition every rung below needs — two distinct
/// descriptions to compare — is supplied by this file rather than borrowed from
/// whatever the normalizer currently declines to converge.
///
/// Both halves of the postcondition matter and neither implies the other:
///
/// * the two spellings must differ, or the comparison is a description against
///   itself and any rung it reports is unearned; and
/// * the result must have been reached over those same two strings, which is
///   what proves nothing normalized them out from under the assertion. A future
///   change routing this entry point back through the normalizer would fail
///   here instead of silently restoring the fragility.
fn denotational(provider: &MockProvider, first: &str, second: &str) -> EquivalenceResult {
    assert_ne!(
        first, second,
        "a denotational rung sits between two distinct descriptions; comparing one with itself \
         measures nothing"
    );
    let checker = EquivalenceChecker::new(provider.clone());
    let result = checker
        .compare_denotations(&parse_hgvs(first).unwrap(), &parse_hgvs(second).unwrap())
        .unwrap();
    assert_eq!(
        (
            result.normalized_first.as_deref(),
            result.normalized_second.as_deref()
        ),
        (Some(first), Some(second)),
        "the verdict was reached over a different pair than the one constructed, so the level \
         below is measuring something other than what it says"
    );
    result
}

// ---------------------------------------------------------------------------
// 1. The counterexample: transcript-equal, genome-different
// ---------------------------------------------------------------------------

/// The reduced `c.3921dup` / `c.3922dup` case. Both duplicate a `T` in the
/// two-base run that the exon junction splits, so they denote one transcript
/// sequence; they duplicate genomic 1010 and genomic 2001 respectively, which
/// are ~990 bp apart in different exons.
///
/// The pair must reach `SequenceMatch` (the transcript axis agrees) and must
/// **not** reach `CrossAxisSequenceMatch` (the genomic axis does not).
///
/// This is the negative control for the cross-axis comparison, so it is driven
/// the same way as the positive one: a `strengthen_across_axes` that stopped
/// consulting the genomic axis and reported `CrossAxisSequenceMatch` for every
/// own-axis agreement would fail *here*, and nowhere else in this file.
#[test]
fn a_junction_straddling_dup_pair_agrees_on_the_transcript_and_not_on_the_genome() {
    let provider = two_exon_provider();
    let result = denotational(
        &provider,
        &format!("{TX}:c.10dup"),
        &format!("{TX}:c.11dup"),
    );
    let observed = result.level;
    assert_eq!(
        observed,
        EquivalenceLevel::SequenceMatch,
        "the transcript axis agrees, so the pair reaches the single-axis rung"
    );
    assert!(
        !observed.is_at_least(EquivalenceLevel::CrossAxisSequenceMatch),
        "the genomic axis disagrees, so the pair must not clear a cross-axis gate"
    );
}

// ---------------------------------------------------------------------------
// 2. The positive case: agreement on every determined axis
// ---------------------------------------------------------------------------

/// #1420-v2's partition pair, written on the transcript: a two-member cis
/// allele against the single `delins` spanning the same block. The two denote
/// one transcript sequence, and they sit wholly inside exon 2, so the genomic
/// axis is a contiguous re-expression of the same edit and agrees there too.
/// The pair therefore clears the gate.
///
/// **The pair is constructed, not normalized (#1800).** It used to be handed to
/// `check`, which required the normalizer to keep the two apart — and #1616
/// merges them, so that reading of the test is gone. What is being exercised is
/// the rung, and the rung's input is two descriptions: `[48dup;52del]` derives
/// two SPDI triples, `49_52delinsATTG` derives one, both are re-expressed on the
/// genome through the exon alignment, and both reconstructions are compared
/// against the reference window. None of that consults the normalizer, so none
/// of it can be undone by making the normalizer more confluent.
#[test]
fn a_pair_agreeing_on_both_determined_axes_reaches_the_cross_axis_rung() {
    let provider = two_exon_provider();
    let result = denotational(
        &provider,
        &format!("{TX}:c.[48dup;52del]"),
        &format!("{TX}:c.49_52delinsATTG"),
    );
    let observed = result.level;
    assert_eq!(observed, EquivalenceLevel::CrossAxisSequenceMatch);
    assert!(observed.is_at_least(EquivalenceLevel::CrossAxisSequenceMatch));
    assert!(
        result
            .notes
            .iter()
            .any(|note| note.contains("every determined axis")),
        "the note must record that the second axis was consulted and agreed, not merely that the \
         own axis did: {:?}",
        result.notes
    );
}

/// The same shape, one base of payload apart, must come back a decided
/// negative.
///
/// Without it, a `same_resulting_sequence` that agreed unconditionally would
/// satisfy every positive assertion in this file. The control is on the
/// transcript rather than the genome so it also covers the derivation the
/// positive above depends on.
#[test]
fn a_partition_that_denotes_different_bases_is_a_decided_negative() {
    let provider = two_exon_provider();
    let result = denotational(
        &provider,
        &format!("{TX}:c.[48dup;52del]"),
        &format!("{TX}:c.49_52delinsATTC"),
    );
    assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
    assert!(result.level.is_decided());
}

/// A genomic description determines exactly one axis — the genome. Apply-
/// equality on it therefore *is* the whole relation, so a genomic pair that
/// matches by resulting sequence reaches the cross-axis rung rather than
/// stopping at the single-axis one. Without this a confluence gate stated as
/// "at least `CrossAxisSequenceMatch`" would reject every `g.` pair, which is
/// most of what a caller has.
#[test]
fn a_genomic_pair_determines_one_axis_and_so_reaches_the_cross_axis_rung() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(CONTIG, format!("{}{EXON2_TX}", "A".repeat(1000)));
    // The same #1420-v2 pair, shifted onto the contig at g.1001 + n.
    //
    // This test used to pin the pair's *normalized* forms, and they were not
    // the two written here: the span spelling re-derived to
    // `g.[1039T>A;1041_1042delinsTG]`, because a genomic axis carries no reading
    // frame and so `delins.md:47`'s payload-coincidence carve-out does not reach
    // it. That asymmetry is `delins-payload-coincidence-carve-out-is-coding-dna-scoped`
    // and is still true; it is simply no longer this test's business, because
    // the rung is now driven over the constructed pair.
    let result = denotational(
        &provider,
        &format!("{CONTIG}:g.[1038dup;1042del]"),
        &format!("{CONTIG}:g.1039_1042delinsATTG"),
    );
    assert_eq!(result.level, EquivalenceLevel::CrossAxisSequenceMatch);
    assert!(
        result
            .notes
            .iter()
            .any(|note| note.contains("the one axis they determine")),
        "a genomic pair reaches the rung by determining one axis, not by corroboration: {:?}",
        result.notes
    );
}

/// The short-circuit that made the old fixture strategy fragile, pinned as
/// behaviour so the restructure above is anchored to something rather than
/// merely asserted in a comment.
///
/// `check` answers `NormalizedMatch` for a converging pair, and
/// `NormalizedMatch` is deliberately off the denotational order — so
/// `is_at_least(CrossAxisSequenceMatch)` *rejects* a pair whose two spellings
/// normalize alike, and no denotational rung is reached for it at all. Every
/// test in this file that went through `check` therefore needed the normalizer
/// to keep its pair apart, which is the dependency #1800 removes.
///
/// The pair here is chosen to converge and to keep converging: `general.md:56`
/// ranks a substitution above a one-base `delins` of the same base, so
/// `g.1005delinsG` canonicalizes onto `g.1005C>G`. This is the one test in the
/// file that wants convergence, and it is the only one a *further* confluence
/// improvement can strengthen rather than break.
#[test]
fn check_stops_at_normalized_match_and_so_cannot_reach_a_denotational_rung() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(CONTIG, format!("{}{EXON2_TX}", "A".repeat(1000)));

    let converged = check(
        &provider,
        &format!("{CONTIG}:g.1005C>G"),
        &format!("{CONTIG}:g.1005delinsG"),
    );
    assert_eq!(converged.level, EquivalenceLevel::NormalizedMatch);
    assert!(
        !converged
            .level
            .is_at_least(EquivalenceLevel::CrossAxisSequenceMatch),
        "the whole hazard in one assertion: `check` reports a rung that a confluence gate must \
         reject, for a pair that is genuinely one variant"
    );

    // The same pair, same provider, through the denotational entry point: the
    // rung is reached and answers.
    let denoted = denotational(
        &provider,
        &format!("{CONTIG}:g.1005C>G"),
        &format!("{CONTIG}:g.1005delinsG"),
    );
    assert_eq!(denoted.level, EquivalenceLevel::CrossAxisSequenceMatch);
    assert!(denoted
        .level
        .is_at_least(EquivalenceLevel::CrossAxisSequenceMatch));
}

// ---------------------------------------------------------------------------
// 3. Indeterminate: "cannot tell" stops being reported as "no"
// ---------------------------------------------------------------------------

/// A non-literal inserted payload denotes no sequence at all, so neither side
/// can be applied. Today the checker reports `NotEquivalent` — a positive claim
/// that the two descriptions denote different sequences, made without either
/// sequence ever having been computed.
#[test]
fn a_non_literal_payload_is_indeterminate_rather_than_a_negative_verdict() {
    let provider = two_exon_provider();
    let observed = level(
        &provider,
        &format!("{CONTIG}:g.1001_1002insN[10]"),
        &format!("{CONTIG}:g.1001_1002insN[12]"),
    );
    assert_eq!(observed, EquivalenceLevel::Indeterminate);
    assert!(
        !observed.is_decided(),
        "undecidable is not a decided verdict"
    );
    assert!(
        !observed.is_equivalent(),
        "undecidable is not a positive one"
    );
}

/// `ins6` — a count of unspecified bases — is the same shape spelled without a
/// repeat unit.
#[test]
fn an_unspecified_insertion_count_is_indeterminate() {
    let provider = two_exon_provider();
    assert_eq!(
        level(
            &provider,
            &format!("{CONTIG}:g.1001_1002ins6"),
            &format!("{CONTIG}:g.1001_1002ins8"),
        ),
        EquivalenceLevel::Indeterminate
    );
}

/// No input that denotes no sequence may reach a **positive denotational**
/// rung, by whichever path — and through whichever **entry point** — the
/// denotational rungs are reached.
///
/// This is stated as an invariant over the rungs rather than as a regression
/// for one measured pair, because the hazard is precisely that the rule was
/// applied on one path and not another. There are two ways to report agreement
/// — the `SequenceVerdict::Same` branch and the fall-through — and
/// normalization *expands* `insN[10]` into a literal `N` run, so an
/// indeterminate payload arrives at the sequence rung looking definite.
/// Asserting the property over a table, rather than pinning one path's answer,
/// is what makes a future third path fail here instead of silently reporting a
/// positive.
///
/// **Both entry points are driven, and that is not redundant.** `#1800` added
/// `compare_denotations`, which reaches the same rungs without normalizing —
/// and therefore without `check`'s string-identity rung in front of it. A pair
/// of byte-equal indeterminate payloads reconstructs to byte-equal windows, so
/// the sequence rung answers `Same` and returns *before* the fall-through
/// indeterminacy check: `g.1001_1002insNNN` against itself reported
/// `CrossAxisSequenceMatch`, which `is_at_least(CrossAxisSequenceMatch)`
/// accepts. `check` could never show it, because string identity answers
/// `Identical` first. Driving one entry point would leave the other's copy of
/// this invariant unmeasured.
///
/// Every row below is apply-*shaped* like a match — same accession, overlapping
/// spans, payloads that expand to comparable runs — so a checker that consulted
/// only the normalized forms would have something to agree about.
#[test]
fn an_indeterminate_input_never_wins_a_decided_denotational_rung() {
    let provider = two_exon_provider();
    for (a, b) in [
        // Two non-literal repeat payloads at one position.
        (
            format!("{CONTIG}:g.1001_1002insN[10]"),
            format!("{CONTIG}:g.1001_1002insN[12]"),
        ),
        // The same payload expanded by hand on one side, so the two normalized
        // forms are byte-comparable runs of `N`.
        (
            format!("{CONTIG}:g.1001_1002insN[3]"),
            format!("{CONTIG}:g.1002_1003insNNN"),
        ),
        // An unspecified count, which never expands at all.
        (
            format!("{CONTIG}:g.1005_1006ins6"),
            format!("{CONTIG}:g.1005delinsA6"),
        ),
        // On the transcript, where a second determined axis exists and so a
        // second chance to report a cross-axis positive.
        (
            format!("{TX}:c.20_21insN[3]"),
            format!("{TX}:c.20delinsANNN"),
        ),
        // Two byte-equal literal `N` payloads: the sequence rung reconstructs
        // identical windows and answers `Same`, so this is the row that reaches
        // the rung *without* needing normalization to expand anything. It is
        // invisible through `check` — string identity answers `Identical` — and
        // is the row that caught `compare_denotations` skipping the hoisted
        // indeterminacy guard.
        (
            format!("{CONTIG}:g.1001_1002insNNN"),
            format!("{CONTIG}:g.1001_1002insNNN"),
        ),
    ] {
        // `compare_denotations` is called directly rather than through
        // `denotational`, because that helper refuses an identical pair — which
        // is exactly the row that broke this invariant, so it must not be the
        // one row exempted from it.
        let via_denotations = EquivalenceChecker::new(provider.clone())
            .compare_denotations(&parse_hgvs(&a).unwrap(), &parse_hgvs(&b).unwrap())
            .unwrap()
            .level;
        // `check` short-circuits on string identity, so it has nothing to say
        // about the identical row; everywhere it does answer, the two entry
        // points must agree, or the invariant holds on one and not the other.
        let observed = if a == b {
            via_denotations
        } else {
            let via_check = level(&provider, &a, &b);
            assert_eq!(
                via_check, via_denotations,
                "{a} vs {b}: the two entry points disagree about an indeterminate pair, so the \
                 invariant holds on one and not the other"
            );
            via_check
        };
        assert!(
            !observed.is_at_least(EquivalenceLevel::SequenceMatch),
            "{a} vs {b}: reported {observed:?}, a positive denotational verdict about a pair \
             where at least one side denotes no sequence"
        );
        // …and the *exact* rung, not merely "not positive". `NotEquivalent` is
        // also off the denotational order, so the assertion above alone would
        // pass on a regression that answered a decided **negative** for every
        // row — which is precisely the conflation `Indeterminate` was added to
        // remove. Pinning the level is what makes this test able to see it.
        assert_eq!(
            observed,
            EquivalenceLevel::Indeterminate,
            "{a} vs {b}: an undecidable pair must be reported undecidable, not as a decided \
             negative"
        );
        assert!(
            !observed.is_decided(),
            "{a} vs {b}: reported {observed:?}, which `is_decided()` calls true"
        );
        assert!(
            !observed.is_equivalent(),
            "{a} vs {b}: reported {observed:?}, which `is_equivalent()` calls true"
        );
    }
}

/// A trans allele denotes two molecules, so it denotes no single sequence that
/// `apply` could produce. That is undecidable for this relation, not negative.
#[test]
fn a_multi_molecule_allele_is_indeterminate() {
    let provider = two_exon_provider();
    assert_eq!(
        level(
            &provider,
            &format!("{TX}:c.[13T>A];[14G>T]"),
            &format!("{TX}:c.[13T>G];[14G>T]"),
        ),
        EquivalenceLevel::Indeterminate
    );
}

/// The reference window the sequence rung reconstructs is capped at 100 kb. The
/// cap's own comment claims declining past it "only ever leaves the pre-existing
/// `NotEquivalent` verdict in place" — which is exactly the conflation this rung
/// exists to remove, and is wrong for a wide cis allele whose two spellings
/// really are one variant.
///
/// **The contig is a period-4 repeat, not a homopolymer, and that is the whole
/// design of the test.** On 200,000 `A` bases — what this test used to use —
/// each single-base deletion 3'-shuffles as far right as the *normalizer's*
/// fetch window allows, so the verdict was a function of
/// `NormalizeConfig::window_size` rather than of the cap the test names. That
/// is not a theoretical objection: #1697 grew that fetch window, and this test
/// was the case that surfaced the resulting non-idempotent normalization. On
/// `ACGT` repeats every `T` is isolated, so neither deletion can move at any
/// window size, the union span is the only free variable, and
/// `MAX_SEQUENCE_COMPARE_WINDOW` is the only constant that can decide it.
#[test]
fn a_span_past_the_compare_window_cap_is_indeterminate() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(CONTIG, "ACGT".repeat(50_000));

    let far = check(
        &provider,
        &format!("{CONTIG}:g.100del"),
        &format!("{CONTIG}:g.150000del"),
    );
    // Neither spelling moved, so the 149,901-base union span is what exceeds
    // the cap — not a shuffle that walked one of them somewhere else. This test
    // stays on `check` deliberately: the cap is reached through normalization,
    // and pinning the normalized pair is what rules out the alternative
    // explanation.
    assert_eq!(
        (
            far.normalized_first.as_deref(),
            far.normalized_second.as_deref()
        ),
        (
            Some(format!("{CONTIG}:g.100del").as_str()),
            Some(format!("{CONTIG}:g.150000del").as_str())
        ),
        "one of the two deletions shuffled, so the union span is not the variable under test"
    );
    assert_eq!(far.level, EquivalenceLevel::Indeterminate);
    assert!(!far.level.is_decided());

    // The control that makes the span the only variable: the same two shapes on
    // the same contig, close enough to fit inside the cap, are compared and get
    // a decided answer. Without this, an `Indeterminate` from any other cause
    // would read as a pass.
    let near = check(
        &provider,
        &format!("{CONTIG}:g.100del"),
        &format!("{CONTIG}:g.104del"),
    );
    assert_eq!(near.level, EquivalenceLevel::NotEquivalent);
    assert!(near.level.is_decided());
}

/// A second axis that could not be **computed** stops at `SequenceMatch`, and
/// that verdict stays **decided**.
///
/// This is the deliberate line between the two undecidables, and it is pinned
/// here because it is a choice rather than a fallout. `Indeterminate` means *no
/// denotation could be computed for a side* — nothing to compare. Here both
/// sides denote, the transcript comparison ran, and it agreed; what is missing
/// is the second axis's *corroboration*. Reporting `Indeterminate` would
/// discard a fact that was measured, and would answer "cannot tell" about a
/// pair whose own-axis equality is established.
///
/// Nothing is over-claimed by keeping it decided, because the **level** is what
/// a gate reads: a confluence gate is written `is_at_least(
/// CrossAxisSequenceMatch)`, and `SequenceMatch` fails it. The pair is
/// correctly excluded from the gate either way; the only difference is whether
/// the true own-axis fact survives in the answer.
///
/// What the level does *not* carry is the difference between "the second axis
/// disagreed" and "the second axis could not be asked" — that lives in the note,
/// asserted below so it cannot be dropped silently. Splitting it into the level
/// needs a rung of its own rather than a reuse of `Indeterminate`; until a
/// consumer is measured needing it, the note is the record.
#[test]
fn a_second_axis_that_cannot_be_computed_stops_at_the_single_axis_rung() {
    let provider = two_exon_provider_without_a_contig();
    let result = denotational(
        &provider,
        &format!("{TX}:c.[48dup;52del]"),
        &format!("{TX}:c.49_52delinsATTG"),
    );
    assert_eq!(
        result.level,
        EquivalenceLevel::SequenceMatch,
        "own-axis equality held and was measured, so the verdict is that and not `Indeterminate`"
    );
    assert!(
        result.level.is_decided(),
        "a measured own-axis agreement is a decided verdict"
    );
    assert!(
        !result
            .level
            .is_at_least(EquivalenceLevel::CrossAxisSequenceMatch),
        "an uncorroborated pair must not clear a cross-axis gate"
    );
    assert!(
        result
            .notes
            .iter()
            .any(|note| note.contains("could not be computed")),
        "the note is the only place the uncomparable second axis is recorded, so it must say \
         so: {:?}",
        result.notes
    );

    // The control: with the contig restored, the very same pair reaches the
    // cross-axis rung. Without it, a `SequenceMatch` from any other cause would
    // read as a pass.
    assert_eq!(
        denotational(
            &two_exon_provider(),
            &format!("{TX}:c.[48dup;52del]"),
            &format!("{TX}:c.49_52delinsATTG"),
        )
        .level,
        EquivalenceLevel::CrossAxisSequenceMatch
    );
}

// ---------------------------------------------------------------------------
// 4. The order, and what is deliberately not on it
// ---------------------------------------------------------------------------

#[test]
fn the_denotational_rungs_are_ordered_by_strength() {
    use EquivalenceLevel::*;
    assert!(Identical.is_at_least(CrossAxisSequenceMatch));
    assert!(Identical.is_at_least(SequenceMatch));
    assert!(CrossAxisSequenceMatch.is_at_least(SequenceMatch));
    assert!(!CrossAxisSequenceMatch.is_at_least(Identical));
    assert!(!SequenceMatch.is_at_least(CrossAxisSequenceMatch));
    // Reflexive on the order.
    assert!(SequenceMatch.is_at_least(SequenceMatch));
    assert!(CrossAxisSequenceMatch.is_at_least(CrossAxisSequenceMatch));
}

/// `NormalizedMatch` consults the normalizer, so gating on it would restore the
/// circularity the whole relation exists to remove. It is off the order in both
/// directions: it never satisfies a floor, and it can never *be* a floor.
#[test]
fn normalized_match_can_neither_satisfy_nor_be_a_gate() {
    use EquivalenceLevel::*;
    assert!(!NormalizedMatch.is_at_least(SequenceMatch));
    assert!(!Identical.is_at_least(NormalizedMatch));
    assert!(!CrossAxisSequenceMatch.is_at_least(NormalizedMatch));
    assert!(!NormalizedMatch.is_at_least(NormalizedMatch));
}

/// `AccessionVersionDifference` is orthogonal rather than stronger or weaker:
/// it is a claim about two *different* reference sequences, and the relation's
/// first conjunct is "same accession". Apply-equality is not even defined
/// across two versions, so the rung is off the order — while staying a
/// perfectly good, decided, positive verdict.
#[test]
fn the_accession_version_rung_is_off_the_order_but_still_decided() {
    use EquivalenceLevel::*;
    assert!(!AccessionVersionDifference.is_at_least(SequenceMatch));
    assert!(!Identical.is_at_least(AccessionVersionDifference));
    assert!(AccessionVersionDifference.is_equivalent());
    assert!(AccessionVersionDifference.is_decided());
}

// ---------------------------------------------------------------------------
// 5. The same counterexample at the spec's own coordinates, on the minus strand
// ---------------------------------------------------------------------------

const DMD_TX: &str = "NM_004006.2";
const DMD_CONTIG: &str = "NC_000023.11";

/// The reduced fixture above is on the **plus** strand, where `tx_to_genomic`
/// adds the in-exon offset. The spec's actual case is `NM_004006.2`, which is on
/// the **minus** strand, where the offset is *subtracted* from the exon's high
/// end — a different arm of the mapper, and the arm on which the junction
/// counterexample was originally observed. So this repeats the pin at the real
/// coordinates rather than trusting the reduced one to cover both.
///
/// Geometry, as `derived_transcript_placements` carries it:
///
/// * exon 27 — cDNA 4031..=4165 at genomic 32441180..=32441314
/// * exon 28 — cDNA 4166..=4315 at genomic 32438241..=32438390
///
/// Minus strand, so cDNA increases as genomic decreases: cDNA 4165 (the **last**
/// base of exon 27) sits at `g.32441180`, and cDNA 4166 (the **first** base of
/// exon 28) at `g.32438390`. The two are 2,790 bp apart, in different exons —
/// which is what `DNA/duplication.md:148` means by "the wrong nucleotide, in the
/// wrong exon".
///
/// `start_codon` 244 puts `c.1` at cDNA 245, so `c.3921` is cDNA 4165 and
/// `c.3922` is cDNA 4166: the spec's own pair, exactly.
///
/// A leading exon carries cDNA 1..=4030 so the transcript is covered end to end;
/// it sits upstream on the genome (higher coordinates, minus strand) and plays
/// no part in the comparison.
struct DmdFixture {
    provider: MockProvider,
}

impl DmdFixture {
    const EXON27_TX: (u64, u64) = (4031, 4165);
    const EXON27_G: (u64, u64) = (32_441_180, 32_441_314);
    const EXON28_TX: (u64, u64) = (4166, 4315);
    const EXON28_G: (u64, u64) = (32_438_241, 32_438_390);
    const LEAD_TX: (u64, u64) = (1, 4030);
    const LEAD_G: (u64, u64) = (32_445_000, 32_449_029);
    const TX_LEN: usize = 4315;
    const CDS_START: u64 = 245;

    fn new() -> Self {
        // A T-free filler, so the only `T` run in the transcript is the one the
        // junction splits.
        let mut tx: Vec<u8> = "ACGCACGC".bytes().cycle().take(Self::TX_LEN).collect();
        tx[4164] = b'T'; // cDNA 4165 — last base of exon 27
        tx[4165] = b'T'; // cDNA 4166 — first base of exon 28
        assert!(
            tx[4163] != b'T' && tx[4166] != b'T',
            "the `T` run must be exactly the two bases the junction splits"
        );

        let tx_sequence = String::from_utf8(tx).expect("ascii");
        let exons = vec![
            Exon::with_genomic(
                1,
                Self::LEAD_TX.0,
                Self::LEAD_TX.1,
                Self::LEAD_G.0,
                Self::LEAD_G.1,
            ),
            Exon::with_genomic(
                27,
                Self::EXON27_TX.0,
                Self::EXON27_TX.1,
                Self::EXON27_G.0,
                Self::EXON27_G.1,
            ),
            Exon::with_genomic(
                28,
                Self::EXON28_TX.0,
                Self::EXON28_TX.1,
                Self::EXON28_G.0,
                Self::EXON28_G.1,
            ),
        ];

        // The contig: `C` filler, with each exon's genomic window written as the
        // reverse complement of the cDNA segment it carries.
        //
        // The filler must be `C`, not `A`. Read in transcript orientation an
        // intronic genomic `A` is a `T`, which would splice the intron onto the
        // two-base `T` run the junction splits — and the normalizer then shifts
        // `c.3921dup` into the intron as `c.3922-1dup`, which has no SPDI and so
        // reports `NotEquivalent` instead of the single-axis match. `C` reads as
        // `G`, so the run stays exactly two bases. (The plus-strand fixture above
        // uses a `G` intron for the same reason; on the minus strand the choice
        // is complemented, which is easy to get backwards.)
        let mut genomic = vec![b'C'; (Self::LEAD_G.1 + 1_000) as usize];
        for (tx_span, g_span) in [
            (Self::LEAD_TX, Self::LEAD_G),
            (Self::EXON27_TX, Self::EXON27_G),
            (Self::EXON28_TX, Self::EXON28_G),
        ] {
            for tx_pos in tx_span.0..=tx_span.1 {
                let g_pos = g_span.1 - (tx_pos - tx_span.0);
                let base = tx_sequence.as_bytes()[(tx_pos - 1) as usize];
                genomic[(g_pos - 1) as usize] = complement(base);
            }
        }

        let transcript = Transcript::new(
            DMD_TX.to_string(),
            Some("DMD".to_string()),
            Strand::Minus,
            tx_sequence,
            Some(Self::CDS_START),
            Some(Self::TX_LEN as u64),
            exons,
            Some(DMD_CONTIG.to_string()),
            Some(Self::EXON28_G.0),
            Some(Self::LEAD_G.1),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );

        let mut provider = MockProvider::new();
        provider.add_genomic_sequence(DMD_CONTIG, String::from_utf8(genomic).expect("ascii"));
        provider.add_transcript(transcript);
        Self { provider }
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// `c.3921dup` and `c.3922dup` on the real minus-strand geometry: apply-equal on
/// the transcript, and *not* on the genome.
///
/// Driven over the constructed pair for the same reason as its plus-strand twin
/// — and here the hazard is not hypothetical:
/// `rulings[exon-junction-dup-converge-from-the-far-side]` is **decided** and
/// says `c.3922dup` normalizes to `c.3921dup`, so the day it is implemented a
/// `check`-driven version of this test would answer `NormalizedMatch` and stop
/// exercising the mapper's minus-strand arm without saying so.
#[test]
fn the_spec_pair_at_its_own_coordinates_is_single_axis_only() {
    let fixture = DmdFixture::new();
    let observed = denotational(
        &fixture.provider,
        &format!("{DMD_TX}:c.3921dup"),
        &format!("{DMD_TX}:c.3922dup"),
    )
    .level;
    assert_eq!(
        observed,
        EquivalenceLevel::SequenceMatch,
        "both positions carry `T`, so the transcript axis agrees"
    );
    assert!(
        !observed.is_at_least(EquivalenceLevel::CrossAxisSequenceMatch),
        "g.32441180 against g.32438390 — 2,790 bp apart, in different exons \
         (`DNA/duplication.md:148`), so the pair must not clear a cross-axis gate"
    );
}

#[test]
fn every_rung_but_indeterminate_is_decided() {
    use EquivalenceLevel::*;
    for level in [
        Identical,
        NormalizedMatch,
        CrossAxisSequenceMatch,
        SequenceMatch,
        AccessionVersionDifference,
        NotEquivalent,
    ] {
        assert!(level.is_decided(), "{level:?} must be a decided verdict");
    }
    assert!(!Indeterminate.is_decided());
    assert!(!Indeterminate.is_equivalent());
}
