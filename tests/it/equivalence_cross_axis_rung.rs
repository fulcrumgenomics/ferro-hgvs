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
/// The rest of exon 2 is the `reported_partition_verdicts` contig, so a
/// partition pair that is known to stay in two distinct canonical forms —
/// #1420-v2's `[37dup;41del]` against `38_41delinsATTG` — can be written on
/// this transcript. Transcript position = that contig's position + 11.
///
/// #1419's `[19_23del;27_33del]` / `19_33delinsCGG` was the pair used here
/// until #1649 converged it; do not reach for it again without re-measuring.
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

/// Assert the two normalized forms a verdict was reached over.
///
/// Every denotational rung below sits *between* two normalized strings, so a
/// test that pins only the level cannot tell "the rung fired" from "the two
/// spellings converged and the rung was never reached". That is not a
/// hypothetical: these tests were first written over #1419's partition pair,
/// #1649 converged it, and both of them silently started asserting
/// `NormalizedMatch` against a rung they no longer exercised. Pinning the pair
/// of strings makes that failure name itself.
fn assert_normalized(result: &EquivalenceResult, first: &str, second: &str) {
    assert_eq!(
        (
            result.normalized_first.as_deref(),
            result.normalized_second.as_deref()
        ),
        (Some(first), Some(second)),
        "the normalized forms this verdict was reached over have moved; the level below is \
         measuring something other than what it says"
    );
    assert_ne!(
        first, second,
        "the two spellings converged, so `check` returned at the NormalizedMatch rung and no \
         denotational rung was exercised"
    );
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
#[test]
fn a_junction_straddling_dup_pair_agrees_on_the_transcript_and_not_on_the_genome() {
    let provider = two_exon_provider();
    let result = check(
        &provider,
        &format!("{TX}:c.10dup"),
        &format!("{TX}:c.11dup"),
    );
    // Neither spelling moves: the exon boundary bounds the duplication on the
    // transcript, so the two stay on opposite sides of it and the rung below is
    // reached over a genuinely distinct pair.
    assert_normalized(&result, &format!("{TX}:c.10dup"), &format!("{TX}:c.11dup"));
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

/// #1420-v2's partition pair, written on the transcript. The two spellings stay
/// in distinct canonical forms — that is the open non-confluence the issue is
/// about — so they reach a sequence rung rather than `NormalizedMatch`; and
/// they sit wholly inside exon 2, so the genomic axis is a contiguous
/// re-expression of the same edit and agrees. The pair therefore clears the
/// gate.
///
/// **#1419's pair used to stand here and no longer can.** #1649 gave the
/// splitter the two-deletion shape, which converged
/// `c.[30_34del;38_44del]` with `c.30_44delinsCGG`; both tests then asserted
/// against `NormalizedMatch` without exercising any denotational rung. #1420-v2
/// is the replacement because it is a *dup*-plus-*del* partition, which the
/// splitter does not reach — see `reported_partition_verdicts`, where its two
/// rows are still pinned as divergent. `assert_normalized` is what makes a
/// repeat of that convergence fail loudly instead of quietly.
#[test]
fn a_pair_agreeing_on_both_determined_axes_reaches_the_cross_axis_rung() {
    let provider = two_exon_provider();
    let result = check(
        &provider,
        &format!("{TX}:c.[48dup;52del]"),
        &format!("{TX}:c.49_52delinsATTG"),
    );
    assert_normalized(
        &result,
        &format!("{TX}:c.[48dup;52del]"),
        &format!("{TX}:c.49_52delinsATTG"),
    );
    let observed = result.level;
    assert_eq!(observed, EquivalenceLevel::CrossAxisSequenceMatch);
    assert!(observed.is_at_least(EquivalenceLevel::CrossAxisSequenceMatch));
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
    // The same #1420-v2 pair, shifted onto the contig at g.1001 + n. The span
    // spelling re-derives here — a genomic axis carries no reading frame, so
    // `delins.md:47`'s payload-coincidence carve-out does not reach it and the
    // members stay individual — which is exactly why the pinned pair below is
    // not the pair pinned on the transcript.
    let result = check(
        &provider,
        &format!("{CONTIG}:g.[1038dup;1042del]"),
        &format!("{CONTIG}:g.1039_1042delinsATTG"),
    );
    assert_normalized(
        &result,
        &format!("{CONTIG}:g.[1038dup;1042del]"),
        &format!("{CONTIG}:g.[1039T>A;1041_1042delinsTG]"),
    );
    assert_eq!(result.level, EquivalenceLevel::CrossAxisSequenceMatch);
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
/// rung, by whichever path `check` happens to take.
///
/// This is stated as an invariant over `check` rather than as a regression for
/// one measured pair, because the hazard is precisely that the rule was applied
/// on one path and not another. `check` has two ways to report agreement — the
/// `SequenceVerdict::Same` branch and the fall-through — and normalization
/// *expands* `insN[10]` into a literal `N` run, so an indeterminate payload
/// arrives at the sequence rung looking definite. Asserting the property over a
/// table, rather than pinning one path's answer, is what makes a future third
/// path fail here instead of silently reporting a positive.
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
    ] {
        let observed = level(&provider, &a, &b);
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
    // the cap — not a shuffle that walked one of them somewhere else.
    assert_normalized(
        &far,
        &format!("{CONTIG}:g.100del"),
        &format!("{CONTIG}:g.150000del"),
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
    let result = check(
        &provider,
        &format!("{TX}:c.[48dup;52del]"),
        &format!("{TX}:c.49_52delinsATTG"),
    );
    assert_normalized(
        &result,
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
        level(
            &two_exon_provider(),
            &format!("{TX}:c.[48dup;52del]"),
            &format!("{TX}:c.49_52delinsATTG"),
        ),
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
#[test]
fn the_spec_pair_at_its_own_coordinates_is_single_axis_only() {
    let fixture = DmdFixture::new();
    let observed = level(
        &fixture.provider,
        &format!("{DMD_TX}:c.3921dup"),
        &format!("{DMD_TX}:c.3922dup"),
    );
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
