//! #1578 follow-up (item 4): the equivalence checker must not return a verdict
//! about a variant that **no rung ever examined**.
//!
//! `EquivalenceChecker::check` decides in rungs: string identity, then
//! accession-version difference, then `NormalizedMatch` on the normalized
//! forms, then `SequenceMatch` on the SPDI-reconstructed edited sequence
//! (#1158). `NotEquivalent` is the fall-through.
//!
//! For a genome ring — and for the `sup`-wrapped ring — **both** of the two
//! deciding rungs are structurally blind:
//!
//! - `Normalizer::normalize` returns a ring unchanged (`src/normalize/mod.rs`,
//!   the `HV::GenomeRing` / `HV::Supernumerary` pass-through arms), so
//!   `NormalizedMatch` degenerates into the string comparison the `Identical`
//!   rung already made; and
//! - `hgvs_to_spdi` declines a ring (`UnsupportedVariantType`), so
//!   `edit_triples` yields `None` and the `SequenceMatch` rung cannot fire
//!   either.
//!
//! So `check` fell through to `NotEquivalent` — a *positive claim that two
//! descriptions denote different variants* — having evaluated neither. That is
//! strictly worse than declining: the other three modules that hand-roll
//! `Allele` recursion (`src/spdi/apply.rs`, `src/vcf/from_hgvs.rs`,
//! `src/project/projector.rs`) all decline a ring, and their declines are
//! pinned in `issue_1578_followup_ring_declines.rs`.
//!
//! `Circular` (`o.`) is the near miss that shows the defect is the *conjunction*
//! of the two blindnesses rather than the pass-through alone: `o.` also passes
//! through normalization untouched, but `hgvs_to_spdi` supports it, so the
//! `SequenceMatch` rung still reaches a sound verdict. That case is pinned
//! below so a future change to either rung cannot quietly move it into the
//! declining set.
//!
//! The two textual rungs are sound for a ring and are deliberately kept:
//! `Identical` and `AccessionVersionDifference` compare strings and need no
//! normalization. Only the fall-through is replaced by a decline.

use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::MockProvider;

/// A contig with two homopolymer runs, so that a single-base `del` inside a run
/// has a 3'-most normalized form that differs from the spelling given:
///
/// - 1-based 100..109 is `AAAAAAAAAA`, so `g.100del`, `g.105del` and `g.109del`
///   all normalize to `g.109del`;
/// - 1-based 200..209 is `TTTTTTTTTT`, so `g.200del` and `g.205del` both
///   normalize to `g.209del`.
///
/// That is what makes the ring rows below *the same ring spelled two ways*
/// rather than two different rings — the property the checker got wrong.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let mut sequence = String::new();
    sequence.push_str(&"C".repeat(99)); // 1..99
    sequence.push_str("AAAAAAAAAA"); // 100..109
    sequence.push_str(&"G".repeat(90)); // 110..199
    sequence.push_str("TTTTTTTTTT"); // 200..209
    sequence.push_str(&"C".repeat(291)); // 210..500
    provider.add_genomic_sequence("NC_000022.11", &sequence);
    provider.add_genomic_sequence("NC_000022.10", &sequence);
    provider.add_genomic_sequence("NC_001422.1", &sequence);
    provider
}

/// The control that establishes the two ring rows below describe one variant:
/// spelled as bare `g.` descriptions, the very same two positions converge
/// under normalization. Any change that breaks this control invalidates the
/// ring assertions rather than merely failing alongside them.
#[test]
fn the_bare_spellings_of_the_ring_segments_are_equivalent() {
    let checker = EquivalenceChecker::new(provider());
    for (left, right) in [
        ("NC_000022.11:g.100del", "NC_000022.11:g.105del"),
        ("NC_000022.11:g.200del", "NC_000022.11:g.205del"),
    ] {
        let result = checker
            .check(&parse_hgvs(left).unwrap(), &parse_hgvs(right).unwrap())
            .expect("a bare genomic pair is decidable");
        assert_eq!(
            result.level,
            EquivalenceLevel::NormalizedMatch,
            "`{left}` and `{right}` must converge under normalization — \
             the ring assertions below depend on it",
        );
    }
}

/// The defect. Two spellings of one ring differ only in a segment's
/// pre-normalization position, and the checker used to answer `NotEquivalent`:
/// a claim that they denote different variants, reached without evaluating
/// either side. It must decline instead.
#[test]
fn a_ring_pair_neither_rung_can_evaluate_is_declined() {
    let checker = EquivalenceChecker::new(provider());
    for (left, right) in [
        // the first segment spelled two ways
        (
            "NC_000022.11:g.100del::200del",
            "NC_000022.11:g.105del::200del",
        ),
        // the second segment spelled two ways
        (
            "NC_000022.11:g.100del::200del",
            "NC_000022.11:g.100del::205del",
        ),
        // and the `sup`-wrapped ring, which is blind through the same arms
        (
            "NC_000022.11:g.[100del::200del]sup",
            "NC_000022.11:g.[105del::200del]sup",
        ),
    ] {
        let error = checker
            .check(&parse_hgvs(left).unwrap(), &parse_hgvs(right).unwrap())
            .expect_err(&format!(
                "`{left}` vs `{right}`: no rung can evaluate a ring, so the \
                 checker must decline rather than return a verdict"
            ));
        let rendered = error.to_string();
        assert!(
            rendered.contains("equivalence"),
            "the decline must say it is an equivalence verdict it withheld, \
             got: {rendered}"
        );
    }
}

/// `Identical` is decided by string comparison before any rung that needs a
/// provider, so it stays sound for a ring and must keep working. A decline here
/// would be a regression in the opposite direction — refusing a question that
/// has an obvious answer.
#[test]
fn two_identical_rings_are_still_identical() {
    let checker = EquivalenceChecker::new(provider());
    for input in [
        "NC_000022.11:g.100del::200del",
        "NC_000022.11:g.[100del::200del]sup",
    ] {
        let variant = parse_hgvs(input).unwrap();
        let result = checker
            .check(&variant, &variant)
            .expect("string identity needs no normalization");
        assert_eq!(
            result.level,
            EquivalenceLevel::Identical,
            "`{input}` compared with itself is `Identical`",
        );
    }
}

/// `AccessionVersionDifference` is likewise purely textual — it compares the
/// accession's prefix/number/version and then the rendered variant part — so it
/// remains sound for a ring and must not be swept into the decline.
#[test]
fn one_ring_on_two_accession_versions_is_still_a_version_difference() {
    let checker = EquivalenceChecker::new(provider());
    let result = checker
        .check(
            &parse_hgvs("NC_000022.10:g.100del::200del").unwrap(),
            &parse_hgvs("NC_000022.11:g.100del::200del").unwrap(),
        )
        .expect("an accession-version difference is decided textually");
    assert_eq!(
        result.level,
        EquivalenceLevel::AccessionVersionDifference,
        "the same ring on two versions is a version difference, not a decline",
    );
}

/// Two rings on plainly disjoint spans are **also** declined, and that is the
/// deliberate reading rather than an over-broad fix.
///
/// It is tempting to answer `NotEquivalent` here, since the two are obviously
/// different. But the checker has no way to distinguish this pair from the
/// respelled pair above: in both cases the ring reaches the fall-through with
/// its segments un-normalized and SPDI-unrepresentable, so the *same* evidence
/// (none) underlies both. Answering `NotEquivalent` on this pair would be right
/// by luck while being wrong on the other, and a rung that is right by luck is
/// the thing this change removes. Declining on both is the position the
/// available evidence actually supports.
///
/// The cost is bounded and visible: a caller that needs a verdict for a ring
/// gets an explicit decline naming why, not a plausible wrong answer.
#[test]
fn two_rings_on_disjoint_spans_are_declined_rather_than_answered_by_luck() {
    let checker = EquivalenceChecker::new(provider());
    let error = checker
        .check(
            &parse_hgvs("NC_000022.11:g.100del::200del").unwrap(),
            &parse_hgvs("NC_000022.11:g.400del::450del").unwrap(),
        )
        .expect_err("the checker cannot tell this pair from a respelled one");
    assert!(
        error.to_string().contains("cannot decide equivalence"),
        "the decline must name itself: {error}"
    );
}

/// The scope boundary, pinned. An RNA fusion (`::` between two transcripts)
/// shares both blindnesses — normalization passes it through and SPDI cannot
/// represent it — so the same argument would reach it. It is deliberately left
/// answering, because no pair of fusion spellings denoting *one* fusion has been
/// exhibited: widening the decline there would refuse pairs the checker answers
/// today without evidence that any of those answers is wrong.
///
/// This test exists so that boundary is a recorded decision rather than an
/// oversight. If a wrong fusion answer is ever measured, this test is the one to
/// change, and changing it should require stating the pair that measures it.
#[test]
fn an_rna_fusion_pair_is_still_answered_a_deliberate_scope_boundary() {
    let checker = EquivalenceChecker::new(provider());
    let result = checker
        .check(
            &parse_hgvs("NM_000088.3:r.1_100::NM_000089.3:r.200_300").unwrap(),
            &parse_hgvs("NM_000088.3:r.1_100::NM_000089.3:r.201_300").unwrap(),
        )
        .expect("fusions are outside this change's scope and still answer");
    assert_eq!(
        result.level,
        EquivalenceLevel::NotEquivalent,
        "two different fusions still answer `NotEquivalent` — the decline is \
         scoped to the ring shapes a defect was measured on",
    );
}

/// The near miss that defines the boundary. A circular (`o.`) variant also
/// passes through normalization unchanged, but `hgvs_to_spdi` supports it, so
/// the `SequenceMatch` rung still evaluates the pair and reaches a sound
/// verdict. It must therefore stay *decidable* — the decline keys on "no rung
/// could examine this", not on "normalization was a pass-through".
#[test]
fn a_circular_pair_is_still_decided_by_the_sequence_rung() {
    let checker = EquivalenceChecker::new(provider());
    let result = checker
        .check(
            &parse_hgvs("NC_001422.1:o.100del").unwrap(),
            &parse_hgvs("NC_001422.1:o.105del").unwrap(),
        )
        .expect("`o.` is SPDI-representable, so the sequence rung can decide it");
    assert_eq!(
        result.level,
        EquivalenceLevel::SequenceMatch,
        "`o.` normalization is a pass-through, but the SPDI rung still decides \
         the pair — this is why the decline must not key on the pass-through",
    );
}
