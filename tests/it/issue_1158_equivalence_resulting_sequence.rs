//! Regression test for issue #1158: `EquivalenceChecker::check` reported
//! `NotEquivalent` for two variants with the **same resulting sequence** when
//! they used different (but equivalent) HGVS encodings of a complex indel.
//!
//! Root cause: the checker only ever compared normalized HGVS **strings**; it
//! never projected either variant to its resulting reference sequence. Two
//! encodings that ferro legitimately keeps in distinct canonical forms — a
//! single length-changing `delins` vs a decomposed cis allele of the same edit
//! — therefore compared unequal even though they produce byte-identical
//! sequence against the reference.
//!
//! Fix: a sequence-level rung (`EquivalenceLevel::SequenceMatch`) that
//! reconstructs each variant's edited window from SPDI triples and compares the
//! resulting sequences. It only ever upgrades a `NotEquivalent` verdict — the
//! Identical / NormalizedMatch / AccessionVersionDifference rungs run first.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};
use ferro_hgvs::parse_hgvs;

fn level(core: &str, a: &str, b: &str) -> EquivalenceLevel {
    let checker = EquivalenceChecker::new(SyntheticBuilder::genomic(core).build());
    checker
        .check(&parse_hgvs(a).unwrap(), &parse_hgvs(b).unwrap())
        .unwrap()
        .level
}

/// Core `AGTCAGT` at HGVS 257..=263. Replacing all seven bases with `GATTA`
/// (`g.257_263delinsGATTA`) yields the same sequence as the decomposed cis
/// allele `[257A>G;258G>A;260C>T;262_263del]`. The two normalize to **different**
/// canonical forms and are caught by the sequence rung. This is issue #1158's
/// case A, and the rung it exercises is the one #1158 asked for.
const CASE_A_CORE: &str = "AGTCAGT";
const CASE_A_DELINS: &str = "NC_TEST.1:g.257_263delinsGATTA";
const CASE_A_ALLELE: &str = "NC_TEST.1:g.[257A>G;258G>A;260C>T;262_263del]";

#[test]
fn delins_and_decomposed_allele_with_same_sequence_are_equivalent() {
    // **Re-blessed back to `SequenceMatch` under the partition model**, which is
    // the rung this file's own header describes: "Two encodings that ferro
    // legitimately keeps in distinct canonical forms — a single length-changing
    // `delins` vs a decomposed cis allele of the same edit". Between #1235 and the
    // partition model the two encodings converged, because a partition was
    // re-derived from `(reference, resulting sequence)` and both spellings gave
    // the same one; the checker then reached `NormalizedMatch` and case A stopped
    // covering the rung it was written for.
    //
    // The two forms are now distinct again, and legitimately so, because the
    // encodings assert different partitions and each is preserved:
    //
    // - `g.257_263delinsGATTA` asserts **one** member. It is length-changing
    //   (7 reference nucleotides against a 5-base payload), so there is no
    //   column-wise correspondence in which an interior unchanged run could be
    //   identified without inventing an alignment, and `general.md:34` — which
    //   speaks of "two variants separated by one or more nucleotides" — has no
    //   purchase on a lone member. It survives whole.
    // - `g.[257A>G;258G>A;260C>T;262_263del]` asserts **four**. The first two abut
    //   (separation 0) and merge unconditionally per `DNA/substitution.md:32`
    //   ("changes involving two or more consecutive nucleotides are described as
    //   deletion-insertion (delins)"); the remaining gaps are one unchanged
    //   nucleotide each, which `general.md:34` keeps individual, so it comes out
    //   as `[257_258delinsGA;260C>T;262_263del]`.
    //
    // Two partitions, two canonical strings, one resulting sequence — which is
    // exactly the shape `SequenceMatch` exists to rate. `general.md:34` cannot
    // adjudicate between them: `DNA/delins.md:79-84` gives the spec's own reason
    // the input's spelling is the authority ("the two variants may have been
    // reported (or might occur) individually"), and `DNA/delins.md:47` recommends
    // the spanning form for the shape the first encoding writes.
    assert_eq!(
        level(CASE_A_CORE, CASE_A_DELINS, CASE_A_ALLELE),
        EquivalenceLevel::SequenceMatch,
    );
}

#[test]
fn sequence_match_is_equivalent() {
    let checker = EquivalenceChecker::new(SyntheticBuilder::genomic(CASE_A_CORE).build());
    let result = checker
        .check(
            &parse_hgvs(CASE_A_DELINS).unwrap(),
            &parse_hgvs(CASE_A_ALLELE).unwrap(),
        )
        .unwrap();
    assert!(result.is_equivalent());
    assert!(EquivalenceLevel::SequenceMatch.is_equivalent());
}

/// Case B from the issue: a `delins` that reduces to a pure deletion vs the
/// direct deletion of the same bases. Core `TAAAAAG`: `g.257_259delinsT`
/// (delete `TAA`, insert `T`) reduces to deleting `AA`, the same edit as
/// `g.258_259del`.
///
/// The issue's requirement is that the checker report these **equivalent**; it
/// does not fix which rung gets there, and the answer legitimately moved. Before
/// the companion normalizer fix (#1157) the reduced deletion skipped the
/// 3'-shift, so the two spellings normalized differently and only this PR's
/// sequence-level rung unified them. #1157 now shifts the reduction, so they
/// normalize identically and the earlier `NormalizedMatch` rung — which runs
/// first by design — catches them. Assert the contract (equivalent), not the
/// rung, so neither fix's landing order can flip this test red. Case A above
/// stays non-confluent by design and remains the sequence rung's own coverage.
#[test]
fn delins_reducing_to_deletion_matches_direct_deletion() {
    let level = level(
        "TAAAAAG",
        "NC_TEST.1:g.257_259delinsT",
        "NC_TEST.1:g.258_259del",
    );
    assert!(
        level.is_equivalent(),
        "a delins reducing to a deletion must be equivalent to the direct \
         deletion of the same bases, got {level:?}"
    );
}

/// The rung must not over-match: two delins that differ in a single inserted
/// base produce different sequences and must stay `NotEquivalent`.
#[test]
fn delins_with_different_resulting_sequence_stays_not_equivalent() {
    assert_eq!(
        level(
            CASE_A_CORE,
            "NC_TEST.1:g.257_263delinsGATTA",
            "NC_TEST.1:g.257_263delinsGATTC",
        ),
        EquivalenceLevel::NotEquivalent,
    );
}

/// A decomposed allele whose edits yield a *different* sequence than the delins
/// must also stay `NotEquivalent` (guards the allele reconstruction path).
#[test]
fn decomposed_allele_with_different_sequence_stays_not_equivalent() {
    // Same as CASE_A_ALLELE but the first substitution is 257A>C (not A>G),
    // so the resulting first base is C, not G — a different sequence.
    assert_eq!(
        level(
            CASE_A_CORE,
            CASE_A_DELINS,
            "NC_TEST.1:g.[257A>C;258G>A;260C>T;262_263del]",
        ),
        EquivalenceLevel::NotEquivalent,
    );
}

// -- Existing behavior must be preserved (the new rung runs only after these) --

#[test]
fn identical_variants_still_identical() {
    assert_eq!(
        level(CASE_A_CORE, "NC_TEST.1:g.258A>G", "NC_TEST.1:g.258A>G"),
        EquivalenceLevel::Identical,
    );
}

#[test]
fn substitution_written_two_ways_still_normalized_match() {
    // A substitution written as a 1-base delins normalizes to the same string,
    // so it is caught by the NormalizedMatch rung before sequence comparison.
    assert_eq!(
        level(CASE_A_CORE, "NC_TEST.1:g.257A>G", "NC_TEST.1:g.257delinsG"),
        EquivalenceLevel::NormalizedMatch,
    );
}
