//! Regression test for issue #1157: normalization of a `delins` that reduces
//! to a pure deletion (or a duplication) was **not idempotent** and **not
//! confluent** with the same edit written directly as a `del`/`dup`.
//!
//! Root cause: in `normalize_na_edit`, when `canonicalize_delins` trimmed a
//! `delins` down to a `DelinsCanonical::Deletion` (or `::Duplication`), the
//! arm returned the reduced edit **directly**, bypassing the 3'-shift
//! `shuffle()` that every genuine `del`/`dup` input receives. The sibling
//! `DelinsCanonical::Insertion` arm already recurses into `normalize_na_edit`
//! for exactly this reason; the del/dup arms did not.
//!
//! Symptom (from the issue):
//!   norm(g.10_20delinsTC) = g.12_20del        (unshifted — WRONG)
//!   norm(g.12_20del)      = g.13_21del         (shifted — the true canonical)
//!   ⇒ norm(norm(x)) != norm(x)   (idempotence violated)
//!
//! These tests exercise the same defect on the synthetic homopolymer/tandem
//! fixtures where the 3'-shift target is unambiguous.

use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};

/// Normalize `input` against a genomic provider whose core sequence is `core`
/// (padded on both sides; the core's first base is at 1-based HGVS position
/// `PAD_OFFSET + 1`).
fn norm(core: &str, input: &str) -> String {
    normalize_to_string(SyntheticBuilder::genomic(core).build(), input)
}

// --------------------------------------------------------------------------
// delins → deletion
// --------------------------------------------------------------------------

/// Core `T AAAAA G`: a `T` at core position 1, an A-tract at core 2..=6, a `G`
/// at core 7. `PAD_OFFSET` = 256, so core position k is HGVS position 256 + k.
///
/// `g.257_259delinsT` deletes ref `TAA` and inserts `T`: shared-affix trimming
/// strips the leading `T`, leaving a pure deletion of `AA` at 258_259 — inside
/// the A-tract 258..=262. A genuine deletion of two A's from that tract 3'-shifts
/// and collapses to the repeat form `g.258_262A[3]`.
const DEL_CORE: &str = "TAAAAAG";

#[test]
fn delins_reducing_to_deletion_is_three_prime_shifted() {
    // The bug returned the unshifted `g.258_259del`; the canonical form is the
    // 3'-shifted repeat contraction.
    assert_eq!(
        norm(DEL_CORE, "NC_TEST.1:g.257_259delinsT"),
        "NC_TEST.1:g.258_262A[3]",
    );
}

#[test]
fn delins_reducing_to_deletion_is_idempotent() {
    let once = norm(DEL_CORE, "NC_TEST.1:g.257_259delinsT");
    let twice = norm(DEL_CORE, &once);
    assert_eq!(once, twice, "norm(norm(x)) must equal norm(x)");
}

#[test]
fn delins_reducing_to_deletion_matches_direct_deletion() {
    // Encoding invariance: the delins spelling and the equivalent direct
    // deletion of the same two A's must normalize to the same canonical form.
    assert_eq!(
        norm(DEL_CORE, "NC_TEST.1:g.257_259delinsT"),
        norm(DEL_CORE, "NC_TEST.1:g.258_259del"),
    );
}

// --------------------------------------------------------------------------
// delins → duplication (sibling defect in the same function)
// --------------------------------------------------------------------------

/// Core `A GGGGG T`: a `G`-tract at core 2..=6. `g.258delinsGG` deletes ref `G`
/// and inserts `GG` — a single-copy duplication of the reference `G` at 258.
/// A genuine duplication inside the G-tract 3'-shifts to the tract's 3' end,
/// `g.262dup`.
const DUP_CORE: &str = "AGGGGGT";

#[test]
fn delins_reducing_to_duplication_is_three_prime_shifted() {
    // The bug returned the unshifted `g.258dup`; the canonical form is `g.262dup`.
    assert_eq!(
        norm(DUP_CORE, "NC_TEST.1:g.258delinsGG"),
        "NC_TEST.1:g.262dup",
    );
}

#[test]
fn delins_reducing_to_duplication_is_idempotent() {
    let once = norm(DUP_CORE, "NC_TEST.1:g.258delinsGG");
    let twice = norm(DUP_CORE, &once);
    assert_eq!(once, twice, "norm(norm(x)) must equal norm(x)");
}

#[test]
fn delins_reducing_to_duplication_matches_direct_duplication() {
    assert_eq!(
        norm(DUP_CORE, "NC_TEST.1:g.258delinsGG"),
        norm(DUP_CORE, "NC_TEST.1:g.258dup"),
    );
}

// --------------------------------------------------------------------------
// Scope boundary: issue #1157 "case A". Two spellings, two canonical forms,
// each selected by its own clause.
// --------------------------------------------------------------------------
//
// A single length-changing `delins` and a decomposed cis allele of the *same*
// edit produce the same resulting sequence and are kept in DISTINCT canonical
// forms. That is not a shrug — a different clause governs each side, and both
// clauses are the spec speaking directly at the shape in front of it:
//
//   * the ALLELE asserts three changed blocks separated by one unchanged
//     nucleotide each (259 and 261). `general.md:34` — "two variants separated
//     by one or more nucleotides should be described individually and **not**
//     as a "delins"" — keeps them individual. `general.md:35`'s codon
//     exception is the only thing that would merge them, and it cannot reach a
//     `g.` description: there is no reading frame, so "together affecting one
//     amino acid" is unanswerable.
//
//   * the DELINS asserts ONE changed block. Splitting it requires inventing an
//     alignment of a 7-nt reference against a 5-nt payload and then reading the
//     columns that happen to agree as unchanged runs. That is exactly the
//     construction `DNA/delins.md:46` performs ("parts of the inserted sequence
//     "align" with the reference sequence, giving an alternative description
//     like `c.[850_869del;874_881del;887_897del;901_902insG]`") and exactly the
//     one `:47` rejects one line later: "**The "delins" format is
//     recommended**: it is simpler and prevents software tools making incorrect
//     predictions for the consequences on protein level."
//
// The ruling record `delins-merge-vs-individual-gap-two-or-more` (DECIDED,
// 2026-08-07) settles that conflict for `:47` and scopes the ruling to precisely
// this class — "a MINIMAL single `delins` that would be split because payload
// bases coincide with reference bases". `g.257_263delinsGATTA` is that shape,
// so the delins is retained.
//
// DISCLOSURE. Under #1235 these two spellings were made to converge on the
// three-member form; they no longer do. That is a representation change and it
// is deliberate: the convergence was bought by splitting the single delins,
// which the ruling above adjudicates as the wrong direction. Nothing here
// claims the two spellings *should* diverge as a matter of taste — it records
// that two rank-(1) clauses reach them separately and disagree, which is the
// spec's own provenance contradiction (`DNA/delins.md:83`: "the two variants
// may have been reported (or might occur) individually").
//
// Both forms are — and must remain — idempotent.
//
// Core `AGTCAGT` at HGVS 257..=263: replacing all seven bases with `GATTA`
// (`g.257_263delinsGATTA`) yields the same sequence as
// `g.[257A>G;258G>A;260C>T;262_263del]`.
const CASE_A_CORE: &str = "AGTCAGT";
const CASE_A_DELINS: &str = "NC_TEST.1:g.257_263delinsGATTA";
const CASE_A_ALLELE: &str = "NC_TEST.1:g.[257A>G;258G>A;260C>T;262_263del]";

#[test]
fn decomposed_cis_allele_is_not_collapsed_into_a_spanning_delins() {
    // Adjacent members 257A>G;258G>A merge to 257_258delinsGA, but the allele
    // keeps three members: the unchanged bases at 259 and 261 are NOT absorbed
    // into a single spanning delins over 257..=263.
    assert_eq!(
        norm(CASE_A_CORE, CASE_A_ALLELE),
        "NC_TEST.1:g.[257_258delinsGA;260C>T;262_263del]",
    );
}

/// NOTE ON THE NAME: it says "decomposed" and the assertion says the opposite.
/// The name is pinned by `msto_regression_corpus::MSTO_ISSUES` (#1157) and
/// `every_cataloged_test_exists_in_its_source_file` fails if it moves, so it is
/// left alone rather than renamed here; renaming it and the catalog entry
/// together is a follow-up. Read the assertion, not the name.
///
/// `g.257_263delinsGATTA` is a MINIMAL single `delins`: one asserted changed
/// block, 7 reference nt replaced by 5. It is retained.
///
/// Under #1235 this decomposed to `g.[257_258delinsGA;260C>T;262_263del]`,
/// citing `delins.md:17` ("two variants separated by one or more nucleotides
/// should be described individually and **not** as a "delins""). That citation
/// does not survive reading three lines further. `delins.md:46` builds an
/// alignment-driven split of its own worked example —
/// `c.[850_869del;874_881del;887_897del;901_902insG]`, whose gaps are 4, 5 and
/// 3, every one of them "one or more" — and `delins.md:47` answers "**The
/// "delins" format is recommended**". So `:17` read as reaching an
/// alignment-driven split demands precisely the description `:47` advises
/// against, in the spec's own example. The ruling record
/// `delins-merge-vs-individual-gap-two-or-more` (DECIDED 2026-08-07) resolves
/// that for `:47`, scoped to "a MINIMAL single `delins` that would be split
/// because payload bases coincide with reference bases".
///
/// This block is that shape. There are no unchanged bases in the input's
/// assertion at all: 259 and 261 are "unchanged" only under one particular
/// alignment of `AGTCAGT` against `GATTA`, and that alignment is not unique —
/// non-uniqueness is ground (4) of the ruling.
#[test]
fn single_delins_is_decomposed_at_its_unchanged_bases() {
    assert_eq!(norm(CASE_A_CORE, CASE_A_DELINS), CASE_A_DELINS);
    // …and it is a fixed point, so this is a canonical form rather than a
    // pass-through that a second normalization would move.
    assert_eq!(
        norm(CASE_A_CORE, &norm(CASE_A_CORE, CASE_A_DELINS)),
        CASE_A_DELINS,
    );
}

/// NOTE ON THE NAME: as above — it says "normalize_equal" and the assertion
/// says they differ. The name is pinned by `msto_regression_corpus` (#1157).
///
/// The two spellings of case A denote one resulting sequence and have two
/// canonical forms, because two different clauses select them:
/// `general.md:34` keeps the allele's three blocks individual (separations of
/// one nucleotide, on a `g.` axis where `general.md:35`'s codon exception
/// cannot apply), and `DNA/delins.md:47` keeps the single delins spanning.
///
/// This assertion was `assert_ne!` under #1160, `assert_eq!` under #1235, and
/// is `assert_ne!` again — so it is pinned here as *both exact strings plus*
/// the inequality, which is strictly stronger than any of the three previous
/// forms. A future change that converges them fails on the exact strings and
/// names which side moved, instead of just flipping a boolean.
#[test]
fn sequence_identical_delins_and_allele_normalize_equal() {
    let from_delins = norm(CASE_A_CORE, CASE_A_DELINS);
    let from_allele = norm(CASE_A_CORE, CASE_A_ALLELE);

    assert_eq!(from_delins, "NC_TEST.1:g.257_263delinsGATTA");
    assert_eq!(
        from_allele,
        "NC_TEST.1:g.[257_258delinsGA;260C>T;262_263del]"
    );
    assert_ne!(
        from_delins, from_allele,
        "the two spellings assert different partitions and are selected by \
         different clauses (delins.md:47 vs general.md:34); if they converge, \
         say which clause was overruled"
    );
}

#[test]
fn decomposed_cis_allele_is_idempotent() {
    let once = norm(CASE_A_CORE, CASE_A_ALLELE);
    assert_eq!(
        once,
        norm(CASE_A_CORE, &once),
        "norm(norm(x)) must equal norm(x)"
    );
}
