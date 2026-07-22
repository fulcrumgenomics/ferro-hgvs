//! Scope the protein cis-allele coalesce's three remaining declines to the
//! affected members/run instead of vetoing the whole allele — #1130, completing
//! what #1125 did for the to-`Ter` decline.
//!
//! # Spec
//!
//! `protein/substitution.md:23` / `protein/delins.md:18`:
//!
//! > changes involving **two or more consecutive amino acids** are described as
//! > a deletion/insertion variant (delins) […] the description
//! > `p.Arg76_Cys77delinsSerTrp` is correct, the description
//! > `p.[Arg76Ser;Cys77Trp]` is not correct.
//!
//! That obligation attaches to a **run** of consecutive changed residues. An
//! unrelated member elsewhere in the allele — an edit kind this narrow rule does
//! not merge, a second edit at one residue, or a `( )` predicted member — is a
//! reason to leave *that* member alone, not to leave a well-formed run of
//! consecutive changes in the bracket shape the spec calls "not correct".
//!
//! The three declines and what now happens instead:
//!
//! | decline | before | now |
//! |---------|--------|-----|
//! | member is not a single-residue sub/del | whole allele | that member is opaque: emitted verbatim, never merged, and it blocks any member it overlaps |
//! | two members occupy one residue | whole allele | those members are emitted verbatim; every other run merges |
//! | a run mixes predicted `( )` and certain members | whole allele | the run splits at each certainty boundary; each same-certainty sub-run of ≥2 merges |
//!
//! These rewrites are reference-independent (all residues are named in the
//! AST), so an empty `MockProvider` exercises them.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Pin the canonical form of `input` and require it to be a normalization fixed
/// point. Idempotency matters here because a partially-coalesced allele is
/// re-parsed as a bracket whose surviving blocked members must not then change
/// the outcome of a second pass.
fn canonicalizes_to(input: &str, expected: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
    let normalized = Normalizer::new(MockProvider::new())
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("normalize {input:?}: {e}"));
    assert_eq!(normalized.to_string(), expected, "for {input:?}");

    let reparsed =
        parse_hgvs(expected).unwrap_or_else(|e| panic!("parse expected {expected:?}: {e}"));
    let renormalized = Normalizer::new(MockProvider::new())
        .normalize(&reparsed)
        .unwrap_or_else(|e| panic!("re-normalize {expected:?}: {e}"));
    assert_eq!(
        renormalized.to_string(),
        expected,
        "canonical form {expected:?} is not a normalization fixed point",
    );
}

// ---------------------------------------------------------------------------
// Decline 1: a member whose edit kind this rule does not merge
// ---------------------------------------------------------------------------

/// A single-residue `dup` is out of scope for the coalesce, but it must not
/// hold back an unrelated adjacent run.
#[test]
fn an_unmergeable_edit_kind_no_longer_blocks_an_unrelated_run() {
    canonicalizes_to(
        "NP_003997.1:p.[Cys100Ser;Asp101dup;Gly200Ala;Ala201Val]",
        "NP_003997.1:p.[Cys100Ser;Asp101dup;Gly200_Ala201delinsAlaVal]",
    );
}

/// A **multi-residue** member is opaque the same way — it occupies its whole
/// span, and a run clear of that span still merges.
#[test]
fn a_multi_residue_member_no_longer_blocks_an_unrelated_run() {
    canonicalizes_to(
        "NP_003997.1:p.[Cys100_Asp105del;Gly200Ala;Ala201Val]",
        "NP_003997.1:p.[Cys100_Asp105del;Gly200_Ala201delinsAlaVal]",
    );
}

/// An opaque member breaks adjacency at the residues it occupies: a sub at 106
/// overlaps the `100_106del` span, so it is emitted verbatim and cannot pair
/// with the sub at 107.
#[test]
fn an_opaque_member_blocks_a_run_it_overlaps() {
    canonicalizes_to(
        "NP_003997.1:p.[Cys100_Gly106del;Gly106Ala;Ala107Val]",
        "NP_003997.1:p.[Cys100_Gly106del;Gly106Ala;Ala107Val]",
    );
}

/// The opaque member itself is never merged into an adjacent run — a `dup` at
/// 101 sits between subs at 100 and 102, and none of the three coalesce.
#[test]
fn an_opaque_member_is_never_merged_into_a_run() {
    canonicalizes_to(
        "NP_003997.1:p.[Cys100Ser;Asp101dup;Gly102Ala]",
        "NP_003997.1:p.[Cys100Ser;Asp101dup;Gly102Ala]",
    );
}

// ---------------------------------------------------------------------------
// Decline 2: two members at one residue
// ---------------------------------------------------------------------------

/// Two edits at one residue are a contradiction, not a delins — but they must
/// not hold back an unrelated adjacent run.
#[test]
fn a_duplicate_residue_no_longer_blocks_an_unrelated_run() {
    canonicalizes_to(
        "NP_003997.1:p.[Arg76Cys;Arg76Ser;Gly200Ala;Ala201Val]",
        "NP_003997.1:p.[Arg76Cys;Arg76Ser;Gly200_Ala201delinsAlaVal]",
    );
}

/// The contradictory pair itself is still never coalesced (#1116).
#[test]
fn a_duplicate_residue_pair_is_still_never_coalesced() {
    canonicalizes_to(
        "NP_003997.1:p.[Arg76Cys;Arg76Ser]",
        "NP_003997.1:p.[Arg76Cys;Arg76Ser]",
    );
}

/// A residue carrying two edits also breaks a run through it: the subs at 100
/// and 102 are each adjacent to 101, but 101 is contradictory, so nothing
/// merges.
#[test]
fn a_duplicate_residue_breaks_a_run_through_it() {
    canonicalizes_to(
        "NP_003997.1:p.[Cys100Ser;Asp101Gly;Asp101Trp;Arg102Ala]",
        "NP_003997.1:p.[Cys100Ser;Asp101Gly;Asp101Trp;Arg102Ala]",
    );
}

// ---------------------------------------------------------------------------
// Decline 3: a run mixing predicted `( )` and certain members
// ---------------------------------------------------------------------------

/// A mixed predicted/certain run is still not merged (merging would fabricate a
/// certainty the members do not jointly assert), but it no longer vetoes an
/// unrelated run elsewhere in the allele.
#[test]
fn a_mixed_certainty_run_no_longer_blocks_an_unrelated_run() {
    canonicalizes_to(
        "NP_003997.1:p.[(Cys100Ser);Asp101Gly;Gly200Ala;Ala201Val]",
        "NP_003997.1:p.[(Cys100Ser);Asp101Gly;Gly200_Ala201delinsAlaVal]",
    );
}

/// A mixed run splits at the certainty boundary: each same-certainty stretch of
/// ≥2 members merges on its own, so a `[certain, certain, predicted, predicted]`
/// run becomes two delins rather than staying a four-member bracket.
#[test]
fn a_mixed_certainty_run_splits_at_the_boundary() {
    canonicalizes_to(
        "NP_003997.1:p.[Cys100Ser;Asp101Gly;(Gly102Ala);(Ala103Val)]",
        "NP_003997.1:p.[Cys100_Asp101delinsSerGly;(Gly102_Ala103delinsAlaVal)]",
    );
}

/// A two-member mixed run has no same-certainty stretch of ≥2, so neither
/// member merges — the #1116 expectation, unchanged.
#[test]
fn a_two_member_mixed_certainty_run_is_still_not_coalesced() {
    canonicalizes_to(
        "NP_003997.1:p.[(Arg76Ser);Cys77Trp]",
        "NP_003997.1:p.[(Arg76Ser);Cys77Trp]",
    );
}

// ---------------------------------------------------------------------------
// Combined
// ---------------------------------------------------------------------------

/// All three blocked shapes in one allele, plus a to-`Ter` member (#1125), must
/// still leave the one clean run free to coalesce.
#[test]
fn every_decline_at_once_still_lets_a_clean_run_coalesce() {
    canonicalizes_to(
        "NP_003997.1:p.[Arg76Cys;Arg76Ser;(Cys100Ser);Asp101Gly;Trp150dup;Gly200Ala;Ala201Val;Arg300Ter]",
        "NP_003997.1:p.[Arg76Cys;Arg76Ser;(Cys100Ser);Asp101Gly;Trp150dup;Gly200_Ala201delinsAlaVal;Arg300Ter]",
    );
}
