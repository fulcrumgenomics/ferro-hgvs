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
// Scope boundary: issue #1157 "case A" is deliberate behavior, not a bug.
// --------------------------------------------------------------------------
//
// A single length-changing `delins` and a decomposed cis allele of the *same*
// edit produce the same resulting sequence but are kept in DISTINCT canonical
// forms: ferro does not merge non-adjacent cis-allele members across unchanged
// reference bases into one spanning `delins` (that would be a coarser, lossy
// description the HGVS spec does not require). "Same resulting sequence" is the
// EquivalenceChecker's job (issue #1158), not the normalizer's.
//
// These are characterization tests: they lock the deliberate non-confluence so
// a future change that accidentally collapses alleles into spanning delins (or
// vice versa) is caught. Both forms are — and must remain — idempotent.
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

#[test]
fn single_delins_is_not_decomposed_into_an_allele() {
    // The length-changing delins stays a single delins; it is not rewritten as
    // the decomposed allele. Together with the test above, this pins that the
    // two sequence-identical encodings remain distinct after normalization.
    assert_eq!(norm(CASE_A_CORE, CASE_A_DELINS), CASE_A_DELINS);
}

#[test]
fn sequence_identical_delins_and_allele_do_not_normalize_equal() {
    // The scope-decision point of #1157 in one assertion: same resulting
    // sequence, different canonical forms — by design.
    assert_ne!(
        norm(CASE_A_CORE, CASE_A_DELINS),
        norm(CASE_A_CORE, CASE_A_ALLELE),
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
