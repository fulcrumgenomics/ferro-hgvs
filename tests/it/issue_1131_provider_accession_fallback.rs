//! A provider that carries protein data for *some* accessions must not block
//! the AST-derivable `delins` minimization for an accession it does **not**
//! carry — #1131.
//!
//! # Scope
//!
//! `try_protein_delins_canonicalize` picks between its reference-backed and
//! reference-free paths on `ReferenceProvider::has_protein_data()`, which is a
//! **provider-wide** capability flag. The reference-backed branch then bails out
//! of the whole function when the fetch fails for the **specific** accession at
//! hand — so the reference-free minimizations of #1119 and #1126 never ran for
//! an accession a real provider does not carry, even though they need no
//! reference at all.
//!
//! A failed fetch means "no reference for this accession", not "cannot
//! canonicalize": fall back to the residues the `delins` location names, exactly
//! as the no-provider path does. The reference is still required for a range
//! spanning ≥ 3 residues (unnamed interior residues) and for the ins→dup
//! refinement, both of which keep declining.
//!
//! Each test builds a provider whose protein set is non-empty (so
//! `has_protein_data()` is true) but does not contain `NP_003997.1`.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// The accession every case below is written against — deliberately absent from
/// the providers built here.
const MISSING: &str = "NP_003997.1";

/// A provider with protein data for an unrelated accession, so
/// `has_protein_data()` is true while `MISSING` cannot be fetched.
fn provider_without_the_accession() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_protein("NP_999999.9", "MAAAAKGSTVE");
    provider
}

/// A protein sequence for `MISSING` positioned so that residue 559 is `V` and
/// residue 560 is `E` (1-based), matching the endpoints the inputs below name.
fn sequence_for_the_accession() -> String {
    let mut seq = String::from("M");
    seq.push_str(&"A".repeat(556)); // residues 2..=557
    seq.push_str("SVE"); // residues 558, 559, 560
    seq.push_str(&"A".repeat(10));
    seq
}

/// Normalize `input` with `provider` and assert the rendered result, then
/// require that result to be a fixed point under the same provider.
fn canonicalizes_with(provider: MockProvider, input: &str, expected: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
    let normalized = Normalizer::new(provider.clone())
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("normalize {input:?}: {e}"));
    assert_eq!(normalized.to_string(), expected, "for {input:?}");

    let reparsed =
        parse_hgvs(expected).unwrap_or_else(|e| panic!("parse expected {expected:?}: {e}"));
    let renormalized = Normalizer::new(provider)
        .normalize(&reparsed)
        .unwrap_or_else(|e| panic!("re-normalize {expected:?}: {e}"));
    assert_eq!(
        renormalized.to_string(),
        expected,
        "canonical form {expected:?} is not a normalization fixed point",
    );
}

/// #1119's single-residue `delins` → substitution reduction must still apply
/// when the provider has protein data for a different accession.
#[test]
fn single_residue_delins_reduces_for_an_accession_the_provider_lacks() {
    canonicalizes_with(
        provider_without_the_accession(),
        &format!("{MISSING}:p.Val559delinsGly"),
        &format!("{MISSING}:p.Val559Gly"),
    );
}

/// #1119's two-residue affix trim → substitution.
#[test]
fn two_residue_delins_trims_for_an_accession_the_provider_lacks() {
    canonicalizes_with(
        provider_without_the_accession(),
        &format!("{MISSING}:p.Val559_Glu560delinsSerGlu"),
        &format!("{MISSING}:p.Val559Ser"),
    );
}

/// #1119's affix trim that leaves a pure deletion.
#[test]
fn two_residue_delins_becomes_deletion_for_an_accession_the_provider_lacks() {
    canonicalizes_with(
        provider_without_the_accession(),
        &format!("{MISSING}:p.Val559_Glu560delinsVal"),
        &format!("{MISSING}:p.Glu560del"),
    );
}

/// #1126's named-flank pure insertion.
#[test]
fn named_flank_insertion_trims_for_an_accession_the_provider_lacks() {
    canonicalizes_with(
        provider_without_the_accession(),
        &format!("{MISSING}:p.Val559_Glu560delinsValSerGlu"),
        &format!("{MISSING}:p.Val559_Glu560insSer"),
    );
}

/// A ≥3-residue range has unnamed interior residues, so the fallback must NOT
/// fabricate them — it still declines, exactly as the no-provider path does.
#[test]
fn three_residue_delins_still_declines_for_an_accession_the_provider_lacks() {
    canonicalizes_with(
        provider_without_the_accession(),
        &format!("{MISSING}:p.Val559_Trp561delinsAlaGlyTrp"),
        &format!("{MISSING}:p.Val559_Trp561delinsAlaGlyTrp"),
    );
}

/// The fallback must not shadow a **successful** reference fetch: with the
/// accession present, the ins→dup refinement the reference enables still wins
/// over the plain `ins` the AST alone would produce.
#[test]
fn a_present_accession_still_takes_the_reference_backed_path() {
    let input = format!("{MISSING}:p.Val559_Glu560delinsValValGlu");

    // Reference-free (the AST-only result): a plain insertion.
    canonicalizes_with(
        provider_without_the_accession(),
        &input,
        &format!("{MISSING}:p.Val559_Glu560insVal"),
    );

    // With the protein sequence, the inserted `Val` is recognized as a tandem
    // duplication of `Val559` — a reduction only the reference can establish.
    let mut with_protein = MockProvider::new();
    with_protein.add_protein(MISSING, sequence_for_the_accession());
    canonicalizes_with(with_protein, &input, &format!("{MISSING}:p.Val559dup"));
}
