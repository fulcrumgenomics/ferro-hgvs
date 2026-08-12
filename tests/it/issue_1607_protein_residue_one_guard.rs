//! Issue #1607 — a protein `delins` whose affix trim lands on the initiation
//! codon must not be canonicalized into a start-loss substitution.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/1607>.
//!
//! `parse_hgvs` refuses `p.Met1Val` (`validate_no_start_loss_substitution`): a
//! variant that changes the translation initiation codon is described neither
//! as a substitution nor as an extension (`protein/substitution.md:49`,
//! `checklist.md:65`, `protein/extension.md:28`). The input
//! `p.Met1_Ala3delinsValAlaAla` parses because that guard keys on a single
//! certain position and the input is a range — so the illegality is *created*
//! by `try_protein_delins_canonicalize`'s affix trim, and the guard belongs in
//! the producer too.
//!
//! Ferro must not pick among `p.0`, `p.0?`, `p.(Met1?)` and the
//! upstream-initiation insertion form: which is correct is a claim about the
//! *consequence*, and two residue lists do not settle it. The only move
//! available is to decline and leave the input as authored.
//!
//! **Scope.** This fix closes exactly the reparse hole: the producer refuses
//! to emit what the parser refuses to accept. Residuals at residue 1 that are
//! a smaller `delins` or a pure `del` are legal, reparseable descriptions, so
//! they are *not* touched here — whether they too should be declined is a
//! representation question for the operator, and the current answers are
//! pinned below so a later decision is a visible edit rather than a drift.
//!
//! **Two entry arms, both covered.** `try_protein_delins_canonicalize` reaches
//! `residual_del` either from the fetched protein window or — with no protein
//! data for the accession — from the residues named in the input's endpoints
//! (#1119/#1131). The guard is downstream of both, so each arm gets its own
//! positive case and its own negative control.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// `MAAGGG…` — Met at residue 1, Ala at 2 and 3, then a run of Gly.
fn reference_backed_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_protein("NP_003997.1", format!("MAA{}", "G".repeat(20)));
    provider
}

/// A provider carrying **no** protein sequence at all, so
/// `ReferenceProvider::has_protein_data` is false and
/// `try_protein_delins_canonicalize` takes its named-endpoint fallback.
fn reference_free_provider() -> MockProvider {
    MockProvider::new()
}

/// Normalize `input` against `provider` and return the emitted string.
///
/// The provider is a parameter rather than a constant because the guard has two
/// distinct entry arms — see
/// [`a_reference_free_delins_trimming_to_a_residue_one_substitution_is_declined`].
fn normalized_with(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let parsed = parse_hgvs(input).expect("input must parse");
    normalizer
        .normalize(&parsed)
        .expect("normalization must succeed")
        .to_string()
}

/// Normalize `input` against the reference-backed provider.
fn normalized(input: &str) -> String {
    normalized_with(reference_backed_provider(), input)
}

/// Every output this file produces must survive a round trip. That is the
/// property the whole issue is about, so assert it by name rather than only
/// asserting the exact string.
fn assert_reparses(out: &str) {
    assert!(
        parse_hgvs(out).is_ok(),
        "normalization emitted a description parse_hgvs refuses: {out:?}",
    );
}

/// The reported case: the common suffix `AlaAla` trims away, the residual is
/// `Met` → `Val` at residue 1, and the `1 -> 1` residual would take the
/// `sub > delins` branch. Declining leaves the input as authored.
#[test]
fn a_delins_trimming_to_a_residue_one_substitution_is_declined() {
    let out = normalized("NP_003997.1:p.Met1_Ala3delinsValAlaAla");
    assert_eq!(out, "NP_003997.1:p.Met1_Ala3delinsValAlaAla");
    assert_reparses(&out);
    assert_ne!(
        out, "NP_003997.1:p.Met1Val",
        "the start-loss substitution is exactly what must not be emitted",
    );
}

/// The same hole reached without any affix to trim. A single-residue `delins`
/// at the initiator parses — the parser's guard inspects the *edit*, and this
/// one is a `Delins`, not a `Substitution` — and the `1 -> 1` residual then
/// takes the same branch. Pinned separately because it exercises the
/// `lcp == 0 && lcs == 0` fall-through rather than the trim.
#[test]
fn a_bare_residue_one_delins_is_declined() {
    let out = normalized("NP_003997.1:p.Met1delinsVal");
    assert_eq!(out, "NP_003997.1:p.Met1delinsVal");
    assert_reparses(&out);
}

/// A **distinct positive arm**, not a variation of the two cases above.
///
/// `try_protein_delins_canonicalize` sources `residual_del` from two places.
/// When `has_protein_data()` holds it is the fetched protein window; otherwise
/// it is the residues **named in the input's own endpoints**, a fallback that is
/// valid precisely because a ≤ 2-residue range names every residue it deletes
/// (#1119, widened to the per-accession case by #1131). The guard keys on
/// `residual_del[0]`, so it sits downstream of both — and every other case in
/// this file is reference-backed, which would leave the fallback arm carrying
/// the fix untested.
///
/// The arm is reachable and the hole was real. Against a bare provider,
/// `p.Met1_Ala2delinsValAla` trims its common `Ala` suffix to a `Met` → `Val`
/// residual at residue 1, and before this fix that emitted `p.Met1Val` — which
/// `parse_hgvs` refuses. Measured A/B, same input, same bare provider:
///
/// ```text
/// origin/main:  NP_003997.1:p.Met1_Ala2delinsValAla -> NP_003997.1:p.Met1Val   reparses=false
/// this branch:  NP_003997.1:p.Met1_Ala2delinsValAla -> (input unchanged)       reparses=true
/// ```
#[test]
fn a_reference_free_delins_trimming_to_a_residue_one_substitution_is_declined() {
    let out = normalized_with(
        reference_free_provider(),
        "NP_003997.1:p.Met1_Ala2delinsValAla",
    );
    assert_eq!(out, "NP_003997.1:p.Met1_Ala2delinsValAla");
    assert_reparses(&out);
    assert_ne!(
        out, "NP_003997.1:p.Met1Val",
        "the start-loss substitution is exactly what must not be emitted, and \
         the reference-free arm reaches the same branch",
    );
}

/// The reference-free arm's negative control: the same fallback, the same
/// residue-1 span, but the initiator is retained — so the guard must stay out
/// of the way and the trim must still happen.
#[test]
fn a_reference_free_residual_that_keeps_the_initiator_still_canonicalizes() {
    // Met-Ala -> Met-Cys: common prefix `Met`, residual is Ala2 -> Cys.
    let out = normalized_with(
        reference_free_provider(),
        "NP_003997.1:p.Met1_Ala2delinsMetCys",
    );
    assert_eq!(out, "NP_003997.1:p.Ala2Cys");
    assert_reparses(&out);
}

/// The guard keys on the initiator *changing*, not merely on the span
/// touching residue 1. A trim whose residual retains `Met` at residue 1 is
/// not a start loss, so canonicalization proceeds.
#[test]
fn a_residual_that_keeps_the_initiator_still_canonicalizes() {
    // Met-Ala-Ala -> Met-Cys-Ala: common prefix `Met`, common suffix `Ala`,
    // residual is Ala2 -> Cys. Residue 1 is untouched.
    let out = normalized("NP_003997.1:p.Met1_Ala3delinsMetCysAla");
    assert_eq!(out, "NP_003997.1:p.Ala2Cys");
    assert_reparses(&out);
}

/// A trim landing wholly downstream of the initiation codon is unaffected.
#[test]
fn a_delins_away_from_the_initiator_still_canonicalizes() {
    // Ala3-Gly4-Gly5 -> Cys-Gly-Gly: common suffix `GlyGly`, residual
    // Ala3 -> Cys.
    let out = normalized("NP_003997.1:p.Ala3_Gly5delinsCysGlyGly");
    assert_eq!(out, "NP_003997.1:p.Ala3Cys");
    assert_reparses(&out);
}

/// Out of scope, pinned: a residual that is still a `delins` beginning at
/// residue 1. `p.Met1_Ala2delinsValCys` is a legal, reparseable description —
/// the parser's start-loss guard does not apply to a range — so the reparse
/// hole this issue reports is not present here and nothing is changed. Whether
/// the initiator being replaced should nonetheless decline the trim is a
/// representation question left to the operator.
#[test]
fn a_residual_delins_at_residue_one_is_left_as_it_is_today() {
    let out = normalized("NP_003997.1:p.Met1_Gly4delinsValCysAlaGly");
    assert_eq!(out, "NP_003997.1:p.Met1_Ala2delinsValCys");
    assert_reparses(&out);
}

/// Out of scope, pinned: a residual that is a pure deletion covering the
/// initiator. `p.Met1_Ala2del` is legal and reparseable, so — as above — this
/// fix leaves it alone.
///
/// Two reasons the parser does not refuse it, and only the second is the one
/// to cite. `validate_no_start_loss_substitution` requires a **single certain
/// position** (`s == e`), which a range is not; and its `is_residue_swap`
/// matches only `Substitution`, `SubstitutionAlternatives` and an N-terminal
/// `Extension`, so a `Deletion` edit falls to the `_ => false` arm whatever
/// its span. Its doc comment's remark that "a downstream one is a deletion"
/// is **not** that rule: it classifies the *consequence* of a start-codon
/// variant that activates a downstream initiation site, and says nothing
/// about a deletion description spanning residue 1.
#[test]
fn a_residual_deletion_at_residue_one_is_left_as_it_is_today() {
    let out = normalized("NP_003997.1:p.Met1_Ala3delinsAla");
    assert_eq!(out, "NP_003997.1:p.Met1_Ala2del");
    assert_reparses(&out);
}
