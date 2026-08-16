//! #1989: a sequence comparison that could not be reconstructed **for want of
//! reference data** must report [`EquivalenceLevel::Indeterminate`], not a
//! decided negative.
//!
//! The sequence rung (`same_resulting_sequence`) reconstructs each side's edited
//! window and compares the bases. When the two variants share an accession and
//! both convert to SPDI triples, but the reference window itself cannot be
//! obtained from the provider, **nothing was compared**. That is the same
//! situation the [`SequenceVerdict::WindowTooWide`] arm already handles — the
//! window that would settle the pair could not be reconstructed — and it
//! resolves to `Indeterminate`. Its sibling, the *reference-unavailable* decline,
//! used to fall through with `SequenceVerdict::Different` to the ladder's
//! `NotEquivalent` tail, turning "I could not tell" into "they are different".
//!
//! # The measured wrong answer
//!
//! `g.2C>A` and `g.2delCinsA` are the **same edit** — both replace base 2's `C`
//! with an `A`, i.e. identical SPDI triples — so the pair is equivalent. Both
//! convert to SPDI without consulting the provider (their bases are stated), but
//! on a provider that serves no bases for the accession the comparison window
//! cannot be fetched. Before this fix the checker answered `NotEquivalent`: a
//! decided negative asserted for two descriptions that denote the *same*
//! sequence, reached without either sequence ever having been reconstructed.
//!
//! # Scope
//!
//! This targets exactly the `WindowTooWide` sibling: a bounded, same-accession,
//! SPDI-representable window that could not be *reconstructed*. It deliberately
//! does **not** touch the other declines the checker routes to `NotEquivalent`
//! on purpose and pins elsewhere — an RNA fusion
//! (`issue_1578_followup_equivalence_rungs`), a genuine cross-accession pair
//! (a decided negative by construction), or an overlapping cis allele
//! (`issue_1244_equivalence_overlap_panic`). Those are representation- and
//! geometry-level declines, not "the reference could not be read".

use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};
use ferro_hgvs::{parse_hgvs, MockProvider};

use crate::common::synthetic::SyntheticBuilder;

/// An accession the provider serves no bases for, so the comparison window can
/// never be reconstructed.
const ABSENT: &str = "NC_ABSENT.1";

/// The primary reproducer: two spellings of one edit, compared as written, on a
/// provider that cannot serve the reference. The comparison never ran, so the
/// verdict must be `Indeterminate` rather than a decided negative.
#[test]
fn an_equivalent_pair_whose_reference_is_unavailable_is_indeterminate() {
    let checker = EquivalenceChecker::new(MockProvider::new());

    let a = parse_hgvs(&format!("{ABSENT}:g.2C>A")).unwrap();
    let b = parse_hgvs(&format!("{ABSENT}:g.2delCinsA")).unwrap();

    // `compare_denotations` reaches the sequence rung over the pair as written,
    // without normalizing — so the assertion cannot silently depend on whether
    // the normalizer happens to converge the two spellings.
    let result = checker.compare_denotations(&a, &b).unwrap();
    assert_eq!(
        result.level,
        EquivalenceLevel::Indeterminate,
        "a reference-unavailable decline is 'could not tell', not a decided negative; got {:?}",
        result.level,
    );
    assert!(
        !result.level.is_decided(),
        "an unreconstructable comparison is not a decided verdict",
    );
    assert!(
        !result.level.is_equivalent(),
        "and it is not a positive one either",
    );
}

/// The same rung reached through the other entry point. `check` normalizes both
/// sides first; on a provider with no bases the two spellings do not converge,
/// so it reaches the sequence rung and must report the same `Indeterminate`.
///
/// Two genuinely different substitutions are used here, so the point is not that
/// the pair is equivalent but that the checker may not *assert a negative* about
/// a pair it never compared.
#[test]
fn check_reports_indeterminate_when_the_reference_cannot_be_reconstructed() {
    let checker = EquivalenceChecker::new(MockProvider::new());

    let a = parse_hgvs(&format!("{ABSENT}:g.1A>G")).unwrap();
    let b = parse_hgvs(&format!("{ABSENT}:g.2C>G")).unwrap();

    let result = checker.check(&a, &b).unwrap();
    assert_eq!(
        result.level,
        EquivalenceLevel::Indeterminate,
        "the comparison window could not be reconstructed, so the verdict is not decided; got {:?}",
        result.level,
    );
}

/// The discriminating control, and the reason this fix does not over-map. When
/// the reference **is** served, two variants that genuinely produce different
/// sequences are reconstructed and compared, and that decided negative must
/// survive unchanged.
#[test]
fn a_served_genuinely_different_pair_stays_a_decided_negative() {
    let checker = EquivalenceChecker::new(SyntheticBuilder::genomic("ACGTACGTACGT").build());

    let a = parse_hgvs("NC_TEST.1:g.1A>G").unwrap();
    let b = parse_hgvs("NC_TEST.1:g.5A>C").unwrap();

    // Both entry points reach and run the sequence rung, which finds the two
    // reconstructed windows differ.
    assert_eq!(
        checker.check(&a, &b).unwrap().level,
        EquivalenceLevel::NotEquivalent,
        "a reconstructed, genuinely different pair is still `NotEquivalent`",
    );
    assert_eq!(
        checker.compare_denotations(&a, &b).unwrap().level,
        EquivalenceLevel::NotEquivalent,
    );
}

/// The other half of the control: when the reference is served, an equivalent
/// pair still reaches a decided **positive**. The fix moves only the case where
/// the window could not be reconstructed, never a case where it could.
#[test]
fn a_served_equivalent_pair_stays_a_decided_positive() {
    let checker = EquivalenceChecker::new(SyntheticBuilder::genomic("ACGTACGTACGT").build());

    let a = parse_hgvs("NC_TEST.1:g.2C>A").unwrap();
    let b = parse_hgvs("NC_TEST.1:g.2delCinsA").unwrap();

    let result = checker.compare_denotations(&a, &b).unwrap();
    assert!(
        result.level.is_equivalent(),
        "a reconstructed, equivalent pair still reaches a decided positive; got {:?}",
        result.level,
    );
}
