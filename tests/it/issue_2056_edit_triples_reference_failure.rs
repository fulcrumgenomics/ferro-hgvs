//! #2056: a reference-data failure **inside `edit_triples`** must report
//! [`EquivalenceLevel::Indeterminate`], not a decided negative.
//!
//! This is the sibling of #1989. #1989 (and its fix, #2053) split a
//! `ReferenceUnavailable` decline out of the sequence rung's **window fetch**
//! (`compare_triples`): the two sides convert to SPDI, but the reference window
//! that would compare them cannot be read. #2056 is the *other* reference route,
//! one step earlier: a member's HGVS→SPDI conversion itself needs to read the
//! reference and cannot.
//!
//! The short forms are the whole point. HGVS routinely omits the deleted or
//! duplicated bases — `g.2del`, `g.2dup`, `c.10_12del` are the ordinary
//! spellings — and every one of those must read the reference to build its SPDI
//! (to learn which bases it deletes / duplicates). On a provider that cannot
//! serve them, `hgvs_to_spdi` returns [`ConversionError::MissingReferenceData`].
//! `edit_triples` used to swallow that with `.ok()?`, turning the reference
//! failure into `None` → `SequenceVerdict::Declined` → the ladder's
//! `NotEquivalent` tail: a decided negative asserted for a pair that was never
//! compared, and which is in fact *equivalent* when the reference is served.
//!
//! # Scope
//!
//! This targets exactly the swallowed HGVS→SPDI reference failure. It does
//! **not** touch the declines the checker routes to `NotEquivalent` on purpose
//! and pins elsewhere — an RNA fusion, a genuine cross-accession pair, an
//! overlapping cis allele (each an *edit-level* decline where the reference was
//! never the obstacle) — nor the window-fetch route #2053 already handles.

use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};
use ferro_hgvs::spdi::{hgvs_to_spdi, ConversionError};
use ferro_hgvs::{parse_hgvs, MockProvider};

use crate::common::synthetic::SyntheticBuilder;

/// An accession the provider serves no bases for, so no short form's SPDI can be
/// built.
const ABSENT: &str = "NC_ABSENT.1";

/// The primary reproducer, and the one that surfaces through **both** entry
/// points. `g.2dup` needs the reference to name its duplicated base and
/// `g.2_3insC` states its own; on an absent provider the two do not converge
/// under normalization, so `check` reaches the sequence rung too. Both must
/// report `Indeterminate`, not the decided negative the swallowed conversion
/// error used to produce.
#[test]
fn a_short_form_needing_an_unavailable_reference_is_indeterminate() {
    let checker = EquivalenceChecker::new(MockProvider::new());

    let a = parse_hgvs(&format!("{ABSENT}:g.2dup")).unwrap();
    let b = parse_hgvs(&format!("{ABSENT}:g.2_3insC")).unwrap();

    for (label, result) in [
        (
            "compare_denotations",
            checker.compare_denotations(&a, &b).unwrap(),
        ),
        ("check", checker.check(&a, &b).unwrap()),
    ] {
        assert_eq!(
            result.level,
            EquivalenceLevel::Indeterminate,
            "{label}: a swallowed reference failure is 'could not tell', not a decided negative; \
             got {:?}",
            result.level,
        );
        assert!(
            !result.level.is_decided(),
            "{label}: an unreconstructable comparison is not a decided verdict",
        );
    }
}

/// The second short form, reached through the direct entry point. `g.2del` needs
/// the reference to name its deleted base; `g.2delC` states it. Via `check` the
/// two happen to converge under normalization (a different, earlier rung), so
/// `compare_denotations` — which does not normalize — is where the swallowed
/// failure is visible.
#[test]
fn a_deletion_short_form_reference_failure_is_indeterminate() {
    let checker = EquivalenceChecker::new(MockProvider::new());

    let a = parse_hgvs(&format!("{ABSENT}:g.2del")).unwrap();
    let b = parse_hgvs(&format!("{ABSENT}:g.2delC")).unwrap();

    let result = checker.compare_denotations(&a, &b).unwrap();
    assert_eq!(
        result.level,
        EquivalenceLevel::Indeterminate,
        "the deletion's base could not be read, so the verdict is not decided; got {:?}",
        result.level,
    );
}

/// The discriminating control, and why this fix does not over-map. When the
/// reference **is** served, the short forms convert and the pair reconstructs to
/// the same sequence — a decided **positive** that must survive unchanged.
#[test]
fn served_short_forms_still_reach_a_decided_positive() {
    let checker = EquivalenceChecker::new(SyntheticBuilder::genomic("ACGTACGTACGT").build());

    for (a, b) in [
        ("NC_TEST.1:g.2dup", "NC_TEST.1:g.2_3insC"),
        ("NC_TEST.1:g.2del", "NC_TEST.1:g.2delC"),
    ] {
        let va = parse_hgvs(a).unwrap();
        let vb = parse_hgvs(b).unwrap();
        let result = checker.compare_denotations(&va, &vb).unwrap();
        assert!(
            result.level.is_equivalent(),
            "{a} vs {b}: a reconstructed, equivalent pair still reaches a decided positive; \
             got {:?}",
            result.level,
        );
    }
}

/// The other half of the control: with a served reference, two variants that
/// genuinely produce different sequences must still be a decided negative. The
/// fix moves only the case where the reference could not be read, never a
/// genuine disagreement.
#[test]
fn served_genuinely_different_short_forms_stay_not_equivalent() {
    let checker = EquivalenceChecker::new(SyntheticBuilder::genomic("ACGTACGTACGT").build());

    // Two deletions at different positions on a served contig — a real
    // disagreement the sequence rung reconstructs.
    //
    // Mind what `g.2` and `g.6` actually address. `SyntheticBuilder::genomic`
    // wraps its core in `PAD_OFFSET` (256) bases of `ACGT…` on each side, so
    // `NC_TEST.1` here is **524** bases and the core `ACGTACGTACGT` starts at
    // `g.257`. Positions 2 and 6 are therefore in the 5' **pad**, not in the
    // core; they happen to read `C` either way only because the pad repeats the
    // same period-4 motif the core does. That padding is also why an
    // out-of-bounds probe on this builder is not out of bounds — `g.20`/`g.21`
    // are served bases, and a stated-ref-base contradiction there declines
    // through a third route (`apply_triples` -> `compare_triples`'s catch-all),
    // which is neither this route nor #2053's window fetch. That third route is
    // #2075 and is still open.
    let a = parse_hgvs("NC_TEST.1:g.2del").unwrap();
    let b = parse_hgvs("NC_TEST.1:g.6del").unwrap();

    assert_eq!(
        checker.compare_denotations(&a, &b).unwrap().level,
        EquivalenceLevel::NotEquivalent,
        "a reconstructed, genuinely different pair is still `NotEquivalent`",
    );
}

/// The transcript `MockProvider::with_test_data` serves completely — bases and
/// all — so nothing on it can decline for want of reference data.
const SERVED_TX: &str = "NM_000088.3";

/// The premise every guard below rests on, asserted rather than assumed: the
/// provider serves this accession's **bases**, not merely its metadata.
///
/// A short form is the instrument. `c.10del` does not state the base it
/// deletes, so its SPDI cannot be built without reading the reference — exactly
/// the dependency this module's reproducers exploit. If this converts, a
/// decline anywhere below is not a reference-availability decline.
#[test]
fn the_probe_transcript_is_genuinely_served() {
    let provider = MockProvider::with_test_data();
    let short_form = parse_hgvs(&format!("{SERVED_TX}:c.10del")).unwrap();
    assert!(
        hgvs_to_spdi(&short_form, &provider).is_ok(),
        "{SERVED_TX} must serve its bases for the intronic guards below to mean anything: a \
         deletion short form has to read the reference to name what it deletes",
    );
}

/// An intronic offset is a limit of **SPDI**, not of the reference, and the
/// conversion error must say so.
///
/// `MissingReferenceData` was the wrong carrier: it is raised here *before* any
/// provider call (`c.`) or on the position's spelling alone (`n.`, `r.`), so it
/// asserts something false — that the reference was consulted and came up
/// short. Consumers classify on that error, and this one classified it as
/// "could not tell", which is what
/// [`intronic_positions_on_a_served_reference_stay_not_equivalent`] measures the
/// consequence of.
#[test]
fn an_intronic_offset_declines_as_a_representation_limit_not_a_reference_failure() {
    let provider = MockProvider::with_test_data();

    for descriptor in [
        format!("{SERVED_TX}:c.10+5A>G"),
        format!("{SERVED_TX}:n.10+5A>G"),
        format!("{SERVED_TX}:r.10+5a>g"),
    ] {
        let variant = parse_hgvs(&descriptor).unwrap();
        let err = hgvs_to_spdi(&variant, &provider)
            .expect_err("SPDI has no offset notation, so an intronic position cannot convert");
        assert!(
            matches!(err, ConversionError::UnrepresentableInSpdi { .. }),
            "{descriptor}: SPDI cannot address an intronic offset however completely the \
             reference is served, so this is `UnrepresentableInSpdi` — a `MissingReferenceData` \
             here is a false statement about the provider, got {err:?}",
        );
    }
}

/// The guard the classification change exists for: on a **fully served**
/// reference, two genuinely different intronic variants stay a decided
/// negative.
///
/// These are ordinary splice-region descriptions — different positions,
/// different bases — and `c.####+/-N` is one of the largest classes in real
/// HGVS, so the blast radius of getting this wrong is wide. Nothing here is
/// unavailable: [`the_probe_transcript_is_genuinely_served`] pins that. The
/// pair declines because SPDI cannot address an intronic offset, which is a
/// representation limit and therefore a decided negative — the same treatment
/// every other unrepresentable shape gets.
///
/// Both entry points are checked because both moved when the classification
/// was wrong.
#[test]
fn intronic_positions_on_a_served_reference_stay_not_equivalent() {
    let checker = EquivalenceChecker::new(MockProvider::with_test_data());

    for (a, b) in [
        // Two different splice-donor variants, one intron apart.
        (
            format!("{SERVED_TX}:c.10+5A>G"),
            format!("{SERVED_TX}:c.20+5T>C"),
        ),
        // An intronic variant against an exonic one at the same anchor.
        (
            format!("{SERVED_TX}:c.10+5A>G"),
            format!("{SERVED_TX}:c.10A>G"),
        ),
        // The same shape on the sibling axes the retype covers.
        (
            format!("{SERVED_TX}:n.10+5A>G"),
            format!("{SERVED_TX}:n.20+5T>C"),
        ),
        (
            format!("{SERVED_TX}:r.10+5a>g"),
            format!("{SERVED_TX}:r.20+5u>c"),
        ),
    ] {
        let va = parse_hgvs(&a).unwrap();
        let vb = parse_hgvs(&b).unwrap();

        for (label, result) in [
            ("compare_denotations", checker.compare_denotations(&va, &vb)),
            ("check", checker.check(&va, &vb)),
        ] {
            let result = result.unwrap();
            assert_eq!(
                result.level,
                EquivalenceLevel::NotEquivalent,
                "{label}: {a} vs {b} are different variants on a served reference; SPDI cannot \
                 address an intronic offset, which is a representation limit and not a reference \
                 failure, so the verdict stays decided. Got {:?} ({:?})",
                result.level,
                result.notes,
            );
            assert!(
                !result.notes.iter().any(|n| n.contains("could not be read")),
                "{label}: {a} vs {b}: the reference was never consulted for an intronic offset, \
                 so no note may claim it could not be read. Got {:?}",
                result.notes,
            );
        }
    }
}
