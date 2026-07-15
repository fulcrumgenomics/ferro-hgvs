//! Issue #1012 item 2 — warn-and-degrade when a genome-requiring
//! normalization step cannot run because the reference provider carries no
//! genomic sequence data.
//!
//! Before this change the behavior was inconsistent: some genome-requiring
//! enhancements were *silently* skipped (a less-normalized result with no
//! warning), while the intronic / boundary-spanning genomic paths returned a
//! hard `Err`. The agreed design surfaces a *uniform* signal on every
//! genome-requiring step: lenient/silent return the best-effort result AND a
//! structured `ReducedCapabilityNoGenome` warning, while strict mode promotes
//! that warning to `FerroError::ReducedReferenceCapability` — so a degraded
//! result is never mistaken for a fully-normalized one in any mode.
//!
//! These tests build a `MockProvider` that knows a transcript with genomic
//! exon coordinates and a chromosome (so intronic/boundary dispatch is
//! reachable) but *no genomic sequence* — hence
//! `ReferenceProvider::has_genomic_data() == false`. They assert that the
//! relevant normalization now returns `Ok` with a best-effort variant and a
//! `REDUCED_CAPABILITY_NO_GENOME` warning, covering both a former-silent path
//! (exon/intron junction 3'-shuffle) and a former-error path (intronic
//! genomic normalization).

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, FerroError, NormalizeConfig, Normalizer, ReferenceProvider};

/// 2-exon plus-strand coding transcript with explicit genomic coordinates on
/// each exon and a chromosome, but with NO genomic sequence registered on the
/// provider. The genomic exon coords make the intronic / boundary-spanning
/// dispatch reachable; the absent genomic sequence makes
/// `has_genomic_data()` false, so the genome-requiring steps must
/// warn-and-degrade.
///
/// Transcript layout (identical shape to issue_392's fixture):
/// ```text
///   tx:   1234567890 ⇢⇢ 1234567890 12
///   seq:  AAAATGCCCC    GGGGTAGAAT AA
///   exon 1 genomic: 100-109   (10 bp)
///   intron 1:       110-129   (20 bp)
///   exon 2 genomic: 130-141   (12 bp)
/// ```
/// `cds_start = 1`, `cds_end = 22`, so `c.<N>` maps 1:1 to tx `<N>`. The
/// `CCCC` homopolymer sits at c.7..=c.10 (the 3' edge of exon 1).
fn provider_no_genomic_sequence() -> MockProvider {
    let mut provider = MockProvider::new();
    let transcript = Transcript::new(
        "NM_INT.1".to_string(),
        Some("INT".to_string()),
        Strand::Plus,
        Some("AAAATGCCCCGGGGTAGAATAA".to_string()),
        Some(1),
        Some(22),
        vec![
            Exon::with_genomic(1, 1, 10, 100, 109),
            Exon::with_genomic(2, 11, 22, 130, 141),
        ],
        Some("chr_int".to_string()),
        Some(100),
        Some(141),
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    // Deliberately do NOT call `add_genomic_sequence` — has_genomic_data() == false.
    provider
}

#[test]
fn provider_reports_no_genomic_data() {
    // Guards the premise of every test below: the fixture genuinely lacks
    // genomic sequence data even though its exons carry genomic coordinates.
    assert!(
        !provider_no_genomic_sequence().has_genomic_data(),
        "fixture must not carry genomic sequence data",
    );
}

/// Former-ERROR path: an intronic variant used to return
/// `Err(FerroError::IntronicVariant)` when the provider had no genomic data.
/// It must now succeed with the best-effort (unchanged) variant plus the
/// reduced-capability warning. The genomic-context (`NG_(NM_)`) form is used
/// so the input is spec-valid and does not also raise the EINTRONIC
/// bare-transcript warning — isolating the capability signal.
#[test]
fn intronic_variant_warns_and_degrades_without_genome() {
    let normalizer =
        Normalizer::with_config(provider_no_genomic_sequence(), NormalizeConfig::lenient());
    let variant = parse_hgvs("NG_012337.1(NM_INT.1):c.10+5del").expect("parse");

    let diag = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("intronic normalization without a genome must not error (warn-and-degrade)");

    // Best-effort result: the input is echoed UNCHANGED (cannot shuffle without
    // a genome). Compare against the parsed input's own Display — an independent
    // oracle — rather than a substring/literal.
    assert_eq!(
        diag.result.to_string(),
        variant.to_string(),
        "best-effort output must preserve the input exactly",
    );

    // The reduced-capability warning is surfaced.
    let reduced = diag
        .warnings
        .iter()
        .find(|w| w.code() == "REDUCED_CAPABILITY_NO_GENOME");
    assert!(
        reduced.is_some(),
        "expected REDUCED_CAPABILITY_NO_GENOME warning, got {:?}",
        diag.warnings,
    );
    // The message is human-readable and names the missing capability.
    let message = reduced.unwrap().message();
    assert!(
        message.contains("genomic sequence data") && message.contains("intronic normalization"),
        "warning message should describe the missing capability, got {message:?}",
    );

    // This spec-valid genomic-context form must NOT raise EINTRONIC.
    assert!(
        !diag
            .warnings
            .iter()
            .any(|w| w.code() == "INTRONIC_ON_BARE_TRANSCRIPT"),
        "genomic-context intronic form must not raise EINTRONIC, got {:?}",
        diag.warnings,
    );
}

/// Former-SILENT path: a purely-exonic deletion that comes to rest at an
/// exon's 3' boundary triggers the #670 exon/intron junction 3'-shuffle, which
/// requires a genome. Previously that enhancement was skipped with NO warning
/// when the provider lacked genomic data, so the caller could not tell the
/// exon-confined result apart from a genuinely-final one. It must now emit the
/// reduced-capability warning while still returning the exon-confined result.
#[test]
fn exon_boundary_del_warns_and_degrades_without_genome() {
    let normalizer =
        Normalizer::with_config(provider_no_genomic_sequence(), NormalizeConfig::lenient());
    // `c.10del` deletes the 3'-most base of the `CCCC` homopolymer (c.7..=c.10),
    // which is also the last base of exon 1 — the exon 3' boundary landing that
    // arms the junction-crossing 3'-shuffle.
    let variant = parse_hgvs("NM_INT.1:c.10del").expect("parse");

    let diag = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("exonic normalization must succeed");

    // Best-effort result is the exon-confined deletion — the input echoed
    // unchanged (compare against the parsed input's own Display, not a literal).
    assert_eq!(
        diag.result.to_string(),
        variant.to_string(),
        "best-effort output must be the exon-confined deletion, unchanged",
    );

    assert!(
        diag.warnings
            .iter()
            .any(|w| w.code() == "REDUCED_CAPABILITY_NO_GENOME"),
        "expected REDUCED_CAPABILITY_NO_GENOME warning on the exon-boundary landing, got {:?}",
        diag.warnings,
    );
}

/// A plain exonic deletion that does NOT land at an exon boundary must not be
/// flagged — the reduced-capability warning is scoped to genome-requiring
/// steps that were actually reached, not every transcript-only normalization.
#[test]
fn interior_exonic_del_is_not_flagged() {
    let normalizer =
        Normalizer::with_config(provider_no_genomic_sequence(), NormalizeConfig::lenient());
    // `c.5del` (the `T` at position 5) sits in the interior of exon 1, far from
    // the exon 3' boundary — no junction shuffle is armed.
    let variant = parse_hgvs("NM_INT.1:c.5del").expect("parse");

    let diag = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("interior exonic normalization must succeed");

    assert!(
        !diag
            .warnings
            .iter()
            .any(|w| w.code() == "REDUCED_CAPABILITY_NO_GENOME"),
        "interior exonic deletion must not raise REDUCED_CAPABILITY_NO_GENOME, got {:?}",
        diag.warnings,
    );
}

/// Strict mode REJECTS a reduced-capability degradation rather than returning a
/// best-effort result. The signal is uniform across every genome-requiring
/// step: `ReducedCapabilityNoGenome` is promoted to
/// `FerroError::ReducedReferenceCapability` in strict mode (like the sibling
/// `CanonicalSplitSkipped`→`VariantExceedsReference` promotion), so a strict
/// caller — which typically inspects only the `Result`, not the warnings vec —
/// never silently accepts a knowingly-degraded variant.
///
/// This covers BOTH categories the warning spans: the former-error intronic
/// path (rejected in strict pre-#1012 too — behavior preserved) and the
/// former-silent exon-junction path (which strict used to accept silently —
/// now uniformly rejected, so strict never accepts a degraded result).
#[test]
fn strict_mode_rejects_reduced_capability() {
    let normalizer =
        Normalizer::with_config(provider_no_genomic_sequence(), NormalizeConfig::strict());

    // Former-ERROR path: an intronic variant with no genome.
    let intronic = parse_hgvs("NG_012337.1(NM_INT.1):c.10+5del").expect("parse");
    assert!(
        matches!(
            normalizer.normalize(&intronic),
            Err(FerroError::ReducedReferenceCapability { .. })
        ),
        "strict mode must reject the intronic reduced-capability path",
    );

    // Former-SILENT path: an exon-3'-boundary deletion that arms the
    // junction-crossing shuffle. Strict used to accept this silently; it is now
    // rejected too, so strict never returns a degraded result.
    let boundary = parse_hgvs("NM_INT.1:c.10del").expect("parse");
    assert!(
        matches!(
            normalizer.normalize(&boundary),
            Err(FerroError::ReducedReferenceCapability { .. })
        ),
        "strict mode must reject the exon-boundary reduced-capability path",
    );
}
