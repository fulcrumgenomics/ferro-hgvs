//! Issue #429 — provider-aware W3016 length-mismatch detection for
//! mixed-shape intronic-endpoint ranges.
//!
//! PR #406 (#390 item 5) added the same-anchor + same-sign intronic
//! detection but explicitly skipped mixed-shape inputs because
//! computing the actual span requires intron-length data only
//! available through a `ReferenceProvider`. This file pins the
//! provider-aware extension:
//!
//!   - `c.100+5_200-3del<refseq>`: both intronic, different anchors.
//!   - `c.100+5_100-3del<refseq>`: same anchor, opposite signs.
//!   - `c.100_200+5del<refseq>`: one exonic + one intronic.
//!
//! The no-provider variant (`detect_length_mismatch`) continues to
//! silently skip these shapes — confirmed by the no-provider parity
//! tests below.
//!
//! # Test fixture shape
//!
//! Two-exon transcript on the plus strand. Genomic layout:
//!
//!   ```text
//!   genomic:  1000 ........ 1099  | 2000 ........ 2099
//!             [-- exon 1 ---]      [-- exon 2 ---]
//!   tx:       1 ........ 100       101 ........ 200
//!   ```
//!
//! Intron between genomic [1100, 1999] (length 900). CDS spans the
//! full transcript (cds_start=1, cds_end=200).

use ferro_hgvs::error_handling::corrections::{
    detect_length_mismatch, detect_length_mismatch_with_provider,
};
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{
    Exon, GenomeBuild, ManeStatus, Strand as TxStrand, Transcript,
};

const ACC: &str = "NM_TEST.1";

/// Build a two-exon plus-strand transcript with a 900-bp intron.
/// CDS spans the full transcript. Used by all provider-aware tests.
fn two_exon_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    // Exon 1: tx 1..100 ↔ genomic 1000..1099.
    let e1 = Exon::with_genomic(1, 1, 100, 1000, 1099);
    // Exon 2: tx 101..200 ↔ genomic 2000..2099.
    let e2 = Exon::with_genomic(2, 101, 200, 2000, 2099);
    provider.add_transcript(Transcript::new(
        ACC.to_string(),
        Some("TEST".to_string()),
        TxStrand::Plus,
        None::<String>,
        Some(1),
        Some(200),
        vec![e1, e2],
        Some("chr1".to_string()),
        Some(1000),
        Some(2099),
        GenomeBuild::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    provider
}

// =============================================================================
// No-provider parity: mixed-shape inputs continue to silently skip
// =============================================================================

#[test]
fn no_provider_mixed_anchor_different_offsets_silently_skips() {
    // `c.100+5_200-3delAAAA` — without a provider, the detector skips
    // (per #390 item 5's "deliberately skipped" contract).
    let input = format!("{ACC}:c.100+5_200-3delAAAA");
    let hits = detect_length_mismatch(&input);
    assert!(
        hits.is_empty(),
        "no-provider mixed-shape input must silently skip; got {hits:?}",
    );
}

#[test]
fn no_provider_mixed_axis_exonic_intronic_silently_skips() {
    let input = format!("{ACC}:c.100_200+5delAAAA");
    let hits = detect_length_mismatch(&input);
    assert!(
        hits.is_empty(),
        "no-provider exonic-intronic input must silently skip; got {hits:?}",
    );
}

// =============================================================================
// Provider-aware detection: each mixed-shape pattern resolves cleanly
// =============================================================================

/// `c.100+5_101-3delAAAA`: both endpoints in the SAME intron
/// (the intron between exon 1 and exon 2), referenced from different
/// boundaries.
///
/// - `100+5`: 5 bp into intron from the 5' boundary (after exon 1).
///   Genomic = intron_5'_start + 5 - 1 = 1100 + 5 - 1 = 1104.
/// - `101-3`: 3 bp before the 3' boundary (before exon 2). Genomic =
///   intron_3'_end - 3 + 1 = 1999 - 3 + 1 = 1997.
///
/// Span = 1997 - 1104 + 1 = 894 bp. Declared ref seq `AAAA` = 4 bp →
/// mismatch (894 != 4).
#[test]
fn provider_mixed_anchor_different_intronic_offsets_detects_mismatch() {
    let provider = two_exon_provider();
    let input = format!("{ACC}:c.100+5_101-3delAAAA");
    let hits = detect_length_mismatch_with_provider(&input, &provider);
    assert_eq!(
        hits.len(),
        1,
        "provider-aware mixed-shape mismatch must surface exactly one hit; got {hits:?}",
    );
    // The hit's `original` field captures the span from just after the
    // `c.` marker through the end of the ref seq (matches the existing
    // detector's span convention).
    assert!(
        hits[0].original.contains("100+5_101-3delAAAA"),
        "hit should cover the mismatched range; got original={:?}",
        hits[0].original,
    );
}

/// Same input but the ref seq matches the resolved span: no hit.
/// The intron is 900 bp; the declared range `100+5_101-3` is
/// genomic 1104..1997 = 894 bp. A 894-bp explicit `delN...` ref would
/// match (impractical to spell out — use a much smaller in-range
/// example to confirm the no-mismatch path doesn't fire).
#[test]
fn provider_intron_internal_same_intron_consistent_no_hit() {
    let provider = two_exon_provider();
    // `c.100+1_100+3delAAA`: both intronic, same anchor (100), same
    // sign (+). Span = 3 bp. Ref = 3 bp. Match — no hit. This shape
    // is already handled by the same-anchor + same-sign branch, so
    // the no-provider variant also returns empty. Asserted as a
    // baseline that the provider path doesn't introduce false
    // positives on the already-handled shape.
    let input = format!("{ACC}:c.100+1_100+3delAAA");
    let hits = detect_length_mismatch_with_provider(&input, &provider);
    assert!(
        hits.is_empty(),
        "consistent length on already-handled shape must not produce a hit; got {hits:?}",
    );
}

/// One exonic + one intronic: `c.100_100+5delAAAA`.
/// Start at exonic c.100 → genomic 1099 (end of exon 1).
/// End at intronic 100+5 → genomic 1100 + 5 - 1 = 1104.
/// Span = 1104 - 1099 + 1 = 6 bp. Ref AAAA = 4 bp → mismatch.
#[test]
fn provider_exonic_to_intronic_detects_mismatch() {
    let provider = two_exon_provider();
    let input = format!("{ACC}:c.100_100+5delAAAA");
    let hits = detect_length_mismatch_with_provider(&input, &provider);
    assert_eq!(
        hits.len(),
        1,
        "exonic-to-intronic mixed-shape must surface a hit when ref length disagrees; got {hits:?}",
    );
}

/// Exonic-to-intronic with a matching ref seq: no hit.
#[test]
fn provider_exonic_to_intronic_consistent_no_hit() {
    let provider = two_exon_provider();
    // `c.100_100+5delAAAAAA`: 6 bp ref = 6 bp span (genomic 1099..1104).
    let input = format!("{ACC}:c.100_100+5delAAAAAA");
    let hits = detect_length_mismatch_with_provider(&input, &provider);
    assert!(
        hits.is_empty(),
        "matching ref length on exonic-to-intronic mixed-shape must not produce a hit; got {hits:?}",
    );
}

/// Gene-in-parentheses transcript selector (`NG_...(NM_...):c.`): the
/// backward accession scan stops at `(` but must not leave a trailing
/// `)` on the captured id, or `get_transcript("NM_TEST.1)")` fails and
/// the mixed-shape range is silently skipped. With the `)` trimmed, the
/// provider lookup resolves and the mismatch is detected — same hit as
/// the bare `NM_TEST.1:c.100+5_101-3delAAAA` form above.
///
/// Regression for the PR #436 CR finding on accession normalization.
#[test]
fn provider_gene_in_parens_selector_resolves_transcript() {
    let provider = two_exon_provider();
    let input = format!("NG_016465.4({ACC}):c.100+5_101-3delAAAA");
    let hits = detect_length_mismatch_with_provider(&input, &provider);
    assert_eq!(
        hits.len(),
        1,
        "gene-in-parens selector must resolve the inner transcript and detect the mismatch; got {hits:?}",
    );
}

// =============================================================================
// Existing shapes continue to work under the provider-aware entry point
// =============================================================================

/// All-exonic `c.100_105delAAAA`: declared 4 bp, span is 6 bp → mismatch.
/// The all-exonic path is handled by both no-provider and provider-aware
/// variants (same branch). Pin that the provider entry point doesn't
/// regress the existing shapes.
#[test]
fn provider_all_exonic_same_axis_detects_mismatch() {
    let provider = two_exon_provider();
    let input = format!("{ACC}:c.100_105delAAAA");
    let hits = detect_length_mismatch_with_provider(&input, &provider);
    assert_eq!(
        hits.len(),
        1,
        "all-exonic mismatch must still surface; got {hits:?}"
    );
}

/// Same-anchor same-sign intronic `c.100+5_100+10delAAAA`: declared
/// 4 bp, span is 6 bp → mismatch. Handled by the (b) branch in the
/// no-provider path; pin that provider-aware path also handles it.
#[test]
fn provider_same_anchor_same_sign_intronic_still_detects_mismatch() {
    let provider = two_exon_provider();
    let input = format!("{ACC}:c.100+5_100+10delAAAA");
    let hits = detect_length_mismatch_with_provider(&input, &provider);
    assert_eq!(
        hits.len(),
        1,
        "same-anchor same-sign intronic mismatch must still surface; got {hits:?}",
    );
}

// =============================================================================
// Graceful fallback: unknown accession / no transcript → silent skip
// =============================================================================

#[test]
fn provider_unknown_accession_silently_skips_mixed_shape() {
    let provider = two_exon_provider();
    // Accession not in provider — `get_transcript` returns Err →
    // `resolve_cds_endpoint_to_genomic` returns None → mixed-shape
    // path returns None → no hit. Crucially: no panic, no false
    // positive.
    let input = "NM_UNKNOWN.1:c.100+5_200-3delAAAA";
    let hits = detect_length_mismatch_with_provider(input, &provider);
    assert!(
        hits.is_empty(),
        "unknown accession must silently skip the mixed-shape path; got {hits:?}",
    );
}

#[test]
fn provider_g_axis_mixed_offsets_skipped_per_scope() {
    // The provider-aware path is scoped to the `c.` axis. A g.-axis
    // input with intronic-looking offsets (which would be malformed
    // HGVS in practice but exercises the scope guard) silently skips.
    let provider = two_exon_provider();
    let input = "NC_000001.11:g.100+5_200-3delAAAA";
    let hits = detect_length_mismatch_with_provider(input, &provider);
    assert!(
        hits.is_empty(),
        "non-`c.` axis must skip the provider path; got {hits:?}",
    );
}

/// Compound-allele inner member axis tracking: a `g.[X;Y]` compound's
/// inner member must NOT be routed through the c.-axis provider path
/// just because the outer accession is non-empty. Pre-fix, the inner
/// dispatch always used a hardcoded `b'c'` axis byte, creating a
/// latent false-positive risk on `g.` compounds. This test pins the
/// fix that propagates the outer axis to inner members.
#[test]
fn provider_g_axis_compound_inner_member_skipped_per_scope() {
    let provider = two_exon_provider();
    // Outer axis is `g.`; the inner member `100+5_200-3delAAAA` would
    // otherwise be a mixed-shape intronic range that the provider path
    // could erroneously process. The axis guard (now correctly using
    // the outer `g.` axis) must keep it out of the provider path.
    let input = "NC_000001.11:g.[100A>G;100+5_200-3delAAAA]";
    let hits = detect_length_mismatch_with_provider(input, &provider);
    // The inner mixed-shape intronic range must not be processed via
    // the provider path. The hit list may be empty or contain only the
    // outer-axis exonic hit (none here). Crucially, NO mixed-shape hit
    // appears that would imply c.-axis resolution against `NC_000001.11`.
    for hit in &hits {
        assert!(
            !hit.original.contains("100+5_200-3delAAAA"),
            "g.-compound inner member must not be c.-axis-resolved; \
             got hit on the mixed-shape inner: {:?}",
            hit.original,
        );
    }
}

// =============================================================================
// Minus-strand resolution
// =============================================================================

/// Build a minus-strand two-exon transcript. Exon ordering follows
/// transcript direction, but genomic coordinates are reversed:
///
///   ```text
///   genomic:  1000 ........ 1099  | 2000 ........ 2099
///                                     [-- exon 1 ---]      <- minus
///             [-- exon 2 ---]                              <- minus
///   tx:                                1 ........ 100
///                                                  101 ........ 200
///   ```
///
/// Equivalently: on minus strand, transcript 5' end maps to the
/// HIGHER genomic position. Exon 1 (tx 1..100) maps to genomic
/// 2000..2099 reversed; exon 2 (tx 101..200) maps to genomic
/// 1000..1099 reversed.
fn minus_strand_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let e1 = Exon::with_genomic(1, 1, 100, 2000, 2099);
    let e2 = Exon::with_genomic(2, 101, 200, 1000, 1099);
    provider.add_transcript(Transcript::new(
        "NM_MINUS.1".to_string(),
        Some("MTEST".to_string()),
        TxStrand::Minus,
        None::<String>,
        Some(1),
        Some(200),
        vec![e1, e2],
        Some("chr1".to_string()),
        Some(1000),
        Some(2099),
        GenomeBuild::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    provider
}

/// Minus-strand exonic-to-intronic span: detection must not break on
/// the genomic-order-reversed case. The `|g_end - g_start| + 1`
/// absolute-value formula handles minus strand transparently because
/// the W3016 length check only cares about span magnitude.
#[test]
fn provider_minus_strand_exonic_to_intronic_detects_mismatch() {
    let provider = minus_strand_provider();
    let input = "NM_MINUS.1:c.100_100+5delAAAA";
    let hits = detect_length_mismatch_with_provider(input, &provider);
    // Whether the hit fires depends on the resolved span. The point
    // here is that the call does not panic and produces a consistent
    // shape (either a hit or an empty list, but never an error). The
    // resolved span should be the same magnitude as the plus-strand
    // case (6 bp) — declared 4 bp → mismatch.
    assert_eq!(
        hits.len(),
        1,
        "minus-strand exonic-to-intronic mismatch must surface a hit; got {hits:?}",
    );
}
