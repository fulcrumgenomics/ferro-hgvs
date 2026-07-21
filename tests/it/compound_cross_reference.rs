//! Cross-reference / cross-coordinate-system compound alleles (#81 H2).
//!
//! HGVS allows compound alleles whose sub-variants reference different
//! sequences and/or coordinate systems, e.g.
//!
//!     [NM_000088.3:c.1A>G;NC_000017.11:g.2A>G]
//!     [NM_000088.3:c.1A>G];[NC_000017.11:g.2A>G]
//!
//! The spec is silent on a single canonical shape for these — there is no
//! shared coordinate system to pivot through, no coordinate range syntax
//! that crosses references, and the merge pass introduced in PR #80
//! (issue #72) explicitly excludes them as a barrier (see
//! `tests/merge_consecutive_edits_tests.rs::test_no_merge_different_accessions`
//! and `::test_no_merge_different_variant_types`).
//!
//! These tests pin the canonicalization that ferro-hgvs has chosen:
//!
//! 1. Parse round-trip — the input string is preserved verbatim, including
//!    the relative order of sub-variants.
//! 2. Phase is preserved (cis `[a;b]`, trans `[a];[b]`, mosaic `a/b`,
//!    chimeric `a//b`, unknown `[a(;)b]`).
//! 3. Display falls back to the *expanded* form (each sub-variant prints
//!    its own `ACC:type.` prefix) whenever the sub-variants do not share
//!    both accession and coordinate type. The compact form
//!    (`ACC:type.[edit;edit]`) is a strict optimization for the common case.
//! 4. Each sub-variant normalizes independently against its own reference;
//!    a c. dup with a transcript provider 3'-shifts even when its
//!    cis-companion is a g. variant on a different (un-resolvable) accession.
//! 5. No merging happens across references or coordinate systems — the PR
//!    #80 contract holds.

use ferro_hgvs::hgvs::variant::{AllelePhase, AlleleVariant};
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, Normalizer};

// =====================================================================
// Helpers
// =====================================================================

/// Parse `input`, normalize through a `MockProvider` (no reference data),
/// and return the rendered string. Mirrors the helper used in
/// `merge_consecutive_edits_tests.rs`.
fn normalize_to_string(input: &str) -> String {
    let normalizer = Normalizer::new(MockProvider::new());
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

/// Parse `input` and assert that it round-trips through `Display` unchanged.
/// Returns the parsed variant for further inspection by the caller.
fn assert_parse_roundtrip(input: &str) -> HgvsVariant {
    let parsed = parse_hgvs(input).expect("parse failed");
    assert_eq!(
        format!("{}", parsed),
        input,
        "round-trip mismatch for {}",
        input,
    );
    parsed
}

/// Assert the parsed variant is an `Allele` with the expected phase and
/// sub-variant count, and return the inner `AlleleVariant`.
fn expect_allele(
    variant: &HgvsVariant,
    expected_phase: AllelePhase,
    expected_len: usize,
) -> &AlleleVariant {
    let allele = match variant {
        HgvsVariant::Allele(a) => a,
        other => panic!("expected Allele, got {:?}", other),
    };
    assert_eq!(
        allele.phase, expected_phase,
        "phase mismatch: got {:?}",
        allele.phase
    );
    assert_eq!(
        allele.variants.len(),
        expected_len,
        "sub-variant count mismatch: got {}",
        allele.variants.len(),
    );
    allele
}

/// Build a 60-base transcript provider keyed by `NM_TEST.1` so that c./n./r.
/// sub-variants can normalize against a real sequence. Genomic accessions
/// (e.g. `NC_000001.11`) remain unresolved — the matching test cases assert
/// that those still pass through unchanged while the c. partner normalizes.
fn provider_with_simple_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence: String =
        "ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT".to_string();
    let len = sequence.len() as u64;
    let exons = vec![Exon::new(1, 1, len)];
    let transcript = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(1),
        Some(60),
        exons,
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

// =====================================================================
// 1. Cis cross-reference round-trip and order preservation.
// =====================================================================

#[test]
fn test_cis_cross_reference_round_trip() {
    // The canonical example from the issue. Cross-reference cis alleles
    // emit the expanded form because the sub-variants do not share an
    // accession; the input string must round-trip verbatim.
    let input = "[NM_000088.3:c.1A>G;NC_000017.11:g.2A>G]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert!(matches!(allele.variants[0], HgvsVariant::Cds(_)));
    assert!(matches!(allele.variants[1], HgvsVariant::Genome(_)));
}

#[test]
fn test_cis_cross_reference_preserves_input_order() {
    // Pin the canonicalization decision: ferro-hgvs preserves the input
    // ordering of sub-variants. It does NOT sort them by coordinate
    // system, accession, or position. The HGVS spec has no canonical
    // ordering for cross-reference compound alleles, so reordering would
    // be an unfounded invention.
    let forward = "[NM_000088.3:c.1A>G;NC_000017.11:g.2A>G]";
    let reverse = "[NC_000017.11:g.2A>G;NM_000088.3:c.1A>G]";
    let forward_parsed = assert_parse_roundtrip(forward);
    let reverse_parsed = assert_parse_roundtrip(reverse);
    // The two parses must NOT be equal — order is significant.
    assert_ne!(
        forward_parsed, reverse_parsed,
        "reversed cross-ref alleles must not compare equal — order is preserved",
    );
}

// =====================================================================
// 2. Trans cross-reference round-trip.
// =====================================================================

#[test]
fn test_trans_cross_reference_round_trip() {
    // Trans phase: each variant lives on its own bracketed group.
    let input = "[NM_000088.3:c.1A>G];[NC_000017.11:g.2A>G]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Trans, 2);
    assert!(matches!(allele.variants[0], HgvsVariant::Cds(_)));
    assert!(matches!(allele.variants[1], HgvsVariant::Genome(_)));
}

// =====================================================================
// 3. Mixed coordinate systems within the same compound allele.
//
// The spec lists c., g., n., r., p., m., o. as distinct coordinate
// systems. Cross-coord compounds are valid HGVS but undefined under
// canonicalization; ferro-hgvs preserves them as-is.
// =====================================================================

#[test]
fn test_cis_mixed_c_and_g() {
    // c. + g. — already covered by the round-trip above; this case
    // pins the explicit "different coordinate types" assertion.
    let input = "[NM_000088.3:c.1A>G;NC_000017.11:g.2A>G]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert_eq!(allele.variants[0].variant_type(), "c");
    assert_eq!(allele.variants[1].variant_type(), "g");
}

#[test]
fn test_cis_mixed_c_and_r() {
    // c. + r. on different accessions (NM_/NR_).
    let input = "[NM_000088.3:c.1A>G;NR_000088.3:r.2a>g]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert_eq!(allele.variants[0].variant_type(), "c");
    assert_eq!(allele.variants[1].variant_type(), "r");
}

#[test]
fn test_cis_mixed_c_and_p() {
    // c. + p. on the matched mRNA / protein pair.
    let input = "[NM_000088.3:c.1A>G;NP_000079.2:p.Met2Val]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert_eq!(allele.variants[0].variant_type(), "c");
    assert_eq!(allele.variants[1].variant_type(), "p");
}

#[test]
fn test_cis_mixed_c_and_n_same_accession() {
    // Same accession, different coord types — still NOT compactable
    // because the type prefix differs, and not mergeable for the same
    // reason.
    let input = "[NM_000088.3:c.1A>G;NM_000088.3:n.2A>G]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert_eq!(allele.variants[0].variant_type(), "c");
    assert_eq!(allele.variants[1].variant_type(), "n");
}

#[test]
fn test_cis_three_different_coord_systems() {
    // Three sub-variants, three different coord systems (c., g., p.).
    // Pins that the expanded form scales beyond two variants.
    let input = "[NM_000088.3:c.1A>G;NC_000017.11:g.2A>G;NP_000079.2:p.Met2Val]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 3);
    assert_eq!(allele.variants[0].variant_type(), "c");
    assert_eq!(allele.variants[1].variant_type(), "g");
    assert_eq!(allele.variants[2].variant_type(), "p");
}

// =====================================================================
// 4. Other phases across coordinate systems.
// =====================================================================

#[test]
fn test_mosaic_cross_reference_round_trip() {
    let input = "NM_000088.3:c.1A>G/NC_000017.11:g.2A>G";
    let parsed = assert_parse_roundtrip(input);
    expect_allele(&parsed, AllelePhase::Mosaic, 2);
}

#[test]
fn test_chimeric_cross_reference_round_trip() {
    let input = "NM_000088.3:c.1A>G//NC_000017.11:g.2A>G";
    let parsed = assert_parse_roundtrip(input);
    expect_allele(&parsed, AllelePhase::Chimeric, 2);
}

#[test]
fn test_unknown_phase_cross_reference_round_trip() {
    let input = "[NM_000088.3:c.1A>G(;)NC_000017.11:g.2A>G]";
    let parsed = assert_parse_roundtrip(input);
    expect_allele(&parsed, AllelePhase::Unknown, 2);
}

// =====================================================================
// 5. Per-variant independent normalization against each variant's own
//    reference.
//
// The spec requires each sub-variant in a compound allele to be
// validated/normalized against the reference its accession names. With
// `provider_with_simple_transcript`, only `NM_TEST.1` is registered, so
// a c. variant on that transcript can normalize while its g. companion
// on `NC_000001.11` (no genome data) passes through unchanged.
// =====================================================================

#[test]
fn test_each_variant_normalizes_against_own_ref() {
    // c.10dupA at the start of a `CCCCC` run normalizes to c.14dup
    // (drop the mismatched stated-ref base AND 3'-shift to the
    // rightmost copy of the homopolymer, both in a single pass —
    // see #219). The g. companion has no resolvable ref data and must
    // round-trip unchanged. Pins both:
    //   * per-variant normalization is independent
    //   * the un-resolvable companion does not block normalization of
    //     its sibling
    let normalizer = Normalizer::new(provider_with_simple_transcript());
    let parsed = parse_hgvs("[NM_TEST.1:c.10dupA;NC_000001.11:g.100A>G]").expect("parse");
    let normalized = normalizer.normalize(&parsed).expect("normalize");
    assert_eq!(
        format!("{}", normalized),
        "[NM_TEST.1:c.14dup;NC_000001.11:g.100A>G]",
    );
}

#[test]
fn test_unresolved_companion_does_not_break_resolved_variant_normalize() {
    // Mirror the previous case in trans phase. Each bracket group
    // normalizes independently; c.10dupA 3'-shifts through the CCCCC
    // homopolymer to c.14dup (single pass per #219), the g. variant
    // rides through.
    let normalizer = Normalizer::new(provider_with_simple_transcript());
    let parsed = parse_hgvs("[NM_TEST.1:c.10dupA];[NC_000001.11:g.100A>G]").expect("parse");
    let normalized = normalizer.normalize(&parsed).expect("normalize");
    assert_eq!(
        format!("{}", normalized),
        "[NM_TEST.1:c.14dup];[NC_000001.11:g.100A>G]",
    );
}

// =====================================================================
// 6. No-merge contract (PR #80).
//
// PR #80 introduced consecutive-edit merging in cis alleles, but only
// for sub-variants sharing both accession and coordinate type
// (`HgvsVariant::all_share_accession_and_type`). These tests pin that
// the merge pass is correctly off-limits for cross-reference and
// cross-coord alleles, even when the positions are numerically
// adjacent.
// =====================================================================

#[test]
fn test_no_merge_across_references() {
    // Two genomic SNVs on adjacent positions but different accessions —
    // covered by the existing merge tests, but re-pinned here in the
    // cross-reference context.
    let result = normalize_to_string("[NC_000001.11:g.100G>A;NC_000002.11:g.101A>C]");
    assert!(result.contains("NC_000001.11"), "got {}", result);
    assert!(result.contains("NC_000002.11"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_across_coord_systems_same_accession() {
    // Same accession (`NM_000088.3`), different coord types (c. vs n.),
    // numerically adjacent positions. The merge gate is shared accession
    // *and* shared coord type — this case must NOT collapse.
    let result = normalize_to_string("[NM_000088.3:c.1A>G;NM_000088.3:n.2A>G]");
    assert!(result.contains("c.1A>G"), "got {}", result);
    assert!(result.contains("n.2A>G"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_across_coord_systems_different_accession() {
    // c. + g. with numerically adjacent positions on different
    // accessions. Two reasons to refuse merge (different accession,
    // different coord system) — verify it stays separated.
    let result = normalize_to_string("[NM_000088.3:c.1A>G;NC_000017.11:g.2A>G]");
    assert!(result.contains("NM_000088.3:c.1A>G"), "got {}", result);
    assert!(result.contains("NC_000017.11:g.2A>G"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_across_coord_systems_three_variant() {
    // Three-variant chain mixing c./g./p. — even though c.1 and (notional)
    // g.2 / p.Met2Val are numerically "adjacent" in a numeric sense, the
    // merge pass cannot bridge the coord-system boundary.
    let result =
        normalize_to_string("[NM_000088.3:c.1A>G;NC_000017.11:g.2A>G;NP_000079.2:p.Met2Val]");
    assert!(result.contains("NM_000088.3:c.1A>G"), "got {}", result);
    assert!(result.contains("NC_000017.11:g.2A>G"), "got {}", result);
    assert!(result.contains("NP_000079.2:p.Met2Val"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_normalize_preserves_allele_wrapper_for_cross_reference() {
    // Sanity check the post-normalize variant kind: a cross-reference
    // cis allele is non-mergeable, so the `Allele` wrapper must
    // survive normalization (cf. `test_singleton_cis_allele_preserves_wrapper`
    // in `merge_consecutive_edits_tests.rs` for the singleton case).
    let normalizer = Normalizer::new(MockProvider::new());
    let parsed = parse_hgvs("[NM_000088.3:c.1A>G;NC_000017.11:g.2A>G]").expect("parse");
    let normalized = normalizer.normalize(&parsed).expect("normalize");
    assert!(
        matches!(normalized, HgvsVariant::Allele(_)),
        "expected Allele preserved after cross-ref normalize, got {:?}",
        normalized,
    );
}

// =====================================================================
// 7. Compact-vs-expanded display contract.
//
// Pins that `Display` only emits the compact `ACC:type.[edit;edit]`
// form when every sub-variant shares both accession and coord type;
// otherwise it falls back to the per-variant expanded form. This is
// the canonicalization decision the issue calls "undefined" — record
// it so future refactors don't silently change shape.
// =====================================================================

#[test]
fn test_display_uses_compact_form_only_when_shared() {
    // Same accession + same type -> compact.
    let compact = format!(
        "{}",
        parse_hgvs("[NM_000088.3:c.1A>G;NM_000088.3:c.2A>G]").expect("parse"),
    );
    assert_eq!(compact, "NM_000088.3:c.[1A>G;2A>G]");

    // Different accession (any types) -> expanded.
    let expanded_acc = format!(
        "{}",
        parse_hgvs("[NM_000088.3:c.1A>G;NM_000089.3:c.2A>G]").expect("parse"),
    );
    assert_eq!(expanded_acc, "[NM_000088.3:c.1A>G;NM_000089.3:c.2A>G]");

    // Same accession, different type -> expanded.
    let expanded_type = format!(
        "{}",
        parse_hgvs("[NM_000088.3:c.1A>G;NM_000088.3:n.2A>G]").expect("parse"),
    );
    assert_eq!(expanded_type, "[NM_000088.3:c.1A>G;NM_000088.3:n.2A>G]");

    // Cross-reference (different accession AND different type) -> expanded.
    let expanded_xref = format!(
        "{}",
        parse_hgvs("[NM_000088.3:c.1A>G;NC_000017.11:g.2A>G]").expect("parse"),
    );
    assert_eq!(expanded_xref, "[NM_000088.3:c.1A>G;NC_000017.11:g.2A>G]");
}

// =====================================================================
// 8. Coord-system pair coverage matrix.
//
// Section 3 covers c.+g., c.+r., c.+p., c.+n. (same accession). The
// HGVS spec defines seven coord systems (g., c., n., r., p., m., o.).
// This section pins the remaining cross-coord pairs that ferro accepts,
// so the audit covers every pair the parser sees in the wild.
// =====================================================================

#[test]
fn test_cis_mixed_g_and_m() {
    // g. + m. (genomic + mitochondrial). NC_012920.1 is the canonical
    // mtDNA accession; ferro keeps the m. coord type even though the
    // accession could theoretically be reused for a nuclear g. variant.
    let input = "[NC_000001.11:g.100A>G;NC_012920.1:m.100A>G]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert_eq!(allele.variants[0].variant_type(), "g");
    assert_eq!(allele.variants[1].variant_type(), "m");
}

#[test]
fn test_cis_mixed_c_and_m() {
    let input = "[NM_000088.3:c.1A>G;NC_012920.1:m.100A>G]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert_eq!(allele.variants[0].variant_type(), "c");
    assert_eq!(allele.variants[1].variant_type(), "m");
}

#[test]
fn test_cis_mixed_g_and_o() {
    // g. + o. (genomic + circular DNA, SVD-WG006).
    let input = "[NC_000001.11:g.100A>G;NC_001802.1:o.100A>G]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert_eq!(allele.variants[0].variant_type(), "g");
    assert_eq!(allele.variants[1].variant_type(), "o");
}

#[test]
fn test_cis_mixed_c_and_o() {
    let input = "[NM_000088.3:c.1A>G;NC_001802.1:o.100A>G]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert_eq!(allele.variants[0].variant_type(), "c");
    assert_eq!(allele.variants[1].variant_type(), "o");
}

#[test]
fn test_cis_mixed_n_and_r() {
    // n. + r. (non-coding transcript + RNA on different accessions).
    let input = "[NR_000001.1:n.100A>G;NM_000088.3:r.50a>g]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert_eq!(allele.variants[0].variant_type(), "n");
    assert_eq!(allele.variants[1].variant_type(), "r");
}

#[test]
fn test_cis_mixed_n_and_p() {
    let input = "[NR_000001.1:n.100A>G;NP_000079.2:p.Met2Val]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert_eq!(allele.variants[0].variant_type(), "n");
    assert_eq!(allele.variants[1].variant_type(), "p");
}

#[test]
fn test_cis_mixed_r_and_p() {
    let input = "[NM_000088.3:r.50a>g;NP_000079.2:p.Met2Val]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert_eq!(allele.variants[0].variant_type(), "r");
    assert_eq!(allele.variants[1].variant_type(), "p");
}

#[test]
fn test_cis_mixed_m_and_r() {
    // Rare but spec-allowed: m. + r. (mtDNA + RNA on different
    // accessions). Pin that ferro doesn't trip on the type combination.
    let input = "[NC_012920.1:m.100A>G;NM_000088.3:r.50a>g]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert_eq!(allele.variants[0].variant_type(), "m");
    assert_eq!(allele.variants[1].variant_type(), "r");
}

#[test]
fn test_cis_mixed_m_and_p() {
    let input = "[NC_012920.1:m.100A>G;NP_000079.2:p.Met2Val]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 2);
    assert_eq!(allele.variants[0].variant_type(), "m");
    assert_eq!(allele.variants[1].variant_type(), "p");
}

// =====================================================================
// 9. Idempotency under multi-pass normalize.
//
// Normalize-then-renormalize must be a fixed point for cross-reference
// compounds. This guards against a class of round-trip bugs where the
// first pass shifts a variant and the second pass shifts it back (or
// further), changing shape on every iteration. Pins both the rendered
// string and the AlleleVariant equality.
// =====================================================================

#[test]
fn test_normalize_idempotent_cross_reference_cis() {
    let normalizer = Normalizer::new(MockProvider::new());
    let input = "[NM_000088.3:c.1A>G;NC_000017.11:g.2A>G]";
    let n1 = normalizer
        .normalize(&parse_hgvs(input).expect("parse"))
        .expect("normalize");
    let s1 = format!("{}", n1);
    let n2 = normalizer
        .normalize(&parse_hgvs(&s1).expect("re-parse"))
        .expect("re-normalize");
    let s2 = format!("{}", n2);
    assert_eq!(s1, s2, "normalize is not idempotent: {} -> {}", s1, s2);
    assert_eq!(
        n1, n2,
        "normalize produced unequal AlleleVariants on re-pass",
    );
}

#[test]
fn test_normalize_idempotent_cross_reference_trans() {
    let normalizer = Normalizer::new(MockProvider::new());
    let input = "[NM_000088.3:c.1A>G];[NC_000017.11:g.2A>G]";
    let n1 = normalizer
        .normalize(&parse_hgvs(input).expect("parse"))
        .expect("normalize");
    let s1 = format!("{}", n1);
    let n2 = normalizer
        .normalize(&parse_hgvs(&s1).expect("re-parse"))
        .expect("re-normalize");
    assert_eq!(s1, format!("{}", n2));
    assert_eq!(n1, n2);
}

#[test]
fn test_normalize_idempotent_with_resolvable_companion() {
    // Per-variant normalize 3'-shifts the c. dup all the way to the end
    // of the `CCCCC` run at positions 10..14 of the test transcript, so
    // a single full normalization pass produces `c.14dup`. From there
    // both subsequent passes are fixed points. Pins idempotency in the
    // presence of mixed resolvable / un-resolvable accessions.
    //
    // (Using `c.10dupA` here would NOT be idempotent: the first pass
    // canonicalizes the explicit-base form to position-only `c.10dup`,
    // and the second pass then 3'-shifts to `c.14dup`. That two-pass
    // gap is a single-variant normalization concern, orthogonal to the
    // H2 cross-reference contract this test is auditing.)
    let normalizer = Normalizer::new(provider_with_simple_transcript());
    let input = "[NM_TEST.1:c.10dup;NC_000001.11:g.100A>G]";
    let n1 = normalizer
        .normalize(&parse_hgvs(input).expect("parse"))
        .expect("normalize");
    let s1 = format!("{}", n1);
    assert_eq!(s1, "[NM_TEST.1:c.14dup;NC_000001.11:g.100A>G]");
    let n2 = normalizer
        .normalize(&parse_hgvs(&s1).expect("re-parse"))
        .expect("re-normalize");
    assert_eq!(s1, format!("{}", n2));
    assert_eq!(n1, n2);
}

// =====================================================================
// 10. Hash and Eq stability.
//
// Two parses of the same cross-reference compound string must produce
// equal AlleleVariant values that hash equal. Reordered or
// phase-different variants must NOT hash equal — order and phase are
// significant per HGVS, and `HgvsVariant` is used as a HashMap key
// in normalization caches.
// =====================================================================

#[test]
fn test_hash_eq_same_input() {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    let s = "[NM_000088.3:c.1A>G;NC_000017.11:g.2A>G]";
    let a = parse_hgvs(s).expect("parse a");
    let b = parse_hgvs(s).expect("parse b");
    assert_eq!(a, b, "two parses of the same string must compare equal");

    let mut ha = DefaultHasher::new();
    a.hash(&mut ha);
    let mut hb = DefaultHasher::new();
    b.hash(&mut hb);
    assert_eq!(
        ha.finish(),
        hb.finish(),
        "Hash must agree with Eq for identical parses",
    );
}

#[test]
fn test_hash_differs_on_reorder() {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    let forward = parse_hgvs("[NM_000088.3:c.1A>G;NC_000017.11:g.2A>G]").expect("parse");
    let reverse = parse_hgvs("[NC_000017.11:g.2A>G;NM_000088.3:c.1A>G]").expect("parse");
    let mut hf = DefaultHasher::new();
    forward.hash(&mut hf);
    let mut hr = DefaultHasher::new();
    reverse.hash(&mut hr);
    // Hash collisions are theoretically allowed, but for this small
    // structural reordering they must differ in practice; if this ever
    // flakes the parser changed shape and we should re-pin.
    assert_ne!(
        hf.finish(),
        hr.finish(),
        "reordered cross-ref alleles must hash differently",
    );
}

#[test]
fn test_cis_and_trans_not_equal() {
    // Same sub-variants, different phase -> different HgvsVariant value.
    let cis = parse_hgvs("[NM_000088.3:c.1A>G;NC_000017.11:g.2A>G]").expect("parse cis");
    let trans = parse_hgvs("[NM_000088.3:c.1A>G];[NC_000017.11:g.2A>G]").expect("parse trans");
    assert_ne!(cis, trans, "cis and trans must not compare equal");
}

// =====================================================================
// 11. Edit-type matrix across references.
//
// Existing tests are SNV-heavy. Pin that the cross-reference shape is
// preserved across the broader edit-type set ferro supports (del, dup,
// inv, delins, repeat). Each test mixes coord systems to keep the
// cross-reference contract in scope.
// =====================================================================

#[test]
fn test_cross_reference_del_and_dup() {
    let input = "[NM_000088.3:c.10_15del;NC_000017.11:g.20dup]";
    let parsed = assert_parse_roundtrip(input);
    expect_allele(&parsed, AllelePhase::Cis, 2);
}

#[test]
fn test_cross_reference_inv_and_delins() {
    let input = "[NM_000088.3:c.10_15inv;NC_000017.11:g.20_25delinsAA]";
    let parsed = assert_parse_roundtrip(input);
    expect_allele(&parsed, AllelePhase::Cis, 2);
}

#[test]
fn test_cross_reference_repeat_edits() {
    // Repeat notation (`A[5]`, `CTG[3]`) on each side of a cross-ref
    // pair. The trailing `]` on a repeat is inside the sub-variant; the
    // outer allele bracket terminates only after both repeat blocks.
    let input = "[NM_000088.3:c.5_7CTG[3];NC_000017.11:g.10A[5]]";
    let parsed = assert_parse_roundtrip(input);
    expect_allele(&parsed, AllelePhase::Cis, 2);
}

#[test]
fn test_cross_reference_normalize_preserves_edit_types() {
    // Round-trip through the normalizer for a mixed edit-type pair; no
    // edit type may be silently rewritten to another (e.g., dup must
    // stay dup, not become a delins, when normalization isn't merging).
    let result = normalize_to_string("[NM_000088.3:c.10_15del;NC_000017.11:g.20dup]");
    assert!(result.contains("c.10_15del"), "got {}", result);
    assert!(result.contains("g.20dup"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

// =====================================================================
// 12. Multi-reference haplotype + trans display fidelity.
//
// The DNA/alleles spec gives an example of a multi-reference haplotype
// with four different genomic accessions all sharing coord type `g.`
// (`docs/recommendations/DNA/alleles.md` line 84). Pin that the expanded
// form is emitted (because accessions differ) and that trans cross-ref
// renders as `[a];[b]`, never `[a;b]`.
// =====================================================================

#[test]
fn test_multi_reference_haplotype_round_trip() {
    // Plain SNVs (the spec's full example uses repeat notation; we strip
    // that to keep the test focused on cross-reference shape).
    let input = "[M59228.1:g.250G>C;AF209160.1:g.572C>T;Z11861.1:g.61T>C;Z16803.1:g.114T>A]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Cis, 4);
    // All sub-variants share coord type `g.` but accessions differ -> expanded form.
    for v in &allele.variants {
        assert_eq!(v.variant_type(), "g");
    }
    let display = format!("{}", parsed);
    assert!(
        !display.starts_with("M59228.1:g.["),
        "expected expanded form, got compact: {}",
        display,
    );
    assert!(display.starts_with('['), "got {}", display);
}

#[test]
fn test_trans_cross_reference_display_uses_split_brackets() {
    // Pin that trans phase emits `[a];[b]` (split brackets), never the
    // cis shape `[a;b]`. Distinguishes phase visually for downstream
    // tools.
    let input = "[NM_000088.3:c.1A>G];[NC_000017.11:g.2A>G]";
    let display = format!("{}", parse_hgvs(input).expect("parse"));
    assert_eq!(display, input);
    assert!(display.contains("];["), "got {}", display);
}

#[test]
fn test_trans_cross_reference_normalize_preserves_split_brackets() {
    // Same as above, but post-normalize. Pins that the trans wrapper
    // survives the normalize pass for cross-reference compounds (which
    // are not merge-eligible — the merge guard requires
    // shared accession + shared type).
    let normalized = normalize_to_string("[NM_000088.3:c.1A>G];[NC_000017.11:g.2A>G]");
    assert_eq!(normalized, "[NM_000088.3:c.1A>G];[NC_000017.11:g.2A>G]");
    assert!(normalized.contains("];["), "got {}", normalized);
}

// =====================================================================
// 13. Compact-form rejection and spec-discouraged acceptance.
//
// The `ACC:type.[edit;edit]` compact form is a strict optimization for
// the shared-accession-and-type case. When a sub-variant inside the
// brackets carries its own accession (e.g.
// `NM_X:c.[1A>G;NC_Y:g.2A>G]`), parse must fail — the compact form
// inherits the prefix and a sub-variant cannot redefine it.
//
// Conversely, the spec (DNA/alleles.md line 22) flags
// `c.[76A>C];g.[10091C>G]` as discouraged. The expanded equivalent
// `[NM_X:c.76A>C];[NC_Y:g.10091C>G]` is syntactically valid HGVS and
// ferro accepts it. Pin both behaviors so the policy is recorded.
// =====================================================================

#[test]
fn test_compact_form_rejects_mixed_prefix_inside_brackets() {
    let input = "NM_000088.3:c.[1A>G;NC_000017.11:g.2A>G]";
    assert!(
        parse_hgvs(input).is_err(),
        "expected parse error for compact form with mixed prefix, parsed OK: {}",
        input,
    );
}

#[test]
fn test_compact_form_rejects_bare_with_mixed_prefix() {
    // Bare compact form (no leading accession) with mixed prefix is
    // also invalid — the parser cannot infer an accession.
    let input = "c.[1A>G;NC_000017.11:g.2A>G]";
    assert!(
        parse_hgvs(input).is_err(),
        "expected parse error for bare compact form with mixed prefix, parsed OK: {}",
        input,
    );
}

#[test]
fn test_spec_discouraged_cross_reference_trans_accepted() {
    // `[NM_X:c.76A>C];[NC_Y:g.10091C>G]` is the spec-discouraged shape
    // in a syntactically valid wrapper. Pin that ferro accepts it and
    // does NOT silently merge / canonicalize / reject it — round-trips
    // verbatim.
    let input = "[NM_000088.3:c.76A>C];[NC_000017.11:g.10091C>G]";
    let parsed = assert_parse_roundtrip(input);
    expect_allele(&parsed, AllelePhase::Trans, 2);
}

// =====================================================================
// 14. NullAllele / UnknownAllele in cross-reference compounds and
//     3-way unknown-phase mixes.
//
// `[a];[0]` is hemizygous (one X-chromosome present), `[a];[?]` is one
// allele identified, the other expected but not yet found. These are
// spec-defined markers (DNA/alleles.md notes on second-allele
// representation). They must work when paired with a cross-reference
// companion variant.
// =====================================================================

#[test]
fn test_trans_cross_reference_with_null_allele() {
    let input = "[NM_000088.3:c.1A>G];[0]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Trans, 2);
    assert!(matches!(allele.variants[0], HgvsVariant::Cds(_)));
    assert!(matches!(allele.variants[1], HgvsVariant::NullAllele));
}

#[test]
fn test_trans_cross_reference_with_unknown_allele() {
    let input = "[NM_000088.3:c.1A>G];[?]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Trans, 2);
    assert!(matches!(allele.variants[0], HgvsVariant::Cds(_)));
    assert!(matches!(allele.variants[1], HgvsVariant::UnknownAllele));
}

#[test]
fn test_trans_cross_reference_marker_normalize_preserves_shape() {
    // Markers ride through normalize untouched; the c. partner
    // normalizes (or no-ops) but the marker leg stays as `[0]` / `[?]`.
    let null_normalized = normalize_to_string("[NM_000088.3:c.1A>G];[0]");
    assert_eq!(null_normalized, "[NM_000088.3:c.1A>G];[0]");
    let unk_normalized = normalize_to_string("[NM_000088.3:c.1A>G];[?]");
    assert_eq!(unk_normalized, "[NM_000088.3:c.1A>G];[?]");
}

#[test]
fn test_unknown_phase_three_way_cross_coord_round_trip() {
    // Three sub-variants with three coord systems in unknown phase.
    let input = "[NM_000088.3:c.1A>G(;)NC_000017.11:g.2A>G(;)NP_000079.2:p.Met2Val]";
    let parsed = assert_parse_roundtrip(input);
    let allele = expect_allele(&parsed, AllelePhase::Unknown, 3);
    assert_eq!(allele.variants[0].variant_type(), "c");
    assert_eq!(allele.variants[1].variant_type(), "g");
    assert_eq!(allele.variants[2].variant_type(), "p");
}
