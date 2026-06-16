//! Tests for issue #165 / tracking issue #81 item A10 (sub-only branch).
//!
//! A `delins` whose post-trim span contains two or more independent
//! single-base mismatches separated by at least one unchanged nucleotide
//! must decompose to the individual `[X>A; B>Y]` form per the HGVS
//! edit-priority rule (`general.md:56`: substitution > deletion >
//! inversion > duplication > insertion; `delins` is the residual when no
//! higher-priority form applies). Adjacent (no-gap) mismatches stay as a
//! delins per `substitution.md` ("two or more consecutive nucleotides are
//! described as deletion/insertion") — that run-merge is handled by the
//! `build_split_variants` adjacency rule introduced in #182.
//!
//! Spec exception (`general.md:35-38`): in a coding sequence, two
//! variants separated by exactly one nucleotide that together affect a
//! single codon are described as a `delins`. The decomposition preserves
//! that exact `[Sub; Identity; Sub]` triplet as a 3-base delins even when
//! it is embedded inside a longer span — the codon-frame merge from #79
//! survives A10 round-trips.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Test transcript used by every test in this file.
///
/// ```text
///   c. axis: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ...
///   base:    A T G C A A A A A  C  C  C  C  C  G  G  G  G  G  T ...
/// ```
///
/// CDS = c.1..=c.60. Codons (1-indexed, formula `(base - 1) / 3 + 1`):
/// codon 4 = c.10..c.12, codon 5 = c.13..c.15. The codon-frame exception
/// is exercised against this transcript by placing a pair of SNVs inside
/// codon 4 (c.10, c.12) and across the codon-4 / codon-5 boundary (c.12,
/// c.14) for the negative case.
fn provider_simple() -> MockProvider {
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

fn normalize_with(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

// =============================================================================
// Genomic (`g.`) — no codon frame.
// =============================================================================

#[test]
fn genomic_delins_with_gap_decomposes() {
    // Genomic variants have no codon frame, so any `[Sub; Identity; Sub]`
    // decomposition emits two separate subs per `general.md:34`. Build a
    // synthetic chromosome where g.100..102 resolves to `CAC`.
    let mut provider = MockProvider::new();
    // 0-based slicing reads g.100..102 from bytes[99..102]; the leading
    // filler keeps the substantive bytes at known indices.
    let mut seq = String::from_utf8(vec![b'N'; 99]).unwrap();
    seq.push_str("CAC");
    seq.push_str(&String::from_utf8(vec![b'N'; 100]).unwrap());
    provider.add_genomic_sequence("NC_000001.11", seq);
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NC_000001.11:g.100_102delinsTAG").expect("parse failed");
    let result = format!(
        "{}",
        normalizer.normalize(&variant).expect("normalize failed")
    );
    assert_eq!(result, "NC_000001.11:g.[100C>T;102C>G]");
}

// =============================================================================
// Tx (`n.`) — no codon frame, gap-of-1 delins must decompose.
// =============================================================================

#[test]
fn tx_delins_with_gap_decomposes_to_individual_subs() {
    // Literal 3-base `n.` delins whose middle base equals the reference.
    // The user may have entered this either directly or via cis-allele
    // merging (`n.[10C>A; 12C>A]` merges to `n.10_12delinsACA`). Either
    // way, the post-merge canonical form per `general.md:34` (variants
    // separated by ≥1 unchanged nt are described individually) is two
    // separate subs.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:n.10_12delinsACA"),
        "NM_TEST.1:n.[10C>A;12C>A]",
    );
}

#[test]
fn tx_cis_allele_with_gap_decomposes() {
    // Bracketed-cis form should produce the same canonical decomposition
    // as the literal delins above. The cis-allele merge does not fire
    // (gap-of-1 is not strictly-adjacent and `n.` has no codon-frame
    // exception), so each sub passes through normalize independently.
    assert_eq!(
        normalize_with(provider_simple(), "[NM_TEST.1:n.10C>A;NM_TEST.1:n.12C>A]",),
        "NM_TEST.1:n.[10C>A;12C>A]",
    );
}

// =============================================================================
// RNA (`r.`) — codon-frame aware (issue #275 item 1), lowercase preserved.
// =============================================================================

#[test]
fn rna_delins_with_gap_preserves_codon_frame_triplet_lowercase() {
    // r. shares the CDS-relative axis with c. and follows the same
    // codon-frame exception (issue #275 item 1). A `[Sub; Identity; Sub]`
    // triplet whose endpoints share a codon stays as a single 3-base
    // delins under the spec's `general.md:35-38` carve-out. Codon 4
    // covers r.10..r.12, so the input round-trips. RNA display preserves
    // lowercase on the user-declared alt bases.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:r.10_12delinsaca"),
        "NM_TEST.1:r.10_12delinsaca",
    );
}

#[test]
fn rna_delins_with_gap_decomposes_when_pair_straddles_codon_boundary() {
    // r.12 sits in codon 4 and r.14 in codon 5, so the codon-frame
    // exception does NOT apply. The `[Sub; Identity; Sub]` triplet must
    // decompose into two separate subs per `general.md:34` (variants
    // separated by ≥1 unchanged nt are described individually). Ref at
    // r.12..r.14 is `CCC`; insert `aca` produces `c>a` at 12 and `c>a`
    // at 14 with an unchanged middle.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:r.12_14delinsaca"),
        "NM_TEST.1:r.[12c>a;14c>a]",
    );
}

// =============================================================================
// Mitochondrial (`m.`) — no codon frame.
// =============================================================================

#[test]
fn mt_delins_with_gap_decomposes() {
    // The mitochondrial coordinate system has no codon-frame exception;
    // `m.` follows the genomic / `n.` decomposition rule.
    //
    // Build a 12-byte synthetic mt sequence so m.10..12 resolves to
    // ref="CAC". `add_genomic_sequence` stores the contig as-is and
    // `get_genomic_sequence` reads a 0-based half-open slice — m.10..12
    // is `bytes[9..12]`.
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_012920.1", "NNNNNNNNNCAC");
    assert_eq!(
        normalize_with(provider, "NC_012920.1:m.10_12delinsTAG"),
        "NC_012920.1:m.[10C>T;12C>G]",
    );
}

// =============================================================================
// Coding (`c.`) — codon-frame exception preserves the [Sub; Identity; Sub] pair.
// =============================================================================

#[test]
fn cds_codon_frame_pair_preserved_as_delins() {
    // c.10 and c.12 are in codon 4 (`(base - 1) / 3` = 3 for both). The
    // spec exception (`general.md:35-38`) keeps the pair as a 3-base
    // delins. Equivalent to issue #79's codon-frame merge — A10 must not
    // regress this.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.10_12delinsTCA"),
        "NM_TEST.1:c.10_12delinsTCA",
    );
}

#[test]
fn cds_codon_frame_pair_via_cis_allele_preserved() {
    // The merge stage of `normalize_allele` codon-frame-merges
    // `c.[10C>T;12C>A]` into `c.10_12delinsTCA`. A10's decomposition
    // recognises the [Sub; Identity; Sub] triplet as the codon-frame
    // exception and keeps it as a single delins.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.[10C>T;12C>A]"),
        "NM_TEST.1:c.10_12delinsTCA",
    );
}

#[test]
fn cds_cross_codon_pair_decomposes() {
    // c.12 is in codon 4, c.14 is in codon 5 — the pair does not "together
    // affect one amino acid", so the codon-frame exception does not
    // apply. Per `general.md:34` the canonical form is two separate subs.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.12_14delinsACA"),
        "NM_TEST.1:c.[12C>A;14C>A]",
    );
}

#[test]
fn cds_embedded_codon_frame_triplet_preserved_within_longer_span() {
    // c.10..c.13 spans codon 4 (10-12) and the first base of codon 5
    // (13). After merging `c.[10C>T;12C>A;13C>G]` the cis-allele path
    // produces `c.10_13delinsTCAG` (codon-frame merge on the first pair
    // then strict-adjacent chain with the third sub).
    //
    // The principled decomposition: the (10, 12) pair satisfies the
    // codon-frame exception and emits as a 3-base delins; c.13 is in a
    // different codon and emits as a separate sub. The chained 4-base
    // delins is non-canonical (spec exception is defined for a *pair*,
    // not a chain).
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.[10C>T;12C>A;13C>G]"),
        "NM_TEST.1:c.[10_12delinsTCA;13C>G]",
    );
}

#[test]
fn cds_embedded_codon_frame_triplet_from_literal_delins() {
    // Same shape as the cis-allele case above but entered as a literal
    // delins. ref c.10..c.13 = CCCC. alt = TCAG matches the [Sub;
    // Identity; Sub; Sub] decomposition: the (10, 12) pair preserves as
    // a codon-frame delins and c.13 splits off as a single sub.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.10_13delinsTCAG"),
        "NM_TEST.1:c.[10_12delinsTCA;13C>G]",
    );
}

#[test]
fn cds_two_codon_frame_triplets_back_to_back() {
    // ref c.10..c.15 = CCCCCG (codon 4 = c.10..c.12, codon 5 = c.13..c.15).
    // alt = TCAGCT decomposes to:
    //   pos 10 C>T   (sub)      — codon 4
    //   pos 11 C=C   (identity) — codon 4
    //   pos 12 C>A   (sub)      — codon 4
    //   pos 13 C>G   (sub)      — codon 5
    //   pos 14 C=C   (identity) — codon 5
    //   pos 15 G>T   (sub)      — codon 5
    //
    // Both (10, 12) and (13, 15) match the codon-frame exception, so
    // each triplet emits as its own 3-base delins. The two delins are
    // strictly adjacent (12 + 1 == 13) but do not re-merge: each is
    // emitted independently by `build_split_variants`'s triplet
    // lookahead, and the outer cis-allele wrapper preserves them.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.10_15delinsTCAGCT"),
        "NM_TEST.1:c.[10_12delinsTCA;13_15delinsGCT]",
    );
}

// =============================================================================
// Adjacent (no-gap) substitutions must stay as a single delins (#182).
// A10's loosened gate must NOT decompose pure-substitution runs.
// =============================================================================

#[test]
fn cds_adjacent_sub_pair_stays_as_delins() {
    // c.[10C>T;11C>A]: strictly-adjacent merge produces
    // `c.10_11delinsTA`. The post-merge decomposition would emit
    // `[Sub@0; Sub@1]` (no Identity, no Inversion); per the
    // `(has_inv || has_identity)` guard `decompose_delins` returns
    // `None`, so `build_split_variants` is never called and the delins
    // stays intact — matching `substitution.md`'s rule that consecutive
    // nucleotide changes are described as `delins` (#182).
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.[10C>T;11C>A]"),
        "NM_TEST.1:c.10_11delinsTA",
    );
}

#[test]
fn cds_three_adjacent_subs_stay_as_delins() {
    // Same principle as above for a 3-base adjacent run.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.10_12delinsTAG"),
        "NM_TEST.1:c.10_12delinsTAG",
    );
}

// =============================================================================
// Round-trip stability — re-normalizing the decomposed form must reproduce it.
// =============================================================================

#[test]
fn round_trip_decomposed_form_is_stable() {
    // n. decomposed pair re-parses as `n.[10C>A;12C>A]` and re-normalizes
    // to itself. The strict-adjacent merge does not fire (gap-of-1), so
    // no merge round-trips back to a delins.
    let once = normalize_with(provider_simple(), "NM_TEST.1:n.10_12delinsACA");
    let twice = normalize_with(provider_simple(), &once);
    assert_eq!(once, twice);
    assert_eq!(once, "NM_TEST.1:n.[10C>A;12C>A]");
}

#[test]
fn round_trip_codon_frame_pair_is_stable() {
    // c. codon-frame pair re-parses and re-normalizes to itself.
    let once = normalize_with(provider_simple(), "NM_TEST.1:c.10_12delinsTCA");
    let twice = normalize_with(provider_simple(), &once);
    assert_eq!(once, twice);
    assert_eq!(once, "NM_TEST.1:c.10_12delinsTCA");
}

#[test]
fn round_trip_embedded_triplet_is_stable() {
    // c.[10_12delinsTCA;13C>G] re-merges via the cis-allele path
    // (strict-adjacent merge of the trailing sub onto the delins
    // produces `c.10_13delinsTCAG`) and then re-decomposes to the same
    // canonical form.
    let once = normalize_with(provider_simple(), "NM_TEST.1:c.10_13delinsTCAG");
    let twice = normalize_with(provider_simple(), &once);
    assert_eq!(once, twice);
    assert_eq!(once, "NM_TEST.1:c.[10_12delinsTCA;13C>G]");
}
