//! Regression coverage for issue #382: an `Insertion`/`Delins` whose
//! inserted sequence matches a contiguous downstream reference window
//! must be canonicalised to `Duplication` of that window per the HGVS
//! `dup > ins/delins` priority rule, with the resulting dup 3'-shifted.
//!
//! Anchor case from the biocommons-normalize corpus (#324):
//!
//!   NC_000001.10:g.1647893delinsCTTTCTT
//!     biocommons: g.1647895_1647900dup
//!     ferro (pre-#356): g.1647893_1647894insTTTCTT
//!     ferro (current):  g.1647895_1647900dup  ✓
//!
//! The fix landed in PR #356 (the delins → ins recursion in
//! `normalize_na_edit`, which routes the shared-affix-trimmed insertion
//! back through the ins-canonicalisation pipeline including
//! `insertion_to_duplication`). PR #346 contributed the direction-aware
//! `InsToDupResult` plumbing that the (c) Single-copy fallback in
//! `mod.rs:~3596` relies on.
//!
//! This file pins the resulting behavior with synthetic MockProvider
//! tests so future refactors of the canonicalisation chain don't
//! silently regress the spec rule.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

// =============================================================================
// Genomic axis — the issue's anchor case shape
// =============================================================================

/// Synthetic genomic contig "TEST.1", padded to 300 bases so the
/// normalize window-fetch (default `window_size = 100`) succeeds even
/// when the variant sits near the start of the contig.
///
/// The substantive bases at 1-based positions 105..119:
///
///   pos: 105 106 107 108 109 110 111 112 113 114
///   seq:  C   T   T   T   C   T   T   T   A   G
///         ^   \---- ref ----/                ^
///         |   matches alt "TTTCTT" at        |
///         delins  positions 106..111         shift-
///         base    plus a trailing "T" at 112 stopper
///                 enables one 3' unit-shift  (≠ T)
///                 of the dup from [106..111]
///                 to [107..112].
///
/// `g.105delinsCTTTCTT` has shared prefix "C" with ref[105], so after
/// shared-affix trim it is `ins TTTCTT` between g.105 and g.106. The
/// 6-base alt "TTTCTT" matches ref[106..111] = "TTTCTT", so the
/// spec-canonical form is `dup` of that window. The 3'-rule unit-shift
/// proceeds: ref[106] = T = ref[112] permits the dup to walk one base
/// right to `g.107_112dup`; ref[107] = T ≠ ref[113] = A stops the walk
/// there. Final canonical output: `TEST.1:g.107_112dup`. This mirrors
/// `NC_000001.10:g.1647893delinsCTTTCTT` → `g.1647895_1647900dup` in
/// shape.
fn genomic_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    // Bytes 0..99 are filler. Bytes 99..119 hold the substantive
    // pattern "GGGGGCTTTCTTTAGGGGGG" (1-based positions 100..119).
    // Pad to length 300 so the default 100-base normalize window
    // always lands inside the contig.
    let mut seq = String::from_utf8(vec![b'G'; 99]).unwrap();
    seq.push_str("GGGGGCTTTCTTTAGGGGGG");
    while seq.len() < 300 {
        seq.push('G');
    }
    provider.add_genomic_sequence("TEST.1", seq);
    provider
}

fn normalize_genomic(input: &str) -> String {
    let normalizer = Normalizer::with_config(
        genomic_provider(),
        NormalizeConfig::default()
            .with_direction(ShuffleDirection::ThreePrime)
            .allow_crossing_boundaries(),
    );
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{}", normalized)
}

#[test]
fn delins_with_alt_matching_downstream_tract_emits_dup() {
    // The anchor case for #382. Equivalent edit shape to
    // NC_000001.10:g.1647893delinsCTTTCTT → g.1647895_1647900dup.
    // Exercises the canonicalize_delins → DelinsCanonical::Insertion →
    // recurse → insertion_to_duplication path.
    assert_eq!(
        normalize_genomic("TEST.1:g.105delinsCTTTCTT"),
        "TEST.1:g.107_112dup",
    );
}

#[test]
fn insertion_with_alt_matching_downstream_tract_emits_dup() {
    // Same edit re-expressed as a literal Insertion (the post-delins-
    // trim shape). Isolates the ins → dup canon path from the
    // delins → ins recursion at mod.rs:3153.
    assert_eq!(
        normalize_genomic("TEST.1:g.105_106insTTTCTT"),
        "TEST.1:g.107_112dup",
    );
}

#[test]
fn insertion_explicit_range_delins_form_emits_dup() {
    // Same edit expressed as an explicit single-base-range delins, to
    // pin parser-shape invariance through canonicalize_delins.
    assert_eq!(
        normalize_genomic("TEST.1:g.105_105delinsCTTTCTT"),
        "TEST.1:g.107_112dup",
    );
}

#[test]
fn insertion_with_no_adjacent_tract_match_stays_as_ins() {
    // Negative regression: an insertion whose alt does not match any
    // upstream OR downstream window in the reference must stay as an
    // Insertion. Synthetic contig provides only "GGGGG...G" around the
    // ins point, so the alt "TTTCTT" has no adjacent tract.
    let mut provider = MockProvider::new();
    let mut seq = String::from_utf8(vec![b'G'; 200]).unwrap();
    while seq.len() < 300 {
        seq.push('G');
    }
    provider.add_genomic_sequence("TEST_NEG.1", seq);
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default()
            .with_direction(ShuffleDirection::ThreePrime)
            .allow_crossing_boundaries(),
    );
    let variant = parse_hgvs("TEST_NEG.1:g.105_106insTTTCTT").expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    assert_eq!(format!("{}", normalized), "TEST_NEG.1:g.105_106insTTTCTT");
}

// =============================================================================
// c.-axis — the issue's "Tests" bullet point shape
// =============================================================================

/// Synthetic transcript "NM_TEST382.1", 200 bases CDS-only.
///
/// The substantive bases at c. positions 5..15:
///
///   c.:  5 6 7 8 9 10 11 12 13 14 15
///   seq: C T T T C T  T  T  A  G  G
///
/// `c.5_6insTTTCTT` (a literal insertion shape) has the same
/// "alt matches a downstream window" structure: ref[c.6..c.11] =
/// "TTTCTT" matches the alt, and ref[c.6] = T = ref[c.12] permits one
/// 3' unit-shift to `c.7_12dup`; ref[c.13] = A stops further shifting.
fn cds_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    // Build a 200-base CDS-only transcript with the substantive
    // pattern at c.5..c.14. Padding around lets the normalize window
    // (default 100) stay inside the sequence regardless of where the
    // shift terminates.
    let mut seq = String::from("ATGG"); // c.1..c.4 = ATGG (Met start codon spillover into a benign next base)
    seq.push_str("CTTTCTTTAG"); // c.5..c.14
    while seq.len() < 200 {
        seq.push('G');
    }
    let len = seq.len() as u64;
    let transcript = Transcript::new(
        "NM_TEST382.1".to_string(),
        Some("TEST382".to_string()),
        Strand::Plus,
        seq,
        Some(1),
        Some(200),
        vec![Exon::new(1, 1, len)],
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

fn normalize_cds(input: &str) -> String {
    let normalizer = Normalizer::with_config(
        cds_provider(),
        NormalizeConfig::default()
            .with_direction(ShuffleDirection::ThreePrime)
            .allow_crossing_boundaries(),
    );
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{}", normalized)
}

#[test]
fn cds_insertion_with_downstream_tract_match_emits_dup() {
    // Spec-canonical form for c.5_6insTTTCTT on a transcript whose
    // c.5..c.14 = "CTTTCTTTAG" is `c.7_12dup`. Mirrors the genomic
    // anchor case on the c. axis.
    assert_eq!(
        normalize_cds("NM_TEST382.1:c.5_6insTTTCTT"),
        "NM_TEST382.1:c.7_12dup",
    );
}

#[test]
fn cds_delins_with_downstream_tract_match_emits_dup() {
    // Delins form of the same edit (shared-affix trim consumes the
    // leading "C") — exercises the canonicalize_delins → recurse path
    // on the c. axis.
    assert_eq!(
        normalize_cds("NM_TEST382.1:c.5delinsCTTTCTT"),
        "NM_TEST382.1:c.7_12dup",
    );
}
