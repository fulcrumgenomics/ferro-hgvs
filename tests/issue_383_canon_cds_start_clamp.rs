//! Regression test for issue #383: 5'-direction canonicalisation
//! rewrite on a c.-axis Insertion/Delins must clamp at the CDS-start
//! boundary instead of silently moving the variant strictly into 5'UTR.
//!
//! Companion to PR #343 (Pattern C) which closed the shuffle-path
//! CDS-start clamp. The remaining divergences live in the
//! canonicalisation-rewrite path — ferro's pre-shift rewrite
//! (`canonicalize_delins` → `DelinsCanonical::Insertion` recursion in
//! `mod.rs:~3153`) already lands the variant in 5'UTR before the
//! shift-path clamp sees it.
//!
//! Biocommons cases (5prime+cross, from `tests/fixtures/biocommons-
//! normalize/cases.json` and the current `/tmp/ferro-xfail/biocommons-
//! normalized.tsv` snapshot):
//!
//!   NM_212556.2:c.1_2insCA    bio: c.1delinsACA  ferro: c.-2_-1dup
//!   NM_212556.2:c.1delinsCA   bio: c.1delinsCA   ferro: c.-1_1insC
//!   NM_212556.2:c.2_3insCAT   bio: c.1delinsATCA ferro: c.-1_1insCAT
//!
//! Unifying spec rule (from `tests/fixtures/biocommons-normalize/
//! failure-patterns.md` and HGVS §general "3'-rule applies to ALL
//! descriptions" combined with the spec's per-axis coordinate
//! treatment):
//!
//!   A 5'-direction canonicalisation rewrite on a c.-axis variant may
//!   not land *strictly* inside 5'UTR. When the rewrite cannot stay
//!   anchored at c.1 or later, emit as `c.1delins<…>` (absorbing the
//!   boundary base into the alt) for Insertion inputs, or suppress the
//!   rewrite entirely for Delins inputs whose canonicalisation would
//!   move them past the boundary.
//!
//! For NM_212556.2 the real reference at c.1..c.3 is "ATG…". The
//! synthetic transcript here reproduces those bases so the test pins
//! exactly the biocommons expectations.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Synthetic transcript whose CDS-start is at c.1 and whose first three
/// CDS bases match NM_212556.2's known c.1..c.3 = "ATG…" context plus
/// a tail that puts `c.2_3insCAT` exactly on the spec-clamp path:
///
///   axis: c.-3 c.-2 c.-1 | c.1 c.2 c.3 c.4 c.5 ...
///   base:  G    G    G   | A   T   G   C   ...
///
/// (HGVS skips c.0 — c.-1 is immediately followed by c.1.)
///
/// CDS spans c.1..c.200 (synthetic — the spec rule fires only on the
/// CDS-start side regardless of CDS-end). The 5'UTR is short (c.-3..c.-1
/// = 3 bases of G filler) and the 3'UTR is empty.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    // 3 bytes 5'UTR + 200 bytes CDS. Total tx length = 203.
    // Tx bytes 0..3 = 5'UTR (positions c.-3..c.-1).
    // Tx byte 3 = c.1 (CDS start).
    // Tx byte 4 = c.2, byte 5 = c.3, …
    let mut seq = String::from("GGG"); // 5'UTR
    seq.push_str("ATGC"); // c.1..c.4 = A, T, G, C
    while seq.len() < 203 {
        seq.push('G');
    }
    let len = seq.len() as u64;
    // cds_start / cds_end in transcript coords (1-based byte positions).
    // c.1 is at byte index 3 (0-based) = byte position 4 (1-based).
    // c.200 is at byte index 202 (0-based) = byte position 203 (1-based).
    let transcript = Transcript::new(
        "NM_TEST383.1".to_string(),
        Some("TEST383".to_string()),
        Strand::Plus,
        seq,
        Some(4),
        Some(203),
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

fn normalize_with(direction: ShuffleDirection, input: &str) -> String {
    let normalizer = Normalizer::with_config(
        provider(),
        NormalizeConfig::default()
            .with_direction(direction)
            .allow_crossing_boundaries(),
    );
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{}", normalized)
}

#[test]
fn five_prime_insertion_at_cds_start_clamps_to_delins_at_c1() {
    // c.1_2insCA on a transcript whose c.1=A: biocommons emits
    // c.1delinsACA (the 5'-shift would land the ins at c.-1_1,
    // which is strictly inside 5'UTR; clamp emits delins-at-c.1
    // absorbing the boundary base "A" into the alt).
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST383.1:c.1_2insCA"),
        "NM_TEST383.1:c.1delinsACA",
    );
}

#[test]
fn five_prime_delins_at_cds_start_suppresses_rewrite() {
    // c.1delinsCA: shared-affix trim would produce an Insertion of "C"
    // landing at c.-1_1 (strictly inside 5'UTR); the clamp suppresses
    // the rewrite and keeps the original delins form.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST383.1:c.1delinsCA"),
        "NM_TEST383.1:c.1delinsCA",
    );
}

#[test]
fn three_prime_delins_at_cds_start_suppresses_rewrite() {
    // Same suppression must apply in the 3prime direction since
    // biocommons emits c.1delinsCA for both directions (the
    // canonicalisation-rewrite path is direction-agnostic for this
    // shape).
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, "NM_TEST383.1:c.1delinsCA"),
        "NM_TEST383.1:c.1delinsCA",
    );
}

#[test]
fn five_prime_insertion_two_bases_in_clamps_to_delins_at_c1() {
    // c.2_3insCAT: 5'-shift walks left two positions (alt rotates
    // CAT → TCA → ATC), the second step would cross into c.0_1 /
    // c.-1_1 (5'UTR); clamp emits delins-at-c.1 absorbing ref[c.1]="A"
    // into the alt: c.1delinsATCA (= "ATC" rotated alt + "A" boundary
    // base).
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST383.1:c.2_3insCAT"),
        "NM_TEST383.1:c.1delinsATCA",
    );
}

#[test]
fn five_prime_cds_interior_rewrite_unaffected_by_clamp() {
    // Negative: a c.-axis insertion whose canonical rewrite stays
    // entirely inside CDS proper must NOT be clamped. Insert "C"
    // between c.10 and c.11 — far from any boundary; the C does not
    // match the surrounding G-tract so no left-shift occurs and the
    // input is preserved verbatim. Pin the exact output so a stray
    // clamp landing at a positive anchor (e.g. `c.1delins…`) would
    // fail the test, not just a sign flip into 5'UTR.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST383.1:c.10_11insC"),
        "NM_TEST383.1:c.10_11insC",
    );
}

#[test]
fn five_prime_utr_interior_rewrite_unaffected_by_clamp() {
    // Negative: a 5'UTR-interior insertion whose canonical rewrite
    // stays inside 5'UTR must NOT be clamped — the CDS-start clamp is
    // about CDS↔5'UTR crossings, not about 5'UTR-resident variants.
    //
    // c.-3_-2insG on the 5'UTR (ref = "GGG") stays in 5'UTR and 5'-
    // shifts within the G-tract to the canonical `c.-3dup`. Pin the
    // exact output so a misplaced clamp that snaps to `c.1…` or
    // produces any other axis token fails the test.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, "NM_TEST383.1:c.-3_-2insG"),
        "NM_TEST383.1:c.-3dup",
    );
}

/// Provider for the long-shift homopolymer regression: c.1..c.20 = 20
/// A's, then G filler, so a CDS-interior insertion of an A-homopolymer
/// can 5'-shift through many CDS positions and land strictly inside
/// 5'UTR before the clamp gets a chance to fire.
///
///   axis: c.-3 c.-2 c.-1 | c.1 c.2 ... c.20 c.21 ...
///   base:  G    G    G   | A   A   ...  A    G   ...
fn provider_homopolymer() -> MockProvider {
    let mut provider = MockProvider::new();
    let mut seq = String::from("GGG"); // 5'UTR (c.-3..c.-1)
    for _ in 0..20 {
        seq.push('A'); // CDS c.1..c.20
    }
    while seq.len() < 203 {
        seq.push('G');
    }
    let len = seq.len() as u64;
    let transcript = Transcript::new(
        "NM_TEST383HP.1".to_string(),
        Some("TEST383HP".to_string()),
        Strand::Plus,
        seq,
        Some(4),
        Some(203),
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

fn normalize_with_provider(
    provider: MockProvider,
    direction: ShuffleDirection,
    input: &str,
) -> String {
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default()
            .with_direction(direction)
            .allow_crossing_boundaries(),
    );
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{}", normalized)
}

#[test]
fn five_prime_insertion_long_homopolymer_shift_clamps_at_cds_start() {
    // CodeRabbit-flagged regression: when a CDS-interior literal
    // insertion 5'-shifts across a long enough homopolymer that
    // `x = tx_start - cds_start + 1 > alt.len() + 1`, the original
    // clamp branch silently no-op'd and the result landed at
    // `c.-1_1insAA` — strictly inside 5'UTR even though the input was
    // CDS-interior.
    //
    // Spec-canonical clamp for the homopolymer case is a delins anchored
    // at c.1 that absorbs the boundary positions: `c.1_3delinsAAAAA`
    // (= delete c.1..c.3 = "AAA", insert 5 A's). Applying this to the
    // ref `GGG | AAAAAAAAAAAAAAAAAAAA | G...` yields `GGG + AAAAA +
    // ref[c.4..]` = `GGG + 22*A + G...`, the same total as the input
    // `c.5_6insAA` (insert 2 A's into the 20-A tract).
    //
    // The exact form here is the *implementation's* spec-canonical
    // emission; assert it precisely so any future regression to
    // c.-N_… or to a different delete window is caught immediately.
    assert_eq!(
        normalize_with_provider(
            provider_homopolymer(),
            ShuffleDirection::FivePrime,
            "NM_TEST383HP.1:c.5_6insAA",
        ),
        "NM_TEST383HP.1:c.1_3delinsAAAAA",
    );
}

#[test]
fn five_prime_insertion_far_interior_homopolymer_shift_clamps() {
    // Same long-shift homopolymer pattern but starting much deeper
    // inside the CDS, with a longer alt. tx_start - cds_start + 1 = 20
    // and alt.len() = 3, so x = 20 >> alt.len() + 1 = 4 — the buggy
    // branch was silently no-op'ing on inputs of this shape.
    //
    // Clamp form: delete c.1..c.{x-L} = c.1..c.17 (17 A's), insert
    // x A's = 20 A's. Equivalent to extending the 20-A tract by 3
    // bases (the alt length).
    assert_eq!(
        normalize_with_provider(
            provider_homopolymer(),
            ShuffleDirection::FivePrime,
            "NM_TEST383HP.1:c.20_21insAAA",
        ),
        "NM_TEST383HP.1:c.1_17delinsAAAAAAAAAAAAAAAAAAAA",
    );
}
