//! Regression test for issue #387: 3'-direction canonicalisation
//! rewrite on a c.-axis Insertion/Delins must clamp at the CDS-end
//! boundary instead of silently moving the variant strictly onto the
//! 3'UTR axis (`c.*<N>` coords).
//!
//! CDS-**end** mirror of #383 (CDS-start clamp). Companion to PR #343
//! which closed the shuffle-path side of the CDS↔UTR clamp; this
//! covers the canonicalisation-rewrite-path side on the 3'-direction.
//!
//! Biocommons case (3prime+cross, 3prime+no-cross — same shape both
//! directions per `tests/fixtures/biocommons-normalize/cases.json`):
//!
//!   NM_212556.2:c.1400_1401insAC   bio: c.1401delinsACA   ferro: c.1401_*1insCA
//!
//! Spec basis (HGVS DNA §general "3'-rule applies to ALL descriptions",
//! combined with the per-axis coordinate treatment): the CDS (`c.<N>`)
//! and 3'UTR (`c.*<N>`) are different sub-axes. A canonicalisation
//! rewrite that silently moves a CDS-interior input onto the 3'UTR
//! axis violates the user's axis intent. Unifying rule:
//!
//!   A 3'-direction canonicalisation rewrite on a c.-axis variant
//!   may not land strictly on the 3'UTR axis. When the rewrite cannot
//!   stay anchored at `c.<cds_end>` or earlier, emit as
//!   `c.<cds_end>delins<…>` (absorbing the boundary base into the
//!   alt) for Insertion inputs, or restore the input form unchanged
//!   for Delins inputs whose canonicalisation would move them past
//!   the boundary.
//!
//! The synthetic transcript here reproduces NM_212556.2's c.<cds_end-1>
//! and c.<cds_end> context (= "AA", consistent with the genotype-
//! equivalence derivation `r[c.1401] = A` from `c.1400_1401insAC` →
//! `c.1401delinsACA`).

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Synthetic CDS-only transcript whose c.4 = c.5 = A (= the
/// boundary-absorbed base), c.cds_end = c.5, c.*1 = G (≠ A so the
/// 3'-shift would terminate after one step past the boundary, just
/// like the real NM_212556.2 input where the unrotated alt won't
/// keep walking into 3'UTR after the first step).
///
///   axis: c.-3 c.-2 c.-1 | c.1 c.2 c.3 c.4 c.5 | c.*1 c.*2 c.*3 ...
///   base:  G    G    G   | A   T   G   A   A   | G    G    G    ...
///                          ^                ^   ^
///                          cds_start        cds_end (c.5)
///                                               first 3'UTR base
///
/// CDS spans c.1..c.5 (synthetic — short for testability). 5'UTR is
/// 3 bases of G filler (so `start_axis` for CDS-interior inputs
/// resolves to `AxisRegion::Cds`); 3'UTR is padded out to the full
/// transcript length.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    // tx bytes 0..3 = c.-3..c.-1 (5'UTR filler)
    // tx bytes 3..8 = c.1..c.5  (CDS: "ATGAA")
    // tx bytes 8..  = c.*1..    (3'UTR: "G…")
    let mut seq = String::from("GGG"); // 5'UTR
    seq.push_str("ATGAA"); // c.1..c.5: A T G A A
    while seq.len() < 203 {
        seq.push('G');
    }
    let len = seq.len() as u64;
    // 1-based tx coords: c.1 at byte 3 → tx 4; c.5 at byte 7 → tx 8.
    let transcript = Transcript::new(
        "NM_TEST387.1".to_string(),
        Some("TEST387".to_string()),
        Strand::Plus,
        seq,
        Some(4),
        Some(8),
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
fn three_prime_insertion_at_cds_end_clamps_to_delins_at_cds_end() {
    // Anchor case for #387, mirror of `NM_212556.2:c.1400_1401insAC` →
    // `c.1401delinsACA`. Input ins "AC" between c.4 and c.5; ref[c.5]
    // = A; 3'-shift would advance the ins to `c.5_*1insCA`, crossing
    // the CDS↔3'UTR axis boundary. Clamp emits `c.5delinsACA` =
    // alt "AC" ++ ref[c.5] = "A".
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, "NM_TEST387.1:c.4_5insAC"),
        "NM_TEST387.1:c.5delinsACA",
    );
}

#[test]
fn no_cross_three_prime_insertion_at_cds_end_clamps_same_way() {
    // The biocommons fixture lists this row for BOTH `3prime+cross`
    // and `3prime+no-cross` directions on `NM_212556.2:c.1400_1401insAC`.
    // The clamp must fire identically regardless of `cross_boundaries`
    // because it is about the per-axis coordinate treatment, not
    // about the exon vs. full-transcript range.
    let normalizer = Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    );
    let variant = parse_hgvs("NM_TEST387.1:c.4_5insAC").expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    assert_eq!(format!("{}", normalized), "NM_TEST387.1:c.5delinsACA");
}

#[test]
fn three_prime_delins_at_cds_end_suppresses_rewrite() {
    // Delins-input mirror of the Delins suppression case in #383's
    // tests. `c.5delinsAC`: shared-affix trim (prefix "A") would
    // leave Insertion of "C" between c.5 and c.6 (= c.*1) — strictly
    // on the 3'UTR axis. The clamp suppresses the rewrite and keeps
    // the original delins form.
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, "NM_TEST387.1:c.5delinsAC"),
        "NM_TEST387.1:c.5delinsAC",
    );
}

#[test]
fn three_prime_cds_interior_rewrite_unaffected_by_clamp() {
    // Negative: a c.-axis insertion whose canonical rewrite stays
    // entirely inside CDS proper must NOT be clamped. `c.1_2insTG` on
    // ref[c.1..c.5] = ATGAA: trim and/or 3'-shift never crosses
    // c.cds_end (= c.5). Whatever ferro emits must stay in positive
    // c. coords (no `c.*`).
    let out = normalize_with(ShuffleDirection::ThreePrime, "NM_TEST387.1:c.1_2insTG");
    assert!(
        !out.contains("c.*"),
        "interior canon was clamped or shifted into 3'UTR: {out}"
    );
    assert!(
        out.starts_with("NM_TEST387.1:c."),
        "interior canon produced unexpected accession/axis: {out}"
    );
}

#[test]
fn three_prime_utr_interior_rewrite_unaffected_by_clamp() {
    // Negative: a 3'UTR-interior insertion whose canonical rewrite
    // stays inside 3'UTR must NOT be clamped — the CDS-end clamp is
    // about CDS↔3'UTR crossings, not about 3'UTR-resident variants.
    let out = normalize_with(ShuffleDirection::ThreePrime, "NM_TEST387.1:c.*1_*2insG");
    assert!(
        out.contains("c.*"),
        "3'UTR-interior canon escaped into CDS via clamp: {out}"
    );
}
