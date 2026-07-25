//! Regression test for issue #340 subgroup (i): 5'-direction shuffle of
//! a single-base insertion in a homopolymer tract picks the 3'-most
//! position instead of the 5'-most. Biocommons cases (from #324):
//!
//!   NM_000051.3:c.14_15insT   bio: c.14dup   ferro: c.15dup   (5prime)
//!   NM_000051.3:c.-4_-3insA   bio: c.-4dup   ferro: c.-3dup   (5prime)
//!   NM_000051.3:c.9171_*1insA bio: c.9171dup ferro: c.9170dup (5prime)
//!   NM_000051.3:c.*4_*5insT   bio: c.*3dup   ferro: c.*4dup   (5prime)
//!
//! Self-contained synthetic transcript pins the pattern without needing
//! the real NM_000051.3 manifest. The (ii) delins → dup and (iii)
//! CDS-start clamp subgroups in the issue body are left for follow-up
//! PRs per the issue author's own scoping note.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Synthetic transcript:
///   tx positions 1-20  (CDS-only, no 5'UTR / 3'UTR — avoids overlap
///                       with the axis-clamp work in #337)
///   c.1..c.20 == "ATGAATTAAATAGAAAAATA"
///                       ↑   ↑↑      — c.5=A, c.6=T, c.7=T, c.8=A
///
/// `c.6_7insT` inserts another T into the 2-T tract at c.6-c.7,
/// growing it to a 3-T tract at c.6, c.7, c.8 (or 5, 6, 7 depending
/// on shift direction). Under 5prime, the canonical dup is c.6dup;
/// under 3prime, c.7dup. ferro currently picks 3prime regardless.
///
/// We use a homopolymer of length 2 specifically so the shift has
/// exactly one position to cover — making the off-by-one isolable.
fn provider_with_t_tract_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    //                       12345678901234567890
    let sequence = String::from("ATGAATTAAATAGAAAAATA");
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(1),
        Some(20),
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

fn normalize_with(direction: ShuffleDirection, cross: bool, input: &str) -> String {
    let mut config = NormalizeConfig::default().with_direction(direction);
    if cross {
        config = config.allow_crossing_boundaries();
    }
    let normalizer = Normalizer::with_config(provider_with_t_tract_transcript(), config);
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{}", normalized)
}

#[test]
fn five_prime_insertion_in_homopolymer_picks_leftmost_dup() {
    // c.6_7insT into the c.6-c.7 T-tract (c.5=A, c.6=T, c.7=T, c.8=A).
    // 5prime-most dup is c.6dup. ferro pre-fix emits c.7dup.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, false, "NM_TEST.1:c.6_7insT"),
        "NM_TEST.1:c.6dup",
    );
}

#[test]
fn three_prime_insertion_in_homopolymer_picks_rightmost_dup() {
    // Sanity check that the 3prime direction still picks c.7dup.
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, false, "NM_TEST.1:c.6_7insT"),
        "NM_TEST.1:c.7dup",
    );
}

#[test]
fn five_prime_cross_insertion_in_homopolymer_picks_leftmost_dup() {
    // cross_boundaries=true must NOT change the homopolymer dup position;
    // it only affects exon-intron crossing. Single-exon transcript here
    // has no introns, but the test pins behavior under both cross modes.
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, true, "NM_TEST.1:c.6_7insT"),
        "NM_TEST.1:c.6dup",
    );
}

/// Synthetic transcript with a longer (4-base) homopolymer tract at
/// positions c.10-c.13 — used to verify the 5'-most slot is picked
/// across multiple candidate positions, not just "the other endpoint
/// of a 2-slot tract" (Opus + Sonnet review feedback).
fn provider_with_long_t_tract_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    //                              1234567890123456789012345
    let sequence = String::from("ATGAGAGAGTTTTGAAAAAAAACCC");
    //                                       ^^^^ c.10-c.13 = TTTT
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_LONG.1".to_string(),
        Some("LONG".to_string()),
        Strand::Plus,
        sequence,
        Some(1),
        Some(25),
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

#[test]
fn five_prime_insertion_in_4base_homopolymer_picks_leftmost_dup() {
    // c.10-c.13 = T,T,T,T. Inserting another T anywhere inside the
    // tract — the 5'-most canonical dup is c.10dup regardless of how
    // the insertion was originally phrased. Pins that the fix
    // generalizes beyond 2-base tracts.
    let normalizer = Normalizer::with_config(
        provider_with_long_t_tract_transcript(),
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let variant = parse_hgvs("NM_LONG.1:c.11_12insT").expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    assert_eq!(format!("{}", normalized), "NM_LONG.1:c.10dup");
}

#[test]
fn three_prime_insertion_in_4base_homopolymer_picks_rightmost_dup() {
    // Same fixture, 3prime direction: rightmost dup is c.13dup.
    let normalizer = Normalizer::with_config(
        provider_with_long_t_tract_transcript(),
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    );
    let variant = parse_hgvs("NM_LONG.1:c.11_12insT").expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    assert_eq!(format!("{}", normalized), "NM_LONG.1:c.13dup");
}

/// Follow-up #357: audit `insertion_to_repeat` for direction-awareness.
/// The output of repeat notation (`A[N]`) covers the full tract, so it
/// should be invariant under `ShuffleDirection`. The audit uses an
/// `n.` (non-coding) context so the codon-frame gate doesn't suppress
/// repeat notation on 1-base units.
fn provider_with_aaa_tract_noncoding() -> MockProvider {
    let mut provider = MockProvider::new();
    //                          123456789012345
    let sequence = "CGTAAAGAGAGGCAA".to_string();
    //                  ^^^ n.4-n.6 = AAA (3-copy A tract)
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NR_REPAUDIT.1".to_string(),
        Some("REPAUDIT".to_string()),
        Strand::Plus,
        sequence,
        None,
        None,
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

/// Follow-up #356: a single-position `delinsXXX` where the post-trim
/// inserted bases duplicate a contiguous span of nearby reference
/// bases should canonicalize to `dup` over that span. Biocommons case
/// from #324: `NC_000001.10:g.1647893delinsCTTTCTT` → `g.1647895_1647900dup`.
///
/// Synthetic reproducer: 13-base genomic window where position 4
/// holds `C` and positions 5-10 hold `TTTCTT`. The delins consumes
/// position 4 (the `C`) and inserts `CTTTCTT`; trim leaves
/// `insTTTCTT` at the 4_5 boundary, which matches ref[5..10] verbatim
/// and should re-classify as `5_10dup`.
#[test]
fn delins_to_dup_canonicalizes_multibase_dup_on_genome() {
    // Build a 250-byte genomic window so the default normalize window
    // size (100 bp on each side) is satisfied. Place the variant at
    // position 150 with the dup-target tract at positions 151-156.
    //
    //   pos:  ... 149 150 151 152 153 154 155 156 157 ...
    //   ref:  ...  A   C   T   T   T   C   T   T   G  ...
    let mut filler = "A".repeat(149);
    filler.push_str("CTTTCTTGGGAAA"); // positions 150..162
    filler.push_str(&"A".repeat(250 - filler.len()));
    let mut provider = ferro_hgvs::MockProvider::new();
    provider.add_genomic_sequence("NC_000001.10", &filler);
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    );
    let variant = parse_hgvs("NC_000001.10:g.150delinsCTTTCTT").expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    assert_eq!(
        format!("{}", normalized),
        "NC_000001.10:g.151_156dup",
        "delinsCTTTCTT at pos 150 (where ref=C and ref[151..156]=TTTCTT) must \
         canonicalize to 151_156dup",
    );
}

/// Normalize `input` against `provider` in `dir`, rendered.
fn norm_with(provider: MockProvider, input: &str, dir: ShuffleDirection) -> String {
    let normalizer =
        Normalizer::with_config(provider, NormalizeConfig::default().with_direction(dir));
    let variant = parse_hgvs(input).expect("parse");
    format!("{}", normalizer.normalize(&variant).expect("normalize"))
}

/// #357: `insertion_to_repeat` on an **unambiguously-phased** tract is
/// direction-invariant.
///
/// 2-copy insertion `AA` extends the n.4-n.6 `AAA` tract to 5 copies. The
/// canonical repeat form covers the full tract, and for a homopolymer the
/// 3'-most and 5'-most phases are the same window, so both directions must
/// render identically.
///
/// Scope, stated precisely because the original name (`..._is_direction_invariant_...`)
/// over-promised: this pins invariance ONLY where the tract's phase is
/// unambiguous. It is not a general direction-invariance guarantee — a tract
/// whose flanking base continues the unit's rotation phases to the
/// direction-appropriate end, which the companion test below pins. Before that
/// companion existed, a homopolymer was the only case under test, so this test
/// would have kept passing even if multi-base-unit direction handling broke
/// entirely.
#[test]
fn insertion_to_repeat_is_direction_invariant_for_unambiguous_tract_on_n_axis() {
    let input = "NR_REPAUDIT.1:n.6_7insAA";
    let five_prime = norm_with(
        provider_with_aaa_tract_noncoding(),
        input,
        ShuffleDirection::FivePrime,
    );
    let three_prime = norm_with(
        provider_with_aaa_tract_noncoding(),
        input,
        ShuffleDirection::ThreePrime,
    );
    // Assert the canonical value on each direction independently, not merely that
    // the two agree: two identically-wrong renderings would satisfy a mutual-
    // equality check. The window names the REFERENCE tract (n.4-n.6, 3 copies of
    // `A`), and the count is the post-edit total (3 ref + 2 inserted = 5).
    const CANONICAL: &str = "NR_REPAUDIT.1:n.4_6A[5]";
    assert_eq!(three_prime, CANONICAL, "3' canonical form");
    assert_eq!(five_prime, CANONICAL, "5' canonical form");
    assert_eq!(
        five_prime, three_prime,
        "a homopolymer tract has one phase, so insertion_to_repeat must render \
         identically in both directions: 5prime={five_prime} 3prime={three_prime}",
    );
}

/// An `n.` transcript whose n.3-n.11 is the odd alternating run `GAGAGAGAG`.
///
/// Odd length is what makes the phase ambiguous: read as `GA` the run starts at
/// n.3, read as `AG` it starts at n.4, and each flanking base continues the
/// other reading's rotation. That is the shape a homopolymer fixture cannot
/// express, and the one where direction actually decides the window.
fn provider_with_ambiguous_ag_tract_noncoding() -> MockProvider {
    let mut provider = MockProvider::new();
    //                       1234567890123
    let sequence = "CCGAGAGAGAGCC".to_string();
    //                ^^^^^^^^^ n.3-n.11 = GAGAGAGAG
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NR_REPAUDIT.1".to_string(),
        Some("REPAUDIT".to_string()),
        Strand::Plus,
        sequence,
        None,
        None,
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

/// #357 companion: on an **ambiguously-phased** tract the repeat form phases to
/// the direction-canonical end, so the two directions deliberately differ.
///
/// This is the case the homopolymer test above cannot reach. Inserting `AGAG`
/// against the odd `GAGAGAGAG` run yields one haplotype with two equally-valid
/// spellings; 3' names the 3'-most phase (`AG[6]` at n.4) and 5' the 5'-most
/// (`GA[6]` at n.3). Asserting the *difference* is what keeps a future change
/// from silently collapsing 5' back onto the 3' phase — the defect class fixed
/// in the repeat/dup direction-awareness work.
#[test]
fn insertion_to_repeat_phases_to_direction_canonical_end_for_ambiguous_tract() {
    let input = "NR_REPAUDIT.1:n.11_12insAGAG";
    let three_prime = norm_with(
        provider_with_ambiguous_ag_tract_noncoding(),
        input,
        ShuffleDirection::ThreePrime,
    );
    let five_prime = norm_with(
        provider_with_ambiguous_ag_tract_noncoding(),
        input,
        ShuffleDirection::FivePrime,
    );
    assert_eq!(three_prime, "NR_REPAUDIT.1:n.4_11AG[6]", "3'-most phase");
    assert_eq!(five_prime, "NR_REPAUDIT.1:n.3_10GA[6]", "5'-most phase");
    assert_ne!(
        three_prime, five_prime,
        "an ambiguous tract must phase to the direction-canonical end",
    );
}

/// Both directions stay confluent on the ambiguous tract: other spellings of the
/// same haplotype (a 5'-side insertion, and the equivalent `dup`) reach the same
/// direction-canonical form as the 3'-side insertion above.
#[test]
fn ambiguous_tract_spellings_are_confluent_within_each_direction() {
    for input in [
        "NR_REPAUDIT.1:n.2_3insGAGA",
        "NR_REPAUDIT.1:n.3_6dup",
        "NR_REPAUDIT.1:n.11_12insAGAG",
    ] {
        assert_eq!(
            norm_with(
                provider_with_ambiguous_ag_tract_noncoding(),
                input,
                ShuffleDirection::ThreePrime
            ),
            "NR_REPAUDIT.1:n.4_11AG[6]",
            "3' confluence for {input}",
        );
        assert_eq!(
            norm_with(
                provider_with_ambiguous_ag_tract_noncoding(),
                input,
                ShuffleDirection::FivePrime
            ),
            "NR_REPAUDIT.1:n.3_10GA[6]",
            "5' confluence for {input}",
        );
    }
}
