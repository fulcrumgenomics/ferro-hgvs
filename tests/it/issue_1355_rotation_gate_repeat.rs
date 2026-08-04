//! Regression for #1355, at the `Normalizer` level.
//!
//! `insertion_to_repeat` scored every cyclic rotation of the inserted unit by the
//! size of the reference tandem tract it abuts, kept the single
//! highest-`ref_count` winner, and only then applied the #882 preservation gate.
//! So a rotation abutting a longer but *out-of-phase* tract discarded a valid
//! in-phase repeat and took the whole result down with it — the rotation that
//! lost on score was never reconsidered.
//!
//! `NM_TEST1355.1:c.2_3insTCCTCC` is the repro. Rotation `TCC` abuts a one-copy
//! in-phase tract at c.-1_2 (the correct answer); rotation `CCT` abuts a
//! four-copy tract further 3' that is out of phase with the insertion, and wins
//! on count. Before the fix the 3' direction emitted the un-promoted
//! `c.2_3insTCCTCC` while the 5' direction emitted `c.1delinsCCTCCTC`, so the two
//! directions also *disagreed* about an edit that is direction-independent. Both
//! now reach `c.-1_2TCC[3]`.
//!
//! This is byte-for-byte the defect #1349 corrected in `insertion_to_duplication`
//! one function over, so this file is deliberately the sibling of
//! `issue_1349_rotation_fallback_dup.rs` and shares its synthetic reference.
//!
//! The unit tests in `rules.rs` call `insertion_to_repeat` directly. These go
//! through `Normalizer::normalize` — and therefore `normalize_core_checked`,
//! where the `FERRO_ASSERT_IDEMPOTENT=1` / `FERRO_ASSERT_REPARSE=1` oracles sit —
//! so the emitted repeat is judged as a whole-variant answer rather than as a
//! tuple of five fields.

use ferro_hgvs::normalize::rules::insertion_to_repeat;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// 5'UTR `GTACGT`, then a CDS opening `CCCCTCCTCCTCCTCCCGGGG` (the #1349
/// sequence, so the sibling defects share a reference), then a CDS-internal `TA`
/// tract for the codon-gate case below.
///
/// ```text
/// axis: c.-6 … c.-1 | c.1 c.2 c.3 … c.21 c.22 c.23 … c.30
/// base:  G  …   T   |  C   C   C  …   G    G    T  …   A
/// ```
///
/// The two tracts that compete for `c.2_3insTCCTCC`:
///   - `TCC` at c.-1_2 (one copy) — abuts the insertion point, in phase;
///   - `CCT` over c.3_14 (four copies, 0-based `[8, 20)`) — out of phase with the
///     insertion, and it outscores `TCC` four to one.
const UTR5: &str = "GTACGT";
const CDS_OPEN: &str = "CCCCTCCTCCTCCTCCCGGGG";
const CDS_TA_TRACT: &str = "GTATATATA";
const UTR3: &str = "AAAAAAAAAA";

/// 1-based position of `c.1`: `UTR5` is six bases, so `c.1` is the seventh.
const CDS_START: u64 = 7;

fn transcript_sequence() -> String {
    format!("{UTR5}{CDS_OPEN}{CDS_TA_TRACT}{UTR3}")
}

/// 1-based inclusive CDS end — everything but the 3'UTR. The CDS is 30 bases
/// (ten codons), so the codon gate is not tripped by a ragged CDS length.
fn cds_end() -> u64 {
    (UTR5.len() + CDS_OPEN.len() + CDS_TA_TRACT.len()) as u64
}

fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = transcript_sequence();
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_TEST1355.1".to_string(),
        Some("TEST1355".to_string()),
        Strand::Plus,
        sequence,
        Some(CDS_START),
        Some(cds_end()),
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
    format!("{normalized}")
}

const REPRO: &str = "NM_TEST1355.1:c.2_3insTCCTCC";
const EXPECTED: &str = "NM_TEST1355.1:c.-1_2TCC[3]";

/// The in-phase `TCC` tract must be found even though the out-of-phase `CCT`
/// rotation outscores it. Before the fix this direction emitted the input back
/// unchanged, having discarded the only valid candidate.
#[test]
fn three_prime_insertion_promotes_to_the_in_phase_repeat() {
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, REPRO),
        EXPECTED
    );
}

/// Same answer in the other direction. `find_tandem_extent` scores the rotations
/// identically either way, so the loser-takes-all bug was direction-blind — but
/// its *symptoms* were not: 5' fell through to the CDS-start clamp and emitted
/// `c.1delinsCCTCCTC`, so the two directions disagreed about an edit whose
/// canonical form is direction-independent.
#[test]
fn five_prime_insertion_promotes_to_the_in_phase_repeat() {
    assert_eq!(normalize_with(ShuffleDirection::FivePrime, REPRO), EXPECTED);
}

/// The directions agree, pinned as its own assertion so a future change that
/// moves the canonical answer still has to move it in both directions together.
#[test]
fn both_directions_agree() {
    assert_eq!(
        normalize_with(ShuffleDirection::ThreePrime, REPRO),
        normalize_with(ShuffleDirection::FivePrime, REPRO),
    );
}

/// `norm(norm(x)) == norm(x)` in both directions, which is what the
/// `FERRO_ASSERT_IDEMPOTENT=1` oracle asserts for every normalize call in the
/// suite; pinned here too so the guard holds without the env var set.
#[test]
fn normalization_is_a_fixed_point_in_both_directions() {
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        let once = normalize_with(direction, REPRO);
        let twice = normalize_with(direction, &once);
        assert_eq!(
            once, twice,
            "{direction:?}: normalization must be idempotent"
        );
    }
}

/// The emitted repeat must describe the *same haplotype* as the input insertion.
///
/// Checked by splicing both edits into the reference by hand and comparing bytes,
/// deliberately **not** through `tandem_edit_preserves_insertion` — that is the
/// production predicate the fix relies on, so using it here would only assert
/// that the code agrees with itself.
///
/// Coordinates, both derived from the pinned strings above:
///   - input `c.2_3ins`: `c.2` is 1-based position 8, so the insertion splices in
///     at 0-based index 8;
///   - emitted `c.-1_2TCC[3]`: `c.-1` is 1-based position 6 and `c.2` is 8, so
///     the repeat replaces 0-based `[5, 8)` with three copies of `TCC`.
#[test]
fn the_emitted_repeat_reconstructs_the_inserted_haplotype() {
    // Tie the hardcoded coordinates below to what normalization actually emits.
    // Without this the test would assert only a property of the reference and
    // would stay green if the emission ever moved.
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_eq!(
            normalize_with(direction, REPRO),
            EXPECTED,
            "{direction:?}: the coordinates below are derived from this emission"
        );
    }

    let reference = transcript_sequence();

    let mut from_input = String::new();
    from_input.push_str(&reference[..8]);
    from_input.push_str("TCCTCC");
    from_input.push_str(&reference[8..]);

    let mut from_repeat = String::new();
    from_repeat.push_str(&reference[..5]);
    from_repeat.push_str(&"TCC".repeat(3));
    from_repeat.push_str(&reference[8..]);

    assert_eq!(
        from_input, from_repeat,
        "`{EXPECTED}` must spell the same haplotype as `{REPRO}`"
    );
}

/// A codon-gated tract declines the repeat and the pipeline still lands a valid
/// answer.
///
/// In a `c.` context a CDS-resident tract with a non-codon-length unit may not be
/// written as repeat notation (#1210), which inside the rotation loop arrives as
/// `RepeatNormResult::Insertion`. `c.24_25insTATA` sits in the `TA` tract at
/// c.23_30, wholly inside the CDS, so `insertion_to_repeat` declines — proven to be
/// the *gate*, not a failed tract search, by the sibling test below. The pipeline
/// keeps going and promotes the edit to a duplication; `dup` carries no codon
/// restriction, since the spec's exception forbids only `[N]`. Which copy is named
/// is the direction rule's call, so the two directions differ here by design.
///
/// **What this does not pin.** #1355 turned that gate's `return None` into a
/// `continue`, and this case cannot tell the two apart: *every* rotation here
/// shares the one CDS-resident tract, so both spellings decline the function and
/// this test passes against `main` as well. The two are observably different only
/// when the highest-scoring rotation is codon-gated while a lower-scoring one is
/// not — which needs two tracts straddling `cds_end` so they get different
/// residency answers, and is not constructed here. The value below is still worth
/// pinning: it is the fixed point the codon-gated path has to keep reaching.
#[test]
fn a_codon_gated_tract_still_reaches_a_fixed_point() {
    for (direction, expected) in [
        (ShuffleDirection::ThreePrime, "NM_TEST1355.1:c.27_30dup"),
        (ShuffleDirection::FivePrime, "NM_TEST1355.1:c.23_26dup"),
    ] {
        let once = normalize_with(direction, "NM_TEST1355.1:c.24_25insTATA");
        assert_eq!(
            once, expected,
            "{direction:?}: a codon-gated repeat must fall through to the dup"
        );
        let twice = normalize_with(direction, &once);
        assert_eq!(once, twice, "{direction:?}: and must be a fixed point");
    }
}

/// It is the *codon gate* that declines the tract above, not a failure to find it.
///
/// Called directly because that is the only way to hold the tract search fixed
/// while flipping the gate: with CDS context every rotation is rejected and the
/// function declines, and with the gate lifted the very same tract is emitted as
/// `TA[6]` over c.23_30 (1-based 29_36).
#[test]
fn the_codon_gate_is_what_declines_the_tract() {
    let reference = transcript_sequence();
    let reference = reference.as_bytes();
    // `c.24` is 1-based position 30 (`c.n` is position `n + 6`), and
    // `insertion_to_repeat` takes the 0-based index of the base to the *left* of
    // the insertion point, so `c.24_25ins` is index 29.
    let insertion_point = CDS_START + 24 - 1 - 1;

    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert!(
            insertion_to_repeat(
                reference,
                insertion_point,
                b"TATA",
                true,
                Some((CDS_START, cds_end())),
                direction,
            )
            .is_none(),
            "{direction:?}: a CDS-resident `TA` unit must be codon-gated"
        );

        let (_, count, start, end, unit) =
            insertion_to_repeat(reference, insertion_point, b"TATA", false, None, direction)
                .unwrap_or_else(|| {
                    panic!("{direction:?}: the same tract must be found once the gate is lifted")
                });
        assert_eq!(
            (unit.as_slice(), count, start, end),
            (b"TA".as_slice(), 6, 29, 36),
            "{direction:?}: expected `TA[6]` over the 1-based tract 29_36"
        );
    }
}
