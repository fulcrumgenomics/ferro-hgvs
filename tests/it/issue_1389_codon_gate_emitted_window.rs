//! Regression for #1389, at the `Normalizer` level.
//!
//! `insertion_to_repeat` answered the #1210 codon gate about the tract
//! `find_tandem_extent` *discovered*, but `normalize_repeat` then phase-aligns
//! that tract before emitting, and the two windows can differ. A repeat whose
//! **emitted** window is wholly CDS-resident could therefore be written with a
//! non-codon-length unit, which `repeated.md`'s codon rule forbids.
//!
//! The reference below puts a `TG` tract across the CDS start, so the 3'-aligned
//! window lands inside the CDS while the raw tract does not:
//!
//! ```text
//! axis:  c.-5 c.-4 c.-3 c.-2 c.-1 | c.1 c.2 c.3 c.4 c.5 c.6 …
//! base:   A    A    A    A    T   |  G   T   G   T   G   T  …
//!                            └────── raw TG tract starts here, in the 5'UTR
//!                                     └──── 3'-aligned GT window starts here
//! ```
//!
//! On `NM_T1389.1:c.3_4insTGTG` under 3', the discovered tract starts at `c.-1`
//! and so reads as non-resident, switching the gate off — but the emitted window
//! is `c.1_6`, wholly inside the CDS, and it was written `c.1_6GT[5]` with a
//! two-base unit. Measured against this file's own fixture:
//!
//! | direction | before | after |
//! |---|---|---|
//! | 3' | `c.1_6GT[5]` (invalid: 2-base unit, wholly CDS-resident) | `c.3_6dup` |
//! | 5' | `c.-1_5TG[5]` | `c.-1_5TG[5]` (unchanged) |
//!
//! The 5' answer must not move: its window genuinely includes `c.-1`, so gating
//! the codon rule off there is correct. That the two directions still emit
//! different forms is not a defect — they align to different windows, and the
//! residency answer legitimately differs between them. What was wrong was
//! computing the verdict from a window neither of them emitted.
//!
//! Pre-existing rather than introduced by #1355: the same 3' output appears on
//! `main` and on #1368's head. #1368 moved the gate to run per rotation but did
//! not change *which window* it asks about.
//!
//! The unit tests in `rules.rs` call `insertion_to_repeat` directly. These go
//! through `Normalizer::normalize`, and therefore `normalize_core_checked` where
//! the `FERRO_ASSERT_IDEMPOTENT` / `REPARSE` / `IN_BOUNDS` oracles sit, so the
//! emitted repeat is judged as a whole-variant answer rather than as a tuple.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// 5'UTR ending in `T`, so the `TG` tract begins one base before `c.1`. That
/// single-base straddle is the whole point of the fixture: it is what makes the
/// raw tract non-resident while the 3'-aligned window is resident.
const UTR5: &str = "AAAAT";
/// CDS opening on `GTGTGT` — the rest of the straddling tract — then filler, a
/// four-copy `CAG` tract for the codon-length control at the bottom of this file,
/// and enough tail that the control's expanded six-copy tract still ends inside
/// the CDS. Thirty-six bases, a whole number of codons, so a ragged CDS length
/// cannot be mistaken for the cause of a codon-gate decision.
const CDS: &str = "GTGTGTCCCGGGCAGCAGCAGCAGTTTTTTTTTTTT";
const UTR3: &str = "AAAAAAAAAA";

/// 1-based position of `c.1`: `UTR5` is five bases, so `c.1` is the sixth.
const CDS_START: u64 = 6;

fn cds_end() -> u64 {
    (UTR5.len() + CDS.len()) as u64
}

fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = format!("{UTR5}{CDS}{UTR3}");
    let len = sequence.len() as u64;
    provider.add_transcript(Transcript::new(
        "NM_T1389.1".to_string(),
        Some("T1389".to_string()),
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
    ));
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

const REPRO: &str = "NM_T1389.1:c.3_4insTGTG";

/// The defect. The 3'-aligned window `c.1_6` is wholly CDS-resident, so a
/// two-base unit may not be spelled as repeat notation there; the edit falls
/// through to the duplication, which carries no codon constraint.
///
/// Asserted as an inequality against the old output as well as an equality on the
/// new one, because the equality alone would pass for the wrong reason if some
/// later change moved the canonical answer to a third form.
#[test]
fn three_prime_gates_a_repeat_whose_emitted_window_is_cds_resident() {
    let got = normalize_with(ShuffleDirection::ThreePrime, REPRO);
    assert_ne!(
        got, "NM_T1389.1:c.1_6GT[5]",
        "c.1_6 is wholly inside the CDS, so a two-base repeat unit is forbidden \
         there (repeated.md); the codon gate was answered about the raw tract"
    );
    assert_eq!(got, "NM_T1389.1:c.3_6dup");
}

/// The 5' direction must not move. Its window starts at `c.-1`, in the 5'UTR, so
/// the tract is genuinely not CDS-resident and the codon rule correctly does not
/// apply — a fix that gated this too would be over-reaching.
#[test]
fn five_prime_window_outside_the_cds_still_emits_the_repeat() {
    assert_eq!(
        normalize_with(ShuffleDirection::FivePrime, REPRO),
        "NM_T1389.1:c.-1_5TG[5]"
    );
}

/// Whatever each direction emits must be a fixed point, which is what
/// `FERRO_ASSERT_IDEMPOTENT=1` asserts for every normalize call in the suite —
/// pinned here too so the guard holds without the env var set.
///
/// This is the assertion that would have caught the defect independently of the
/// expected strings above: `c.1_6GT[5]` is a description whose own re-derivation
/// has to agree with it.
#[test]
fn both_directions_reach_a_fixed_point() {
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        let once = normalize_with(direction, REPRO);
        let twice = normalize_with(direction, &once);
        assert_eq!(
            once, twice,
            "{direction:?}: normalization is not idempotent"
        );
    }
}

/// A repeat whose emitted window is wholly CDS-resident *and* whose unit is a
/// whole codon must still be spelled as a repeat. The gate keys on unit length,
/// not on residency alone, so this pins that the fix did not turn CDS-residency
/// itself into a refusal.
#[test]
fn a_codon_length_unit_inside_the_cds_is_still_a_repeat() {
    // The four-copy `CAG` tract sits at c.13_24; adding two copies expands it to
    // c.13_24CAG[6], which ends at c.24 — still inside the 36-base CDS. So this
    // window is wholly CDS-resident, exactly like the gated case above, and
    // differs from it only in that the unit is a whole codon.
    //
    // The unit has to be genuinely three-periodic. `smallest_repeat_unit` reduces
    // `CCCCCC` to the one-base `C`, so a homopolymer would be gated on its unit
    // length and would assert the opposite of what this test claims — a mistake
    // this test was written with before being corrected.
    let got = normalize_with(ShuffleDirection::ThreePrime, "NM_T1389.1:c.15_16insCAGCAG");
    assert_eq!(
        got, "NM_T1389.1:c.13_24CAG[6]",
        "a codon-length unit inside the CDS must still reach repeat notation"
    );
}
