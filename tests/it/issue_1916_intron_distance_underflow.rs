//! Issue #1916 — `IntronPosition::distance_from_donor` / `distance_from_acceptor`
//! must decline an offset they cannot measure, rather than underflow on it.
//!
//! Both functions compute the *cross-boundary* distance as
//! `intron_length - offset.unsigned_abs()` on two `u64`s. When the magnitude
//! exceeds `intron_length` that underflows: `attempt to subtract with overflow`
//! under `overflow-checks` (the `dev`, `test` and `soak` profiles all enable
//! it), and in release a wrap to a value near `u64::MAX` returned as a distance
//! in bases from a splice site.
//!
//! This is the residue of #1767. That PR made this struct's five *predicates*
//! total by reading `offset` through `unsigned_abs` rather than `abs`, and the
//! block comment it left on `src/reference/transcript.rs` names these two
//! functions as the part it did not finish — noting that they "ALREADY return
//! `Option<u64>`, so they have a decline channel and do not use it".
//! `issue_1767_unknown_offset_splice_classifiers::the_intron_position_predicates_are_total`
//! builds the very record that reproduces this (`offset: i64::MIN`,
//! `intron_length: 1_000`) and calls only the five predicates.
//!
//! **Two regimes reach the subtraction, and the answer to both is a refusal.**
//!
//! * `offset.unsigned_abs() > intron_length` — the record claims a position
//!   inside an intron of length `intron_length` sitting further than that from
//!   one of its own boundaries. No position in that intron does. The record is
//!   self-inconsistent and no distance to *either* boundary follows from it.
//! * `offset` is an unknown-offset sentinel (`c.N+?` / `c.N-?`). Per
//!   `is_unknown_offset`'s own contract these "denote a position unbounded in
//!   one direction, from which no distance can be derived at all" — the ruling
//!   #1767 applied to `IntronicConsequence::from_cds_pos` and #1841 to
//!   `IntronicRegion::from_offset` and `EffectPredictor::classify_splice_variant`.
//!
//! **Why not `saturating_sub`.** A clamp would answer `Some(0)`, and `0` is the
//! single worst value available: it reads as *at the splice site*, so every
//! downstream band (`is_canonical_splice_site`, `SpliceSiteType::DonorCanonical`,
//! `IntronicConsequence::SpliceDonorVariant`, HIGH impact) would promote an
//! unmeasurable offset to the most clinically significant answer there is. It
//! would also only fix the arm that panics: the *same*-boundary arm returns
//! `Some(offset.unsigned_abs())` unconditionally and never underflows, so
//! clamping the other one leaves the two functions disagreeing about one
//! position — `distance_from_donor() == None` beside
//! `distance_from_acceptor() == Some(300)`. Hence
//! [`the_two_distances_never_disagree_about_resolvability`], which is the
//! assertion that forbids a one-armed fix.
//!
//! The refusal must not reach a position the library can actually produce, so
//! [`every_derived_intronic_position_still_answers`] walks every genomic base of
//! a real two-intron transcript through `find_intron_at_genomic` and requires an
//! answer at each one. That is the negative control for the whole change: it
//! fails if the guard is drawn even one base too tight.

use ferro_hgvs::hgvs::parser::position::{OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE};
use ferro_hgvs::reference::transcript::{
    Exon, GenomeBuild, IntronBoundary, IntronPosition, ManeStatus, Strand, Transcript,
};

/// Both unknown-offset sentinels, so `+?` and `-?` are covered symmetrically —
/// #1767's own lesson that a guard keyed on `i64::MIN` alone leaves its
/// neighbour passing by accident rather than by decision.
const SENTINELS: [i64; 2] = [OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE];

fn intron_position(boundary: IntronBoundary, offset: i64, intron_length: u64) -> IntronPosition {
    IntronPosition {
        intron_number: 1,
        boundary,
        offset,
        tx_boundary_pos: 100,
        intron_length,
    }
}

/// A two-intron transcript with genomic coordinates on every exon, built
/// programmatically so no fixture is committed. Whichever strand is asked for,
/// the two genomic gaps between the exons are 1101..=1200 (100 bases) and
/// 1301..=1600 (300 bases) — deliberately different lengths, so a guard
/// accidentally keyed on one of them fails on the other.
///
/// The exon layout is strand-specific because `compute_introns` derives an
/// intron's genomic span differently per strand: on `Minus` it reads
/// `downstream.genomic_end + 1 ..= upstream.genomic_start - 1`, so transcript
/// order must run *down* the genome or the span comes out inverted and is
/// dropped (`ge >= gs` fails), leaving `find_intron_at_genomic` with nothing to
/// find.
fn two_intron_transcript(strand: Strand) -> Transcript {
    let exons = match strand {
        // Transcript order ascends the genome.
        Strand::Plus => vec![
            Exon::with_genomic(1, 1, 100, 1001, 1100),
            Exon::with_genomic(2, 101, 200, 1201, 1300),
            Exon::with_genomic(3, 201, 300, 1601, 1700),
        ],
        // Transcript order descends the genome; the gaps are the same two.
        _ => vec![
            Exon::with_genomic(1, 1, 100, 1601, 1700),
            Exon::with_genomic(2, 101, 200, 1201, 1300),
            Exon::with_genomic(3, 201, 300, 1001, 1100),
        ],
    };
    Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        strand,
        "A".repeat(300),
        Some(1),
        Some(300),
        exons,
        Some("chrTEST".to_string()),
        Some(1001),
        Some(1700),
        GenomeBuild::GRCh38,
        ManeStatus::default(),
        None,
        None,
    )
}

/// The two genomic ranges the transcript above leaves between its exons, with
/// the length each intron therefore has.
const INTRON_SPANS: [(u64, u64, u64); 2] = [(1101, 1200, 100), (1301, 1600, 300)];

/// The defect as reported: a *measured* offset larger than the intron it names.
/// Panicked at `transcript.rs:1029` before the fix (`100 - 300`); in release it
/// returned `Some(18446744073709551516)`.
#[test]
fn a_measured_offset_larger_than_the_intron_is_declined() {
    let pos = intron_position(IntronBoundary::FivePrime, 300, 100);
    assert_eq!(
        pos.distance_from_acceptor(),
        None,
        "a position 300 bases 3' of the donor is not inside a 100-base intron, \
         so it has no distance to that intron's acceptor"
    );

    let pos = intron_position(IntronBoundary::ThreePrime, -300, 100);
    assert_eq!(
        pos.distance_from_donor(),
        None,
        "the mirror image: 300 bases 5' of the acceptor of a 100-base intron"
    );
}

/// The boundary of the guard, stated as three adjacent magnitudes so a fix that
/// is off by one is caught.
///
/// Offsets are 1-based from their boundary, so a magnitude of exactly
/// `intron_length` names the intron's far-edge base — still inside the intron,
/// and still answerable. The `Some(0)` it yields is the shipped arithmetic's
/// own answer for that base, not a clamp: a magnitude of `intron_length + 1`
/// names the first base of the flanking exon and declines.
#[test]
fn the_guard_fires_exactly_where_the_magnitude_leaves_the_intron() {
    let length = 100;
    for (offset, expected) in [(99u64, Some(1u64)), (100, Some(0)), (101, None)] {
        let pos = intron_position(IntronBoundary::FivePrime, offset as i64, length);
        assert_eq!(
            pos.distance_from_acceptor(),
            expected,
            "offset {offset} against an intron of length {length}"
        );
    }
}

/// An unknown-offset sentinel states no distance in *either* direction, so both
/// arms of both functions must decline — including the same-boundary arm, which
/// never underflowed and answered `Some(9223372036854775807)`.
#[test]
fn an_unknown_offset_sentinel_is_declined_from_either_boundary() {
    for offset in SENTINELS {
        for boundary in [IntronBoundary::FivePrime, IntronBoundary::ThreePrime] {
            let pos = intron_position(boundary, offset, 1_000);
            assert_eq!(
                pos.distance_from_donor(),
                None,
                "offset {offset} from {boundary:?} states no distance to the donor"
            );
            assert_eq!(
                pos.distance_from_acceptor(),
                None,
                "offset {offset} from {boundary:?} states no distance to the acceptor"
            );
        }
    }
}

/// The assertion a one-armed fix cannot satisfy. The two functions describe one
/// position, so they must agree on whether that position is resolvable at all —
/// clamping only the arm that panics leaves `None` beside `Some(300)` for the
/// same record.
#[test]
fn the_two_distances_never_disagree_about_resolvability() {
    let offsets = [
        i64::MIN,
        i64::MIN + 1,
        OFFSET_UNKNOWN_NEGATIVE,
        -1_000_000,
        -301,
        -300,
        -100,
        -1,
        0,
        1,
        100,
        300,
        301,
        1_000_000,
        i64::MAX - 1,
        OFFSET_UNKNOWN_POSITIVE,
    ];
    for boundary in [IntronBoundary::FivePrime, IntronBoundary::ThreePrime] {
        for offset in offsets {
            for length in [0u64, 1, 100, 300, u64::MAX] {
                let pos = intron_position(boundary, offset, length);
                // Totality first: neither call may panic for any `i64`.
                let donor = pos.distance_from_donor();
                let acceptor = pos.distance_from_acceptor();
                assert_eq!(
                    donor.is_some(),
                    acceptor.is_some(),
                    "{boundary:?} offset {offset} length {length}: one distance \
                     resolved and the other did not, for the same position"
                );
            }
        }
    }
}

/// The negative control, and the reason this change moves no in-tree behaviour:
/// **every** position the library can actually derive must still answer.
///
/// `find_intron_at_genomic` maintains `dist_to_5prime + dist_to_3prime ==
/// intron_length + 1` with both terms at least 1, so a derived `IntronPosition`
/// always satisfies `offset.unsigned_abs() <= intron_length`. This walks every
/// genomic base of both introns, on both strands, and requires `Some` from both
/// functions at each one — so a guard drawn one base too tight fails here rather
/// than silently narrowing the API.
#[test]
fn every_derived_intronic_position_still_answers() {
    for strand in [Strand::Plus, Strand::Minus] {
        let transcript = two_intron_transcript(strand);
        let mut checked = 0usize;

        for (g_start, g_end, length) in INTRON_SPANS {
            for genomic_pos in g_start..=g_end {
                let (_, pos) = transcript
                    .find_intron_at_genomic(genomic_pos)
                    .unwrap_or_else(|| {
                        panic!("{genomic_pos} is inside intron {g_start}..={g_end}")
                    });
                assert_eq!(pos.intron_length, length, "at {genomic_pos} on {strand:?}");

                let donor = pos.distance_from_donor().unwrap_or_else(|| {
                    panic!(
                        "a derived position ({genomic_pos} on {strand:?}, offset {}) must have \
                         a donor distance",
                        pos.offset
                    )
                });
                let acceptor = pos.distance_from_acceptor().unwrap_or_else(|| {
                    panic!(
                        "a derived position ({genomic_pos} on {strand:?}, offset {}) must have \
                         an acceptor distance",
                        pos.offset
                    )
                });

                // The shipped arithmetic, restated as an invariant rather than
                // re-implemented: the two distances partition the intron.
                assert_eq!(
                    donor + acceptor,
                    length,
                    "at {genomic_pos} on {strand:?}, offset {}",
                    pos.offset
                );
                checked += 1;
            }
        }

        // A zero here would make every assertion above vacuous.
        assert_eq!(
            checked,
            INTRON_SPANS
                .iter()
                .map(|(_, _, l)| *l as usize)
                .sum::<usize>(),
            "every intronic base on {strand:?} must have been visited"
        );
    }
}

/// The shipped semantics of the answering path are unchanged. Pinned here as
/// well as in `transcript.rs`'s own unit test, because the fix rewrites the
/// bodies of both functions and this is the arithmetic it must not move.
#[test]
fn the_answering_path_keeps_its_shipped_arithmetic() {
    let pos = intron_position(IntronBoundary::FivePrime, 10, 500);
    assert_eq!(pos.distance_from_donor(), Some(10));
    assert_eq!(pos.distance_from_acceptor(), Some(490));

    let pos = intron_position(IntronBoundary::ThreePrime, -10, 500);
    assert_eq!(pos.distance_from_acceptor(), Some(10));
    assert_eq!(pos.distance_from_donor(), Some(490));
}
