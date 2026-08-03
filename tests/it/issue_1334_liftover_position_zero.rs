//! #1334 — liftover must refuse position `0` rather than panic on it.
//!
//! `Liftover::lift` and `Liftover::lift_interval` take raw `u64` positions on a
//! `pub` API. The web service reaches them with an unvalidated value —
//! `parse_genomic_position` in `src/service/handlers/liftover.rs` parses with a
//! bare `pos_part.parse::<u64>()` and applies no lower bound, so `chr7:0` yields
//! `Ok(("chr7", 0))` — and the `ferro liftover` CLI reaches `lift` the same way.
//!
//! Since #1282 made `hgvs_pos_to_index` assert on `0`, that input panicked: in
//! the service it killed the handler task and dropped the connection, where
//! before #1282 it had wrapped to an index near `usize::MAX` in release (the
//! release profile sets no `overflow-checks`), found no chains, and rendered a
//! graceful error body. Neither is acceptable; #1282 asked for the guard:
//!
//! > Guard the conversion (return `None`/`Err` for `pos < 1`) rather than
//! > relying on callers never producing a zero.
//!
//! This file pins the guarded behaviour from *outside* the crate, which is where
//! the offending callers live — the in-module unit tests in `src/liftover/lift.rs`
//! exercise the same paths but cannot show that the failure is reachable and
//! catchable across the public API boundary.
//!
//! There is no oracle movement to report: liftover is chain-file arithmetic, not
//! an HGVS description, so no reference oracle (biocommons hgvs, Mutalyzer,
//! VariantValidator) has an opinion about it. For every position these entry
//! points previously answered — every `pos >= 1` — the answer is byte-identical;
//! the only inputs whose behaviour moves are the ones that used to panic.

use ferro_hgvs::error::FerroError;
use ferro_hgvs::liftover::{ChainFile, Liftover};
use ferro_hgvs::reference::transcript::GenomeBuild;

/// A single 1:1 chain covering target `chr1:10001..20000` (1-based inclusive),
/// mapping onto query `chr1:10051..20150` with a 100 bp offset introduced by the
/// gap between the two 5000 bp blocks.
///
/// Built inline rather than read from a fixture so the corpus stays generated
/// rather than committed.
const CHAIN: &str =
    "chain 1000 chr1 1000000 + 10000 20000 chr1 1000100 + 10050 20150 1\n5000\n5000\n\n";

fn liftover() -> Liftover {
    Liftover::one_way(ChainFile::parse(CHAIN.as_bytes()).expect("chain parses"))
}

/// The reported reproducer: `lift` with `pos == 0`.
///
/// Pins the *new* output — a returned `InvalidCoordinates` naming the offending
/// coordinate and contig — rather than merely "is an error", so a future change
/// that swaps the variant or drops the coordinate from the message fails here.
#[test]
fn lift_returns_invalid_coordinates_for_position_zero() {
    let err = liftover()
        .lift(GenomeBuild::GRCh37, GenomeBuild::GRCh38, "chr7", 0)
        .expect_err("position 0 must be an error, not a panic");

    let FerroError::InvalidCoordinates { msg } = &err else {
        panic!("expected InvalidCoordinates, got {err:?}");
    };
    assert!(
        msg.contains("chr7"),
        "error must name the contig the caller passed: {msg}"
    );
    assert!(
        msg.contains("1-based"),
        "error must say what was wrong with the position: {msg}"
    );
}

/// Both interval endpoints carry the same guard.
///
/// This is the case a partial fix loses: guarding only `start` leaves
/// `lift_interval(.., 10001, 0)` converting `end` with the panicking helper, and
/// every test that only ever zeroes the start still passes.
#[test]
fn lift_interval_returns_invalid_coordinates_for_a_zero_at_either_end() {
    for (start, end) in [(0u64, 10100u64), (10001, 0), (0, 0)] {
        let err = liftover()
            .lift_interval(GenomeBuild::GRCh37, GenomeBuild::GRCh38, "chr1", start, end)
            .expect_err("a zero endpoint must be an error, not a panic");
        assert!(
            matches!(err, FerroError::InvalidCoordinates { .. }),
            "expected InvalidCoordinates for ({start}, {end}), got {err:?}"
        );
    }
}

/// The negative control: positions the guard must not touch still lift to the
/// coordinates they lifted to before it existed.
///
/// Without this, a guard that rejected far more than `0` — say `pos <= 1`, or
/// every position — would satisfy both tests above and silently break liftover.
#[test]
fn valid_positions_still_lift_to_their_previous_coordinates() {
    let liftover = liftover();

    // First base of the chain's first block: 1:1, no gap offset yet.
    let first = liftover
        .lift(GenomeBuild::GRCh37, GenomeBuild::GRCh38, "chr1", 10001)
        .expect("the first covered base lifts");
    assert_eq!(first.source_pos, 10001);
    assert_eq!(first.target_contig, "chr1");
    assert_eq!(first.target_pos, 10051);
    assert!(!first.in_gap);

    // An interval wholly inside that block lifts both endpoints the same way.
    let interval = liftover
        .lift_interval(
            GenomeBuild::GRCh37,
            GenomeBuild::GRCh38,
            "chr1",
            10001,
            10100,
        )
        .expect("a covered interval lifts");
    assert_eq!(interval.source_interval, (10001, 10100));
    assert_eq!(interval.target_interval, Some((10051, 10150)));

    // And position 1 — the smallest position that exists — still reaches the
    // chain lookup and fails there as "no chain covers this", not at the guard.
    //
    // Both failures are spelled `InvalidCoordinates`, so the discriminator has
    // to be the message: the uncovered path names the chain lookup, the guard
    // names the 1-based convention. That shared variant is pre-existing and out
    // of scope here; what matters is that the guard did not swallow position 1.
    let err = liftover
        .lift(GenomeBuild::GRCh37, GenomeBuild::GRCh38, "chr1", 1)
        .expect_err("position 1 is outside the chain");
    let FerroError::InvalidCoordinates { msg } = &err else {
        panic!("expected InvalidCoordinates, got {err:?}");
    };
    assert!(
        msg.contains("No chain found"),
        "position 1 is valid; it must fail at the chain lookup, not at the \
         1-based guard: {msg}"
    );
}
