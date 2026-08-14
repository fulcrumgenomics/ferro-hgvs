//! `Normalizer::sequence_normalize`: the "one canonical description per variant"
//! round trip, and its per-side widening loop.
//!
//! `sequence_normalize` expresses a variant as a padded reference/alternate
//! window (`to_sequences`), derives a description from those bases alone
//! (`from_sequences`), and — while a member still rests on a window edge that
//! can still move — doubles the pad and retries. Two spellings of one variant
//! therefore reach one description, decided by the observed bases.
//!
//! The loop reads the two per-side flags apart, so a placement pinned to the
//! sequence's own start (5') or end (3') is recognised as settled rather than
//! chased: a variant at position 1 of a long contig does not widen 5' forever
//! against an edge that cannot move.
//!
//! Every test builds a genomic `JsonProvider` over one synthetic contig, since
//! `sequence_normalize` reads the reference (unlike the free `from_sequences`).

use ferro_hgvs::{parse_hgvs, FromSequencesOptions, JsonProvider, Normalizer, ShuffleDirection};
use std::io::Write;

/// A genome-capable provider over one contig named `NC_TEST.1` carrying `seq`.
///
/// The name is genomic on purpose: `from_sequences` (which this method calls)
/// emits `g.` and refuses a transcript/protein accession, so a non-genomic
/// contig name would be rejected before any derivation.
fn provider(seq: &str) -> JsonProvider {
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { "NC_TEST.1": seq },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

/// Run `sequence_normalize` over a parsed description, returning the string.
fn seqnorm(
    nz: &Normalizer<JsonProvider>,
    desc: &str,
    direction: ShuffleDirection,
    normalize: bool,
) -> String {
    let variant = parse_hgvs(desc).unwrap_or_else(|e| panic!("parse {desc}: {e}"));
    let options = FromSequencesOptions::default().with_direction(direction);
    nz.sequence_normalize(&variant, &options, normalize)
        .unwrap_or_else(|e| panic!("sequence_normalize {desc}: {e}"))
        .to_string()
}

/// Three spellings of one deletion in a 300-base run reach one description, and
/// reaching it requires the loop to widen: the run is far longer than the
/// 128-base start pad (`merge::CANONICAL_PAD`), so a single fetch cannot contain
/// the whole roll.
#[test]
fn confluence_and_widening_over_a_long_interior_run() {
    // unique(10) + A*300 (positions 11..=310) + unique(12).
    let seq = format!("GCTAGCTAGC{}GCTAGCTAGCTA", "A".repeat(300));
    let nz = Normalizer::new(provider(&seq));

    // Delete one A from anywhere in the run: one variant, three spellings. All
    // must land on the 3'-most placement, 310.
    for spelling in [
        "NC_TEST.1:g.11del",
        "NC_TEST.1:g.160del",
        "NC_TEST.1:g.310del",
    ] {
        assert_eq!(
            seqnorm(&nz, spelling, ShuffleDirection::ThreePrime, false),
            "NC_TEST.1:g.310del",
            "{spelling} did not converge on the 3'-most placement",
        );
    }
}

/// A duplication and the insertion spelling of the same variant converge — the
/// case a narrow first window gets silently wrong.
///
/// `dup` typing reads the reference bases immediately 5' of the insertion point
/// (`duplication.md:18`), so a window that does not reach them derives an `ins`
/// instead. That mis-typed `ins` rests on **neither** window edge, so both
/// `bounded_*` flags are clear and the loop returns it without ever widening —
/// there is no flag for "the answer is wrong for a reason the edges cannot
/// show". Starting the loop at `merge::CANONICAL_PAD` rather than at a bare
/// couple of bases is what makes the two spellings agree; with a 1-base start
/// pad `g.14_15insACGTACGT` comes back unchanged while `g.7_14dup` does not.
#[test]
fn a_duplication_and_its_insertion_spelling_converge() {
    // unique(6) + ACGTACGT (positions 7..=14) + unique(6).
    let seq = "GCTAGC".to_string() + "ACGTACGT" + "TCGATC";
    let nz = Normalizer::new(provider(&seq));

    for spelling in ["NC_TEST.1:g.7_14dup", "NC_TEST.1:g.14_15insACGTACGT"] {
        assert_eq!(
            seqnorm(&nz, spelling, ShuffleDirection::ThreePrime, false),
            "NC_TEST.1:g.7_14dup",
            "{spelling} did not converge on the duplication spelling",
        );
    }
}

/// A run reaching the contig's 3' end settles by the sequence, not the window:
/// the 3' edge is the last base there is, so the loop returns rather than
/// declining or widening past it.
#[test]
fn a_run_reaching_the_three_prime_end_settles_by_the_sequence() {
    // unique(10) + A*40 (positions 11..=50, the contig ends inside the run).
    let seq = format!("GCTAGCTAGC{}", "A".repeat(40));
    let nz = Normalizer::new(provider(&seq));

    assert_eq!(
        seqnorm(
            &nz,
            "NC_TEST.1:g.15del",
            ShuffleDirection::ThreePrime,
            false
        ),
        "NC_TEST.1:g.50del",
        "a deletion in a contig-terminal run must roll to the last base",
    );
}

/// The case the per-side split exists for: a variant whose placement is pinned
/// to position 1. Under 5'-shuffle a deletion in a run at the contig start rolls
/// to base 1, resting on the 5' edge of a window that already begins there — a
/// stall, not an overrun. The loop must return `g.1del`, not chase an edge that
/// cannot move nor decline it as unbounded.
#[test]
fn a_run_at_the_contig_start_settles_at_position_one_under_five_prime() {
    // A*40 (positions 1..=40, the contig begins inside the run) + unique(10).
    let seq = format!("{}GCTAGCTAGC", "A".repeat(40));
    let nz = Normalizer::new(provider(&seq));

    assert_eq!(
        seqnorm(&nz, "NC_TEST.1:g.20del", ShuffleDirection::FivePrime, false),
        "NC_TEST.1:g.1del",
        "a 5'-shuffled deletion in a contig-initial run must settle at base 1",
    );
}

/// `normalize = true` routes the derived description through `normalize`. For a
/// placement already at its 3'-most, reference-anchored base that is a no-op, so
/// the result matches the `normalize = false` derivation — what this pins is
/// that the routing runs without error and returns the same settled variant.
#[test]
fn normalize_true_routes_through_normalize() {
    let seq = format!("GCTAGCTAGC{}GCTAGCTAGCTA", "A".repeat(300));
    let nz = Normalizer::new(provider(&seq));

    assert_eq!(
        seqnorm(&nz, "NC_TEST.1:g.11del", ShuffleDirection::ThreePrime, true),
        "NC_TEST.1:g.310del",
    );
}

/// An interior tract longer than the widest window the loop will read is
/// declined, not answered over a truncated window: where the change shifts to
/// would depend on how much reference was read, which is exactly the
/// non-confluence the widening removes. The message names the offending side.
#[test]
fn an_unbounded_interior_tract_is_declined() {
    // A run longer than MAX_SHIFT_TRACT (32768), flanked by unique sequence on
    // both sides so neither edge is ever the contig's own — the loop can never
    // settle the 3' side, and reaches the pad ceiling.
    let seq = format!("GCTAGCTAGC{}GCTAGCTAGCTA", "A".repeat(40_000));
    let nz = Normalizer::new(provider(&seq));

    let variant = parse_hgvs("NC_TEST.1:g.15del").unwrap();
    let options = FromSequencesOptions::default().with_direction(ShuffleDirection::ThreePrime);
    let err = nz
        .sequence_normalize(&variant, &options, false)
        .expect_err("an unbounded interior tract must be declined");
    let msg = err.to_string();
    assert!(
        msg.contains("3' edge") && msg.contains("window-settled"),
        "the refusal must name the unsettled side: {msg}",
    );
}
