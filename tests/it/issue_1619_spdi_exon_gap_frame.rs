//! #1619 — `hgvs_to_spdi` must index the same base the normalizer reads.
//!
//! `Exon::start`/`end` are documented as **transcript** coordinates, and a
//! transcript sequence is the concatenation of its exons, so it has no holes.
//! `CoordinateMapper` nonetheless walked the exon list whenever those
//! coordinates showed a gap, skipping each gap as though it were an intron —
//! and the normalizer does no such thing: `region_sequence_delta` maps `c.N` to
//! the flat offset `cds_start + N - 1`. So the two halves of ferro named
//! different bases of one accession, and the gap grew by one per "intron"
//! crossed.
//!
//! Measured on `main`, on the committed
//! `tests/fixtures/sequences/normalization_transcripts.json`, whose
//! `NM_000492.4` record carries the real 6,070-base CFTR transcript,
//! `cds_start: 71`, and 27 exons spelled with a one-base gap between each pair:
//! `c.1520` is flat offset 1589 and `hgvs_to_spdi` placed it at **1599**, the
//! ten gaps its walk crossed.
//!
//! The visible symptom was a correct normalization reading as a corruption.
//! `NM_000492.4:c.1520_1522del` -> `c.1521_1523del` is the documented CFTR
//! 3'-shift — `c.1515..1530` reads `TATCATCTTTGGTGTT`, so deleting either `TCT`
//! or `CTT` leaves `TATCATTGGTGTT` over that window — but read through the walk the pair landed on
//! `[1599, 1603) = TTCC`, where the two deletions leave `C` and `T`. That is
//! the row the denoted-sequence oracle (#1615) fires on, and one of the two
//! keeping `FERRO_ASSERT_SEQUENCE` out of CI's `test-oracle` job.
//!
//! The gaps are a property of that hand-written fixture, not of real data: no
//! cdot-derived transcript has them, because `MultiFastaProvider` builds each
//! exon's transcript interval from cdot's own tiling. 22 of the fixture's 29
//! records carry them.

use std::path::Path;

use ferro_hgvs::reference::mock::JsonProvider;
use ferro_hgvs::reference::provider::ReferenceProvider;
use ferro_hgvs::spdi::convert::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, Normalizer};

const FIXTURE: &str = "tests/fixtures/sequences/normalization_transcripts.json";

/// The fixture's records all carry their sequence; a missing one is a fixture
/// problem, not a case to skip.
fn sequence_of(transcript: &ferro_hgvs::reference::transcript::Transcript) -> &str {
    transcript
        .sequence
        .as_deref()
        .expect("the fixture record carries its sequence")
}

fn provider() -> JsonProvider {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(FIXTURE);
    JsonProvider::from_json(&path).expect("load normalization_transcripts.json")
}

/// The fixture really does carry transcript-coordinate gaps, so the test below
/// is exercising the shape it claims to. Asserted rather than assumed: if the
/// fixture is ever regenerated contiguously, this test stops proving anything
/// and should say so loudly rather than pass for the wrong reason.
#[test]
fn the_fixture_still_carries_the_gaps_that_made_the_two_halves_disagree() {
    let transcript = provider()
        .get_transcript("NM_000492.4")
        .expect("the fixture serves NM_000492.4");
    let mut exons: Vec<_> = transcript.exons.iter().collect();
    exons.sort_by_key(|e| e.start);
    let gaps = exons
        .windows(2)
        .filter(|w| w[0].end + 1 != w[1].start)
        .count();
    assert_eq!(
        gaps, 26,
        "NM_000492.4's 27 exons no longer show 26 transcript-coordinate gaps; \
         this test's premise has moved"
    );
    let spanned: u64 = exons.iter().map(|e| e.end - e.start + 1).sum();
    assert_eq!(
        (spanned, sequence_of(&transcript).len()),
        (6043, 6070),
        "the fixture's internal inconsistency — exon spans that do not sum to \
         the sequence length — has moved"
    );
}

/// `c.1520` must resolve to the base the normalizer reads, flat offset 1589.
///
/// Asserted through the SPDI position rather than through an internal helper,
/// because the position `to_spdi` publishes is what every caller downstream of
/// it — `apply_to_reference`, `canonical_spdi`, `compare_denoted_sequences` —
/// actually indexes with.
#[test]
fn a_cds_position_resolves_onto_the_flat_transcript() {
    let provider = provider();
    let variant = parse_hgvs("NM_000492.4:c.1520_1522del").expect("fixture must parse");
    let spdi = hgvs_to_spdi(&variant, &provider).expect("converts against the fixture");

    // cds_start - 1 + 1520 - 1 = 70 + 1519 = 1589, 0-based.
    assert_eq!(
        spdi.position, 1589,
        "c.1520 must land on the flat transcript offset the normalizer reads; \
         1599 is the exon walk crossing ten synthetic gaps"
    );

    let transcript = provider
        .get_transcript("NM_000492.4")
        .expect("the fixture serves NM_000492.4");
    assert_eq!(
        &sequence_of(&transcript)[1589..1592],
        "TCT",
        "the deleted bases the description names"
    );
    assert_eq!(spdi.deletion, "TCT");
}

/// The invariant behind the row, and the reason this is filed as a defect in
/// the applier rather than in the normalizer: the two halves of ferro must
/// agree on which bases a description denotes.
///
/// `c.1520_1522del` 3'-shifts to `c.1521_1523del`, and both spellings delete
/// one base from the same `TCTTT` run. Applied over a window containing both,
/// they must leave the same sequence — which is exactly what the walk broke,
/// and it broke it while the normalization itself was right.
#[test]
fn the_shifted_pair_denotes_one_sequence_to_the_applier() {
    let provider = provider();
    let input = parse_hgvs("NM_000492.4:c.1520_1522del").expect("fixture must parse");
    let normalized = Normalizer::new(provider.clone())
        .normalize(&input)
        .expect("normalizes against the fixture");
    assert_eq!(
        normalized.to_string(),
        "NM_000492.4:c.1521_1523del",
        "the documented CFTR 3'-shift"
    );

    let before = hgvs_to_spdi(&input, &provider).expect("input converts");
    let after = hgvs_to_spdi(&normalized, &provider).expect("output converts");

    let transcript = provider
        .get_transcript("NM_000492.4")
        .expect("the fixture serves NM_000492.4");
    let sequence = sequence_of(&transcript);
    // The window is `c.1515`..`c.1530` on the flat transcript
    // (`cds_start - 1 + N - 1`), wide enough to contain both spellings — a
    // per-description window would give each its own frame and report every
    // 3'-shift as a difference.
    const WINDOW_START: usize = 1584;
    const WINDOW_END: usize = 1600;
    assert_eq!(
        &sequence[WINDOW_START..WINDOW_END],
        "TATCATCTTTGGTGTT",
        "the CFTR homopolymer the 3'-shift runs in"
    );
    let applied = |spdi: &ferro_hgvs::spdi::SpdiVariant| {
        let start = spdi.position as usize;
        let end = start + spdi.deletion.len();
        format!(
            "{}{}{}",
            &sequence[WINDOW_START..start],
            spdi.insertion,
            &sequence[end..WINDOW_END]
        )
    };
    assert_eq!(
        applied(&before),
        applied(&after),
        "the input and its 3'-shifted normalization denote different bases to \
         the applier; the normalization is correct, so the applier is wrong"
    );
    assert_eq!(
        applied(&before),
        "TATCATTGGTGTT",
        "and both leave one base short of the `TTT` run, over the window above"
    );
}
