//! The two coordinate frames must **disagree** on a gapped transcript.
//!
//! `CoordinateMapper` carries two families of conversion that look like
//! duplicates and are not:
//!
//! - the **sequence frame** — `cds_to_tx` / `tx_to_cds`. A `c.`/`n.` position
//!   names a base of the transcript *reference sequence*, so `c.N` is
//!   `cds_start + N - 1` on the flat transcript, with no exon walk.
//!   `ReferenceProvider::get_sequence` on a transcript accession serves that
//!   flat sequence, so this is the only offset that can index it.
//! - the **genome frame** — `genomic_to_tx` / `tx_to_genomic` and the projector
//!   above them. Mapping between a contig and a transcript genuinely is exon-
//!   and CIGAR-aware, and a transcript-coordinate gap there correctly means
//!   "this transcript base aligns to nothing".
//!
//! On a transcript whose cdot alignment carries a transcript-coordinate gap the
//! two answers **differ, by construction**, for every position 3' of the gap.
//! That divergence is the design, not a defect: #1619 was caused by making the
//! sequence frame walk the exon table, which moved every `c.` position 3' of
//! `NM_033517.1`'s real 39-base hole while `get_sequence` kept serving the flat
//! transcript — a correct normalization then read as a corruption. PR #1735
//! settled it: each family is correct for its own axis and neither may borrow
//! the other's arithmetic.
//!
//! Nothing in the tree previously *pinned* that, so the two families read as
//! duplication inviting a "cleanup" that collapses them. This file is that pin.
//! If it goes red, a change has made the two frames agree — read the assertion
//! message before concluding the test is wrong.

use ferro_hgvs::convert::CoordinateMapper;
use ferro_hgvs::hgvs::location::{CdsPos, TxPos};
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};

/// Flat length of the transcript reference sequence, in bases.
const TX_LEN: usize = 400;
/// 1-based transcript position of `c.1`.
const CDS_START: u64 = 51;
/// 1-based transcript position of the last CDS base.
const CDS_END: u64 = 350;
/// First transcript base of the transcript-coordinate gap (aligns to no exon).
const GAP_TX_START: u64 = 201;
/// Last transcript base of the transcript-coordinate gap.
const GAP_TX_END: u64 = 239;
/// Width of the gap, in transcript bases. `NM_033517.1` carries a real hole of
/// exactly this size, inside a CDS that RefSeq annotates as contiguous across
/// it.
const GAP_WIDTH: u64 = GAP_TX_END - GAP_TX_START + 1;

/// Build a gapped transcript with the geometry of `NM_033517.1`: a
/// transcript-coordinate gap *inside* the CDS, with the CDS annotated
/// contiguously across it.
///
/// Constructed programmatically on purpose — this guard must assert on every
/// run, so it may not depend on a prepared reference or a downloaded cdot that
/// could be absent and turn the test into a silent skip.
///
/// Geometry (1-based, inclusive):
///
/// | span | transcript | genome |
/// |---|---|---|
/// | exon 1 | 1..=200 | 10_001..=10_200 |
/// | **gap** | 201..=239 | *aligns to nothing* |
/// | exon 2 | 240..=400 | 12_001..=12_161 |
///
/// The two exons are separated by an intron as well, so "advance one aligned
/// base" is genuinely an exon-table walk and not plain genomic arithmetic.
fn gapped_transcript() -> Transcript {
    // Background 'A', with two distinguishable marker bases planted at the two
    // transcript indices the two frames name for `c.220` (see
    // `MARKED_CDS_POS`). Their identity is what proves the frames name
    // different *bases*, not merely different integers.
    let mut sequence = vec![b'A'; TX_LEN];
    sequence[(SEQUENCE_FRAME_TX - 1) as usize] = b'C';
    sequence[(GENOME_FRAME_WALK_TX - 1) as usize] = b'G';

    Transcript::new(
        "NM_GAPPED.1".to_string(),
        Some("GAPTEST".to_string()),
        Strand::Plus,
        String::from_utf8(sequence).expect("marker sequence is valid UTF-8"),
        Some(CDS_START),
        Some(CDS_END),
        vec![
            Exon::with_genomic(1, 1, 200, 10_001, 10_200),
            Exon::with_genomic(2, 240, 400, 12_001, 12_161),
        ],
        Some("chr1".to_string()),
        Some(10_001),
        Some(12_161),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    )
}

/// A CDS position 3' of the gap, used as the worked example throughout.
const MARKED_CDS_POS: i64 = 220;
/// What the **sequence frame** answers for `c.220`: `cds_start + 220 - 1`.
const SEQUENCE_FRAME_TX: u64 = CDS_START + MARKED_CDS_POS as u64 - 1; // 270
/// What a **genome-frame** walk answers for `c.220`: the same count of
/// *aligned* bases, which steps over the 39-base hole and lands 39 further on.
const GENOME_FRAME_WALK_TX: u64 = SEQUENCE_FRAME_TX + GAP_WIDTH; // 309
/// A CDS position 5' of the gap, where the two frames must still agree.
const UNGAPPED_CDS_POS: i64 = 100;

/// Place `c.<cds_pos>` by walking the **genome frame**: start at `c.1` and
/// advance `cds_pos - 1` transcript bases that are actually aligned to the
/// genome, skipping any that are not.
///
/// This is the exact operation #1619 was: placing a `c.` position through the
/// exon table. It is written entirely in terms of the public genome-frame API
/// (`tx_to_genomic`), so unifying the two frames in *either* direction changes
/// what it returns and reddens the assertions below.
fn place_by_genome_frame_walk(mapper: &CoordinateMapper<'_>, cds_pos: i64) -> u64 {
    assert!(cds_pos >= 1, "walk is defined for coding positions only");
    let mut tx = CDS_START;
    let mut remaining = cds_pos as u64 - 1;
    while remaining > 0 {
        tx += 1;
        assert!(
            tx <= TX_LEN as u64,
            "walked off the transcript placing c.{cds_pos}",
        );
        let aligned = mapper
            .tx_to_genomic(&TxPos::new(tx as i64))
            .expect("transcript has genomic coordinates")
            .is_some();
        if aligned {
            remaining -= 1;
        }
    }
    tx
}

/// The gap is real: its transcript bases exist in the sequence frame and have
/// `c.` numbers, and they align to nothing in the genome frame.
///
/// This is the whole disagreement in one fixture. If either half of it stops
/// holding, the two frames have been collapsed onto one.
#[test]
fn a_transcript_coordinate_gap_exists_in_one_frame_and_not_the_other() {
    let tx = gapped_transcript();
    let mapper = CoordinateMapper::new(&tx);

    for gap_tx in [GAP_TX_START, GAP_TX_START + 1, GAP_TX_END] {
        // Genome frame: nothing there.
        assert_eq!(
            mapper
                .tx_to_genomic(&TxPos::new(gap_tx as i64))
                .expect("transcript has genomic coordinates"),
            None,
            "n.{gap_tx} is inside the transcript-coordinate gap, so the genome \
             frame must report no aligned genomic base. A Some(..) here means \
             tx_to_genomic has stopped consulting the exon table.",
        );

        // Sequence frame: a real base, with a real CDS number.
        let cds = mapper
            .tx_to_cds(&TxPos::new(gap_tx as i64))
            .expect("flat conversion cannot fail inside the transcript");
        assert_eq!(
            cds.base,
            (gap_tx - CDS_START + 1) as i64,
            "n.{gap_tx} is a base of the transcript reference sequence, so it \
             has a flat c. number even though it aligns to nothing. Only the \
             genome has a hole here; the sequence does not.",
        );
        assert!(
            !cds.utr3,
            "n.{gap_tx} is inside the CDS ({CDS_START}..={CDS_END})",
        );
    }
}

/// The two frames must place a coding position 3' of the gap on **different**
/// transcript bases, and the difference must be exactly the gap width.
///
/// This is the guard the m8 finding asks for. A change that "unifies" the CDS↔tx
/// conversion with the exon-aware mapping — in either direction — makes these
/// two numbers equal and turns this red.
#[test]
fn sequence_frame_and_genome_frame_disagree_past_a_gap() {
    let tx = gapped_transcript();
    let mapper = CoordinateMapper::new(&tx);

    let sequence_frame = mapper
        .cds_to_tx(&CdsPos::new(MARKED_CDS_POS))
        .expect("flat conversion cannot fail inside the CDS")
        .base as u64;
    let genome_frame = place_by_genome_frame_walk(&mapper, MARKED_CDS_POS);

    assert_eq!(
        sequence_frame, SEQUENCE_FRAME_TX,
        "c.{MARKED_CDS_POS} must be cds_start + N - 1 on the flat transcript \
         ({CDS_START} + {MARKED_CDS_POS} - 1 = {SEQUENCE_FRAME_TX}). \
         `get_sequence` on a transcript accession serves the flat sequence, so \
         an exon-walked offset cannot index it — that was #1619.",
    );
    assert_eq!(
        genome_frame,
        GENOME_FRAME_WALK_TX,
        "walking {} aligned bases from c.1 must step over the {GAP_WIDTH}-base \
         hole and land on n.{GENOME_FRAME_WALK_TX}. A different answer means \
         tx_to_genomic no longer reports the gap.",
        MARKED_CDS_POS - 1,
    );
    assert_ne!(
        sequence_frame, genome_frame,
        "The sequence frame and the genome frame MUST disagree here, and this \
         disagreement is correct — they answer different questions. \
         `cds_to_tx` (src/convert/mapper.rs) places c.N on the transcript's own \
         bases, which is what `get_sequence` and the normalizer index; \
         `tx_to_genomic` walks the exon alignment, where a \
         transcript-coordinate gap really is a hole. If you got here by \
         collapsing the two onto one implementation, you have reintroduced \
         #1619 (adjudicated by PR #1735): the c. form and the base it names \
         come apart on the {GAP_WIDTH}-base hole real transcripts such as \
         NM_033517.1 carry. Note also that gaps are not malformed fixture data \
         — 58 of 474,818 GRCh38 builds and 159 of 190,754 GRCh37 builds carry \
         one, which is why PR #1665 was closed.",
    );
    assert_eq!(
        genome_frame - sequence_frame,
        GAP_WIDTH,
        "the two frames must diverge by exactly the width of the gap between \
         them, not by some other amount",
    );
}

/// …and they must **agree** 5' of the gap. The disagreement above is caused by
/// the gap, not by the two families being unrelated — which is what makes
/// collapsing them look safe on any ungapped fixture.
#[test]
fn sequence_frame_and_genome_frame_agree_before_the_gap() {
    let tx = gapped_transcript();
    let mapper = CoordinateMapper::new(&tx);

    let sequence_frame = mapper
        .cds_to_tx(&CdsPos::new(UNGAPPED_CDS_POS))
        .expect("flat conversion cannot fail inside the CDS")
        .base as u64;
    let genome_frame = place_by_genome_frame_walk(&mapper, UNGAPPED_CDS_POS);

    assert_eq!(
        sequence_frame,
        CDS_START + UNGAPPED_CDS_POS as u64 - 1,
        "c.{UNGAPPED_CDS_POS} is 5' of the gap",
    );
    assert_eq!(
        sequence_frame, genome_frame,
        "5' of the gap the two frames coincide. Every ungapped transcript looks \
         like this, which is exactly why a collapse of the two passes casual \
         review — see the sibling test for the case that catches it.",
    );
}

/// The two frames do not merely produce different integers: they name different
/// **bases** of the served transcript sequence.
///
/// `ReferenceProvider::get_sequence` on a transcript accession serves the flat
/// sequence, so the sequence-frame offset reads the base the `c.` form denotes
/// and the genome-frame offset reads an unrelated one.
#[test]
fn the_two_frames_name_different_bases_of_the_served_transcript() {
    let tx = gapped_transcript();
    let mapper = CoordinateMapper::new(&tx);
    let sequence = tx
        .sequence
        .as_deref()
        .expect("fixture carries a flat transcript sequence")
        .as_bytes();
    assert_eq!(sequence.len(), TX_LEN, "flat sequence has no gap in it");

    let sequence_frame = mapper
        .cds_to_tx(&CdsPos::new(MARKED_CDS_POS))
        .expect("flat conversion cannot fail inside the CDS")
        .base as u64;
    let genome_frame = place_by_genome_frame_walk(&mapper, MARKED_CDS_POS);

    assert_eq!(
        sequence[(sequence_frame - 1) as usize],
        b'C',
        "the sequence-frame offset must index the base c.{MARKED_CDS_POS} \
         denotes on the transcript the reference provider actually serves",
    );
    assert_eq!(
        sequence[(genome_frame - 1) as usize],
        b'G',
        "the genome-frame walk lands {GAP_WIDTH} bases further on, on a \
         different base entirely — which is why an exon-walked offset must \
         never be used to index the flat transcript (#1619)",
    );
}

/// The flat round trip closes across the gap. `tx_to_cds` is the exact inverse
/// of `cds_to_tx`, on the same axis, and the gap changes nothing about that.
#[test]
fn the_sequence_frame_round_trips_across_the_gap() {
    let tx = gapped_transcript();
    let mapper = CoordinateMapper::new(&tx);

    // Positions before, inside, and after the gap.
    for cds_pos in [1, 100, 150, 151, 170, 189, 190, MARKED_CDS_POS, 300] {
        let tx_pos = mapper
            .cds_to_tx(&CdsPos::new(cds_pos))
            .expect("flat conversion cannot fail inside the CDS");
        assert_eq!(
            tx_pos.base,
            CDS_START as i64 + cds_pos - 1,
            "c.{cds_pos} is a flat offset from cds_start",
        );
        let back = mapper
            .tx_to_cds(&tx_pos)
            .expect("flat inverse cannot fail inside the CDS");
        assert_eq!(
            back.base, cds_pos,
            "c.{cds_pos} must round-trip through the sequence frame; the gap is \
             invisible to it",
        );
    }
}
