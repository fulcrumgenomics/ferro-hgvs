//! Hermetic pin on the **two** cdot coordinate layers, and on which is which.
//!
//! # Why this exists
//!
//! `cdot` coordinates appear in two layers with the *same field names*, and the
//! two disagree:
//!
//! | layer | genomic | transcript |
//! |---|---|---|
//! | raw cdot JSON, on disk | 0-based half-open | **1-based inclusive** |
//! | ferro's in-memory `CdotTranscript` (post-#742) | 1-based `[incl, excl)` | **0-based half-open** |
//!
//! `coords::cdot_tx_coords` documents the *raw* layer; the basis table at the
//! top of `data::cdot` documents the *internal* one. Two tables, two layers, one
//! set of names — which is how two independent audits of that ambiguity each
//! resolved it the wrong way.
//!
//! # Why it is here rather than beside the code
//!
//! The test that is cited as *the* pin on the internal transcript basis —
//! `ferro_hgvs::convert::mapper`'s `test_cdot_tx_coordinates_are_0_based` —
//! opens by looking for `benchmark-output/cdot/cdot-0.2.32.refseq.GRCh38.json`
//! and `return`ing if it is absent. That file is a multi-gigabyte download
//! present on approximately no default checkout, so on almost every machine and
//! almost every CI job that test asserts **nothing** and reports green. A skip
//! that reads as a pass is the failure mode `tests/fixtures/CORPUS_LAYOUT.md`
//! exists to remove, and the convention this repo uses elsewhere is to pair any
//! skipping tier with a hermetic tier that never skips.
//!
//! **This module is that hermetic tier.** Its fixture is a synthetic cdot
//! payload in this file, so it downloads nothing, reads no LFS object, consults
//! no environment variable, and therefore has no way to skip. Fixing the
//! skipping sibling in place is a separate change to `src/convert/mapper.rs`;
//! this guard is what makes the invariant asserted on every run in the meantime.

use ferro_hgvs::coords::{cdot_genomic_to_closed, cdot_tx_coords};
use ferro_hgvs::data::cdot::CdotMapper;

/// Raw cdot JSON: one plus-strand transcript, two exons.
///
/// Exon rows are `[genome_start, genome_end, exon_number, tx_start, tx_end, cigar]`
/// in cdot's own conventions — genome 0-based half-open, transcript **1-based
/// inclusive**. So exon 1 covers 100 genomic bases at 0-based `[1000, 1100)` and
/// transcript bases 1..=100; exon 2 covers 50 at `[2000, 2050)` and 101..=150.
const RAW_CDOT_JSON: &str = r#"
{
    "transcripts": {
        "NM_LAYER.1": {
            "gene_name": "LAYERTEST",
            "genome_builds": {
                "GRCh38": {
                    "contig": "NC_000001.11",
                    "strand": "+",
                    "exons": [
                        [1000, 1100, 0, 1, 100, "M100"],
                        [2000, 2050, 1, 101, 150, "M50"]
                    ]
                }
            },
            "start_codon": 30,
            "stop_codon": 120
        }
    }
}
"#;

/// The raw values written into [`RAW_CDOT_JSON`], named so the assertions below
/// can state the before/after rather than restating literals.
mod raw {
    /// Exon 1: 0-based inclusive genomic start / 0-based exclusive genomic end.
    pub const EXON1_GENOME: (u64, u64) = (1000, 1100);
    /// Exon 1: 1-based **inclusive** transcript start / end.
    pub const EXON1_TX: (u64, u64) = (1, 100);
    /// Exon 2, same conventions.
    pub const EXON2_GENOME: (u64, u64) = (2000, 2050);
    /// Exon 2, same conventions.
    pub const EXON2_TX: (u64, u64) = (101, 150);
}

fn load() -> ferro_hgvs::data::cdot::CdotTranscript {
    let mapper = CdotMapper::from_reader_with_build(RAW_CDOT_JSON.as_bytes(), "GRCh38")
        .expect("the synthetic cdot fixture must parse");
    mapper
        .get_transcript("NM_LAYER.1")
        .expect("the synthetic cdot fixture must contain NM_LAYER.1")
        .clone()
}

/// The fixture is real, so nothing downstream can be vacuously green.
///
/// Asserted before the basis checks, in the house style, so an empty or
/// silently-dropped exon table cannot make them pass over zero rows.
#[test]
fn the_synthetic_fixture_loads_and_is_not_empty() {
    let tx = load();
    assert_eq!(
        tx.exons.len(),
        2,
        "the fixture must load both exon rows; the ingestion path silently \
         drops rows with fewer than five fields, so a shortfall here would make \
         every basis assertion below a claim about a corpus that is not there"
    );
}

/// The claim `test_cdot_tx_coordinates_are_0_based` makes, asserted hermetically.
///
/// Raw cdot spells the first exon's transcript bounds `1..=100`; after loading,
/// the first exon must start at transcript **0**, not 1.
#[test]
fn raw_one_based_transcript_bounds_load_as_zero_based_half_open() {
    let tx = load();

    assert_eq!(
        raw::EXON1_TX.0,
        1,
        "precondition: the fixture must spell the first exon's transcript start \
         as raw cdot does (1-based), or this test proves nothing about the \
         conversion"
    );

    assert_eq!(
        tx.exons[0][2], 0,
        "first exon tx_start must be 0 (0-based half-open) after loading, not \
         the raw JSON's 1"
    );
    assert_eq!(
        (tx.exons[0][2], tx.exons[0][3]),
        (raw::EXON1_TX.0 - 1, raw::EXON1_TX.1),
        "tx bounds convert 1-based inclusive -> 0-based half-open: start -1, \
         end unchanged"
    );
    assert_eq!(
        (tx.exons[1][2], tx.exons[1][3]),
        (raw::EXON2_TX.0 - 1, raw::EXON2_TX.1),
    );

    // Both exons are exactly as long on the transcript axis as their genomic
    // spans, so the conversion cannot have shifted one end without the other.
    assert_eq!(tx.exons[0][3] - tx.exons[0][2], 100);
    assert_eq!(tx.exons[1][3] - tx.exons[1][2], 50);
}

/// The genomic axis converts in the opposite direction on the same rows, which
/// is the mixed-basis trap the two tables exist to warn about.
#[test]
fn raw_zero_based_genomic_bounds_load_as_one_based() {
    let tx = load();
    assert_eq!(
        (tx.exons[0][0], tx.exons[0][1]),
        (raw::EXON1_GENOME.0 + 1, raw::EXON1_GENOME.1 + 1),
        "genome bounds convert 0-based half-open -> 1-based [incl, excl)"
    );
    assert_eq!(
        (tx.exons[1][0], tx.exons[1][1]),
        (raw::EXON2_GENOME.0 + 1, raw::EXON2_GENOME.1 + 1),
    );
}

/// The CDS bounds ride through untouched **on a contiguous exon table**, so a
/// reader cannot infer that *everything* on the record was shifted.
///
/// The qualifier is load-bearing rather than a hedge, and it is a third
/// coordinate space rather than a third layer of the two above. cdot's
/// `start_codon`/`stop_codon` are offsets into the gap-COLLAPSED transcript —
/// exon bases only — while the exon table's own `tx_start`/`tx_end` are offsets
/// into the flat sequence. `data::cdot::collapsed_to_flat_tx_pos` maps the
/// former onto the latter (#1619), and that mapping is the **identity** exactly
/// when the exon table is contiguous and starts at 0.
///
/// [`RAW_CDOT_JSON`] is both: exon 1 covers tx `1..=100` and exon 2 `101..=150`,
/// so the internal rows abut at 100 with no hole between them. That is why
/// `start_codon` 30 and `stop_codon` 120 arrive here unchanged — not because the
/// loader leaves the CDS alone.
///
/// Give the fixture a transcript-coordinate gap and the identity stops holding:
/// the bounds would then be shifted by the enclosed hole while the exon rows
/// around them are not, which is this module's mixed-basis trap in its third
/// space. So do not read this test as "the loader never touches the CDS" — read
/// it as "on this geometry, the mapping has nothing to do".
#[test]
fn cds_bounds_survive_loading_unchanged_on_a_contiguous_exon_table() {
    let tx = load();
    assert_eq!(tx.cds_start, Some(30));
    assert_eq!(tx.cds_end, Some(120));
}

/// `coords::cdot_tx_coords` describes the RAW layer, and this pins that the raw
/// and internal layers genuinely differ — so pointing the helper at an
/// in-memory `CdotTranscript` is a real error, not a harmless no-op.
#[test]
fn the_coords_helpers_describe_the_raw_layer_not_the_internal_one() {
    let tx = load();

    // On the raw layer the helper is correctly a no-op: raw tx bounds are
    // already 1-based.
    assert_eq!(
        cdot_tx_coords(raw::EXON1_TX.0, raw::EXON1_TX.1),
        raw::EXON1_TX
    );

    // On the internal layer it is wrong, because the internal values are no
    // longer 1-based. Stated as an inequality rather than as an expected
    // wrong answer: the point is that the two layers are distinguishable.
    assert_ne!(
        cdot_tx_coords(tx.exons[0][2], tx.exons[0][3]),
        raw::EXON1_TX,
        "if these were equal the two layers would be indistinguishable and the \
         helper's 'already 1-based' claim would be safe to apply to either — it \
         is not"
    );

    // Same asymmetry on the genomic axis, in the other direction: the raw
    // helper's input is the raw value, and the loaded value has already moved.
    let (raw_g_start, _raw_g_end) =
        cdot_genomic_to_closed(raw::EXON1_GENOME.0, raw::EXON1_GENOME.1);
    assert_eq!(
        raw_g_start, tx.exons[0][0],
        "the genomic start agrees across layers (both 1-based inclusive) — it is \
         the END that differs: the helper returns 1-based inclusive, ingestion \
         stores 1-based exclusive"
    );
    assert_ne!(
        cdot_genomic_to_closed(raw::EXON1_GENOME.0, raw::EXON1_GENOME.1).1,
        tx.exons[0][1],
        "1-based inclusive end vs 1-based exclusive end must not coincide"
    );
}
