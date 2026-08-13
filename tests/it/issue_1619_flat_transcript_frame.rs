//! Issue #1619: a `c.`/`n.` position is an offset on the FLAT transcript
//! sequence, so it must not be resolved by walking the exon table.
//!
//! # The two frames, and which one a `c.` number lives in
//!
//! `docs/background/numbering.md:21` anchors `c.1` on "the **`A`** of the `ATG`
//! translation initiation (start) codon" of the *reference sequence*, and
//! `:52` numbers `n.` "from the first to the last nucleotide of the reference
//! sequence". Neither sentence mentions a genome, an alignment or an exon: for
//! a transcript accession the reference sequence **is** the spliced transcript,
//! and a position on it is a count of its own bases. `:40` says so from the
//! other side — a coding DNA reference sequence "**does not contain** intron or
//! 5' and 3' gene flanking sequences", so there is nothing in that axis for an
//! exon table to skip over.
//!
//! That is the **sequence frame**, and it is the frame `ReferenceProvider::
//! get_sequence` serves and the normalizer indexes. The **genome frame** —
//! `genomic_to_tx` / `tx_to_genomic` and everything `VariantProjector` builds on
//! them — is a different question with a different answer, and it stays
//! exon- and CIGAR-aware. This module pins both halves, because the defect was
//! one frame's arithmetic leaking into the other.
//!
//! # Why the accession's own numbering is flat, measured rather than assumed
//!
//! `NM_033517.1` (SHANK3) is one of the 58 GRCh38 cdot builds carrying a real
//! transcript-coordinate gap: exon 10 ends at tx 1302, exon 11 starts at tx
//! 1342, so 39 transcript bases align to no genomic exon. The gap is genuine
//! upstream annotation (`tests/it/normalization_transcripts_exon_contract.rs`
//! exempts it by name and size), but it is a property of cdot's **genome
//! alignment**, not of the accession's numbering:
//!
//! - `NP_277052.1`, the protein `NM_033517.1` codes for, is **1731 aa** (NCBI
//!   esummary `slen`; GenPept `LOCUS NP_277052 1731 aa`).
//! - The GenBank record for `NM_033517.1` annotates `CDS 1..5196` and
//!   `/coded_by="NM_033517.1:1..5196"`. `5196 / 3 = 1732` codons = 1731 residues
//!   plus the stop.
//! - That CDS span is contiguous over the mRNA and 39 bases longer than the
//!   exon-covered span, so RefSeq's own CDS annotation counts the unaligned
//!   bases. RefSeq's own exon table for the accession tiles `1..63`, `64..267`,
//!   … `4605..7096` with no hole at all.
//!
//! So the accession's native coordinate space is flat and contiguous, and a
//! `c.` number written against that accession is a flat offset. The fixture
//! carries the flat `cds_end: 5196`, matching GenBank.
//!
//! # What was wrong
//!
//! `CoordinateMapper::cds_to_tx` walked the exon list whenever it saw a
//! junction that did not abut, so every `c.` position 3' of the hole resolved
//! 39 bases 3' of the base it names, while `get_sequence` (and therefore the
//! normalizer) indexed the flat transcript. `hgvs_to_spdi` and the normalizer
//! then named different bases of the same accession — the denoted-sequence
//! oracle's report on `NM_033517.1:c.4818dupC`:
//!
//! ```text
//! FERRO_ASSERT_SEQUENCE: normalization changed the sequence the description denotes
//!   input: NM_033517.1:c.4818dupC   output: NM_033517.1:c.4818dup
//!   over NM_033517.1 [4877, 4877): input applies C / output applies T
//! ```

use std::collections::HashMap;
use std::sync::OnceLock;

use ferro_hgvs::convert::CoordinateMapper;
use ferro_hgvs::hgvs::location::{CdsPos, TxPos};
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};
use serde::Deserialize;

/// The gapped record this module is built on, and the flat facts about it.
const ACCESSION: &str = "NM_033517.1";
/// cdot's hole: exon 10 ends here …
const GAP_LAST_TX_BASE_BEFORE: u64 = 1302;
/// … and exon 11 resumes here, 39 transcript bases later.
const GAP_FIRST_TX_BASE_AFTER: u64 = 1342;
/// The size of that hole — what the exon walk used to add to every position 3'
/// of it.
const GAP_SIZE: u64 = GAP_FIRST_TX_BASE_AFTER - GAP_LAST_TX_BASE_BEFORE - 1;

#[derive(Debug, Deserialize)]
struct ExonData {
    number: u32,
    start: u64,
    end: u64,
}

#[derive(Debug, Deserialize)]
struct TranscriptData {
    id: String,
    gene_symbol: String,
    strand: String,
    sequence: String,
    cds_start: u64,
    cds_end: u64,
    exons: Vec<ExonData>,
}

fn fixture_records() -> &'static HashMap<String, TranscriptData> {
    static RECORDS: OnceLock<HashMap<String, TranscriptData>> = OnceLock::new();
    RECORDS.get_or_init(|| {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/sequences/normalization_transcripts.json"
        );
        let json = std::fs::read_to_string(path).expect("read normalization_transcripts.json");
        let records: Vec<TranscriptData> =
            serde_json::from_str(&json).expect("parse normalization_transcripts.json");
        records.into_iter().map(|r| (r.id.clone(), r)).collect()
    })
}

fn shank3() -> Transcript {
    let data = fixture_records()
        .get(ACCESSION)
        .unwrap_or_else(|| panic!("{ACCESSION} present in the fixture"));
    Transcript::new(
        data.id.clone(),
        Some(data.gene_symbol.clone()),
        if data.strand == "+" {
            Strand::Plus
        } else {
            Strand::Minus
        },
        data.sequence.clone(),
        Some(data.cds_start),
        Some(data.cds_end),
        data.exons
            .iter()
            .map(|e| Exon::new(e.number, e.start, e.end))
            .collect(),
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    )
}

fn shank3_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_transcript(shank3());
    provider
}

/// The flat base at a 1-based transcript coordinate, read straight out of the
/// fixture's stored sequence — the same bytes `get_sequence` serves.
fn flat_base_at(tx_base: u64) -> char {
    let data = fixture_records().get(ACCESSION).expect("record");
    data.sequence
        .as_bytes()
        .get((tx_base - 1) as usize)
        .map(|b| *b as char)
        .unwrap_or_else(|| panic!("tx {tx_base} within the stored sequence"))
}

/// The gap this module rests on is really in the fixture, at the size the
/// module doc claims. A test that silently lost its geometry would pass
/// everything below for the wrong reason.
#[test]
fn the_record_carries_the_real_cdot_gap_and_the_flat_genbank_cds() {
    let tx = shank3();
    let junction = tx
        .exons
        .windows(2)
        .find(|w| w[0].end + 1 != w[1].start)
        .expect("NM_033517.1 carries cdot's transcript-coordinate gap");
    assert_eq!(junction[0].end, GAP_LAST_TX_BASE_BEFORE);
    assert_eq!(junction[1].start, GAP_FIRST_TX_BASE_AFTER);
    assert_eq!(GAP_SIZE, 39);

    // GenBank: `CDS 1..5196`, `/coded_by="NM_033517.1:1..5196"`, and
    // NP_277052.1 is 1731 aa == (5196 / 3) - 1 residues plus the stop codon.
    assert_eq!(tx.cds_start, Some(1));
    assert_eq!(tx.cds_end, Some(5196));
    assert_eq!((5196 / 3) - 1, 1731);
    // The exon table covers 39 fewer bases than that CDS, which is the whole
    // point: cdot's alignment, not the accession's numbering.
    assert_eq!(5196 - GAP_SIZE, 5157);
}

/// The sequence frame: `c.N` is `cds_start + N - 1` on the flat transcript,
/// including for positions 3' of a gap in the genome alignment.
#[test]
fn a_coding_position_resolves_on_the_flat_transcript_axis() {
    let tx = shank3();
    let mapper = CoordinateMapper::new(&tx);
    let cds_start = tx.cds_start.expect("coding");

    // 5' of the hole — unaffected either way, so it is the control.
    for c in [1_i64, 100, 1302] {
        assert_eq!(
            mapper.cds_to_tx(&CdsPos::new(c)).unwrap().base,
            cds_start as i64 + c - 1,
            "c.{c} (5' of the gap) must be the flat offset",
        );
    }

    // 3' of the hole — the walk added GAP_SIZE here.
    for c in [1342_i64, 4818, 5196] {
        assert_eq!(
            mapper.cds_to_tx(&CdsPos::new(c)).unwrap().base,
            cds_start as i64 + c - 1,
            "c.{c} (3' of the gap) must be the flat offset, not +{GAP_SIZE}",
        );
    }

    // And the inverse, so the two directions cannot drift apart.
    for c in [1_i64, 1302, 1342, 4818, 5196] {
        let tx_pos = mapper.cds_to_tx(&CdsPos::new(c)).unwrap();
        assert_eq!(
            mapper.tx_to_cds(&tx_pos).unwrap().base,
            c,
            "round trip c.{c}"
        );
    }
}

/// `c.4818` and `n.4818` name the same base of `NM_033517.1` — `cds_start` is
/// 1, so the two spellings are the same offset — and `hgvs_to_spdi` must agree
/// with the stored sequence about which base that is. Before the fix the `c.`
/// arm walked the exon table and the `n.` arm did not, so the two disagreed by
/// exactly the gap.
#[test]
fn the_spdi_for_a_coding_position_names_the_base_the_transcript_has_there() {
    let provider = shank3_provider();
    let coding = hgvs_to_spdi(
        &parse_hgvs(&format!("{ACCESSION}:c.4818del")).expect("parse c."),
        &provider,
    )
    .expect("resolve c.4818del");
    let noncoding = hgvs_to_spdi(
        &parse_hgvs(&format!("{ACCESSION}:n.4818del")).expect("parse n."),
        &provider,
    )
    .expect("resolve n.4818del");

    assert_eq!(
        coding, noncoding,
        "c.4818 and n.4818 are the same base of {ACCESSION} (cds_start == 1)",
    );
    // Interbase, so the 1-based transcript base 4818 starts at 4817.
    assert_eq!(coding.position, 4817);
    assert_eq!(coding.deletion, flat_base_at(4818).to_string());
    assert_eq!(coding.deletion, "C");
}

/// The reproducer #1619 was filed on, expressed without the oracle flag: the
/// input and the normalized output must denote the same base, and that base is
/// the flat transcript's.
#[test]
fn the_dup_reproducer_denotes_one_base_before_and_after_normalization() {
    let provider = shank3_provider();
    let input = parse_hgvs(&format!("{ACCESSION}:c.4818dupC")).expect("parse input");
    let normalized = Normalizer::new(shank3_provider())
        .normalize(&input)
        .expect("normalize");
    assert_eq!(normalized.to_string(), format!("{ACCESSION}:c.4818dup"));

    let before = hgvs_to_spdi(&input, &provider).expect("resolve input");
    let after = hgvs_to_spdi(&normalized, &provider).expect("resolve output");
    assert_eq!(
        before.insertion, after.insertion,
        "the duplicated base must not change identity across normalization",
    );
    assert_eq!(before.insertion, "C");
    assert_eq!(before, after);
}

/// A gapped transcript that also carries genomic coordinates. The gap's 39
/// bases exist in the transcript sequence and align to nothing, which is what
/// makes it the right shape for the genome-frame half below.
fn gapped_transcript_with_genomic_coords() -> Transcript {
    // tx 1..10 -> g.1001..1010, then 10 unaligned transcript bases (tx 11..20),
    // then tx 21..30 -> g.2001..2010.
    Transcript::new(
        "NM_GENOMEFRAME.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        "A".repeat(30),
        Some(1),
        Some(30),
        vec![
            Exon::with_genomic(1, 1, 10, 1001, 1010),
            Exon::with_genomic(2, 21, 30, 2001, 2010),
        ],
        Some("chr1".to_string()),
        Some(1001),
        Some(2010),
        Default::default(),
        ManeStatus::None,
        None,
        None,
    )
}

/// The genome frame stays exon-aware. This is the half that must NOT move:
/// `genomic_to_tx` / `tx_to_genomic` read the exon table, an unaligned
/// transcript base has no genomic position, and a genomic position 3' of the
/// hole still maps to the flat transcript base the exon table gives it.
#[test]
fn the_genome_frame_remains_exon_aware() {
    let tx = gapped_transcript_with_genomic_coords();
    let mapper = CoordinateMapper::new(&tx);

    // Exon 1: genomic and transcript run together.
    assert_eq!(mapper.genomic_to_tx(1001).unwrap().unwrap().base, 1);
    assert_eq!(mapper.genomic_to_tx(1010).unwrap().unwrap().base, 10);
    // The intervening genome is intron — no transcript position at all.
    assert!(mapper.genomic_to_tx(1500).unwrap().is_none());
    // Exon 2 resumes at transcript base 21, NOT at 11: the exon table, not a
    // flat count, decides where the genome lands on the transcript.
    assert_eq!(mapper.genomic_to_tx(2001).unwrap().unwrap().base, 21);
    assert_eq!(mapper.genomic_to_tx(2010).unwrap().unwrap().base, 30);

    // And back, exon-aware in the same way.
    assert_eq!(mapper.tx_to_genomic(&TxPos::new(10)).unwrap(), Some(1010));
    assert_eq!(mapper.tx_to_genomic(&TxPos::new(21)).unwrap(), Some(2001));
    // A transcript base inside the hole aligns to nothing, and says so.
    assert_eq!(mapper.tx_to_genomic(&TxPos::new(15)).unwrap(), None);

    // Composed: g.2001 is transcript base 21, which with cds_start == 1 is
    // c.21 on the flat axis. The exon-aware step is genome<->tx; the CDS step
    // is the flat shift.
    assert_eq!(mapper.genomic_to_cds(2001).unwrap().unwrap().base, 21);
    assert_eq!(
        mapper
            .cds_to_genomic(&CdsPos::new(21))
            .unwrap()
            .expect("c.21 aligns"),
        2001
    );
    // A CDS position inside the hole is a real transcript base with no genomic
    // counterpart — `None`, not a silently shifted coordinate.
    assert_eq!(mapper.cds_to_genomic(&CdsPos::new(15)).unwrap(), None);
}
