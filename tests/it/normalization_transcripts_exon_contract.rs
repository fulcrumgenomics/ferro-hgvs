//! Structural guard for `tests/fixtures/sequences/normalization_transcripts.json`.
//!
//! That fixture is the transcript table behind `normalize_property_tests`,
//! `normalize_tests` and `idempotency_tests`: each record is loaded into a
//! `MockProvider` as a [`Transcript`], and its `exons` list is what every
//! `c.`/`n.` position in those suites is resolved through. Nothing checked that
//! the exon list was well formed, and it was not — see the contract below.
//!
//! # The `Exon` contract
//!
//! `src/reference/transcript.rs` documents `Exon.start`/`Exon.end` as
//! "position in **transcript** coordinates (1-based, inclusive)". Transcript
//! coordinates address the *spliced* transcript — introns have already been
//! removed — so the exons must **tile a contiguous run** of transcript space
//! with no holes and no overlaps: consecutive exons abut. A transcript
//! coordinate that falls between two exons is not an intronic position, it is
//! an unaddressable one.
//!
//! What the contract does **not** say is where that run begins or ends. Neither
//! endpoint is anchored — `1` at the 5' end and `sequence.len()` at the 3' end
//! are both left free, each for its own measured reason, in the two sections
//! below. The invariant is the junctions.
//!
//! That is also **almost** exactly what cdot, the upstream producer, emits: its
//! `cds_start_i`/`cds_end_i` are 1-based inclusive transcript bounds, and
//! `CdotTranscript` consumes them as such (`src/data/cdot.rs`, `from_genome_build`).
//! "Almost", because 58 of its 474,818 multi-exon builds do carry a real gap —
//! so the contract above is what cdot means, not a property every cdot record
//! has. One such record is in this fixture, exempted by name in
//! [`CDOT_GAP_JUNCTIONS`]; see the section on it below.
//!
//! # Why this file exists
//!
//! Measured against cdot-0.2.32.refseq.GRCh38 (482,519 exon-bearing builds,
//! 474,818 of them multi-exon): **58** builds carry any transcript-coordinate
//! gap at all (0.012%), the smallest such gap is **23** bases, and **zero**
//! builds have a one-base gap. The fixture used to carry a one-base gap at
//! *every* junction of *every* multi-exon record — a shape that does not occur
//! in the upstream data even once.
//!
//! The effect was not cosmetic. On `NM_000492.4` (27 exons, 26 synthetic holes)
//! a consumer that resolved `c.` positions by walking the exon list landed
//! **10 bases** 3' of the flat-transcript answer by `c.1520`, which is the
//! `hgvs_to_spdi`-vs-normalizer disagreement recorded against #1619.
//!
//! # Coverage is deliberately a bound, not an equality
//!
//! `exons.last().end == sequence.len()` is **not** asserted. Two records
//! (`NR_046285.1`, `NR_153405.1`) have real cdot exon tables that stop short of
//! the true RefSeq transcript length (2599/2619 and 4299/4337) because their 3'
//! bases do not align to the genome. That is genuine upstream data, so the
//! assertion here is the one-sided `<=`: the exon list may under-span the stored
//! sequence, but it may never address a base the sequence does not have.
//!
//! # The 5' anchor is not asserted either, and for a measured reason
//!
//! `exons[0].start == 1` looks like the natural companion to the bound above and
//! is **not** checked. An earlier revision of this file asserted it, on no
//! measurement: over cdot-0.2.32.refseq.GRCh38, **109** builds start their first
//! exon past transcript coordinate 1 — `NM_000280.5` at 8, `NM_000910.2` at 7,
//! `NM_001100631.1` at 870 — which is **1.9×** the 58-build gapped population
//! this module's whole argument rests on. A hard 5' anchor would therefore
//! reject genuine upstream data almost twice as often as the synthetic shape it
//! was meant to catch, and it is the same under-span at the other end.
//!
//! No consumer asks for it, and since #1619 fewer consumers ask for anything
//! here. `CoordinateMapper` used to choose between flat and exon-walking
//! arithmetic on the junction predicate `w[0].end + 1 != w[1].start`; it no
//! longer does, because `cds_to_tx`/`tx_to_cds` are flat sequence-axis
//! conversions that never read the exon table (`src/convert/mapper.rs`, "the
//! two frames"). What still reads it is the **genome** frame —
//! `genomic_to_tx`/`tx_to_genomic`, which scan exons for the one spanning a
//! position — and nothing under `src/` reads `exons[0].start` as an anchor.
//! Every record in this fixture does begin at 1; that is an observation about
//! the fixture, not a contract the data owes.
//!
//! # Where the exon tables come from, and the one that is not cdot's
//!
//! 28 of the 29 records carry their accession's cdot-0.2.32.refseq.GRCh38 exon
//! table exactly, in transcript order. The one exception is `NM_001127687.1`,
//! which cdot does not carry at all; it is named in [`FLATTENED_RECORDS`] and
//! [`no_unnamed_record_is_flattened`] fails if a second single-exon record ever
//! appears. That guard exists because a single-exon record satisfies every
//! junction assertion below **for free**, so a flattening is invisible here —
//! six of them were, until they were restored.
//!
//! # One record has a REAL gap, and it is exempted by name
//!
//! `NM_033517.1` is one of the 58: cdot gives it a genuine 39-base
//! transcript-coordinate hole between exon 10 (ends 1302) and exon 11 (starts
//! 1342). That is upstream annotation, not a fixture defect, so it is pinned in
//! [`CDOT_GAP_JUNCTIONS`] — **by accession, junction and exact size** — rather
//! than accommodated by loosening the contiguity rule. A general relaxation
//! would re-admit the one-base synthetic holes this whole file exists to keep
//! out; an exemption that names its one instance and its measurement cannot.
//!
//! Carrying it is also what closed #1619, and it is the reason this record was
//! nearly left flattened. While `hgvs_to_spdi` resolved a `c.` position by
//! walking the exon list across the hole and the normalizer indexed the flat
//! transcript, the existing case `NM_033517.1:c.4818dupC` -> `c.4818dup` fired
//! under `FERRO_ASSERT_SEQUENCE=1`: the input applied `C`, the output `T`, at
//! transcript position 4877. That was the reproducer #1619 had not otherwise
//! had — it is why the issue could be settled against real cdot geometry
//! instead of a shape nobody had measured — and it is now the regression guard
//! for the ruling that closed it
//! (`c-and-n-positions-are-flat-transcript-offsets`, in
//! `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`). Flattening
//! this record would delete that guard, which is the second reason not to,
//! independent of the contract above. `tests/it/issue_1619_flat_transcript_frame.rs`
//! pins the frame directly.

use serde::Deserialize;
use std::collections::BTreeSet;

/// Records whose `exons` list is deliberately **not** the cdot table for that
/// accession, with the reason.
///
/// A single-exon record has no junction, so it satisfies [`check_exon_contract`]
/// vacuously — which is the whole hazard, and why the one that remains is
/// written down rather than left to be rediscovered.
const FLATTENED_RECORDS: &[(&str, &str)] = &[(
    "NM_001127687.1",
    "absent from cdot-0.2.32.refseq.GRCh38 entirely — there is no upstream exon \
     table to restore, so its single exon is a placeholder rather than a \
     flattening of anything",
)];

/// Junctions permitted to carry a transcript-coordinate gap, because cdot
/// itself puts one there.
///
/// Each entry is `(accession, exon number on the 5' side of the junction, exact
/// gap in bases)` and every field is checked, so this cannot drift into a
/// blanket: a different junction, a different size, or a gap on any accession
/// not listed is still a violation. Measured against
/// cdot-0.2.32.refseq.GRCh38, where 58 of 474,818 multi-exon builds carry a
/// gap, the smallest is 23 bases and **none** is one base.
///
/// `NM_033517.1` is this fixture's only member of that population. See the
/// module docs for what carrying it bought (#1619) and why it must stay.
const CDOT_GAP_JUNCTIONS: &[(&str, u32, i64)] = &[("NM_033517.1", 10, 39)];

#[derive(Debug, Deserialize)]
struct TranscriptRecord {
    id: String,
    sequence: String,
    exons: Vec<ExonRecord>,
}

#[derive(Debug, Deserialize)]
struct ExonRecord {
    number: u32,
    start: u64,
    end: u64,
}

fn load_records() -> Vec<TranscriptRecord> {
    let path = format!(
        "{}/tests/fixtures/sequences/normalization_transcripts.json",
        env!("CARGO_MANIFEST_DIR")
    );
    let json =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    serde_json::from_str(&json).unwrap_or_else(|e| panic!("failed to parse {path}: {e}"))
}

/// The fixture must actually be able to *fail* the contiguity assertion below.
///
/// This is the non-vacuity guard, and it is not ceremony. Contiguity is a
/// property of exon **junctions**, so a single-exon record has nothing to get
/// wrong and passes for free. When this fixture was first audited, the seven
/// records that satisfied contiguity were precisely its seven single-exon ones —
/// every record with a junction violated it, so "7 records pass" was a
/// structural artifact rather than evidence of anything. Six of those seven
/// have since had their real cdot tables restored; the seventh is named in
/// [`FLATTENED_RECORDS`] because cdot has no table for it.
///
/// Asserting on the junction count keeps that from recurring silently: if the
/// fixture is ever emptied, truncated, or flattened to one exon per record,
/// this fails loudly instead of letting `check_exon_contract` sail through with
/// nothing to check.
#[test]
fn fixture_exercises_the_contiguity_property() {
    let records = load_records();
    assert!(
        records.len() >= 29,
        "fixture shrank to {} records; it had 29 when this guard was written",
        records.len()
    );

    let junctions: usize = records
        .iter()
        .map(|r| r.exons.len().saturating_sub(1))
        .sum();
    let multi_exon = records.iter().filter(|r| r.exons.len() > 1).count();
    assert!(
        multi_exon >= 28 && junctions >= 517,
        "contiguity is a property of exon junctions, so it is only checkable on \
         multi-exon records: found {multi_exon} multi-exon record(s) spanning \
         {junctions} junction(s), expected at least 28 and 517. A pass from \
         `check_exon_contract` with fewer than this is a structural zero, not a \
         clean fixture."
    );
}

/// Every pinned gap exemption must still describe a gap the fixture has.
///
/// [`CDOT_GAP_JUNCTIONS`] suppresses a violation, so a stale entry is the one
/// kind of error [`check_exon_contract`] cannot report: if the gap it names were
/// edited away, the exemption would go quietly dead and the guard would still
/// pass. Requiring each entry to match a junction that is really there keeps the
/// list **shrink-only** — restoring contiguity means deleting the row, and
/// deleting the row is a visible diff.
#[test]
fn every_pinned_gap_exemption_is_still_load_bearing() {
    let records = load_records();
    for &(accession, number, size) in CDOT_GAP_JUNCTIONS {
        let record = records
            .iter()
            .find(|r| r.id == accession)
            .unwrap_or_else(|| panic!("{accession} is pinned in CDOT_GAP_JUNCTIONS but absent"));
        let found = record.exons.windows(2).any(|pair| {
            pair[0].number == number
                && pair[1].start as i64 - pair[0].end as i64 - 1 == size
                && size != 0
        });
        assert!(
            found,
            "CDOT_GAP_JUNCTIONS pins a {size}-base gap after exon {number} of {accession}, but \
             the fixture has no such junction. If the table was restored to contiguity, delete \
             the row — an exemption that matches nothing silently weakens this guard."
        );
    }
}

/// No record may be flattened to a single exon without saying so.
///
/// [`check_exon_contract`] cannot see a flattening — a single-exon record has no
/// junction to violate — so a real multi-exon transcript rewritten as one exon
/// passes it silently, and 111 exons across six records did exactly that before
/// [`FLATTENED_RECORDS`] existed. This test makes the set of single-exon records
/// an exact match against that list, so adding one costs a line and a stated
/// reason, and restoring one costs deleting its line.
#[test]
fn no_unnamed_record_is_flattened() {
    let records = load_records();
    let single_exon: BTreeSet<&str> = records
        .iter()
        .filter(|r| r.exons.len() == 1)
        .map(|r| r.id.as_str())
        .collect();
    let named: BTreeSet<&str> = FLATTENED_RECORDS.iter().map(|(id, _)| *id).collect();

    let roster = FLATTENED_RECORDS
        .iter()
        .map(|(id, reason)| format!("  {id}: {reason}"))
        .collect::<Vec<_>>()
        .join("\n");
    assert_eq!(
        single_exon,
        named,
        "the fixture's single-exon records must be exactly those named in \
         FLATTENED_RECORDS. Unnamed and single-exon: {:?}. Named but no longer \
         single-exon: {:?}. A single-exon record passes `check_exon_contract` \
         vacuously, so an unnamed one is coverage this file does not have.\n\
         FLATTENED_RECORDS is:\n{roster}",
        single_exon.difference(&named).collect::<Vec<_>>(),
        named.difference(&single_exon).collect::<Vec<_>>(),
    );
}

/// Every record must satisfy the `Exon` contract stated in the module docs.
#[test]
fn check_exon_contract() {
    let records = load_records();
    assert!(!records.is_empty(), "fixture is empty");

    let mut violations: Vec<String> = Vec::new();
    for record in &records {
        let id = &record.id;
        let exons = &record.exons;

        if exons.is_empty() {
            violations.push(format!("{id}: has no exons"));
            continue;
        }

        // Exon numbers are 1-based and consecutive in transcript order.
        for (index, exon) in exons.iter().enumerate() {
            let expected = index as u32 + 1;
            if exon.number != expected {
                violations.push(format!(
                    "{id}: exon at index {index} is numbered {}, expected {expected}",
                    exon.number
                ));
            }
            // 1-based inclusive bounds describing at least one base.
            if exon.start == 0 {
                violations.push(format!("{id}: exon {expected} starts at 0 (1-based)"));
            }
            if exon.end < exon.start {
                violations.push(format!(
                    "{id}: exon {expected} is inverted ({}..={})",
                    exon.start, exon.end
                ));
            }
        }

        // Deliberately NOT asserted here: `exons[0].start == 1`. 109 real cdot
        // GRCh38 builds begin past transcript coordinate 1 and no consumer reads
        // the field — see "The 5' anchor is not asserted either" in the module
        // docs for the measurement and the reasoning.

        // The load-bearing invariant: exons tile transcript space with no holes
        // and no overlaps, because transcript coordinates address the spliced
        // sequence.
        for pair in exons.windows(2) {
            let (left, right) = (&pair[0], &pair[1]);
            if right.start != left.end + 1 {
                let gap = right.start as i64 - left.end as i64 - 1;
                // A gap cdot itself publishes here, of exactly this size, is
                // upstream annotation rather than a defect. Every field must
                // match — see `CDOT_GAP_JUNCTIONS`.
                if CDOT_GAP_JUNCTIONS
                    .iter()
                    .any(|&(acc, number, size)| acc == id && number == left.number && size == gap)
                {
                    continue;
                }
                violations.push(format!(
                    "{id}: exon {} ends at {} but exon {} starts at {} ({} of {} transcript \
                     base(s)) — transcript coordinates must tile the spliced sequence",
                    left.number,
                    left.end,
                    right.number,
                    right.start,
                    if gap > 0 { "hole" } else { "overlap" },
                    gap.abs(),
                ));
            }
        }

        // The exon list may under-span the stored sequence (see module docs) but
        // must never address a base beyond it.
        let last_end = exons[exons.len() - 1].end;
        if last_end > record.sequence.len() as u64 {
            violations.push(format!(
                "{id}: last exon ends at transcript position {last_end}, past the {}-base sequence",
                record.sequence.len()
            ));
        }
    }

    assert!(
        violations.is_empty(),
        "{} record(s) of {} violate the Exon contract:\n  {}",
        violations
            .iter()
            .filter_map(|v| v.split(':').next())
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
        records.len(),
        violations.join("\n  ")
    );
}
