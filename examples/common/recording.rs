//! A [`ReferenceProvider`] that wraps a real manifest-backed provider and
//! records every reference access made through it, plus the window-assembly
//! that turns those recordings into a committed
//! [`WindowFixture`](ferro_hgvs::conformance::reference_window::WindowFixture).
//!
//! This is how a hermetic fixture is captured without guessing what a pass
//! needs: run the pass for real against the manifest, record exactly the
//! transcripts, base ranges and genomic-parent placements it touched, and
//! commit those (padded and merged) instead of the multi-GB reference.
//!
//! Mirrors the recorder embedded in `extract_mutalyzer_genomic_windows.rs`,
//! lifted here so a second extractor does not restate it.

#![allow(dead_code)]

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;
use std::sync::Arc;

use ferro_hgvs::conformance::reference_window::{GenomicWindow, WindowFixture};
use ferro_hgvs::hgvs::variant::Accession;
use ferro_hgvs::reference::derived_placement::DerivedPlacement;
use ferro_hgvs::reference::provider::GenomicPlacement;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{FerroError, HgvsVariant, MultiFastaProvider, ReferenceProvider};

/// Safety margin (bp) added to each side of a captured window beyond the
/// realized access range. Generous relative to any HGVS shuffle distance.
pub const PAD: u64 = 64;
/// Recorded ranges on one accession closer than this gap are merged into a
/// single window (keeps the fixture small while still splitting distant loci).
pub const MERGE_GAP: u64 = 2 * PAD;

/// Everything the wrapped provider was asked for during a pass.
#[derive(Default)]
pub struct Recorder {
    /// Every transcript resolved during the pass, keyed by id (deduped).
    pub transcripts: BTreeMap<String, Transcript>,
    /// Recorded genomic ranges per contig, 0-based half-open.
    pub contig_ranges: BTreeMap<String, Vec<(u64, u64)>>,
    /// `get_sequence` ranges per id; an id that is not a resolved transcript
    /// becomes its own window accession.
    pub seq_ranges: Vec<(String, u64, u64)>,
    /// Every `NG_`/`LRG_` placement consulted, keyed by parent accession.
    pub placements: BTreeMap<String, GenomicPlacement>,
}

/// Wraps the real manifest provider, forwarding every call and recording what
/// the pass touches. `Clone` shares the recorder (`Rc<RefCell<_>>`) so a fresh
/// clone handed to each per-case projector accumulates into one record.
#[derive(Clone)]
pub struct RecordingProvider {
    inner: Arc<MultiFastaProvider>,
    rec: Rc<RefCell<Recorder>>,
}

impl RecordingProvider {
    pub fn new(inner: Arc<MultiFastaProvider>) -> Self {
        Self {
            inner,
            rec: Rc::new(RefCell::new(Recorder::default())),
        }
    }

    /// The wrapped manifest provider (used to read captured windows back out).
    pub fn inner(&self) -> &MultiFastaProvider {
        &self.inner
    }

    fn capture(&self, tx: &Arc<Transcript>) {
        self.rec
            .borrow_mut()
            .transcripts
            .entry(tx.id.clone())
            .or_insert_with(|| (**tx).clone());
    }
}

impl ReferenceProvider for RecordingProvider {
    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError> {
        let tx = self.inner.get_transcript(id)?;
        self.capture(&tx);
        Ok(tx)
    }

    fn get_transcript_for_variant(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Arc<Transcript>, FerroError> {
        let tx = self.inner.get_transcript_for_variant(variant)?;
        self.capture(&tx);
        Ok(tx)
    }

    fn get_transcript_for_accession(
        &self,
        accession: &Accession,
    ) -> Result<Arc<Transcript>, FerroError> {
        let tx = self.inner.get_transcript_for_accession(accession)?;
        self.capture(&tx);
        Ok(tx)
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        let out = self.inner.get_sequence(id, start, end)?;
        self.rec
            .borrow_mut()
            .seq_ranges
            .push((id.to_string(), start, end));
        Ok(out)
    }

    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        let out = self.inner.get_genomic_sequence(contig, start, end)?;
        self.rec
            .borrow_mut()
            .contig_ranges
            .entry(contig.to_string())
            .or_default()
            .push((start, end));
        Ok(out)
    }

    fn genomic_placement(&self, parent: &Accession) -> Option<GenomicPlacement> {
        let p = self.inner.genomic_placement(parent)?;
        self.rec
            .borrow_mut()
            .placements
            .entry(parent.full())
            .or_insert_with(|| p.clone());
        Some(p)
    }

    fn genomic_placement_on_build(
        &self,
        parent: &Accession,
        build: Option<&str>,
    ) -> Option<GenomicPlacement> {
        let p = self.inner.genomic_placement_on_build(parent, build)?;
        self.rec
            .borrow_mut()
            .placements
            .entry(parent.full())
            .or_insert_with(|| p.clone());
        Some(p)
    }

    fn resolve_legacy_gene_selector(
        &self,
        selector: &str,
        ng_parent: Option<&Accession>,
    ) -> Option<String> {
        self.inner.resolve_legacy_gene_selector(selector, ng_parent)
    }

    fn has_transcript(&self, id: &str) -> bool {
        self.inner.has_transcript(id)
    }

    fn has_transcript_version_exact(&self, id: &str) -> bool {
        self.inner.has_transcript_version_exact(id)
    }

    fn has_genomic_data(&self) -> bool {
        self.inner.has_genomic_data()
    }

    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        self.inner.get_protein_sequence(accession, start, end)
    }

    fn has_protein_data(&self) -> bool {
        self.inner.has_protein_data()
    }

    fn get_sequence_length(&self, id: &str) -> Result<u64, FerroError> {
        self.inner.get_sequence_length(id)
    }
}

/// Merge `(start, end)` ranges that overlap or sit within [`MERGE_GAP`] into
/// padded, disjoint intervals (each clamped to `[0, len]`).
pub fn merged_windows(mut ranges: Vec<(u64, u64)>, len: Option<u64>) -> Vec<(u64, u64)> {
    ranges.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::new();
    for (s, e) in ranges {
        let s = s.saturating_sub(PAD);
        let e = match len {
            Some(l) => (e + PAD).min(l),
            None => e + PAD,
        };
        match out.last_mut() {
            Some(last) if s <= last.1 + MERGE_GAP => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

fn strand_str(s: Strand) -> &'static str {
    match s {
        Strand::Plus => "+",
        Strand::Minus => "-",
        _ => "?",
    }
}

/// Assemble the committed fixture from a completed recording.
///
/// `keep_transcript` decides which resolved transcripts are kept: an
/// `NG_`/`LRG_` resolved as a "transcript" carries a multi-hundred-KB genomic
/// sequence and is served through its placement instead, so extractors filter
/// those out to keep the committed file small.
pub fn build_window_fixture(
    provider: &RecordingProvider,
    captured_from: String,
    description: String,
    keep_transcript: impl Fn(&Transcript) -> bool,
) -> WindowFixture {
    let rec = provider.rec.borrow();

    // Cluster every read range by accession. A `get_sequence` id that is not a
    // resolved transcript is itself a windowed accession (a chromosome read via
    // `get_sequence`, or a synthesized `NG_`).
    let mut acc_ranges: BTreeMap<String, Vec<(u64, u64)>> = rec.contig_ranges.clone();
    for (id, start, end) in &rec.seq_ranges {
        if !rec.transcripts.contains_key(id) {
            acc_ranges
                .entry(id.clone())
                .or_default()
                .push((*start, *end));
        }
    }

    let transcripts: Vec<Transcript> = rec
        .transcripts
        .values()
        .filter(|tx| keep_transcript(tx))
        .cloned()
        .collect();

    let mut contig_lengths: BTreeMap<String, u64> = BTreeMap::new();
    let mut genomic: Vec<GenomicWindow> = Vec::new();
    for (acc, ranges) in &acc_ranges {
        let len = provider.inner.get_sequence_length(acc).ok();
        if let Some(l) = len {
            contig_lengths.insert(acc.clone(), l);
        }
        for (start, end) in merged_windows(ranges.clone(), len) {
            let bases = provider
                .inner
                .get_sequence(acc, start, end)
                .unwrap_or_else(|e| panic!("capturing window {acc}:{start}-{end}: {e}"));
            genomic.push(GenomicWindow {
                contig: acc.clone(),
                start,
                bases,
            });
        }
    }
    genomic.sort_by(|a, b| (a.contig.as_str(), a.start).cmp(&(b.contig.as_str(), b.start)));

    let mut placements: Vec<DerivedPlacement> = rec
        .placements
        .iter()
        .map(|(parent, p)| DerivedPlacement {
            parent: parent.clone(),
            nc: p.nc.full(),
            nc_start: p.nc_start,
            nc_end: p.nc_end,
            strand: strand_str(p.strand).to_string(),
            anchored_by: String::new(),
            mismatch_fraction: 0.0,
        })
        .collect();
    placements.sort_by(|a, b| a.parent.cmp(&b.parent));

    WindowFixture {
        description,
        captured_from,
        contig_lengths,
        transcripts,
        genomic,
        placements,
    }
}

#[cfg(test)]
mod tests {
    use super::{merged_windows, MERGE_GAP, PAD};

    // Guard the arithmetic below against a future retuning of the constants.
    const _: () = assert!(PAD == 64 && MERGE_GAP == 128);

    #[test]
    fn start_padding_saturates_at_zero() {
        // 10 - PAD underflows u64 → clamps to 0; end grows by PAD.
        assert_eq!(merged_windows(vec![(10, 20)], None), vec![(0, 20 + PAD)]);
    }

    #[test]
    fn end_padding_clamps_to_len() {
        // 20 + PAD = 84 would exceed len 50 → clamped to len.
        assert_eq!(merged_windows(vec![(10, 20)], Some(50)), vec![(0, 50)]);
    }

    #[test]
    fn overlapping_padded_ranges_merge() {
        // (100,110)→(36,174) and (150,160)→(86,224) overlap → one interval.
        assert_eq!(
            merged_windows(vec![(100, 110), (150, 160)], None),
            vec![(36, 224)]
        );
    }

    #[test]
    fn ranges_within_merge_gap_merge() {
        // (100,110)→(36,174) and (300,310)→(236,374): a 62 bp gap ≤ MERGE_GAP
        // (128) → still merged into one interval.
        assert_eq!(
            merged_windows(vec![(100, 110), (300, 310)], None),
            vec![(36, 374)]
        );
    }

    #[test]
    fn disjoint_ranges_stay_separate() {
        // (100,110)→(36,174) and (1000,1010)→(936,1074): the 762 bp gap exceeds
        // MERGE_GAP → two intervals.
        assert_eq!(
            merged_windows(vec![(100, 110), (1000, 1010)], None),
            vec![(36, 174), (936, 1074)]
        );
    }

    #[test]
    fn input_is_sorted_before_merging() {
        // Reversed input yields the same merged result as the sorted order.
        assert_eq!(
            merged_windows(vec![(150, 160), (100, 110)], None),
            vec![(36, 224)]
        );
    }
}
