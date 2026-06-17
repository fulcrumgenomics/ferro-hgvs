//! Hermetic, fixture-backed reference windows for the biocommons conformance
//! gate (#325, #478 pillar 4).
//!
//! The manifest-mode `axis_normalized` test needs a multi-GB reference manifest
//! that exists only on a few machines, so CI skips it and the merge gate never
//! runs. This module is the fixture half of the fix: a [`WindowFixture`] is a
//! small, committed snapshot of *exactly* the reference bases the biocommons
//! corpus's normalize pass touches — every transcript it resolves (stored whole
//! via [`Transcript`]'s own serde) plus a window of each genomic contig it reads
//! — and [`WindowProvider`] serves them through the [`ReferenceProvider`] trait
//! with no out-of-band data.
//!
//! The fixture is *generated* by `examples/extract_biocommons_windows.rs` (with
//! a `--check` mode, like the spec-fixture and conformance-summary generators):
//! it wraps the real manifest provider in a recording shim, runs the same
//! normalize loop the test does, and captures whatever bases were requested. So
//! a window is never hand-sized — it is precisely the access set of a real run,
//! padded by a safety margin. Serving is then a pure re-base + slice, with no
//! assumptions about the genomic coordinate convention: the extractor captured
//! each window's bases by asking the manifest provider for the same
//! `[start, end)` range the normalizer requested, so byte `i` of `bases` is the
//! base at global coordinate `start + i` by construction.

use std::collections::{BTreeMap, HashMap};
use std::path::Path;
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::reference::derived_placement::{DerivedPlacement, DerivedPlacements};
use crate::reference::provider::{select_placement_for_build, GenomicPlacement, ReferenceProvider};
use crate::reference::transcript::Transcript;
use crate::FerroError;

/// A captured slice of one genomic contig: the bases the normalize pass read,
/// keyed by the global coordinate of the first base so serve-time lookups can
/// re-base into the window. `start` and the indices of `bases` use the exact
/// convention the [`ReferenceProvider::get_genomic_sequence`] caller used when
/// the extractor recorded the access (0-based, half-open), so serving is a pure
/// subtraction and slice.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenomicWindow {
    /// Contig/chromosome accession (e.g. `"NC_000006.11"`).
    pub contig: String,
    /// Global coordinate of `bases[0]` (0-based), in the provider's
    /// genomic-coordinate convention.
    pub start: u64,
    /// The captured bases for `[start, start + bases.len())`.
    pub bases: String,
}

impl GenomicWindow {
    /// Global end coordinate (exclusive) covered by this window.
    fn end(&self) -> u64 {
        self.start + self.bases.len() as u64
    }

    /// Bases for `[start, end)` if fully contained, else `None`.
    fn slice(&self, start: u64, end: u64) -> Option<String> {
        if start < self.start || end > self.end() || start > end {
            return None;
        }
        let lo = (start - self.start) as usize;
        let hi = (end - self.start) as usize;
        Some(self.bases[lo..hi].to_string())
    }
}

/// The committed hermetic fixture: every transcript and genomic window the
/// biocommons corpus's normalize pass touches, captured from a manifest run by
/// `examples/extract_biocommons_windows.rs`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WindowFixture {
    /// Human-facing provenance note (schema version, generator command).
    #[serde(default)]
    pub description: String,
    /// Manifest `prepared_at` the windows were captured from, for traceability.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub captured_from: String,
    /// True length of each windowed contig, captured from the manifest. Lets
    /// [`WindowProvider`] reproduce `MultiFastaProvider`'s clamping of an
    /// over-range read to the contig end (the `CanonicalSplitSkipped` / W5003
    /// path for variants that exceed the reference) instead of erroring, while
    /// still erroring on a read that is *within* the contig but outside a
    /// captured window — a genuine extraction defect, not a clamp.
    #[serde(default)]
    pub contig_lengths: BTreeMap<String, u64>,
    /// Every transcript resolved during the pass, stored whole.
    pub transcripts: Vec<Transcript>,
    /// Per-contig genomic windows (one or more disjoint slices per contig).
    #[serde(default)]
    pub genomic: Vec<GenomicWindow>,
    /// Version-independent `NG_`/`LRG_` → chromosome placements (#728) the pass
    /// used to re-anchor a transcript-coordinate variant into its genomic
    /// parent's own frame (#480). Empty for corpora that never reference an
    /// `NG_`/`LRG_` parent (e.g. biocommons); the mutalyzer genomic gate
    /// populates it so [`WindowProvider`] can serve `genomic_placement`.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub placements: Vec<DerivedPlacement>,
}

impl WindowFixture {
    /// Load a fixture from a JSON file.
    pub fn from_json_path(path: &Path) -> Result<Self, FerroError> {
        let content = std::fs::read_to_string(path)?;
        let fixture: WindowFixture = serde_json::from_str(&content)?;
        Ok(fixture)
    }

    /// Serialize to pretty JSON with a trailing newline (matches the other
    /// generated fixtures, so `--check` byte-compares cleanly).
    pub fn to_json(&self) -> Result<String, FerroError> {
        let mut s = serde_json::to_string_pretty(self)?;
        s.push('\n');
        Ok(s)
    }

    /// Build a [`WindowProvider`] serving this fixture's transcripts and windows.
    pub fn to_provider(&self) -> WindowProvider {
        let transcripts: HashMap<String, Arc<Transcript>> = self
            .transcripts
            .iter()
            .map(|tx| (tx.id.clone(), Arc::new(tx.clone())))
            .collect();
        let mut genomic: HashMap<String, Vec<GenomicWindow>> = HashMap::new();
        for w in &self.genomic {
            genomic.entry(w.contig.clone()).or_default().push(w.clone());
        }
        // Resolve the serialized `DerivedPlacement` records into the runtime
        // `GenomicPlacement` map, grouped per parent accession (one entry per
        // genome build, mirroring `MultiFastaProvider::refseqgene_placements`).
        // A record with an unparseable `nc` or a non-`+`/`-` strand is dropped
        // by `to_placements` rather than silently mis-placed.
        let mut placements: HashMap<String, Vec<GenomicPlacement>> = HashMap::new();
        let resolved = DerivedPlacements {
            description: String::new(),
            placements: self.placements.clone(),
        };
        for (parent, placement) in resolved.to_placements() {
            placements.entry(parent).or_default().push(placement);
        }
        WindowProvider {
            transcripts,
            genomic,
            contig_lengths: self
                .contig_lengths
                .iter()
                .map(|(k, v)| (k.clone(), *v))
                .collect(),
            placements,
        }
    }
}

/// A [`ReferenceProvider`] backed entirely by committed [`WindowFixture`] data —
/// no manifest, no out-of-band files. Drives the per-PR hermetic conformance
/// gate. Transcript and genomic dispatch mirror
/// [`MockProvider`](crate::reference::mock::MockProvider): a contig is anything
/// registered in `genomic` or named by a transcript's `chromosome`, and a bare
/// (unversioned) id falls back to a base-accession match.
#[derive(Debug, Clone)]
pub struct WindowProvider {
    transcripts: HashMap<String, Arc<Transcript>>,
    genomic: HashMap<String, Vec<GenomicWindow>>,
    contig_lengths: HashMap<String, u64>,
    /// `NG_`/`LRG_` parent → its chromosome placement(s), one per genome build.
    /// Empty unless the fixture carried `placements` (the mutalyzer genomic
    /// gate). Served via [`ReferenceProvider::genomic_placement_on_build`].
    placements: HashMap<String, Vec<GenomicPlacement>>,
}

impl WindowProvider {
    /// Whether `id` names a genomic contig (a `genomic` key or a transcript's
    /// `chromosome`) rather than a transcript accession. Mirrors
    /// `MockProvider::is_known_contig` so `get_sequence` dispatches identically.
    fn is_known_contig(&self, id: &str) -> bool {
        if self.genomic.contains_key(id) {
            return true;
        }
        self.transcripts
            .values()
            .any(|tx| tx.chromosome.as_deref() == Some(id))
    }
}

impl ReferenceProvider for WindowProvider {
    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError> {
        if let Some(tx) = self.transcripts.get(id) {
            return Ok(Arc::clone(tx));
        }
        // Unversioned fallback only (a versioned miss must not cross-version
        // substitute) — same rule as MockProvider::get_transcript (#311).
        if !id.contains('.') {
            for (key, tx) in &self.transcripts {
                if key.split('.').next().unwrap_or(key) == id {
                    return Ok(Arc::clone(tx));
                }
            }
        }
        Err(FerroError::ReferenceNotFound { id: id.to_string() })
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        // Contigs resolve strictly against the genomic path; otherwise treat
        // the id as a transcript accession (mirrors MockProvider dispatch).
        if self.is_known_contig(id) {
            return self.get_genomic_sequence(id, start, end);
        }
        let transcript = self.get_transcript(id)?;
        transcript
            .get_sequence(start, end)
            .map(|s| s.to_string())
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: format!(
                    "transcript {id} has no bases for {start}-{end} in the window fixture"
                ),
            })
    }

    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        let windows =
            self.genomic
                .get(contig)
                .ok_or_else(|| FerroError::GenomicReferenceNotAvailable {
                    contig: contig.to_string(),
                    start,
                    end,
                })?;
        // Reproduce MultiFastaProvider: a read whose end runs past the contig is
        // clamped to the contig end and returns the available prefix (this is
        // what lets normalize see a short read and emit CanonicalSplitSkipped /
        // W5003 for a variant that exceeds the reference). Clamping is keyed on
        // the *true* contig length, so a read that ends within the contig but
        // outside a captured window still errors as an extraction defect.
        let end = match self.contig_lengths.get(contig) {
            Some(&len) => end.min(len),
            None => end,
        };
        for w in windows {
            if let Some(s) = w.slice(start, end) {
                return Ok(s);
            }
        }
        // A request inside the contig but outside every captured window means
        // the fixture is missing bases the normalizer needs — a clear
        // extraction defect, not a missing-data skip.
        Err(FerroError::InvalidCoordinates {
            msg: format!(
                "genomic request {contig}:{start}-{end} is not covered by any committed window \
                 (regenerate reference-windows.json via \
                 `cargo run --features dev --example extract_biocommons_windows`)"
            ),
        })
    }

    fn has_genomic_data(&self) -> bool {
        !self.genomic.is_empty()
    }

    fn genomic_placement(
        &self,
        parent: &crate::hgvs::variant::Accession,
    ) -> Option<GenomicPlacement> {
        self.genomic_placement_on_build(parent, None)
    }

    fn genomic_placement_on_build(
        &self,
        parent: &crate::hgvs::variant::Accession,
        build: Option<&str>,
    ) -> Option<GenomicPlacement> {
        // Mirror `MultiFastaProvider`: select the placement matching the
        // resolved build (or GRCh38-preferred with no hint), declining rather
        // than mis-anchoring onto another build's placement.
        let list = self.placements.get(&parent.full())?;
        select_placement_for_build(list, build)
    }

    fn get_sequence_length(&self, id: &str) -> Result<u64, FerroError> {
        // Captured true contig length is authoritative for windowed contigs.
        if let Some(&len) = self.contig_lengths.get(id) {
            return Ok(len);
        }
        // Resolve through `get_transcript` so the unversioned base-id fallback
        // (NM_TEST → NM_TEST.1) applies consistently with the transcript path;
        // a raw map `get` here would spuriously miss an id `get_transcript`
        // resolves. An `Err` (no such transcript) falls through to the
        // genomic-window and error branches below.
        if let Ok(tx) = self.get_transcript(id) {
            if let Some(seq) = &tx.sequence {
                return Ok(seq.len() as u64);
            }
        }
        // Fallback for a windowed contig with no captured length: the highest
        // captured coordinate keeps a bound-probing caller inside the window.
        if let Some(windows) = self.genomic.get(id) {
            if let Some(max_end) = windows.iter().map(GenomicWindow::end).max() {
                return Ok(max_end);
            }
        }
        Err(FerroError::ReferenceNotFound { id: id.to_string() })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::{Exon, ManeStatus, Strand};

    fn tx(id: &str, seq: &str) -> Transcript {
        let len = seq.len() as u64;
        Transcript::new(
            id.to_string(),
            Some("GENE".to_string()),
            Strand::Plus,
            seq.to_string(),
            Some(1),
            Some(len),
            vec![Exon::new(1, 1, len)],
            Some("NC_000001.11".to_string()),
            None,
            None,
            Default::default(),
            ManeStatus::None,
            None,
            None,
        )
    }

    fn fixture() -> WindowFixture {
        WindowFixture {
            description: "test".to_string(),
            captured_from: String::new(),
            // NC_000006.12: window reaches the contig end (1012) — over-range
            // reads clamp. NC_DEFECT.1: contig is longer (5000) than the window
            // (ends at 1012) — a read inside the contig but past the window is
            // an extraction defect and must error, not clamp.
            contig_lengths: [
                ("NC_000006.12".to_string(), 1012u64),
                ("NC_DEFECT.1".to_string(), 5000u64),
            ]
            .into_iter()
            .collect(),
            transcripts: vec![tx("NM_TEST.1", "ACGTACGTACGT")],
            genomic: vec![
                GenomicWindow {
                    contig: "NC_000006.12".to_string(),
                    start: 1000,
                    bases: "AAACCCGGGTTT".to_string(),
                },
                GenomicWindow {
                    contig: "NC_DEFECT.1".to_string(),
                    start: 1000,
                    bases: "AAACCCGGGTTT".to_string(),
                },
            ],
            placements: Vec::new(),
        }
    }

    #[test]
    fn serves_ng_placement_from_fixture_and_round_trips() {
        use crate::hgvs::variant::Accession;
        use crate::reference::Strand;

        let mut f = fixture();
        f.placements = vec![DerivedPlacement {
            parent: "NG_012337.1".to_string(),
            nc: "NC_000011.10".to_string(),
            nc_start: 112_081_847,
            nc_end: 112_097_794,
            strand: "+".to_string(),
            anchored_by: String::new(),
            mismatch_fraction: 0.0,
        }];

        // The placements field survives a JSON round trip (committed-fixture path).
        let reloaded: WindowFixture = serde_json::from_str(&f.to_json().unwrap()).unwrap();
        assert_eq!(reloaded.placements, f.placements);

        let p = reloaded.to_provider();
        let ng = Accession::new("NG", "012337", Some(1));
        let placement = p
            .genomic_placement(&ng)
            .expect("WindowProvider serves the fixture's NG_ placement");
        assert_eq!(placement.nc.full(), "NC_000011.10");
        assert_eq!(placement.parent_start, 1);
        assert_eq!(placement.nc_start, 112_081_847);
        assert_eq!(placement.nc_end, 112_097_794);
        assert_eq!(placement.strand, Strand::Plus);

        // An unplaced parent declines.
        assert!(p
            .genomic_placement(&Accession::new("NG", "999999", Some(9)))
            .is_none());
    }

    #[test]
    fn transcript_exact_and_unversioned_lookup() {
        let p = fixture().to_provider();
        assert!(p.get_transcript("NM_TEST.1").is_ok());
        // unversioned base-id fallback
        assert!(p.get_transcript("NM_TEST").is_ok());
        // a versioned miss must NOT cross-version substitute
        assert!(p.get_transcript("NM_TEST.2").is_err());
        assert!(p.get_transcript("NM_OTHER.1").is_err());
    }

    #[test]
    fn genomic_window_rebase_and_slice() {
        let p = fixture().to_provider();
        // full window
        assert_eq!(
            p.get_genomic_sequence("NC_000006.12", 1000, 1012).unwrap(),
            "AAACCCGGGTTT"
        );
        // interior slice re-bases by start
        assert_eq!(
            p.get_genomic_sequence("NC_000006.12", 1003, 1006).unwrap(),
            "CCC"
        );
    }

    #[test]
    fn over_range_read_clamps_to_contig_end() {
        // Window reaches the contig end (1012): a read past it returns the
        // available prefix, mirroring MultiFastaProvider (the short read that
        // drives CanonicalSplitSkipped / W5003).
        let p = fixture().to_provider();
        assert_eq!(
            p.get_genomic_sequence("NC_000006.12", 1006, 9999).unwrap(),
            "GGGTTT"
        );
    }

    #[test]
    fn within_contig_but_outside_window_is_extraction_defect_error() {
        let p = fixture().to_provider();
        // start before the window
        assert!(p.get_genomic_sequence("NC_000006.12", 999, 1005).is_err());
        // NC_DEFECT.1 contig is longer (5000) than its window (ends 1012): a
        // read that stays within the contig but past the window must error —
        // not clamp — so a too-narrow extraction is caught loudly.
        assert!(p.get_genomic_sequence("NC_DEFECT.1", 1006, 1500).is_err());
        // unknown contig is unavailable, not out-of-range
        assert!(matches!(
            p.get_genomic_sequence("NC_000009.12", 1, 2),
            Err(FerroError::GenomicReferenceNotAvailable { .. })
        ));
    }

    #[test]
    fn get_sequence_length_reports_true_contig_length() {
        let p = fixture().to_provider();
        assert_eq!(p.get_sequence_length("NC_DEFECT.1").unwrap(), 5000);
        assert_eq!(p.get_sequence_length("NM_TEST.1").unwrap(), 12);
        // Unversioned id resolves via the same base-accession fallback as
        // `get_transcript`, so the two lookups never disagree on the same id.
        assert_eq!(p.get_sequence_length("NM_TEST").unwrap(), 12);
    }

    #[test]
    fn get_sequence_dispatches_contig_vs_transcript() {
        let p = fixture().to_provider();
        // contig id → genomic path
        assert_eq!(p.get_sequence("NC_000006.12", 1000, 1003).unwrap(), "AAA");
        // transcript id → transcript-relative bases (0-based, half-open)
        assert_eq!(p.get_sequence("NM_TEST.1", 1, 4).unwrap(), "CGT");
    }

    #[test]
    fn roundtrips_through_json() {
        let json = fixture().to_json().unwrap();
        let back = serde_json::from_str::<WindowFixture>(&json).unwrap();
        let p = back.to_provider();
        assert_eq!(
            p.get_genomic_sequence("NC_000006.12", 1000, 1012).unwrap(),
            "AAACCCGGGTTT"
        );
        assert!(json.ends_with('\n'));
    }

    #[test]
    fn has_genomic_data_reflects_windows() {
        assert!(fixture().to_provider().has_genomic_data());
        let empty = WindowFixture {
            description: String::new(),
            captured_from: String::new(),
            contig_lengths: BTreeMap::new(),
            transcripts: vec![tx("NM_TEST.1", "ACGT")],
            genomic: vec![],
            placements: Vec::new(),
        };
        assert!(!empty.to_provider().has_genomic_data());
    }
}
