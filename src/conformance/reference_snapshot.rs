//! Pinned, version-exact transcript reference snapshot for a conformance corpus.
//!
//! A conformance corpus pins exact `accession.version`s; to judge it
//! *reproducibly* the reference oracle must be pinned too. This module is the
//! committed, reviewable form of that oracle for the **transcript** axis: for
//! every transcript accession the corpus references (see
//! [`accession_inventory`](super::accession_inventory)), it carries the
//! transcript's bases at the exact pinned version plus that version's **own**
//! CDS bounds — never a fuzzy sibling version's frame (#714).
//!
//! The snapshot is two committed files in a corpus's fixture directory:
//!
//! - `transcripts.fna` — a plain FASTA, one record per versioned accession,
//!   carrying the exact bases. Standard format so the bytes are reviewable and
//!   greppable.
//! - `transcripts.metadata.json` — this [`TranscriptSnapshot`]: per-accession
//!   CDS bounds, gene/protein labels, and provenance (source + a
//!   `sha256` of the bases) so the committed sequence is auditable against its
//!   authoritative source.
//!
//! The builder (`examples/build_conformance_snapshot.rs`) harvests both files
//! from a prepared reference using the version-exact / cross-build-exact cdot
//! lookups (#717/#720); the harness (a later increment) loads them to serve the
//! corpus hermetically, with no out-of-band reference data.
//!
//! **CDS coordinate convention:** `cds_start` / `cds_end` are **1-based,
//! inclusive** transcript coordinates — the GenBank/HGVS-natural convention,
//! matching the existing `supplemental_cds` metadata (so a loader can feed this
//! snapshot straight into that path). `c.1` is the base at `cds_start`. This
//! differs from cdot's internal 0-based representation; the builder converts on
//! harvest (`cds_start = cdot.cds_start + 1`, `cds_end = cdot.cds_end`).

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::FerroError;

/// File name of the committed transcript FASTA within a snapshot directory.
pub const TRANSCRIPTS_FASTA: &str = "transcripts.fna";
/// File name of the committed transcript metadata within a snapshot directory.
pub const TRANSCRIPTS_METADATA: &str = "transcripts.metadata.json";

/// Where a snapshot entry's bytes came from, recorded so the committed sequence
/// is auditable against its authoritative source.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Provenance {
    /// How the entry was sourced. The builder emits
    /// `"manifest:transcript-fasta+cdot"` — the bases from the transcript FASTA
    /// manifest, paired with the cdot build that supplied the CDS.
    pub source: String,
    /// Lowercase hex SHA-256 of the entry's bases (uppercase, no newlines),
    /// for auditing the committed FASTA against the source.
    pub sha256: String,
}

/// One transcript's pinned metadata: its version-exact CDS bounds plus labels
/// and provenance. The bases themselves live in the companion FASTA, keyed by
/// the same versioned accession.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TranscriptEntry {
    /// CDS start, **1-based inclusive** transcript coordinate (`c.1` is here).
    /// `None` for a non-coding transcript (`NR_`) or when CDS is unknown.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cds_start: Option<u64>,
    /// CDS end, **1-based inclusive** transcript coordinate. `None` as above.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cds_end: Option<u64>,
    /// Gene symbol, if the source recorded one.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_symbol: Option<String>,
    /// Protein accession for the CDS (e.g. `"NP_002993.1"`), if recorded.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub protein: Option<String>,
    /// Transcript length in bases (must equal the companion FASTA record's
    /// length; the loader validates this).
    pub length: u64,
    /// Source + hash of the bases.
    pub provenance: Provenance,
}

/// The committed transcript snapshot metadata: every transcript accession the
/// corpus references, keyed by its versioned accession (e.g. `"NM_003002.2"`),
/// carrying CDS bounds + provenance. Mirrors
/// [`WindowFixture`](super::reference_window::WindowFixture)'s shape: a human
/// provenance note, the manifest it was captured from, and deterministic
/// ordering so `--check` byte-compares cleanly.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TranscriptSnapshot {
    /// Human-facing provenance note (schema version, generator command).
    #[serde(default)]
    pub description: String,
    /// Manifest `prepared_at` (or path) the snapshot was captured from, for
    /// traceability.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub captured_from: String,
    /// Per-accession metadata, ordered by versioned accession.
    pub transcripts: BTreeMap<String, TranscriptEntry>,
}

impl TranscriptSnapshot {
    /// Load the metadata from a JSON file.
    pub fn from_json_path<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        let content = std::fs::read_to_string(path.as_ref())?;
        let snapshot: TranscriptSnapshot = serde_json::from_str(&content)?;
        Ok(snapshot)
    }

    /// Serialize to pretty JSON with a trailing newline (matches the other
    /// generated fixtures so `--check` byte-compares cleanly).
    pub fn to_json(&self) -> Result<String, FerroError> {
        let mut s = serde_json::to_string_pretty(self)?;
        s.push('\n');
        Ok(s)
    }
}

/// Render `sequences` (keyed by versioned accession) as a FASTA string: records
/// in key order, one unwrapped sequence line each, trailing newline. The fixed
/// ordering and single-line bodies keep the committed file stable for `--check`.
pub fn render_fasta(sequences: &BTreeMap<String, String>) -> String {
    let mut out = String::new();
    for (accession, bases) in sequences {
        out.push('>');
        out.push_str(accession);
        out.push('\n');
        out.push_str(bases);
        out.push('\n');
    }
    out
}

/// Parse a FASTA string written by [`render_fasta`] into a map of versioned
/// accession → bases. Tolerant of wrapped sequence lines (joins them) so a
/// hand-inspected/edited file still loads; the header is taken up to the first
/// whitespace.
pub fn parse_fasta(text: &str) -> BTreeMap<String, String> {
    let mut map = BTreeMap::new();
    let mut current: Option<(String, String)> = None;
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some((acc, seq)) = current.take() {
                map.insert(acc, seq);
            }
            let acc = header.split_whitespace().next().unwrap_or("").to_string();
            current = Some((acc, String::new()));
        } else if let Some((_, seq)) = current.as_mut() {
            seq.push_str(line.trim());
        }
    }
    if let Some((acc, seq)) = current.take() {
        map.insert(acc, seq);
    }
    map
}

/// Read the FASTA companion (`transcripts.fna`) in a snapshot directory.
pub fn load_sequences<P: AsRef<Path>>(dir: P) -> Result<BTreeMap<String, String>, FerroError> {
    let path = dir.as_ref().join(TRANSCRIPTS_FASTA);
    let text = std::fs::read_to_string(&path)?;
    Ok(parse_fasta(&text))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(length: u64, cds: Option<(u64, u64)>) -> TranscriptEntry {
        TranscriptEntry {
            cds_start: cds.map(|c| c.0),
            cds_end: cds.map(|c| c.1),
            gene_symbol: Some("SDHD".to_string()),
            protein: Some("NP_002993.1".to_string()),
            length,
            provenance: Provenance {
                source: "manifest:transcript-fasta+cdot(GRCh37)".to_string(),
                sha256: "deadbeef".to_string(),
            },
        }
    }

    #[test]
    fn fasta_round_trips() {
        let mut seqs = BTreeMap::new();
        seqs.insert("NM_003002.2".to_string(), "ACGTACGT".to_string());
        seqs.insert("NR_024274.1".to_string(), "TTTT".to_string());
        let rendered = render_fasta(&seqs);
        // Records emitted in key order, one unwrapped line each.
        assert_eq!(rendered, ">NM_003002.2\nACGTACGT\n>NR_024274.1\nTTTT\n");
        assert_eq!(parse_fasta(&rendered), seqs);
    }

    #[test]
    fn parse_fasta_joins_wrapped_lines() {
        let text = ">NM_003002.2 SDHD\nACGT\nACGT\n";
        let parsed = parse_fasta(text);
        // Header trimmed at whitespace; wrapped body lines joined.
        assert_eq!(parsed.get("NM_003002.2"), Some(&"ACGTACGT".to_string()));
    }

    #[test]
    fn metadata_json_round_trips() {
        let mut snapshot = TranscriptSnapshot {
            description: "test snapshot".to_string(),
            captured_from: "2026-06-14".to_string(),
            transcripts: BTreeMap::new(),
        };
        snapshot
            .transcripts
            .insert("NM_003002.2".to_string(), entry(1382, Some((62, 541))));
        snapshot
            .transcripts
            .insert("NR_024274.1".to_string(), entry(2000, None));

        let json = snapshot.to_json().expect("serializes");
        assert!(json.ends_with('\n'), "trailing newline for stable --check");
        let parsed: TranscriptSnapshot = serde_json::from_str(&json).expect("round-trips");
        assert_eq!(parsed, snapshot);
    }

    #[test]
    fn noncoding_entry_omits_cds_fields() {
        // An `NR_` entry has no CDS; the fields are skipped, not serialized null.
        let mut snapshot = TranscriptSnapshot::default();
        snapshot
            .transcripts
            .insert("NR_024274.1".to_string(), entry(2000, None));
        let json = snapshot.to_json().expect("serializes");
        assert!(!json.contains("cds_start"), "cds_start omitted for NR_");
        assert!(!json.contains("cds_end"), "cds_end omitted for NR_");
    }

    #[test]
    fn unknown_field_is_rejected() {
        // `deny_unknown_fields` guards against silent schema drift.
        let json = r#"{"description":"x","transcripts":{},"bogus":1}"#;
        assert!(serde_json::from_str::<TranscriptSnapshot>(json).is_err());
    }
}
