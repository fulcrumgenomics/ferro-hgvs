//! Authoritative transcript records and the canonical-overrides schema.
//!
//! The canonical-record validation effort (issue #520) needs the *authoritative*
//! RefSeq record for an exact accession version — its CDS coordinates, the
//! paired protein accession, and the transcript length — to detect (and later
//! correct) transcript records whose cdot-derived metadata disagrees with the
//! canonical record (e.g. `NM_012459.2` cdot-paired with `NP_036591.3` instead
//! of the canonical `NP_036591.2`).
//!
//! This module provides:
//!
//! - [`AuthoritativeRecord`] — the fields we compare against, parsed from a
//!   GenBank flat-file record ([`parse_authoritative_genbank`]). Parsing is
//!   pure and network-free; the actual NCBI efetch that produces the GenBank
//!   text lives in the `prepare` pipeline.
//! - [`CanonicalOverrides`] — a serialisable accession→record map, the
//!   "authoritative-overrides" file written by `ferro prepare` and consumed by
//!   the later validation / correction phases.
//!
//! All coordinates here are **1-based inclusive** (as written in the GenBank
//! `CDS` feature). cdot stores CDS as 0-based; callers comparing the two must
//! convert.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// The authoritative facts about a transcript at an exact accession version,
/// taken from its canonical RefSeq (GenBank) record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct AuthoritativeRecord {
    /// The exact versioned accession from the record's `VERSION` line
    /// (e.g. `NM_012459.2`).
    pub accession: String,
    /// Transcript length in nucleotides (the length of the `ORIGIN` sequence).
    pub tx_length: u64,
    /// CDS start, 1-based inclusive, in transcript coordinates. `None` for a
    /// non-coding record.
    pub cds_start: Option<u64>,
    /// CDS end, 1-based inclusive. `None` for a non-coding record.
    pub cds_end: Option<u64>,
    /// The protein accession paired with the CDS (`/protein_id`, e.g.
    /// `NP_036591.2`). `None` for a non-coding record or one without the
    /// qualifier.
    pub protein_id: Option<String>,
}

/// Parse the authoritative facts out of a single GenBank flat-file record.
///
/// Extracts the `VERSION` accession, the transcript length (counted from the
/// `ORIGIN` sequence), and — from the first `CDS` feature — the 1-based CDS
/// coordinates and the `/protein_id`. Returns `None` if the record has no
/// `VERSION` line or no `ORIGIN` sequence (i.e. it is not a usable record).
///
/// Location parsing handles the common RefSeq mRNA shapes: a plain `a..b`
/// range and a `join(a..b,...)` (including one that wraps onto continuation
/// lines), `complement(...)`, and the `<`/`>` partial-boundary markers.
/// Compound multi-segment CDS spans collapse to their outer bounds (first
/// segment start, last segment end), which is sufficient for length/pairing
/// validation.
///
/// Parses the **first** GenBank record only (it stops at the first `//`), so
/// callers should efetch one accession per call rather than a concatenated
/// multi-record batch.
pub fn parse_authoritative_genbank(genbank: &str) -> Option<AuthoritativeRecord> {
    let mut accession: Option<String> = None;
    let mut cds_start: Option<u64> = None;
    let mut cds_end: Option<u64> = None;
    let mut protein_id: Option<String> = None;
    let mut tx_length: u64 = 0;

    let mut in_origin = false;
    // True while reading the `CDS` feature's qualifier lines (so a stray
    // `/protein_id` under another feature isn't captured).
    let mut in_cds_feature = false;
    // The CDS location string, accumulated across any continuation lines (a
    // long `join(...)` wraps onto qualifier-indented lines) until the first
    // `/qualifier`; parsed once after the loop.
    let mut cds_loc = String::new();
    let mut collecting_loc = false;
    // Only the first CDS feature is used (its bounds and protein_id). A second
    // CDS (rare in mRNA records) must not overwrite the first's bounds while
    // keeping the first's protein_id.
    let mut seen_cds = false;

    for line in genbank.lines() {
        if line.starts_with("ORIGIN") {
            in_origin = true;
            in_cds_feature = false;
            collecting_loc = false;
            continue;
        }

        if in_origin {
            if line.starts_with("//") {
                break;
            }
            // Sequence lines: "        1 acgtacgt acgt...". Count alphabetic
            // bases only (the leading number and spaces are not sequence).
            tx_length += line.bytes().filter(|b| b.is_ascii_alphabetic()).count() as u64;
            continue;
        }

        if accession.is_none() {
            if let Some(rest) = line.strip_prefix("VERSION") {
                accession = rest.split_whitespace().next().map(str::to_string);
            }
        }

        // Feature keys sit at column 6 (exactly 5 leading spaces); qualifier
        // and location-continuation lines are indented further (~21 spaces).
        let trimmed = line.trim_start();
        let is_feature_key_line = line.starts_with("     ") && !line.starts_with("      ");
        if is_feature_key_line {
            // Accept only the first CDS feature; later features (CDS or not)
            // end the CDS block.
            let is_cds = !seen_cds && trimmed.split_whitespace().next() == Some("CDS");
            in_cds_feature = is_cds;
            collecting_loc = is_cds;
            if is_cds {
                seen_cds = true;
                cds_loc.clear();
                if let Some(loc) = trimmed.split_whitespace().nth(1) {
                    cds_loc.push_str(loc);
                }
            }
        } else if collecting_loc {
            // Inside the CDS feature, before the first qualifier: either a
            // wrapped continuation of the location, or the first `/qualifier`.
            if trimmed.starts_with('/') {
                collecting_loc = false;
                if protein_id.is_none() {
                    if let Some(pid) = extract_qualifier(trimmed, "protein_id") {
                        protein_id = Some(pid);
                    }
                }
            } else {
                cds_loc.push_str(trimmed); // wrapped location segment
            }
        } else if in_cds_feature && protein_id.is_none() {
            if let Some(pid) = extract_qualifier(trimmed, "protein_id") {
                protein_id = Some(pid);
            }
        }
    }

    if !cds_loc.is_empty() {
        if let Some((s, e)) = parse_cds_bounds(&cds_loc) {
            cds_start = Some(s);
            cds_end = Some(e);
        }
    }

    let accession = accession?;
    if tx_length == 0 {
        return None;
    }
    Some(AuthoritativeRecord {
        accession,
        tx_length,
        cds_start,
        cds_end,
        protein_id,
    })
}

/// Parse the outer 1-based bounds of a GenBank CDS location string.
///
/// Handles `a..b` and `join(a..b,...)`; tolerates `<`/`>` partial markers and
/// `complement(...)` wrappers. Returns `None` if no `..` range is present.
fn parse_cds_bounds(loc: &str) -> Option<(u64, u64)> {
    let inner = loc
        .strip_prefix("complement(")
        .map(|s| s.trim_end_matches(')'))
        .unwrap_or(loc);
    let inner = inner
        .strip_prefix("join(")
        .map(|s| s.trim_end_matches(')'))
        .unwrap_or(inner);

    // For join(...), the outer CDS bounds are the start of the first segment
    // and the end of the last segment. Each segment is a single `x..y` range.
    let first = inner.split(',').next()?;
    let last = inner.rsplit(',').next()?; // char pattern → reverse-iterable
    let start = first.split_once("..")?.0;
    let end = last.split_once("..")?.1;

    let start = start
        .trim_matches(|c: char| !c.is_ascii_digit())
        .parse()
        .ok()?;
    let end = end
        .trim_matches(|c: char| !c.is_ascii_digit())
        .parse()
        .ok()?;
    Some((start, end))
}

/// Extract a quoted GenBank qualifier value, e.g. `/protein_id="NP_036591.2"`
/// → `NP_036591.2`.
fn extract_qualifier(line: &str, name: &str) -> Option<String> {
    let needle = format!("/{name}=\"");
    let start = line.find(&needle)? + needle.len();
    let rest = &line[start..];
    let end = rest.find('"')?;
    Some(rest[..end].to_string())
}

/// Current on-disk schema version for [`CanonicalOverrides`].
pub const CANONICAL_OVERRIDES_SCHEMA_VERSION: u32 = 1;

/// An accession→[`AuthoritativeRecord`] map: the "authoritative-overrides" file
/// written by `ferro prepare --validate-canonical` and consumed by the
/// validation / correction phases. Keyed by exact versioned accession.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct CanonicalOverrides {
    /// On-disk schema version, so a future reader can detect a format change.
    /// Defaults to 0 when reading a file written before the field existed.
    #[serde(default)]
    pub schema_version: u32,
    /// Authoritative records keyed by exact versioned accession.
    pub records: BTreeMap<String, AuthoritativeRecord>,
}

impl Default for CanonicalOverrides {
    fn default() -> Self {
        Self {
            schema_version: CANONICAL_OVERRIDES_SCHEMA_VERSION,
            records: BTreeMap::new(),
        }
    }
}

impl CanonicalOverrides {
    /// Look up the authoritative record for an exact versioned accession.
    pub fn get(&self, accession: &str) -> Option<&AuthoritativeRecord> {
        self.records.get(accession)
    }

    /// Insert (or replace) a record, keyed by its own accession.
    pub fn insert(&mut self, record: AuthoritativeRecord) {
        self.records.insert(record.accession.clone(), record);
    }

    /// Number of records.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether the map is empty.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Serialise to pretty JSON.
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    /// Parse from JSON (the on-disk overrides file).
    pub fn from_json(json: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(json)
    }
}

/// Build a [`CanonicalOverrides`] for `accessions` by fetching and parsing each
/// one's GenBank record.
///
/// `fetch_gb` is the GenBank-text source: given an exact versioned accession it
/// returns the record text, or `None` if the fetch failed. Injecting it keeps
/// this orchestration network-free and testable — the `ferro prepare` caller
/// supplies an NCBI-efetch implementation, tests supply fixtures.
///
/// Returns the assembled overrides plus the list of accessions that could not
/// be fetched or parsed (so the caller can report them). A record whose parsed
/// `accession` differs from the requested one is still keyed by its own
/// `VERSION` accession (the authoritative truth), and the requested accession is
/// recorded as unresolved.
pub fn build_canonical_overrides<F>(
    accessions: &[String],
    mut fetch_gb: F,
) -> (CanonicalOverrides, Vec<String>)
where
    F: FnMut(&str) -> Option<String>,
{
    let mut overrides = CanonicalOverrides::default();
    let mut failed = Vec::new();
    for accession in accessions {
        match fetch_gb(accession)
            .as_deref()
            .and_then(parse_authoritative_genbank)
        {
            // The fetched record's VERSION must match what we asked for; a
            // mismatch means efetch resolved to a different version and the
            // record is not authoritative for `accession`.
            Some(record) if record.accession == *accession => overrides.insert(record),
            _ => failed.push(accession.clone()),
        }
    }
    (overrides, failed)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A minimal but realistic GenBank record: NM_012459.2-shaped, with the
    /// canonical CDS `31..327` paired with `NP_036591.2`, total length 30
    /// (sequence below) — the values matter, not biological realism.
    const FIXTURE: &str = "\
LOCUS       NM_012459                 30 bp    mRNA    linear   PRI 01-JAN-2020
DEFINITION  Homo sapiens translocase (TIMM8B), mRNA.
ACCESSION   NM_012459
VERSION     NM_012459.2
FEATURES             Location/Qualifiers
     source          1..30
                     /organism=\"Homo sapiens\"
     gene            1..30
                     /gene=\"TIMM8B\"
     CDS             4..30
                     /gene=\"TIMM8B\"
                     /codon_start=1
                     /product=\"mitochondrial import inner membrane translocase\"
                     /protein_id=\"NP_036591.2\"
                     /translation=\"MAEL\"
ORIGIN
        1 acgatggcag agctgtaagg ccacgagcct
//
";

    #[test]
    fn parses_version_length_cds_and_protein_id() {
        let rec = parse_authoritative_genbank(FIXTURE).expect("record should parse");
        assert_eq!(rec.accession, "NM_012459.2");
        assert_eq!(rec.tx_length, 30, "should count the 30 sequence bases");
        assert_eq!(rec.cds_start, Some(4));
        assert_eq!(rec.cds_end, Some(30));
        assert_eq!(rec.protein_id.as_deref(), Some("NP_036591.2"));
    }

    #[test]
    fn parses_join_location_outer_bounds() {
        let gb = "\
VERSION     NM_TEST.1
FEATURES             Location/Qualifiers
     CDS             join(10..20,31..45)
                     /protein_id=\"NP_TEST.1\"
ORIGIN
        1 acgtacgtac
//
";
        let rec = parse_authoritative_genbank(gb).expect("record should parse");
        assert_eq!(rec.cds_start, Some(10));
        assert_eq!(rec.cds_end, Some(45), "outer end of the last join segment");
        assert_eq!(rec.protein_id.as_deref(), Some("NP_TEST.1"));
    }

    #[test]
    fn parses_wrapped_join_location() {
        // A join() that wraps onto a continuation line (the CDS key line ends
        // with a trailing comma) must still yield outer bounds, not None.
        let gb = "\
VERSION     NM_WRAP.1
FEATURES             Location/Qualifiers
     CDS             join(1..100,200..300,
                     400..567)
                     /protein_id=\"NP_WRAP.1\"
ORIGIN
        1 acgtacgtac
//
";
        let rec = parse_authoritative_genbank(gb).expect("record should parse");
        assert_eq!(rec.cds_start, Some(1));
        assert_eq!(
            rec.cds_end,
            Some(567),
            "outer end of the wrapped last segment"
        );
        assert_eq!(rec.protein_id.as_deref(), Some("NP_WRAP.1"));
    }

    #[test]
    fn parses_complement_join_location() {
        let gb = "\
VERSION     NM_CJ.1
     CDS             complement(join(10..20,31..45))
                     /protein_id=\"NP_CJ.1\"
ORIGIN
        1 acgtacgtac
//
";
        let rec = parse_authoritative_genbank(gb).expect("record should parse");
        assert_eq!(rec.cds_start, Some(10));
        assert_eq!(rec.cds_end, Some(45));
    }

    #[test]
    fn first_cds_wins_over_a_second_cds() {
        // The first CDS's bounds AND protein_id are kept; a second CDS must not
        // overwrite the bounds.
        let gb = "\
VERSION     NM_TWO.1
     CDS             1..9
                     /protein_id=\"NP_FIRST.1\"
     CDS             10..30
                     /protein_id=\"NP_SECOND.2\"
ORIGIN
        1 acgtacgtac acgtacgtac acgtacgtac
//
";
        let rec = parse_authoritative_genbank(gb).expect("record should parse");
        assert_eq!(rec.cds_start, Some(1));
        assert_eq!(
            rec.cds_end,
            Some(9),
            "second CDS must not overwrite the first's bounds"
        );
        assert_eq!(rec.protein_id.as_deref(), Some("NP_FIRST.1"));
    }

    #[test]
    fn legacy_overrides_without_schema_version_default_to_zero() {
        // A file written before the schema_version field existed reads as 0.
        let ov = CanonicalOverrides::from_json(r#"{"records":{}}"#).unwrap();
        assert_eq!(ov.schema_version, 0);
        assert!(ov.is_empty());
    }

    #[test]
    fn tolerates_partial_boundary_markers() {
        let gb = "\
VERSION     NM_PARTIAL.1
     CDS             <1..>27
                     /protein_id=\"NP_PARTIAL.1\"
ORIGIN
        1 acgtacgtac acgtacgtac acgtacg
//
";
        let rec = parse_authoritative_genbank(gb).expect("record should parse");
        assert_eq!(rec.cds_start, Some(1));
        assert_eq!(rec.cds_end, Some(27));
    }

    #[test]
    fn noncoding_record_has_no_cds_or_protein() {
        let gb = "\
VERSION     NR_TEST.1
FEATURES             Location/Qualifiers
     gene            1..10
                     /gene=\"LINC\"
ORIGIN
        1 acgtacgtac
//
";
        let rec = parse_authoritative_genbank(gb).expect("record should parse");
        assert_eq!(rec.accession, "NR_TEST.1");
        assert_eq!(rec.tx_length, 10);
        assert_eq!(rec.cds_start, None);
        assert_eq!(rec.cds_end, None);
        assert_eq!(rec.protein_id, None);
    }

    #[test]
    fn returns_none_without_version_or_sequence() {
        assert!(parse_authoritative_genbank("LOCUS only, no version or origin\n").is_none());
        let no_seq = "VERSION     NM_X.1\nORIGIN\n//\n";
        assert!(
            parse_authoritative_genbank(no_seq).is_none(),
            "a record with no sequence bases is not usable"
        );
    }

    #[test]
    fn does_not_capture_protein_id_outside_the_cds_feature() {
        // A /protein_id appearing under a non-CDS feature must be ignored.
        let gb = "\
VERSION     NM_TEST.1
     misc_feature    1..10
                     /protein_id=\"NP_WRONG.9\"
     CDS             1..9
                     /protein_id=\"NP_RIGHT.1\"
ORIGIN
        1 acgtacgtac
//
";
        let rec = parse_authoritative_genbank(gb).expect("record should parse");
        assert_eq!(rec.protein_id.as_deref(), Some("NP_RIGHT.1"));
    }

    #[test]
    fn canonical_overrides_roundtrips_through_json() {
        let mut ov = CanonicalOverrides::default();
        assert_eq!(ov.schema_version, CANONICAL_OVERRIDES_SCHEMA_VERSION);
        ov.insert(AuthoritativeRecord {
            accession: "NM_012459.2".to_string(),
            tx_length: 327,
            cds_start: Some(31),
            cds_end: Some(327),
            protein_id: Some("NP_036591.2".to_string()),
        });
        assert_eq!(ov.len(), 1);
        let json = ov.to_json().unwrap();
        let back = CanonicalOverrides::from_json(&json).unwrap();
        assert_eq!(back, ov);
        assert_eq!(
            back.get("NM_012459.2")
                .and_then(|r| r.protein_id.as_deref()),
            Some("NP_036591.2")
        );
    }

    #[test]
    fn build_overrides_collects_parsed_records() {
        let accs = vec!["NM_012459.2".to_string()];
        let (ov, failed) = build_canonical_overrides(&accs, |acc| {
            assert_eq!(acc, "NM_012459.2");
            Some(FIXTURE.to_string())
        });
        assert!(failed.is_empty());
        assert_eq!(ov.len(), 1);
        assert_eq!(
            ov.get("NM_012459.2").and_then(|r| r.protein_id.as_deref()),
            Some("NP_036591.2")
        );
    }

    #[test]
    fn build_overrides_records_fetch_failures() {
        let accs = vec!["NM_NOPE.1".to_string()];
        let (ov, failed) = build_canonical_overrides(&accs, |_| None);
        assert!(ov.is_empty());
        assert_eq!(failed, vec!["NM_NOPE.1".to_string()]);
    }

    #[test]
    fn build_overrides_rejects_version_mismatch() {
        // Requested .9 but efetch returned the .2 record → not authoritative
        // for the requested version, so it is recorded as failed, not inserted.
        let accs = vec!["NM_012459.9".to_string()];
        let (ov, failed) = build_canonical_overrides(&accs, |_| Some(FIXTURE.to_string()));
        assert!(ov.is_empty());
        assert_eq!(failed, vec!["NM_012459.9".to_string()]);
    }
}
