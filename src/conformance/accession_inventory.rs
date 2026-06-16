//! Corpus accession inventory — the version-exact set of reference accessions a
//! conformance corpus refers to.
//!
//! A corpus row pins specific `accession.version`s in its input (and expected
//! outputs). To judge that corpus reproducibly, the reference oracle must be
//! *version-complete* for it: every accession any row references, carried at the
//! exact pinned version. This module answers "what does the corpus need?" — it
//! scans the corpus text and returns the referenced accessions, classified by
//! the reference source that must supply each (transcript FASTA + CDS, genomic
//! FASTA, protein FASTA, Ensembl, LRG).
//!
//! It is the foundation of the pinned-snapshot builder (which fetches each
//! inventoried accession by exact version) and the systematic form of the
//! reference-availability triage behind the `reference_unavailable` disposition:
//! [`Inventory::missing_against`] diffs the corpus's needs against the set of
//! accessions a prepared reference actually carries.
//!
//! Extraction is **regex-based, not parser-based**, on purpose: corpus error
//! cases are deliberately unparseable, and cross-reference accessions are
//! embedded inside `delins` payloads (e.g. `…delinsNM_003002.4:100_102`) where a
//! variant parse would not surface them. A lexical scan over the well-defined
//! accession grammar catches both.

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::sync::LazyLock;

use regex::Regex;

use super::mutalyzer::Fixture;

/// The reference source class an accession belongs to — i.e. which kind of
/// reference data a version-complete snapshot must carry for it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum AccessionClass {
    /// RefSeq transcript (`NM_`/`NR_`/`XM_`/`XR_`) — needs the transcript
    /// sequence and its CDS bounds (the data missing for `NM_003002.2`).
    RefSeqTranscript,
    /// RefSeq genomic (`NG_`/`NC_`/`NW_`/`NT_`/`NZ_`) — needs the genomic
    /// sequence and, for projection, its placement.
    RefSeqGenomic,
    /// RefSeq protein (`NP_`/`XP_`) — needs the protein sequence (protein axes).
    RefSeqProtein,
    /// Ensembl (`ENST`/`ENSG`/`ENSP`) — absent from a RefSeq-only reference.
    Ensembl,
    /// LRG (`LRG_<n>`) — needs the LRG record (XML/FASTA).
    Lrg,
}

impl AccessionClass {
    /// Stable snake_case identifier for grouping/printing.
    pub fn as_str(self) -> &'static str {
        match self {
            AccessionClass::RefSeqTranscript => "refseq_transcript",
            AccessionClass::RefSeqGenomic => "refseq_genomic",
            AccessionClass::RefSeqProtein => "refseq_protein",
            AccessionClass::Ensembl => "ensembl",
            AccessionClass::Lrg => "lrg",
        }
    }

    /// Classify by accession prefix. Returns `None` for an unrecognized prefix.
    fn from_accession(accession: &str) -> Option<AccessionClass> {
        // RefSeq: the two-letter prefix before `_` determines the molecule type.
        if let Some(prefix) = accession.split('_').next() {
            match prefix {
                "NM" | "NR" | "XM" | "XR" => return Some(AccessionClass::RefSeqTranscript),
                "NG" | "NC" | "NW" | "NT" | "NZ" => return Some(AccessionClass::RefSeqGenomic),
                "NP" | "XP" => return Some(AccessionClass::RefSeqProtein),
                "LRG" => return Some(AccessionClass::Lrg),
                _ => {}
            }
        }
        if accession.starts_with("ENST")
            || accession.starts_with("ENSG")
            || accession.starts_with("ENSP")
        {
            return Some(AccessionClass::Ensembl);
        }
        None
    }
}

/// One accession reference found in corpus text: the bare accession, the pinned
/// version (if any), and its class.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AccessionRef {
    /// Bare accession without version, e.g. `"NM_003002"`, `"LRG_199"`.
    pub accession: String,
    /// The pinned version, or `None` when the corpus referenced it unversioned
    /// (e.g. `ENST00000375549:c.100del`).
    pub version: Option<u32>,
    /// Reference source class.
    pub class: AccessionClass,
}

impl AccessionRef {
    /// The versioned string as it would key a reference index — `"NM_003002.2"`,
    /// or the bare accession when unversioned.
    pub fn versioned(&self) -> String {
        match self.version {
            Some(v) => format!("{}.{}", self.accession, v),
            None => self.accession.clone(),
        }
    }
}

/// Matches an accession token (RefSeq / Ensembl / LRG) with an optional
/// `.version`. Distinctive prefixes mean coordinate digits never false-match.
/// LRG is intentionally versionless: `LRG_24t1` matches `LRG_24` and stops
/// before the `t1` transcript selector.
static ACCESSION_RE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(
        r"(?P<acc>(?:NM_|NR_|XM_|XR_|NG_|NC_|NW_|NT_|NZ_|NP_|XP_)\d+|ENS[TGP]\d+|LRG_\d+)(?:\.(?P<ver>\d+))?",
    )
    .expect("accession regex is valid")
});

/// Extract every accession reference in `text`, in order of appearance
/// (duplicates preserved). Robust to unparseable inputs and embedded
/// cross-reference accessions.
pub fn extract_accessions(text: &str) -> Vec<AccessionRef> {
    ACCESSION_RE
        .captures_iter(text)
        .filter_map(|caps| {
            let acc = caps.name("acc")?.as_str().to_string();
            let class = AccessionClass::from_accession(&acc)?;
            // LRG carries no `.version`; ignore any captured version for it.
            let version = if class == AccessionClass::Lrg {
                None
            } else {
                caps.name("ver")
                    .and_then(|m| m.as_str().parse::<u32>().ok())
            };
            Some(AccessionRef {
                accession: acc,
                version,
                class,
            })
        })
        .collect()
}

/// One inventory entry: an accession reference plus the set of corpus inputs
/// that reference it (for traceability when provisioning or triaging).
#[derive(Debug, Clone)]
pub struct InventoryEntry {
    /// The referenced accession (versioned as the corpus pinned it).
    pub acc_ref: AccessionRef,
    /// The `case.input` strings that reference this accession (sorted, unique).
    pub referencing_inputs: BTreeSet<String>,
}

/// The version-complete accession inventory of a corpus, keyed by the
/// versioned-accession string so a pinned `.2` and an unversioned reference to
/// the same accession are tracked as distinct needs.
#[derive(Debug, Clone, Default)]
pub struct Inventory {
    entries: BTreeMap<String, InventoryEntry>,
}

impl Inventory {
    /// Build the inventory for `fixture`: scan every case's input and expected
    /// outputs for referenced accessions. The input side is what ferro must
    /// resolve to avoid a no-op; the expected-output side adds accessions the
    /// corpus judges ferro's rendering/translation against (e.g. `NP_` proteins).
    pub fn from_fixture(fixture: &Fixture) -> Inventory {
        let mut inv = Inventory::default();
        for case in &fixture.cases {
            for text in case_reference_texts(case) {
                for acc_ref in extract_accessions(&text) {
                    inv.entries
                        .entry(acc_ref.versioned())
                        .or_insert_with(|| InventoryEntry {
                            acc_ref: acc_ref.clone(),
                            referencing_inputs: BTreeSet::new(),
                        })
                        .referencing_inputs
                        .insert(case.input.clone());
                }
            }
        }
        inv
    }

    /// All entries, ordered by versioned accession.
    pub fn entries(&self) -> impl Iterator<Item = &InventoryEntry> {
        self.entries.values()
    }

    /// Number of distinct (accession, version) references.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the inventory is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Count of distinct references in each class.
    pub fn counts_by_class(&self) -> BTreeMap<AccessionClass, usize> {
        let mut counts = BTreeMap::new();
        for entry in self.entries.values() {
            *counts.entry(entry.acc_ref.class).or_insert(0) += 1;
        }
        counts
    }

    /// The corpus references the snapshot does NOT yet cover: entries whose
    /// versioned accession is absent from `available` (e.g. the union of a
    /// prepared reference's FASTA-index names). This is the systematic form of
    /// the `reference_unavailable` triage — every returned entry is a row that
    /// no-ops for lack of a reference, not a ferro bug. Ordered by versioned
    /// accession.
    pub fn missing_against<'a>(&'a self, available: &HashSet<String>) -> Vec<&'a InventoryEntry> {
        self.entries
            .iter()
            .filter(|(versioned, _)| !available.contains(*versioned))
            .map(|(_, entry)| entry)
            .collect()
    }
}

/// The corpus text fields that may carry reference accessions, for one case:
/// the input plus the per-axis expected outputs (including `noncoding`, which
/// holds full `NG_(NR_):n.` descriptions, not bare codes). `errors`/`infos`
/// carry mutalyzer codes, not accessions, so they are excluded.
fn case_reference_texts(case: &super::mutalyzer::Case) -> Vec<String> {
    let mut texts = vec![case.input.clone()];
    texts.extend(case.normalized.clone());
    texts.extend(case.genomic.clone());
    texts.extend(case.protein_description.clone());
    texts.extend(case.rna_description.clone());
    if let Some(noncoding) = &case.noncoding {
        texts.extend(noncoding.iter().cloned());
    }
    if let Some(pairs) = &case.coding_protein_descriptions {
        for pair in pairs {
            texts.extend(pair.iter().cloned());
        }
    }
    texts
}

#[cfg(test)]
mod tests {
    use super::*;

    fn refs(text: &str) -> Vec<(String, Option<u32>, AccessionClass)> {
        extract_accessions(text)
            .into_iter()
            .map(|r| (r.accession, r.version, r.class))
            .collect()
    }

    #[test]
    fn extracts_versioned_refseq_transcript() {
        assert_eq!(
            refs("NM_003002.2:c.273del"),
            vec![(
                "NM_003002".to_string(),
                Some(2),
                AccessionClass::RefSeqTranscript
            )]
        );
    }

    #[test]
    fn extracts_both_accessions_from_genomic_context() {
        // NG_(NM_) compound: the genomic parent and the transcript are distinct
        // needs.
        assert_eq!(
            refs("NG_012337.1(NM_003002.2):c.273del"),
            vec![
                (
                    "NG_012337".to_string(),
                    Some(1),
                    AccessionClass::RefSeqGenomic
                ),
                (
                    "NM_003002".to_string(),
                    Some(2),
                    AccessionClass::RefSeqTranscript
                ),
            ]
        );
    }

    #[test]
    fn extracts_embedded_delins_cross_reference() {
        // The cross-reference accession inside a delins payload is a real
        // reference need a variant parse would not surface.
        let got = refs("NG_012337.1(NM_012459.2):c.274delinsNG_012337.3(NM_003002.4):c.200_203");
        assert!(got.contains(&(
            "NM_003002".to_string(),
            Some(4),
            AccessionClass::RefSeqTranscript
        )));
        assert!(got.contains(&(
            "NG_012337".to_string(),
            Some(3),
            AccessionClass::RefSeqGenomic
        )));
    }

    #[test]
    fn unversioned_ensembl_has_no_version() {
        assert_eq!(
            refs("ENST00000375549:c.100del"),
            vec![("ENST00000375549".to_string(), None, AccessionClass::Ensembl)]
        );
    }

    #[test]
    fn versioned_ensembl_gene_and_transcript() {
        assert_eq!(
            refs("ENSG00000204370.13(ENST00000375549.8):c.100del"),
            vec![
                (
                    "ENSG00000204370".to_string(),
                    Some(13),
                    AccessionClass::Ensembl
                ),
                (
                    "ENST00000375549".to_string(),
                    Some(8),
                    AccessionClass::Ensembl
                ),
            ]
        );
    }

    #[test]
    fn lrg_is_versionless_and_stops_before_transcript_selector() {
        // `LRG_199t1` → accession `LRG_199` (the `t1` selector is not a version).
        assert_eq!(
            refs("LRG_199t1:c.11del"),
            vec![("LRG_199".to_string(), None, AccessionClass::Lrg)]
        );
    }

    #[test]
    fn protein_accession_classified() {
        assert_eq!(
            refs("NP_036591.2:p.(Arg2Val)"),
            vec![(
                "NP_036591".to_string(),
                Some(2),
                AccessionClass::RefSeqProtein
            )]
        );
    }

    #[test]
    fn unparseable_text_still_scans_known_tokens() {
        // Error-case-style junk: no valid HGVS, but an accession token is still
        // surfaced (and non-accession noise is ignored).
        assert_eq!(
            refs("garbage NM_000193.2 !!! 273del ???"),
            vec![(
                "NM_000193".to_string(),
                Some(2),
                AccessionClass::RefSeqTranscript
            )]
        );
    }

    #[test]
    fn no_false_match_from_coordinates() {
        // Bare coordinates / edits carry no accession prefix → nothing matched.
        assert!(refs("c.273del").is_empty());
        assert!(refs("g.5525_5532delinsGTGGGAATTGT").is_empty());
    }

    fn fixture_from_cases(cases_json: &str) -> Fixture {
        let json = format!(
            r#"{{"description":"t","source":"t","source_commit":"t","license":"t","refreshed_at":"t","cases":[{cases_json}]}}"#
        );
        serde_json::from_str(&json).expect("fixture deserializes")
    }

    #[test]
    fn inventory_collects_distinct_versions_and_tracks_inputs() {
        let fixture = fixture_from_cases(
            r#"
            {"input":"NM_003002.2:c.273del","normalized":"NM_003002.2:c.274del"},
            {"input":"NM_003002.4:c.273del","normalized":"NM_003002.4:c.274del"},
            {"input":"ENST00000375549:c.100del"}
            "#,
        );
        let inv = Inventory::from_fixture(&fixture);
        // Three distinct (accession, version) needs: NM_003002.2, .4, ENST…(none).
        assert_eq!(inv.len(), 3);

        let keys: Vec<String> = inv.entries().map(|e| e.acc_ref.versioned()).collect();
        assert!(keys.contains(&"NM_003002.2".to_string()));
        assert!(keys.contains(&"NM_003002.4".to_string()));
        assert!(keys.contains(&"ENST00000375549".to_string()));

        // .2 is referenced only by the .2 input.
        let v2 = inv
            .entries()
            .find(|e| e.acc_ref.versioned() == "NM_003002.2")
            .unwrap();
        assert_eq!(
            v2.referencing_inputs.iter().cloned().collect::<Vec<_>>(),
            vec!["NM_003002.2:c.273del".to_string()]
        );
    }

    #[test]
    fn inventory_scans_noncoding_axis() {
        // Corpus `noncoding` rows are accession-bearing `NG_(NR_):n.`
        // descriptions, not bare codes, so both accessions must be captured and
        // linked to the case input.
        let fixture = fixture_from_cases(
            r#"
            {"input":"NG_007485.1(NM_000077.4):c.151_153del","noncoding":["NG_007485.1(NR_024274.1):n.616+3443_616+3444insGAT"]}
            "#,
        );
        let inv = Inventory::from_fixture(&fixture);

        let ng = inv
            .entries()
            .find(|e| e.acc_ref.versioned() == "NG_007485.1")
            .expect("NG_ from noncoding axis is captured");
        assert!(ng
            .referencing_inputs
            .contains("NG_007485.1(NM_000077.4):c.151_153del"));

        let nr = inv
            .entries()
            .find(|e| e.acc_ref.versioned() == "NR_024274.1")
            .expect("NR_ from noncoding axis is captured");
        assert_eq!(nr.acc_ref.class, AccessionClass::RefSeqTranscript);
        assert!(nr
            .referencing_inputs
            .contains("NG_007485.1(NM_000077.4):c.151_153del"));
    }

    #[test]
    fn missing_against_reports_only_uncovered() {
        let fixture = fixture_from_cases(
            r#"
            {"input":"NM_003002.2:c.273del"},
            {"input":"NM_003002.4:c.273del"}
            "#,
        );
        let inv = Inventory::from_fixture(&fixture);
        let mut available = HashSet::new();
        available.insert("NM_003002.4".to_string()); // .4 present, .2 absent

        let missing = inv.missing_against(&available);
        assert_eq!(missing.len(), 1);
        assert_eq!(missing[0].acc_ref.versioned(), "NM_003002.2");
    }

    #[test]
    fn counts_by_class_groups_references() {
        let fixture = fixture_from_cases(
            r#"
            {"input":"NG_012337.1(NM_003002.2):c.273del","protein_description":"NP_002993.1:p.(Asp92ThrfsTer43)"}
            "#,
        );
        let inv = Inventory::from_fixture(&fixture);
        let counts = inv.counts_by_class();
        assert_eq!(counts.get(&AccessionClass::RefSeqGenomic), Some(&1));
        assert_eq!(counts.get(&AccessionClass::RefSeqTranscript), Some(&1));
        assert_eq!(counts.get(&AccessionClass::RefSeqProtein), Some(&1));
    }
}
