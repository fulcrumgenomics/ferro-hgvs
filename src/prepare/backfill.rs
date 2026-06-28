//! Version-aware transcript backfill for `ferro prepare` (issue #842).
//!
//! `ferro prepare` populates the transcript FASTA from NCBI's **bulk** RefSeq
//! RNA release, a current-snapshot product that carries roughly one version per
//! accession. cdot, downloaded separately, carries exon alignments for many
//! **historical** versions. When a query names a superseded version whose bases
//! are absent from the prepared FASTA but whose alignment is in cdot, the
//! provider falls back to reconstructing bases from the genome via the cdot exon
//! CIGAR — a lossy path that cannot see deposited-transcript substitutions
//! (#807/#471/#400).
//!
//! This module closes that reference-coverage gap at prepare time: it diffs the
//! `accession.version` keys present in cdot against those present in the
//! transcript FASTA(s), then fetches the deposited sequence for the missing ones
//! and appends them to a backfill FASTA so the primary FASTA path serves them
//! and synthesis is never invoked for those accessions.
//!
//! The network fetch is injected (see [`execute_backfill`]), so the diff +
//! append logic is unit-testable with a mocked fetcher and CI stays offline. The
//! design mirrors [`crate::reference::authoritative::build_canonical_overrides`].

use crate::FerroError;
use std::collections::HashSet;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::path::Path;

/// The result of diffing a requested accession set against cdot and the
/// transcript FASTA(s).
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct BackfillPlan {
    /// Requested accessions that are in cdot but absent from the FASTA — these
    /// are the genuine gap to fetch. Sorted and deduplicated.
    pub to_fetch: Vec<String>,
    /// Requested accessions already present in the transcript FASTA(s) (which
    /// includes any prior backfill output) — skipped as an incremental-cache
    /// hit. Sorted and deduplicated.
    pub already_present: Vec<String>,
    /// Requested accessions absent from cdot — fetched bases would still serve
    /// the primary path, but with no exon alignment they fall outside the
    /// cdot-vs-FASTA gap this feature targets, so they are reported and skipped.
    /// Sorted and deduplicated.
    pub not_in_cdot: Vec<String>,
}

/// Diff `requested` against the `accession.version` keys present in cdot and in
/// the transcript FASTA(s), classifying each requested accession.
///
/// Classification is order-independent and deterministic: an accession already
/// present in the FASTA is always `already_present` (so a re-run never re-fetches
/// it); otherwise an accession in cdot is `to_fetch` and one absent from cdot is
/// `not_in_cdot`. Each output list is sorted and deduplicated.
pub fn plan_backfill(
    requested: &[String],
    cdot_versions: &HashSet<String>,
    present_versions: &HashSet<String>,
) -> BackfillPlan {
    let mut to_fetch = Vec::new();
    let mut already_present = Vec::new();
    let mut not_in_cdot = Vec::new();

    for acc in requested {
        if present_versions.contains(acc) {
            already_present.push(acc.clone());
        } else if cdot_versions.contains(acc) {
            to_fetch.push(acc.clone());
        } else {
            not_in_cdot.push(acc.clone());
        }
    }

    for v in [&mut to_fetch, &mut already_present, &mut not_in_cdot] {
        v.sort();
        v.dedup();
    }

    BackfillPlan {
        to_fetch,
        already_present,
        not_in_cdot,
    }
}

/// The outcome of fetching and appending the planned accessions.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct BackfillOutcome {
    /// Accessions whose deposited sequence was fetched and appended.
    pub fetched: Vec<String>,
    /// Accessions whose fetch or parse failed (warned, but not fatal).
    pub failed: Vec<String>,
}

/// Extract the concatenated, upper-cased sequence from a single FASTA record,
/// validating that the record's header accession equals `expected`.
///
/// The header token compared is the first whitespace-delimited field after the
/// leading `>` (e.g. `>NM_000088.3 Homo sapiens ...` → `NM_000088.3`). A header
/// whose accession differs from `expected` is rejected (returns `None`): NCBI
/// efetch resolves a bare accession to its current version, so this guard keeps
/// us from writing a *different* version's bases under the requested key (the
/// same requested-vs-fetched check `build_canonical_overrides` applies). An empty
/// sequence also returns `None`.
pub fn parse_fasta_record(text: &str, expected: &str) -> Option<String> {
    let mut lines = text.lines();
    let header = lines.next()?;
    let accession = header.strip_prefix('>')?.split_whitespace().next()?;
    if accession != expected {
        return None;
    }
    // Keep only the first record's bases: stop at the next record header so a
    // multi-record efetch response (or trailing text) can't be folded in.
    let sequence: String = lines
        .take_while(|line| !line.starts_with('>'))
        .flat_map(str::chars)
        .filter(|c| c.is_ascii_alphabetic())
        .map(|c| c.to_ascii_uppercase())
        .collect();
    if sequence.is_empty() {
        None
    } else {
        Some(sequence)
    }
}

/// Fetch each accession in `to_fetch` via the injected `fetch` and append its
/// sequence to `backfill_fasta`, returning which accessions were written vs
/// failed.
///
/// `fetch` is the FASTA-text source: given an exact versioned accession it
/// returns the efetch FASTA record text, or `None` on a fetch failure. Injecting
/// it keeps this orchestration network-free and testable.
///
/// The output is opened in **append** mode, and only **lazily on the first
/// successful record**, so an incremental re-run adds only the newly planned
/// records and a run where every fetch fails (or `to_fetch` is empty) leaves no
/// orphan 0-byte file behind — mirroring the no-orphan-empty-file rule in
/// `fetch_canonical_overrides`. The caller reindexes the whole file afterward.
/// Records are written `>{accession}\n` followed by the sequence wrapped at 70
/// columns, so the resulting `.fai` key is exactly the requested
/// `accession.version`.
///
/// A per-accession fetch, parse, or version-mismatch failure is collected into
/// `BackfillOutcome::failed` and does not abort; only an I/O error on the output
/// file returns `Err`.
pub fn execute_backfill<F>(
    to_fetch: &[String],
    backfill_fasta: &Path,
    mut fetch: F,
) -> Result<BackfillOutcome, FerroError>
where
    F: FnMut(&str) -> Option<String>,
{
    let mut outcome = BackfillOutcome::default();
    // Created lazily on the first record actually written, so a fully-failed or
    // empty run never orphans a 0-byte FASTA.
    let mut writer: Option<BufWriter<std::fs::File>> = None;

    for acc in to_fetch {
        let Some(sequence) = fetch(acc)
            .as_deref()
            .and_then(|t| parse_fasta_record(t, acc))
        else {
            outcome.failed.push(acc.clone());
            continue;
        };

        if writer.is_none() {
            let file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(backfill_fasta)
                .map_err(|e| FerroError::Io {
                    msg: format!(
                        "Failed to open backfill FASTA {}: {}",
                        backfill_fasta.display(),
                        e
                    ),
                })?;
            writer = Some(BufWriter::new(file));
        }
        let w = writer.as_mut().expect("writer initialized above");

        writeln!(w, ">{}", acc).map_err(|e| FerroError::Io {
            msg: format!("Failed to write backfill record {}: {}", acc, e),
        })?;
        for chunk in sequence.as_bytes().chunks(70) {
            writeln!(w, "{}", std::str::from_utf8(chunk).unwrap_or("")).map_err(|e| {
                FerroError::Io {
                    msg: format!("Failed to write backfill sequence {}: {}", acc, e),
                }
            })?;
        }
        outcome.fetched.push(acc.clone());
    }

    if let Some(mut w) = writer {
        w.flush().map_err(|e| FerroError::Io {
            msg: format!("Failed to flush backfill FASTA: {}", e),
        })?;
    }

    Ok(outcome)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;
    use std::fs;

    fn set(items: &[&str]) -> HashSet<String> {
        items.iter().map(|s| s.to_string()).collect()
    }

    fn vec_of(items: &[&str]) -> Vec<String> {
        items.iter().map(|s| s.to_string()).collect()
    }

    // ---- plan_backfill ----

    #[test]
    fn plan_identifies_exactly_the_missing() {
        let plan = plan_backfill(
            &vec_of(&["NM_A.1", "NM_B.2", "NM_C.3"]),
            &set(&["NM_A.1", "NM_B.2", "NM_C.3", "NM_D.4"]),
            &set(&["NM_A.1"]),
        );
        assert_eq!(plan.to_fetch, vec_of(&["NM_B.2", "NM_C.3"]));
        assert_eq!(plan.already_present, vec_of(&["NM_A.1"]));
        assert!(plan.not_in_cdot.is_empty());
    }

    #[test]
    fn plan_routes_not_in_cdot_never_to_fetch() {
        let plan = plan_backfill(&vec_of(&["NM_GHOST.9"]), &set(&["NM_A.1"]), &set(&[]));
        assert!(plan.to_fetch.is_empty());
        assert_eq!(plan.not_in_cdot, vec_of(&["NM_GHOST.9"]));
    }

    #[test]
    fn plan_present_wins_for_incremental_rerun() {
        // An accession already written by a prior backfill (so in present_versions)
        // must be skipped even though it is also in cdot.
        let plan = plan_backfill(&vec_of(&["NM_B.2"]), &set(&["NM_B.2"]), &set(&["NM_B.2"]));
        assert!(plan.to_fetch.is_empty());
        assert_eq!(plan.already_present, vec_of(&["NM_B.2"]));
    }

    #[test]
    fn plan_sorts_and_deduplicates() {
        let plan = plan_backfill(
            &vec_of(&["NM_C.3", "NM_B.2", "NM_B.2", "NM_A.1"]),
            &set(&["NM_A.1", "NM_B.2", "NM_C.3"]),
            &set(&[]),
        );
        assert_eq!(plan.to_fetch, vec_of(&["NM_A.1", "NM_B.2", "NM_C.3"]));
    }

    #[test]
    fn plan_empty_requested_is_empty() {
        let plan = plan_backfill(&[], &set(&["NM_A.1"]), &set(&[]));
        assert_eq!(plan, BackfillPlan::default());
    }

    // ---- parse_fasta_record ----

    #[test]
    fn parse_extracts_sequence_uppercased() {
        let text = ">NM_000088.3 Homo sapiens\nacgt\nACGT\n";
        assert_eq!(
            parse_fasta_record(text, "NM_000088.3"),
            Some("ACGTACGT".to_string())
        );
    }

    #[test]
    fn parse_rejects_header_accession_mismatch() {
        // efetch resolved the bare accession to a newer version: reject it.
        let text = ">NM_000088.4 newer version\nACGT\n";
        assert_eq!(parse_fasta_record(text, "NM_000088.3"), None);
    }

    #[test]
    fn parse_keeps_only_first_record_bases() {
        // A multi-record response must not fold the second record's bases in.
        let text = ">NM_A.1 first\nACGT\n>NM_B.2 second\nTTTT\n";
        assert_eq!(parse_fasta_record(text, "NM_A.1"), Some("ACGT".to_string()));
    }

    #[test]
    fn parse_rejects_empty_sequence() {
        let text = ">NM_000088.3 header only\n";
        assert_eq!(parse_fasta_record(text, "NM_000088.3"), None);
    }

    // ---- execute_backfill ----

    #[test]
    fn execute_appends_fetched_records_with_accession_headers() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("backfill.fna");
        let outcome = execute_backfill(&vec_of(&["NM_A.1", "NM_B.2"]), &fasta, |acc| {
            Some(format!(">{} desc\nACGTACGT\n", acc))
        })
        .unwrap();
        assert_eq!(outcome.fetched, vec_of(&["NM_A.1", "NM_B.2"]));
        assert!(outcome.failed.is_empty());

        let body = fs::read_to_string(&fasta).unwrap();
        assert!(body.contains(">NM_A.1\n"));
        assert!(body.contains(">NM_B.2\n"));
        assert_eq!(body.matches('>').count(), 2);
    }

    #[test]
    fn execute_continues_on_fetch_failure() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("backfill.fna");
        let outcome = execute_backfill(&vec_of(&["NM_OK.1", "NM_BAD.2"]), &fasta, |acc| {
            (acc == "NM_OK.1").then(|| format!(">{} desc\nACGT\n", acc))
        })
        .unwrap();
        assert_eq!(outcome.fetched, vec_of(&["NM_OK.1"]));
        assert_eq!(outcome.failed, vec_of(&["NM_BAD.2"]));
    }

    #[test]
    fn execute_rejects_version_mismatch_as_failed() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("backfill.fna");
        let outcome = execute_backfill(&vec_of(&["NM_X.1"]), &fasta, |_acc| {
            // efetch returns the *current* version, not the requested .1.
            Some(">NM_X.2 newer\nACGT\n".to_string())
        })
        .unwrap();
        assert!(outcome.fetched.is_empty());
        assert_eq!(outcome.failed, vec_of(&["NM_X.1"]));
        // Nothing written for the mismatched record.
        let body = fs::read_to_string(&fasta).unwrap_or_default();
        assert!(!body.contains('>'));
    }

    #[test]
    fn execute_appends_incrementally_without_clobbering() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("backfill.fna");
        execute_backfill(&vec_of(&["NM_A.1"]), &fasta, |acc| {
            Some(format!(">{} desc\nACGT\n", acc))
        })
        .unwrap();
        execute_backfill(&vec_of(&["NM_B.2"]), &fasta, |acc| {
            Some(format!(">{} desc\nGGGG\n", acc))
        })
        .unwrap();
        let body = fs::read_to_string(&fasta).unwrap();
        assert!(body.contains(">NM_A.1\n"), "first record preserved");
        assert!(body.contains(">NM_B.2\n"), "second record appended");
    }

    #[test]
    fn execute_leaves_no_orphan_file_when_all_fetches_fail() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("backfill.fna");
        let outcome = execute_backfill(&vec_of(&["NM_A.1", "NM_B.2"]), &fasta, |_| None).unwrap();
        assert!(outcome.fetched.is_empty());
        assert_eq!(outcome.failed, vec_of(&["NM_A.1", "NM_B.2"]));
        assert!(
            !fasta.exists(),
            "a fully-failed run must not orphan a 0-byte FASTA"
        );
    }

    #[test]
    fn execute_empty_plan_writes_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("backfill.fna");
        let outcome = execute_backfill(&[], &fasta, |_| panic!("must not fetch")).unwrap();
        assert_eq!(outcome, BackfillOutcome::default());
        assert!(!fasta.exists(), "no file created for an empty plan");
    }
}
