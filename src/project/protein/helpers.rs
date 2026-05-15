//! Pure helper functions for protein consequence prediction.
//!
//! All functions in this module are side-effect-free and easily unit-tested
//! independently of the rest of the projector machinery.

use crate::error::FerroError;
use crate::hgvs::edit::{InsertedSequence, NaEdit};
use crate::hgvs::location::AminoAcid;
use crate::reference::transcript::Transcript;
use crate::sequence::reverse_complement;

use super::substitution::translate;

// ── CDS sequence readers ──────────────────────────────────────────────────────

/// Read the full CDS as an uppercase `String` from a transcript.
///
/// `cds_start` and `cds_end` are the 1-based inclusive positions of the CDS in
/// transcript space (i.e. `Transcript.cds_start` / `Transcript.cds_end`).
pub(crate) fn read_full_cds(transcript: &Transcript) -> Result<String, FerroError> {
    let cds_start = transcript
        .cds_start
        .ok_or_else(|| FerroError::ConversionError {
            msg: format!("transcript {} has no CDS start", transcript.id),
        })? as usize;
    let cds_end = transcript
        .cds_end
        .ok_or_else(|| FerroError::ConversionError {
            msg: format!("transcript {} has no CDS end", transcript.id),
        })? as usize;
    let seq =
        transcript
            .sequence
            .as_deref()
            .ok_or_else(|| FerroError::ProteinSequenceUnavailable {
                accession: transcript.id.clone(),
            })?;
    if cds_start < 1 || cds_end > seq.len() || cds_start > cds_end + 1 {
        return Err(FerroError::ProteinSequenceUnavailable {
            accession: transcript.id.clone(),
        });
    }
    // Transcript.cds_start is 1-based → 0-based index = cds_start - 1.
    // Transcript.cds_end is 1-based inclusive → exclusive index = cds_end.
    Ok(seq[cds_start - 1..cds_end].to_uppercase())
}

/// Build the mutated CDS sequence by applying `edit` at a given CDS coordinate.
///
/// `cds_start` and `cds_end` are **1-based CDS positions** (as recorded in a c. variant),
/// not transcript-space positions.  For a point edit (substitution, point deletion) the
/// caller sets `cds_start == cds_end`.
///
/// The returned string is the full CDS after the edit, in uppercase.
pub(crate) fn build_mutated_cds(
    transcript: &Transcript,
    cds_pos_start: i64, // 1-based CDS position of edit start
    cds_pos_end: i64,   // 1-based CDS position of edit end (inclusive)
    edit: &NaEdit,
) -> Result<String, FerroError> {
    let tx_cds_start = transcript
        .cds_start
        .ok_or_else(|| FerroError::ConversionError {
            msg: format!("transcript {} has no CDS", transcript.id),
        })? as usize;

    let full_cds = read_full_cds(transcript)?;

    // Convert 1-based CDS positions to 0-based CDS-string indices.
    let idx_start = (cds_pos_start - 1) as usize;
    let idx_end = cds_pos_end as usize; // exclusive upper bound

    if idx_end > full_cds.len() {
        return Err(FerroError::ProteinSequenceUnavailable {
            accession: transcript.id.clone(),
        });
    }

    let before = &full_cds[..idx_start];
    let affected = &full_cds[idx_start..idx_end];
    let after = &full_cds[idx_end..];

    let mutated = match edit {
        NaEdit::Deletion { .. } => {
            // Remove the affected bases.
            format!("{}{}", before, after)
        }
        NaEdit::Insertion { sequence } => {
            // HGVS insertion c.N_N+1insX: the inserted bases go between cds_pos_start
            // and cds_pos_start+1, i.e. at idx_start (after `before`).
            // The `affected` region here should be 0 bases wide (cds_pos_start == cds_pos_end
            // for most insertions, or the caller can pass both as the same boundary base).
            // We treat `before` as everything up through and including cds_pos_start.
            let inserted = extract_literal_sequence(sequence, transcript)?;
            // For HGVS insertions, idx_end == idx_start (the interval covers the gap),
            // so before == full_cds[..idx_start] and after == full_cds[idx_start..].
            // We want to insert between position idx_start-1 and idx_start.
            // The caller passes cds_pos_start as the last unaffected base's c. position,
            // so idx_end == idx_start and we insert at that junction.
            format!("{}{}{}{}", before, affected, inserted, after)
        }
        NaEdit::Duplication { .. } => {
            // Duplicate: append a copy of `affected` right after its end.
            format!("{}{}{}{}", before, affected, affected, after)
        }
        NaEdit::Delins { sequence, .. } => {
            let inserted = extract_literal_sequence(sequence, transcript)?;
            format!("{}{}{}", before, inserted, after)
        }
        NaEdit::Inversion { .. } => {
            let rc = reverse_complement(affected);
            format!("{}{}{}", before, rc, after)
        }
        NaEdit::Substitution { alternative, .. } => {
            // Single-base substitution.
            let alt_char = alternative.to_u8() as char;
            format!("{}{}{}", before, alt_char, after)
        }
        _ => {
            return Err(FerroError::UnsupportedProjection {
                reason: format!("build_mutated_cds does not support edit type: {:?}", edit),
            })
        }
    };

    let _ = tx_cds_start; // currently unused but kept for clarity
    Ok(mutated.to_uppercase())
}

/// Extract a literal nucleotide sequence from an `InsertedSequence`.
///
/// Returns an error if the sequence is not a literal (we don't have sequence
/// data to resolve count / range / named etc. without a reference lookup).
fn extract_literal_sequence(
    seq: &InsertedSequence,
    transcript: &Transcript,
) -> Result<String, FerroError> {
    match seq {
        InsertedSequence::Literal(s) => Ok(s.to_string()),
        _ => Err(FerroError::UnsupportedProjection {
            reason: format!(
                "protein prediction requires a literal inserted sequence; \
                 got non-literal for transcript {}",
                transcript.id
            ),
        }),
    }
}

// ── Translation helpers ───────────────────────────────────────────────────────

/// Translate a CDS sequence to a protein, stopping before the first stop codon.
///
/// Incomplete codons at the end are silently dropped.
pub(crate) fn translate_full_cds(cds: &str) -> Vec<AminoAcid> {
    cds.as_bytes()
        .chunks_exact(3)
        .filter_map(|c| std::str::from_utf8(c).ok().and_then(translate))
        .take_while(|aa| *aa != AminoAcid::Ter)
        .collect()
}

/// Translate a CDS sequence to a protein, including the termination codon.
///
/// Stops translation immediately after the first `Ter` is encountered.
pub(crate) fn translate_full_cds_with_stop(cds: &str) -> Vec<AminoAcid> {
    let mut result = Vec::new();
    for chunk in cds.as_bytes().chunks_exact(3) {
        if let Ok(s) = std::str::from_utf8(chunk) {
            if let Some(aa) = translate(s) {
                result.push(aa);
                if aa == AminoAcid::Ter {
                    break;
                }
            }
        }
    }
    result
}

// ── Diff helpers ─────────────────────────────────────────────────────────────

/// Return the 0-based index of the first amino acid that differs between `ref_prot`
/// and `alt_prot`.  If they are identical up to the shorter length, returns that length.
pub(crate) fn first_diff_position(ref_prot: &[AminoAcid], alt_prot: &[AminoAcid]) -> usize {
    ref_prot
        .iter()
        .zip(alt_prot.iter())
        .position(|(r, a)| r != a)
        .unwrap_or(ref_prot.len().min(alt_prot.len()))
}

/// Compute the net number of nucleotides added (positive) or removed (negative) by `edit`.
///
/// Returns `None` when the net change cannot be determined (e.g. non-literal inserted
/// sequences or uncertain-length edits).
pub(crate) fn net_length_change(edit: &NaEdit, del_len: usize) -> Option<i64> {
    match edit {
        NaEdit::Deletion { .. } => Some(-(del_len as i64)),
        NaEdit::Insertion { sequence } => sequence.len().map(|n| n as i64),
        NaEdit::Duplication { .. } => Some(del_len as i64),
        NaEdit::Delins { sequence, .. } => sequence.len().map(|n| n as i64 - del_len as i64),
        NaEdit::Inversion { .. } => Some(0),
        NaEdit::Substitution { .. } => Some(0),
        _ => None,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::edit::{InsertedSequence, Sequence};
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use std::sync::OnceLock;

    fn tx(seq: &str, cds_start: u64, cds_end: u64) -> Transcript {
        Transcript {
            id: "NM_TEST.1".to_string(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: Some(seq.to_string()),
            cds_start: Some(cds_start),
            cds_end: Some(cds_end),
            exons: vec![Exon::new(1, 1, seq.len() as u64)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        }
    }

    // ── read_full_cds ─────────────────────────────────────────────────────────

    #[test]
    fn read_full_cds_no_utr() {
        // cds_start=1, cds_end=9 (whole transcript = CDS "ATGCGCTAA")
        let t = tx("ATGCGCTAA", 1, 9);
        assert_eq!(read_full_cds(&t).unwrap(), "ATGCGCTAA");
    }

    #[test]
    fn read_full_cds_with_5utr() {
        // 2-base 5' UTR "AA" + CDS "ATGCCCTAG" (9 bp).
        let t = tx("AAATGCCCTAG", 3, 11);
        assert_eq!(read_full_cds(&t).unwrap(), "ATGCCCTAG");
    }

    // ── build_mutated_cds ─────────────────────────────────────────────────────

    #[test]
    fn mutated_cds_deletion_single_base() {
        // CDS "ATGCGCTAA"; delete c.4 (C, first base of Arg codon).
        // Mutated CDS = "ATG" + "GCTAA" = "ATGGCTAA"
        let t = tx("ATGCGCTAA", 1, 9);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = build_mutated_cds(&t, 4, 4, &edit).unwrap();
        assert_eq!(result, "ATGGCTAA");
    }

    #[test]
    fn mutated_cds_deletion_three_bases() {
        // CDS "ATGCGCTAA"; delete c.4_6 (CGC = Arg codon).
        // Mutated CDS = "ATG" + "TAA" = "ATGTAA"
        let t = tx("ATGCGCTAA", 1, 9);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = build_mutated_cds(&t, 4, 6, &edit).unwrap();
        assert_eq!(result, "ATGTAA");
    }

    #[test]
    fn mutated_cds_insertion() {
        // CDS "ATGCGCTAA"; insert "GGG" at c.3_4 (between positions 3 and 4).
        // build_mutated_cds receives cds_pos_start=3, cds_pos_end=3 (length=0 gap),
        // or start=3, end=3 giving affected=""...
        // Actually for HGVS insertion c.3_4ins the caller should pass start=3, end=3
        // so affected = "" and we insert between.
        let t = tx("ATGCGCTAA", 1, 9);
        let seq: Sequence = "GGG".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: InsertedSequence::Literal(seq),
        };
        // Insert at position end-of-codon1 / start-of-codon2: cds 3 to 3
        // before = "ATG", affected = "", insert = "GGG", after = "CGCTAA"
        let result = build_mutated_cds(&t, 3, 3, &edit).unwrap();
        assert_eq!(result, "ATGGGGCGCTAA");
    }

    #[test]
    fn mutated_cds_duplication() {
        // CDS "ATGCGCTAA"; dup c.4_6 (CGC).
        // Mutated CDS = "ATG" + "CGC" + "CGC" + "TAA" = "ATGCGCCGCTAA"
        let t = tx("ATGCGCTAA", 1, 9);
        let edit = NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: None,
        };
        let result = build_mutated_cds(&t, 4, 6, &edit).unwrap();
        assert_eq!(result, "ATGCGCCGCTAA");
    }

    #[test]
    fn mutated_cds_delins() {
        // CDS "ATGCGCTAA"; replace c.4_6 (CGC) with "TCC".
        // Mutated CDS = "ATG" + "TCC" + "TAA" = "ATGTCCTAA"
        let t = tx("ATGCGCTAA", 1, 9);
        let seq: Sequence = "TCC".parse().unwrap();
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(seq),
            deleted: None,
            deleted_length: None,
        };
        let result = build_mutated_cds(&t, 4, 6, &edit).unwrap();
        assert_eq!(result, "ATGTCCTAA");
    }

    #[test]
    fn mutated_cds_inversion() {
        // CDS "ATGCGCTAA"; invert c.4_6 (CGC → revcomp = GCG).
        // Mutated CDS = "ATG" + "GCG" + "TAA" = "ATGGCGTAA"
        let t = tx("ATGCGCTAA", 1, 9);
        let edit = NaEdit::Inversion {
            sequence: None,
            length: None,
        };
        let result = build_mutated_cds(&t, 4, 6, &edit).unwrap();
        assert_eq!(result, "ATGGCGTAA");
    }

    // ── translate helpers ─────────────────────────────────────────────────────

    #[test]
    fn translate_full_cds_met_arg() {
        // "ATGCGCTAA" → [Met, Arg] (stops before Ter)
        let aas = translate_full_cds("ATGCGCTAA");
        assert_eq!(aas, vec![AminoAcid::Met, AminoAcid::Arg]);
    }

    #[test]
    fn translate_full_cds_with_stop_includes_ter() {
        let aas = translate_full_cds_with_stop("ATGCGCTAA");
        assert_eq!(aas, vec![AminoAcid::Met, AminoAcid::Arg, AminoAcid::Ter]);
    }

    #[test]
    fn translate_full_cds_incomplete_codon_dropped() {
        // 7 bases → 2 complete codons + 1 incomplete (dropped).
        let aas = translate_full_cds("ATGCGCTA");
        assert_eq!(aas, vec![AminoAcid::Met, AminoAcid::Arg]);
    }

    // ── first_diff_position ───────────────────────────────────────────────────

    #[test]
    fn first_diff_identical() {
        let r = vec![AminoAcid::Met, AminoAcid::Arg];
        let a = vec![AminoAcid::Met, AminoAcid::Arg];
        assert_eq!(first_diff_position(&r, &a), 2);
    }

    #[test]
    fn first_diff_at_start() {
        let r = vec![AminoAcid::Met, AminoAcid::Arg];
        let a = vec![AminoAcid::Val, AminoAcid::Arg];
        assert_eq!(first_diff_position(&r, &a), 0);
    }

    #[test]
    fn first_diff_second_position() {
        let r = vec![AminoAcid::Met, AminoAcid::Arg];
        let a = vec![AminoAcid::Met, AminoAcid::Ser];
        assert_eq!(first_diff_position(&r, &a), 1);
    }

    // ── net_length_change ─────────────────────────────────────────────────────

    #[test]
    fn net_change_deletion() {
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        assert_eq!(net_length_change(&edit, 3), Some(-3));
    }

    #[test]
    fn net_change_insertion() {
        let seq: Sequence = "GGG".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: InsertedSequence::Literal(seq),
        };
        assert_eq!(net_length_change(&edit, 0), Some(3));
    }

    #[test]
    fn net_change_duplication() {
        let edit = NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: None,
        };
        assert_eq!(net_length_change(&edit, 3), Some(3));
    }

    #[test]
    fn net_change_inversion() {
        let edit = NaEdit::Inversion {
            sequence: None,
            length: None,
        };
        assert_eq!(net_length_change(&edit, 6), Some(0));
    }
}
