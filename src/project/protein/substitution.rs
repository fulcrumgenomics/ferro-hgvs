//! Protein consequence prediction for CDS substitutions.

use crate::backtranslate::{Codon, CodonTable};
use crate::error::FerroError;
use crate::hgvs::edit::{Base, NaEdit, ProteinEdit};
use crate::hgvs::interval::ProtInterval;
use crate::hgvs::location::{AminoAcid, ProtPos};
use crate::hgvs::variant::{HgvsVariant, LocEdit, ProteinVariant};
use crate::project::accession::parse_accession;
use crate::reference::transcript::Transcript;

/// Read the codon (3 bases) covering a given CDS position from a transcript's CDS sequence.
///
/// Returns the codon bases (uppercase) plus the 0/1/2 frame index (which base
/// within the codon the CDS position corresponds to).
pub(crate) fn read_ref_codon(
    transcript: &Transcript,
    cds_pos: i64,
) -> Result<(String, u8), FerroError> {
    let cds_start = transcript
        .cds_start
        .ok_or_else(|| FerroError::ConversionError {
            msg: format!("transcript {} has no CDS", transcript.id),
        })?;
    let cds_end = transcript
        .cds_end
        .ok_or_else(|| FerroError::ConversionError {
            msg: format!("transcript {} has no CDS", transcript.id),
        })?;
    if cds_pos < 1 {
        return Err(FerroError::ConversionError {
            msg: format!("CDS position {} is not in the CDS", cds_pos),
        });
    }
    // Transcript.cds_start/cds_end are 1-based inclusive transcript positions
    // (see `src/reference/transcript.rs`), so the CDS spans cds_end - cds_start + 1
    // bases and cds_pos is bounded by that length.
    let cds_len = cds_end.saturating_sub(cds_start) + 1;
    if (cds_pos as u64) > cds_len {
        return Err(FerroError::ConversionError {
            msg: format!(
                "CDS position {} is outside the CDS of transcript {} (length {})",
                cds_pos, transcript.id, cds_len
            ),
        });
    }
    let codon_index = ((cds_pos - 1) / 3) as u64;
    let frame = ((cds_pos - 1) % 3) as u8;

    let tx_offset = (cds_start as i64 - 1 + codon_index as i64 * 3) as usize;
    let seq =
        transcript
            .sequence
            .as_deref()
            .ok_or_else(|| FerroError::ProteinSequenceUnavailable {
                accession: transcript.id.clone(),
            })?;
    if tx_offset + 3 > seq.len() {
        return Err(FerroError::ProteinSequenceUnavailable {
            accession: transcript.id.clone(),
        });
    }
    let codon = seq[tx_offset..tx_offset + 3].to_uppercase();
    Ok((codon, frame))
}

/// Apply a single-base substitution to a codon at the given 0/1/2 frame index.
pub(crate) fn apply_substitution(codon: &str, frame: u8, alt: Base) -> String {
    let mut bytes = codon.as_bytes().to_vec();
    bytes[frame as usize] = alt.to_u8();
    String::from_utf8(bytes).expect("alt base is ASCII")
}

/// Translate a 3-character codon to an `AminoAcid`. Returns `Some(Ter)` for stop codons,
/// `Some(Xaa)` if the codon is unrecognized; `None` only if the input is not 3 ASCII bases.
pub(crate) fn translate(codon: &str) -> Option<AminoAcid> {
    let parsed = Codon::parse(codon)?;
    let table = CodonTable::standard();
    if table.is_stop(&parsed) {
        return Some(AminoAcid::Ter);
    }
    Some(*table.amino_acid_for(&parsed).unwrap_or(&AminoAcid::Xaa))
}

/// Build a predicted `p.(...)` ProteinVariant for a CDS substitution.
///
/// Returns `Err(ProteinSequenceUnavailable)` if the transcript's CDS sequence
/// cannot supply the codon at `cds_pos`. Returns `Err(UnsupportedProjection)`
/// if `edit` is not a `NaEdit::Substitution`.
pub(crate) fn predict_substitution_protein(
    transcript: &Transcript,
    cds_pos: i64,
    edit: &NaEdit,
    protein_accession: &str,
) -> Result<HgvsVariant, FerroError> {
    let alt = match edit {
        NaEdit::Substitution { alternative, .. } => *alternative,
        _ => {
            return Err(FerroError::UnsupportedProjection {
                reason: "predict_substitution_protein called with non-substitution edit"
                    .to_string(),
            })
        }
    };
    let (ref_codon, frame) = read_ref_codon(transcript, cds_pos)?;
    let alt_codon = apply_substitution(&ref_codon, frame, alt);
    let ref_aa = translate(&ref_codon).unwrap_or(AminoAcid::Xaa);
    let alt_aa = translate(&alt_codon).unwrap_or(AminoAcid::Xaa);
    let aa_number = ((cds_pos - 1) / 3 + 1) as u64;

    // For the identity case, use predicted: false — the outer Mu::Uncertain wrapper
    // (from LocEdit::new_predicted) already provides the p.(...) predicted wrapping.
    // Using predicted: true here would produce double parens: p.(Arg1(=)).
    let protein_edit = if ref_aa == alt_aa {
        ProteinEdit::Identity {
            predicted: false,
            whole_protein: false,
        }
    } else {
        ProteinEdit::Substitution {
            reference: ref_aa,
            alternative: alt_aa,
        }
    };

    let accession = parse_accession(protein_accession);
    let loc = ProtInterval::point(ProtPos::new(ref_aa, aa_number));
    let variant = ProteinVariant {
        accession,
        gene_symbol: transcript.gene_symbol.clone(),
        loc_edit: LocEdit::new_predicted(loc, protein_edit),
    };
    Ok(HgvsVariant::Protein(variant))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use std::sync::OnceLock;

    fn tx_with_seq(seq: &str, cds_start: u64, cds_end: u64) -> Transcript {
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

    #[test]
    fn read_ref_codon_first_codon_first_base() {
        // 5' UTR "AA" + CDS "ATGCCCTAG" (Met-Pro-Stop). cds_start=3 (1-based).
        let tx = tx_with_seq("AAATGCCCTAG", 3, 11);
        let (codon, frame) = read_ref_codon(&tx, 1).unwrap();
        assert_eq!(codon, "ATG");
        assert_eq!(frame, 0);
    }

    #[test]
    fn read_ref_codon_second_codon_middle_base() {
        let tx = tx_with_seq("AAATGCCCTAG", 3, 11);
        let (codon, frame) = read_ref_codon(&tx, 5).unwrap();
        assert_eq!(codon, "CCC");
        assert_eq!(frame, 1);
    }

    #[test]
    fn read_ref_codon_rejects_pos_below_one() {
        let tx = tx_with_seq("AAATGCCCTAG", 3, 11);
        let err = read_ref_codon(&tx, 0).expect_err("cds_pos < 1 must error");
        assert!(matches!(err, FerroError::ConversionError { .. }), "{err:?}");
    }

    #[test]
    fn read_ref_codon_rejects_pos_past_cds_end() {
        // CDS spans c.1..c.9 (cds_start=3, cds_end=11; CDS length = 11 - 3 + 1 = 9).
        let tx = tx_with_seq("AAATGCCCTAG", 3, 11);
        // c.9 is the last CDS position and still valid.
        assert!(read_ref_codon(&tx, 9).is_ok());
        // c.10 is past cds_end → must error rather than read 3'UTR bytes.
        let err = read_ref_codon(&tx, 10).expect_err("cds_pos past cds_end must error");
        match &err {
            FerroError::ConversionError { msg } => {
                assert!(
                    msg.contains("outside the CDS"),
                    "expected message to mention 'outside the CDS', got: {msg}"
                );
            }
            other => panic!("expected ConversionError, got: {other:?}"),
        }
    }

    #[test]
    fn apply_substitution_first_base() {
        assert_eq!(apply_substitution("ATG", 0, Base::C), "CTG");
    }

    #[test]
    fn apply_substitution_middle_base() {
        assert_eq!(apply_substitution("ATG", 1, Base::A), "AAG");
    }

    #[test]
    fn translate_known_codons() {
        assert_eq!(translate("ATG"), Some(AminoAcid::Met));
        assert_eq!(translate("TAA"), Some(AminoAcid::Ter));
        assert_eq!(translate("CGG"), Some(AminoAcid::Arg));
        assert_eq!(translate("GGG"), Some(AminoAcid::Gly));
    }

    #[test]
    fn substitution_missense_arg_to_gly() {
        // CDS "CGGGGG": codon 1 = CGG = Arg. c.1C>G → codon 1 = GGG = Gly.
        let tx = tx_with_seq("CGGGGG", 1, 6);
        let edit = NaEdit::Substitution {
            reference: Base::C,
            alternative: Base::G,
        };
        let pv = predict_substitution_protein(&tx, 1, &edit, "NP_000288.1").unwrap();
        let s = match &pv {
            HgvsVariant::Protein(p) => p.to_string(),
            _ => panic!("expected Protein variant"),
        };
        assert_eq!(s, "NP_000288.1:p.(Arg1Gly)");
    }

    #[test]
    fn substitution_synonymous_arg_arg() {
        // CGG (Arg) → CGC (Arg). Frame 2 of codon 1, alt C.
        let tx = tx_with_seq("CGGGGG", 1, 6);
        let edit = NaEdit::Substitution {
            reference: Base::G,
            alternative: Base::C,
        };
        let pv = predict_substitution_protein(&tx, 3, &edit, "NP_000288.1").unwrap();
        let s = match &pv {
            HgvsVariant::Protein(p) => p.to_string(),
            _ => panic!("expected Protein variant"),
        };
        assert_eq!(s, "NP_000288.1:p.(Arg1=)");
    }

    #[test]
    fn substitution_nonsense_arg_ter() {
        // CGA (Arg) → TGA (Ter). Frame 0 of codon 1, alt T.
        let tx = tx_with_seq("CGAAAA", 1, 6);
        let edit = NaEdit::Substitution {
            reference: Base::C,
            alternative: Base::T,
        };
        let pv = predict_substitution_protein(&tx, 1, &edit, "NP_000288.1").unwrap();
        let s = match &pv {
            HgvsVariant::Protein(p) => p.to_string(),
            _ => panic!("expected Protein variant"),
        };
        assert_eq!(s, "NP_000288.1:p.(Arg1Ter)");
    }
}
