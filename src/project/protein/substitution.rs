//! Protein consequence prediction for CDS substitutions.

use crate::backtranslate::{Codon, CodonTable};
use crate::error::FerroError;
use crate::hgvs::edit::{Base, NaEdit, ProteinEdit};
use crate::hgvs::interval::ProtInterval;
use crate::hgvs::location::{AminoAcid, ProtPos};
use crate::hgvs::variant::{HgvsVariant, LocEdit, ProteinVariant};
use crate::project::accession::parse_accession;
use crate::project::protein::helpers::{affects_initiation_codon, build_initiator_unknown};
use crate::reference::transcript::Transcript;
use once_cell::sync::Lazy;

/// Standard genetic-code table, built once and reused. Constructing a fresh
/// `CodonTable` per call dominates protein-prediction CPU because
/// `translate_full_cds` invokes `translate` once per codon of every CDS.
static STANDARD_CODON_TABLE: Lazy<CodonTable> = Lazy::new(CodonTable::standard);

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
    translate_bytes(codon.as_bytes())
}

/// Byte-slice variant of [`translate`] for hot callers that have already
/// produced a `&[u8]` (e.g. `chunks_exact(3)` over a CDS sequence). Skips the
/// `from_utf8` step that the `&str`-taking entry point would imply for
/// `chunks_exact` output. samply showed that detour at ~18% of indel CPU
/// after the prior round of perf wins.
pub(crate) fn translate_bytes(bytes: &[u8]) -> Option<AminoAcid> {
    let parsed = Codon::parse_bytes(bytes)?;
    let table = &*STANDARD_CODON_TABLE;
    if table.is_stop(&parsed) {
        return Some(AminoAcid::Ter);
    }
    Some(table.amino_acid_for(&parsed).unwrap_or(AminoAcid::Xaa))
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
    // A substitution in the translation initiation codon (CDS 1–3) has an
    // unpredictable protein consequence — report p.(Met1?) rather than a
    // concrete missense, which the spec disallows (#498).
    if affects_initiation_codon(edit, cds_pos, cds_pos) {
        return Ok(build_initiator_unknown(protein_accession, transcript));
    }
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
            protein_id: None,
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

    /// Exhaustive: every 3-base ACGT codon must round-trip to the canonical
    /// amino acid. Used as the safety net for perf changes in this file —
    /// callers depend on these mappings being byte-stable across releases.
    #[test]
    fn translate_all_64_codons() {
        let table: &[(&str, AminoAcid)] = &[
            ("TTT", AminoAcid::Phe),
            ("TTC", AminoAcid::Phe),
            ("TTA", AminoAcid::Leu),
            ("TTG", AminoAcid::Leu),
            ("CTT", AminoAcid::Leu),
            ("CTC", AminoAcid::Leu),
            ("CTA", AminoAcid::Leu),
            ("CTG", AminoAcid::Leu),
            ("ATT", AminoAcid::Ile),
            ("ATC", AminoAcid::Ile),
            ("ATA", AminoAcid::Ile),
            ("ATG", AminoAcid::Met),
            ("GTT", AminoAcid::Val),
            ("GTC", AminoAcid::Val),
            ("GTA", AminoAcid::Val),
            ("GTG", AminoAcid::Val),
            ("TCT", AminoAcid::Ser),
            ("TCC", AminoAcid::Ser),
            ("TCA", AminoAcid::Ser),
            ("TCG", AminoAcid::Ser),
            ("AGT", AminoAcid::Ser),
            ("AGC", AminoAcid::Ser),
            ("CCT", AminoAcid::Pro),
            ("CCC", AminoAcid::Pro),
            ("CCA", AminoAcid::Pro),
            ("CCG", AminoAcid::Pro),
            ("ACT", AminoAcid::Thr),
            ("ACC", AminoAcid::Thr),
            ("ACA", AminoAcid::Thr),
            ("ACG", AminoAcid::Thr),
            ("GCT", AminoAcid::Ala),
            ("GCC", AminoAcid::Ala),
            ("GCA", AminoAcid::Ala),
            ("GCG", AminoAcid::Ala),
            ("TAT", AminoAcid::Tyr),
            ("TAC", AminoAcid::Tyr),
            ("TAA", AminoAcid::Ter),
            ("TAG", AminoAcid::Ter),
            ("TGA", AminoAcid::Ter),
            ("CAT", AminoAcid::His),
            ("CAC", AminoAcid::His),
            ("CAA", AminoAcid::Gln),
            ("CAG", AminoAcid::Gln),
            ("AAT", AminoAcid::Asn),
            ("AAC", AminoAcid::Asn),
            ("AAA", AminoAcid::Lys),
            ("AAG", AminoAcid::Lys),
            ("GAT", AminoAcid::Asp),
            ("GAC", AminoAcid::Asp),
            ("GAA", AminoAcid::Glu),
            ("GAG", AminoAcid::Glu),
            ("TGT", AminoAcid::Cys),
            ("TGC", AminoAcid::Cys),
            ("TGG", AminoAcid::Trp),
            ("CGT", AminoAcid::Arg),
            ("CGC", AminoAcid::Arg),
            ("CGA", AminoAcid::Arg),
            ("CGG", AminoAcid::Arg),
            ("AGA", AminoAcid::Arg),
            ("AGG", AminoAcid::Arg),
            ("GGT", AminoAcid::Gly),
            ("GGC", AminoAcid::Gly),
            ("GGA", AminoAcid::Gly),
            ("GGG", AminoAcid::Gly),
        ];
        assert_eq!(table.len(), 64, "must enumerate all 64 codons");
        for (codon, aa) in table {
            assert_eq!(translate(codon), Some(*aa), "codon {codon}");
            // Lowercase and U-substitution must yield the same result.
            assert_eq!(
                translate(&codon.to_ascii_lowercase()),
                Some(*aa),
                "lowercase codon {codon}"
            );
            let u_codon: String = codon
                .chars()
                .map(|c| if c == 'T' { 'U' } else { c })
                .collect();
            assert_eq!(translate(&u_codon), Some(*aa), "RNA codon {u_codon}");
        }
    }

    #[test]
    fn translate_rejects_bad_input() {
        assert_eq!(translate(""), None, "empty");
        assert_eq!(translate("AT"), None, "two bases");
        assert_eq!(translate("ATGC"), None, "four bases");
        assert_eq!(
            translate("ATN"),
            None,
            "ambiguous N must not silently translate"
        );
        assert_eq!(translate("AT*"), None, "non-base char");
    }

    #[test]
    fn substitution_missense_arg_to_gly() {
        // CDS "ATGCGGGGG": codon 2 = CGG = Arg. c.4C>G → codon 2 = GGG = Gly.
        // (Codon 2, not 1, so the initiation-codon guard does not apply.)
        let tx = tx_with_seq("ATGCGGGGG", 1, 9);
        let edit = NaEdit::Substitution {
            reference: Base::C,
            alternative: Base::G,
        };
        let pv = predict_substitution_protein(&tx, 4, &edit, "NP_000288.1").unwrap();
        let s = match &pv {
            HgvsVariant::Protein(p) => p.to_string(),
            _ => panic!("expected Protein variant"),
        };
        assert_eq!(s, "NP_000288.1:p.(Arg2Gly)");
    }

    #[test]
    fn substitution_synonymous_arg_arg() {
        // codon 2 = CGG (Arg) → CGC (Arg). Frame 2 of codon 2 (c.6), alt C.
        let tx = tx_with_seq("ATGCGGGGG", 1, 9);
        let edit = NaEdit::Substitution {
            reference: Base::G,
            alternative: Base::C,
        };
        let pv = predict_substitution_protein(&tx, 6, &edit, "NP_000288.1").unwrap();
        let s = match &pv {
            HgvsVariant::Protein(p) => p.to_string(),
            _ => panic!("expected Protein variant"),
        };
        assert_eq!(s, "NP_000288.1:p.(Arg2=)");
    }

    #[test]
    fn substitution_in_initiation_codon_is_met1_unknown() {
        // A substitution anywhere in the initiation codon (CDS 1–3) has an
        // unpredictable protein consequence → p.(Met1?), NOT a concrete
        // missense (HGVS substitution.md:51; #498). Sequence codon 1 = CTG
        // (Leu) — protein residue 1 is still the initiator Met.
        let tx = tx_with_seq("CTGGGGTAA", 1, 9);
        // Reference base differs per position in codon CTG, so build a
        // position-specific substitution for each rather than reusing one
        // edit — every case is a true substitution against the actual base.
        for (pos, reference, alternative) in [
            (1, Base::C, Base::G),
            (2, Base::T, Base::G),
            (3, Base::G, Base::A),
        ] {
            let edit = NaEdit::Substitution {
                reference,
                alternative,
            };
            let pv = predict_substitution_protein(&tx, pos, &edit, "NP_TEST.1").unwrap();
            let s = match &pv {
                HgvsVariant::Protein(p) => p.to_string(),
                _ => panic!("expected Protein variant"),
            };
            assert_eq!(s, "NP_TEST.1:p.(Met1?)", "cds_pos {pos}");
        }
    }

    #[test]
    fn substitution_nonsense_arg_ter() {
        // codon 2 = CGA (Arg) → TGA (Ter). Frame 0 of codon 2 (c.4), alt T.
        let tx = tx_with_seq("ATGCGAAAA", 1, 9);
        let edit = NaEdit::Substitution {
            reference: Base::C,
            alternative: Base::T,
        };
        let pv = predict_substitution_protein(&tx, 4, &edit, "NP_000288.1").unwrap();
        let s = match &pv {
            HgvsVariant::Protein(p) => p.to_string(),
            _ => panic!("expected Protein variant"),
        };
        assert_eq!(s, "NP_000288.1:p.(Arg2Ter)");
    }
}
