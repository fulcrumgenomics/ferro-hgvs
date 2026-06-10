//! The `c↔n` axis: derive the `n.` (transcript-relative) representation of a
//! coding (`c.`) variant, genome-free, via CDS-offset arithmetic.
//!
//! `c→n` is pure transcript-local arithmetic (CDS position → transcript
//! position through [`crate::convert::coding::cds_to_transcript_pos`]), so it is
//! defined even for intronic positions — unlike `c→p`, which needs the spliced
//! CDS.
#![forbid(unsafe_code)]

use crate::convert::coding::cds_to_transcript_pos;
use crate::error::FerroError;
use crate::hgvs::edit::NaEdit;
use crate::hgvs::interval::TxInterval;
use crate::hgvs::location::{CdsPos, TxPos};
use crate::hgvs::variant::{HgvsVariant, LocEdit, TxVariant};
use crate::project::accession::parse_accession;
use crate::reference::transcript::Transcript;

/// Convert one `c.` position to its `n.` (transcript-relative) [`TxPos`].
///
/// `n.base` is the transcript-frame position from `cds_to_transcript_pos`
/// (where `c.1 → cds_start`), and any intronic `CdsPos.offset` is carried over
/// verbatim. A coding 3'UTR base (`c.*N`) maps to a plain positive `n.`
/// coordinate (not `n.*N`), matching `cds_to_transcript_pos`'s `cds_end + N`.
fn cds_pos_to_tx_pos(pos: &CdsPos, transcript: &Transcript) -> Result<TxPos, FerroError> {
    let base = cds_to_transcript_pos(pos, transcript)? as i64;
    Ok(match pos.offset {
        Some(offset) => TxPos::with_offset(base, offset),
        None => TxPos::new(base),
    })
}

/// Derive the `n.` form ([`HgvsVariant::Tx`]) of a coding variant from its
/// resolved `c.` start/end CDS positions and transcript-axis edit.
///
/// The edit is identical to the one carried in the `c.` form — only the
/// coordinate frame reframes from CDS-relative to transcript-relative.
pub(crate) fn noncoding_from_coding(
    cds_start: &CdsPos,
    cds_end: &CdsPos,
    transcript: &Transcript,
    c_edit: &NaEdit,
    transcript_id: &str,
    gene_symbol: Option<String>,
) -> Result<HgvsVariant, FerroError> {
    let n_start = cds_pos_to_tx_pos(cds_start, transcript)?;
    let n_end = cds_pos_to_tx_pos(cds_end, transcript)?;
    Ok(HgvsVariant::Tx(TxVariant {
        accession: parse_accession(transcript_id),
        gene_symbol,
        loc_edit: LocEdit::new(TxInterval::new(n_start, n_end), c_edit.clone()),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::edit::{Base, NaEdit};
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use std::sync::OnceLock;

    fn make_test_transcript() -> Transcript {
        Transcript {
            id: "NM_TEST.1".to_string(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: Some("AAAAATGCCCAAAGGGTTTTAAAAAA".to_string()), // 26 bases
            cds_start: Some(6),
            cds_end: Some(20),
            exons: vec![Exon::new(1, 1, 26)],
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

    fn sub_a_g() -> NaEdit {
        // A>G substitution; exact bases are irrelevant to coordinate reframing.
        // `Base` is an enum (src/hgvs/edit.rs) — construct variants directly.
        NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        }
    }

    #[test]
    fn cds_first_base_maps_to_cds_start_offset() {
        // c.1 (first CDS base) → n.6 (cds_start), per cds_to_transcript_pos.
        let tx = make_test_transcript();
        let p = cds_pos_to_tx_pos(&CdsPos::new(1), &tx).unwrap();
        assert_eq!(p.base, 6);
        assert_eq!(p.offset, None);
        assert!(!p.downstream);
    }

    #[test]
    fn five_utr_maps_below_cds_start() {
        // c.-1 → n.5 (cds_start - 1).
        let tx = make_test_transcript();
        let p = cds_pos_to_tx_pos(&CdsPos::new(-1), &tx).unwrap();
        assert_eq!(p.base, 5);
    }

    #[test]
    fn three_utr_maps_past_cds_end_as_plain_positive() {
        // c.*1 → n.21 (cds_end + 1), plain positive (NOT n.* notation).
        let tx = make_test_transcript();
        let p = cds_pos_to_tx_pos(&CdsPos::utr3(1), &tx).unwrap();
        assert_eq!(p.base, 21);
        assert!(!p.downstream);
    }

    #[test]
    fn intronic_offset_is_carried() {
        // c.3+5 → n.8+5 : base from cds_to_transcript_pos(c.3)=8, offset +5 kept.
        let tx = make_test_transcript();
        let mut pos = CdsPos::new(3);
        pos.offset = Some(5);
        let p = cds_pos_to_tx_pos(&pos, &tx).unwrap();
        assert_eq!(p.base, 8);
        assert_eq!(p.offset, Some(5));
    }

    #[test]
    fn builds_tx_variant_rendering_as_n_dot() {
        // c.1A>G on NM_TEST.1 → n.6A>G.
        let tx = make_test_transcript();
        let v = noncoding_from_coding(
            &CdsPos::new(1),
            &CdsPos::new(1),
            &tx,
            &sub_a_g(),
            "NM_TEST.1",
            None,
        )
        .unwrap();
        assert!(matches!(v, HgvsVariant::Tx(_)));
        assert_eq!(v.to_string(), "NM_TEST.1:n.6A>G");
    }
}
