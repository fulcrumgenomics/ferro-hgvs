//! The `c↔n` axis: derive the `n.` (transcript-relative) representation of a
//! coding (`c.`) variant.
//!
//! The CDS→transcript conversion is routed through the projector's authoritative
//! [`CoordinateMapper::cds_to_tx`] — the **same** exon- and CIGAR-aware mapping
//! that `genome_to_cds` / `project_single_inner` already trust. A flat
//! CDS-offset shift (`cds_start + base - 1`) would silently diverge on
//! multi-exon transcripts with tx-coordinate gaps or on transcripts whose
//! genome alignment carries CIGAR insertions, producing an `n.` base that
//! disagrees with the `c.` form's true transcript position. Routing through the
//! mapper keeps the two axes consistent and makes `c→n` defined even for
//! intronic positions (the mapper carries the offset).
#![forbid(unsafe_code)]

use crate::convert::mapper::CoordinateMapper;
use crate::error::FerroError;
use crate::hgvs::edit::NaEdit;
use crate::hgvs::interval::TxInterval;
use crate::hgvs::location::CdsPos;
use crate::hgvs::variant::{HgvsVariant, LocEdit, TxVariant};
use crate::project::accession::parse_accession;
use crate::reference::transcript::Transcript;

/// Derive the `n.` form ([`HgvsVariant::Tx`]) of a coding variant from its
/// resolved `c.` start/end CDS positions and transcript-axis edit.
///
/// The edit is identical to the one carried in the `c.` form — only the
/// coordinate frame reframes from CDS-relative to transcript-relative, via the
/// exon/CIGAR-aware [`CoordinateMapper::cds_to_tx`]. Intronic offsets and the
/// 3'UTR (`c.*N` → a plain positive `n.`, not `n.*N`) are handled by the mapper.
pub(crate) fn noncoding_from_coding(
    cds_start: &CdsPos,
    cds_end: &CdsPos,
    transcript: &Transcript,
    c_edit: &NaEdit,
    transcript_id: &str,
    gene_symbol: Option<String>,
) -> Result<HgvsVariant, FerroError> {
    let mapper = CoordinateMapper::new(transcript);
    let n_start = mapper.cds_to_tx(cds_start)?;
    let n_end = mapper.cds_to_tx(cds_end)?;
    Ok(HgvsVariant::Tx(TxVariant {
        accession: parse_accession(transcript_id),
        gene_symbol,
        loc_edit: LocEdit::new(TxInterval::new(n_start, n_end), c_edit.clone()),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::cdot::CigarOp;
    use crate::hgvs::edit::Base;
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use std::sync::OnceLock;

    /// Single-exon coding transcript: `cds_start = 6`, `cds_end = 20`, no CIGAR.
    /// `n.` coordinates equal plain transcript offsets and are hand-checkable.
    fn make_test_transcript() -> Transcript {
        Transcript {
            cds_start_incomplete: false,
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

    /// Render the `n.` form of a single-position coding variant for assertions.
    fn n_string(cds: CdsPos, tx: &Transcript, id: &str) -> String {
        noncoding_from_coding(&cds, &cds, tx, &sub_a_g(), id, None)
            .unwrap()
            .to_string()
    }

    #[test]
    fn cds_first_base_renders_as_cds_start_offset() {
        // c.1 (first CDS base) → n.6 (cds_start).
        let tx = make_test_transcript();
        assert_eq!(
            n_string(CdsPos::new(1), &tx, "NM_TEST.1"),
            "NM_TEST.1:n.6A>G"
        );
    }

    #[test]
    fn five_utr_renders_below_cds_start() {
        // c.-1 → n.5 (cds_start - 1).
        let tx = make_test_transcript();
        assert_eq!(
            n_string(CdsPos::new(-1), &tx, "NM_TEST.1"),
            "NM_TEST.1:n.5A>G"
        );
    }

    #[test]
    fn three_utr_renders_plain_positive_not_star() {
        // c.*1 → n.21 (cds_end + 1), plain positive — NOT n.*1 (which would mean
        // past the transcript 3' end).
        let tx = make_test_transcript();
        assert_eq!(
            n_string(CdsPos::utr3(1), &tx, "NM_TEST.1"),
            "NM_TEST.1:n.21A>G"
        );
    }

    #[test]
    fn intronic_offset_is_carried() {
        // c.3+5 → n.8+5 : base 8 (cds_start 6 + 2), intronic offset +5 carried.
        let tx = make_test_transcript();
        let mut pos = CdsPos::new(3);
        pos.offset = Some(5);
        assert_eq!(n_string(pos, &tx, "NM_TEST.1"), "NM_TEST.1:n.8+5A>G");
    }

    #[test]
    fn builds_tx_variant_not_cds() {
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
        assert!(
            matches!(v, HgvsVariant::Tx(_)),
            "expected an n. (Tx) variant"
        );
    }

    #[test]
    fn cigar_insertion_n_coordinate_is_transcript_native() {
        // #944: HGVS c./n. numbering is transcript-native — a CIGAR insertion adds
        // transcript bases that BOTH the CDS and n. coordinate already count, so
        // c.6 maps to n.6 with NO shift (the 2 inserted bases are c.4/c.5). The n.
        // derivation must still route through the same CIGAR-aware `cds_to_tx` as
        // the c. form (c/n consistency — this test's original #592 intent); after
        // the #944 fix both agree at n.6.
        let tx = Transcript {
            cds_start_incomplete: false,
            id: "NM_CIGAR.1".to_string(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: Some("A".repeat(20)),
            cds_start: Some(1),
            cds_end: Some(20),
            exons: vec![Exon::new(1, 1, 20)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            // 3 matched, 2 inserted (tx bases not in genome), then 15 matched.
            exon_cigars: vec![Some(vec![
                CigarOp::Match(3),
                CigarOp::Insertion(2),
                CigarOp::Match(15),
            ])],
            cached_introns: OnceLock::new(),
        };
        // Transcript-native: c.6 -> n.6.
        assert_eq!(
            n_string(CdsPos::new(6), &tx, "NM_CIGAR.1"),
            "NM_CIGAR.1:n.6A>G"
        );
        // Consistency: the mapper agrees.
        let mapper = CoordinateMapper::new(&tx);
        assert_eq!(mapper.cds_to_tx(&CdsPos::new(6)).unwrap().base, 6);
    }
}
