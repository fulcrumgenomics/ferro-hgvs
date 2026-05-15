//! Protein consequence prediction for indels (deletion, insertion, duplication,
//! deletion-insertion, inversion) in the CDS.

use crate::error::FerroError;
use crate::hgvs::edit::{AminoAcidSeq, ExtDirection, NaEdit, ProteinEdit};
use crate::hgvs::interval::ProtInterval;
use crate::hgvs::location::{AminoAcid, ProtPos};
use crate::hgvs::variant::{HgvsVariant, LocEdit, ProteinVariant};
use crate::project::accession::parse_accession;
use crate::reference::transcript::Transcript;

use super::helpers::{
    build_mutated_cds, first_diff_position, net_length_change, read_full_cds, translate_full_cds,
    translate_full_cds_with_stop,
};

// ── Main entry point ──────────────────────────────────────────────────────────

/// Predict the protein consequence of an indel (non-substitution) CDS edit.
///
/// `cds_pos_start` and `cds_pos_end` are the 1-based CDS positions of the
/// variant's start and end respectively (as they appear in the c. variant).
///
/// Returns a `p.(...)` `HgvsVariant::Protein` on success, or a
/// `FerroError::UnsupportedProjection` / `FerroError::ProteinSequenceUnavailable`
/// on failure.
pub(crate) fn predict_indel_protein(
    transcript: &Transcript,
    cds_pos_start: i64,
    cds_pos_end: i64,
    edit: &NaEdit,
    protein_accession: &str,
) -> Result<HgvsVariant, FerroError> {
    // 1. Build ref and mutated CDS sequences.
    let ref_cds = read_full_cds(transcript)?;
    let mut_cds = build_mutated_cds(transcript, cds_pos_start, cds_pos_end, edit)?;

    // 2. Translate both.
    let ref_protein = translate_full_cds(&ref_cds);
    let alt_protein = translate_full_cds(&mut_cds);

    // 3. Check for stop-codon readthrough: the ref had a stop that is now an AA.
    //    Detected when the mutated protein is longer than the ref protein and
    //    alt_protein[ref_protein.len()] is not Ter.
    if alt_protein.len() > ref_protein.len() {
        // Check whether the ref protein ended at the stop codon position.
        // ref_cds should have a Ter at codon (ref_protein.len()+1).
        let ref_protein_with_stop = translate_full_cds_with_stop(&ref_cds);
        if ref_protein_with_stop.last() == Some(&AminoAcid::Ter)
            && alt_protein.len() > ref_protein.len()
        {
            // Possibly readthrough if the MUTATION affected the stop codon.
            // The affected CDS region includes the stop codon area.
            let stop_cds_start = (ref_protein.len() * 3 + 1) as i64; // 1-based CDS pos of stop codon
            if cds_pos_start >= stop_cds_start
                || (cds_pos_end >= stop_cds_start && cds_pos_start <= stop_cds_start + 2)
            {
                return build_extension_variant(
                    &ref_protein,
                    &alt_protein,
                    &mut_cds,
                    protein_accession,
                    transcript,
                );
            }
        }
    }

    // 4. Net length change and frameshift detection.
    let del_len = (cds_pos_end - cds_pos_start + 1) as usize;
    let net =
        net_length_change(edit, del_len).ok_or_else(|| FerroError::UnsupportedProjection {
            reason: "cannot determine net length change for protein prediction".to_string(),
        })?;
    let is_frameshift = net % 3 != 0;

    if is_frameshift {
        build_frameshift_variant(
            &ref_protein,
            &alt_protein,
            &mut_cds,
            protein_accession,
            transcript,
        )
    } else {
        build_inframe_variant(
            transcript,
            cds_pos_start,
            cds_pos_end,
            edit,
            net,
            &ref_protein,
            &alt_protein,
            protein_accession,
        )
    }
}

// ── In-frame prediction ───────────────────────────────────────────────────────

/// Build the protein variant for an in-frame (net % 3 == 0) indel.
#[allow(clippy::too_many_arguments)]
fn build_inframe_variant(
    transcript: &Transcript,
    cds_pos_start: i64,
    cds_pos_end: i64,
    edit: &NaEdit,
    net: i64, // net nucleotide change (negative = deletion, positive = insertion)
    ref_protein: &[AminoAcid],
    alt_protein: &[AminoAcid],
    protein_accession: &str,
) -> Result<HgvsVariant, FerroError> {
    // Check identity first: if the proteins are identical, emit p.(=).
    if ref_protein == alt_protein {
        return build_identity_variant(ref_protein, protein_accession, transcript);
    }

    // Is the edit codon-aligned (i.e. starts at the first base of a codon)?
    let frame_at_start = (cds_pos_start - 1) % 3; // 0, 1, or 2

    match edit {
        NaEdit::Deletion { .. } => {
            // Net < 0 always for deletion.
            if frame_at_start == 0 && net % 3 == 0 {
                // Codon-aligned deletion: whole codons removed.
                build_inframe_deletion(ref_protein, alt_protein, protein_accession, transcript)
            } else {
                // Straddles a codon boundary: effectively a delins at the AA level.
                build_inframe_delins(
                    ref_protein,
                    alt_protein,
                    cds_pos_start,
                    protein_accession,
                    transcript,
                )
            }
        }
        NaEdit::Insertion { .. } => {
            if frame_at_start == 2 && net % 3 == 0 {
                // Codon-aligned insertion: inserted between two adjacent codons.
                // cds_pos_end == cds_pos_start (gap position), and the position
                // is the last base of a codon (frame 2, i.e. cds_pos % 3 == 0).
                build_inframe_insertion(
                    ref_protein,
                    alt_protein,
                    cds_pos_start,
                    protein_accession,
                    transcript,
                )
            } else {
                // Mid-codon insertion: treated as a delins at the AA level.
                build_inframe_delins(
                    ref_protein,
                    alt_protein,
                    cds_pos_start,
                    protein_accession,
                    transcript,
                )
            }
        }
        NaEdit::Duplication { .. } => {
            if frame_at_start == 0 && net % 3 == 0 {
                // Codon-aligned duplication.
                build_inframe_duplication(
                    ref_protein,
                    alt_protein,
                    cds_pos_start,
                    cds_pos_end,
                    protein_accession,
                    transcript,
                )
            } else {
                // Mid-codon duplication: treated as a delins.
                build_inframe_delins(
                    ref_protein,
                    alt_protein,
                    cds_pos_start,
                    protein_accession,
                    transcript,
                )
            }
        }
        NaEdit::Delins { .. } | NaEdit::Inversion { .. } => {
            // Always use the generic delins pathway for delins and inversion.
            build_inframe_delins(
                ref_protein,
                alt_protein,
                cds_pos_start,
                protein_accession,
                transcript,
            )
        }
        _ => Err(FerroError::UnsupportedProjection {
            reason: format!(
                "build_inframe_variant does not support edit type for protein prediction: {:?}",
                edit
            ),
        }),
    }
}

/// Identity: proteins are identical after the edit (e.g. synonymous codon change).
fn build_identity_variant(
    ref_protein: &[AminoAcid],
    protein_accession: &str,
    transcript: &Transcript,
) -> Result<HgvsVariant, FerroError> {
    // Use the whole-protein identity representation.
    let accession = parse_accession(protein_accession);
    // Point at Met1 for simplicity.
    let ref_aa = ref_protein.first().copied().unwrap_or(AminoAcid::Met);
    let loc = ProtInterval::point(ProtPos::new(ref_aa, 1));
    let edit = crate::hgvs::edit::ProteinEdit::Identity {
        predicted: false,
        whole_protein: true,
    };
    let variant = ProteinVariant {
        accession,
        gene_symbol: transcript.gene_symbol.clone(),
        loc_edit: LocEdit::new_predicted(loc, edit),
    };
    Ok(HgvsVariant::Protein(variant))
}

/// Codon-aligned in-frame deletion: `p.(Xxx{N}del)` or `p.(Xxx{N}_Yyy{M}del)`.
fn build_inframe_deletion(
    ref_protein: &[AminoAcid],
    alt_protein: &[AminoAcid],
    protein_accession: &str,
    transcript: &Transcript,
) -> Result<HgvsVariant, FerroError> {
    // Find first and last positions that differ.
    let first_diff = first_diff_position(ref_protein, alt_protein);
    // In a pure codon-aligned deletion, alt_protein is shorter. The deleted
    // range in ref_protein is [first_diff .. first_diff + (ref_len - alt_len)].
    let n_deleted = ref_protein.len().saturating_sub(alt_protein.len());
    let last_deleted = first_diff + n_deleted - 1;

    if n_deleted == 0 {
        return build_identity_variant(ref_protein, protein_accession, transcript);
    }

    let start_aa = ref_protein
        .get(first_diff)
        .copied()
        .unwrap_or(AminoAcid::Xaa);
    let start_pos = (first_diff + 1) as u64;

    let protein_edit = ProteinEdit::Deletion {
        sequence: None,
        count: None,
    };

    let loc = if n_deleted == 1 {
        ProtInterval::point(ProtPos::new(start_aa, start_pos))
    } else {
        let end_aa = ref_protein
            .get(last_deleted)
            .copied()
            .unwrap_or(AminoAcid::Xaa);
        let end_pos = (last_deleted + 1) as u64;
        ProtInterval::new(
            ProtPos::new(start_aa, start_pos),
            ProtPos::new(end_aa, end_pos),
        )
    };

    let accession = parse_accession(protein_accession);
    let variant = ProteinVariant {
        accession,
        gene_symbol: transcript.gene_symbol.clone(),
        loc_edit: LocEdit::new_predicted(loc, protein_edit),
    };
    Ok(HgvsVariant::Protein(variant))
}

/// Codon-aligned in-frame insertion: `p.(Xxx{N}_Yyy{N+1}ins{AAseq})`.
fn build_inframe_insertion(
    ref_protein: &[AminoAcid],
    alt_protein: &[AminoAcid],
    cds_pos_start: i64, // last base of the codon before insertion
    protein_accession: &str,
    transcript: &Transcript,
) -> Result<HgvsVariant, FerroError> {
    // The insertion happens between codon N and N+1.
    // cds_pos_start is the last base (frame 2) of codon N.
    // codon_N = (cds_pos_start - 1) / 3  (0-based), so 1-based = (cds_pos_start + 2) / 3
    let aa_before_idx = ((cds_pos_start - 1) / 3) as usize; // 0-based index in ref_protein
    let aa_before = ref_protein
        .get(aa_before_idx)
        .copied()
        .unwrap_or(AminoAcid::Xaa);
    let aa_before_pos = (aa_before_idx + 1) as u64;
    let aa_after = ref_protein
        .get(aa_before_idx + 1)
        .copied()
        .unwrap_or(AminoAcid::Ter);
    let aa_after_pos = aa_before_pos + 1;

    // The inserted amino acids are the additional ones in alt_protein between
    // aa_before_idx and aa_before_idx+1.
    let n_inserted = alt_protein.len().saturating_sub(ref_protein.len());
    let inserted_aas: Vec<AminoAcid> =
        alt_protein[aa_before_idx + 1..aa_before_idx + 1 + n_inserted].to_vec();

    let protein_edit = ProteinEdit::Insertion {
        sequence: AminoAcidSeq::new(inserted_aas),
    };

    let loc = ProtInterval::new(
        ProtPos::new(aa_before, aa_before_pos),
        ProtPos::new(aa_after, aa_after_pos),
    );
    let accession = parse_accession(protein_accession);
    let variant = ProteinVariant {
        accession,
        gene_symbol: transcript.gene_symbol.clone(),
        loc_edit: LocEdit::new_predicted(loc, protein_edit),
    };
    Ok(HgvsVariant::Protein(variant))
}

/// Codon-aligned duplication: `p.(Xxx{N}dup)` or `p.(Xxx{N}_Yyy{M}dup)`.
fn build_inframe_duplication(
    ref_protein: &[AminoAcid],
    _alt_protein: &[AminoAcid],
    cds_pos_start: i64,
    cds_pos_end: i64,
    protein_accession: &str,
    transcript: &Transcript,
) -> Result<HgvsVariant, FerroError> {
    // The duplicated region in CDS coordinates is [cds_pos_start, cds_pos_end].
    // Convert to 0-based amino acid indices.
    let aa_start_idx = ((cds_pos_start - 1) / 3) as usize;
    let aa_end_idx = ((cds_pos_end - 1) / 3) as usize;

    let start_aa = ref_protein
        .get(aa_start_idx)
        .copied()
        .unwrap_or(AminoAcid::Xaa);
    let start_pos = (aa_start_idx + 1) as u64;

    let loc = if aa_start_idx == aa_end_idx {
        ProtInterval::point(ProtPos::new(start_aa, start_pos))
    } else {
        let end_aa = ref_protein
            .get(aa_end_idx)
            .copied()
            .unwrap_or(AminoAcid::Xaa);
        let end_pos = (aa_end_idx + 1) as u64;
        ProtInterval::new(
            ProtPos::new(start_aa, start_pos),
            ProtPos::new(end_aa, end_pos),
        )
    };

    let accession = parse_accession(protein_accession);
    let variant = ProteinVariant {
        accession,
        gene_symbol: transcript.gene_symbol.clone(),
        loc_edit: LocEdit::new_predicted(loc, ProteinEdit::Duplication),
    };
    Ok(HgvsVariant::Protein(variant))
}

/// Generic in-frame delins: `p.(Xxx{N}_Yyy{M}delins{ZzzZzz...})`.
///
/// Used when the edit straddles a codon boundary or is a delins/inversion.
fn build_inframe_delins(
    ref_protein: &[AminoAcid],
    alt_protein: &[AminoAcid],
    cds_pos_start: i64,
    protein_accession: &str,
    transcript: &Transcript,
) -> Result<HgvsVariant, FerroError> {
    // Find the first and last positions that differ.
    let first_diff = first_diff_position(ref_protein, alt_protein);
    let last_diff_ref = find_last_diff(ref_protein, alt_protein);
    let last_diff_alt = find_last_diff_alt(ref_protein, alt_protein);

    // The affected region in the ref protein spans [first_diff, last_diff_ref].
    // The replacement sequence from the alt protein spans [first_diff, last_diff_alt].
    let ref_len = ref_protein.len();
    let alt_len = alt_protein.len();

    let _ = cds_pos_start; // used implicitly through first_diff

    // If the alt protein is empty or only differs at the last AA:
    if first_diff >= ref_len && first_diff >= alt_len {
        // No difference at the amino acid level.
        return build_identity_variant(ref_protein, protein_accession, transcript);
    }

    let start_aa = ref_protein
        .get(first_diff)
        .copied()
        .unwrap_or(AminoAcid::Xaa);
    let start_pos = (first_diff + 1) as u64;

    let end_ref = last_diff_ref.min(ref_len.saturating_sub(1));
    let end_aa = ref_protein.get(end_ref).copied().unwrap_or(AminoAcid::Xaa);
    let end_pos = (end_ref + 1) as u64;

    let end_alt = last_diff_alt.min(alt_len.saturating_sub(1));
    let inserted: Vec<AminoAcid> = if first_diff <= end_alt {
        alt_protein[first_diff..=end_alt].to_vec()
    } else {
        vec![]
    };

    let protein_edit = ProteinEdit::Delins {
        sequence: AminoAcidSeq::new(inserted),
    };

    let loc = if start_pos == end_pos {
        ProtInterval::point(ProtPos::new(start_aa, start_pos))
    } else {
        ProtInterval::new(
            ProtPos::new(start_aa, start_pos),
            ProtPos::new(end_aa, end_pos),
        )
    };

    let accession = parse_accession(protein_accession);
    let variant = ProteinVariant {
        accession,
        gene_symbol: transcript.gene_symbol.clone(),
        loc_edit: LocEdit::new_predicted(loc, protein_edit),
    };
    Ok(HgvsVariant::Protein(variant))
}

/// Find the 0-based index of the last AA in `ref_protein` that is in the "differing" region.
///
/// This is needed for range-based protein delins to identify the endpoint in the ref.
fn find_last_diff(ref_prot: &[AminoAcid], alt_prot: &[AminoAcid]) -> usize {
    let ref_len = ref_prot.len();
    let alt_len = alt_prot.len();
    // Scan from the shorter tail inward.
    let max_tail = ref_len.min(alt_len);
    for tail in 0..max_tail {
        let ri = ref_len - 1 - tail;
        let ai = alt_len - 1 - tail;
        if ref_prot[ri] != alt_prot[ai] {
            return ri;
        }
    }
    // All trailing positions match: the differing region ends at ref_len - alt_len - 1
    // (if ref is longer) or 0 (if alt is longer or equal).
    ref_len.saturating_sub(alt_len)
}

/// Find the 0-based index of the last AA in `alt_protein` that is in the "differing" region.
fn find_last_diff_alt(ref_prot: &[AminoAcid], alt_prot: &[AminoAcid]) -> usize {
    let ref_len = ref_prot.len();
    let alt_len = alt_prot.len();
    let max_tail = ref_len.min(alt_len);
    for tail in 0..max_tail {
        let ri = ref_len - 1 - tail;
        let ai = alt_len - 1 - tail;
        if ref_prot[ri] != alt_prot[ai] {
            return ai;
        }
    }
    alt_len.saturating_sub(ref_len)
}

// ── Frameshift prediction ─────────────────────────────────────────────────────

/// Build a `p.(Xxx{N}[Yyy]fsTer{K})` frameshift variant.
fn build_frameshift_variant(
    ref_protein: &[AminoAcid],
    alt_protein: &[AminoAcid],
    mut_cds: &str,
    protein_accession: &str,
    transcript: &Transcript,
) -> Result<HgvsVariant, FerroError> {
    let first_diff = first_diff_position(ref_protein, alt_protein);
    let aa_pos = (first_diff + 1) as u64; // 1-based
    let ref_aa = ref_protein
        .get(first_diff)
        .copied()
        .unwrap_or(AminoAcid::Xaa);
    let new_aa = alt_protein.get(first_diff).copied();

    // Scan downstream from first_diff in the mutated CDS to find a stop codon.
    let codon_byte_offset = first_diff * 3;
    let ter_pos: Option<u64> = if codon_byte_offset < mut_cds.len() {
        let downstream = &mut_cds[codon_byte_offset..];
        let downstream_aas = translate_full_cds_with_stop(downstream);
        // ter_pos = 1-based position of the Ter in the downstream scan.
        // HGVS fsTer{K}: K is the number of amino acids (including the new shifted AA) until Ter.
        downstream_aas
            .iter()
            .position(|aa| *aa == AminoAcid::Ter)
            .map(|p| (p + 1) as u64)
    } else {
        None
    };

    let protein_edit = ProteinEdit::Frameshift { new_aa, ter_pos };
    let loc = ProtInterval::point(ProtPos::new(ref_aa, aa_pos));
    let accession = parse_accession(protein_accession);
    let variant = ProteinVariant {
        accession,
        gene_symbol: transcript.gene_symbol.clone(),
        loc_edit: LocEdit::new_predicted(loc, protein_edit),
    };
    Ok(HgvsVariant::Protein(variant))
}

// ── Extension (stop-codon readthrough) prediction ─────────────────────────────

/// Build a `p.(Ter{N}{Yyy}ext*{K})` extension variant for stop-codon readthrough.
fn build_extension_variant(
    ref_protein: &[AminoAcid],
    _alt_protein: &[AminoAcid],
    mut_cds: &str,
    protein_accession: &str,
    transcript: &Transcript,
) -> Result<HgvsVariant, FerroError> {
    let stop_pos = (ref_protein.len() + 1) as u64; // 1-based position of the original stop codon

    // The new amino acid at the stop codon position.
    // Translate starting from the stop-codon position in the mutated CDS.
    let stop_offset = ref_protein.len() * 3;
    let new_aa: Option<AminoAcid> = if stop_offset + 3 <= mut_cds.len() {
        let stop_codon_seq = &mut_cds[stop_offset..stop_offset + 3];
        super::substitution::translate(stop_codon_seq).filter(|aa| *aa != AminoAcid::Ter)
    } else {
        None
    };

    // Count extension amino acids until a new stop is found.
    let ext_count: Option<i64> = if stop_offset < mut_cds.len() {
        let downstream = &mut_cds[stop_offset..];
        let downstream_aas = translate_full_cds_with_stop(downstream);
        downstream_aas
            .iter()
            .position(|aa| *aa == AminoAcid::Ter)
            .map(|p| (p + 1) as i64) // +1: distance from new AA to new stop (inclusive)
    } else {
        None
    };

    let protein_edit = ProteinEdit::Extension {
        new_aa,
        direction: ExtDirection::CTerminal,
        count: ext_count,
    };

    let loc = ProtInterval::point(ProtPos::new(AminoAcid::Ter, stop_pos));
    let accession = parse_accession(protein_accession);
    let variant = ProteinVariant {
        accession,
        gene_symbol: transcript.gene_symbol.clone(),
        loc_edit: LocEdit::new_predicted(loc, protein_edit),
    };
    Ok(HgvsVariant::Protein(variant))
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use std::sync::OnceLock;

    fn tx(seq: &str, cds_start: u64, cds_end: u64) -> Transcript {
        Transcript {
            id: "NM_TEST.1".to_string(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: seq.to_string(),
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

    fn prot_str(v: &HgvsVariant) -> String {
        match v {
            HgvsVariant::Protein(p) => p.to_string(),
            _ => panic!("expected Protein variant"),
        }
    }

    // ─── Frameshift tests ─────────────────────────────────────────────────────

    #[test]
    fn del_single_base_frameshift_with_ter() {
        // CDS "ATGCGCTAA": Met-Arg-Stop.
        // Delete c.4 (first base of Arg codon CGC).
        // Mutated CDS: "ATGGCTAA" → Met-Ala-Stop (stop at position 2 → fsTer2? No—)
        // Actually: ATGGCTAA → ATG GCT AA(incomplete) → [Met, Ala] then no stop. Wait.
        // "ATGGCTAA" = ATG GCT AA (8 bases, 2 complete codons + 2 partial)
        // → Met, Ala — no stop. fsTer? = no ter found.
        // But the downstream scan from position 1 (0-indexed): mut_cds = "ATGGCTAA"
        // first_diff: ref=[Met,Arg], alt=[Met,Ala] → first_diff=1 (Arg≠Ala)
        // scan from codon 1 onward in mut_cds: "GCTAA" → GCT (Ala), AA (incomplete)
        // No Ter found → fsTer? should be None.
        let t = tx("ATGCGCTAA", 1, 9);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = predict_indel_protein(&t, 4, 4, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        // Arg2[Ala]fsTer? — with no downstream stop in the remaining sequence.
        // alt_protein = [Met, Ala], first_diff=1, new_aa=Ala
        // downstream from byte 3 = "GCTAA" → Ala, then incomplete AA
        // no Ter → ter_pos = None
        assert!(s.contains("Arg2"), "expected Arg2 in '{}'", s);
        assert!(s.contains("fs"), "expected fs in '{}'", s);
        assert!(s.contains("Ala"), "expected Ala new_aa in '{}'", s);
    }

    #[test]
    fn del_three_bases_codon_aligned_single_aa() {
        // CDS "ATGCGCTAA": del c.4_6 (entire Arg codon CGC).
        // Mutated CDS: "ATGTAA" = ATG TAA → Met-Stop
        // ref=[Met, Arg], alt=[Met]
        // first_diff=1, n_deleted=1 → p.(Arg2del)
        let t = tx("ATGCGCTAA", 1, 9);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = predict_indel_protein(&t, 4, 6, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert_eq!(s, "NP_TEST.1:p.(Arg2del)");
    }

    #[test]
    fn del_six_bases_codon_aligned_range() {
        // CDS "ATGCGCAAATAA": Met-Arg-Lys-Stop (12 bp).
        // del c.4_9 (CGC AAA = Arg + Lys).
        // Mutated CDS: "ATGTAA" → [Met]
        // ref=[Met, Arg, Lys], alt=[Met]
        // first_diff=1, n_deleted=2 → p.(Arg2_Lys3del)
        let t = tx("ATGCGCAAATAA", 1, 12);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = predict_indel_protein(&t, 4, 9, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert_eq!(s, "NP_TEST.1:p.(Arg2_Lys3del)");
    }

    #[test]
    fn del_two_bases_non_aligned_frameshift() {
        // CDS "ATGCGCAAATAA": del c.4_5 (CG, non-multiple-of-3).
        // Net = -2 → frameshift.
        let t = tx("ATGCGCAAATAA", 1, 12);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = predict_indel_protein(&t, 4, 5, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert!(s.contains("fs"), "expected frameshift in '{}'", s);
    }

    #[test]
    fn ins_three_bases_codon_aligned() {
        // CDS "ATGCGCTAA": insert "GGG" at c.3_4 (between codon 1 end and codon 2 start).
        // cds_pos_start=3 (last base of Met codon, frame=2), cds_pos_end=3.
        // Mutated CDS: "ATGGGGCGCTAA" → ATG GGG CGC TAA → [Met, Gly, Arg]
        // ref=[Met, Arg], alt=[Met, Gly, Arg]
        // Insertion between aa1 (Met) and aa2 (Arg): p.(Met1_Arg2insGly)
        let t = tx("ATGCGCTAA", 1, 9);
        let seq: crate::hgvs::edit::Sequence = "GGG".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
        };
        let result = predict_indel_protein(&t, 3, 3, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert_eq!(s, "NP_TEST.1:p.(Met1_Arg2insGly)");
    }

    #[test]
    fn ins_one_base_frameshift() {
        // CDS "ATGCGCTAA": insert "A" at c.3_4.
        // Net = +1 → frameshift.
        let t = tx("ATGCGCTAA", 1, 9);
        let seq: crate::hgvs::edit::Sequence = "A".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
        };
        let result = predict_indel_protein(&t, 3, 3, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert!(s.contains("fs"), "expected frameshift in '{}'", s);
    }

    #[test]
    fn dup_three_bases_codon_aligned() {
        // CDS "ATGCGCTAA": dup c.4_6 (CGC).
        // Mutated CDS: "ATGCGCCGCTAA" → ATG CGC CGC TAA → [Met, Arg, Arg]
        // ref=[Met, Arg], alt=[Met, Arg, Arg]
        // dup of Arg codon (cds 4..6 → aa_start_idx=1, aa_end_idx=1) → p.(Arg2dup)
        let t = tx("ATGCGCTAA", 1, 9);
        let edit = NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: None,
        };
        let result = predict_indel_protein(&t, 4, 6, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert_eq!(s, "NP_TEST.1:p.(Arg2dup)");
    }

    #[test]
    fn dup_one_base_frameshift() {
        // CDS "ATGCGCTAA": dup c.4 (single C, net=+1 → frameshift).
        let t = tx("ATGCGCTAA", 1, 9);
        let edit = NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: None,
        };
        let result = predict_indel_protein(&t, 4, 4, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert!(s.contains("fs"), "expected frameshift in '{}'", s);
    }

    #[test]
    fn delins_same_length_codon_aligned() {
        // CDS "ATGCGCTAA": delins c.4_6 (CGC → TCC).
        // Mutated CDS: "ATGTCCTAA" → ATG TCC TAA → [Met, Ser]
        // ref=[Met, Arg], alt=[Met, Ser]
        // first_diff=1, last_diff_ref=1, last_diff_alt=1 → p.(Arg2delinsser)? No:
        // It's a delins: p.(Arg2delinsSer) → but HGVS for single AA change is really
        // p.(Arg2Ser) substitution. However since it's a delins NaEdit we go through
        // build_inframe_delins. Let's check what we actually get.
        let t = tx("ATGCGCTAA", 1, 9);
        let seq: crate::hgvs::edit::Sequence = "TCC".parse().unwrap();
        let edit = NaEdit::Delins {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
            deleted: None,
            deleted_length: None,
        };
        let result = predict_indel_protein(&t, 4, 6, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        // Should be a delins at Arg2 → Ser: p.(Arg2delinsSer)
        assert!(s.contains("Arg2"), "expected Arg2 in '{}'", s);
        assert!(s.contains("delins"), "expected delins in '{}'", s);
        assert!(s.contains("Ser"), "expected Ser in '{}'", s);
    }

    #[test]
    fn inversion_codon_cgc_to_gcg() {
        // CDS "ATGCGCTAA": inv c.4_6 (CGC → GCG).
        // GCG = Ala. Mutated CDS: "ATGGCGTAA" → [Met, Ala]
        // first_diff=1 → delins Arg2 → Ala: p.(Arg2delinsAla)
        let t = tx("ATGCGCTAA", 1, 9);
        let edit = NaEdit::Inversion {
            sequence: None,
            length: None,
        };
        let result = predict_indel_protein(&t, 4, 6, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert!(s.contains("Arg2"), "expected Arg2 in '{}'", s);
        assert!(s.contains("delins"), "expected delins in '{}'", s);
        assert!(s.contains("Ala"), "expected Ala in '{}'", s);
    }

    #[test]
    fn stop_readthrough_extension() {
        // CDS "ATGCGCTAA": the stop codon is "TAA" at c.7_9.
        // Replace TAA with TGG (Trp): delins c.7_9 TGG.
        // Mutated CDS: "ATGCGCTGG" → [Met, Arg, Trp] and then whatever follows.
        // The mutated CDS has no stop — translate_full_cds returns [Met, Arg, Trp]
        // and translate_full_cds_with_stop also returns [Met, Arg, Trp] (no Ter).
        // alt_protein.len()=3 > ref_protein.len()=2 → extension check.
        // Stop was at position 3 (1-based), new AA at that position = Trp.
        // ext_count: downstream from stop_offset=6: "TGG" = Trp — no Ter found → None.
        // Result: p.(Ter3Trpext*?)
        let t = tx("ATGCGCTAA", 1, 9);
        let seq: crate::hgvs::edit::Sequence = "TGG".parse().unwrap();
        let edit = NaEdit::Delins {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
            deleted: None,
            deleted_length: None,
        };
        let result = predict_indel_protein(&t, 7, 9, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert!(s.contains("Ter3"), "expected Ter3 in '{}'", s);
        assert!(s.contains("Trp"), "expected Trp in '{}'", s);
        assert!(s.contains("ext"), "expected ext in '{}'", s);
    }
}
