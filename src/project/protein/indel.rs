//! Protein consequence prediction for indels (deletion, insertion, duplication,
//! deletion-insertion, inversion) in the CDS.

use crate::error::FerroError;
use crate::hgvs::edit::{AminoAcidSeq, FrameshiftTer, NaEdit, ProteinEdit};
use crate::hgvs::interval::ProtInterval;
use crate::hgvs::location::{AminoAcid, ProtPos};
use crate::hgvs::variant::{HgvsVariant, LocEdit, ProteinVariant};
use crate::project::accession::parse_accession;
use crate::reference::transcript::Transcript;

use super::helpers::{
    affects_initiation_codon, build_cterminal_extension, build_initiator_unknown,
    build_mutated_cds_with_ref, cds_has_recognized_start, first_diff_position, force_initiator_met,
    mut_cds_with_3utr, net_length_change, translate_full_cds_with_stop, translate_mutated_cds,
    translate_mutated_cds_inframe, RefProteinBundle,
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
    ref_bundle: &RefProteinBundle,
    cds_pos_start: i64,
    cds_pos_end: i64,
    edit: &NaEdit,
    protein_accession: &str,
) -> Result<HgvsVariant, FerroError> {
    // A change touching the translation initiation codon (CDS 1–3) has an
    // unpredictable protein consequence — report p.(Met1?) up front rather
    // than a concrete del/ins/delins, which the spec disallows (#498).
    if affects_initiation_codon(edit, cds_pos_start, cds_pos_end) {
        return Ok(build_initiator_unknown(protein_accession, transcript));
    }

    // An HGVS insertion `c.A_(A+1)insX` is anchored at A, with the inserted bases
    // in the gap immediately 3' of A. Callers pass the literal `A_(A+1)` interval,
    // but the indel machinery treats an insertion's span as the single anchor
    // position A (the gap after it). Collapse the end to the start so the inserted
    // bases land between A and A+1; without this, `build_mutated_cds_with_ref`
    // splices them one position too far 3' (after A+1), spuriously preserving the
    // flanking codon's residue and yielding a clean `ins` where a mid-codon
    // insertion must be a `delins` (#511). Other edits keep their `[A,B]` span.
    let cds_pos_end = if matches!(edit, NaEdit::Insertion { .. }) {
        cds_pos_start
    } else {
        cds_pos_end
    };
    // 1. The reference CDS + translation are pre-computed by the caller (and
    //    cached on `VariantProjector` so fan-out across many indels on the
    //    same transcript doesn't re-translate (or re-read+re-uppercase) the
    //    full CDS each time). The mutated CDS is built by splicing the cached
    //    `ref_cds` directly — no `read_full_cds` per call.
    let RefProteinBundle {
        ref_cds,
        ref_protein,
        ref_protein_with_stop,
    } = ref_bundle;
    let mut_cds =
        build_mutated_cds_with_ref(ref_cds, transcript, cds_pos_start, cds_pos_end, edit)?;

    // 2. Compute net length change up-front so we can pick the right alt
    //    translation strategy. (We also use `net` again below for the
    //    frameshift / in-frame branch — moving it here just lifts the
    //    existing call site, no semantic change.)
    let del_len = (cds_pos_end - cds_pos_start + 1) as usize;
    let net =
        net_length_change(edit, del_len).ok_or_else(|| FerroError::UnsupportedProjection {
            reason: "cannot determine net length change for protein prediction".to_string(),
        })?;
    let is_frameshift = net % 3 != 0;

    // 3. Translate the mutated CDS. We can always skip codons whose three
    //    bytes are byte-identical to the reference (the unchanged prefix
    //    before the edit). For in-frame edits we can also skip the unchanged
    //    suffix after the edit: those bytes are a codon-aligned shift of
    //    `ref_cds`, so their AAs are already in `ref_protein`. For
    //    frameshift edits the suffix can't be lifted (frame shift breaks the
    //    codon alignment between mut and ref) — only the prefix lift applies.
    let alt_protein = if is_frameshift {
        let unchanged_prefix_codons = ((cds_pos_start - 1).max(0) / 3) as usize;
        translate_mutated_cds(ref_protein, &mut_cds, unchanged_prefix_codons)
    } else {
        translate_mutated_cds_inframe(ref_protein, &mut_cds, cds_pos_start, cds_pos_end, net)
    };

    // 3. Check for stop-codon readthrough (stop-loss extension). The reference
    //    ends in a terminator and the edit disrupts that terminator codon while
    //    leaving every sense residue intact, so the former stop now codes an
    //    amino acid and translation reads through the 3'UTR to the next stop:
    //    `p.(Ter<pos><aa>extTer<k>)` (extension.md:30; #498, #615).
    //
    //    Detection must read the mutated CDS *through the 3'UTR*: a small del
    //    that disrupts the stop shortens the CDS-only `alt_protein` so that
    //    `alt_protein.len()` is NOT longer than `ref_protein.len()` — the old
    //    length-based check missed it and the variant fell through to the
    //    frameshift path, which then read the (nonexistent) reference residue
    //    at the stop slot and emitted the invalid `Xaa` token (#615). Instead
    //    we translate the full read-through sequence once and route to the
    //    extension builder iff the protein is unchanged up to and including the
    //    last sense residue and the former stop now codes a non-stop AA. A
    //    frameshift *before* the stop diverges earlier and is excluded.
    if ref_protein_with_stop.last() == Some(&AminoAcid::Ter)
        && is_stop_loss_readthrough(
            ref_protein,
            &mut_cds,
            transcript,
            cds_pos_end,
            is_frameshift,
        )?
    {
        return build_extension_variant(
            ref_protein,
            &alt_protein,
            &mut_cds,
            protein_accession,
            transcript,
        );
    }

    // 4. Build the protein variant (net + is_frameshift were computed above).
    if is_frameshift {
        build_frameshift_variant(
            ref_protein,
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
            ref_protein,
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

    // Did the edit introduce a premature termination codon? An in-frame edit
    // that does not truncate yields `ref_protein.len() + net/3` residues; a
    // shorter `alt_protein` means translation halted at a new stop before the
    // reference terminus (the truncating `Ter` is dropped by translation). Such
    // a change must be rendered with the terminating `Ter`, capped at the change
    // site — never as a C-terminal-spanning delins (#911, insertion.md:37-41,
    // delins.md:27-28). The codon-aligned-insertion arm handles its own stop in
    // `build_inframe_insertion`; every other path routes through the
    // stop-aware delins builder below.
    let expected_alt_len = (ref_protein.len() as i64 + net / 3).max(0) as usize;
    let premature_stop = alt_protein.len() < expected_alt_len;

    // Dispatch to the stop-aware or plain delins builder as appropriate.
    let delins = |cds_start: i64| -> Result<HgvsVariant, FerroError> {
        if premature_stop {
            build_inframe_stop_delins(
                ref_protein,
                alt_protein,
                cds_pos_end,
                protein_accession,
                transcript,
            )
        } else {
            build_inframe_delins(
                ref_protein,
                alt_protein,
                cds_start,
                protein_accession,
                transcript,
            )
        }
    };

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
                delins(cds_pos_start)
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
                    premature_stop,
                    protein_accession,
                    transcript,
                )
            } else {
                // Mid-codon insertion: treated as a delins at the AA level.
                delins(cds_pos_start)
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
                delins(cds_pos_start)
            }
        }
        NaEdit::Delins { .. } | NaEdit::Inversion { .. } => {
            // Always use the generic delins pathway for delins and inversion.
            delins(cds_pos_start)
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

    if n_deleted == 0 {
        // Guard before subtracting: alt_protein.len() >= ref_protein.len() means the
        // "deletion" didn't shorten the protein (e.g. a stop-disrupting deletion whose
        // extension detection missed an edge case). Avoid the underflow in the
        // last_deleted computation by short-circuiting to an identity variant.
        return build_identity_variant(ref_protein, protein_accession, transcript);
    }

    let last_deleted = first_diff + n_deleted - 1;

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
    cds_pos_start: i64,   // last base of the codon before insertion
    premature_stop: bool, // did translation halt at a new stop before the reference terminus?
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

    // The inserted amino acids are the residues in `alt_protein` after the
    // unchanged prefix `ref_protein[..=aa_before_idx]`. Decide the branch from
    // the caller's `premature_stop` flag, NOT `alt.len() > ref.len()`: a
    // codon-aligned insertion near the C-terminus can introduce a premature stop
    // while `alt_protein` is still longer than `ref_protein` (when the inserted
    // sense residues before the new stop outnumber the reference residues after
    // the insertion point), and a length comparison would misclassify it as a
    // plain insertion and drop the terminating `Ter` (#911).
    let inserted_aas: Vec<AminoAcid> = if !premature_stop {
        // No premature stop: the inserted residues are the extra ones in alt.
        let n_inserted = alt_protein.len() - ref_protein.len();
        let insert_end = (aa_before_idx + 1 + n_inserted).min(alt_protein.len());
        alt_protein[aa_before_idx + 1..insert_end].to_vec()
    } else {
        // The inserted nucleotides translate to a premature stop, so
        // translation truncated `alt_protein` at the new (dropped) Ter. Per
        // `insertion.md` L37-41 (`p.(Met3_His4insGlyTer)`,
        // `p.(Pro46_Asn47insSerSerTer)`), this is described as an **insertion**
        // of the translated residues up to and **including** the terminating
        // `Ter` — NOT a deletion-insertion replacing the entire C-terminal
        // sequence (explicitly forbidden), and residues after the stop are not
        // listed. The inserted residues are the alt tail after the prefix; the
        // stop that truncated translation is re-attached as `Ter` (#911).
        let tail_start = (aa_before_idx + 1).min(alt_protein.len());
        let mut ins = alt_protein[tail_start..].to_vec();
        ins.push(AminoAcid::Ter);
        ins
    };

    // An **immediate** stop insertion — the inserted sequence is a bare `Ter`
    // with no sense residue before it — is a nonsense variant, NOT an insertion
    // of `Ter` (insertion.md:25; the `ins…Ter` form of insertion.md:26 is for a
    // stop *encoded within* a sense-carrying insert). The stop lands at the
    // residue immediately 3' of the insertion point (#911).
    if inserted_aas.as_slice() == [AminoAcid::Ter] && aa_before_idx + 1 < ref_protein.len() {
        let variant = ProteinVariant {
            accession: parse_accession(protein_accession),
            gene_symbol: transcript.gene_symbol.clone(),
            loc_edit: LocEdit::new_predicted(
                ProtInterval::point(ProtPos::new(aa_after, aa_after_pos)),
                ProteinEdit::Substitution {
                    reference: aa_after,
                    alternative: AminoAcid::Ter,
                },
            ),
        };
        return Ok(HgvsVariant::Protein(variant));
    }

    // A directly-C-terminal tandem copy is a duplication, not an insertion
    // (duplication.md). The insertion sits before `ref_protein[aa_before_idx + 1]`.
    // A stop-terminated insertion is never a plain tandem duplication.
    if inserted_aas.last() != Some(&AminoAcid::Ter) {
        if let Some((dup_start, dup_end)) =
            protein_tandem_dup_range(ref_protein, aa_before_idx + 1, &inserted_aas)
        {
            return Ok(build_dup_variant(
                ref_protein,
                dup_start,
                dup_end,
                protein_accession,
                transcript,
            ));
        }
    }

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

/// Detect whether a pure protein insertion is a tandem duplication.
///
/// Per HGVS `duplication.md`, an inserted residue block that is a copy of the
/// residues **directly C-terminal** of (i.e. immediately 5' of, in tandem with)
/// the insertion point is a *duplication*, not an insertion. The non-adjacent
/// case — where the inserted sequence is present elsewhere in the protein but not
/// directly flanking — stays an insertion (the `p.His7_Gln8insGly4_Ser6` example
/// in duplication.md's discussion). The protein 3' rule ("the most C-terminal
/// residue is arbitrarily assigned to have been duplicated") is applied by
/// sliding the insertion point C-terminal across any equal-residue run before the
/// tandem test.
///
/// `point` is the 0-based index in `ref_protein` *before which* `inserted` is
/// added (so `inserted` sits between `ref_protein[point-1]` and
/// `ref_protein[point]`). Returns the 0-based inclusive `[start, end]` range of
/// the duplicated residues in `ref_protein`, or `None` when the insertion is not
/// a directly-C-terminal tandem copy (a genuine insertion). Mirrors the
/// DNA-level `normalize::rules::insertion_is_duplication` semantics at the
/// amino-acid level.
fn protein_tandem_dup_range(
    ref_protein: &[AminoAcid],
    point: usize,
    inserted: &[AminoAcid],
) -> Option<(usize, usize)> {
    let unit = inserted.len();
    if unit == 0 || point > ref_protein.len() {
        return None;
    }
    // Apply the protein 3' rule: slide the insertion point C-terminal across any
    // run where the residue at the point matches the head of the (rotating)
    // inserted unit. This places a homopolymer/tandem insertion at its most
    // C-terminal equivalent position before the tandem test.
    let mut point = point;
    let mut rotated: Vec<AminoAcid> = inserted.to_vec();
    while point < ref_protein.len() && ref_protein[point] == rotated[0] {
        rotated.rotate_left(1);
        point += 1;
    }
    // Tandem test: the `unit` residues immediately 5' of the (shifted) insertion
    // point must equal the (rotated) inserted unit — a copy lying directly
    // C-terminal of the original.
    let start = point.checked_sub(unit)?;
    if ref_protein[start..point] == rotated[..] {
        Some((start, point - 1))
    } else {
        None
    }
}

/// Build a `p.(Xxx{N}dup)` / `p.(Xxx{N}_Yyy{M}dup)` variant for the 0-based
/// inclusive residue range `[start_idx, end_idx]` in `ref_protein`. Shares the
/// location shape with [`build_inframe_duplication`] (point when single residue,
/// interval otherwise).
fn build_dup_variant(
    ref_protein: &[AminoAcid],
    start_idx: usize,
    end_idx: usize,
    protein_accession: &str,
    transcript: &Transcript,
) -> HgvsVariant {
    let start_aa = ref_protein
        .get(start_idx)
        .copied()
        .unwrap_or(AminoAcid::Xaa);
    let start_pos = (start_idx + 1) as u64;
    let loc = if start_idx == end_idx {
        ProtInterval::point(ProtPos::new(start_aa, start_pos))
    } else {
        let end_aa = ref_protein.get(end_idx).copied().unwrap_or(AminoAcid::Xaa);
        ProtInterval::new(
            ProtPos::new(start_aa, start_pos),
            ProtPos::new(end_aa, (end_idx + 1) as u64),
        )
    };
    HgvsVariant::Protein(ProteinVariant {
        accession: parse_accession(protein_accession),
        gene_symbol: transcript.gene_symbol.clone(),
        loc_edit: LocEdit::new_predicted(loc, ProteinEdit::Duplication),
    })
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

    // Pure-deletion guard (#847): when the last differing ref index falls
    // *before* `first_diff` and the alt protein is *shorter* than the ref, no
    // ref residue is replaced — the net effect is a pure deletion, not a
    // delins. This arises for a codon-straddling in-frame deletion inside a run
    // of identical residues (e.g. one Leu removed from a Leu run): the forward
    // `first_diff` and backward `find_last_diff` scans cross over because every
    // residue in the run is identical. Emit a `del` over the 3'-most deleted
    // residues, mirroring `build_inframe_deletion`'s arithmetic. Without this,
    // the pure-insertion guard below fires with `n_inserted == 0` and renders
    // the invalid `p.(Xaa_Yaa ins)` (empty inserted sequence).
    if last_diff_ref < first_diff && alt_len < ref_len && first_diff < ref_len {
        let n_deleted = ref_len - alt_len;
        let first_del = first_diff;
        let last_del = (first_diff + n_deleted - 1).min(ref_len - 1);
        let start_aa = ref_protein
            .get(first_del)
            .copied()
            .unwrap_or(AminoAcid::Xaa);
        let start_pos = (first_del + 1) as u64;
        let loc = if first_del == last_del {
            ProtInterval::point(ProtPos::new(start_aa, start_pos))
        } else {
            let end_aa = ref_protein.get(last_del).copied().unwrap_or(AminoAcid::Xaa);
            let end_pos = (last_del + 1) as u64;
            ProtInterval::new(
                ProtPos::new(start_aa, start_pos),
                ProtPos::new(end_aa, end_pos),
            )
        };
        let variant = ProteinVariant {
            accession: parse_accession(protein_accession),
            gene_symbol: transcript.gene_symbol.clone(),
            loc_edit: LocEdit::new_predicted(
                loc,
                ProteinEdit::Deletion {
                    sequence: None,
                    count: None,
                },
            ),
        };
        return Ok(HgvsVariant::Protein(variant));
    }

    // Pure-insertion guard (#498): when the last differing ref index falls
    // *before* `first_diff` and the alt protein is *longer* than the ref, the
    // ref differing region is empty — no ref residues are replaced, so this is
    // an insertion, not a delins. Rendering it as a delins would cross the
    // positions (`start > end`) into malformed HGVS. Emit an insertion between
    // the two flanking ref residues with ascending positions. (Requires an
    // interior insertion point with a residue on each side; the N-terminal edge
    // `first_diff == 0` falls through to the clamped delins path below, which
    // stays ascending.)
    if last_diff_ref < first_diff && alt_len > ref_len && first_diff >= 1 && first_diff < ref_len {
        let n_inserted = alt_len.saturating_sub(ref_len);
        let left_idx = first_diff - 1;
        let left_aa = ref_protein.get(left_idx).copied().unwrap_or(AminoAcid::Xaa);
        let right_aa = ref_protein
            .get(first_diff)
            .copied()
            .unwrap_or(AminoAcid::Ter);
        let inserted: Vec<AminoAcid> = alt_protein
            .get(first_diff..first_diff + n_inserted)
            .map(|s| s.to_vec())
            .unwrap_or_default();
        // A directly-C-terminal tandem copy is a duplication, not an insertion
        // (duplication.md). `first_diff` is the index in `ref_protein` before
        // which the residues are inserted.
        if let Some((dup_start, dup_end)) =
            protein_tandem_dup_range(ref_protein, first_diff, &inserted)
        {
            return Ok(build_dup_variant(
                ref_protein,
                dup_start,
                dup_end,
                protein_accession,
                transcript,
            ));
        }
        let loc = ProtInterval::new(
            ProtPos::new(left_aa, (left_idx + 1) as u64),
            ProtPos::new(right_aa, (first_diff + 1) as u64),
        );
        let variant = ProteinVariant {
            accession: parse_accession(protein_accession),
            gene_symbol: transcript.gene_symbol.clone(),
            loc_edit: LocEdit::new_predicted(
                loc,
                ProteinEdit::Insertion {
                    sequence: AminoAcidSeq::new(inserted),
                },
            ),
        };
        return Ok(HgvsVariant::Protein(variant));
    }

    let start_aa = ref_protein
        .get(first_diff)
        .copied()
        .unwrap_or(AminoAcid::Xaa);
    let start_pos = (first_diff + 1) as u64;

    // Clamp the end to never fall below the start: `find_last_diff` can return
    // an index below `first_diff` for the degenerate cases not handled by the
    // pure-insertion guard above, which would otherwise yield a descending
    // (malformed) range.
    let end_ref = last_diff_ref.max(first_diff).min(ref_len.saturating_sub(1));
    let end_aa = ref_protein.get(end_ref).copied().unwrap_or(AminoAcid::Xaa);
    let end_pos = (end_ref + 1) as u64;

    let end_alt = last_diff_alt.min(alt_len.saturating_sub(1));
    let inserted: Vec<AminoAcid> = if first_diff <= end_alt {
        alt_protein[first_diff..=end_alt].to_vec()
    } else {
        vec![]
    };

    // One amino acid replaced by one (non-stop) amino acid is, by definition, a
    // substitution, not a deletion-insertion (delins.md: "when **one** amino acid
    // is replaced by **one** other amino acid, the change is a substitution").
    // This applies to every edit routed here, including a whole-codon inversion
    // that changes a single residue. A premature-stop edit never reaches this
    // function — it is routed to `build_inframe_stop_delins` upstream (#911) —
    // so the `inserted[0] != Ter` guard below only ever excludes a spurious Ter.
    if start_pos == end_pos && inserted.len() == 1 && inserted[0] != AminoAcid::Ter {
        let protein_edit = ProteinEdit::Substitution {
            reference: start_aa,
            alternative: inserted[0],
        };
        let loc = ProtInterval::point(ProtPos::new(start_aa, start_pos));
        let variant = ProteinVariant {
            accession: parse_accession(protein_accession),
            gene_symbol: transcript.gene_symbol.clone(),
            loc_edit: LocEdit::new_predicted(loc, protein_edit),
        };
        return Ok(HgvsVariant::Protein(variant));
    }

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

/// Build the protein consequence for an in-frame edit that introduces a
/// **premature stop** codon (#911).
///
/// Translation truncated `alt_protein` at the new stop (the `Ter` was dropped),
/// so `alt_protein` = the unchanged prefix + the new residues up to the stop.
/// Per `insertion.md:27-28` / `delins.md:27-28`, the change is described capped
/// at the edit's site with the terminating `Ter`, never as a deletion-insertion
/// replacing the entire C-terminal sequence (the forbidden form), and residues
/// after the stop are not listed:
/// - an **immediate** stop (the inserted sequence is a bare `Ter`, no sense
///   residue before it) is a nonsense **substitution** (`p.Xaa{N}Ter`),
///   regardless of how many codons the edit deletes (delins.md:27,
///   deletion.md:22-24);
/// - otherwise a `delins` over the codons the edit deletes, whose inserted
///   sequence is the new residues followed by `Ter`
///   (`p.(Asn47delinsSerSerTer)`, `p.(Pro578_Lys579delinsLeuTer)`).
fn build_inframe_stop_delins(
    ref_protein: &[AminoAcid],
    alt_protein: &[AminoAcid],
    cds_pos_end: i64,
    protein_accession: &str,
    transcript: &Transcript,
) -> Result<HgvsVariant, FerroError> {
    let ref_len = ref_protein.len();
    // Range **start**: the first codon whose residue changed (protein-derived).
    let first_diff = first_diff_position(ref_protein, alt_protein);
    let start_idx = first_diff.min(ref_len.saturating_sub(1));
    // Range **end**: the last codon the edit deletes, `codon(cds_pos_end)`
    // (DNA-derived), never past the reference terminus. Start and end are on
    // different bases on purpose; `.max(start_idx)` keeps the range ascending
    // even in the (contrived) case where the leading deleted codons are
    // protein-silent so `first_diff` runs past the deleted span.
    let end_codon = ((cds_pos_end - 1).max(0) / 3) as usize;
    let end_idx = end_codon.max(start_idx).min(ref_len.saturating_sub(1));

    // Inserted residues: the new residues from the first change up to the stop,
    // plus the terminating `Ter` (dropped by translation).
    let mut inserted: Vec<AminoAcid> = alt_protein.get(start_idx..).unwrap_or(&[]).to_vec();
    inserted.push(AminoAcid::Ter);

    let start_aa = ref_protein
        .get(start_idx)
        .copied()
        .unwrap_or(AminoAcid::Xaa);
    let start_pos = (start_idx + 1) as u64;
    let accession = parse_accession(protein_accession);

    // An **immediate** stop (the inserted sequence is a bare `Ter` — alt ends
    // exactly at the first changed codon) is a nonsense substitution anchored at
    // that codon, regardless of how many codons the edit deletes (delins.md:27):
    // never the C-terminal-spanning `p.(X_YdelinsTer)` range form.
    if inserted.as_slice() == [AminoAcid::Ter] {
        let variant = ProteinVariant {
            accession,
            gene_symbol: transcript.gene_symbol.clone(),
            loc_edit: LocEdit::new_predicted(
                ProtInterval::point(ProtPos::new(start_aa, start_pos)),
                ProteinEdit::Substitution {
                    reference: start_aa,
                    alternative: AminoAcid::Ter,
                },
            ),
        };
        return Ok(HgvsVariant::Protein(variant));
    }

    let loc = if start_idx == end_idx {
        ProtInterval::point(ProtPos::new(start_aa, start_pos))
    } else {
        let end_aa = ref_protein.get(end_idx).copied().unwrap_or(AminoAcid::Xaa);
        ProtInterval::new(
            ProtPos::new(start_aa, start_pos),
            ProtPos::new(end_aa, (end_idx + 1) as u64),
        )
    };
    let variant = ProteinVariant {
        accession,
        gene_symbol: transcript.gene_symbol.clone(),
        loc_edit: LocEdit::new_predicted(
            loc,
            ProteinEdit::Delins {
                sequence: AminoAcidSeq::new(inserted),
            },
        ),
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

    // Defense-in-depth (#615): a frameshift whose first changed residue is at
    // or beyond the reference protein's end has no indexable reference residue,
    // so the legacy `ref_protein.get(first_diff).unwrap_or(Xaa)` below would
    // emit the invalid `Xaa` token. That situation is a stop-region
    // readthrough: the divergence is the former terminator (or past it),
    // which is a C-terminal extension, not a frameshift. Route it to the
    // extension builder, which anchors at the stop position and reads the new
    // residue from the mutated sequence — never `Xaa`. (Normal frameshifts
    // diverge mid-protein with `first_diff < ref_protein.len()` and skip this.)
    if first_diff >= ref_protein.len() {
        return build_extension_variant(
            ref_protein,
            alt_protein,
            mut_cds,
            protein_accession,
            transcript,
        );
    }

    let aa_pos = (first_diff + 1) as u64; // 1-based
    let ref_aa = ref_protein
        .get(first_diff)
        .copied()
        .unwrap_or(AminoAcid::Xaa);
    let new_aa = alt_protein.get(first_diff).copied();

    // Scan downstream from first_diff in the mutated CDS to find a stop codon.
    // A computed frameshift always carries an explicit termination state: the
    // new stop position when found (`fsTer{K}`), else `fsTer?` (no new stop in
    // the available sequence). The predictor never emits the bare short form
    // `fs` — that is reserved for human-authored input with no detail. The
    // degenerate `else` (frameshift starts at/beyond the CDS end) is likewise
    // an unknown-termination case (frameshift.md:44; spec example p.Ile327Argfs*?).
    // Scan for the new stop over the mutated CDS *plus* the unchanged 3'UTR:
    // a frameshift's new in-frame stop commonly lies past the annotated CDS
    // end. Scanning only the CDS-truncated `mut_cds` mis-emits `fsTer?`.
    let scan_seq = mut_cds_with_3utr(mut_cds, transcript)?;
    let codon_byte_offset = first_diff * 3;
    let ter: FrameshiftTer = if codon_byte_offset < scan_seq.len() {
        let downstream = &scan_seq[codon_byte_offset..];
        let downstream_aas = translate_full_cds_with_stop(downstream);
        // fsTer{K}: K is the 1-based position of the Ter in the downstream
        // scan (counting the new shifted amino acid as position 1).
        downstream_aas
            .iter()
            .position(|aa| *aa == AminoAcid::Ter)
            .map(|p| FrameshiftTer::At((p + 1) as u64))
            .unwrap_or(FrameshiftTer::Unknown)
    } else {
        FrameshiftTer::Unknown
    };

    // Spec (frameshift.md:22, 36-39): `fsTer1` is illegal — the minimum
    // frameshift count is `fsTer2`. A shifted frame whose first changed codon
    // is immediately a stop is a NONSENSE substitution `p.(<aa><pos>Ter)`, not
    // a frameshift. Emit the nonsense form (mirrors `predict_substitution`'s
    // `Ter` substitution) rather than the invalid `fsTer1`.
    if matches!(ter, FrameshiftTer::At(1)) {
        let nonsense_edit = ProteinEdit::Substitution {
            reference: ref_aa,
            alternative: AminoAcid::Ter,
        };
        let loc = ProtInterval::point(ProtPos::new(ref_aa, aa_pos));
        let accession = parse_accession(protein_accession);
        return Ok(HgvsVariant::Protein(ProteinVariant {
            accession,
            gene_symbol: transcript.gene_symbol.clone(),
            loc_edit: LocEdit::new_predicted(loc, nonsense_edit),
        }));
    }

    let protein_edit = ProteinEdit::Frameshift { new_aa, ter };
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

/// Does the edit produce a stop-codon readthrough (stop-loss extension)?
///
/// A stop-loss extension occurs when the edit disrupts the terminator codon
/// while leaving every sense residue intact, so the former stop now codes an
/// amino acid and translation continues into the 3'UTR. We detect it by
/// translating the *full* mutated read-through sequence (CDS + 3'UTR) and
/// requiring all of:
///
/// 1. the terminator is genuinely destroyed — either the edit is a frameshift
///    (`net % 3 != 0`, which scrambles the reading frame through the stop) or
///    its affected CDS span overlaps the terminator codon
///    (`cds_pos_end >= stop_cds_start`). Without this, an in-frame duplication
///    of a sense codon *before* the stop (e.g. `c.4_6dup` on Met-Arg-Ter) also
///    leaves the sense residues intact and pushes a residue into the
///    former-stop slot, yet the terminator is untouched (just relocated by one
///    codon) — that is a plain `dup`, not a stop-loss.
/// 2. the first `ref_protein.len()` translated residues are byte-identical to
///    `ref_protein` (no sense residue changed — a frameshift that begins
///    *before* the last sense codon diverges earlier and is excluded), and
/// 3. there is a residue at the former-stop slot (index `ref_protein.len()`)
///    that is not a terminator (the stop was genuinely read through).
///
/// Reading through the 3'UTR is essential: a small deletion/duplication that
/// disrupts the stop shortens (or doesn't lengthen) the CDS-only translation,
/// so a length-based check would miss the readthrough and misroute to the
/// frameshift path — which then reads the (nonexistent) reference residue at
/// the empty stop slot and emits the invalid `Xaa` token (#615).
fn is_stop_loss_readthrough(
    ref_protein: &[AminoAcid],
    mut_cds: &str,
    transcript: &Transcript,
    cds_pos_end: i64,
    is_frameshift: bool,
) -> Result<bool, FerroError> {
    // 1-based CDS position of the (reference) terminator codon.
    let stop_cds_start = (ref_protein.len() * 3 + 1) as i64;
    // A frameshift scrambles the frame through the stop even when its anchor
    // lies a base or two upstream of the terminator codon; an in-frame edit
    // only disrupts the stop if its span actually overlaps it.
    if !is_frameshift && cds_pos_end < stop_cds_start {
        return Ok(false);
    }

    let scan_seq = mut_cds_with_3utr(mut_cds, transcript)?;
    let mut readthrough = translate_full_cds_with_stop(&scan_seq);
    // `ref_protein` had residue 1 forced to the initiator `Met` for a recognized
    // initiator (`RefProteinBundle::from_transcript`, #801); apply the same
    // normalization to the freshly re-translated read-through so the sense-prefix
    // comparison below stays consistent. `scan_seq` shares the CDS start codon
    // (a downstream stop-loss edit never touches codon 0), so the gate matches.
    if cds_has_recognized_start(&scan_seq) {
        force_initiator_met(&mut readthrough);
    }
    let stop_idx = ref_protein.len();
    // The former-stop slot must hold a non-terminator residue (the stop was
    // read through) and every sense residue before it must be unchanged.
    match readthrough.get(stop_idx) {
        Some(&aa) if aa != AminoAcid::Ter => Ok(readthrough[..stop_idx] == *ref_protein),
        _ => Ok(false),
    }
}

/// Build a `p.(Ter{N}{Yyy}ext*{K})` extension variant for stop-codon readthrough.
fn build_extension_variant(
    ref_protein: &[AminoAcid],
    _alt_protein: &[AminoAcid],
    mut_cds: &str,
    protein_accession: &str,
    transcript: &Transcript,
) -> Result<HgvsVariant, FerroError> {
    // 1-based protein position of the original stop codon; the read-through
    // tail is computed by the shared C-terminal extension builder.
    let stop_pos = (ref_protein.len() + 1) as u64;

    let scan_seq = mut_cds_with_3utr(mut_cds, transcript)?;
    build_cterminal_extension(stop_pos, &scan_seq, protein_accession, transcript)
}

/// Predict a C-terminal extension for a **deletion that spans the CDS→3'UTR
/// boundary** (5' end in the CDS, 3' end a `*N` position). Unlike the in-CDS
/// stop-loss path (`build_extension_variant` / `is_stop_loss_readthrough`, which
/// append the *unchanged* 3'UTR via `mut_cds_with_3utr`), the deleted bases here
/// include 3'UTR positions — so we splice the edit against the **extended
/// reference** (CDS ++ 3'UTR) ONCE and translate that directly. Routing through
/// the in-CDS helpers would double-append the 3'UTR and reintroduce the deleted
/// `*N` base (#857).
///
/// `cds_pos_start` is the 1-based CDS position of the 5' end; `cds_end_utr3_n`
/// is the `N` of the `*N` 3' end. Returns `Ok(Some(p.(Ter<pos><aa>extTer<k>)))`
/// on a clean stop-loss read-through (sense residues intact, former stop now
/// codes an amino acid), else `Ok(None)`.
pub(crate) fn predict_stop_region_extension(
    transcript: &Transcript,
    ref_bundle: &RefProteinBundle,
    cds_pos_start: i64,
    cds_end_utr3_n: i64,
    edit: &NaEdit,
    protein_accession: &str,
) -> Result<Option<HgvsVariant>, FerroError> {
    // Precondition: the reference CDS must end in a terminator (mirrors the
    // `ref_protein_with_stop.last() == Ter` gate the in-CDS path applies).
    if ref_bundle.ref_protein_with_stop.last() != Some(&AminoAcid::Ter) {
        return Ok(None);
    }
    // Extended reference: CDS start through end of transcript (CDS ++ 3'UTR), in
    // transcript orientation, uppercase. `cds_start` is 1-based → 0-based index
    // `cds_start - 1` (same convention as `read_full_cds`).
    let seq =
        transcript
            .sequence
            .as_deref()
            .ok_or_else(|| FerroError::ProteinSequenceUnavailable {
                accession: transcript.id.clone(),
            })?;
    let cds_start = transcript
        .cds_start
        .ok_or_else(|| FerroError::ConversionError {
            msg: format!("transcript {} has no CDS start", transcript.id),
        })? as usize;
    if cds_start < 1 || cds_start > seq.len() {
        return Ok(None);
    }
    let extended = seq[cds_start - 1..].to_ascii_uppercase();

    // Convert the `*N` 3' end to an absolute 1-based CDS-relative position:
    // `ref_cds` covers the CDS incl. stop, so `*N` lands at `ref_cds.len() + N`.
    let abs_cds_pos_end = ref_bundle.ref_cds.len() as i64 + cds_end_utr3_n;

    // Splice the edit against the extended reference (no overrun, no double-3'UTR).
    let mut_full =
        build_mutated_cds_with_ref(&extended, transcript, cds_pos_start, abs_cds_pos_end, edit)?;

    // Stop-loss check directly on `mut_full` (replicates `is_stop_loss_readthrough`
    // without re-appending the 3'UTR). `force_initiator_met` keeps the sense-prefix
    // compare consistent with `ref_protein` (residue 0 was forced to Met for a
    // recognized initiator, #801).
    let mut readthrough = translate_full_cds_with_stop(&mut_full);
    if cds_has_recognized_start(&mut_full) {
        force_initiator_met(&mut readthrough);
    }
    let stop_idx = ref_bundle.ref_protein.len();
    // Former-stop slot holds a non-Ter residue AND every sense residue before it
    // is unchanged. (The `.get(stop_idx)` match already bounds-checks.)
    let is_readthrough = match readthrough.get(stop_idx) {
        Some(&aa) if aa != AminoAcid::Ter => readthrough[..stop_idx] == *ref_bundle.ref_protein,
        _ => false,
    };
    if !is_readthrough {
        return Ok(None);
    }

    let stop_pos = stop_idx as u64 + 1;
    Ok(Some(build_cterminal_extension(
        stop_pos,
        &mut_full,
        protein_accession,
        transcript,
    )?))
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

    fn prot_str(v: &HgvsVariant) -> String {
        match v {
            HgvsVariant::Protein(p) => p.to_string(),
            _ => panic!("expected Protein variant"),
        }
    }

    /// Regression (#498): `build_inframe_delins` must never emit a descending
    /// (start > end) protein range. A pure insertion mis-routed here — alt is
    /// `ref` with residues inserted between two unchanged residues — used to
    /// produce a crossed range like `p.(Gly3_Ala2delinsIle)` because the
    /// forward `first_diff` scan and the backward `find_last_diff` scan cross
    /// over (`last_diff_ref < first_diff`). The spec-correct rendering is an
    /// insertion with ascending positions.
    #[test]
    fn inframe_delins_renders_pure_insertion_as_ascending_ins() {
        use crate::hgvs::location::AminoAcid::{Ala, Gly, Ile, Met};
        let t = tx("ATGGCTGGTTAA", 1, 12); // sequence unused by the position logic
        let ref_protein = [Met, Ala, Gly];
        let alt_protein = [Met, Ala, Ile, Gly]; // Ile inserted between Ala(2) and Gly(3)
        let v = build_inframe_delins(&ref_protein, &alt_protein, 5, "NP_TEST.1", &t).unwrap();
        assert_eq!(prot_str(&v), "NP_TEST.1:p.(Ala2_Gly3insIle)");
    }

    /// #855 (delins.md): one residue replaced by one residue is a substitution,
    /// not a delins — even when routed through `build_inframe_delins`.
    #[test]
    fn inframe_delins_single_residue_renders_substitution() {
        use crate::hgvs::location::AminoAcid::{Arg, Gly, Met, Val};
        let t = tx("ATGCGCGGTTAA", 1, 12);
        let ref_protein = [Met, Arg, Gly];
        let alt_protein = [Met, Val, Gly]; // Arg2 → Val (1 AA → 1 AA)
        let v = build_inframe_delins(&ref_protein, &alt_protein, 4, "NP_TEST.1", &t).unwrap();
        assert_eq!(prot_str(&v), "NP_TEST.1:p.(Arg2Val)");
    }

    /// #855 (duplication.md): a pure insertion that is a tandem copy of the
    /// residue directly C-terminal (immediately 5') of the insertion point is a
    /// duplication, not an insertion. Mirrors the corpus `c.15_17dup`
    /// → `p.(Cys6dup)` shape (mid-codon dup routed through the delins guard).
    #[test]
    fn inframe_delins_tandem_single_residue_renders_dup() {
        use crate::hgvs::location::AminoAcid::{Arg, Cys, Met};
        let t = tx("ATGTGCCGCTAA", 1, 12);
        let ref_protein = [Met, Cys, Arg];
        let alt_protein = [Met, Cys, Cys, Arg]; // extra Cys inserted after Cys2
        let v = build_inframe_delins(&ref_protein, &alt_protein, 6, "NP_TEST.1", &t).unwrap();
        assert_eq!(prot_str(&v), "NP_TEST.1:p.(Cys2dup)");
    }

    /// #855 (duplication.md 3' rule): for a tandem run the **most C-terminal**
    /// residue is assigned as duplicated. Inserting an extra Ala into an Ala run
    /// must report the 3'-most Ala, not the 5'-most.
    #[test]
    fn inframe_delins_tandem_run_picks_most_cterminal() {
        use crate::hgvs::location::AminoAcid::{Ala, Gly, Met};
        let t = tx("ATGGCAGCAGGTTAA", 1, 15);
        let ref_protein = [Met, Ala, Ala, Gly]; // Ala at positions 2 and 3
        let alt_protein = [Met, Ala, Ala, Ala, Gly]; // one extra Ala in the run
        let v = build_inframe_delins(&ref_protein, &alt_protein, 5, "NP_TEST.1", &t).unwrap();
        assert_eq!(prot_str(&v), "NP_TEST.1:p.(Ala3dup)");
    }

    /// #855 regression guard: a multi-residue tandem copy renders as a ranged dup.
    #[test]
    fn inframe_delins_tandem_two_residues_renders_ranged_dup() {
        use crate::hgvs::location::AminoAcid::{Gly, His, Met, Ser};
        let t = tx("ATGGGCTCCCACTAA", 1, 15);
        let ref_protein = [Met, Gly, Ser, His];
        let alt_protein = [Met, Gly, Ser, Gly, Ser, His]; // GlySer copied in tandem
        let v = build_inframe_delins(&ref_protein, &alt_protein, 8, "NP_TEST.1", &t).unwrap();
        assert_eq!(prot_str(&v), "NP_TEST.1:p.(Gly2_Ser3dup)");
    }

    /// #855 regression guard: a non-tandem insertion (inserted residue does not
    /// match the directly-5' residue) must stay an insertion, not become a dup.
    #[test]
    fn inframe_delins_non_tandem_insertion_stays_insertion() {
        use crate::hgvs::location::AminoAcid::{Ala, Gly, Ile, Met};
        let t = tx("ATGGCTGGTTAA", 1, 12);
        let ref_protein = [Met, Ala, Gly];
        let alt_protein = [Met, Ala, Ile, Gly]; // Ile ≠ Ala(2) → genuine insertion
        let v = build_inframe_delins(&ref_protein, &alt_protein, 5, "NP_TEST.1", &t).unwrap();
        assert_eq!(prot_str(&v), "NP_TEST.1:p.(Ala2_Gly3insIle)");
    }

    /// #855 bounds guard: an insertion near the N-terminus (fewer preceding
    /// residues than the inserted unit) must not panic and must stay an insertion.
    #[test]
    fn protein_tandem_dup_range_n_terminal_no_panic() {
        use crate::hgvs::location::AminoAcid::{Ala, Gly, Met};
        // Insert [Gly, Ala] at index 1 (only one residue, Met, precedes it).
        let ref_protein = [Met, Ala];
        assert_eq!(protein_tandem_dup_range(&ref_protein, 1, &[Gly, Ala]), None);
        // Empty insertion and out-of-range point are also None.
        assert_eq!(protein_tandem_dup_range(&ref_protein, 1, &[]), None);
        assert_eq!(protein_tandem_dup_range(&ref_protein, 99, &[Met]), None);
    }

    /// Regression (#847): `build_inframe_delins` called directly with a `ref`/
    /// `alt` pair whose backward (`find_last_diff`) and forward (`first_diff`)
    /// scans cross over (`last_diff_ref < first_diff`) AND whose `alt` is
    /// *shorter* than `ref` must render a `del`, never an empty `ins`. This
    /// locks the pure-deletion guard's invariant: when it fires, `alt` is
    /// provably `ref` with a contiguous block removed, so `del` is exact. Here
    /// one Leu is removed from a two-Leu run between Ala and Gly →
    /// `p.(Leu4del)`. Calling the builder directly (bypassing edit-type routing)
    /// proves the guard, not just the `Deletion` entry path.
    #[test]
    fn inframe_delins_renders_pure_deletion_as_del() {
        use crate::hgvs::location::AminoAcid::{Ala, Gly, Leu, Met};
        let t = tx("ATGGCTCTTCTTGGTTAA", 1, 18); // sequence unused by the position logic
        let ref_protein = [Met, Ala, Leu, Leu, Gly];
        let alt_protein = [Met, Ala, Leu, Gly]; // one Leu removed from the Leu-Leu run
        let v = build_inframe_delins(&ref_protein, &alt_protein, 8, "NP_TEST.1", &t).unwrap();
        assert_eq!(prot_str(&v), "NP_TEST.1:p.(Leu4del)");
    }

    /// Test convenience wrapper: build the `RefProteinBundle` on demand so
    /// existing tests can keep their original `predict_indel_protein` calls.
    fn predict_indel(
        t: &Transcript,
        cds_pos_start: i64,
        cds_pos_end: i64,
        edit: &NaEdit,
        protein_accession: &str,
    ) -> Result<HgvsVariant, FerroError> {
        let bundle = RefProteinBundle::from_transcript(t)?;
        predict_indel_protein(
            t,
            &bundle,
            cds_pos_start,
            cds_pos_end,
            edit,
            protein_accession,
        )
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
        let result = predict_indel(&t, 4, 4, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        // alt_protein = [Met, Ala], first_diff=1, new_aa=Ala. The downstream
        // scan from byte 3 ("GCTAA" → Ala, then an incomplete codon) finds no
        // new stop, so the shifted frame's termination is unknown: the
        // predictor must emit the explicit `fsTer?` marker (not bare `fs`),
        // matching the spec's no-stop example p.Ile327Argfs*? (frameshift.md:44).
        assert_eq!(s, "NP_TEST.1:p.(Arg2AlafsTer?)");
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
        let result = predict_indel(&t, 4, 6, &edit, "NP_TEST.1").unwrap();
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
        let result = predict_indel(&t, 4, 9, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert_eq!(s, "NP_TEST.1:p.(Arg2_Lys3del)");
    }

    /// Regression (#847): an in-frame single-codon deletion that straddles a
    /// codon boundary (`frame_at_start != 0`) inside a run of identical
    /// residues must render as a single-residue `del`, NOT a malformed `ins`
    /// with an empty inserted sequence.
    ///
    /// CDS "ATG" + "CTG"×5 + "GGG" + "TAA" = Met-Leu-Leu-Leu-Leu-Leu-Gly-Ter.
    /// Delete c.6_8 (3 bp, codon-straddling: `(6-1)%3 == 2`) which removes the
    /// "G C T" spanning the 1st Leu codon's last base and the 2nd Leu codon's
    /// first two bases — net effect is one Leu codon removed from the run:
    ///   mutated CDS "ATGCTGCTGCTGCTGGGGTAA" = Met-Leu-Leu-Leu-Leu-Gly-Ter.
    /// ref = [Met,Leu,Leu,Leu,Leu,Leu,Gly], alt = [Met,Leu,Leu,Leu,Leu,Gly].
    /// The 3'-most Leu is deleted ⇒ p.(Leu6del). Before the fix this misrouted
    /// through `build_inframe_delins`'s pure-insertion guard and emitted the
    /// invalid `p.(Leu5_Leu6ins)` (empty inserted sequence).
    #[test]
    fn del_single_codon_in_homopolymer_run_renders_del_not_empty_ins() {
        let t = tx("ATGCTGCTGCTGCTGCTGGGGTAA", 1, 24);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = predict_indel(&t, 6, 8, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert_eq!(s, "NP_TEST.1:p.(Leu6del)");
        assert!(
            !s.contains("ins"),
            "single-codon in-frame deletion must not render as an insertion: '{}'",
            s
        );
    }

    /// Regression (#847): a two-codon in-frame deletion straddling a codon
    /// boundary inside a homopolymer run must render as a range `del`, not an
    /// empty `ins`/`delins`.
    ///
    /// Same fixture as the single-codon case but delete c.6_11 (6 bp,
    /// codon-straddling): removes two Leu codons from the run.
    ///   mutated CDS "ATGCTGCTGCTGGGGTAA" = Met-Leu-Leu-Leu-Gly-Ter.
    /// ref = [Met,Leu,Leu,Leu,Leu,Leu,Gly], alt = [Met,Leu,Leu,Leu,Gly].
    /// The two 3'-most Leu are deleted ⇒ p.(Leu5_Leu6del).
    #[test]
    fn del_two_codons_in_homopolymer_run_renders_range_del() {
        let t = tx("ATGCTGCTGCTGCTGCTGGGGTAA", 1, 24);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = predict_indel(&t, 6, 11, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert_eq!(s, "NP_TEST.1:p.(Leu5_Leu6del)");
        assert!(
            !s.contains("ins"),
            "two-codon in-frame deletion must not render as an insertion: '{}'",
            s
        );
    }

    /// Companion to the #847 fix: a codon-straddling in-frame deletion that is
    /// NOT inside a run still changes a residue at the AA level and must remain
    /// a `delins` (its differing ref region is non-empty), unaffected by the
    /// pure-deletion guard.
    ///
    /// CDS "ATGCGCAAAGGGTAA" = Met-Arg-Lys-Gly-Ter. Delete c.5_7 ("GCA",
    /// codon-straddling): mutated CDS "ATGCAAGGGTAA"? Recompute bases —
    ///   ATG C G C A A A G G G T A A (positions 1..15)
    ///   delete 5,6,7 (G,C,A) → ATG C | A A G G G T A A = "ATGCAAGGGTAA"
    ///   = ATG CAA GGG TAA = Met-Gln-Gly-Ter.
    /// ref = [Met,Arg,Lys,Gly], alt = [Met,Gln,Gly]: Arg2/Lys3 collapse to a
    /// single Gln ⇒ a genuine delins, p.(Arg2_Lys3delinsGln).
    #[test]
    fn del_codon_straddling_not_in_run_stays_delins() {
        let t = tx("ATGCGCAAAGGGTAA", 1, 15);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = predict_indel(&t, 5, 7, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert_eq!(s, "NP_TEST.1:p.(Arg2_Lys3delinsGln)");
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
        let result = predict_indel(&t, 4, 5, &edit, "NP_TEST.1").unwrap();
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
        let result = predict_indel(&t, 3, 3, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert_eq!(s, "NP_TEST.1:p.(Met1_Arg2insGly)");
    }

    /// A mid-codon in-frame insertion that disrupts the flanking codon's residue
    /// must be a delins, not a clean insertion (#511). Uses the real caller's
    /// `A_(A+1)` interval convention (`cds_pos_end == cds_pos_start + 1`).
    ///
    /// CDS "ATGATGGGCTAA" = Met-Met-Gly-Stop. Insert "ATC" at c.5_6 (between the
    /// 2nd and 3rd base of Met2's codon): ATG AT|ATC|G GGC TAA → ATG ATA TCG GGC
    /// TAA → Met-Ile-Ser-Gly. So Met2 changes to Ile and Ser is added (Gly3
    /// preserved) → p.(Met2delinsIleSer), NOT the spurious p.(Met2_Gly3insIle).
    #[test]
    fn ins_mid_codon_disrupting_flanking_residue_is_delins() {
        let t = tx("ATGATGGGCTAA", 1, 12);
        let seq: crate::hgvs::edit::Sequence = "ATC".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
        };
        let result = predict_indel(&t, 5, 6, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Met2delinsIleSer)");
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
        let result = predict_indel(&t, 3, 3, &edit, "NP_TEST.1").unwrap();
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
        let result = predict_indel(&t, 4, 6, &edit, "NP_TEST.1").unwrap();
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
        let result = predict_indel(&t, 4, 4, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert!(s.contains("fs"), "expected frameshift in '{}'", s);
    }

    #[test]
    fn delins_same_length_codon_aligned() {
        // CDS "ATGCGCTAA": delins c.4_6 (CGC → TCC).
        // Mutated CDS: "ATGTCCTAA" → ATG TCC TAA → [Met, Ser]
        // ref=[Met, Arg], alt=[Met, Ser].
        // One amino acid replaced by one other amino acid is a SUBSTITUTION, not a
        // delins (delins.md:17), even though the DNA-level edit is a delins → so the
        // spec-correct protein consequence is p.(Arg2Ser), not p.(Arg2delinsSer).
        let t = tx("ATGCGCTAA", 1, 9);
        let seq: crate::hgvs::edit::Sequence = "TCC".parse().unwrap();
        let edit = NaEdit::Delins {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
            deleted: None,
            deleted_length: None,
        };
        let result = predict_indel(&t, 4, 6, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Arg2Ser)");
    }

    #[test]
    fn inversion_codon_cgc_to_gcg() {
        // CDS "ATGCGCTAA": inv c.4_6 (CGC → GCG).
        // GCG = Ala. Mutated CDS: "ATGGCGTAA" → [Met, Ala].
        // A whole-codon inversion that changes a single residue is a SUBSTITUTION
        // per delins.md:17 (1 AA → 1 AA), not a delins → p.(Arg2Ala). (The spec
        // says inversions are rendered as delins on the protein level only when
        // they span 2+ residues.)
        let t = tx("ATGCGCTAA", 1, 9);
        let edit = NaEdit::Inversion {
            sequence: None,
            length: None,
        };
        let result = predict_indel(&t, 4, 6, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Arg2Ala)");
    }

    // ─── Regression tests for CodeRabbit-flagged edge cases ───────────────────

    #[test]
    fn inframe_deletion_with_zero_diff_returns_identity_not_panic() {
        // Regression for the underflow in build_inframe_deletion: if it is ever
        // called with alt_protein.len() >= ref_protein.len(), `n_deleted` is 0
        // and `first_diff + n_deleted - 1` would underflow usize. The guard
        // must run before the subtraction. Calling directly bypasses upstream
        // detection so we are testing the function's own safety.
        let t = tx("ATGCGCTAA", 1, 9);
        let ref_protein = [AminoAcid::Met, AminoAcid::Arg];
        let alt_protein = [AminoAcid::Met, AminoAcid::Arg];
        let result =
            build_inframe_deletion(&ref_protein, &alt_protein, "NP_TEST.1", &t).expect("no panic");
        let s = prot_str(&result);
        // Identity variant -- whole-protein "=" representation.
        assert!(s.contains("(="), "expected identity '(=' in '{}'", s);
    }

    #[test]
    fn ins_codon_aligned_immediate_stop_renders_nonsense_substitution() {
        // #911: insert the stop codon "TAA" between c.3 (last base of Met codon)
        // and c.4 (first base of Arg codon). The inserted codon is a bare stop
        // (no sense residue before it), so translation truncates alt_protein to
        // [Met]. ref_protein = [Met, Arg, Lys] (CDS "ATGCGCAAATAA"). Per
        // insertion.md:25 an *immediate* stop is a nonsense variant — the stop
        // lands at the residue 3' of the insertion point: `p.(Arg2Ter)`. NOT the
        // `ins…Ter` form (reserved for a sense-carrying insert, insertion.md:26)
        // and NOT the spec-forbidden C-terminal delins.
        let t = tx("ATGCGCAAATAA", 1, 12);
        let seq: crate::hgvs::edit::Sequence = "TAA".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
        };
        let result = predict_indel(&t, 3, 3, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Arg2Ter)");
    }

    #[test]
    fn ins_codon_aligned_residues_then_stop_renders_insertion_with_ter() {
        // #911: insert "GCGTAA" (Ala + stop) between c.3 and c.4. New CDS
        // "ATG GCG TAA CGCAAATAA" = Met-Ala-Stop; alt truncates to [Met, Ala].
        // ref = [Met, Arg, Lys]. Per insertion.md:41 (`p.(Pro46_Asn47insSerSerTer)`)
        // the inserted residues up to and including the Ter are listed, using the
        // reference flanking positions: `p.(Met1_Arg2insAlaTer)`.
        let t = tx("ATGCGCAAATAA", 1, 12);
        let seq: crate::hgvs::edit::Sequence = "GCGTAA".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
        };
        let result = predict_indel(&t, 3, 3, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Met1_Arg2insAlaTer)");
    }

    #[test]
    fn ins_codon_aligned_late_stop_longer_than_ref_tail_keeps_ter() {
        // #911 regression (CodeRabbit MAJOR): a codon-aligned insertion can
        // introduce a premature stop while the truncated `alt_protein` is still
        // LONGER than `ref_protein` — when the inserted sense residues before the
        // new stop outnumber the reference residues after the insertion point.
        // Deciding the stop branch from `alt.len() > ref.len()` misclassifies this
        // as a plain insertion and drops the terminating `Ter`; the caller's
        // `premature_stop` flag must be used instead.
        //
        // CDS "ATGCGCTAA" = Met-Arg-Ter; ref_protein = [Met, Arg] (len 2). Insert
        // "GCGGCGTAA" (Ala-Ala-Stop) between c.3 and c.4 → new CDS
        // "ATG GCGGCGTAA CGCTAA"; alt truncates to [Met, Ala, Ala] (len 3 > 2).
        // Per insertion.md:37-41 the inserted residues up to and including the Ter
        // are listed at the reference flanking positions: `p.(Met1_Arg2insAlaAlaTer)`.
        let t = tx("ATGCGCTAA", 1, 9);
        let seq: crate::hgvs::edit::Sequence = "GCGGCGTAA".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
        };
        let result = predict_indel(&t, 3, 3, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Met1_Arg2insAlaAlaTer)");
    }

    #[test]
    fn delins_immediate_stop_renders_nonsense_substitution() {
        // #911: a delins whose inserted sequence starts with a stop is a nonsense
        // substitution, not a C-terminal-spanning (or empty) delins. CDS
        // "ATGAAAGGGCCCTAA" = Met-Lys-Gly-Pro-Stop; `c.4_6delinsTAAGGG` replaces
        // codon 2 (Lys) → new CDS "ATG TAAGGG GGGCCCTAA" = Met-Stop.
        // delins.md:27 / deletion.md:22-24: `p.(Lys2Ter)`, not
        // `p.(Lys2_Pro4delins)`.
        let t = tx("ATGAAAGGGCCCTAA", 1, 15);
        let seq: crate::hgvs::edit::Sequence = "TAAGGG".parse().unwrap();
        let edit = NaEdit::Delins {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
            deleted: None,
            deleted_length: None,
        };
        let result = predict_indel(&t, 4, 6, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Lys2Ter)");
    }

    #[test]
    fn delins_multi_codon_immediate_stop_renders_nonsense_substitution() {
        // #911 (MAJOR-1 guard): an immediate stop whose DNA edit deletes TWO
        // codons must still collapse to a nonsense substitution, not a range
        // delins `p.(Lys2_Gly3delinsTer)`. CDS "ATGAAAGGGCCCTAA"; `c.4_9delinsTAAAAA`
        // replaces codons 2-3 → new CDS "ATG TAAAAA CCCTAA" = Met-Stop.
        // delins.md:27: `p.(Lys2Ter)`, anchored at the first changed codon.
        let t = tx("ATGAAAGGGCCCTAA", 1, 15);
        let seq: crate::hgvs::edit::Sequence = "TAAAAA".parse().unwrap();
        let edit = NaEdit::Delins {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
            deleted: None,
            deleted_length: None,
        };
        let result = predict_indel(&t, 4, 9, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Lys2Ter)");
    }

    #[test]
    fn delins_residues_then_stop_caps_at_change_site() {
        // #911: a delins whose inserted sequence encodes residues then a stop is
        // a delins capped at the changed codon(s), ending in Ter — not spanning
        // the C-terminal (delins.md:27-28, `p.(Asn47delinsSerSerTer)`). CDS
        // "ATGAAAGGGCCCTAA" = Met-Lys-Gly-Pro-Stop; `c.7_9delinsAGCTAA` replaces
        // codon 3 (Gly) → new CDS "ATGAAA AGCTAA CCCTAA" = Met-Lys-Ser-Stop.
        let t = tx("ATGAAAGGGCCCTAA", 1, 15);
        let seq: crate::hgvs::edit::Sequence = "AGCTAA".parse().unwrap();
        let edit = NaEdit::Delins {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
            deleted: None,
            deleted_length: None,
        };
        let result = predict_indel(&t, 7, 9, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Gly3delinsSerTer)");
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
        let result = predict_indel(&t, 7, 9, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert!(s.contains("Ter3"), "expected Ter3 in '{}'", s);
        assert!(s.contains("Trp"), "expected Trp in '{}'", s);
        assert!(s.contains("ext"), "expected ext in '{}'", s);
    }

    #[test]
    fn stop_readthrough_extension_on_non_aug_start() {
        // #801 regression: stop-loss readthrough on a non-AUG-initiation
        // transcript must still route to the C-terminal extension path. The
        // RefProteinBundle now forces residue 1 to Met for a recognized
        // initiator (CTG here), so `is_stop_loss_readthrough` must apply the
        // same Met1 normalization to the re-translated read-through sequence
        // before comparing — otherwise readthrough[0] (Leu, from the literal
        // CTG) would not equal ref_protein[0] (Met) and the variant would
        // misroute out of the extension path.
        //
        // CDS "CTGCGCTAA": Met(forced)-Arg-Ter; stop "TAA" at c.7_9.
        // delins c.7_9 TGG converts the stop to Trp → readthrough, no further
        // stop in CDS → p.(Ter3Trpext*?).
        let t = tx("CTGCGCTAA", 1, 9);
        let seq: crate::hgvs::edit::Sequence = "TGG".parse().unwrap();
        let edit = NaEdit::Delins {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
            deleted: None,
            deleted_length: None,
        };
        let result = predict_indel(&t, 7, 9, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert!(
            s.contains("Ter3") && s.contains("Trp") && s.contains("ext"),
            "expected a C-terminal extension (Ter3...ext) on a CTG-start \
             transcript, got '{}'",
            s
        );
    }

    #[test]
    fn stop_readthrough_extension_into_3utr() {
        // Stop-loss read-through whose new stop lies in the 3'UTR, not the
        // annotated CDS. The mutated CDS alone has no further stop, so the
        // extension count can only be found by reading through into the
        // downstream transcript sequence. Regression for passing a CDS-only
        // `mut_cds` into `build_cterminal_extension` (degraded to `extTer?`).
        //
        // Transcript "ATGCGCTAAGCATAA": CDS c.1_9 = "ATGCGCTAA"
        // ([Met, Arg, Ter]); 3'UTR = "GCATAA". delins c.7_9 TGG converts the
        // stop "TAA" to "TGG" (Trp), so translation reads through into the
        // 3'UTR: "TGG"(Trp) "GCA"(Ala) "TAA"(Ter) → the new stop is the 2nd
        // added residue → extTer2 (the converted-stop residue is not counted).
        let t = tx("ATGCGCTAAGCATAA", 1, 9);
        let seq: crate::hgvs::edit::Sequence = "TGG".parse().unwrap();
        let edit = NaEdit::Delins {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
            deleted: None,
            deleted_length: None,
        };
        let result = predict_indel(&t, 7, 9, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert!(
            s.contains("Ter3Trp") && s.contains("extTer2"),
            "expected p.(Ter3TrpextTer2) reading into the 3'UTR, got '{}'",
            s
        );
    }

    // ─── Structural diversity tests (c17) ─────────────────────────────────────
    //
    // These cases pin behaviour at the structural edges of the CDS that the
    // existing tests don't reach: edits that remove the start codon, edits
    // that delete the entire stop, indels larger than 6 bases, mid-codon
    // insertions that introduce a new stop, and frame-restoring compound
    // edits.

    /// Initiator Met loss: deleting c.1_3 removes the entire start codon
    /// (ATG). HGVS recommendation: report as `p.(Met1?)` because we can't
    /// predict where translation actually starts on the altered transcript.
    ///
    /// Deleting the initiation codon yields an unpredictable protein
    /// consequence, reported as `p.(Met1?)` (HGVS deletion.md:62; #498) — not
    /// a concrete `Met1del`.
    #[test]
    fn del_start_codon_initiator_met_loss() {
        // CDS "ATGCGCTAA"; del c.1_3 removes the entire start codon (ATG).
        let t = tx("ATGCGCTAA", 1, 9);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = predict_indel(&t, 1, 3, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Met1?)");
    }

    /// Initiator Met loss via delins: replacing c.1_3 (the start codon) with
    /// any sequence has an unpredictable protein consequence, reported as
    /// `p.(Met1?)` rather than a concrete `Met1delins…` (#498). The
    /// `[start,end]` span overlaps CDS 1–3 so the guard fires for delins.
    #[test]
    fn delins_start_codon_initiator_met_loss() {
        // CDS "ATGCGCTAA"; delins c.1_3 (ATG → TCC) overwrites the start codon.
        let t = tx("ATGCGCTAA", 1, 9);
        let seq: crate::hgvs::edit::Sequence = "TCC".parse().unwrap();
        let edit = NaEdit::Delins {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
            deleted: None,
            deleted_length: None,
        };
        let result = predict_indel(&t, 1, 3, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Met1?)");
    }

    /// Initiator Met loss via duplication: duplicating c.1_3 (the start codon)
    /// disrupts initiation, reported as `p.(Met1?)` rather than a concrete
    /// `Met1dup` (#498). The `[start,end]` span overlaps CDS 1–3.
    #[test]
    fn dup_start_codon_initiator_met_loss() {
        // CDS "ATGCGCTAA"; dup c.1_3 (ATG) overlaps the start codon.
        let t = tx("ATGCGCTAA", 1, 9);
        let edit = NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: None,
        };
        let result = predict_indel(&t, 1, 3, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Met1?)");
    }

    /// Initiator Met loss via inversion: inverting c.1_3 (the start codon)
    /// disrupts initiation, reported as `p.(Met1?)` rather than a concrete
    /// `Met1delins…` (#498). The `[start,end]` span overlaps CDS 1–3.
    #[test]
    fn inversion_start_codon_initiator_met_loss() {
        // CDS "ATGCGCTAA"; inv c.1_3 (ATG → CAT) overwrites the start codon.
        let t = tx("ATGCGCTAA", 1, 9);
        let edit = NaEdit::Inversion {
            sequence: None,
            length: None,
        };
        let result = predict_indel(&t, 1, 3, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Met1?)");
    }

    /// Initiator Met loss via insertion *inside* the start codon: an insertion
    /// disrupts initiation only when the gap is strictly inside CDS 1–3, i.e.
    /// `cds_pos_start ∈ {1,2}`. An insertion at c.1_2 (start=1) splits the
    /// start codon and is reported as `p.(Met1?)` (#498).
    #[test]
    fn ins_within_start_codon_initiator_met_loss() {
        // CDS "ATGCGCTAA"; insert "GGG" at c.1_2 (gap between c.1 and c.2,
        // strictly inside the start codon ATG).
        let t = tx("ATGCGCTAA", 1, 9);
        let seq: crate::hgvs::edit::Sequence = "GGG".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
        };
        let result = predict_indel(&t, 1, 1, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Met1?)");
    }

    /// Initiator Met loss via insertion at the *second* interior gap of the
    /// start codon (gap between c.2 and c.3, `cds_pos_start = 2`). The helper
    /// treats both interior gaps (`cds_pos_start ∈ {1,2}`) as disrupting
    /// initiation, so this also reports `p.(Met1?)` (#498). Pins the
    /// `cds_pos_start == 2` branch that `ins_within_start_codon` (start=1) and
    /// `ins_after_start_codon` (start=3) leave uncovered, guarding against an
    /// off-by-one regression on the interior-gap range.
    #[test]
    fn ins_second_gap_within_start_codon_initiator_met_loss() {
        // CDS "ATGCGCTAA"; insert "GGG" at c.2_3 (gap between c.2 and c.3,
        // strictly inside the start codon ATG).
        let t = tx("ATGCGCTAA", 1, 9);
        let seq: crate::hgvs::edit::Sequence = "GGG".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
        };
        let result = predict_indel(&t, 2, 2, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Met1?)");
    }

    /// Insertion boundary: an insertion *after* the start codon (gap between
    /// c.3 and c.4, `cds_pos_start = 3`) does NOT disrupt initiation, so the
    /// guard must NOT fire — the normal in-frame insertion consequence is
    /// reported instead of `p.(Met1?)`. Pins the insertion-specific boundary
    /// of `affects_initiation_codon` (#498).
    #[test]
    fn ins_after_start_codon_does_not_trigger_initiator_unknown() {
        // CDS "ATGCGCTAA"; insert "GGG" at c.3_4 (gap just past the start
        // codon). Same fixture as ins_three_bases_codon_aligned.
        let t = tx("ATGCGCTAA", 1, 9);
        let seq: crate::hgvs::edit::Sequence = "GGG".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
        };
        let result = predict_indel(&t, 3, 3, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert_ne!(s, "NP_TEST.1:p.(Met1?)", "guard must not fire past CDS 3");
        assert_eq!(s, "NP_TEST.1:p.(Met1_Arg2insGly)");
    }

    /// Large multi-codon in-frame deletion (4 codons, 12 bases) — covers
    /// the codon-aligned `inframe_deletion` builder with N > 2.
    #[test]
    fn del_twelve_bases_four_codons_aligned() {
        // CDS for Met-Arg-Lys-Gly-Phe-Pro-Stop (21 bases).
        // Delete c.4_15 = "CGCAAAGGGTTT" (codons 2-5).
        // Mutated CDS: "ATG" + "CCATAA" — wait, c.16_21 = "CCATAA"? Let's
        // recompute. Original = "ATG CGC AAA GGG TTT CCA TAA" = 21 bases,
        // 7 codons (Met Arg Lys Gly Phe Pro Stop). Delete c.4_15 removes
        // 12 bases → leaves "ATG" + "CCATAA" = "ATGCCATAA" = Met-Pro-Stop.
        let t = tx("ATGCGCAAAGGGTTTCCATAA", 1, 21);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = predict_indel(&t, 4, 15, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert!(s.contains("del"), "expected del in '{}'", s);
        // Range deletion notation: Arg2_Phe5del or similar.
        assert!(s.contains("Arg2"), "expected Arg2 in '{}'", s);
        assert!(
            s.contains("Phe5") || s.contains("Phe") && s.contains("_"),
            "expected Phe5 (range end) in '{}'",
            s
        );
    }

    /// Mid-codon in-frame insertion that does NOT truncate the protein — the
    /// inserted "TAA" is distributed across codons so no premature stop appears
    /// before the reference terminus. Covers the plain mid-codon
    /// `inframe_delins` path (`premature_stop == false`), i.e.
    /// `build_inframe_delins`, not the stop-aware builder. The truncating
    /// counterpart is `ins_mid_codon_premature_stop_caps_at_change_site`.
    #[test]
    fn ins_mid_codon_no_premature_stop_renders_plain_delins() {
        // CDS "ATGCGCTGCAAATAA" = Met-Arg-Cys-Lys-Stop (15 bases, 5 codons).
        // Insert "TAA" at c.4_5 (mid-codon, between first and second base of
        // codon 2 / Arg = "CGC"). Mutated CDS:
        //   "ATG" + "C" + "TAA" + "GCTGCAAATAA" = "ATGCTAAGCTGCAAATAA"
        //   = ATG CTA AGC TGC AAA TAA = Met-Leu-Ser-Cys-Lys-Stop.
        // No premature stop in the codons before the original stop.
        // (Choosing this fixture so the inserted "TAA" doesn't sit in
        // codon position; it gets distributed across codons.)
        let t = tx("ATGCGCTGCAAATAA", 1, 15);
        let seq: crate::hgvs::edit::Sequence = "TAA".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
        };
        let result = predict_indel(&t, 4, 4, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        // Mid-codon ins of multiple of 3 bases is still in-frame; the
        // protein changes via a delins-at-the-AA-level rewrite. No stop was
        // introduced, so there must be no terminating `Ter` in the output.
        assert!(
            s.contains("delins") || s.contains("ins"),
            "expected delins/ins, got '{}'",
            s
        );
        assert!(
            !s.contains("Ter"),
            "no premature stop was introduced, so no Ter expected, got '{}'",
            s
        );
    }

    /// Mid-codon in-frame insertion that DOES introduce a premature stop,
    /// truncating `alt_protein` before the reference terminus. This is the
    /// route the previous test's name promised but never exercised: a
    /// `NaEdit::Insertion` at a non-codon-aligned position (`frame_at_start !=
    /// 2`) reaches `build_inframe_variant`'s mid-codon delins arm with
    /// `premature_stop == true`, so it must dispatch to
    /// `build_inframe_stop_delins` and be capped at the first changed codon
    /// ending in `Ter` — never a C-terminal-spanning delins (#911).
    #[test]
    fn ins_mid_codon_premature_stop_caps_at_change_site() {
        // CDS "ATGCGCTGCAAATAA" = Met-Arg-Cys-Lys-Stop (15 bases, 5 codons).
        // Insert "GCGAATAAG" (9 bases, in-frame) at c.4_5 (mid-codon, after the
        // first base of codon 2 / Arg = "CGC"). Mutated CDS:
        //   "ATG" + "C" + "GCGAATAAG" + "GCTGCAAATAA" = "ATGCGCGAATAAGGCTGCAAATAA"
        //   = ATG CGC GAA TAA GGC TGC AAA TAA = Met-Arg-Glu-Stop.
        // Translation truncates at the new stop: alt = [Met, Arg, Glu],
        // ref = [Met, Arg, Cys, Lys]. The first changed codon is Cys3, which
        // becomes Glu then an immediate stop ⇒ a delins capped at Cys3 ending
        // in Ter: p.(Cys3delinsGluTer). NOT a C-terminal-spanning delins.
        let t = tx("ATGCGCTGCAAATAA", 1, 15);
        let seq: crate::hgvs::edit::Sequence = "GCGAATAAG".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
        };
        let result = predict_indel(&t, 4, 4, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Cys3delinsGluTer)");
    }

    /// Frame-restoring compound: an insertion of net 0 mod 3 that crosses
    /// a codon boundary. Specifically a 6-base insertion mid-codon — still
    /// in-frame at the protein level, but every codon downstream of the
    /// insertion shifts by 6 bases.
    #[test]
    fn ins_six_bases_mid_codon_inframe() {
        // CDS "ATGCGCTGCAAATAA" (Met-Arg-Cys-Lys-Stop).
        // Insert "TTTCCC" at c.4_5 (mid-codon, between first and second
        // base of codon 2). Mutated CDS:
        //   "ATG" + "C" + "TTTCCC" + "GCTGCAAATAA"
        //   = "ATGCTTTCCCGCTGCAAATAA"
        //   = ATG CTT TCC CGC TGC AAA TAA
        //   = Met-Leu-Ser-Arg-Cys-Lys-Stop.
        // alt = [Met, Leu, Ser, Arg, Cys, Lys], ref = [Met, Arg, Cys, Lys].
        // At the AA level this is a 2-AA insertion + 1-AA substitution =
        // an inframe delins.
        let t = tx("ATGCGCTGCAAATAA", 1, 15);
        let seq: crate::hgvs::edit::Sequence = "TTTCCC".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
        };
        let result = predict_indel(&t, 4, 4, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        // Mid-codon in-frame edit emits delins at AA level; just assert
        // it's not a frameshift and not silently identity.
        assert!(
            !s.contains('='),
            "expected non-identity prediction, got '{}'",
            s
        );
        assert!(
            !s.contains("fs"),
            "in-frame insertion should not be frameshift, got '{}'",
            s
        );
    }

    /// Regression: a frameshift whose new in-frame stop lies in the 3'UTR
    /// (past `cds_end`) must emit `fsTer{K}` (counting the new stop's
    /// position from the first shifted residue), NOT `fsTer?`. Before the
    /// 3'UTR-extension fix, `build_frameshift_variant` scanned only the
    /// CDS-truncated `mut_cds`, so the downstream stop was never found.
    ///
    /// Construction (worked out codon-by-codon):
    ///   ref CDS  "ATGAAGAAGAAGTGA" = Met Lys Lys Lys Ter   (cds 1..=15)
    ///   3'UTR    "CTAACCC"                                  (seq 16..=22)
    ///   delete c.4 (first A of codon-2 AAG) →
    ///     mutated CDS  "ATGAGAAGAAGTGA"  (CDS-only: Met Arg Arg Ser … no stop)
    ///     + 3'UTR      "ATGAGAAGAAGTGACTAACCC"
    ///       → ATG(Met) AGA(Arg) AGA(Arg) AGT(Ser) GAC(Asp) TAA(*) — stop at
    ///         downstream position 5 (Arg=1, Arg=2, Ser=3, Asp=4, Ter=5).
    ///   first changed AA: Lys2 → Arg ⇒ p.(Lys2ArgfsTer5).
    #[test]
    fn frameshift_new_stop_in_3utr_emits_fster_count() {
        let t = tx("ATGAAGAAGAAGTGACTAACCC", 1, 15);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = predict_indel(&t, 4, 4, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Lys2ArgfsTer5)");
    }

    /// Spec (frameshift.md:22, 36-39): a frameshift whose first changed residue
    /// is immediately a stop is a NONSENSE substitution (`p.(<aa><pos>Ter)`),
    /// never `fsTer1` — `fsTer1` is an illegal count (the minimum is `fsTer2`).
    /// `c.3_4insT` on Met-Lys-Ter shifts the frame so codon 2 reads TAA: the
    /// mutated CDS "ATGTAAATAA" translates Met then an immediate stop at the
    /// first changed codon ⇒ `p.(Lys2Ter)`, not `p.(Lys2fsTer1)`.
    #[test]
    fn frameshift_immediate_stop_is_nonsense_not_fster1() {
        let t = tx("ATGAAATAA", 1, 9);
        let seq: crate::hgvs::edit::Sequence = "T".parse().unwrap();
        let edit = NaEdit::Insertion {
            sequence: crate::hgvs::edit::InsertedSequence::Literal(seq),
        };
        let result = predict_indel(&t, 3, 3, &edit, "NP_TEST.1").unwrap();
        assert_eq!(prot_str(&result), "NP_TEST.1:p.(Lys2Ter)");
    }

    // ─── Stop-loss indel extension tests (#615) ──────────────────────────────
    //
    // A small del/dup that disrupts the terminator codon is a stop-loss: the
    // protein is intact up to the last sense residue and the former stop now
    // codes an amino acid, so translation reads through the 3'UTR to the next
    // stop, yielding `p.(Ter<pos><aa>extTer<k>)` (or `extTer?` with no
    // downstream stop). Such a variant must NEVER emit the `Xaa` undetermined
    // residue nor be misrouted to the frameshift/delins path (#615).

    /// (a) Single-base deletion *inside* the terminator codon with a new
    /// in-frame stop in the 3'UTR.
    ///
    /// Transcript "ATGAAATAATTAAGGG": CDS c.1_9 = "ATGAAATAA" (Met-Lys-Ter,
    /// stop TAA at c.7_9); 3'UTR = "TTAAGGG". Delete c.7 (the first T of the
    /// stop). The CDS becomes "ATGAAAAA"; reading through the 3'UTR gives
    /// "ATGAAAAATTAAGGG" → Met-Lys-Asn-Ter, so the former stop (Ter3) now
    /// codes Asn and the very next codon is the new stop ⇒ extTer1. Crucially:
    /// no `Xaa`, and anchored at Ter3 (not the frameshift form).
    #[test]
    fn del_in_stop_codon_is_extension_with_downstream_stop() {
        let t = tx("ATGAAATAATTAAGGG", 1, 9);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = predict_indel(&t, 7, 7, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert_eq!(s, "NP_TEST.1:p.(Ter3AsnextTer1)");
        assert!(!s.contains("Xaa"), "stop-loss must not emit Xaa: '{}'", s);
    }

    /// (b) Single-base duplication that disrupts the terminator codon with a
    /// new in-frame stop in the 3'UTR.
    ///
    /// Transcript "ATGAAATAATTTAAGGG": CDS c.1_9 = "ATGAAATAA" (Met-Lys-Ter);
    /// 3'UTR = "TTTAAGGG". Duplicate c.7 (the first T of the stop). The CDS
    /// becomes "ATGAAATTAA"; reading through gives "ATGAAATTAATTTAAGGG" →
    /// Met-Lys-Leu-Ile-Ter. The former stop (Ter3) now codes Leu; the new stop
    /// is the 2nd added residue ⇒ extTer2. No `Xaa`.
    #[test]
    fn dup_disrupting_stop_codon_is_extension_with_downstream_stop() {
        let t = tx("ATGAAATAATTTAAGGG", 1, 9);
        let edit = NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: None,
        };
        let result = predict_indel(&t, 7, 7, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert_eq!(s, "NP_TEST.1:p.(Ter3LeuextTer2)");
        assert!(!s.contains("Xaa"), "stop-loss must not emit Xaa: '{}'", s);
    }

    /// (c) Stop-disrupting deletion with NO in-frame stop in the available
    /// downstream sequence → the extension length is unknown, reported as
    /// `extTer?` (extension.md:33,50). Never `Xaa`.
    ///
    /// Transcript "ATGAAATAATTTGGG": CDS c.1_9 = "ATGAAATAA"; 3'UTR = "TTTGGG"
    /// (no in-frame stop). Delete c.7 → CDS "ATGAAAAA"; read-through
    /// "ATGAAAAATTTGGG" → Met-Lys-Asn-Leu with no further stop ⇒ extTer?.
    #[test]
    fn del_in_stop_codon_without_downstream_stop_is_ext_ter_unknown() {
        let t = tx("ATGAAATAATTTGGG", 1, 9);
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let result = predict_indel(&t, 7, 7, &edit, "NP_TEST.1").unwrap();
        let s = prot_str(&result);
        assert_eq!(s, "NP_TEST.1:p.(Ter3AsnextTer?)");
        assert!(!s.contains("Xaa"), "stop-loss must not emit Xaa: '{}'", s);
    }

    // ─── Boundary-spanning (CDS→3'UTR) stop-loss extension tests (#857) ───────
    //
    // A deletion whose 5' end is in the CDS and whose 3' end is a 3'UTR (`*N`)
    // position disrupts the stop codon by removing it across the boundary. The
    // splice must run against the extended (CDS + 3'UTR) reference, not the
    // CDS-only string, and must NOT re-append the 3'UTR (no double-append).

    fn predict_stop_ext(
        t: &Transcript,
        cds_pos_start: i64,
        cds_end_utr3_n: i64,
        edit: &NaEdit,
        acc: &str,
    ) -> Result<Option<HgvsVariant>, FerroError> {
        let bundle = RefProteinBundle::from_transcript(t)?;
        predict_stop_region_extension(t, &bundle, cds_pos_start, cds_end_utr3_n, edit, acc)
    }

    /// CDS "ATGAAATAA" (Met-Lys-Ter, stop c.7_9) + 3'UTR "TCTAA". Delete c.9
    /// (last stop base 'A') + *1 ('T'). mut_full = "ATGAAATACTAA" → Met-Lys-Tyr-Ter
    /// ⇒ former stop (Ter3) codes Tyr, new stop is the 1st added residue ⇒ extTer1.
    #[test]
    fn boundary_del_into_3utr_is_extension_with_downstream_stop() {
        let t = tx("ATGAAATAATCTAA", 1, 9);
        let del = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let v = predict_stop_ext(&t, 9, 1, &del, "NP_TEST.1")
            .unwrap()
            .unwrap();
        let s = prot_str(&v);
        assert_eq!(s, "NP_TEST.1:p.(Ter3TyrextTer1)");
        assert!(!s.contains("Xaa"));
    }

    /// Same shape, 3'UTR "TCGGG" → no downstream in-frame stop ⇒ extTer?.
    #[test]
    fn boundary_del_into_3utr_without_downstream_stop_is_ext_ter_unknown() {
        let t = tx("ATGAAATAATCGGG", 1, 9);
        let del = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let v = predict_stop_ext(&t, 9, 1, &del, "NP_TEST.1")
            .unwrap()
            .unwrap();
        assert_eq!(prot_str(&v), "NP_TEST.1:p.(Ter3TyrextTer?)");
    }

    /// A boundary del that also changes a SENSE residue (5' end well inside the
    /// CDS) is not a clean stop-loss ⇒ None (prefix-equality fails).
    #[test]
    fn boundary_del_changing_sense_residue_returns_none() {
        let t = tx("ATGAAATAATCTAA", 1, 9);
        let del = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        assert!(predict_stop_ext(&t, 2, 1, &del, "NP_TEST.1")
            .unwrap()
            .is_none());
    }

    /// CDS not terminated by a stop (cds_end before the real stop) ⇒ guard ⇒ None.
    #[test]
    fn boundary_del_non_stop_terminated_cds_returns_none() {
        // CDS c.1_6 = "ATGAAA" (Met-Lys, no terminal stop).
        let t = tx("ATGAAATCTAA", 1, 6);
        let del = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        assert!(predict_stop_ext(&t, 6, 1, &del, "NP_TEST.1")
            .unwrap()
            .is_none());
    }

    /// Non-ATG (CTG) recognized-initiator start: force_initiator_met keeps the
    /// sense-prefix compare consistent so the extension still fires.
    #[test]
    fn boundary_del_non_atg_start_still_extends() {
        let t = tx("CTGAAATAATCTAA", 1, 9);
        let del = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let v = predict_stop_ext(&t, 9, 1, &del, "NP_TEST.1")
            .unwrap()
            .unwrap();
        assert_eq!(prot_str(&v), "NP_TEST.1:p.(Ter3TyrextTer1)");
    }
}
