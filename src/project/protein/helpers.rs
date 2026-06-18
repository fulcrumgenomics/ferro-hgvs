//! Pure helper functions for protein consequence prediction.
//!
//! All functions in this module are side-effect-free and easily unit-tested
//! independently of the rest of the projector machinery.

use crate::error::FerroError;
use crate::hgvs::edit::{ExtDirection, InsertedSequence, NaEdit, ProteinEdit};
use crate::hgvs::interval::ProtInterval;
use crate::hgvs::location::{AminoAcid, CdsPos, ProtPos};
use crate::hgvs::variant::{HgvsVariant, LocEdit, ProteinVariant};
use crate::project::accession::parse_accession;
use crate::reference::transcript::Transcript;
use crate::sequence::reverse_complement;

use super::substitution::{translate, translate_bytes};

/// Does `edit` modify a base of the translation initiation codon (1-based CDS
/// positions 1–3)?
///
/// An insertion places bases *between* `cds_pos_start` and the next position,
/// so it disrupts the start codon only when the gap is strictly inside it
/// (between c.1/c.2 or c.2/c.3, i.e. `cds_pos_start ∈ {1,2}`). Every other edit
/// (substitution, deletion, delins, duplication, inversion) modifies the
/// `[start,end]` span, so it disrupts the start codon whenever that span
/// overlaps `[1,3]`.
pub(crate) fn affects_initiation_codon(
    edit: &NaEdit,
    cds_pos_start: i64,
    cds_pos_end: i64,
) -> bool {
    match edit {
        NaEdit::Insertion { .. } => (1..=2).contains(&cds_pos_start),
        _ => cds_pos_start <= 3 && cds_pos_end >= 1,
    }
}

/// Does `edit` reach the translation initiation codon on the coding axis,
/// taking the full CDS-position context (offsets, intronic state) into account?
///
/// True only for an **exonic, non-offset** edit whose span overlaps CDS 1–3
/// ([`affects_initiation_codon`]). An intronic edit, or one whose start/end
/// carries an intronic offset, never reaches the start codon and returns
/// `false`. This is the exact gate the single-variant protein path uses to
/// short-circuit to the initiation-codon-unknown form (`p.(Met1?)` on an `ATG`
/// start, `p.?` on a non-AUG-initiation transcript; #771). It is exposed so the
/// allele-protein path can detect an initiation-disrupting cis member and
/// collapse the whole-allele consequence to that unknown form.
pub(crate) fn edit_reaches_initiation_codon(
    edit: &NaEdit,
    cds_start: &CdsPos,
    cds_end: &CdsPos,
    is_intronic: bool,
) -> bool {
    !is_intronic
        && cds_start.offset.is_none()
        && cds_end.offset.is_none()
        && affects_initiation_codon(edit, cds_start.base, cds_end.base)
}

/// Build a predicted `p.(Met1?)` variant for a change affecting the translation
/// initiation codon. Per HGVS (`recommendations/protein/substitution.md:51`,
/// `deletion.md:62`) the protein consequence of an initiation-codon variant
/// cannot be predicted, so it is reported as `Met1?`; emitting a concrete
/// missense (`p.Met1Thr`-style) is explicitly disallowed by the spec. Protein
/// residue 1 is the initiator `Met` regardless of the reference codon.
pub(crate) fn build_initiator_unknown(
    protein_accession: &str,
    transcript: &Transcript,
) -> HgvsVariant {
    let loc = ProtInterval::point(ProtPos::new(AminoAcid::Met, 1));
    HgvsVariant::Protein(ProteinVariant {
        accession: parse_accession(protein_accession),
        gene_symbol: transcript.gene_symbol.clone(),
        loc_edit: LocEdit::new_predicted(loc, ProteinEdit::position_unknown()),
    })
}

/// Build a whole-protein-unknown `p.?` variant.
///
/// Used for an initiation-codon-affecting edit on a **non-AUG-initiation**
/// transcript: the consequence is unknown, but the `p.(Met1?)` form
/// [`build_initiator_unknown`] emits is wrong because residue 1 is not `Met`
/// when the start codon is not `ATG`. The whole-protein `p.?` ("an effect is
/// expected but cannot be reliably predicted", `recommendations/uncertain.md`)
/// is the spec-correct form and matches mutalyzer. The location is ignored when
/// the edit is a whole-protein unknown, so a placeholder `Met1` point is fine.
pub(crate) fn build_whole_protein_unknown(
    protein_accession: &str,
    transcript: &Transcript,
) -> HgvsVariant {
    let loc = ProtInterval::point(ProtPos::new(AminoAcid::Met, 1));
    HgvsVariant::Protein(ProteinVariant {
        accession: parse_accession(protein_accession),
        gene_symbol: transcript.gene_symbol.clone(),
        loc_edit: LocEdit::new(loc, ProteinEdit::whole_protein_unknown()),
    })
}

/// Does the reference CDS begin with a canonical `ATG` start codon?
///
/// The protein-consequence path translates a transcript's CDS bases directly,
/// trusting that `Transcript.cds_start` points at the first base of the
/// initiation codon. When the cdot CDS annotation is inconsistent with the
/// transcript FASTA (cds_start lands mid-frame, not on an `ATG`), translation
/// reads the wrong frame and fabricates a garbage protein (issue #625). A
/// non-`ATG` first codon is the dominant signal of that inconsistency, so
/// callers decline protein prediction when this returns `false`.
///
/// This is a deliberate over-approximation: RefSeq does annotate a small number
/// of transcripts with experimentally-supported non-AUG initiation (CUG/GUG
/// starts carried with a `transl_except` exception), whose `cds_start` *is*
/// consistent with the FASTA. Keying on `ATG` declines protein on those too. We
/// accept that rare false-decline rather than translate the common
/// inconsistent-annotation case from the wrong frame.
///
/// The check keys on the **start codon only**. It deliberately does *not*
/// inspect internal codons: selenoproteins (e.g. SELENON / NM_020451.2)
/// legitimately carry an in-frame internal `TGA` recoded as selenocysteine,
/// and they begin with `ATG` — rejecting internal stops would wrongly decline
/// every selenoprotein.
///
/// `ref_cds` is expected uppercase (the form [`read_full_cds`] /
/// [`read_cds_start_codon`] return); the comparison is therefore a plain byte
/// compare of the first three bases. A CDS shorter than three bases returns
/// `false` (no start codon is possible), so callers correctly decline on it.
pub(crate) fn cds_has_valid_start(ref_cds: &str) -> bool {
    ref_cds.as_bytes().get(..3) == Some(b"ATG")
}

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

/// Read **only** the three start-codon bases of the CDS, uppercase.
///
/// Same coordinate validation as [`read_full_cds`] (returns
/// [`FerroError::ProteinSequenceUnavailable`] on a missing sequence or
/// degenerate coords, including a CDS shorter than three bases), but slices and
/// uppercases just the first codon instead of the whole 1–10 kbp CDS. Used by
/// the issue-#625 CDS-start sanity guard, which needs nothing more than the
/// first codon to decide whether the annotation begins with an `ATG`.
pub(crate) fn read_cds_start_codon(transcript: &Transcript) -> Result<String, FerroError> {
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
    // Reject the same degenerate coords as `read_full_cds`, plus require at
    // least three CDS bases so the start codon exists.
    if cds_start < 1 || cds_end > seq.len() || cds_start > cds_end + 1 || cds_end < cds_start + 2 {
        return Err(FerroError::ProteinSequenceUnavailable {
            accession: transcript.id.clone(),
        });
    }
    // Transcript.cds_start is 1-based → 0-based index = cds_start - 1.
    Ok(seq[cds_start - 1..cds_start + 2].to_uppercase())
}

/// Build the mutated CDS sequence by applying `edit` at a given CDS coordinate.
///
/// `cds_start` and `cds_end` are **1-based CDS positions** (as recorded in a c. variant),
/// not transcript-space positions.  For a point edit (substitution, point deletion) the
/// caller sets `cds_start == cds_end`.
///
/// The returned string is the full CDS after the edit, in uppercase.
#[cfg(test)]
pub(crate) fn build_mutated_cds(
    transcript: &Transcript,
    cds_pos_start: i64, // 1-based CDS position of edit start
    cds_pos_end: i64,   // 1-based CDS position of edit end (inclusive)
    edit: &NaEdit,
) -> Result<String, FerroError> {
    let ref_cds = read_full_cds(transcript)?;
    build_mutated_cds_with_ref(&ref_cds, transcript, cds_pos_start, cds_pos_end, edit)
}

/// Fast path of [`build_mutated_cds`] for callers that already hold the
/// reference CDS (e.g. via `RefProteinBundle::ref_cds`). Skips the per-call
/// `read_full_cds` + `to_uppercase` that would otherwise re-uppercase a
/// 1-10 kbp string on every indel projection.
///
/// The `ref_cds` argument MUST be uppercase ACGT (the form that
/// `read_full_cds` returns). Inserted / reverse-complemented slices are
/// uppercased individually before splicing so the result is still uniformly
/// uppercase without a final whole-string pass.
pub(crate) fn build_mutated_cds_with_ref(
    ref_cds: &str,
    transcript: &Transcript,
    cds_pos_start: i64,
    cds_pos_end: i64,
    edit: &NaEdit,
) -> Result<String, FerroError> {
    // Convert 1-based CDS positions to 0-based CDS-string indices.
    let idx_start = (cds_pos_start - 1) as usize;
    let idx_end = cds_pos_end as usize; // exclusive upper bound

    if idx_end > ref_cds.len() {
        return Err(FerroError::ProteinSequenceUnavailable {
            accession: transcript.id.clone(),
        });
    }

    let before = &ref_cds[..idx_start];
    let affected = &ref_cds[idx_start..idx_end];
    let after = &ref_cds[idx_end..];

    let mutated = match edit {
        NaEdit::Deletion { .. } => {
            // Remove the affected bases.
            let mut out = String::with_capacity(before.len() + after.len());
            out.push_str(before);
            out.push_str(after);
            out
        }
        NaEdit::Insertion { sequence } => {
            // HGVS insertion c.N_N+1insX: the inserted bases go between
            // cds_pos_start and cds_pos_start+1, i.e. at idx_start (after
            // `before`). Caller passes idx_end == idx_start so `affected` is
            // empty; preserved in the concat to make this explicit.
            let inserted = extract_literal_sequence(sequence, transcript)?;
            let mut out =
                String::with_capacity(before.len() + affected.len() + inserted.len() + after.len());
            out.push_str(before);
            out.push_str(affected);
            push_ascii_uppercase(&mut out, &inserted);
            out.push_str(after);
            out
        }
        NaEdit::Duplication { .. } => {
            let mut out = String::with_capacity(before.len() + 2 * affected.len() + after.len());
            out.push_str(before);
            out.push_str(affected);
            out.push_str(affected);
            out.push_str(after);
            out
        }
        NaEdit::Delins { sequence, .. } => {
            let inserted = extract_literal_sequence(sequence, transcript)?;
            let mut out = String::with_capacity(before.len() + inserted.len() + after.len());
            out.push_str(before);
            push_ascii_uppercase(&mut out, &inserted);
            out.push_str(after);
            out
        }
        NaEdit::Inversion { .. } => {
            // `affected` is already uppercase; reverse_complement preserves
            // case, so no additional uppercase pass is needed.
            let rc = reverse_complement(affected);
            let mut out = String::with_capacity(before.len() + rc.len() + after.len());
            out.push_str(before);
            out.push_str(&rc);
            out.push_str(after);
            out
        }
        NaEdit::Substitution { alternative, .. } => {
            // Single-base substitution. `alternative.to_u8()` is ASCII upper.
            let mut out = String::with_capacity(before.len() + 1 + after.len());
            out.push_str(before);
            out.push(alternative.to_u8() as char);
            out.push_str(after);
            out
        }
        _ => {
            return Err(FerroError::UnsupportedProjection {
                reason: format!("build_mutated_cds does not support edit type: {:?}", edit),
            })
        }
    };

    Ok(mutated)
}

/// Push `src` onto `dst`, uppercasing ASCII letters as we go. Avoids the
/// allocation of `src.to_uppercase()` for the common case where the caller
/// only has a small inserted slice (typically 1-100 bp).
fn push_ascii_uppercase(dst: &mut String, src: &str) {
    for &b in src.as_bytes() {
        dst.push(b.to_ascii_uppercase() as char);
    }
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
        .filter_map(translate_bytes)
        .take_while(|aa| *aa != AminoAcid::Ter)
        .collect()
}

/// Pre-translated reference CDS for a transcript.
///
/// Cached on [`VariantProjector`] so that fan-out across many indel variants
/// on the same transcript translates the reference CDS exactly once instead
/// of once per variant. The fields are exactly the values that
/// `predict_indel_protein` needs from the reference side of the comparison.
#[derive(Debug)]
pub(crate) struct RefProteinBundle {
    /// Uppercase CDS bytes. Cached so `predict_indel_protein` doesn't
    /// re-read + re-uppercase the same string for every variant on the same
    /// transcript; the alt-side CDS is built by slicing this with
    /// `build_mutated_cds_with_ref`.
    pub ref_cds: String,
    pub ref_protein: Vec<AminoAcid>,
    pub ref_protein_with_stop: Vec<AminoAcid>,
}

impl RefProteinBundle {
    /// Build a bundle by reading and translating the CDS of `transcript`.
    pub(crate) fn from_transcript(transcript: &Transcript) -> Result<Self, FerroError> {
        let ref_cds = read_full_cds(transcript)?;
        let ref_protein = translate_full_cds(&ref_cds);
        let ref_protein_with_stop = translate_full_cds_with_stop(&ref_cds);
        Ok(Self {
            ref_cds,
            ref_protein,
            ref_protein_with_stop,
        })
    }
}

/// Translate a *mutated* CDS to a protein, reusing the already-translated
/// reference prefix wherever the alt and ref CDS bytes are identical.
///
/// The first `unchanged_prefix_codons` codons of `mut_cds` (positions
/// `0..unchanged_prefix_codons * 3`) are byte-identical to `ref_cds` by
/// construction — they're the unedited prefix that `build_mutated_cds_with_ref`
/// preserved verbatim. We can lift their amino-acid translations straight out
/// of `ref_protein` and only translate the remaining tail of `mut_cds`.
///
/// Correctness notes:
/// * `ref_protein` may be shorter than `unchanged_prefix_codons` if the
///   reference's first stop falls inside the unchanged prefix. In that case
///   the alt protein also stops there (same bytes, same `take_while`), so we
///   return the ref protein verbatim.
/// * Stop-codon readthrough (an edit that removes the ref's terminator) still
///   works because we only short-circuit when the prefix already contained
///   the stop. When the edit is at or past the ref stop's codon we fall into
///   the normal translation path.
/// * Output matches `translate_full_cds(mut_cds)` exactly — verified by
///   `translate_mutated_cds_matches_full_translation` over substitution /
///   deletion / insertion / duplication / delins / inversion cases.
pub(crate) fn translate_mutated_cds(
    ref_protein: &[AminoAcid],
    mut_cds: &str,
    unchanged_prefix_codons: usize,
) -> Vec<AminoAcid> {
    let kept = unchanged_prefix_codons.min(ref_protein.len());
    let mut out: Vec<AminoAcid> = ref_protein[..kept].to_vec();
    if kept < unchanged_prefix_codons {
        // The reference already terminated inside the unchanged prefix, so
        // the alt protein terminates at the same point (same bytes).
        return out;
    }
    let start_byte = unchanged_prefix_codons * 3;
    if start_byte >= mut_cds.len() {
        return out;
    }
    out.extend(
        mut_cds.as_bytes()[start_byte..]
            .chunks_exact(3)
            .filter_map(translate_bytes)
            .take_while(|aa| *aa != AminoAcid::Ter),
    );
    out
}

/// In-frame variant of [`translate_mutated_cds`] that also lifts the
/// post-edit *suffix* from `ref_protein` whenever the mutated and reference
/// CDS bytes have re-converged.
///
/// For an in-frame indel (`net % 3 == 0`), the bytes of `mut_cds` after the
/// modification region are byte-identical to a (codon-aligned) slice of
/// `ref_cds`, so we can skip translating them and copy the corresponding
/// `ref_protein` AAs directly. Only the small edit window — bytes from the
/// first edited codon's start up to the first codon boundary at-or-past the
/// last edited byte — is translated normally.
///
/// Parameters:
/// * `cds_pos_start` / `cds_pos_end` — the 1-based CDS positions of the edit
///   (same values predict_indel_protein passes in).
/// * `net` — net base-length change (`mut_cds.len() - ref_cds.len()`). Must
///   be a multiple of 3 (the in-frame contract).
///
/// Correctness invariants:
/// * Output matches `translate_full_cds(mut_cds)` byte-for-byte; verified by
///   `translate_mutated_cds_inframe_matches_full_translation` over deletion,
///   insertion, duplication, delins, inversion, and readthrough cases.
/// * If a stop appears inside the edit window, translation terminates there
///   and the suffix is not lifted (same as `take_while(!= Ter)`).
/// * If the alt sequence extends past the ref's stop (readthrough), the
///   tail of `mut_cds` after the suffix-join point is translated normally.
pub(crate) fn translate_mutated_cds_inframe(
    ref_protein: &[AminoAcid],
    mut_cds: &str,
    cds_pos_start: i64,
    cds_pos_end: i64,
    net: i64,
) -> Vec<AminoAcid> {
    debug_assert_eq!(net % 3, 0, "in-frame fast path requires net % 3 == 0");

    let unchanged_prefix_codons = ((cds_pos_start - 1).max(0) / 3) as usize;
    let kept = unchanged_prefix_codons.min(ref_protein.len());
    let mut out: Vec<AminoAcid> = ref_protein[..kept].to_vec();
    if kept < unchanged_prefix_codons {
        // Reference already stopped inside the unchanged prefix; alt does too.
        return out;
    }

    // Locate the first mut-CDS byte at-or-past the modification (this is the
    // first byte whose value comes from the unchanged ref suffix, shifted by
    // `net` in ref-space). For all edit categories this is
    // `idx_start + mut_edit_length` where mut_edit_length = ref_edit_length + net.
    let idx_start = (cds_pos_start - 1).max(0) as usize;
    let idx_end = cds_pos_end as usize;
    let ref_edit_length = idx_end.saturating_sub(idx_start) as i64;
    let mut_edit_length = (ref_edit_length + net).max(0) as usize;
    let bx = idx_start + mut_edit_length;
    let bx_first_full_codon = bx.div_ceil(3) * 3;

    // Translate the edit window: codons starting at `start_byte` up to the
    // first codon boundary on/after the modification's end.
    let start_byte = unchanged_prefix_codons * 3;
    let edit_window_end = bx_first_full_codon.min(mut_cds.len());
    if start_byte < edit_window_end {
        for chunk in mut_cds.as_bytes()[start_byte..edit_window_end].chunks_exact(3) {
            if let Some(aa) = translate_bytes(chunk) {
                if aa == AminoAcid::Ter {
                    return out;
                }
                out.push(aa);
            }
        }
    }

    // The remaining mut suffix (bytes from `bx_first_full_codon` onward) is
    // byte-identical to ref bytes from `bx_first_full_codon - net` onward.
    // Lift the corresponding amino acids from `ref_protein`.
    let ref_byte_for_join = (bx_first_full_codon as i64 - net).max(0) as usize;
    let ref_codon_idx = ref_byte_for_join / 3;
    if ref_codon_idx <= ref_protein.len() {
        out.extend_from_slice(&ref_protein[ref_codon_idx..]);
        return out;
    }

    // Readthrough: alt has codons past ref's first stop. Translate the rest
    // of mut_cds normally (the join point is past ref_protein's end).
    if bx_first_full_codon < mut_cds.len() {
        out.extend(
            mut_cds.as_bytes()[bx_first_full_codon..]
                .chunks_exact(3)
                .filter_map(translate_bytes)
                .take_while(|aa| *aa != AminoAcid::Ter),
        );
    }
    out
}

/// Translate a CDS sequence to a protein, including the termination codon.
///
/// Stops translation immediately after the first `Ter` is encountered.
pub(crate) fn translate_full_cds_with_stop(cds: &str) -> Vec<AminoAcid> {
    let mut result = Vec::new();
    for chunk in cds.as_bytes().chunks_exact(3) {
        if let Some(aa) = translate_bytes(chunk) {
            result.push(aa);
            if aa == AminoAcid::Ter {
                break;
            }
        }
    }
    result
}

// ── C-terminal extension ──────────────────────────────────────────────────────

/// Build a C-terminal extension protein variant for a stop-loss (no-stop)
/// change.
///
/// The reference stop codon at 1-based protein position `stop_pos` now codes
/// for an amino acid, so translation reads through `mut_cds` — the mutated CDS
/// concatenated with the downstream 3'UTR, in transcript orientation — until a
/// new stop codon. The result is `p.(TerNNN<aa>extTer<k>)`, where `k` is the
/// 1-based residue offset of the new stop within the added tail, or `extTer?`
/// when no new stop is reached within the available sequence
/// (`recommendations/protein/extension.md:30,50`).
///
/// Shared by the substitution and indel stop-loss paths so both render
/// extensions identically.
pub(crate) fn build_cterminal_extension(
    stop_pos: u64,
    mut_cds: &str,
    protein_accession: &str,
    transcript: &Transcript,
) -> Result<HgvsVariant, FerroError> {
    let stop_offset = (stop_pos.saturating_sub(1) as usize) * 3;

    // The amino acid now coded by the former stop codon (None if the read-through
    // codon is itself a stop, or the sequence is too short).
    let new_aa: Option<AminoAcid> = if stop_offset + 3 <= mut_cds.len() {
        translate(&mut_cds[stop_offset..stop_offset + 3]).filter(|aa| *aa != AminoAcid::Ter)
    } else {
        None
    };

    // `extTer{K}`: K is the position of the new stop in the *added* sequence —
    // the amino acids added **beyond** the original Ter (`docs/syntax.yaml:78`,
    // extension_length). The converted-stop residue occupies the original Ter's
    // slot (the C-terminal mirror of `Met1`, the N-terminal `ext-N` anchor) and
    // is NOT counted. The downstream scan starts at the original stop codon, so
    // `downstream[0]` is that converted-stop residue and `downstream[k]` is the
    // k-th added residue; K is therefore the 0-based scan index of the new Ter
    // (do NOT add 1 — that would count the converted-stop residue). `None` when
    // the available sequence has no further stop (`extTer?`). NB: the
    // frameshift `fsTer` count differs — there the first *changed* residue is a
    // genuinely-new residue and is position 1 — so do not unify the two.
    let ext_count: Option<i64> = if stop_offset < mut_cds.len() {
        translate_full_cds_with_stop(&mut_cds[stop_offset..])
            .iter()
            .position(|aa| *aa == AminoAcid::Ter)
            .map(|p| p as i64)
    } else {
        None
    };

    let protein_edit = ProteinEdit::Extension {
        new_aa,
        direction: ExtDirection::CTerminal,
        count: ext_count,
    };
    let loc = ProtInterval::point(ProtPos::new(AminoAcid::Ter, stop_pos));
    Ok(HgvsVariant::Protein(ProteinVariant {
        accession: parse_accession(protein_accession),
        gene_symbol: transcript.gene_symbol.clone(),
        loc_edit: LocEdit::new_predicted(loc, protein_edit),
    }))
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

/// Append the transcript's unchanged 3'UTR (`seq[cds_end..]`) to an
/// already-built mutated CDS, so a shifted/read-through reading frame can find
/// a stop codon that lies past the annotated CDS end. Both the frameshift and
/// the stop-loss/extension predictors need this: without it a downstream stop
/// in the 3'UTR is missed and the prediction degrades to `fsTer?` / `extTer?`.
///
/// `mut_cds` must be uppercase. The edit always lies within the CDS, so the
/// 3'UTR slice is unchanged by it.
pub(crate) fn mut_cds_with_3utr(
    mut_cds: &str,
    transcript: &Transcript,
) -> Result<String, FerroError> {
    let seq =
        transcript
            .sequence
            .as_deref()
            .ok_or_else(|| FerroError::ProteinSequenceUnavailable {
                accession: transcript.id.clone(),
            })?;
    let cds_end = transcript
        .cds_end
        .ok_or_else(|| FerroError::ConversionError {
            msg: format!("transcript {} has no CDS end", transcript.id),
        })? as usize;
    // Guard the slice (mirrors `read_full_cds`): in the normal flow the caller
    // has already validated the transcript via `RefProteinBundle::from_transcript`,
    // but this is `pub(crate)` — return a clean error rather than panic if a
    // future caller passes a transcript whose `cds_end` exceeds its sequence.
    // A `cds_end` of 0 is equally degenerate: `seq[0..]` would append the whole
    // transcript as if it were 3'UTR, so reject it here too.
    if cds_end == 0 || cds_end > seq.len() {
        return Err(FerroError::ProteinSequenceUnavailable {
            accession: transcript.id.clone(),
        });
    }
    Ok(format!("{mut_cds}{}", seq[cds_end..].to_ascii_uppercase()))
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

/// Clamp a whole-exon deletion to the exonic CDS bases it removes, but only for
/// the **well-formed canonical bracketing** form — the one configuration whose
/// removed exonic span is unambiguous.
///
/// Returns `Some((start.base, end.base))` — the inclusive 1-based CDS span
/// removed — when ALL of these hold:
///   * `edit` is a plain `Deletion`;
///   * both endpoints are intronic (`offset.is_some()`);
///   * the interval is well-formed ascending: `start.base <= end.base`;
///   * the 5' endpoint's offset is `< 0` (it lies in the intron *before*
///     `start.base`, so `start.base` is the first exonic base removed);
///   * the 3' endpoint's offset is `> 0` (it lies in the intron *after*
///     `end.base`, so `end.base` is the last exonic base removed);
///   * the span is entirely within the CDS (`start.base > 0`).
///
/// Returns `None` for everything else.
///
/// Why only the canonical form: an intronic-flanked deletion's removed exonic
/// span is trustworthy only when the endpoints bracket the exon(s) from outside
/// (5' offset `< 0`, 3' offset `> 0`) in a well-formed ascending interval. Other
/// shapes — a descending/crossed-offset interval (e.g. ferro's reverse-strand
/// mis-normalization `c.704-18_677+65del`, #762), a 5' offset `> 0`, an
/// exonic endpoint, or a pure-/deep-intronic span — are ambiguous or known to
/// come from unreliable normalization, so predicting from them risks a *wrong*
/// protein. We decline rather than guess; cases blocked only by #762 begin
/// predicting automatically once that normalization is corrected. See the
/// design spec and HGVS `protein/deletion.md` "one or more exons".
pub(crate) fn whole_exon_deletion_span(
    edit: &NaEdit,
    cds_start: &CdsPos,
    cds_end: &CdsPos,
) -> Option<(i64, i64)> {
    if !matches!(edit, NaEdit::Deletion { .. }) {
        return None;
    }
    // Both endpoints must be intronic.
    let start_off = cds_start.offset?;
    let end_off = cds_end.offset?;
    // Reject any UTR endpoint. A 3'UTR position (`c.*N`) carries `utr3: true`
    // with a *positive* `base`, so it would otherwise pass the `base > 0`
    // check below and be mis-treated as a CDS coordinate — splicing 3'UTR
    // offsets into the CDS and emitting a wrong protein. (5'UTR is `base <= 0`,
    // already rejected below, but guard it here too for clarity.)
    if cds_start.utr3 || cds_end.utr3 {
        return None;
    }
    // Canonical bracketing only: well-formed ascending interval whose 5' end
    // sits in the intron before its base and 3' end in the intron after its
    // base. Anything else (incl. the reverse-strand crossed form) is declined.
    if cds_start.base > cds_end.base || start_off >= 0 || end_off <= 0 {
        return None;
    }
    let lo = cds_start.base;
    let hi = cds_end.base;
    // Span must sit inside the CDS, not the 5'UTR (`base <= 0`). The ascending
    // `base` order (`lo <= hi`) is already guaranteed by the `cds_start.base >
    // cds_end.base` rejection above.
    if lo <= 0 {
        return None;
    }
    Some((lo, hi))
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
            protein_id: None,
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

    /// Locks the per-codon behaviour of `translate_full_cds` end-to-end across
    /// every ACGT codon. Concatenating all 64 codons gives a 192-base CDS;
    /// adding three trailing Ns checks that an invalid codon at the end is
    /// silently skipped (filter_map convention) rather than aborting.
    #[test]
    fn translate_full_cds_all_64_codons() {
        // Build a 192-base CDS that contains every codon exactly once, in a
        // fixed canonical order. None of the codons before the first stop
        // matter — the take_while-before-Ter contract means we'd cut off at
        // the first stop. So instead we use a stop-free synthetic CDS by
        // omitting TAA/TAG/TGA from the body and asserting via the
        // `_with_stop` variant for the stop-codon path.
        let codons = [
            "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA", "ATG",
            "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "CCT", "CCC",
            "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG", "TAT", "TAC",
            "CAT", "CAC", "CAA", "CAG", "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
            "TGT", "TGC", "TGG", "CGT", "CGC", "CGA", "CGG", "AGA", "AGG", "GGT", "GGC", "GGA",
            "GGG",
        ];
        // 61 sense codons (no stops). Concat + 2 trailing bases (incomplete).
        let cds: String = codons.iter().copied().collect::<String>() + "NN";
        let aas = translate_full_cds(&cds);
        assert_eq!(
            aas.len(),
            61,
            "61 sense codons translated, incomplete tail dropped"
        );
        // First codon TTT → Phe, last codon GGG → Gly.
        assert_eq!(aas[0], AminoAcid::Phe);
        assert_eq!(aas[60], AminoAcid::Gly);

        // _with_stop variant: prepend a stop right after the first codon and
        // confirm translation halts immediately after Ter.
        let cds_with_stop = String::from("ATGTAA") + &cds;
        let aas = translate_full_cds_with_stop(&cds_with_stop);
        assert_eq!(aas, vec![AminoAcid::Met, AminoAcid::Ter]);
    }

    /// `translate_mutated_cds` must produce byte-identical results to
    /// `translate_full_cds(mut_cds)` for every prefix offset between 0 and
    /// the codon count, across every indel category exercised below. This is
    /// the safety net for the localised-translation perf change: any
    /// off-by-one in the prefix-keep arithmetic falls out here.
    #[test]
    fn translate_mutated_cds_matches_full_translation() {
        // Reference CDS: 6 sense codons + stop = 21 bases.
        // Met-Arg-Pro-Thr-Ala-Cys-Ter
        let ref_cds = "ATGCGGCCCACCGCGTGCTAA";
        let ref_protein = translate_full_cds(ref_cds);
        assert_eq!(ref_protein.len(), 6, "Met..Cys (stops before Ter)");

        // Every (mut_cds, unchanged_prefix_codons) pair below was constructed
        // so that the first `unchanged_prefix_codons * 3` bytes of mut_cds
        // are identical to `ref_cds[..unchanged_prefix_codons * 3]` — i.e.
        // it preserves the precondition of `translate_mutated_cds`.
        let cases: &[(&str, &str, usize)] = &[
            // edit at codon 0 — no prefix kept; substitution of first base.
            ("subst codon 0", "GTGCGGCCCACCGCGTGCTAA", 0),
            // edit at codon 1 — 1 codon kept; mid-codon deletion.
            ("del codon 1", "ATGCGGCCACCGCGTGCTAA", 1),
            // edit at codon 2 — 2 codons kept; codon-aligned 3-base deletion.
            ("del codon 2 whole", "ATGCGGACCGCGTGCTAA", 2),
            // edit at codon 3 — 3 codons kept; in-frame insertion of 3 bases.
            ("ins after codon 3", "ATGCGGCCCACCGGGGCGTGCTAA", 3),
            // edit at codon 4 — 4 codons kept; frameshift insertion (+1 base).
            ("fs ins codon 4", "ATGCGGCCCACCGCGCTGCTAA", 4),
            // edit at codon 5 — 5 codons kept; substitution that introduces a
            // new stop one codon earlier than the ref's stop.
            ("nonsense codon 5", "ATGCGGCCCACCGCGTAATGCTAA", 5),
            // edit at codon 6 (past first stop) — 6 codons kept; alt also
            // stopped at codon 6 in the ref translation, so kept >= ref_len.
            ("readthrough", "ATGCGGCCCACCGCGTGCCCATAA", 6),
        ];

        for (name, mut_cds, prefix) in cases {
            let from_full = translate_full_cds(mut_cds);
            let from_partial = translate_mutated_cds(&ref_protein, mut_cds, *prefix);
            assert_eq!(
                from_full, from_partial,
                "case '{}': translate_mutated_cds (prefix={}) must match \
                 translate_full_cds(mut_cds); full={:?}, partial={:?}",
                name, prefix, from_full, from_partial
            );
        }
    }

    /// Byte-for-byte equivalence between `translate_mutated_cds_inframe`
    /// (the prefix + suffix lift) and `translate_full_cds` (the unchanged
    /// reference implementation) across every in-frame indel category. Any
    /// off-by-one in the suffix-join arithmetic falls out here.
    #[test]
    fn translate_mutated_cds_inframe_matches_full_translation() {
        // Reference CDS: 6 sense codons + stop = 21 bases.
        // Met-Arg-Pro-Thr-Ala-Cys-Ter
        let ref_cds = "ATGCGGCCCACCGCGTGCTAA";
        let ref_protein = translate_full_cds(ref_cds);
        assert_eq!(ref_protein.len(), 6);

        // (name, mut_cds, cds_pos_start, cds_pos_end, net)
        // Each mut_cds was constructed by applying the named edit to ref_cds
        // so that its prefix/suffix invariants match `build_mutated_cds_with_ref`.
        let cases: &[(&str, &str, i64, i64, i64)] = &[
            // Codon-aligned in-frame deletion of CCC at c.7-9.
            ("del whole codon 3", "ATGCGGACCGCGTGCTAA", 7, 9, -3),
            // 6-base codon-aligned deletion at c.7-12.
            ("del two whole codons", "ATGCGGGCGTGCTAA", 7, 12, -6),
            // Mid-codon in-frame deletion (6 bases straddling codon boundaries
            // c.5-10 — deletes GGCCCA, leaves ATGC + CCGCGTGCTAA = ATGCCCGCGTGCTAA).
            ("mid-codon 6-base del", "ATGCCCGCGTGCTAA", 5, 10, -6),
            // Codon-aligned in-frame insertion of GGG after c.3 (between Met
            // and Arg codons).
            (
                "codon-aligned ins after c.3",
                "ATGGGGCGGCCCACCGCGTGCTAA",
                3,
                3,
                3,
            ),
            // Codon-aligned in-frame insertion that also introduces a stop.
            ("ins introducing stop", "ATGTAACGGCCCACCGCGTGCTAA", 3, 3, 3),
            // In-frame delins: replace TGC (Cys) with TAATGC (Ter + Cys) at
            // c.16-18. Length: +3. Introduces a premature stop.
            (
                "delins introducing stop",
                "ATGCGGCCCACCGCGTAATGCTAA",
                16,
                18,
                3,
            ),
            // Inversion of codon 2 (CCC → GGG via revcomp).
            ("inversion of codon 2", "ATGCGGGGGACCGCGTGCTAA", 7, 9, 0),
            // Same-length delins (net = 0) replacing one codon with another.
            ("same-length delins", "ATGCGGAAAACCGCGTGCTAA", 7, 9, 0),
            // Readthrough: delete the entire stop codon, alt extends.
            // mut = ref minus TAA + 3 new bases after CDS, but for simplicity
            // we model a same-length stop-substituting delins (TAA → CCG +
            // 3 new bases that include a downstream stop).
            (
                "readthrough w/ downstream stop",
                "ATGCGGCCCACCGCGTGCCCGCCATAA",
                19,
                21,
                6,
            ),
        ];

        for (name, mut_cds, start, end, net) in cases {
            let from_full = translate_full_cds(mut_cds);
            let from_inframe =
                translate_mutated_cds_inframe(&ref_protein, mut_cds, *start, *end, *net);
            assert_eq!(
                from_full, from_inframe,
                "case '{}' (start={}, end={}, net={}): inframe path must \
                 match full translate; full={:?}, inframe={:?}",
                name, start, end, net, from_full, from_inframe
            );
        }
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

    // ── whole_exon_deletion_span ──────────────────────────────────────────

    #[test]
    fn whole_exon_span_canonical_bracketing() {
        // c.677-18_704+65del: 5' offset <0, 3' offset >0 → clamp [677, 704].
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let start = CdsPos::with_offset(677, -18);
        let end = CdsPos::with_offset(704, 65);
        assert_eq!(
            whole_exon_deletion_span(&edit, &start, &end),
            Some((677, 704))
        );
    }

    #[test]
    fn whole_exon_span_inframe_codon_aligned() {
        // DMD c.961-1_1149+3del → [961, 1149] (189 bases).
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let start = CdsPos::with_offset(961, -1);
        let end = CdsPos::with_offset(1149, 3);
        assert_eq!(
            whole_exon_deletion_span(&edit, &start, &end),
            Some((961, 1149))
        );
    }

    #[test]
    fn whole_exon_span_noncanonical_offset_positive_start_declines() {
        // c.960+1_1149+3del: 5' offset > 0 is non-canonical. We decline rather
        // than guess (an offset>0 start can't be safely disambiguated from
        // unreliable normalization). → None.
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let start = CdsPos::with_offset(960, 1);
        let end = CdsPos::with_offset(1149, 3);
        assert_eq!(whole_exon_deletion_span(&edit, &start, &end), None);
    }

    #[test]
    fn whole_exon_span_pure_intron_declines() {
        // c.681+1_682-1del: lo = 682, hi = 681 → lo > hi → None.
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let start = CdsPos::with_offset(681, 1);
        let end = CdsPos::with_offset(682, -1);
        assert_eq!(whole_exon_deletion_span(&edit, &start, &end), None);
    }

    #[test]
    fn whole_exon_span_deep_intron_declines() {
        // c.961-50_961-10del: both in the same intron, no exonic base → None.
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let start = CdsPos::with_offset(961, -50);
        let end = CdsPos::with_offset(961, -10);
        assert_eq!(whole_exon_deletion_span(&edit, &start, &end), None);
    }

    #[test]
    fn whole_exon_span_requires_both_intronic() {
        // One exonic endpoint (partial-exon) → None.
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let start = CdsPos::with_offset(632, -5);
        let end = CdsPos::new(670); // exonic, offset None
        assert_eq!(whole_exon_deletion_span(&edit, &start, &end), None);
    }

    #[test]
    fn whole_exon_span_only_plain_deletions() {
        // A substitution with offsets is never a whole-exon deletion.
        let edit = NaEdit::Substitution {
            reference: crate::hgvs::edit::Base::A,
            alternative: crate::hgvs::edit::Base::G,
        };
        let start = CdsPos::with_offset(677, -18);
        let end = CdsPos::with_offset(704, 65);
        assert_eq!(whole_exon_deletion_span(&edit, &start, &end), None);
    }

    #[test]
    fn whole_exon_span_reverse_strand_crossed_form_declines() {
        // The exact #762 reverse-strand mis-normalization for VSIR:
        // cds_start = (704, -18), cds_end = (677, +65). Descending base order
        // with crossed offsets — declining here is what prevents emitting a
        // wrong protein (e.g. p.(Ala227GlnfsTer6) instead of the correct
        // p.(Arg226ProfsTer102)). Predicts correctly once #762 is fixed and
        // the endpoints arrive canonical as (677, -18)/(704, +65).
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let start = CdsPos::with_offset(704, -18);
        let end = CdsPos::with_offset(677, 65);
        assert_eq!(whole_exon_deletion_span(&edit, &start, &end), None);
    }

    #[test]
    fn whole_exon_span_utr_only_declines() {
        // Exonic overlap entirely in the 5'UTR (CDS base <= 0) → None.
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let start = CdsPos::with_offset(-5, -2);
        let end = CdsPos::with_offset(-1, 2);
        assert_eq!(whole_exon_deletion_span(&edit, &start, &end), None);
    }

    #[test]
    fn whole_exon_span_utr3_intronic_declines() {
        // 3'UTR positions (`c.*N`) carry utr3=true with a *positive* base, so
        // they'd pass the `base > 0` check; the explicit utr3 guard must reject
        // them rather than splice 3'UTR offsets into the CDS. (c.*100-5_*200+5del)
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let start = CdsPos {
            utr3: true,
            ..CdsPos::with_offset(100, -5)
        };
        let end = CdsPos {
            utr3: true,
            ..CdsPos::with_offset(200, 5)
        };
        assert_eq!(whole_exon_deletion_span(&edit, &start, &end), None);
    }
}
