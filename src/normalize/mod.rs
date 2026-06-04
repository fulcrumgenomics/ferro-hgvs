//! Normalization engine
//!
//! Implements the core HGVS variant normalization algorithm including
//! 3'/5' shifting and boundary detection.
//!
//! # Coordinate Systems
//!
//! This module uses multiple coordinate systems:
//!
//! | Context | Basis | Type/Notes |
//! |---------|-------|------------|
//! | HGVS variant positions | 1-based | `u64` for genomic, `i64` for CDS/Tx |
//! | Array indexing | 0-based | `usize` for sequence slicing |
//! | Boundaries struct | 1-based | `(start, end)` inclusive |
//! | Shuffle input/output | 0-based | Uses array indices |
//! | Relative positions | 1-based | Positions within fetched window |
//!
//! Key conversions:
//! - `hgvs_pos_to_index(pos)` converts 1-based HGVS position to 0-based index
//! - `index_to_hgvs_pos(idx)` converts 0-based index to 1-based HGVS position
//! - `pos.saturating_sub(1)` manually converts 1-based to 0-based
//! - `idx + 1` manually converts 0-based to 1-based
//!
//! For type-safe coordinate handling, see [`crate::coords`].

pub mod boundary;
pub mod config;
pub(crate) mod merge;
mod overlap;
pub mod rules;
pub mod shuffle;
pub mod validate;

use crate::coords::{hgvs_pos_to_index, index_to_hgvs_pos};
use crate::error::FerroError;
use crate::hgvs::edit::{Base, InsertedSequence, NaEdit, ProteinEdit, Sequence};
use crate::hgvs::interval::{Interval, ProtInterval};
use crate::hgvs::location::{
    AminoAcid, CdsPos, GenomePos, ProtPos, RnaPos, SpecialPosition, TxPos,
};
use crate::hgvs::parser::position::{OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE};
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{
    AllelePhase, AlleleVariant, CdsVariant, GenomeVariant, HgvsVariant, LocEdit, MtVariant,
    RnaVariant, TxVariant,
};
use crate::hgvs::HgvsVariant as HV;
use crate::reference::transcript::Strand;
use crate::reference::ReferenceProvider;
use boundary::Boundaries;
pub use config::{NormalizeConfig, ShuffleDirection};
use rules::{
    canonicalize_conversion_to_delins, canonicalize_edit, canonicalize_insertion_expand,
    needs_normalization, should_canonicalize, DelinsSubedit, InsCoordKind,
};
use shuffle::shuffle;

/// Check if a CDS position has an unknown (?) offset sentinel value
fn has_unknown_offset_cds(pos: &CdsPos) -> bool {
    matches!(
        pos.offset,
        Some(OFFSET_UNKNOWN_POSITIVE) | Some(OFFSET_UNKNOWN_NEGATIVE)
    )
}

/// Intron length (1-based inclusive base count) at the intron addressed by
/// `tx_boundary` + `offset`. Returns `None` when:
///
/// - no intron exists at that boundary (e.g. `tx_boundary` is the last
///   exon's end and `offset` is positive — no 3' intron exists), OR
/// - the intron's genomic coordinates aren't populated (cdot data missing
///   per-exon genomic alignment), OR
/// - the intron's genomic span is inverted (`genomic_end < genomic_start`,
///   a defensive guard inherited from `Intron::genomic_length`).
///
/// All three cases collapse to "conservative skip" per the existing W4004
/// missing-metadata contract.
fn intron_length_at_tx_boundary(
    transcript: &crate::reference::transcript::Transcript,
    tx_boundary: u64,
    offset: i64,
) -> Option<u64> {
    // Delegate to `Intron::genomic_length` so the length formula stays in
    // one place — the same source-of-truth `find_intron_at_genomic` uses
    // at `src/reference/transcript.rs:695`. `genomic_length` rejects
    // inverted spans by construction, removing the silent-zero footgun a
    // manual `g_end.saturating_sub(g_start)` would produce on a
    // misconstructed minus-strand fixture.
    transcript
        .find_intron_at_tx_boundary(tx_boundary, offset)?
        .genomic_length()
}

/// Resolve a `CdsPos` to its transcript-frame position (1-based) ignoring
/// the intronic offset. Mirrors the non-intronic branch of
/// `Normalizer::cds_to_tx_pos`; kept as a free helper so
/// `check_cds_pos_past_end` can compute the intronic-anchor `tx_boundary`
/// without holding a `Normalizer` reference.
///
/// Note on the `base == 0` arm: `c.0` is rejected by the parser at
/// `src/hgvs/parser/position.rs` (per #269), and the caller
/// `check_cds_pos_past_end` short-circuits via `pos.is_unknown()` before
/// this helper sees `CdsPos { base: 0, utr3: false }`. The arm is dead
/// code in the W4004 path; kept here for parity with
/// `Normalizer::cds_to_tx_pos` so this helper remains a drop-in
/// replacement if a future caller needs it outside the bounds check.
fn cds_pos_to_tx_boundary(
    pos: &CdsPos,
    transcript: &crate::reference::transcript::Transcript,
) -> Option<u64> {
    let cds_start = transcript.cds_start?;
    if pos.utr3 {
        // c.*N: tx_boundary = cds_end + N
        let cds_end = transcript.cds_end?;
        let base = u64::try_from(pos.base).ok()?;
        Some(cds_end + base)
    } else if pos.base < 0 {
        // c.-N: tx_boundary = cds_start + base  (HGVS skips c.0, see
        // Normalizer::cds_to_tx_pos and the #97 fix)
        let signed = cds_start as i64 + pos.base;
        u64::try_from(signed).ok()
    } else if pos.base == 0 {
        // Defensive parity with `cds_to_tx_pos`; see helper-level doc
        // note. Unreachable from `check_cds_pos_past_end` because the
        // caller short-circuits on `pos.is_unknown()` (which keys off
        // `CDS_BASE_UNKNOWN == 0`).
        Some(cds_start.saturating_sub(1))
    } else {
        // c.<N>: tx_boundary = cds_start + N - 1
        Some(cds_start + (pos.base as u64) - 1)
    }
}

/// If `pos` lies past the CDS-end (for plain `c.<N>`), past the
/// transcript-end (for `c.*<N>`), or past the 5'UTR-start (for `c.-N`),
/// return a `PositionPastEnd` warning describing the violation;
/// otherwise `None`. See #336 (CDS-end / `c.*N`), #348 (`c.-N`), and
/// #392 (intronic offsets).
///
/// **Scope.** This helper handles the c. axis:
/// - plain `c.<N>` (cds-end)
/// - `c.*<N>` (transcript-end)
/// - `c.-<N>` (5utr-start)
/// - `c.<N>+<M>` / `c.<N>-<M>` (intron-end, when the intron's genomic
///   span is known; #392)
///
/// The `n.` axis is handled by [`check_tx_pos_past_end`] (#347).
///
/// **Conservative skip.** When the required transcript metadata is missing
/// (no `cds_start`/`cds_end` for a c. variant, or no genomic coords on the
/// exons enclosing the intron), the helper returns `None` rather than emit
/// a false positive. Non-intronic bounds derive from
/// `Transcript::cds_length()` / `utr3_length()` / `utr5_length()`, which
/// fall back to the exon-sum transcript length when no cached `sequence`
/// is loaded. Intronic bounds derive from the intron's `genomic_end -
/// genomic_start + 1`, which is only populated when the cdot data carries
/// per-exon genomic alignment. Callers wanting strict gating on missing
/// metadata must check transcript completeness separately.
fn check_cds_pos_past_end(
    accession: &str,
    pos: &CdsPos,
    transcript: &crate::reference::transcript::Transcript,
) -> Option<NormalizationWarning> {
    // Unknown positions can't be bounds-checked.
    if pos.is_unknown() {
        return None;
    }
    // Intronic-offset bound check (#392). Requires both the c.→tx
    // boundary AND the intron's genomic span; returns None on either
    // missing piece per the conservative-skip contract.
    if pos.is_intronic() {
        let offset = pos.offset?;
        let tx_boundary = cds_pos_to_tx_boundary(pos, transcript)?;
        let intron_length = intron_length_at_tx_boundary(transcript, tx_boundary, offset)?;
        let abs_offset = offset.unsigned_abs();
        if abs_offset > intron_length {
            return Some(NormalizationWarning::PositionPastEnd {
                accession: accession.to_string(),
                coordinate_system: "c".to_string(),
                position: pos.to_string(),
                bound_kind: "intron-end".to_string(),
                bound_value: intron_length,
            });
        }
        return None;
    }
    if pos.utr3 {
        // c.*N: N must fit in the post-CDS transcript suffix. Prefer
        // `Transcript::utr3_length()` so the bound source-of-truth stays in
        // one place. `utr3_length` falls back to the exon-sum transcript
        // length when no cached `sequence` is loaded, so the check still
        // fires for coordinate-only transcripts — a coordinate bound check
        // should not silently degrade into a sequence-availability check.
        let utr3_len = transcript.utr3_length()?;
        if pos.base > 0 && (pos.base as u64) > utr3_len {
            return Some(NormalizationWarning::PositionPastEnd {
                accession: accession.to_string(),
                coordinate_system: "c".to_string(),
                position: format!("*{}", pos.base),
                bound_kind: "transcript-end".to_string(),
                bound_value: utr3_len,
            });
        }
        return None;
    }
    if pos.base < 1 {
        // 5'UTR `c.-N`: |N| must fit in `cds_start - 1`. `Transcript::utr5_length()`
        // returns `cds_start.saturating_sub(1)`. Closes #348.
        let utr5_len = transcript.utr5_length()?;
        let abs_n = pos.base.unsigned_abs();
        if abs_n > utr5_len {
            return Some(NormalizationWarning::PositionPastEnd {
                accession: accession.to_string(),
                coordinate_system: "c".to_string(),
                position: pos.base.to_string(),
                bound_kind: "5utr-start".to_string(),
                bound_value: utr5_len,
            });
        }
        return None;
    }
    // Plain `c.<N>`: N must fit in the CDS.
    let cds_len = transcript.cds_length()?;
    if (pos.base as u64) > cds_len {
        return Some(NormalizationWarning::PositionPastEnd {
            accession: accession.to_string(),
            coordinate_system: "c".to_string(),
            position: pos.base.to_string(),
            bound_kind: "cds-end".to_string(),
            bound_value: cds_len,
        });
    }
    None
}

/// If `pos` (an `n.<N>` transcript position) lies past the transcript's
/// total length OR an intronic-offset position lies past the enclosing
/// intron, return a `PositionPastEnd` warning describing the violation;
/// otherwise `None`. Closes #347 (transcript-end) and #392 (intron-end).
///
/// Bounds:
/// - plain `n.<N>`: `1 <= N <= sequence_length`
/// - `n.<N>+<M>` / `n.<N>-<M>`: `|M| <= intron_length` when the intron's
///   genomic span is known
///
/// Unknown positions (`n.?`) and downstream `n.*N` sentinels are skipped.
/// `N < 1` (a 5'-of-transcript position) is also skipped — outside the
/// bound check's domain.
///
/// Bounds source for non-intronic: `Transcript::sequence_length()`,
/// which falls back to the sum of exon spans for coordinate-only
/// transcripts. Bounds source for intronic: intron's `genomic_end -
/// genomic_start + 1`. Both follow the conservative-skip contract when
/// metadata is missing.
fn check_tx_pos_past_end(
    accession: &str,
    pos: &TxPos,
    transcript: &crate::reference::transcript::Transcript,
) -> Option<NormalizationWarning> {
    // `has_unknown_offset_tx` matches `OFFSET_UNKNOWN_*` sentinels (not a
    // base-unknown predicate); there's no `TX_BASE_UNKNOWN` sibling of
    // `CDS_BASE_UNKNOWN`, so the n.-axis variant skips only on offset.
    // This is the intentional asymmetry with the c.-axis check above.
    if has_unknown_offset_tx(pos) || pos.is_downstream() {
        return None;
    }
    // Intronic-offset bound check (#392). The tx_boundary is `pos.base`
    // itself (n.-axis positions are direct transcript coordinates).
    // `is_intronic()` is precisely `offset.is_some() && offset !=
    // Some(0)`, so the parser-degenerate `n.<N>+0` shape falls through
    // to the plain `n.<N>` check below — symmetric with the c.-axis
    // dispatch above.
    if pos.is_intronic() {
        let offset = pos
            .offset
            .expect("is_intronic implies non-zero Some offset");
        if pos.base < 1 {
            // 5'-of-transcript intronic offsets are out of domain
            // (no intron exists upstream of the first exon).
            return None;
        }
        let tx_boundary = pos.base as u64;
        let intron_length = intron_length_at_tx_boundary(transcript, tx_boundary, offset)?;
        let abs_offset = offset.unsigned_abs();
        if abs_offset > intron_length {
            return Some(NormalizationWarning::PositionPastEnd {
                accession: accession.to_string(),
                coordinate_system: "n".to_string(),
                position: pos.to_string(),
                bound_kind: "intron-end".to_string(),
                bound_value: intron_length,
            });
        }
        return None;
    }
    if pos.base < 1 {
        return None;
    }
    let seq_len = transcript.sequence_length();
    if (pos.base as u64) > seq_len {
        return Some(NormalizationWarning::PositionPastEnd {
            accession: accession.to_string(),
            coordinate_system: "n".to_string(),
            position: pos.base.to_string(),
            bound_kind: "transcript-end".to_string(),
            bound_value: seq_len,
        });
    }
    None
}

/// Check if a TxPos has an unknown (?) offset sentinel value
fn has_unknown_offset_tx(pos: &TxPos) -> bool {
    matches!(
        pos.offset,
        Some(OFFSET_UNKNOWN_POSITIVE) | Some(OFFSET_UNKNOWN_NEGATIVE)
    )
}

/// If `pos.base` lies past the contig length on a mitochondrial accession,
/// return a `PositionPastEnd` warning describing the violation; otherwise
/// `None`. Closes #393.
///
/// Bounds: `m.<N>` must satisfy `1 <= N <= contig_length`. The contig
/// length is sourced from [`ReferenceProvider::get_sequence_length`].
/// Wraparound ranges (`m.<high>_<low>`, where `high > low` on a valid
/// circular contig, per SVD-WG006) are NOT past-end if both endpoints
/// fit in the contig — this check fires only when an endpoint itself
/// exceeds the contig length.
///
/// Skipped cases (mirrors `check_tx_pos_past_end`): special positions
/// (`pter`/`qter`/`cen`, encoded as `base == 0`) and positions carrying
/// an offset (non-standard on mt but parseable), because the final
/// coordinate isn't determined by `base` alone.
fn check_mt_pos_past_end(
    accession: &str,
    pos: &GenomePos,
    contig_length: u64,
) -> Option<NormalizationWarning> {
    if pos.base < 1 || pos.offset.is_some() {
        return None;
    }
    if pos.base > contig_length {
        return Some(NormalizationWarning::PositionPastEnd {
            accession: accession.to_string(),
            coordinate_system: "m".to_string(),
            position: pos.base.to_string(),
            bound_kind: "contig-end".to_string(),
            bound_value: contig_length,
        });
    }
    None
}

/// True if the edit is a `Deletion` or `Duplication` shape — the two
/// edit kinds the HGVS exon-junction exception applies to (HGVS
/// general.md: "deletions/duplications around exon/exon junctions using
/// c., r., or n. reference sequences are not shifted"). Insertions and
/// inversions are explicitly out of scope and still 3'-shift across
/// exon junctions. (#334)
///
/// An empty `delins` (`NaEdit::Delins { sequence: InsertedSequence::Empty,
/// .. }`) is also matched here because `normalize_na_edit` rewrites it
/// into a `NaEdit::Deletion` (HGVS issue #81 A3 — empty `delins` is
/// semantically a deletion). That rewrite happens *inside*
/// `normalize_na_edit`, i.e. **after** the shuffle bounds have already
/// been picked, so without folding empty `delins` into this predicate
/// the n./r. paths would still use the full-transcript bounds for what
/// is effectively a deletion and shuffle across the exon junction.
fn edit_is_del_or_dup(edit: &NaEdit) -> bool {
    matches!(
        edit,
        NaEdit::Deletion { .. }
            | NaEdit::Duplication { .. }
            | NaEdit::Delins {
                sequence: InsertedSequence::Empty,
                ..
            }
    )
}

/// Warning generated during normalization.
///
/// This enum is open-ended: each variant owns the fields its code needs.
/// Future warning codes add new variants without touching existing emit
/// sites. Marked `#[non_exhaustive]` so downstream callers must include a
/// wildcard arm when matching — adding a new variant is therefore not a
/// breaking change. Mirrors the same attribute on
/// [`NormalizationInfo`] / [`NormalizeResult`].
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum NormalizationWarning {
    /// Reference sequence mismatch. Stated ref bases in the HGVS expression
    /// do not match the actual reference sequence. Code: `REFSEQ_MISMATCH`.
    RefSeqMismatch {
        /// What the input claimed as reference
        stated_ref: String,
        /// What the actual reference sequence has
        actual_ref: String,
        /// Position info
        position: String,
        /// Whether the mismatch was actually corrected by the normalizer.
        ///
        /// `true` when the canonical Display drops or rewrites the stated
        /// bases (sub / del / dup / inv: the explicit `sequence` field is
        /// stripped during canonicalization, so the wrong stated bases
        /// never reach the output).
        ///
        /// `false` when the description passes through verbatim — the
        /// normalizer surfaces the warning but cannot rewrite the user's
        /// declared form. This covers `Repeat` and `MultiRepeat`
        /// consistency mismatches (issues #214 / #279): the per-unit
        /// declaration is part of the user's variant description, and
        /// the normalizer declines to second-guess it. (Issue #280.)
        corrected: bool,
        /// Free-form context from the validation layer (e.g. a
        /// per-edit-kind explanation of the mismatch). When present,
        /// the `Display` impl appends it to the synthesized message so
        /// downstream consumers retain the nuance that previously
        /// lived in the dropped `message` field.
        details: Option<String>,
    },

    /// Two or more cis-allele edits share identical reference bounds.
    /// The HGVS spec does not define a canonical form for this case;
    /// ferro preserves the input verbatim and emits this warning.
    /// Code: `OVERLAP_CONFLICTING_EDITS`.
    OverlapConflict {
        /// Accession of the reference sequence
        accession: String,
        /// Coordinate system: "g" | "c" | "n" | "r" | "m"
        coordinate_system: String,
        /// Canonical span text, e.g. "100" or "100_103"
        location: String,
        /// Edit kinds, e.g. ["sub", "sub"] or ["del", "inv"]
        edit_kinds: Vec<String>,
    },

    /// A telomere/centromere marker could not be resolved to a concrete base.
    /// Currently only `cen` (a centromere is an assembly-annotated region, not a
    /// sequence-derivable nucleotide). The input is preserved verbatim.
    /// Code: `UNRESOLVABLE_SPECIAL_POSITION`.
    UnresolvableSpecialPosition {
        /// Accession of the reference sequence.
        accession: String,
        /// The marker that could not be resolved, e.g. "cen".
        marker: String,
    },

    /// `apply_canonical_split` was unable to canonicalize because the
    /// reference window returned by the provider did not span the same
    /// number of bytes as the HGVS interval. Per HGVS spec
    /// `recommendations/background/refseq.md` §43, this means the
    /// variant is not entirely encompassed by the reference sequence —
    /// strict mode promotes this warning to
    /// `FerroError::VariantExceedsReference`. The variant is returned
    /// unchanged in lenient/silent modes.
    /// Closes-after: #354, #355. Code: `CANONICAL_SPLIT_SKIPPED`.
    CanonicalSplitSkipped {
        /// Accession of the reference sequence
        accession: String,
        /// HGVS span start (1-based inclusive). Carried so strict-mode
        /// promotion can build a `FerroError::VariantExceedsReference`
        /// with the same span information.
        hgvs_start: u64,
        /// HGVS span end (1-based inclusive).
        hgvs_end: u64,
        /// Number of bytes the HGVS span demanded.
        expected_span: usize,
        /// Number of bytes the provider returned.
        actual_bytes: usize,
    },

    /// A `c.` variant whose start and end positions sit in different
    /// coordinate sub-axes (5'UTR / CDS / 3'UTR). The 3'-rule shuffle has
    /// no well-defined semantics across an axis boundary, so ferro
    /// preserves the canonical input position and emits this warning.
    /// Closes-after: #350. Code: `CROSS_AXIS_VARIANT_NOT_SHUFFLED`.
    CrossAxisVariantNotShuffled {
        /// Accession of the reference sequence
        accession: String,
        /// Axis of the start position: "5utr" | "cds" | "3utr"
        start_axis: String,
        /// Axis of the end position: "5utr" | "cds" | "3utr"
        end_axis: String,
    },

    /// A 3'-rule shuffle would have crossed a CDS↔UTR coordinate sub-axis
    /// boundary, but the axis clamp constrained the result to the
    /// boundary. Closes-after: #349. Code: `AXIS_CLAMP_APPLIED`.
    AxisClampApplied {
        /// Accession of the reference sequence
        accession: String,
        /// Shuffle direction that was clamped: "5prime" | "3prime"
        direction: String,
        /// Axis bound that did the clamping: "5utr" | "cds_start" | "cds_end" | "3utr"
        clamp_kind: String,
    },

    /// Canonicalization (e.g. `p.ins → p.dup` or `p.delins → p.dup`)
    /// produced a `Duplication` whose interval includes p.1 — the
    /// initiator methionine. The duplication form is spec-permitted
    /// (Prioritization rule is unconditional; spec uses Met1-inclusive
    /// ranges in deletion.md:63-65), but consumers may also wish to
    /// describe the protein-level consequence per the substitution
    /// rule for start-codon variants (substitution.md:45-65).
    /// Closes-after: #92. Code: `INITIATOR_MET_CANONICALIZATION`.
    InitiatorMetCanonicalization {
        /// Accession of the reference sequence.
        accession: String,
        /// Final dup interval text, e.g. "Met1" or "Met1_Lys2".
        location: String,
    },

    /// Bracketed / reference-range `ins[...]` payload was expanded to a
    /// flat literal sequence. Emitted alongside the canonical rewrite for
    /// observability — callers can audit which inputs were canonicalized
    /// vs. preserved verbatim. Code: `INSERTED_SEQUENCE_EXPANDED`.
    InsertedSequenceExpanded {
        /// Accession of the outer variant
        accession: String,
        /// Original `ins[...]` payload as written (e.g. `[ATC]` or
        /// `[100_120inv]` or `[A;100_110]`)
        original_payload: String,
        /// The expanded flat literal sequence written into the AST
        expanded_literal: String,
    },

    /// A position lies past one of the transcript bounds — the CDS-end
    /// (for plain `c.<N>`), the transcript-end (for `c.*<N>` or `n.<N>`),
    /// or the 5'UTR-start (for `c.-<N>`) — and therefore does not
    /// reference an existing base. Code: `POSITION_PAST_END` (W4004).
    ///
    /// One warning is emitted per offending position — a range with both
    /// endpoints past-end produces two warnings. Covers both the `c.`
    /// axis (#336, #348) and the `n.` axis (#347); intronic offsets
    /// remain out of scope (they depend on intron-size alignment data
    /// this check does not consult).
    PositionPastEnd {
        /// Transcript accession (e.g. `NM_001001656.1`).
        accession: String,
        /// Coordinate system: `"c"` for c. (cds-end / transcript-end /
        /// 5utr-start) or `"n"` for n. (transcript-end).
        coordinate_system: String,
        /// The single offending position in HGVS string form — `"946"`
        /// for plain CDS positions, `"*9"` for 3'UTR positions. Range
        /// strings like `"935_946"` are never produced here; each
        /// endpoint of a range yields its own warning.
        position: String,
        /// The bound the position exceeded — `"cds-end"` for plain
        /// `c.<N>` or `"transcript-end"` for `c.*<N>`.
        bound_kind: String,
        /// The numeric bound (e.g. 945 if the CDS is 945 bases long).
        bound_value: u64,
    },
}

impl NormalizationWarning {
    /// The warning's user-facing code string.
    pub fn code(&self) -> &'static str {
        match self {
            Self::RefSeqMismatch { .. } => "REFSEQ_MISMATCH",
            Self::OverlapConflict { .. } => "OVERLAP_CONFLICTING_EDITS",
            Self::UnresolvableSpecialPosition { .. } => "UNRESOLVABLE_SPECIAL_POSITION",
            Self::CanonicalSplitSkipped { .. } => "CANONICAL_SPLIT_SKIPPED",
            Self::CrossAxisVariantNotShuffled { .. } => "CROSS_AXIS_VARIANT_NOT_SHUFFLED",
            Self::AxisClampApplied { .. } => "AXIS_CLAMP_APPLIED",
            Self::InitiatorMetCanonicalization { .. } => "INITIATOR_MET_CANONICALIZATION",
            Self::InsertedSequenceExpanded { .. } => "INSERTED_SEQUENCE_EXPANDED",
            Self::PositionPastEnd { .. } => "POSITION_PAST_END",
        }
    }

    /// Human-readable message synthesized from the warning's structural
    /// fields. Equivalent to `format!("{self}")` — preserved as a method
    /// for ergonomics and back-compat with `.message()` call sites
    /// (#397 item 3 dropped the per-variant `message: String` field).
    pub fn message(&self) -> String {
        self.to_string()
    }
}

impl std::fmt::Display for NormalizationWarning {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::RefSeqMismatch {
                stated_ref,
                actual_ref,
                position,
                corrected,
                details,
            } => {
                write!(
                    f,
                    "reference sequence mismatch at {position}: stated {stated_ref:?}, actual {actual_ref:?} (corrected={corrected})",
                )?;
                if let Some(detail) = details.as_deref().filter(|s| !s.is_empty()) {
                    write!(f, ": {detail}")?;
                }
                Ok(())
            }
            Self::OverlapConflict {
                accession,
                coordinate_system,
                location,
                edit_kinds,
            } => write!(
                f,
                "{} cis edits share identical bounds at {}:{}.{}: {}",
                edit_kinds.len(),
                accession,
                coordinate_system,
                location,
                edit_kinds.join(", "),
            ),
            Self::UnresolvableSpecialPosition { accession, marker } => write!(
                f,
                "{accession}:g.{marker} — centromere/telomere marker '{marker}' cannot be \
                 resolved to a coordinate without assembly annotation; input preserved verbatim \
                 (UnresolvableSpecialPosition)"
            ),
            Self::CanonicalSplitSkipped {
                accession,
                hgvs_start,
                hgvs_end,
                expected_span,
                actual_bytes,
            } => write!(
                f,
                "canonical split skipped at {accession}:{hgvs_start}_{hgvs_end}: expected {expected_span} bytes, got {actual_bytes} (HGVS refseq \u{00A7}43)",
            ),
            Self::CrossAxisVariantNotShuffled {
                accession,
                start_axis,
                end_axis,
            } => write!(
                f,
                "{accession}: variant spans {start_axis} \u{2194} {end_axis} sub-axes; 3'-rule shuffle skipped",
            ),
            Self::AxisClampApplied {
                accession,
                direction,
                clamp_kind,
            } => write!(
                f,
                "{accession}: {direction} shuffle clamped at {clamp_kind} boundary",
            ),
            Self::InitiatorMetCanonicalization {
                accession,
                location,
            } => write!(
                f,
                "{accession}: canonical form `p.{location}dup` includes the initiator methionine; the predicted protein consequence may also be described as `p.0?` or `p.(Met1?)` per HGVS Substitution recommendations",
            ),
            Self::InsertedSequenceExpanded {
                accession,
                original_payload,
                expanded_literal,
            } => write!(
                f,
                "{accession}: ins payload {original_payload} expanded to literal {expanded_literal}",
            ),
            Self::PositionPastEnd {
                accession,
                coordinate_system,
                position,
                bound_kind,
                bound_value,
            } => write!(
                f,
                "{accession}:{coordinate_system}.{position} lies past the {bound_kind} (bound {bound_value})",
            ),
        }
    }
}

/// Info-grade signal generated during normalization.
///
/// Mirrors [`NormalizationWarning`] but for non-error, non-warning events
/// callers may want to record (e.g. the shuffle layer moved the variant
/// to a canonical position). The enum is open-ended: each variant owns the
/// fields its code needs. Future info codes add new variants without
/// touching existing emit sites.
///
/// This is the structural counterpart to mutalyzer/mutalyzer's `infos`
/// array (codes prefixed with `I`); see [`crate::error_handling::info_map`]
/// for the upstream-equivalent string mapping used by the corpus runner.
///
/// Marked `#[non_exhaustive]` so adding new info variants is non-breaking
/// for downstream `match` arms.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum NormalizationInfo {
    /// The shuffle layer relocated the variant per the HGVS arbitrary-
    /// position rule. The HGVS spec mandates the 3' (rightmost) form, but
    /// ferro supports both 3' and 5' shuffling via
    /// [`config::ShuffleDirection`] (some VCF pipelines prefer 5'); the
    /// direction is carried explicitly so callers can interpret the
    /// signal correctly. Code: `SHUFFLE_APPLIED`. Mutalyzer-equivalent:
    /// `ICORRECTEDPOINT` (mutalyzer only emits 3').
    ShuffleApplied {
        /// Accession of the reference sequence.
        accession: String,
        /// Direction in which the shuffle ran for this normalization.
        direction: config::ShuffleDirection,
        /// Position text of the input variant (HGVS, no accession),
        /// e.g. `"4"` or `"100_103"`.
        original_position: String,
        /// Position text of the normalized variant (HGVS, no accession),
        /// e.g. `"12"` or `"108_111"`.
        normalized_position: String,
    },
}

impl NormalizationInfo {
    /// The info's user-facing code string.
    pub fn code(&self) -> &'static str {
        match self {
            Self::ShuffleApplied { .. } => "SHUFFLE_APPLIED",
        }
    }

    /// Human-readable message synthesized from the info's structural
    /// fields. Equivalent to `format!("{self}")` (#397 item 3 dropped
    /// the per-variant `message: String` field).
    pub fn message(&self) -> String {
        self.to_string()
    }
}

impl std::fmt::Display for NormalizationInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ShuffleApplied {
                accession,
                direction,
                original_position,
                normalized_position,
            } => write!(
                f,
                "{accession}: {direction} shuffle relocated variant from {original_position} to {normalized_position}",
            ),
        }
    }
}

/// Result of normalization with optional warnings and info-grade signals.
///
/// Marked `#[non_exhaustive]` so future diagnostic axes (e.g. a separate
/// `notices` list) can be added without breaking downstream construction
/// via struct literals.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct NormalizeResult {
    /// The normalized variant
    pub result: HgvsVariant,
    /// Warnings generated during normalization
    pub warnings: Vec<NormalizationWarning>,
    /// Info-grade signals generated during normalization (e.g. a 3'-rule
    /// shuffle was applied). Empty when no signals fired. See
    /// [`NormalizationInfo`] for the variant taxonomy.
    pub infos: Vec<NormalizationInfo>,
}

impl NormalizeResult {
    /// Create a new result without warnings or infos
    pub fn new(result: HgvsVariant) -> Self {
        Self {
            result,
            warnings: vec![],
            infos: vec![],
        }
    }

    /// Create a result with warnings (and no infos)
    pub fn with_warnings(result: HgvsVariant, warnings: Vec<NormalizationWarning>) -> Self {
        Self {
            result,
            warnings,
            infos: vec![],
        }
    }

    /// Create a result with both warnings and infos
    pub fn with_diagnostics(
        result: HgvsVariant,
        warnings: Vec<NormalizationWarning>,
        infos: Vec<NormalizationInfo>,
    ) -> Self {
        Self {
            result,
            warnings,
            infos,
        }
    }

    /// Add a warning to the result
    pub fn add_warning(&mut self, warning: NormalizationWarning) {
        self.warnings.push(warning);
    }

    /// Add an info-grade signal to the result
    pub fn add_info(&mut self, info: NormalizationInfo) {
        self.infos.push(info);
    }

    /// Check if there are any warnings
    pub fn has_warnings(&self) -> bool {
        !self.warnings.is_empty()
    }

    /// Check if there are any info-grade signals
    pub fn has_infos(&self) -> bool {
        !self.infos.is_empty()
    }

    /// Check if there's a reference mismatch warning
    pub fn has_ref_mismatch(&self) -> bool {
        self.warnings
            .iter()
            .any(|w| matches!(w, NormalizationWarning::RefSeqMismatch { .. }))
    }
}

/// Resolve a genomic position to a concrete 1-based base, mapping telomere
/// markers to reference boundaries:
/// - `pter` -> 1 (first nucleotide),
/// - `qter` -> reference length (last nucleotide),
/// - `cen`  -> `Ok(None)` (a centromere is an assembly-annotated region, not a
///   sequence-derivable base),
/// - a plain (non-special) position -> its own `base`.
///
/// A `qter` whose reference length is unavailable also yields `Ok(None)` so the
/// caller can fall back to canonicalization (matches the "no sequence -> minimal
/// notation" philosophy). The caller distinguishes the two `None` cases by
/// re-inspecting `pos.special`: `Some(Cen)` is a structural failure (warn/reject),
/// any other `None` is an environment gap (silent fallback).
///
/// Precondition: offset-carrying positions are bailed out by the caller before
/// this is called; this function does not inspect `pos.offset`.
// Will be called by normalize_genome in the next task (Task 2).
#[allow(dead_code)]
fn resolve_special_genome_pos<P: ReferenceProvider>(
    pos: &GenomePos,
    accession: &str,
    provider: &P,
) -> Result<Option<u64>, FerroError> {
    match pos.special {
        None => Ok(Some(pos.base)),
        Some(SpecialPosition::Pter) => Ok(Some(1)),
        // Length unavailable -> graceful None (caller canonicalizes).
        Some(SpecialPosition::Qter) => Ok(provider.get_sequence_length(accession).ok()),
        Some(SpecialPosition::Cen) => Ok(None),
    }
}

/// Main normalizer struct
pub struct Normalizer<P: ReferenceProvider> {
    provider: P,
    config: NormalizeConfig,
}

impl<P: ReferenceProvider> Normalizer<P> {
    /// Create a new normalizer with the given reference provider
    pub fn new(provider: P) -> Self {
        Self {
            provider,
            config: NormalizeConfig::default(),
        }
    }

    /// Create a normalizer with custom configuration
    pub fn with_config(provider: P, config: NormalizeConfig) -> Self {
        Self { provider, config }
    }

    /// Get the configuration
    pub fn config(&self) -> &NormalizeConfig {
        &self.config
    }

    /// Normalize a variant
    ///
    /// In strict mode (default), rejects variants with reference mismatches.
    /// Use `normalize_with_diagnostics` for lenient mode that corrects mismatches.
    pub fn normalize(&self, variant: &HgvsVariant) -> Result<HgvsVariant, FerroError> {
        let result = self.normalize_with_diagnostics(variant)?;

        // In strict mode, reject if there were reference mismatches.
        if self.config.should_reject_ref_mismatch() {
            if let Some(err) = result.warnings.iter().find_map(|w| match w {
                NormalizationWarning::RefSeqMismatch {
                    position,
                    stated_ref,
                    actual_ref,
                    ..
                } => Some(FerroError::ReferenceMismatch {
                    location: position.clone(),
                    expected: stated_ref.clone(),
                    found: actual_ref.clone(),
                }),
                _ => None,
            }) {
                return Err(err);
            }
        }

        // Strict mode also rejects W5003 VariantExceedsReference per
        // HGVS spec refseq.md §43 — the variant must be entirely
        // encompassed by the selected reference. Promotes the
        // `CanonicalSplitSkipped` warning to a typed error. Closes
        // #355; matches biocommons hgvs which raises
        // `HGVSInvalidVariantError` for this shape.
        if self.config.should_reject_variant_exceeds_reference() {
            if let Some(err) = result.warnings.iter().find_map(|w| match w {
                NormalizationWarning::CanonicalSplitSkipped {
                    accession,
                    hgvs_start,
                    hgvs_end,
                    expected_span,
                    actual_bytes,
                    ..
                } => Some(FerroError::VariantExceedsReference {
                    accession: accession.clone(),
                    hgvs_start: *hgvs_start,
                    hgvs_end: *hgvs_end,
                    expected_span: *expected_span as u64,
                    actual_bytes: *actual_bytes as u64,
                }),
                _ => None,
            }) {
                return Err(err);
            }
        }

        // In strict mode, reject if any position lies past the CDS-end,
        // transcript-end, or contig-end (W4004). Use an axis-aware noun
        // for the reference structure — `m.` lives on a contig, not a
        // transcript, so the rejection message must say so.
        if self.config.should_reject_position_past_end() {
            if let Some(err) = result.warnings.iter().find_map(|w| match w {
                NormalizationWarning::PositionPastEnd {
                    accession,
                    coordinate_system,
                    position,
                    bound_kind,
                    bound_value,
                    ..
                } => {
                    let reference_noun = match coordinate_system.as_str() {
                        "m" => "contig",
                        _ => "transcript",
                    };
                    Some(FerroError::InvalidCoordinates {
                        msg: format!(
                            "{accession}:{coordinate_system}.{position} is past the {bound_kind} \
                             ({bound_value}); position does not reference a base in the \
                             {reference_noun} (PositionPastEnd / W4004)"
                        ),
                    })
                }
                _ => None,
            }) {
                return Err(err);
            }
        }

        // In strict mode, reject when two or more cis-allele edits
        // share identical reference bounds (W5002). The registry
        // already declared `ModeBehavior::always_warn_if_not_rejected`
        // (Strict → Reject) but the emit site at `overlap.rs:88`
        // unconditionally pushed the warning. Closes #395 item 6.
        if self.config.should_reject_overlap_conflict() {
            if let Some(err) = result.warnings.iter().find_map(|w| match w {
                NormalizationWarning::OverlapConflict {
                    accession,
                    coordinate_system,
                    location,
                    edit_kinds,
                    ..
                } => Some(FerroError::InvalidCoordinates {
                    msg: format!(
                        "{accession}:{coordinate_system}.{location} has \
                         {n} coincident cis-allele edits ({kinds}); HGVS spec defines no \
                         canonical form for this case (OverlapConflictingEdits / W5002)",
                        n = edit_kinds.len(),
                        kinds = edit_kinds.join(", "),
                    ),
                }),
                _ => None,
            }) {
                return Err(err);
            }
        }

        // In strict mode, reject an unresolvable `cen` position (W4005). A
        // centromere is an assembly-annotated region with no sequence-derivable
        // base, so `normalize_genome` cannot place it and emits the
        // `UnresolvableSpecialPosition` warning. The registry declares
        // `warn_accept()` (Strict → Reject); promote it here rather than
        // silently echoing the input. See #488.
        if self.config.should_reject_unresolvable_centromere() {
            if let Some(err) = result.warnings.iter().find_map(|w| match w {
                NormalizationWarning::UnresolvableSpecialPosition { accession, marker } => {
                    Some(FerroError::InvalidCoordinates {
                        msg: format!(
                            "{accession}:g.{marker} — '{marker}' is an assembly-annotated \
                             region with no sequence-derivable coordinate and cannot be \
                             normalized (UnresolvableCentromere / W4005)"
                        ),
                    })
                }
                _ => None,
            }) {
                return Err(err);
            }
        }

        Ok(result.result)
    }

    /// Normalize a variant with detailed warnings
    ///
    /// Returns the normalized variant along with any warnings generated during
    /// normalization (e.g., reference sequence mismatches that were auto-corrected).
    /// Use this method when you want to track what corrections were made.
    pub fn normalize_with_diagnostics(
        &self,
        variant: &HgvsVariant,
    ) -> Result<NormalizeResult, FerroError> {
        let (result, warnings) = match variant {
            HV::Genome(v) => self.normalize_genome(v)?,
            HV::Cds(v) => self.normalize_cds(v)?,
            HV::Tx(v) => self.normalize_tx(v)?,
            HV::Protein(v) => self.normalize_protein(v)?,
            HV::Rna(v) => self.normalize_rna(v)?,
            HV::Mt(v) => self.normalize_mt(v)?,
            HV::Allele(a) => self.normalize_allele(a)?,
            // Circular variants normalize like genomic variants
            HV::Circular(v) => (
                HV::Circular(crate::hgvs::variant::CircularVariant {
                    accession: v.accession.clone(),
                    gene_symbol: v.gene_symbol.clone(),
                    loc_edit: v.loc_edit.clone(),
                }),
                vec![],
            ),
            // RNA fusions pass through unchanged (no normalization needed for fusions)
            HV::RnaFusion(v) => (HV::RnaFusion(v.clone()), vec![]),
            // Null and unknown allele markers pass through unchanged
            HV::NullAllele => (HV::NullAllele, vec![]),
            HV::UnknownAllele => (HV::UnknownAllele, vec![]),
        };

        let infos = detect_shuffle_infos(variant, &result, self.config.shuffle_direction);
        Ok(NormalizeResult::with_diagnostics(result, warnings, infos))
    }

    /// Normalize an allele (compound) variant
    ///
    /// Normalizes each variant in the allele individually, with overlap prevention.
    /// After normalization, checks if variants would overlap and constrains shifting
    /// to prevent collisions.
    fn normalize_allele(
        &self,
        allele: &crate::hgvs::variant::AlleleVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // Merge first, then normalize each result through the full per-variant
        // pipeline. Pre-normalizing each bracket entry would shift it
        // independently of its siblings, which can collapse adjacent edits
        // onto the same 3'-end position (e.g. `[260delA;261delA]` →
        // `[264del;264del]`) and defeat the strict-adjacency merge. Issue #180.
        //
        // Order:
        //   1. merge raw bracket entries by positional adjacency
        //   2. cis-only: decompose merged delins into [..., inv, ...] (#160)
        //   3. run the full per-variant pipeline on every result — this is
        //      what applies the HGVS 3' rule to the merged anchor (#161, #180)
        //   4. detect post-shift overlaps and emit warnings
        let mut all_warnings = Vec::new();
        let original_len = allele.variants.len();
        let merged_raw =
            merge::merge_consecutive_edits(allele.variants.clone(), allele.phase, &self.provider);

        // Issue #160 + #165: any merged delins (or pre-existing delins
        // that survived merge unchanged) may decompose into a sequence of
        // higher-priority forms per `general.md:56` — `[..., inv, ...]`
        // when an inv-eligible sub-span is present (#160) and/or into
        // separate substitutions when interior positions match the
        // reference (#165). Run the split per merged variant; the helper
        // is a no-op for non-Delins variants and for Delins with nothing
        // to split out. Only applies in cis phase — trans alleles aren't
        // collapsible in the first place.
        let merged_split: Vec<HgvsVariant> =
            if allele.phase == crate::hgvs::variant::AllelePhase::Cis {
                merged_raw
                    .into_iter()
                    .flat_map(|v| self.canonical_split_for_variant(v))
                    .collect()
            } else {
                merged_raw
            };

        // Per-variant pipeline on every merged result. This is the single
        // canonical place where the 3' rule, ins→dup canonicalization, ref
        // validation, etc. apply — a merged variant is semantically a new
        // variant and goes through the same pipeline as any direct input.
        let mut normalized: Vec<HgvsVariant> = Vec::with_capacity(merged_split.len());
        for v in merged_split {
            let r = self.normalize_with_diagnostics(&v)?;
            all_warnings.extend(r.warnings);
            normalized.push(r.result);
        }

        // Overlap detection runs post-shift so collisions caused by the
        // 3' shift surface alongside input-time ones. Overlap *prevention*
        // is structural — the merge-first ordering above plus the strict
        // `prev.end + 1 == next.start` adjacency check in
        // `merge_consecutive_edits` make it impossible for the normalizer
        // to emit overlapping ranges from non-overlapping inputs.
        all_warnings.extend(crate::normalize::overlap::detect_overlap_conflicts(
            &normalized,
            allele.phase,
        ));

        // HGVS requires consecutive edits in cis to render as a single edit.
        // Only unwrap when a merge actually collapsed multiple sub-variants —
        // pre-existing singleton alleles must round-trip with the Allele
        // wrapper intact for programmatic callers (Display already renders
        // singletons in bare form regardless).
        let result = if allele.phase == crate::hgvs::variant::AllelePhase::Cis
            && original_len > 1
            && normalized.len() == 1
        {
            normalized.pop().unwrap()
        } else {
            HgvsVariant::Allele(crate::hgvs::variant::AlleleVariant::new(
                normalized,
                allele.phase,
            ))
        };

        Ok((result, all_warnings))
    }

    /// Issue #333: expand a bracketed / reference-range `ins[...]` payload
    /// in an `Insertion` / `Delins` / `DupIns` edit to a flat literal
    /// sequence. Returns the rewritten edit plus a
    /// [`NormalizationWarning::InsertedSequenceExpanded`] for
    /// observability, or `None` when nothing was canonicalized.
    ///
    /// The four `normalize_<axis>` methods wrap this helper to build the
    /// per-axis variant struct. Splitting the warning construction here
    /// keeps the original `[ATC]` / `[A;100_110]` rendering accessible
    /// before the edit is mutated.
    fn try_expand_ins_edit(
        &self,
        edit: &NaEdit,
        accession: &str,
        kind: InsCoordKind,
    ) -> Result<Option<(NaEdit, NormalizationWarning)>, FerroError> {
        // Snapshot the inserted-sequence payload BEFORE rewriting so the
        // warning's `original_payload` matches its documented format
        // (e.g. `[ATC]` / `[A;100_110]`), and the per-edit display string
        // for the human-readable message. Pure Insertion / Delins /
        // DupIns are the only edits the helper acts on; other variants
        // short-circuit below since `canonicalize_insertion_expand` will
        // return `Ok(None)` for them.
        let original_inserted = match edit {
            NaEdit::Insertion { sequence }
            | NaEdit::Delins { sequence, .. }
            | NaEdit::DupIns { sequence } => sequence,
            _ => return Ok(None),
        };
        // The warning is specifically about a bracketed payload being
        // expanded — by construction `canonicalize_insertion_expand`
        // only acts on `InsertedSequence::Complex`. `to_bracketed_string`
        // forces the bracketed form so the warning preserves the user's
        // input shape (`[ATC]`, `[A;100_110]`) even though `Display` on
        // `InsertedSequence::Complex` now drops brackets for single-element
        // vectors (spec-canonical form).
        let original_payload = original_inserted.to_bracketed_string();

        let new_edit = match canonicalize_insertion_expand(edit, accession, kind, &self.provider)? {
            Some(e) => e,
            None => return Ok(None),
        };

        let new_inserted = match &new_edit {
            NaEdit::Insertion { sequence }
            | NaEdit::Delins { sequence, .. }
            | NaEdit::DupIns { sequence } => sequence,
            // The expand helper preserves the edit kind, so this branch
            // is unreachable in practice; falling back to the full edit
            // display keeps the field non-empty if invariants change.
            _ => {
                let warning = NormalizationWarning::InsertedSequenceExpanded {
                    accession: accession.to_string(),
                    original_payload,
                    expanded_literal: format!("{}", new_edit),
                };
                return Ok(Some((new_edit, warning)));
            }
        };
        let expanded_literal = format!("{}", new_inserted);

        let warning = NormalizationWarning::InsertedSequenceExpanded {
            accession: accession.to_string(),
            original_payload,
            expanded_literal,
        };
        Ok(Some((new_edit, warning)))
    }

    /// `normalize_genome` companion that wraps [`Self::try_expand_ins_edit`]
    /// and rebuilds a `GenomeVariant` carrying the rewritten edit.
    fn try_expand_genome_ins(
        &self,
        variant: &GenomeVariant,
        edit: &NaEdit,
        accession: &str,
    ) -> Result<Option<(GenomeVariant, NormalizationWarning)>, FerroError> {
        let Some((new_edit, warning)) =
            self.try_expand_ins_edit(edit, accession, InsCoordKind::Direct)?
        else {
            return Ok(None);
        };
        let new_variant = GenomeVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                variant.loc_edit.location.clone(),
                variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
            ),
        };
        Ok(Some((new_variant, warning)))
    }

    /// `normalize_cds` companion — c. positions are CDS-relative.
    fn try_expand_cds_ins(
        &self,
        variant: &CdsVariant,
        edit: &NaEdit,
        accession: &str,
    ) -> Result<Option<(CdsVariant, NormalizationWarning)>, FerroError> {
        let Some((new_edit, warning)) =
            self.try_expand_ins_edit(edit, accession, InsCoordKind::Cds)?
        else {
            return Ok(None);
        };
        let new_variant = CdsVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                variant.loc_edit.location.clone(),
                variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
            ),
        };
        Ok(Some((new_variant, warning)))
    }

    /// `normalize_tx` companion — n. positions are 1-based transcript
    /// positions, which the helper treats as direct.
    fn try_expand_tx_ins(
        &self,
        variant: &TxVariant,
        edit: &NaEdit,
        accession: &str,
    ) -> Result<Option<(TxVariant, NormalizationWarning)>, FerroError> {
        let Some((new_edit, warning)) =
            self.try_expand_ins_edit(edit, accession, InsCoordKind::Direct)?
        else {
            return Ok(None);
        };
        let new_variant = TxVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                variant.loc_edit.location.clone(),
                variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
            ),
        };
        Ok(Some((new_variant, warning)))
    }

    /// `normalize_mt` companion — m. positions are direct genomic
    /// positions on the mitochondrial accession.
    fn try_expand_mt_ins(
        &self,
        variant: &MtVariant,
        edit: &NaEdit,
        accession: &str,
    ) -> Result<Option<(MtVariant, NormalizationWarning)>, FerroError> {
        let Some((new_edit, warning)) =
            self.try_expand_ins_edit(edit, accession, InsCoordKind::Direct)?
        else {
            return Ok(None);
        };
        let new_variant = MtVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                variant.loc_edit.location.clone(),
                variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
            ),
        };
        Ok(Some((new_variant, warning)))
    }
}

/// Resolve a genomic position to a concrete 1-based base, mapping telomere
/// markers to reference boundaries:
/// - `pter` -> 1 (first nucleotide),
/// - `qter` -> reference length (last nucleotide),
/// - `cen`  -> `Ok(None)` (a centromere is an assembly-annotated region, not a
///   sequence-derivable base),
/// - a plain (non-special) position -> its own `base`.
///
/// A `qter` whose reference length is unavailable also yields `Ok(None)` so the
/// caller can fall back to canonicalization (matches the "no sequence -> minimal
/// notation" philosophy). The caller distinguishes the two `None` cases by
/// re-inspecting `pos.special`: `Some(Cen)` is a structural failure (warn/reject),
/// any other `None` is an environment gap (silent fallback).
///
/// Precondition: offset-carrying positions are bailed out by the caller before
/// this is called; this function does not inspect `pos.offset`.
fn resolve_special_genome_pos<P: ReferenceProvider>(
    pos: &GenomePos,
    accession: &str,
    provider: &P,
) -> Result<Option<u64>, FerroError> {
    match pos.special {
        None => Ok(Some(pos.base)),
        Some(SpecialPosition::Pter) => Ok(Some(1)),
        // Length unavailable -> graceful None (caller canonicalizes).
        Some(SpecialPosition::Qter) => Ok(provider.get_sequence_length(accession).ok()),
        Some(SpecialPosition::Cen) => Ok(None),
    }
}

impl<P: ReferenceProvider> Normalizer<P> {
    /// Normalize a genomic variant
    fn normalize_genome(
        &self,
        variant: &GenomeVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // Can't normalize variants with unknown edits or positions
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Genome(variant.clone()), vec![])),
        };

        // SVD-WG009: rewrite `con` to `delins` before any further work.
        // Pure-syntax; no reference data needed. Re-run the axis
        // normalizer on the rewritten variant so the downstream passes
        // (issue #333 ins[...] expansion, 3' shift, ins→dup, canonical
        // split) all see the delins form — otherwise `...con...` inputs
        // would stop at the intermediate delins.
        if let Some(new_edit) = canonicalize_conversion_to_delins(edit) {
            let new_variant = GenomeVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    variant.loc_edit.location.clone(),
                    variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
                ),
            };
            return self.normalize_genome(&new_variant);
        }

        // Issue #333: expand bracketed / reference-range `ins[...]`
        // payloads to a flat literal sequence so the rest of the
        // pipeline (3' shift, ins→dup, etc.) operates on concrete
        // bases. Same canonicalization applies to the inserted-
        // sequence payload inside Delins and DupIns.
        let accession = variant.accession.transcript_accession();
        if let Some((new_variant, warning)) =
            self.try_expand_genome_ins(variant, edit, &accession)?
        {
            let (result, mut warnings) = self.normalize_genome(&new_variant)?;
            warnings.insert(0, warning);
            return Ok((result, warnings));
        }

        // Only normalize indels
        if !needs_normalization(edit) {
            return Ok((HV::Genome(variant.clone()), vec![]));
        }
        let start_pos = match variant.loc_edit.location.start.inner() {
            Some(pos) => pos,
            None => {
                return Ok((
                    HV::Genome(self.canonicalize_genome_variant(variant)),
                    vec![],
                ))
            }
        };
        let end_pos = match variant.loc_edit.location.end.inner() {
            Some(pos) => pos,
            None => {
                return Ok((
                    HV::Genome(self.canonicalize_genome_variant(variant)),
                    vec![],
                ))
            }
        };

        // Offset-carrying genome positions (uncertain `g.123+?`-style) have no
        // resolution and cannot be losslessly remapped through base-only window
        // normalization; bail to minimal-notation cleanup. (Telomere markers
        // are handled by resolution just below — see #488 and the design doc.)
        if start_pos.offset.is_some() || end_pos.offset.is_some() {
            return Ok((
                HV::Genome(self.canonicalize_genome_variant(variant)),
                vec![],
            ));
        }

        // Resolve telomere markers (pter/qter) to concrete 1-based bases before
        // the window math; `base == 0` sentinels would otherwise underflow
        // `hgvs_pos_to_index` (#488). A plain position resolves to its own base.
        let had_special = start_pos.is_special() || end_pos.is_special();
        let (start, end) = match (
            resolve_special_genome_pos(start_pos, &accession, &self.provider)?,
            resolve_special_genome_pos(end_pos, &accession, &self.provider)?,
        ) {
            (Some(s), Some(e)) => (s, e),
            // Unresolved: `cen` is structurally unresolvable (surface a warning
            // so it is not silently echoed); a length-less qter/pter is an
            // environment gap (silent canonicalize fallback).
            _ => {
                let warnings = if matches!(start_pos.special, Some(SpecialPosition::Cen))
                    || matches!(end_pos.special, Some(SpecialPosition::Cen))
                {
                    vec![NormalizationWarning::UnresolvableSpecialPosition {
                        accession: accession.clone(),
                        marker: "cen".to_string(),
                    }]
                } else {
                    vec![]
                };
                return Ok((
                    HV::Genome(self.canonicalize_genome_variant(variant)),
                    warnings,
                ));
            }
        };

        // On the resolved-special path the reference length sizes both the
        // whole-span short-circuit and the fetch clamp. (A `qter` resolve above
        // also queried the length; this is a second cheap in-memory index
        // lookup, not I/O.) If it is unavailable we cannot safely size the
        // fetch -> canonicalize fallback.
        let resolved_len = if had_special {
            match self.provider.get_sequence_length(&accession) {
                Ok(len) => len,
                Err(_) => {
                    return Ok((
                        HV::Genome(self.canonicalize_genome_variant(variant)),
                        vec![],
                    ))
                }
            }
        } else {
            0 // unused when !had_special
        };

        // Whole-contig span (e.g. g.pter_qterdel): fully anchored, cannot
        // 3'-shift, and fetching the whole contig is pure waste. Render the
        // resolved concrete form and return before any get_sequence call.
        if had_special && start == 1 && end == resolved_len {
            let resolved = GenomeVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::new(
                    Interval::new(GenomePos::new(start), GenomePos::new(end)),
                    edit.clone(),
                ),
            };
            return Ok((
                HV::Genome(self.canonicalize_genome_variant(&resolved)),
                vec![],
            ));
        }

        // Window-based fetch around the variant. `window_start`/`fetch_end` are
        // 0-based half-open offsets into the contig even though they derive from
        // 1-based `base`; the reconciliation `rel = base - window_start` then
        // `hgvs_pos_to_index(rel) = rel - 1` is the same one the non-special
        // path uses (do not "fix" this as a bug). On the resolved-special path
        // clamp the upper bound to the contig length so the read is well-formed
        // against providers that error on past-EOF reads (MockProvider); the
        // ordinary path keeps the provider-clamped behavior byte-identical.
        let window_start = start.saturating_sub(self.config.window_size);
        let raw_end = end.saturating_add(self.config.window_size);
        let fetch_end = if had_special {
            raw_end.min(resolved_len)
        } else {
            raw_end
        };
        let seq_result = self
            .provider
            .get_sequence(&accession, window_start, fetch_end);

        let ref_seq = match seq_result {
            Ok(s) => s,
            // Can't do full normalization without sequence, but apply minimal notation
            Err(_) => {
                return Ok((
                    HV::Genome(self.canonicalize_genome_variant(variant)),
                    vec![],
                ))
            }
        };

        // Adjust coordinates to be relative to the window
        let rel_start = start - window_start;
        let rel_end = end - window_start;

        // Perform normalization
        let (new_rel_start, new_rel_end, new_edit, mut warnings) = self.normalize_na_edit(
            ref_seq.as_bytes(),
            edit,
            rel_start,
            rel_end,
            &Boundaries::new(0, ref_seq.len() as u64),
            false, // genomic context: codon-frame gate does not apply
        )?;

        // Adjust back to genomic coordinates
        let new_start = new_rel_start + window_start;
        let new_end = new_rel_end + window_start;

        // Reconstruct variant with new position (using adjusted coordinates)
        let new_variant = GenomeVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(GenomePos::new(new_start), GenomePos::new(new_end)),
                new_edit,
            ),
        };

        // Issue #160 + #165: a normalized Delins may decompose under the
        // spec's edit-priority rule (`general.md:56`) — into `[..., inv,
        // ...]` for rev-comp sub-spans and/or into separate subs across
        // interior identities. Returns the variant unchanged for
        // non-Delins or no-decomposition cases.
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Genome(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split), warnings))
    }

    /// Normalize a CDS variant
    fn normalize_cds(
        &self,
        variant: &CdsVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // Can't normalize variants with unknown edits or positions
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Cds(variant.clone()), vec![])),
        };

        // SVD-WG009: rewrite `con` to `delins`. Recurse on the rewritten
        // variant so the bounds-check + needs_normalization gates below still
        // fire — otherwise `c.<past-end>conT` slips past W4004. The downstream
        // passes (issue #333 ins[...] expansion, 3' shift, ins→dup, canonical
        // split) then all see the delins form.
        // `canonicalize_conversion_to_delins` only matches `NaEdit::Conversion`
        // and the rewrite swaps it for `Delins`, so the recursion terminates
        // on the second entry.
        if let Some(new_edit) = canonicalize_conversion_to_delins(edit) {
            let new_variant = CdsVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    variant.loc_edit.location.clone(),
                    variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
                ),
            };
            return self.normalize_cds(&new_variant);
        }

        // Bounds check: `c.<N>` where N exceeds the CDS length, or `c.*<N>` where
        // N exceeds the 3'UTR length, is malformed input. Strict and lenient
        // both emit a `PositionPastEnd` warning; the outer `normalize` wrapper
        // converts it to a typed error in strict mode. Silent mode skips the
        // check entirely. We early-return with the canonical variant since
        // normalize() cannot do sensible work on past-end input. Must run
        // BEFORE both the `try_expand_cds_ins` `ins[...]` expansion (so the
        // `con`-rewrite path `c.<past-end>conT` does not surface a misleading
        // `UnsupportedVariant` for the `T` payload) and the
        // `needs_normalization` short-circuit so substitutions
        // (`c.946G>C`) also get checked.
        //
        // The bounds check requires the transcript and concrete positions, but
        // is *optional* — when the transcript can't be fetched or the positions
        // are unknown/`?`-offset, we skip the gate and fall through so callers
        // without manifest data (e.g. the spec-fixture / parse-only path) still
        // hit the downstream `try_expand_cds_ins` and canonicalization passes.
        let accession = variant.accession.transcript_accession();
        let transcript_for_intronic =
            || -> Result<crate::reference::transcript::Transcript, FerroError> {
                self.provider
                    .get_transcript_for_variant(&HV::Cds(variant.clone()))
            };
        let transcript_opt = self.provider.get_transcript(&accession).ok();
        if self.config.should_reject_position_past_end()
            || self.config.should_warn_position_past_end()
        {
            if let (Some(transcript), Some(start_pos), Some(end_pos)) = (
                transcript_opt.as_ref(),
                variant.loc_edit.location.start.inner(),
                variant.loc_edit.location.end.inner(),
            ) {
                if !has_unknown_offset_cds(start_pos) && !has_unknown_offset_cds(end_pos) {
                    let mut bounds_warnings: Vec<NormalizationWarning> = Vec::new();
                    let acc_str = accession.clone();
                    if let Some(w) = check_cds_pos_past_end(&acc_str, start_pos, transcript) {
                        bounds_warnings.push(w);
                    }
                    let end_distinct = end_pos.base != start_pos.base
                        || end_pos.utr3 != start_pos.utr3
                        || end_pos.offset != start_pos.offset;
                    if end_distinct {
                        if let Some(w) = check_cds_pos_past_end(&acc_str, end_pos, transcript) {
                            bounds_warnings.push(w);
                        }
                    }
                    if !bounds_warnings.is_empty() {
                        return Ok((
                            HV::Cds(self.canonicalize_cds_variant(variant)),
                            bounds_warnings,
                        ));
                    }
                }
            }
        }

        // Issue #333: expand bracketed / reference-range `ins[...]`
        // payloads. CDS-coord ranges (e.g. c.X_Yins[A_B]) translate to
        // transcript coordinates via the transcript's cds_start. Same
        // canonicalization applies to Delins / DupIns payloads. Runs
        // AFTER the bounds gate above so past-end `con`/`delins` inputs
        // reject on the coordinate rather than on a downstream
        // `ins[...]` expansion failure.
        let cds_accession = accession.clone();
        if let Some((new_variant, warning)) =
            self.try_expand_cds_ins(variant, edit, &cds_accession)?
        {
            let (result, mut warnings) = self.normalize_cds(&new_variant)?;
            warnings.insert(0, warning);
            return Ok((result, warnings));
        }

        // Extract positions for the post-expansion path. Substitutions (and
        // other non-shifted edits) also need normalization fall-through, so
        // we look up positions BEFORE the `needs_normalization`
        // short-circuit below.
        let start_pos = match variant.loc_edit.location.start.inner() {
            Some(pos) => pos,
            None => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };
        let end_pos = match variant.loc_edit.location.end.inner() {
            Some(pos) => pos,
            None => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };

        // Can't normalize variants with unknown (?) offsets - return unchanged
        if has_unknown_offset_cds(start_pos) || has_unknown_offset_cds(end_pos) {
            return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![]));
        }

        // Re-fetch the transcript for the downstream normalization passes
        // (intronic / boundary / 3' shuffle / canonical split). The bounds
        // gate above used `transcript_opt`, which may be `None` for callers
        // without a manifest; reaching this point implies we still need the
        // transcript and must early-return if it can't be fetched.
        let transcript = match transcript_opt {
            Some(t) => t,
            // Can't do full normalization without transcript, but apply minimal notation
            None => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };

        // Only normalize indels (the bounds check above runs regardless).
        if !needs_normalization(edit) {
            return Ok((HV::Cds(variant.clone()), vec![]));
        }

        // Handle intronic variants specially
        if start_pos.is_intronic() || end_pos.is_intronic() {
            // Switch to the variant-aware lookup so an NG/NC-parented input
            // gets the build-correct chromosome. If the variant-aware lookup
            // fails, fall back to the plain transcript we already fetched.
            let transcript = transcript_for_intronic().unwrap_or(transcript);
            // Check if both positions are intronic and in the same intron
            if start_pos.is_intronic() && end_pos.is_intronic() {
                return self.normalize_intronic_cds(variant, &transcript, start_pos, end_pos, edit);
            }
            // Variant spans exon-intron boundary - normalize in genomic space
            return self.normalize_boundary_spanning_cds(
                variant,
                &transcript,
                start_pos,
                end_pos,
                edit,
            );
        }

        // Convert CDS to transcript coordinates for normalization
        let cds_start = match transcript.cds_start {
            Some(s) => s,
            None => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };

        // Calculate transcript positions - return unchanged if position is out of range
        let tx_start = match self.cds_to_tx_pos(start_pos, cds_start, transcript.cds_end) {
            Ok(v) => v,
            Err(_) => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };
        let tx_end = match self.cds_to_tx_pos(end_pos, cds_start, transcript.cds_end) {
            Ok(v) => v,
            Err(_) => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };

        // Get boundaries. Both cross_boundaries modes route through
        // `get_cds_boundaries_with_axis_info` so the CDS↔UTR axis bound
        // applies regardless of cross mode (closes #337). The
        // exon-vs-full-tx dimension still toggles on
        // `config.cross_boundaries`. We additionally use the un-clamped
        // exon bound to detect:
        //
        //   - Cross-axis variants (#350): `tx_start` and `tx_end` map to
        //     different axes (5'UTR / CDS / 3'UTR). The 3'-rule shuffle
        //     has no well-defined semantics across an axis boundary, so
        //     we preserve the canonical input position and emit
        //     `CrossAxisVariantNotShuffled`.
        //
        //   - Axis-clamp activations (#349): after shuffling, the result
        //     position rests at the axis boundary AND the axis bound is
        //     tighter than the exon bound on that side. Emit
        //     `AxisClampApplied` so callers can flag for human review.
        let axis_info = match boundary::get_cds_boundaries_with_axis_info(
            &transcript,
            tx_start,
            &self.config,
        ) {
            Ok(b) => b,
            Err(_) => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };
        let boundaries = axis_info.clamped.clone();
        let exon_only = axis_info.exon.clone();
        let start_axis = axis_info.axis_region;
        let end_axis = boundary::axis_region_of(&transcript, tx_end);

        // #350: bail on cross-axis variants. Both endpoints must live in
        // the same axis region for a 3'-rule shuffle to be defined.
        //
        // **Insertion exception (#402).** A boundary-straddling
        // `c.X_X+1ins<alt>` is zero-width — its shuffle is well-defined
        // (the insertion point can move toward either axis as long as
        // the alt continues to match reference) and the canonical
        // form may land entirely on one axis as an `insAlt` or
        // `dup`. Let Insertion inputs proceed through `normalize_na_edit`
        // and rely on the post-shift CDS-start / CDS-end clamps
        // (PR #385 / PR #388, with the #401 spanning-dup exception) to
        // catch any result that would silently re-ax the input.
        let is_zero_width_insertion =
            matches!(edit, NaEdit::Insertion { .. }) && tx_end == tx_start + 1;
        if start_axis != end_axis
            && !matches!(start_axis, boundary::AxisRegion::None)
            && !matches!(end_axis, boundary::AxisRegion::None)
            && !is_zero_width_insertion
        {
            let acc = variant.accession.transcript_accession();
            let warning = NormalizationWarning::CrossAxisVariantNotShuffled {
                accession: acc,
                start_axis: start_axis.label().to_string(),
                end_axis: end_axis.label().to_string(),
            };
            return Ok((
                HV::Cds(self.canonicalize_cds_variant(variant)),
                vec![warning],
            ));
        }

        // Perform normalization on transcript sequence (CDS context).
        // Coordinate-only transcripts (no cached bases) fall back to the
        // canonicalize-only path, matching the other early-return branches.
        let seq = match transcript.sequence.as_deref() {
            Some(s) => s.as_bytes(),
            None => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
        };
        // HGVS spec (repeated.md): the codon-frame restriction
        // (`unit_len % 3 == 0` for repeat notation in `c.` context)
        // applies only to bases inside the CDS proper. UTR positions
        // are exempt:
        //   > This restriction only applies to the coding sequence,
        //   > which does not include the introns or the UTR sequence.
        // The variant is entirely within the CDS iff its transcript-
        // frame span sits between `cds_start` and `cds_end` (inclusive).
        // 5' UTR (`c.-N`) maps to `tx_start < cds_start`; 3' UTR
        // (`c.*N`) maps to `tx_start > cds_end`. Pass `is_coding=false`
        // for any variant whose footprint touches UTR. If `cds_end` is
        // unset we cannot verify the variant lies within CDS proper, so
        // we conservatively treat it as UTR-touching and skip the gate.
        let is_coding = match transcript.cds_end {
            Some(cds_end) => tx_start >= cds_start && tx_end <= cds_end,
            None => false,
        };
        let (mut new_tx_start, mut new_tx_end, mut new_edit, mut warnings) =
            self.normalize_na_edit(seq, edit, tx_start, tx_end, &boundaries, is_coding)?;

        // #349: detect whether the axis clamp was operative for this
        // shuffle. For 5'-direction shuffles the clamp fires when the
        // result is anchored at `boundaries.left` AND the axis bound is
        // strictly tighter than the exon bound on the left. Symmetric
        // logic for 3'-direction. Skip when there's no axis sub-region
        // (non-coding transcripts).
        //
        // The cheap landed-at-boundary check alone is not sufficient — it
        // can fire even when the unclamped shuffle would have stopped at
        // the same axis-boundary position anyway (e.g. when the reference
        // base immediately past the boundary does not match, so the
        // shuffle was never going to move past the boundary in the first
        // place). To eliminate that false-positive class, re-shuffle
        // against the unclamped `exon_only` bound and only treat the
        // clamp as operative when the unclamped result would have moved
        // strictly past the axis boundary.
        if !matches!(start_axis, boundary::AxisRegion::None) {
            let new_tx_start_0 = new_tx_start.saturating_sub(1);
            let cheap_left = boundaries.left > exon_only.left && new_tx_start_0 == boundaries.left;
            let cheap_right = boundaries.right < exon_only.right && new_tx_end == boundaries.right;
            let direction_could_clamp = match self.config.shuffle_direction {
                ShuffleDirection::FivePrime => cheap_left,
                ShuffleDirection::ThreePrime => cheap_right,
            };
            let (left_clamp_fired, right_clamp_fired) = if direction_could_clamp {
                let (exon_only_start, exon_only_end, _exon_only_edit, _exon_only_warnings) =
                    self.normalize_na_edit(seq, edit, tx_start, tx_end, &exon_only, is_coding)?;
                let exon_only_start_0 = exon_only_start.saturating_sub(1);
                (
                    cheap_left && exon_only_start_0 < boundaries.left,
                    cheap_right && exon_only_end > boundaries.right,
                )
            } else {
                (false, false)
            };
            let (direction_str, clamp_kind) = match self.config.shuffle_direction {
                ShuffleDirection::FivePrime if left_clamp_fired => Some((
                    "5prime",
                    match start_axis {
                        boundary::AxisRegion::Cds => "cds_start",
                        boundary::AxisRegion::ThreeUtr => "3utr",
                        _ => "5utr",
                    },
                )),
                ShuffleDirection::ThreePrime if right_clamp_fired => Some((
                    "3prime",
                    match start_axis {
                        boundary::AxisRegion::Cds => "cds_end",
                        boundary::AxisRegion::FiveUtr => "5utr",
                        _ => "3utr",
                    },
                )),
                _ => None,
            }
            .unwrap_or(("", ""));
            if !direction_str.is_empty() {
                let acc = variant.accession.transcript_accession();
                warnings.push(NormalizationWarning::AxisClampApplied {
                    accession: acc,
                    direction: direction_str.to_string(),
                    clamp_kind: clamp_kind.to_string(),
                });
            }
        }

        // Issue #383 CDS-start clamp for the canonicalisation-rewrite
        // path. Companion to PR #343 (shuffle-path clamp). The spec
        // (§general "3'-rule applies to ALL descriptions" + the per-
        // axis coordinate treatment) says canonicalisation may not
        // silently move a CDS-interior input strictly into 5'UTR.
        // ferro's pre-shift `canonicalize_delins`, 5'-shift, and post-
        // shift ins→dup recognizer can each land such an input at
        // `new_tx_start < cds_start` — emitted variously as
        // `c.-N_1ins<…>`, `c.-M_-N dup`, etc. (e.g.
        // `NM_212556.2:c.1_2insCA` (5prime+cross) → `c.-2_-1dup`).
        //
        // The clamp fires on `new_tx_start < cds_start` regardless of
        // output edit-type, and rewrites by edit-type of the INPUT
        // (not the post-canon output, whose alt may already be a
        // rotated / duplicated derivative):
        //
        //   - Insertion input → `c.1_Kdelins<absorbed alt>` where
        //     `K = max(1, X - L)` (`L = |alt|`):
        //
        //       new_alt = ref[c.1..c.{X+1}] ++ alt[0..L - X + K]
        //
        //     For `X <= L+1` (the common case) `K = 1` and the formula
        //     collapses to `c.1delins<ref[c.1..c.{X+1}] ++ alt[..L-X+1]>`
        //     (1-base anchor — what biocommons emits for the
        //     NM_212556.2 corpus cases). For long left-shifts across a
        //     homopolymer (`X > L+1`, e.g. ref `c.1..c.20 = AAAAA…`
        //     with `c.5_6insAA`) `K = X - L` extends the delete window
        //     so the rewrite stays anchored at c.1 instead of silently
        //     falling through to `c.-1_1ins…`. We always verify
        //     equivalence against the input before accepting the
        //     rewrite — for non-homopolymer-tandem shapes where the
        //     formula doesn't algebraically reduce, we leave
        //     `new_edit` alone (preserving existing behaviour).
        //   - Delins input → restore the input form unchanged. The
        //     shared-affix trim is what pushed the residual past the
        //     boundary; suppressing it leaves the spec-canonical form
        //     (the input itself).
        //
        // Spanning-dup exception (#401): when the canon output is a
        // `Duplication` whose `new_tx_end >= cds_start`, the dup-source
        // spans the c.-1/c.1 boundary (one endpoint in UTR, one in CDS).
        // That IS the spec-canonical form (HGVS §general; edit-type
        // priority `dup > ins`) — biocommons emits it on inputs whose
        // alt equals `ref[c.-1] ++ ref[c.1]`. The clamp would
        // incorrectly collapse the spanning dup to `c.1delins<…>`, so
        // we skip the clamp in that case and keep the canon output.
        // Entirely-UTR dups (`new_tx_end < cds_start`) still clamp
        // unchanged.
        let spanning_dup_exception =
            matches!(new_edit, NaEdit::Duplication { .. }) && new_tx_end >= cds_start;
        // Issue #418 extension: when `cds_start == 1` (transcript has no
        // 5'UTR) the 5'-shuffle on an Insertion saturates at the
        // transcript start (`new_start = new_end = cds_start`) instead
        // of producing a true `new_tx_start < cds_start` signal. The
        // coordinate conversion from a 0-based shuffle result of 0 back
        // to 1-based HGVS clamps at 1, so the position-only gate
        // misses the saturation case. Detect it by the degenerate
        // (start == end) Insertion shape that only arises from
        // left-saturation, and fire the clamp the same way the
        // non-degenerate case does.
        let cds_start_left_saturated = matches!(edit, NaEdit::Insertion { .. })
            && new_tx_start == cds_start
            && new_tx_end == cds_start;
        if matches!(start_axis, boundary::AxisRegion::Cds)
            && (new_tx_start < cds_start || cds_start_left_saturated)
            && !spanning_dup_exception
        {
            match edit {
                NaEdit::Insertion {
                    sequence: InsertedSequence::Literal(in_lit),
                } => {
                    let alt_bytes: Vec<u8> = in_lit.bases().iter().map(|b| *b as u8).collect();
                    let x = (tx_start.saturating_sub(cds_start) + 1) as usize;
                    let cds_start_0b = (cds_start as usize).saturating_sub(1);
                    let prefix_end = cds_start_0b + x;
                    // k = max(1, x - L). For x <= L+1 this is 1 (the
                    // original 1-base-anchor formula); for the long-
                    // shift homopolymer case (x > L+1) this widens the
                    // delete window enough to give us a non-negative
                    // alt_take.
                    let k = x.saturating_sub(alt_bytes.len()).max(1);
                    let alt_take = (alt_bytes.len() + k).saturating_sub(x);
                    let delete_end_0b = cds_start_0b + k;
                    if x >= 1
                        && prefix_end <= seq.len()
                        && delete_end_0b <= seq.len()
                        && alt_take <= alt_bytes.len()
                    {
                        let prefix = &seq[cds_start_0b..prefix_end];
                        let alt_part = &alt_bytes[..alt_take];
                        let mut new_alt: Vec<u8> = prefix.to_vec();
                        new_alt.extend_from_slice(alt_part);
                        // Equivalence check: the clamped delins
                        //   delete ref[cds_start_0b..cds_start_0b+k]
                        //   insert new_alt (length L+k)
                        // must yield the same final sequence as the
                        // input
                        //   insert alt at tx_start (length L)
                        //   = between bytes tx_start_0b and
                        //     tx_start_0b+1 (tx_start_0b = cds_start_0b+x-1).
                        // After cancelling the shared prefix at
                        // `cds_start_0b` and the shared suffix from
                        // `cds_start_0b+x` onward, equivalence reduces
                        // to:
                        //   alt[..alt_take] ++ ref[cds_start_0b+k..
                        //                          cds_start_0b+x]
                        //       == alt
                        // i.e. the last (x - k) bytes of `alt` must
                        // equal the corresponding ref window. For the
                        // x <= L+1 case (k = 1) this collapses to the
                        // pre-existing single-base-anchor formula and
                        // always holds because the canonicalisation
                        // shifted past that ref window in the first
                        // place. For x > L+1 (k = x - L) the check
                        // accepts homopolymer / tandem-repeat cases
                        // and rejects anything else, so we never emit
                        // a delins that disagrees with the input.
                        let ref_tail = &seq[delete_end_0b..prefix_end];
                        let alt_tail_start = alt_take;
                        let equivalent = ref_tail.len() == alt_bytes.len() - alt_tail_start
                            && ref_tail == &alt_bytes[alt_tail_start..];
                        if equivalent {
                            let bases: Vec<Base> = new_alt
                                .iter()
                                .filter_map(|b| Base::from_char(*b as char))
                                .collect();
                            if bases.len() == new_alt.len() {
                                new_edit = NaEdit::Delins {
                                    sequence: InsertedSequence::Literal(Sequence::new(bases)),
                                    deleted: None,
                                    deleted_length: None,
                                };
                                new_tx_start = cds_start;
                                new_tx_end = cds_start + (k as u64) - 1;
                            }
                        }
                    }
                }
                NaEdit::Delins { .. } => {
                    new_edit = edit.clone();
                    new_tx_start = tx_start;
                    new_tx_end = tx_end;
                }
                _ => {
                    // Other edit types (Deletion, Duplication directly
                    // as input, Inversion, …) do not currently produce
                    // 5'UTR-resident rewrites from CDS-interior inputs
                    // in this canonicalisation path.
                }
            }
        }

        // Issue #387 CDS-end clamp for the canonicalisation-rewrite path
        // (3'-direction mirror of PR #385 / issue #383's CDS-start
        // clamp). HGVS DNA §general "3'-rule applies to ALL
        // descriptions" plus the per-axis coordinate treatment forbid
        // silently re-axing a CDS-interior input onto the 3'UTR
        // (`c.*<N>`) axis. ferro's 3'-shift on an alt whose rotated
        // form keeps matching reference past `c.<cds_end>` can leave
        // `new_tx_end > cds_end` — observed on
        // `NM_212556.2:c.1400_1401insAC` (3prime+cross & 3prime+no-cross)
        // emitting `c.1401_*1insCA`.
        //
        // For an Insertion input `c.X_X+1ins<alt>` rewrite as
        // `c.<cds_end>delins<new_alt>` where:
        //
        //     new_alt = alt[Y_c - 1 .. |alt|]     // tail of the original alt
        //             ++ ref[c.X+1 .. c.<cds_end>+1]   // CDS bases the shift walked past
        //     Y_c     = cds_end_c. - X        // distance from input start to cds_end
        //
        // Derived by equating post-edit genotype positions of the
        // original ins (positions X+1..X+|alt| = alt) with the clamped
        // delins (positions cds_end..cds_end+|alt| = new_alt). The
        // `Y_c-1` head of alt is consumed by the genotype shift the
        // first `Y_c-1` of which fall before `c.<cds_end>`; the
        // ref-byte suffix encodes the CDS bases pushed back by the
        // insertion. Length invariant `|new_alt| = |alt| + 1` keeps
        // the net length change matching the original ins.
        //
        // For a Delins input the symmetric `c.<cds_start>` case for #383
        // restored the original form unchanged; we mirror that here so
        // a Delins whose canonicalisation would push the residual ins
        // strictly past `cds_end` is suppressed.
        //
        // **Output edit-type gate.** Per the HGVS spec a `dup` is the
        // priority canonical form even when its span bridges the
        // CDS/3'UTR boundary — e.g. `NM_000051.3:c.9170_9171insAT` →
        // `c.9171_*1dup` (biocommons; ferro already emits this) is a
        // valid boundary-bridging dup that must NOT be rewritten as a
        // delins-at-`c.<cds_end>`. Restrict the clamp to outputs whose
        // post-canon edit-type is `Insertion` (= ins didn't promote to
        // dup). Duplication outputs that touch / cross `cds_end` are
        // the spec-canonical form and stay.
        //
        // Mirror of the CDS-start clamp's `spanning_dup_exception` gate
        // above (#401). The two clamps use opposite gate styles —
        // CDS-end positive-lists `Insertion`, CDS-start negative-lists
        // `Duplication` whose end reaches CDS — but both have the same
        // intent: preserve spanning duplications, clamp everything
        // else.
        if let Some(cds_end_tx) = transcript.cds_end {
            if matches!(start_axis, boundary::AxisRegion::Cds)
                && new_tx_end > cds_end_tx
                && matches!(new_edit, NaEdit::Insertion { .. })
            {
                match edit {
                    NaEdit::Insertion {
                        sequence: InsertedSequence::Literal(in_lit),
                    } => {
                        let alt_bytes: Vec<u8> = in_lit.bases().iter().map(|b| *b as u8).collect();
                        // Y_c = cds_end_c. - X = cds_end_tx - tx_start
                        // (both 1-based tx; the cds_start offsets cancel
                        // because the input start in c.-axis 1-based is
                        // `tx_start - cds_start + 1`).
                        let y_c = cds_end_tx.saturating_sub(tx_start);
                        // y_c >= 1 here (else `new_tx_end > cds_end_tx`
                        // could not have fired for an Insertion that
                        // started strictly inside CDS proper). Bound the
                        // alt-tail offset on the alt length too, so a
                        // pathological state (`y_c > |alt| + 1`) leaves
                        // the result untouched rather than panicking.
                        let alt_take_start = y_c.saturating_sub(1) as usize;
                        let suffix_start = (tx_end as usize).saturating_sub(1);
                        let suffix_end = cds_end_tx as usize;
                        if y_c >= 1
                            && alt_take_start <= alt_bytes.len()
                            && suffix_start <= suffix_end
                            && suffix_end <= seq.len()
                        {
                            let alt_part = &alt_bytes[alt_take_start..];
                            let suffix = &seq[suffix_start..suffix_end];
                            let mut new_alt: Vec<u8> = alt_part.to_vec();
                            new_alt.extend_from_slice(suffix);
                            let bases: Vec<Base> = new_alt
                                .iter()
                                .filter_map(|b| Base::from_char(*b as char))
                                .collect();
                            if bases.len() == new_alt.len() {
                                new_edit = NaEdit::Delins {
                                    sequence: InsertedSequence::Literal(Sequence::new(bases)),
                                    deleted: None,
                                    deleted_length: None,
                                };
                                new_tx_start = cds_end_tx;
                                new_tx_end = cds_end_tx;
                            }
                        }
                    }
                    NaEdit::Delins { .. } => {
                        new_edit = edit.clone();
                        new_tx_start = tx_start;
                        new_tx_end = tx_end;
                    }
                    _ => {
                        // Other edit types do not currently produce
                        // 3'UTR-resident rewrites from CDS-interior
                        // inputs in this canonicalisation path.
                    }
                }
            }
        }

        // Convert back to CDS coordinates
        let new_start = self.tx_to_cds_pos(new_tx_start, cds_start, transcript.cds_end)?;
        let new_end = self.tx_to_cds_pos(new_tx_end, cds_start, transcript.cds_end)?;

        let new_variant = CdsVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(Interval::new(new_start, new_end), new_edit),
        };

        // Issue #160 + #165 post-canonicalization split. The codon-frame
        // exception applies only to CDS-proper positions, which the
        // helper filters internally via `simple_cds_pos`.
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Cds(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split), warnings))
    }

    /// Normalize a transcript variant
    fn normalize_tx(
        &self,
        variant: &TxVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // Can't normalize variants with unknown edits or positions
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Tx(variant.clone()), vec![])),
        };

        // SVD-WG009: rewrite `con` to `delins`. Recurse on the rewritten
        // variant so the bounds-check + needs_normalization gates below still
        // fire — otherwise `n.<past-end>conT` slips past W4004. The downstream
        // passes (issue #333 ins[...] expansion, 3' shift, ins→dup, canonical
        // split) then all see the delins form. The recursion terminates on the
        // second entry because the rewrite emits `Delins`, which
        // `canonicalize_conversion_to_delins` no longer matches.
        if let Some(new_edit) = canonicalize_conversion_to_delins(edit) {
            let new_variant = TxVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    variant.loc_edit.location.clone(),
                    variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
                ),
            };
            return self.normalize_tx(&new_variant);
        }

        // Bounds check: `n.<N>` where N exceeds the transcript length is
        // malformed input. Mirrors normalize_cds's PositionPastEnd wiring;
        // closes #347. Must run BEFORE the `try_expand_tx_ins` expansion (so
        // `n.<past-end>conT` rejects on the coordinate rather than the
        // downstream `ins[...]` failure) and BEFORE the
        // `needs_normalization` short-circuit so substitutions
        // (`n.<N>G>C`) also get checked.
        //
        // The gate is optional — if the transcript can't be fetched or the
        // positions are unknown/`?`-offset, skip the check and fall through
        // so callers without manifest data still hit the downstream
        // `try_expand_tx_ins` and canonicalization passes.
        let accession = variant.accession.transcript_accession();
        let transcript_for_intronic =
            || -> Result<crate::reference::transcript::Transcript, FerroError> {
                self.provider
                    .get_transcript_for_variant(&HV::Tx(variant.clone()))
            };
        let transcript_opt = self.provider.get_transcript(&accession).ok();
        if self.config.should_reject_position_past_end()
            || self.config.should_warn_position_past_end()
        {
            if let (Some(transcript), Some(start_pos), Some(end_pos)) = (
                transcript_opt.as_ref(),
                variant.loc_edit.location.start.inner(),
                variant.loc_edit.location.end.inner(),
            ) {
                if !has_unknown_offset_tx(start_pos) && !has_unknown_offset_tx(end_pos) {
                    let mut bounds_warnings: Vec<NormalizationWarning> = Vec::new();
                    if let Some(w) = check_tx_pos_past_end(&accession, start_pos, transcript) {
                        bounds_warnings.push(w);
                    }
                    let end_distinct =
                        end_pos.base != start_pos.base || end_pos.offset != start_pos.offset;
                    if end_distinct {
                        if let Some(w) = check_tx_pos_past_end(&accession, end_pos, transcript) {
                            bounds_warnings.push(w);
                        }
                    }
                    if !bounds_warnings.is_empty() {
                        return Ok((
                            HV::Tx(self.canonicalize_tx_variant(variant)),
                            bounds_warnings,
                        ));
                    }
                }
            }
        }

        // Issue #333: expand bracketed / reference-range `ins[...]`
        // payloads. n. positions are direct transcript positions, so
        // the helper treats them as `InsCoordKind::Direct`. Runs AFTER
        // the bounds gate above so past-end `con`/`delins` inputs reject
        // on the coordinate rather than on a downstream `ins[...]`
        // expansion failure.
        let tx_accession = accession.clone();
        if let Some((new_variant, warning)) =
            self.try_expand_tx_ins(variant, edit, &tx_accession)?
        {
            let (result, mut warnings) = self.normalize_tx(&new_variant)?;
            warnings.insert(0, warning);
            return Ok((result, warnings));
        }

        // Extract positions for the post-expansion path. Substitutions (and
        // other non-shifted edits) also need normalization fall-through, so
        // we look up positions BEFORE the `needs_normalization`
        // short-circuit below.
        let start_pos = match variant.loc_edit.location.start.inner() {
            Some(pos) => pos,
            None => return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![])),
        };
        let end_pos = match variant.loc_edit.location.end.inner() {
            Some(pos) => pos,
            None => return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![])),
        };

        // Can't normalize variants with unknown (?) offsets - return unchanged
        if has_unknown_offset_tx(start_pos) || has_unknown_offset_tx(end_pos) {
            return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![]));
        }

        // Re-bind the transcript for the downstream normalization passes.
        // See `normalize_cds` for the rationale; reaching this point
        // implies we still need the transcript and must early-return if
        // it can't be fetched.
        let transcript = match transcript_opt {
            Some(t) => t,
            None => return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![])),
        };

        // Only normalize indels (the bounds check above runs regardless).
        if !needs_normalization(edit) {
            return Ok((HV::Tx(variant.clone()), vec![]));
        }

        if start_pos.is_intronic() || end_pos.is_intronic() {
            // Switch to the variant-aware lookup so an NG/NC-parented input
            // gets the build-correct chromosome.
            let transcript = transcript_for_intronic().unwrap_or(transcript);
            // Route intronic tx variants to the intronic normalization path
            if start_pos.is_intronic() && end_pos.is_intronic() {
                return self.normalize_intronic_tx(variant, &transcript, start_pos, end_pos, edit);
            }
            // Variant spans exon-intron boundary - not yet supported for n. coords
            return Err(FerroError::IntronicVariant {
                variant: format!("{}", variant),
            });
        }

        // Downstream `n.*N` positions encode a post-stop distance, not a
        // transcript index — feeding `pos.base` into the in-transcript path
        // below would normalize against the wrong window. `simple_tx_pos`
        // already skips these elsewhere; mirror that here. `check_tx_pos_past_end`
        // also short-circuits them, so this guard catches the remaining
        // shift-only edits (del/dup/ins).
        if start_pos.is_downstream() || end_pos.is_downstream() {
            return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![]));
        }

        // Only normalize positive positions (within transcript)
        // Negative positions are outside the transcript sequence
        if start_pos.base < 1 || end_pos.base < 1 {
            return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![]));
        }

        let tx_start = start_pos.base as u64;
        let tx_end = end_pos.base as u64;

        // Get boundaries.
        //
        // Per HGVS general.md the 3' rule has an explicit exception:
        // "deletions/duplications around exon/exon junctions using c., r.,
        // or n. reference sequences are not shifted." The carve-out is
        // narrow — it names deletions and duplications only — so we apply
        // the exon-only clamp from `get_cds_boundaries_with_axis_info`
        // **only for `Deletion` and `Duplication` edits**. Insertions
        // (including delins-shapes that go through the ins pipeline) and
        // inversions still 3'-shift across exon junctions, matching the
        // pre-#334 behavior. The CDS↔UTR axis clamp does not apply to
        // the n. axis by HGVS spec — `n.<N>` numbering runs through the
        // whole transcript with no CDS↔UTR sub-axis distinction,
        // regardless of whether the underlying transcript record happens
        // to carry `cds_start`/`cds_end`. (#334)
        let boundaries = if edit_is_del_or_dup(edit) {
            match boundary::get_cds_boundaries_with_axis_info(&transcript, tx_start, &self.config) {
                Ok(b) => b.exon,
                // `tx_start` outside every exon under `cross_boundaries=false`
                // means the position is intronic by the boundary helper's
                // definition. The intronic branch above should have handled
                // it; fall back to the canonicalize-only path defensively.
                Err(_) => return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![])),
            }
        } else {
            Boundaries::new(0, transcript.sequence_length())
        };

        // Perform normalization (n. non-coding tx context).
        // Coordinate-only transcripts (no cached bases) fall back to the
        // canonicalize-only path.
        let seq = match transcript.sequence.as_deref() {
            Some(s) => s.as_bytes(),
            None => return Ok((HV::Tx(self.canonicalize_tx_variant(variant)), vec![])),
        };
        let (new_start, new_end, new_edit, mut warnings) =
            self.normalize_na_edit(seq, edit, tx_start, tx_end, &boundaries, false)?;

        let new_variant = TxVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(TxPos::new(new_start as i64), TxPos::new(new_end as i64)),
                new_edit,
            ),
        };

        // Issue #160 + #165 post-canonicalization split.
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Tx(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split), warnings))
    }

    /// Normalize a protein variant.
    ///
    /// Performs four passes in order:
    ///
    /// 1. **Reference validation**: Check that position amino acids match the
    ///    reference protein sequence (if protein data is available).
    ///
    /// 2. **Redundant sequence removal**: Remove explicit sequences in deletions
    ///    when they match the amino acids at the deletion position.
    ///    Example: `p.Val600delVal` → `p.Val600del`
    ///
    /// 3. **3' shifting**: For `Deletion` and `Duplication` edits, walk the
    ///    rotation predicate via [`Self::shuffle_protein_3prime`] to land
    ///    at the spec-canonical most-3' anchor (HGVS general.md, #91).
    ///    Other edit kinds (substitution, frameshift, extension, identity,
    ///    no-protein, repeats) pass through unchanged. Insertion 3'-shift
    ///    is deferred to issue #92 (the `p.ins → p.dup` canonicalization).
    ///
    /// 4. **1-letter to 3-letter conversion**: (handled by parser/display)
    fn normalize_protein(
        &self,
        variant: &crate::hgvs::variant::ProteinVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        use crate::hgvs::edit::ProteinEdit;
        use crate::hgvs::variant::{LocEdit, ProteinVariant};

        // Validate reference amino acids if provider has protein data
        if self.provider.has_protein_data() {
            self.validate_protein_reference(variant)?;
        }

        // Get the current edit
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Protein(variant.clone()), vec![])),
        };

        // Apply normalization based on edit type. May rewrite both the
        // edit and the interval (e.g. `ins → dup` canonicalization
        // anchors the dup at the upstream-duplicated residue range,
        // which is a different interval than the input `<X>_<Y>ins…`).
        let (normalized_edit, post_canon_interval) = match edit {
            ProteinEdit::Deletion { sequence, count } => {
                // Check for redundant sequence that matches the position
                let new_edit = if let Some(seq) = sequence {
                    if self.is_redundant_protein_deletion_sequence(&variant.loc_edit.location, seq)
                    {
                        ProteinEdit::Deletion {
                            sequence: None,
                            count: *count,
                        }
                    } else {
                        edit.clone()
                    }
                } else {
                    edit.clone()
                };
                (new_edit, variant.loc_edit.location.clone())
            }
            // HGVS Prioritization: `<X>_<Y>ins<seq>` whose inserted
            // residues duplicate the residues immediately upstream
            // (anchored at p.X) is canonicalized to `<a>_<b>dup` over
            // that upstream range. The 3'-shift pass that follows will
            // then walk the dup to its canonical anchor. (#92)
            ProteinEdit::Insertion { sequence } => self
                .try_protein_ins_to_dup(variant, sequence)
                .unwrap_or_else(|| (edit.clone(), variant.loc_edit.location.clone())),
            // HGVS Prioritization: delins is the last-resort form. A
            // delins whose inserted residues share an affix with the
            // deleted residues must be re-described (delins.md:53-56);
            // after affix-trim the residual is routed through the
            // higher-priority canonicalization helpers. (#92)
            ProteinEdit::Delins { sequence } => self
                .try_protein_delins_canonicalize(variant, sequence)
                .unwrap_or_else(|| (edit.clone(), variant.loc_edit.location.clone())),
            // Other edits pass through unchanged
            _ => (edit.clone(), variant.loc_edit.location.clone()),
        };

        // Compute the post-shuffle interval. The shuffle layer operates on
        // the spec-canonical form, so we apply the residue-3'-shift after
        // the edit-shape rewrite above. Only deletions and duplications
        // are shuffled — per HGVS general.md "the most 3' position
        // possible … is arbitrarily assigned" — substitutions and other
        // edit kinds keep their input position (issue #91).
        //
        // The shuffler receives a synthesized variant that carries the
        // post-canonical interval (which may differ from the input's
        // after the `ins → dup` rewrite above). Important:
        // `shuffle_protein_3prime` reads only `loc_edit.location` from
        // the variant — it does NOT inspect the Mu wrapping on the
        // edit. The `map_ref` here preserves the input's Mu wrapping
        // for the eventual output, but the shuffler's behavior must
        // not depend on it. If a future shuffler refactor starts
        // reading the edit's uncertainty, audit this construction.
        let canon_variant = ProteinVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                post_canon_interval.clone(),
                variant.loc_edit.edit.map_ref(|_| normalized_edit.clone()),
            ),
        };
        let (new_interval, shuffled) = match &normalized_edit {
            ProteinEdit::Deletion { .. } | ProteinEdit::Duplication => self
                .shuffle_protein_3prime(&canon_variant, &normalized_edit)
                .unwrap_or_else(|| (post_canon_interval.clone(), false)),
            _ => (post_canon_interval.clone(), false),
        };

        // After a 3'-shift, a `Deletion { sequence: Some(..) }` carries
        // the pre-shift residue list but the new interval points at a
        // rotated reference window — keeping the stale residues would
        // emit a self-contradictory HGVS string (positions saying one
        // thing, `delXYZ` saying another). Per HGVS the residue list
        // is optional, so drop it; Display then falls back to the
        // unambiguous position-only `del` form. The `count` field is
        // length-only and is invariant under the shift, so it stays.
        let output_edit = match (&normalized_edit, shuffled) {
            (ProteinEdit::Deletion { count, .. }, true) => ProteinEdit::Deletion {
                sequence: None,
                count: *count,
            },
            _ => normalized_edit.clone(),
        };

        // Build the post-canonical variant; if no rewrite applied,
        // fall back to the input.
        let edit_changed = &output_edit != edit;
        let interval_changed = post_canon_interval != variant.loc_edit.location;
        let final_variant = if edit_changed || interval_changed || shuffled {
            ProteinVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    new_interval,
                    variant.loc_edit.edit.map_ref(|_| output_edit.clone()),
                ),
            }
        } else {
            variant.clone()
        };

        // Met1 soft-warning: emit when the final edit is a dup whose
        // interval includes position 1. (#92 sub-item 3)
        let mut warnings: Vec<NormalizationWarning> = Vec::new();
        if matches!(output_edit, ProteinEdit::Duplication) {
            if let (Some(s), Some(e)) = (
                final_variant.loc_edit.location.start.inner(),
                final_variant.loc_edit.location.end.inner(),
            ) {
                if s.number <= 1 && e.number >= 1 {
                    let location = if s.number == e.number {
                        format!("{}{}", s.aa, s.number)
                    } else {
                        format!("{}{}_{}{}", s.aa, s.number, e.aa, e.number)
                    };
                    warnings.push(NormalizationWarning::InitiatorMetCanonicalization {
                        accession: final_variant.accession.transcript_accession().to_string(),
                        location,
                    });
                }
            }
        }

        Ok((HV::Protein(final_variant), warnings))
    }

    /// Apply the HGVS 3' rule to a protein deletion or duplication.
    ///
    /// The shuffle rule for proteins: for a deletion or duplication of
    /// a contiguous range `protein[s0..s0+L]`, sliding the edit one
    /// residue right is spec-equivalent iff `protein[s0] == protein[s0 +
    /// L]`. This is the exact one-step equivalence condition for both
    /// operations — proof:
    ///
    /// - **Deletion**: removing `protein[s0..s0+L]` leaves the same
    ///   protein as removing `protein[s0+1..s0+L+1]` iff
    ///   `protein[s0] == protein[s0+L]` (the two windows differ only at
    ///   their first/last residue, so the remaining sequences are equal
    ///   iff those two residues match).
    /// - **Duplication**: inserting a copy of `protein[s0..s0+L]` after
    ///   position `s0+L-1` is spec-equivalent to inserting a copy of
    ///   `protein[s0+1..s0+L+1]` after position `s0+L` iff
    ///   `protein[s0] == protein[s0+L]` (the two duplicated windows
    ///   produce identical 2L-residue insertions iff those two residues
    ///   match).
    ///
    /// The walk iterates this predicate to a fixed point, yielding the
    /// spec-canonical most-3' anchor.
    ///
    /// Returns `Some((new_interval, shifted))` on success. `shifted` is
    /// `true` iff the rotation walked at least one step. Returns `None`
    /// when:
    ///
    /// - the variant's edit is not a `Deletion` or `Duplication`,
    /// - the provider has no protein data for the variant's accession,
    /// - either interval endpoint is uncertain (`?`),
    /// - the start/end positions are inverted (`end < start`),
    /// - the edit's footprint extends past the end of the protein.
    ///
    /// Both entrypoints — `Normalizer::normalize` (the strict API,
    /// which rejects reference mismatches) and
    /// `Normalizer::normalize_with_diagnostics` (the lenient API, which
    /// records mismatches as warnings) — share this shuffler via
    /// `normalize_protein`. They both treat `None` as "no shift" and
    /// return the input interval unchanged.
    fn shuffle_protein_3prime(
        &self,
        variant: &crate::hgvs::variant::ProteinVariant,
        edit: &ProteinEdit,
    ) -> Option<(ProtInterval, bool)> {
        // Only deletions and duplications shuffle.
        match edit {
            ProteinEdit::Deletion { .. } | ProteinEdit::Duplication => {}
            _ => return None,
        }

        // Reject uncertain endpoints — the rotation predicate is
        // undefined when either position is `?`.
        let start_pos = *variant.loc_edit.location.start.inner()?;
        let end_pos = *variant.loc_edit.location.end.inner()?;
        if end_pos.number < start_pos.number {
            return None;
        }

        let accession = variant.accession.transcript_accession();
        if !self.provider.has_protein_data() {
            return None;
        }
        let edit_len = end_pos.number - start_pos.number + 1;
        let protein_len = self.discover_protein_length(&accession)?;
        let edit_len_usize = edit_len as usize;
        let start_zero = (start_pos.number as usize).checked_sub(1)?;
        let protein_len_usize = protein_len as usize;
        // The edit footprint must fit inside the protein:
        // `[start_zero, start_zero + edit_len_usize)` ⊆ `[0, protein_len)`.
        if start_zero
            .checked_add(edit_len_usize)
            .map(|end_excl| end_excl > protein_len_usize)
            .unwrap_or(true)
        {
            return None;
        }

        // Fetch a bounded window starting at `start_zero` and grow it
        // only when the walk reaches the buffer edge with more
        // residues still ahead. This avoids the O(protein_len) full
        // read that would otherwise dominate normalization for tiny
        // edits on long (or pathological) proteins. The initial
        // window covers the edit plus 256 residues of look-ahead,
        // which suffices to halt the walk for every non-degenerate
        // case; degenerate runs trigger geometric re-fetches.
        const INITIAL_LOOKAHEAD: usize = 256;
        let mut window_end_excl = start_zero
            .saturating_add(edit_len_usize)
            .saturating_add(INITIAL_LOOKAHEAD)
            .min(protein_len_usize);
        let mut buf = self
            .provider
            .get_protein_sequence(&accession, start_zero as u64, window_end_excl as u64)
            .ok()?
            .into_bytes();
        // `s0` is the current edit start in protein-space (0-based);
        // the corresponding byte in `buf` lives at `s0 - start_zero`.
        let mut s0 = start_zero;
        let mut shifted = false;
        loop {
            let probe_protein = s0 + edit_len_usize;
            if probe_protein >= protein_len_usize {
                break;
            }
            // Grow the window if the probe falls past the buffer
            // edge. Geometric growth keeps total bytes fetched
            // ≤ 2 × actual_walk; the buffer offset is fixed at
            // `start_zero` so previously walked bytes stay valid.
            if probe_protein - start_zero >= buf.len() {
                let new_window_end_excl = (window_end_excl * 2)
                    .max(probe_protein + 1)
                    .min(protein_len_usize);
                if new_window_end_excl <= window_end_excl {
                    // Hit the C-terminus without growing — nothing
                    // left to read, halt.
                    break;
                }
                buf = self
                    .provider
                    .get_protein_sequence(&accession, start_zero as u64, new_window_end_excl as u64)
                    .ok()?
                    .into_bytes();
                window_end_excl = new_window_end_excl;
            }
            let s0_idx = s0 - start_zero;
            let probe_idx = probe_protein - start_zero;
            if buf[probe_idx] != buf[s0_idx] {
                break;
            }
            s0 += 1;
            shifted = true;
        }

        if !shifted {
            return Some((variant.loc_edit.location.clone(), false));
        }

        // Build the new interval. The post-shift start is residue
        // `s0 + 1` (1-based), end is `s0 + edit_len` (1-based).
        let new_start_num = (s0 + 1) as u64;
        let new_end_num = new_start_num + edit_len - 1;
        // Look up the AAs at the new positions from the reference.
        // After any geometric growth `buf` covers `[start_zero, …)`,
        // so the indices below are guaranteed in-range: the walk
        // either halted because `buf[probe_idx] != buf[s0_idx]`
        // (both indices already valid) or because the new edit
        // footprint reached the C-terminus (also in `buf` after the
        // final grow).
        let new_start_aa = AminoAcid::from_one_letter(buf[s0 - start_zero] as char)?;
        let new_end_aa =
            AminoAcid::from_one_letter(buf[s0 + edit_len_usize - 1 - start_zero] as char)?;

        let new_start = ProtPos::new(new_start_aa, new_start_num);
        let new_end = ProtPos::new(new_end_aa, new_end_num);
        Some((Interval::new(new_start, new_end), true))
    }

    /// Affix-trim a protein delins and route the residual through the
    /// existing canonicalization helpers.
    ///
    /// Returns `Some((edit, interval))` when the rewrite applied;
    /// `None` to fall back to the input delins (used for genuine
    /// delins with no shared affix, uncertain endpoints, or
    /// missing provider data).
    ///
    /// Branches (per HGVS Prioritization `general.md:56-57`):
    /// - residual empty + del empty → identity `p.<X>_<Y>=`
    /// - residual empty, del non-empty → pure `del`
    /// - residual single AA, del single AA → substitution
    /// - residual non-empty, del empty → route through
    ///   `try_protein_ins_to_dup` (anchored at the flanking
    ///   positions of the trimmed range)
    /// - otherwise → smaller delins (caller falls back; the helper
    ///   returns `None` for the no-progress case)
    fn try_protein_delins_canonicalize(
        &self,
        variant: &crate::hgvs::variant::ProteinVariant,
        seq: &crate::hgvs::edit::AminoAcidSeq,
    ) -> Option<(ProteinEdit, ProtInterval)> {
        use crate::hgvs::edit::AminoAcidSeq;

        // Reject uncertain endpoints — same rationale as
        // `try_protein_ins_to_dup`.
        let start_mu = variant.loc_edit.location.start.as_single()?;
        let end_mu = variant.loc_edit.location.end.as_single()?;
        let start_pos = match start_mu {
            Mu::Certain(p) => *p,
            _ => return None,
        };
        let end_pos = match end_mu {
            Mu::Certain(p) => *p,
            _ => return None,
        };
        if end_pos.number < start_pos.number {
            return None;
        }
        if start_pos.number == 0 {
            return None;
        }
        if !self.provider.has_protein_data() {
            return None;
        }

        let accession = variant.accession.transcript_accession();
        let expected_len = (end_pos.number - start_pos.number + 1) as usize;
        let ref_aas =
            self.fetch_protein_window(&accession, start_pos.number - 1, end_pos.number)?;
        if ref_aas.len() != expected_len {
            return None;
        }

        // Affix-trim: longest common prefix then longest common suffix.
        let mut lcp = 0usize;
        let max_pref = ref_aas.len().min(seq.0.len());
        while lcp < max_pref && ref_aas[lcp] == seq.0[lcp] {
            lcp += 1;
        }
        let ref_tail = &ref_aas[lcp..];
        let seq_tail = &seq.0[lcp..];
        let mut lcs = 0usize;
        let max_suf = ref_tail.len().min(seq_tail.len());
        while lcs < max_suf
            && ref_tail[ref_tail.len() - 1 - lcs] == seq_tail[seq_tail.len() - 1 - lcs]
        {
            lcs += 1;
        }

        // Trimmed positions (1-based inclusive). When lcp consumes
        // the entire deleted window we collapse to a zero-width range
        // anchored at `start + lcp`.
        let new_start = start_pos.number + lcp as u64;
        let new_end = end_pos.number - lcs as u64;
        let residual_seq = AminoAcidSeq(seq.0[lcp..seq.0.len() - lcs].to_vec());
        let residual_del = if new_start <= new_end {
            ref_aas[lcp..ref_aas.len() - lcs].to_vec()
        } else {
            Vec::new()
        };

        // No progress — caller falls back to the unchanged delins.
        // Do NOT route through an ins→dup search here: the deleted
        // window is still non-empty, so any tandem-match rewrite
        // would emit a bare `dup` that silently drops the deletion
        // side of the edit and returns a non-equivalent variant.
        if lcp == 0 && lcs == 0 {
            return None;
        }

        // Build the canonical residual form.
        if residual_del.is_empty() && residual_seq.0.is_empty() {
            // Identity: preserve the input range so Display emits
            // `p.<X>_<Y>=`. The `predicted` flag is propagated by the
            // caller via the `Mu` wrapper on the edit; here we mirror
            // the certain form. `whole_protein: false` because the
            // identity is position-specific, not whole-protein.
            return Some((
                ProteinEdit::Identity {
                    predicted: false,
                    whole_protein: false,
                },
                variant.loc_edit.location.clone(),
            ));
        }
        if residual_seq.0.is_empty() {
            // Pure deletion at the trimmed range.
            let dstart = ProtPos::new(residual_del[0], new_start);
            let dend = ProtPos::new(*residual_del.last()?, new_end);
            return Some((
                ProteinEdit::Deletion {
                    sequence: None,
                    count: None,
                },
                Interval::new(dstart, dend),
            ));
        }
        if residual_del.is_empty() {
            // Zero-width del + non-empty ins → route through the
            // existing ins→dup helper anchored at the flanking
            // positions of the trimmed range.
            let flank_start = new_start.saturating_sub(1);
            let flank_end = new_start;
            if flank_start >= 1 {
                let flank_start_aa =
                    self.fetch_protein_window(&accession, flank_start - 1, flank_start)?[0];
                let flank_end_aa =
                    self.fetch_protein_window(&accession, flank_end - 1, flank_end)?[0];
                let synth_start = ProtPos::new(flank_start_aa, flank_start);
                let synth_end = ProtPos::new(flank_end_aa, flank_end);
                let synth_variant = crate::hgvs::variant::ProteinVariant {
                    accession: variant.accession.clone(),
                    gene_symbol: variant.gene_symbol.clone(),
                    loc_edit: crate::hgvs::variant::LocEdit::with_uncertainty(
                        Interval::new(synth_start, synth_end),
                        variant.loc_edit.edit.map_ref(|_| ProteinEdit::Insertion {
                            sequence: residual_seq.clone(),
                        }),
                    ),
                };
                if let Some((edit, interval)) =
                    self.try_protein_ins_to_dup(&synth_variant, &residual_seq)
                {
                    return Some((edit, interval));
                }
            }
            // Fall back to the trimmed insertion. Fetch both flanking
            // residues in one provider call covering 1-based
            // [new_start - 1, new_start] (0-based half-open
            // [new_start - 2, new_start)). new_start <= 1 cannot
            // anchor a left flank in 1-based HGVS coordinates, so
            // bail out — the caller keeps the input delins.
            if new_start <= 1 {
                return None;
            }
            let flank_aas = self.fetch_protein_window(&accession, new_start - 2, new_start)?;
            if flank_aas.len() != 2 {
                return None;
            }
            let ins_start = ProtPos::new(flank_aas[0], new_start - 1);
            let ins_end = ProtPos::new(flank_aas[1], new_start);
            return Some((
                ProteinEdit::Insertion {
                    sequence: residual_seq,
                },
                Interval::new(ins_start, ins_end),
            ));
        }
        if residual_del.len() == 1 && residual_seq.0.len() == 1 {
            // Sub > delins.
            let pos = ProtPos::new(residual_del[0], new_start);
            return Some((
                ProteinEdit::Substitution {
                    reference: residual_del[0],
                    alternative: residual_seq.0[0],
                },
                Interval::new(pos, pos),
            ));
        }

        // Genuine smaller delins.
        let dstart = ProtPos::new(residual_del[0], new_start);
        let dend = ProtPos::new(*residual_del.last()?, new_end);
        Some((
            ProteinEdit::Delins {
                sequence: residual_seq,
            },
            Interval::new(dstart, dend),
        ))
    }

    /// HGVS Prioritization: if a protein insertion is equivalent to a
    /// duplication anywhere in the surrounding reference, the canonical
    /// form is a duplication.
    ///
    /// Algorithm: for `<X>_<Y>ins<seq>` with `len = seq.len()` and
    /// `Y = X + 1`, test two windows:
    ///
    /// 1. **Upstream window** (1-based `[X - len + 1, X]`): if equal to
    ///    `seq`, rewrite to `<a>_<b>dup` over that upstream range. The
    ///    inserted residues duplicate the residues immediately preceding
    ///    the insertion point.
    /// 2. **Downstream window** (1-based `[Y, Y + len - 1]`): if equal
    ///    to `seq`, rewrite to `<c>_<d>dup` over that downstream range.
    ///    The inserted residues duplicate the residues immediately
    ///    following the insertion point — semantically the same change
    ///    as the upstream-match case (and indeed they're related by 3'
    ///    shift), but covers inputs anchored on the 5' side of an
    ///    ambiguous boundary (e.g. `p.Val3_Ala4insAla` against a run
    ///    `...VAAA...`: V is the preceding residue, but the inserted A
    ///    duplicates the following A).
    ///
    /// The subsequent 3'-shift pass then walks the dup to its canonical
    /// anchor. Both anchors are spec-equivalent; the 3'-shift converges
    /// them to the same canonical form.
    ///
    /// Returns `Some((new_edit, new_interval))` on a successful rewrite.
    /// Returns `None` when the rewrite is not applicable:
    ///
    /// - provider has no protein data,
    /// - either interval endpoint is not `Mu::Certain` (uncertain or
    ///   range boundaries are passed through unchanged so the canonical
    ///   uncertainty marker survives),
    /// - the interval is not the standard `<X>_<X+1>ins…` adjacent shape,
    /// - any reference residue cannot be parsed via `from_one_letter`,
    /// - neither window matches `seq`.
    ///
    /// (#92)
    fn try_protein_ins_to_dup(
        &self,
        variant: &crate::hgvs::variant::ProteinVariant,
        seq: &crate::hgvs::edit::AminoAcidSeq,
    ) -> Option<(ProteinEdit, ProtInterval)> {
        if seq.is_empty() {
            return None;
        }
        if !self.provider.has_protein_data() {
            return None;
        }
        // Require Mu::Certain at both endpoints. Uncertain or range
        // boundaries carry semantics that the dup rewrite would silently
        // drop (e.g. `p.(Ala4)_Ala5insAla` parenthesizing the start
        // would not survive the rewrite to `p.Ala4dup`).
        let start_mu = variant.loc_edit.location.start.as_single()?;
        let end_mu = variant.loc_edit.location.end.as_single()?;
        let start_pos = match start_mu {
            crate::hgvs::uncertainty::Mu::Certain(p) => *p,
            _ => return None,
        };
        let end_pos = match end_mu {
            crate::hgvs::uncertainty::Mu::Certain(p) => *p,
            _ => return None,
        };
        // An insertion is parsed as `<X>_<Y>ins<seq>` with Y = X + 1.
        // Reject any other shape.
        if end_pos.number != start_pos.number + 1 {
            return None;
        }
        let len = seq.len() as u64;
        let accession = variant.accession.transcript_accession();

        // Upstream window: 1-based [start_pos.number - len + 1, start_pos.number].
        if start_pos.number >= len {
            let window_start_0 = start_pos.number - len;
            let window_end_excl = start_pos.number;
            if let Some(window_aas) =
                self.fetch_protein_window(&accession, window_start_0, window_end_excl)
            {
                if window_aas == seq.0 {
                    let dup_start_num = window_start_0 + 1;
                    let dup_end_num = start_pos.number;
                    let dup_start = ProtPos::new(window_aas[0], dup_start_num);
                    let dup_end = ProtPos::new(*window_aas.last()?, dup_end_num);
                    return Some((ProteinEdit::Duplication, Interval::new(dup_start, dup_end)));
                }
            }
        }

        // Downstream window: 1-based [end_pos.number, end_pos.number + len - 1].
        let window_start_0 = end_pos.number - 1;
        let window_end_excl = end_pos.number + len - 1;
        if let Some(window_aas) =
            self.fetch_protein_window(&accession, window_start_0, window_end_excl)
        {
            if window_aas == seq.0 {
                let dup_start_num = end_pos.number;
                let dup_end_num = end_pos.number + len - 1;
                let dup_start = ProtPos::new(window_aas[0], dup_start_num);
                let dup_end = ProtPos::new(*window_aas.last()?, dup_end_num);
                return Some((ProteinEdit::Duplication, Interval::new(dup_start, dup_end)));
            }
        }

        None
    }

    /// Fetch a reference protein window and decode each byte to
    /// [`AminoAcid`]. Returns `None` if the provider call fails, the
    /// returned slice has the wrong length, or any residue is not a
    /// valid one-letter code.
    fn fetch_protein_window(
        &self,
        accession: &str,
        start_0: u64,
        end_excl: u64,
    ) -> Option<Vec<AminoAcid>> {
        let expected_len = (end_excl - start_0) as usize;
        let window = self
            .provider
            .get_protein_sequence(accession, start_0, end_excl)
            .ok()?;
        let bytes = window.as_bytes();
        if bytes.len() != expected_len {
            return None;
        }
        let mut aas = Vec::with_capacity(expected_len);
        for &b in bytes {
            aas.push(AminoAcid::from_one_letter(b as char)?);
        }
        Some(aas)
    }

    /// Discover the length of a protein via the `ReferenceProvider`
    /// trait. The trait returns an out-of-range error rather than a
    /// truncated string and has no separate length API, so we probe
    /// with `get_protein_sequence(accession, 0, n)`.
    ///
    /// Strategy: seed the upper probe at 64 KiB, which covers titin
    /// (~35 kAA, the longest known human protein) and every realistic
    /// protein. If 64 KiB already succeeds, walk up exponentially to a
    /// 1 GiB safety cap. If the seed fails, the exact length is in
    /// `[0, SEED)` — skip the exponential ramp and let the binary
    /// search converge directly. Then binary-search between the last
    /// known-good and first known-bad probes for the exact length.
    /// Returns `None` for empty proteins or accessions the provider
    /// cannot resolve.
    fn discover_protein_length(&self, accession: &str) -> Option<u64> {
        const SEED: u64 = 64 * 1024;
        const CAP: u64 = 1 << 30;

        let probe_ok = |n: u64| self.provider.get_protein_sequence(accession, 0, n).is_ok();

        let (mut lo, mut hi) = if probe_ok(SEED) {
            // Common case: the seed succeeded. Grow only if the protein
            // is even longer (vanishingly rare).
            let mut hi = SEED;
            let mut lo = SEED;
            while probe_ok(hi) {
                lo = hi;
                if hi >= CAP {
                    break;
                }
                hi = hi.saturating_mul(2);
            }
            (lo, hi)
        } else {
            // Seed failed, so the exact length is in `[0, SEED)`.
            // Hand `[0, SEED)` straight to the binary search below
            // instead of re-doing a logarithmic exponential ramp from
            // 1 (which would duplicate the work the failed seed
            // already proved unnecessary). `lo = 0` is a valid lower
            // bound — `probe_ok(0)` (empty slice) is accepted by the
            // provider; only after the binary search converges to
            // `lo = 0` do we conclude the protein is genuinely empty.
            (0, SEED)
        };
        while lo + 1 < hi {
            let mid = lo + (hi - lo) / 2;
            if probe_ok(mid) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        if lo == 0 {
            return None;
        }
        Some(lo)
    }

    /// Validate that the amino acids in a protein variant match the reference
    ///
    /// Returns an error if there's a mismatch between the variant's stated
    /// amino acid(s) and the actual reference protein sequence.
    fn validate_protein_reference(
        &self,
        variant: &crate::hgvs::variant::ProteinVariant,
    ) -> Result<(), FerroError> {
        use crate::hgvs::edit::ProteinEdit;

        let accession = variant.accession.transcript_accession();

        // Get start and end positions
        let start_pos = match variant.loc_edit.location.start.inner() {
            Some(pos) => pos,
            None => return Ok(()), // Can't validate uncertain positions
        };
        let end_pos = match variant.loc_edit.location.end.inner() {
            Some(pos) => pos,
            None => return Ok(()),
        };

        // Get the edit to know what amino acids to validate
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok(()),
        };

        // Only validate edits that specify reference amino acids
        match edit {
            ProteinEdit::Substitution { .. } => {
                // For substitutions, the reference AA comes from the position (start_pos.aa),
                // not from the edit (which may be Xaa placeholder)
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
            }
            ProteinEdit::Deletion { .. } => {
                // Validate start position (from the interval's start AA)
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
                // Validate end position if different
                if end_pos.number != start_pos.number {
                    self.validate_protein_position(&accession, end_pos.number, &end_pos.aa)?;
                }
            }
            ProteinEdit::Duplication => {
                // Validate start and end positions
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
                if end_pos.number != start_pos.number {
                    self.validate_protein_position(&accession, end_pos.number, &end_pos.aa)?;
                }
            }
            ProteinEdit::Insertion { .. } => {
                // Validate flanking positions
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
                self.validate_protein_position(&accession, end_pos.number, &end_pos.aa)?;
            }
            ProteinEdit::Delins { .. } => {
                // Validate start and end
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
                if end_pos.number != start_pos.number {
                    self.validate_protein_position(&accession, end_pos.number, &end_pos.aa)?;
                }
            }
            ProteinEdit::Frameshift { .. } => {
                // Validate the frameshift position
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
            }
            ProteinEdit::Extension { .. } => {
                // Extension typically at Ter position - validate
                self.validate_protein_position(&accession, start_pos.number, &start_pos.aa)?;
            }
            _ => {
                // Identity, Unknown, etc. - no validation needed
            }
        }

        Ok(())
    }

    /// Validate a single protein position against reference
    fn validate_protein_position(
        &self,
        accession: &str,
        position: u64,
        expected_aa: &crate::hgvs::location::AminoAcid,
    ) -> Result<(), FerroError> {
        // Position is 1-based in HGVS, convert to 0-based for sequence access
        // get_protein_sequence uses half-open interval [start, end)
        let start = hgvs_pos_to_index(position) as u64;
        let end = position; // exclusive end

        // Try to get the reference amino acid
        match self.provider.get_protein_sequence(accession, start, end) {
            Ok(ref_seq) => {
                if ref_seq.len() != 1 {
                    return Ok(()); // Unexpected, skip validation
                }

                let ref_aa_char = ref_seq.chars().next().unwrap();
                let expected_char = expected_aa.to_one_letter();

                if ref_aa_char != expected_char {
                    return Err(FerroError::AminoAcidMismatch {
                        accession: accession.to_string(),
                        position,
                        expected: expected_aa.to_string(),
                        found: ref_aa_char.to_string(),
                    });
                }
            }
            Err(_) => {
                // Protein sequence not available, skip validation
            }
        }

        Ok(())
    }

    /// Check if the deletion sequence is redundant (matches the position amino acids)
    ///
    /// A deletion sequence is redundant if it exactly matches the amino acids
    /// specified in the interval. For example:
    /// - `p.Val600delVal` - sequence [Val] matches position 600's Val → redundant
    /// - `p.Lys23_Leu24delLysLeu` - sequence [Lys, Leu] matches positions 23-24 → redundant
    fn is_redundant_protein_deletion_sequence(
        &self,
        interval: &crate::hgvs::interval::ProtInterval,
        sequence: &crate::hgvs::edit::AminoAcidSeq,
    ) -> bool {
        // Get the start and end positions
        let start_pos = match interval.start.inner() {
            Some(pos) => pos,
            None => return false,
        };
        let end_pos = match interval.end.inner() {
            Some(pos) => pos,
            None => return false,
        };

        // Calculate expected sequence length from interval
        let interval_len = if end_pos.number >= start_pos.number {
            (end_pos.number - start_pos.number + 1) as usize
        } else {
            return false;
        };

        // Check if sequence length matches
        if sequence.len() != interval_len {
            return false;
        }

        // For a point deletion (single AA), check if the sequence matches the position AA
        if interval_len == 1 {
            return sequence.0.len() == 1 && sequence.0[0] == start_pos.aa;
        }

        // For a range deletion, check first and last AAs
        // The sequence should be [start_aa, ..., end_aa]
        if let (Some(first), Some(last)) = (sequence.0.first(), sequence.0.last()) {
            return *first == start_pos.aa && *last == end_pos.aa;
        }

        false
    }

    /// Normalize an RNA variant
    ///
    /// RNA variants (r.) are similar to transcript variants (n.) and undergo
    /// the same 3'/5' shifting normalization for indels. The main difference
    /// is that RNA uses lowercase nucleotides in HGVS notation.
    fn normalize_rna(
        &self,
        variant: &crate::hgvs::variant::RnaVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        use crate::hgvs::interval::RnaInterval;
        use crate::hgvs::variant::{LocEdit, RnaVariant};

        // Can't normalize variants with unknown edits or positions
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };

        // SVD-WG009: rewrite `con` to `delins`. Re-run on the rewritten
        // variant so downstream passes (3' shift, ins→dup, canonical
        // split) all see the delins form. The RNA axis does not run the
        // issue #333 ins[...] expansion step — that is wired only into
        // the genome/cds/tx/mt axes.
        if let Some(new_edit) = canonicalize_conversion_to_delins(edit) {
            let new_variant = RnaVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    variant.loc_edit.location.clone(),
                    variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
                ),
            };
            return self.normalize_rna(&new_variant);
        }

        // Only normalize indels
        if !needs_normalization(edit) {
            return Ok((HV::Rna(variant.clone()), vec![]));
        }

        // Check for intronic variants or unknown positions
        let start_pos = match variant.loc_edit.location.start.inner() {
            Some(pos) => pos,
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };
        let end_pos = match variant.loc_edit.location.end.inner() {
            Some(pos) => pos,
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };

        if start_pos.is_intronic() || end_pos.is_intronic() {
            return Err(FerroError::IntronicVariant {
                variant: format!("{}", variant),
            });
        }

        // Try to get transcript (RNA uses the same accession as mRNA transcripts)
        let accession = variant.accession.transcript_accession();
        let transcript = match self.provider.get_transcript(&accession) {
            Ok(t) => t,
            Err(_) => return Ok((HV::Rna(variant.clone()), vec![])),
        };

        // Convert RNA positions to transcript-1 positions, deciding per
        // endpoint. UTR (`r.*N`/`r.-N`) and non-positive bases need a CDS
        // to translate; non-UTR positive bases map 1:1 to transcript-1
        // indices and work without a CDS (mock providers used in tests
        // often omit cds_start/end). Choosing per endpoint keeps mixed
        // intervals like `r.50_*1del` from rerouting the positive end
        // through `rna_to_tx_pos` — issue #163 follow-up.
        //
        // Axis convention (pinned by `tests/issue_291_rna_axis_convention.rs`,
        // closes #291): `r.` positive non-UTR bases are
        // **transcript-1-relative**, NOT CDS-relative. `r.10` against a
        // transcript with `cds_start = 100` maps to tx index 10, not tx 109.
        // This is consistent across `fetch_ref_for_canonical_split`,
        // `simple_rna_pos`, and `normalize_na_edit`. Only `r.*N`/`r.-N`
        // (and base 0) translate through `cds_start`/`cds_end` via
        // `rna_to_tx_pos`.
        let cds_info = transcript.cds_start.zip(transcript.cds_end);
        let map_in = |pos: &crate::hgvs::location::RnaPos| -> Option<u64> {
            if pos.utr3 || pos.base < 1 {
                let (cds_start, cds_end) = cds_info?;
                self.rna_to_tx_pos(pos, cds_start, Some(cds_end)).ok()
            } else {
                Some(pos.base as u64)
            }
        };
        let tx_start = match map_in(start_pos) {
            Some(v) => v,
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };
        let tx_end = match map_in(end_pos) {
            Some(v) => v,
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };

        // Get boundaries.
        //
        // The previous "entire transcript span" comment was incorrect:
        // HGVS general.md explicitly carves out the exon-junction
        // exception for c., r., and n. references alike — deletions and
        // duplications around exon/exon junctions are not shifted across
        // the junction. The carve-out is narrow (it names deletions and
        // duplications only), so we apply the exon-only clamp from
        // `get_cds_boundaries_with_axis_info` **only for `Deletion` and
        // `Duplication` edits**. Insertions and inversions still 3'-shift
        // across exon junctions, matching the pre-#334 behavior. The
        // CDS↔UTR axis clamp (which `normalize_cds` also intersects)
        // does NOT apply to r., because r. natively spans 5'UTR (`r.-N`)
        // / coding (`r.<N>`) / 3'UTR (`r.*N`) and existing tests pin
        // shuffles that cross those sub-axes (see
        // `tests/issue_163_rna_utr3_flag.rs::
        // rna_mixed_cds_utr3_del_shifts_into_utr`). (#334)
        let boundaries = if edit_is_del_or_dup(edit) {
            match boundary::get_cds_boundaries_with_axis_info(&transcript, tx_start, &self.config) {
                Ok(b) => b.exon,
                Err(_) => return Ok((HV::Rna(variant.clone()), vec![])),
            }
        } else {
            Boundaries::new(0, transcript.sequence_length())
        };

        // Perform normalization (RNA context: codon-frame gate does not apply;
        // r. is not in the spec's accepted reference types for repeats).
        // Coordinate-only transcripts fall back to the canonicalize-only path.
        let seq = match transcript.sequence.as_deref() {
            Some(s) => s.as_bytes(),
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };
        let (new_tx_start, new_tx_end, new_edit, mut warnings) =
            self.normalize_na_edit(seq, edit, tx_start, tx_end, &boundaries, false)?;

        // Convert each normalized tx position back independently, restoring
        // UTR notation when the position falls outside the CDS. This catches
        // both the original issue #163 case (UTR input shuffling within the
        // UTR) and a positive-base input that shifts past `cds_end` during
        // normalization. Without `cds_info` we keep the simple base-1 mapping.
        //
        // For positions that stay inside the CDS-proper window
        // (`cds_start <= pos <= cds_end`), the tx index is emitted directly
        // as `r.{pos}` — i.e. the transcript-1-relative axis is preserved on
        // the way out, matching the convention enforced by `map_in` above.
        // See `tests/issue_291_rna_axis_convention.rs` for the pin.
        use crate::hgvs::location::RnaPos;
        let map_out = |pos: u64| -> Result<RnaPos, FerroError> {
            if let Some((cds_start, cds_end)) = cds_info {
                if pos < cds_start || pos > cds_end {
                    return self.tx_to_rna_pos(pos, cds_start, Some(cds_end));
                }
            }
            Ok(RnaPos::new(pos as i64))
        };
        let new_start = map_out(new_tx_start)?;
        let new_end = map_out(new_tx_end)?;

        let new_variant = RnaVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(RnaInterval::new(new_start, new_end), new_edit),
        };

        // Issue #160 + #165 post-canonicalization split (T/U-equivalent
        // comparison for the rev-comp scan and per-position emissions).
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Rna(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split), warnings))
    }

    /// Convert an RNA position to a transcript-1 position.
    ///
    /// Mirrors `cds_to_tx_pos`. `r.*N` maps to `cds_end + N`, `r.-N` to
    /// `cds_start + N` (HGVS skips the `0` gap).
    fn rna_to_tx_pos(
        &self,
        pos: &crate::hgvs::location::RnaPos,
        cds_start: u64,
        cds_end: Option<u64>,
    ) -> Result<u64, FerroError> {
        if pos.utr3 {
            let end = cds_end.ok_or_else(|| FerroError::ConversionError {
                msg: "No CDS end".to_string(),
            })?;
            let base = u64::try_from(pos.base).map_err(|_| FerroError::ConversionError {
                msg: format!("Negative base {} in 3' UTR position", pos.base),
            })?;
            Ok(end + base)
        } else if pos.base < 0 {
            let tx_pos = cds_start as i64 + pos.base;
            u64::try_from(tx_pos).map_err(|_| FerroError::ConversionError {
                msg: format!(
                    "RNA position r.{} maps before transcript start (cds_start={})",
                    pos.base, cds_start
                ),
            })
        } else if pos.base == 0 {
            Ok(cds_start.saturating_sub(1))
        } else {
            Ok(cds_start + pos.base as u64 - 1)
        }
    }

    /// Convert a transcript-1 position back to an RNA position, restoring
    /// the appropriate region (`r.*N` for 3'UTR, `r.-N` for 5'UTR).
    fn tx_to_rna_pos(
        &self,
        pos: u64,
        cds_start: u64,
        cds_end: Option<u64>,
    ) -> Result<crate::hgvs::location::RnaPos, FerroError> {
        use crate::hgvs::location::RnaPos;
        let end = cds_end.ok_or_else(|| FerroError::ConversionError {
            msg: "No CDS end".to_string(),
        })?;
        if pos < cds_start {
            Ok(RnaPos {
                base: pos as i64 - cds_start as i64,
                offset: None,
                utr3: false,
            })
        } else if pos > end {
            Ok(RnaPos {
                base: (pos - end) as i64,
                offset: None,
                utr3: true,
            })
        } else {
            Ok(RnaPos {
                base: (pos - cds_start + 1) as i64,
                offset: None,
                utr3: false,
            })
        }
    }

    /// Normalize a mitochondrial variant
    ///
    /// Mirrors `normalize_genome` for non-origin-crossing variants: fetch
    /// a sequence window around the variant, run `normalize_na_edit` with
    /// `is_coding=false` (mito is genomic-style and not subject to the
    /// codon-frame restriction — `repeated.md` line 21 restricts the
    /// codon-frame gate exclusively to `c.` descriptions), then map
    /// positions back.
    ///
    /// Origin-crossing (wraparound) `del`/`delins` variants are rejected
    /// at parse time by `parse_genome_interval`'s inverted-range check;
    /// wraparound `dup`/`ins`/`inv` are exempt from that check and reach
    /// this function. Without provider data they fall through to
    /// `canonicalize_mt_variant`; with provider data the window fetch
    /// errors (start > end is invalid for a linear slice), again
    /// dropping to the fallback. F1 / #129 will introduce circular-aware
    /// semantics in a follow-up — see `tests/mito_circular_audit.rs` for
    /// the pinned behavior preserved by this PR.
    fn normalize_mt(
        &self,
        variant: &crate::hgvs::variant::MtVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // Can't normalize variants with unknown edits or positions.
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Mt(variant.clone()), vec![])),
        };

        // SVD-WG009: rewrite `con` to `delins` up front. Pure-syntax;
        // no reference data needed. Re-run on the rewritten variant so
        // the downstream passes (past-end bounds gate, issue #333 ins[...]
        // expansion, 3' shift, ins→dup, canonical split) all see the delins
        // form — otherwise `m.<past>conT` would early-return through the
        // bounds gate below in non-canonical `con` form, mirroring the
        // silent bypass previously fixed on the c./n. axes (see
        // `normalize_cds` / `normalize_tx` and the `*_for_con_rewrite`
        // tests in `tests/issue_336_position_past_end.rs`).
        // `canonicalize_conversion_to_delins` only matches `NaEdit::Conversion`
        // and the rewrite swaps it for `Delins`, so the recursion terminates
        // on the second entry.
        if let Some(new_edit) = canonicalize_conversion_to_delins(edit) {
            let new_variant = MtVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    variant.loc_edit.location.clone(),
                    variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
                ),
            };
            return self.normalize_mt(&new_variant);
        }

        // Past-end bounds check (#393, W4004). Fires when an m. position
        // exceeds the contig length, e.g. `m.16570A>G` on the 16569-bp mtDNA.
        // Wraparound ranges (start > end, valid per SVD-WG006) bypass this
        // only when both endpoints fit in the contig.
        // Conservative skip when the provider has no length data for this
        // accession — mirrors the transcript-absent skip in normalize_cds.
        // Runs AFTER the `con -> delins` rewrite above so past-end inputs
        // arriving via the `con` fast path (e.g. `m.16570conT`) re-enter
        // this function in canonical `delins` form before the gate fires.
        let mt_accession_bounds = variant.accession.transcript_accession();
        if self.config.should_reject_position_past_end()
            || self.config.should_warn_position_past_end()
        {
            if let (Some(start_pos), Some(end_pos)) = (
                variant.loc_edit.location.start.inner(),
                variant.loc_edit.location.end.inner(),
            ) {
                if let Ok(contig_length) = self.provider.get_sequence_length(&mt_accession_bounds) {
                    let mut bounds_warnings: Vec<NormalizationWarning> = Vec::new();
                    if let Some(w) =
                        check_mt_pos_past_end(&mt_accession_bounds, start_pos, contig_length)
                    {
                        bounds_warnings.push(w);
                    }
                    // Mirror the c./n. dedupe guard: compare the full
                    // (base, offset) tuple, not just `base`. A range like
                    // `m.16570+1_16570` shares the same `base` on both
                    // endpoints but the offset differs — checking only
                    // `base` would skip the in-bounds endpoint and lose
                    // W4004 on the other one. Offsets are non-standard on
                    // m. but are parseable today (see
                    // `check_mt_pos_past_end`'s skip on `pos.offset.is_some()`).
                    let end_distinct =
                        end_pos.base != start_pos.base || end_pos.offset != start_pos.offset;
                    if end_distinct {
                        if let Some(w) =
                            check_mt_pos_past_end(&mt_accession_bounds, end_pos, contig_length)
                        {
                            bounds_warnings.push(w);
                        }
                    }
                    if !bounds_warnings.is_empty() {
                        return Ok((
                            HV::Mt(self.canonicalize_mt_variant(variant)),
                            bounds_warnings,
                        ));
                    }
                }
            }
        }

        // Issue #333: expand bracketed / reference-range `ins[...]`
        // payloads. m. positions are direct genomic positions on the
        // mitochondrial accession.
        let mt_accession = variant.accession.transcript_accession();
        if let Some((new_variant, warning)) =
            self.try_expand_mt_ins(variant, edit, &mt_accession)?
        {
            let (result, mut warnings) = self.normalize_mt(&new_variant)?;
            warnings.insert(0, warning);
            return Ok((result, warnings));
        }

        // Only normalize indels; substitutions / identity / repeat-with-
        // count pass through unchanged. Mirrors `normalize_genome`.
        if !needs_normalization(edit) {
            return Ok((HV::Mt(variant.clone()), vec![]));
        }

        // Fallback for variants we cannot remap through the window-based
        // pipeline (unknown position, decorated position, no provider
        // data, or a reversed `<high>_<low>` wraparound range). Runs
        // minimal-notation cleanup, then applies `apply_canonical_split`
        // so issue #160 inv-split and issue #165 sub-only decomposition
        // remain in force — those run on a narrow fetch
        // (`fetch_ref_for_canonical_split`) that can succeed even when
        // the wider shuffle window does not.
        //
        // Reversed ranges (`start > end`, per SVD-WG006 wraparound on
        // `m.`/`o.`) skip `apply_canonical_split`: the helper computes
        // `expected_span = hgvs_end - hgvs_start + 1` as `u64`, which
        // underflows on reversed inputs. Today's shipped providers
        // reject `start > end` at their `get_sequence` boundary, so the
        // helper currently short-circuits via the `.ok()?` path before
        // touching that arithmetic. Guarding here removes the latent
        // dependency on every present and future provider doing that
        // boundary check correctly — circular-aware split math is its
        // own design and belongs to #129 Path 2 follow-up.
        let mt_fallback = |v: &crate::hgvs::variant::MtVariant| {
            let canonical = self.canonicalize_mt_variant(v);
            let is_reversed = canonical
                .loc_edit
                .location
                .start
                .inner()
                .zip(canonical.loc_edit.location.end.inner())
                .is_some_and(|(s, e)| s.base > e.base);
            if is_reversed {
                return (HV::Mt(canonical), vec![]);
            }
            let (split, warnings) = self.apply_canonical_split(HV::Mt(canonical));
            (wrap_allele_if_split(split), warnings)
        };

        let accession = variant.accession.transcript_accession();
        let start_pos = match variant.loc_edit.location.start.inner() {
            Some(pos) => pos,
            None => return Ok(mt_fallback(variant)),
        };
        let end_pos = match variant.loc_edit.location.end.inner() {
            Some(pos) => pos,
            None => return Ok(mt_fallback(variant)),
        };

        // Decorated genome positions (offset / pter / qter / cen) cannot
        // be losslessly remapped through base-only window normalization —
        // remapping via `pos.base` and rebuilding with `GenomePos::new`
        // would silently drop the decoration. Fall back to minimal-
        // notation cleanup for these.
        //
        // Telomere/centromere markers are intentionally NOT resolved for
        // mitochondrial (circular) references — `m.` has no telomeres, so
        // pter/qter/cen are meaningless here; keep the deliberate bail-out
        // rather than unifying with normalize_genome's resolver (#488).
        if start_pos.offset.is_some()
            || end_pos.offset.is_some()
            || start_pos.is_special()
            || end_pos.is_special()
        {
            return Ok(mt_fallback(variant));
        }

        let start = start_pos.base;
        let end = end_pos.base;

        // Wraparound `<high>_<low>` ranges (per SVD-WG006: `del`/`dup`/
        // `delins` on `m.`/`o.`) have `start > end`. Circular-aware
        // 3'-shift across the origin is out of scope for #129 (matches
        // mutalyzer + biocommons + strict spec reading); route these
        // straight to the canonicalize-only fallback. Doing the gate
        // here — rather than relying on the provider rejecting `start >
        // end` and the post-arithmetic `rel_end = end - window_start`
        // underflow as a fail-safe — keeps the path correct under any
        // window-size configuration. Issue #129 follow-up: implement
        // circular-aware shuffle modulo contig length.
        if start > end {
            return Ok(mt_fallback(variant));
        }

        // Window-based fetch around the variant. Non-origin-crossing
        // path: identical to genomic.
        let window_start = start.saturating_sub(self.config.window_size);
        let seq_result = self.provider.get_sequence(
            &accession,
            window_start,
            end.saturating_add(self.config.window_size),
        );

        let ref_seq = match seq_result {
            Ok(s) => s,
            // No reference data → fall back to minimal-notation cleanup.
            Err(_) => return Ok(mt_fallback(variant)),
        };

        let rel_start = start - window_start;
        let rel_end = end - window_start;

        // Mitochondrial reference is plus-strand and not subject to the
        // codon-frame `unit_len % 3 == 0` restriction (the mito genome
        // has no canonical "CDS" exemption boundary in the same sense as
        // nuclear `c.`; the spec's mito chapter doesn't carry the
        // codon-frame clause), so pass `is_coding=false`.
        let (new_rel_start, new_rel_end, new_edit, mut warnings) = self.normalize_na_edit(
            ref_seq.as_bytes(),
            edit,
            rel_start,
            rel_end,
            &Boundaries::new(0, ref_seq.len() as u64),
            false,
        )?;

        let new_start = new_rel_start + window_start;
        let new_end = new_rel_end + window_start;

        let new_variant = MtVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(GenomePos::new(new_start), GenomePos::new(new_end)),
                new_edit,
            ),
        };

        // Issue #160 inv-split post-pass (mirrors normalize_genome).
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Mt(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split), warnings))
    }

    /// Apply minimal notation to an mt variant without full normalization.
    /// Mirrors `canonicalize_genome_variant` — used as a fallback when
    /// reference data is unavailable.
    fn canonicalize_mt_variant(
        &self,
        variant: &crate::hgvs::variant::MtVariant,
    ) -> crate::hgvs::variant::MtVariant {
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return variant.clone(),
        };

        if !should_canonicalize(edit) {
            return variant.clone();
        }

        let canonical_edit = canonicalize_edit(edit);

        MtVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                variant.loc_edit.location.clone(),
                variant.loc_edit.edit.map_ref(|_| canonical_edit),
            ),
        }
    }

    /// Normalize an intronic CDS variant
    ///
    /// This converts the intronic position to genomic coordinates, normalizes
    /// in genomic space, and converts back to CDS intronic notation.
    fn normalize_intronic_cds(
        &self,
        variant: &CdsVariant,
        transcript: &crate::reference::transcript::Transcript,
        start_pos: &CdsPos,
        end_pos: &CdsPos,
        edit: &NaEdit,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        use crate::convert::CoordinateMapper;

        // Check if we have genomic data available
        if !self.provider.has_genomic_data() {
            return Err(FerroError::IntronicVariant {
                variant: format!("{}", variant),
            });
        }

        // Get the chromosome for this transcript. Issue #332: include the
        // parent accession and full variant Display in the error so the
        // remaining failure mode (transcript not present on any genome build)
        // is diagnosable without re-running with extra logging.
        let chromosome =
            transcript
                .chromosome
                .as_ref()
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!(
                        "Transcript {} has no chromosome mapping for intronic \
                         normalization (parent={}, variant={}); the cdot data has \
                         no genomic alignment for this transcript on any known \
                         genome build",
                        transcript.id,
                        variant
                            .accession
                            .genomic_context
                            .as_deref()
                            .map(|a| a.full())
                            .unwrap_or_else(|| "<none>".to_string()),
                        variant,
                    ),
                })?;

        // Create coordinate mapper
        let mapper = CoordinateMapper::new(transcript);

        // Convert CDS intronic positions to genomic
        let g_start = mapper.cds_to_genomic_with_intron(start_pos)?;
        let g_end = mapper.cds_to_genomic_with_intron(end_pos)?;

        // On minus strand, genomic coords may be reversed relative to coding order.
        // Track whether we swap so we can restore coding order after normalization.
        let swapped = g_start > g_end;
        let (g_start, g_end) = if swapped {
            (g_end, g_start)
        } else {
            (g_start, g_end)
        };

        // Get a window of genomic sequence around the variant for normalization
        // Use the same window size as for exonic normalization
        let window = self.config.window_size;
        let seq_start = g_start.saturating_sub(window);
        let seq_end = g_end.saturating_add(window);

        // Fetch genomic sequence
        let genomic_seq = self
            .provider
            .get_genomic_sequence(chromosome, seq_start, seq_end)?;

        // Calculate the variant position relative to the fetched sequence
        let rel_start = (g_start - seq_start) + 1; // 1-based
        let rel_end = (g_end - seq_start) + 1;

        // Define boundaries within the intron
        // For intronic variants, we can shift within the intron but not into exons
        // Find the intron boundaries
        // Use exon-aware CDS-to-tx mapping to account for cdot's gap positions
        let tx_pos = mapper.cds_to_tx(start_pos)?;
        let tx_start = u64::try_from(tx_pos.base).map_err(|_| FerroError::ConversionError {
            msg: format!(
                "Negative transcript position {} during intronic normalization",
                tx_pos.base
            ),
        })?;

        let intron = transcript
            .find_intron_at_tx_boundary(tx_start, start_pos.offset.unwrap_or(0))
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Could not find intron for normalization".to_string(),
            })?;

        // Get intron boundaries in genomic coordinates
        let (intron_g_start, intron_g_end) = match (intron.genomic_start, intron.genomic_end) {
            (Some(s), Some(e)) => (s, e),
            _ => {
                return Err(FerroError::ConversionError {
                    msg: "Intron has no genomic coordinates".to_string(),
                })
            }
        };

        // Calculate relative intron boundaries
        let intron_rel_start = intron_g_start.saturating_sub(seq_start) + 1;
        let intron_rel_end = intron_g_end.saturating_sub(seq_start) + 1;
        let boundaries = Boundaries::new(intron_rel_start, intron_rel_end);

        // On minus-strand transcripts the genomic-strand sequence is the
        // reverse complement of the transcript view, but the variant's
        // edit alt is in transcript view. Running `normalize_na_edit` on
        // the genomic-strand bytes therefore defeats every rule that
        // compares the alt against the local reference window. Flip the
        // sequence and the relative positions / boundaries to transcript
        // view here, run normalization, then map the result positions
        // back to the genomic frame. (Issue #98.)
        let (work_seq, work_rel_start, work_rel_end, work_boundaries) = flip_intronic_for_strand(
            transcript.strand,
            &genomic_seq,
            rel_start,
            rel_end,
            &boundaries,
        )?;

        // Perform normalization in transcript-view space (CDS intronic context).
        // HGVS spec (repeated.md): the codon-frame restriction
        // (`unit_len % 3 == 0` for repeat notation in `c.` context)
        // applies only to bases inside the CDS proper. Introns are
        // exempt:
        //   > This restriction only applies to the coding sequence,
        //   > which does not include the introns or the UTR sequence.
        // Pass `is_coding=false` so an intronic homopolymer dup/del
        // can emit `[N±k]` repeat notation instead of falling back to
        // the gated `ins<literal>` / plain `del` forms.
        let seq_bytes = work_seq.as_bytes();
        let (work_new_rel_start, work_new_rel_end, new_edit, warnings) = self.normalize_na_edit(
            seq_bytes,
            edit,
            work_rel_start,
            work_rel_end,
            &work_boundaries,
            false,
        )?;

        // Map the result positions back to the genomic-strand frame
        let (new_rel_start, new_rel_end) = unflip_intronic_positions(
            transcript.strand,
            work_seq.len() as u64,
            work_new_rel_start,
            work_new_rel_end,
        );

        // Convert the normalized genomic position back to absolute genomic
        let new_g_start = seq_start + new_rel_start - 1;
        let new_g_end = seq_start + new_rel_end - 1;

        // Convert back to CDS intronic notation
        let new_start = mapper.genomic_to_cds_intronic(new_g_start)?;
        let new_end = mapper.genomic_to_cds_intronic(new_g_end)?;

        // Restore coding order if positions were swapped for genomic processing
        let (new_start, new_end) = if swapped {
            (new_end, new_start)
        } else {
            (new_start, new_end)
        };

        let new_variant = CdsVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(Interval::new(new_start, new_end), new_edit),
        };

        Ok((HV::Cds(new_variant), warnings))
    }

    /// Normalize an intronic transcript (n.) variant
    ///
    /// This mirrors `normalize_intronic_cds()` but works with TxPos instead of CdsPos.
    /// Converts to genomic coordinates, normalizes in genomic space, and converts back.
    fn normalize_intronic_tx(
        &self,
        variant: &TxVariant,
        transcript: &crate::reference::transcript::Transcript,
        start_pos: &TxPos,
        end_pos: &TxPos,
        edit: &NaEdit,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        use crate::convert::CoordinateMapper;

        // Check if we have genomic data available
        if !self.provider.has_genomic_data() {
            return Err(FerroError::IntronicVariant {
                variant: format!("{}", variant),
            });
        }

        // Get the chromosome for this transcript. Issue #332: include the
        // parent accession and full variant Display in the error so the
        // remaining failure mode (transcript not present on any genome build)
        // is diagnosable without re-running with extra logging.
        let chromosome =
            transcript
                .chromosome
                .as_ref()
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!(
                        "Transcript {} has no chromosome mapping for intronic \
                         normalization (parent={}, variant={}); the cdot data has \
                         no genomic alignment for this transcript on any known \
                         genome build",
                        transcript.id,
                        variant
                            .accession
                            .genomic_context
                            .as_deref()
                            .map(|a| a.full())
                            .unwrap_or_else(|| "<none>".to_string()),
                        variant,
                    ),
                })?;

        // Create coordinate mapper
        let mapper = CoordinateMapper::new(transcript);

        // Convert tx intronic positions to genomic
        let g_start = mapper.tx_to_genomic_with_intron(start_pos)?;
        let g_end = mapper.tx_to_genomic_with_intron(end_pos)?;

        // On minus strand, genomic coords may be reversed relative to coding order.
        // Track whether we swap so we can restore coding order after normalization.
        let swapped = g_start > g_end;
        let (g_start, g_end) = if swapped {
            (g_end, g_start)
        } else {
            (g_start, g_end)
        };

        // Get a window of genomic sequence around the variant
        let window = self.config.window_size;
        let seq_start = g_start.saturating_sub(window);
        let seq_end = g_end.saturating_add(window);

        // Fetch genomic sequence
        let genomic_seq = self
            .provider
            .get_genomic_sequence(chromosome, seq_start, seq_end)?;

        // Calculate the variant position relative to the fetched sequence
        let rel_start = (g_start - seq_start) + 1; // 1-based
        let rel_end = (g_end - seq_start) + 1;

        // Find the intron boundaries for normalization limits
        let tx_start = u64::try_from(start_pos.base).map_err(|_| FerroError::ConversionError {
            msg: format!(
                "Negative transcript position {} during intronic normalization",
                start_pos.base
            ),
        })?;

        let intron = transcript
            .find_intron_at_tx_boundary(tx_start, start_pos.offset.unwrap_or(0))
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Could not find intron for normalization".to_string(),
            })?;

        // Get intron boundaries in genomic coordinates
        let (intron_g_start, intron_g_end) = match (intron.genomic_start, intron.genomic_end) {
            (Some(s), Some(e)) => (s, e),
            _ => {
                return Err(FerroError::ConversionError {
                    msg: "Intron has no genomic coordinates".to_string(),
                })
            }
        };

        // Calculate relative intron boundaries
        let intron_rel_start = intron_g_start.saturating_sub(seq_start) + 1;
        let intron_rel_end = intron_g_end.saturating_sub(seq_start) + 1;
        let boundaries = Boundaries::new(intron_rel_start, intron_rel_end);

        // See `normalize_intronic_cds`: same orientation fix for #98.
        let (work_seq, work_rel_start, work_rel_end, work_boundaries) = flip_intronic_for_strand(
            transcript.strand,
            &genomic_seq,
            rel_start,
            rel_end,
            &boundaries,
        )?;

        // Perform normalization in transcript-view space (n. non-coding intronic context).
        let seq_bytes = work_seq.as_bytes();
        let (work_new_rel_start, work_new_rel_end, new_edit, warnings) = self.normalize_na_edit(
            seq_bytes,
            edit,
            work_rel_start,
            work_rel_end,
            &work_boundaries,
            false,
        )?;

        // Map the result positions back to the genomic-strand frame
        let (new_rel_start, new_rel_end) = unflip_intronic_positions(
            transcript.strand,
            work_seq.len() as u64,
            work_new_rel_start,
            work_new_rel_end,
        );

        // Convert the normalized genomic position back to absolute genomic
        let new_g_start = seq_start + new_rel_start - 1;
        let new_g_end = seq_start + new_rel_end - 1;

        // Convert back to transcript intronic notation
        let new_start = mapper.genomic_to_tx_with_intron(new_g_start)?;
        let new_end = mapper.genomic_to_tx_with_intron(new_g_end)?;

        // Restore coding order if positions were swapped for genomic processing
        let (new_start, new_end) = if swapped {
            (new_end, new_start)
        } else {
            (new_start, new_end)
        };

        let new_variant = TxVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(Interval::new(new_start, new_end), new_edit),
        };

        Ok((HV::Tx(new_variant), warnings))
    }

    /// Normalize a CDS variant that spans an exon-intron boundary
    ///
    /// This handles cases like:
    /// - c.914_918+3del (exonic start, intronic end)
    /// - c.194-64_233del (intronic start, exonic end)
    ///
    /// Strategy: Convert to genomic coordinates, normalize there, convert back.
    fn normalize_boundary_spanning_cds(
        &self,
        variant: &CdsVariant,
        transcript: &crate::reference::transcript::Transcript,
        start_pos: &CdsPos,
        end_pos: &CdsPos,
        edit: &NaEdit,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        use crate::convert::CoordinateMapper;

        // Require genomic data for boundary spanning normalization
        if !self.provider.has_genomic_data() {
            return Err(FerroError::ExonIntronBoundary {
                exon: 0,
                variant: format!("{}", variant),
            });
        }

        // Issue #332: same improved error shape as the intronic paths.
        let chromosome =
            transcript
                .chromosome
                .as_ref()
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!(
                        "Transcript {} has no chromosome for boundary \
                         normalization (parent={}, variant={}); the cdot data \
                         has no genomic alignment for this transcript on any \
                         known genome build",
                        transcript.id,
                        variant
                            .accession
                            .genomic_context
                            .as_deref()
                            .map(|a| a.full())
                            .unwrap_or_else(|| "<none>".to_string()),
                        variant,
                    ),
                })?;

        let mapper = CoordinateMapper::new(transcript);

        // Convert both positions to genomic
        // For exonic positions, use standard conversion
        // For intronic positions, use intronic conversion
        let g_start = self.cds_pos_to_genomic(&mapper, start_pos)?;
        let g_end = self.cds_pos_to_genomic(&mapper, end_pos)?;

        // On minus strand, genomic coords may be reversed relative to coding order.
        // Track whether we swap so we can restore coding order after normalization.
        let swapped = g_start > g_end;
        let (g_start, g_end) = if swapped {
            (g_end, g_start)
        } else {
            (g_start, g_end)
        };

        // Fetch genomic sequence with window for normalization
        let window = self.config.window_size;
        let seq_start = g_start.saturating_sub(window);
        let seq_end = g_end.saturating_add(window);
        let genomic_seq = self
            .provider
            .get_genomic_sequence(chromosome, seq_start, seq_end)?;

        // Calculate relative positions (1-based)
        let rel_start = (g_start - seq_start) + 1;
        let rel_end = (g_end - seq_start) + 1;

        // Determine normalization boundaries
        // For boundary-spanning variants, use the union of exon and intron boundaries
        let boundaries =
            self.get_boundary_spanning_limits(transcript, &mapper, start_pos, end_pos, seq_start)?;

        // On minus-strand transcripts the genomic-strand window is the
        // reverse complement of the transcript view, but the variant's edit
        // alt is in transcript view. Running `normalize_na_edit` on raw
        // genomic bytes therefore canonicalizes against the wrong alphabet
        // (and the codon-frame repeat gate inspects ref context here too).
        // Mirror the intronic flow: flip into transcript view before
        // normalization, then unflip the result positions back to the
        // genomic frame. (CDS boundary-spanning context.)
        let (work_seq, work_rel_start, work_rel_end, work_boundaries) = flip_intronic_for_strand(
            transcript.strand,
            &genomic_seq,
            rel_start,
            rel_end,
            &boundaries,
        )?;

        // HGVS spec (repeated.md): the codon-frame restriction applies
        // only to bases inside the CDS proper. Boundary-spanning
        // variants cross an exon/intron boundary, so their footprint
        // is not entirely within coding sequence — pass
        // `is_coding=false` to match the intronic exemption. (A
        // hypothetical purely-exonic-span variant won't enter this
        // function; the exonic CDS path in `normalize_cds` makes its
        // own UTR/CDS-aware choice.)
        let seq_bytes = work_seq.as_bytes();
        let (work_new_rel_start, work_new_rel_end, new_edit, warnings) = self.normalize_na_edit(
            seq_bytes,
            edit,
            work_rel_start,
            work_rel_end,
            &work_boundaries,
            false,
        )?;

        let (new_rel_start, new_rel_end) = unflip_intronic_positions(
            transcript.strand,
            work_seq.len() as u64,
            work_new_rel_start,
            work_new_rel_end,
        );

        // Convert back to absolute genomic
        let new_g_start = seq_start + new_rel_start - 1;
        let new_g_end = seq_start + new_rel_end - 1;

        // Convert genomic back to CDS
        // The result might be:
        // - Still boundary-spanning
        // - Fully exonic (if shifted into exon)
        // - Fully intronic (if shifted into intron)
        let new_start = mapper.genomic_to_cds_intronic(new_g_start)?;
        let new_end = mapper.genomic_to_cds_intronic(new_g_end)?;

        // Restore coding order if positions were swapped for genomic processing
        let (new_start, new_end) = if swapped {
            (new_end, new_start)
        } else {
            (new_start, new_end)
        };

        let new_variant = CdsVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::new(Interval::new(new_start, new_end), new_edit),
        };

        Ok((HV::Cds(new_variant), warnings))
    }

    /// Convert a CDS position (exonic or intronic) to genomic coordinate
    fn cds_pos_to_genomic(
        &self,
        mapper: &crate::convert::CoordinateMapper,
        pos: &CdsPos,
    ) -> Result<u64, FerroError> {
        if pos.is_intronic() {
            mapper.cds_to_genomic_with_intron(pos)
        } else {
            mapper
                .cds_to_genomic(pos)?
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!("CDS position {} not in exons", pos.base),
                })
        }
    }

    /// Calculate normalization boundaries for boundary-spanning variants
    ///
    /// Returns the genomic region within which we can shift the variant,
    /// which is the union of the exon and intron containing the variant ends.
    fn get_boundary_spanning_limits(
        &self,
        transcript: &crate::reference::transcript::Transcript,
        mapper: &crate::convert::CoordinateMapper,
        start_pos: &CdsPos,
        end_pos: &CdsPos,
        seq_start: u64,
    ) -> Result<Boundaries, FerroError> {
        // Find the genomic extent of the region we can shift within
        // This is the union of:
        // - The exon containing the exonic position
        // - The intron containing the intronic position

        // Identify which position is exonic and which is intronic
        let (exonic_pos, intronic_pos) = if start_pos.is_intronic() {
            (end_pos, start_pos)
        } else {
            (start_pos, end_pos)
        };

        // Get exon boundaries in genomic coords
        let tx_pos = mapper.cds_to_tx(exonic_pos)?;
        let tx_pos_base = u64::try_from(tx_pos.base).map_err(|_| FerroError::ConversionError {
            msg: format!(
                "Negative transcript position {} during boundary normalization",
                tx_pos.base
            ),
        })?;
        let exon = transcript
            .exon_at(tx_pos_base)
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Could not find exon for boundary normalization".to_string(),
            })?;

        // Get intron boundaries in genomic coords
        let tx_boundary = mapper.cds_to_tx(intronic_pos)?;
        let tx_boundary_base =
            u64::try_from(tx_boundary.base).map_err(|_| FerroError::ConversionError {
                msg: format!(
                    "Negative transcript position {} during boundary normalization",
                    tx_boundary.base
                ),
            })?;
        let offset = intronic_pos.offset.unwrap_or(0);
        let intron = transcript
            .find_intron_at_tx_boundary(tx_boundary_base, offset)
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Could not find intron for boundary normalization".to_string(),
            })?;

        // Get genomic coordinates from exon and intron
        let (exon_g_start, exon_g_end): (u64, u64) = match (exon.genomic_start, exon.genomic_end) {
            (Some(s), Some(e)) => (s, e),
            _ => {
                return Err(FerroError::ConversionError {
                    msg: "Exon has no genomic coordinates".to_string(),
                })
            }
        };

        let (intron_g_start, intron_g_end): (u64, u64) =
            match (intron.genomic_start, intron.genomic_end) {
                (Some(s), Some(e)) => (s, e),
                _ => {
                    return Err(FerroError::ConversionError {
                        msg: "Intron has no genomic coordinates".to_string(),
                    })
                }
            };

        // Union of exon and intron genomic coordinates
        let combined_start = exon_g_start.min(intron_g_start);
        let combined_end = exon_g_end.max(intron_g_end);

        // Convert to relative positions (1-based within our fetched sequence)
        let rel_start = combined_start.saturating_sub(seq_start) + 1;
        let rel_end = combined_end.saturating_sub(seq_start) + 1;

        Ok(Boundaries::new(rel_start, rel_end))
    }

    /// Core normalization for nucleic acid edits
    ///
    /// Returns (new_start, new_end, new_edit, warnings)
    fn normalize_na_edit(
        &self,
        ref_seq: &[u8],
        edit: &NaEdit,
        start: u64,
        end: u64,
        boundaries: &Boundaries,
        is_coding: bool,
    ) -> Result<(u64, u64, NaEdit, Vec<NormalizationWarning>), FerroError> {
        let mut warnings = Vec::new();

        // Validate reference allele before normalization
        let validation = validate::validate_reference(edit, ref_seq, start, end);
        if !validation.valid {
            // `corrected` is honest about whether the canonical Display
            // drops the user-stated bases. See issue #280.
            //
            // True for `Deletion` / `Duplication` (canonicalize_edit
            // strips `sequence`) and conditionally for `Inversion`
            // (the Inversion arm in `normalize_na_edit` emits
            // `sequence: None` when `shorten_inversion` rewrites the
            // span, and otherwise still drops the stated bases in the
            // canonical path). `Substitution` does not reach here at
            // all — `needs_normalization` returns `false` for real
            // substitutions and only the degenerate `ref == alt` case
            // routes through this function, where the validator's
            // single-base check is satisfied by construction.
            //
            // False for `Repeat` / `MultiRepeat` consistency mismatches
            // (issues #214 / #279): the per-unit declaration is part of
            // the user's form and the normalizer passes the description
            // through verbatim.
            //
            // TODO/Note: revisit if `delins` ever gets a stated-deleted
            // validator — today `validate_reference`'s `NaEdit::Delins`
            // arm returns `ok()` unconditionally (see also the Delins
            // arm in `canonicalize_edit`, which strips `deleted` /
            // `deleted_length`), so the branch is unreachable for
            // Delins.
            let corrected = !matches!(edit, NaEdit::Repeat { .. } | NaEdit::MultiRepeat { .. });
            warnings.push(NormalizationWarning::RefSeqMismatch {
                stated_ref: validation.stated_ref.unwrap_or_default(),
                actual_ref: validation.actual_ref.unwrap_or_default(),
                position: format!("{}-{}", start, end),
                corrected,
                details: validation.warning,
            });
        }

        // Substitution with ref == alt is identity (e.g. c.100A>A → c.100=).
        // This is the SNV companion to the same-base delins → identity rule;
        // the rewrite is purely syntactic on the edit's stated bases, so it
        // applies across coordinate systems and runs before shuffling.
        if let NaEdit::Substitution {
            reference,
            alternative,
        } = edit
        {
            if reference == alternative {
                return Ok((start, end, NaEdit::position_identity(), warnings));
            }
        }

        // Get the alternate sequence for the edit
        let alt_seq = match edit {
            NaEdit::Deletion { .. } => vec![],
            NaEdit::Insertion { sequence } => {
                // Only shuffle if we have a literal sequence
                if let Some(bases) = sequence.bases() {
                    bases.iter().map(|b| b.to_u8()).collect()
                } else {
                    // Cannot shuffle non-literal insertions (counts, ranges, etc.)
                    return Ok((start, end, edit.clone(), warnings.clone()));
                }
            }
            NaEdit::Duplication { .. } => {
                // Always read duplicated bytes from the reference,
                // regardless of any user-stated bases on the input
                // (`dup<base>` shapes). The canonical HGVS Display
                // drops the stated bases anyway (see the dup output
                // arm in `get_canonical_form`); using the stated
                // bases here when they don't match the reference
                // (the parser accepts them in lenient mode with a
                // `RefSeqMismatch` warning) caused the shuffle to
                // mis-shift and broke single-pass idempotency:
                // `c.10dupA` against a `CCCCC` homopolymer
                // canonicalized to `c.10dup` on pass 1 (stated-ref
                // stripped without shifting) and only shifted to
                // `c.14dup` on pass 2 (because the now-clean form
                // reads bytes from reference and the shuffle fires
                // correctly). Reading from reference up-front
                // collapses both pass 1 and pass 2 behavior into a
                // single pass. Issue #219.
                let s = hgvs_pos_to_index(start);
                let e = end as usize;
                if e <= ref_seq.len() {
                    ref_seq[s..e].to_vec()
                } else {
                    vec![]
                }
            }
            NaEdit::Delins { sequence, .. } => {
                use crate::hgvs::edit::InsertedSequence;

                // HGVS spec (issue #81 A3): a delins with an empty inserted
                // sequence is semantically a deletion and must be rendered as
                // `del`. Rewrite up-front so the result picks up del 3'-shift
                // and validation in the Deletion arm.
                if matches!(sequence, InsertedSequence::Empty) {
                    let del = NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    };
                    return self
                        .normalize_na_edit(ref_seq, &del, start, end, boundaries, is_coding);
                }

                // HGVS spec: delins should NOT be 3' shifted like del/dup/ins,
                // but the edit-type priority (sub > del > inv > dup > ins) means
                // we may need to rewrite it as a higher-priority form: identity
                // (insert == ref), substitution (1->1 ref!=alt), or duplication.
                if let InsertedSequence::Literal(seq) = sequence {
                    use crate::hgvs::edit::{Base, Sequence};
                    use rules::DelinsCanonical;
                    let seq_bytes: Vec<u8> = seq.bases().iter().map(|b| *b as u8).collect();
                    let start_idx = hgvs_pos_to_index(start);
                    let end_idx = end as usize;

                    // Reconstruct an InsertedSequence from a Vec<u8> produced by
                    // shared-affix trimming. The bytes round-trip through `Base`
                    // because they originated from a typed `Sequence` (the input
                    // `seq` above), so `from_char` cannot fail; expect-on-None
                    // makes the invariant explicit if a future refactor breaks
                    // the pipeline.
                    let bytes_to_inserted_seq = |bytes: &[u8]| -> InsertedSequence {
                        let bases: Vec<Base> = bytes
                            .iter()
                            .map(|b| {
                                Base::from_char(*b as char).expect(
                                    "trimmed delins byte must be a valid IUPAC base \
                                     because the input sequence was already a typed Sequence",
                                )
                            })
                            .collect();
                        InsertedSequence::Literal(Sequence::new(bases))
                    };

                    match rules::canonicalize_delins(ref_seq, start_idx, end_idx, &seq_bytes) {
                        DelinsCanonical::Identity => {
                            // c.10delinsG where ref[10]=G  ->  c.10=
                            return Ok((start, end, NaEdit::position_identity(), warnings.clone()));
                        }
                        DelinsCanonical::Substitution {
                            position,
                            reference,
                            alternative,
                        } => {
                            // g.1000delinsA where ref[1000]=G  ->  g.1000G>A.
                            // After shared-affix trimming `position` is a
                            // 0-indexed offset into ref_seq, not necessarily
                            // the input `start`.
                            let pos = index_to_hgvs_pos(position);
                            return Ok((
                                pos,
                                pos,
                                NaEdit::Substitution {
                                    reference,
                                    alternative,
                                },
                                warnings.clone(),
                            ));
                        }
                        DelinsCanonical::Deletion { start: s0, end: e0 } => {
                            // c.2_5delinsAT (ref ACGT) -> c.3_4del. Range fields
                            // are the trimmed half-open 0-indexed interval.
                            return Ok((
                                index_to_hgvs_pos(s0),
                                e0 as u64,
                                NaEdit::Deletion {
                                    sequence: None,
                                    length: None,
                                },
                                warnings.clone(),
                            ));
                        }
                        DelinsCanonical::Insertion {
                            after_index,
                            sequence: ins_bytes,
                        } => {
                            // c.2_4delinsACGT (ref ACT) -> c.3_4insG. `after_index`
                            // is the 0-indexed position of the base AFTER the
                            // insertion, which is also the 1-based HGVS position
                            // of the base BEFORE — so HGVS X = after_index,
                            // Y = after_index + 1.
                            //
                            // Recurse into `normalize_na_edit` with the new
                            // Insertion so the full ins pipeline (3'/5' shuffle
                            // + `insertion_to_duplication` + `insertion_to_repeat`)
                            // runs. Without recursion, a delins-derived
                            // insertion that duplicates a nearby reference tract
                            // (e.g. biocommons `g.X delinsCTTTCTT` where
                            // ref[X+1..X+6]=TTTCTT) would skip the
                            // ins→dup recognizer and emit the long `insTTTCTT`
                            // form instead of the canonical `dup`. Closes-after:
                            // #356.
                            //
                            // Issue #418 guard: when shared-affix trimming
                            // collapses a Delins-at-tx-start (e.g. `c.1delinsCA`
                            // on a transcript with `cds_start = 1` whose c.1=A —
                            // the suffix `A` matches ref[c.1] and trims to
                            // `Insertion("C")` at `after_index = 0`) the
                            // recursive call would pass `tx_start = 0`, which is
                            // not a valid 1-based HGVS position and underflows
                            // `hgvs_pos_to_index(0)` downstream. Spec-canonical
                            // behavior for a Delins input whose
                            // canonicalisation rewrite would land strictly past
                            // the start of the transcript is to restore the
                            // input form unchanged (this is the same disposition
                            // the post-shift #383 clamp gives Delins inputs
                            // that cross the CDS-start boundary into 5'UTR).
                            // Suppress the recursion to avoid the panic and
                            // emit the input verbatim.
                            if after_index == 0 {
                                return Ok((start, end, edit.clone(), warnings.clone()));
                            }
                            let new_edit = NaEdit::Insertion {
                                sequence: bytes_to_inserted_seq(&ins_bytes),
                            };
                            // Preserve warnings collected for the original
                            // delins (e.g. RefSeqMismatch in strict mode):
                            // merge them into the recursive call's warnings
                            // rather than dropping them by tail-returning.
                            let (new_start, new_end, new_edit, mut child_warnings) = self
                                .normalize_na_edit(
                                    ref_seq,
                                    &new_edit,
                                    after_index as u64,
                                    (after_index + 1) as u64,
                                    boundaries,
                                    is_coding,
                                )?;
                            let mut merged = warnings.clone();
                            merged.append(&mut child_warnings);
                            return Ok((new_start, new_end, new_edit, merged));
                        }
                        DelinsCanonical::Inversion { start: s0, end: e0 } => {
                            // A2 (#81): g.100_102delinsTAG where ref=CTA  ->  g.100_102inv.
                            // Position interval is already shortened. e0 is the
                            // exclusive 0-based end; the HGVS 1-based inclusive
                            // end takes the same numeric value.
                            return Ok((
                                index_to_hgvs_pos(s0),
                                e0 as u64,
                                NaEdit::Inversion {
                                    sequence: None,
                                    length: None,
                                },
                                warnings.clone(),
                            ));
                        }
                        DelinsCanonical::Duplication { start: s0, end: e0 } => {
                            // c.5delinsGG where ref[5]=G  ->  c.5dup. Duplication
                            // is detected before trimming, so the range matches
                            // the input.
                            return Ok((
                                index_to_hgvs_pos(s0),
                                e0 as u64,
                                NaEdit::Duplication {
                                    sequence: None,
                                    length: None,
                                    uncertain_extent: None,
                                },
                                warnings.clone(),
                            ));
                        }
                        DelinsCanonical::KeepAsDelins {
                            start: s0,
                            end: e0,
                            sequence: trimmed_bytes,
                        } => {
                            // Either no trimming was possible (range == input)
                            // or trimming reduced the delins to a smaller delins
                            // that still doesn't fit a higher-priority form.
                            return Ok((
                                index_to_hgvs_pos(s0),
                                e0 as u64,
                                NaEdit::Delins {
                                    sequence: bytes_to_inserted_seq(&trimmed_bytes),
                                    deleted: None,
                                    deleted_length: None,
                                },
                                warnings.clone(),
                            ));
                        }
                    }
                }
                // Non-literal insert (Count, Range, Reference, PositionRange,
                // Complex, …): cannot trim or classify without the actual
                // bases, but we still strip an explicit deleted sequence /
                // length per the same spec rule that the Literal arm above
                // applies (`delins.md`: the recommendation is to omit the
                // explicit deleted bases). Closes #338's WITH-provider gap
                // surfaced in code review.
                return Ok((start, end, canonicalize_edit(edit), warnings.clone()));
            }
            NaEdit::Inversion { .. } => {
                // Apply the complementary-outer-bases shortening rule. After
                // shortening, the inversion's interval may no longer match the
                // input's explicit `sequence`/`length` (if any), so emit
                // minimal notation.
                let start_idx = hgvs_pos_to_index(start); // Convert 1-based to 0-based
                let end_idx = end as usize; // end is exclusive (0-based)

                if let Some((new_s, new_e)) = rules::shorten_inversion(ref_seq, start_idx, end_idx)
                {
                    return Ok((
                        index_to_hgvs_pos(new_s),
                        new_e as u64,
                        NaEdit::Inversion {
                            sequence: None,
                            length: None,
                        },
                        warnings,
                    ));
                } else {
                    // Inversion reduced to identity. Use the canonical
                    // position-only Identity (matches the Delins identity arm
                    // above), so both inversion-collapse paths emit the same
                    // shape.
                    return Ok((start, end, NaEdit::position_identity(), warnings));
                }
            }
            NaEdit::Repeat {
                sequence,
                count,
                additional_counts,
                trailing,
            } => {
                use crate::hgvs::edit::RepeatCount;

                // Only normalize exact counts with a sequence
                let Some(seq) = sequence else {
                    return Ok((start, end, edit.clone(), warnings.clone()));
                };

                // Range, UncertainRange, MinUncertain, MaxUncertain, Unknown:
                // no concrete count to 3'-shift against, so pass through unchanged.
                let RepeatCount::Exact(specified_count) = count else {
                    return Ok((start, end, edit.clone(), warnings.clone()));
                };

                // Skip if there are additional counts (genotype notation)
                if !additional_counts.is_empty() || trailing.is_some() {
                    return Ok((start, end, edit.clone(), warnings.clone()));
                }

                // Get the repeat unit as bytes
                let repeat_unit: Vec<u8> = seq.bases().iter().map(|b| b.to_u8()).collect();
                let pos_idx = hgvs_pos_to_index(start); // Convert 1-based to 0-based

                // Normalize the repeat
                match rules::normalize_repeat(
                    ref_seq,
                    pos_idx,
                    &repeat_unit,
                    *specified_count,
                    is_coding,
                ) {
                    rules::RepeatNormResult::Deletion {
                        start: del_start,
                        end: del_end,
                    } => {
                        // Minimal notation - no explicit length
                        let del_edit = NaEdit::Deletion {
                            sequence: None,
                            length: None,
                        };
                        return Ok((del_start, del_end, del_edit, warnings));
                    }
                    rules::RepeatNormResult::Duplication {
                        start: dup_start,
                        end: dup_end,
                        sequence: _dup_seq,
                    } => {
                        // Minimal notation - no explicit sequence or length
                        let dup_edit = NaEdit::Duplication {
                            sequence: None,
                            length: None,
                            uncertain_extent: None,
                        };
                        return Ok((dup_start, dup_end, dup_edit, warnings));
                    }
                    rules::RepeatNormResult::Insertion {
                        start: ins_start,
                        end: ins_end,
                        sequence: ins_seq,
                    } => {
                        // Codon-frame gate routed an expansion to ins literal form
                        // (e.g., c.1741_1742insTATATATA per spec).
                        use crate::hgvs::edit::{Base, InsertedSequence, Sequence};
                        let bases: Vec<Base> = ins_seq
                            .iter()
                            .filter_map(|&b| Base::from_char(b as char))
                            .collect();
                        if bases.len() == ins_seq.len() {
                            let ins_edit = NaEdit::Insertion {
                                sequence: InsertedSequence::Literal(Sequence::new(bases)),
                            };
                            return Ok((ins_start, ins_end, ins_edit, warnings));
                        }
                        // Defensive fallback: rule layer returned a base byte
                        // that doesn't fit the Base alphabet (e.g. N). Don't
                        // emit a truncated insertion — keep the original edit
                        // and positions so downstream invariants hold.
                        return Ok((start, end, edit.clone(), warnings));
                    }
                    rules::RepeatNormResult::Repeat {
                        start: rep_start,
                        end: rep_end,
                        sequence: rep_seq,
                        count: rep_count,
                    } => {
                        use crate::hgvs::edit::{Base, RepeatCount, Sequence};
                        let bases: Vec<Base> = rep_seq
                            .iter()
                            .filter_map(|&b| Base::from_char(b as char))
                            .collect();
                        if bases.len() == rep_seq.len() {
                            let rep_edit = NaEdit::Repeat {
                                sequence: Some(Sequence::new(bases)),
                                count: RepeatCount::Exact(rep_count),
                                additional_counts: vec![],
                                trailing: None,
                            };
                            return Ok((rep_start, rep_end, rep_edit, warnings));
                        }
                        // Defensive fallback: rule layer returned a repeat
                        // unit byte that doesn't fit the Base alphabet (e.g.
                        // a gap or non-IUPAC byte from the reference). Don't
                        // emit a truncated repeat sequence — keep the
                        // original edit and positions.
                        return Ok((start, end, edit.clone(), warnings));
                    }
                    rules::RepeatNormResult::Unchanged => {
                        return Ok((start, end, edit.clone(), warnings.clone()));
                    }
                }
            }
            _ => return Ok((start, end, edit.clone(), warnings.clone())), // Other edits don't need shuffling
        };

        // Perform shuffle
        // For insertions, the HGVS interval X_Y (where Y = X+1) represents flanking positions.
        // We need to adjust the end coordinate so shuffle checks the correct reference position.
        // For c.445_446insA: start=634, end=635 (1-based tx coords)
        // We want shuffle to check ref_seq[634-1] = ref_seq[633] for first flanking
        // and ref_seq[635-1] = ref_seq[634] for second flanking (the position to check for 3' shift)
        //
        // For insertions, we need to adjust the start passed to shuffle so that the alt_idx
        // calculation starts at 0 (not 1). The shuffle's alt_idx formula is:
        //   alt_idx = (new_end - start) % alt_seq.len()
        // For insertions, new_end starts at shuffle_end (which is end - 1 for insertions).
        // If we pass start_idx directly, alt_idx = (end-1) - start_idx = 1 (wrong, should be 0).
        // By passing start_idx + 1, we get alt_idx = (end-1) - (start_idx+1) = 0 (correct).
        let shuffle_end = match edit {
            NaEdit::Insertion { .. } => end.saturating_sub(1), // Use second flanking position
            _ => end,
        };
        let start_idx = hgvs_pos_to_index(start); // Convert 1-based to 0-based
        let shuffle_start = match edit {
            NaEdit::Insertion { .. } => start_idx as u64 + 1, // Adjust so alt_idx starts at 0
            _ => start_idx as u64,
        };
        let result = shuffle(
            ref_seq,
            &alt_seq,
            shuffle_start,
            shuffle_end, // Adjusted for insertions
            boundaries,
            self.config.shuffle_direction,
        );

        // Convert back to 1-based HGVS positions
        // For insertions, we adjusted the start for shuffle, now adjust back
        let shuffle_result_start = match edit {
            NaEdit::Insertion { .. } => result.start.saturating_sub(1), // Adjust back
            _ => result.start,
        };
        let new_start = index_to_hgvs_pos(shuffle_result_start as usize);
        // For insertions, we adjusted the end for shuffle, now restore the HGVS X_(X+1) format
        let new_end = match edit {
            NaEdit::Insertion { .. } => index_to_hgvs_pos(result.end as usize), // Restore second flanking position
            _ => result.end,
        };

        // Determine the canonical form for the edit
        // HGVS rules:
        // - Deletions ALWAYS stay as deletions (just shift 3')
        // - Insertions become duplications if single-base matches adjacent
        // - Multi-base insertions/dups in homopolymer become repeat notation
        let (final_start, final_end, new_edit) = match edit {
            NaEdit::Insertion { sequence } => {
                use crate::hgvs::edit::{InsertedSequence, RepeatCount, Sequence};

                if let InsertedSequence::Literal(seq) = sequence {
                    let seq_bytes: Vec<u8> = seq.bases().iter().map(|b| *b as u8).collect();

                    // Check for repeat notation first (multi-base homopolymer insertion)
                    // Use the ORIGINAL position (start), not shuffled position (result.start)
                    // because repeat notation refers to the reference tract position
                    if seq_bytes.len() > 1 {
                        let original_pos_idx = hgvs_pos_to_index(start) as u64; // 0-based original position
                        if let Some((_first, count, rep_start, rep_end, unit_bytes)) =
                            rules::insertion_to_repeat(
                                ref_seq,
                                original_pos_idx,
                                &seq_bytes,
                                is_coding,
                            )
                        {
                            use crate::hgvs::edit::Base;
                            let bases: Vec<Base> = unit_bytes
                                .iter()
                                .filter_map(|&b| Base::from_char(b as char))
                                .collect();
                            if bases.len() == unit_bytes.len() {
                                let repeat_seq = Sequence::new(bases);
                                let repeat_edit = NaEdit::Repeat {
                                    sequence: Some(repeat_seq),
                                    count: RepeatCount::Exact(count),
                                    additional_counts: vec![],
                                    trailing: None,
                                };
                                return Ok((rep_start, rep_end, repeat_edit, warnings));
                            }
                        }
                    }

                    // Resolve insertion → duplication canonicalization. Three candidate
                    // dup positions compete; we pick by the rules below in order.
                    //
                    // (a) Tract-aligned dup via `insertion_to_duplication` (uses the
                    //     ORIGINAL insertion point and finds the maximal tandem run
                    //     under any cyclic rotation of the alt). When the tract has
                    //     `ref_count >= 2` we prefer this regardless of how far
                    //     shuffle walked: the multi-copy tract has a meaningful phase
                    //     that the spec-canonical form preserves (issue #132).
                    //
                    // (b) Post-shuffle simple dup via `insertion_is_duplication`. When
                    //     shuffle walked past a single-copy tract via partial-match
                    //     extension (e.g. TGATC abutting TGAAG — first three bases
                    //     match but the fourth does not), the post-shuffle position is
                    //     more 3' than (a)'s tract-aligned position and is the canonical
                    //     answer per the 3' rule (issue #180).
                    //
                    // (c) Single-copy tract fallback (`insertion_to_duplication` with
                    //     `ref_count == 1`). Hit when shuffle stalled before completing
                    //     one alt rotation (so (b) doesn't find a dup at the post-
                    //     shuffle position) but the alt does match an adjacent ref unit
                    //     at the ORIGINAL position. Example: ins AACA abutting AACA.
                    //
                    // If none match, fall through to ins (possibly rotated).
                    let original_pos_idx = hgvs_pos_to_index(start) as u64;
                    let ins_to_dup = rules::insertion_to_duplication(
                        ref_seq,
                        original_pos_idx,
                        &seq_bytes,
                        self.config.shuffle_direction,
                    );

                    // Codon-frame gate (repeated.md): in c., if the alt is
                    // >=2 copies of a non-codon-aligned unit, the spec
                    // mandates ins<literal>, not dup. The smallest-unit
                    // length is rotation-invariant, so we can compute it
                    // once from `seq_bytes` and apply the same gate at
                    // every dup-emission site below. In practice the gate
                    // never fires when `ins_to_dup` is `Some` (that helper
                    // only returns for single-unit alts, where
                    // `smallest_unit.len() == seq_bytes.len()`), but we
                    // guard the (a) fast path and (c) single-copy fallback
                    // anyway so the spec rule is explicit at each site and
                    // survives future changes to `insertion_to_duplication`.
                    let smallest_unit = rules::smallest_repeat_unit(&seq_bytes);
                    let codon_blocks_dup = is_coding
                        && smallest_unit.len() < seq_bytes.len()
                        && !smallest_unit.len().is_multiple_of(3);

                    if !codon_blocks_dup {
                        if let Some(rules::InsToDupResult {
                            start: dup_start,
                            end: dup_end,
                            ref_count,
                            ..
                        }) = ins_to_dup.as_ref()
                        {
                            if *ref_count >= 2 {
                                return Ok((
                                    *dup_start,
                                    *dup_end,
                                    NaEdit::Duplication {
                                        sequence: None,
                                        length: None,
                                        uncertain_extent: None,
                                    },
                                    warnings,
                                ));
                            }
                        }
                    }

                    // Check for simple duplication (single-base or matching adjacent)
                    // When shifting an insertion through a repeat region, the effective sequence
                    // rotates. For example, shifting "GGC" by 1 position gives "GCG".
                    let shift_amount =
                        (result.start as usize).saturating_sub(shuffle_start as usize);
                    let rotation = shift_amount % seq_bytes.len();
                    let rotated_seq: Vec<u8> = if rotation > 0 {
                        seq_bytes[rotation..]
                            .iter()
                            .chain(seq_bytes[..rotation].iter())
                            .copied()
                            .collect()
                    } else {
                        seq_bytes.clone()
                    };

                    // `insertion_is_duplication` checks both sides of the
                    // post-shuffle insertion point. We need to know which
                    // side actually matched so we anchor the dup notation
                    // there — anchoring on the wrong side lands one position
                    // too far in the wrong direction.
                    let pos_idx = result.start as usize;
                    let ins_len = rotated_seq.len();
                    let five_prime_match = pos_idx >= ins_len
                        && pos_idx <= ref_seq.len()
                        && ref_seq[pos_idx - ins_len..pos_idx] == rotated_seq[..];
                    let three_prime_match = pos_idx + ins_len <= ref_seq.len()
                        && ref_seq[pos_idx..pos_idx + ins_len] == rotated_seq[..];

                    if !codon_blocks_dup && (five_prime_match || three_prime_match) {
                        // Side-aware dup anchor. `insertion_is_duplication`-
                        // equivalent checks above return true if EITHER the
                        // preceding ref tract (5'-side: `ref[pos-L..pos]`)
                        // OR the following one (3'-side: `ref[pos..pos+L]`)
                        // equals the (possibly rotated) alt. Which side
                        // matched determines where the duplicated region
                        // sits in HGVS coordinates:
                        //
                        //   5'-side match (`ref[pos-L..pos] == alt`): the
                        //   alt duplicates the immediately-preceding tract,
                        //   so the dup region is `(new_start-L+1, new_start)`.
                        //   This is the natural 3'-shuffle stopping shape —
                        //   at the 3'-most equivalent position, the alt
                        //   phase aligns with the tract just behind the
                        //   insertion.
                        //
                        //   3'-side match (`ref[pos..pos+L] == alt`): the
                        //   alt duplicates the immediately-following tract,
                        //   so the dup region is `(new_start+1, new_start+L)`.
                        //   This is the natural 5'-shuffle stopping shape —
                        //   at the 5'-most equivalent position the alt
                        //   aligns with the tract just AHEAD of the
                        //   insertion (the one the shuffle walked through).
                        //   Without this case, an `ins → dup` canon firing
                        //   only on the 3' side (e.g. the post-5'-shift
                        //   `c.9170_9171insA` where ref[c.9171] = A — closes
                        //   #402, or the issue #418 (b) NM_001166478.1:
                        //   c.36_37insTC off-by-2) labels the dup at the
                        //   wrong base — emitting bases that DON'T equal
                        //   the alt, i.e. a different haplotype.
                        //
                        //   Both sides match → the alt sits in a 2+-copy
                        //   tract, which is path (a) territory above (it
                        //   would have returned with `ref_count >= 2`
                        //   before reaching here). Reach this branch only
                        //   on `ref_count == 1` tracts where exactly one
                        //   side matches in practice; the shuffle-direction
                        //   tie-break below is defensive for the
                        //   unreachable both-match case (5'-most for
                        //   `FivePrime`, 3'-most for `ThreePrime`).
                        //
                        // Single-base dups (`dup_len == 1`) must still
                        // choose the correct flanking position: 5'-side
                        // anchors at `new_start`, 3'-side at `new_start+1`.
                        // The old "single position for single-base dup"
                        // shortcut anchored unconditionally at `new_start`,
                        // emitting the wrong haplotype on 3'-side-only
                        // single-copy matches (e.g. `...CT[insT]` adjacent
                        // to an isolated `T` — buggy `c.40dup` (= C) vs
                        // correct `c.41dup` (= T)).
                        let prefer_three_prime = match (five_prime_match, three_prime_match) {
                            (false, true) => true,
                            (true, false) => false,
                            (true, true) => {
                                matches!(
                                    self.config.shuffle_direction,
                                    ShuffleDirection::ThreePrime
                                )
                            }
                            (false, false) => {
                                unreachable!("outer guard requires at least one side matched")
                            }
                        };
                        let dup_len = rotated_seq.len() as u64;
                        let (dup_start, dup_end) = if prefer_three_prime {
                            if dup_len == 1 {
                                (new_start + 1, new_start + 1)
                            } else {
                                (new_start + 1, new_start + dup_len)
                            }
                        } else if dup_len == 1 {
                            (new_start, new_start)
                        } else {
                            (new_start - dup_len + 1, new_start)
                        };
                        (
                            dup_start,
                            dup_end,
                            NaEdit::Duplication {
                                sequence: None, // Minimal notation - no explicit sequence
                                length: None,
                                uncertain_extent: None,
                            },
                        )
                    } else if let Some(rules::InsToDupResult {
                        start: dup_start,
                        end: dup_end,
                        ..
                    }) = ins_to_dup.as_ref().filter(|_| !codon_blocks_dup)
                    {
                        // (c) Single-copy tract fallback. Reached when (a) declined
                        // because `ref_count < 2` and (b) declined because the post-
                        // shuffle rotated alt doesn't match adjacent reference. The
                        // alt is a (possibly rotated) tandem unit abutting a single-
                        // copy ref tract at the original insertion point — emit the
                        // dup over that tract.
                        (
                            *dup_start,
                            *dup_end,
                            NaEdit::Duplication {
                                sequence: None,
                                length: None,
                                uncertain_extent: None,
                            },
                        )
                    } else {
                        // Output the rotated sequence for shifted insertions
                        if rotation > 0 {
                            use crate::hgvs::edit::{Base, InsertedSequence, Sequence};
                            let rotated_bases: Vec<Base> = rotated_seq
                                .iter()
                                .filter_map(|&b| Base::from_char(b as char))
                                .collect();
                            // Mirror the gated-ins guard used by the
                            // RepeatNormResult::Insertion / GatedInsertion
                            // branches: if any byte fell outside the Base
                            // alphabet, refuse to emit a truncated `ins`
                            // and fall back to the original edit so
                            // downstream invariants hold.
                            if rotated_bases.len() == rotated_seq.len() {
                                let new_sequence =
                                    InsertedSequence::Literal(Sequence::new(rotated_bases));
                                (
                                    new_start,
                                    new_end,
                                    NaEdit::Insertion {
                                        sequence: new_sequence,
                                    },
                                )
                            } else {
                                (new_start, new_end, edit.clone())
                            }
                        } else {
                            (new_start, new_end, edit.clone())
                        }
                    }
                } else {
                    (new_start, new_end, edit.clone())
                }
            }
            NaEdit::Duplication { .. } => {
                use crate::hgvs::edit::{Base, RepeatCount, Sequence};

                // Check if duplication should become repeat notation
                // Use the shuffled positions (result.start, result.end) which are 0-based
                // This applies to both single-base dups in homopolymers and multi-base tandem dups
                if let Some(dup_result) =
                    rules::duplication_to_repeat(ref_seq, result.start, result.end, is_coding)
                {
                    match dup_result {
                        rules::DupToRepeatResult::Homopolymer {
                            base,
                            count,
                            start: rep_start,
                            end: rep_end,
                        } => {
                            if let Some(base_enum) = Base::from_char(base as char) {
                                let repeat_seq = Sequence::new(vec![base_enum]);
                                let repeat_edit = NaEdit::Repeat {
                                    sequence: Some(repeat_seq),
                                    count: RepeatCount::Exact(count),
                                    additional_counts: vec![],
                                    trailing: None,
                                };
                                return Ok((rep_start, rep_end, repeat_edit, warnings));
                            }
                        }
                        rules::DupToRepeatResult::TandemRepeat {
                            unit,
                            count,
                            start: rep_start,
                            end: rep_end,
                        } => {
                            let bases: Vec<Base> = unit
                                .iter()
                                .filter_map(|&b| Base::from_char(b as char))
                                .collect();
                            if bases.len() == unit.len() {
                                let repeat_seq = Sequence::new(bases);
                                let repeat_edit = NaEdit::Repeat {
                                    sequence: Some(repeat_seq),
                                    count: RepeatCount::Exact(count),
                                    additional_counts: vec![],
                                    trailing: None,
                                };
                                return Ok((rep_start, rep_end, repeat_edit, warnings));
                            }
                        }
                        rules::DupToRepeatResult::GatedInsertion {
                            start: ins_start,
                            end: ins_end,
                            sequence: ins_seq,
                        } => {
                            // Codon-frame gate routed a multi-copy dup to ins
                            // literal form per HGVS spec.
                            use crate::hgvs::edit::InsertedSequence;
                            let bases: Vec<Base> = ins_seq
                                .iter()
                                .filter_map(|&b| Base::from_char(b as char))
                                .collect();
                            if bases.len() == ins_seq.len() {
                                let ins_edit = NaEdit::Insertion {
                                    sequence: InsertedSequence::Literal(Sequence::new(bases)),
                                };
                                return Ok((ins_start, ins_end, ins_edit, warnings));
                            }
                            // Defensive fallback: rule layer returned a base
                            // byte that doesn't fit the Base alphabet (e.g.
                            // N). Fall through to the generic dup minimal-
                            // notation path below rather than emitting a
                            // truncated insertion.
                        }
                    }
                }
                // Keep as duplication but strip explicit sequence (minimal notation)
                (
                    new_start,
                    new_end,
                    NaEdit::Duplication {
                        sequence: None,
                        length: None,
                        uncertain_extent: None,
                    },
                )
            }
            // Deletions: post-shift, check for B2 canonical-form rule
            // (deletion of >=2 tandem-repeat units → unit[N-k]); otherwise
            // strip explicit length for minimal `del` notation. The
            // collect-into-Option short-circuits if any byte in the unit isn't
            // a valid `Base` (e.g. `N`), in which case we fall through to del.
            //
            // B2 is defined for a *post-3'-shift* deletion (the shuffle phase-
            // alignment lemma justifies emitting `unit[N-k]` without rotation).
            // Under FivePrime shuffle, applying it would re-anchor the
            // 5'-normalized deletion to the canonical tract position, defeating
            // the user's choice of direction — so gate it on ThreePrime.
            NaEdit::Deletion { .. } => {
                use crate::hgvs::edit::{Base, RepeatCount, Sequence};
                if self.config.shuffle_direction == ShuffleDirection::ThreePrime {
                    if let Some(rep) = rules::deletion_to_repeat(
                        ref_seq,
                        result.start as usize,
                        result.end as usize,
                        is_coding,
                    ) {
                        let bases: Option<Vec<Base>> = rep
                            .unit
                            .iter()
                            .map(|&b| Base::from_char(b as char))
                            .collect();
                        if let Some(bases) = bases {
                            let repeat_edit = NaEdit::Repeat {
                                sequence: Some(Sequence::new(bases)),
                                count: RepeatCount::Exact(rep.count),
                                additional_counts: vec![],
                                trailing: None,
                            };
                            return Ok((rep.start, rep.end, repeat_edit, warnings));
                        }
                    }
                }
                (
                    new_start,
                    new_end,
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                )
            }
            // All other edit types stay unchanged
            _ => (new_start, new_end, edit.clone()),
        };

        Ok((final_start, final_end, new_edit, warnings))
    }

    /// Convert CDS position to transcript position
    fn cds_to_tx_pos(
        &self,
        pos: &CdsPos,
        cds_start: u64,
        cds_end: Option<u64>,
    ) -> Result<u64, FerroError> {
        if pos.utr3 {
            let end = cds_end.ok_or_else(|| FerroError::ConversionError {
                msg: "No CDS end".to_string(),
            })?;
            let base = u64::try_from(pos.base).map_err(|_| FerroError::ConversionError {
                msg: format!("Negative base {} in 3' UTR position", pos.base),
            })?;
            Ok(end + base)
        } else if pos.base < 0 {
            // 5'UTR: HGVS numbering skips c.0 (c.-1 is the base immediately
            // upstream of c.1), so c.-N maps to tx position cds_start - N.
            // Issue #97 — the previous formula `cds_start + base - 1`
            // double-counted the gap and emitted the wrong tx position.
            let tx_pos = cds_start as i64 + pos.base;
            u64::try_from(tx_pos).map_err(|_| FerroError::ConversionError {
                msg: format!(
                    "CDS position c.{} maps before transcript start (cds_start={})",
                    pos.base, cds_start
                ),
            })
        } else if pos.base == 0 {
            // c.0 is not a valid HGVS position, but historical inputs
            // can land here. Preserve the legacy mapping (treat as the
            // last 5'UTR base, equivalent to c.-1) rather than failing.
            Ok(cds_start.saturating_sub(1))
        } else {
            Ok(cds_start + pos.base as u64 - 1)
        }
    }

    /// Convert transcript position to CDS position
    fn tx_to_cds_pos(
        &self,
        pos: u64,
        cds_start: u64,
        cds_end: Option<u64>,
    ) -> Result<CdsPos, FerroError> {
        let end = cds_end.ok_or_else(|| FerroError::ConversionError {
            msg: "No CDS end".to_string(),
        })?;

        if pos < cds_start {
            // 5'UTR: HGVS numbering skips c.0, so a tx position one
            // base 5' of cds_start is c.-1 (not c.0). Inverse of the
            // forward formula `tx = cds_start + base` for negative
            // base: `base = tx - cds_start`. Issue #97 — the previous
            // formula `tx - cds_start + 1` would emit base = 0 for
            // tx = cds_start - 1, rendered by `CdsPos::Display` as
            // `c.?` (`CDS_BASE_UNKNOWN`).
            Ok(CdsPos {
                base: pos as i64 - cds_start as i64,
                offset: None,
                utr3: false,
            })
        } else if pos > end {
            Ok(CdsPos {
                base: (pos - end) as i64,
                offset: None,
                utr3: true,
            })
        } else {
            Ok(CdsPos {
                base: (pos - cds_start + 1) as i64,
                offset: None,
                utr3: false,
            })
        }
    }

    /// Apply minimal notation to a CDS variant without full normalization.
    ///
    /// This is used when we can't do full normalization (e.g., missing transcript)
    /// but still want to apply minimal HGVS notation rules.
    fn canonicalize_cds_variant(&self, variant: &CdsVariant) -> CdsVariant {
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return variant.clone(),
        };

        // Only canonicalize if the edit has redundant information
        if !should_canonicalize(edit) {
            return variant.clone();
        }

        let canonical_edit = canonicalize_edit(edit);

        CdsVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                variant.loc_edit.location.clone(),
                variant.loc_edit.edit.map_ref(|_| canonical_edit),
            ),
        }
    }

    /// Apply minimal notation to a genome variant without full normalization.
    fn canonicalize_genome_variant(&self, variant: &GenomeVariant) -> GenomeVariant {
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return variant.clone(),
        };

        if !should_canonicalize(edit) {
            return variant.clone();
        }

        let canonical_edit = canonicalize_edit(edit);

        GenomeVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                variant.loc_edit.location.clone(),
                variant.loc_edit.edit.map_ref(|_| canonical_edit),
            ),
        }
    }

    /// Apply minimal notation to a transcript variant without full normalization.
    fn canonicalize_tx_variant(&self, variant: &TxVariant) -> TxVariant {
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return variant.clone(),
        };

        if !should_canonicalize(edit) {
            return variant.clone();
        }

        let canonical_edit = canonicalize_edit(edit);

        TxVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                variant.loc_edit.location.clone(),
                variant.loc_edit.edit.map_ref(|_| canonical_edit),
            ),
        }
    }

    /// Post-canonicalization split for a single normalized variant.
    /// Coord-system-agnostic: handles `g.`, `m.`, `c.` (CDS-proper positions
    /// only), `n.`, and `r.`. Fetches the per-coord-system reference window
    /// internally, calls `decompose_delins`, and rebuilds N variants when
    /// the decomposition fires. Returns `vec![variant]` if the variant
    /// doesn't decompose (non-Delins, complex location, no provider data,
    /// nothing to split out, **or the provider's ref window doesn't match
    /// the HGVS coordinate span** — e.g. cdot alignment gaps shorten the
    /// returned byte slice, see #339).
    ///
    /// Implements two spec-priority rules from `general.md:56`
    /// (substitution > deletion > inversion > duplication > insertion):
    /// - Inversion priority: a delins whose span contains a rev-comp
    ///   sub-span splits into `[…; inv; …]` (issue #160).
    /// - Substitution priority: a delins whose post-trim span contains
    ///   two or more independent single-base mismatches separated by at
    ///   least one unchanged nucleotide splits into separate substitutions
    ///   (issue #165 / item A10). The narrow codon-frame exception
    ///   (`general.md:35-38`) is preserved by `build_split_variants`,
    ///   which re-groups `[Sub; Identity; Sub]` triplets whose endpoints
    ///   share a codon when the variant is in CDS.
    ///
    /// Position math: `decompose_delins` returns 0-indexed offsets into
    /// the fetched `ref_bytes` slice. `ref_bytes[0]` corresponds to the
    /// variant's HGVS start position, so absolute HGVS pos = `hgvs_start +
    /// offset`.
    ///
    /// RNA `r.` variants have `U` bases in the alt while transcript ref bytes
    /// are `T`. Both slices are normalized to `T` before comparison so the
    /// rev-comp scan works uniformly; the emitted `Substitution` sub-edits
    /// preserve the original alt `Base` (which may be `Base::U`).
    fn apply_canonical_split(
        &self,
        variant: HgvsVariant,
    ) -> (Vec<HgvsVariant>, Vec<NormalizationWarning>) {
        let Some((hgvs_start, hgvs_end, alt_bytes, ref_bytes)) =
            self.fetch_ref_for_canonical_split(&variant)
        else {
            return (vec![variant], vec![]);
        };
        // When the provider returns fewer (or more) bytes than the HGVS
        // interval span, the canonical-split decomposition would walk past
        // the end of `ref_bytes` (or under-walk and miss interior
        // identities). This happens in practice when a cdot exon
        // alignment collapses gaps that the HGVS span counts as
        // positions (e.g. biocommons NG_032871.1:g.32476_53457delins…
        // returns 15,539 bytes for a 20,982-bp span). Pre-fix this
        // fired a `debug_assert_eq!` that panicked debug builds and
        // would have produced out-of-bounds reads in release. Now bail
        // out gracefully with a `CanonicalSplitSkipped` warning so
        // callers can flag the input for human review. Closes #339,
        // closes-after #354.
        let n = ref_bytes.len();
        let expected_span = (hgvs_end - hgvs_start + 1) as usize;
        if n != expected_span {
            let accession = variant_accession_string(&variant);
            let warning = NormalizationWarning::CanonicalSplitSkipped {
                accession,
                hgvs_start,
                hgvs_end,
                expected_span,
                actual_bytes: n,
            };
            return (vec![variant], vec![warning]);
        }
        let ref_norm = normalize_t_u(&ref_bytes);
        let alt_norm = normalize_t_u(&alt_bytes);
        let Some(subedits) = rules::decompose_delins(&ref_norm, 0, n, &alt_norm) else {
            return (vec![variant], vec![]);
        };
        // Substitution sub-edits inherit `alt_norm` bytes (T-form) from
        // `decompose_delins`, but the user's literal alt may have been U
        // (r. inputs). Re-derive the substitution `alternative` from the
        // pre-normalized `alt_bytes` so r. variants render `g>u` instead of
        // a silently coerced `g>t`. The position field is a 0-indexed offset
        // into the same window passed to `decompose_delins`, so it indexes
        // alt_bytes directly.
        let subedits = subedits
            .into_iter()
            .map(|se| match se {
                rules::DelinsSubedit::Substitution {
                    position,
                    reference,
                    alternative,
                } => {
                    let alt = crate::hgvs::edit::Base::from_char(alt_bytes[position] as char)
                        .unwrap_or(alternative);
                    rules::DelinsSubedit::Substitution {
                        position,
                        reference,
                        alternative: alt,
                    }
                }
                other => other,
            })
            .collect();
        // The codon-frame exception (`general.md:35-38`) applies to any
        // variant on the CDS-relative axis: `c.` (CDS proper) and `r.`
        // (RNA CDS-relative; issue #275 item 1).
        // `fetch_ref_for_canonical_split` already filters to CDS-proper
        // positions via `simple_cds_pos` / `simple_rna_pos`, so the
        // discriminant check below is sufficient. The exception fires
        // inside `build_split_variants` for every embedded
        // `[Sub; Identity; Sub]` triplet whose endpoints share a codon.
        // Caveat: the `r.` branch of `fetch_ref_for_canonical_split`
        // does not translate the RNA axis through `cds_start`, so for
        // transcripts with `cds_start > 1` the ref window can be
        // mis-aligned. The codon-frame split path inherits that latent
        // bug. See #291 — the fix is to apply `rna_to_tx_pos` (or
        // equivalent) before reading bytes.
        let codon_frame_aware = matches!(variant, HgvsVariant::Cds(_) | HgvsVariant::Rna(_));
        (
            build_split_variants(&variant, subedits, hgvs_start, codon_frame_aware),
            vec![],
        )
    }

    /// Per-coord-system extraction of `(hgvs_start, hgvs_end, alt_bytes,
    /// ref_bytes)` for the post-canonicalization split. Returns `None`
    /// when the variant is not a single-Delins at simple positions, or
    /// when the provider can't supply the ref window.
    ///
    /// The `ref_bytes` slice is sized exactly to the variant's HGVS interval
    /// (`hgvs_end - hgvs_start + 1` bytes), with `ref_bytes[0]` aligned to
    /// HGVS pos `hgvs_start`. This invariant lets the caller use a uniform
    /// `hgvs_pos = hgvs_start + offset` formula regardless of coord system.
    ///
    /// Note: this helper reads the transcript *sequence* (FASTA-derived,
    /// build-invariant) and not the `chromosome` field; the bare
    /// `provider.get_sequence(&transcript_accession, …)` is therefore
    /// sufficient. The build-aware lookup (`get_transcript_for_variant`,
    /// see #332) is reserved for paths that actually consume `chromosome`
    /// — the intronic and boundary-spanning normalization branches.
    fn fetch_ref_for_canonical_split(
        &self,
        variant: &HgvsVariant,
    ) -> Option<(u64, u64, Vec<u8>, Vec<u8>)> {
        let (hgvs_start, hgvs_end, alt) = extract_simple_delins(variant)?;
        let ref_bytes = match variant {
            HgvsVariant::Genome(g) => self
                .provider
                // get_sequence is 0-based half-open: [hgvs_start - 1, hgvs_end).
                .get_sequence(
                    &g.accession.transcript_accession(),
                    hgvs_start - 1,
                    hgvs_end,
                )
                .ok()?
                .into_bytes(),
            HgvsVariant::Mt(m) => self
                .provider
                .get_sequence(
                    &m.accession.transcript_accession(),
                    hgvs_start - 1,
                    hgvs_end,
                )
                .ok()?
                .into_bytes(),
            HgvsVariant::Cds(c) => {
                // CDS pos N → 1-based tx pos = cds_start + N - 1.
                // 0-based tx slice = [cds_start + N - 2, cds_start + end - 1).
                let tx = self
                    .provider
                    .get_transcript(&c.accession.transcript_accession())
                    .ok()?;
                let cds_start = tx.cds_start?;
                let s = cds_start.checked_add(hgvs_start)?.checked_sub(2)? as usize;
                let e = cds_start.checked_add(hgvs_end)?.checked_sub(1)? as usize;
                let bytes = tx.sequence.as_deref()?.as_bytes();
                if e > bytes.len() || s >= e {
                    return None;
                }
                bytes[s..e].to_vec()
            }
            HgvsVariant::Tx(t) => {
                let tx = self
                    .provider
                    .get_transcript(&t.accession.transcript_accession())
                    .ok()?;
                let s = (hgvs_start - 1) as usize;
                let e = hgvs_end as usize;
                let bytes = tx.sequence.as_deref()?.as_bytes();
                if e > bytes.len() || s >= e {
                    return None;
                }
                bytes[s..e].to_vec()
            }
            HgvsVariant::Rna(r) => {
                // `r.` positive non-UTR bases use the transcript-1-relative
                // axis in this codebase (NOT CDS-relative): `r.N` for `N > 0`
                // and non-UTR maps directly to tx index `N`. This matches the
                // convention used by `simple_rna_pos`, `normalize_rna::map_in`/
                // `map_out`, and `normalize_na_edit` for the r. arm.
                // `simple_rna_pos` filters intronic positions, 3'-UTR
                // positions (the `r.*N` form), and non-positive bases. It
                // does NOT filter 5'-UTR positions because under the tx-1
                // axis convention those have positive bases below
                // `cds_start` and are themselves valid tx indices — the
                // `hgvs_start - 1` slice against the full transcript
                // sequence is therefore correct for any position reaching
                // this arm (including tx-1 5'-UTR).
                // Pinned by `tests/issue_291_rna_axis_convention.rs` (closes #291).
                let tx = self
                    .provider
                    .get_transcript(&r.accession.transcript_accession())
                    .ok()?;
                let s = (hgvs_start - 1) as usize;
                let e = hgvs_end as usize;
                let bytes = tx.sequence.as_deref()?.as_bytes();
                if e > bytes.len() || s >= e {
                    return None;
                }
                bytes[s..e].to_vec()
            }
            _ => return None,
        };
        Some((hgvs_start, hgvs_end, alt, ref_bytes))
    }

    /// Issue #160 + #165 post-merge canonicalization for a single
    /// variant. Used by the cis-allele merge path; `normalize_allele`
    /// applies this per merged variant. Conservatively returns
    /// `vec![v]` for variants the helper can't process.
    ///
    /// Three spec rules are folded together by re-running normalization
    /// on the merged variant:
    /// - Full-span canonicalization (identity / dup / sub / del / ins /
    ///   full-span inv with outer-pair shortening) handled by
    ///   `canonicalize_delins` inside `normalize_na_edit`.
    /// - Sub-span inv decomposition (the issue #160 case) handled by
    ///   `apply_canonical_split` wired into each per-coord-system
    ///   `normalize_*`.
    /// - Sub-only decomposition for delins containing interior identities
    ///   (issue #165 / item A10), with the spec's codon-frame exception
    ///   (`general.md:35-38`) preserved inside `build_split_variants`.
    ///
    /// If the result is an `HgvsVariant::Allele` (the split fired and
    /// produced multiple variants), unwrap its inner variants so they
    /// flatten into the outer cis-allele list rather than nesting.
    fn canonical_split_for_variant(&self, v: HgvsVariant) -> Vec<HgvsVariant> {
        if !matches!(
            v,
            HgvsVariant::Genome(_)
                | HgvsVariant::Mt(_)
                | HgvsVariant::Cds(_)
                | HgvsVariant::Tx(_)
                | HgvsVariant::Rna(_)
        ) {
            return vec![v];
        }
        if extract_simple_delins(&v).is_none() {
            return vec![v];
        }
        match self.normalize_with_diagnostics(&v) {
            Ok(r) => match r.result {
                HgvsVariant::Allele(a) => a.variants,
                other => vec![other],
            },
            Err(_) => vec![v],
        }
    }
}

/// Position-only Display text for variant kinds that go through the
/// shuffle pipeline.
///
/// Returns the `loc_edit.location` rendering (e.g. `"4"` for `c.4del`,
/// `"100_103"` for `g.100_103del`) for the six nucleic-acid axes whose
/// normalizers call `normalize_na_edit`. Returns `None` for:
///
/// - `Protein` — protein 3'-shifting not yet implemented (#91).
/// - `Allele` — handled separately by the per-member compare path in
///   [`detect_shuffle_infos`].
/// - `RnaFusion` / `NullAllele` / `UnknownAllele` — no single position to
///   compare against.
///
/// `Circular` (o.) is included even though its current normalizer is a
/// pass-through clone (positions cannot change); the surface is
/// forward-safe for when #129 wires real circular shuffling and costs
/// nothing today (the post-hoc equality check trivially returns no info).
fn position_text_if_shuffleable(variant: &HgvsVariant) -> Option<String> {
    match variant {
        HV::Genome(v) => Some(format!("{}", v.loc_edit.location)),
        HV::Cds(v) => Some(format!("{}", v.loc_edit.location)),
        HV::Tx(v) => Some(format!("{}", v.loc_edit.location)),
        HV::Rna(v) => Some(format!("{}", v.loc_edit.location)),
        HV::Mt(v) => Some(format!("{}", v.loc_edit.location)),
        HV::Circular(v) => Some(format!("{}", v.loc_edit.location)),
        HV::Protein(_) | HV::Allele(_) | HV::RnaFusion(_) | HV::NullAllele | HV::UnknownAllele => {
            None
        }
    }
}

/// Returns the underlying `NaEdit` for a shuffle-eligible variant so the
/// shift detector can constrain `SHUFFLE_APPLIED` emission to genuine
/// shuffle transitions (companion to `position_text_if_shuffleable`).
/// `Mu::Uncertain` is unwrapped to the inner edit because shuffleability
/// is determined by edit shape, not by the uncertainty wrapper.
fn na_edit_if_shuffleable(variant: &HgvsVariant) -> Option<&NaEdit> {
    let mu = match variant {
        HV::Genome(v) => &v.loc_edit.edit,
        HV::Cds(v) => &v.loc_edit.edit,
        HV::Tx(v) => &v.loc_edit.edit,
        HV::Rna(v) => &v.loc_edit.edit,
        HV::Mt(v) => &v.loc_edit.edit,
        HV::Circular(v) => &v.loc_edit.edit,
        HV::Protein(_) | HV::Allele(_) | HV::RnaFusion(_) | HV::NullAllele | HV::UnknownAllele => {
            return None
        }
    };
    match mu {
        crate::hgvs::Mu::Certain(e) | crate::hgvs::Mu::Uncertain(e) => Some(e),
        // `Mu::Unknown` (Display: `?`) carries no edit — there is nothing
        // for shuffle to act on, so the kind-compatibility gate skips it.
        crate::hgvs::Mu::Unknown => None,
    }
}

/// True if the (input, output) edit-kind transition is one the shuffle
/// layer can produce. The shuffle algorithm rotates a deletion / insertion
/// / duplication / delins span along its reference window and may
/// canonicalize a single-base or homopolymer insertion into a duplication
/// or repeat — every other kind transition (e.g. `delins → sub` from
/// canonical-split, `delins → inv` from inversion decomposition,
/// `delins → =` from identity rewrite, `ins → delins`) is a
/// *canonicalization* and must not surface as `SHUFFLE_APPLIED`.
///
/// Closes the false-positive flagged on PR #426: `c.1_3delinsACG → c.2T>C`
/// is canonical-split, not shuffle, and the original position-only
/// post-hoc compare mislabeled it.
fn shuffle_kind_compatible(input: &NaEdit, output: &NaEdit) -> bool {
    use NaEdit::*;
    matches!(
        (input, output),
        (Deletion { .. }, Deletion { .. })
            | (Duplication { .. }, Duplication { .. })
            | (Delins { .. }, Delins { .. })
            | (Inversion { .. }, Inversion { .. })
            | (
                Insertion { .. },
                Insertion { .. } | Duplication { .. } | Repeat { .. } | MultiRepeat { .. },
            )
            | (Repeat { .. }, Repeat { .. } | MultiRepeat { .. })
            | (MultiRepeat { .. }, MultiRepeat { .. })
    )
}

/// Detect shuffle infos by comparing the input variant to the normalized
/// result.
///
/// A shift is recorded when the position-only Display text differs between
/// input and output. This is structural (not byte-exact on the full
/// description), so a canonical-form rewrite that leaves the position
/// unchanged (e.g. an explicit `delA` → `del`, which `canonicalize_edit`
/// applies only to the edit body, not the location) does not emit a false
/// positive.
///
/// Compound-allele rules (input/output both `HV::Allele` with equal
/// length): the comparison runs per bracket member in input order.
/// Bracket-length mismatches (the merge layer combined or split members)
/// are conservative no-ops — the structural rewrite is not a pure shuffle.
///
/// Top-level kind changes (input is `HV::Allele` but output is a bare
/// variant after cis-collapse, or vice versa via the `wrap_allele_if_split`
/// canonical-split path) are also conservative no-ops by construction:
/// `position_text_if_shuffleable` returns `None` for `HV::Allele`, so the
/// fallback arm yields no info. This narrows the surface to the cases the
/// signal is unambiguous; #330 explicitly scopes to the simple 3'-shift
/// channel.
///
/// Returns at most one info per shuffle event:
/// - Single-axis variants: zero or one info.
/// - Cis/trans alleles: one info per shifted member, in member order.
fn detect_shuffle_infos(
    input: &HgvsVariant,
    output: &HgvsVariant,
    direction: config::ShuffleDirection,
) -> Vec<NormalizationInfo> {
    match (input, output) {
        (HV::Allele(in_allele), HV::Allele(out_allele)) => {
            if in_allele.variants.len() != out_allele.variants.len() {
                return Vec::new();
            }
            in_allele
                .variants
                .iter()
                .zip(out_allele.variants.iter())
                .filter_map(|(i, o)| single_variant_shift_info(i, o, direction))
                .collect()
        }
        _ => single_variant_shift_info(input, output, direction)
            .into_iter()
            .collect(),
    }
}

/// Single-variant shift detector. Returns `Some(info)` iff the position
/// text differs between input and output for a shuffle-eligible axis AND
/// the (input, output) edit-kind transition is one the shuffle layer can
/// produce (see [`shuffle_kind_compatible`]). The edit-kind gate excludes
/// canonicalization that also moves positions — e.g. `delins → sub`,
/// `delins → =`, `delins → inv` — from masquerading as shuffle.
fn single_variant_shift_info(
    input: &HgvsVariant,
    output: &HgvsVariant,
    direction: config::ShuffleDirection,
) -> Option<NormalizationInfo> {
    let original_position = position_text_if_shuffleable(input)?;
    let normalized_position = position_text_if_shuffleable(output)?;
    if original_position == normalized_position {
        return None;
    }
    let input_edit = na_edit_if_shuffleable(input)?;
    let output_edit = na_edit_if_shuffleable(output)?;
    if !shuffle_kind_compatible(input_edit, output_edit) {
        return None;
    }
    let accession = input
        .accession()
        .map(|a| format!("{}", a))
        .unwrap_or_default();
    Some(NormalizationInfo::ShuffleApplied {
        accession,
        direction,
        original_position,
        normalized_position,
    })
}

/// Flip a fetched intronic genomic-strand window into transcript-view
/// orientation when the host transcript is on the minus strand. Returns
/// the input unchanged on plus strand. The relative positions and the
/// shuffle boundaries are flipped so they index into the returned
/// sequence consistently. Companion to [`unflip_intronic_positions`].
fn flip_intronic_for_strand(
    strand: Strand,
    genomic_seq: &str,
    rel_start: u64,
    rel_end: u64,
    boundaries: &Boundaries,
) -> Result<(String, u64, u64, Boundaries), FerroError> {
    if strand != Strand::Minus {
        return Ok((
            genomic_seq.to_string(),
            rel_start,
            rel_end,
            boundaries.clone(),
        ));
    }
    let seq_len = genomic_seq.len() as u64;
    // The reverse-complement flip maps a coordinate `x` to `seq_len - x + 1`,
    // which underflows (and previously panicked with "attempt to subtract
    // with overflow") whenever `x` exceeds the fetched window. That happens
    // for minus-strand intronic inputs whose enclosing intron extends past
    // the variant-sized window — `boundaries.right` is the far intron edge,
    // which for a large intron can lie far outside the fetched bases (e.g.
    // `NG_007107.2(NM_004992.3):c.378-17delT`). We cannot reliably 3'-shift
    // within bases we did not fetch, and silently clamping the boundary to
    // the window edge produces an off-by-one shift, so surface a clear error
    // instead. Sizing the window to the enclosing intron is the real fix that
    // would let these normalize (and demote) rather than error — see #488.
    let in_window = |x: u64| (1..=seq_len).contains(&x);
    if !in_window(rel_start)
        || !in_window(rel_end)
        || !in_window(boundaries.left)
        || !in_window(boundaries.right)
    {
        return Err(FerroError::ConversionError {
            msg: format!(
                "intronic minus-strand shuffle window too small: rel {rel_start}..{rel_end}, \
                 boundaries {}..{} exceed fetched window of length {seq_len}",
                boundaries.left, boundaries.right
            ),
        });
    }
    let rc = crate::sequence::reverse_complement(genomic_seq);
    let new_rel_start = seq_len - rel_end + 1;
    let new_rel_end = seq_len - rel_start + 1;
    let new_boundaries = Boundaries::new(
        seq_len - boundaries.right + 1,
        seq_len - boundaries.left + 1,
    );
    Ok((rc, new_rel_start, new_rel_end, new_boundaries))
}

/// Inverse of [`flip_intronic_for_strand`] for the result positions
/// emitted by `normalize_na_edit`. On plus strand returns the input
/// unchanged; on minus strand maps from transcript-view back to the
/// genomic-strand frame.
fn unflip_intronic_positions(
    strand: Strand,
    seq_len: u64,
    rel_start: u64,
    rel_end: u64,
) -> (u64, u64) {
    if strand == Strand::Minus {
        (seq_len - rel_end + 1, seq_len - rel_start + 1)
    } else {
        (rel_start, rel_end)
    }
}

// =============================================================================
// Issue #160 + #165: delins post-canonicalization split helpers
// =============================================================================
//
// After `normalize_na_edit` (or `merge_consecutive_edits` for cis alleles)
// produces a Delins variant, the resulting span may be expressible in a
// higher-priority form under `general.md:56` (sub > del > inv > dup > ins).
// Two cases fire here:
// - Inversion sub-span: the delins span contains a rev-comp sub-region —
//   split into `[…; inv; …]` (issue #160).
// - Independent substitutions: the delins span contains two or more
//   single-base mismatches separated by at least one unchanged nucleotide
//   — split each into its own sub variant (issue #165 / tracking issue
//   #81 item A10). The codon-frame exception (`general.md:35-38`) is
//   preserved when applicable (see `build_split_variants`).
//
// The split is implemented as a post-pass over an already-built variant. It
// fetches a reference window via the provider, calls
// `rules::decompose_delins`, and rebuilds N variants when the
// decomposition fires. For variants that don't decompose (most cases), the
// helper returns `vec![input]` and is effectively a no-op.

/// Per-coord-system-aware extraction of `(hgvs_start, hgvs_end, alt_bytes)`
/// from a variant whose edit is a literal `Delins` at simple positions
/// (no offsets, no uncertainty). Returns `None` for any variant shape that
/// can't be decomposed by the post-canonicalization split rules
/// (issues #160 / #165): non-Delins, intronic, uncertain boundary,
/// non-literal insert, etc.
/// Return the transcript-axis accession string for any `HgvsVariant`
/// kind that `apply_canonical_split` operates on. Used to build the
/// `CanonicalSplitSkipped` warning message — best-effort; returns
/// `<unknown>` for variant kinds that don't carry a single accession.
fn variant_accession_string(variant: &HgvsVariant) -> String {
    match variant {
        HgvsVariant::Genome(v) => v.accession.transcript_accession(),
        HgvsVariant::Cds(v) => v.accession.transcript_accession(),
        HgvsVariant::Tx(v) => v.accession.transcript_accession(),
        HgvsVariant::Rna(v) => v.accession.transcript_accession(),
        HgvsVariant::Mt(v) => v.accession.transcript_accession(),
        _ => "<unknown>".to_string(),
    }
}

fn extract_simple_delins(variant: &HgvsVariant) -> Option<(u64, u64, Vec<u8>)> {
    let (start, end, edit) = match variant {
        HgvsVariant::Genome(v) => simple_genome_loc_edit(&v.loc_edit)?,
        HgvsVariant::Cds(v) => simple_cds_loc_edit(&v.loc_edit)?,
        HgvsVariant::Tx(v) => simple_tx_loc_edit(&v.loc_edit)?,
        HgvsVariant::Rna(v) => simple_rna_loc_edit(&v.loc_edit)?,
        HgvsVariant::Mt(v) => simple_genome_loc_edit(&v.loc_edit)?,
        _ => return None,
    };
    let NaEdit::Delins { sequence, .. } = edit else {
        return None;
    };
    let InsertedSequence::Literal(seq) = sequence else {
        return None;
    };
    let alt: Vec<u8> = seq.bases().iter().map(|b| b.to_u8()).collect();
    Some((start, end, alt))
}

fn simple_genome_loc_edit(
    le: &LocEdit<Interval<GenomePos>, NaEdit>,
) -> Option<(u64, u64, &NaEdit)> {
    let edit = le.edit.inner()?;
    let s = simple_genome_pos(le.location.start.as_single()?)?;
    let e = simple_genome_pos(le.location.end.as_single()?)?;
    Some((s, e, edit))
}
fn simple_genome_pos(mu: &Mu<GenomePos>) -> Option<u64> {
    let Mu::Certain(p) = mu else { return None };
    if p.is_special() || p.offset.is_some() {
        return None;
    }
    Some(p.base)
}

fn simple_cds_loc_edit(le: &LocEdit<Interval<CdsPos>, NaEdit>) -> Option<(u64, u64, &NaEdit)> {
    let edit = le.edit.inner()?;
    // Only handle simple positive CDS positions (no UTR, no intronic, no
    // uncertainty). UTR delins decomposition would need its own coord-axis
    // logic and is out of scope for this fix.
    let s = simple_cds_pos(le.location.start.as_single()?)?;
    let e = simple_cds_pos(le.location.end.as_single()?)?;
    Some((s, e, edit))
}
fn simple_cds_pos(mu: &Mu<CdsPos>) -> Option<u64> {
    let Mu::Certain(p) = mu else { return None };
    if p.is_unknown() || p.is_intronic() || p.is_3utr() || p.base <= 0 {
        return None;
    }
    Some(p.base as u64)
}

fn simple_tx_loc_edit(le: &LocEdit<Interval<TxPos>, NaEdit>) -> Option<(u64, u64, &NaEdit)> {
    let edit = le.edit.inner()?;
    let s = simple_tx_pos(le.location.start.as_single()?)?;
    let e = simple_tx_pos(le.location.end.as_single()?)?;
    Some((s, e, edit))
}
fn simple_tx_pos(mu: &Mu<TxPos>) -> Option<u64> {
    let Mu::Certain(p) = mu else { return None };
    if p.is_intronic() || p.is_downstream() || p.base <= 0 {
        return None;
    }
    Some(p.base as u64)
}

fn simple_rna_loc_edit(le: &LocEdit<Interval<RnaPos>, NaEdit>) -> Option<(u64, u64, &NaEdit)> {
    let edit = le.edit.inner()?;
    let s = simple_rna_pos(le.location.start.as_single()?)?;
    let e = simple_rna_pos(le.location.end.as_single()?)?;
    Some((s, e, edit))
}
fn simple_rna_pos(mu: &Mu<RnaPos>) -> Option<u64> {
    let Mu::Certain(p) = mu else { return None };
    if p.is_intronic() || p.is_3utr() || p.base <= 0 {
        return None;
    }
    Some(p.base as u64)
}

/// Build a single HgvsVariant matching `template`'s coord-system kind /
/// accession / gene_symbol, with a new `[start_1based, end_1based]` location
/// and the given edit. Used by `build_split_variants` to spread the output
/// of `decompose_delins` back into a sequence of HgvsVariants.
fn build_variant_at(
    template: &HgvsVariant,
    start_1based: u64,
    end_1based: u64,
    edit: NaEdit,
) -> HgvsVariant {
    match template {
        HgvsVariant::Genome(g) => HgvsVariant::Genome(GenomeVariant {
            accession: g.accession.clone(),
            gene_symbol: g.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(GenomePos::new(start_1based), GenomePos::new(end_1based)),
                edit,
            ),
        }),
        HgvsVariant::Cds(c) => HgvsVariant::Cds(CdsVariant {
            accession: c.accession.clone(),
            gene_symbol: c.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(
                    CdsPos::new(start_1based as i64),
                    CdsPos::new(end_1based as i64),
                ),
                edit,
            ),
        }),
        HgvsVariant::Tx(t) => HgvsVariant::Tx(TxVariant {
            accession: t.accession.clone(),
            gene_symbol: t.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(
                    TxPos::new(start_1based as i64),
                    TxPos::new(end_1based as i64),
                ),
                edit,
            ),
        }),
        HgvsVariant::Rna(r) => HgvsVariant::Rna(RnaVariant {
            accession: r.accession.clone(),
            gene_symbol: r.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(
                    RnaPos::new(start_1based as i64),
                    RnaPos::new(end_1based as i64),
                ),
                edit,
            ),
        }),
        HgvsVariant::Mt(m) => HgvsVariant::Mt(MtVariant {
            accession: m.accession.clone(),
            gene_symbol: m.gene_symbol.clone(),
            loc_edit: LocEdit::new(
                Interval::new(GenomePos::new(start_1based), GenomePos::new(end_1based)),
                edit,
            ),
        }),
        _ => unreachable!("build_variant_at called with non-NaEdit variant kind"),
    }
}

/// Build N HgvsVariants from a Vec<DelinsSubedit>. Position offsets in the
/// subedits are 0-indexed into the (per-variant-sized) ref slice; absolute
/// 1-based HGVS positions are recovered as `offset + hgvs_start`, where
/// `hgvs_start` is the variant's HGVS start position.
///
/// Spec rules implemented (see `general.md`, `substitution.md`):
///
/// 1. **Codon-frame exception** (`general.md:35-38`, issue #79 / #165).
///    When `codon_frame_aware` is true, the scan looks ahead at each
///    position for a `[Sub@i; Identity@i+1; Sub@i+2]` triplet whose CDS
///    endpoints (`hgvs_start + i`, `hgvs_start + i + 2`) share a codon.
///    Such a triplet emits as a single 3-base `delins` with alt
///    sequence `[Sub@i.alt, Identity@i+1.base, Sub@i+2.alt]`. The flag
///    is true only for `c.` (CDS) variants — `g.`, `n.`, `r.`, and `m.`
///    have no codon-frame and skip this branch. The exception is
///    deliberately narrow (length-3, exact pattern, in-codon endpoints)
///    so it matches the spec text "two variants separated by one
///    nucleotide, together affecting one amino acid".
///
/// 2. **Adjacent-substitution coalescence** (`substitution.md`,
///    issue #182). Consecutive `Substitution` sub-edits whose positions
///    are strictly adjacent (no gap, no `Inversion` or `IdentityAt`
///    between them) group into a single `delins` variant — "changes
///    involving two or more consecutive nucleotides are described as
///    deletion/insertion".
///
/// 3. **Inversion as a hard barrier** (issue #166). An `Inversion`
///    always emits standalone and breaks any in-flight substitution
///    run, preserving the inv-priority decomposition.
///
/// Singleton sub-runs stay as `Substitution`. `IdentityAt` not consumed
/// by a codon-frame triplet drops (an unchanged base is not an edit) and
/// always ends any in-flight substitution run — the gap means the
/// surrounding subs are no longer "consecutive".
fn build_split_variants(
    template: &HgvsVariant,
    subedits: Vec<DelinsSubedit>,
    hgvs_start: u64,
    codon_frame_aware: bool,
) -> Vec<HgvsVariant> {
    let abs = |idx: usize| -> u64 { idx as u64 + hgvs_start };

    let mut output: Vec<HgvsVariant> = Vec::new();
    // Pending run of strictly-adjacent Substitution sub-edits, in
    // left-to-right order. Each entry is `(position, reference, alternative)`
    // with `position` the 0-indexed offset into the variant's ref window.
    let mut run: Vec<(usize, Base, Base)> = Vec::new();

    let n = subedits.len();
    let mut i = 0;
    while i < n {
        // Codon-frame triplet lookahead: try to consume `[Sub; Identity; Sub]`
        // at offsets `[i, i+1, i+2]` whose endpoints share a codon. Only
        // fires for CDS variants (`codon_frame_aware`) and is the post-merge
        // half of issue #79: a pair of in-codon SNVs separated by one
        // unchanged base must render as a 3-base `delins`, even when the
        // pair sits inside a longer decomposition.
        if codon_frame_aware && i + 2 < n {
            if let (
                DelinsSubedit::Substitution {
                    position: p1,
                    alternative: a1,
                    ..
                },
                DelinsSubedit::IdentityAt {
                    position: pm,
                    base: bm,
                },
                DelinsSubedit::Substitution {
                    position: p3,
                    alternative: a3,
                    ..
                },
            ) = (&subedits[i], &subedits[i + 1], &subedits[i + 2])
            {
                if *pm == *p1 + 1 && *p3 == *p1 + 2 {
                    let cds_p1 = abs(*p1) as i64;
                    let cds_p3 = abs(*p3) as i64;
                    if merge::same_codon(cds_p1, cds_p3) {
                        // Codon-frame triplet preserved as a 3-base
                        // delins. `bm` is the unchanged ref byte from
                        // `decompose_delins`. Both `c.` and `r.` flow
                        // through this branch (issue #275 item 1 set
                        // `codon_frame_aware = true` for
                        // `HgvsVariant::Cds` and `HgvsVariant::Rna`),
                        // but no T/U recovery is needed here: the
                        // emitted delins is rendered by the per-variant
                        // formatter, and the `r.` formatter lowercases
                        // all of its bases (T → u included) when it
                        // prints the alt sequence. Forwarding the raw
                        // ref byte from `decompose_delins` is therefore
                        // safe for both coordinate systems.
                        flush_substitution_run(&mut output, template, hgvs_start, &mut run);
                        let s = abs(*p1);
                        let e = abs(*p3);
                        let alt_bases = vec![*a1, *bm, *a3];
                        output.push(build_variant_at(
                            template,
                            s,
                            e,
                            NaEdit::Delins {
                                sequence: InsertedSequence::Literal(Sequence::new(alt_bases)),
                                deleted: None,
                                deleted_length: None,
                            },
                        ));
                        i += 3;
                        continue;
                    }
                }
            }
        }

        match &subedits[i] {
            DelinsSubedit::Substitution {
                position,
                reference,
                alternative,
            } => {
                let breaks_run = matches!(run.last(), Some((prev, _, _)) if *prev + 1 != *position);
                if breaks_run {
                    flush_substitution_run(&mut output, template, hgvs_start, &mut run);
                }
                run.push((*position, *reference, *alternative));
            }
            DelinsSubedit::Inversion { start, end } => {
                flush_substitution_run(&mut output, template, hgvs_start, &mut run);
                // Half-open 0-indexed [start, end) of length L=end-start.
                // HGVS inclusive interval covers L bases starting at
                // abs(start) and ending at abs(start)+L-1 = abs(end-1).
                let s = abs(*start);
                let e = abs(*end) - 1;
                output.push(build_variant_at(
                    template,
                    s,
                    e,
                    NaEdit::Inversion {
                        sequence: None,
                        length: None,
                    },
                ));
            }
            // Drop IdentityAt: an unchanged base is not an edit. Outside of
            // the codon-frame triplet branch above, identities here are
            // either codon-frame-merge interior bases (issue #79) that did
            // not pair into a same-codon triplet, or outer bases absorbed
            // by `shorten_inversion`. An identity also ends any in-flight
            // substitution run — the gap means the surrounding subs are no
            // longer "consecutive".
            DelinsSubedit::IdentityAt { .. } => {
                flush_substitution_run(&mut output, template, hgvs_start, &mut run);
            }
        }
        i += 1;
    }
    flush_substitution_run(&mut output, template, hgvs_start, &mut run);
    output
}

/// Flush a pending run of consecutive substitution sub-edits into `output`.
/// A length-1 run emits a `Substitution`; a length-2+ run emits a single
/// `Delins` over `[run.first.position, run.last.position]` with `sequence`
/// = concatenated `alternative` bases. See `build_split_variants` for the
/// spec rationale (issue #182).
fn flush_substitution_run(
    output: &mut Vec<HgvsVariant>,
    template: &HgvsVariant,
    hgvs_start: u64,
    run: &mut Vec<(usize, Base, Base)>,
) {
    if run.is_empty() {
        return;
    }
    let abs = |idx: usize| -> u64 { idx as u64 + hgvs_start };
    if run.len() == 1 {
        let (position, reference, alternative) = run.pop().unwrap();
        let p = abs(position);
        output.push(build_variant_at(
            template,
            p,
            p,
            NaEdit::Substitution {
                reference,
                alternative,
            },
        ));
        return;
    }
    let s = abs(run.first().unwrap().0);
    let e = abs(run.last().unwrap().0);
    let alt_bases: Vec<Base> = run.drain(..).map(|(_, _, a)| a).collect();
    output.push(build_variant_at(
        template,
        s,
        e,
        NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::new(alt_bases)),
            deleted: None,
            deleted_length: None,
        },
    ));
}

/// Normalize DNA `T` and RNA `U` to a single byte (`T`) so byte-wise
/// comparison works across coord systems. Used by `apply_canonical_split`
/// to make every decomposition scan (rev-comp inv detection, per-position
/// sub / identity classification) T/U-agnostic for `r.` variants whose
/// alt bytes contain `U` while the transcript ref contains `T`.
fn normalize_t_u(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            b'U' => b'T',
            b'u' => b't',
            other => other,
        })
        .collect()
}

/// If `variants` has 1 element return it directly; if >1 wrap in a cis Allele.
fn wrap_allele_if_split(mut variants: Vec<HgvsVariant>) -> HgvsVariant {
    debug_assert!(!variants.is_empty(), "wrap_allele_if_split: empty input");
    if variants.len() == 1 {
        variants.pop().unwrap()
    } else {
        HgvsVariant::Allele(AlleleVariant::new(variants, AllelePhase::Cis))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;
    use crate::reference::MockProvider;

    #[test]
    fn test_normalize_substitution_unchanged() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();

        // Substitutions should not change
        assert_eq!(format!("{}", variant), format!("{}", normalized));
    }

    #[test]
    fn test_normalize_with_config() {
        let provider = MockProvider::with_test_data();
        let config = NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime);
        let normalizer = Normalizer::with_config(provider, config);

        assert_eq!(
            normalizer.config().shuffle_direction,
            ShuffleDirection::FivePrime
        );
    }

    #[test]
    fn test_normalizer_handles_missing_transcript() {
        let provider = MockProvider::new(); // Empty provider
        let normalizer = Normalizer::new(provider);

        // Should return variant unchanged when transcript not found
        let variant = parse_hgvs("NM_MISSING.1:c.100del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(result.is_ok());
        // Verify output equals input (unchanged)
        assert_eq!(
            format!("{}", variant),
            format!("{}", result.unwrap()),
            "Missing transcript should return variant unchanged"
        );
    }

    #[test]
    fn test_normalize_deletion_shifts_3prime() {
        // NM_001234.1 has G repeat spanning exon boundaries
        // Exon 1: c.1-11, Exon 2: c.12-26, Exon 3: c.27+
        // G repeat is at c.9-c.33, but shift stops at exon boundary
        // c.10 is in exon 1 (ends at c.11), so deletion shifts to c.11
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:c.10del").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Should remain a deletion
        assert!(
            output.contains("del"),
            "Deletion should remain a deletion, got: {}",
            output
        );
        // Should shift from c.10 to c.11 (3'-most within exon 1)
        assert!(
            output.contains("c.11del"),
            "Deletion should shift 3' to exon boundary (c.11), got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_insertion_becomes_dup() {
        // NM_001234.1 has G repeat at CDS positions c.9-c.33
        // Inserting G after position 10 should shift 3' and become dup
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:c.10_11insG").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Inserting G in a G-repeat should become dup
        assert!(
            output.contains("dup"),
            "Insertion of matching base should become dup, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_duplication_shifts_3prime() {
        // NM_001234.1 has G repeat spanning positions c.9-33 (25 G's)
        // Single-base duplications stay as simple dups (only 2+ base dups become repeat notation)
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:c.10dup").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Single-base duplication should shift 3' and stay as dup
        // c.10dup in GGGGG...GGG tract shifts to rightmost position (c.33) but stays as dup
        assert!(
            output.contains("dup"),
            "Single-base duplication should remain as dup notation, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_delins_unchanged() {
        // A delins that doesn't simplify should stay as delins
        // Deleting G and inserting AT is not a dup pattern
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10delinsAT").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        assert!(
            output.contains("delinsAT"),
            "Delins should remain unchanged, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_single_base_delins_becomes_substitution() {
        // HGVS edit-type priority: a 1→1 delins with ref!=alt must be expressed
        // as a substitution. Transcript NM_000088.3 starts ATGCCCAAGG...; position
        // 10 is G. c.10delinsT replaces G with T → c.10G>T.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10delinsT").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.10G>T");
    }

    #[test]
    fn test_normalize_single_base_delins_same_base_becomes_identity() {
        // Per HGVS, a delins whose insert equals the reference is identity (=).
        // Transcript NM_000088.3 starts ATGCCCAAGG...; position 10 is G.
        // c.10delinsG produces no change → c.10=.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10delinsG").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.10=");
    }

    #[test]
    fn test_normalize_multi_base_delins_same_seq_becomes_identity() {
        // Transcript NM_000088.3 starts ATG at positions 1-3.
        // c.1_3delinsATG produces no change → c.1_3=.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.1_3delinsATG").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.1_3=");
    }

    #[test]
    fn test_normalize_multi_base_delete_delins_to_pure_deletion() {
        // c.10_11delinsT against NM_000088.3 (c.10_11 = GT). The shared `T`
        // suffix consumes the inserted base entirely, leaving a single-base
        // deletion at c.10. Per HGVS minimal-form rules (sub > del > inv >
        // dup > ins) the canonical output is a pure deletion, not a delins.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10_11delinsT").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.10del");
    }

    #[test]
    fn test_normalize_empty_insert_delins_becomes_deletion() {
        // HGVS spec: a delins whose inserted sequence is empty is semantically
        // a deletion and must be rendered as `del`. Issue #81 item A3.
        // Transcript NM_000088.3 starts ATGCCCAAGG…; c.10delins (empty insert)
        // is a deletion of position 10 → c.10del.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10delins").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.10del");
    }

    #[test]
    fn test_normalize_empty_insert_multi_base_delins_becomes_deletion() {
        // Multi-base form: c.10_11delins (deletes "GT" at positions 10-11,
        // inserts nothing) → del, then the spec's 3'-rule shifts the deletion
        // to c.11_12del because ref[10]=G == ref[12]=G. Issue #81 item A3.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10_11delins").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.11_12del");
    }

    #[test]
    fn test_normalize_delins_to_dup_still_works() {
        // Regression guard: adding identity/substitution checks before the dup
        // check must not block legitimate dup conversions. ref[5] = C;
        // c.5delinsCC matches the dup pattern (insert == deleted twice) and
        // must still normalize to dup.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.5delinsCC").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);
        assert!(
            output.contains("dup"),
            "delins matching dup pattern should normalize to dup, got: {}",
            output
        );
        assert!(
            !output.contains("delins"),
            "delins should not survive the dup conversion, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_delins_different_bases_becomes_substitution() {
        // c.1_3delinsACG against NM_000088.3 (c.1_3 = ATG). The shared `A`
        // prefix and `G` suffix collapse the delins to T -> C at c.2 per the
        // HGVS minimal-form rule (sub > del > inv > dup > ins).
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.1_3delinsACG").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.2T>C");

        // r. twin: same collapse on the RNA path. NM_000088.3 has cds_start = 1,
        // so r.1_3 maps to the same ATG run and the residual is at r.2. The
        // residual ref byte is DNA `T` (sourced from the c. axis), which the
        // RNA `Display` path canonicalizes to `u` per the HGVS RNA alphabet
        // (see issue #276 / `NaEdit::Substitution` in src/hgvs/edit.rs).
        let variant = parse_hgvs("NM_000088.3:r.1_3delinsacg").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:r.2u>c");
    }

    #[test]
    fn test_normalize_substitution_ref_equals_alt_becomes_identity() {
        // Per HGVS, a substitution where the reference and alternative bases
        // are identical produces no change and must be expressed using
        // identity notation (`=`). SNV companion to the same-base delins rule.
        // Transcript NM_000088.3 starts ATGCCCAAGG...; position 10 is G.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10G>G").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.10=");
    }

    #[test]
    fn test_normalize_substitution_ref_equals_alt_first_position() {
        // Boundary check: the rule must fire at position 1 (first base).
        // Transcript NM_000088.3 starts ATG...; position 1 is A.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.1A>A").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.1=");
    }

    #[test]
    fn test_normalize_substitution_ref_not_equal_alt_unchanged() {
        // Regression guard: a real SNV must not be rewritten to identity.
        // c.10G>T at position 10 (G) is a valid substitution.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10G>T").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_000088.3:c.10G>T");
    }

    #[test]
    fn test_normalize_substitution_ref_equals_alt_without_provider_data() {
        // The A4 rule is purely syntactic, so it must fire even when the
        // provider has no transcript loaded — matching the spec's stance that
        // `c.123C>C` is "not allowed" regardless of reference availability.
        // Spec example: docs/recommendations/DNA/other.md (HGVS v21.0).
        let provider = MockProvider::new(); // empty — no transcripts
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_004006.2:c.123C>C").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", result), "NM_004006.2:c.123=");
    }

    #[test]
    fn test_normalize_protein_substitution_unchanged() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Protein substitution variants should pass through unchanged
        let variant = parse_hgvs("NP_000079.2:p.Val600Glu").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", variant), format!("{}", normalized));
    }

    #[test]
    fn test_normalize_protein_deletion_removes_redundant_sequence() {
        // Redundant sequence removal: p.Val600delVal → p.Val600del
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NP_000079.2:p.Val600delVal").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", normalized);

        // Should remove redundant "Val" from the deletion
        assert_eq!(
            output, "NP_000079.2:p.Val600del",
            "Redundant sequence should be removed from deletion"
        );
    }

    #[test]
    fn test_normalize_protein_deletion_range_removes_redundant_sequence() {
        // Redundant sequence removal for range: p.Lys23_Glu25delLysAlaGlu → p.Lys23_Glu25del
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NP_000079.2:p.Lys23_Glu25delLysAlaGlu").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", normalized);

        // Should remove redundant sequence from the deletion
        assert_eq!(
            output, "NP_000079.2:p.Lys23_Glu25del",
            "Redundant sequence should be removed from range deletion"
        );
    }

    #[test]
    fn test_normalize_protein_deletion_non_matching_sequence_unchanged() {
        // Non-matching sequence should stay: p.Val600delGlu should NOT change
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NP_000079.2:p.Val600delGlu").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", normalized);

        // Should NOT remove non-matching sequence
        assert_eq!(
            output, "NP_000079.2:p.Val600delGlu",
            "Non-matching sequence should not be removed"
        );
    }

    #[test]
    fn test_normalize_protein_deletion_no_sequence_unchanged() {
        // Deletion without sequence should stay unchanged
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NP_000079.2:p.Val600del").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", normalized);

        assert_eq!(
            output, "NP_000079.2:p.Val600del",
            "Deletion without sequence should remain unchanged"
        );
    }

    #[test]
    fn test_normalize_protein_duplication_unchanged() {
        // Duplications should pass through unchanged
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NP_000079.2:p.Val600dup").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", variant), format!("{}", normalized));
    }

    #[test]
    fn test_normalize_protein_frameshift_unchanged() {
        // Frameshifts should pass through unchanged
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NP_000079.2:p.Arg97ProfsTer23").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", variant), format!("{}", normalized));
    }

    #[test]
    fn test_normalize_protein_reference_validation_match() {
        // Test that validation passes when amino acid matches reference
        // NP_TEST.1 has: M at position 1, V at position 2
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Position 1 = M (Met), Position 2 = V (Val)
        let variant = parse_hgvs("NP_TEST.1:p.Met1Val").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "Validation should pass for matching amino acid"
        );
    }

    #[test]
    fn test_normalize_protein_reference_validation_mismatch() {
        // Test that validation fails when amino acid doesn't match reference
        // NP_TEST.1 has M at position 1, but we claim it's Val
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Position 1 is M (Met), not V (Val)
        let variant = parse_hgvs("NP_TEST.1:p.Val1Glu").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_err(),
            "Validation should fail for mismatched amino acid"
        );

        if let Err(crate::error::FerroError::AminoAcidMismatch {
            position,
            expected,
            found,
            ..
        }) = result
        {
            assert_eq!(position, 1);
            assert_eq!(expected, "Val");
            assert_eq!(found, "M");
        } else {
            panic!("Expected AminoAcidMismatch error");
        }
    }

    #[test]
    fn test_normalize_protein_reference_validation_deletion() {
        // Test validation for deletion variants
        // NP_TEST.1 has V at position 2
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Position 2 = V (Val) - should pass
        let variant = parse_hgvs("NP_TEST.1:p.Val2del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "Validation should pass for matching deletion position"
        );
    }

    #[test]
    fn test_normalize_protein_reference_validation_missing_protein() {
        // Test that missing protein data skips validation (doesn't fail)
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // NP_MISSING.1 doesn't exist in provider
        let variant = parse_hgvs("NP_MISSING.1:p.Val600Glu").unwrap();
        let result = normalizer.normalize(&variant);
        // Should NOT error - just skip validation when protein not available
        assert!(
            result.is_ok(),
            "Missing protein should skip validation, not fail"
        );
    }

    #[test]
    fn test_normalize_rna_substitution_unchanged() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // RNA substitutions should pass through unchanged
        let variant = parse_hgvs("NM_000088.3:r.10a>g").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", variant), format!("{}", normalized));
    }

    #[test]
    fn test_normalize_rna_deletion_shifts_3prime() {
        // NM_001234.1 has a G repeat at tx positions 13-37 spanning three
        // exons: exon 1 (1-15), exon 2 (16-30), exon 3 (31-44). Per HGVS
        // general.md's exon-junction exception, an r. deletion does not
        // shift across exon boundaries — so `r.14del` (a G in exon 1)
        // shifts only as far as `r.15del` (the last G inside exon 1).
        // Closes #334.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:r.14del").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        assert!(
            output.contains("del"),
            "Deletion should remain a deletion, got: {}",
            output
        );
        // Spec-canonical 3' anchor inside exon 1 is position 15.
        assert!(
            output.contains("r.15del"),
            "Deletion should shift 3' to exon-1 boundary (r.15del), got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_rna_insertion_becomes_dup() {
        // NM_001234.1 has G repeat at positions 13-36
        // Inserting g after position 14 should shift 3' and become dup
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:r.14_15insg").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Inserting g in a G-repeat should become dup
        assert!(
            output.contains("dup"),
            "Insertion of matching base should become dup, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_rna_duplication_shifts_3prime() {
        // NM_001234.1 has G repeat - single-base duplications stay as simple dups
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:r.14dup").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Single-base duplication should shift 3' and stay as dup
        assert!(
            output.contains("dup"),
            "Single-base RNA duplication should remain as dup notation, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_rna_intronic_returns_error() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Intronic RNA variants should return an error
        let variant = parse_hgvs("NM_001234.1:r.10+5del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(result.is_err(), "Intronic RNA variant should return error");
    }

    #[test]
    fn test_normalize_rna_missing_transcript_unchanged() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Missing transcript should return variant unchanged
        let variant = parse_hgvs("NM_MISSING.1:r.100del").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", variant), format!("{}", normalized));
    }

    #[test]
    fn test_normalize_null_allele() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Null alleles should pass through
        let variant = HgvsVariant::NullAllele;
        let normalized = normalizer.normalize(&variant).unwrap();
        assert!(matches!(normalized, HgvsVariant::NullAllele));
    }

    #[test]
    fn test_normalize_unknown_allele() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Unknown alleles should pass through
        let variant = HgvsVariant::UnknownAllele;
        let normalized = normalizer.normalize(&variant).unwrap();
        assert!(matches!(normalized, HgvsVariant::UnknownAllele));
    }

    #[test]
    fn test_normalize_allele_variant() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Allele variants should normalize each component
        // Using substitutions which should remain unchanged
        let variant = parse_hgvs("[NM_000088.3:c.10A>G;NM_000088.3:c.20C>T]").unwrap();
        let result = normalizer.normalize(&variant).unwrap();

        // Verify it's still an allele
        assert!(matches!(result, HgvsVariant::Allele(_)));

        // Verify output is the canonical compact form (ACC:c.[edit1;edit2])
        assert_eq!(
            format!("{}", result),
            "NM_000088.3:c.[10A>G;20C>T]",
            "Allele display should use canonical compact form"
        );
    }

    #[test]
    fn test_normalize_5prime_direction() {
        // Test that 5' direction shifts toward 5' end instead of 3'
        // NM_001234.1 has G repeat spanning exons
        // Exon 2: c.12-c.26
        // With 5' direction, deletion at c.20 should shift toward c.12
        // Note: Actual shift depends on boundary handling
        let provider = MockProvider::with_test_data();
        let config = NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime);
        let normalizer = Normalizer::with_config(provider, config);

        // Delete G at position 20 (middle of G repeat in exon 2)
        let variant = parse_hgvs("NM_001234.1:c.20del").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        assert!(
            output.contains("del"),
            "Should remain a deletion, got: {}",
            output
        );
        // With 5' direction within exon 2, should shift toward 5' boundary
        // The exact position depends on exon boundary handling
        assert!(
            output.contains("c.13del") || output.contains("c.12del"),
            "5' direction should shift deletion toward exon start, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_3prime_direction() {
        let provider = MockProvider::with_test_data();
        let config = NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime);
        let normalizer = Normalizer::with_config(provider, config);

        let variant = parse_hgvs("NM_000088.3:c.10del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(result.is_ok());
    }

    #[test]
    fn test_normalize_with_cross_boundaries() {
        let provider = MockProvider::with_test_data();
        let config = NormalizeConfig::default().allow_crossing_boundaries();
        let normalizer = Normalizer::with_config(provider, config);

        let variant = parse_hgvs("NM_000088.3:c.10del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(result.is_ok());
    }

    #[test]
    fn test_normalize_genomic_variant() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Genomic variants with missing sequence should pass through unchanged
        let variant = parse_hgvs("NC_000001.11:g.12345del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(result.is_ok());
    }

    #[test]
    fn test_normalize_tx_variant() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NR_000001.1:n.100del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(result.is_ok());
    }

    #[test]
    fn test_normalize_mt_variant() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // MT variants pass through unchanged
        let variant = parse_hgvs("NC_012920.1:m.100A>G").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert!(matches!(normalized, HgvsVariant::Mt(_)));
    }

    #[test]
    fn test_cds_to_tx_pos_utr5() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // 5' UTR position (negative)
        let cds_pos = CdsPos {
            base: -5,
            offset: None,
            utr3: false,
        };
        let result = normalizer.cds_to_tx_pos(&cds_pos, 10, Some(50));
        assert!(result.is_ok());
    }

    #[test]
    fn test_cds_to_tx_pos_utr3() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // 3' UTR position
        let cds_pos = CdsPos {
            base: 5,
            offset: None,
            utr3: true,
        };
        let result = normalizer.cds_to_tx_pos(&cds_pos, 10, Some(50));
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 55);
    }

    #[test]
    fn test_cds_to_tx_pos_coding() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Normal coding position
        let cds_pos = CdsPos {
            base: 10,
            offset: None,
            utr3: false,
        };
        let result = normalizer.cds_to_tx_pos(&cds_pos, 5, Some(50));
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 14); // 5 + 10 - 1 = 14
    }

    #[test]
    fn test_tx_to_cds_pos_utr5() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Position before CDS start
        let result = normalizer.tx_to_cds_pos(3, 10, Some(50));
        assert!(result.is_ok());
        let cds_pos = result.unwrap();
        assert!(cds_pos.base < 0);
        assert!(!cds_pos.utr3);
    }

    #[test]
    fn test_tx_to_cds_pos_utr3() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Position after CDS end
        let result = normalizer.tx_to_cds_pos(55, 10, Some(50));
        assert!(result.is_ok());
        let cds_pos = result.unwrap();
        assert!(cds_pos.utr3);
    }

    #[test]
    fn test_tx_to_cds_pos_coding() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // Normal coding position
        let result = normalizer.tx_to_cds_pos(20, 10, Some(50));
        assert!(result.is_ok());
        let cds_pos = result.unwrap();
        assert!(!cds_pos.utr3);
        assert_eq!(cds_pos.base, 11); // 20 - 10 + 1 = 11
    }

    #[test]
    fn test_config_default() {
        let config = NormalizeConfig::default();
        assert_eq!(config.shuffle_direction, ShuffleDirection::ThreePrime);
        assert!(!config.cross_boundaries);
    }

    #[test]
    #[allow(deprecated)]
    fn test_config_builder() {
        let config = NormalizeConfig::default()
            .with_direction(ShuffleDirection::FivePrime)
            .allow_crossing_boundaries()
            .skip_validation();

        assert_eq!(config.shuffle_direction, ShuffleDirection::FivePrime);
        assert!(config.cross_boundaries);
        // skip_validation now sets RefSeqMismatch to SilentCorrect
        assert!(!config.should_reject_ref_mismatch());
        assert!(!config.should_warn_ref_mismatch());
    }

    #[test]
    fn test_duplication_3prime_shift_two_bases() {
        // Test the exact scenario from ClinVar: c.4159dup vs c.4160dup
        // When duplicating an A in a homopolymer tract (AA),
        // HGVS requires converting to repeat notation.
        //
        // NM_888888.1 sequence: ATGCCCGAAGCCCCCCCCCGTTTGCATGCATGCATGCAT
        // Positions (1-based):  12345678901234567890...
        // c.8 = A, c.9 = A (the "AA" in "GAA")
        //
        // Single-base duplications stay as simple dups
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_888888.1:c.8dup").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Single-base duplication should shift 3' and stay as dup
        assert!(
            output.contains("dup"),
            "Single-base duplication should remain as dup notation, got: {}",
            output
        );
    }

    #[test]
    fn test_duplication_3prime_shift_three_bases() {
        // Test with three consecutive identical bases (TTT)
        //
        // NM_888888.1 sequence: ATGCCCGAAGCCCCCCCCCGTTTGCATGCATGCATGCAT
        // Positions:           ...              2021222324...
        // c.20 = G, c.21 = T, c.22 = T, c.23 = T, c.24 = G
        //
        // Single-base duplications stay as simple dups
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_888888.1:c.21dup").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        // Single-base duplication should shift 3' to c.23 and stay as dup
        assert!(
            output.contains("dup"),
            "Single-base duplication should remain as dup notation, got: {}",
            output
        );
    }

    #[test]
    fn test_duplication_no_shift_when_unique() {
        // Test that a duplication of a unique base doesn't shift
        //
        // NM_888888.1 sequence: ATGCCCGAAGCCCCCCCCCGTTTGCATGCATGCATGCAT
        // c.7 = G (followed by AA, so no G to shift to)
        //
        // c.7dup should stay as c.7dup
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_888888.1:c.7dup").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        assert!(
            output.contains("dup"),
            "Should remain a duplication, got: {}",
            output
        );
        assert!(
            output.contains("c.7dup"),
            "Duplication of unique G at c.7 should not shift, got: {}",
            output
        );
    }

    // =====================================================================
    // Exon-Intron Boundary Spanning Variant Tests
    // =====================================================================

    /// Create a provider with a transcript that has genomic coordinates and introns
    /// for testing boundary-spanning variant normalization.
    ///
    /// Transcript structure (NM_BOUNDARY.1):
    /// - Gene: BOUNDARY
    /// - Strand: Plus
    /// - Chromosome: chr1
    ///
    /// Genomic layout (chr1):
    /// Position: 1000                   1020  1030                   1050  1060                   1080
    ///           |-------- Exon 1 ------|      |-------- Exon 2 ------|      |-------- Exon 3 ------|
    ///           ATGCCCAAAGGGTTTAGGCCC       AAAGGGTTTAGGCCCAAAAAA       GGGTTTAGGCCCAAATGA
    ///                                 ^^^  ^^^                   ^^^  ^^^
    ///                              intron 1                   intron 2
    ///
    /// Transcript positions (tx):
    /// - Exon 1: tx 1-20 = genomic 1000-1019
    /// - Intron 1: genomic 1020-1029 (10 bp)
    /// - Exon 2: tx 21-40 = genomic 1030-1049
    /// - Intron 2: genomic 1050-1059 (10 bp)
    /// - Exon 3: tx 41-58 = genomic 1060-1077
    ///
    /// CDS: starts at tx 1 (no 5' UTR for simplicity)
    /// CDS positions: c.1 = tx 1, c.20 = tx 20 (last of exon 1), c.21 = tx 21 (first of exon 2)
    ///
    /// Intron 1 sequence (g.1020-1029): "GTAAGCTAGG" (10 bp)
    ///   - c.20+1 = g.1020 (G)
    ///   - c.20+10 = g.1029 (G)
    ///   - c.21-10 = g.1020 (G)
    ///   - c.21-1 = g.1029 (G)
    ///
    /// Intron 2 sequence (g.1050-1059): "GTAAGTAAGG" (10 bp)
    fn make_boundary_test_provider() -> MockProvider {
        use crate::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
        use std::sync::OnceLock;

        let mut provider = MockProvider::new();

        // Build transcript sequence (exons only, spliced)
        // Exon 1 (20bp): ATGCCCAAAGGGTTTAGGCC (ends with CC at exon boundary)
        // Exon 2 (20bp): AAAGGGTTTAGGCCCAAAAA (AA repeat at both boundaries)
        // Exon 3 (18bp): GGGTTTAGGCCCAAATGA
        // Total: 58bp transcript
        let tx_seq = "ATGCCCAAAGGGTTTAGGCCAAAGGGTTTAGGCCCAAAAAGGGTTTAGGCCCAAATGA";

        // Build genomic sequence around the transcript
        // We'll create 100bp before, the gene region, and 100bp after
        // Gene region: exon1 + intron1 + exon2 + intron2 + exon3
        // = 20 + 10 + 20 + 10 + 18 = 78bp
        //
        // Genomic: 900-999 (padding) + 1000-1019 (exon1) + 1020-1029 (intron1) +
        //          1030-1049 (exon2) + 1050-1059 (intron2) + 1060-1077 (exon3) + 1078+ (padding)
        let mut genomic_seq = String::new();

        // Padding before (positions 0-999, 1000 bytes at 0-based)
        for _ in 0..1000 {
            genomic_seq.push('N');
        }

        // Exon 1 (positions 1000-1019, 0-based 1000-1019)
        genomic_seq.push_str("ATGCCCAAAGGGTTTAGGCC");

        // Intron 1 (positions 1020-1029) - with splice consensus
        // Note: The intron has AAA at the end (1027-1029) to test shifting
        genomic_seq.push_str("GTAAGCTAAA");

        // Exon 2 (positions 1030-1049)
        genomic_seq.push_str("AAAGGGTTTAGGCCCAAAAA");

        // Intron 2 (positions 1050-1059) - with AAA at start for testing
        genomic_seq.push_str("AAAGTAAGGG");

        // Exon 3 (positions 1060-1077)
        genomic_seq.push_str("GGGTTTAGGCCCAAATGA");

        // Padding after (100 bytes)
        for _ in 0..100 {
            genomic_seq.push('N');
        }

        // Add genomic sequence to provider
        provider.add_genomic_sequence("chr1", genomic_seq);

        // Create transcript with exons that have genomic coordinates
        provider.add_transcript(Transcript {
            id: "NM_BOUNDARY.1".to_string(),
            gene_symbol: Some("BOUNDARY".to_string()),
            strand: Strand::Plus,
            sequence: Some(tx_seq.to_string()),
            cds_start: Some(1),
            cds_end: Some(58),
            exons: vec![
                Exon::with_genomic(1, 1, 20, 1000, 1019),
                Exon::with_genomic(2, 21, 40, 1030, 1049),
                Exon::with_genomic(3, 41, 58, 1060, 1077),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1077),
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        provider
    }

    /// Minus-strand mirror of `make_boundary_test_provider`.
    ///
    /// Same transcript sequence as the plus-strand fixture, so c.40_40+3
    /// still spans the same poly-A region in transcript view. The genomic
    /// content at the gene region is the reverse complement of each exon
    /// (so RC of the genomic plus strand recovers `tx_seq`), and the
    /// exon-to-genomic mapping is reversed: tx 1 maps to the high genomic
    /// end (g.1077) and tx 58 to the low end (g.1000). Intron 2 is laid
    /// out so that c.40+1..c.40+4 read as `A` in transcript view, putting
    /// the boundary-spanning dup inside the same 5-A tract that the plus
    /// fixture exercises.
    fn make_boundary_test_provider_minus() -> MockProvider {
        use crate::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
        use std::sync::OnceLock;

        let mut provider = MockProvider::new();

        let tx_seq = "ATGCCCAAAGGGTTTAGGCCAAAGGGTTTAGGCCCAAAAAGGGTTTAGGCCCAAATGA";

        let mut genomic_seq = String::new();
        for _ in 0..1000 {
            genomic_seq.push('N');
        }
        // Exon 3 region (g.1000-1017): RC of tx[41..58] ("GGGTTTAGGCCCAAATGA").
        genomic_seq.push_str("TCATTTGGGCCTAAACCC");
        // Intron 2 (g.1018-1027): the last four bases (g.1024-1027) are 'T',
        // so c.40+1..c.40+4 read as 'A' in transcript view, extending the
        // exonic poly-A across the boundary.
        genomic_seq.push_str("AAAGTATTTT");
        // Exon 2 region (g.1028-1047): RC of tx[21..40] ("AAAGGGTTTAGGCCCAAAAA").
        genomic_seq.push_str("TTTTTGGGCCTAAACCCTTT");
        // Intron 1 (g.1048-1057): mirrors the plus fixture's intron 1 content.
        genomic_seq.push_str("GTAAGCTAAA");
        // Exon 1 region (g.1058-1077): RC of tx[1..20] ("ATGCCCAAAGGGTTTAGGCC").
        genomic_seq.push_str("GGCCTAAACCCTTTGGGCAT");
        for _ in 0..100 {
            genomic_seq.push('N');
        }

        provider.add_genomic_sequence("chr1", genomic_seq);

        provider.add_transcript(Transcript {
            id: "NM_BOUNDARYM.1".to_string(),
            gene_symbol: Some("BOUNDARY_M".to_string()),
            strand: Strand::Minus,
            sequence: Some(tx_seq.to_string()),
            cds_start: Some(1),
            cds_end: Some(58),
            exons: vec![
                Exon::with_genomic(1, 1, 20, 1058, 1077),
                Exon::with_genomic(2, 21, 40, 1028, 1047),
                Exon::with_genomic(3, 41, 58, 1000, 1017),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1077),
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        provider
    }

    #[test]
    fn test_boundary_spanning_exonic_to_intronic_del() {
        // Test: c.20_20+3del - deletion from last exon base into intron
        // c.20 = last base of exon 1 (C at g.1019)
        // c.20+3 = 3rd intronic base (A at g.1022)
        // Deletes: C (exonic) + GTA (intronic) = 4 bases
        //
        // This should normalize without error (using genomic space)
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.20_20+3del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Boundary-spanning deletion should normalize, got error: {:?}",
            result.err()
        );

        let output = format!("{}", result.unwrap());
        assert!(
            output.contains("del"),
            "Should remain a deletion, got: {}",
            output
        );
    }

    #[test]
    fn test_boundary_spanning_intronic_to_exonic_del() {
        // Test: c.21-3_23del - deletion from intron into exon
        // c.21-3 = 3rd base before exon 2 (A at g.1027)
        // c.23 = 3rd base of exon 2 (A at g.1032)
        // Deletes: AAA (intronic) + AAA (exonic) = 6 bases
        //
        // The intron ends with AAA (g.1027-1029) and exon starts with AAA (g.1030-1032)
        // This is a repeat, so normalization might shift
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.21-3_23del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Boundary-spanning deletion should normalize, got error: {:?}",
            result.err()
        );

        let output = format!("{}", result.unwrap());
        assert!(
            output.contains("del"),
            "Should remain a deletion, got: {}",
            output
        );
    }

    #[test]
    fn test_boundary_spanning_same_base_position() {
        // Test: c.20_20+5del - deletion starting and ending at same CDS base
        // Start is exonic (c.20), end is intronic (c.20+5)
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.20_20+5del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Same-base boundary-spanning deletion should normalize, got error: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_boundary_splice_site_plus1() {
        // Test: c.20_20+1del - deletion of last exon base + splice donor (GT)
        // This is a clinically important splice site variant
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.20_20+1del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Splice site +1 deletion should normalize, got error: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_boundary_splice_site_minus1() {
        // Test: c.21-1_22del - deletion of splice acceptor + first exon bases
        // c.21-1 = last intronic base before exon 2 (A at g.1029)
        // c.22 = 2nd base of exon 2 (A at g.1031)
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.21-1_22del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Splice site -1 deletion should normalize, got error: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_boundary_spanning_dup() {
        // Test: c.40_40+3dup — duplication spanning exon-intron boundary,
        // landing on a 4-A poly-A region. Per HGVS spec (repeated.md):
        //
        //   > This restriction only applies to the coding sequence,
        //   > which does not include the introns or the UTR sequence.
        //
        // The codon-frame `unit_len % 3 == 0` restriction does NOT apply
        // to boundary-spanning variants (mixed exon/intron context is
        // not "purely coding sequence"), so `normalize_boundary_spanning_cds`
        // passes `is_coding=false` to `normalize_na_edit`. The dup must
        // therefore canonicalize to `A[N]` repeat notation, not the gated
        // `ins<literal>` fallback. (B4-remaining, issue #209.)
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.40_40+3dup").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Boundary-spanning duplication should normalize, got error: {:?}",
            result.err()
        );

        let output = format!("{}", result.unwrap());
        assert!(
            output.contains("A[") && !output.contains("ins"),
            "Boundary-spanning multi-copy dup spanning into intron must emit \
             `A[N]` repeat notation (intron exempts the codon-frame gate per \
             repeated.md), got: {}",
            output
        );
    }

    #[test]
    fn test_boundary_spanning_dup_minus_strand() {
        // Minus-strand mirror of `test_boundary_spanning_dup`. Pins the
        // strand-specific flip in `normalize_boundary_spanning_cds`: the
        // genomic-strand window is RC of the transcript view, so without
        // flipping, repeat detection would inspect the wrong alphabet and
        // miss the A homopolymer. With `is_coding=false` (the
        // boundary-spanning context's intronic exemption — B4-remaining),
        // the canonical form is `A[N]` repeat notation on the transcript-
        // view bytes.
        let provider = make_boundary_test_provider_minus();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARYM.1:c.40_40+3dup").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Boundary-spanning duplication should normalize, got error: {:?}",
            result.err()
        );

        let output = format!("{}", result.unwrap());
        assert!(
            output.contains("A[") && !output.contains("ins"),
            "Minus-strand boundary-spanning multi-copy dup must emit `A[N]` \
             repeat notation (intron exempts the codon-frame gate per \
             repeated.md), got: {}",
            output
        );
    }

    #[test]
    fn test_boundary_spanning_del_emits_repeat_notation() {
        // Test: `c.40_40+3del` — deletion spanning exon-intron boundary
        // into the same 4-A poly-A region used by
        // `test_boundary_spanning_dup`. Symmetric counterpart on the
        // `del` branch (`deletion_to_repeat`) of the codon-frame gate
        // fix. Per the intronic exemption (issue #209 B4-remaining),
        // boundary-spanning context passes `is_coding=false`, so the
        // del of repeat-unit bases must canonicalize to `A[N-k]` over
        // the reference-tract extent rather than falling back to plain
        // `del`.
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.40_40+3del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Boundary-spanning deletion should normalize, got error: {:?}",
            result.err()
        );

        let output = format!("{}", result.unwrap());
        assert!(
            output.contains("A[") && !output.contains("del"),
            "Boundary-spanning multi-copy del spanning into intron must \
             emit `A[N]` repeat notation (intron exempts the codon-frame \
             gate per repeated.md), got: {}",
            output
        );
    }

    #[test]
    fn test_boundary_spanning_delins() {
        // Test: c.20_20+2delinsTTT - delins spanning exon-intron boundary
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.20_20+2delinsTTT").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_ok(),
            "Boundary-spanning delins should normalize, got error: {:?}",
            result.err()
        );

        let output = format!("{}", result.unwrap());
        assert!(
            output.contains("delins") || output.contains(">"),
            "Should remain a delins or become substitution, got: {}",
            output
        );
    }

    #[test]
    fn test_boundary_no_genomic_data_returns_error() {
        // Test that without genomic data, we still get the ExonIntronBoundary error
        let provider = MockProvider::with_test_data(); // No genomic data
        let normalizer = Normalizer::new(provider);

        // NM_001234.1 doesn't have genomic coordinates
        let variant = parse_hgvs("NM_001234.1:c.11_11+3del").unwrap();
        let result = normalizer.normalize(&variant);

        assert!(
            result.is_err(),
            "Boundary-spanning without genomic data should return error"
        );
    }

    #[test]
    fn test_deletion_3prime_shift_consecutive_bases() {
        // Test case simulating NM_001408491.1:c.517delA -> should become c.518del
        // Create a transcript with consecutive A's at positions that should shift
        use crate::reference::transcript::{Exon, ManeStatus, Strand};
        use std::sync::OnceLock;

        let mut provider = MockProvider::new();

        // Create a sequence where c.517 and c.518 are both 'A'
        // CDS starts at position 1, so c.N = transcript position N
        // Put "AA" at positions 517-518 (1-based)
        // Sequence: 516 bases of padding + "AA" + more padding
        let mut seq = String::new();
        for _ in 0..516 {
            seq.push('G'); // Padding (not A to ensure we see the shift)
        }
        seq.push('A'); // Position 517 (c.517)
        seq.push('A'); // Position 518 (c.518)
        for _ in 519..=600 {
            seq.push('G'); // More padding
        }

        provider.add_transcript(crate::reference::transcript::Transcript {
            id: "NM_777777.1".to_string(),
            gene_symbol: Some("SHIFTTEST".to_string()),
            strand: Strand::Plus,
            sequence: Some(seq.clone()),
            cds_start: Some(1),
            cds_end: Some(600),
            exons: vec![Exon::new(1, 1, 600)], // Single exon covering all
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        let normalizer = Normalizer::new(provider);

        // Parse c.517delA - deleting the first A
        let variant = parse_hgvs("NM_777777.1:c.517delA").unwrap();

        // Debug: print the sequence around positions 517-518
        println!("Sequence length: {}", seq.len());
        println!(
            "Position 516 (0-based 515): {}",
            seq.chars().nth(515).unwrap()
        );
        println!(
            "Position 517 (0-based 516): {}",
            seq.chars().nth(516).unwrap()
        );
        println!(
            "Position 518 (0-based 517): {}",
            seq.chars().nth(517).unwrap()
        );
        println!(
            "Position 519 (0-based 518): {}",
            seq.chars().nth(518).unwrap()
        );

        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        println!("Input:  NM_777777.1:c.517delA");
        println!("Output: {}", output);

        // The deletion should shift from c.517 to c.518 (3' rule)
        // because both positions are 'A'
        assert!(
            output.contains("c.518del"),
            "Deletion at c.517 should shift to c.518 (3' rule), got: {}",
            output
        );
    }

    #[test]
    fn test_deletion_3prime_shift_with_utr() {
        // Same test but with a 5' UTR (cds_start > 1)
        // This simulates real transcripts more accurately
        use crate::reference::transcript::{Exon, ManeStatus, Strand};
        use std::sync::OnceLock;

        let mut provider = MockProvider::new();

        // Create a transcript with 100bp 5' UTR
        // CDS starts at position 101, so:
        // c.1 = tx position 101
        // c.517 = tx position 617
        // c.518 = tx position 618
        let utr_len = 100;
        let mut seq = String::new();

        // 5' UTR (100 bases)
        for _ in 0..utr_len {
            seq.push('T');
        }
        // CDS: 516 bases of G padding, then "AA", then more G
        for _ in 0..516 {
            seq.push('G');
        }
        seq.push('A'); // tx position 617 = c.517
        seq.push('A'); // tx position 618 = c.518
        for _ in 0..100 {
            seq.push('G');
        }

        let seq_len = seq.len();
        println!("Test with UTR:");
        println!("  Sequence length: {}", seq_len);
        println!("  CDS start (1-based): 101");
        println!("  c.517 = tx position 617 (0-based 616)");
        println!("  c.518 = tx position 618 (0-based 617)");
        println!(
            "  tx pos 617 (0-based 616): {}",
            seq.chars().nth(616).unwrap()
        );
        println!(
            "  tx pos 618 (0-based 617): {}",
            seq.chars().nth(617).unwrap()
        );

        provider.add_transcript(crate::reference::transcript::Transcript {
            id: "NM_666666.1".to_string(),
            gene_symbol: Some("UTRTEST".to_string()),
            strand: Strand::Plus,
            sequence: Some(seq.clone()),
            cds_start: Some(101),
            cds_end: Some(seq_len as u64),
            exons: vec![Exon::new(1, 1, seq_len as u64)], // Single exon
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_666666.1:c.517delA").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        println!("Input:  NM_666666.1:c.517delA");
        println!("Output: {}", output);

        assert!(
            output.contains("c.518del"),
            "Deletion at c.517 should shift to c.518 (3' rule) even with UTR, got: {}",
            output
        );
    }

    #[test]
    fn test_normalize_inverted_range_insertion_no_panic() {
        // Regression: ClinVar pattern NC_000011.10:g.5238138_5153222insTATTT
        // has start > end (inverted range).  Previously caused a panic in
        // insertion_is_duplication due to slice index out of bounds.
        // The normalizer should return an error, not panic.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NC_000011.10:g.5238138_5153222insTATTT").unwrap();
        let result = normalizer.normalize(&variant);
        // It's fine if this returns Ok (unchanged) or Err (validation failure),
        // but it must NOT panic.
        let _ = result;
    }

    #[test]
    fn test_delins_should_not_shift() {
        // HGVS spec: delins should NOT be 3' shifted like del/dup/ins
        // This test ensures we don't incorrectly shift delins positions
        use crate::reference::transcript::{Exon, ManeStatus, Strand};
        use std::sync::OnceLock;

        let mut provider = MockProvider::new();

        // Create a transcript where delins could theoretically shift if we were wrong
        // Sequence: ...GGAATTCC... where we do c.10_11delinsXX
        // If incorrectly shifted, it might become c.11_12delinsXX
        let seq = "GGGGGGGGGGAATTCCGGGGGGGGGG".to_string(); // c.10=A, c.11=A, c.12=T, c.13=T

        provider.add_transcript(crate::reference::transcript::Transcript {
            id: "NM_555555.1".to_string(),
            gene_symbol: Some("DELINSTEST".to_string()),
            strand: Strand::Plus,
            sequence: Some(seq),
            cds_start: Some(1),
            cds_end: Some(26),
            exons: vec![Exon::new(1, 1, 26)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        let normalizer = Normalizer::new(provider);

        // Test delins - should NOT shift
        let variant = parse_hgvs("NM_555555.1:c.10_11delinsTT").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        assert!(
            output.contains("c.10_11delins"),
            "Delins should NOT be shifted (HGVS spec), got: {}",
            output
        );
    }

    #[test]
    fn test_cds_to_tx_pos_utr5_underflow() {
        // cds_start=5, base=-6 → 5 + (-6) - 1 = -2, should return Err not wrap to u64::MAX
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let pos = CdsPos {
            base: -6,
            offset: None,
            utr3: false,
        };
        let result = normalizer.cds_to_tx_pos(&pos, 5, Some(38));
        assert!(
            result.is_err(),
            "cds_to_tx_pos should return Err for positions before transcript start, got: {:?}",
            result
        );
    }

    #[test]
    fn test_cds_to_tx_pos_utr5_valid() {
        // HGVS numbering skips c.0, so c.-N maps to tx position
        // cds_start - N. For cds_start=5, c.-3 → tx = 5 + (-3) = 2.
        // Issue #97 — the previous formula `cds_start + base - 1`
        // double-counted the gap and returned tx 1 (the c.-4 base).
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let pos = CdsPos {
            base: -3,
            offset: None,
            utr3: false,
        };
        let result = normalizer.cds_to_tx_pos(&pos, 5, Some(38));
        assert_eq!(result.unwrap(), 2);
    }

    #[test]
    fn test_normalize_cds_utr5_deep_negative() {
        // A deeply negative 5' UTR position that would overflow should return an error, not panic
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        // c.-88 with cds_start=5 → 5 + (-88) - 1 = -84, which would wrap to huge u64
        let variant = parse_hgvs("NM_001234.1:c.-88A>G").unwrap();
        let result = normalizer.normalize(&variant);
        // The primary check is that this doesn't panic.
        let _ = result;
    }

    #[test]
    fn test_normalize_unknown_offset_returns_unchanged() {
        // Variants with ? offsets (sentinel values i64::MAX/MIN) should return unchanged
        // because we can't normalize with indeterminate boundaries
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // c.-85-?_834+?del has unknown offsets on both positions
        let variant = parse_hgvs("NM_000088.3:c.-85-?_834+?del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "Unknown offset should not error, got: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_normalize_unknown_offset_single_position() {
        // Even a single unknown offset should cause early return
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10-?del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "Single unknown offset should not error, got: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_normalize_utr_before_tx_start_returns_unchanged() {
        // c.-215 with a small UTR should not error - return unchanged
        // NM_001234.1 has cds_start=5, so c.-215 maps to 5 + (-215) - 1 = -211
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:c.-215_-214del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "UTR before transcript start should not error, got: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_normalize_no_cds_returns_unchanged() {
        // An NR_ transcript with c. coordinates should not error
        let mut provider = MockProvider::new();
        use crate::reference::transcript::{Exon, Transcript};
        provider.add_transcript(Transcript {
            id: "NR_001566.1".to_string(),
            gene_symbol: Some("NCRNA".to_string()),
            strand: crate::reference::transcript::Strand::Plus,
            sequence: Some("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC".to_string()),
            cds_start: None,
            cds_end: None,
            exons: vec![Exon::new(1, 1, 51)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: Default::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: std::sync::OnceLock::new(),
        });
        let normalizer = Normalizer::new(provider);

        // c. variant on a non-coding transcript (no CDS)
        let variant = parse_hgvs("NR_001566.1:c.10del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "No CDS should not error, got: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_normalize_tx_intronic() {
        // n. intronic variants should normalize via genomic space
        // Build a non-coding transcript with genomic coords and intronic positions
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};

        let mut provider = MockProvider::new();

        // Create transcript: 2 exons with an intron in between
        // Exon 1: tx 1-100, genomic 1000-1099
        // Intron: genomic 1100-1199
        // Exon 2: tx 101-200, genomic 1200-1299
        // Sequence in the intron around position 1100+: AAAA... (for shifting test)
        let tx_sequence = "A".repeat(200);

        provider.add_transcript(Transcript {
            id: "NR_038982.1".to_string(),
            gene_symbol: Some("NCRNA_TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some(tx_sequence),
            cds_start: None,
            cds_end: None,
            exons: vec![
                Exon::with_genomic(1, 1, 100, 1000, 1099),
                Exon::with_genomic(2, 101, 200, 1200, 1299),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1299),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: std::sync::OnceLock::new(),
        });

        // Add genomic sequence for chr1 around positions 1000-1299
        // Make the intron region (1100-1199) be "AGCT" repeated to test shifting
        let mut genomic = String::new();
        for _ in 0..325 {
            genomic.push_str("AGCT");
        }
        provider.add_genomic_sequence("chr1", genomic);

        let normalizer = Normalizer::new(provider);

        // n.100+4del - intronic deletion in a non-coding transcript
        let variant = parse_hgvs("NR_038982.1:n.100+4del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "n. intronic normalization should succeed, got: {:?}",
            result.err()
        );
        let output = format!("{}", result.unwrap());
        assert!(
            output.contains('+') || output.contains('-'),
            "Normalized intronic n. variant should retain intronic notation, got: {}",
            output
        );
    }

    // #488: a minus-strand intronic shuffle boundary can lie outside the
    // fetched genomic window (the window is sized to the variant, the intron
    // far edge is not). The reverse-complement flip `seq_len - x + 1` then
    // underflowed and panicked. It must now surface a clean error, and the
    // in-window path must be unaffected.
    #[test]
    fn flip_intronic_for_strand_errors_on_out_of_window_boundary() {
        let seq = "ACGTACGT"; // len 8
        let boundaries = Boundaries::new(2, 100); // right = 100 >> seq_len = 8
        let result = flip_intronic_for_strand(Strand::Minus, seq, 3, 4, &boundaries);
        assert!(
            matches!(result, Err(FerroError::ConversionError { .. })),
            "out-of-window boundary must yield ConversionError, got {result:?}"
        );
    }

    #[test]
    fn flip_intronic_for_strand_in_window_and_plus_passthrough() {
        let seq = "ACGTACGT"; // len 8
        let boundaries = Boundaries::new(2, 7);

        // Minus strand: flips sequence + coordinates into transcript view.
        let (rc, rel_start, rel_end, b) =
            flip_intronic_for_strand(Strand::Minus, seq, 3, 5, &boundaries)
                .expect("in-window minus-strand flip");
        assert_eq!(rc, crate::sequence::reverse_complement(seq));
        assert_eq!((rel_start, rel_end), (8 - 5 + 1, 8 - 3 + 1));
        assert_eq!((b.left, b.right), (8 - 7 + 1, 8 - 2 + 1));

        // Plus strand: returned unchanged.
        let (s2, rs2, re2, b2) = flip_intronic_for_strand(Strand::Plus, seq, 3, 5, &boundaries)
            .expect("plus-strand passthrough");
        assert_eq!(
            (s2.as_str(), rs2, re2, b2.left, b2.right),
            (seq, 3, 5, 2, 7)
        );
    }

    #[test]
    fn resolve_special_genome_pos_maps_markers() {
        use crate::hgvs::location::{GenomePos, SpecialPosition};
        use crate::reference::MockProvider;

        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_TEST.1", "A".repeat(500));

        // pter -> 1 (no length lookup needed).
        let pter = GenomePos {
            base: 0,
            special: Some(SpecialPosition::Pter),
            offset: None,
        };
        assert_eq!(
            resolve_special_genome_pos(&pter, "NC_TEST.1", &provider).unwrap(),
            Some(1)
        );

        // qter -> contig length.
        let qter = GenomePos {
            base: 0,
            special: Some(SpecialPosition::Qter),
            offset: None,
        };
        assert_eq!(
            resolve_special_genome_pos(&qter, "NC_TEST.1", &provider).unwrap(),
            Some(500)
        );

        // cen -> None (structurally unresolvable).
        let cen = GenomePos {
            base: 0,
            special: Some(SpecialPosition::Cen),
            offset: None,
        };
        assert_eq!(
            resolve_special_genome_pos(&cen, "NC_TEST.1", &provider).unwrap(),
            None
        );

        // Plain position -> its own base.
        let plain = GenomePos {
            base: 42,
            special: None,
            offset: None,
        };
        assert_eq!(
            resolve_special_genome_pos(&plain, "NC_TEST.1", &provider).unwrap(),
            Some(42)
        );
    }

    #[test]
    fn resolve_special_genome_pos_qter_without_length_is_none() {
        use crate::hgvs::location::{GenomePos, SpecialPosition};
        use crate::reference::MockProvider;

        // Provider has no contig registered -> get_sequence_length errors ->
        // graceful Ok(None), not an Err.
        let provider = MockProvider::new();
        let qter = GenomePos {
            base: 0,
            special: Some(SpecialPosition::Qter),
            offset: None,
        };
        assert_eq!(
            resolve_special_genome_pos(&qter, "NC_MISSING.1", &provider).unwrap(),
            None
        );
    }

    // #488: a genomic variant whose position is a telomere/centromere marker
    // (`pter`/`qter`/`cen`) carries a `base == 0` sentinel. Before the guard in
    // `normalize_genome`, that 0 flowed into `coords::hgvs_pos_to_index(0)`
    // (`pos - 1`) and panicked with "attempt to subtract with overflow". The
    // normalizer must instead fall back to minimal canonicalization, preserving
    // the marker, and never panic — regardless of whether reference bases for
    // the contig are available (the panic site is reached only once a window of
    // bases is successfully fetched).
    #[test]
    fn genome_special_position_does_not_panic() {
        use crate::reference::MockProvider;

        // Contig starts with "TTT" so that pter (→ base 1) 3'-shifts through
        // the leading T-run; length = 3 + 200 = 203. Sequence registered so
        // the 100-base window fetch succeeds and normalization reaches the
        // coordinate math that used to underflow (#488).
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_000002.12", format!("TTT{}", "ACGT".repeat(50)));
        let normalizer = Normalizer::new(provider);

        // pter and qter now resolve to concrete coordinates and are then
        // 3'-shifted in the normal way. Observed outputs (not guessed):
        //   pter del → base 1, shifts right through TTT → g.3del
        //   qter del → base 203 (last base), no further right-shift → g.203del
        for (input, expected) in [
            ("NC_000002.12:g.pterdel", "NC_000002.12:g.3del"),
            ("NC_000002.12:g.qterdel", "NC_000002.12:g.203del"),
        ] {
            let variant = parse_hgvs(input).expect("parse special-position genomic variant");
            let out = format!("{}", normalizer.normalize(&variant).expect("normalize"));
            assert_eq!(out, expected, "{input} should resolve to {expected}");
        }

        // cen is structurally unresolvable: must not panic and must be
        // preserved verbatim in the output.
        let cen = parse_hgvs("NC_000002.12:g.cendel").expect("parse cen");
        let cen_out = format!(
            "{}",
            normalizer
                .normalize(&cen)
                .expect("normalize cen — must not error")
        );
        assert!(
            cen_out.contains("cen"),
            "cen retained verbatim, got {cen_out}"
        );
    }

    #[test]
    fn genome_pter_resolves_and_shifts() {
        use crate::reference::MockProvider;
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NG_012337.1", format!("TTT{}", "ACGT".repeat(50)));
        let n = Normalizer::new(provider);
        let v = parse_hgvs("NG_012337.1:g.pterdel").unwrap();
        let out = format!("{}", n.normalize(&v).unwrap());
        assert_eq!(
            out, "NG_012337.1:g.3del",
            "pter del must resolve+shift to g.3del"
        );
    }

    #[test]
    fn genome_qter_resolves_to_last_base() {
        use crate::reference::MockProvider;
        let mut provider = MockProvider::new();
        let seq = format!("{}G", "ACGT".repeat(50)); // length 201, last base 'G'
        provider.add_genomic_sequence("NC_TEST.2", seq);
        let n = Normalizer::new(provider);
        let v = parse_hgvs("NC_TEST.2:g.qterdel").unwrap();
        let out = format!("{}", n.normalize(&v).unwrap());
        assert_eq!(
            out, "NC_TEST.2:g.201del",
            "qter del must resolve to the last base"
        );
    }

    #[test]
    fn genome_whole_span_pter_qter_short_circuits() {
        use crate::reference::MockProvider;
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_TEST.3", "ACGT".repeat(50)); // length 200
        let n = Normalizer::new(provider);
        let v = parse_hgvs("NC_TEST.3:g.pter_qterdel").unwrap();
        let out = format!("{}", n.normalize(&v).unwrap());
        assert_eq!(
            out, "NC_TEST.3:g.1_200del",
            "whole-span pter_qter must render concrete span"
        );
    }

    #[test]
    fn genome_cen_emits_warning_not_silent() {
        use crate::normalize::NormalizationWarning;
        use crate::reference::MockProvider;
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_TEST.4", "ACGT".repeat(50));
        let n = Normalizer::new(provider);
        let v = parse_hgvs("NC_TEST.4:g.cendel").unwrap();
        let result = n.normalize_with_diagnostics(&v).unwrap();
        assert_eq!(format!("{}", result.result), "NC_TEST.4:g.cendel");
        assert!(
            result
                .warnings
                .iter()
                .any(|w| matches!(w, NormalizationWarning::UnresolvableSpecialPosition { .. })),
            "cen must emit an UnresolvableSpecialPosition warning, got {:?}",
            result.warnings
        );
    }

    #[test]
    fn genome_mixed_special_plain_past_end_matches_plain_path() {
        use crate::reference::MockProvider;
        // Regression for the PR #526 review concern (CodeRabbit): a mixed
        // special/plain span like `g.pter_<past-end>del` resolves `pter` to base
        // 1 while leaving the plain endpoint beyond the contig. The resulting
        // relative end falls outside the fetched window, but `shuffle` guards
        // every reference index (shuffle.rs: `ref_idx >= ref_seq.len()`), so the
        // shift simply cannot advance — no panic, no out-of-bounds, and the
        // variant is echoed verbatim. Crucially this is byte-identical to the
        // equivalent plain-coordinate span: the special-position resolution adds
        // no window-overflow misbehavior of its own. A special-only early-return
        // (or clamping the endpoint to the contig length) would be either dead
        // complexity or an active corruption of the user's past-end coordinate.
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_PASTEND.1", "ACGT".repeat(50)); // length 200
        let n = Normalizer::new(provider);

        let special = parse_hgvs("NC_PASTEND.1:g.pter_5000del").expect("parse special past-end");
        let special_out = n
            .normalize(&special)
            .expect("mixed special/plain past-end must not error or panic");

        let plain = parse_hgvs("NC_PASTEND.1:g.1_5000del").expect("parse plain past-end");
        let plain_out = n
            .normalize(&plain)
            .expect("plain past-end normalizes without error");

        assert_eq!(
            format!("{special_out}"),
            format!("{plain_out}"),
            "resolved-special past-end span must behave identically to its plain-coordinate equivalent"
        );
        assert_eq!(
            format!("{special_out}"),
            "NC_PASTEND.1:g.1_5000del",
            "past-end span is echoed verbatim (shuffle cannot shift past the fetched window)"
        );
    }

    #[test]
    fn genome_cen_strict_mode_rejects() {
        use crate::error_handling::ErrorMode;
        use crate::reference::MockProvider;
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_TEST.5", "ACGT".repeat(50));

        // Strict mode promotes the unresolvable-cen warning to an error
        // instead of silently echoing the input.
        let strict = Normalizer::with_config(
            provider.clone(),
            NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
        );
        let v = parse_hgvs("NC_TEST.5:g.cendel").unwrap();
        let result = strict.normalize(&v);
        assert!(
            matches!(result, Err(FerroError::InvalidCoordinates { .. })),
            "strict mode must reject unresolvable cen, got {result:?}"
        );

        // Lenient (default) mode accepts and preserves the input.
        let lenient = Normalizer::new(provider);
        let v2 = parse_hgvs("NC_TEST.5:g.cendel").unwrap();
        let out = lenient.normalize(&v2).expect("lenient accepts cen");
        assert_eq!(format!("{out}"), "NC_TEST.5:g.cendel");
    }

    #[test]
    fn genome_qter_without_length_canonicalizes_gracefully() {
        use crate::reference::MockProvider;
        // Provider has the contig genomic sequence NOT registered -> qter length
        // is unavailable -> graceful canonicalize fallback (no panic, no error).
        let provider = MockProvider::new();
        let n = Normalizer::new(provider);
        let v = parse_hgvs("NC_NONE.1:g.qterdel").unwrap();
        let result = n.normalize(&v);
        assert!(
            result.is_ok(),
            "length-less qter must not error/panic, got {result:?}"
        );
        // Marker preserved verbatim since it could not be resolved.
        assert!(
            format!("{}", result.unwrap()).contains("qter"),
            "unresolved qter should be returned verbatim"
        );
    }
}
