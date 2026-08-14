//! Normalization engine
//!
//! Implements the core HGVS variant normalization algorithm including
//! 3' shifting and boundary detection. (A 5' arm exists as an internal
//! differential oracle and is not reachable from any entry point — see
//! [`config::ShuffleDirection`] and `README.md` rule 6.)
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
// `pub` only under the `dev` feature (test/bench builds), so the
// `seqfirst_align` criterion benchmark — an external crate — can reach
// `AlignmentDag::build`/`edit_distance`; `pub(crate)` otherwise, since there is
// no public API here yet.
//
// Not "unwired into `normalize_allele`": it IS compiled in and run there, as a
// shadow comparison behind `FERRO_SEQFIRST_SHADOW=1` that never affects output.
// What remains is promoting it from shadow to authoritative. (Same correction as
// the module's own header — this comment is its sibling and said the same wrong
// thing.)
pub mod from_sequences;
#[cfg(feature = "dev")]
pub mod seqfirst;
#[cfg(not(feature = "dev"))]
pub(crate) mod seqfirst;
pub mod sequence_pair;
pub mod shuffle;
pub mod validate;

// The three block partitioners, re-exported for the `dump_partitions` example.
// Dev-only measurement surface, not API — see the module's own doc.
#[cfg(feature = "dev")]
pub use merge::dev_partitioners;

// How often a sequence-first partitioner declined and `partition_block` answered
// under its name. Not gated on `dev`: a bake-off is run from whatever build the
// measurement uses, and a census that exists only in some builds is one a run can
// forget to read. See `PartitionDeclineCounts`.
pub use merge::{partition_decline_counts, PartitionDeclineCounts};

// The denominator beside that census: blocks cut on EVERY arm, `Live` included.
//
// This one IS a census that exists only in some builds, which is exactly what the
// comment above declines to do — so the trade is stated rather than left to be
// read off the `cfg`. The reason it goes the other way is that the shipped `Live`
// path must not grow a per-block atomic to serve a measurement, and `Live` is the
// arm this counter exists to cover. Its only consumer,
// `examples/measure_spec_conformance_per_arm.rs`, refuses to run in a build where
// it is absent rather than reporting a zero it cannot distinguish from a corpus
// that never reached the partitioner — see `partition_blocks_cut`'s own doc.
#[cfg(debug_assertions)]
pub use merge::partition_blocks_cut;

// The refusal a binary should print and exit on when `FERRO_PARTITION` named an
// arm this build has not got. A library cannot refuse from
// `canonicalize_from_sequence` -- it is infallible -- so the refusal is offered
// to entry points that can return a failure.
pub use merge::partition_switch_startup_error;

use crate::coords::{hgvs_pos_to_index, index_to_hgvs_pos};
use crate::error::FerroError;
use crate::hgvs::edit::{
    Base, InsertedSequence, NaEdit, ProteinEdit, ProteinInsSeq, RepeatCount, Sequence,
};
use crate::hgvs::interval::{CdsInterval, Interval, ProtInterval};
use crate::hgvs::location::{
    AminoAcid, CdsPos, GenomePos, ProtPos, RnaPos, SpecialPosition, TxPos,
};
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{
    Accession, AllelePhase, AlleleVariant, CdsVariant, GenomeVariant, HgvsVariant, LocEdit,
    MtVariant, RnaVariant, TxVariant,
};
use crate::hgvs::HgvsVariant as HV;
use crate::reference::transcript::Strand;
use crate::reference::ReferenceProvider;
use boundary::Boundaries;
pub use config::NormalizeConfig;
#[doc(hidden)]
pub use config::ShuffleDirection;
#[cfg(debug_assertions)]
use merge::OutOfBoundsCoordinate;
use rules::{
    canonicalize_conversion_to_delins, canonicalize_edit, canonicalize_insertion_expand,
    is_real_substitution, is_uncertain_real_substitution, needs_normalization,
    rna_uracil_to_thymine, should_canonicalize, DelinsSubedit, InsCoordKind,
};
use shuffle::shuffle;

/// Check if a CDS position has an unknown (?) offset sentinel value
fn has_unknown_offset_cds(pos: &CdsPos) -> bool {
    pos.has_unknown_offset()
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
/// this helper sees `CdsPos { base: 0, utr3: false }`. Special positions
/// (pter/qter/cen, which also carry `base == 0`) are similarly skipped by
/// the caller before reaching this helper. The arm is dead code in the
/// W4004 path; kept here for parity with `Normalizer::cds_to_tx_pos` so
/// this helper remains a drop-in replacement if a future caller needs it
/// outside the bounds check.
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
    // Unknown positions and special markers (pter/qter/cen) can't be bounds-checked.
    if pos.is_unknown() || pos.is_special() {
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

/// Canonicalize a plain past-CDS coordinate to its 3'UTR `*N` form for output.
///
/// A plain `c.<N>` with `N > cds_len` is out-of-scheme HGVS: coding-DNA
/// numbering ends at the last base of the stop codon, and 3'UTR positions use
/// the `*` prefix (`background/numbering.md` L21/L30) — so the nucleotide `N -
/// cds_len` bases past the stop is *named* `c.*(N - cds_len)`, not `c.<N>`.
/// When `N` maps into the 3'UTR (`1 <= N - cds_len <= utr3_len`) this returns
/// that canonical `*N` position.
///
/// This is a rendering canonicalization only — the same nucleotide, in valid
/// syntax — so callers still emit the `PositionPastEnd` (W4004) warning to flag
/// that the input used out-of-scheme numbering, and strict mode still rejects.
/// A position past the whole transcript (`N - cds_len > utr3_len`) has no valid
/// `*N` form and is returned unchanged, as are 3'UTR (`*N`), 5'UTR (`-N`),
/// intronic-offset, special, and unknown positions. See #920/#336.
fn canonicalize_pastcds_pos_to_utr3(
    pos: &CdsPos,
    transcript: &crate::reference::transcript::Transcript,
) -> CdsPos {
    if pos.is_unknown() || pos.utr3 || pos.special.is_some() || pos.offset.is_some() || pos.base < 1
    {
        return *pos;
    }
    let (Some(cds_len), Some(utr3_len)) = (transcript.cds_length(), transcript.utr3_length())
    else {
        return *pos;
    };
    let base = pos.base as u64;
    if base > cds_len {
        let utr3_base = base - cds_len;
        if utr3_base >= 1 && utr3_base <= utr3_len {
            return CdsPos::utr3(utr3_base as i64);
        }
    }
    *pos
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

/// Check if a TxPos has an unknown (?) offset sentinel value.
///
/// Delegates to [`TxPos::has_unknown_offset`] so the predicate has one
/// definition. It used to be spelled out here as well, and a rule written in
/// several places is how they drift apart.
fn has_unknown_offset_tx(pos: &TxPos) -> bool {
    pos.has_unknown_offset()
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

/// Everything the codon-frame gate needs for one [`Normalizer::normalize_na_edit`]
/// call: whether the axis has a reading frame at all, and if so where the CDS is.
///
/// # The rule
///
/// `DNA/repeated.md` L21-24 (and `RNA/repeated.md` L24-27 in the same terms):
///
/// > **exception**: using a coding DNA reference sequence ("c." description), a
/// > repeated sequence variant description can be used only for repeat units with
/// > a length which is a multiple of 3 […] This restriction only applies to the
/// > coding sequence, which does not include the introns or the UTR sequence.
///
/// So the gate needs two things: does this axis have a reading frame, and does the
/// span under consideration lie inside the CDS proper.
///
/// # Why one type instead of two arguments
///
/// This replaces the pair `is_coding: bool` + `cds_span: Option<(u64, u64)>`
/// (#1206). They encoded the same question at different scopes and had to agree
/// per axis, which the type system did not enforce — the invalid third
/// combination, "a coding verdict with no CDS bounds", was representable and
/// silently disabled the gate for any site asking about a span other than the
/// input's. Passing `None` on `r.` compiled, passed clippy, and reverted #1192
/// (verified during #1185's rebase: `r.4_9a[8]`, the spec-forbidden form, instead
/// of the gated answer). The parameter's own doc had already gone stale claiming
/// `r.` should pass `None`.
///
/// Now an axis either has a frame — and therefore has bounds — or it does not.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CodonGate {
    /// The axis has no reading frame, so the multiple-of-3 rule never applies:
    /// `g.`, `m.`/`o.` (genomic-style), `n.`, and every intronic / boundary-
    /// spanning context that the spec's own carve-out exempts.
    NotApplicable,
    /// A coding axis — `c.`, and `r.` on a coding transcript, which IS the `c.`
    /// axis (`background/numbering.md` L58/L61, #469).
    Coding {
        /// Precomputed verdict for the span the INPUT edit occupies: `true` when
        /// it lies wholly inside the CDS proper.
        input_span_is_coding: bool,
        /// CDS bounds, 1-based inclusive, in the same (transcript) frame as the
        /// positions handed to `normalize_na_edit` — so a site deciding about a
        /// span *other* than the input's can re-answer the question for that span
        /// (#1185) instead of reusing a verdict computed for the wrong one.
        cds: (u64, u64),
    },
}

impl CodonGate {
    /// The gate for an edit occupying `input_start..=input_end` (1-based,
    /// inclusive, transcript frame) on an axis whose CDS is `cds`.
    ///
    /// `cds == None` means the record carries no CDS end, so we cannot verify that
    /// anything lies inside the CDS proper; the pre-#1206 code conservatively
    /// treated that as UTR-touching and skipped the gate, which is exactly
    /// [`CodonGate::NotApplicable`].
    ///
    /// `normalize_cds` and `normalize_rna` computed this verdict with two separate
    /// copies of the same expression, which is how #1192's divergence became
    /// possible in the first place; there is now one.
    fn for_input_span(cds: Option<(u64, u64)>, input_start: u64, input_end: u64) -> Self {
        match cds {
            Some((cds_start, cds_end)) => CodonGate::Coding {
                input_span_is_coding: input_start >= cds_start && input_end <= cds_end,
                cds: (cds_start, cds_end),
            },
            None => CodonGate::NotApplicable,
        }
    }

    /// Does the gate apply to the span the input edit occupies?
    ///
    /// This is the verdict the `rules` helpers take as their `is_coding` argument.
    fn input_span_is_coding(self) -> bool {
        match self {
            CodonGate::NotApplicable => false,
            CodonGate::Coding {
                input_span_is_coding,
                ..
            } => input_span_is_coding,
        }
    }

    /// The CDS bounds, or `None` on an axis with no reading frame.
    ///
    /// Prefer [`Self::span_is_coding`] — it answers the question rather than
    /// handing out the raw bounds to re-pair with a separate verdict, which is the
    /// shape this type exists to eliminate.
    ///
    /// The one caller is `rules::insertion_to_repeat` (#1210/#1212), which cannot
    /// use `span_is_coding`: it *discovers* the tandem tract internally, so there
    /// is no span to ask about until it is already inside the helper. Its
    /// signature still takes the verdict and the bounds separately, so the pairing
    /// is reconstructed here — but from a gate that is internally consistent by
    /// construction, so the invalid combination cannot be produced. Pushing
    /// `CodonGate` down into `rules` would close that last seam; it is left out of
    /// this refactor because several of that helper\'s unit tests deliberately pass
    /// a coding verdict with no bounds, so the conversion is not mechanical and
    /// would change what they test.
    fn cds_bounds(self) -> Option<(u64, u64)> {
        match self {
            CodonGate::NotApplicable => None,
            CodonGate::Coding { cds, .. } => Some(cds),
        }
    }

    /// Does the gate apply to `start..=end` (1-based, inclusive, transcript
    /// frame)?
    ///
    /// Use this — not [`Self::input_span_is_coding`] — at any site deciding about
    /// a span the input edit does not occupy, such as a tract a shuffle walked
    /// into. A 3'-shifted tract that runs out of the CDS into the 3'UTR is exempt
    /// by the spec's carve-out even though the input span was CDS-resident
    /// (#1185).
    fn span_is_coding(self, start: u64, end: u64) -> bool {
        match self {
            CodonGate::NotApplicable => false,
            CodonGate::Coding {
                cds: (cds_start, cds_end),
                ..
            } => start >= cds_start && end <= cds_end,
        }
    }
}

/// Which side of `anchor` a shuffled insertion payload came to rest on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum BoundarySide {
    /// Payload rests immediately 5' of `anchor` (`…ins` before `ref[anchor]`).
    FivePrime,
    /// Payload rests immediately 3' of `anchor` (`…ins` after `ref[anchor]`).
    ThreePrime,
}

/// The coordinate identity shared by every boundary clamp in this module:
/// "insert `A'` immediately outside `anchor`" is *identically* "delete
/// `ref[anchor]`, insert the concatenation".
///
/// - [`BoundarySide::FivePrime`] → `delins(A' ++ ref[anchor])`
/// - [`BoundarySide::ThreePrime`] → `delins(ref[anchor] ++ A')`
///
/// It merely moves the delete boundary one base, so it is exact for **any**
/// `A'` and needs no equivalence check. That argument is what PR #1170 (5', at
/// `cds_start`), issue #387 (3', at `cds_end`) and issue #1202 (both, at the
/// transcript bounds) each rely on — this function is the single place it is
/// implemented, so the proof lives once rather than once per boundary.
///
/// Only the *rewrite* is shared. Each caller keeps its own **detection**,
/// because the two boundary kinds differ in kind, not degree: `cds_start` /
/// `cds_end` are *axis-change* boundaries whose far side (`c.-1`, `c.*1`) is
/// perfectly representable, so clamping there is a policy with carve-outs
/// (#401 spanning dups, #387's edit-type allow-list, #383's Delins restore);
/// the transcript bounds are *representability* boundaries with no far side at
/// all, so that clamp is unconditional. Folding the detections together would
/// mean threading those carve-outs through a single predicate for no gain.
///
/// Returns `None` — leaving the caller's values untouched — when `anchor` lies
/// outside `seq` or the reference byte there is not a base we can spell.
pub(crate) fn insertion_to_boundary_delins(
    seq: &[u8],
    a_prime: &Sequence,
    anchor: u64,
    side: BoundarySide,
) -> Option<NaEdit> {
    let bases = boundary_delins_bases(seq, a_prime.bases(), anchor, side)?;
    Some(NaEdit::Delins {
        sequence: InsertedSequence::Literal(Sequence::new(bases)),
        deleted: None,
        deleted_length: None,
        substitution_reference: None,
    })
}

/// The replacement bases the boundary identity produces: `payload` concatenated
/// with `ref[anchor]`, on whichever side the payload rests.
///
/// The whole of [`insertion_to_boundary_delins`]'s arithmetic, split out because
/// the sequence-first derivation in [`merge`] needs the same rewrite in a
/// different currency: it types its derived pieces straight into `Vec<Base>` and
/// only later hands them to a builder that decides the `NaEdit`, so routing
/// through the `NaEdit`-returning form above would mean constructing a `Delins`
/// solely to take it apart again. Both callers therefore share one
/// implementation of the identity rather than two copies that can drift.
///
/// Returns `None` when `anchor` lies outside `seq` or the reference byte there
/// is not a base we can spell — the same two refusals the caller above documents.
pub(crate) fn boundary_delins_bases(
    seq: &[u8],
    payload: &[Base],
    anchor: u64,
    side: BoundarySide,
) -> Option<Vec<Base>> {
    let anchor_0b = (anchor as usize).saturating_sub(1);
    let ref_base = Base::from_char(*seq.get(anchor_0b)? as char)?;
    let mut bases = Vec::with_capacity(payload.len() + 1);
    match side {
        BoundarySide::FivePrime => {
            bases.extend_from_slice(payload);
            bases.push(ref_base);
        }
        BoundarySide::ThreePrime => {
            bases.push(ref_base);
            bases.extend_from_slice(payload);
        }
    }
    Some(bases)
}

/// Which ends of the `seq` slice handed to [`clamp_insertion_at_sequence_bounds`]
/// are true boundaries of the underlying entity.
///
/// The transcript-axis callers pass the whole transcript, so both ends are real
/// ([`SequenceEnds::WHOLE`]). `normalize_genome` passes a **window** into the
/// contig, where an end is real only when the window happens to be flush with the
/// contig there. Clamping at a window edge that is merely where the fetch stopped
/// would assert a boundary that does not exist and rewrite a perfectly valid
/// interior insertion into a `delins` at the wrong place.
/// The sequence-first canonicalization in [`merge`] fetches a window too, and
/// applies the same clamp to the pieces it derives, so this is `pub(crate)`
/// rather than private to this module.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SequenceEnds {
    /// `seq[0]` is the entity's first base.
    pub(crate) five_prime: bool,
    /// `seq[len-1]` is the entity's last base.
    pub(crate) three_prime: bool,
}

impl SequenceEnds {
    /// The slice is the entire entity, so both of its ends are real boundaries.
    const WHOLE: Self = Self {
        five_prime: true,
        three_prime: true,
    };

    /// Neither end of the slice is known to be a real boundary, so no clamp may
    /// fire. The conservative answer whenever the entity's length is unknown or
    /// the window plainly sits inside it.
    pub(crate) const INTERIOR: Self = Self {
        five_prime: false,
        three_prime: false,
    };
}

/// Rewrite in place a post-shuffle insertion that has saturated at a **sequence**
/// boundary, into the equivalent `delins` — the only spelling of that haplotype
/// the `n.` / `r.` axes, the `c.` axis outside the CDS, and the `g.` axis can
/// express. Leaves `edit`/`start`/`end` untouched for anything that is not a
/// boundary-saturated literal insertion.
///
/// Named for the *sequence* rather than the transcript since #1205: nothing in it
/// is transcript-specific — it enforces `[1, seq.len()]`, which is equally the
/// contig bounds — and `normalize_genome` is now a caller.
///
/// `n.` and `r.` number the transcript directly: there is no position `0` and
/// no `*`-suffixed overflow axis, so an insertion driven to rest immediately
/// 5' of base 1 or immediately 3' of the last base has no valid pair of
/// adjacent anchors. Left alone it escapes in one of two broken shapes (#1202):
///
/// - 5': the 0-based shuffle result of `0` clamps back through
///   `saturating_sub(1)` to the degenerate `start == end == 1`, emitted as
///   `n.1ins<A'>` — a single-position insertion, which ferro's own parser
///   rejects (`DNA/insertion.md:95-101` requires two adjacent anchors).
/// - 3': the shuffle rests against `boundaries.right`, emitted as
///   `n.<len>_<len+1>ins<A'>`, whose second anchor is past the transcript end
///   and so raises W4004 `PositionPastEnd` when re-normalized in strict mode.
///
/// Both are repaired by applying [`insertion_to_boundary_delins`] at the
/// transcript bounds — the same identity PR #1170 and issue #387 apply at
/// `cds_start` / `cds_end` — giving `n.1delins<A' ++ ref[1]>` at the 5' bound
/// and `n.<len>delins<ref[len] ++ A'>` at the 3' bound.
///
/// This is deliberately keyed on the *post-shuffle* edit rather than the input,
/// so every spelling that shuffles to the same resting payload (a directly
/// written `ins`, a `dup` routed through the shuffle, any start position)
/// collapses to the same minimal, idempotent output.
///
/// Note this is **not** a non-coding-transcript special case: the `n.` axis
/// numbers the whole transcript with no CDS sub-axis regardless of whether the
/// record carries `cds_start`/`cds_end` (#334), so an `n.` variant on a coding
/// transcript reaches this the same way.
///
/// **Relationship to the `c.`-axis clamps.** `normalize_cds` clamps at
/// `cds_start` (#1170) and `cds_end` (#387), gated on `AxisRegion::Cds`. Those
/// are a different boundary, not a superset: a `c.` insertion in a **UTR**
/// region (`AxisRegion::FiveUtr` / `ThreeUtr`) is outside that gate and can
/// still saturate against the transcript bounds. So `normalize_cds` calls this
/// helper too, after its own clamps. All three callers therefore hold the same
/// invariant — a post-shuffle insertion never rests outside `[1, len]` — while
/// each axis keeps its own axis-specific boundary behavior:
///
/// | axis | boundary the axis clamps at | clamped by |
/// |---|---|---|
/// | `c.` (CDS region) | `cds_start` / `cds_end` | `normalize_cds` |
/// | `c.` (UTR regions) | transcript `1` / `len` | this helper |
/// | `n.`, `r.` | transcript `1` / `len` | this helper |
/// | `g.` | contig `1` / `len` | `normalize_genome` (#1205) |
/// | `m.` | contig `1` / `len` | `normalize_mt` (#1217) |
/// | `o.` | n/a — returned unchanged, no shuffle runs | n/a |
///
/// `normalize_genome` became the fourth caller in #1205, which closes the last
/// non-circular gap: before it, a 5'-saturated insertion at a contig start
/// escaped as `g.1ins<A\'>` — the same unparseable shape, on the axis the spec
/// uses for its own counter-example (`DNA/insertion.md:95-101` rejects
/// `g.123insG`) — and a 3'-saturated one as `g.<len>_<len+1>ins<A\'>`, whose
/// second anchor is past the contig end.
///
/// `normalize_mt` became the fifth in #1217. #1205 had left it out on the
/// grounds that a circular reference wraps, so position 1 HAS a valid 5\'
/// neighbour (the last base) and the correct output there is a wraparound
/// description rather than this `delins` rewrite. That does not hold: #129
/// established, and `issue_129_mt_circular_wraparound.rs` pins, that ferro
/// **rejects** `m.<high>_<low>ins`. The spec's reversed-range exception
/// (`deletion.md:17`, SVD-WG006) is granted to `del`/`dup` — `delins` inherits
/// it by composition — while `insertion.md` is silent, so the general 5\'→3\'
/// rule applies, which is also how mutalyzer and biocommons read it. The
/// wraparound `ins` is therefore not an available answer and the
/// single-position `delins`, which needs no reversed range at all, is. (#129
/// separately disabled 3\'-rule shifting across the origin, so coming to rest
/// at the terminus is already the intended behavior; only its rendering was
/// wrong.)
///
/// `o.` stays out: `normalize_core` returns circular `o.` variants unchanged, so
/// no shuffle runs and nothing can saturate. A genuine circular normalizer is
/// #466's circular candidate.
///
/// The genomic caller passes a **window** rather than the whole contig, so it
/// must say which of its ends are real — see [`SequenceEnds`].
fn clamp_insertion_at_sequence_bounds(
    seq: &[u8],
    edit: &mut NaEdit,
    start: &mut u64,
    end: &mut u64,
    ends: SequenceEnds,
) {
    // Cheapest discriminating guard first: only two resting places qualify, and
    // both are two register compares against values already in hand. `seq_len ==
    // 0` needs no separate guard — `insertion_to_boundary_delins` bails on an
    // anchor outside `seq`, which an empty `seq` always is.
    let seq_len = seq.len() as u64;
    let (anchor, side) = if ends.five_prime && *start == 1 && *end == 1 {
        (1, BoundarySide::FivePrime)
    } else if ends.three_prime && *start == seq_len && *end == seq_len + 1 {
        (seq_len, BoundarySide::ThreePrime)
    } else {
        return;
    };
    let NaEdit::Insertion {
        sequence: InsertedSequence::Literal(a_prime),
    } = &*edit
    else {
        return;
    };
    let Some(delins) = insertion_to_boundary_delins(seq, a_prime, anchor, side) else {
        return;
    };
    *edit = delins;
    *start = anchor;
    *end = anchor;
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

    /// Normalization returned FEWER cis members than the input described: two
    /// or more separately-reported variants were coalesced into one.
    /// Code: `MEMBERS_COALESCED_FROM_REPORTED_FORM` (W5005).
    ///
    /// Purely informational — the normalized string is identical whether or not
    /// this warning is produced. It carries the *provenance* the canonical form
    /// deliberately does not: see `ErrorType::MembersCoalescedFromReportedForm`
    /// for why that separation is the point rather than a compromise.
    MembersCoalesced {
        /// Accession the coalesced members sit on.
        accession: String,
        /// Coordinate system: "g" | "c" | "n" | "r" | "m"
        coordinate_system: String,
        /// How many cis members the INPUT described.
        reported_members: usize,
        /// How many the normalized form describes.
        normalized_members: usize,
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

    /// A telomere marker on a genomic-reference `c.` description denotes a
    /// transcript-flank position not numberable in `c.` (#488 Phase 2b).
    /// Code: `TRANSCRIPT_FLANK_NOT_DESCRIBABLE`.
    TranscriptFlankNotDescribable {
        /// Accession of the reference (the NG/LRG/NC-parent c. accession).
        accession: String,
        /// The marker, "pter" or "qter".
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

    /// An intronic offset (`c.<N>±<M>` / `n.<N>±<M>`) appears on a bare
    /// transcript reference (`NM_` c. / `NR_`/`XR_` n. with
    /// `genomic_context: None`) — a spec-invalid description form. Code:
    /// `INTRONIC_ON_BARE_TRANSCRIPT` (W4007). Strict mode (or the errors-axis
    /// override) escalates this to `FerroError::IntronicVariant`; lenient
    /// surfaces it and returns the existing value unchanged. See #486.
    IntronicOnBareTranscript {
        /// Full variant Display, e.g. `"NM_003002.2:c.274+20C>T"`.
        variant: String,
        /// Coordinate system: `"c"` (NM_ coding) or `"n"` (NR_/XR_ non-coding).
        coordinate_system: String,
    },

    /// A `c.`/`p.`/`r.` variant is described against a transcript whose 5'
    /// CDS is annotated incomplete (`cds_start_NF`): no confirmed `ATG`
    /// initiation codon means `c.1`/`p.1` are undefined relative to it, so
    /// it is not an HGVS-recommended reference for a coding-axis
    /// description. The coordinate is preserved verbatim (no re-numbering).
    /// Input-side counterpart to Task 4's projection-side decline. Code:
    /// `INCOMPLETE_CDS_START_REFERENCE` (W5004). See #972.
    IncompleteCdsStartReference {
        /// Transcript accession (e.g. `ENST00000011700.10`).
        accession: String,
        /// Coordinate system of the input variant: `"c"`, `"p"`, or `"r"`.
        coordinate_system: String,
    },

    /// A genome-requiring normalization step could not run because the
    /// configured reference provider carries no genomic sequence data
    /// (`ReferenceProvider::has_genomic_data() == false`, e.g. a
    /// transcripts-only provider). The best-effort result is returned
    /// UNCHANGED from the point the step would have refined it, and this
    /// warning marks the output as degraded so a reduced-capability result
    /// is never mistaken for a fully-normalized one. This is an
    /// environmental limitation, not a defect in the variant: it means a
    /// genome-requiring step was reached but could not be run. The same input
    /// against a genome-backed reference *may* normalize further — the
    /// intronic / boundary-spanning paths genuinely could; the exon-junction
    /// 3'-shuffle landing only *might*, since whether the pattern extends into
    /// the intron is exactly what could not be checked without a genome.
    /// Code: `REDUCED_CAPABILITY_NO_GENOME`.
    ///
    /// Emitted by both the formerly-silent exon/intron junction 3'-shuffle
    /// enhancement (`#670`/`#704`) and the formerly-erroring intronic /
    /// boundary-spanning genomic normalization paths — item 2 of the #1012
    /// warn-and-degrade unification.
    ///
    /// Mode interaction: lenient/silent return the best-effort variant with
    /// this advisory; strict mode promotes it to
    /// `FerroError::ReducedReferenceCapability` (via
    /// `NormalizeConfig::should_reject_reduced_capability`) so a strict caller
    /// never receives a knowingly-degraded result.
    ReducedCapabilityNoGenome {
        /// Full variant Display, e.g. `"NM_INT.1:c.10+5del"`.
        variant: String,
        /// Short description of the genome-requiring step that was skipped,
        /// e.g. `"intronic normalization"` or
        /// `"exon/intron junction 3'-shuffle"`.
        capability: String,
    },
}

/// True iff any concrete position reachable from `boundary` is intronic —
/// covering both a `Single` position and either endpoint of a `Range`
/// boundary (uncertain breakpoints like `c.(100+1_101-1)_(200+1_201-1)del`).
/// `Unknown` (`?`) and otherwise-absent inner positions contribute `false`.
fn boundary_has_intronic<T>(
    boundary: &crate::hgvs::interval::UncertainBoundary<T>,
    is_intronic: impl Fn(&T) -> bool,
) -> bool {
    use crate::hgvs::interval::UncertainBoundary;
    let mu_intronic = |mu: &crate::hgvs::uncertainty::Mu<T>| mu.inner().is_some_and(&is_intronic);
    match boundary {
        UncertainBoundary::Single(mu) => mu_intronic(mu),
        UncertainBoundary::Range { start, end } => mu_intronic(start) || mu_intronic(end),
    }
}

/// EINTRONIC (#486): detect an intronic offset on a **bare transcript
/// reference** (no `NG_(…)`/`NC_(…)` genomic context). Returns the warning
/// to emit, or `None`.
///
/// In scope: a bare coding transcript (`NM_`/`XM_` or LRG `LRG_<N>t<k>`) used
/// with `c.` and a bare non-coding transcript (`NR_`/`XR_` via
/// `is_noncoding_rna()`, or LRG) used with `n.`, with
/// `Accession.genomic_context == None`. The spec's "a (non-)coding DNA
/// reference sequence does not contain introns" rule applies equally to curated
/// (`NM_`/`NR_`), predicted (`XM_`/`XR_`), and LRG transcript references — an
/// LRG transcript is itself a bare reference with no `NG_`/`NC_` genomic
/// context — so all are covered on each axis (#834). Out of scope:
/// genomic-context forms (`genomic_context: Some`), `NG_`/`NC_` references
/// (which never reach the c./n. transcript path), Ensembl `ENST`, and the r.
/// axis. Both `Single`
/// and `Range` (uncertain-breakpoint) position boundaries are inspected; an
/// unknown (`?`) offset still counts as intronic (`CdsPos::is_intronic` treats
/// the unknown-offset sentinel as intronic), which is correct — it is an
/// intronic position whose exact offset is unspecified.
fn intronic_on_bare_transcript_warning(variant: &HgvsVariant) -> Option<NormalizationWarning> {
    intronic_on_bare_transcript_axis(variant).map(|coordinate_system| {
        NormalizationWarning::IntronicOnBareTranscript {
            variant: format!("{variant}"),
            coordinate_system: coordinate_system.to_string(),
        }
    })
}

/// The predicate half of [`intronic_on_bare_transcript_warning`]: the coordinate
/// system (`"c"` or `"n"`) whose bare transcript reference `variant` names an
/// intronic position on, or `None`.
///
/// **One rule, two callers, and the split is deliberate.** The scope prose above
/// documents this function as much as the warning constructor — the clause has
/// exactly one reading in this crate, and that is load-bearing in both
/// directions: the same question decides whether strict mode *refuses an input*
/// and whether `#1704` must *re-parent an output*. Two readings of one clause is
/// how ferro came to refuse a description in strict mode while manufacturing the
/// identical description in lenient. So the warning is derived from this
/// predicate rather than the two being written out separately, which is what
/// makes them unable to drift.
///
/// The reason it is a separate function rather than
/// `…_warning(v).is_some()` is cost, not taste. Every `normalize()` asks this
/// question at least twice — once on the strict ladder and once at
/// [`Normalizer::reparent_junction_exit`] — and the `is_some()` form built a
/// `NormalizationWarning` carrying `format!("{variant}")`, rendering the whole
/// description to a `String` only to drop it. The authored-bare-intronic class
/// is common enough for that to be a real per-call allocation on a path that
/// wants a `bool`.
fn intronic_on_bare_transcript_axis(variant: &HgvsVariant) -> Option<&'static str> {
    match variant {
        HV::Cds(v) => {
            // A bare coding-DNA reference does not contain introns, so an
            // intronic offset on it is a spec-invalid form regardless of which
            // transcript namespace addresses it. Curated/predicted RefSeq
            // (`NM_`/`XM_`) and LRG transcript references (`LRG_<N>t<k>`, which
            // carry no `NG_`/`NC_` genomic context — the LRG *is* the reference)
            // are all bare coding transcripts; treat them uniformly (#834).
            let is_bare_coding_transcript =
                matches!(&*v.accession.prefix, "NM" | "XM") || v.accession.is_lrg();
            if !is_bare_coding_transcript || v.accession.genomic_context.is_some() {
                return None;
            }
            let intronic = boundary_has_intronic(&v.loc_edit.location.start, CdsPos::is_intronic)
                || boundary_has_intronic(&v.loc_edit.location.end, CdsPos::is_intronic);
            intronic.then_some("c")
        }
        HV::Tx(v) => {
            // Same rule on the non-coding axis: a bare `NR_`/`XR_` or LRG
            // non-coding transcript reference used with `n.` has no introns
            // (#834 extends the LRG coverage to match the `c.` arm).
            let is_bare_noncoding_transcript =
                v.accession.is_noncoding_rna() || v.accession.is_lrg();
            if !is_bare_noncoding_transcript || v.accession.genomic_context.is_some() {
                return None;
            }
            let intronic = boundary_has_intronic(&v.loc_edit.location.start, TxPos::is_intronic)
                || boundary_has_intronic(&v.loc_edit.location.end, TxPos::is_intronic);
            intronic.then_some("n")
        }
        _ => None,
    }
}

/// One intronic offset ferro **manufactured**, recorded at the `#670` junction
/// gate that made it (`#1723`).
///
/// # Why this is carried rather than re-derived
///
/// `#1704` asked the provenance question — "did ferro manufacture this intronic
/// offset, or did the author write it?" — at the top of the pipeline, by
/// comparing the whole input description against the whole output. That cannot
/// be made per-leaf, and per-leaf is what the question is: normalization
/// reorders, merges and splits members, so no identity map from an output leaf
/// back to an input leaf survives to the top. The consequence was measured in
/// two directions, both pinned in `defect_371_transcript_exit`: one authored
/// intronic member vetoed the wrapper for a sibling ferro had moved itself, and
/// a second bare accession needing a wrapper was never looked at.
///
/// # The fact is free at the gate
///
/// `normalize_cds` diverts **every** input with an intronic endpoint to
/// `normalize_intronic_cds` / `normalize_boundary_spanning_cds` before the gate
/// is reached (`start_pos.is_intronic() || end_pos.is_intronic()`), and the
/// mirrors in `normalize_tx` / `normalize_rna` do the same. So reaching the gate
/// *is* the statement that both endpoints were exonic, and the gate's own
/// success condition is `crossed_into_intron`. Nothing is compared: manufacture
/// is what the branch means.
///
/// The contig is free there too. The crossing is computed by fetching genomic
/// bases from `boundary_transcript.chromosome`, so a resolved genomic reference
/// is already in hand — `#1704` re-derived it at the top with a second provider
/// lookup that could disagree with the one the crossing actually used.
///
/// # What a record IS evidence of, and what it deliberately is not
///
/// **A record says only that *some* gate produced this exact leaf.** It is not
/// evidence that the pass which produced it survived. `manufactured` is a
/// `&mut Vec` threaded through `normalize_core`, and at least two callers throw
/// their *results* away while keeping the records: `normalize_allele`'s
/// warnings-only loop on the raw-conflict branch, which returns the authored
/// members verbatim, and `canonicalize_from_sequence`'s alternation loop, whose
/// non-converged intermediates are discarded. Records minted on those paths
/// still reach the seam.
///
/// That is safe because matching is by **leaf value**, and value equality is a
/// strong statement: the output leaf must be byte-identical to one a `#670`
/// junction gate produced, on the same accession. So an authored leaf can only
/// be repaired by coinciding exactly with a manufactured one — at which point
/// the two descriptions are the same string and the provenance question has no
/// observable answer.
///
/// **Not evidenced.** No input was constructed in which an *authored* intronic
/// leaf coincides with a manufactured one on the same accession, so this is an
/// unverified residual rather than a demonstrated hazard. It is stated here
/// because the sibling half of the doc — "an unmatched record leaves the output
/// exactly as `#1704` left it" — addresses only the *moved-leaf* direction and
/// reads as though it covered both.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ManufacturedJunctionExit {
    /// The leaf exactly as the gate produced it.
    ///
    /// Identity is by value rather than by index because there is no stable
    /// index: the sibling passes below `normalize_allele` may reorder, merge or
    /// drop members after the gate has run. A leaf that a later pass *moved* no
    /// longer matches, and that is the safe direction — an unmatched record
    /// leaves the output exactly as `#1704` left it rather than wrapping
    /// something whose provenance is no longer known.
    leaf: HgvsVariant,
    /// The bare transcript accession the leaf is rendered against.
    bare: Accession,
    /// `bare` carrying the genomic context the crossing resolved —
    /// `NC_…(NM_…)`, the form `checklist.md:20` names.
    wrapper: Accession,
}

/// Every leaf that is the `checklist.md:20` form — an intronic position named on
/// a **bare** transcript reference — in document order.
///
/// Delegates to [`intronic_on_bare_transcript_axis`] per leaf rather than
/// re-stating its scope, so the clause has exactly one reading in this crate —
/// the same reading [`intronic_on_bare_transcript_warning`] reports, since that
/// warning is derived from this predicate rather than written out beside it.
/// That matters in both directions: the same predicate decides whether strict
/// mode *refuses an input* (the W4007 rung of `normalize_core_checked`'s ladder)
/// and whether `#1704` must *re-parent an output*. Two readings of one clause is
/// how ferro came to refuse a description in strict mode while manufacturing the
/// identical description in lenient.
///
/// **Every leaf, not the first.** `#1704` reached for the offending leaf with a
/// `find_map`, so it could only ever repair one accession — and it additionally
/// asked an ANY-leaf existence question of the *input* as a provenance proxy,
/// which is the granularity mismatch `#1723` removes. Its own doc comment named
/// the first-leaf consequence for a `products` allele and called it silent; it is
/// not confined to `products` — an ordinary two-accession allele reaches it, and
/// `a_second_bare_accession_is_repaired_too` pins that shape.
fn bare_transcript_intronic_leaves(variant: &HgvsVariant) -> Vec<&HgvsVariant> {
    let mut found = Vec::new();
    collect_bare_transcript_intronic_leaves(variant, &mut found);
    found
}

fn collect_bare_transcript_intronic_leaves<'a>(
    variant: &'a HgvsVariant,
    found: &mut Vec<&'a HgvsVariant>,
) {
    match variant {
        HV::Allele(allele) => {
            for member in &allele.variants {
                collect_bare_transcript_intronic_leaves(member, found);
            }
        }
        HV::Supernumerary(inner) => collect_bare_transcript_intronic_leaves(inner, found),
        // The predicate, not the warning constructor: this runs on every
        // `normalize()` and only wants the answer, while building the warning
        // would render the description into a `String` and drop it.
        leaf => {
            if intronic_on_bare_transcript_axis(leaf).is_some() {
                found.push(leaf);
            }
        }
    }
}

/// Re-parent every `c.`/`n.` leaf of `variant` whose accession is `bare` onto
/// `wrapper`, which is `bare` carrying a genomic context.
///
/// **Every leaf on that accession, not only the intronic ones.** An allele
/// renders in its compact form (`ACC:c.[a;b]`) only when its members share an
/// accession — `use_compact_form` / `all_share_accession_and_type` — so wrapping
/// the crossing member alone would drop `c.[10+2del;18del]` into the expanded
/// `[NC_…(NM_…):c.10+2del;NM_…:c.18del]`. That is a far larger representation
/// change than the defect it repairs, and no clause asks for it: a genomic
/// wrapper on an *exonic* position is unremarkable
/// (`NC_000023.10(NM_004006.2):c.94del`), so lifting it to the whole description
/// is the cheap side of the trade.
///
/// **That argument is about EXONIC members and does not settle the mixed case.**
/// Where the sibling is an *authored* intronic position, lifting re-spells a leaf
/// the `bare-transcript-intronic-position` ruling says to leave as authored, so
/// the two answers are a genuine policy choice rather than a cheap side. See
/// [`Normalizer::reparent_junction_exit`] and the `undecided`
/// `junction-exit-wrapper-scope-in-a-mixed-allele` record.
fn reparent_leaves(variant: HgvsVariant, bare: &Accession, wrapper: &Accession) -> HgvsVariant {
    match variant {
        HV::Cds(mut v) if v.accession == *bare => {
            v.accession = wrapper.clone();
            HV::Cds(v)
        }
        HV::Tx(mut v) if v.accession == *bare => {
            v.accession = wrapper.clone();
            HV::Tx(v)
        }
        HV::Allele(allele) => HV::Allele(AlleleVariant {
            variants: allele
                .variants
                .into_iter()
                .map(|member| reparent_leaves(member, bare, wrapper))
                .collect(),
            ..allele
        }),
        HV::Supernumerary(inner) => {
            HV::Supernumerary(Box::new(reparent_leaves(*inner, bare, wrapper)))
        }
        other => other,
    }
}

/// Record every `checklist.md:20` leaf the `#670` junction gate just produced,
/// against the contig the crossing was computed on.
///
/// Called from all three copies of that gate — `normalize_cds`, `normalize_tx`
/// and `normalize_rna` — at the point the crossed answer is adopted. Reaching
/// that point is itself the proof of manufacture: an input naming an intronic
/// endpoint is diverted to `normalize_intronic_*` / `normalize_boundary_spanning_*`
/// well before the gate, so both endpoints were exonic and the intronic offset
/// in `produced` is one ferro made. See [`ManufacturedJunctionExit`].
///
/// **Scope is delegated, not restated.** The `r.` gate calls this too and it
/// records nothing there, because `intronic_on_bare_transcript_warning` — the
/// single reading of the clause in this crate — does not reach the `r.` axis
/// (`#486`/`#834`). Wiring the call anyway is deliberate: if that scope is ever
/// widened, the provenance follows it instead of silently staying behind.
fn record_manufactured_junction_exits(
    produced: &HgvsVariant,
    chromosome: Option<&str>,
    manufactured: &mut Vec<ManufacturedJunctionExit>,
) {
    let Some(chromosome) = chromosome else {
        return;
    };
    for leaf in bare_transcript_intronic_leaves(produced) {
        let Some(bare) = leaf.accession() else {
            continue;
        };
        let Some(wrapper) = genomic_wrapper_for(bare, chromosome) else {
            continue;
        };
        manufactured.push(ManufacturedJunctionExit {
            leaf: leaf.clone(),
            bare: bare.clone(),
            wrapper,
        });
    }
}

/// Build `NC_…(NM_…)` from the contig name the `#670` crossing resolved, or
/// decline.
///
/// # Guards on data, not on a re-derivation
///
/// `#1704` ran these three conditions at the top of the pipeline against a
/// *second* provider lookup. They now run at the gate, against the very
/// `Transcript` whose `chromosome` field the genomic re-shuffle fetched bases
/// from — so the reference this returns is the one the coordinate was computed
/// in, not one re-derived and hoped to agree.
///
/// # What each condition actually excludes — measured, not assumed
///
/// - **The bare parse.** `chromosome` is a provider-supplied string, not
///   necessarily an accession. A SAM-style `chr17` or an assembly name `GRCh38`
///   is refused *here*, by `parse_accession` failing on the bare string. It is
///   **not** refused by the pairing rule below: `is_valid_compound_outer`
///   deliberately admits an unclassifiable custom accession
///   (`inferred_variant_type().is_none()`, `#1146`) so a custom reference can
///   still carry a specification, and ferro's own parser reads
///   `chr17(NM_TEST.1):c.10+2del` back happily. Pinned by
///   `a_sam_style_chromosome_is_excluded_by_the_bare_accession_parse`.
/// - **`rest.is_empty()`.** `parse_accession` is a `nom` parser and stops at the
///   first byte it cannot use, so `NC_SYNTH.1junk` parses to `NC_SYNTH.1` with a
///   trailing remainder. Without this, a garbage suffix would be silently
///   discarded rather than declined.
/// - **`is_valid_compound_outer`.** Refuses a *transcript* named as the outer
///   reference — `NM_OTHER.1(NM_TEST.1)` is backwards and `parse_hgvs` will not
///   read it back, which is the property `FERRO_ASSERT_REPARSE` checks at this
///   same seam.
/// - **`outer == *bare`.** Refuses a self-referential wrapper, and it is **not**
///   subsumed by the condition above. That was assumed once and it is false: the
///   `checklist.md:20` predicate admits LRG by prefix (`Accession::is_lrg`, so a
///   bare `LRG_<N>` as well as `LRG_<N>t<M>`), while `is_valid_compound_outer`
///   keys off `inferred_variant_type`, which reads a bare `LRG_<N>` as
///   **genomic** and admits it. For `bare == LRG_5` the two conditions above both
///   pass and only this one stops `LRG_5(LRG_5):c.…`. Nothing downstream would —
///   that string re-parses, so `FERRO_ASSERT_REPARSE` is blind to it — and the
///   provider shape is not contrived, since an LRG record *is* its own genomic
///   reference and naming itself in `Transcript::chromosome` is reasonable.
///   Pinned by `a_self_referential_wrapper_is_declined` and, on this function
///   directly, by `the_genomic_wrapper_builder_declines_what_it_cannot_justify`.
fn genomic_wrapper_for(bare: &Accession, chromosome: &str) -> Option<Accession> {
    let (rest, outer) = crate::hgvs::parser::accession::parse_accession(chromosome).ok()?;
    if !rest.is_empty()
        || !crate::hgvs::parser::accession::is_valid_compound_outer(&outer)
        || outer == *bare
    {
        return None;
    }
    Some(bare.clone().with_genomic_context(outer))
}

/// The number of cis members a variant describes.
///
/// A cis allele reports its member count; everything else is one member. A
/// `trans`/unphased group is deliberately NOT counted — its members are on
/// different molecules, so they were never candidates for coalescing and a
/// change in their number would mean something else entirely. An `uncertain`
/// group is excluded for the same reason: `(...)` marks a predicted grouping
/// rather than a set of separately-reported observations.
fn cis_member_count(variant: &HgvsVariant) -> usize {
    match variant {
        HgvsVariant::Allele(allele)
            if allele.phase == crate::hgvs::variant::AllelePhase::Cis && !allele.uncertain =>
        {
            allele.variants.len()
        }
        _ => 1,
    }
}

/// The accession and coordinate-system letter to report a coalesce against.
///
/// Reads the first leaf carrying an accession, so a cis allele reports the
/// accession its members share rather than `None`.
fn accession_and_axis(variant: &HgvsVariant) -> Option<(String, String)> {
    let leaf = crate::hgvs::variant::first_leaf_with_accession(variant)?;
    let accession = leaf.accession()?.to_string();
    let axis = match leaf {
        HgvsVariant::Genome(_) => "g",
        HgvsVariant::Cds(_) => "c",
        HgvsVariant::Tx(_) => "n",
        HgvsVariant::Rna(_) => "r",
        HgvsVariant::Mt(_) => "m",
        // Protein too: a `p.` cis allele coalesces just as the nucleotide axes
        // do — `p.[Ala100Gly;Ala101Gly]` normalizes to
        // `p.Ala100_Ala101delinsGlyGly`, two reported members becoming one —
        // and the provenance `DNA/delins.md:79-84` is about is exactly as
        // unrecoverable from that string. Omitting it dropped the warning
        // silently on the axis where a database query is most likely to be
        // looking for the individually reported form.
        HgvsVariant::Protein(_) => "p",
        _ => return None,
    };
    Some((accession, axis.to_string()))
}

impl NormalizationWarning {
    /// The warning's user-facing code string.
    pub fn code(&self) -> &'static str {
        match self {
            Self::RefSeqMismatch { .. } => "REFSEQ_MISMATCH",
            Self::OverlapConflict { .. } => "OVERLAP_CONFLICTING_EDITS",
            Self::MembersCoalesced { .. } => "MEMBERS_COALESCED_FROM_REPORTED_FORM",
            Self::UnresolvableSpecialPosition { .. } => "UNRESOLVABLE_SPECIAL_POSITION",
            Self::TranscriptFlankNotDescribable { .. } => "TRANSCRIPT_FLANK_NOT_DESCRIBABLE",
            Self::CanonicalSplitSkipped { .. } => "CANONICAL_SPLIT_SKIPPED",
            Self::CrossAxisVariantNotShuffled { .. } => "CROSS_AXIS_VARIANT_NOT_SHUFFLED",
            Self::AxisClampApplied { .. } => "AXIS_CLAMP_APPLIED",
            Self::InitiatorMetCanonicalization { .. } => "INITIATOR_MET_CANONICALIZATION",
            Self::InsertedSequenceExpanded { .. } => "INSERTED_SEQUENCE_EXPANDED",
            Self::PositionPastEnd { .. } => "POSITION_PAST_END",
            Self::IntronicOnBareTranscript { .. } => "INTRONIC_ON_BARE_TRANSCRIPT",
            Self::IncompleteCdsStartReference { .. } => "INCOMPLETE_CDS_START_REFERENCE",
            Self::ReducedCapabilityNoGenome { .. } => "REDUCED_CAPABILITY_NO_GENOME",
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
            Self::MembersCoalesced {
                accession,
                coordinate_system,
                reported_members,
                normalized_members,
            } => write!(
                f,
                "{accession}:{coordinate_system}. — input described {reported_members} cis \
                 members, normalized form describes {normalized_members}; the individually \
                 reported form is not recoverable from the normalized string \
                 (DNA/delins.md:79-84)"
            ),
            Self::UnresolvableSpecialPosition { accession, marker } => write!(
                f,
                "{accession}: '{marker}' — centromere/telomere marker cannot be \
                 resolved to a coordinate without assembly annotation; input preserved verbatim \
                 (UnresolvableSpecialPosition)"
            ),
            Self::TranscriptFlankNotDescribable { accession, marker } => write!(
                f,
                "{accession}: '{marker}' on a genomic-reference c. description denotes a 5'/3' \
                 transcript-flank position, which HGVS does not permit numbering in c. coordinates; \
                 use the genomic g. form (TranscriptFlankNotDescribable / W4006)"
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
            Self::IntronicOnBareTranscript {
                variant,
                coordinate_system,
            } => {
                // Show a parent-reference example matching the input axis: a
                // coding `NM_` parent for c., a non-coding `NR_` parent for n.
                let parent_example = if coordinate_system == "n" {
                    "NG_(NR_)/NC_(NR_)"
                } else {
                    "NG_(NM_)/NC_(NM_)"
                };
                write!(
                    f,
                    "{variant}: intronic offset on a bare {coordinate_system}. transcript \
                     reference; a genomic reference sequence is required (e.g. {parent_example}) \
                     (IntronicOnBareTranscript / W4007 / EINTRONIC)",
                )
            }
            Self::IncompleteCdsStartReference {
                accession,
                coordinate_system,
            } => write!(
                f,
                "{accession}: transcript has an incomplete (unconfirmed ATG) CDS start \
                 (cds_start_NF); not an HGVS-recommended reference for {coordinate_system}. \
                 description — use the genomic (g.) or non-coding (n.) representation instead \
                 (IncompleteCdsStartReference / W5004)",
            ),
            Self::ReducedCapabilityNoGenome {
                variant,
                capability,
            } => write!(
                f,
                "{variant}: {capability} requires genomic sequence data, which the configured \
                 reference provider does not carry; returning the best-effort result unchanged \
                 (ReducedCapabilityNoGenome)",
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
    /// position rule. Every shipped path shuffles 3' (the rightmost form the
    /// HGVS recommendations mandate); the direction is still carried
    /// explicitly because ferro's own tests drive a 5' arm as a differential
    /// oracle — see [`config::ShuffleDirection`]. Code: `SHUFFLE_APPLIED`.
    /// Mutalyzer-equivalent: `ICORRECTEDPOINT` (mutalyzer only emits 3').
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

/// Main normalizer struct
pub struct Normalizer<P: ReferenceProvider> {
    provider: P,
    config: NormalizeConfig,
}

/// Resolve a special `c.` position to a concrete [`CdsPos`] using the transcript's
/// CDS structure. `pter` maps to transcript position 1 projected into c. coords;
/// `qter` maps to the last transcript position projected into c. coords; `cen` returns
/// `Ok(None)` (a centromere has no coordinate on a transcript); a plain (non-special)
/// position is returned as-is inside `Some`.
///
/// Returns `Ok(None)` when the projection cannot be computed (non-coding / missing
/// CDS metadata) so the caller can fall back to canonicalization. The caller can
/// distinguish `cen` (warn/reject) from a projection gap by re-inspecting
/// `pos.special`.
fn resolve_special_cds_pos(
    pos: &CdsPos,
    transcript: &crate::reference::transcript::Transcript,
) -> Result<Option<CdsPos>, FerroError> {
    use crate::convert::mapper::CoordinateMapper;
    use crate::hgvs::location::{SpecialPosition, TxPos};

    let special = match pos.special {
        None => return Ok(Some(*pos)),
        Some(s) => s,
    };
    let tx_pos = match special {
        SpecialPosition::Pter => TxPos::new(1),
        // NOTE: `sequence_length()` is the transcript-coordinate space that
        // `tx_to_cds` consumes — it counts mRNA bases (exon-sum fallback when
        // no cached sequence) so qter maps to the last exonic base, not past it.
        SpecialPosition::Qter => TxPos::new(transcript.sequence_length() as i64),
        SpecialPosition::Cen => return Ok(None),
    };
    let mapper = CoordinateMapper::new(transcript);
    // A projection error (non-coding / missing CDS) is a graceful fallback.
    Ok(mapper.tx_to_cds(&tx_pos).ok())
}

/// Sort cis-allele members into genomic (coordinate) order using the total
/// order `cis_member_order_key` (#1098). A no-op unless every member shares a
/// single accession — a mixed-accession bracketed allele (`[NC_…;NM_…]`,
/// #218/#219) has no cross-molecule genomic order to canonicalize to, so those
/// are left in authored order. The key is total (the rendered descriptor is the
/// final tie-break), so the result never depends on input order even when two
/// members share a start point.
///
/// Callers gate on the cis / not-uncertain / not-overlap-conflicting conditions
/// (see `normalize_allele`): this is used both *before* the merge, so
/// `merge_consecutive_edits` fires deterministically regardless of input order
/// (#1103), and *after* it, to render the surviving members in genomic order
/// (#1098/#1101).
fn sort_cis_members_by_genomic_order(members: &mut [HgvsVariant]) {
    let first_accession = members
        .first()
        .and_then(|m| m.accession().map(|a| a.full()));
    let single_accession = first_accession.is_some()
        && members
            .iter()
            .all(|m| m.accession().map(|a| a.full()) == first_accession);
    if single_accession {
        members.sort_by(|a, b| {
            crate::hgvs::variant::cis_member_order_key(a)
                .cmp(&crate::hgvs::variant::cis_member_order_key(b))
        });
    }
}

/// Whether the opt-in normalization idempotency self-check is active.
///
/// Enabled when the `FERRO_ASSERT_IDEMPOTENT` environment variable is set to
/// anything other than `0`/empty. Read once and cached, so the per-normalization
/// cost when disabled is a single relaxed atomic load. Only compiled in debug
/// builds; see the call site in [`Normalizer::normalize_core_checked`].
#[cfg(debug_assertions)]
fn idempotency_self_check_enabled() -> bool {
    use std::sync::OnceLock;
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| {
        std::env::var("FERRO_ASSERT_IDEMPOTENT")
            .map(|v| !v.is_empty() && v != "0")
            .unwrap_or(false)
    })
}

/// Whether the opt-in re-parse self-check is active.
///
/// Enabled when `FERRO_ASSERT_REPARSE` is set to anything other than
/// `0`/empty. Read once and cached, like the idempotency gate beside it.
///
/// Separate from `FERRO_ASSERT_IDEMPOTENT` on purpose. The two oracles overlap
/// but neither subsumes the other, and the idempotency one has a blind spot the
/// re-parse one covers: it verifies by *re-normalizing* its own output, which it
/// cannot do for an output that fails to parse. An unparseable result is
/// therefore invisible to it, so folding this in under the same flag would hide
/// which invariant actually broke. Both are on in CI.
#[cfg(debug_assertions)]
fn reparse_self_check_enabled() -> bool {
    use std::sync::OnceLock;
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| {
        std::env::var("FERRO_ASSERT_REPARSE")
            .map(|v| !v.is_empty() && v != "0")
            .unwrap_or(false)
    })
}

/// Whether the opt-in in-bounds self-check is active.
///
/// Enabled when `FERRO_ASSERT_IN_BOUNDS` is set to anything other than
/// `0`/empty. Read once and cached, like the two gates above it.
///
/// The third of a set, and separate from both for the same reason they are
/// separate from each other: neither existing oracle asks this question
/// *directly*. `parse_hgvs` has no provider and so cannot bounds-check, which
/// makes an out-of-range output perfectly re-parseable; and while the
/// idempotency oracle does happen to catch some instances — the #1307 output is
/// not a fixed point — that is incidental, not the invariant. An overrun that
/// re-normalizes to itself would pass both.
///
/// Cheaper than either, too: one length lookup per member, no re-normalization
/// and no parse.
#[cfg(debug_assertions)]
fn in_bounds_self_check_enabled() -> bool {
    use std::sync::OnceLock;
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| {
        std::env::var("FERRO_ASSERT_IN_BOUNDS")
            .map(|v| !v.is_empty() && v != "0")
            .unwrap_or(false)
    })
}

/// Whether the opt-in denoted-sequence self-check is active.
///
/// Enabled when `FERRO_ASSERT_SEQUENCE` is set to anything other than
/// `0`/empty. Read once and cached, like the three gates above it.
///
/// The fourth at the same seam (#1615), and the only one that asks what the
/// output *means* rather than how it is written. The three before it are all
/// well-formedness questions — is it a fixed point, does it parse, is it in
/// bounds — and a description denoting entirely different bases can satisfy all
/// three at once. That is not hypothetical: #1592 and #1600 were both live on
/// `main`, both emitted different bases than their input, and both issues record
/// all three oracles passing on the reproducer.
///
/// The most expensive of the four when enabled — a provider fetch of the union
/// window plus two splices per normalization — which is why it runs last in
/// [`Normalizer::assert_seam_oracles`] and why it is its own flag rather than
/// being folded into one of the others.
#[cfg(debug_assertions)]
fn denoted_sequence_self_check_enabled() -> bool {
    use std::sync::OnceLock;
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| {
        std::env::var("FERRO_ASSERT_SEQUENCE")
            .map(|v| !v.is_empty() && v != "0")
            .unwrap_or(false)
    })
}

/// Normalizations the denoted-sequence oracle looked at but could not compare.
///
/// A skip that reads as a pass is the exact failure mode a sequence oracle
/// exists to remove, so the skips are *counted* rather than being silent. A run
/// that reports zero comparisons made is a run that checked nothing, however
/// green it looks. See [`denoted_sequence_oracle_counts`].
#[cfg(debug_assertions)]
static DENOTED_SEQUENCE_SKIPPED: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0);

/// Normalizations the denoted-sequence oracle compared and found to agree.
#[cfg(debug_assertions)]
static DENOTED_SEQUENCE_COMPARED: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0);

/// `(compared, skipped)` for the denoted-sequence oracle, process-wide.
///
/// `compared` counts the normalizations whose input and output sequences were
/// actually derived and found equal; `skipped` counts those the comparison
/// declined, for any of the [`crate::spdi::NotComparable`] reasons. A mismatch
/// panics rather than being counted, so there is no third figure.
///
/// Exposed so a caller can assert its run made comparisons at all. Both counters
/// stay at zero when `FERRO_ASSERT_SEQUENCE` is unset — the oracle returns before
/// touching them — so a nonzero `compared` also witnesses that the flag is on.
///
/// Compiled out in release, like the oracle itself.
#[cfg(debug_assertions)]
pub fn denoted_sequence_oracle_counts() -> (u64, u64) {
    use std::sync::atomic::Ordering;
    (
        DENOTED_SEQUENCE_COMPARED.load(Ordering::Relaxed),
        DENOTED_SEQUENCE_SKIPPED.load(Ordering::Relaxed),
    )
}

#[cfg(debug_assertions)]
thread_local! {
    /// Re-entrancy guard for the idempotency self-check's verification pass.
    ///
    /// The check lives in `normalize_core_checked` so it covers the projector
    /// too, but its verification pass re-enters that very method — without this
    /// flag each check would recurse forever. Thread-local (not a global) so
    /// concurrent normalizations on other threads stay checked.
    ///
    /// Read by **two** oracles, for different reasons. `assert_idempotent` must
    /// honour it or it recurses; `assert_denoted_sequence` honours it because the
    /// re-entrant pass asks a question the outer call already asked, at the price
    /// of a second provider fetch and a counter the outer call did not ask for.
    /// `assert_reparseable` and `assert_in_bounds` deliberately do not: they are
    /// pure predicates over one description, hold no counters, and re-running
    /// them on the verification pass costs nothing worth naming.
    static IN_IDEMPOTENCY_CHECK: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
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

    /// Get a reference to the underlying reference provider.
    pub fn provider(&self) -> &P {
        &self.provider
    }

    /// Apply this variant to the reference and return the window before and
    /// after (#1159).
    ///
    /// The ground truth for equivalence: two descriptions denote the same edit
    /// exactly when they produce the same resulting bases. Use
    /// [`canonical_spdi`](Self::canonical_spdi) when a compact comparable key is
    /// wanted instead of the bases themselves.
    ///
    /// # Errors
    ///
    /// Declines rather than guessing for anything without one well-defined
    /// resulting sequence — a multi-molecule or null allele, members on different
    /// accessions, members that overlap, an edit SPDI cannot represent, or a span
    /// wider than [`crate::spdi::MAX_APPLY_WINDOW`]. See
    /// [`crate::spdi::apply_to_reference`].
    pub fn apply_to_reference(
        &self,
        variant: &HgvsVariant,
    ) -> Result<crate::spdi::AppliedVariant, FerroError> {
        crate::spdi::apply_to_reference(variant, &self.provider)
    }

    /// The reference/alternate pair this variant denotes, in the shape
    /// [`crate::from_sequences`] consumes.
    ///
    /// The inverse of [`Self::from_sequences`]. Where
    /// [`Self::apply_to_reference`] answers "what bases does this describe", in
    /// SPDI's 0-based frame, this answers the same question 1-based — the
    /// convention `VcfRecord`, HGVS and `from_sequences` all use. The two are
    /// otherwise the same computation.
    ///
    /// `pad` extra reference bases are fetched past the members' own span, to
    /// give a subsequent derivation room for its 3' placement. **128 is the
    /// value the normalizer's own derivation uses** (`merge::CANONICAL_PAD`);
    /// with `pad = 0` you get the bare changed block, and every
    /// `from_sequences` call on it will report `placement_bounded_by_window`, because
    /// with no flank every member is against an edge.
    ///
    /// **Both sides are padded by `pad`**, so the window a caller sizes from
    /// this doc is `span + 2 * pad`, not `span + pad`. This paragraph used to
    /// say the opposite — "only the 3' side is padded, and that is deliberate
    /// rather than an oversight" — four lines above a body that calls
    /// [`Self::prepend_five_prime_flank`]. The 3'-only argument it was quoting
    /// belongs to [`crate::spdi::apply_to_reference`] and is sound *there*:
    /// prepending `n` bases of 5' flank adds exactly `n` to the common prefix,
    /// so it cannot move an SPDI key. It does not carry to this method, whose
    /// consumer is `dup` typing rather than a trim — see the comment in the body
    /// for the ten corpus classes that measured the difference.
    ///
    /// The bases are returned **upper-cased**. `fetch_window` does not
    /// case-fold, so a soft-masked region arrives lower-case while the
    /// prepended flank was upper-cased — manufacturing a mixed-case pair that
    /// no caller supplied. Folding the whole window keeps
    /// `to_sequences` -> `from_sequences` an inverse on a masked region.
    ///
    /// Declines whatever [`Self::apply_to_reference`] declines: a
    /// multi-molecule or null allele, members on different accessions, members
    /// that overlap, an edit SPDI cannot represent, or a span wider than
    /// [`crate::spdi::MAX_APPLY_WINDOW`].
    pub fn to_sequences(
        &self,
        variant: &HgvsVariant,
        pad: u64,
    ) -> Result<crate::normalize::sequence_pair::SequencePair, FerroError> {
        let (applied, window_is_final) =
            crate::spdi::apply::apply_to_reference_padded(variant, &self.provider, pad)?;

        // `apply_to_reference_padded` pads only the 3' side, on the argument
        // that 5' flank adds equally to both strings' common prefix and so
        // cannot change what a *trim-based* consumer sees. That argument is
        // sound for the trim and wrong for the consumer this method serves:
        // `dup` typing reads the reference bases immediately 5' of an
        // insertion point (`duplication.md:18`), so a member flush with the
        // window's 5' edge cannot be recognised as a duplication and comes back
        // as an `ins` instead.
        //
        // Measured: over the cis confluence corpus, ten `all-dup` classes
        // reached both `g.[10_17dup;…]` and `g.[17_18insAATATATT;…]` — one
        // variant, two spellings — purely because the two inputs' member spans
        // gave them different amounts of 5' flank. Padding both sides closes it.
        //
        // Prepending the same bases to both strings cannot change what the pair
        // denotes, so this widens the window without touching the answer.
        let (reference, resulting, start) = self.prepend_five_prime_flank(
            &applied.accession,
            applied.start,
            applied.reference,
            applied.resulting,
            pad,
        );

        Ok(crate::normalize::sequence_pair::SequencePair {
            accession: applied.accession,
            // `AppliedVariant::start` is 0-based and every public surface here
            // is 1-based. This `+ 1` is the entire reason `SequencePair` exists
            // rather than handing `AppliedVariant` straight to `from_sequences`.
            position: start + 1,
            // Folded whole. `fetch_window` serves the provider's own case, so a
            // soft-masked region comes back lower-case, while the flank above is
            // upper-cased — the pair would otherwise be mixed-case in a way no
            // caller wrote. `from_sequences` folds too, so this is belt and
            // braces for the derivation; what it actually fixes is the *pair* a
            // caller reads, stores and compares.
            reference: reference.to_ascii_uppercase(),
            alternate: resulting.to_ascii_uppercase(),
            window_is_final,
        })
    }

    /// Move a window to `[start, end]`, padding from the reference or trimming,
    /// whichever each side needs.
    ///
    /// The reference-holding half of re-anchoring;
    /// [`crate::SequencePair::trim_to`] is the pure half and can only narrow.
    /// `None` leaves that edge where it is. Both bounds are **1-based
    /// inclusive**.
    ///
    /// # It moves a window's edges; it does not relocate the window
    ///
    /// The requested window must **overlap the pair's own**, and the overlap
    /// must still hold the bases the two sequences disagree on. Each edge may
    /// move outwards (padded from the reference) or inwards (trimmed), in either
    /// combination — but a window disjoint from the pair's is refused, not
    /// fetched. The changed bases exist only in the pair, so there is nothing to
    /// carry to a region the pair does not cover; a caller who wants the
    /// reference over some other interval wants `get_sequence`, not this.
    ///
    /// The bases come back **upper-cased**, for the reason
    /// [`Self::to_sequences`] folds its window: the flank is read from the
    /// provider and the rest is the caller's, so a soft-masked region would
    /// otherwise return a mixed-case pair that no caller wrote. `trim_to`, which
    /// fetches nothing and so cannot mix, leaves case alone.
    ///
    /// # What this is for
    ///
    /// Holding a derivation inside a region it must not leave — a target region,
    /// an amplicon, a tiling window — and doing so from whatever raw window each
    /// caller happens to have, **provided every such window overlaps the
    /// region**. Anchoring every input to one window makes them agree, because
    /// `from_sequences` is a pure function of the window it is given.
    ///
    /// # What this is NOT for
    ///
    /// Making heterogeneous raw pairs agree *in general*. That is already
    /// available, twice, and both are better answers:
    /// [`Self::from_sequences`] with `normalize = true`, or a round trip through
    /// [`Self::to_sequences`]. Both reach the **reference-anchored** placement,
    /// which can shift as far as the sequence allows rather than as far as a
    /// chosen window allows. Measured on three reads over one homopolymer
    /// deletion, both converge where the raw derivations do not.
    ///
    /// So reach for this when the bound is a *requirement*, not when it is a
    /// convenience. Anchoring to a window that cuts an ambiguous run makes every
    /// caller using that window agree with each other and disagree with the
    /// reference — legitimate as a stated contract, misleading as a default. The
    /// derivation reports it as `placement_bounded_by_window`.
    ///
    /// # Errors
    ///
    /// Refuses rather than clamping, always:
    ///
    /// * a bound outside the sequence — `start` of 0, or `end` past the last
    ///   base. A window is not silently pulled back to the contig, because a
    ///   caller who asked for bases that do not exist has a bug upstream and a
    ///   clamped window would hide it;
    /// * a **requested window disjoint from the pair's own** — see the section
    ///   above. This is the refusal a caller is most likely to meet, because
    ///   "re-anchor to my target region" reads like it should work from any
    ///   pair, and it works only from one the region overlaps;
    /// * a bound the provider cannot serve;
    /// * anything [`crate::SequencePair::trim_to`] refuses, for the narrowing
    ///   side — a cut through a base the two sequences disagree on, an emptied
    ///   reference, `start` past `end`.
    pub fn reanchor(
        &self,
        pair: &crate::normalize::sequence_pair::SequencePair,
        start: Option<u64>,
        end: Option<u64>,
    ) -> Result<crate::normalize::sequence_pair::SequencePair, FerroError> {
        let start = start.unwrap_or(pair.position);
        let end = end.unwrap_or(pair.end());
        let invalid = |msg: String| FerroError::InvalidCoordinates { msg };
        if start == 0 {
            return Err(invalid(
                "reanchor: start is 1-based; 0 does not name a base".to_string(),
            ));
        }
        if start > end {
            return Err(invalid(format!(
                "reanchor({start}, {end}) on {}: start is past end",
                pair.accession
            )));
        }

        let length = self.provider.get_sequence_length(&pair.accession)?;
        if end > length {
            return Err(FerroError::VariantExceedsReference {
                accession: pair.accession.clone(),
                hgvs_start: start,
                hgvs_end: end,
                expected_span: end - start + 1,
                actual_bytes: length,
            });
        }

        // Narrow first, then widen. Doing it in this order means the trim is
        // always over the caller's own bases — bases they observed — rather than
        // over bases this method just fetched, so a disagreement between the
        // supplied window and the reference cannot be trimmed away silently.
        //
        // Which makes the narrowing target the *intersection* of the requested
        // window and the pair's, and an empty intersection is a refusal this
        // method owes its own message. Delegating it to `trim_to` produced
        // `trim_to(1000, 16) on chr1: start is past end` — a sibling method's
        // name and two coordinates the caller never passed, for the mistake
        // callers make most often here.
        let (overlap_start, overlap_end) = (start.max(pair.position), end.min(pair.end()));
        if overlap_start > overlap_end {
            return Err(invalid(format!(
                "reanchor({start}, {end}) on {}: the requested window does not overlap the \
                 pair's own [{}, {}]. `reanchor` moves a window's edges over the reference; \
                 the changed bases exist only in the pair, so there is nothing to carry to a \
                 region the pair does not cover",
                pair.accession,
                pair.position,
                pair.end()
            )));
        }
        let narrowed = pair.trim_to(Some(overlap_start), Some(overlap_end))?;

        // Both edges are read before either string moves — `end()` is derived
        // from `reference`, so consuming it first would make the 3' test consult
        // a value that no longer exists.
        let (narrowed_start, narrowed_end) = (narrowed.position, narrowed.end());
        let mut reference = narrowed.reference;
        let mut alternate = narrowed.alternate;
        if start < narrowed_start {
            let flank = self.fetch_flank(&pair.accession, start, narrowed_start)?;
            reference.insert_str(0, &flank);
            alternate.insert_str(0, &flank);
        }
        if end > narrowed_end {
            let flank = self.fetch_flank(&pair.accession, narrowed_end + 1, end + 1)?;
            reference.push_str(&flank);
            alternate.push_str(&flank);
        }

        Ok(crate::normalize::sequence_pair::SequencePair {
            accession: pair.accession.clone(),
            position: start,
            // Folded for the same reason `to_sequences` folds: `fetch_flank`
            // upper-cases what it reads from the provider and the middle is
            // whatever the caller supplied, so a soft-masked window widened here
            // would come back mixed-case — a pair nobody wrote, and one that
            // compares unequal to both of its own halves.
            reference: reference.to_ascii_uppercase(),
            alternate: alternate.to_ascii_uppercase(),
            // Reaching the sequence's own end settles the 3' edge. So does
            // *leaving the pair's 3' edge where it was*, when that edge was
            // already settled — which is the only way a `true` ever enters this
            // method, since `to_sequences` is what produces one. Recomputing
            // `end == length` unconditionally threw that away, so an identity
            // `reanchor(pair, None, None)` and a 5'-only widen both downgraded a
            // settled window to `false`. The carry-through is `trim_to`'s
            // (`self.window_is_final && tail == 0`), stated over coordinates
            // because this method can widen as well as narrow.
            window_is_final: end == length || (end == pair.end() && pair.window_is_final),
        })
    }

    /// Reference bases over the half-open range `[from, to)`, 1-based.
    ///
    /// Refuses when the provider cannot serve them or serves the wrong number,
    /// because [`Self::reanchor`] promises an exact window and a short read here
    /// would silently produce a different one.
    fn fetch_flank(&self, accession: &str, from: u64, to: u64) -> Result<String, FerroError> {
        let wanted = (to - from) as usize;
        let bases = self.provider.get_sequence(accession, from - 1, to - 1)?;
        if bases.len() != wanted {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "reanchor: the provider served {} of the {wanted} bases requested over \
                     {accession}:{from}-{}",
                    bases.len(),
                    to - 1
                ),
            });
        }
        Ok(bases.to_ascii_uppercase())
    }

    /// Widen a window 5' by `pad` bases, returning the two strings and the new
    /// 0-based start.
    ///
    /// Best effort by design: near a sequence start there are simply fewer bases
    /// to prepend, and a short read there is not a failure — the window is still
    /// a correct description of its own extent. Returns the inputs unchanged
    /// when the provider cannot serve the flank, so a reduced-capability
    /// provider degrades to today's 3'-only behaviour rather than erroring.
    ///
    /// **The flank must be served in full or not at all**, which is the same
    /// discipline [`crate::spdi::apply::fetch_window`] applies and for the same
    /// reason. `get_sequence` is asked for `[start - wanted, start)`; a provider
    /// that serves fewer bases has served the **front** of that range, but the
    /// only way to accept a short read is to re-label those bases as the *last*
    /// `served` of it and report a new start of `start - served`. That silently
    /// shifts every derived coordinate by `wanted - served`, and the pair still
    /// looks well-formed, so nothing downstream can notice. Over-serving is the
    /// same fault mirrored: `start - served` underflows, and `[profile.release]`
    /// sets no `overflow-checks`, so a release build wraps to a colossal
    /// position instead of panicking.
    ///
    /// Declining an inexact read costs only the 5' flank on a provider that
    /// cannot serve it, which is precisely the degradation this function was
    /// already written to accept.
    fn prepend_five_prime_flank(
        &self,
        accession: &str,
        start: u64,
        reference: String,
        resulting: String,
        pad: u64,
    ) -> (String, String, u64) {
        let wanted = pad.min(start);
        if wanted == 0 {
            return (reference, resulting, start);
        }
        let flank_start = start - wanted;
        match self.provider.get_sequence(accession, flank_start, start) {
            Ok(flank) if flank.len() as u64 == wanted => {
                let flank = flank.to_ascii_uppercase();
                (
                    format!("{flank}{reference}"),
                    format!("{flank}{resulting}"),
                    flank_start,
                )
            }
            _ => (reference, resulting, start),
        }
    }

    /// [`crate::from_sequences`], against this normalizer's reference.
    ///
    /// The free function is a pure function of its arguments and therefore
    /// cannot range-check `position` — doing so needs the reference, which would
    /// make the provider a hidden input. This one holds a provider, so it does
    /// check, and refuses an interval running past the end of the sequence. It
    /// also offers `normalize`, which the free function cannot: normalizing
    /// needs the reference too.
    ///
    /// **Prefer `normalize = true` unless you have a reason not to.** Over a
    /// 6,000-shape sweep `normalize` moved **8.6%** of derived descriptions:
    /// repeat notation (`g.27_28insAAA` -> `g.27A[4]`), reference-anchored
    /// member re-derivation, and an inversion spread across several members
    /// (`g.[17C>A;19T>A;21T>G]` -> `g.17_21inv`, which the alignment DAG
    /// partitions before anything can see it, since it minimises edit distance
    /// and an inversion is not in that cost model).
    ///
    /// All three are rule 2 and rule 3 — the pair this design assigns to
    /// `normalize` — so `false` still yields a conformant, deterministic
    /// description. It is simply not the recommended form as often as one might
    /// assume. This doc previously said the flag was "a no-op on measured
    /// material"; that was measured on too narrow a corpus and is withdrawn.
    ///
    /// # Which seam oracles this reaches, stated because it is not all of them
    ///
    /// `normalize = true` routes the derived description through
    /// [`Self::normalize`] and therefore through `normalize_core_checked`, the
    /// single exit `assert_seam_oracles` runs from — so all four oracles
    /// (`FERRO_ASSERT_IDEMPOTENT`, `_REPARSE`, `_IN_BOUNDS`, `_SEQUENCE`) fire
    /// on it, and this is the **only** entry to the derivation that reaches
    /// them.
    ///
    /// `normalize = false`, the free [`crate::from_sequences`] and
    /// [`crate::SequencePair::derive`] reach none of them, and that is a
    /// deliberate limit rather than an omission to be closed later:
    ///
    /// - The three form oracles compare an output against **its input
    ///   description**. There is no input description here — the input is bases
    ///   — so there is nothing to be idempotent with respect to, nothing to
    ///   re-parse a "before" of, and no prior spelling whose coordinates could
    ///   be re-checked. (`_IN_BOUNDS` is the one that could be asked
    ///   independently, and this method answers it directly instead: the
    ///   range check above refuses an interval past the sequence end **before**
    ///   deriving, which the free function cannot do at all.)
    /// - `_SEQUENCE` compares the bases an input and an output denote. The
    ///   caller's `reference` is taken on trust and is explicitly not required
    ///   to be the reference, so comparing the derived description against the
    ///   held reference would report a caller's own choice as a fault.
    ///
    /// What the derivation has instead is `from_sequences`'
    /// `verify_round_trip`, and its doc says plainly what that is worth: it
    /// shares an applier with the derivation, so it is a self-consistency check.
    /// The independent comparison is made in the corpus and multi-member axes,
    /// both of which key the derived side through `hgvs_to_spdi`.
    pub fn from_sequences(
        &self,
        accession: &str,
        position: u64,
        reference: &str,
        alternate: &str,
        options: &crate::normalize::from_sequences::FromSequencesOptions,
        normalize: bool,
    ) -> Result<HgvsVariant, FerroError> {
        let length = self.provider.get_sequence_length(accession)?;
        // Checked before the derivation, not after: a description of bases that
        // are not in the sequence is invalid however well it is derived, and
        // saying so up front names the real fault.
        let end = position
            .checked_add(reference.len() as u64)
            .and_then(|end| end.checked_sub(1))
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: format!("position {position} plus the reference length overflows"),
            })?;
        if position > length || end > length {
            return Err(FerroError::VariantExceedsReference {
                accession: accession.to_string(),
                hgvs_start: position,
                hgvs_end: end,
                expected_span: reference.len() as u64,
                actual_bytes: length,
            });
        }
        let variant = crate::normalize::from_sequences::from_sequences(
            accession, position, reference, alternate, options,
        )?;
        if normalize {
            self.normalize(&variant)
        } else {
            Ok(variant)
        }
    }

    /// An encoding-invariant SPDI key for this variant, derived from the bases it
    /// results in rather than from how it was written (#1159).
    ///
    /// Two descriptions on one accession with the same resulting sequence give the
    /// same triple, whatever their spelling, member count or member order — which
    /// is what [`crate::spdi::hgvs_to_spdi`] cannot do, because it transliterates
    /// the caller's own partitioning. `g.8_14delinsGATTA` and
    /// `g.[8A>G;9G>A;11C>T;13_14del]` transliterate to one triple and four; they
    /// canonicalize to the same one.
    ///
    /// See [`crate::spdi::apply`]'s module docs for the exact guarantee, and in
    /// particular for what "canonical" does *not* claim here.
    ///
    /// # Errors
    ///
    /// As [`apply_to_reference`](Self::apply_to_reference).
    pub fn canonical_spdi(
        &self,
        variant: &HgvsVariant,
    ) -> Result<crate::spdi::SpdiVariant, FerroError> {
        crate::spdi::canonical_spdi(variant, &self.provider)
    }

    /// The reference-aware SPDI transliteration of this variant.
    ///
    /// Unlike [`crate::spdi::hgvs_to_spdi_simple`], which is what the Python
    /// `hgvs_to_spdi` used to be limited to, this resolves the edits that need the
    /// reference to know their bases — `del`, `delins`, `inv`, `dup` — instead of
    /// failing on them (#1159).
    ///
    /// This preserves the input's partitioning, so it is **not** an
    /// encoding-invariant key; a cis allele is not one triple at all and is
    /// refused here. Use [`canonical_spdi`](Self::canonical_spdi) for a key.
    ///
    /// # Errors
    ///
    /// Whatever [`crate::spdi::hgvs_to_spdi`] reports, as a [`FerroError`].
    pub fn to_spdi(&self, variant: &HgvsVariant) -> Result<crate::spdi::SpdiVariant, FerroError> {
        crate::spdi::hgvs_to_spdi(variant, &self.provider).map_err(|e| {
            FerroError::UnsupportedVariant {
                variant_type: format!("{variant}: cannot convert to SPDI — {e}"),
            }
        })
    }

    /// Normalize a variant — **the quiet exit**.
    ///
    /// This returns the normalized variant and **discards every
    /// [`NormalizationWarning`] the normalization raised**. The core produces
    /// them (`normalize_core_checked` returns a `(variant, warnings)` pair); the
    /// `.0` below throws the second half away. That is a deliberate signature
    /// choice — the return type has nowhere to put them — but it means a caller
    /// on this exit cannot learn that a repair happened.
    ///
    /// Several of those repairs are lossy in a way the output string does not
    /// record: `MEMBERS_COALESCED_FROM_REPORTED_FORM` (the individually reported
    /// members are gone), `INSERTED_SEQUENCE_EXPANDED` (the `[100_110]` payload
    /// is now a literal), `REFSEQ_MISMATCH` with `corrected=true` (the caller's
    /// stated base was wrong and was silently replaced). On this exit each of
    /// those is indistinguishable from a clean normalization.
    ///
    /// **Use [`normalize_with_diagnostics`](Self::normalize_with_diagnostics)
    /// unless you specifically want the quiet behaviour.** It routes through the
    /// same `normalize_core_checked`, so the normalized variant is
    /// byte-identical; it simply also hands back the warnings, plus the
    /// info-grade shuffle signals. `ferro normalize`, `ferro project` and the
    /// Python `normalize_with_warnings` bindings are all on that exit.
    ///
    /// `error_mode` is orthogonal and does **not** substitute for reading the
    /// warnings: strict mode promotes a specific ladder of conditions to hard
    /// errors, and every warning outside that ladder is raised identically in
    /// strict and lenient — so a strict `Normalizer` on this exit drops exactly
    /// the same diagnostics a lenient one does.
    ///
    /// # Errors
    ///
    /// In strict mode, rejects variants with reference mismatches and the rest
    /// of the `should_reject_*` ladder. In the default lenient mode those
    /// conditions become warnings — which this exit then discards.
    pub fn normalize(&self, variant: &HgvsVariant) -> Result<HgvsVariant, FerroError> {
        // Thin wrapper: `normalize()` returns only the variant. Both the core
        // normalization AND the strict-mode rejection ladder live in
        // `normalize_core_checked`, so a strict config rejects identically
        // whether a variant is normalized directly or through the projector.
        // The opt-in idempotency self-check lives in `normalize_core_checked`
        // (below), so it covers this method AND the projector.
        Ok(self.normalize_core_checked(variant)?.0)
    }

    /// Opt-in idempotency self-check (issue #1157 follow-up): a normalizer must
    /// satisfy `norm(norm(x)) == norm(x)`. Re-normalizes `normalized` and panics
    /// on any drift, so every test that normalizes becomes an idempotency oracle
    /// for the cost of one gated check.
    ///
    /// Called from [`Normalizer::normalize_core_checked`] — the shared core
    /// behind both `normalize()` and `VariantProjector` — so the projection-driven
    /// axes (genomic/coding/protein) are covered too. It previously sat in
    /// `normalize()` alone, which left every projector-only path unchecked.
    ///
    /// The verification pass re-enters `normalize_core_checked`, so
    /// `IN_IDEMPOTENCY_CHECK` breaks the recursion: the inner call skips its own
    /// check. The flag is cleared BEFORE asserting so a panic cannot leave it
    /// stuck set for later normalizations on this thread.
    ///
    /// Compiled out entirely in release (`debug_assertions`), so production is
    /// untouched and zero-cost.
    #[cfg(debug_assertions)]
    fn assert_idempotent(&self, variant: &HgvsVariant, normalized: &HgvsVariant) {
        if !idempotency_self_check_enabled() || IN_IDEMPOTENCY_CHECK.with(|f| f.get()) {
            return;
        }
        // RAII, not a bare `set(false)` after the call: if the verification pass
        // panics (or unwinds for any reason) a manual reset is skipped, leaving
        // the flag stuck `true` — which silently disables the oracle for every
        // later normalization on this thread. A self-check that can quietly turn
        // itself off is worse than one with a documented limit. The guard drops
        // at the end of the block, i.e. BEFORE the assert below, so the flag is
        // also clear on the ordinary non-idempotent-panic path.
        struct ReentrancyGuard;
        impl Drop for ReentrancyGuard {
            fn drop(&mut self) {
                IN_IDEMPOTENCY_CHECK.with(|f| f.set(false));
            }
        }

        let second = {
            IN_IDEMPOTENCY_CHECK.with(|f| f.set(true));
            let _guard = ReentrancyGuard;
            self.normalize_core_checked(normalized)
        };

        match second {
            Ok((again, _)) => assert_eq!(
                normalized.to_string(),
                again.to_string(),
                "FERRO_ASSERT_IDEMPOTENT: normalize is not idempotent\n  \
                 input: {variant}\n  once:  {normalized}\n  twice: {again}",
            ),
            Err(e) => panic!(
                "FERRO_ASSERT_IDEMPOTENT: normalized output failed to re-normalize\n  \
                 input: {variant}\n  once:  {normalized}\n  error: {e}",
            ),
        }
    }

    /// Assert that a normalized description is one `parse_hgvs` accepts.
    ///
    /// Not a canonicalization question. Whatever the right canonical form is,
    /// `normalize` is documented to return a valid HGVS description, and every
    /// consumer that chains a second call depends on that — including the
    /// idempotency oracle above, which re-normalizes its own output and so
    /// cannot even run on a string that will not parse. An unparseable result is
    /// invisible to that oracle, which is why this one is separate.
    ///
    /// The case that motivated it — a deletion and an
    /// insertion-turned-duplication claiming one position, which ferro emitted
    /// as `g.[261_262del;262dup;…]` and its own parser rejected as a
    /// self-cancelling allele — is fixed by `respell_colliding_duplications` in
    /// this change, and the reduced form is pinned by
    /// `tests/it/normalize_reparse_invariant.rs`. CI therefore sets
    /// `FERRO_ASSERT_REPARSE` alongside `FERRO_ASSERT_IDEMPOTENT`; the whole
    /// suite passes with it on.
    ///
    /// Compiled out in release, and cheaper than the idempotency check it sits
    /// beside — one parse rather than a whole re-normalization.
    #[cfg(debug_assertions)]
    fn assert_reparseable(&self, variant: &HgvsVariant, normalized: &HgvsVariant) {
        if !reparse_self_check_enabled() {
            return;
        }
        // `0` and `?` are legal whole-allele descriptions but are not
        // self-describing: `parse_hgvs` wants an accession, so it rejects them
        // standalone. That is a limit of the entry point, not a malformed
        // output, and exempting them here keeps the oracle's failures meaningful.
        if matches!(
            normalized,
            HgvsVariant::NullAllele | HgvsVariant::UnknownAllele
        ) {
            return;
        }
        let rendered = normalized.to_string();
        let Err(e) = crate::hgvs::parser::parse_hgvs(&rendered) else {
            return;
        };
        // Only fire when normalization *introduced* the breakage — but say
        // exactly which inputs are passed over, rather than exempting anything
        // that fails to re-parse (#1264).
        //
        // The blanket version ("skip when the input does not itself re-parse")
        // was doing far more work than its comment claimed. Instrumenting it to
        // report instead of return silently, and re-running the suite, found 18
        // hits in four shapes — including two live projector defects that were
        // building descriptions ferro could not read back: the RNA-only `spl`
        // edit carried onto the `g.`/`c.` axes, and an insertion whose anchors
        // straddle a splice junction. Both are fixed in this change; a blanket
        // exemption would have re-absorbed them, and any future defect of the
        // same class, without a sound.
        //
        // What remains is a closed list. Each entry is a shape ferro constructs
        // deliberately and does not claim is renderable:
        //
        // * An **empty allele** (`[]`). The parser refuses empty brackets
        //   everywhere they can appear, so this is reachable only by direct
        //   construction — which the projector's own tests do on purpose, to
        //   pin that it declines them. `AlleleVariant::checked` now refuses it
        //   at the boundaries where a member list is assembled from data.
        // * A **non-flanking genomic insertion**. This is the projection
        //   *pivot*: its coordinates are sound and the downstream cdot
        //   derivation of the c./n./p. axes needs it, but its spelling is not
        //   one HGVS admits, so the projector withholds the *reported* genomic
        //   axis (see `non_flanking_genomic_insertion_anchor`). Normalizing the pivot
        //   is deliberate; reporting it is what would be wrong.
        // * A **non-coding downstream position** (`n.*N`), which #1748 refuses at
        //   parse in every mode — `background/numbering.md:52` numbers that axis
        //   from the first nucleotide to the last and puts no `*` zone on it —
        //   while `TxPos::downstream` stays public API. So the spelling is gone
        //   and the AST is not, and a Rust library caller can still build one and
        //   hand it here. Normalization does not *produce* the shape: the two
        //   sites that can emit it — the `is_downstream` short-circuit below,
        //   which returns the input through `canonicalize_tx_variant` (it clones
        //   the location, edit only), and `build_tx_merged`'s `Region::TxDownstream`
        //   arm, whose region is minted at exactly one place gated on
        //   `pos.is_downstream()` — are both strictly flag-preserving, and no
        //   string entry point (parser, CLI, Python, VCF, SPDI, projector) can
        //   reach them. What arrives is the caller's own construction, so ferro
        //   does not claim it is renderable and the oracle has nothing to say
        //   about it.
        //
        // Anything else that arrives unparseable is now a failure, which is the
        // point.
        //
        // Note the third entry is keyed on the AST, via the same predicate the
        // parser refuses on, and never on the rendered string. That matters for
        // one specific reason: `noncoding_zone_marker` matches `HgvsVariant::Tx`
        // exhaustively and yields `None` for `Cds`, so `c.*N` — anchored to the
        // CDS, still legal, still parsing — cannot be swept in by it. A string
        // test for `*` would not have that property.
        let input_is_a_deliberate_non_renderable = matches!(
            variant,
            HgvsVariant::Allele(allele) if allele.variants.is_empty()
        )
            || crate::hgvs::variant::non_flanking_genomic_insertion_anchor(variant).is_some()
            || crate::hgvs::noncoding_zones::noncoding_zone_marker(
                variant,
                crate::hgvs::noncoding_zones::NonCodingZone::Downstream,
            )
            .is_some();
        if input_is_a_deliberate_non_renderable {
            return;
        }
        debug_assert!(
            crate::hgvs::parser::parse_hgvs(&variant.to_string()).is_ok(),
            "FERRO_ASSERT_REPARSE: normalization was handed an input it cannot re-parse, and \
             that input is not one of the deliberate non-renderable shapes\n  input: {variant}",
        );
        panic!(
            "FERRO_ASSERT_REPARSE: normalization produced a description ferro cannot re-parse\n  \
             input:  {variant}\n  output: {rendered}\n  error:  {e}",
        );
    }

    /// Assert that no coordinate in a normalized description names a position
    /// its own sequence does not have.
    ///
    /// The third seam oracle (#1353), and the one that asks the question
    /// directly. It exists because the class kept being found by hand, one shape
    /// at a time — #1274 (`T:g.[8_9insA;10del]` → `g.10_11=` on ten bases),
    /// #1343 (`c.[*10dup;*11dup]` → `c.*11_*12insAA`) and #1307
    /// (`g.[24dup;24C>G]` → `g.[24C>G;24_25insC]`) — each filed, fixed and
    /// regression-tested separately. Three instances of one defect class is the
    /// argument for an invariant at the seam rather than a fourth per-shape test.
    ///
    /// Neither existing oracle covers it. `FERRO_ASSERT_REPARSE` cannot:
    /// `parse_hgvs` holds no provider, so `g.24_25insC` is well-formed to it.
    /// `FERRO_ASSERT_IDEMPOTENT` catches some instances incidentally — the #1307
    /// output re-normalizes to `g.23_24insG`, so it is not a fixed point — but an
    /// out-of-range output that *is* a fixed point passes it, which is exactly
    /// the #1327 shape (`m.16569_16570insAA`).
    ///
    /// See [`merge::first_out_of_bounds_coordinate`] for what is checked, what is
    /// deliberately not, and why a reversed circular range is not an error.
    ///
    /// Compiled out in release, like the two beside it.
    #[cfg(debug_assertions)]
    fn assert_in_bounds(&self, variant: &HgvsVariant, normalized: &HgvsVariant) {
        if !in_bounds_self_check_enabled() {
            return;
        }
        let Some(overrun) = merge::first_out_of_bounds_coordinate(normalized, &self.provider)
        else {
            return;
        };
        // Only fire when normalization *introduced* it, matching
        // `assert_reparseable`'s discipline: an input that already names a
        // position off the end is a caller error, and W4004 `PositionPastEnd` is
        // the mechanism that reports it. Blaming normalization for preserving it
        // would make the oracle's failures mean two different things.
        if merge::first_out_of_bounds_coordinate(variant, &self.provider).is_some() {
            return;
        }
        let OutOfBoundsCoordinate {
            accession,
            position,
            length,
            member,
        } = overrun;
        panic!(
            "FERRO_ASSERT_IN_BOUNDS: normalization produced a coordinate the sequence does not \
             have\n  input:  {variant}\n  output: {normalized}\n  member: {member}\n  names \
             position {position} on {accession}, which has {length} bases",
        );
    }

    /// Assert that a normalized description denotes the same bases its input
    /// does.
    ///
    /// The fourth seam oracle (#1615), and the first that asks about *meaning*
    /// rather than form. Normalization may re-spell a variant however the
    /// canonical form requires; what it may never do is change the sequence the
    /// description resolves to. Nothing else at this seam checks that:
    ///
    /// * `FERRO_ASSERT_IDEMPOTENT` asks whether the output is a fixed point. A
    ///   wrong sequence normalizes to itself perfectly well — #1592, #1600,
    ///   #1304, #1308 and #1312 are all fixed points.
    /// * `FERRO_ASSERT_REPARSE` asks whether it parses. `parse_hgvs` holds no
    ///   provider, so it cannot know what any description denotes.
    /// * `FERRO_ASSERT_IN_BOUNDS` asks whether its coordinates exist. Every one
    ///   of the eight defects below names positions that do exist.
    ///
    /// The class had been found by hand six times — #1254, #1281, #1290, #1304,
    /// #1308, #1312 — each filed, fixed and regression-tested separately, before
    /// #1592 and #1600 made it eight. That is the same argument #1353 made for
    /// the in-bounds oracle after three.
    ///
    /// # The applier is not the normalizer
    ///
    /// [`crate::spdi::compare_denoted_sequences`] reaches the bases through
    /// `hgvs_to_spdi` and an SPDI splice, with no normalization anywhere in the
    /// path — otherwise the check would agree with the output merely because
    /// normalization produced it. `EquivalenceChecker` is *not* usable here for
    /// exactly that reason.
    ///
    /// # Skipping is counted, and an inapplicable output is not a skip
    ///
    /// Whenever the input itself denotes no single sequence — a trans allele, a
    /// `REFSEQ_MISMATCH` whose stated base is wrong, an edit SPDI cannot carry —
    /// there is no baseline, and the comparison declines. That is the same
    /// discipline `assert_reparseable` and `assert_in_bounds` apply, and it is
    /// deliberately *counted* ([`denoted_sequence_oracle_counts`]) rather than
    /// silent: a run that compared nothing looks identical to a run that
    /// compared everything and found no fault.
    ///
    /// An output that is **self-contradictory** while its input denotes a
    /// sequence is the opposite case and fires. #1281's `g.[1del;9_10delinsA]`
    /// normalized to `g.[1del;1del]`, two members claiming one base: the
    /// description has no resulting sequence at all, which is worse than having
    /// the wrong one. Treating it as a skip would file the more severe defect
    /// under the milder one's exemption.
    ///
    /// That fire is deliberately narrower than "the output produced no triples".
    /// An output the applier merely cannot *transliterate* is a skip, because
    /// the two sides do not need the reference equally: `g.1000delC` states its
    /// deleted base and converts with no provider at all, while the `g.1000del`
    /// normalization emits for it must read one. Reading that asymmetry as a
    /// fault raised 328 false alarms across the suite in a single measured run,
    /// against 16 genuine ones.
    ///
    /// # It does not run on the idempotency oracle's verification pass
    ///
    /// `assert_idempotent` re-enters `normalize_core_checked`, so this seam is
    /// reached a second time with `(normalized, again)` — a pair the outer call
    /// has already asked about, since the outer comparison covers
    /// `(input, normalized)` and the idempotency assertion covers
    /// `normalized == again`. Honouring `IN_IDEMPOTENCY_CHECK` therefore costs no
    /// coverage and buys three things: it halves the provider fetches when both
    /// flags are set (`sweeps` sets both), it keeps
    /// [`denoted_sequence_oracle_counts`] a count of *normalizations the caller
    /// asked for* rather than a mixture of those and re-entrant verification
    /// passes, and it stops a non-idempotent output being reported as a sequence
    /// fault by the inner call before `assert_idempotent` can name it as the
    /// non-idempotency it is.
    ///
    /// The return is before either counter, like the flag-off path: a re-entrant
    /// pass is not a skip the oracle declined, so counting it as one would put
    /// the same inflation into `skipped` that it removes from `compared`.
    ///
    /// Compiled out in release, like the three beside it, and the most expensive
    /// of the four when on — one provider fetch and two splices per
    /// normalization — so it runs last.
    #[cfg(debug_assertions)]
    fn assert_denoted_sequence(
        &self,
        variant: &HgvsVariant,
        normalized: &HgvsVariant,
        warnings: &[NormalizationWarning],
    ) {
        use crate::spdi::DenotedSequenceComparison as Outcome;
        use std::sync::atomic::Ordering;

        if !denoted_sequence_self_check_enabled() || IN_IDEMPOTENCY_CHECK.with(|f| f.get()) {
            return;
        }
        // A corrected `REFSEQ_MISMATCH` is a sequence change the normalizer is
        // *supposed* to make, and the only one there is. The input stated a base
        // the reference does not have — `c.10dupA` where the reference reads `C`
        // — and canonicalization drops the wrong bases, so the two descriptions
        // legitimately denote different sequences. Counted as a skip rather than
        // exempted silently, and kept narrow: only the `corrected` half, since a
        // mismatch the normalizer declines to rewrite leaves both sides denoting
        // whatever the input said.
        if warnings.iter().any(|w| {
            matches!(
                w,
                NormalizationWarning::RefSeqMismatch {
                    corrected: true,
                    ..
                }
            )
        }) {
            DENOTED_SEQUENCE_SKIPPED.fetch_add(1, Ordering::Relaxed);
            return;
        }
        match crate::spdi::compare_denoted_sequences(variant, normalized, &self.provider) {
            Outcome::Agree => {
                DENOTED_SEQUENCE_COMPARED.fetch_add(1, Ordering::Relaxed);
            }
            Outcome::NotComparable(_) => {
                DENOTED_SEQUENCE_SKIPPED.fetch_add(1, Ordering::Relaxed);
            }
            Outcome::OutputContradictsItself => panic!(
                "FERRO_ASSERT_SEQUENCE: normalization produced a description that denotes no \
                 sequence at all, from an input that denotes one\n  input:  {variant}\n  output: \
                 {normalized}\n  its members claim the same base, or two of its insertions share \
                 one interbase with no stated order",
            ),
            Outcome::Differ {
                accession,
                start,
                reference,
                from_input,
                from_output,
            } => panic!(
                "FERRO_ASSERT_SEQUENCE: normalization changed the sequence the description \
                 denotes\n  input:  {variant}\n  output: {normalized}\n  over {accession} \
                 [{start}, {}):\n    reference       {reference}\n    input applies   \
                 {from_input}\n    output applies  {from_output}",
                start + reference.len() as u64,
            ),
        }
    }

    /// Whether a lone (non-allele) member is a shape the sequence-first pass
    /// can improve.
    ///
    /// Any edit type may enter: the derivation is a function of (reference,
    /// denoted sequence), so restricting by the input's *spelling* — as this
    /// used to do, admitting only `delins`/`inv` — is exactly the
    /// input-relativity #1235 is about. A lone `g.263_264insGAAA` and the
    /// two-member `g.[261_262insGA;263_264insAA]` denote the same variant, and
    /// until this widened, only one of them was allowed to be re-derived.
    ///
    /// The **axis** gate is open to the transcript axes as of #1235: `c.`, `n.`
    /// and `r.` join `g.`/`m.` in the match below. #1237 had narrowed it to
    /// `g.`/`m.` on the grounds that `merge::canonicalize_from_sequence` refused
    /// the transcript axes anyway; a prior PR removed that refusal for
    /// **alleles**, which left a lone `c.`/`n.`/`r.` member as the one shape
    /// still held back on axis grounds. The asymmetry was real — a lone spelling
    /// that the allele spelling of the same variant would re-partition could
    /// stay unsplit — and closing it is what the widening does.
    ///
    /// # What widening the axis cost, measured rather than predicted
    ///
    /// An earlier revision of this comment held the axis shut behind four tests
    /// it predicted would move. **Exactly one of the four moved anything
    /// committed, and it was fixed rather than re-blessed** — so the conformance
    /// decision the comment was waiting on turned out not to exist.
    /// Measured twice: once with only the widening applied (this branch's first
    /// commit, `feat(normalize): widen the axis gate to c./n./r.`) and again at
    /// the branch tip. Commit subjects rather than hashes, because rebasing
    /// rewrites the hashes and a citation nobody can resolve is worse than none.
    ///
    /// * `normalize_tests::test_normalize_inversion_unchanged` and
    ///   `rna_coding_consistency::parity_safe_regime::inversion_parity` — **did
    ///   not move.** Both pass unmodified with only the widening applied and at
    ///   the tip; the predicted lone-`inv` split does not arise on the `main`
    ///   this branch sits on. Neither file is touched by this branch.
    /// * `spec_enumeration_tests::enumeration_replays_recorded_behavior` — the
    ///   gitignored recording does regenerate, but that test compares ferro
    ///   against itself; the two **committed** census guards are what judge
    ///   behaviour, and they do not move. `DIVERGENCE_BUDGET` and
    ///   `PASSING_CENSUS` are byte-identical to `main`, so the diff there is
    ///   commentary only. Per this branch's own two measurements
    ///   `ProjectionSplitsSingleMember` went 9 -> 13 and back to 9 inside the
    ///   branch; only the net is asserted anywhere.
    /// * `mutalyzer_normalize_tests::gate_normalized_snapshot` — the prediction
    ///   that held, and at exactly the stated count. With only the widening
    ///   applied it fails with **3 hermetic divergences**, one variant in three
    ///   spellings:
    ///   `NM_000143.3:c.44_47delinsATC` (also `…delTGCGinsATC`, `…del4insATC`)
    ///   rendered as `c.[44_45delinsAT;49del]`.
    ///
    /// **How those three were adjudicated: they were not.** The divergence is
    /// `general.md:35`'s coding exception — two changes one base apart within
    /// one reading frame — which the live path could not apply because it had no
    /// axis to apply it on. `merge::coalesce_coding_frame_separation` supplies
    /// it, and that function's own doc names this same `c.44_47delinsATC` as the
    /// case it was written for. At the tip `gate_normalized_snapshot` passes all
    /// **26** covered rows (0 fail, 0 known_bug, 0 improvement) with **no
    /// Mutalyzer snapshot, `cases.json` row or oracle baseline changed** — this
    /// branch touches no file under `tests/fixtures/`. So ferro agrees with the
    /// oracle again and there was no conformance decision left to take.
    ///
    /// **CORRECTION: they have since been adjudicated, the other way.** The
    /// paragraph above reads "the exception applies here" off the *distance*
    /// half of `general.md:35` alone. The clause has two conjuncts, and the
    /// second — "together affecting one amino acid" — is not met: `c.44_47` is
    /// four positions, spanning codon 15 (`c.43_45`) and codon 16 (`c.46_48`),
    /// and no codon holds four consecutive positions at any phase. Once
    /// `coalesce_coding_frame_separation` tests that conjunct, the three rows go
    /// back to `c.[44_45delinsAT;49del]` permanently and they now carry a
    /// `spec_citation` annotation in `cases.json` recording it as a deliberate
    /// spec-over-Mutalyzer divergence. Precedence is **spec-explicit >
    /// Mutalyzer**, and Mutalyzer has no separation rule at all — it minimises a
    /// weighted description length — so its merge here is two objectives meeting
    /// rather than evidence about the clause. Do not read this branch's
    /// "0 fail" as still current.
    fn is_splittable_single_member(variant: &HgvsVariant) -> bool {
        let edit = match variant {
            HV::Genome(g) => &g.loc_edit.edit,
            HV::Mt(m) => &m.loc_edit.edit,
            HV::Cds(c) => &c.loc_edit.edit,
            HV::Tx(t) => &t.loc_edit.edit,
            HV::Rna(r) => &r.loc_edit.edit,
            _ => return false,
        };
        edit.inner().is_some()
    }

    /// Re-derive the variant from the sequence it produces (#1229-#1235).
    ///
    /// The per-member pipeline normalizes each cis-allele member in isolation,
    /// which makes the canonical form depend on the input spelling and lets one
    /// member's 3'-shift run over a sibling. This pass runs once, at the top
    /// level, over the already-normalized variant: it applies the whole thing to
    /// the reference, re-derives the edit set from the (reference, result) pair,
    /// and partitions it globally. See `merge::canonicalize_from_sequence`.
    ///
    /// Canonicalizing can expose further work for the per-member pipeline (a
    /// derived member may itself 3'-shift), so the two alternate to a fixed
    /// point. The loop is bounded: each pass either reaches a fixed point or is
    /// abandoned in favour of the last stable value, so this can never diverge.
    fn canonicalize_from_sequence(
        &self,
        variant: HgvsVariant,
        warnings: Vec<NormalizationWarning>,
        manufactured: &mut Vec<ManufacturedJunctionExit>,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        /// Alternations between the sequence-first pass and `normalize_core`.
        /// Two is enough for every shape seen so far; the bound only exists so
        /// a pathological input cannot spin.
        const MAX_PASSES: usize = 4;

        let mut current = variant;
        // The value the loop started from, kept for the exhaustion fallback
        // below — that branch needs the pre-loop value rather than whichever
        // way-point the loop stopped on.
        //
        // Cloned lazily. Exhaustion requires the loop to have advanced at least
        // once, so a variant the pass declines on the first line never needs the
        // copy — and that is the overwhelming majority, since this method runs on
        // every `normalize()` and every projection while `sequence_first_pass`
        // serves only cis alleles and lone `delins`/`inv` members. Cloning
        // eagerly put a heap allocation on that whole path.
        let mut original: Option<HgvsVariant> = None;
        // Warnings belonging to the pass whose result is the one we return.
        //
        // Only the settled pass's warnings are kept, exactly as the sibling
        // loop in `normalize_allele` does with its own `settled_warnings`.
        // Accumulating across passes is wrong twice over: extending *before*
        // the `renormalized == current` check duplicates every warning the
        // re-derivation reproduces, and extending on a pass whose result a
        // later pass overwrites attaches warnings describing a discarded
        // intermediate state to the final variant. `assert_idempotent` cannot
        // catch either, because it compares `to_string()` of the variant and
        // never looks at the warning set.
        let mut settled_warnings: Vec<NormalizationWarning> = Vec::new();
        // Every `break` below is a *proof* that `current` is a fixed point of
        // this loop: the pass declined, or re-deriving reproduced `current`.
        // Running out of iterations proves nothing, so the two must be told
        // apart — see the exhaustion branch below.
        let mut converged = false;
        for _ in 0..MAX_PASSES {
            let Some(canonical) = self.sequence_first_pass(&current) else {
                converged = true;
                break;
            };
            if canonical == current {
                converged = true;
                break;
            }
            // The re-derived form is expressed in raw coordinates; hand it back
            // to the per-member pipeline so axis-specific canonicalization
            // (position shapes, boundary clamps, warnings) is applied to it.
            let (renormalized, pass_warnings) = self.normalize_core(&canonical, manufactured)?;
            if renormalized == current {
                converged = true;
                break;
            }
            original.get_or_insert_with(|| current.clone());
            settled_warnings = pass_warnings;
            current = renormalized;
        }

        // Exhausting the bound means the last pass's output was never shown to
        // be stable, and returning it breaks idempotency: `normalize` would
        // hand back pass 4's result, and normalizing *that* runs the same
        // alternation again and reaches pass 5. Since the loop is a pure
        // function of its input, returning the value it started from is
        // idempotent by construction — the second call re-enters with the same
        // value, exhausts identically, and returns it again.
        //
        // A missed canonicalization is the cost, and it is the right one: this
        // is the same trade every refusal in `merge` makes, and the pre-loop
        // value is the per-member pipeline's own output rather than an
        // arbitrary way-point. `assert_idempotent` is `#[cfg(debug_assertions)]`,
        // so nothing would catch the alternative in the builds that process
        // real data.
        //
        // Not reachable on any shape measured so far — two alternations settle
        // every case seen — but "not currently reachable" is exactly the
        // circumstance in which a silent wrong answer survives.
        if !converged {
            // `original` is always `Some` here: reaching the exhaustion branch
            // means every iteration ran, and each one sets it.
            return Ok((original.unwrap_or(current), warnings));
        }

        let mut warnings = warnings;
        warnings.extend(settled_warnings);
        Ok((current, warnings))
    }

    /// One sequence-first re-derivation. `None` when the variant is not a shape
    /// the pass serves (protein, trans alleles, fusions, and so on).
    fn sequence_first_pass(&self, variant: &HgvsVariant) -> Option<HgvsVariant> {
        // Only shapes where a partition decision actually exists. A lone
        // substitution / deletion / insertion / duplication already has exactly
        // one canonical partition, so re-deriving it can only lose what the
        // per-member pipeline added (gene symbol, RNA `U` alphabet, boundary and
        // exon-junction clamps). A lone `delins` or `inv` is different: those are
        // precisely the single-member forms that may need splitting (#1230,
        // #1232).
        //
        // Borrowed, not cloned: `canonicalize_from_sequence` refuses far more
        // often than it rebuilds, and every refusal used to cost a full copy of
        // the member list (a thousand members for a large allele) first.
        let members: &[HgvsVariant] = match variant {
            HV::Allele(allele) => {
                if allele.uncertain || allele.phase != AllelePhase::Cis || allele.variants.len() < 2
                {
                    return None;
                }
                &allele.variants
            }
            single if Self::is_splittable_single_member(single) => std::slice::from_ref(variant),
            _ => return None,
        };
        // When strict mode declares an input has no canonical form, the
        // derivation must not manufacture one.
        //
        // This is a deliberate narrowing of "the derivation is authoritative":
        // it is authoritative over *how* to spell a variant, not over *whether*
        // an ill-formed input has a spelling at all. An allele carrying an
        // overlap conflict is rejected by strict mode as `OverlapConflictingEdits
        // / W5002` precisely because the spec defines no canonical form for it;
        // re-deriving it from the sequence hands lenient mode a single tidy
        // member that strict mode then *accepts*, which launders the conflict out
        // of existence instead of leaving it visible to be rejected.
        //
        // The check has to be here, on `detect_overlap_conflicts`, because the
        // two refusals in play answer different questions in different coordinate
        // spaces and disagree for exactly this shape:
        //
        // * `merge::apply_edits_to_window` refuses on **interbase** geometry, so
        //   `24dup` (the zero-length junction `[24, 24]`) and `24C>G` (the span
        //   `[23, 24]`) are flush, not overlapping, and it admits the pair.
        // * `detect_overlap_conflicts` works in **HGVS-coordinate** space, where
        //   both members name base 24, and reports the conflict — and it is this
        //   one that strict mode and the #395 / #1276 / #1235 guards use.
        //
        // #1307 is the instance that showed the gap: `g.[24dup;24C>G]` on a
        // 24-base contig derived to `g.23_24insG`, which is in range and denotes
        // the same bases, yet turned a strict-rejected input into a
        // strict-accepted output. The lone-pure-insertion bail in
        // `canonicalize_from_sequence` had been masking it; removing that bail
        // made it reachable rather than creating it.
        //
        // Recomputed from the member list on every pass rather than read off a
        // carried `OverlapConflict` warning. That distinction is what keeps the
        // gate idempotent: the warning does not survive re-normalization, so
        // gating on it would refuse the first pass and admit the second.
        if members.len() > 1
            && crate::normalize::overlap::detect_overlap_conflicts(members, AllelePhase::Cis)
                .iter()
                .any(|w| matches!(w, NormalizationWarning::OverlapConflict { .. }))
        {
            return None;
        }
        // `phase` is `Cis` on both arms — the allele arm refuses anything else
        // above — and `uncertain` is likewise always false, which is exactly what
        // `AlleleVariant::new` builds.
        let rebuilt = merge::canonicalize_from_sequence(
            members,
            AllelePhase::Cis,
            &self.provider,
            self.config.shuffle_direction,
        )?;
        match rebuilt.len() {
            0 => None,
            1 => rebuilt.into_iter().next(),
            _ => Some(HV::Allele(crate::hgvs::variant::AlleleVariant::new(
                rebuilt,
                AllelePhase::Cis,
            ))),
        }
    }

    /// `normalize_core` followed by the sequence-first canonicalization.
    ///
    /// This is the seam every public exit must share. The pass cannot live in
    /// `normalize_core` itself — `normalize_allele` re-enters that per member, so
    /// it would recurse — but it also must not be wired into only one of the two
    /// exits: `normalize()` and `normalize_with_diagnostics()` would then return
    /// *different canonical strings* for exactly the shapes #1229-#1235 exist to
    /// canonicalize, and `normalize_with_diagnostics` is what the Python bindings
    /// and the web service call.
    ///
    /// An allele whose members overlap has no well-defined resulting sequence,
    /// so canonicalizing one would be canonicalizing a corrupted form. The
    /// sibling clamp (#1234) keeps the members the per-member pipeline *shifts*
    /// disjoint, but it does not cover a conflict that was in the input all
    /// along — an insertion interior to another member's span is reported
    /// (`OverlapConflict`) rather than resolved, and so survives into here.
    ///
    /// **This used to say `merge::apply_edits_to_window` refuses those on the
    /// edit geometry, "which is where the check has to live". It does not, and
    /// it is not.** That refusal is interbase-geometric, so it admits a pair that
    /// is flush in interbase space while coincident in HGVS-coordinate space —
    /// a `dup` and a substitution on one base, which is #1307. The gate that
    /// actually holds this line is in [`Self::sequence_first_pass`], on
    /// `detect_overlap_conflicts`; see the reasoning there, including why
    /// recomputing it per pass (rather than gating on the carried warning, the
    /// idempotency cost this note was warning about) keeps the pass a fixed
    /// point. Do not restore the older claim — leaving it in place is how the
    /// hole gets re-opened.
    fn normalize_core_canonical(
        &self,
        variant: &HgvsVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // The one place the per-leaf junction-exit provenance is collected. It is
        // created here rather than threaded in because this is the single seam
        // every public normalization exits through, and it is the only consumer.
        // `Vec::new()` does not allocate, so a normalization that manufactures
        // nothing — essentially all of them — pays nothing for carrying it.
        let mut manufactured: Vec<ManufacturedJunctionExit> = Vec::new();
        let (normalized, warnings) = self.normalize_core(variant, &mut manufactured)?;
        let (normalized, warnings) =
            self.canonicalize_from_sequence(normalized, warnings, &mut manufactured)?;
        let normalized = self.split_protein_separation(normalized);
        Ok((
            self.reparent_junction_exit(normalized, &manufactured),
            warnings,
        ))
    }

    /// Render a junction-crossing output against the genomic reference the
    /// crossing was computed against (#1704).
    ///
    /// # The defect
    ///
    /// `#670` (and its `normalize_tx`/`normalize_rna` mirrors) applies the 3'
    /// rule across an exon/intron junction: when the exon-confined shuffle lands
    /// exactly on an exon's last base, the shuffle is re-run in genomic space and
    /// the answer adopted **only if it crossed into the intron**. The coordinate
    /// that produces is correct — `background/numbering.md:26` explicitly
    /// withholds `general.md:44`'s exon/exon exception from an exon/intron
    /// junction, so the shift is required rather than merely permitted — but the
    /// *accession* it was rendered against was not. `checklist.md:20`: an `NM_`
    /// "can only be used to describe variants in introns using a `c.` prefix when
    /// a genomic reference sequence is given". So ferro manufactured 371
    /// descriptions (3' direction; 0 at 5', which has no mirror of the gate) that
    /// its own strict mode refuses as `IntronicOnBareTranscript` / W4007 — the
    /// reverse of the laundering `#1406` legislated against, and the residue the
    /// `bare-transcript-intronic-position` ruling recorded as unexcused.
    ///
    /// # Why re-parenting rather than declining, and rather than clamping
    ///
    /// **Not clamping.** Converging these back on the in-exon answer would
    /// implement the exception `numbering.md:26` withholds, and would silently
    /// revert `#670`, whose own fixtures pin the intronic coordinate as correct.
    /// The coordinate is not what is wrong.
    ///
    /// **Not declining.** The counter is a *lenient*-mode figure, so refusing in
    /// strict mode leaves every one of the 371 exactly where it is; and the input
    /// (`c.7del`) is a perfectly valid description, so refusing it because ferro's
    /// own chosen output form is inexpressible fails rule 1 in the other
    /// direction. It would also make the answer depend on the error mode, which
    /// rule 4 does not admit.
    ///
    /// **Re-parenting closes the class.** The genomic reference `checklist.md:20`
    /// asks for is not merely obtainable — the crossing could not have been
    /// computed without it, since the genomic re-shuffle fetches bases from a
    /// named contig and `#670` declines outright when `transcript.chromosome` is
    /// `None`. So the wrapper is already in hand at the moment the answer is
    /// adopted, and `NC_…(NM_…):c.10+2del` is the exact form the clause names.
    ///
    /// # Scope, and the two deliberate limits
    ///
    /// The trigger is *ferro manufactured the offset*, and since `#1723` that is
    /// a **per-leaf** fact carried from the `#670` gate in
    /// [`ManufacturedJunctionExit`] rather than a whole-description guess made
    /// here. An authored intronic position is left exactly as authored — that
    /// class is settled by the `bare-transcript-intronic-position` ruling (strict
    /// refuses, lenient accepts as written), and re-spelling it here would
    /// overturn a decided record as a side effect of fixing a different one.
    ///
    /// The `r.` axis is **not** covered, matching
    /// [`intronic_on_bare_transcript_warning`]'s scope exactly (`#486`/`#834`
    /// put it out of scope, and `checklist.md:20` names the `c.` prefix). One
    /// predicate governs input refusal and output rendering; widening only this
    /// half would put the two back out of step, which is the defect being fixed.
    ///
    /// Declining is still possible and is the safe direction: when no usable
    /// genomic accession can be built the output is returned untouched, which is
    /// the pre-`#1704` string rather than a worse one.
    ///
    /// # The one question this deliberately does NOT answer
    ///
    /// `ACC:c.[a;b]` renders compactly only when its members share an accession
    /// (`use_compact_form` / `all_share_accession_and_type`). So when one leaf on
    /// an accession needs a wrapper and a sibling on the same accession does not,
    /// there is a choice — **lift** the wrapper to the whole description, or
    /// **expand** to per-member accessions — and it is a representation-policy
    /// choice, not a code one.
    ///
    /// It is answered here only where it is not really a choice. Where every
    /// bare-intronic leaf on the accession is one ferro manufactured, lifting
    /// re-spells only *exonic* siblings, which [`reparent_leaves`] argues is the
    /// cheap side. Where a sibling is an **authored** intronic position, lifting
    /// would re-spell a leaf the `bare-transcript-intronic-position` ruling says
    /// to leave alone, and expanding would change the description's shape — so
    /// this declines, which is byte-for-byte what `#1704` shipped, and the
    /// question is on the record as the `undecided`
    /// `junction-exit-wrapper-scope-in-a-mixed-allele` ruling.
    ///
    /// The point of `#1723` is that both answers are now *expressible*: the
    /// decision reads a per-leaf classification and could take either branch.
    /// `#1704` could not represent the question at all — its input predicate was
    /// an any-leaf existence test, so one authored member silently vetoed the
    /// repair for every other member, including one ferro had moved itself.
    fn reparent_junction_exit(
        &self,
        output: HgvsVariant,
        manufactured: &[ManufacturedJunctionExit],
    ) -> HgvsVariant {
        // Cheapest question first, and the one that is false for essentially
        // every normalization: did any gate manufacture an offset at all?
        if manufactured.is_empty() {
            return output;
        }
        // Group the recorded manufactures by the accession they landed on, so a
        // description carrying two offending accessions repairs both. `#1704`
        // took the FIRST offending leaf and matched on accession equality with no
        // loop, so a second accession was never considered.
        //
        // **Two records for one bare accession must agree on the wrapper, and
        // disagreement DECLINES rather than picking one.** Keeping the first
        // silently would be the same shape as the `find_map` this pass removes,
        // one level in: the whole thesis here is that the contig must be the one
        // the crossing was computed against, so a second contig arriving for the
        // same transcript means the two crossings resolved differently and
        // neither is entitled to speak for the accession. `None` is the declined
        // state and is sticky — a third agreeing record does not revive it.
        //
        // **Not evidenced, deliberately stated anyway.** No provider in this tree
        // reaches it: `Transcript::chromosome` is constant per accession in every
        // fixture, and no input was constructed that makes two gates on one
        // accession resolve different contigs. So this is a latent invariant
        // written where it is relied on, not a fixed bug — the same status
        // `outer == *bare` had until LRG turned out to reach it.
        let mut by_accession: Vec<(&Accession, Option<&Accession>)> = Vec::new();
        for record in manufactured {
            match by_accession
                .iter_mut()
                .find(|(bare, _)| *bare == &record.bare)
            {
                Some((_, wrapper)) => {
                    if *wrapper != Some(&record.wrapper) {
                        *wrapper = None;
                    }
                }
                None => by_accession.push((&record.bare, Some(&record.wrapper))),
            }
        }

        let mut output = output;
        for (bare, wrapper) in by_accession {
            // The disagreement case above: two contigs for one transcript, so
            // this accession is left exactly as `#1704` left it.
            let Some(wrapper) = wrapper else {
                continue;
            };
            // Every bare-intronic leaf still standing on this accession, and
            // whether each is one a gate recorded. A leaf a later sibling pass
            // moved no longer matches any record, so it counts as unexplained and
            // the accession declines — the safe direction.
            let verdict = {
                let leaves: Vec<&HgvsVariant> = bare_transcript_intronic_leaves(&output)
                    .into_iter()
                    .filter(|leaf| leaf.accession() == Some(bare))
                    .collect();
                // An empty set is NOT "all manufactured". The gate recorded a
                // leaf, but nothing on this accession still names a bare intronic
                // position — a later sibling pass moved or dropped it — so there
                // is nothing to repair, and wrapping the exonic leaves that
                // remain would be a re-spelling no clause asks for.
                !leaves.is_empty()
                    && leaves
                        .iter()
                        .all(|leaf| manufactured.iter().any(|record| &record.leaf == *leaf))
            };
            if !verdict {
                // Either nothing to repair, or the mixed case. Declining
                // reproduces `#1704` exactly; see the `undecided` ruling named
                // above for what the mixed half leaves open.
                continue;
            }
            output = reparent_leaves(output, bare, wrapper);
        }
        output
    }

    /// The protein-axis partition pass, run once at the top level.
    ///
    /// [`Self::try_protein_split_delins`] holds the rule and its citations; this
    /// is only about **where** it runs. It runs here, beside
    /// [`Self::canonicalize_from_sequence`] and for the same reason: the move
    /// turns one member into several, so it belongs at the seam that sees the
    /// whole description rather than in
    /// [`normalize_protein`](Self::normalize_protein), which `normalize_allele`
    /// re-enters once per member.
    ///
    /// Wiring it into `normalize_protein` was tried, and produced two defects
    /// that every bare-description test passes straight over:
    ///
    /// ```text
    ///   p.[Ser44_Trp46delinsArgLeuArg;Ala60Gly]
    ///     -> [NP_003997.1:p.Ala60Gly;NP_003997.1:p.[Ser44Arg;Trp46Arg]]
    ///   p.[Ser44_Trp46delinsArgLeuArg];[Ala60Gly]
    ///     -> NP_003997.1:p.[Ser44Arg;Trp46Arg];[Ala60Gly]     (does not re-parse)
    /// ```
    ///
    /// The first is a nested bracket carrying a repeated accession, which is not
    /// a description at all. The second *is* legal HGVS — the DNA axis emits
    /// `g.[A;B];[C]` and reads it back — but ferro's **protein** allele grammar
    /// accepts only single-member arms, so normalization would emit a string
    /// `parse_hgvs` refuses. That is exactly what `FERRO_ASSERT_REPARSE` exists to
    /// catch, and its exemption list is closed on purpose.
    ///
    /// So the move applies to a whole `p.` description and **declines for a
    /// member of one**: `p.[<separated delins>;…]` keeps its delins. Closing that
    /// needs two decisions this pass deliberately does not take — where an inner
    /// allele's `( )` marker goes once its members flatten into the outer bracket
    /// (`protein/alleles.md:34` puts it per member, `:90` puts it round the group)
    /// and a protein allele grammar that reads `p.[A;B];[C]`. Pinned by
    /// `protein_axis_split_move::a_delins_inside_a_cis_allele_is_left_alone` and
    /// its trans sibling.
    fn split_protein_separation(&self, variant: HgvsVariant) -> HgvsVariant {
        let HV::Protein(protein) = &variant else {
            return variant;
        };
        match self.try_protein_split_delins(protein) {
            Some(members) => {
                let uncertain = protein.loc_edit.edit.is_uncertain();
                wrap_allele_if_split(members.into_iter().map(HV::Protein).collect(), uncertain)
            }
            None => variant,
        }
    }

    /// Normalize a variant, apply the strict-mode rejection ladder, and return
    /// the normalized variant together with any warnings that were NOT promoted
    /// to hard errors.
    ///
    /// This is the shared core behind both `normalize()` (which discards the
    /// warnings) and `VariantProjector`, which surfaces them on each projection.
    /// Routing the projector through here — rather than the raw `normalize_core`
    /// — is what keeps a strict-configured projector rejecting the same inputs
    /// (`EINTRONIC`, `RefSeqMismatch`, W5002/W5003/W5004/W4004/W4005/W4006, …) it
    /// would reject via `normalize()`. In the default lenient config every
    /// `should_reject_*` is false, so the ladder is a no-op and the warnings pass
    /// straight through.
    pub(crate) fn normalize_core_checked(
        &self,
        variant: &HgvsVariant,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // `standards.md:39` (#1627): a member stating an alignment-only symbol
        // cannot be normalized, so this refuses BEFORE any work — and, unlike
        // every rung of the ladder below, in EVERY mode.
        //
        // The mode gate belongs at parse, not here. Per the decided
        // `rulings[absolute-prohibition-enforcement-stage]`, it governs whether
        // the INPUT is judged; it "does not, and cannot, govern whether the
        // output conforms", because rule 1 of the README ruleset — "Output
        // follows the HGVS recommendations. Absolute — never traded." — is about
        // output. `X` is a masked nucleotide: the base it stands for was not
        // resolved, so there is no sequence to shuffle, no denotation to derive,
        // and nothing a lenient mode could be lenient *toward*. Lenient
        // therefore fails here on exactly the ground the ruling gives it — it
        // cannot normalize — rather than on an input-conformance check it does
        // not run.
        //
        // Before this, all three modes emitted the offending member back
        // verbatim with an EMPTY warning vector: normalization was not
        // impossible, it was VACUOUS, which is worse, because the output looks
        // normalized while carrying a spelling the recommendations prohibit.
        if let Some(found) = crate::hgvs::alignment_symbols::alignment_only_symbol(variant) {
            // `InvalidCoordinates` is this function's existing carrier for a
            // normalize-stage conformance refusal (W5004 and W3022 below both
            // use it); the fault is in the stated bases rather than a position,
            // and the `[W3028]` tag in the message is what names it.
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "[{}] cannot normalize `{}`: it states `{}` ({}), which {} lists as used \
                     in alignment only. A masked base names no nucleotide, so the description \
                     denotes no sequence and there is nothing to normalize. State the resolved \
                     bases instead (`{}` is the IUPAC symbol for an unknown base).",
                    crate::error_handling::ErrorType::AlignmentOnlySymbolInDescription.code(),
                    found.stated,
                    found.symbol,
                    found.meaning(),
                    found.clause(),
                    found.unknown_base(),
                ),
            });
        }

        // `background/numbering.md:6`/`:8`/`:11` (#1628): a `g.`/`o.`/`m.`
        // position may not carry a `+`/`-` offset, so a description that states
        // one cannot be normalized — and, like the rung above and unlike every
        // rung of the ladder below, this refuses in EVERY mode.
        //
        // The mode gate belongs at parse (`apply_genomic_offset_rule`), for the
        // reason `rulings[absolute-prohibition-enforcement-stage]` gives: it
        // governs whether the INPUT is judged and "does not, and cannot, govern
        // whether the output conforms". A genomic accession carries no exon
        // table, so the offset is measured from nothing and the position names
        // no nucleotide; there is no sequence to shuffle and nothing a lenient
        // mode could be lenient toward. Lenient therefore fails here on exactly
        // the ground the ruling gives it — it cannot normalize — rather than on
        // an input-conformance check it does not run.
        //
        // Before this, all three modes returned the offending description
        // BYTE-IDENTICALLY with an EMPTY warning vector, from the
        // offset-carrying short circuit in `normalize_genome_variant`:
        // normalization was not impossible, it was VACUOUS. `hgvs_to_spdi`
        // meanwhile DROPPED the offset and answered for `g.266del`, so two
        // halves of ferro called one description two different variants — the
        // confluence half of #1628, closed separately by #1641 and #1734.
        if let Some(found) = crate::hgvs::genomic_offsets::genomic_axis_offset(variant) {
            // `InvalidCoordinates` is this function's existing carrier for a
            // normalize-stage conformance refusal (W5004, W3022 and W3028 all
            // use it); here the fault genuinely is in a position, and the
            // `[W4009]` tag in the message is what names it.
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "[{}] cannot normalize `{}`: {}.",
                    crate::error_handling::ErrorType::GenomicPositionOffset.code(),
                    variant,
                    found.refusal(),
                ),
            });
        }

        // Call the canonical core directly (skipping the per-call
        // `detect_shuffle_infos` work `normalize_with_diagnostics` does) and wrap
        // with empty infos; the ladder below only inspects warnings.
        let (normalized, warnings) = self.normalize_core_canonical(variant)?;
        let result = NormalizeResult::with_warnings(normalized, warnings);

        // EINTRONIC (#486): reject an intronic offset on a bare transcript
        // reference. The reference *form* being invalid (the spec marks
        // `NM_:c.357+1` "not correct") is more fundamental than any
        // offset-magnitude error, so this is checked BEFORE PositionPastEnd —
        // a bare intronic-and-past-intron-end input rejects as EINTRONIC, not
        // W4004. Reuses `FerroError::IntronicVariant` because only that variant
        // carries the EINTRONIC tag in mutalyzer_map.rs (the runner
        // substring-matches the `IntronicVariant` variant name). This is a
        // spec-form rejection and is distinct from the reduced-capability path
        // in `normalize_intronic_cds`/`_tx` (no genomic data), which since
        // #1012 no longer errors — it warn-and-degrades with
        // `ReducedCapabilityNoGenome`.
        if self.config.should_reject_intronic_bare_transcript() {
            if let Some(err) = result.warnings.iter().find_map(|w| match w {
                NormalizationWarning::IntronicOnBareTranscript {
                    variant,
                    coordinate_system,
                } => Some(FerroError::IntronicVariant {
                    // Keep `variant` the raw HGVS string (machine-readable).
                    variant: variant.clone(),
                    // Carry the distinguishing clarifier in `detail` instead, so
                    // this spec-form rejection is not confused (in logs) with the
                    // capability-failure `IntronicVariant` from
                    // `normalize_intronic_cds` ("no genomic data") — without
                    // polluting the `variant` field. The EINTRONIC tag in
                    // mutalyzer_map.rs keys off the variant *name*, not this
                    // string, so the clarifier does not affect conformance
                    // mapping.
                    detail: Some({
                        // Parent-reference example matching the input axis: a
                        // coding `NM_` parent for c., a non-coding `NR_` for n.
                        let parent_example = if coordinate_system == "n" {
                            "NG_(NR_)/NC_(NR_)"
                        } else {
                            "NG_(NM_)/NC_(NM_)"
                        };
                        format!(
                            "intronic offset on a bare {coordinate_system}. transcript \
                             reference; a genomic reference is required, e.g. {parent_example} — \
                             IntronicOnBareTranscript / W4007 / EINTRONIC"
                        )
                    }),
                }),
                _ => None,
            }) {
                return Err(err);
            }
        }

        // Strict mode also rejects W5004 IncompleteCdsStartReference: a
        // `c.`/`p.`/`r.` variant described against a transcript whose 5' CDS
        // is annotated incomplete (`cds_start_NF`) has no confirmed ATG, so
        // `c.1`/`p.1`/`r.1` are undefined against it — not an
        // HGVS-recommended reference for that axis. `normalize_cds` /
        // `normalize_protein` / `normalize_rna` already decline to
        // re-number the coordinate (verbatim pass-through) in every mode;
        // this only promotes the warning to a hard reject in strict mode.
        // Input-side counterpart to Task 4's projection-side decline (#972).
        if self.config.should_reject_incomplete_cds_start() {
            if let Some(err) = result.warnings.iter().find_map(|w| match w {
                NormalizationWarning::IncompleteCdsStartReference {
                    accession,
                    coordinate_system,
                } => Some(FerroError::InvalidCoordinates {
                    msg: format!(
                        "{accession}: transcript has an incomplete (unconfirmed ATG) CDS start \
                         (cds_start_NF); not an HGVS-recommended reference for \
                         {coordinate_system}. description — use the genomic (g.) or non-coding \
                         (n.) representation instead (IncompleteCdsStartReference / W5004)"
                    ),
                }),
                _ => None,
            }) {
                return Err(err);
            }
        }

        // W3022 InitiatorMetCanonicalization is an advisory about ferro's own
        // canonical output, so — unlike every other rung of this ladder — the
        // base mode does not promote it; only an explicit `--reject W3022`
        // does. See `NormalizeConfig::should_reject_initiator_met_canonicalization`
        // for why, and #1196 for why the `--reject` direction has to do
        // *something* rather than be silently inert.
        if self.config.should_reject_initiator_met_canonicalization() {
            if let Some(err) = result.warnings.iter().find_map(|w| match w {
                NormalizationWarning::InitiatorMetCanonicalization {
                    accession,
                    location,
                } => Some(FerroError::InvalidCoordinates {
                    msg: format!(
                        "{accession}: canonical form `p.{location}dup` includes the initiator \
                         methionine, and this code was explicitly rejected \
                         (InitiatorMetCanonicalization / W3022). The duplication is \
                         HGVS-canonical; drop `--reject W3022` to accept it, or describe the \
                         predicted consequence as `p.0?` / `p.(Met1?)` instead"
                    ),
                }),
                _ => None,
            }) {
                return Err(err);
            }
        }

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
                            "{accession}: '{marker}' is an assembly-annotated region with no \
                             sequence-derivable coordinate and cannot be normalized \
                             (UnresolvableCentromere / W4005)"
                        ),
                    })
                }
                _ => None,
            }) {
                return Err(err);
            }
        }

        // In strict mode, reject a telomere marker that resolves to a
        // transcript-flank position on a genomic-reference c. (W4006). The
        // flank is not numberable in c. per HGVS numbering.md; #488 Phase 2b.
        if self.config.should_reject_transcript_flank() {
            if let Some(err) = result.warnings.iter().find_map(|w| match w {
                NormalizationWarning::TranscriptFlankNotDescribable { accession, marker } => {
                    Some(FerroError::InvalidCoordinates {
                        msg: format!(
                            "{accession}: '{marker}' on a genomic-reference c. denotes a 5'/3' \
                             transcript-flank position, not numberable in c. (HGVS numbering.md \
                             transcript-flanking); use the genomic g. form \
                             (TranscriptFlankNotDescribable / W4006)"
                        ),
                    })
                }
                _ => None,
            }) {
                return Err(err);
            }
        }

        // In strict mode, reject a reduced-capability (no-genomic-data)
        // degradation. Lenient/silent return the best-effort variant with the
        // `ReducedCapabilityNoGenome` advisory, but strict mode must not hand
        // back a knowingly-degraded result, so promote it to a typed error —
        // the same reject-in-strict contract the intronic/boundary genomic
        // paths carried before #1012 (now unified across every genome-requiring
        // step). #1012 item 2.
        //
        // Ordering: this runs LAST, after every spec-validity rejection above.
        // A variant that is both genuinely invalid (e.g. past the CDS end,
        // W4004) *and* degraded must report the input defect first — the
        // reduced-capability limit is environmental, so it is the lowest-
        // priority strict rejection.
        if self.config.should_reject_reduced_capability() {
            if let Some(err) = result.warnings.iter().find_map(|w| match w {
                NormalizationWarning::ReducedCapabilityNoGenome {
                    variant,
                    capability,
                } => Some(FerroError::ReducedReferenceCapability {
                    variant: variant.clone(),
                    capability: capability.clone(),
                }),
                _ => None,
            }) {
                return Err(err);
            }
        }

        self.assert_seam_oracles(variant, &result.result, &result.warnings);

        // Provenance (#1235, `DNA/delins.md:79-84`). If the input described more
        // cis members than the normalized form does, two or more separately
        // reported variants have been coalesced, and that fact is not
        // recoverable from the output string — by design, since making the
        // string depend on the input's spelling is precisely the
        // non-confluence #1235 removes.
        //
        // Emitted here rather than at any single coalescing site on purpose:
        // several passes can reduce the member count (the live rule's adjacency
        // coalesce, the prioritisation collapse, the sequence-first
        // re-derivation), and the caller cares that it happened, not which pass
        // did it. Comparing the counts at the one exit catches all of them and
        // cannot drift as those passes change.
        //
        // It never rejects and never alters `result.result`.
        let mut result = result;
        let reported_members = cis_member_count(variant);
        let normalized_members = cis_member_count(&result.result);
        if normalized_members < reported_members {
            if let Some((accession, coordinate_system)) = accession_and_axis(variant) {
                result
                    .warnings
                    .push(NormalizationWarning::MembersCoalesced {
                        accession,
                        coordinate_system,
                        reported_members,
                        normalized_members,
                    });
            }
        }

        Ok((result.result, result.warnings))
    }

    /// The four seam oracles, run on a normalized result before it is returned.
    ///
    /// Ordered cheapest-and-most-specific first, so the failure a run reports is
    /// the most precise one available: a bad coordinate is a bad coordinate
    /// whether or not the string also re-parses or re-normalizes, and
    /// `assert_idempotent` would otherwise report an out-of-range output as a
    /// non-idempotency, which is the symptom rather than the fault. See
    /// `assert_idempotent` for the re-entrancy guard.
    ///
    /// `assert_denoted_sequence` (#1615) is last for the same reason: it is the
    /// only one that reads the reference, and an output that is out of bounds or
    /// unparseable should be reported as *that* rather than as a sequence it
    /// could not be applied to. It is also the most expensive.
    ///
    /// A method rather than an inline block because there is more than one public
    /// exit to cover. `normalize_core_checked` is the seam for `normalize()` and
    /// every `VariantProjector` path, but `normalize_with_diagnostics` reaches
    /// `normalize_core_canonical` directly and so returned a normalized variant
    /// that no oracle had ever inspected — for all three checks, not just the new
    /// one. Calling this from both is what makes "the oracles cover every
    /// normalized output ferro hands back" true rather than nearly true.
    ///
    /// Each call stays individually `#[cfg(debug_assertions)]`-gated, so the whole
    /// body compiles away in release exactly as the inline block did.
    fn assert_seam_oracles(
        &self,
        variant: &HgvsVariant,
        normalized: &HgvsVariant,
        warnings: &[NormalizationWarning],
    ) {
        #[cfg(debug_assertions)]
        self.assert_in_bounds(variant, normalized);
        #[cfg(debug_assertions)]
        self.assert_reparseable(variant, normalized);
        #[cfg(debug_assertions)]
        self.assert_idempotent(variant, normalized);
        #[cfg(debug_assertions)]
        self.assert_denoted_sequence(variant, normalized, warnings);
        // Release builds read none of the parameters once the four calls above
        // are compiled out.
        #[cfg(not(debug_assertions))]
        let _ = (variant, normalized, warnings);
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
        // Through `normalize_core_checked`, not `normalize_core_canonical`: this
        // is a *public* exit, so it owes the caller the same strict-mode contract
        // `normalize()` does (#1382). It used to call the canonical core directly
        // and so applied none of the `should_reject_*` ladder, which meant a
        // strict-configured `Normalizer` rejected a variant through `normalize()`
        // and accepted the very same one through here. The asymmetry was invisible
        // in the default lenient config, where every rung is a no-op.
        //
        // `normalize_core_checked` routes through `normalize_core_canonical`
        // itself, so the variant this returns is the same canonical one as before
        // and `infos` still describes the shift from the input to *that* variant —
        // which is why the ladder can be added without moving the output.
        //
        // The warnings it returns are the ones NOT promoted to hard errors, so a
        // lenient caller still sees everything this entry point exists to report.
        //
        // This also subsumes #1366's separate `assert_seam_oracles` call here.
        // That call existed precisely because this exit bypassed
        // `normalize_core_checked`; now that it does not, the oracles run once at
        // that method's own exit. Calling them here as well would re-run all three
        // on every diagnostics normalization — including `assert_idempotent`'s full
        // re-normalization, the most expensive of them.
        let (result, warnings) = self.normalize_core_checked(variant)?;
        let infos = detect_shuffle_infos(variant, &result, self.config.shuffle_direction);
        Ok(NormalizeResult::with_diagnostics(result, warnings, infos))
    }

    /// The requested transcript accession a transcript-coordinate variant would
    /// be SILENTLY SUBSTITUTED for, if any — the #785 trigger.
    ///
    /// Returns `Some(requested_id)` only when ALL of the following hold:
    /// - the variant is a transcript-coordinate axis (`c.`/`n.`/`r.`);
    /// - the caller named an EXPLICIT version (bare references intentionally
    ///   keep lenient "latest version" resolution and are not gated; LRG
    ///   accessions never carry a version, so they are likewise unaffected);
    /// - the named accession is a transcript ([`is_transcript_accession`], so a
    ///   genomic/gene wrapper like an unrewritten `NG_(GENE):c.…` selector is
    ///   not mistaken for a versioned transcript);
    /// - the provider would NOT serve the bases of that exact version
    ///   ([`ReferenceProvider::has_transcript_version_exact`] is false); and
    /// - the provider nonetheless resolves the request — i.e. a sibling version
    ///   is actually substituted.
    ///
    /// The version-exact predicate is the load-bearing signal: it answers "do
    /// the bases the read path would serve correspond to the exact requested
    /// version?" and so catches BOTH substitution paths — a FASTA version-strip
    /// fallback (which serves the sibling's accession) AND a cdot-genome
    /// reconstruction off a sibling record (which serves the *requested*
    /// accession with the sibling's frame). A naive "served accession differs"
    /// check would miss the latter. The `get_transcript(...).is_ok()` guard then
    /// distinguishes a genuine substitution (a sibling is served) from a
    /// transcript that is simply absent with no sibling to fall back to — the
    /// latter is left to the existing missing-transcript handling (echoed
    /// unchanged) rather than gated. The genomic-context wrapper is stripped, so
    /// a compound `NC_…(NM_….3):c.…` is gated on its inner `NM_….3`.
    fn silently_substituted_transcript(&self, variant: &HV) -> Option<String> {
        let accession = match variant {
            HV::Cds(v) => &v.accession,
            HV::Tx(v) => &v.accession,
            HV::Rna(v) => &v.accession,
            _ => return None,
        };
        accession.version?; // explicit version only
        let requested = accession.transcript_accession();
        if !is_transcript_accession(&requested) {
            return None;
        }
        if self.provider.has_transcript_version_exact(&requested) {
            return None; // the exact version's bases are served — no substitution
        }
        // Not version-exact: a sibling would be served (substitution) — gate it —
        // unless the transcript is simply absent with no sibling (Err), which is
        // left to the existing missing-transcript path.
        self.provider
            .get_transcript(&requested)
            .is_ok()
            .then_some(requested)
    }

    /// Core normalization: dispatch by variant kind and return the
    /// normalized variant plus warnings, WITHOUT computing the diagnostic
    /// `infos` axis. `normalize()` discards infos, so it calls this
    /// directly to skip the per-call `detect_shuffle_infos` cost on the
    /// hot path; `normalize_with_diagnostics` layers infos on top. The
    /// (variant, warnings) output is identical to the previous inline match.
    pub(crate) fn normalize_core(
        &self,
        variant: &HgvsVariant,
        manufactured: &mut Vec<ManufacturedJunctionExit>,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // Rewrite a legacy gene-model selector (`NG_(GENE_v001):c.…`) to the
        // spec-preferred transcript form (`NG_(NM_):c.…`) up front, so the rest
        // of normalization — and projection, which normalizes first (#637) —
        // operates on the resolved transcript. Unresolvable selectors are left
        // unchanged (#500).
        let rewritten = self.rewrite_legacy_gene_selector(variant);
        let variant = rewritten.as_ref().unwrap_or(variant);

        // Silent version-substitution gate (#785). A c./n./r. description is
        // defined only against its named `accession.version`: per HGVS, `c.1` is
        // the first base of *that* version's start codon and a versioned
        // identifier is required precisely because versions differ in length/CDS
        // offset (background/numbering.md, background/refseq.md). When the caller
        // names an EXPLICIT transcript version the reference does not carry, a
        // lenient provider silently falls back to a *sibling* version and serves
        // its frame/sequence — yielding a spec-invalid result whose stated
        // reference and coordinate frame disagree, with no error. Decline with
        // `TranscriptVersionNotExact` (a clean reference-unavailable miss the
        // conformance mapping and `ferro project` CLI already treat as soft)
        // rather than substitute. Projection inherits this because it normalizes
        // first (#637), and `normalize_allele` recurses through `normalize_core`
        // so each allele member is checked.
        if let Some(requested) = self.silently_substituted_transcript(variant) {
            return Err(FerroError::TranscriptVersionNotExact { requested });
        }

        // EINTRONIC (#486): an intronic offset on a bare transcript reference
        // (NM_ c. / NR_ n., no NG_(…)/NC_(…) context) is a spec-invalid
        // description form. Detect it up front — before any per-axis
        // short-circuit, so substitutions are covered, and once per allele
        // member (normalize_allele recurses through normalize_core).
        let eintronic_warning = if self.config.should_reject_intronic_bare_transcript()
            || self.config.should_warn_intronic_bare_transcript()
        {
            intronic_on_bare_transcript_warning(variant)
        } else {
            None
        };
        // In reject mode, short-circuit BEFORE normalization: the form is
        // invalid regardless of whether the intronic position would resolve,
        // and running normalize_cds/normalize_tx first can surface an
        // incidental ConversionError (or capability-failure IntronicVariant)
        // instead of the EINTRONIC form error. The echoed result is discarded
        // when `normalize()` escalates the warning to FerroError::IntronicVariant.
        if self.config.should_reject_intronic_bare_transcript() {
            if let Some(w) = &eintronic_warning {
                return Ok((variant.clone(), vec![w.clone()]));
            }
        }

        // Dispatch normalization. On success, warn-only mode attaches W4007 to the
        // fully normalized output below (3'-shift, W4004, repeat notation are all
        // preserved). On failure, warn-only mode must NOT propagate the error for
        // an intronic-on-bare-transcript form: the spec-invalidity of that form is
        // independent of provider capability, so when the intronic position cannot
        // resolve (e.g. an indel with no genomic sequence to anchor it) we recover
        // by echoing the input with W4007 instead of dropping the warning and
        // surfacing a capability IntronicVariant/ConversionError (#682). A
        // substitution never reaches the resolution pass, so this only changes the
        // resolve-failure path. Reject mode already short-circuited above; silent
        // mode (no eintronic_warning) propagates the error unchanged.
        let (result, mut warnings) = match self.normalize_dispatch(variant, manufactured) {
            Ok(resolved) => resolved,
            Err(e) => match &eintronic_warning {
                Some(w) => return Ok((variant.clone(), vec![w.clone()])),
                None => return Err(e),
            },
        };

        // Warn-only: attach the EINTRONIC warning to the resolved output. (Reject
        // mode already short-circuited above; silent mode has no warning.)
        if let Some(w) = eintronic_warning {
            warnings.push(w);
        }

        // A cis allele's members are merged *before* any of them reaches the
        // per-member reference validator, and a merge keeps only their alt
        // bases — so a member consumed by one has its reference assertions
        // dropped without ever being checked (#1543). Ask the question here, on
        // the authored input, where the assertions still exist.
        //
        // Deliberately outside `normalize_allele` rather than inside it: that
        // function has several early returns, and the point of this guard is
        // that no rearrangement of the pipeline below can make it vacuous. It
        // adds warnings and changes nothing else, so the merged form a
        // correctly-spelled allele produces is untouched.
        if let HV::Allele(allele) = variant {
            for warning in merge::authored_member_reference_mismatches(
                &allele.variants,
                allele.phase,
                &self.provider,
            ) {
                // The per-member pipeline reports the same finding for any
                // member a merge did *not* consume, and it gets there first.
                // Position is the whole key: two `RefSeqMismatch` warnings over
                // one sequence-axis span are one mismatch reported twice.
                let NormalizationWarning::RefSeqMismatch { position, .. } = &warning else {
                    continue;
                };
                let already = warnings.iter().any(|w| {
                    matches!(
                        w,
                        NormalizationWarning::RefSeqMismatch { position: seen, .. }
                            if seen == position
                    )
                });
                if !already {
                    warnings.push(warning);
                }
            }
        }

        Ok((result, warnings))
    }

    /// Dispatch normalization by variant kind, returning the normalized variant
    /// plus warnings. Extracted from [`Self::normalize_core`] so the caller can
    /// recover a warn-only EINTRONIC resolve-failure (#682) without the per-arm
    /// `?` propagating the error and discarding the W4007 warning; this method
    /// itself performs no EINTRONIC handling.
    fn normalize_dispatch(
        &self,
        variant: &HgvsVariant,
        manufactured: &mut Vec<ManufacturedJunctionExit>,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        Ok(match variant {
            // A `g.` variant on a mitochondrial reference (NC_012920 / NC_001807)
            // must be described with `m.` (#487). Coerce it to an `MtVariant`
            // — same accession/location/edit, only the coordinate system label
            // differs — and normalize on the mito path. Parsing is unchanged
            // (still yields `Genome`); the coercion happens only here.
            HV::Genome(v) if v.accession.is_mitochondrial() => {
                let mt = crate::hgvs::variant::MtVariant {
                    accession: v.accession.clone(),
                    gene_symbol: v.gene_symbol.clone(),
                    loc_edit: v.loc_edit.clone(),
                };
                self.normalize_mt(&mt)?
            }
            HV::Genome(v) => self.normalize_genome(v)?,
            HV::Cds(v) => self.normalize_cds(v, manufactured)?,
            HV::Tx(v) => self.normalize_tx(v, manufactured)?,
            HV::Protein(v) => self.normalize_protein(v)?,
            HV::Rna(v) => self.normalize_rna(v, manufactured)?,
            HV::Mt(v) => self.normalize_mt(v)?,
            HV::Allele(a) => self.normalize_allele(a, manufactured)?,
            // Circular (`o.`, SVD-WG006) variants are returned unchanged: no
            // 3'-shift runs and no warning is emitted. A genuine circular
            // normalizer would 3'-shift with origin-wraparound semantics (cf.
            // the mitochondrial wraparound handling in `normalize_mt`), which is
            // not yet implemented.
            //
            // Tracked by #466's circular candidate, which asks SVD-WG to state
            // the origin-wraparound shift semantics first — "the most 3'
            // position possible" is ill-defined on a ring, so there is nothing
            // to implement against yet.
            //
            // NOT #951, and not #129 before it: both are closed. This pointer
            // has now gone stale twice (#129 -> #951 -> here), and #951's own
            // body was filed *because* the comment then cited the closed #129.
            // If #466's candidate is ever closed too, repoint this rather than
            // leaving the third generation of the same dangling reference.
            // (The earlier "normalize like genomic variants" note overstated
            // this pass-through clone.)
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
            // Genome rings pass through unchanged (not normalized; #546)
            HV::GenomeRing(g) => (HV::GenomeRing(g.clone()), vec![]),
            // Supernumerary markers pass through unchanged (not normalized; #546)
            HV::Supernumerary(inner) => (HV::Supernumerary(inner.clone()), vec![]),
            // Null and unknown allele markers pass through unchanged
            HV::NullAllele => (HV::NullAllele, vec![]),
            HV::UnknownAllele => (HV::UnknownAllele, vec![]),
        })
    }

    /// Normalize an allele (compound) variant
    ///
    /// Normalizes each variant in the allele individually, with overlap prevention.
    /// After normalization, checks if variants would overlap and constrains shifting
    /// to prevent collisions.
    fn normalize_allele(
        &self,
        allele: &crate::hgvs::variant::AlleleVariant,
        manufactured: &mut Vec<ManufacturedJunctionExit>,
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
        // Detect insertion overlaps on the *raw* members before merge collapses
        // them: two insertions sharing a junction, or an insertion interior to a
        // span edit, are mutalyzer `EOVERLAP` cases (#486). The merge below
        // would otherwise fold them into one combined edit, hiding the overlap
        // from the post-shift `detect_overlap_conflicts` pass. Strict mode then
        // promotes the W5002 warning to a typed error.
        all_warnings.extend(crate::normalize::overlap::detect_insertion_overlaps(
            &allele.variants,
            allele.phase,
        ));

        // Coincident *bounds* among the raw members, for exactly the reason the
        // insertion detector above runs pre-merge — and the half that was
        // missing (#1406).
        //
        // `detect_overlap_conflicts` also runs post-shift, on the normalized
        // members, which meant W5002 reported whether the repair had *succeeded*
        // rather than whether the input conflicted. Two inputs of the same shape
        // got opposite verdicts:
        //
        // ```text
        // reference  TTTTTTTTTAATATATTTTAATAC     24 bases
        // g.[23dup;23A>G]  -> g.22_23insG          no warning, strict ACCEPTS
        // g.[24dup;24C>G]  -> unchanged            W5002, strict REJECTS
        // ```
        //
        // Both are a `dup` of one base beside a substitution of that same base.
        // The only difference is that at 23 the repair collapses the pair into a
        // single member, so by the time the post-shift detector looks there is
        // nothing left to conflict; at 24 the terminal decline (#1307) leaves
        // the pair intact and it is reported. Proximity to the end of the contig
        // decided whether a conflicting description was rejected.
        //
        // The same mechanism let lenient mode launder a conflict outright:
        // `g.[11_12inv;11_12insAA]` is strict-rejected, and its own lenient
        // output `g.[11_12=;10_11A[4]]` is strict-*accepted*.
        //
        // So the verdict is taken from the input, where the conflict lives, and
        // a conflicting allele is preserved as authored — the same answer
        // `has_same_gap_insertions` gives just below, and ferro's standing one
        // since #395/#486/#1004: report the conflict and leave the description
        // alone rather than pick a winner among orderings the spec does not
        // rank.
        // Scoped to *coincident bounds*, and deliberately not widened to the
        // insertion-interior-to-a-span conflicts `detect_insertion_overlaps`
        // reports. Those alleles are still expected to have each member
        // normalized on its own — `[4_10inv;5_6insAA]` settles as
        // `[5_9inv;6_7insAA]`, the inversion 3'-shifted and the insertion
        // respelled at its own junction — and preserving them verbatim would
        // drop the per-member 3' rule that #1276 and #1235's transcript axes
        // pin. What those alleles do not get is whole-allele re-derivation,
        // which they already did not.
        let raw_conflicts =
            crate::normalize::overlap::detect_overlap_conflicts(&allele.variants, allele.phase);
        if !raw_conflicts.is_empty() {
            // Harvest the per-member diagnostics before preserving. Returning
            // here skips the per-member loop below, which is where every
            // *member-local* finding is raised — so without this a conflicting
            // allele that ALSO misstates a reference base reports only the
            // conflict. Measured on `g.[13G>A;13A>C]` (base 13 is `A`): the
            // warning set went from `[RefSeqMismatch, OverlapConflict]` to
            // `[OverlapConflict]`, and strict mode's error changed from
            // "Reference mismatch at 13-13: expected G" to the coincidence
            // message. A reference mismatch is a data-integrity signal about the
            // caller's input, independent of whether the members also collide;
            // losing it to a sibling's conflict is a strictly worse report.
            //
            // Each member is normalized in isolation purely to collect what it
            // says — the results are discarded, because preserving the authored
            // members is the whole point of this branch.
            //
            // Errors propagate rather than being swallowed. An earlier cut used
            // `if let Ok(..)`, which quietly made this branch *more* permissive
            // than the ordinary per-member path: a member that fails hard —
            // `TranscriptVersionNotExact`, say — surfaces to the caller there
            // and would have been discarded here, so an allele would go from
            // erroring to returning a preserved value merely because a sibling
            // conflicted with it. Whether the members collide is unrelated to
            // whether one of them is resolvable at all.
            for member in &allele.variants {
                let (_, member_warnings) = self.normalize_core(member, manufactured)?;
                all_warnings.extend(member_warnings);
            }
            all_warnings.extend(raw_conflicts);
            let mut preserved =
                crate::hgvs::variant::AlleleVariant::new(allele.variants.clone(), allele.phase);
            preserved.uncertain = allele.uncertain;
            return Ok((HgvsVariant::Allele(preserved), all_warnings));
        }

        // Two separate insertions at the *same* junction (`g.[4_5insT;4_5insA]`)
        // are an order-ambiguous overlap conflict, not a normalizable form:
        // there is no canonical order for the two inserted sequences, so the
        // HGVS spec expresses "insert both" as a single ordered compound payload
        // (`ins[T;A]`), not two members (general.md:79). Every reference agrees —
        // mutalyzer rejects it (`EOVERLAP`), VariantValidator rejects it
        // (`AlleleSyntaxError: … ranges overlap`), and ferro's own strict mode
        // rejects it via the W5002 warning emitted just above (#486).
        //
        // In the non-strict modes the strict reject is not applied, so *preserve*
        // the allele as authored here rather than falling through to the
        // collapse/merge/shift pipeline below. That pipeline fabricated an
        // order-dependent merged insertion (`insAC`) which no other tool
        // produces and which then 3'-shifted flush against a neighbour and
        // collapsed to a `delins` only on re-normalization — leaving
        // `normalize` non-idempotent (#1004). Declining to canonicalize an
        // unresolvable overlap is both spec-consistent and idempotent, and it
        // still carries the W5002 warning so strict mode rejects unchanged.
        if merge::has_same_gap_insertions(&allele.variants, allele.phase) {
            let mut preserved =
                crate::hgvs::variant::AlleleVariant::new(allele.variants.clone(), allele.phase);
            preserved.uncertain = allele.uncertain;
            return Ok((HgvsVariant::Allele(preserved), all_warnings));
        }

        // Protein-axis adjacency (protein/substitution.md:23, delins.md:18):
        // two or more consecutive-residue changes (substitutions and/or
        // single-residue deletions) in one cis allele describe a single delins
        // (`p.[Arg76Ser;Cys77Trp]` → `p.Arg76_Cys77delinsSerTrp`;
        // `p.[Arg76Ser;Cys77del]` → `p.Arg76_Cys77delinsSer`), or a pure range
        // deletion when every member is a deletion. The positional nucleotide
        // merge below never fires on protein members, so
        // coalesce the protein axis here and re-dispatch the result (a bare
        // `Protein` when fully merged, else a smaller `Allele`) through the
        // normal pipeline. The helper is a no-op once no adjacent run remains,
        // so the re-dispatch cannot loop.
        if let Some(coalesced) = merge::coalesce_protein_adjacent_changes(allele) {
            let (normalized, mut warnings) = self.normalize_dispatch(&coalesced, manufactured)?;
            all_warnings.append(&mut warnings);
            return Ok((normalized, all_warnings));
        }

        // Collapse → merge → split → per-member normalize, iterated to a fixed
        // point. The per-member 3'-shift (inside `normalize_core`) can move an
        // insertion — or the `dup` it canonicalises to — flush against an
        // adjacent substitution/deletion. That *shift-created* adjacency is only
        // visible to `collapse_overlapping_cis_edits`, which runs on its input
        // members, on the *next* iteration; running a single pass left the two
        // as separate members (#999) and made `normalize` non-idempotent (#1000).
        //
        // Termination: each iteration's collapse/merge either strictly reduces
        // the member count or the per-member results already equal the input, so
        // the loop settles in at most a few passes. `max_passes` is a defensive
        // backstop against pathological oscillation, never the normal exit.
        //   - collapse: insertions flanking a deletion/sub at one locus → one
        //     delins; `merge_consecutive_edits` only handles strictly-consecutive
        //     non-overlapping edits (#487).
        //   - split (#160 + #165, cis only): a merged/pre-existing delins may
        //     decompose into higher-priority forms per `general.md:56` —
        //     `[..., inv, ...]` when a whole maximal contiguous run is a
        //     reverse complement (#160, corrected by #1034; sub-runs are not
        //     carved out) and/or into separate substitutions where interior
        //     positions match the reference (#165). The helper is a no-op
        //     otherwise.
        //   - per-member pipeline: the single canonical place where the 3' rule,
        //     ins→dup canonicalization, ref validation, etc. apply — a merged
        //     variant is semantically a new variant and goes through the same
        //     pipeline as any direct input. `normalize()` discards the diagnostic
        //     `infos` axis, so members use `normalize_core` (skipping the
        //     per-member `detect_shuffle_infos` cost).
        let is_cis = allele.phase == crate::hgvs::variant::AllelePhase::Cis;
        // Iterate the pipeline to a fixed point in cis; non-cis alleles are never
        // collapsible, so a single pass suffices. Iterating is sound here because
        // the one input that made it unsound — two insertion-like members at the
        // same gap, which are order-ambiguous and must not collapse (#487) and
        // whose per-member shift could move them apart to hide the collision — is
        // preserved-and-returned above as an overlap conflict (#1004) and so never
        // reaches this loop. `max_passes` is a defensive backstop.
        let max_passes = if is_cis {
            allele.variants.len().max(1) + 2
        } else {
            1
        };
        // #1103: sort cis members into genomic order *before* the merge, so
        // `merge_consecutive_edits` sees a canonical member order and fires the
        // same merges regardless of input order. #1101 sorts *after* the merge,
        // which restores order-independence for disjoint members but cannot undo
        // an order-dependent merge: `merge_consecutive_edits` only combines pairs
        // that are adjacent *in the input list*, so two permutations of one
        // adjacency-mergeable member set collapsed to different strings — and
        // different member counts (e.g. `g.[104_105insA;105del;105_106insC]`
        // merged all three into `g.105delinsAC`, but `g.[105del;104_105insA;
        // 105_106insC]` left them as three members). Feeding the merge a
        // genomic-order member list makes the merge deterministic; the two
        // flanking insertions of a co-located ins/del/ins trio land on distinct
        // junctions, so their genomic order (not input order) fixes the merged
        // `delins` payload unambiguously.
        //
        // Gated exactly like the post-merge sort (cis, certain, single-accession
        // — the latter checked inside the helper), plus: skip when the *input*
        // members already conflict. Reordering is only meaning-preserving for
        // disjoint members; an *input-time* overlap-conflicting allele must
        // preserve its authored order (#395), so it is left intact here and the
        // post-merge overlap guard handles it. (A conflict that emerges only
        // from the per-member 3' shift is not covered by this input-time gate:
        // such members are reordered by this sort and then rendered in the
        // deterministic descriptor tie-break order rather than authored order —
        // input-order-independent, but a narrowing of the #395 verbatim
        // contract for that shift-emergent case.) Same-junction insertions are
        // already returned verbatim above (#1004), so they never reach here.
        let input_has_overlap_conflict = all_warnings
            .iter()
            .any(|w| matches!(w, NormalizationWarning::OverlapConflict { .. }))
            || crate::normalize::overlap::detect_overlap_conflicts(&allele.variants, allele.phase)
                .iter()
                .any(|w| matches!(w, NormalizationWarning::OverlapConflict { .. }));
        let reorder_before_merge = is_cis && !allele.uncertain && !input_has_overlap_conflict;

        let mut current = allele.variants.clone();
        let mut normalized: Vec<HgvsVariant>;
        let settled_warnings: Vec<NormalizationWarning>;
        let mut pass = 0;
        loop {
            if reorder_before_merge {
                merge::sort_cis_members_for_merge(&mut current);
            }
            let pre_collapsed = merge::collapse_overlapping_cis_edits(
                current.clone(),
                allele.phase,
                &self.provider,
            );
            let merged_raw =
                merge::merge_consecutive_edits(pre_collapsed, allele.phase, &self.provider);
            // `canonical_split_for_variant` is the sequence-first
            // canonicalization (#1237), and like `normalize_core` below it acts
            // on one member with no sibling context. It can *reposition* a
            // member, not merely re-spell it: on a nine-`T` tract
            // `g.9_10delinsA` reduces to `g.1del` under a 5' shuffle, an
            // eight-position move. The sibling passes take their "before"
            // snapshot to decide what crossed what, so a snapshot taken *after*
            // this step shows no movement and they let the member sit wherever
            // it landed — including on top of a sibling (#1281).
            //
            // So record, per surviving member, the member it came from. A
            // one-for-one rewrite is a repositioning of the same member, and
            // its origin span is what the clamp needs. A genuine split into
            // several pieces has no single origin span per piece — the pieces
            // partition the original between them — so each piece stands for
            // itself, exactly as before.
            let (merged_split, pre_split): (Vec<HgvsVariant>, Vec<HgvsVariant>) = if is_cis {
                let mut split = Vec::with_capacity(merged_raw.len());
                let mut origins = Vec::with_capacity(merged_raw.len());
                for variant in merged_raw {
                    let origin = variant.clone();
                    let pieces = self.canonical_split_for_variant(variant, manufactured);
                    if pieces.len() == 1 {
                        origins.push(origin);
                    } else {
                        origins.extend(pieces.iter().cloned());
                    }
                    split.extend(pieces);
                }
                (split, origins)
            } else {
                // Non-cis has no origin snapshot to build: all three sibling
                // passes below return early on `phase != AllelePhase::Cis`, so
                // cloning the member list here would be pure waste on the
                // normalization hot path. All three also refuse a `before` whose
                // length differs from `after`, so an empty snapshot stays safe
                // if that phase gate ever moves.
                (merged_raw, Vec::new())
            };

            let mut result: Vec<HgvsVariant> = Vec::with_capacity(merged_split.len());
            let mut pass_warnings: Vec<NormalizationWarning> = Vec::new();
            for v in &merged_split {
                let (r, warnings) = self.normalize_core(v, manufactured)?;
                pass_warnings.extend(warnings);
                result.push(r);
            }

            // Four sibling-unaware decisions were taken per member above, and
            // all four are corrected here, in this order. The order matters:
            // each pass leaves the member set in the shape the next one expects.
            //
            // First: a deletion inside a tandem tract is re-spelled as a copy
            // count over the *whole tract* (`g.1_2del` -> `g.1_9T[7]` on a
            // nine-base run), and that widened span can swallow a sibling. Undo
            // it before the clamp runs, so the clamp sees a deletion span it can
            // pull back rather than a repeat it would skip.
            merge::demote_repeats_spanning_siblings(
                &pre_split,
                &mut result,
                allele.phase,
                allele.uncertain,
                &self.provider,
            );
            // Second: each member went to its standalone most-3' position, which
            // can carry it clear over a sibling's bases and make the pair denote
            // a different sequence (#1254). Pull any such member back to the last
            // position before the sibling — still the most 3' placement that
            // keeps the members on disjoint windows, and now merely *adjacent*,
            // so the next pass's merge coalesces the two.
            // `pre_split`, not `merged_split`: a member repositioned by the
            // canonicalization above must be measured from where it started,
            // not from where that step already put it (#1281).
            merge::clamp_sibling_crossing_shifts(
                &pre_split,
                &mut result,
                allele.phase,
                allele.uncertain,
                &self.provider,
            );
            // Third: an insertion or duplication consumes no base, so the clamp
            // above leaves it alone — but its *junction* shifts through a tract
            // too, and carrying it past a base a sibling edits changes what the
            // allele denotes. Bound the junction at the sibling's 5' edge, which
            // is still flush against it, so the #999 collapse keeps firing.
            merge::clamp_sibling_crossing_junctions(
                &pre_split,
                &mut result,
                allele.phase,
                allele.uncertain,
                &self.provider,
            );
            // Finally, a spelling repair rather than a repositioning: a
            // duplication whose span collides with a sibling's bases becomes
            // the equivalent insertion. `Xdup` claims position X, so a deletion
            // covering X makes the pair contradictory and ferro's own parser
            // rejects the result; `X_X+1ins<ref[X]>` denotes the same edit and
            // claims no base.
            //
            // Inside the loop, not after it. Run once at the end this produced
            // a non-idempotent result — the re-spelled allele is not a fixed
            // point, because the del and ins it exposes then cancel further
            // (`[261_262del;262_263insA;…]` reduces to `[261del;…]`). Feeding it
            // back through settles on that reduction instead.
            merge::respell_colliding_duplications(
                &mut result,
                allele.phase,
                allele.uncertain,
                &self.provider,
            );
            // Then, so the coalesce below can see them: two members that each
            // grew the SAME tract are both spelled as a copy count over it, and
            // a repeat carries no junction — so the pair claims one tract twice
            // and the coalesce never considers it (#1316). Re-spell each as the
            // insertion it stands for, which restores the junction. Deliberately
            // after the clamps rather than inside
            // `demote_repeats_spanning_siblings`: an insertion blocks no
            // sibling's shift, so demoting a repeat before the clamps would
            // release a sibling that #1296 clamped out of its tract.
            merge::demote_coincident_tract_repeats(
                &mut result,
                allele.phase,
                allele.uncertain,
                &self.provider,
            );
            // Last: two junction-occupying members that settled on the SAME
            // junction. The clamps above bound a member against a sibling's
            // bases, and an insertion or duplication has none, so a pair can
            // land on one interbase and claim it twice (#1286). Merge them into
            // one insertion carrying both payloads. This must happen here rather
            // than in the sequence-first pass, which derives from the allele's
            // resulting sequence — an overlapping allele has none, so that pass
            // declines exactly the input it would have fixed.
            merge::coalesce_members_at_one_junction(
                &pre_split,
                &mut result,
                allele.phase,
                allele.uncertain,
                &self.provider,
            );
            // Last, once every span above has settled: a member that cancelled
            // to `=` while a sibling grew over the bases it names is a
            // contradiction, and an overlap (#1297).
            merge::drop_identity_members_covered_by_siblings(
                &mut result,
                allele.phase,
                allele.uncertain,
                &self.provider,
            );

            // And last for the protein axis specifically: the coalesce above
            // this loop runs on the *authored* members, so it only ever sees
            // adjacency the input already had. The per-member 3' shift can
            // *create* it — `p.[Gly13del;Gly16Ala]` shifts the deletion down
            // its Gly run to `Gly17del`, which lands flush against `Gly16Ala`
            // — and the nucleotide `merge_consecutive_edits` at the top of the
            // pass never fires on protein members, so nothing closed that gap.
            // The result was that pass one emitted `p.[Gly16Ala;Gly17del]` and
            // pass two merged it to `p.Gly16_Gly17delinsAla`, i.e. `normalize`
            // was not idempotent (found by #1614's protein corpus axis under
            // `FERRO_ASSERT_IDEMPOTENT`).
            //
            // `substitution.md:32` is why the merged form is the right answer
            // rather than a matter of taste: at separation zero the split
            // spelling is marked `class="invalid"`, which the decided
            // `rulings[delins-adjacent-members-when-both-consume-reference]`
            // record scopes to members that both consume reference bases — a
            // substitution and a single-residue deletion both do.
            //
            // `general.md:157-160` is the authority for the other half — the
            // shift that *creates* the adjacency — and it is easy to miss
            // because it is stated as a Q&A rather than as a rule: "protein
            // variant descriptions should be derived from comparing the variant
            // protein sequence with the reference protein sequence. Knowledge on
            // the underlying change on the DNA level should not be used." Its
            // worked example is a 3'-rule application at the protein level: a
            // `Ser` lost from a `SerSer` run is `p.Ser5del`, and the passage
            // says in as many words that placing it at the deleted *codon* —
            // `p.Ser4del` — "is not correct". So carrying `Gly13del` down its
            // Gly run to `Gly17del` is required rather than incidental, and
            // `substitution.md:32` governs what is then made of the result. Both
            // halves are published clauses; neither rests on judgement.
            //
            // Inside the loop for the same reason `respell_colliding_duplications`
            // is: the coalesced allele is a new variant that must go back through
            // the per-member pipeline. Termination is unchanged — a coalesce
            // strictly reduces the member count, and the helper is a no-op once
            // no adjacent run remains.
            // Gated on the helper's *own* preconditions rather than on `is_cis`
            // alone. It declines an allele with fewer than two members or with
            // any non-protein member, so building the candidate under a bare
            // `is_cis` clones the member vector on every pass of every
            // nucleotide cis allele only to have the helper return `None`.
            let coalescible = is_cis
                && result.len() >= 2
                && result.iter().all(|v| matches!(v, HgvsVariant::Protein(_)));
            if coalescible {
                let mut candidate = crate::hgvs::variant::AlleleVariant::new(
                    result.clone(),
                    crate::hgvs::variant::AllelePhase::Cis,
                );
                candidate.uncertain = allele.uncertain;
                if let Some(coalesced) = merge::coalesce_protein_adjacent_changes(&candidate) {
                    result = match coalesced {
                        HgvsVariant::Allele(inner) => inner.variants,
                        single => vec![single],
                    };
                }
            }

            pass += 1;
            // Stable once a full pass leaves the member set unchanged.
            // `max_passes == 1` (non-cis) forces the original single-pass
            // behavior; otherwise iterate until the fixed point, with
            // `max_passes` as a defensive backstop.
            if result == current || pass >= max_passes {
                normalized = result;
                settled_warnings = pass_warnings;
                break;
            }
            current = result;
        }
        all_warnings.extend(settled_warnings);

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

        // Render cis-allele members in coordinate order, independent of input
        // order (#1098), so normalization produces a single canonical form:
        // leaving members in input order meant two inputs for the same allele
        // (`[a;b]` vs `[b;a]`) normalized to two different strings, breaking use
        // of the normalized descriptor as a stable key. Sort by a *total* order
        // (`cis_member_order_key`) so the result is deterministic even when two
        // members share a start.
        //
        // Basis: for DNA the spec discusses listing haplotype members "in
        // genomic order" (DNA/alleles.md, a discussion note); for protein it is
        // the machinery that makes the consecutive-residue delins reachable
        // order-independently (protein/substitution.md:23 requires
        // `p.[Arg76Ser;Cys77Trp]` to render as `p.Arg76_Cys77delinsSerTrp`
        // whatever order the two subs are authored in). Either way the aim is a
        // deterministic, input-order-independent canonical string.
        //
        // Only reorder when every member is on the **same accession**. The
        // spec's genomic-order rule is about variants on the same chromosome;
        // for a mixed-accession bracketed allele (`[NC_…;NM_…]`, #218/#219)
        // there is no cross-molecule genomic order to canonicalize to, and an
        // accession-string sort is not it — so those are left in authored order.
        //
        // Skip when the overlap passes above (`detect_insertion_overlaps`,
        // `detect_overlap_conflicts`) flagged an `OverlapConflict` **and** the
        // allele really was handed back untouched. Reordering is only
        // meaning-preserving for *disjoint* members; a conflict means two members
        // collide on one molecule, so they are NOT independent.
        // `merge_consecutive_edits` collapses strictly-adjacent members but does
        // not resolve a genuine overlap — that is exactly what those passes warn
        // about — so overlapping members can reach here. Leaving them in authored
        // order both avoids reordering edits whose combination is already
        // ill-defined and honors the `OverlapConflictingEdits` "preserves input
        // verbatim" contract (#395). Restricted to cis: trans / chimeric / mosaic
        // member order is not a genomic-order question, and uncertain groups
        // (`[(a;b)]`) are left untouched to preserve their authored form.
        // `sort_by` (not `sort_unstable_by`) is a belt-and-suspenders choice — the
        // key is already total, so stability is not relied upon.
        //
        // **The conflict alone does not earn the skip (#1414).** The reasoning
        // above is a contract about output the pipeline did not touch, so it
        // applies only when the output *is* untouched. A conflicting allele whose
        // members were all rewritten has no verbatim form left to preserve, and
        // skipping the sort for it produced two defects at once:
        // out-of-genomic-order members — a direct violation of #1235's criterion
        // 2 — and non-idempotence, since the rewritten output no longer conflicts
        // and so takes the *other* branch of this gate on the second pass.
        // Measured on `origin/main` at `94817bf6`, 41 inputs in four shape
        // families reach it (`[inv;insAA]` 16, `[inv;insA]` 16, `[dup;inv]` 7,
        // `[insA;inv]` 2); #1423 had
        // repaired only the single example the issue cited, because it drops an
        // identity a sibling *overlaps* and a zero-width junction insertion
        // overlaps nothing.
        //
        // So the predicate is the conjunction, not the conflict alone: an
        // allele emitted verbatim keeps its authored order, and everything else
        // is sorted exactly as before.
        let has_overlap_conflict = all_warnings
            .iter()
            .any(|w| matches!(w, NormalizationWarning::OverlapConflict { .. }));
        let emitted_verbatim = normalized == allele.variants;
        if is_cis
            && !allele.uncertain
            && !(has_overlap_conflict && emitted_verbatim)
            && normalized.len() > 1
        {
            sort_cis_members_by_genomic_order(&mut normalized);
        }

        // A conflict the input really had must survive into the output (#1406).
        //
        // #395 says an overlap-conflicting allele has no canonical form, so
        // ferro preserves it verbatim and reports `W5002`, which strict mode
        // promotes to an error. But strict mode re-reads a *description*, and it
        // has no memory of what that description was normalized from — so the
        // contract only holds if the emitted description still looks conflicting.
        //
        // It did not. The per-member pipeline repairs the members one at a time
        // — an inversion whose span is its own reverse complement cancels to an
        // identity, an interior insertion respells as a repeat — and the merge
        // then collapses what is left. The conflict is gone from the output, so
        // strict mode accepts a description it would have rejected had the
        // caller written it directly:
        //
        //     in      g.[11_12inv;11_12insAA]   strict: REJECT (W5002)
        //     lenient g.10_11A[4]               strict: ACCEPT
        //
        // That is the laundering #1406 row 3 reports, and it breaks the
        // invariant `issue_1276_dup_junction_overlap::
        // lenient_output_of_a_conflict_is_still_rejected_by_strict` already
        // pins for the shapes it happens to cover.
        //
        // So: when the input conflicts and the output does not, hand back the
        // input. That is #395's contract stated exactly — preserve it verbatim —
        // rather than approximated by "skip some passes and hope the conflict
        // survives them".
        //
        // Deliberately asked of the *members*, both detectors, on each side.
        // `detect_insertion_overlaps` documents that it must see pre-merge
        // members because the merge collapses the overlap out of view — and
        // that collapse is precisely the erasure being detected here, so asking
        // it of the post-merge output is the right question, not a misuse.
        //
        // Idempotent by construction: the value handed back is the input, whose
        // conflict the next pass detects again, so it is returned unchanged a
        // second time.
        // Uncertain alleles included, unlike the sort gate above. That gate
        // excludes them because member *order* inside `[(a;b)]` is authored
        // presentation, which is not ours to change. This one is about whether
        // the conflict survives into the output at all, and an uncertain allele
        // launders exactly like a certain one — measured:
        //
        //     in  g.[(9_10insA;9_10inv)]   strict: REJECT (W5002)
        //     out g.[(11dup;9_10=)]        strict: ACCEPT
        //
        // Excluding them would leave the hole this gate exists to close open on
        // one spelling of the same input.
        if is_cis {
            let conflicts = |members: &[HgvsVariant]| {
                !crate::normalize::overlap::detect_overlap_conflicts(members, allele.phase)
                    .is_empty()
                    || !crate::normalize::overlap::detect_interior_junction_conflicts(
                        members,
                        allele.phase,
                    )
                    .is_empty()
            };
            // A **plural** output only. Both detectors return nothing for fewer
            // than two members (`overlap.rs`, first statement of each), so for a
            // single-member output `!conflicts(&normalized)` is vacuously true
            // and says nothing about whether the conflict was erased — there is
            // simply nothing left for a second member to collide with.
            //
            // The distinction is real, not defensive. Members *shifted apart*
            // while still plural is the laundering this gate is for. Members
            // *composed into one edit* is the opposite: the colliding writes
            // have been resolved into a single description that denotes a
            // definite sequence, which is what `delins.md:86-89` asks for and
            // what the deletion exemption below already relies on for
            // `g.[2_3del;2_3insAA]` -> `g.2_3delinsAA`.
            //
            // Without this term the gate reverts #1423: `g.[11_12inv;11_12insAA]`
            // collapses to the single member `g.10_11A[4]`, which cannot
            // conflict, so the input was handed back and a merged, shipped form
            // was undone.
            if conflicts(&allele.variants) && normalized.len() > 1 && !conflicts(&normalized) {
                normalized = allele.variants.clone();
            }
        }

        // HGVS requires consecutive edits in cis to render as a single edit.
        // Only unwrap when a merge actually collapsed multiple sub-variants —
        // pre-existing singleton alleles must round-trip with the Allele
        // wrapper intact for programmatic callers (Display already renders
        // singletons in bare form regardless).
        // Never unwrap a predicted/uncertain cis allele (`[(a;b)]`): even if
        // the members merge down to one, the `uncertain` flag and the `[(...)]`
        // notation must be preserved.  Dropping the wrapper would silently
        // change semantics and break round-trip fidelity.
        let result = if allele.phase == crate::hgvs::variant::AllelePhase::Cis
            && !allele.uncertain
            && original_len > 1
            && normalized.len() == 1
        {
            normalized.pop().unwrap()
        } else {
            let mut rebuilt = crate::hgvs::variant::AlleleVariant::new(normalized, allele.phase);
            // Preserve the uncertainty wrapper `(...)` (e.g. the and/or
            // group `c.(370A>C^372C>R)`) — rebuilding via `new` would
            // otherwise drop it and re-render the group expanded.
            rebuilt.uncertain = allele.uncertain;
            HgvsVariant::Allele(rebuilt)
        };

        Ok((result, all_warnings))
    }

    /// Issue #333: expand a bracketed / reference-range `ins[...]` payload
    /// in an `Insertion` / `Delins` / `DupIns` edit to a flat literal
    /// sequence. Returns the rewritten edit plus a
    /// [`NormalizationWarning::InsertedSequenceExpanded`] for
    /// observability, or `None` when nothing was canonicalized.
    ///
    /// The five `normalize_<axis>` methods — genome, cds, tx, rna (#1183), and
    /// mt — wrap this helper to build the per-axis variant struct, one
    /// `try_expand_<axis>_ins` each. Splitting the warning construction here
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
        // only acts on `InsertedSequence::Complex`. Since #856 a single
        // bracketed *literal* (`[ATC]`) is collapsed to a plain `Literal` at
        // parse time, so the surviving `Complex` payloads here are genuinely
        // multi-part (`[A;T]`) or single non-literal (`[100_110]`).
        // `to_bracketed_string` forces the bracketed form so the warning
        // preserves the user's input shape (`[A;100_110]`) even though
        // `Display` on `InsertedSequence::Complex` drops brackets for
        // single-element vectors (spec-canonical form).
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

    /// `normalize_rna` companion — `r.` numbering follows the associated DNA
    /// reference, so it is CDS-relative (== `c.`) on a coding transcript and
    /// transcript-relative (== `n.`) on a non-coding one. That split needs the
    /// provider, so it is carried by `InsCoordKind::Rna` and resolved at fetch
    /// time in `fetch_position_range_bases` — the same kind the `r.`
    /// cross-reference *payload* path already uses (#773). No `U`→`T` step is
    /// needed here: the payload resolves to DNA bases from the provider, and
    /// `RnaVariant`'s Display renders `T`→`u` again on output. #1183.
    fn try_expand_rna_ins(
        &self,
        variant: &crate::hgvs::variant::RnaVariant,
        edit: &NaEdit,
        accession: &str,
    ) -> Result<Option<(crate::hgvs::variant::RnaVariant, NormalizationWarning)>, FerroError> {
        let Some((new_edit, warning)) =
            self.try_expand_ins_edit(edit, accession, InsCoordKind::Rna)?
        else {
            return Ok(None);
        };
        let new_variant = crate::hgvs::variant::RnaVariant {
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

/// Widest **half**-width [`Normalizer::normalize_in_grown_window`] will grow the
/// fetch to before refusing to shift any further.
///
/// **This bounds `window`, not the span read.** `normalize_in_grown_window`
/// fetches `[start - window, end + window]`, so the widest reference span is
/// roughly *twice* this plus the variant's own length — measured, a 60,000-base
/// tract at contig offset 100,000 converges to `g.160000del` off a window
/// covering `[34,465 .. 165,537]`, 131,073 bases. The *travel* bound is the
/// half-width itself: a 65,000-base tract converges and a 66,000-base one is
/// refused. [`MAX_INTRONIC_SHUFFLE_WINDOW`] at its own definition is tested as a
/// **full** span, so the two constants are numerically equal and semantically a
/// factor of two apart; do not read one off the other.
///
/// The value is **reused from** [`MAX_INTRONIC_SHUFFLE_WINDOW`] rather than
/// independently sourced — that constant is this module's existing answer to the
/// same question (how far a shuffle may fetch before the cost stops being worth
/// the answer) and is itself unsourced. What justifies the reuse is that a run
/// long enough to exhaust it is a structural feature (an assembly `N` gap, a
/// synthetic contig) rather than a repeat tract: the pathogenic expansions that
/// do run long — FMR1 `CGG`, RFC1 `AAGGG` — reach thousands of bases, which is
/// still an order of magnitude inside 64 KiB. It is a cost bound picked to sit
/// clear of the biology, not a bound derived from it.
const MAX_SHUFFLE_FETCH_WINDOW: u64 = MAX_INTRONIC_SHUFFLE_WINDOW;

/// Next window to try, or `None` once [`MAX_SHUFFLE_FETCH_WINDOW`] has been
/// tried. Doubling keeps the number of re-fetches logarithmic in the distance
/// the shuffle actually travels, so the overwhelmingly common case — a result
/// that rests well inside the configured window — costs exactly one fetch and
/// never consults this function at all.
///
/// `max(1)` guards a `window_size` of 0, which would otherwise double to itself
/// forever.
fn grow_shuffle_fetch_window(window: u64) -> Option<u64> {
    if window >= MAX_SHUFFLE_FETCH_WINDOW {
        return None;
    }
    Some(
        window
            .max(1)
            .saturating_mul(2)
            .min(MAX_SHUFFLE_FETCH_WINDOW),
    )
}

/// Bases of reference kept either side of an edit's own span by
/// [`Normalizer::canonicalize_without_shifting`].
///
/// Every type rule that survives the growth cap is decided within one payload
/// length of the edit: `duplication.md:5` defines a duplication as a copy
/// inserted "**directly 3'** of the original copy", so seeing the original costs
/// exactly `|payload|` bases of 5' flank, and the mirror check on the 3' side
/// costs the same. Delins trimming, the delins type ladder and the inversion
/// shortening all read inside the span itself.
///
/// **The value must depend on the edit and on nothing else.** It is what bounds
/// how far the capped answer can look, so a flank derived from the reference —
/// a tract length, a window width — would make the answer depend on the very
/// distance the cap declined to travel, and the walk would come straight back.
/// The `max(1)` floor keeps a zero-payload edit (`del`, `inv`) from being handed
/// a slice with no context at all.
fn capped_typing_flank(edit: &NaEdit) -> u64 {
    let payload = match edit {
        NaEdit::Insertion { sequence } | NaEdit::Delins { sequence, .. } => {
            sequence.bases().map_or(0, <[Base]>::len)
        }
        _ => 0,
    };
    (payload as u64).max(1)
}

/// One settled attempt at normalizing an edit against a fetched reference
/// window. See [`Normalizer::normalize_in_grown_window`].
struct WindowedNormalization {
    /// 0-based offset of the window's first base within the contig. Positions
    /// convert back with `contig = window_relative + window_start`.
    window_start: u64,
    /// The fetched bases.
    ref_seq: String,
    /// Window-relative 1-based start of the normalized edit.
    new_rel_start: u64,
    /// Window-relative 1-based end of the normalized edit.
    new_rel_end: u64,
    /// The normalized edit.
    new_edit: NaEdit,
    /// Warnings raised while normalizing it.
    warnings: Vec<NormalizationWarning>,
}

impl<P: ReferenceProvider> Normalizer<P> {
    /// Clamp a window-fetch upper bound to the contig length.
    ///
    /// The `±window_size` fetch window around a variant can run past the
    /// contig 3' end for a variant near the end. Every provider ERRORS on a
    /// past-EOF read rather than clamping, so an unclamped `raw_end` drops the
    /// whole variant into the minimal-notation fallback — skipping 3' shifting
    /// and the `delins -> inv/sub/dup` canonicalization (#1041).
    ///
    /// Clamp to the contig length, but ONLY when the variant itself fits within
    /// the contig (`end <= len`). If the variant *span* runs past the contig
    /// end (`end > len`), clamping would fetch a window shorter than the span
    /// and `normalize_na_edit` would read a truncated reference and
    /// mis-normalize (e.g. `g.99_103inv` on a 100 bp contig collapsing to
    /// `g.99_103=`); keep the raw window so the read errors into the safe
    /// pass-through fallback instead. When the length is unavailable, likewise
    /// keep the raw window.
    ///
    /// Shared by `normalize_genome` (#1042) and the mitochondrial/circular
    /// path (#1044) so the two cannot silently drift — the divergence that
    /// caused #1044 in the first place.
    fn clamp_fetch_end_to_contig(&self, accession: &str, end: u64, raw_end: u64) -> u64 {
        match self.provider.get_sequence_length(accession) {
            Ok(len) if end <= len => raw_end.min(len),
            _ => raw_end,
        }
    }

    /// Normalize an edit against a reference window sized so the shuffle can
    /// actually reach its fixed point (#1691).
    ///
    /// The `±window_size` fetch is a *cost* bound, but `normalize_na_edit` is
    /// handed `Boundaries::new(0, ref_seq.len())` and so treats the window's far
    /// edges as hard limits on the shuffle. Those edges are where the fetch
    /// stopped, not properties of the contig, so inside a repeat tract longer
    /// than the window the shuffle reports the edge it reached rather than the
    /// tract's end. Re-normalizing that answer re-centres the window on it and
    /// advances another `window_size` bases — `g.100del -> g.200del ->
    /// g.300del -> …`, which never terminates and trips
    /// `FERRO_ASSERT_IDEMPOTENT`.
    ///
    /// So the window is grown geometrically until the result rests strictly
    /// inside it, and the caller gets whichever attempt settled. Growth is only
    /// attempted when the edge is *provably* artificial — the contig is known to
    /// continue past the window's 3' end, or the window's 5' end is not the
    /// contig's first base. Where flushness cannot be established (no length
    /// from the provider) the window is left alone, which is the same
    /// conservative choice [`Self::clamp_fetch_end_to_contig`] makes.
    ///
    /// # Returns
    ///
    /// `Ok(None)` when the caller should fall back to minimal notation: **any**
    /// fetch failed — the first one or a grown one (a grown failure is not
    /// allowed to fall back to the truncated attempt that preceded it, which
    /// would reinstate the walk) — or the growth cap was reached by an edit
    /// [`Self::canonicalize_without_shifting`] declines to re-type.
    ///
    /// **Refusing to shift is the only capped answer that is still a fixed
    /// point.** Capping the *shift* instead — returning the far edge of the
    /// largest window we were willing to fetch — just changes the walk's step
    /// size from `window_size` to the cap; the second call re-centres on the new
    /// position and walks again.
    ///
    /// **What is refused at the cap is the travel, and only the travel.** The
    /// edit is still re-typed against the reference by
    /// [`Self::canonicalize_without_shifting`], because the type rules —
    /// `ins -> dup`, `delins -> sub/inv/dup/del`, the inversion shortening — are
    /// decided by the bases flanking the edit rather than by how far a tract
    /// runs, so they are answerable at the cap and `DNA/duplication.md:18` says
    /// they must be answered ("when a variant can be described as a duplication,
    /// it **must** be described as a duplication and not as, e.g., an
    /// insertion"). Refusing them alongside the shift bought idempotency with
    /// conformance, which is the one thing the README ruleset never trades.
    ///
    /// `fetch_end_of` maps a raw `end + window` upper bound to the one the
    /// caller wants read (each axis clamps differently), and is re-consulted on
    /// every attempt because the raw bound changes as the window grows.
    fn normalize_in_grown_window(
        &self,
        accession: &str,
        start: u64,
        end: u64,
        edit: &NaEdit,
        fetch_end_of: impl Fn(u64) -> u64,
    ) -> Result<Option<WindowedNormalization>, FerroError> {
        // Whether the contig continues past a window is a question about the
        // contig, so the length is read once rather than per attempt.
        let contig_len = self.provider.get_sequence_length(accession).ok();
        let mut window = self.config.window_size;

        loop {
            let window_start = start.saturating_sub(window);
            let fetch_end = fetch_end_of(end.saturating_add(window));
            let Ok(ref_seq) = self
                .provider
                .get_sequence(accession, window_start, fetch_end)
            else {
                // A *grown* fetch that fails is reported the same as a failed
                // first fetch rather than falling back to the truncated attempt,
                // which would reinstate the walk this method exists to stop.
                return Ok(None);
            };

            // Window-relative 1-based positions; `hgvs_pos_to_index` inside
            // `normalize_na_edit` takes them the rest of the way.
            let (new_rel_start, new_rel_end, new_edit, warnings) = self.normalize_na_edit(
                ref_seq.as_bytes(),
                edit,
                start - window_start,
                end - window_start,
                &Boundaries::new(0, ref_seq.len() as u64),
                // Neither the genomic nor the mitochondrial axis carries a
                // reading frame; both callers passed this before the growth
                // loop was factored out of them.
                CodonGate::NotApplicable,
            )?;

            let seq_len = ref_seq.len() as u64;
            // The 3' edge is artificial only where the contig is known to
            // continue past it; the 5' edge only where the window did not start
            // at the contig's first base.
            let ran_to_three_prime_edge = new_rel_end >= seq_len
                && contig_len.is_some_and(|len| window_start.saturating_add(seq_len) < len);
            let ran_to_five_prime_edge = new_rel_start <= 1 && window_start > 0;

            if !ran_to_three_prime_edge && !ran_to_five_prime_edge {
                return Ok(Some(WindowedNormalization {
                    window_start,
                    ref_seq,
                    new_rel_start,
                    new_rel_end,
                    new_edit,
                    warnings,
                }));
            }

            match grow_shuffle_fetch_window(window) {
                Some(next) => window = next,
                // Cap reached. Refuse the travel, not the typing: the sequence
                // for the latter is already in hand, so it is re-derived from
                // this last window rather than re-fetched.
                None => {
                    return self.canonicalize_without_shifting(
                        &ref_seq,
                        window_start,
                        start,
                        end,
                        edit,
                    )
                }
            }
        }
    }

    /// Re-type an edit against the reference **without letting it travel**, for
    /// the growth cap in [`Self::normalize_in_grown_window`].
    ///
    /// A normalized description answers three separable questions, and running
    /// out of window costs it only two of them:
    ///
    /// - **Travel** — the 3' rule's resting place. Unanswerable past the cap by
    ///   construction, and refusing it is what makes the capped answer a fixed
    ///   point.
    /// - **Extent** — repeat notation's `<first>_<last>unit[N]`, which states
    ///   where the reference tract begins and ends. Equally unanswerable, and a
    ///   *truncated* extent is worse than none: it is a false claim about the
    ///   reference rather than a merely unshifted one. Declined below.
    /// - **Type** — `ins -> dup`, `delins -> sub/inv/dup/del`, the inversion
    ///   shortening. These read the bases flanking the edit and nothing else, so
    ///   the cap does not touch them and `DNA/duplication.md:18` requires them.
    ///
    /// So the edit is re-normalized against a **slice** of the window whose
    /// width is derived from the edit alone — never from the tract — with the
    /// shuffle pinned by [`Boundaries::pinned`]. Both halves are load-bearing
    /// for idempotency. Pinning alone is not enough: `insertion_to_duplication`
    /// consults no boundaries, so on the full window it anchors the `dup` at the
    /// far end of whatever tract it can see, which is the walk again wearing a
    /// different edit type. Slicing alone is not enough either, for the same
    /// reason the cap exists at all. Together, every base the answer reads sits
    /// at a fixed offset from the edit the caller handed in, so re-normalizing
    /// that answer reads the same relative bases and returns it unchanged.
    ///
    /// The result is returned as an ordinary [`WindowedNormalization`] over the
    /// slice, so the callers' coordinate reconstruction, boundary clamp and
    /// canonical split all apply to it unchanged — the capped answer is a
    /// normalization like any other, not a second output path to keep in sync.
    ///
    /// `Ok(None)` means "no re-typing was available", and the caller falls back
    /// to the syntactic minimal-notation cleanup exactly as it does for a failed
    /// fetch.
    fn canonicalize_without_shifting(
        &self,
        ref_seq: &str,
        window_start: u64,
        start: u64,
        end: u64,
        edit: &NaEdit,
    ) -> Result<Option<WindowedNormalization>, FerroError> {
        // A repeat's canonical form IS an extent claim, so there is no
        // travel-free part of it to salvage: `normalize_repeat` would report the
        // slice's own edges as the tract's. Hand it back untouched.
        if matches!(edit, NaEdit::Repeat { .. } | NaEdit::MultiRepeat { .. }) {
            return Ok(None);
        }

        let flank = capped_typing_flank(edit);
        // Window-relative, 0-based: `p` sits at index `p - window_start - 1`.
        let span_start_idx = (start - window_start).saturating_sub(1);
        let span_end_idx = end - window_start;
        let window_len = ref_seq.len() as u64;
        // The slice has to hold the whole span, or the type rules would read a
        // truncated edit and answer about a variant nobody described. The growth
        // loop normalized this same span against this same window a moment ago,
        // so a span that does not fit is a contradiction rather than a case —
        // decline instead of clamping into a wrong answer.
        if span_end_idx > window_len || span_start_idx >= span_end_idx {
            return Ok(None);
        }
        let lo = span_start_idx.saturating_sub(flank);
        let hi = span_end_idx.saturating_add(flank).min(window_len);
        let slice = &ref_seq[lo as usize..hi as usize];
        // The slice's own 0-based offset within the contig, which is what the
        // caller adds back to the returned relative positions.
        let slice_start = window_start + lo;

        let (new_rel_start, new_rel_end, new_edit, warnings) = self.normalize_na_edit(
            slice.as_bytes(),
            edit,
            start - slice_start,
            end - slice_start,
            &Boundaries::pinned(),
            // Same reasoning as the growth loop's own call: neither axis
            // reaching here carries a reading frame.
            CodonGate::NotApplicable,
        )?;

        // Extent claim, arrived at from the other direction — an `ins` or `dup`
        // the type rules promoted into repeat notation. The promotion itself is
        // right; the tract it would name is the slice's, so decline it here
        // rather than publish a span the reference does not have.
        if matches!(new_edit, NaEdit::Repeat { .. } | NaEdit::MultiRepeat { .. }) {
            return Ok(None);
        }

        Ok(Some(WindowedNormalization {
            window_start: slice_start,
            ref_seq: slice.to_string(),
            new_rel_start,
            new_rel_end,
            new_edit,
            warnings,
        }))
    }

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

        // #1052: an uncertain/predicted-wrapped substitution must stay a
        // silent pass-through — see `is_uncertain_real_substitution`'s doc
        // comment.
        if is_uncertain_real_substitution(edit, &variant.loc_edit.edit) {
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
            // Unresolved. `cen` is a structural *refusal*: preserve the input
            // verbatim (do not rewrite the edit body) and surface W4005 so it is
            // not silently echoed. A length-less qter/pter is instead an
            // environment gap — a genuine resolution attempt that fell back, so
            // it takes the silent canonicalize path.
            _ => {
                if matches!(start_pos.special, Some(SpecialPosition::Cen))
                    || matches!(end_pos.special, Some(SpecialPosition::Cen))
                {
                    return Ok((
                        HV::Genome(variant.clone()),
                        vec![NormalizationWarning::UnresolvableSpecialPosition {
                            accession: accession.clone(),
                            marker: "cen".to_string(),
                        }],
                    ));
                }
                return Ok((
                    HV::Genome(self.canonicalize_genome_variant(variant)),
                    vec![],
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
                loc_edit: LocEdit::with_uncertainty(
                    Interval::new(GenomePos::new(start), GenomePos::new(end)),
                    variant.loc_edit.edit.with_same_certainty(edit.clone()),
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
        // path uses (do not "fix" this as a bug).
        //
        // Clamp the upper bound to the contig length so the read is well-formed
        // against providers that ERROR on a past-EOF read (#1041); see
        // `clamp_fetch_end_to_contig` for the `end <= len` gate and rationale.
        //
        // Special path (unchanged): it already resolved `resolved_len` and
        // clamps unconditionally — a mixed special/plain past-end span like
        // `g.pter_<past-end>del` fetches the whole contig and relies on
        // `shuffle`'s per-index bounds guard to echo the input verbatim
        // (see `genome_mixed_special_plain_past_end_matches_plain_path`).
        //
        // #1691: the window also bounds the SHUFFLE, and its far edges are not
        // contig properties, so a tract longer than the window walks instead of
        // converging. `normalize_in_grown_window` re-fetches until the result
        // settles strictly inside; `None` means either no sequence or a shuffle
        // still running at the growth cap, and both take the minimal-notation
        // fallback that a failed fetch always took.
        let attempt = self.normalize_in_grown_window(&accession, start, end, edit, |raw_end| {
            if had_special {
                raw_end.min(resolved_len)
            } else {
                self.clamp_fetch_end_to_contig(&accession, end, raw_end)
            }
        })?;
        let Some(WindowedNormalization {
            window_start,
            ref_seq,
            mut new_rel_start,
            mut new_rel_end,
            mut new_edit,
            mut warnings,
        }) = attempt
        else {
            return Ok((
                HV::Genome(self.canonicalize_genome_variant(variant)),
                vec![],
            ));
        };

        // #1205: an insertion driven to rest against a contig end has no valid
        // pair of adjacent anchors there, exactly as on the transcript axes
        // (#1202) — 5' it escapes as the single-position `g.1ins<A'>`, which
        // `DNA/insertion.md:95-101` rejects by name and ferro's own parser
        // refuses; 3' as `g.<len>_<len+1>ins<A'>`, whose second anchor is past
        // the contig end. Rewrite both to the equivalent `delins`.
        //
        // `ref_seq` is a WINDOW into the contig, not the whole thing, so its ends
        // are contig boundaries only where the window is flush with the contig —
        // otherwise the "boundary" is just where the fetch stopped, and clamping
        // there would invent one. The 5' side is flush iff the window started at
        // the contig's first base; the 3' side iff the window's last base is the
        // contig's last. With no length available we cannot establish flushness,
        // so we do not clamp — the same conservative choice
        // `clamp_fetch_end_to_contig` makes on the read side.
        //
        // Runs on window-relative coordinates, before the conversion back to
        // genomic ones, mirroring the ordering the three transcript callers hold
        // (clamp in the frame `ref_seq` is indexed in, then convert).
        let contig_len = self.provider.get_sequence_length(&accession).ok();
        clamp_insertion_at_sequence_bounds(
            ref_seq.as_bytes(),
            &mut new_edit,
            &mut new_rel_start,
            &mut new_rel_end,
            SequenceEnds {
                five_prime: window_start == 0,
                three_prime: contig_len
                    .is_some_and(|len| window_start + ref_seq.len() as u64 >= len),
            },
        );

        // Adjust back to genomic coordinates
        let new_start = new_rel_start + window_start;
        let new_end = new_rel_end + window_start;

        // Reconstruct variant with new position (using adjusted coordinates)
        let new_variant = GenomeVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                Interval::new(GenomePos::new(new_start), GenomePos::new(new_end)),
                variant.loc_edit.edit.with_same_certainty(new_edit),
            ),
        };

        // Issue #160 + #165: a normalized Delins may decompose under the
        // spec's edit-priority rule (`general.md:56`) — into `[..., inv,
        // ...]` when a whole maximal contiguous run is a reverse complement
        // (#1034: rev-comp sub-runs are not carved out) and/or into separate
        // subs across interior identities. Returns the variant unchanged for
        // non-Delins or no-decomposition cases.
        let uncertain = new_variant.loc_edit.edit.is_uncertain();
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Genome(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split, uncertain), warnings))
    }

    /// Normalize a CDS variant
    fn normalize_cds(
        &self,
        variant: &CdsVariant,
        manufactured: &mut Vec<ManufacturedJunctionExit>,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        // #972 Task 5: a transcript whose 5' CDS is annotated incomplete
        // (`cds_start_NF`) has no confirmed ATG, so `c.1` is undefined
        // against it and it is not an HGVS-recommended `c.` reference
        // (W5004). Decline to re-number/normalize — pass the input through
        // verbatim — in EVERY mode; only whether the warning surfaces (and
        // so whether `normalize()` promotes it to a hard reject) is
        // mode-gated below. Must run before every other check (bounds,
        // telomere resolution, `ins[...]` expansion, 3' shuffle, canonical
        // split) since none of those can be trusted without a confirmed CDS
        // start. Input-side counterpart to Task 4's projection-side decline.
        if let Ok(transcript) = self
            .provider
            .get_transcript(&variant.accession.transcript_accession())
        {
            if transcript.cds_start_incomplete {
                let mut warnings = Vec::new();
                if self.config.should_warn_incomplete_cds_start()
                    || self.config.should_reject_incomplete_cds_start()
                {
                    warnings.push(NormalizationWarning::IncompleteCdsStartReference {
                        accession: variant.accession.transcript_accession(),
                        coordinate_system: "c".to_string(),
                    });
                }
                return Ok((HV::Cds(variant.clone()), warnings));
            }
        }

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
            return self.normalize_cds(&new_variant, manufactured);
        }

        // Resolve telomere/centromere markers (pter/qter) to concrete CdsPos
        // BEFORE the bounds gate; recurse so the normal pipeline sees concrete
        // positions. `cen` (unresolvable) takes the W4005 path. (#488 Phase 2)
        let start_special = variant
            .loc_edit
            .location
            .start
            .inner()
            .map(|p| p.is_special())
            .unwrap_or(false);
        let end_special = variant
            .loc_edit
            .location
            .end
            .inner()
            .map(|p| p.is_special())
            .unwrap_or(false);
        if start_special || end_special {
            // Phase 2b (#488): on a genomic-reference (NG_/NC_/LRG_) c. description,
            // pter/qter denote the genomic parent's terminus, in the 5'/3' transcript
            // flank — not numberable in c. per HGVS numbering.md (the flank-numbering
            // proposal was rejected, open-issues.md). Refuse: leave verbatim + emit
            // W4006 (strict mode rejects in the normalize() wrapper). Read
            // genomic_context from the ORIGINAL accession before transcript_accession()
            // strips it. Only pter/qter (not cen) hit this branch.
            if variant.accession.genomic_context.is_some() {
                let marker_of = |p: Option<&CdsPos>| match p.and_then(|p| p.special) {
                    Some(SpecialPosition::Pter) => Some("pter"),
                    Some(SpecialPosition::Qter) => Some("qter"),
                    _ => None,
                };
                if let Some(marker) = marker_of(variant.loc_edit.location.start.inner())
                    .or_else(|| marker_of(variant.loc_edit.location.end.inner()))
                {
                    let warnings = vec![NormalizationWarning::TranscriptFlankNotDescribable {
                        accession: format!("{}", variant.accession),
                        marker: marker.to_string(),
                    }];
                    // Refusal is verbatim: preserve the original c. variant
                    // exactly. `canonicalize_cds_variant` rewrites edit bodies
                    // (e.g. an explicit `delA` loses its deleted base), which
                    // would break the W4006 leave-as-is contract.
                    return Ok((HV::Cds(variant.clone()), warnings));
                }
            }
            let acc = variant.accession.transcript_accession();
            // `cen` has no numberable CDS coordinate, so it is unresolvable with
            // or without transcript metadata. Surface W4005 up front, before the
            // transcript lookup — otherwise a missing transcript silently
            // swallows it and strict mode wrongly accepts `c.cendel`. (#488)
            let is_cen = |p: Option<&CdsPos>| {
                matches!(p.and_then(|p| p.special), Some(SpecialPosition::Cen))
            };
            if is_cen(variant.loc_edit.location.start.inner())
                || is_cen(variant.loc_edit.location.end.inner())
            {
                // Structural refusal: preserve the input verbatim (don't rewrite
                // the edit body) and surface W4005. `canonicalize_cds_variant`
                // would normalize the edit (e.g. drop an explicit `delA` base).
                return Ok((
                    HV::Cds(variant.clone()),
                    vec![NormalizationWarning::UnresolvableSpecialPosition {
                        accession: acc.clone(),
                        marker: "cen".to_string(),
                    }],
                ));
            }
            let transcript = match self.provider.get_transcript(&acc).ok() {
                Some(t) => t,
                None => return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![])),
            };
            let start_in = variant.loc_edit.location.start.inner();
            let end_in = variant.loc_edit.location.end.inner();
            let rs = start_in
                .map(|p| resolve_special_cds_pos(p, &transcript))
                .transpose()?;
            let re = end_in
                .map(|p| resolve_special_cds_pos(p, &transcript))
                .transpose()?;
            match (rs, re) {
                (Some(Some(s)), Some(Some(e))) => {
                    let new_variant = CdsVariant {
                        accession: variant.accession.clone(),
                        gene_symbol: variant.gene_symbol.clone(),
                        loc_edit: LocEdit::with_uncertainty(
                            Interval::new(s, e),
                            variant.loc_edit.edit.clone(),
                        ),
                    };
                    return self.normalize_cds(&new_variant, manufactured);
                }
                _ => {
                    // A pter/qter bound did not project (e.g. a non-coding
                    // transcript with no CDS frame). `cen` was already handled
                    // and returned above, so this is a pure projection gap:
                    // preserve the input verbatim with no warning.
                    return Ok((HV::Cds(self.canonicalize_cds_variant(variant)), vec![]));
                }
            }
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
            || -> Result<std::sync::Arc<crate::reference::transcript::Transcript>, FerroError> {
                // Resolve directly from the accession — the build-aware lookup
                // only needs the accession, so we avoid cloning the whole
                // variant just to satisfy the by-variant signature.
                self.provider
                    .get_transcript_for_accession(&variant.accession)
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
                        // Render any plain past-CDS coordinate in canonical
                        // `c.*N` form (the input used out-of-scheme numbering;
                        // the W4004 warning above still flags it, and strict
                        // mode still rejects). Only a simple certain
                        // single/point interval is rewritten; a
                        // genuinely-past-transcript position is left untouched
                        // by the helper. See #920/#336.
                        let new_start = canonicalize_pastcds_pos_to_utr3(start_pos, transcript);
                        let new_end = canonicalize_pastcds_pos_to_utr3(end_pos, transcript);
                        let both_certain = variant
                            .loc_edit
                            .location
                            .start
                            .as_single()
                            .is_some_and(|m| m.is_certain())
                            && variant
                                .loc_edit
                                .location
                                .end
                                .as_single()
                                .is_some_and(|m| m.is_certain());
                        // Only pure positional shuffle edits get the idempotent
                        // full-normalization treatment. These never carry an
                        // `ins[...]`/`con` cross-reference payload, so rewriting
                        // to `c.*N` and recursing cannot re-enter (and fail) the
                        // `try_expand_cds_ins` expansion the bounds gate is
                        // ordered to precede. Substitutions don't shuffle;
                        // ins/delins/dupins/con keep the plain early-return so a
                        // past-end input still rejects on the coordinate (strict
                        // mode -> InvalidCoordinates).
                        let is_shuffle_only_edit = matches!(
                            edit,
                            NaEdit::Deletion { .. }
                                | NaEdit::NPaddedDeletion { .. }
                                | NaEdit::Duplication { .. }
                                | NaEdit::Inversion { .. }
                        );
                        if both_certain
                            && is_shuffle_only_edit
                            && (new_start != *start_pos || new_end != *end_pos)
                        {
                            // The plain past-CDS coordinate maps into the 3'UTR:
                            // rewrite it to the in-bounds `c.*N` form and run
                            // FULL normalization on the result. Running the
                            // shuffle here (rather than returning the merely
                            // position-canonicalized variant) keeps the output
                            // idempotent — a shufflable `del`/`dup` in a 3'UTR
                            // repeat would otherwise render `c.*Ndel` on the
                            // first pass and then 3'-shift again on a second
                            // normalize pass. The rewritten position is in
                            // bounds, so this recursion does not re-enter the
                            // past-end gate. Re-attach the W4004 warning that
                            // flags the out-of-scheme input. See #920/#336.
                            let mut v = variant.clone();
                            v.loc_edit.location = CdsInterval::new(new_start, new_end);
                            let (normalized, mut warns) = self.normalize_cds(&v, manufactured)?;
                            let mut merged = bounds_warnings;
                            merged.append(&mut warns);
                            return Ok((normalized, merged));
                        }
                        // Non-shuffle edit (substitution, `=`, or an
                        // ins/delins/dupins/con that must not be expanded on a
                        // `*N` marker): still render the plain past-CDS
                        // coordinate in its canonical `c.*N` form when one exists
                        // (the PR's core feature), canonicalizing only the edit
                        // body. These edits do not 3'-shift, so no recursion is
                        // needed for idempotency. A genuinely-past-transcript or
                        // uncertain position leaves `out_variant == variant`.
                        let rewritten;
                        let out_variant =
                            if both_certain && (new_start != *start_pos || new_end != *end_pos) {
                                let mut v = variant.clone();
                                v.loc_edit.location = CdsInterval::new(new_start, new_end);
                                rewritten = v;
                                &rewritten
                            } else {
                                variant
                            };
                        return Ok((
                            HV::Cds(self.canonicalize_cds_variant(out_variant)),
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
            let (result, mut warnings) = self.normalize_cds(&new_variant, manufactured)?;
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

        // #1052: an uncertain/predicted-wrapped substitution must stay a
        // silent pass-through — see `is_uncertain_real_substitution`'s doc
        // comment. Must run before the intronic dispatch below too (an
        // uncertain intronic sub must stay silent regardless of the
        // intronic guard there).
        if is_uncertain_real_substitution(edit, &variant.loc_edit.edit) {
            return Ok((HV::Cds(variant.clone()), vec![]));
        }

        // #760: a UTR offset with no enclosing intron — past a transcript
        // terminus (`c.*824+10`, 3' of the final exon → `c.*834`) or on a
        // transcript whose model has no intron there at all — is not intronic.
        // Fold it to a plain (non-intronic) position and recurse so the standard
        // path maps it, instead of the intron-window shuffle that finds no
        // intron to anchor on. The folded position carries no offset, so the
        // recursion folds nothing and terminates.
        {
            use crate::hgvs::interval::UncertainBoundary;
            let mapper = crate::convert::CoordinateMapper::new(&transcript);
            let folded_start = mapper.fold_non_intronic_utr_offset(start_pos);
            let folded_end = mapper.fold_non_intronic_utr_offset(end_pos);
            if folded_start.is_some() || folded_end.is_some() {
                // Replace each folded boundary while preserving its certainty.
                let rebuild =
                    |boundary: &UncertainBoundary<CdsPos>, folded: Option<CdsPos>| match folded {
                        Some(p) => match boundary.as_single() {
                            Some(mu) => UncertainBoundary::Single(mu.map_ref(|_| p)),
                            None => UncertainBoundary::certain(p),
                        },
                        None => boundary.clone(),
                    };
                let new_location = Interval::with_complex_boundaries(
                    rebuild(&variant.loc_edit.location.start, folded_start),
                    rebuild(&variant.loc_edit.location.end, folded_end),
                );
                let new_variant = CdsVariant {
                    accession: variant.accession.clone(),
                    gene_symbol: variant.gene_symbol.clone(),
                    loc_edit: LocEdit::with_uncertainty(
                        new_location,
                        variant.loc_edit.edit.clone(),
                    ),
                };
                return self.normalize_cds(&new_variant, manufactured);
            }
        }

        // Handle intronic variants specially
        if start_pos.is_intronic() || end_pos.is_intronic() {
            // #1052: real substitutions are validated on the plain exonic path
            // only. Intronic-sub validation would require genomic projection
            // (an explicit non-goal); routing a sub through
            // `normalize_intronic_cds` / `normalize_boundary_spanning_cds` would
            // emit `ReducedCapabilityNoGenome` (no genomic data) or a hard
            // `ConversionError` (unmappable intron) on a variant that previously
            // passed through unchanged. Preserve that silent pass-through.
            if is_real_substitution(edit) {
                return Ok((HV::Cds(variant.clone()), vec![]));
            }
            // Switch to the accession-aware lookup so an NG/NC-parented input
            // gets the build-correct chromosome. If the accession-aware lookup
            // fails, fall back to the plain transcript we already fetched.
            let transcript = transcript_for_intronic().unwrap_or(transcript);
            // Both intronic *in the same intron* → the dedicated intronic
            // shuffle (which assumes a single enclosing intron). Both intronic
            // in DIFFERENT introns (#670: a deletion spanning intron→exon→intron)
            // falls through to the genomic boundary path, which sizes its window
            // to the full span — `normalize_intronic_cds` would bound the shuffle
            // to only `start`'s intron and so never shift such a span.
            if start_pos.is_intronic()
                && end_pos.is_intronic()
                && self.same_intron(&transcript, start_pos, end_pos)
            {
                return self.normalize_intronic_cds(variant, &transcript, start_pos, end_pos, edit);
            }
            // One endpoint exon/intron boundary, or a multi-intron span: genomic space.
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
        let mut boundaries = axis_info.clamped.clone();
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
        // **Whole-CDS del/dup exception (#918).** A deletion/duplication that
        // spans the ENTIRE CDS — one endpoint in the 5'UTR (`c.-N`), the other
        // in the 3'UTR (`c.*M`) — has a well-defined 3'-rule shuffle even though
        // its endpoints sit in different axes: the coding DNA reference is a
        // single contiguous spliced string (5'UTR+CDS+3'UTR; refseq.md), so the
        // deleted block simply rolls 3' within that string (e.g.
        // `NM_012459.2:c.-1_*1del` → `c.1_*2del`, matching mutalyzer). The axis
        // labels change but the sequence operation does not. Scoped to
        // 3'-direction del/dup; the result flows through `normalize_na_edit` over
        // the full transcript span below. A non-shifting span is left unchanged
        // (the CDS-start/-end clamps are gated on `start_axis == Cds`, which this
        // is not), so this only ever *adds* the spec-mandated shift.
        let is_whole_cds_span_delup = self.config.shuffle_direction == ShuffleDirection::ThreePrime
            && matches!(edit, NaEdit::Deletion { .. } | NaEdit::Duplication { .. })
            && matches!(start_axis, boundary::AxisRegion::FiveUtr)
            && matches!(end_axis, boundary::AxisRegion::ThreeUtr);
        // **Span-preserving re-typing exception (#1536).** The bail's stated
        // ground is that *the shuffle* has no defined semantics across an axis
        // boundary — and it then throws away the whole per-member pipeline,
        // canonicalization included. Re-typing a `delins` against the bases it
        // denotes (`canonicalize_delins`: `inv`, `del`, `ins`, `dup`, `sub`,
        // plus the shared-affix trim) is not a shuffle. It reads the reference
        // under the member's own span and moves nothing, and the `c.` reference
        // is a single contiguous string (`refseq.md`) whose axis labels change
        // at `cds_start`/`cds_end` while the sequence does not — the argument
        // #918 already makes verbatim for a whole-CDS del/dup.
        //
        // Discarding it is what made `c.9_*1delinsTGTGCATT` and `c.9_*1inv` two
        // fixed points for one variant, discriminated by nothing but where the
        // stop codon falls (#1536). Measured on a 40-mer with a real UTR at both
        // ends: over the 33 placements of one 8-base whole-block reverse
        // complement, the 14 that cross a CDS boundary all diverged and the 19
        // that do not all converged — including the three that sit wholly in a
        // UTR and the one that runs to the transcript end, so the discriminator
        // is the boundary and not the sequence end.
        //
        // `Delins` and `Inversion` together, not `Delins` alone. Confluence is a
        // property of the *pair*: the `delins` spelling reaches `inv` only if the
        // `inv` spelling reaches the same string, and the `inv` spelling is
        // subject to the same shared-affix trim (`c.3_10inv` -> `c.4_9inv`).
        // Letting one through and not the other trades one non-confluence for
        // another.
        //
        // Every other edit kind keeps the bail. For a `del`/`dup` the
        // canonicalization *is* the shuffle, so admitting them would be the
        // cross-axis shift #350 refuses and #918 carves out only for a
        // whole-CDS span; a `sub` cannot straddle a boundary; an `ins` is
        // already carved out by #402.
        let is_span_preserving_retype =
            matches!(edit, NaEdit::Delins { .. } | NaEdit::Inversion { .. });
        // Warnings raised before `normalize_na_edit` runs, so they survive the
        // carve-out below rather than being dropped with the early return.
        let mut pre_warnings: Vec<NormalizationWarning> = Vec::new();
        // Whether the span-preserving re-type carve-out below was actually
        // taken. Recorded explicitly rather than inferred from
        // `!pre_warnings.is_empty()`, even though the two agree today: the
        // recursion gate near the end of this function is a control-flow
        // decision about *this* path, and keying it on whether some vector
        // happens to be non-empty makes it fire for any future warning pushed
        // before `normalize_na_edit` — which the name `pre_warnings` positively
        // invites. That is the same mistake, one level up, as deciding which
        // warnings to drop by what the vector happened to hold; see the note on
        // the gate itself.
        let mut took_span_preserving_retype = false;
        if start_axis != end_axis
            && !matches!(start_axis, boundary::AxisRegion::None)
            && !matches!(end_axis, boundary::AxisRegion::None)
            && !is_zero_width_insertion
            && !is_whole_cds_span_delup
        {
            let acc = variant.accession.transcript_accession();
            let warning = NormalizationWarning::CrossAxisVariantNotShuffled {
                accession: acc,
                start_axis: start_axis.label().to_string(),
                end_axis: end_axis.label().to_string(),
            };
            if !is_span_preserving_retype {
                return Ok((
                    HV::Cds(self.canonicalize_cds_variant(variant)),
                    vec![warning],
                ));
            }
            // The warning still stands, and is still true: nothing below may
            // move this member off the footprint it arrived on. Pinning the
            // shuffle bound to the input span is what enforces that — the
            // re-typing runs, and a derived piece may settle anywhere inside the
            // member's own hull, but no shift can carry it into a region the
            // input did not already occupy. `axis_info.clamped` cannot serve
            // here: it is clamped to the *start* endpoint's region, which for a
            // straddling member is narrower than the member itself.
            pre_warnings.push(warning);
            took_span_preserving_retype = true;
            boundaries = Boundaries::new(tx_start.saturating_sub(1), tx_end);
        }

        // #918: relax the CDS↔3'UTR axis clamp for a CDS-resident del/dup so
        // the 3'-rule shuffle reaches its true most-3' position even when that
        // position lies past `cds_end`, in the 3'UTR (positive `c.<N>` →
        // `c.*<M>`). The HGVS 3' rule is unconditional — "for all descriptions,
        // the most 3' position possible **of the reference sequence** is
        // arbitrarily assigned to have been changed" (general.md; deletion.md
        // L20; duplication.md L24) — and its **only** stated exception is
        // deletions/duplications around exon/exon junctions (general.md;
        // numbering.md#DNAc), which is a genomic-projection concern, not the
        // CDS/UTR coordinate-label transition. The coding DNA reference
        // sequence is a single contiguous string (5'UTR+CDS+3'UTR; refseq.md),
        // so shifting from `c.250` into `c.*189` stays within that one sequence.
        // mutalyzer performs this shift (e.g. `c.250del` → `c.*189del`); ferro's
        // maintainers track the prior no-shift behavior as a bug (#487/#918),
        // not an accepted divergence.
        //
        // Scope (deliberately narrow — refuse anything the spec doesn't clearly
        // govern):
        //   - only when BOTH endpoints are CDS-resident (`start_axis` /
        //     `end_axis` == Cds) so we never reinterpret a variant that already
        //     straddles an axis (those keep the #350 bail above);
        //   - only del/dup/delins/ins (the spec-clear "most-3' of the reference
        //     sequence" shapes);
        //
        //     `Insertion` was originally excluded on the theory that an
        //     insertion's re-axing into the 3'UTR was already governed by the
        //     dedicated #387 CDS-end clamp. It is not, and the exclusion was
        //     itself a non-idempotency (#1209): with the bound left at
        //     `cds_end`, an insertion that could shift *past* the boundary
        //     instead stopped on it, and the #387 clamp then rewrote it to
        //     `c.<cds_end>delins<…>`. That output is a `Delins`, which DOES get
        //     this relaxation, so a second pass shifted it the rest of the way —
        //     the two passes literally took different branches of this `match`.
        //     `c.25_26insGAT` over a CDS ending `GGGG` with an `ACGT…` 3'UTR
        //     gave `c.26delinsGATG` then `c.*1_*2insTGA`.
        //
        //     Relaxing the bound does not disarm the #387 clamp: an insertion
        //     that genuinely *saturates* `cds_end` still has nowhere to shift,
        //     so the shuffle leaves it resting on the boundary and the clamp
        //     fires exactly as before. Only insertions that had somewhere to go
        //     stop reaching it — which is the whole defect.
        //
        //     `Delins` is included because gating on the input's edit-kind
        //     *spelling* made this relaxation miss its own case (#1185). A
        //     `delins` whose ref and alt share a prefix or suffix reduces to a
        //     pure deletion during normalization — `c.7_9delinsA` over ref `AAC`
        //     is a deletion of `AC` at `c.8_9` — and that deletion is exactly
        //     the shape this relaxation exists for. Clamped, it stopped at
        //     `cds_end` on the first pass and shifted only on a second, so
        //     normalization was not idempotent. The `r.` spelling of the
        //     byte-identical edit already shifted in one pass, which is what
        //     showed the clamp (not the shuffle) was at fault.
        //
        //     This also relaxes the bound for a `delins` that does NOT reduce,
        //     which is a slightly wider change than the non-idempotency strictly
        //     required. That is deliberate and consistent: the 3' rule is stated
        //     "for all descriptions", and the single-contiguous-sequence argument
        //     above does not depend on the edit kind. It is also what keeps the
        //     rule from depending on how a caller happened to spell an edit.
        //     The same reasoning is why `Insertion` now joins them rather than
        //     being special-cased to only the non-idempotent shape.
        //   - only the 3'-direction shuffle;
        //   - the bound is relaxed to the **exon** bound, not seq_len, so the
        //     spec's exon/exon exception stays enforced: a 3'UTR intron still
        //     stops the shuffle at the exon edge (the #670 machinery then
        //     governs any genomic-space continuation).
        // `tx_to_cds_pos` re-expresses any resulting `tx > cds_end` as `c.*<M>`.
        //     The `end_axis` test asks whether the variant's *reference span*
        //     straddles the boundary. An insertion has no reference span: its
        //     two endpoints are the flanks of a gap, so an insertion resting on
        //     the CDS/3'UTR junction necessarily reads `end_axis == ThreeUtr`
        //     while covering no 3'UTR base at all. Requiring both flanks
        //     therefore excluded every junction insertion from the relaxation
        //     above — including the ones this relaxation was extended to
        //     `Insertion` for — and reproduced #1209's own alternation one shape
        //     over (#1426):
        //
        //         c.18_*1insACT -> c.18delinsTACT -> c.*1_*2insCTA
        //
        //     The first pass is barred here, stops on the boundary and is
        //     rewritten by the #387 clamp; the resulting `Delins` has both
        //     flanks on `c.18`, so the second pass is admitted and shifts. Two
        //     passes, two branches — the exact defect #1209 records.
        //
        //     So for an insertion the question is asked of the 5' flank alone.
        //     A del/dup/delins still needs both, because for those the span is
        //     real and #350's bail is what keeps a genuinely straddling variant
        //     from being reinterpreted.
        //
        //     `is_zero_width_insertion`, not `matches!(edit, Insertion { .. })`.
        //     The argument above rests on an insertion covering no reference
        //     base, and that is what `tx_end == tx_start + 1` establishes. A
        //     malformed insertion with non-adjacent flanks (`c.10_20insA`) does
        //     name a span, and it is exactly the straddling shape #350's bail
        //     exists to refuse — so it must not be waved through on the strength
        //     of its edit kind. This is the same predicate the 5' mirror below
        //     uses, and the two arms should not disagree about what counts as an
        //     insertion for this purpose.
        if self.config.shuffle_direction == ShuffleDirection::ThreePrime
            && matches!(start_axis, boundary::AxisRegion::Cds)
            && (matches!(end_axis, boundary::AxisRegion::Cds) || is_zero_width_insertion)
            && matches!(
                edit,
                NaEdit::Deletion { .. }
                    | NaEdit::Duplication { .. }
                    | NaEdit::Delins { .. }
                    | NaEdit::Insertion { .. }
            )
        {
            boundaries = Boundaries::new(boundaries.left, exon_only.right);
        }

        // The 5' mirror, for zero-width insertions only (#1426).
        //
        // The relaxation above is 3'-direction-only, deliberately: the 3' rule
        // is what HGVS mandates. But the *bound* it relaxes is not a rule, it is
        // an axis clamp derived from the region the variant currently sits in —
        // and for a zero-width insertion that region flips as the insertion
        // moves, so the bound that stops it is one the next pass no longer
        // applies:
        //
        //     c.*1_*2insCTA  tx=(19,20) start_axis=ThreeUtr clamped.left=18
        //         -> stops at c.18_*1
        //     c.18_*1insACT  tx=(18,19) start_axis=Cds      clamped.left=0
        //         -> continues to c.17_18insTAC
        //
        // Two passes, two bounds, because the first pass's answer changed which
        // clamp the second pass reads. That is the same defect the 3' half of
        // #1426 was, arriving from the other side, and it is not a question of
        // whether the 5' shuffle *should* cross `cds_end` — it already does,
        // just one base per pass.
        //
        // Scoped to a zero-width insertion for the same reason the 3' arm asks
        // only its 5' flank: an insertion covers no reference base, so the axis
        // clamp's purpose — never reinterpret a variant whose span straddles the
        // boundary — has nothing to bite on. A `del`/`dup`/`delins` keeps both
        // bounds; for those the span is real.
        //
        // This reuses `is_zero_width_insertion` rather than re-deriving the same
        // predicate locally. The 3' arm's note above asks that the two arms not
        // disagree about what counts as an insertion here, and two bindings that
        // are equal today are exactly how that promise is broken later — a
        // narrowing applied to one would leave the other admitting a malformed
        // `c.10_20insA` on the strength of its edit kind alone.
        if self.config.shuffle_direction == ShuffleDirection::FivePrime
            && is_zero_width_insertion
            && matches!(
                start_axis,
                boundary::AxisRegion::Cds | boundary::AxisRegion::ThreeUtr
            )
        {
            boundaries = Boundaries::new(exon_only.left, boundaries.right);
        }

        // #918: a whole-CDS-spanning del/dup (5'UTR→3'UTR, carved out of the
        // #350 bail above) shuffles over the full contiguous spliced transcript
        // — its deleted block already crosses every internal CDS exon/intron
        // junction, so the single-exon `exon_only` bound does not apply; the
        // spliced sequence has no introns and the natural bound is the whole
        // transcript. The genomic-projection exon/exon exception is handled
        // downstream (#670), not by this spliced-coordinate bound.
        if is_whole_cds_span_delup {
            boundaries = Boundaries::new(0, transcript.sequence_length());
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
        // (`c.*N`) maps to `tx_start > cds_end`. A footprint touching UTR is
        // therefore not gated. If `cds_end` is unset we cannot verify the variant
        // lies within CDS proper, so we conservatively treat it as UTR-touching
        // and skip the gate — `CodonGate::for_input_span` owns both rules, and
        // the same call in `normalize_rna` is what keeps the two axes in step
        // (#469, #1192).
        let gate = CodonGate::for_input_span(
            transcript.cds_end.map(|cds_end| (cds_start, cds_end)),
            tx_start,
            tx_end,
        );
        let (mut new_tx_start, mut new_tx_end, mut new_edit, mut warnings) =
            self.normalize_na_edit(seq, edit, tx_start, tx_end, &boundaries, gate)?;

        // Substitutions are validated-then-returned-unchanged by
        // `normalize_na_edit`'s `Substitution` arm; nothing downstream (axis
        // clamps, the #670 junction-crossing shuffle continuation, canonical
        // split, `LocEdit` reconstruction) applies to a never-shuffled edit
        // kind, and letting a substitution fall through it anyway causes two
        // known regressions (#1052 follow-up): spurious genome-capability
        // rejections at exon boundaries, and silent `Mu::Uncertain` /
        // allele-member loss in cis alleles. Return the ORIGINAL variant
        // (not a `LocEdit::new`-reconstructed one) so uncertainty and all
        // other input structure survive untouched. Guard on
        // `is_real_substitution` (ref != alt) only — a `ref == alt` sub must
        // still flow through to the identity (`=`) rewrite below.
        if is_real_substitution(edit) {
            return Ok((HV::Cds(variant.clone()), warnings));
        }

        // #670: apply the 3' rule across the exon/intron junction. The
        // exon-confined shuffle above only sees spliced (exon) bases, so a
        // purely-exonic indel that comes to rest at an exon's 3' edge is never
        // given the chance to continue into the following intron — even though
        // the spec's "exception 3' rule" (numbering.md NOTE; deletion.md
        // exon/intron border) requires it. When the shuffle lands exactly at
        // the exon boundary, a downstream intron exists, and we have genomic
        // context, re-run the shuffle in genomic space (which spans the
        // junction naturally) and adopt the result only if it actually crossed
        // into the intron. The trigger is a rare edge landing, so the hot path
        // is untouched; it is 3'-only (5' shuffles stay exon-confined).
        // Perf note: an edge-landing variant whose genomic re-shuffle does NOT
        // cross still pays one `transcript_for_intronic()` lookup + one
        // `normalize_boundary_spanning_cds` fetch, then discards them. That cost
        // is confined to variants ending exactly at an exon's 3' base; if a bulk
        // 3'-direction `c.` workload regresses, this is the place to memoize.
        if self.config.shuffle_direction == ShuffleDirection::ThreePrime
            && new_tx_end == exon_only.right
        {
            if !self.provider.has_genomic_data() {
                // #1012: the junction-crossing 3'-shuffle is a genome-requiring
                // enhancement. With no genomic data we cannot even attempt it,
                // so the exon-confined result flows through as before — but mark
                // it degraded rather than silently returning a less-normalized
                // variant.
                warnings.push(NormalizationWarning::ReducedCapabilityNoGenome {
                    variant: format!("{variant}"),
                    capability: "exon/intron junction 3'-shuffle".to_string(),
                });
            } else {
                // Prefer the accession-aware transcript so an NG/NC-parented input
                // resolves the build-correct chromosome; fall back to the plain
                // transcript. Clone keeps `transcript` available for the fall-through.
                let boundary_transcript =
                    transcript_for_intronic().unwrap_or_else(|_| transcript.clone());
                if boundary_transcript.chromosome.is_some() {
                    // Engine errors (no following intron — e.g. last exon — no
                    // genomic alignment, …) fall through to the exon-confined
                    // result, the safe pre-#670 behavior. The exon/EXON suppression
                    // rule is preserved structurally: the genomic shuffle boundary
                    // is capped at the adjacent intron's far edge (never the next
                    // exon), so a homopolymer spanning an exon/exon junction can
                    // shift into the intron but never bridge into the downstream
                    // exon.
                    if let Ok((boundary_variant, boundary_warnings)) = self
                        .normalize_boundary_spanning_cds(
                            variant,
                            &boundary_transcript,
                            start_pos,
                            end_pos,
                            edit,
                        )
                    {
                        let crossed_into_intron = matches!(
                            &boundary_variant,
                            HV::Cds(cv)
                                if cv.loc_edit.location.start.inner().is_some_and(|p| p.is_intronic())
                                    || cv.loc_edit.location.end.inner().is_some_and(|p| p.is_intronic())
                        );
                        if crossed_into_intron {
                            // #1723: the offset in `boundary_variant` is ferro's
                            // own — an intronic input never reaches this gate —
                            // and `boundary_transcript.chromosome` is the contig
                            // the crossing was computed against. Record both here
                            // rather than guessing at the top of the pipeline.
                            record_manufactured_junction_exits(
                                &boundary_variant,
                                boundary_transcript.chromosome.as_deref(),
                                manufactured,
                            );
                            let mut combined = warnings;
                            combined.extend(boundary_warnings);
                            return Ok((boundary_variant, combined));
                        }
                    }
                }
            }
        }

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
                    self.normalize_na_edit(seq, edit, tx_start, tx_end, &exon_only, gate)?;
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
        //
        // Rewrite via an exact coordinate identity keyed on the POST-shuffle
        // edit (`new_edit`), not the input `edit`. After shuffling, an
        // insertion of the (rotated) sequence `A'` that has come to rest
        // immediately 5' of `cds_start` occupies flanks
        // `(cds_start-1, cds_start)` [`c.-1_1`] when a 5'UTR exists, or the
        // degenerate `(cds_start, cds_start)` when it does not (the 0-based-0
        // shuffle result clamps back to HGVS 1). "Insert `A'` between
        // `cds_start-1` and `cds_start`" is *identically* "delete
        // `ref[cds_start]`, insert `A' ++ ref[cds_start]`" — moving the delete
        // boundary one base, exact for ANY `A'`, no equivalence check needed.
        //
        // Keying on `new_edit` unifies three inputs that all reach here as the
        // same shuffled `A'`: a directly-written insertion, a codon-gated
        // `dup` (routed through the shuffle in `normalize_na_edit`), and any
        // start position — so every spelling of one haplotype collapses to a
        // single minimal, idempotent `c.<cds_start>delins<A' ++ ref[cds_start]>`
        // (confluence; supersedes the earlier width-varying `k`-widened form).
        let insertion_rests_at_cds_start = matches!(new_edit, NaEdit::Insertion { .. })
            && new_tx_end == cds_start
            && (new_tx_start == cds_start || new_tx_start + 1 == cds_start);
        if matches!(start_axis, boundary::AxisRegion::Cds)
            && !spanning_dup_exception
            && (insertion_rests_at_cds_start || new_tx_start < cds_start)
        {
            if matches!(edit, NaEdit::Delins { .. }) {
                // A Delins input whose canonicalisation pushed the residual
                // past the boundary restores to its own form (#383): the
                // shared-affix trim is what drove it past `cds_start`.
                new_edit = edit.clone();
                new_tx_start = tx_start;
                new_tx_end = tx_end;
            } else if insertion_rests_at_cds_start {
                if let NaEdit::Insertion {
                    sequence: InsertedSequence::Literal(a_prime),
                } = &new_edit
                {
                    if let Some(delins) = insertion_to_boundary_delins(
                        seq,
                        a_prime,
                        cds_start,
                        BoundarySide::FivePrime,
                    ) {
                        new_edit = delins;
                        new_tx_start = cds_start;
                        new_tx_end = cds_start;
                    }
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
        //
        // Rewrite via the mirror of the CDS-start identity, keyed on the
        // POST-shuffle `new_edit`. After shuffling, an insertion of the
        // (rotated) sequence `A'` resting immediately 3' of `cds_end` occupies
        // flanks `(cds_end, cds_end+1)` [`c.<cds_end>_*1`]. "Insert `A'` between
        // `cds_end` and `cds_end+1`" is *identically* "delete `ref[cds_end]`,
        // insert `ref[cds_end] ++ A'`" — exact for ANY `A'`. This keeps a
        // CDS-interior edit on the CDS axis instead of drifting onto the 3'UTR
        // (`c.*1`) axis, and (keyed on `new_edit`) unifies dup- and
        // insertion-spelled inputs and every start position into one minimal,
        // idempotent `c.<cds_end>delins<ref[cds_end] ++ A'>` (confluence;
        // resolves the #387 CDS-end saturation latent). Spanning *duplications*
        // (`new_edit` = Duplication, e.g. `c.9171_*1dup`) are the spec-canonical
        // form and are NOT clamped.
        if let Some(cds_end_tx) = transcript.cds_end {
            let insertion_rests_at_cds_end = matches!(new_edit, NaEdit::Insertion { .. })
                && new_tx_start == cds_end_tx
                && (new_tx_end == cds_end_tx || new_tx_end == cds_end_tx + 1);
            if matches!(start_axis, boundary::AxisRegion::Cds) && insertion_rests_at_cds_end {
                if matches!(edit, NaEdit::Delins { .. }) {
                    // Delins input restores to its own form (mirror #383/#387).
                    new_edit = edit.clone();
                    new_tx_start = tx_start;
                    new_tx_end = tx_end;
                } else if let NaEdit::Insertion {
                    sequence: InsertedSequence::Literal(a_prime),
                } = &new_edit
                {
                    if let Some(delins) = insertion_to_boundary_delins(
                        seq,
                        a_prime,
                        cds_end_tx,
                        BoundarySide::ThreePrime,
                    ) {
                        new_edit = delins;
                        new_tx_start = cds_end_tx;
                        new_tx_end = cds_end_tx;
                    }
                }
            }
        }

        // #1202: both clamps above are gated on `AxisRegion::Cds`, so a
        // **UTR-region** insertion is outside them and can still saturate
        // against the transcript bounds (reachable once `cds_start > 1`).
        //
        // Unconditional: the helper is self-gating on shape and position, and
        // either clamp above has already rewritten `new_edit` to a `Delins`, so
        // it is inert after them — no double-clamp. Runs before the CDS
        // conversion so the arithmetic stays in transcript space.
        clamp_insertion_at_sequence_bounds(
            seq,
            &mut new_edit,
            &mut new_tx_start,
            &mut new_tx_end,
            SequenceEnds::WHOLE,
        );

        // Convert back to CDS coordinates
        let new_start = self.tx_to_cds_pos(new_tx_start, cds_start, transcript.cds_end)?;
        let new_end = self.tx_to_cds_pos(new_tx_end, cds_start, transcript.cds_end)?;

        let new_variant = CdsVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                Interval::new(new_start, new_end),
                variant.loc_edit.edit.with_same_certainty(new_edit),
            ),
        };

        // #1536: the span-preserving carve-out above pinned the shuffle to the
        // input's own span, because a straddling member has no defined 3'-most
        // position. Re-typing can *resolve* that — the shared-affix trim moves
        // both endpoints inward, so `c.72_*3delins…` comes back as `c.*1_*4del`,
        // which straddles nothing. At that point the pin has no ground left and
        // the ordinary 3' rule applies again.
        //
        // Not cosmetic. Without this the pass emits an output that is not its
        // own fixed point: `c.*1_*4del` sits wholly in the 3'UTR, so normalizing
        // it a second time shuffles it on to `c.*2_*5del`. Measured by forcing the
        // gate below closed and re-running, on the branch as rebased:
        // `spec_conformance_axis`'s 3' `non_idempotent_outputs` reads **6** instead
        // of **4** — the two extra rows being exactly
        // `s00-c3{m,p}-cds-end-del-del-p2-sep1` — and confluence degrades with it,
        // 3' `split_three`/`split_more` reading 246/55 instead of 248/53. That is a
        // rank-1 conformance regression rather than a representation choice.
        //
        // (An earlier form of this comment quoted "rose from 7 to 9". That was
        // measured against the pre-#1599 baseline of 7; the baseline is now 4, and
        // 4 -> 6 is re-measured rather than the old figure rebased by arithmetic.)
        //
        // Recursion depth is one. The re-entry is gated on the new span sitting
        // in a single axis region, and a single-region member cannot reach the
        // carve-out; it is gated on the variant having actually changed, so a
        // fixed point cannot re-enter either.
        //
        // Only the warnings the recursion *invalidates* are dropped. Which ones
        // those are is decided by what each warning asserts, not by what happens
        // to be in the vector:
        //
        // - `pre_warnings` holds exactly one warning, pushed at exactly one site
        //   above — `CrossAxisVariantNotShuffled`. Invalidated: the member *is*
        //   shuffled, by the recursive call.
        // - `AxisClampApplied` is invalidated too. Its contract is to say *why
        //   the position did not shift further*, and on this path it did shift
        //   further (`c.*1_*4del` -> `c.*2_*5del`), so keeping it would misreport
        //   the result.
        // - **Everything else is carried forward**, because nothing about
        //   re-typing the span makes it untrue. `warnings` comes from
        //   `normalize_na_edit`, which validates the reference and pushes
        //   `RefSeqMismatch`; `NormalizeResult::has_ref_mismatch` is what strict
        //   mode rejects on. Dropping that one would silently stop strict
        //   `normalize()` refusing a straddling delins whose stated deleted bases
        //   do not match the reference. `ReducedCapabilityNoGenome` is the same
        //   shape — a statement about the *provider*, which a re-type cannot fix.
        //
        // **This is a correction, and the measurement that preceded it is why the
        // correction was needed.** An assertion at this line over the whole
        // `--lib --test it` suite fired on exactly one of 9 956 tests, carrying
        // only `AxisClampApplied` — and that was read as licence to drop the lot.
        // It is a claim about the corpus, not about the code: no corpus row builds
        // a straddling retype whose explicit deleted bases mismatch, so
        // `RefSeqMismatch` never reached the line and its loss was invisible.
        // Filter on the predicate, never on what the corpus happened to produce.
        //
        // Carried warnings are prepended without de-duplication. The re-typed edit
        // usually no longer states its deleted bases, so the recursive call cannot
        // re-derive the mismatch — which is exactly why carrying is necessary; and
        // `NormalizationWarning` is not `PartialEq`, while de-duplicating on
        // `code()` would conflate two `RefSeqMismatch` warnings at different
        // positions (position is the whole key — see the merge-time comment that
        // says so). A duplicated warning is strictly preferable to a lost one.
        // The gate reads the explicit `took_span_preserving_retype` flag rather
        // than `!pre_warnings.is_empty()`. The two are equivalent today, since
        // that vector has exactly one push site and it is the carve-out — but
        // equivalent-today is what this whole comment is about: the flag says
        // *the carve-out ran*, which is the fact the recursion depends on, while
        // emptiness is an observable that a second push site would silently
        // decouple from it.
        if took_span_preserving_retype && new_variant != *variant {
            let new_start_axis = boundary::axis_region_of(&transcript, new_tx_start);
            let new_end_axis = boundary::axis_region_of(&transcript, new_tx_end);
            if new_start_axis == new_end_axis
                && !matches!(new_start_axis, boundary::AxisRegion::None)
            {
                let mut carried: Vec<NormalizationWarning> = warnings
                    .into_iter()
                    .filter(|w| !matches!(w, NormalizationWarning::AxisClampApplied { .. }))
                    .collect();
                let (inner, inner_warnings) = self.normalize_cds(&new_variant, manufactured)?;
                carried.extend(inner_warnings);
                return Ok((inner, carried));
            }
        }
        // Prepended, not appended: `pre_warnings` describes a decision taken
        // before the edit was touched, so it reads first — and on the
        // early-return path above it was the only warning, which keeps the two
        // paths' ordering comparable.
        if !pre_warnings.is_empty() {
            pre_warnings.append(&mut warnings);
            warnings = pre_warnings;
        }

        // Issue #160 + #165 post-canonicalization split. The codon-frame
        // exception applies only to CDS-proper positions, which the
        // helper filters internally via `simple_cds_pos`.
        let uncertain = new_variant.loc_edit.edit.is_uncertain();
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Cds(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split, uncertain), warnings))
    }

    /// Normalize a transcript variant
    fn normalize_tx(
        &self,
        variant: &TxVariant,
        manufactured: &mut Vec<ManufacturedJunctionExit>,
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
            return self.normalize_tx(&new_variant, manufactured);
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
            || -> Result<std::sync::Arc<crate::reference::transcript::Transcript>, FerroError> {
                // Resolve directly from the accession — avoids cloning the
                // whole variant just to call the by-variant lookup.
                self.provider
                    .get_transcript_for_accession(&variant.accession)
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
            let (result, mut warnings) = self.normalize_tx(&new_variant, manufactured)?;
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

        // #1052: an uncertain/predicted-wrapped substitution must stay a
        // silent pass-through — see `is_uncertain_real_substitution`'s doc
        // comment. Must run before the intronic dispatch below too.
        if is_uncertain_real_substitution(edit, &variant.loc_edit.edit) {
            return Ok((HV::Tx(variant.clone()), vec![]));
        }

        if start_pos.is_intronic() || end_pos.is_intronic() {
            // #1052: real substitutions are validated on the plain exonic path
            // only. Intronic-sub validation would require genomic projection
            // (an explicit non-goal); routing a sub through
            // `normalize_intronic_tx` / `normalize_boundary_spanning_tx` would
            // emit `ReducedCapabilityNoGenome` (no genomic data) or a hard
            // `ConversionError` (unmappable intron) on a variant that previously
            // passed through unchanged. Preserve that silent pass-through.
            if is_real_substitution(edit) {
                return Ok((HV::Tx(variant.clone()), vec![]));
            }
            // Switch to the accession-aware lookup so an NG/NC-parented input
            // gets the build-correct chromosome.
            let transcript = transcript_for_intronic().unwrap_or(transcript);
            // Mirror of the `c.` dispatch (#670/#704): both intronic *in the
            // same intron* → the dedicated single-intron shuffle; one endpoint
            // on an exon/intron boundary, or a multi-intron span → genomic
            // boundary-spanning space (`normalize_intronic_tx` would bound the
            // shuffle to only `start`'s intron and never shift such a span).
            if start_pos.is_intronic()
                && end_pos.is_intronic()
                && self.same_intron_tx(&transcript, start_pos, end_pos)
            {
                return self.normalize_intronic_tx(variant, &transcript, start_pos, end_pos, edit);
            }
            return self.normalize_boundary_spanning_tx(
                variant,
                &transcript,
                start_pos,
                end_pos,
                edit,
            );
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
        // `mut` bindings are #1202's (its transcript-bound clamp rewrites these
        // in place). `CodonGate::NotApplicable`: `n.` has no
        // reading frame, so there is no CDS span to re-check a shuffled tract
        // against — unlike `r.` below, which does pass real bounds.
        let (mut new_start, mut new_end, mut new_edit, mut warnings) = self.normalize_na_edit(
            seq,
            edit,
            tx_start,
            tx_end,
            &boundaries,
            CodonGate::NotApplicable,
        )?;

        // See the identical guard + comment in `normalize_cds` (#1052 follow-up):
        // substitutions are validated-then-returned-unchanged by
        // `normalize_na_edit`; nothing downstream applies to a never-shuffled
        // edit kind, so return the ORIGINAL variant to preserve `Mu::Uncertain`
        // and all other input structure. `ref == alt` still flows through to
        // the identity (`=`) rewrite below.
        if is_real_substitution(edit) {
            return Ok((HV::Tx(variant.clone()), warnings));
        }

        // #704 sub-problem A (mirror of the `normalize_cds` block, #670): apply
        // the 3' rule across the exon/intron junction for `n.`. The exon-confined
        // shuffle above only sees spliced (exon) bases, so a purely-exonic del/dup
        // that comes to rest at an exon's 3' edge is never given the chance to
        // continue into the following intron — even though the spec's exception
        // 3' rule requires it. When the shuffle lands exactly at the exon
        // boundary (`boundaries.right`, the del/dup exon bound), a downstream
        // intron exists, and we have genomic context, re-run the shuffle in
        // genomic space (which spans the junction naturally) and adopt the result
        // only if it actually crossed into the intron. The trigger is a rare edge
        // landing, so the hot path is untouched; it is 3'-only (5' shuffles
        // stay exon-confined), exactly as the `c.` path.
        //
        // Ordering constraint (#1202, other half documented at the clamp
        // below): a 3'-saturated insertion arrives here with
        // `new_end == boundaries.right + 1`, so it does NOT trip this trigger —
        // and the clamp below then moves it onto exactly `boundaries.right`.
        // Neither block may be reordered past the other.
        if self.config.shuffle_direction == ShuffleDirection::ThreePrime
            && new_end == boundaries.right
        {
            if !self.provider.has_genomic_data() {
                // #1012: mirror of the `normalize_cds` reduced-capability guard.
                // The junction-crossing 3'-shuffle needs a genome; without one
                // the exon-confined result flows through, marked degraded.
                warnings.push(NormalizationWarning::ReducedCapabilityNoGenome {
                    variant: format!("{variant}"),
                    capability: "exon/intron junction 3'-shuffle".to_string(),
                });
            } else {
                // Prefer the accession-aware transcript so an NG/NC-parented input
                // resolves the build-correct chromosome; fall back to the plain
                // transcript. Clone keeps `transcript` available for the fall-through.
                let boundary_transcript =
                    transcript_for_intronic().unwrap_or_else(|_| transcript.clone());
                if boundary_transcript.chromosome.is_some() {
                    // Engine errors (no following intron — e.g. last exon — no genomic
                    // alignment, …) fall through to the exon-confined result, the safe
                    // pre-#704 behavior. The exon/EXON suppression rule is preserved
                    // structurally: the genomic shuffle window is capped at the
                    // adjacent intron's far edge (never the next exon).
                    if let Ok((boundary_variant, boundary_warnings)) = self
                        .normalize_boundary_spanning_tx(
                            variant,
                            &boundary_transcript,
                            start_pos,
                            end_pos,
                            edit,
                        )
                    {
                        let crossed_into_intron = matches!(
                            &boundary_variant,
                            HV::Tx(tv)
                                if tv.loc_edit.location.start.inner().is_some_and(|p| p.is_intronic())
                                    || tv.loc_edit.location.end.inner().is_some_and(|p| p.is_intronic())
                        );
                        if crossed_into_intron {
                            // #1723 — mirror of the `normalize_cds` recording.
                            record_manufactured_junction_exits(
                                &boundary_variant,
                                boundary_transcript.chromosome.as_deref(),
                                manufactured,
                            );
                            let mut combined = warnings;
                            combined.extend(boundary_warnings);
                            return Ok((boundary_variant, combined));
                        }
                    }
                }
            }
        }

        // #1202: rewrite an insertion the shuffle drove off either end of the
        // transcript; see the helper for the identity and its proof.
        //
        // Must stay AFTER the exon/intron junction block above: that block
        // tests `new_end == boundaries.right`, and the 3' rewrite moves
        // `new_end` from `len + 1` onto exactly `len` (== `boundaries.right`
        // for the insertion bounds), so clamping first would divert a saturated
        // insertion into the junction-crossing engine.
        clamp_insertion_at_sequence_bounds(
            seq,
            &mut new_edit,
            &mut new_start,
            &mut new_end,
            SequenceEnds::WHOLE,
        );

        let new_variant = TxVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                Interval::new(TxPos::new(new_start as i64), TxPos::new(new_end as i64)),
                variant.loc_edit.edit.with_same_certainty(new_edit),
            ),
        };

        // Issue #160 + #165 post-canonicalization split.
        let uncertain = new_variant.loc_edit.edit.is_uncertain();
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Tx(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split, uncertain), warnings))
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

        // #972 Task 5: mirrors the `normalize_cds` gate. A `p.` variant
        // whose underlying transcript has an incomplete (unconfirmed ATG)
        // 5' CDS start (`cds_start_NF`) has no defined `p.1`, so decline to
        // validate/re-number — pass through verbatim in every mode — before
        // even attempting the reference-amino-acid check below (that check
        // itself assumes a trustworthy p.-numbering). Best-effort: this
        // resolves the transcript via `accession.transcript_accession()`,
        // which succeeds for a genomic-context-wrapped coding selector
        // (`NG_(NM_):p.…`) but not for a bare NP_/XP_ protein accession with
        // no nucleotide counterpart in the reference data — that shape
        // simply falls through unchecked, same as other reference-dependent
        // checks in this module when the transcript can't be resolved.
        if let Ok(transcript) = self
            .provider
            .get_transcript(&variant.accession.transcript_accession())
        {
            if transcript.cds_start_incomplete {
                let mut warnings = Vec::new();
                if self.config.should_warn_incomplete_cds_start()
                    || self.config.should_reject_incomplete_cds_start()
                {
                    warnings.push(NormalizationWarning::IncompleteCdsStartReference {
                        accession: variant.accession.transcript_accession(),
                        coordinate_system: "p".to_string(),
                    });
                }
                return Ok((HV::Protein(variant.clone()), warnings));
            }
        }

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
            // Only a literal inserted sequence can duplicate upstream residues.
            // `insXaa[n]` / `insTer<n>` have no concrete residues to compare, so
            // they pass through unchanged via the catch-all arm below.
            ProteinEdit::Insertion {
                sequence: ProteinInsSeq::Literal(seq),
            } => self
                .try_protein_ins_to_dup(variant, seq)
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
        //
        // Gated on the error configuration (#1196) — the push used to be
        // unconditional, so `--error-mode silent` and `--ignore W3022` were
        // both inert for this code. See
        // `NormalizeConfig::should_warn_initiator_met_canonicalization` for
        // why an advisory about ferro's own canonical output is surfaced
        // rather than promoted to a rejection in strict mode.
        let mut warnings: Vec<NormalizationWarning> = Vec::new();
        if self.config.should_warn_initiator_met_canonicalization()
            && matches!(output_edit, ProteinEdit::Duplication)
        {
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

        // #502: a `p.` coordinate on a genomic reference whose parenthetical
        // selector is a coding transcript (`NM_`/`XM_`) should carry the
        // protein accession (`NP_`/`XP_`) — the spec-preferred selector for a
        // protein-level statement (background/refseq.md L38-42). Rewrite the
        // selector when the reference data resolves the paired protein;
        // otherwise the input selector is preserved (#121).
        let final_variant = match self.resolve_genomic_protein_selector(&final_variant.accession) {
            Some(protein_accession) => ProteinVariant {
                accession: protein_accession,
                ..final_variant
            },
            None => final_variant,
        };

        // NOTE: the protein-axis split move (`protein/delins.md:21`, `:64`) is
        // deliberately NOT applied here. It turns one member into several, so it
        // runs at the whole-description seam —
        // [`Self::split_protein_separation`], called from
        // `normalize_core_canonical` — which this method cannot be, because
        // `normalize_allele` re-enters it once per member. That helper records
        // what wiring it in here was measured to produce.
        Ok((HV::Protein(final_variant), warnings))
    }

    /// Resolve the spec-preferred protein-accession (`NP_`/`XP_`) selector for a
    /// `p.` variant on a genomic reference (#502).
    ///
    /// For an input like `NG_012337.3(NM_003002.4):p.(Asp92Tyr)` the
    /// parenthetical selector is a coding *transcript* (`NM_`), but a `p.`
    /// coordinate is a protein-level statement; per HGVS
    /// (`background/refseq.md` L38-42) a protein accession (`NP_`) is an
    /// accepted parenthetical specification for a protein coordinate (the
    /// spec-preferred selector here). When the reference
    /// data resolves the transcript's paired protein, this returns that
    /// `NP_`/`XP_` accession with the genomic-context wrapper preserved, to
    /// replace the `NM_`/`XM_` selector.
    ///
    /// Returns `None` (leave the selector unchanged) when there is no genomic
    /// context, the selector is not a coding transcript, the provider has no
    /// such transcript, or it carries no `protein_id`. This is never an error:
    /// a genome-less or data-poor environment (e.g. `MockProvider`) simply
    /// preserves the input `NM_` selector, mirroring the #121 "preserve when
    /// present, don't synthesize" policy.
    fn resolve_genomic_protein_selector(
        &self,
        accession: &crate::hgvs::variant::Accession,
    ) -> Option<crate::hgvs::variant::Accession> {
        use crate::hgvs::parser::accession::parse_accession;

        // Only a genomic-reference compound accession (NG_/LRG_/NC_ wrapper) is
        // in scope; a bare transcript selector is left untouched.
        let context = accession.genomic_context.as_ref()?;
        // Only a coding-transcript selector (NM_/XM_) has a paired NP_/XP_.
        if !matches!(accession.prefix.as_ref(), "NM" | "XM") {
            return None;
        }
        // Resolve the transcript's paired protein accession from reference data.
        let transcript = self
            .provider
            .get_transcript(&accession.transcript_accession())
            .ok()?;
        let protein_id = transcript.protein_id.as_deref()?;
        // Parse the protein accession and re-attach the genomic-context wrapper.
        let (rest, protein_accession) = parse_accession(protein_id).ok()?;
        if !rest.is_empty() {
            // Malformed `protein_id`: preserve the input selector rather than
            // emit a partially-parsed accession.
            return None;
        }
        // The provider's `protein_id` must be an actual protein accession
        // (`NP_`/`XP_`). A non-protein value (e.g. a transcript or gene
        // accession) is semantically malformed for a protein selector, so
        // preserve the input selector rather than emit it as a `p.` selector.
        if !matches!(protein_accession.prefix.as_ref(), "NP" | "XP") {
            return None;
        }
        Some(protein_accession.with_genomic_context((**context).clone()))
    }

    /// Rewrite a legacy LOVD gene-model selector on a genomic reference —
    /// `NG_/NC_/LRG(GENE[_v001]):c.…` — to the spec-preferred transcript-accession
    /// form `NG_/NC_/LRG(NM_):c.…`, when the provider resolves the gene's
    /// transcript (#500/#637). `_v001` and the bare gene name resolve identically;
    /// the `c.` coordinates are unchanged (the gene-model selector *is* that
    /// transcript).
    ///
    /// Resolution is parent-relative-first (#792): the `NG_` parent accession is
    /// passed through (`Some(accession)` below), so when that exact `NG_` version
    /// uniquely hosts the gene, its hosted transcript wins; only otherwise does
    /// resolution fall back to the global reference-standard map.
    ///
    /// With no gene-symbol selector, synthesizes the transcript the `NG_`
    /// uniquely hosts via [`ReferenceProvider::sole_hosted_transcript`] (#923),
    /// returning `Some(NG_(NM_))`; it returns `None` there only when that hosted
    /// mapping is absent or ambiguous. Also returns `None` when the reference is
    /// not genomic (`NM_(GENE)` keeps the #121 behavior), the selector already
    /// names a transcript context, or a present gene-symbol selector is
    /// unresolvable (unknown, higher locus version, or no summary ingested) —
    /// mirroring the #121 "preserve when present, don't synthesize" policy for
    /// the gene-symbol case. Never an error.
    fn rewrite_legacy_gene_selector(
        &self,
        variant: &crate::hgvs::variant::HgvsVariant,
    ) -> Option<crate::hgvs::variant::HgvsVariant> {
        use crate::hgvs::parser::accession::parse_accession;
        use crate::hgvs::variant::{CdsVariant, HgvsVariant};

        // Only the coding axis carries the corpus's legacy selectors; the
        // reference-standard map resolves to `NM_`, which suits `c.` coordinates.
        let HgvsVariant::Cds(c) = variant else {
            return None;
        };
        let accession = &c.accession;
        // Only a genomic-reference selector slot (`NG_`/`NC_`/`LRG_`) is in scope.
        let is_genomic_ref = matches!(accession.prefix.as_ref(), "NG" | "NC") || accession.is_lrg();
        // A selector that already names a transcript context is not a bare
        // gene-model selector — leave it.
        if !is_genomic_ref || accession.genomic_context.is_some() {
            return None;
        }
        // With a gene-symbol selector, resolve it via the legacy gene model
        // (#500/#792). With no selector at all, synthesize the transcript the
        // `NG_` parent uniquely hosts (#923) — declining when ambiguous.
        let nm = match c.gene_symbol.as_deref() {
            Some(gene_symbol) => self
                .provider
                .resolve_legacy_gene_selector(gene_symbol, Some(accession))?,
            None => self.provider.sole_hosted_transcript(accession)?,
        };
        let (rest, nm_accession) = parse_accession(&nm).ok()?;
        // The resolver yields a reference-standard `NM_` coding transcript
        // (`parse_refseqgene_summary` keeps only `NM_` rows); reject anything else
        // rather than emit a malformed or non-coding selector for a `c.` variant.
        if !rest.is_empty() || nm_accession.prefix.as_ref() != "NM" {
            return None;
        }
        // The genomic reference becomes the context wrapper; the resolved
        // transcript becomes the selector; the now-redundant gene symbol is dropped.
        Some(HgvsVariant::Cds(CdsVariant {
            accession: nm_accession.with_genomic_context(accession.clone()),
            gene_symbol: None,
            loc_edit: c.loc_edit.clone(),
        }))
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
    /// delins with no shared affix, uncertain endpoints, or a
    /// ≥3-residue range with no provider data — the ≤2-residue
    /// cases are trimmed reference-free from the named endpoints).
    ///
    /// Branches (per HGVS Prioritization `general.md:56-57`):
    /// - residual empty + del empty → identity `p.<X>_<Y>=`
    /// - residual empty, del non-empty → pure `del`
    /// - residual single AA, del single AA → substitution
    /// - residual non-empty, del empty → route through
    ///   `try_protein_ins_to_dup` (anchored at the flanking
    ///   positions of the trimmed range); reference-free, emit the
    ///   plain `<X>_<Y>ins<residual>` when both flanking residues are
    ///   named by the input's endpoints (#1126)
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
        let expected_len = (end_pos.number - start_pos.number + 1) as usize;

        // Obtain the deleted reference residues. Prefer the provider's protein
        // sequence (works for any range length); failing that, fall back to the
        // residues NAMED in the `delins` location endpoints. That fallback is
        // valid only when *every* deleted residue is named — i.e. the range
        // spans ≤ 2 residues (start and/or end). A ≥ 3-residue range has
        // unnamed interior residues and cannot be trimmed reference-free
        // (#1119). This mirrors the nucleotide axis's reference-driven
        // `canonicalize_delins`, restricted here to the AST-derivable cases.
        //
        // `has_protein_data()` is a **provider-wide** flag, so a provider that
        // carries proteins for other accessions still fails to fetch the one at
        // hand. That is "no reference for this accession", not "cannot
        // canonicalize": take the same named-endpoint fallback rather than
        // declining outright (#1131).
        let fetched: Option<Vec<AminoAcid>> = if self.provider.has_protein_data() {
            let accession = variant.accession.transcript_accession();
            self.fetch_protein_window(&accession, start_pos.number - 1, end_pos.number)
                .filter(|aas| aas.len() == expected_len)
        } else {
            None
        };
        // Whether `ref_aas` came from the reference. Gates the branches below
        // that need more of the protein than the location names — the ins→dup
        // search and its flanking-residue lookups.
        let reference_backed = fetched.is_some();
        let ref_aas: Vec<AminoAcid> = match fetched {
            Some(aas) => aas,
            None => match expected_len {
                1 => vec![start_pos.aa],
                2 => vec![start_pos.aa, end_pos.aa],
                _ => return None,
            },
        };

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

        // No affix trimmed. A 1→1 delins is still a substitution (sub > delins
        // priority, `delins.md`) even with no shared affix, so let that single
        // case fall through to the reduction below; otherwise there is no
        // rewrite to make. Do NOT route the non-1→1 case through an ins→dup
        // search here: the deleted window is still non-empty, so any
        // tandem-match rewrite would emit a bare `dup` that silently drops the
        // deletion side of the edit and returns a non-equivalent variant.
        if lcp == 0 && lcs == 0 && !(residual_del.len() == 1 && residual_seq.0.len() == 1) {
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
            // Trimmed to a pure insertion. The ins→dup search below needs the
            // protein sequence; without it, emit the minimal insertion when —
            // and only when — the AST names the residue on BOTH sides of the
            // insertion point (#1126).
            //
            // `ref_aas[k]` is the residue at 1-based position
            // `start_pos.number + k`, so the insertion point's flanks —
            // `new_start - 1` and `new_start`, i.e. `start_pos.number + lcp - 1`
            // and `start_pos.number + lcp` — are `ref_aas[lcp - 1]` and
            // `ref_aas[lcp]`. Both are in range exactly when
            // `1 <= lcp < ref_aas.len()`: a leading trim that consumed nothing
            // puts the point 5′ of the named window, one that consumed all of it
            // puts the point 3′ of it, and in either case the missing flank
            // residue is only knowable from the reference — so keep the input
            // delins (#1119's conservative default).
            //
            // The ins→dup refinement stays reference-gated: when the reference
            // really backed `ref_aas`, it may reduce this same insertion further
            // to a `dup`, which only the reference can establish. Gate on
            // `reference_backed`, not the provider-wide capability flag — for an
            // accession the provider lacks, the reference path can only bail
            // (#1131).
            if !reference_backed {
                if lcp == 0 || lcp >= ref_aas.len() {
                    return None;
                }
                let ins_start = ProtPos::new(ref_aas[lcp - 1], new_start - 1);
                let ins_end = ProtPos::new(ref_aas[lcp], new_start);
                return Some((
                    ProteinEdit::Insertion {
                        sequence: ProteinInsSeq::Literal(residual_seq),
                    },
                    Interval::new(ins_start, ins_end),
                ));
            }
            let accession = variant.accession.transcript_accession();
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
                            sequence: ProteinInsSeq::Literal(residual_seq.clone()),
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
                    sequence: ProteinInsSeq::Literal(residual_seq),
                },
                Interval::new(ins_start, ins_end),
            ));
        }
        if residual_del.len() == 1 && residual_seq.0.len() == 1 {
            // Sub > delins — except at the translation initiation codon, where
            // the substitution form does not exist.
            //
            // `validate_no_start_loss_substitution` (parser) refuses
            // `p.Met1<Xaa>`: a variant that changes the initiation codon is
            // described neither as a substitution nor as an extension
            // (`protein/substitution.md:49`, `checklist.md:65`,
            // `protein/extension.md:28`). A *range* delins is not refused —
            // that guard keys on a single certain position — so an input like
            // `p.Met1_Ala3delinsValAlaAla` parses, and it is the affix trim
            // here that would manufacture the illegal `p.Met1Val`. The guard
            // therefore belongs in the producer as well (#1607).
            //
            // Declining is the only move available. The spec offers `p.0`,
            // `p.0?`, `p.(Met1?)` and — when an upstream initiation site is
            // activated — the insertion form `p.Met1_Leu2ins…`; which one holds
            // is a claim about the *consequence*, and the reference and payload
            // residues do not settle it. So the input is left as authored
            // rather than canonicalised into a choice ferro cannot make.
            //
            // The condition mirrors **`validate_no_start_loss_substitution`**
            // — the parser's guard, and only that one — evaluated against the
            // substitution this branch is about to build: residue 1, an
            // initiator `Met` as the reference residue, and a different
            // residue replacing it. So the producer cannot emit a
            // substitution that entry point would refuse. A non-`Met` residue
            // 1 and a `Met`→`Met` residual are both accepted there and so are
            // left to canonicalise.
            //
            // It is deliberately **not** parity with the sibling producer
            // guard in `try_protein_delins_split_unchanged_interior`
            // (#1606), which declines on the broader
            // `start_pos.number == 1 && ref_aas[0] != seq.0[0]` — any change
            // at residue 1, `Met` or not. That path emits an allele of
            // members and takes the wider berth; this one is scoped to the
            // reparse hole, so it tracks the parser instead.
            if new_start == 1
                && residual_del[0] == AminoAcid::Met
                && residual_seq.0[0] != AminoAcid::Met
            {
                return None;
            }
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

    /// The protein-axis split move: an **equal-length** `delins` whose interior
    /// leaves at least one residue unchanged describes two or more separate
    /// changes, not one `delins`.
    ///
    /// `recommendations/protein/delins.md:21` states it outright — "two variants
    /// separated by one or more amino acids should be described individually and
    /// not as a `delins`" — and `:64` publishes the worked answer:
    /// `p.[Ser44Arg;Trp46Arg]`, with the explicit note that "the variant is
    /// **not** described as `p.Ser44_Trp46delinsArgLeuArg`". This is the protein
    /// analogue of the nucleotide separation rule (`general.md:34`) that
    /// [`Self::apply_canonical_split`] applies via
    /// [`rules::decompose_delins`](crate::normalize::rules::decompose_delins).
    ///
    /// # Why equal length, and only equal length
    ///
    /// A `delins` whose reference span and payload have the **same** length has
    /// exactly one residue-wise correspondence, so "the middle residue is
    /// unchanged" is a fact about the two sequences rather than a choice. An
    /// unequal-length `delins` has no such correspondence: finding an unchanged
    /// run inside it means first picking an alignment, and *which* alignment is
    /// not determined by the reference and the payload — several are equally
    /// good and they disagree about where the unchanged residue sits. On the
    /// nucleotide axis that choice is settled by applying the edits to the
    /// reference and re-deriving from the resulting sequence
    /// (`canonical-form-choice-when-both-legal`), which is not available here:
    /// there is no apply-to-reference on the protein axis, which is why
    /// [`merge::cis_kind_of`](crate::normalize::merge) returns `None` for
    /// `Protein`. With nothing to derive the answer from, ferro declines rather
    /// than invent an alignment, and those are left as authored.
    ///
    /// The floor is **one** unchanged residue: `delins.md:21` says "separated by
    /// one or more amino acids", and unlike the nucleotide axis there is no
    /// codon carve-out here (`general.md:35-38` is a *reading-frame* exception,
    /// and a protein description has no codons of its own to share).
    ///
    /// # Why this is not routed through the nucleotide canonicalizer
    ///
    /// [`merge::canonicalize_from_sequence`](crate::normalize::merge) applies the
    /// edit set to a reference window and re-derives the partition from the
    /// resulting sequence. There is no apply-to-reference on the protein axis —
    /// [`merge::cis_kind_of`](crate::normalize::merge) returns `None` for
    /// `Protein` for exactly that reason — so the move is implemented here, on
    /// the reference residues themselves.
    ///
    /// # When it declines
    ///
    /// Returning `None` leaves the input `delins` exactly as authored. It
    /// declines when:
    ///
    /// - the edit is not a `Delins`, or either endpoint is uncertain or a range
    ///   (there is then no residue-wise correspondence to read);
    /// - the payload length differs from the span, per the section above;
    /// - the span is under three residues, which has no interior to be
    ///   unchanged;
    /// - **the protein reference is unavailable for this accession.** The
    ///   unchanged middle residue is a claim about the *reference*, and the
    ///   payload cannot testify to it: assuming `payload[i]` matches the
    ///   reference is the guess this rule must never make. `has_protein_data()`
    ///   is provider-wide, so the fetch itself is what decides (#1131);
    /// - the reference span or the payload names an unknown residue (`Xaa`),
    ///   where equality is not a statement about identity;
    /// - **either side carries a `Ter`.** In the *payload*, a split across a stop
    ///   would emit members 3' of it, which `delins.md:45` forbids ("amino acids
    ///   after the translation termination codon are **not** listed"). In the
    ///   *reference span*, a stop replaced by an ordinary residue is a
    ///   **stop-loss**, and `extension.md:18` ranks it — "prioritisation: (1)
    ///   extension, (2) frameshift or deletion-insertion" — so it is described
    ///   with `ProteinEdit::Extension` (`:30`, `p.Ter110GlnextTer17`), never a
    ///   substitution. Splitting `Ser-Leu-Ter` against a payload of
    ///   `Ala-Leu-His` would emit `p.[Ser988Ala;Ter990His]`, which spells an
    ///   extension as a substitution and states the wrong consequence. Both are
    ///   the mirror of the `first_ter_pos` gate in
    ///   [`merge::coalesce_protein_adjacent_changes`](crate::normalize::merge);
    /// - **residue 1 is the changed one.** A member at the translation
    ///   initiation codon renders as `p.Met1Xxx`, which `substitution.md:49`
    ///   forbids ("not described as a substitution or as an extension") and
    ///   which `parse_hgvs` refuses, so the split would emit a description ferro
    ///   cannot read back. The correct description — `p.0`, `p.0?`,
    ///   `p.(Met1?)`, or the insertion form when an upstream initiation site is
    ///   activated — is a claim about the *consequence*, and two protein
    ///   sequences do not settle which. The range input parses only because
    ///   `validate_no_start_loss_substitution` keys on a single certain
    ///   position; the illegality is created by the split, so the guard belongs
    ///   here;
    /// - fewer than two changed runs survive, i.e. the change really is one
    ///   contiguous `delins`.
    fn try_protein_split_delins(
        &self,
        variant: &crate::hgvs::variant::ProteinVariant,
    ) -> Option<Vec<crate::hgvs::variant::ProteinVariant>> {
        use crate::hgvs::edit::AminoAcidSeq;
        use crate::hgvs::variant::{LocEdit, ProteinVariant};

        let seq = match variant.loc_edit.edit.inner()? {
            ProteinEdit::Delins { sequence } => sequence,
            _ => return None,
        };
        // Certain single points only — same rationale as
        // `try_protein_delins_canonicalize`.
        let start_pos = match variant.loc_edit.location.start.as_single()? {
            Mu::Certain(p) => *p,
            _ => return None,
        };
        let end_pos = match variant.loc_edit.location.end.as_single()? {
            Mu::Certain(p) => *p,
            _ => return None,
        };
        if start_pos.number == 0 || end_pos.number < start_pos.number {
            return None;
        }
        let span = (end_pos.number - start_pos.number + 1) as usize;
        // Equal length, and at least one interior residue to be unchanged.
        if span < 3 || span != seq.0.len() {
            return None;
        }
        if seq.0.contains(&AminoAcid::Ter) || seq.0.contains(&AminoAcid::Xaa) {
            return None;
        }

        // The reference residues, or nothing. Never the payload's own middle.
        if !self.provider.has_protein_data() {
            return None;
        }
        let accession = variant.accession.transcript_accession();
        let ref_aas =
            self.fetch_protein_window(&accession, start_pos.number - 1, end_pos.number)?;
        if ref_aas.len() != span
            || ref_aas.contains(&AminoAcid::Xaa)
            || ref_aas.contains(&AminoAcid::Ter)
        {
            return None;
        }

        // The translation initiation codon, changed. A member here would be
        // `p.Met1Xxx`, which `substitution.md:49` forbids and `parse_hgvs`
        // refuses — so the split would emit a description ferro cannot read
        // back. Which description is correct (`p.0`, `p.0?`, `p.(Met1?)`, or
        // the upstream-initiation insertion form) is a claim about the
        // *consequence*, which two protein sequences cannot settle, so decline
        // and leave the input as authored rather than guess.
        if start_pos.number == 1 && ref_aas[0] != seq.0[0] {
            return None;
        }

        // Maximal runs of residue-wise disagreement. Two or more runs means an
        // unchanged residue separates them, which is the separation the clause
        // is about.
        let mut runs: Vec<(usize, usize)> = Vec::new();
        let mut offset = 0usize;
        while offset < span {
            if ref_aas[offset] == seq.0[offset] {
                offset += 1;
                continue;
            }
            let lo = offset;
            while offset < span && ref_aas[offset] != seq.0[offset] {
                offset += 1;
            }
            runs.push((lo, offset - 1));
        }
        if runs.len() < 2 {
            return None;
        }

        // One member per run: a single-residue run is a substitution
        // (`general.md:56` ranks substitution above delins), a longer one the
        // smaller delins it really is. Members are always certain — a predicted
        // `p.(…)` input carries its marker onto the whole allele, which
        // `wrap_allele_if_split` does.
        let members = runs
            .into_iter()
            .map(|(lo, hi)| {
                let member_start = ProtPos::new(ref_aas[lo], start_pos.number + lo as u64);
                let member_end = ProtPos::new(ref_aas[hi], start_pos.number + hi as u64);
                let edit = if lo == hi {
                    ProteinEdit::Substitution {
                        reference: ref_aas[lo],
                        alternative: seq.0[lo],
                    }
                } else {
                    ProteinEdit::Delins {
                        sequence: AminoAcidSeq::new(seq.0[lo..=hi].to_vec()),
                    }
                };
                ProteinVariant {
                    accession: variant.accession.clone(),
                    gene_symbol: variant.gene_symbol.clone(),
                    loc_edit: LocEdit::new(Interval::new(member_start, member_end), edit),
                }
            })
            .collect();
        Some(members)
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
    /// trait's `get_protein_length` API.
    ///
    /// The trait's default implementation derives the length by probing
    /// `get_protein_sequence(accession, 0, n)` and binary-searching for
    /// the largest accepted `n` (length-equals-largest-accepted-`n`
    /// semantics), preserving the historical behavior. Providers that
    /// store length metadata (e.g. `MockProvider`) override it to return
    /// the length directly without cloning the sequence. Either way, an
    /// empty protein or an accession the provider cannot resolve yields a
    /// length of `0`.
    ///
    /// Returns `None` for empty proteins or accessions the provider
    /// cannot resolve.
    fn discover_protein_length(&self, accession: &str) -> Option<u64> {
        // A length of `0` — an empty protein or an accession the provider
        // cannot resolve — maps to `None`, matching the prior behavior. The
        // `Err(_)` arm is defensive: the trait contract maps unresolvable
        // accessions to `Ok(0)`, but a provider may still surface an error.
        match self.provider.get_protein_length(accession) {
            Ok(0) | Err(_) => None,
            Ok(len) => Some(len),
        }
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
        manufactured: &mut Vec<ManufacturedJunctionExit>,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        use crate::hgvs::interval::RnaInterval;
        use crate::hgvs::variant::{LocEdit, RnaVariant};

        // #972 Task 5: mirrors the `normalize_cds` gate. On a coding
        // transcript, `r.` numbering is CDS-relative — identical to `c.`
        // (see `project::rna` module doc) — so it inherits the same
        // "no confirmed ATG" problem when the transcript's 5' CDS is
        // annotated incomplete (`cds_start_NF`). A non-coding transcript
        // never sets `cds_start_incomplete`, so this check is a no-op for
        // the `n.`-equivalent `r.` case and only fires on a genuinely
        // CDS-relative `r.` input. Decline to re-number — pass through
        // verbatim in every mode — before any of the u/T rewrite, con→delins
        // rewrite, or intronic dispatch below.
        if let Ok(transcript) = self
            .provider
            .get_transcript(&variant.accession.transcript_accession())
        {
            if transcript.cds_start_incomplete {
                let mut warnings = Vec::new();
                if self.config.should_warn_incomplete_cds_start()
                    || self.config.should_reject_incomplete_cds_start()
                {
                    warnings.push(NormalizationWarning::IncompleteCdsStartReference {
                        accession: variant.accession.transcript_accession(),
                        coordinate_system: "r".to_string(),
                    });
                }
                return Ok((HV::Rna(variant.clone()), warnings));
            }
        }

        // Can't normalize variants with unknown edits or positions
        let edit = match variant.loc_edit.edit.inner() {
            Some(e) => e,
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };

        // #736: re-express any RNA base `u` in the edit as `T` before
        // normalization. The transcript reference is stored as DNA, and the
        // parser keeps `u` as a distinct `Base::U` (separate from `Base::T`), so
        // a literal `u` in an insertion / delins would never match the reference
        // and the edit would silently fail to canonicalize or 3'-shift. Re-run
        // on the DNA-rewritten variant; `RnaVariant`'s Display renders `T`→`u`
        // again on output, so this is display-neutral.
        if let Some(dna_edit) = rna_uracil_to_thymine(edit) {
            let new_variant = RnaVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    variant.loc_edit.location.clone(),
                    variant.loc_edit.edit.map_ref(|_| dna_edit.clone()),
                ),
            };
            return self.normalize_rna(&new_variant, manufactured);
        }

        // SVD-WG009: rewrite `con` to `delins`. Re-run on the rewritten
        // variant so downstream passes (the #333 `ins[...]` expansion below,
        // 3' shift, ins→dup, canonical split) all see the delins form.
        if let Some(new_edit) = canonicalize_conversion_to_delins(edit) {
            let new_variant = RnaVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    variant.loc_edit.location.clone(),
                    variant.loc_edit.edit.map_ref(|_| new_edit.clone()),
                ),
            };
            return self.normalize_rna(&new_variant, manufactured);
        }

        // Issue #333 / #1183: expand bracketed / reference-range `ins[...]`
        // payloads, mirroring the genome/cds/tx/mt axes. Runs here — after the
        // u/T and con→delins rewrites, before the `needs_normalization` and
        // intronic dispatches below — so the payload is already a literal by the
        // time either the intronic (tx-delegating) or exonic path sees it, and so
        // an out-of-scope payload is *refused* on `r.` exactly as it is on `c.`
        // instead of passing through unvalidated.
        let rna_accession = variant.accession.transcript_accession();
        if let Some((new_variant, warning)) =
            self.try_expand_rna_ins(variant, edit, &rna_accession)?
        {
            let (result, mut warnings) = self.normalize_rna(&new_variant, manufactured)?;
            warnings.insert(0, warning);
            return Ok((result, warnings));
        }

        // Only normalize indels
        if !needs_normalization(edit) {
            return Ok((HV::Rna(variant.clone()), vec![]));
        }

        // #1052: an uncertain/predicted-wrapped substitution must stay a
        // silent pass-through — see `is_uncertain_real_substitution`'s doc
        // comment. Must run before the intronic dispatch below too.
        if is_uncertain_real_substitution(edit, &variant.loc_edit.edit) {
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

        // #704: route intronic / boundary-spanning `r.` through the same
        // genomic-space machinery the `n.` path uses. On a coding transcript
        // `r.` shares `c.`/`n.` numbering for the exonic anchor and the intron
        // offset is identical across all three axes, so convert the `r.`
        // endpoints to `TxPos`, delegate to the tx intronic dispatch (the
        // `same_intron_tx` split → `normalize_intronic_tx`, else
        // `normalize_boundary_spanning_tx`), then convert the result back to
        // `r.`. Pre-#704 this errored for any intronic `r.`.
        if start_pos.is_intronic() || end_pos.is_intronic() {
            // #1052: real substitutions are validated on the plain exonic path
            // only. Intronic-sub validation would require genomic projection
            // (an explicit non-goal); routing a sub through the tx intronic
            // dispatch this block delegates to would emit
            // `ReducedCapabilityNoGenome` (no genomic data) or a hard
            // `ConversionError` (unmappable intron) on a variant that previously
            // passed through unchanged. Preserve that silent pass-through.
            if is_real_substitution(edit) {
                return Ok((HV::Rna(variant.clone()), vec![]));
            }
            use crate::hgvs::variant::{LocEdit, TxVariant};
            // Accession-aware lookup so an NG/NC-parented input resolves the
            // build-correct chromosome; fall back to the plain transcript. A
            // missing transcript preserves the historical "can't normalize
            // intronic r." signal.
            let plain_accession = variant.accession.transcript_accession();
            let transcript = self
                .provider
                .get_transcript_for_accession(&variant.accession)
                .or_else(|_| self.provider.get_transcript(&plain_accession))
                .map_err(|_| FerroError::IntronicVariant {
                    variant: format!("{}", variant),
                    detail: None,
                })?;
            let cds_info = transcript.cds_start.zip(transcript.cds_end);

            let (tx_start, tx_end) = match (
                self.rna_pos_to_txpos(start_pos, cds_info),
                self.rna_pos_to_txpos(end_pos, cds_info),
            ) {
                (Some(s), Some(e)) => (s, e),
                _ => {
                    return Err(FerroError::IntronicVariant {
                        variant: format!("{}", variant),
                        detail: None,
                    })
                }
            };

            let tx_variant = TxVariant {
                accession: variant.accession.clone(),
                gene_symbol: variant.gene_symbol.clone(),
                loc_edit: LocEdit::with_uncertainty(
                    Interval::new(tx_start, tx_end),
                    variant.loc_edit.edit.with_same_certainty(edit.clone()),
                ),
            };

            // Mirror of the `normalize_tx` intronic dispatch (both intronic in the
            // same intron → single-intron shuffle; one endpoint exonic or a
            // multi-intron span → genomic boundary-spanning).
            let (tx_result, tx_warnings) = if tx_start.is_intronic()
                && tx_end.is_intronic()
                && self.same_intron_tx(&transcript, &tx_start, &tx_end)
            {
                self.normalize_intronic_tx(&tx_variant, &transcript, &tx_start, &tx_end, edit)?
            } else {
                self.normalize_boundary_spanning_tx(
                    &tx_variant,
                    &transcript,
                    &tx_start,
                    &tx_end,
                    edit,
                )?
            };

            // The tx engines always return a bare `HV::Tx`; convert it back to
            // `r.` and finish through the same canonical-split tail as the exonic
            // path (T/U-equivalent rev-comp scan).
            let HV::Tx(tv) = tx_result else {
                return Err(FerroError::IntronicVariant {
                    variant: format!("{}", variant),
                    detail: None,
                });
            };
            let rna_variant = self.txvariant_to_rnavariant(&tv, cds_info)?;
            let uncertain = rna_variant.loc_edit.edit.is_uncertain();
            let (split, mut split_warnings) = self.apply_canonical_split(HV::Rna(rna_variant));
            let mut warnings = tx_warnings;
            warnings.append(&mut split_warnings);
            return Ok((wrap_allele_if_split(split, uncertain), warnings));
        }

        // Try to get transcript (RNA uses the same accession as mRNA transcripts)
        let accession = variant.accession.transcript_accession();
        let transcript = match self.provider.get_transcript(&accession) {
            Ok(t) => t,
            Err(_) => return Ok((HV::Rna(variant.clone()), vec![])),
        };

        // Convert RNA positions to transcript-1 positions.
        //
        // Axis convention (#469, HGVS `background/numbering.md` L58/L61): on a
        // **coding** transcript, `r.` shares `c.`'s **CDS-relative** numbering
        // — `r.123` relates to `c.123` — so every region (`r.-N` 5'UTR, `r.N`
        // CDS, `r.*N` 3'UTR) translates through `cds_start`/`cds_end` via
        // `rna_to_tx_pos`, exactly like `c.`. `r.10` against `cds_start = 100`
        // maps to tx 109 (== `c.10`), not tx 10. This supersedes the
        // transcript-1-relative pin PR #304 added when closing #291 (made on
        // internal-consistency grounds without checking the spec).
        //
        // When the transcript carries no CDS (coordinate-only / mock providers
        // that omit cds_start/end), positive non-UTR bases fall back to
        // transcript-1 indices and UTR positions can't be resolved.
        let cds_info = transcript.cds_start.zip(transcript.cds_end);
        let map_in = |pos: &crate::hgvs::location::RnaPos| -> Option<u64> {
            match cds_info {
                Some((cds_start, cds_end)) => {
                    self.rna_to_tx_pos(pos, cds_start, Some(cds_end)).ok()
                }
                None => {
                    if pos.utr3 || pos.base < 1 {
                        None
                    } else {
                        Some(pos.base as u64)
                    }
                }
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

        // Coordinate-only transcripts fall back to the canonicalize-only path.
        let seq = match transcript.sequence.as_deref() {
            Some(s) => s.as_bytes(),
            None => return Ok((HV::Rna(variant.clone()), vec![])),
        };

        // #1192: the codon-frame gate applies on `r.` exactly as it does on
        // `c.`. `RNA/repeated.md` L24-27 states the restriction for an RNA
        // reference in the same terms `DNA/repeated.md` uses:
        //   > using a coding RNA reference sequence, a repeated sequence
        //   > variant description can be used only for repeat units with a
        //   > length which is a multiple of 3 [...] This restriction only
        //   > applies to the coding sequence, which does not include the UTR
        //   > sequence.
        // and gives `NM_024312.4:r.2686a[10]` as an explicit counter-example.
        // (The previous comment here asserted that `r.` is not an accepted
        // reference type for repeats, which that page contradicts.)
        //
        // The condition mirrors `normalize_cds` exactly — deliberately, since
        // `r.` on a coding transcript is CDS-relative, i.e. the same axis as
        // `c.` (#469). `cds_info` is `None` for a non-coding transcript, which
        // has no reading frame to preserve, and a footprint touching either UTR
        // falls outside `cds_start..=cds_end` and is exempt by the carve-out.
        // Built through the SAME constructor `normalize_cds` uses, which is the
        // point: the two axes cannot drift apart again the way they had before
        // #1192. `cds_info` is `None` for a non-coding transcript, which has no
        // reading frame to preserve, so the gate is `NotApplicable` there; a
        // footprint touching either UTR falls outside `cds_start..=cds_end` and is
        // exempt by the carve-out.
        //
        // Passing the real bounds (rather than dropping them) is what lets the
        // shuffled-tract re-check run on `r.` too. Before #1206 that was a
        // separate `Option` argument, and passing `None` for it made the verdict
        // unconditionally "exempt" and silently re-opened the divergence #1192
        // closed; carrying the bounds inside the gate makes that unexpressible.
        let gate = CodonGate::for_input_span(cds_info, tx_start, tx_end);
        // `mut` bindings are #1202's (its transcript-bound clamp rewrites these
        // in place).
        let (mut new_tx_start, mut new_tx_end, mut new_edit, mut warnings) =
            self.normalize_na_edit(seq, edit, tx_start, tx_end, &boundaries, gate)?;

        // See the identical guard + comment in `normalize_cds` (#1052 follow-up):
        // substitutions are validated-then-returned-unchanged by
        // `normalize_na_edit`; nothing downstream applies to a never-shuffled
        // edit kind, so return the ORIGINAL variant to preserve `Mu::Uncertain`
        // and all other input structure. `ref == alt` still flows through to
        // the identity (`=`) rewrite below.
        if is_real_substitution(edit) {
            return Ok((HV::Rna(variant.clone()), warnings));
        }

        // #704 sub-problem A (mirror of the `normalize_cds`/`normalize_tx`
        // post-check, #670): apply the 3' rule across the exon/intron junction
        // for a purely-exonic `r.` del/dup that comes to rest at an exon's 3'
        // edge. The exon-confined shuffle above only sees spliced bases; when it
        // lands exactly at the exon boundary (`boundaries.right`, the del/dup
        // exon bound), a downstream intron exists, and we have genomic context,
        // re-run the shuffle in genomic space (which spans the junction) and
        // adopt the result only if it actually crossed into the intron. The
        // original endpoints are purely exonic here, so they map to plain
        // (offset-less) `TxPos`. 3'-only; the hot path is untouched.
        //
        // Ordering constraint (#1202, other half documented at the clamp
        // below): a 3'-saturated insertion arrives here with
        // `new_tx_end == boundaries.right + 1`, so it does NOT trip this
        // trigger — and the clamp below then moves it onto exactly
        // `boundaries.right`. Neither block may be reordered past the other.
        if self.config.shuffle_direction == ShuffleDirection::ThreePrime
            && new_tx_end == boundaries.right
        {
            if !self.provider.has_genomic_data() {
                // #1012: mirror of the `normalize_cds`/`normalize_tx` reduced-
                // capability guard for the `r.` axis. The junction-crossing
                // 3'-shuffle needs a genome; without one the exon-confined result
                // flows through, marked degraded.
                warnings.push(NormalizationWarning::ReducedCapabilityNoGenome {
                    variant: format!("{variant}"),
                    capability: "exon/intron junction 3'-shuffle".to_string(),
                });
            } else {
                use crate::hgvs::variant::{LocEdit, TxVariant};
                let boundary_transcript = self
                    .provider
                    .get_transcript_for_accession(&variant.accession)
                    .unwrap_or_else(|_| transcript.clone());
                if boundary_transcript.chromosome.is_some() {
                    let bs_start = TxPos::new(tx_start as i64);
                    let bs_end = TxPos::new(tx_end as i64);
                    let bs_variant = TxVariant {
                        accession: variant.accession.clone(),
                        gene_symbol: variant.gene_symbol.clone(),
                        loc_edit: LocEdit::with_uncertainty(
                            Interval::new(bs_start, bs_end),
                            variant.loc_edit.edit.with_same_certainty(edit.clone()),
                        ),
                    };
                    // Engine errors (no following intron — last exon — no genomic
                    // alignment, …) fall through to the exon-confined result, the
                    // safe pre-#704 behavior.
                    if let Ok((HV::Tx(tv), boundary_warnings)) = self
                        .normalize_boundary_spanning_tx(
                            &bs_variant,
                            &boundary_transcript,
                            &bs_start,
                            &bs_end,
                            edit,
                        )
                    {
                        let crossed_into_intron = tv
                            .loc_edit
                            .location
                            .start
                            .inner()
                            .is_some_and(|p| p.is_intronic())
                            || tv
                                .loc_edit
                                .location
                                .end
                                .inner()
                                .is_some_and(|p| p.is_intronic());
                        if crossed_into_intron {
                            if let Ok(rna_variant) = self.txvariant_to_rnavariant(&tv, cds_info) {
                                let uncertain = rna_variant.loc_edit.edit.is_uncertain();
                                let (split, mut split_warnings) =
                                    self.apply_canonical_split(HV::Rna(rna_variant));
                                let mut combined = warnings;
                                combined.extend(boundary_warnings);
                                combined.append(&mut split_warnings);
                                let produced = wrap_allele_if_split(split, uncertain);
                                // #1723 — the third copy of the gate. This records
                                // NOTHING today: the clause's single reading in
                                // this crate does not reach the `r.` axis, so
                                // `bare_transcript_intronic_leaves` finds none
                                // here. Wired anyway so the provenance follows the
                                // scope rather than lagging behind a widening.
                                record_manufactured_junction_exits(
                                    &produced,
                                    boundary_transcript.chromosome.as_deref(),
                                    manufactured,
                                );
                                return Ok((produced, combined));
                            }
                        }
                    }
                }
            }
        }

        // #1207: mirror `normalize_cds`'s CDS-boundary clamps — `cds_start`
        // (#1170) and `cds_end` (#387) — on the `r.` axis.
        //
        // On a coding transcript `r.` IS the `c.` axis (`background/numbering.md`
        // L58/L61, #469), so it needs the same rewrite for the same reason, and
        // #1192 made that non-optional: once the codon-frame gate is active on
        // `r.`, an unclamped saturated insertion is not merely ugly, it is
        // **non-idempotent**. Pass 1 refuses repeat notation because the *input*
        // span is CDS-resident and emits `r.-1_1ins<A'>`; pass 2 sees a span
        // touching the 5'UTR, does not gate, and produces `r.1_4c[6]` instead.
        // The `c.` spelling of the same edit was already stable precisely
        // because this clamp exists there (`c.1delinsCCC`).
        //
        // #1202's `clamp_insertion_at_sequence_bounds` below does NOT cover
        // this: it fires only at transcript `1` / `len`, whereas these saturate
        // at the CDS bounds, where the far side (`r.-1`, `r.*1`) is perfectly
        // representable. Two different boundaries, both needed.
        //
        // Keyed on the POST-shuffle `new_edit` and on a literal `Insertion`, so
        // a spanning *duplication* — the spec-canonical form, deliberately not
        // clamped on `c.` either (#401) — cannot reach this. Same identity as
        // every other caller: moving the delete boundary one base is exact for
        // any `A'`, so no equivalence check is needed.
        if let Some((cds_start, cds_end)) = cds_info {
            // Mirrors `normalize_cds`'s `AxisRegion::Cds` gate: the input's
            // start sits in the CDS proper. Deliberately keyed on the input
            // rather than on the shuffled result, so a UTR-region input keeps
            // the transcript-bound behavior below.
            let input_starts_in_cds = tx_start >= cds_start && tx_start <= cds_end;
            if input_starts_in_cds {
                let rests_at_cds_start = new_tx_end == cds_start
                    && (new_tx_start == cds_start || new_tx_start + 1 == cds_start);
                let rests_at_cds_end =
                    new_tx_start == cds_end && (new_tx_end == cds_end || new_tx_end == cds_end + 1);
                let anchor = if rests_at_cds_start {
                    Some((cds_start, BoundarySide::FivePrime))
                } else if rests_at_cds_end {
                    Some((cds_end, BoundarySide::ThreePrime))
                } else {
                    None
                };
                // A `Delins` **input** whose shared-affix trim pushed the
                // residual past a CDS boundary restores to its own form instead
                // of being clamped — the #383 (5') / #387 (3') carve-out, which
                // `normalize_cds` applies at both boundaries and
                // `issue_387_canon_cds_end_clamp::three_prime_delins_at_cds_end_suppresses_rewrite`
                // pins there (`c.5delinsAC` stays `c.5delinsAC`). Without this
                // arm the clamp above would rewrite the reduced insertion into a
                // boundary `delins`, so `r.` would newly diverge from `c.` for
                // exactly the inputs that carve-out exists for — trading the
                // divergence this PR fixes for a different one.
                let escaped_five_prime = rests_at_cds_start || new_tx_start < cds_start;
                if matches!(edit, NaEdit::Delins { .. }) && (escaped_five_prime || rests_at_cds_end)
                {
                    new_edit = edit.clone();
                    new_tx_start = tx_start;
                    new_tx_end = tx_end;
                } else if let (
                    Some((anchor, side)),
                    NaEdit::Insertion {
                        sequence: InsertedSequence::Literal(a_prime),
                    },
                ) = (anchor, &new_edit)
                {
                    if let Some(delins) = insertion_to_boundary_delins(seq, a_prime, anchor, side) {
                        new_edit = delins;
                        new_tx_start = anchor;
                        new_tx_end = anchor;
                    }
                }
            }
        }

        // #1202: mirror of the `normalize_tx` clamp. Runs in transcript space
        // (before the `r.` mapping below) so it shares the helper verbatim, and
        // must stay AFTER the junction block above for the same
        // `new_tx_end == boundaries.right` reason documented there. Runs after
        // the CDS clamps above, mirroring `normalize_cds`'s ordering; either of
        // those has already rewritten `new_edit` to a `Delins`, so this is inert
        // after them.
        clamp_insertion_at_sequence_bounds(
            seq,
            &mut new_edit,
            &mut new_tx_start,
            &mut new_tx_end,
            SequenceEnds::WHOLE,
        );

        // Convert each normalized tx position back to a CDS-relative `r.`
        // position via `tx_to_rna_pos`, which restores the correct region for
        // every tx index: `r.-N` (5'UTR, `pos < cds_start`), `r.N` (CDS,
        // `pos - cds_start + 1`), `r.*N` (3'UTR, `pos > cds_end`). This keeps
        // the output axis CDS-relative (== `c.`) and symmetric with `map_in`
        // (#469). Without `cds_info` (coordinate-only / mock) we keep the
        // simple transcript-1 mapping.
        use crate::hgvs::location::RnaPos;
        let map_out = |pos: u64| -> Result<RnaPos, FerroError> {
            match cds_info {
                Some((cds_start, cds_end)) => self.tx_to_rna_pos(pos, cds_start, Some(cds_end)),
                None => Ok(RnaPos::new(pos as i64)),
            }
        };
        let new_start = map_out(new_tx_start)?;
        let new_end = map_out(new_tx_end)?;

        let new_variant = RnaVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                RnaInterval::new(new_start, new_end),
                variant.loc_edit.edit.with_same_certainty(new_edit),
            ),
        };

        // Issue #160 + #165 post-canonicalization split (T/U-equivalent
        // comparison for the rev-comp scan and per-position emissions).
        let uncertain = new_variant.loc_edit.edit.is_uncertain();
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Rna(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split, uncertain), warnings))
    }

    /// Convert an `r.` position to a transcript-numbered [`TxPos`], carrying the
    /// intron offset verbatim (#704). On a coding transcript the exonic anchor
    /// maps through CDS-relative numbering (`rna_to_tx_pos`); without a CDS
    /// (coordinate-only / mock providers) positive non-UTR bases fall back to
    /// transcript-1 indices and UTR / non-positive bases are unresolvable
    /// (`None`), mirroring the exonic `map_in` in `normalize_rna`.
    ///
    /// The intron offset is identical across the `c.`/`n.`/`r.` axes — it is
    /// measured from the nearest exon boundary regardless of how the anchor base
    /// is numbered — so it is copied through unchanged.
    fn rna_pos_to_txpos(&self, pos: &RnaPos, cds_info: Option<(u64, u64)>) -> Option<TxPos> {
        let base = match cds_info {
            Some((cds_start, cds_end)) => self.rna_to_tx_pos(pos, cds_start, Some(cds_end)).ok()?,
            None => {
                if pos.utr3 || pos.base < 1 {
                    return None;
                }
                pos.base as u64
            }
        };
        Some(TxPos {
            base: base as i64,
            offset: pos.offset,
            downstream: false,
        })
    }

    /// Convert a transcript-numbered [`TxPos`] back to an `r.` position,
    /// restoring the CDS-relative region and carrying the intron offset (#704).
    /// Inverse of [`Self::rna_pos_to_txpos`].
    fn txpos_to_rnapos(
        &self,
        pos: &TxPos,
        cds_info: Option<(u64, u64)>,
    ) -> Result<RnaPos, FerroError> {
        let mut rna = match cds_info {
            Some((cds_start, cds_end)) => {
                self.tx_to_rna_pos(pos.base as u64, cds_start, Some(cds_end))?
            }
            None => RnaPos::new(pos.base),
        };
        rna.offset = pos.offset;
        Ok(rna)
    }

    /// Re-express a normalized `n.`-axis [`TxVariant`] (the result of the tx
    /// intronic / boundary-spanning engines) as an `r.` [`RnaVariant`] (#704),
    /// mapping both endpoints back to CDS-relative `r.` positions. The U/T
    /// rendering is left to `RnaVariant`'s Display (`to_rna_string`), so the
    /// DNA-base edit produced by the tx engine needs no translation.
    fn txvariant_to_rnavariant(
        &self,
        tv: &crate::hgvs::variant::TxVariant,
        cds_info: Option<(u64, u64)>,
    ) -> Result<crate::hgvs::variant::RnaVariant, FerroError> {
        use crate::hgvs::interval::RnaInterval;
        use crate::hgvs::variant::{LocEdit, RnaVariant};
        let start_tp =
            tv.loc_edit
                .location
                .start
                .inner()
                .ok_or_else(|| FerroError::ConversionError {
                    msg: "boundary-spanning r. result has no start position".to_string(),
                })?;
        let end_tp =
            tv.loc_edit
                .location
                .end
                .inner()
                .ok_or_else(|| FerroError::ConversionError {
                    msg: "boundary-spanning r. result has no end position".to_string(),
                })?;
        let new_start = self.txpos_to_rnapos(start_tp, cds_info)?;
        let new_end = self.txpos_to_rnapos(end_tp, cds_info)?;
        let new_edit =
            tv.loc_edit
                .edit
                .inner()
                .cloned()
                .ok_or_else(|| FerroError::ConversionError {
                    msg: "boundary-spanning r. result has no edit".to_string(),
                })?;
        Ok(RnaVariant {
            accession: tv.accession.clone(),
            gene_symbol: tv.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                RnaInterval::new(new_start, new_end),
                tv.loc_edit.edit.with_same_certainty(new_edit),
            ),
        })
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
    /// `CodonGate::NotApplicable` (mito is genomic-style and not subject to the
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

        // #1052: an uncertain/predicted-wrapped substitution must stay a
        // silent pass-through — see `is_uncertain_real_substitution`'s doc
        // comment.
        if is_uncertain_real_substitution(edit, &variant.loc_edit.edit) {
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
            let uncertain = canonical.loc_edit.edit.is_uncertain();
            let (split, warnings) = self.apply_canonical_split(HV::Mt(canonical));
            (wrap_allele_if_split(split, uncertain), warnings)
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
        // path: identical to genomic, including the contig-length clamp
        // (#1044, mirroring #1042 on `normalize_genome`). Every provider
        // ERRORS on a past-EOF read rather than clamping, so an unclamped
        // `end + window_size` past the contig 3' end made `get_sequence`
        // fail and dropped the variant into `mt_fallback` — skipping 3'
        // shift AND the delins->inv/sub/dup canonicalization. The shared
        // `clamp_fetch_end_to_contig` helper clamps only when the variant
        // fits within the contig (`end <= len`); a span past the end keeps
        // the raw window and errors into the safe fallback. Wraparound
        // (`start > end`) already returned above, so `end` here is the
        // linear span end.
        //
        // #1691: shared with `normalize_genome` for the same reason
        // `clamp_fetch_end_to_contig` is — the window bounds the shuffle as well
        // as the read, so a repeat tract longer than it walks instead of
        // converging, and two copies of the growth loop would drift exactly as
        // the two copies of the clamp did. The mitochondrial genome carries no
        // homopolymer near that long, so this is a structural guarantee rather
        // than a behaviour change for any real `m.` input; it is here so a
        // synthetic contig cannot reach the defect through this door.
        //
        // Mitochondrial reference is plus-strand and not subject to the
        // codon-frame `unit_len % 3 == 0` restriction (the mito genome
        // has no canonical "CDS" exemption boundary in the same sense as
        // nuclear `c.`; the spec's mito chapter doesn't carry the
        // codon-frame clause), so the gate `normalize_in_grown_window` passes
        // (`NotApplicable`) is the one this path wants.
        let attempt = self.normalize_in_grown_window(&accession, start, end, edit, |raw_end| {
            self.clamp_fetch_end_to_contig(&accession, end, raw_end)
        })?;
        let Some(WindowedNormalization {
            window_start,
            ref_seq,
            mut new_rel_start,
            mut new_rel_end,
            mut new_edit,
            mut warnings,
        }) = attempt
        else {
            // No reference data, or a shuffle still running at the growth cap →
            // minimal-notation cleanup.
            return Ok(mt_fallback(variant));
        };

        // #1217: an insertion driven to rest against either mitochondrial
        // terminus has no valid pair of adjacent anchors there, exactly as on
        // the transcript axes (#1202) and the `g.` axis (#1205) — 5' it escapes
        // as the single-position `m.1ins<A'>`, which `DNA/insertion.md:95-101`
        // rejects by name and ferro's own parser refuses; 3' as
        // `m.<len>_<len+1>ins<A'>`, whose second anchor is past the contig end.
        //
        // Circularity does not supply an escape here. #1205 left `m.` out on the
        // grounds that position 1 has a valid 5' neighbour on a circular
        // reference, so the answer ought to be a wraparound description — but
        // #129 established that ferro **rejects** `m.<high>_<low>ins` (the
        // reversed-range exception at `deletion.md:17` / SVD-WG006 covers
        // `del`/`dup`, and `delins` by composition; `insertion.md` is silent, so
        // the general 5'→3' rule applies, as both mutalyzer and biocommons
        // read it). So the wraparound `ins` is not an available answer, and the
        // single-position `delins` — which needs no reversed range — is. #129
        // also disabled 3'-rule shifting across the origin, so resting at the
        // terminus is already the intended behavior; only its rendering was
        // wrong.
        //
        // Same windowed-flushness reasoning as `normalize_genome`: `ref_seq` is
        // a window, so an end is a real contig bound only where the window is
        // flush with it, and with no length available we cannot establish
        // flushness and so do not clamp.
        let contig_len = self.provider.get_sequence_length(&accession).ok();
        clamp_insertion_at_sequence_bounds(
            ref_seq.as_bytes(),
            &mut new_edit,
            &mut new_rel_start,
            &mut new_rel_end,
            SequenceEnds {
                five_prime: window_start == 0,
                three_prime: contig_len
                    .is_some_and(|len| window_start + ref_seq.len() as u64 >= len),
            },
        );

        let new_start = new_rel_start + window_start;
        let new_end = new_rel_end + window_start;

        let new_variant = MtVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                Interval::new(GenomePos::new(new_start), GenomePos::new(new_end)),
                variant.loc_edit.edit.with_same_certainty(new_edit),
            ),
        };

        // Issue #160 inv-split post-pass (mirrors normalize_genome).
        let uncertain = new_variant.loc_edit.edit.is_uncertain();
        let (split, mut split_warnings) = self.apply_canonical_split(HV::Mt(new_variant));
        warnings.append(&mut split_warnings);
        Ok((wrap_allele_if_split(split, uncertain), warnings))
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

    /// Fetch a genomic window addressed by **1-based inclusive** coordinates.
    ///
    /// All genomic coordinates inside the normalizer's intronic / boundary-
    /// spanning paths are 1-based: `Exon::genomic_start/end` are 1-based, and
    /// every coordinate derived from them via `CoordinateMapper` (`cds_to_genomic`,
    /// `cds_to_genomic_with_intron`, `genomic_to_cds_intronic`) is 1-based too.
    /// `ReferenceProvider::get_genomic_sequence`, by contrast, is **0-based
    /// half-open** (`seq[start..end]`). Calling it with a 1-based `start` reads
    /// one base too far in the +forward direction, which silently corrupts the
    /// homopolymer the shuffle inspects at an exon/intron junction (and only
    /// there — non-shuffling variants round-trip correctly because
    /// `new_g = seq_start + rel - 1` cancels the error, which is why this hid
    /// until minus-strand junction cases like `NM_004992.3:c.378-17del`).
    ///
    /// Converting the 1-based inclusive `[start, end]` request to the 0-based
    /// half-open `[start - 1, end)` the provider expects keeps the callers'
    /// `rel = g - seq_start + 1` arithmetic correct: the base at 1-based
    /// `start` lands at index 0 of the returned string.
    fn get_genomic_sequence_1based(
        &self,
        chromosome: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        // Reject an invalid 1-based start. A `start` of 0 has no valid 0-based
        // index (`start - 1` would saturate to 0 and fetch `[0, end)`), and
        // because callers derive the relative variant position as
        // `rel = g - seq_start + 1`, a `seq_start` of 0 would silently shift
        // every relative coordinate by one. Callers clamp their window starts
        // to >= 1; this guards the invariant rather than papering over it.
        if start < 1 {
            return Err(FerroError::ConversionError {
                msg: format!(
                    "invalid 1-based genomic window start {start} on {chromosome} \
                     (must be >= 1)"
                ),
            });
        }
        self.provider
            .get_genomic_sequence(chromosome, start - 1, end)
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

        // #1012: warn-and-degrade when the reference carries no genomic data.
        // Intronic normalization projects the variant into genomic space to
        // shuffle across the intron; without a genome that step cannot run, so
        // return the canonicalized input unchanged and mark the result degraded
        // rather than erroring. (This is an environmental limitation, not a
        // spec violation — the same variant against a genome-backed reference
        // normalizes further. Distinct from the EINTRONIC spec-form rejection
        // in `normalize_core`, which still errors on a bare transcript.)
        if !self.provider.has_genomic_data() {
            return Ok((
                HV::Cds(self.canonicalize_cds_variant(variant)),
                vec![NormalizationWarning::ReducedCapabilityNoGenome {
                    variant: format!("{variant}"),
                    capability: "intronic normalization".to_string(),
                }],
            ));
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

        // Find the enclosing intron boundaries first (genomic coordinates), so
        // the fetched window can be sized to cover them. A variant-only window
        // leaves a large intron's far edge outside the fetched bases, which
        // fails the minus-strand shuffle boundary check (issue #573). The
        // intron lookup needs only the transcript/mapper, not the sequence.
        // `cds_to_tx` is the flat sequence-axis conversion (#1619): the
        // transcript position is `cds_start + N - 1`, and the exon table is
        // consulted only by the genome-frame lookup below.
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

        // Get a window of genomic sequence around the variant, extended to
        // cover the enclosing intron boundaries (capped for huge introns).
        let window = self.config.window_size;
        let (seq_start, seq_end) =
            intronic_window_bounds(g_start, g_end, intron_g_start, intron_g_end, window);

        // Fetch genomic sequence (1-based window; see `get_genomic_sequence_1based`)
        let genomic_seq = self.get_genomic_sequence_1based(chromosome, seq_start, seq_end)?;

        // Calculate the variant position relative to the fetched sequence
        let rel_start = (g_start - seq_start) + 1; // 1-based
        let rel_end = (g_end - seq_start) + 1;

        // Calculate relative intron boundaries (now within the fetched window).
        // Guard against an inverted intron span before any subtraction: the intron
        // builder only sets genomic coords when `genomic_end >= genomic_start`, so
        // a valid intron is never inverted, but defend against malformed boundaries
        // (e.g. injected via the public `Intron` fields) that would otherwise make
        // `intron_g_end - seq_start` underflow (debug panic / release wraparound).
        if intron_g_end < intron_g_start {
            return Err(FerroError::ConversionError {
                msg: format!("inverted intron span ({intron_g_start}..{intron_g_end})"),
            });
        }
        // When the huge-intron cap fires, `intronic_window_bounds` falls back to
        // the variant-only window, which may not include one or both intron edges.
        // Silently using `saturating_sub` in that case clamps the boundary to 1,
        // which passes the `flip_intronic_for_strand` in-window guard on minus-strand
        // transcripts and produces incorrect normalization. Return an error instead.
        if intron_g_start < seq_start || intron_g_end.saturating_add(1) > seq_end {
            return Err(FerroError::ConversionError {
                msg: format!(
                    "intron span ({intron_g_start}..{intron_g_end}) exceeds \
                     intronic shuffle window ({seq_start}..{seq_end})"
                ),
            });
        }
        let intron_rel_start = intron_g_start - seq_start + 1;
        let intron_rel_end = intron_g_end - seq_start + 1;
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
        // Pass `CodonGate::NotApplicable` so an intronic homopolymer dup/del
        // can emit `[N±k]` repeat notation instead of falling back to
        // the gated `ins<literal>` / plain `del` forms.
        let seq_bytes = work_seq.as_bytes();
        let (work_new_rel_start, work_new_rel_end, new_edit, warnings) = self.normalize_na_edit(
            seq_bytes,
            edit,
            work_rel_start,
            work_rel_end,
            &work_boundaries,
            CodonGate::NotApplicable,
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
            loc_edit: LocEdit::with_uncertainty(
                Interval::new(new_start, new_end),
                variant.loc_edit.edit.with_same_certainty(new_edit),
            ),
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

        // #1012: warn-and-degrade when the reference carries no genomic data —
        // the `n.` mirror of the `normalize_intronic_cds` guard. Intronic
        // normalization needs the genome to shuffle in genomic space; without
        // one, return the canonicalized input unchanged, marked degraded.
        if !self.provider.has_genomic_data() {
            return Ok((
                HV::Tx(self.canonicalize_tx_variant(variant)),
                vec![NormalizationWarning::ReducedCapabilityNoGenome {
                    variant: format!("{variant}"),
                    capability: "intronic normalization".to_string(),
                }],
            ));
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

        // Find the enclosing intron boundaries first so the fetched window can
        // be sized to cover them (issue #573); see `normalize_intronic_cds`.
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

        // Get a window of genomic sequence around the variant, extended to
        // cover the enclosing intron boundaries (capped for huge introns).
        let window = self.config.window_size;
        let (seq_start, seq_end) =
            intronic_window_bounds(g_start, g_end, intron_g_start, intron_g_end, window);

        // Fetch genomic sequence (1-based window; see `get_genomic_sequence_1based`)
        let genomic_seq = self.get_genomic_sequence_1based(chromosome, seq_start, seq_end)?;

        // Calculate the variant position relative to the fetched sequence
        let rel_start = (g_start - seq_start) + 1; // 1-based
        let rel_end = (g_end - seq_start) + 1;

        // Calculate relative intron boundaries (now within the fetched window).
        // See `normalize_intronic_cds`: same inverted-span guard applies here,
        // defending against an `intron_g_end - seq_start` underflow.
        if intron_g_end < intron_g_start {
            return Err(FerroError::ConversionError {
                msg: format!("inverted intron span ({intron_g_start}..{intron_g_end})"),
            });
        }
        // See `normalize_intronic_cds`: same cap-fallback guard applies here.
        if intron_g_start < seq_start || intron_g_end.saturating_add(1) > seq_end {
            return Err(FerroError::ConversionError {
                msg: format!(
                    "intron span ({intron_g_start}..{intron_g_end}) exceeds \
                     intronic shuffle window ({seq_start}..{seq_end})"
                ),
            });
        }
        let intron_rel_start = intron_g_start - seq_start + 1;
        let intron_rel_end = intron_g_end - seq_start + 1;
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
            CodonGate::NotApplicable,
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
            loc_edit: LocEdit::with_uncertainty(
                Interval::new(new_start, new_end),
                variant.loc_edit.edit.with_same_certainty(new_edit),
            ),
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

        // #1012: warn-and-degrade when the reference carries no genomic data.
        // Boundary-spanning normalization shuffles the exon∪intron extent in
        // genomic space; without a genome that cannot run, so return the
        // canonicalized input unchanged and mark the result degraded rather
        // than erroring with `ExonIntronBoundary`.
        if !self.provider.has_genomic_data() {
            return Ok((
                HV::Cds(self.canonicalize_cds_variant(variant)),
                vec![NormalizationWarning::ReducedCapabilityNoGenome {
                    variant: format!("{variant}"),
                    capability: "exon/intron boundary normalization".to_string(),
                }],
            ));
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

        // The shuffle boundary is the genomic union of the relevant exon and
        // intron; size the fetch window to cover it (not just the variant ±
        // window), so a long adjacent intron does not push the boundary past
        // the fetched sequence. Mirrors the intronic path's `intronic_window_bounds`.
        let (bound_g_start, bound_g_end) =
            self.get_boundary_spanning_genomic_extent(transcript, &mapper, start_pos, end_pos)?;

        // Fetch genomic sequence (1-based window; see `get_genomic_sequence_1based`).
        let window = self.config.window_size;
        let (seq_start, seq_end) =
            intronic_window_bounds(g_start, g_end, bound_g_start, bound_g_end, window);
        let genomic_seq = self.get_genomic_sequence_1based(chromosome, seq_start, seq_end)?;

        // Calculate relative positions (1-based)
        let rel_start = (g_start - seq_start) + 1;
        let rel_end = (g_end - seq_start) + 1;

        // Boundaries within the fetched window. If the huge-intron cap in
        // `intronic_window_bounds` fired, the boundary may lie outside the
        // window; clamp into range so the shuffle simply can't reach that far
        // (the homopolymer it follows ends well before any realistic cap).
        let seq_len = genomic_seq.len() as u64;
        let boundaries = Boundaries::new(
            bound_g_start.saturating_sub(seq_start) + 1,
            (bound_g_end.saturating_sub(seq_start) + 1).min(seq_len),
        );

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
        // `CodonGate::NotApplicable` to match the intronic exemption. (A
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
            CodonGate::NotApplicable,
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
            loc_edit: LocEdit::with_uncertainty(
                Interval::new(new_start, new_end),
                variant.loc_edit.edit.with_same_certainty(new_edit),
            ),
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

    /// True iff two intronic CDS positions lie in the **same** intron. Used by
    /// the dispatch to tell a normal intronic variant (one intron → the
    /// dedicated intronic shuffle) from a deletion spanning multiple introns
    /// (→ the genomic boundary path, #670). If the intron can't be resolved for
    /// either endpoint, defaults to `true` so the conservative single-intron
    /// path is kept (no behavior change for those).
    fn same_intron(
        &self,
        transcript: &crate::reference::transcript::Transcript,
        a: &CdsPos,
        b: &CdsPos,
    ) -> bool {
        use crate::convert::CoordinateMapper;
        let mapper = CoordinateMapper::new(transcript);
        let intron_number = |p: &CdsPos| -> Option<u32> {
            let tx = mapper.cds_to_tx(p).ok()?;
            let base = u64::try_from(tx.base).ok()?;
            transcript
                .find_intron_at_tx_boundary(base, p.offset.unwrap_or(0))
                .map(|i| i.number)
        };
        match (intron_number(a), intron_number(b)) {
            (Some(x), Some(y)) => x == y,
            _ => true,
        }
    }

    /// Genomic extent within which a boundary-spanning variant may shift: the
    /// union (in genomic coordinates) of the relevant exon and intron.
    ///
    /// Three cases:
    /// - **One endpoint intronic** (classic boundary-spanning): union the exon
    ///   of the exonic endpoint with the intron of the intronic endpoint.
    /// - **Both endpoints exonic** (#670): a purely-exonic indel resting at an
    ///   exon's 3' edge that must be free to shift across the exon/intron
    ///   junction under the 3' rule. Union the exon containing the 3' endpoint
    ///   with the intron *immediately following* it (in transcript order).
    ///   Capping at that single intron — never the next exon — preserves the
    ///   exon/exon suppression rule (numbering.md NOTE; deletion.md).
    /// - **Both endpoints intronic in different introns** (#670): a deletion
    ///   spanning intron→exon→intron. Union the two endpoints' introns; the
    ///   exon(s) between them fall inside the min/max span.
    ///
    /// Returns `(combined_g_start, combined_g_end)` (inclusive, 1-based). The
    /// caller sizes the fetch window to cover this so the shuffle boundary is
    /// always in-window (a long adjacent intron must not push the boundary past
    /// the fetched sequence — issue #573 / bug surfaced by #670).
    fn get_boundary_spanning_genomic_extent(
        &self,
        transcript: &crate::reference::transcript::Transcript,
        mapper: &crate::convert::CoordinateMapper,
        start_pos: &CdsPos,
        end_pos: &CdsPos,
    ) -> Result<(u64, u64), FerroError> {
        let exon_at_pos =
            |pos: &CdsPos| -> Result<&crate::reference::transcript::Exon, FerroError> {
                let tx_pos = mapper.cds_to_tx(pos)?;
                let tx_pos_base =
                    u64::try_from(tx_pos.base).map_err(|_| FerroError::ConversionError {
                        msg: format!(
                            "Negative transcript position {} during boundary normalization",
                            tx_pos.base
                        ),
                    })?;
                transcript
                    .exon_at(tx_pos_base)
                    .ok_or_else(|| FerroError::ConversionError {
                        msg: "Could not find exon for boundary normalization".to_string(),
                    })
            };

        // #670: both endpoints intronic in DIFFERENT introns (a deletion
        // spanning intron→exon→intron). Union the two introns' genomic coords;
        // the exon(s) between them are covered by the min/max span. Bounding at
        // the two outer intron edges keeps the shuffle from sliding out of the
        // spanned region into a flanking exon.
        if start_pos.is_intronic() && end_pos.is_intronic() {
            let intron_genomic = |pos: &CdsPos| -> Result<(u64, u64), FerroError> {
                let tx = mapper.cds_to_tx(pos)?;
                let base = u64::try_from(tx.base).map_err(|_| FerroError::ConversionError {
                    msg: format!(
                        "Negative transcript position {} during boundary normalization",
                        tx.base
                    ),
                })?;
                let intron = transcript
                    .find_intron_at_tx_boundary(base, pos.offset.unwrap_or(0))
                    .ok_or_else(|| FerroError::ConversionError {
                        msg: "Could not find intron for boundary normalization".to_string(),
                    })?;
                match (intron.genomic_start, intron.genomic_end) {
                    (Some(s), Some(e)) => Ok((s, e)),
                    _ => Err(FerroError::ConversionError {
                        msg: "Intron has no genomic coordinates".to_string(),
                    }),
                }
            };
            let (a_start, a_end) = intron_genomic(start_pos)?;
            let (b_start, b_end) = intron_genomic(end_pos)?;
            return Ok((a_start.min(b_start), a_end.max(b_end)));
        }

        let (exon, intron) = if start_pos.is_intronic() || end_pos.is_intronic() {
            // Classic boundary-spanning: one endpoint exonic, one intronic.
            let (exonic_pos, intronic_pos) = if start_pos.is_intronic() {
                (end_pos, start_pos)
            } else {
                (start_pos, end_pos)
            };
            let exon = exon_at_pos(exonic_pos)?;
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
            (exon, intron)
        } else {
            // #670: purely-exonic edge variant. Union the exon containing the 3'
            // endpoint with the intron immediately downstream of it (its 5' tx
            // boundary is the exon's 3' end, `c.{exon.end}+1`).
            let exon = exon_at_pos(end_pos)?;
            let intron = transcript
                .find_intron_at_tx_boundary(exon.end, 1)
                .ok_or_else(|| FerroError::ConversionError {
                    msg: "Could not find following intron for boundary normalization".to_string(),
                })?;
            (exon, intron)
        };

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

        Ok((
            exon_g_start.min(intron_g_start),
            exon_g_end.max(intron_g_end),
        ))
    }

    // ------------------------------------------------------------------
    // `n.` (transcript-coordinate) mirrors of the `c.` exon/intron 3'-rule
    // machinery (#704). These mirror `cds_pos_to_genomic`, `same_intron`,
    // `get_boundary_spanning_genomic_extent`, and `normalize_boundary_spanning_cds`
    // but key off the transcript position directly — a `TxPos.base` IS the
    // 1-based transcript index, so no `cds_to_tx` step is needed. The `c.`
    // path is untouched. (`n.` has no CDS frame, so reuse-via-conversion isn't
    // possible; the codebase already keeps `normalize_intronic_cds`/`_tx` as
    // separate mirrors for the same reason.)
    // ------------------------------------------------------------------

    /// Convert a transcript position (exonic or intronic) to genomic coordinate.
    /// Mirror of [`Self::cds_pos_to_genomic`] for `n.` coordinates.
    fn tx_pos_to_genomic(
        &self,
        mapper: &crate::convert::CoordinateMapper,
        pos: &TxPos,
    ) -> Result<u64, FerroError> {
        if pos.is_intronic() {
            mapper.tx_to_genomic_with_intron(pos)
        } else {
            mapper
                .tx_to_genomic(pos)?
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!("transcript position {} not in exons", pos.base),
                })
        }
    }

    /// True iff two intronic transcript positions lie in the **same** intron.
    /// Mirror of [`Self::same_intron`] for `n.` coordinates.
    fn same_intron_tx(
        &self,
        transcript: &crate::reference::transcript::Transcript,
        a: &TxPos,
        b: &TxPos,
    ) -> bool {
        let intron_number = |p: &TxPos| -> Option<u32> {
            let base = u64::try_from(p.base).ok()?;
            transcript
                .find_intron_at_tx_boundary(base, p.offset.unwrap_or(0))
                .map(|i| i.number)
        };
        match (intron_number(a), intron_number(b)) {
            (Some(x), Some(y)) => x == y,
            _ => true,
        }
    }

    /// Genomic extent within which a boundary-spanning `n.` variant may shift.
    /// Mirror of [`Self::get_boundary_spanning_genomic_extent`] for `n.` coords.
    fn get_boundary_spanning_genomic_extent_tx(
        &self,
        transcript: &crate::reference::transcript::Transcript,
        start_pos: &TxPos,
        end_pos: &TxPos,
    ) -> Result<(u64, u64), FerroError> {
        let exon_at_pos = |pos: &TxPos| -> Result<&crate::reference::transcript::Exon, FerroError> {
            let base = u64::try_from(pos.base).map_err(|_| FerroError::ConversionError {
                msg: format!(
                    "Negative transcript position {} during boundary normalization",
                    pos.base
                ),
            })?;
            transcript
                .exon_at(base)
                .ok_or_else(|| FerroError::ConversionError {
                    msg: "Could not find exon for boundary normalization".to_string(),
                })
        };
        let intron_at = |pos: &TxPos| -> Result<&crate::reference::transcript::Intron, FerroError> {
            let base = u64::try_from(pos.base).map_err(|_| FerroError::ConversionError {
                msg: format!(
                    "Negative transcript position {} during boundary normalization",
                    pos.base
                ),
            })?;
            transcript
                .find_intron_at_tx_boundary(base, pos.offset.unwrap_or(0))
                .ok_or_else(|| FerroError::ConversionError {
                    msg: "Could not find intron for boundary normalization".to_string(),
                })
        };
        let intron_genomic =
            |intron: &crate::reference::transcript::Intron| -> Result<(u64, u64), FerroError> {
                match (intron.genomic_start, intron.genomic_end) {
                    (Some(s), Some(e)) => Ok((s, e)),
                    _ => Err(FerroError::ConversionError {
                        msg: "Intron has no genomic coordinates".to_string(),
                    }),
                }
            };

        // Both intronic in DIFFERENT introns: union the two introns' genomic span.
        if start_pos.is_intronic() && end_pos.is_intronic() {
            let (a_start, a_end) = intron_genomic(intron_at(start_pos)?)?;
            let (b_start, b_end) = intron_genomic(intron_at(end_pos)?)?;
            return Ok((a_start.min(b_start), a_end.max(b_end)));
        }

        let (exon, intron) = if start_pos.is_intronic() || end_pos.is_intronic() {
            // Classic boundary-spanning: one endpoint exonic, one intronic.
            let (exonic_pos, intronic_pos) = if start_pos.is_intronic() {
                (end_pos, start_pos)
            } else {
                (start_pos, end_pos)
            };
            (exon_at_pos(exonic_pos)?, intron_at(intronic_pos)?)
        } else {
            // Purely-exonic edge variant (#670 sub-problem A): union the exon
            // containing the 3' endpoint with the intron immediately downstream.
            let exon = exon_at_pos(end_pos)?;
            let intron = transcript
                .find_intron_at_tx_boundary(exon.end, 1)
                .ok_or_else(|| FerroError::ConversionError {
                    msg: "Could not find following intron for boundary normalization".to_string(),
                })?;
            (exon, intron)
        };

        let (exon_g_start, exon_g_end): (u64, u64) = match (exon.genomic_start, exon.genomic_end) {
            (Some(s), Some(e)) => (s, e),
            _ => {
                return Err(FerroError::ConversionError {
                    msg: "Exon has no genomic coordinates".to_string(),
                })
            }
        };
        let (intron_g_start, intron_g_end) = intron_genomic(intron)?;

        Ok((
            exon_g_start.min(intron_g_start),
            exon_g_end.max(intron_g_end),
        ))
    }

    /// Normalize a transcript (`n.`) variant that spans an exon-intron boundary
    /// (or multiple introns). Mirror of [`Self::normalize_boundary_spanning_cds`]:
    /// convert to genomic, 3'-shuffle there over the exon∪intron extent, convert
    /// the result back to transcript notation.
    fn normalize_boundary_spanning_tx(
        &self,
        variant: &TxVariant,
        transcript: &crate::reference::transcript::Transcript,
        start_pos: &TxPos,
        end_pos: &TxPos,
        edit: &NaEdit,
    ) -> Result<(HgvsVariant, Vec<NormalizationWarning>), FerroError> {
        use crate::convert::CoordinateMapper;

        // #1012: warn-and-degrade when the reference carries no genomic data —
        // the `n.` mirror of the `normalize_boundary_spanning_cds` guard.
        if !self.provider.has_genomic_data() {
            return Ok((
                HV::Tx(self.canonicalize_tx_variant(variant)),
                vec![NormalizationWarning::ReducedCapabilityNoGenome {
                    variant: format!("{variant}"),
                    capability: "exon/intron boundary normalization".to_string(),
                }],
            ));
        }

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

        let g_start = self.tx_pos_to_genomic(&mapper, start_pos)?;
        let g_end = self.tx_pos_to_genomic(&mapper, end_pos)?;

        // On minus strand, genomic coords may be reversed relative to tx order.
        let swapped = g_start > g_end;
        let (g_start, g_end) = if swapped {
            (g_end, g_start)
        } else {
            (g_start, g_end)
        };

        let (bound_g_start, bound_g_end) =
            self.get_boundary_spanning_genomic_extent_tx(transcript, start_pos, end_pos)?;

        let window = self.config.window_size;
        let (seq_start, seq_end) =
            intronic_window_bounds(g_start, g_end, bound_g_start, bound_g_end, window);
        let genomic_seq = self.get_genomic_sequence_1based(chromosome, seq_start, seq_end)?;

        let rel_start = (g_start - seq_start) + 1;
        let rel_end = (g_end - seq_start) + 1;

        let seq_len = genomic_seq.len() as u64;
        let boundaries = Boundaries::new(
            bound_g_start.saturating_sub(seq_start) + 1,
            (bound_g_end.saturating_sub(seq_start) + 1).min(seq_len),
        );

        // Flip to transcript view before normalization (see the cds mirror).
        let (work_seq, work_rel_start, work_rel_end, work_boundaries) = flip_intronic_for_strand(
            transcript.strand,
            &genomic_seq,
            rel_start,
            rel_end,
            &boundaries,
        )?;

        let seq_bytes = work_seq.as_bytes();
        let (work_new_rel_start, work_new_rel_end, new_edit, warnings) = self.normalize_na_edit(
            seq_bytes,
            edit,
            work_rel_start,
            work_rel_end,
            &work_boundaries,
            CodonGate::NotApplicable,
        )?;

        let (new_rel_start, new_rel_end) = unflip_intronic_positions(
            transcript.strand,
            work_seq.len() as u64,
            work_new_rel_start,
            work_new_rel_end,
        );

        let new_g_start = seq_start + new_rel_start - 1;
        let new_g_end = seq_start + new_rel_end - 1;

        let new_start = mapper.genomic_to_tx_with_intron(new_g_start)?;
        let new_end = mapper.genomic_to_tx_with_intron(new_g_end)?;

        let (new_start, new_end) = if swapped {
            (new_end, new_start)
        } else {
            (new_start, new_end)
        };

        let new_variant = TxVariant {
            accession: variant.accession.clone(),
            gene_symbol: variant.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(
                Interval::new(new_start, new_end),
                variant.loc_edit.edit.with_same_certainty(new_edit),
            ),
        };

        Ok((HV::Tx(new_variant), warnings))
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
        // See [`CodonGate`]. Ask it about the span you are deciding about:
        // `input_span_is_coding()` for the input edit's own span,
        // `span_is_coding(a, b)` for any other (a shuffled tract, say).
        gate: CodonGate,
    ) -> Result<(u64, u64, NaEdit, Vec<NormalizationWarning>), FerroError> {
        let mut warnings = Vec::new();

        // An insertion resting at interbase 0 — "before base 1" — reaches here
        // with `start == 0`, which is not a position any HGVS axis has (#1282).
        //
        // It is produced by the cis-allele derivation, not by the shuffle: an
        // allele whose applied sequence is the reference with bases added at the
        // very front (`g.[3_4insT;1T>A]` on a leading `T` run denotes exactly
        // that) has its changed block trimmed to the 3' side, leaving the
        // payload before the first base. That is a real, describable variant —
        // it simply has no `ins` spelling, because `g.0_1ins` does not exist.
        //
        // The answer is the one every other boundary clamp in this module gives:
        // "insert A' immediately 5' of anchor" IS "delete anchor, insert
        // A' ++ ref[anchor]". So rewrite to `1delins<A' ++ ref[1]>` here, before
        // the position ever reaches `hgvs_pos_to_index`.
        //
        // This has to happen *before* the shuffle rather than in
        // `clamp_insertion_at_sequence_bounds`, which is keyed on the
        // post-shuffle resting shape (`start == end == 1`) and so never sees a
        // member that arrives already at interbase 0. The two are the same
        // repair at different stages; #418 patched a third instance of the class
        // inside the delins canonicalisation recursion, one site at a time.
        //
        // The constructed delins goes back through this function rather than
        // being returned as-is. `insertion_to_boundary_delins` spells
        // `1delins<A' ++ ref[1]>` without asking whether that reduces, and a
        // payload sharing an affix with `ref[1]` makes it trimmable (a payload
        // *equal* to it would make it a `dup`). Every other call site of the
        // helper runs after a completed shuffle walk, whose completion property
        // already rules a reducible result out; this one fires on a cis-derived
        // member before any shuffle, so it does not inherit that protection. The
        // recursion is what makes the result canonical rather than merely valid,
        // and so a fixed point — see
        // `a_reducible_boundary_payload_still_lands_canonical`. It cannot
        // re-enter this branch: it is dispatched at `start == 1`.
        if start == 0 {
            if let NaEdit::Insertion {
                sequence: InsertedSequence::Literal(seq),
            } = edit
            {
                return match insertion_to_boundary_delins(ref_seq, seq, 1, BoundarySide::FivePrime)
                {
                    // Names avoid shadowing `edit`: the sibling arm below
                    // interpolates the *outer* one, and two `edit`s a line apart
                    // meaning different things is a trap.
                    Some(delins) => {
                        let (delins_start, delins_end, canonical, delins_warnings) =
                            self.normalize_na_edit(ref_seq, &delins, 1, 1, boundaries, gate)?;
                        warnings.extend(delins_warnings);
                        Ok((delins_start, delins_end, canonical, warnings))
                    }
                    // The rewrite spells the result as `1delins<payload ++ ref[1]>`,
                    // so it needs the entity's first base. That is unavailable
                    // only when the provider handed back an empty sequence or a
                    // first byte that is not a nucleotide — a reference problem,
                    // not a coordinate one. Say which, rather than reusing the
                    // "no such position" wording below and sending the reader
                    // after the coordinate.
                    None => Err(FerroError::InvalidCoordinates {
                        msg: format!(
                            "cannot rewrite the interbase-0 insertion {edit} as a boundary \
                             delins: the reference sequence has no usable first base"
                        ),
                    }),
                };
            }
            // Any *other* edit at position 0 — anything that is not a literal
            // insertion — is a coordinate that should never have been
            // constructed. Fail with an error rather than letting it underflow
            // into a garbage index in release (#1282).
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "normalization produced position 0 for {edit}, which no HGVS axis has"
                ),
            });
        }

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
            // canonical path).
            //
            // False for `Repeat` / `MultiRepeat` consistency mismatches
            // (issues #214 / #279): the per-unit declaration is part of
            // the user's form and the normalizer passes the description
            // through verbatim. Also false for `Substitution` (#1052):
            // real substitutions now DO reach here (routed by
            // `needs_normalization` for validation only, see rules.rs)
            // and are returned unchanged by the `Substitution` arm below
            // — the canonical Display keeps the stated ref base, so
            // reporting `corrected: true` would be dishonest.
            //
            // `delins` with a stated `deleted` sequence is validated by
            // `validate_reference`'s `NaEdit::Delins` arm (#486); the
            // canonical form drops `deleted` / `deleted_length` (see the
            // Delins arm in `canonicalize_edit`), so `corrected` is `true`.
            let corrected = !matches!(
                edit,
                NaEdit::Repeat { .. } | NaEdit::MultiRepeat { .. } | NaEdit::Substitution { .. }
            );
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

        // A bracketed inserted repeat with an exact count — either a multi-base
        // `ins<seq>[N]` (`SequenceRepeat`) or a single-base `ins<base>[N]`
        // (`Repeat`) — is syntactic shorthand the plain-insertion path never
        // canonicalizes: both `bases()` to `None`, so the alt-extraction below
        // returns the edit verbatim. Rewrite the degenerate counts, matching
        // mutalyzer:
        //   [0] — zero copies inserted is a no-op → whole-entity identity (`g.=`).
        //   [N>=1] — N inserted copies of the unit are equivalent to the plain
        //         `ins<unit repeated N times>` literal, which the insertion
        //         path below canonicalizes: 1 copy of the adjacent reference
        //         becomes `dup` (`insertion_to_duplication`, `duplication.md`
        //         L17), and >= 2 copies become `[N]` repeat notation merged
        //         with any abutting reference tract (`insertion_to_repeat`,
        //         `repeated.md`; e.g. `g.4_5insGT[5]` against an adjacent
        //         `GT` tract → `g.1_4GT[7]`). Rewrite to the literal and recurse.
        // The rewrite is purely syntactic on the stated sequence, so it runs
        // before shuffling; `Literal::bases()` is `Some`, so the recursion
        // reaches the alt-extraction below without re-entering this
        // bracketed-count intercept. See #920.
        if let NaEdit::Insertion { sequence } = edit {
            let degenerate = match sequence {
                InsertedSequence::SequenceRepeat {
                    sequence,
                    count: RepeatCount::Exact(n),
                } => Some((*n, sequence.clone())),
                InsertedSequence::Repeat {
                    base,
                    count: RepeatCount::Exact(n),
                } => Some((*n, Sequence::new(vec![*base]))),
                _ => None,
            };
            if let Some((count, unit)) = degenerate {
                if count == 0 {
                    return Ok((start, end, NaEdit::whole_entity_identity(), warnings));
                }
                // Materializing the expansion bounds memory: a pathologically
                // large count (`insGT[1e12]`, or a `u64` that overflows the
                // capacity product) is left to fall through to the verbatim
                // pass-through below — the pre-#920 behavior — rather than
                // allocating gigabytes or panicking. No real or conformance
                // repeat approaches this cap.
                const MAX_INS_REPEAT_EXPANSION: usize = 1 << 20; // 1 MiB of bases
                if let Some(expanded_len) = unit.len().checked_mul(count as usize) {
                    if expanded_len <= MAX_INS_REPEAT_EXPANSION {
                        let mut expanded = Vec::with_capacity(expanded_len);
                        for _ in 0..count {
                            expanded.extend_from_slice(unit.bases());
                        }
                        let literal = NaEdit::Insertion {
                            sequence: InsertedSequence::Literal(Sequence::new(expanded)),
                        };
                        let (new_start, new_end, new_edit, mut child_warnings) = self
                            .normalize_na_edit(ref_seq, &literal, start, end, boundaries, gate)?;
                        warnings.append(&mut child_warnings);
                        return Ok((new_start, new_end, new_edit, warnings));
                    }
                }
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
                    return self.normalize_na_edit(ref_seq, &del, start, end, boundaries, gate);
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
                            //
                            // Recurse into `normalize_na_edit` with the reduced
                            // deletion so the full del pipeline (3'-shift
                            // `shuffle()` + repeat contraction) runs — exactly as
                            // the `Insertion` arm below already does. Returning
                            // the trimmed deletion directly skipped the shift,
                            // making the output non-idempotent and non-confluent
                            // with the same edit written as a plain `del`:
                            // `g.257_259delinsT` emitted the unshifted
                            // `g.258_259del`, while a direct `g.258_259del`
                            // 3'-shifts (and contracts) to `g.258_262A[3]`, so
                            // `norm(norm(x)) != norm(x)`. Issue #1157.
                            let new_edit = NaEdit::Deletion {
                                sequence: None,
                                length: None,
                            };
                            let (new_start, new_end, new_edit, mut child_warnings) = self
                                .normalize_na_edit(
                                    ref_seq,
                                    &new_edit,
                                    index_to_hgvs_pos(s0),
                                    e0 as u64,
                                    boundaries,
                                    gate,
                                )?;
                            let mut merged = warnings.clone();
                            merged.append(&mut child_warnings);
                            return Ok((new_start, new_end, new_edit, merged));
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
                                    gate,
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
                            //
                            // Recurse so the reduced duplication takes the same
                            // 3'-shift `shuffle()` a plain `dup` receives, mirror
                            // of the Deletion/Insertion arms. Direct-return here
                            // skipped the shift, so `g.258delinsGG` in a G-tract
                            // emitted the unshifted `g.258dup` while a direct
                            // `g.258dup` shifts to `g.262dup` — non-idempotent and
                            // non-confluent. Issue #1157.
                            let new_edit = NaEdit::Duplication {
                                sequence: None,
                                length: None,
                                uncertain_extent: None,
                            };
                            let (new_start, new_end, new_edit, mut child_warnings) = self
                                .normalize_na_edit(
                                    ref_seq,
                                    &new_edit,
                                    index_to_hgvs_pos(s0),
                                    e0 as u64,
                                    boundaries,
                                    gate,
                                )?;
                            let mut merged = warnings.clone();
                            merged.append(&mut child_warnings);
                            return Ok((new_start, new_end, new_edit, merged));
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
                                    substitution_reference: None,
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
                    // A one-base residue is a substitution to that base's
                    // complement, not an inversion: `inversion.md:5` defines an
                    // inversion as more than one nucleotide and `:16` names the
                    // substitution as its replacement. #1079 enforces the same
                    // rule in the parser, which can only reject `g.234inv`
                    // because deriving the substitution needs the reference
                    // sequence; here we have it. Emitting identity instead
                    // silently discarded the variant (#1249).
                    if new_e == new_s + 1 {
                        if let Some((reference, alternative)) =
                            rules::complementary_substitution(ref_seq[new_s])
                        {
                            let pos = index_to_hgvs_pos(new_s);
                            return Ok((
                                pos,
                                pos,
                                NaEdit::Substitution {
                                    reference,
                                    alternative,
                                },
                                warnings,
                            ));
                        }
                        // Not a base we can complement into a typed
                        // substitution. **Reachable since #1318**, and the
                        // comment here used to say the opposite.
                        //
                        // The old justification was that an unmodelled byte was
                        // returned unchanged by `complement`, so it read as
                        // self-complementary and collapsed to identity before
                        // reaching here. `complement_base` answers `None` for
                        // such a byte instead, so `is_self_complementary` is
                        // now `false` for it and `shorten_inversion` hands back
                        // a one-base residue rather than `None` — measured:
                        // `shorten_inversion(b"X", 0, 1)` was `None`, is now
                        // `Some((0, 1))`.
                        //
                        // Landing here is the right outcome: leave the
                        // inversion as authored rather than claim an identity
                        // for a byte we cannot complement, which is the #1249
                        // class of silent loss. Emitting nothing typed is a
                        // refusal, not a fallback nobody reaches.
                        return Ok((start, end, canonicalize_edit(edit), warnings));
                    }
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
                let end_idx = hgvs_pos_to_index(end); // 0-based inclusive end of input range

                // Normalize the repeat
                match rules::normalize_repeat(
                    ref_seq,
                    pos_idx,
                    end_idx,
                    &repeat_unit,
                    *specified_count,
                    gate.input_span_is_coding(),
                    self.config.shuffle_direction,
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
                            // Re-enter on the derived insertion rather than
                            // returning it raw (#1204), exactly as the
                            // `DupToRepeatResult::GatedInsertion` arm below does.
                            // All three spellings of a gated expansion — `unit[N]`
                            // here, `dup` there, a directly written `ins` — then
                            // reach the SAME insertion canonicalization, so they
                            // agree on the output form. Returning the literal
                            // straight out left this one path unable to promote to
                            // the `dup` the spec prescribes (`repeated.md` L22),
                            // which made `norm(unit[N])` an `ins` while
                            // `norm(norm(unit[N]))` was a `dup` — a
                            // non-idempotency the moment the insertion path
                            // learned to promote. The rule layer also anchors this
                            // insertion at the tract's 3' flank regardless of
                            // direction, so re-entry is what shuffles it to the
                            // direction-appropriate resting place and subjects it
                            // to the same boundary clamps.
                            //
                            // Terminates: the derived edit is an `Insertion`, and
                            // the `Insertion` arm never re-enters this `Repeat`
                            // arm. The `gate` is threaded unchanged — same
                            // transcript, so the same CDS bounds (#1185).
                            let (rec_start, rec_end, rec_edit, mut rec_warnings) = self
                                .normalize_na_edit(
                                    ref_seq, &ins_edit, ins_start, ins_end, boundaries, gate,
                                )?;
                            warnings.append(&mut rec_warnings);
                            return Ok((rec_start, rec_end, rec_edit, warnings));
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
            // MultiRepeat (compound tandem repeat, e.g. `GT[2]GC[2]…`) is
            // deliberately NOT 3'-shifted: unlike a single-unit `Repeat`, a
            // multi-unit tract has no canonical shuffle target — the spec
            // describes it as reported, and the rules-layer shuffler
            // (`normalize_repeat`) only handles one repeat unit + count. It is
            // still routed here (via `needs_normalization` → `true`) for the
            // reference-tract validation (`validate_multirepeat_tract`, run at the
            // top of this function) — that validation, not shuffling, is why the
            // flag stays `true`. Return the (validated) edit unchanged. Giving it
            // an explicit arm — rather than letting it fall to the generic `_`
            // below — makes this pass-through intentional and documented, not a
            // silent contradiction of the `needs_normalization` flag (#953).
            NaEdit::MultiRepeat { .. } => return Ok((start, end, edit.clone(), warnings.clone())),
            // A real substitution reaches here only via `needs_normalization`'s
            // validation routing (#1052). The reference-base check already ran
            // at the top of this function; substitutions are never shuffled, so
            // return the (validated) edit unchanged. Explicit arm — not the
            // generic `_` — to make the pass-through intentional and documented
            // (mirrors the `MultiRepeat` arm, #953).
            NaEdit::Substitution { .. } => return Ok((start, end, edit.clone(), warnings.clone())),
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
                                gate.input_span_is_coding(),
                                // #1210: the verdict above is the INPUT span's;
                                // the bounds let the helper re-ask about the tract
                                // the repeat would actually occupy, which differs
                                // whenever the edit shifts across `cds_end`. This
                                // helper is the one place that still needs the two
                                // separately — see `CodonGate::cds_bounds`.
                                gate.cds_bounds(),
                                self.config.shuffle_direction,
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

                    // No codon-frame gate on the dup paths below (#1204). The
                    // `repeated.md` exception forbids *repeat notation* for a
                    // non-triplet unit in the CDS, and names `dup` as the
                    // replacement for exactly this case — "use
                    // `NM_024312.4:c.2692_2693dup` and **not**
                    // `NM_024312.4:c.2686A[10]`" (DNA L22, RNA L25). The gate is
                    // therefore applied where repeat notation is decided —
                    // `insertion_to_repeat` above and `duplication_to_repeat` in
                    // the Duplication arm, both of which return `None` / reroute
                    // here when it fires — and must NOT be re-applied here, or
                    // the prescribed `dup` is unreachable and the spec's *other*
                    // replacement (a flat `ins`) is emitted for both cases.
                    //
                    // Which of the two the spec prescribes is decided by the
                    // reference, not by the gate, and the three paths below
                    // already decide it correctly: each emits a `dup` only where
                    // the alt equals an adjacent same-length reference tract, so
                    // the change genuinely IS a duplication. The spec's second
                    // example — four `TA` copies added to a `TA[2]` tract — has
                    // no such tract (the added bases are twice its length), so
                    // none of them fire and it falls through to
                    // `c.1741_1742insTATATATA` as prescribed.
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

                    // Check for simple duplication (single-base or matching adjacent)
                    // When shifting an insertion through a repeat region, the effective sequence
                    // rotates. For example, shifting "GGC" by 1 position gives "GCG".
                    //
                    // The rotation must be direction-aware. A 3' (rightward)
                    // shuffle moves the insertion point up (`result.start >=
                    // shuffle_start`) and rotates the alt LEFT by the shift `k`.
                    // A 5' (leftward) shuffle moves it down (`result.start <
                    // shuffle_start`) and must rotate the alt RIGHT by `k` —
                    // expressed here as an equivalent left-rotation of
                    // `(L - k mod L)` so the shared rotate-left below applies to
                    // both. The previous code computed the shift with
                    // `saturating_sub`, which clamped every leftward shift to 0,
                    // so 5' multi-base insertions were emitted UNROTATED at the
                    // shifted position — a different haplotype than the input,
                    // and non-idempotent (e.g. `g.260_261insCA` in an A-tract ->
                    // `g.259_260insCA` instead of the correct `g.259_260insAC`).
                    // Issue #1157 follow-up (idempotency proptest campaign).
                    //
                    // Empty-alt guard is defence in depth: `build_naedit` no
                    // longer merges a cancelling del+ins into an empty insertion
                    // (#1135), but the type still permits one (the parser accepts
                    // a bare `ins`), so nothing-to-rotate means no rotation.
                    let rotation = if seq_bytes.is_empty() {
                        0
                    } else if result.start >= shuffle_start {
                        (result.start as usize - shuffle_start as usize) % seq_bytes.len()
                    } else {
                        let k = (shuffle_start as usize - result.start as usize) % seq_bytes.len();
                        (seq_bytes.len() - k) % seq_bytes.len()
                    };
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

                    if five_prime_match || three_prime_match {
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
                        // 3'-side anchor is derived from `pos_idx` (the 0-based
                        // matched-tract start `ref[pos_idx..pos_idx+L]`), NOT from
                        // `new_start`: when a 5'-shuffled insertion comes to rest
                        // at the transcript start (`result.start == 0`), the
                        // insertion coordinate adjust-back `saturating_sub(1)`
                        // collapses `new_start` to c.1, breaking the usual
                        // `new_start == pos_idx + 1` relation and mislabelling the
                        // dup one position 3' (e.g. `c.2_3insGA` in `G[A..]` →
                        // `c.2_3dup` = a different haplotype instead of the correct
                        // `c.1_2dup`). `pos_idx` is unaffected by that clamp, and
                        // `index_to_hgvs_pos(pos_idx) == new_start + 1` whenever no
                        // saturation occurred, so this is identical off-boundary.
                        let three_prime_anchor = index_to_hgvs_pos(pos_idx);
                        let (dup_start, dup_end) = if prefer_three_prime {
                            (three_prime_anchor, three_prime_anchor + dup_len - 1)
                        } else if dup_len == 1 {
                            (new_start, new_start)
                        } else {
                            // `saturating_sub` for parity with the 3' arm above: that
                            // arm switched off `new_start` precisely because the
                            // insertion adjust-back clamp can collapse it to c.1. A
                            // 5'-side match needs a tract 5' of the insertion point,
                            // which cannot exist at the transcript start, so this is
                            // believed unreachable — but plain `-` on `u64` would
                            // panic in debug / wrap in release if it ever were, and
                            // the rest of this path guards the same way.
                            (new_start.saturating_sub(dup_len - 1).max(1), new_start)
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
                    }) = ins_to_dup.as_ref()
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
                if let Some(dup_result) = rules::duplication_to_repeat(
                    ref_seq,
                    result.start,
                    result.end,
                    // #1185: this call decides about the SHUFFLED tract
                    // (`result.*`), so the codon-frame gate is asked about that
                    // tract — NOT `gate.input_span_is_coding()`, which answers for
                    // the span the input edit occupied.
                    //
                    // A dup whose 3'-shifted tract runs out of the CDS into
                    // the 3'UTR is exempt from the multiple-of-3 rule:
                    // "This restriction only applies to the coding sequence,
                    // which does not include the introns or the UTR
                    // sequence." (DNA/repeated.md) — the added copies land
                    // past `cds_end` and so cannot move the reading frame.
                    // Keyed on the input span it was gated to an `ins`
                    // literal on pass 1 and collapsed to repeat notation only
                    // on pass 2, once the straddling tract had made the input
                    // span itself non-coding.
                    //
                    // `result.*` are 0-based; the conversion to the 1-based
                    // transcript frame `cds` lives in, and the inclusive-end
                    // `saturating_sub(1)`, stay here — they are this caller's
                    // frame, not the gate's.
                    gate.span_is_coding(
                        index_to_hgvs_pos(result.start as usize),
                        index_to_hgvs_pos(result.end.saturating_sub(1) as usize),
                    ),
                ) {
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
                                // The rule layer anchors the gated insertion at
                                // the tract's 3' flank regardless of shuffle
                                // direction. Re-enter `normalize_na_edit` on the
                                // derived insertion so it shuffles to the
                                // direction-appropriate boundary and reaches the
                                // same CDS-start/-end boundary clamps a directly
                                // written insertion does — otherwise a
                                // `dup`-spelled homopolymer edit and its
                                // insertion-spelled twin normalize to different
                                // (and, at a coding boundary, invalid) forms.
                                // Terminates: the derived edit is an Insertion,
                                // which never re-enters this Duplication arm.
                                // The `gate` is threaded through unchanged: the
                                // recursion stays on the same transcript, so the
                                // CDS bounds are identical, and the derived
                                // insertion covers a DIFFERENT span than the
                                // input — exactly the case the gate carries its
                                // bounds for (#1185). Before #1206 this was a
                                // `cds_span: Option<_>` argument and passing
                                // `None` here would have silently switched the
                                // codon-frame gate off for the derived insertion;
                                // that mistake is no longer expressible.
                                let (rec_start, rec_end, rec_edit, mut rec_warnings) = self
                                    .normalize_na_edit(
                                        ref_seq, &ins_edit, ins_start, ins_end, boundaries, gate,
                                    )?;
                                let mut merged = warnings;
                                merged.append(&mut rec_warnings);
                                return Ok((rec_start, rec_end, rec_edit, merged));
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
                        gate.input_span_is_coding(),
                        // #1270: as on the insertion path above, the verdict is
                        // the INPUT span's; the bounds let the helper re-ask
                        // about the tract the repeat would actually occupy.
                        gate.cds_bounds(),
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

    /// Convert CDS position to transcript position.
    ///
    /// **Sequence frame.** `c.N` is `cds_start + N - 1` on the transcript's own
    /// bases — the same flat arithmetic as
    /// [`crate::convert::mapper::CoordinateMapper::cds_to_tx`] — and that is
    /// what the normalizer needs, because the offsets it produces index the
    /// sequence `ReferenceProvider::get_sequence` serves for a transcript
    /// accession, and that sequence is flat. It deliberately does **not**
    /// consult the exon table.
    ///
    /// **Same frame, but not interchangeable at `c.0`.** The two implementations
    /// disagree on exactly one input, and it is a reachable one, so do not read
    /// "same arithmetic" as "call either". `c.0` is not a valid HGVS position;
    /// the legacy arm below maps it to `cds_start - 1` (i.e. `c.-1`, the formula
    /// continued), while `CoordinateMapper::cds_to_tx` folds `base < 1` into one
    /// branch and answers `cds_start` (i.e. `c.1`). Both arms are live —
    /// `merge::collapse_overlapping_cis_edits` can build a `base == 0` anchor,
    /// and `spdi::convert::resolve_cds_to_tx` routes the `c.?` sentinel
    /// (`CDS_BASE_UNKNOWN == 0`) through the mapper. Which one is right needs
    /// its own ruling; see #1772.
    ///
    /// A second, quieter divergence at the other end of the 5'UTR: for
    /// `pos.base < -(cds_start as i64)` — a `c.-N` further upstream than the
    /// transcript is long — this method **errors**
    /// (`u64::try_from` on a negative), while `CoordinateMapper::cds_to_tx`
    /// returns a `TxPos` with a *negative* `base`.
    /// `CoordinateMapper::fold_non_intronic_utr_offset` relies on being handed
    /// that negative base, so it is deliberate there rather than an oversight.
    ///
    /// Note also that neither of these pairs round-trips through `c.0`, because
    /// HGVS has no `c.0` to round-trip to: here `c.0` goes to `cds_start - 1`
    /// and [`Self::tx_to_cds_pos`] brings it back as `c.-1`; on the mapper,
    /// `c.0` goes to `cds_start` and comes back as `c.1`.
    ///
    /// It is therefore *not* interchangeable with the genome-frame conversions
    /// (`genomic_to_tx` / `tx_to_genomic` and the projector above them), which
    /// are exon- and CIGAR-aware. On a transcript whose alignment carries a
    /// transcript-coordinate gap the two families give different answers, by
    /// design; `tests/it/axis_frame_disagreement.rs` pins that. Making this one
    /// walk the exon table is #1619, adjudicated by PR #1735.
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

    /// Convert transcript position to CDS position.
    ///
    /// **Sequence frame**, the inverse of [`Self::cds_to_tx_pos`] and on the
    /// same flat axis: no exon walk, no CIGAR adjustment. See that method for
    /// why the genome-frame conversions are not a substitute — and for the one
    /// input on which the round trip does not close, `c.0`, which comes back as
    /// `c.-1` because HGVS has no `c.0` for it to come back as.
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
                special: None,
            })
        } else if pos > end {
            Ok(CdsPos {
                base: (pos - end) as i64,
                offset: None,
                utr3: true,
                special: None,
            })
        } else {
            Ok(CdsPos {
                base: (pos - cds_start + 1) as i64,
                offset: None,
                utr3: false,
                special: None,
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
    /// - Inversion priority: a delins whose whole *maximal contiguous run*
    ///   is a reverse complement splits out that run as `inv` (issue #160,
    ///   corrected by issue #1034). A reverse-complement *sub-run* of a
    ///   longer contiguous change is NOT carved out — that change stays a
    ///   single `delins`.
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
        // The `r.` branch of `fetch_ref_for_canonical_split` translates
        // the RNA axis through `cds_start` (`r.N` → 1-based tx pos
        // `cds_start + N - 1`, falling back to transcript-1 when the
        // transcript has no CDS), so its ref window is CDS-aligned and the
        // codon-frame split path reads correctly for `cds_start > 1` (#469,
        // superseding the latent #291 mis-alignment).
        //
        // The `r.` arm asks the transcript, not the discriminant. An `r.`
        // description is an *RNA* description, not necessarily a *coding* one:
        // on a non-coding transcript there is no reading frame, so there is no
        // codon for the exception to preserve. Keying off `HgvsVariant::Rna`
        // alone stamped a frame onto `NR_` transcripts and re-merged triplets
        // the separation rule (`general.md:34`) says must stay apart — #1241,
        // where the `n.` and `r.` spellings of one non-coding change settled on
        // two different strings. `Cds` needs no such check: `c.` positions are
        // CDS-relative by construction, and the `Cds` arm of
        // `fetch_ref_for_canonical_split` has already refused a transcript with
        // no `cds_start`.
        let codon_frame_aware = match &variant {
            HgvsVariant::Cds(_) => true,
            HgvsVariant::Rna(r) => self
                .provider
                .get_transcript(&r.accession.transcript_accession())
                .ok()
                .and_then(|tx| tx.cds_start)
                .is_some(),
            _ => false,
        };
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
                // `r.` shares c.'s CDS-relative numbering on a coding
                // transcript (#469): `r.N` maps to 1-based tx pos
                // `cds_start + N - 1`, exactly like the Cds arm above.
                // `simple_rna_pos` has already filtered UTR (`r.*N`/`r.-N`,
                // now a negative base) and intronic positions, so any position
                // reaching this arm is a positive CDS-relative base. Without a
                // CDS (coordinate-only / mock) fall back to transcript-1.
                let tx = self
                    .provider
                    .get_transcript(&r.accession.transcript_accession())
                    .ok()?;
                let bytes = tx.sequence.as_deref()?.as_bytes();
                let (s, e) = match tx.cds_start {
                    Some(cds_start) => (
                        cds_start.checked_add(hgvs_start)?.checked_sub(2)? as usize,
                        cds_start.checked_add(hgvs_end)?.checked_sub(1)? as usize,
                    ),
                    None => ((hgvs_start - 1) as usize, hgvs_end as usize),
                };
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
    fn canonical_split_for_variant(
        &self,
        v: HgvsVariant,
        manufactured: &mut Vec<ManufacturedJunctionExit>,
    ) -> Vec<HgvsVariant> {
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
        // `normalize_core`, deliberately, NOT `normalize_with_diagnostics`.
        //
        // This helper runs *inside* `normalize_allele`, which runs inside
        // `normalize_core`. Since #1237 the public exits route through
        // `normalize_core_canonical`, which runs the sequence-first pass — and
        // that pass calls `normalize_core` on its own re-derivation. Going
        // through a public exit from here therefore closes the loop
        // `normalize_core → normalize_allele → canonical_split_for_variant →
        // normalize_core_canonical → normalize_core`, which recurses until the
        // stack is gone. That is exactly the recursion `normalize_core_canonical`
        // documents as the reason the pass cannot live in `normalize_core`; the
        // rule is the same here, one level down.
        //
        // The switch also restores this helper's pre-#1237 meaning. It wants the
        // per-member pipeline's split of a simple delins, which is what
        // `normalize_core` is; and it discards everything the diagnostics
        // wrapper adds, so the `detect_shuffle_infos` pass it used to run was
        // wasted work on every call.
        match self.normalize_core(&v, manufactured) {
            Ok((result, _warnings)) => match result {
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
/// forward-safe for when circular shuffling is implemented (tracked by #466's
/// circular candidate; the earlier pointers at #951 and, before it, the
/// MT-scoped #129 are both closed) and costs nothing today (the post-hoc
/// equality check trivially returns no info).
fn position_text_if_shuffleable(variant: &HgvsVariant) -> Option<String> {
    match variant {
        HV::Genome(v) => Some(format!("{}", v.loc_edit.location)),
        HV::Cds(v) => Some(format!("{}", v.loc_edit.location)),
        HV::Tx(v) => Some(format!("{}", v.loc_edit.location)),
        HV::Rna(v) => Some(format!("{}", v.loc_edit.location)),
        HV::Mt(v) => Some(format!("{}", v.loc_edit.location)),
        HV::Circular(v) => Some(format!("{}", v.loc_edit.location)),
        HV::Protein(_)
        | HV::Allele(_)
        | HV::RnaFusion(_)
        | HV::GenomeRing(_)
        | HV::Supernumerary(_)
        | HV::NullAllele
        | HV::UnknownAllele => None,
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
        HV::Protein(_)
        | HV::Allele(_)
        | HV::RnaFusion(_)
        | HV::GenomeRing(_)
        | HV::Supernumerary(_)
        | HV::NullAllele
        | HV::UnknownAllele => return None,
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

/// Maximum genomic span fetched for intronic shuffle normalization. Sizing
/// the window to the enclosing intron (issue #573) lets the minus-strand
/// boundary check see the intron edges; a pathologically large intron is
/// capped here, in which case the variant-sized window is used and the
/// downstream guard returns a clean `Err` rather than fetching megabases.
const MAX_INTRONIC_SHUFFLE_WINDOW: u64 = 64 * 1024;

/// Genomic `[seq_start, seq_end]` window for intronic normalization, sized to
/// cover both the variant (± `window`) and the enclosing intron boundaries so
/// the 3'-shift boundary check has the intron edges in-window. Falls back to
/// the variant-sized window when the intron span exceeds
/// [`MAX_INTRONIC_SHUFFLE_WINDOW`].
fn intronic_window_bounds(
    g_start: u64,
    g_end: u64,
    intron_g_start: u64,
    intron_g_end: u64,
    window: u64,
) -> (u64, u64) {
    let var_start = g_start.saturating_sub(window);
    let var_end = g_end.saturating_add(window);
    let want_start = var_start.min(intron_g_start);
    // The fetch is end-exclusive (relative length = `seq_end - seq_start`), so
    // the far intron edge `intron_g_end` only lands inside the 1-based window
    // when `seq_end > intron_g_end`. Extend one past it.
    let want_end = var_end.max(intron_g_end.saturating_add(1));
    // Clamp the window start to the 1-based minimum. Near the start of a contig
    // `saturating_sub(window)` can drive the start to 0, which has no valid
    // 1-based coordinate and would break the caller's
    // `rel = g - seq_start + 1` arithmetic (and is rejected by
    // `get_genomic_sequence_1based`). Clamping to 1 keeps both correct.
    if want_end.saturating_sub(want_start) <= MAX_INTRONIC_SHUFFLE_WINDOW {
        (want_start.max(1), want_end)
    } else {
        (var_start.max(1), var_end)
    }
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
// - Inversion: a whole maximal contiguous run within the delins span is a
//   reverse complement — split that run out as `inv` (issue #160, corrected
//   by issue #1034; a rev-comp sub-run of a longer contiguous change is not
//   carved out).
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
///    Such a triplet contributes the three positions
///    `[Sub@i.alt, Identity@i+1.base, Sub@i+2.alt]` to the pending run,
///    which renders as a 3-base `delins` when nothing adjoins it. The flag
///    is true only for `c.` (CDS) variants — `g.`, `n.`, `r.`, and `m.`
///    have no codon-frame and skip this branch. The exception is
///    deliberately narrow (length-3, exact pattern, in-codon endpoints)
///    so it matches the spec text "two variants separated by one
///    nucleotide, together affecting one amino acid".
///
/// 2. **Adjacent-change coalescence** (`delins.md:16`, issue #182,
///    issue #1524). Sub-edits occupying strictly adjacent positions —
///    no gap, no `Inversion` between them — group into a single `delins`
///    variant: "changes involving two or more consecutive nucleotides are
///    described as deletion/insertion (delins) variants".
///
///    The codon-frame triplet is inside this rule, not beside it, and that
///    is what #1524 fixed. The triplet used to be emitted as its own
///    member, which put two members on consecutive nucleotides whenever
///    `i - 1` or `i + 3` was also changed. Both edges are now closed, and
///    differently: on the right the triplet stays in the run so `i + 3`
///    joins it, while on the left rule 1 **declines** — a changed `i - 1`
///    means `i` is not one of `general.md:35`'s "two variants" at all.
///
/// 3. **Inversion as a hard barrier** (issue #166). An `Inversion`
///    always emits standalone and breaks any in-flight run,
///    preserving the inv-priority decomposition. It can never *be*
///    adjacent to another member: `decompose_delins` emits one only for a
///    whole **maximal** contiguous mismatch run, so a neighbouring
///    mismatch would have been part of the same run.
///
/// Singleton sub-runs stay as `Substitution`. `IdentityAt` not consumed
/// by a codon-frame triplet drops (an unchanged base is not an edit) and
/// always ends any in-flight run — the gap means the
/// surrounding changes are no longer "consecutive".
fn build_split_variants(
    template: &HgvsVariant,
    subedits: Vec<DelinsSubedit>,
    hgvs_start: u64,
    codon_frame_aware: bool,
) -> Vec<HgvsVariant> {
    let abs = |idx: usize| -> u64 { idx as u64 + hgvs_start };

    let mut output: Vec<HgvsVariant> = Vec::new();
    // Pending run of strictly-adjacent sub-edit positions, in left-to-right
    // order. Each entry is `(position, reference, alternative)` with
    // `position` the 0-indexed offset into the variant's ref window.
    //
    // Entries are mismatches except for the unchanged centre a codon-frame
    // triplet contributes, which enters as `reference == alternative`
    // (#1524). A run therefore always *begins* and *ends* on a mismatch,
    // which is the precondition `push_typed_replacement` relies on.
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
                    reference: r1,
                    alternative: a1,
                },
                DelinsSubedit::IdentityAt {
                    position: pm,
                    base: bm,
                },
                DelinsSubedit::Substitution {
                    position: p3,
                    reference: r3,
                    alternative: a3,
                },
            ) = (&subedits[i], &subedits[i + 1], &subedits[i + 2])
            {
                // `general.md:35` licenses the merge for "two **variants**
                // separated by one nucleotide". When `p1 - 1` is itself changed,
                // `p1` is not a variant: `delins.md:16` makes it part of the
                // `delins` spanning the run it sits in, and the thing separated
                // from `p3` by one nucleotide is that whole run — which reaches
                // beyond the one codon the exception is about. So the pattern
                // the exception describes is not present and the branch declines
                // (#1524).
                //
                // This is the same precondition `merge::apply_coding_codon_exception`
                // has always enforced on its own side of the seam
                // (`is_substitution(&pieces[index - 1])` — the left piece must
                // be a *lone* substitution). The two are implementations of one
                // rule and this one was missing the test, so adding it removes a
                // disagreement between the paths rather than imposing a new
                // restriction.
                //
                // It is also what keeps the joining below from over-merging, and
                // that was measured rather than reasoned: without this guard,
                // `c.9_13delinsACAAC` on `TTTTTTTTTAATATAT…`
                // (`[Sub@9; Sub@10; Identity@11; Sub@12; Sub@13]`, codon 4 =
                // 10-12) collapses to a single `c.9_13delinsACAAC` spanning the
                // unchanged 11 — where the two real variants, `9_10delinsAC` and
                // `12_13delinsAC`, span three codons and are exactly what
                // `general.md:34` says to describe individually. That form cost
                // 44 classes of `cis_confluence_axis`'s 3' census
                // (converged 6633 -> 6589) and 44 of its 5' census, because the
                // multi-member spelling of each still split.
                let left_endpoint_is_a_lone_variant =
                    !matches!(run.last(), Some((previous, _, _)) if *previous + 1 == *p1);
                if *pm == *p1 + 1 && *p3 == *p1 + 2 && left_endpoint_is_a_lone_variant {
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
                        // The triplet **enters the pending run** instead of
                        // being emitted beside it (#1524).
                        //
                        // The run is empty or separated here — the guard above
                        // has already declined the adjacent-left case — so this
                        // is not about the left edge. It is about the right one.
                        // Emitting the triplet directly ended the member at `p3`
                        // and started the next one at `p3 + 1` whenever the
                        // following position was also changed, which puts two
                        // members on *consecutive* nucleotides with nothing
                        // unchanged between them. `delins.md:16` forbids exactly
                        // that — "changes involving two or more consecutive
                        // nucleotides are described as deletion/insertion
                        // (delins) variants" — and `general.md:34` governs
                        // members "separated by one or more" nucleotides, so it
                        // does not reach separation 0: there is no competing
                        // clause and nothing to trade off. Leaving the triplet in
                        // the run lets the ordinary adjacency test below absorb
                        // `p3 + 1`, and the member comes out as one `delins`.
                        //
                        // Measured against the prepared reference:
                        // `NM_000083.3:c.2461_2464delinsCTCC` decomposed to
                        // `[Sub@2461; Identity@2462; Sub@2463; Sub@2464]`, the
                        // triplet fired at 2461 because 2461-2463 is one codon,
                        // and the output was `c.[2461_2463delinsCTC;2464G>C]` —
                        // 2463 and 2464 consecutive, in two members. It is now
                        // `c.2461_2464delinsCTCC`, the input's own spelling: this
                        // is one of several rows the corpus audit found where the
                        // input was already correct and ferro made it wrong.
                        //
                        // The unchanged centre enters the run as its own
                        // ref-equals-alt entry, which is what lets a run carry an
                        // interior identity at all. `flush_substitution_run`
                        // re-emits the whole run through
                        // `push_typed_replacement`, so the triplet alone still
                        // renders exactly as it did — the span `[p1, p3]` with
                        // alt `[a1, bm, a3]`, typed rather than assumed to be a
                        // `delins`. That typing stays correct for a joined run
                        // too: both ends are still mismatches, so `Identity`,
                        // `Substitution`, `Deletion`, `Insertion` and
                        // `Duplication` remain unreachable and only `Inversion`
                        // can come back, exactly as `push_typed_replacement`'s
                        // own note argues. (The triplet reaches an inversion only
                        // when the untouched middle base is its own complement —
                        // `W`, `S`, `N` — which is why the reference in
                        // `codon_triplet_over_an_ambiguous_centre_is_an_inversion`
                        // carries an `N`.)
                        //
                        // The flush is still unconditional: a run that ends at
                        // `p1 - 2` or earlier is separated from the triplet by at
                        // least one unchanged nucleotide, and `general.md:34`
                        // keeps it a member of its own. That boundary is the
                        // discriminating case, pinned by
                        // `issue_1524_adjacent_split_members::
                        // a_run_separated_from_the_triplet_still_splits`.
                        flush_substitution_run(&mut output, template, hgvs_start, &mut run);
                        run.push((*p1, *r1, *a1));
                        run.push((*pm, *bm, *bm));
                        run.push((*p3, *r3, *a3));
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

/// Emit one split member spanning `[start_1based, end_1based]`, typed under the
/// spec's edit priority instead of assumed to be a `delins` (#1454).
///
/// `build_split_variants` re-groups the sub-edits `decompose_delins` reports,
/// and its multi-base grouping — the run of strictly adjacent positions, which
/// since #1524 absorbs the codon-frame triplet rather than sitting beside it —
/// **forms spans that did not exist when `decompose_delins` typed the input**.
/// The inversion test it ran was against
/// the input's own maximal contiguous mismatch runs, and it deliberately does
/// not carve a reverse-complement *sub*-run out of a longer contiguous change
/// (#1034, #1040). But once a re-grouping has made such a sub-run a member in
/// its own right, `general.md:56` applies to it directly: inversion outranks
/// deletion-insertion, and `delins.md:5` defines a delins as a replacement
/// "which is not a substitution or inversion", so `delins` over a
/// reverse-complement span is not a description the spec admits at all. The
/// same argument `merge::coalesce_whole_block_inversion` makes for the
/// derivation's pieces, made here for the split's members.
///
/// Before this, that typing arrived one normalization pass late. The member was
/// emitted as a `delins`, and the *next* pass's per-member pipeline re-typed it
/// through `normalize_na_edit`'s Delins arm — the same
/// [`rules::canonicalize_delins`] called here — so the first output was not a
/// fixed point:
///
/// ```text
/// NM_TEST.1:c.10_15delinsTAATAT
///   pass 1 -> NM_TEST.1:c.[10_12delinsTAA;13_15delinsTAT]
///   pass 2 -> NM_TEST.1:c.[10_12delinsTAA;13_15inv]
/// ```
///
/// The typing authority is [`rules::canonicalize_delins`] rather than an
/// open-coded reverse-complement test, precisely so this pass and the next one
/// cannot disagree: reaching the fixed point in one pass means answering the
/// question the second pass would ask, with the function it would ask it of.
///
/// Only `Inversion` can come back for these inputs, which is what makes the
/// narrow match sound rather than lossy. Both groupings hand over equal-length
/// slices whose first and last positions are mismatches, so `Identity`
/// (needs every position equal), `Duplication` (needs `alt.len() == 2 *
/// ref.len()`), `Insertion`/`Deletion` (need an empty side) and `Substitution`
/// (needs a length-1 post-trim range, but greedy shared-affix trimming cannot
/// consume a mismatched end) are all unreachable. Anything else is therefore a
/// change in that function's contract, and falling through to `delins` — the
/// form this code emitted unconditionally before — is the safe reading of it.
fn push_typed_replacement(
    output: &mut Vec<HgvsVariant>,
    template: &HgvsVariant,
    start_1based: u64,
    end_1based: u64,
    ref_bases: &[Base],
    alt_bases: Vec<Base>,
) {
    debug_assert_eq!(
        ref_bases.len(),
        alt_bases.len(),
        "a split member replaces its span base for base"
    );
    let edit = if spans_a_whole_inversion(ref_bases, &alt_bases) {
        // Canonical form states neither the inverted bases nor the length
        // (`inversion.md`, #352), matching `canonicalize_edit`'s own arm.
        NaEdit::Inversion {
            sequence: None,
            length: None,
        }
    } else {
        NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::new(alt_bases)),
            deleted: None,
            deleted_length: None,
            substitution_reference: None,
        }
    };
    output.push(build_variant_at(template, start_1based, end_1based, edit));
}

/// Whether replacing `ref_bases` with `alt_bases` inverts the **whole** span.
///
/// `U` is folded to `T` on both sides before the comparison, exactly as
/// `apply_canonical_split` does before calling `decompose_delins`: on the `r.`
/// axis the alt bases come from the author's literal (so they may be `U`) while
/// the transcript reference is `T`, and [`rules::is_revcomp`]'s complement table
/// maps `A` to `T`, never to `U`. Without the fold an `r.` inversion reads as an
/// ordinary delins and the axis silently keeps the defect this closes.
///
/// A partial match is rejected rather than emitted narrower: a shortened range
/// would leave the peeled outer bases undescribed, and this function's callers
/// have no way to emit them. It cannot arise for their inputs — a peelable
/// outer pair means `alt[first] == ref[first]`, i.e. an identity where both
/// callers guarantee a mismatch — so the guard costs nothing and is the
/// conservative reading if that ever stops holding.
fn spans_a_whole_inversion(ref_bases: &[Base], alt_bases: &[Base]) -> bool {
    // `inversion.md` requires more than one nucleotide; a one-base complement is
    // a substitution, which `flush_substitution_run` already emits directly.
    if ref_bases.len() < 2 || ref_bases.len() != alt_bases.len() {
        return false;
    }
    let to_bytes =
        |bases: &[Base]| normalize_t_u(&bases.iter().map(|base| *base as u8).collect::<Vec<u8>>());
    let ref_bytes = to_bytes(ref_bases);
    let alt_bytes = to_bytes(alt_bases);
    matches!(
        rules::canonicalize_delins(&ref_bytes, 0, ref_bytes.len(), &alt_bytes),
        rules::DelinsCanonical::Inversion { start, end } if start == 0 && end == ref_bytes.len()
    )
}

/// Flush a pending run of consecutive sub-edit positions into `output`.
/// A length-1 run emits a `Substitution`; a length-2+ run emits a single
/// `Delins` over `[run.first.position, run.last.position]` with `sequence`
/// = concatenated `alternative` bases — or an `Inversion` when that span is a
/// whole reverse complement, see [`push_typed_replacement`]. See
/// `build_split_variants` for the spec rationale (issue #182, issue #1524).
///
/// A length-1 run is always a genuine substitution: the only run entry whose
/// `reference` equals its `alternative` is the unchanged centre of a
/// codon-frame triplet, and that arrives flanked by the triplet's two
/// mismatches, so it can never be alone in the run.
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
    let (ref_bases, alt_bases): (Vec<Base>, Vec<Base>) =
        run.drain(..).map(|(_, r, a)| (r, a)).unzip();
    push_typed_replacement(output, template, s, e, &ref_bases, alt_bases);
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

/// Whether `id` is a versioned **transcript** accession subject to the
/// base→latest version substitution the #785 gate guards against.
///
/// RefSeq (`NM_`/`NR_`/`XM_`/`XR_`) and Ensembl (`ENST`) transcripts qualify.
/// Genomic references (`NG_`/`NC_`/`NW_`) and Ensembl *gene* identifiers
/// (`ENSG`) do **not**: they are not transcript reading frames, and a `c.`/`n.`
/// input that resolves to one of those (e.g. an unrewritten gene-model selector
/// `NG_(GENE):c.…`) must not be gated as if it named a transcript version.
fn is_transcript_accession(id: &str) -> bool {
    id.starts_with("NM_")
        || id.starts_with("NR_")
        || id.starts_with("XM_")
        || id.starts_with("XR_")
        || id.starts_with("ENST")
}

/// If `variants` has 1 element return it directly; if >1 wrap in a cis
/// Allele. `uncertain` is the pre-split variant's edit certainty (#1063):
/// a single predicted/uncertain edit that canonically decomposes into a
/// multi-member allele carries its `(...)` wrapper onto the whole allele
/// (`g.[(200C>T;202C>G)]`), never onto the individual members — members
/// built by `build_variant_at` are always `Mu::Certain`. Ignored when the
/// split doesn't fire (single-element `variants`).
fn wrap_allele_if_split(mut variants: Vec<HgvsVariant>, uncertain: bool) -> HgvsVariant {
    debug_assert!(!variants.is_empty(), "wrap_allele_if_split: empty input");
    if variants.len() == 1 {
        variants.pop().unwrap()
    } else if uncertain {
        HgvsVariant::Allele(AlleleVariant::new_uncertain(variants, AllelePhase::Cis))
    } else {
        HgvsVariant::Allele(AlleleVariant::new(variants, AllelePhase::Cis))
    }
}

#[cfg(test)]
mod codon_gate_tests {
    use super::CodonGate;

    /// A record with no CDS end yields the non-applicable gate, matching the
    /// pre-#1206 `(is_coding = false, cds_span = None)` pairing: without bounds we
    /// cannot verify anything lies inside the CDS proper, so we do not gate.
    #[test]
    fn no_cds_bounds_is_not_applicable() {
        assert_eq!(
            CodonGate::for_input_span(None, 10, 20),
            CodonGate::NotApplicable
        );
    }

    /// The input-span verdict is the inclusive containment test both
    /// `normalize_cds` and `normalize_rna` used to compute separately.
    #[test]
    fn input_span_verdict_is_inclusive_containment() {
        let cds = Some((10, 20));
        // wholly inside, including flush against either bound
        assert!(CodonGate::for_input_span(cds, 10, 20).input_span_is_coding());
        assert!(CodonGate::for_input_span(cds, 12, 15).input_span_is_coding());
        // one base into either UTR
        assert!(!CodonGate::for_input_span(cds, 9, 20).input_span_is_coding());
        assert!(!CodonGate::for_input_span(cds, 10, 21).input_span_is_coding());
    }

    /// The gate keeps its bounds, so a site deciding about a span the input edit
    /// does not occupy can re-answer for that span (#1185) — the capability that
    /// the old `Option` argument made optional, and therefore losable.
    #[test]
    fn a_coding_gate_can_re_answer_for_another_span() {
        // Input span inside the CDS, so the input verdict is `true`…
        let gate = CodonGate::for_input_span(Some((10, 20)), 12, 13);
        assert!(gate.input_span_is_coding());
        // …but a tract the shuffle walked into the 3'UTR is exempt.
        assert!(gate.span_is_coding(12, 20));
        assert!(!gate.span_is_coding(12, 21));
        assert!(!gate.span_is_coding(9, 12));
    }

    /// Both questions are answered "not coding" when there is no reading frame,
    /// so a `NotApplicable` axis can never be gated by either accessor.
    #[test]
    fn not_applicable_is_never_coding() {
        assert!(!CodonGate::NotApplicable.input_span_is_coding());
        assert!(!CodonGate::NotApplicable.span_is_coding(1, 1));
        assert!(!CodonGate::NotApplicable.span_is_coding(u64::MIN, u64::MAX));
    }

    /// The property the type exists for: a coding gate always carries bounds.
    ///
    /// This cannot fail at runtime — `Coding` has no constructor that omits `cds`,
    /// which is the whole point — so what this really pins is that the
    /// constructor never produces the pairing the old two-argument form allowed
    /// ("coding verdict, no bounds"). If a future change adds a way to build one,
    /// this is where it should be caught.
    #[test]
    fn a_coding_gate_always_has_bounds() {
        for (cds, start, end) in [
            (Some((1, 3)), 1, 3),
            (Some((1, 3)), 4, 9),
            (Some((100, 100)), 100, 100),
        ] {
            match CodonGate::for_input_span(cds, start, end) {
                CodonGate::Coding { cds: bounds, .. } => assert_eq!(Some(bounds), cds),
                CodonGate::NotApplicable => panic!("{cds:?} should have produced a coding gate"),
            }
        }
    }
}

#[cfg(test)]
mod sequence_bounds_clamp_tests {
    use super::*;
    use crate::hgvs::edit::{Base, InsertedSequence, Sequence};

    fn literal_insertion(bases: &str) -> NaEdit {
        let bases: Vec<Base> = bases
            .chars()
            .map(|c| Base::from_char(c).expect("test bases must be A/C/G/T"))
            .collect();
        NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::new(bases)),
        }
    }

    /// The 5'-saturated shape is rewritten when the slice's first base really is
    /// the entity's first base.
    #[test]
    fn clamps_at_a_real_five_prime_bound() {
        let (mut edit, mut start, mut end) = (literal_insertion("CG"), 1, 1);
        clamp_insertion_at_sequence_bounds(
            b"CCCCAAAA",
            &mut edit,
            &mut start,
            &mut end,
            SequenceEnds::WHOLE,
        );
        assert!(matches!(edit, NaEdit::Delins { .. }), "{edit:?}");
        assert_eq!((start, end), (1, 1));
    }

    /// …and NOT when that end of the slice is merely where a fetch stopped.
    ///
    /// This is the guard `normalize_genome` needs and the transcript callers do
    /// not: it passes a window into the contig, so an interior window edge is not
    /// a contig bound, and clamping there would rewrite a valid insertion into a
    /// `delins` asserting a boundary that does not exist. Unit-tested directly
    /// because the end-to-end route to it is not currently constructible — on the
    /// `g.` axis an alt has to be tandem-compatible with the reference to shuffle
    /// at all, and such an alt becomes repeat notation or a `dup` long before it
    /// could walk a whole window's width. The guard is therefore correctness by
    /// construction rather than a fix for an observed output.
    #[test]
    fn does_not_clamp_at_a_window_edge_that_is_not_a_bound() {
        for ends in [
            SequenceEnds {
                five_prime: false,
                three_prime: true,
            },
            SequenceEnds {
                five_prime: false,
                three_prime: false,
            },
        ] {
            let (mut edit, mut start, mut end) = (literal_insertion("CG"), 1, 1);
            clamp_insertion_at_sequence_bounds(b"CCCCAAAA", &mut edit, &mut start, &mut end, ends);
            assert!(
                matches!(edit, NaEdit::Insertion { .. }),
                "{ends:?}: {edit:?}"
            );
            assert_eq!((start, end), (1, 1));
        }
    }

    /// The 3' side, both ways round.
    #[test]
    fn three_prime_bound_is_gated_the_same_way() {
        let seq = b"CCCCAAAA"; // len 8, so the saturated shape is (8, 9)
        let (mut edit, mut start, mut end) = (literal_insertion("CG"), 8, 9);
        clamp_insertion_at_sequence_bounds(
            seq,
            &mut edit,
            &mut start,
            &mut end,
            SequenceEnds::WHOLE,
        );
        assert!(matches!(edit, NaEdit::Delins { .. }), "{edit:?}");
        assert_eq!((start, end), (8, 8));

        let (mut edit, mut start, mut end) = (literal_insertion("CG"), 8, 9);
        clamp_insertion_at_sequence_bounds(
            seq,
            &mut edit,
            &mut start,
            &mut end,
            SequenceEnds {
                five_prime: true,
                three_prime: false,
            },
        );
        assert!(matches!(edit, NaEdit::Insertion { .. }), "{edit:?}");
        assert_eq!((start, end), (8, 9));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;
    use crate::reference::MockProvider;

    // -----------------------------------------------------------------------
    // #1723: the junction-exit fold, exercised directly
    //
    // `reparent_junction_exit` is driven at the seam by whatever the `#670`
    // gates recorded, and two of its branches are not reachable from a fixture
    // provider — a record whose leaf a later sibling pass moved, and a record
    // for an accession that is no longer in the output at all. Constructing the
    // records here is the only way to exercise those, and it is honest about
    // which they are: the end-to-end shapes live in
    // `tests/it/defect_371_transcript_exit.rs`.
    // -----------------------------------------------------------------------

    /// A normalizer over an empty provider. The fold never consults the provider
    /// — that is the point of carrying the contig from the gate — so an empty one
    /// is sufficient and proves the absence of a second lookup.
    fn fold_only_normalizer() -> Normalizer<MockProvider> {
        Normalizer::new(MockProvider::new())
    }

    fn manufactured_record(output_leaf: &str) -> ManufacturedJunctionExit {
        let leaf = parse_hgvs(output_leaf).expect("fixture leaf parses");
        let bare = leaf.accession().expect("a leaf has an accession").clone();
        let wrapper = genomic_wrapper_for(&bare, "NC_SYNTH.1").expect("a usable contig");
        ManufacturedJunctionExit {
            leaf,
            bare,
            wrapper,
        }
    }

    /// The empty fast path: no gate fired, so the output is returned untouched
    /// and no accession is even looked at.
    #[test]
    fn the_junction_exit_fold_is_a_no_op_with_no_records() {
        let output = parse_hgvs("NM_TEST.1:c.10+2del").unwrap();
        let folded = fold_only_normalizer().reparent_junction_exit(output.clone(), &[]);
        assert_eq!(folded, output, "no record means no repair");
    }

    /// The recorded leaf is the one in the output: it is repaired.
    #[test]
    fn the_junction_exit_fold_repairs_a_recorded_leaf() {
        let records = [manufactured_record("NM_TEST.1:c.10+2del")];
        let output = parse_hgvs("NM_TEST.1:c.10+2del").unwrap();
        assert_eq!(
            fold_only_normalizer()
                .reparent_junction_exit(output, &records)
                .to_string(),
            "NC_SYNTH.1(NM_TEST.1):c.10+2del"
        );
    }

    /// **The moved-leaf branch.** A record survives to the seam but the leaf it
    /// names no longer matches anything in the output — a sibling pass shifted
    /// it. There is then no leaf whose provenance is known, so the accession
    /// declines rather than being repaired on the strength of a record that no
    /// longer describes it.
    ///
    /// Not reachable from a fixture provider, which is why this is a unit test.
    /// Dropping the `all(...)` conjunct turns this assertion into
    /// `NC_SYNTH.1(NM_TEST.1):c.10+9del`.
    #[test]
    fn the_junction_exit_fold_declines_a_leaf_no_record_explains() {
        let records = [manufactured_record("NM_TEST.1:c.10+2del")];
        let output = parse_hgvs("NM_TEST.1:c.10+9del").unwrap();
        assert_eq!(
            fold_only_normalizer()
                .reparent_junction_exit(output, &records)
                .to_string(),
            "NM_TEST.1:c.10+9del",
            "a leaf that no record explains is left exactly as it stands"
        );
    }

    /// **The empty-leaf-set branch.** A record survives to the seam but nothing
    /// in the output names a bare intronic position on that accession any more.
    /// `all()` over an empty set is vacuously true, so without the explicit
    /// `!leaves.is_empty()` the exonic leaves that remain would be wrapped for no
    /// reason. Dropping that conjunct turns this assertion into
    /// `NC_SYNTH.1(NM_TEST.1):c.18del`.
    #[test]
    fn the_junction_exit_fold_wraps_nothing_when_no_intronic_leaf_survives() {
        let records = [manufactured_record("NM_TEST.1:c.10+2del")];
        let output = parse_hgvs("NM_TEST.1:c.18del").unwrap();
        assert_eq!(
            fold_only_normalizer()
                .reparent_junction_exit(output, &records)
                .to_string(),
            "NM_TEST.1:c.18del",
            "an exonic-only output is not wrapped merely because a gate once fired"
        );
    }

    /// **The mixed case, isolated from the pipeline.** One recorded leaf and one
    /// the author spelled, on one accession: the accession declines as a whole.
    /// This is the residue the `undecided`
    /// `junction-exit-wrapper-scope-in-a-mixed-allele` ruling records, and
    /// pinning it here as well as end-to-end separates "the policy is decline"
    /// from "the pipeline happened not to reach it".
    #[test]
    fn the_junction_exit_fold_declines_a_mixed_accession() {
        let records = [manufactured_record("NM_TEST.1:c.10+2del")];
        let output = parse_hgvs("NM_TEST.1:c.[10+2del;30+5del]").unwrap();
        assert_eq!(
            fold_only_normalizer()
                .reparent_junction_exit(output, &records)
                .to_string(),
            "NM_TEST.1:c.[10+2del;30+5del]",
            "an authored sibling on the same accession is the undecided case; decline"
        );
    }

    /// **Two accessions, both recorded, both repaired.** `#1704` matched on
    /// accession equality with no loop after taking the *first* offending leaf,
    /// so the second accession was never reached. Reverting the fold to a single
    /// accession leaves `NM_OTHER.1:c.10+2del` bare here.
    #[test]
    fn the_junction_exit_fold_repairs_every_recorded_accession() {
        let records = [
            manufactured_record("NM_TEST.1:c.10+2del"),
            manufactured_record("NM_OTHER.1:c.10+2del"),
        ];
        let output = parse_hgvs("[NM_TEST.1:c.10+2del;NM_OTHER.1:c.10+2del]").unwrap();
        assert_eq!(
            fold_only_normalizer()
                .reparent_junction_exit(output, &records)
                .to_string(),
            "[NC_SYNTH.1(NM_TEST.1):c.10+2del;NC_SYNTH.1(NM_OTHER.1):c.10+2del]"
        );
    }

    /// The wrapper builder's own declines, on the strings that reach it. Each
    /// case names the conjunct that refuses it — the SAM-style and assembly names
    /// die on the **bare parse**, not on the pairing rule, which is the half a
    /// reader attributes wrongly.
    #[test]
    fn the_genomic_wrapper_builder_declines_what_it_cannot_justify() {
        let bare = parse_hgvs("NM_TEST.1:c.10+2del")
            .unwrap()
            .accession()
            .unwrap()
            .clone();
        assert_eq!(
            genomic_wrapper_for(&bare, "NC_SYNTH.1")
                .expect("a genomic accession is usable")
                .to_string(),
            "NC_SYNTH.1(NM_TEST.1)"
        );
        for refused in [
            "chr17",          // bare parse
            "GRCh38",         // bare parse
            "NC_SYNTH.1junk", // trailing remainder
        ] {
            assert!(
                genomic_wrapper_for(&bare, refused).is_none(),
                "{refused:?} must not become an outer reference"
            );
        }
        // The pairing rule, which is a different mechanism from the two above.
        assert!(
            crate::hgvs::parser::accession::parse_accession("NM_OTHER.1").is_ok(),
            "the bare parse SUCCEEDS here, so the decline below is the pairing rule"
        );
        assert!(
            genomic_wrapper_for(&bare, "NM_OTHER.1").is_none(),
            "a transcript named as the outer reference is backwards"
        );

        // And the third conjunct, which is a THIRD mechanism and is reachable
        // only through LRG. `a_self_referential_wrapper_is_declined` pins it
        // end-to-end; this pins the two premises that make it reachable, so a
        // future reader can see why `outer == *bare` is not subsumed by the
        // pairing rule rather than having to re-derive it.
        let lrg = parse_hgvs("LRG_5:c.10+2del")
            .expect("a bare LRG record parses")
            .accession()
            .expect("and carries an accession")
            .clone();
        let (rest, outer) = crate::hgvs::parser::accession::parse_accession("LRG_5")
            .expect("premise 1: the bare parse SUCCEEDS for an LRG record");
        assert!(rest.is_empty(), "premise 1: and consumes the whole string");
        assert!(
            crate::hgvs::parser::accession::is_valid_compound_outer(&outer),
            "premise 2: the pairing rule ADMITS a bare LRG, because \
             `inferred_variant_type` reads it as genomic — so it cannot be what declines"
        );
        assert!(
            genomic_wrapper_for(&lrg, "LRG_5").is_none(),
            "so only `outer == *bare` stops `LRG_5(LRG_5):c.…`, which re-parses and is \
             therefore invisible to FERRO_ASSERT_REPARSE"
        );
    }

    /// The growth ladder must reach the cap and then stop, from any starting
    /// window. #1691's fix depends on both halves: never terminating would spin
    /// the fetch loop, and stopping early would leave a shorter, position-
    /// relative bound in place — which is the defect itself.
    #[test]
    fn shuffle_fetch_window_growth_reaches_the_cap_and_then_stops() {
        for start in [1_u64, 100, 1000, MAX_SHUFFLE_FETCH_WINDOW - 1] {
            let mut window = start;
            let mut steps = 0;
            while let Some(next) = grow_shuffle_fetch_window(window) {
                assert!(next > window, "growth from {window} did not advance");
                assert!(next <= MAX_SHUFFLE_FETCH_WINDOW, "growth overshot the cap");
                window = next;
                steps += 1;
                assert!(steps < 64, "growth from {start} did not terminate");
            }
            assert_eq!(
                window, MAX_SHUFFLE_FETCH_WINDOW,
                "growth from {start} stopped short of the cap",
            );
        }
        assert_eq!(grow_shuffle_fetch_window(MAX_SHUFFLE_FETCH_WINDOW), None);
    }

    /// A configured `window_size` of 0 must still advance. Doubling alone leaves
    /// it at 0 forever, and the fetch loop would never terminate — the one input
    /// that turns a bounded retry into a hang.
    #[test]
    fn shuffle_fetch_window_growth_escapes_a_zero_window() {
        assert_eq!(grow_shuffle_fetch_window(0), Some(2));
    }

    /// `Normalizer` has exactly two public normalizing methods — `normalize` and
    /// `normalize_with_diagnostics` — and only the first went through
    /// `normalize_core_checked`. So this one returned a normalized variant none of
    /// the three seam oracles had inspected: `FERRO_ASSERT_IN_BOUNDS`,
    /// `FERRO_ASSERT_REPARSE` and `FERRO_ASSERT_IDEMPOTENT` alike. #1366 covered it
    /// by calling `assert_seam_oracles` from this exit directly; #1382 then routed
    /// the exit through `normalize_core_checked`, so the coverage now comes from the
    /// same single call site as every other public normalization — and the strict
    /// ladder comes with it.
    ///
    /// The flags themselves are process-wide `OnceLock`s read from the
    /// environment, so a unit test cannot toggle one; asserting the *invariant*
    /// with the same predicate the oracle uses is the env-independent equivalent,
    /// and it fails if this path ever starts emitting an out-of-range coordinate
    /// again. In CI, where all three flags are set, this path also gets the
    /// live assertions.
    ///
    /// The inputs are the shapes the in-bounds class was found in (#1307's
    /// terminal-dup collision and #1274's two-base identity), each authored at the
    /// end of a 24-base contig where an overrun is one step away.
    #[test]
    fn the_diagnostics_exit_returns_in_bounds_coordinates_too() {
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_TEST.1", "TTTTTTTTTAATATATTTTAATAC".to_string());
        let normalizer = Normalizer::new(provider);
        for descriptor in [
            "NC_TEST.1:g.[24dup;24C>G]",
            "NC_TEST.1:g.[24dup;24delinsGG]",
            "NC_TEST.1:g.24dup",
            "NC_TEST.1:g.23_24insC",
            "NC_TEST.1:g.24C>G",
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            let diagnosed = normalizer
                .normalize_with_diagnostics(&variant)
                .unwrap_or_else(|e| panic!("`{descriptor}` must normalize: {e}"));
            let overrun =
                merge::first_out_of_bounds_coordinate(&diagnosed.result, normalizer.provider());
            assert!(
                overrun.is_none(),
                "`{descriptor}` -> `{}` through normalize_with_diagnostics names a \
                 coordinate off the end: {overrun:?}",
                diagnosed.result
            );
        }
    }

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

    /// Build a single-exon coding transcript whose poly-A run starts inside
    /// the CDS (`c.5`) and continues into the 3'UTR, ending at `c.*4`:
    ///
    /// ```text
    ///   tx pos:  1  2  3  4  5  6  7  8 | 9 10 11 12 |13 14 15
    ///   base:    A  T  G  C  A  A  A  A | A  A  A  A | G  C  T
    ///   c.:      1  2  3  4  5  6  7  8 |*1 *2 *3 *4 |*5 *6 *7
    ///                        \___ CDS ___/  \_ 3'UTR run _/
    /// ```
    ///
    /// `cds_end = 8`, so the A-run spans the CDS↔3'UTR boundary. A single-base
    /// deletion/duplication anywhere in the run has its most-3' equivalent at
    /// `c.*4` (tx 12), which is in the 3'UTR. Single exon so the exon/exon
    /// exception never applies.
    fn provider_with_boundary_homopolymer() -> MockProvider {
        use crate::reference::transcript::{Exon, ManeStatus, Strand};
        let mut provider = MockProvider::new();
        provider.add_transcript(crate::reference::transcript::Transcript {
            cds_start_incomplete: false,
            id: "NM_TEST918.1".to_string(),
            gene_symbol: Some("BND918".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGCAAAAAAAAGCT".to_string()),
            cds_start: Some(1),
            cds_end: Some(8),
            exons: vec![Exon::new(1, 1, 15)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: std::sync::OnceLock::new(),
        });
        provider
    }

    #[test]
    fn cds_deletion_3prime_shifts_across_cds_utr_boundary() {
        // #918: the HGVS 3' rule (general.md; deletion.md L20) mandates the
        // most-3' position "of the reference sequence" even when it lies in the
        // 3'UTR (positive c.N -> c.*M). A CDS-resident del in a poly-A run that
        // continues past cds_end must shift to c.*4del, not clamp at cds_end.
        let normalizer = Normalizer::new(provider_with_boundary_homopolymer());
        let variant = parse_hgvs("NM_TEST918.1:c.5del").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(
            format!("{}", normalized),
            "NM_TEST918.1:c.*4del",
            "a CDS del whose most-3' equivalent is in the 3'UTR must re-express as c.*Ndel"
        );
    }

    #[test]
    fn cds_duplication_3prime_shifts_across_cds_utr_boundary() {
        // Mirror of the deletion case for a duplication (duplication.md L24).
        let normalizer = Normalizer::new(provider_with_boundary_homopolymer());
        let variant = parse_hgvs("NM_TEST918.1:c.5dup").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(
            format!("{}", normalized),
            "NM_TEST918.1:c.*4dup",
            "a CDS dup whose most-3' equivalent is in the 3'UTR must re-express as c.*Ndup"
        );
    }

    #[test]
    fn cds_utr_boundary_shift_is_idempotent() {
        // Normalizing the already-shifted 3'UTR form is a fixed point.
        let normalizer = Normalizer::new(provider_with_boundary_homopolymer());
        let once = normalizer
            .normalize(&parse_hgvs("NM_TEST918.1:c.5del").unwrap())
            .unwrap();
        let twice = normalizer.normalize(&once).unwrap();
        assert_eq!(
            format!("{}", once),
            format!("{}", twice),
            "re-normalizing the shifted c.*4del must be a fixed point"
        );
        assert_eq!(format!("{}", twice), "NM_TEST918.1:c.*4del");
    }

    /// Transcript for the whole-CDS-spanning deletion case (#918 cross-axis).
    ///
    /// ```text
    ///   tx:   1  2  3  4  5  6  7  8  9 10
    ///   seq:  C  A  T  G  A  A  A  G  C  T
    ///   c.:  -1  1  2  3  4  5  6 *1 *2 *3
    ///           \______ CDS ______/ \_3'UTR_/
    /// ```
    ///
    /// `cds_start = 2`, `cds_end = 7`. A deletion spanning the entire CDS
    /// (`c.-1_*1del`, 5'UTR→3'UTR) 3'-shifts exactly one base to `c.1_*2del`
    /// because `seq[c.-1] == seq[c.*2] == C`, then stops (`seq[c.1]=A !=
    /// seq[c.*3]=T`). Single exon, so the exon/exon exception never applies.
    fn provider_with_whole_cds_span() -> MockProvider {
        use crate::reference::transcript::{Exon, ManeStatus, Strand};
        let mut provider = MockProvider::new();
        provider.add_transcript(crate::reference::transcript::Transcript {
            cds_start_incomplete: false,
            id: "NM_TESTCDS.1".to_string(),
            gene_symbol: Some("SPAN918".to_string()),
            strand: Strand::Plus,
            sequence: Some("CATGAAAGCT".to_string()),
            cds_start: Some(2),
            cds_end: Some(7),
            exons: vec![Exon::new(1, 1, 10)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: std::sync::OnceLock::new(),
        });
        provider
    }

    #[test]
    fn whole_cds_span_deletion_3prime_shifts_across_both_axes() {
        // #918: a deletion spanning the ENTIRE CDS (one endpoint in the 5'UTR,
        // the other in the 3'UTR) is a well-defined 3'-rule shuffle on the
        // contiguous spliced transcript, even though its endpoints sit in
        // different axes. `c.-1_*1del` must 3'-shift to `c.1_*2del` (matching
        // mutalyzer), not echo the input via the #350 cross-axis bail.
        let normalizer = Normalizer::new(provider_with_whole_cds_span());
        let variant = parse_hgvs("NM_TESTCDS.1:c.-1_*1del").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(
            format!("{}", normalized),
            "NM_TESTCDS.1:c.1_*2del",
            "a whole-CDS-spanning deletion must 3'-shift across the CDS/UTR axes"
        );
    }

    #[test]
    fn whole_cds_span_deletion_shift_is_idempotent() {
        // Re-normalizing the already-shifted `c.1_*2del` is a fixed point.
        let normalizer = Normalizer::new(provider_with_whole_cds_span());
        let once = normalizer
            .normalize(&parse_hgvs("NM_TESTCDS.1:c.-1_*1del").unwrap())
            .unwrap();
        let twice = normalizer.normalize(&once).unwrap();
        assert_eq!(format!("{}", once), format!("{}", twice));
        assert_eq!(format!("{}", twice), "NM_TESTCDS.1:c.1_*2del");
    }

    #[test]
    fn whole_cds_span_deletion_that_cannot_shift_is_unchanged() {
        // A whole-CDS-spanning deletion already at its most-3' position must be
        // left unchanged — the carve-out only ever ADDS the spec-mandated shift,
        // never invents one. `c.-1_*2del` cannot roll right (`seq[c.-1]=C !=
        // seq[c.*3]=T`), so it stays put.
        let normalizer = Normalizer::new(provider_with_whole_cds_span());
        let variant = parse_hgvs("NM_TESTCDS.1:c.-1_*2del").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        assert_eq!(format!("{}", normalized), "NM_TESTCDS.1:c.-1_*2del");
    }

    #[test]
    fn normalize_rejects_versioned_request_on_silent_substitution() {
        // #785: a c./n./r. description is defined only against its named
        // accession.version (background/numbering.md, refseq.md). When the
        // request names an EXPLICIT version the reference lacks, a lenient
        // provider silently falls back to a sibling version and normalizes
        // against *its* frame/sequence — a confidently-wrong result attributed
        // to the requested version. Decline with TranscriptVersionNotExact
        // instead — the clean "reference unavailable" miss the conformance
        // mapping already understands.
        let mut provider = MockProvider::with_test_data();
        // The reference carries only NM_000088.3; a request for .2 is silently
        // served the .3 bases/frame.
        provider.mark_version_substitution("NM_000088.2", "NM_000088.3");
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.2:c.10del").unwrap();
        let err = normalizer.normalize(&variant).expect_err(
            "an explicitly-versioned request the provider would silently \
             substitute must be rejected, not normalized against the sibling",
        );
        match err {
            FerroError::TranscriptVersionNotExact { requested } => {
                assert_eq!(requested, "NM_000088.2");
            }
            other => panic!("expected TranscriptVersionNotExact, got {other:?}"),
        }
    }

    #[test]
    fn normalize_allows_bare_request_even_if_resolved_to_sibling() {
        // #785: the gate applies ONLY to explicitly-versioned requests. A bare
        // (unversioned) accession intentionally keeps lenient "latest version"
        // resolution, so it must NOT be gated even though it resolves to a
        // specific version under the hood.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088:c.10A>G").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "a bare (unversioned) request must not be gated, got {result:?}"
        );
    }

    #[test]
    fn normalize_allows_versioned_request_served_at_exact_version() {
        // #785: a versioned request the reference serves at the EXACT version
        // (no substitution) must pass the gate unchanged (no false positive).
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        assert!(
            normalizer.normalize(&variant).is_ok(),
            "an exact-version request must not be gated"
        );
    }

    #[test]
    fn normalize_allows_versioned_request_absent_without_sibling() {
        // #785: the gate fires only on a silent SUBSTITUTION (a sibling is
        // actually served), never on a genuinely-absent transcript with no
        // sibling to fall back to. Even when the provider reports the version as
        // not exact, an absent transcript (`get_transcript` errors) falls through
        // to the existing missing-transcript path (echoed unchanged), not a
        // TranscriptVersionNotExact rejection — this guards against over-gating
        // absent Ensembl/RefSeq references.
        let mut provider = MockProvider::with_test_data();
        // Not version-exact AND not registered → no sibling is served.
        provider.mark_non_version_exact("NM_ABSENT.2");
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_ABSENT.2:c.100del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "an absent versioned transcript with no sibling must not be gated, got {result:?}"
        );
        assert_eq!(
            format!("{}", variant),
            format!("{}", result.unwrap()),
            "an absent transcript should be echoed unchanged"
        );
    }

    #[test]
    fn normalize_rejects_versioned_noncoding_request_on_silent_substitution() {
        // #785 parity: the gate covers the non-coding (`n.`) axis, not just `c.`.
        let mut provider = MockProvider::with_test_data();
        provider.mark_version_substitution("NM_000088.2", "NM_000088.3");
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.2:n.10del").unwrap();
        let err = normalizer
            .normalize(&variant)
            .expect_err("a versioned n. silent substitution must be rejected");
        assert!(
            matches!(err, FerroError::TranscriptVersionNotExact { ref requested } if requested == "NM_000088.2"),
            "expected TranscriptVersionNotExact for NM_000088.2, got {err:?}"
        );
    }

    #[test]
    fn normalize_rejects_versioned_rna_request_on_silent_substitution() {
        // #785 parity: the gate covers the RNA (`r.`) axis too.
        let mut provider = MockProvider::with_test_data();
        provider.mark_version_substitution("NM_000088.2", "NM_000088.3");
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.2:r.10del").unwrap();
        let err = normalizer
            .normalize(&variant)
            .expect_err("a versioned r. silent substitution must be rejected");
        assert!(
            matches!(err, FerroError::TranscriptVersionNotExact { ref requested } if requested == "NM_000088.2"),
            "expected TranscriptVersionNotExact for NM_000088.2, got {err:?}"
        );
    }

    #[test]
    fn normalize_rejects_versioned_parented_request_on_silent_substitution() {
        // #785: a genomic-context-wrapped request (`NC_…(NM_….v):c.…`) is gated on
        // its inner transcript — `transcript_accession()` strips the wrapper, so
        // the substituted inner NM is still caught.
        let mut provider = MockProvider::with_test_data();
        provider.mark_version_substitution("NM_000088.2", "NM_000088.3");
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NC_000017.11(NM_000088.2):c.10del").unwrap();
        let err = normalizer
            .normalize(&variant)
            .expect_err("a versioned parented silent substitution must be rejected");
        assert!(
            matches!(err, FerroError::TranscriptVersionNotExact { ref requested } if requested == "NM_000088.2"),
            "expected TranscriptVersionNotExact for inner NM_000088.2, got {err:?}"
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

    // Real #882 locus GACACACATAGGT padded with 120-base neutral flanks so the
    // +/-100 normalize window stays in-contig. Motif G is at 1-based g.121, so
    // the A|C cut of the ACACACA tract is between g.122 (A) and g.123 (C).
    fn phase_gate_882_provider() -> MockProvider {
        let flank = "ACGT".repeat(30); // 120 bases
        let seq = format!("{flank}GACACACATAGGT{flank}"); // len 253
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_000001.11", &seq);
        provider
    }

    #[test]
    fn normalize_out_of_phase_ins_stays_ins_882() {
        // g.122_123insAC sits at the A|C cut, out of phase with the ACACACA
        // tract, so it must NOT become a dup (that was the silent-corruption bug).
        let normalizer = Normalizer::new(phase_gate_882_provider());
        let variant = parse_hgvs("NC_000001.11:g.122_123insAC").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", normalized);
        assert!(
            output.contains("ins") && !output.contains("dup"),
            "out-of-phase insAC must stay an insertion, got: {output}"
        );
        // Re-normalizing the result must not flip the phase gate on a second pass.
        let again = format!("{}", normalizer.normalize(&normalized).unwrap());
        assert_eq!(output, again, "normalize must be idempotent, got: {again}");
    }

    #[test]
    fn normalize_in_phase_ins_becomes_dup_882() {
        // g.122_123insCA at the same cut IS a legitimate tandem dup.
        let normalizer = Normalizer::new(phase_gate_882_provider());
        let variant = parse_hgvs("NC_000001.11:g.122_123insCA").unwrap();
        let normalized = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", normalized);
        assert!(
            output.contains("dup"),
            "in-phase insCA should become a dup, got: {output}"
        );
        // The emitted dup must stay a dup when normalized again.
        let again = format!("{}", normalizer.normalize(&normalized).unwrap());
        assert_eq!(output, again, "normalize must be idempotent, got: {again}");
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

        // Position 1 = M (Met), Position 2 = V (Val). Position 2 is used
        // because a substitution at the initiator Met is spec-forbidden
        // (protein/substitution.md:49).
        let variant = parse_hgvs("NP_TEST.1:p.Val2Leu").unwrap();
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
        // NM_001234.1 has cds_start = 5 and a G repeat at tx positions 13-37
        // spanning three exons: exon 1 (1-15), exon 2 (16-30), exon 3 (31-44).
        // `r.` is CDS-relative (== `c.`, #469), so `r.10` = tx 14 — a G in
        // exon 1. Per HGVS general.md's exon-junction exception, an r.
        // deletion does not shift across exon boundaries, so `r.10del` shifts
        // only as far as `r.11del` (= tx 15, the last G inside exon 1).
        // Closes #334.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:r.10del").unwrap();
        let result = normalizer.normalize(&variant).unwrap();
        let output = format!("{}", result);

        assert!(
            output.contains("del"),
            "Deletion should remain a deletion, got: {}",
            output
        );
        // Spec-canonical 3' anchor inside exon 1 is tx 15 = r.11 (CDS-relative).
        assert!(
            output.contains("r.11del"),
            "Deletion should shift 3' to exon-1 boundary (r.11del), got: {}",
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
    fn test_normalize_rna_intronic_warns_and_degrades_without_genome() {
        // #1012: an intronic RNA variant used to error when the provider had
        // no genomic data (`with_test_data` registers no genomic sequences).
        // It now warn-and-degrades — returning the best-effort (unchanged)
        // variant plus a `ReducedCapabilityNoGenome` warning — so a degraded
        // result is distinguishable from a genuinely-final one.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:r.10+5del").unwrap();
        let diag = normalizer
            .normalize_with_diagnostics(&variant)
            .expect("intronic RNA without a genome must warn-and-degrade, not error");
        assert!(
            diag.warnings
                .iter()
                .any(|w| w.code() == "REDUCED_CAPABILITY_NO_GENOME"),
            "expected REDUCED_CAPABILITY_NO_GENOME warning, got {:?}",
            diag.warnings,
        );
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
    fn test_normalize_cis_allele_reorders_to_genomic_order() {
        // #1098: cis-allele members must be rendered in genomic (coordinate)
        // order regardless of the order they were written. A descending-order
        // input normalizes to ascending genomic order.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("[NM_000088.3:c.20C>T;NM_000088.3:c.10A>G]").unwrap();
        let result = normalizer.normalize(&variant).unwrap();

        assert_eq!(
            format!("{}", result),
            "NM_000088.3:c.[10A>G;20C>T]",
            "Descending-order cis members must be reordered into ascending genomic order"
        );
    }

    #[test]
    fn test_normalize_cis_allele_order_independent() {
        // #1098: two inputs that denote the same allele in different member
        // orders must normalize to the same canonical string.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let ascending = normalizer
            .normalize(&parse_hgvs("[NM_000088.3:c.10A>G;NM_000088.3:c.20C>T]").unwrap())
            .unwrap();
        let descending = normalizer
            .normalize(&parse_hgvs("[NM_000088.3:c.20C>T;NM_000088.3:c.10A>G]").unwrap())
            .unwrap();

        assert_eq!(
            format!("{}", ascending),
            format!("{}", descending),
            "Ascending and descending inputs must normalize to the same canonical string"
        );
    }

    #[test]
    fn test_normalize_cis_allele_reorders_three_members() {
        // #1098: a real sort, not just a two-member swap — three members given
        // out of order settle into ascending genomic order.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant =
            parse_hgvs("[NM_000088.3:c.30G>A;NM_000088.3:c.10A>G;NM_000088.3:c.20C>T]").unwrap();
        let result = normalizer.normalize(&variant).unwrap();

        assert_eq!(
            format!("{}", result),
            "NM_000088.3:c.[10A>G;20C>T;30G>A]",
            "Three shuffled cis members must be reordered into ascending genomic order"
        );
    }

    #[test]
    fn test_cis_member_order_key_is_total_order_with_content_tiebreak() {
        // #1098: the ordering key must be a *total* order — two members that
        // share a start position get a deterministic content tie-break (the
        // rendered descriptor), so member order never falls back to input
        // order. This is why a total order beats a stable sort keyed on start
        // alone.
        use crate::hgvs::variant::cis_member_order_key;

        let sub = parse_hgvs("NC_000001.11:g.10A>G").unwrap();
        let del = parse_hgvs("NC_000001.11:g.10del").unwrap();

        let key_sub = cis_member_order_key(&sub);
        let key_del = cis_member_order_key(&del);

        // Same accession, same start point, and — since both are one base wide
        // — the same end point too, so every positional field of the key ties
        // and only the descriptor can separate them (#1261 added the end).
        assert_eq!(
            (&key_sub.0, key_sub.1, key_sub.2),
            (&key_del.0, key_del.1, key_del.2),
            "the two members must share the whole positional portion of the key"
        );
        // ...but distinct keys overall, via the rendered-descriptor tie-break.
        assert_ne!(
            key_sub, key_del,
            "members sharing a start must still get distinct keys (total order)"
        );
    }

    #[test]
    fn test_normalize_preserves_uncertain_and_or_group() {
        // An uncertainty-wrapped and/or group `(a^b)` must keep its `(...)`
        // wrapper (and compact form) through normalization — normalize must
        // not drop the `uncertain` flag when it rebuilds the allele.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_004006.2:c.(370A>C^372C>R)").unwrap();
        let result = normalizer.normalize(&variant).unwrap();

        assert!(matches!(result, HgvsVariant::Allele(_)));
        assert_eq!(
            format!("{}", result),
            "NM_004006.2:c.(370A>C^372C>R)",
            "Uncertain and/or group must round-trip with its parens after normalize"
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
            special: None,
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
            special: None,
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
            special: None,
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
    /// Genomic layout (chr1). Genomic coordinates are 1-based (matching real
    /// cdot data); the gene region begins at g.1001 (after 1000 bp of padding):
    ///
    /// Transcript positions (tx):
    /// - Exon 1: tx 1-20 = genomic 1001-1020
    /// - Intron 1: genomic 1021-1030 (10 bp)
    /// - Exon 2: tx 21-40 = genomic 1031-1050
    /// - Intron 2: genomic 1051-1060 (10 bp)
    /// - Exon 3: tx 41-58 = genomic 1061-1078
    ///
    /// CDS: starts at tx 1 (no 5' UTR for simplicity)
    /// CDS positions: c.1 = tx 1, c.20 = tx 20 (last of exon 1), c.21 = tx 21 (first of exon 2)
    ///
    /// Intron 1 (g.1021-1030): "GTAAGCTAAA" (10 bp)
    ///   - c.20+1 = g.1021, c.20+10 = c.21-1 = g.1030
    ///
    /// Intron 2 (g.1051-1060): "AAAGTAAGGG" (10 bp)
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
        // Genomic coordinates are 1-based (matching real cdot data). The 1000-byte
        // padding occupies 1-based positions 1-1000 (0-based string indices 0-999),
        // so the gene region begins at 1-based g.1001 (string index 1000):
        //   exon1 g.1001-1020, intron1 g.1021-1030, exon2 g.1031-1050,
        //   intron2 g.1051-1060, exon3 g.1061-1078, then trailing padding.
        let mut genomic_seq = String::new();

        // Padding before (1-based g.1-1000, 1000 bytes at 0-based indices 0-999)
        for _ in 0..1000 {
            genomic_seq.push('N');
        }

        // Exon 1 (g.1001-1020)
        genomic_seq.push_str("ATGCCCAAAGGGTTTAGGCC");

        // Intron 1 (g.1021-1030) - with splice consensus
        // Note: The intron has AAA at the end (g.1028-1030) to test shifting
        genomic_seq.push_str("GTAAGCTAAA");

        // Exon 2 (g.1031-1050)
        genomic_seq.push_str("AAAGGGTTTAGGCCCAAAAA");

        // Intron 2 (g.1051-1060) - with AAA at start for testing
        genomic_seq.push_str("AAAGTAAGGG");

        // Exon 3 (g.1061-1078)
        genomic_seq.push_str("GGGTTTAGGCCCAAATGA");

        // Padding after (100 bytes)
        for _ in 0..100 {
            genomic_seq.push('N');
        }

        // Add genomic sequence to provider
        provider.add_genomic_sequence("chr1", genomic_seq);

        // Create transcript with exons that have genomic coordinates (1-based)
        provider.add_transcript(Transcript {
            cds_start_incomplete: false,
            id: "NM_BOUNDARY.1".to_string(),
            gene_symbol: Some("BOUNDARY".to_string()),
            strand: Strand::Plus,
            sequence: Some(tx_seq.to_string()),
            cds_start: Some(1),
            cds_end: Some(58),
            exons: vec![
                Exon::with_genomic(1, 1, 20, 1001, 1020),
                Exon::with_genomic(2, 21, 40, 1031, 1050),
                Exon::with_genomic(3, 41, 58, 1061, 1078),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1001),
            genomic_end: Some(1078),
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
    /// end (g.1078) and tx 58 to the low end (g.1001). Genomic coordinates are
    /// 1-based (the 1000-byte padding occupies g.1-1000). Intron 2 is laid
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
        // Exon 3 region (g.1001-1018): RC of tx[41..58] ("GGGTTTAGGCCCAAATGA").
        genomic_seq.push_str("TCATTTGGGCCTAAACCC");
        // Intron 2 (g.1019-1028): the last four bases (g.1025-1028) are 'T',
        // so c.40+1..c.40+4 read as 'A' in transcript view, extending the
        // exonic poly-A across the boundary.
        genomic_seq.push_str("AAAGTATTTT");
        // Exon 2 region (g.1029-1048): RC of tx[21..40] ("AAAGGGTTTAGGCCCAAAAA").
        genomic_seq.push_str("TTTTTGGGCCTAAACCCTTT");
        // Intron 1 (g.1049-1058): mirrors the plus fixture's intron 1 content.
        genomic_seq.push_str("GTAAGCTAAA");
        // Exon 1 region (g.1059-1078): RC of tx[1..20] ("ATGCCCAAAGGGTTTAGGCC").
        genomic_seq.push_str("GGCCTAAACCCTTTGGGCAT");
        for _ in 0..100 {
            genomic_seq.push('N');
        }

        provider.add_genomic_sequence("chr1", genomic_seq);

        provider.add_transcript(Transcript {
            cds_start_incomplete: false,
            id: "NM_BOUNDARYM.1".to_string(),
            gene_symbol: Some("BOUNDARY_M".to_string()),
            strand: Strand::Minus,
            sequence: Some(tx_seq.to_string()),
            cds_start: Some(1),
            cds_end: Some(58),
            exons: vec![
                Exon::with_genomic(1, 1, 20, 1059, 1078),
                Exon::with_genomic(2, 21, 40, 1029, 1048),
                Exon::with_genomic(3, 41, 58, 1001, 1018),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1001),
            genomic_end: Some(1078),
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
    fn test_get_genomic_sequence_1based_addresses_one_based_inclusive() {
        // Regression guard for the intronic/boundary-spanning off-by-one:
        // `get_genomic_sequence_1based` must treat its `[start, end]` as 1-based
        // inclusive, mapping the base at 1-based position `p` to 0-based index
        // `p - 1` of the underlying 0-based-half-open provider. Dropping the
        // `-1` conversion (the original bug) silently reads one base too far in
        // the +forward direction, corrupting the homopolymer the shuffle
        // inspects at an exon/intron junction (e.g. NM_004992.3:c.378-17del).
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("chrX", "ACGTACGTTT");
        let normalizer = Normalizer::new(provider);

        // 1-based position 1 is the first base.
        assert_eq!(
            normalizer
                .get_genomic_sequence_1based("chrX", 1, 1)
                .unwrap(),
            "A"
        );
        // 1-based inclusive [1, 4] spans the first four bases.
        assert_eq!(
            normalizer
                .get_genomic_sequence_1based("chrX", 1, 4)
                .unwrap(),
            "ACGT"
        );
        // An interior window: 1-based inclusive [5, 7].
        assert_eq!(
            normalizer
                .get_genomic_sequence_1based("chrX", 5, 7)
                .unwrap(),
            "ACG"
        );
        // A single interior base: 1-based position 4 is 'T'.
        assert_eq!(
            normalizer
                .get_genomic_sequence_1based("chrX", 4, 4)
                .unwrap(),
            "T"
        );

        // A 1-based start of 0 is invalid and must be rejected rather than
        // silently fetching `[0, end)` from index 0. Near a contig start the
        // window math can drive `seq_start` to 0; callers clamp to 1, and this
        // guard ensures a stray 0 surfaces a clear error instead of shifting
        // every `rel = g - seq_start + 1` coordinate by one.
        assert!(
            matches!(
                normalizer.get_genomic_sequence_1based("chrX", 0, 4),
                Err(FerroError::ConversionError { .. })
            ),
            "1-based start of 0 must be rejected"
        );
    }

    #[test]
    fn test_boundary_spanning_exonic_to_intronic_del() {
        // Test: c.20_20+3del - deletion from last exon base into intron
        // c.20 = last base of exon 1 (C at g.1020)
        // c.20+3 = 3rd intronic base (A at g.1023)
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
        // c.21-3 = 3rd base before exon 2 (A at g.1028)
        // c.23 = 3rd base of exon 2 (A at g.1033)
        // Deletes: AAA (intronic) + AAA (exonic) = 6 bases
        //
        // The intron ends with AAA (g.1028-1030) and exon starts with AAA (g.1031-1033)
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
        // c.21-1 = last intronic base before exon 2 (A at g.1030)
        // c.22 = 2nd base of exon 2 (A at g.1032)
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

    /// #670 exon/EXON-suppression fixture (plus strand, 1-based genomic).
    ///
    /// A poly-A spans the exon1/exon2 junction in the spliced transcript
    /// (c.9,c.10 in exon 1; c.11,c.12 in exon 2 = `AAAA`), but the intron
    /// between them starts with a non-A base, so the genomic-space run is
    /// broken at the junction. A `c.10del` must therefore stay at `c.10`
    /// (the upstream exon's last base): it may not shift into the non-matching
    /// intron, and it may not bridge across the exon/exon junction into exon 2.
    fn make_exon_exon_homopolymer_provider() -> MockProvider {
        use crate::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
        use std::sync::OnceLock;

        let mut provider = MockProvider::new();

        // tx (spliced): exon1 "ATGCACGTAA" (c.10=A, c.9=A) + exon2 "AACCGGTTAC"
        let tx_seq = "ATGCACGTAAAACCGGTTAC"; // 20 bp
        let mut genomic_seq = String::new();
        for _ in 0..1000 {
            genomic_seq.push('N'); // padding g.1-1000
        }
        genomic_seq.push_str("ATGCACGTAA"); // exon 1, g.1001-1010
        genomic_seq.push_str("GTAAGTACAG"); // intron 1, g.1011-1020 (starts non-A)
        genomic_seq.push_str("AACCGGTTAC"); // exon 2, g.1021-1030
        for _ in 0..100 {
            genomic_seq.push('N');
        }
        provider.add_genomic_sequence("chr1", genomic_seq);

        provider.add_transcript(Transcript {
            cds_start_incomplete: false,
            id: "NM_EXEX.1".to_string(),
            gene_symbol: Some("EXEX".to_string()),
            strand: Strand::Plus,
            sequence: Some(tx_seq.to_string()),
            cds_start: Some(1),
            cds_end: Some(20),
            exons: vec![
                Exon::with_genomic(1, 1, 10, 1001, 1010),
                Exon::with_genomic(2, 11, 20, 1021, 1030),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1001),
            genomic_end: Some(1030),
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
    fn test_purely_exonic_dup_shifts_across_exon_intron_boundary() {
        // #670 edit-type coverage: a purely-exonic dup at the exon's 3' edge
        // must also apply the 3' rule across the junction (same poly-A run as
        // the del test: exon 2 ends AAAAA, intron 2 starts AAA).
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs("NM_BOUNDARY.1:c.40dup").unwrap();
        let output = format!("{}", normalizer.normalize(&variant).unwrap());
        assert!(
            output.contains(":c.40+") && (output.contains("dup") || output.contains("A[")),
            "purely-exonic dup must shift across the exon/intron junction into the intron, got: {}",
            output
        );
    }

    #[test]
    fn test_purely_exonic_ins_canonicalizes_and_crosses_boundary() {
        // #670 edit-type coverage: an insertion of `A` into the exon-2 poly-A
        // (c.38_39insA) canonicalizes to a dup; in the junction-spanning A-run
        // the 3'-most dup lands intronic. Exercises the ins→dup path through the
        // genomic re-shuffle.
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs("NM_BOUNDARY.1:c.38_39insA").unwrap();
        let output = format!("{}", normalizer.normalize(&variant).unwrap());
        assert_eq!(
            output, "NM_BOUNDARY.1:c.40+3dup",
            "purely-exonic insertion into the junction-spanning A-run must \
             canonicalize to a dup that shifts into the intron, got: {}",
            output
        );
    }

    #[test]
    fn test_minus_strand_purely_exonic_del_shifts_across_boundary() {
        // #670 minus-strand mirror: c.40del on NM_BOUNDARYM.1. The poly-A run is
        // c.36-40 (exon 2, AAAAA) continuing into the intron (c.40+1.. read as A
        // in transcript view), so the 3'-most deletion lands intronic.
        let provider = make_boundary_test_provider_minus();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs("NM_BOUNDARYM.1:c.40del").unwrap();
        let output = format!("{}", normalizer.normalize(&variant).unwrap());
        assert!(
            output.contains(":c.40+") && output.contains("del"),
            "minus-strand purely-exonic del must shift across the exon/intron \
             junction into the intron, got: {}",
            output
        );
    }

    #[test]
    fn test_purely_exonic_del_5prime_direction_stays_exon_confined() {
        // #670 is 3'-only. Under the internal 5' shuffle direction, a purely-exonic del
        // must NOT cross into the intron — it shifts toward the 5' end within
        // the exon, never emitting an intronic offset.
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::with_config(
            provider,
            NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
        );
        let variant = parse_hgvs("NM_BOUNDARY.1:c.40del").unwrap();
        let output = format!("{}", normalizer.normalize(&variant).unwrap());
        assert!(
            !output.contains('+') && !output.contains('-'),
            "5'-direction exonic del must stay exon-confined (no intronic offset), got: {}",
            output
        );
    }

    #[test]
    fn test_exon_exon_junction_homopolymer_del_not_bridged() {
        // #670 exon/EXON suppression: a poly-A spanning the exon1/exon2 junction
        // (AAAA across c.9..c.12 in the spliced transcript), but the intron
        // between the exons starts with a non-A base. `c.10del` must stay at
        // `c.10` — not shift into the non-matching intron, and not bridge across
        // the exon/exon junction into exon 2 (which would give `c.11del`).
        let provider = make_exon_exon_homopolymer_provider();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs("NM_EXEX.1:c.10del").unwrap();
        let output = format!("{}", normalizer.normalize(&variant).unwrap());
        assert_eq!(
            output, "NM_EXEX.1:c.10del",
            "exon/exon-junction-spanning homopolymer must not bridge into the \
             next exon or shift into a non-matching intron, got: {}",
            output
        );
    }

    /// #670 sub-problem B fixture (plus strand, 1-based genomic): a poly-A run
    /// spans intron 1 + exon 2 + intron 2 — `g.1017-1034` (= c.11-4 .. c.20+4),
    /// 18 A's. A deletion with BOTH endpoints intronic but in DIFFERENT introns
    /// exercises the multi-intron-span path (`normalize_intronic_cds` assumes a
    /// single intron and would not shift it).
    fn make_multi_intron_homopolymer_provider() -> MockProvider {
        use crate::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
        use std::sync::OnceLock;

        let mut provider = MockProvider::new();
        let tx_seq = "ATGCACGTACAAAAAAAAAACGTACGTACG"; // exon1 + 10·A (exon2) + exon3
        let mut g = String::new();
        for _ in 0..1000 {
            g.push('N'); // padding g.1-1000
        }
        g.push_str("ATGCACGTAC"); // exon 1   g.1001-1010
        g.push_str("GTAAGCAAAA"); // intron 1 g.1011-1020 (last 4 = A)
        g.push_str("AAAAAAAAAA"); // exon 2   g.1021-1030 (10 A)
        g.push_str("AAAAGTACGT"); // intron 2 g.1031-1040 (first 4 = A)
        g.push_str("CGTACGTACG"); // exon 3   g.1041-1050
        for _ in 0..100 {
            g.push('N');
        }
        provider.add_genomic_sequence("chr1", g);
        provider.add_transcript(Transcript {
            cds_start_incomplete: false,
            id: "NM_MI.1".to_string(),
            gene_symbol: Some("MI".to_string()),
            strand: Strand::Plus,
            sequence: Some(tx_seq.to_string()),
            cds_start: Some(1),
            cds_end: Some(30),
            exons: vec![
                Exon::with_genomic(1, 1, 10, 1001, 1010),
                Exon::with_genomic(2, 11, 20, 1021, 1030),
                Exon::with_genomic(3, 21, 30, 1041, 1050),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1001),
            genomic_end: Some(1050),
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
    fn test_multi_intron_spanning_del_normalizes_across_span() {
        // #670 sub-problem B: a deletion whose endpoints are intronic but in
        // DIFFERENT introns (c.11-2 in intron 1, c.20+2 in intron 2, spanning
        // all of exon 2). Dispatch must route it to the genomic boundary path
        // (not `normalize_intronic_cds`, which assumes one intron and would
        // leave it unshifted). 14 A's removed from the 18-A run (c.11-4..c.20+4)
        // canonicalizes to repeat notation `A[4]` over the full run extent.
        let normalizer = Normalizer::new(make_multi_intron_homopolymer_provider());
        let variant = parse_hgvs("NM_MI.1:c.11-2_20+2del").unwrap();
        let output = format!("{}", normalizer.normalize(&variant).unwrap());
        assert_eq!(
            output, "NM_MI.1:c.11-4_20+4A[4]",
            "multi-intron-spanning deletion in a homopolymer must normalize \
             across the full span, got: {}",
            output
        );

        // Idempotent.
        let reparsed = parse_hgvs(&output).unwrap();
        let again = format!(
            "{}",
            Normalizer::new(make_multi_intron_homopolymer_provider())
                .normalize(&reparsed)
                .unwrap()
        );
        assert_eq!(
            again, output,
            "multi-intron-span normalization must be idempotent"
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
    fn test_purely_exonic_del_shifts_across_exon_intron_boundary() {
        // #670: a purely-exonic deletion at the last base of an exon must apply
        // the 3' rule ACROSS the exon/intron junction (numbering.md:22-26 NOTE;
        // deletion.md exon/intron border). In the boundary fixture (1-based
        // genomic), exon 2 ends with AAAAA (c.36-40, g.1046-1050) and intron 2
        // starts with AAA (c.40+1..+3, g.1051-1053): an 8-A run straddling the
        // c.40 junction. Deleting one A (c.40del) is 3'-most at the last A of
        // the run, which is intronic: c.40+3del. ferro previously kept it
        // exon-confined (c.40del).
        let provider = make_boundary_test_provider();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_BOUNDARY.1:c.40del").unwrap();
        let output = format!("{}", normalizer.normalize(&variant).unwrap());
        assert_eq!(
            output, "NM_BOUNDARY.1:c.40+3del",
            "exonic deletion in a junction-spanning homopolymer must shift into the intron"
        );
    }

    #[test]
    fn test_purely_exonic_del_stays_exon_confined_without_genomic_context() {
        // #670 guard: with NO genomic context (bare NM_, no chromosome), there is
        // no intron to shift into, so a purely-exonic deletion must remain
        // exon-confined. MockProvider::with_test_data's NM_001234.1 has no
        // genomic coordinates.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        // A plain exonic deletion resolves and normalizes within the transcript;
        // it must not attempt (or error on) a genomic-space boundary shift.
        let variant = parse_hgvs("NM_001234.1:c.5del").unwrap();
        let result = normalizer.normalize(&variant);
        assert!(
            result.is_ok(),
            "bare-NM exonic del should normalize exon-confined, got: {:?}",
            result.err()
        );
        let output = format!("{}", result.unwrap());
        assert!(
            !output.contains('+') && !output.contains('-'),
            "bare-NM exonic del must stay exon-confined (no intronic offset), got: {}",
            output
        );
    }

    #[test]
    fn test_boundary_bare_intronic_no_genomic_data_warns_w4007() {
        // #682: `c.11_11+3del` spans the exon/intron boundary on a bare `NM_`
        // (no NG_(…)/NC_(…) context) — itself a spec-invalid intronic form. With
        // no genomic data the intronic position cannot resolve. Default (warn)
        // mode previously propagated that capability error; now the form-level
        // invalidity (W4007) takes precedence: the input is echoed and the
        // warning surfaced rather than hard-failing. (Strict mode still rejects
        // it as EINTRONIC — see tests/issue_486_eintronic.rs.)
        let provider = MockProvider::with_test_data(); // No genomic data
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_001234.1:c.11_11+3del").unwrap();
        // Warn-only must not propagate the capability error.
        let out = normalizer
            .normalize(&variant)
            .expect("bare intronic form must warn, not error, in default mode (#682)");
        assert_eq!(out.to_string(), "NM_001234.1:c.11_11+3del");
        // W4007 is surfaced via diagnostics rather than silently dropped.
        let diag = normalizer
            .normalize_with_diagnostics(&variant)
            .expect("diagnostics should succeed in warn-only mode");
        assert!(
            diag.warnings
                .iter()
                .any(|w| w.code() == "INTRONIC_ON_BARE_TRANSCRIPT"),
            "expected INTRONIC_ON_BARE_TRANSCRIPT (W4007), got {:?}",
            diag.warnings
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
            cds_start_incomplete: false,
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
            cds_start_incomplete: false,
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
        //
        // Since #1079 the description never reaches the normalizer at all:
        // an insertion anchor MUST name two flanking positions listed 5' to
        // 3' (`DNA/insertion.md:15-16`), so the parser refuses it up front.
        // That is a strictly stronger guarantee than "does not panic", but
        // keep the normalizer half exercised on a legal inverted-range edit
        // so the original slice-indexing regression stays covered.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        assert!(
            parse_hgvs("NC_000011.10:g.5238138_5153222insTATTT").is_err(),
            "non-flanking insertion anchor must be refused at parse time"
        );

        let variant = parse_hgvs("NC_000011.10:g.5238138_5153222dup").unwrap();
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
            cds_start_incomplete: false,
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
            special: None,
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
            special: None,
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
        // A coding (`c.`) variant whose provider transcript record carries no
        // CDS bounds (cds_start/cds_end = None) must normalize without error.
        // Uses an NM_ accession so the c. coordinate system is valid (a c.
        // description on a non-coding NR_/XR_ reference is rejected at parse,
        // #486); the missing-CDS path is exercised via the provider record.
        let mut provider = MockProvider::new();
        use crate::reference::transcript::{Exon, Transcript};
        provider.add_transcript(Transcript {
            cds_start_incomplete: false,
            id: "NM_001566.1".to_string(),
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

        // c. variant whose transcript record carries no CDS bounds
        let variant = parse_hgvs("NM_001566.1:c.10del").unwrap();
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

        // Create transcript: 2 exons with an intron in between (1-based genomic)
        // Exon 1: tx 1-100, genomic 1001-1100
        // Intron: genomic 1101-1200
        // Exon 2: tx 101-200, genomic 1201-1300
        // Sequence in the intron around position 1101+: AAAA... (for shifting test)
        let tx_sequence = "A".repeat(200);

        provider.add_transcript(Transcript {
            cds_start_incomplete: false,
            id: "NR_038982.1".to_string(),
            gene_symbol: Some("NCRNA_TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some(tx_sequence),
            cds_start: None,
            cds_end: None,
            exons: vec![
                Exon::with_genomic(1, 1, 100, 1001, 1100),
                Exon::with_genomic(2, 101, 200, 1201, 1300),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1001),
            genomic_end: Some(1300),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: std::sync::OnceLock::new(),
        });

        // Add genomic sequence for chr1 around 1-based positions 1001-1300
        // Make the intron region (g.1101-1200) be "AGCT" repeated to test shifting
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
    fn intronic_window_bounds_covers_intron_and_caps() {
        // Common case: the intron fits inside the variant window, so the
        // window is unchanged (the fix is a no-op here).
        assert_eq!(
            intronic_window_bounds(1000, 1010, 990, 1020, 100),
            (900, 1110)
        );

        // Large intron: the window extends to cover the far intron edge, plus
        // one (the fetch is end-exclusive, so the edge must be < seq_end to
        // land inside the 1-based window) — issue #573.
        let (s, e) = intronic_window_bounds(1000, 1010, 200, 5000, 100);
        assert_eq!((s, e), (200, 5001));

        // Asymmetric extension: only the start edge is outside the variant
        // window (intron_g_start < var_start, intron_g_end is within reach).
        // This is the shape closest to the #573 bug — variant near the far end
        // of a moderately large intron, where only one edge needs extending.
        let (s2, e2) = intronic_window_bounds(1000, 1010, 850, 1050, 100);
        assert_eq!((s2, e2), (850, 1110));

        // Pathologically large intron beyond the cap: fall back to the
        // variant-sized window (the downstream guard then returns a clean Err
        // rather than forcing a huge fetch).
        assert_eq!(
            intronic_window_bounds(1_000_000, 1_000_010, 1, 5_000_000, 100),
            (999_900, 1_000_110)
        );

        // Near the start of a contig the window start is clamped to the 1-based
        // minimum: `g_start - window` (and the intron-extended `want_start`)
        // would otherwise saturate to 0, which has no valid 1-based coordinate
        // and breaks the caller's `rel = g - seq_start + 1` arithmetic.
        // Variant near the contig start, intron entirely to the right.
        let (s_near, _) = intronic_window_bounds(50, 60, 40, 200, 100);
        assert_eq!(s_near, 1, "window start near contig start must clamp to 1");
        // Cap-fallback path also clamps: a huge intron forces the variant-sized
        // window, whose start (`g_start - window`) still underflows to 0 here.
        let (s_cap, _) = intronic_window_bounds(50, 60, 1, 5_000_000, 100);
        assert_eq!(
            s_cap, 1,
            "cap-fallback window start near contig start must clamp to 1"
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
    fn genome_cen_refusal_preserves_edit_body_verbatim() {
        use crate::normalize::NormalizationWarning;
        use crate::reference::MockProvider;
        // The W4005 cen refusal must echo the input byte-for-byte, including an
        // explicit deleted base. `canonicalize_genome_variant` would rewrite the
        // edit body, so the refusal path must preserve the original variant.
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_TEST.4", "ACGT".repeat(50));
        let n = Normalizer::new(provider);
        let input = "NC_TEST.4:g.cendelA";
        let v = parse_hgvs(input).unwrap();
        let result = n.normalize_with_diagnostics(&v).unwrap();
        assert_eq!(
            format!("{}", result.result),
            input,
            "cen refusal must preserve the explicit deleted base verbatim"
        );
        assert!(
            result
                .warnings
                .iter()
                .any(|w| matches!(w, NormalizationWarning::UnresolvableSpecialPosition { .. })),
            "must still emit W4005, got {:?}",
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

    #[test]
    fn resolve_special_cds_pos_maps_markers() {
        use crate::hgvs::location::CdsPos;
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
        // 5'UTR = 61 (cds_start=62), CDS 62..=124, sequence length 200 -> qter c.*(200-124)=c.*76.
        let tx = Transcript::new(
            "NM_TEST.1".into(),
            Some("T".into()),
            Strand::Plus,
            "ACGT".repeat(50),
            Some(62),
            Some(124),
            vec![Exon::new(1, 1, 200)],
            None,
            None,
            None,
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );
        let pter = resolve_special_cds_pos(&CdsPos::pter(), &tx)
            .unwrap()
            .unwrap();
        assert_eq!(format!("{pter}"), "-61");
        let qter = resolve_special_cds_pos(&CdsPos::qter(), &tx)
            .unwrap()
            .unwrap();
        assert_eq!(format!("{qter}"), "*76");
        assert_eq!(resolve_special_cds_pos(&CdsPos::cen(), &tx).unwrap(), None);
        let plain = CdsPos::new(42);
        assert_eq!(resolve_special_cds_pos(&plain, &tx).unwrap(), Some(plain));
    }

    #[test]
    fn resolve_special_cds_pos_boundary_branches() {
        use crate::hgvs::location::CdsPos;
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
        // No 5'UTR (cds_start=1) -> pter -> c.1; no 3'UTR (cds_end=L=100) -> qter -> c.100.
        let tx = Transcript::new(
            "NM_NOUTR.1".into(),
            Some("T".into()),
            Strand::Plus,
            "ACGT".repeat(25),
            Some(1),
            Some(100),
            vec![Exon::new(1, 1, 100)],
            None,
            None,
            None,
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );
        assert_eq!(
            format!(
                "{}",
                resolve_special_cds_pos(&CdsPos::pter(), &tx)
                    .unwrap()
                    .unwrap()
            ),
            "1"
        );
        assert_eq!(
            format!(
                "{}",
                resolve_special_cds_pos(&CdsPos::qter(), &tx)
                    .unwrap()
                    .unwrap()
            ),
            "100"
        );
    }

    #[test]
    fn cds_pter_qter_resolve_and_normalize() {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
        let mut provider = crate::reference::MockProvider::new();
        provider.add_transcript(Transcript::new(
            "NM_TEST.1".into(),
            Some("T".into()),
            Strand::Plus,
            "ACGT".repeat(50),
            Some(62),
            Some(124),
            vec![Exon::new(1, 1, 200)],
            None,
            None,
            None,
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        ));
        let n = Normalizer::new(provider);
        let out = format!(
            "{}",
            n.normalize(&parse_hgvs("NM_TEST.1:c.pterdel").unwrap())
                .unwrap()
        );
        assert!(
            out.contains("c.-"),
            "pter del resolves into the 5'UTR, got {out}"
        );
        let out2 = format!(
            "{}",
            n.normalize(&parse_hgvs("NM_TEST.1:c.qterdel").unwrap())
                .unwrap()
        );
        assert!(
            out2.contains("c.*"),
            "qter del resolves into the 3'UTR, got {out2}"
        );
    }

    #[test]
    fn cds_cen_strict_rejects_lenient_accepts() {
        use crate::error_handling::ErrorMode;
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
        let mut provider = crate::reference::MockProvider::new();
        provider.add_transcript(Transcript::new(
            "NM_TEST.2".into(),
            Some("T".into()),
            Strand::Plus,
            "ACGT".repeat(50),
            Some(62),
            Some(124),
            vec![Exon::new(1, 1, 200)],
            None,
            None,
            None,
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        ));
        let strict = Normalizer::with_config(
            provider.clone(),
            NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
        );
        assert!(
            matches!(
                strict.normalize(&parse_hgvs("NM_TEST.2:c.cendel").unwrap()),
                Err(FerroError::InvalidCoordinates { .. })
            ),
            "strict rejects c.cen",
        );
        let lenient = Normalizer::new(provider);
        let out = format!(
            "{}",
            lenient
                .normalize(&parse_hgvs("NM_TEST.2:c.cendel").unwrap())
                .unwrap()
        );
        assert!(out.contains("cen"), "lenient preserves c.cen, got {out}");
    }

    #[test]
    fn cds_cen_emits_warning_even_without_transcript() {
        use crate::error_handling::ErrorMode;
        use crate::normalize::NormalizationWarning;
        // The provider has NO record for this accession, so `get_transcript`
        // fails. `cen` is structurally unresolvable regardless of transcript
        // availability, so the W4005 warning must still surface (and strict mode
        // reject) rather than silently canonicalizing the input away.
        let provider = crate::reference::MockProvider::new();
        let lenient = Normalizer::new(provider.clone());
        let diag = lenient
            .normalize_with_diagnostics(&parse_hgvs("NM_ABSENT.1:c.cendel").unwrap())
            .expect("lenient must not error on an absent-transcript c.cen");
        assert!(
            diag.warnings
                .iter()
                .any(|w| matches!(w, NormalizationWarning::UnresolvableSpecialPosition { .. })),
            "absent-transcript c.cen must still emit W4005, got {:?}",
            diag.warnings
        );
        let strict = Normalizer::with_config(
            provider,
            NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
        );
        assert!(
            matches!(
                strict.normalize(&parse_hgvs("NM_ABSENT.1:c.cendel").unwrap()),
                Err(FerroError::InvalidCoordinates { .. })
            ),
            "strict mode must reject an absent-transcript c.cen",
        );
    }

    #[test]
    fn cds_cen_refusal_preserves_edit_body_verbatim() {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
        // The c. cen refusal must echo the input byte-for-byte, including an
        // explicit deleted base — `canonicalize_cds_variant` would rewrite the
        // edit body, so the refusal path must preserve the original variant.
        let mut provider = crate::reference::MockProvider::new();
        provider.add_transcript(Transcript::new(
            "NM_TEST.2".into(),
            Some("T".into()),
            Strand::Plus,
            "ACGT".repeat(50),
            Some(62),
            Some(124),
            vec![Exon::new(1, 1, 200)],
            None,
            None,
            None,
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        ));
        let input = "NM_TEST.2:c.cendelA";
        let v = parse_hgvs(input).unwrap();
        let out = format!("{}", Normalizer::new(provider).normalize(&v).unwrap());
        assert_eq!(
            out, input,
            "c. cen refusal must preserve the explicit deleted base verbatim"
        );
    }

    #[test]
    fn cds_pter_minus_strand_resolves() {
        // `c.` coordinates are transcript-intrinsic, so `pter` maps to transcript
        // position 1 regardless of strand; this guards against any future
        // strand-dependent regression in the resolve path.
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
        let mut provider = crate::reference::MockProvider::new();
        provider.add_transcript(Transcript::new(
            "NM_MINUS.1".into(),
            Some("T".into()),
            Strand::Minus,
            "ACGT".repeat(50),
            Some(62),
            Some(124),
            vec![Exon::new(1, 1, 200)],
            None,
            None,
            None,
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        ));
        let n = Normalizer::new(provider);
        let out = format!(
            "{}",
            n.normalize(&parse_hgvs("NM_MINUS.1:c.pterdel").unwrap())
                .unwrap()
        );
        assert!(
            out.contains("c.-"),
            "minus-strand pter still resolves into the 5'UTR, got {out}"
        );
    }

    #[test]
    fn cds_ng_parent_pter_refused_not_misresolved() {
        use crate::error_handling::ErrorMode;
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
        let mut provider = crate::reference::MockProvider::new();
        provider.add_transcript(Transcript::new(
            "NM_003002.2".into(),
            Some("SDHD".into()),
            Strand::Plus,
            "ACGT".repeat(50),
            Some(62),
            Some(124),
            vec![Exon::new(1, 1, 200)],
            None,
            None,
            None,
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        ));
        let v = parse_hgvs("NG_012337.1(NM_003002.2):c.pterdel").unwrap();

        // Lenient: left verbatim (NOT c.-61), with the W4006 warning.
        let lenient = Normalizer::new(provider.clone());
        let r = lenient.normalize_with_diagnostics(&v).unwrap();
        assert_eq!(
            format!("{}", r.result),
            "NG_012337.1(NM_003002.2):c.pterdel",
            "NG-parent pter must stay verbatim, not c.-61"
        );
        assert!(
            r.warnings.iter().any(|w| matches!(
                w,
                NormalizationWarning::TranscriptFlankNotDescribable { .. }
            )),
            "must emit W4006, got {:?}",
            r.warnings
        );

        // Strict: rejects.
        let strict = Normalizer::with_config(
            provider,
            NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
        );
        assert!(
            matches!(
                strict.normalize(&v),
                Err(FerroError::InvalidCoordinates { .. })
            ),
            "strict rejects NG-parent flank pter"
        );
    }

    #[test]
    fn cds_ng_parent_pter_refusal_preserves_edit_body_verbatim() {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
        // The W4006 refusal must echo the ORIGINAL c. variant byte-for-byte,
        // including an explicit deleted base. `canonicalize_cds_variant` would
        // rewrite the edit body (dropping the `A` from `delA`), so the refusal
        // path must clone the input, not canonicalize it.
        let mut provider = crate::reference::MockProvider::new();
        provider.add_transcript(Transcript::new(
            "NM_003002.2".into(),
            Some("SDHD".into()),
            Strand::Plus,
            "ACGT".repeat(50),
            Some(62),
            Some(124),
            vec![Exon::new(1, 1, 200)],
            None,
            None,
            None,
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        ));
        let input = "NG_012337.1(NM_003002.2):c.pterdelA";
        let v = parse_hgvs(input).unwrap();
        let lenient = Normalizer::new(provider);
        let r = lenient.normalize_with_diagnostics(&v).unwrap();
        assert_eq!(
            format!("{}", r.result),
            input,
            "refusal must preserve the explicit deleted base verbatim"
        );
        assert!(
            r.warnings.iter().any(|w| matches!(
                w,
                NormalizationWarning::TranscriptFlankNotDescribable { .. }
            )),
            "must still emit W4006, got {:?}",
            r.warnings
        );
    }

    #[test]
    fn cds_bare_transcript_pter_still_resolves() {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
        let mut provider = crate::reference::MockProvider::new();
        provider.add_transcript(Transcript::new(
            "NM_003002.2".into(),
            Some("SDHD".into()),
            Strand::Plus,
            "ACGT".repeat(50),
            Some(62),
            Some(124),
            vec![Exon::new(1, 1, 200)],
            None,
            None,
            None,
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        ));
        let n = Normalizer::new(provider);
        let out = format!(
            "{}",
            n.normalize(&parse_hgvs("NM_003002.2:c.pterdel").unwrap())
                .unwrap()
        );
        assert_eq!(
            out, "NM_003002.2:c.-61del",
            "bare NM_ pter still resolves to the transcript boundary"
        );
    }
}

#[cfg(test)]
mod eintronic_helper_tests {
    use super::*;
    use crate::parse_hgvs;

    #[test]
    fn bare_nm_intronic_sub_yields_warning() {
        let v = parse_hgvs("NM_004006.2:c.357+1G>A").unwrap();
        let w = intronic_on_bare_transcript_warning(&v).expect("bare NM intronic must warn");
        assert_eq!(w.code(), "INTRONIC_ON_BARE_TRANSCRIPT");
    }

    #[test]
    fn bare_nr_intronic_del_yields_warning() {
        let v = parse_hgvs("NR_038420.1:n.100+10del").unwrap();
        assert!(intronic_on_bare_transcript_warning(&v).is_some());
    }

    #[test]
    fn ng_nm_intronic_sub_no_warning() {
        let v = parse_hgvs("NG_012337.1(NM_004006.2):c.357+1G>A").unwrap();
        assert!(intronic_on_bare_transcript_warning(&v).is_none());
    }

    #[test]
    fn bare_nm_exonic_no_warning() {
        let v = parse_hgvs("NM_004006.2:c.357G>A").unwrap();
        assert!(intronic_on_bare_transcript_warning(&v).is_none());
    }

    #[test]
    fn bare_nm_minus_offset_intronic_yields_warning() {
        // The check is offset-sign-agnostic: a `-` (3'-of-intron) offset on a
        // bare transcript is just as spec-invalid as a `+` offset.
        let v = parse_hgvs("NM_004006.2:c.357-1G>A").unwrap();
        assert!(intronic_on_bare_transcript_warning(&v).is_some());
    }

    #[test]
    fn bare_xm_intronic_yields_warning() {
        // Predicted coding mRNA (XM_) has the same "no introns on a coding
        // reference" property as NM_; the bare intronic form is spec-invalid.
        let v = parse_hgvs("XM_011535484.2:c.357+1G>A").unwrap();
        assert!(intronic_on_bare_transcript_warning(&v).is_some());
    }

    #[test]
    fn bare_nm_uncertain_breakpoint_range_yields_warning() {
        // Uncertain-breakpoint exon deletion: both endpoints are `Range`
        // boundaries whose inner positions are intronic. `.inner()` skips
        // `Range`, so the boundary-aware check must inspect the endpoints.
        let v = parse_hgvs("NM_004006.2:c.(4185+1_4186-1)_(4357+1_4358-1)del").unwrap();
        assert!(intronic_on_bare_transcript_warning(&v).is_some());
    }

    #[test]
    fn ng_nm_uncertain_breakpoint_range_no_warning() {
        // Same uncertain-breakpoint form on the genomic-context reference is
        // the spec-valid description — no warning.
        let v =
            parse_hgvs("NG_012337.1(NM_004006.2):c.(4185+1_4186-1)_(4357+1_4358-1)del").unwrap();
        assert!(intronic_on_bare_transcript_warning(&v).is_none());
    }

    #[test]
    fn bare_lrg_coding_intronic_yields_warning() {
        // #834: an LRG transcript (`LRG_<N>t<k>`) is a bare coding reference
        // with no NG_/NC_ genomic context — same "no introns on a coding
        // reference" property as NM_, so a bare intronic c. form is spec-invalid
        // and must warn, exactly like its NM-addressed equivalent.
        let v = parse_hgvs("LRG_741t1:c.7007+30del").unwrap();
        let w = intronic_on_bare_transcript_warning(&v).expect("bare LRG intronic must warn");
        assert_eq!(w.code(), "INTRONIC_ON_BARE_TRANSCRIPT");
    }

    #[test]
    fn bare_lrg_coding_minus_offset_intronic_yields_warning() {
        // Offset-sign-agnostic, matching the NM behavior.
        let v = parse_hgvs("LRG_763t1:c.53-6del").unwrap();
        assert!(intronic_on_bare_transcript_warning(&v).is_some());
    }

    #[test]
    fn bare_lrg_coding_exonic_no_warning() {
        // A purely exonic LRG c. variant is a valid form — no warning.
        let v = parse_hgvs("LRG_741t1:c.7007del").unwrap();
        assert!(intronic_on_bare_transcript_warning(&v).is_none());
    }

    #[test]
    fn bare_lrg_noncoding_intronic_yields_warning() {
        // #834: the `HV::Tx` (`n.`) arm was widened to treat an LRG transcript
        // as a bare non-coding reference, exactly like `NR_`/`XR_`. An LRG
        // `n.<pos>+<offset>` intronic form has no introns on its bare reference,
        // so it must warn just like the `NR_` form above.
        let v = parse_hgvs("LRG_741t1:n.100+10del").unwrap();
        let w = intronic_on_bare_transcript_warning(&v).expect("bare LRG n. intronic must warn");
        assert_eq!(w.code(), "INTRONIC_ON_BARE_TRANSCRIPT");
    }

    #[test]
    fn bare_lrg_noncoding_minus_offset_intronic_yields_warning() {
        // Offset-sign-agnostic on the `n.` axis too, matching the `NR_` behavior.
        let v = parse_hgvs("LRG_741t1:n.100-6del").unwrap();
        assert!(intronic_on_bare_transcript_warning(&v).is_some());
    }

    #[test]
    fn bare_lrg_noncoding_exonic_no_warning() {
        // A purely exonic LRG n. variant is a valid form — no warning.
        let v = parse_hgvs("LRG_741t1:n.100del").unwrap();
        assert!(intronic_on_bare_transcript_warning(&v).is_none());
    }
}
