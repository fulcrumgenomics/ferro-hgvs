//! HGVS to SPDI conversion.
//!
//! This module provides conversion functions between HGVS and SPDI formats.
//!
//! # Coordinate Systems
//!
//! - HGVS uses 1-based coordinates (see [`OneBasedPos`])
//! - SPDI uses 0-based interbase coordinates (see [`ZeroBasedPos`])
//!
//! For a variant at genomic position 12345 (1-based):
//! - HGVS: `NC_000001.11:g.12345A>G`
//! - SPDI: `NC_000001.11:12344:A:G` (0-based)
//!
//! # Supported Conversions
//!
//! | HGVS | SPDI | Provider |
//! |------|------|----------|
//! | Substitution `g.12345A>G` | `seq:12344:A:G` | not required |
//! | Deletion `g.100_102del` | `seq:99:ATG:` | required (short form) |
//! | Deletion `g.100_102delATG` | `seq:99:ATG:` | not required (explicit form) |
//! | Insertion `g.100_101insATG` | `seq:100::ATG` | not required |
//! | Delins `g.100_102delinsATG` | `seq:99:ATG:ATG` | required |
//! | Duplication `g.100_102dup` | `seq:102::ATG` | required (short form) |
//! | Duplication `g.100_102dupATG` | `seq:102::ATG` | not required (explicit form) |
//! | Inversion `g.100_102inv` | `seq:99:ATG:CAT` | required (short form) |
//! | Inversion `g.100_102invATG` | `seq:99:ATG:CAT` | not required (explicit form) |
//! | Repeat `g.100_105AT[5]` | `seq:99:ATATAT:ATATATATAT` | required |
//!
//! Short-form deletions, delins, duplications, inversions, and repeats use
//! [`hgvs_to_spdi`] with a [`ReferenceProvider`] to fetch reference bases for
//! SPDI's mandatory `del` (or, for duplication, `ins`) field. Explicit forms
//! do not consult the provider.
//!
//! SPDI has no native `inv` or repeat shape; both convert to `delins`. The
//! inverse direction ([`spdi_to_hgvs`]) is therefore lossy for these edits —
//! an SPDI built from `g.100_102inv` round-trips back to
//! `g.100_102delinsCAT`, not to `inv`. Detecting reverse-complement and
//! repeat structure on input SPDI is tracked in #81 (items A2, B1).
//!
//! # Coordinate-system support
//!
//! [`hgvs_to_spdi_simple`] accepts coordinate systems whose positions are
//! resolvable without provider data:
//!
//! - `g.` (genomic) — direct, for a position that names a coordinate. A
//!   `+`/`-` offset (#1628) or a `pter`/`qter`/`cen` (#1643) is refused
//!   wherever it is written, a complex `(a_b)` boundary included, and an end
//!   boundary that names no single coordinate — `(a_b)` or `?` — is refused
//!   rather than collapsed onto the start (#1795).
//! - `m.` (mito) — the mito accession is genomic; same path as `g.`.
//!   Wraparound ranges (`start > end`) on circular references are rejected
//!   per issue #399 — SPDI has no native wraparound representation.
//! - `o.` (circular) — same path as `g.`/`m.`; wraparound ranges rejected
//!   as above.
//! - `n.` (non-coding tx) — exonic, positive base; SPDI on the transcript
//!   accession.
//! - `r.` (RNA) — exonic, positive base; `u`/`U` rewritten to `T` for
//!   SPDI's DNA alphabet convention.
//!
//! [`hgvs_to_spdi`] additionally handles `c.` (CDS) and UTR-style
//! `n.`/`r.` positions by consulting a [`ReferenceProvider`] for transcript
//! metadata, and uses the same provider to fetch reference bases for
//! short-form `Deletion` / `Duplication` / `Delins` edits across all
//! coordinate systems. Intronic `n.`/`r.` positions remain unsupported
//! (SPDI has no offset notation; genomic projection is future work). Per
//! the [SPDI spec], the SPDI accession matches the HGVS accession (NCBI
//! Variation Services emits SPDI on transcript accessions the same way).
//!
//! `p.` (protein) variants are not representable in SPDI and are rejected.
//!
//! [`ReferenceProvider`]: crate::reference::provider::ReferenceProvider
//! [SPDI spec]: https://www.ncbi.nlm.nih.gov/variation/notation/
//!
//! [`OneBasedPos`]: crate::coords::OneBasedPos
//! [`ZeroBasedPos`]: crate::coords::ZeroBasedPos

use super::SpdiVariant;
use crate::convert::CoordinateMapper;
use crate::coords::{OneBasedPos, ZeroBasedPos};
use crate::error::FerroError;
use crate::hgvs::edit::{InsertedSequence, NaEdit, RepeatCount, Sequence};
use crate::hgvs::interval::{Interval, UncertainBoundary};
use crate::hgvs::location::{CdsPos, GenomePos, RnaPos, TxPos};
use crate::hgvs::parser::accession::parse_accession;
use crate::hgvs::variant::{
    Accession, CdsVariant, CircularVariant, GenomeVariant, HgvsVariant, LocEdit, MtVariant,
    RnaVariant, TxVariant,
};
use crate::reference::provider::ReferenceProvider;
use crate::reference::transcript::Transcript;
use crate::sequence::reverse_complement;

/// Maximum number of bases allowed in an SPDI `ins` string emitted from a
/// repeat expansion. The repeat count is user-controlled (HGVS `RepeatCount`
/// is a `u64`), so an unbounded `unit.repeat(count)` can be forced into a
/// huge allocation. 100 KB is well above any biological tandem-repeat tract
/// we'd plausibly emit as a single SPDI.
const MAX_REPEAT_EXPANSION_BASES: usize = 100_000;

/// Widest reference window [`resolve_repeat_tract_span`] will read while looking
/// for the physical extent of a start-only repeat anchor.
///
/// Named separately from [`MAX_REPEAT_EXPANSION_BASES`] because it bounds a
/// different thing: that one caps the `ins` string a repeat *count* can
/// generate, this one caps how much reference is *read* to find the tract the
/// anchor names. Reusing the expansion bound for both read as if one number
/// governed both, and the growth loop compared it against a half-width, so the
/// window it declined at was twice the figure its own error message quoted.
///
/// 200 KB is far above any biological tandem repeat — the largest known
/// pathogenic expansions run to a few tens of kilobases — so a tract still
/// growing at this width is a degenerate reference, not a variant.
const MAX_REPEAT_SEARCH_BASES: u64 = 200_000;

/// Error type for conversion failures.
///
/// Marked `#[non_exhaustive]` so a new conversion-failure case is additive
/// rather than breaking for downstream crates matching on it — the same
/// treatment [`crate::FerroError`] and [`crate::error::ErrorCode`] received in
/// #1033.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ConversionError {
    /// The variant type is not supported for conversion.
    UnsupportedVariantType {
        /// Description of the unsupported type.
        description: String,
    },
    /// Missing reference sequence data needed for conversion.
    MissingReferenceData {
        /// Description of what data is missing.
        description: String,
    },
    /// The edit type is not supported for conversion.
    UnsupportedEditType {
        /// Description of the unsupported edit.
        description: String,
    },
    /// Invalid position or interval.
    InvalidPosition {
        /// Description of the position error.
        description: String,
    },
    /// Invalid accession format.
    InvalidAccession {
        /// Description of the accession error.
        description: String,
    },
    /// A reference provider is required to perform this conversion, but none
    /// was supplied. Distinct from [`MissingReferenceData`], which means a
    /// provider was supplied but does not have data for the requested region.
    ///
    /// [`MissingReferenceData`]: ConversionError::MissingReferenceData
    ProviderRequired {
        /// HGVS coordinate-system letter that triggered the error (`c`, `n`, `r`, ...).
        variant_type: String,
        /// Why a provider is needed for this variant.
        reason: String,
    },
}

impl std::fmt::Display for ConversionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConversionError::UnsupportedVariantType { description } => {
                write!(
                    f,
                    "unsupported variant type for conversion: {}",
                    description
                )
            }
            ConversionError::MissingReferenceData { description } => {
                write!(f, "missing reference data: {}", description)
            }
            ConversionError::UnsupportedEditType { description } => {
                write!(f, "unsupported edit type for conversion: {}", description)
            }
            ConversionError::InvalidPosition { description } => {
                write!(f, "invalid position: {}", description)
            }
            ConversionError::InvalidAccession { description } => {
                write!(f, "invalid accession: {}", description)
            }
            ConversionError::ProviderRequired {
                variant_type,
                reason,
            } => {
                write!(
                    f,
                    "reference provider required to convert {}. variant: {}",
                    variant_type, reason
                )
            }
        }
    }
}

impl std::error::Error for ConversionError {}

impl From<ConversionError> for FerroError {
    fn from(err: ConversionError) -> Self {
        FerroError::ConversionError {
            msg: err.to_string(),
        }
    }
}

/// Helper to convert a Sequence to a String.
fn sequence_to_string(seq: &Sequence) -> String {
    seq.to_string()
}

/// Helper to convert an InsertedSequence to a String (for literal sequences only).
fn inserted_sequence_to_string(seq: &InsertedSequence) -> Option<String> {
    match seq {
        InsertedSequence::Literal(s) => Some(s.to_string()),
        _ => None,
    }
}

/// A same-reference position-range insert, decomposed to `(start, end, invert)`
/// 1-based inclusive coordinates, or `None` when `seq` is not one.
///
/// Recognises the bare `PositionRange`/`PositionRangeInv` variants and the
/// single-part `Complex([PositionRange | PositionRangeInv])` shape the parser
/// produces for range inserts — both denote one contiguous slice of the same
/// accession. Multi-part `Complex` payloads, external references and every
/// non-range variant return `None`, so the caller keeps its existing handling
/// for them.
fn position_range_insert(seq: &InsertedSequence) -> Option<(u64, u64, bool)> {
    use crate::hgvs::edit::InsertedPart;
    match seq {
        InsertedSequence::PositionRange { start, end } => Some((*start, *end, false)),
        InsertedSequence::PositionRangeInv { start, end } => Some((*start, *end, true)),
        InsertedSequence::Complex(parts) => match parts.as_slice() {
            [InsertedPart::PositionRange { start, end }] => Some((*start, *end, false)),
            [InsertedPart::PositionRangeInv { start, end }] => Some((*start, *end, true)),
            _ => None,
        },
        _ => None,
    }
}

/// Read a same-reference span `[start, end]` (1-based inclusive) to its literal
/// bases, reverse-complementing when `invert`, and fold the result through the
/// output alphabet — the one place the fetch-invert-fold logic lives.
///
/// Both the top-level `ins50_57`/`ins50_57inv` shape and a `[…;213_271;…]`
/// bracket part resolve through here, so a coordinate span reads identically
/// however it was spelled. Declines with `MissingReferenceData` when no provider
/// is available to read the span, propagating any fetch failure otherwise.
fn read_reference_span<P>(
    accession: &str,
    start: u64,
    end: u64,
    invert: bool,
    alphabet: AlphabetMode,
    provider: Option<&P>,
) -> Result<String, ConversionError>
where
    P: ReferenceProvider + ?Sized,
{
    let provider = provider.ok_or_else(|| ConversionError::MissingReferenceData {
        description: format!(
            "position-range insertion {start}_{end} names bases in the same reference; \
             a provider is required to read them"
        ),
    })?;
    let bases = fetch_reference_bases(provider, accession, start, end)?;
    let bases = if invert {
        reverse_complement(&bases)
    } else {
        bases
    };
    Ok(apply_alphabet(&bases, alphabet))
}

/// Expand a repeat `unit` `count` times, folded through the output alphabet —
/// the one place the overflow and [`MAX_REPEAT_EXPANSION_BASES`] cap checks
/// live, so a top-level `insA[10]` and a `[…;A[10];…]` bracket part expand
/// identically. `unit` is the literal unit spelled once.
fn expand_repeat_unit(
    unit: &str,
    count: u64,
    alphabet: AlphabetMode,
) -> Result<String, ConversionError> {
    let unit = apply_alphabet(unit, alphabet);
    let count = count as usize;
    let expansion_bases =
        unit.len()
            .checked_mul(count)
            .ok_or_else(|| ConversionError::UnsupportedEditType {
                description: format!(
                    "repeat expansion {} x {} overflows usize",
                    unit.len(),
                    count
                ),
            })?;
    if expansion_bases > MAX_REPEAT_EXPANSION_BASES {
        return Err(ConversionError::UnsupportedEditType {
            description: format!(
                "repeat expansion {} bases exceeds SPDI ins-string cap of {} bases",
                expansion_bases, MAX_REPEAT_EXPANSION_BASES
            ),
        });
    }
    Ok(unit.repeat(count))
}

/// Resolve a same-reference position-range insert (`ins50_57`, `ins50_57inv`)
/// to its literal bases by reading them from the reference, mirroring how the
/// omitted bases of `del`/`dup`/`delins` are fetched.
///
/// Returns `Ok(None)` when `seq` is not a position-range insert, so the caller
/// can fall back to literal handling. The read itself is [`read_reference_span`],
/// shared with the bracket-part path so both spell a span identically.
fn resolve_position_range_insert<P>(
    seq: &InsertedSequence,
    accession: &str,
    alphabet: AlphabetMode,
    provider: Option<&P>,
) -> Result<Option<String>, ConversionError>
where
    P: ReferenceProvider + ?Sized,
{
    let Some((start, end, invert)) = position_range_insert(seq) else {
        return Ok(None);
    };
    read_reference_span(accession, start, end, invert, alphabet, provider).map(Some)
}

/// An inserted sequence whose bases are an exact tandem repeat, decomposed to
/// `(unit, count)` — the literal unit spelled once and the exact number of
/// copies — or `None` when `seq` is not one.
///
/// Recognises the single-base `Repeat` (`insA[10]`), the multi-base
/// `SequenceRepeat` (`insAT[3]`), and the single-part `Complex([Repeat])` shape
/// the parser produces for a lone repeat unit. Only an [`RepeatCount::Exact`]
/// count is a determinate string of bases; every uncertain or range count
/// (`insA[10_15]`, `insA[?]`, …) names no single expansion and returns `None`,
/// so the caller keeps declining it. Multi-part `Complex` payloads and every
/// non-repeat variant also return `None`.
fn exact_repeat_insert(seq: &InsertedSequence) -> Option<(String, u64)> {
    use crate::hgvs::edit::InsertedPart;
    let exact = |count: &RepeatCount| match count {
        RepeatCount::Exact(n) => Some(*n),
        _ => None,
    };
    match seq {
        InsertedSequence::Repeat { base, count } => {
            Some((base.to_char().to_string(), exact(count)?))
        }
        InsertedSequence::SequenceRepeat { sequence, count } => {
            Some((sequence_to_string(sequence), exact(count)?))
        }
        InsertedSequence::Complex(parts) => match parts.as_slice() {
            [InsertedPart::Repeat { base, count }] => {
                Some((base.to_char().to_string(), exact(count)?))
            }
            _ => None,
        },
        _ => None,
    }
}

/// Resolve an exact tandem-repeat insert (`insA[10]`, `insAT[3]`) to its literal
/// bases by expanding the unit the named number of times, mirroring how the
/// short-form `Repeat` edit arm expands its own explicit unit and exact count.
///
/// Returns `Ok(None)` when `seq` is not an exact repeat, so the caller can fall
/// back to literal handling. The expansion itself is [`expand_repeat_unit`],
/// shared with the bracket-part path so both apply the same overflow and cap.
fn resolve_exact_repeat_insert(
    seq: &InsertedSequence,
    alphabet: AlphabetMode,
) -> Result<Option<String>, ConversionError> {
    let Some((unit, count)) = exact_repeat_insert(seq) else {
        return Ok(None);
    };
    expand_repeat_unit(&unit, count, alphabet).map(Some)
}

/// Resolve one part of a compound insert to its literal bases, dispatching each
/// determinate shape to the same leaf helper the single-part path uses.
///
/// A `PositionRange`/`PositionRangeInv` reads through [`read_reference_span`], an
/// exact `Repeat` expands through [`expand_repeat_unit`], and a `Literal` is
/// rendered directly — the same three leaves that back
/// [`resolve_position_range_insert`] and [`resolve_exact_repeat_insert`], so a
/// part spells its bases identically whether it stands alone or inside a bracket.
/// The remaining shapes name bases this path cannot yet read out and decline as
/// `UnsupportedEditType`: a non-exact repeat count (`N[2800]`, `A[10_15]`), a CDS
/// position range carrying intronic offsets, or an external reference SPDI
/// cannot dereference.
fn resolve_inserted_part<P>(
    part: &crate::hgvs::edit::InsertedPart,
    accession: &str,
    alphabet: AlphabetMode,
    provider: Option<&P>,
) -> Result<String, ConversionError>
where
    P: ReferenceProvider + ?Sized,
{
    use crate::hgvs::edit::InsertedPart;

    match part {
        InsertedPart::Literal(seq) => Ok(apply_alphabet(&seq.to_string(), alphabet)),
        InsertedPart::PositionRange { start, end } => {
            read_reference_span(accession, *start, *end, false, alphabet, provider)
        }
        InsertedPart::PositionRangeInv { start, end } => {
            read_reference_span(accession, *start, *end, true, alphabet, provider)
        }
        InsertedPart::Repeat {
            base,
            count: RepeatCount::Exact(n),
        } => expand_repeat_unit(&base.to_char().to_string(), *n, alphabet),
        // A non-exact `Repeat` count (`N[2800]`, `A[10_15]`) names no determinate
        // expansion.
        InsertedPart::Repeat { base, count } => Err(ConversionError::UnsupportedEditType {
            description: format!(
                "compound insert part {base}{count} has no determinate expansion; \
                 this shape is not yet encodable as SPDI"
            ),
        }),
        InsertedPart::CdsPositionRange(range) => Err(ConversionError::UnsupportedEditType {
            description: format!(
                "compound insert part {range} names CDS positions with intronic offsets; \
                 this shape is not yet encodable as SPDI"
            ),
        }),
        InsertedPart::ExternalRef(reference) => Err(ConversionError::UnsupportedEditType {
            description: format!(
                "compound insert part {reference} references an external sequence SPDI \
                 cannot dereference; this shape is not yet encodable as SPDI"
            ),
        }),
    }
}

/// Resolve a multi-part `Complex` insert whose every part names determinate
/// bases — a literal, an exact tandem repeat, or a same-reference position
/// range (optionally inverted) — to the single string of bases they spell in
/// order, mirroring how `del`/`dup` resolve their omitted bases.
///
/// This is the bracketed-composition sibling of [`resolve_position_range_insert`]
/// and [`resolve_exact_repeat_insert`], which each handle only a *single*-part
/// `Complex`. A compound insert that mixes shapes — `delins[T;213_271]`, a
/// literal followed by a reference span — is one contiguous inserted sequence
/// once each part is read out, so it encodes to SPDI just as a plain literal
/// does. Each part is resolved by [`resolve_inserted_part`], which shares its
/// leaf helpers with the single-part path, so no shape logic is duplicated.
///
/// Returns `Ok(None)` when `seq` is not a `Complex`, so the caller falls back to
/// its single-shape resolvers and literal handling. Any part whose bases are
/// genuinely undetermined declines the whole insert as `UnsupportedEditType`.
fn resolve_complex_insert<P>(
    seq: &InsertedSequence,
    accession: &str,
    alphabet: AlphabetMode,
    provider: Option<&P>,
) -> Result<Option<String>, ConversionError>
where
    P: ReferenceProvider + ?Sized,
{
    let InsertedSequence::Complex(parts) = seq else {
        return Ok(None);
    };

    let mut resolved = String::new();
    for part in parts {
        resolved.push_str(&resolve_inserted_part(part, accession, alphabet, provider)?);
    }
    Ok(Some(resolved))
}

/// Helper to get start position from an interval.
fn get_start_pos(interval: &Interval<GenomePos>) -> Option<u64> {
    interval.start.inner().map(|p| p.base)
}

/// Helper to get end position from an interval.
fn get_end_pos(interval: &Interval<GenomePos>) -> Option<u64> {
    interval.end.inner().map(|p| p.base)
}

/// Every [`GenomePos`] a boundary actually **writes down**, whether or not any
/// one of them is a coordinate the conversion could use.
///
/// [`UncertainBoundary::inner`] answers a deliberately narrower question — "is
/// this one position I can resolve?" — and so returns `None` for a `Range`,
/// hiding *both* its endpoints. That is the right answer for a caller asking
/// for a coordinate and the wrong one for a caller asking what the description
/// says, and [`reject_unresolvable_genomic_position`] is the second kind: a
/// prohibited offset is prohibited wherever it is written, so
/// `g.10_(20+1_30)delACGT` carried one straight past the guard and converted
/// as though it read `g.10delACGT` (#1795).
fn stated_positions(boundary: &UncertainBoundary<GenomePos>) -> [Option<&GenomePos>; 2] {
    match boundary {
        UncertainBoundary::Single(mu) => [mu.inner(), None],
        UncertainBoundary::Range { start, end } => [start.inner(), end.inner()],
    }
}

/// Resolve a genomic interval's **end** coordinate, declining when the boundary
/// names no single position.
///
/// The call sites previously wrote `get_end_pos(..).unwrap_or(start_pos)`,
/// which is the same silent drop [`reject_unresolvable_genomic_position`] was
/// written to stop, one field over: an end that resolves to nothing — a
/// `(a_b)` range, or a bare `?` — was replaced by the start, so the interval
/// collapsed to a point and the boundary vanished from the answer. That is what
/// made `g.10_(20+1_30)delACGT`, `g.10_(20_30)delACGT` and `g.10delACGT` one
/// triple; the offset is a symptom of the dropped boundary, and closing only
/// the guard hole would have left the last two of those three still colliding.
///
/// Declining is not a new contract. The **start** side has always declined
/// (`get_start_pos(..).ok_or_else(..)`), and the VCF sibling
/// (`vcf::from_hgvs`) declines on both endpoints with no fallback at all — the
/// fallback here was the asymmetry, not the refusal. An SPDI position is an
/// exact interbase coordinate, so there is nothing an uncertain region or an
/// unknown position could be rendered as.
///
/// Note this is about a boundary that names *no* position, never about one the
/// author flagged as approximate: `(13)` is `Mu::Uncertain`, `inner()` returns
/// 13, and it resolves here exactly as a bare `13` does.
///
/// **Genomic axes only.** The `c.`/`n.`/`r.` paths carried the same fallback at
/// five further sites and collapsed the same way — measured, and deferred to
/// #1804 rather than swept in here, because those axes decide the *neighbouring*
/// question differently: a `+`/`-` offset inside a `Range` is *legitimate*
/// there (`c.(4185+1_4186-1)_(4357+1_4358-1)del`), so
/// [`reject_unresolvable_genomic_position`] cannot simply be pointed at them.
///
/// #1804 has since been closed by #1823, whose
/// [`resolve_transcript_end_boundary`] is this function's transcript-axis
/// counterpart. The two stay separate for a typing reason rather than a policy
/// one — this resolves a [`GenomePos`] to the `u64` its call sites need, while
/// that one hands the axis' own position type back for a further conversion —
/// and they agree on the verdict: both refuse with
/// [`ConversionError::InvalidPosition`]. See that function's doc for why it
/// takes that verdict from its own axis rather than inheriting it from here.
fn resolve_genomic_end_pos(
    interval: &Interval<GenomePos>,
    coord: &str,
) -> Result<u64, ConversionError> {
    get_end_pos(interval).ok_or_else(|| ConversionError::InvalidPosition {
        description: format!(
            "{coord}. end boundary `{}` names no single coordinate: an SPDI position is an \
             exact interbase coordinate, and an uncertain range or unknown (`?`) boundary \
             supplies none. Collapsing it onto the start position would make distinct \
             descriptions share one triple. Give the end a definite position, or keep the \
             description in HGVS",
            interval.end
        ),
    })
}

/// Refuse a genomic interval whose start or end holds something SPDI cannot be
/// given a coordinate for: a `+`/`-` offset (#1628) or a `pter`/`qter`/`cen`
/// special position (#1643).
///
/// The two are the same defect in two fields of one struct. `GenomePos` can
/// hold either because the parser accepts both, and neither half of the
/// conversion can honour them:
///
/// - an **offset** has nothing on a genomic accession to be measured against —
///   there is no exon structure, and `checklist.md:16` prohibits an offset on a
///   genomic position outright — and SPDI has no offset notation;
/// - a **special position** names a landmark of the assembled chromosome
///   (`pter`, `qter`) or of its centromere annotation (`cen`), none of which a
///   sequence accession carries. `GenomePos` stores `base: 0` beside the
///   `special` marker precisely because there is no coordinate to store.
///
/// What the conversion did instead was drop both, which is the one answer worse
/// than either honest option, and it produced the same confluence failure
/// twice. `g.266+2del`, `g.266-268del` and `g.266del` flattened onto one triple
/// while `normalize` keeps them as three distinct descriptions. So did
/// `g.10_qterdelACGTACGTAC`, `g.10_cendelACGTACGTAC` and
/// `g.10_19delACGTACGTAC` — a deletion to the q-arm telomere, one to the
/// centromere and a literal ten-base deletion, all `NC_000001.11:9:ACGTACGTAC:`.
///
/// **The special-position half was reachable only when the description spells
/// its own bases**, which is why it survived #1641 and why the guard is what
/// closes it rather than the checks that appear to. Every other shape happened
/// to be refused *incidentally*, by a check asking a different question:
/// `g.pter_10del` by the 1-based position check ("position 0 is not valid in
/// HGVS"), `g.10_qterdel` by the reference fetch (`invalid 1-based interval
/// [10, 0]`), and the `m.`/`o.` forms by the circular-wraparound check, which
/// reads `start > end` and reports a wraparound that is not there. Those
/// messages diagnose the wrong thing and none of them is a guarantee — a
/// description carrying its own bases needs no fetch and no position check, and
/// so converted.
///
/// **Both endpoints of a complex `(a_b)` boundary are walked, not just the
/// boundary as a whole** (#1795). The guard originally read each side through
/// [`UncertainBoundary::inner`], which reports `None` for a `Range`, so a
/// prohibited offset or special position written *inside* the parentheses was
/// never inspected: `g.10_(20+1_30)delACGT` and `g.10_(20_qter)delACGT` both
/// converted. See [`stated_positions`], which is the walk this uses and is
/// deliberately a different question from `inner()`.
///
/// Refusing is also what the other axes already do with the offsets that *are*
/// legitimate there — see [`resolve_cds_to_tx`], [`resolve_tx_pos`] and
/// [`require_simple_tx_pos`], each of which declines an intronic `c.`/`n.`/`r.`
/// position for want of an SPDI representation. Those decline because the
/// position cannot be projected onto the transcript accession; this one declines
/// because there is nothing on a genomic accession to project it against.
fn reject_unresolvable_genomic_position(
    interval: &Interval<GenomePos>,
    coord: &str,
) -> Result<(), ConversionError> {
    for endpoint in [&interval.start, &interval.end]
        .into_iter()
        .flat_map(stated_positions)
        .flatten()
    {
        if endpoint.offset.is_some() {
            return Err(ConversionError::InvalidPosition {
                description: format!(
                    "{coord}. position {endpoint} carries a +/- offset: a genomic position \
                     cannot have one and SPDI has no offset notation to express it. Drop the \
                     offset, or describe the variant on a transcript accession where an offset \
                     is meaningful"
                ),
            });
        }
        if let Some(special) = endpoint.special {
            return Err(ConversionError::InvalidPosition {
                description: format!(
                    "{coord}. position `{special}` names no numeric coordinate: `pter`, `qter` \
                     and `cen` are landmarks of the assembled chromosome, not positions on the \
                     sequence, and SPDI has no notation for them. Give the position a number, \
                     or keep the description in HGVS"
                ),
            });
        }
    }
    Ok(())
}

/// Convert an HGVS variant to SPDI format without consulting a reference provider.
///
/// This is the "simple" conversion path: the SPDI is emitted on the same
/// accession the HGVS variant uses (no genomic projection of transcript
/// variants). It accepts coordinate systems whose positions can be resolved
/// without transcript metadata:
///
/// | HGVS coord | Supported here | Notes |
/// |------------|----------------|-------|
/// | `g.` (genomic) | yes | direct 1→0-based conversion; a `+`/`-` offset (#1628) or a `pter`/`qter`/`cen` (#1643) is refused wherever it is written, inside a complex `(a_b)` boundary included, and an end boundary naming no single coordinate — `(a_b)` or `?` — is refused rather than collapsed onto the start (#1795) |
/// | `m.` (mito) | yes | mito accession is genomic; same as `g.`; wraparound rejected |
/// | `o.` (circular) | yes | same path as `g.`/`m.`; wraparound rejected |
/// | `n.` (non-coding tx) | exonic, positive base | SPDI sits on the transcript accession |
/// | `r.` (RNA) | exonic, positive base | `u`/`U` rewritten to `T`; SPDI uses DNA alphabet |
/// | `c.` (CDS) | NO — needs CDS start | use [`hgvs_to_spdi`] |
/// | `p.` (protein) | NO | not representable in SPDI |
///
/// For `c.`, UTR-style `n.`/`r.`, or deletion/dup variants without an
/// explicit deleted sequence, use [`hgvs_to_spdi`] which consults a
/// [`ReferenceProvider`]. Intronic `n.`/`r.` positions are not supported
/// by either entry point — SPDI has no offset notation, and genomic
/// projection is future work; both functions return
/// [`ConversionError::MissingReferenceData`]. A `g.`/`m.`/`o.` position
/// carrying an offset (#1628) or a `pter`/`qter`/`cen` special position
/// (#1643) is refused outright by both entry points with
/// [`ConversionError::InvalidPosition`]: neither has anything on a genomic
/// accession to resolve it against, and dropping either made distinct
/// descriptions share one triple. That holds **wherever the position is
/// written**, including inside a complex `(a_b)` boundary (#1795).
///
/// A `g.`/`m.`/`o.` **end** boundary that names no single coordinate — a
/// `(a_b)` range or a bare `?` — is likewise refused with
/// [`ConversionError::InvalidPosition`] rather than collapsed onto the start
/// position, which is what the start side has always done. An uncertain but
/// numeric position (`(13)`) names a coordinate and still converts. On the
/// **transcript** axes, a `c.`/`n.`/`r.` interval whose **end** names no single
/// coordinate is refused with the same variant (#1804); see
/// [`resolve_transcript_end_boundary`] for why those axes need their own
/// resolver rather than the genomic guard.
///
/// # Arguments
///
/// * `variant` - The HGVS variant to convert.
///
/// # Returns
///
/// * `Ok(SpdiVariant)` - Successfully converted variant.
/// * `Err(ConversionError::ProviderRequired)` - The variant requires a
///   provider (typically `c.` or intronic `n.`/`r.`).
/// * `Err(ConversionError)` - Conversion failed for another reason.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::spdi::convert::hgvs_to_spdi_simple;
/// use ferro_hgvs::parse_hgvs;
///
/// // Genomic variant
/// let hgvs = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
/// let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
/// assert_eq!(spdi.to_string(), "NC_000001.11:12344:A:G");
///
/// // Non-coding transcript: SPDI emitted on the transcript accession
/// let hgvs = parse_hgvs("NR_046018.2:n.5C>G").unwrap();
/// let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
/// assert_eq!(spdi.to_string(), "NR_046018.2:4:C:G");
/// ```
///
/// [`ReferenceProvider`]: crate::reference::provider::ReferenceProvider
pub fn hgvs_to_spdi_simple(variant: &HgvsVariant) -> Result<SpdiVariant, ConversionError> {
    match variant {
        HgvsVariant::Genome(g) => genome_to_spdi_simple(g),
        HgvsVariant::Mt(m) => mt_to_spdi_simple(m),
        HgvsVariant::Circular(o) => circular_to_spdi_simple(o),
        HgvsVariant::Tx(n) => tx_to_spdi_simple(n),
        HgvsVariant::Rna(r) => rna_to_spdi_simple(r),
        HgvsVariant::Cds(_) => Err(ConversionError::ProviderRequired {
            variant_type: "c".to_string(),
            reason:
                "CDS positions need transcript metadata (CDS start) to resolve to a transcript \
                 position; call hgvs_to_spdi with a ReferenceProvider"
                    .to_string(),
        }),
        HgvsVariant::Protein(_) => Err(ConversionError::UnsupportedVariantType {
            description: "protein variants cannot be represented in SPDI; SPDI describes \
                          nucleotide variants on a sequence accession"
                .to_string(),
        }),
        _ => Err(ConversionError::UnsupportedVariantType {
            description: format!(
                "variant type {} cannot be converted to SPDI",
                variant.variant_type()
            ),
        }),
    }
}

/// Convert an HGVS variant to SPDI, consulting a reference provider for
/// transcript metadata (CDS start, exon coordinates) and for the deleted /
/// duplicated bases of short-form `Deletion`, `Duplication`, and `Delins`
/// edits.
///
/// This is the provider-aware companion to [`hgvs_to_spdi_simple`]. It
/// handles `c.` (CDS) variants, `n.`/`r.` variants with UTR-style positions
/// (`*N` downstream), and populates SPDI's mandatory `del` field for
/// short-form deletions / delins / identities (and the symmetric `ins` field
/// for short-form duplications) by fetching the reference bases for the
/// variant's interval. Explicit-form input (`g.100_102delATG`,
/// `g.100_102dupATG`, `g.100A=`, etc.) emits the user-supplied bases as-is
/// and does not consult the provider. Intronic `n.`/`r.` positions remain
/// unsupported (SPDI has no offset notation) and return
/// [`ConversionError::MissingReferenceData`]. An offset or a `pter`/`qter`/`cen`
/// special position on a `g.`/`m.`/`o.` position is refused with
/// [`ConversionError::InvalidPosition`] — see
/// [`reject_unresolvable_genomic_position`], which reads inside a complex
/// `(a_b)` boundary as well as a simple one. A genomic end boundary that names
/// no single coordinate is refused by [`resolve_genomic_end_pos`] rather than
/// collapsed onto the start. A `c.`/`n.`/`r.` interval whose **end** names no
/// single coordinate — a range boundary `(20_30)` or an unknown `?` — is
/// refused with the same variant (#1804); see
/// [`resolve_transcript_end_boundary`].
///
/// An **unspelled identity** (`g.100=`, `g.100_102=`) therefore needs a
/// provider: on [`hgvs_to_spdi_simple`] it returns
/// [`ConversionError::MissingReferenceData`] rather than a zero-width triple
/// that would claim no bases. A **whole-entity** identity (`g.=`) names no
/// interval at all and returns [`ConversionError::UnsupportedEditType`].
///
/// The resulting SPDI uses the **same accession** as the HGVS variant — for
/// `NM_000088.3:c.1A>G` the SPDI sequence is `NM_000088.3`, not the
/// underlying genomic accession. This matches NCBI Variation Services'
/// behavior: SPDI is positional on whichever accession is provided. No axis is
/// projected onto another.
///
/// The one thing that is *not* carried over verbatim is a **compound
/// reference's genomic-context wrapper** on a transcript axis: for
/// `NG_008939.1(NM_000532.5):c.156A>G` the SPDI sequence is `NM_000532.5`. That
/// is not a projection either — it is the accession this function already
/// resolved the position against, and naming the wrapper instead pointed every
/// downstream fetch at the genomic parent. See [`transcript_axis_sequence`]. A
/// `g.`/`m.`/`o.` description keeps its accession exactly as written.
///
/// # Arguments
///
/// * `variant` - The HGVS variant to convert.
/// * `provider` - A reference provider used both for transcript metadata
///   (`get_transcript`) and reference-base fetches for short-form edits
///   (`get_genomic_sequence` first, then `get_sequence` as fallback).
///
/// # Errors
///
/// * [`ConversionError::MissingReferenceData`] if the provider does not
///   have the transcript or the requested reference interval, or if the
///   position cannot be resolved (e.g. intronic without exon data).
/// * Other `ConversionError` variants for the same reasons as
///   [`hgvs_to_spdi_simple`].
///
/// # Examples
///
/// Short-form deletion with provider-backed ref fetch:
///
/// ```
/// use ferro_hgvs::spdi::convert::hgvs_to_spdi;
/// use ferro_hgvs::reference::mock::MockProvider;
/// use ferro_hgvs::parse_hgvs;
///
/// let mut provider = MockProvider::new();
/// // Build a contig where 1-based 100..102 = "ATG".
/// let mut seq = "N".repeat(99);
/// seq.push_str("ATG");
/// provider.add_genomic_sequence("NC_000001.11", &seq);
///
/// let hgvs = parse_hgvs("NC_000001.11:g.100_102del").unwrap();
/// let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
/// assert_eq!(spdi.to_string(), "NC_000001.11:99:ATG:");
/// ```
pub fn hgvs_to_spdi<P: ReferenceProvider + ?Sized>(
    variant: &HgvsVariant,
    provider: &P,
) -> Result<SpdiVariant, ConversionError> {
    match variant {
        HgvsVariant::Genome(g) => genome_to_spdi_with_provider(g, provider),
        HgvsVariant::Mt(m) => mt_to_spdi_with_provider(m, provider),
        HgvsVariant::Circular(o) => circular_to_spdi_with_provider(o, provider),
        HgvsVariant::Tx(n) => tx_to_spdi_with_provider(n, provider),
        HgvsVariant::Rna(r) => rna_to_spdi_with_provider(r, provider),
        HgvsVariant::Cds(c) => cds_to_spdi_with_provider(c, provider),
        HgvsVariant::Protein(_) => Err(ConversionError::UnsupportedVariantType {
            description: "protein variants cannot be represented in SPDI; SPDI describes \
                          nucleotide variants on a sequence accession"
                .to_string(),
        }),
        _ => Err(ConversionError::UnsupportedVariantType {
            description: format!(
                "variant type {} cannot be converted to SPDI",
                variant.variant_type()
            ),
        }),
    }
}

/// Convert a genomic variant to SPDI (simple conversion).
fn genome_to_spdi_simple(variant: &GenomeVariant) -> Result<SpdiVariant, ConversionError> {
    reject_unresolvable_genomic_position(&variant.loc_edit.location, "g")?;
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_pos = get_start_pos(&variant.loc_edit.location).ok_or_else(|| {
        ConversionError::InvalidPosition {
            description: "cannot convert variant with unknown start position".to_string(),
        }
    })?;
    let end_pos = resolve_genomic_end_pos(&variant.loc_edit.location, "g")?;
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_pos,
        end_pos,
        edit,
        AlphabetMode::Dna,
        None::<&dyn ReferenceProvider>,
    )
}

/// Convert a mitochondrial variant to SPDI (simple conversion).
///
/// Mitochondrial accessions (e.g. `NC_012920.1`) are themselves genomic
/// accessions, so the conversion is identical to the `g.` path with a
/// different coordinate prefix on the HGVS side.
///
/// Wraparound variants (start > end, per SVD-WG006) are rejected: SPDI is a
/// single-edit format with no native representation for circular-contig
/// wraparound.
fn mt_to_spdi_simple(variant: &MtVariant) -> Result<SpdiVariant, ConversionError> {
    reject_unresolvable_genomic_position(&variant.loc_edit.location, "m")?;
    if let (Some(s), Some(e)) = (
        get_start_pos(&variant.loc_edit.location),
        get_end_pos(&variant.loc_edit.location),
    ) {
        if s > e {
            return Err(ConversionError::InvalidPosition {
                description: format!(
                    "Cannot convert wraparound m. variant to SPDI: SPDI is a single edit and \
                     has no representation for circular-contig wraparound. Variant accession: {}",
                    variant.accession
                ),
            });
        }
    }
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_pos = get_start_pos(&variant.loc_edit.location).ok_or_else(|| {
        ConversionError::InvalidPosition {
            description: "cannot convert variant with unknown start position".to_string(),
        }
    })?;
    let end_pos = resolve_genomic_end_pos(&variant.loc_edit.location, "m")?;
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_pos,
        end_pos,
        edit,
        AlphabetMode::Dna,
        None::<&dyn ReferenceProvider>,
    )
}

/// The accession a **transcript-axis** (`c.`/`n.`/`r.`) SPDI is emitted on: the
/// inner transcript, with any `NG_`/`NC_`/`LRG_` genomic-context wrapper
/// stripped.
///
/// This is the same string the coordinate resolution already uses — see
/// [`resolve_cds_to_tx`], which calls `accession.transcript_accession()`, reads
/// that transcript's `cds_start`, and bounds the result against *that*
/// transcript's length. So the number `emit_spdi_for_edit` receives is an offset
/// on the transcript, and this makes the accession beside it say so.
///
/// **Why not `Display`.** `Accession`'s `Display` renders the compound form
/// (`NG_008939.1(NM_000532.5)`), and that string is both the SPDI's `sequence`
/// field and the key every short-form edit's reference fetch is made against.
/// `MultiFastaProvider::resolve_name` reaches a record for it through the
/// version-strip fallback — `"NG_008939.1(NM_000532.5)".split('.').next()` is
/// `"NG_008939"` — so a **transcript** offset was read out of the **genomic
/// parent**. Measured on the prepared reference:
///
/// ```text
/// get_sequence("NG_008939.1",              191, 197) = CACACA   (the parent)
/// get_sequence("NG_008939.1(NM_000532.5)", 191, 197) = CACACA   <-- same record
/// get_sequence("NM_000532.5",              191, 197) = ACGCCG   <-- the 191's own frame
/// ```
///
/// Splitting the label from the fetch key would not close it: the SPDI's
/// `sequence` field *is* the fetch key for every downstream consumer
/// ([`crate::spdi::apply_to_reference`], and `fetch_window` inside
/// [`crate::spdi::compare_denoted_sequences`]), so a compound label just moves
/// the wrong fetch one hop out. The compound spelling is also not a valid SPDI
/// sequence in the first place — it is HGVS syntax, not an accession, and SPDI's
/// `sequence` must be the accession the position is an offset on.
///
/// **Genomic axes deliberately keep `Display`.** On a `g.`/`m.`/`o.` description
/// the coordinates are the *parent's*, so stripping the wrapper here would name
/// the transcript for a genomic offset — the same defect with the frames
/// swapped. Those paths are untouched.
fn transcript_axis_sequence(accession: &Accession) -> String {
    accession.transcript_accession()
}

/// Convert a non-coding transcript (`n.`) variant to SPDI without consulting
/// a provider. The SPDI is emitted on the transcript accession.
///
/// Returns `MissingReferenceData` for cases that need provider-backed
/// metadata: intronic offsets, downstream (`*N`) positions, and non-positive
/// bases (5' UTR). Use [`hgvs_to_spdi`] with a provider for those.
fn tx_to_spdi_simple(variant: &TxVariant) -> Result<SpdiVariant, ConversionError> {
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_tx = tx_pos_for_simple_path(&variant.loc_edit.location, "n")?;
    let end_tx = tx_end_for_simple_path(&variant.loc_edit.location, "n")?;
    emit_spdi_for_edit(
        transcript_axis_sequence(&variant.accession),
        start_tx,
        end_tx,
        edit,
        AlphabetMode::Dna,
        None::<&dyn ReferenceProvider>,
    )
}

/// Convert an RNA (`r.`) variant to SPDI without consulting a provider.
///
/// Identical to [`tx_to_spdi_simple`] in terms of position handling. The
/// edit's deletion/insertion sequences are rewritten with `u`/`U → T` so the
/// output uses the DNA alphabet that SPDI uses by convention (RefSeq stores
/// transcript sequences as DNA even on `NR_*` and `NM_*` accessions).
fn rna_to_spdi_simple(variant: &RnaVariant) -> Result<SpdiVariant, ConversionError> {
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_pos = rna_pos_for_simple_path(&variant.loc_edit.location, "r")?;
    let end_pos = rna_end_for_simple_path(&variant.loc_edit.location, "r")?;
    emit_spdi_for_edit(
        transcript_axis_sequence(&variant.accession),
        start_pos,
        end_pos,
        edit,
        AlphabetMode::Rna,
        None::<&dyn ReferenceProvider>,
    )
}

/// Convert a genomic variant to SPDI with provider-backed reference fetch.
///
/// Same as [`genome_to_spdi_simple`] for substitution / insertion / identity,
/// but uses the provider to populate the `del` field for short-form
/// `Deletion`, `Duplication`, and `Delins` edits.
fn genome_to_spdi_with_provider<P: ReferenceProvider + ?Sized>(
    variant: &GenomeVariant,
    provider: &P,
) -> Result<SpdiVariant, ConversionError> {
    reject_unresolvable_genomic_position(&variant.loc_edit.location, "g")?;
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_pos = get_start_pos(&variant.loc_edit.location).ok_or_else(|| {
        ConversionError::InvalidPosition {
            description: "cannot convert variant with unknown start position".to_string(),
        }
    })?;
    let end_pos = resolve_genomic_end_pos(&variant.loc_edit.location, "g")?;
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_pos,
        end_pos,
        edit,
        AlphabetMode::Dna,
        Some(provider),
    )
}

/// Convert a mitochondrial variant to SPDI with provider-backed reference fetch.
///
/// The mito accession is genomic, so the path mirrors
/// [`genome_to_spdi_with_provider`].
///
/// Wraparound variants (start > end, per SVD-WG006) are rejected: SPDI is a
/// single-edit format with no native representation for circular-contig
/// wraparound.
fn mt_to_spdi_with_provider<P: ReferenceProvider + ?Sized>(
    variant: &MtVariant,
    provider: &P,
) -> Result<SpdiVariant, ConversionError> {
    reject_unresolvable_genomic_position(&variant.loc_edit.location, "m")?;
    if let (Some(s), Some(e)) = (
        get_start_pos(&variant.loc_edit.location),
        get_end_pos(&variant.loc_edit.location),
    ) {
        if s > e {
            return Err(ConversionError::InvalidPosition {
                description: format!(
                    "Cannot convert wraparound m. variant to SPDI: SPDI is a single edit and \
                     has no representation for circular-contig wraparound. Variant accession: {}",
                    variant.accession
                ),
            });
        }
    }
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_pos = get_start_pos(&variant.loc_edit.location).ok_or_else(|| {
        ConversionError::InvalidPosition {
            description: "cannot convert variant with unknown start position".to_string(),
        }
    })?;
    let end_pos = resolve_genomic_end_pos(&variant.loc_edit.location, "m")?;
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_pos,
        end_pos,
        edit,
        AlphabetMode::Dna,
        Some(provider),
    )
}

/// Convert a circular (`o.`) variant to SPDI (simple conversion).
///
/// Circular accessions follow the same coordinate layout as genomic accessions,
/// so the conversion mirrors the `g.`/`m.` path.
///
/// Wraparound variants (start > end, per SVD-WG006) are rejected: SPDI is a
/// single-edit format with no native representation for circular-contig
/// wraparound.
fn circular_to_spdi_simple(variant: &CircularVariant) -> Result<SpdiVariant, ConversionError> {
    reject_unresolvable_genomic_position(&variant.loc_edit.location, "o")?;
    if let (Some(s), Some(e)) = (
        get_start_pos(&variant.loc_edit.location),
        get_end_pos(&variant.loc_edit.location),
    ) {
        if s > e {
            return Err(ConversionError::InvalidPosition {
                description: format!(
                    "Cannot convert wraparound o. variant to SPDI: SPDI is a single edit and \
                     has no representation for circular-contig wraparound. Variant accession: {}",
                    variant.accession
                ),
            });
        }
    }
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_pos = get_start_pos(&variant.loc_edit.location).ok_or_else(|| {
        ConversionError::InvalidPosition {
            description: "cannot convert variant with unknown start position".to_string(),
        }
    })?;
    let end_pos = resolve_genomic_end_pos(&variant.loc_edit.location, "o")?;
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_pos,
        end_pos,
        edit,
        AlphabetMode::Dna,
        None::<&dyn ReferenceProvider>,
    )
}

/// Convert a circular (`o.`) variant to SPDI with provider-backed reference fetch.
///
/// The circular accession is genomic, so the path mirrors
/// [`genome_to_spdi_with_provider`].
///
/// Wraparound variants (start > end, per SVD-WG006) are rejected: SPDI is a
/// single-edit format with no native representation for circular-contig
/// wraparound.
fn circular_to_spdi_with_provider<P: ReferenceProvider + ?Sized>(
    variant: &CircularVariant,
    provider: &P,
) -> Result<SpdiVariant, ConversionError> {
    reject_unresolvable_genomic_position(&variant.loc_edit.location, "o")?;
    if let (Some(s), Some(e)) = (
        get_start_pos(&variant.loc_edit.location),
        get_end_pos(&variant.loc_edit.location),
    ) {
        if s > e {
            return Err(ConversionError::InvalidPosition {
                description: format!(
                    "Cannot convert wraparound o. variant to SPDI: SPDI is a single edit and \
                     has no representation for circular-contig wraparound. Variant accession: {}",
                    variant.accession
                ),
            });
        }
    }
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_pos = get_start_pos(&variant.loc_edit.location).ok_or_else(|| {
        ConversionError::InvalidPosition {
            description: "cannot convert variant with unknown start position".to_string(),
        }
    })?;
    let end_pos = resolve_genomic_end_pos(&variant.loc_edit.location, "o")?;
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_pos,
        end_pos,
        edit,
        AlphabetMode::Dna,
        Some(provider),
    )
}

/// Convert a CDS (`c.`) variant to SPDI by resolving CDS coordinates to
/// transcript positions through the supplied provider.
///
/// The resulting SPDI uses the variant's transcript accession (e.g.
/// `NM_000088.3`), matching NCBI Variation Services' convention. Short-form
/// `Deletion` / `Duplication` / `Delins` edits trigger a provider fetch on
/// the transcript accession to populate SPDI's `del` field.
fn cds_to_spdi_with_provider<P: ReferenceProvider + ?Sized>(
    variant: &CdsVariant,
    provider: &P,
) -> Result<SpdiVariant, ConversionError> {
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_cds = variant.loc_edit.location.start.inner().ok_or_else(|| {
        ConversionError::InvalidPosition {
            description: "cannot convert c. variant with unknown start position".to_string(),
        }
    })?;
    let end_cds = resolve_transcript_end_boundary(&variant.loc_edit.location, "c")?;
    let (start_tx, end_tx) = resolve_cds_to_tx(&variant.accession, start_cds, end_cds, provider)?;
    emit_spdi_for_edit(
        transcript_axis_sequence(&variant.accession),
        start_tx,
        end_tx,
        edit,
        AlphabetMode::Dna,
        Some(provider),
    )
}

/// Convert an `n.` variant with provider-backed reference fetch and
/// transcript-aware position resolution (intronic offsets, downstream
/// `*N`, non-positive base).
fn tx_to_spdi_with_provider<P: ReferenceProvider + ?Sized>(
    variant: &TxVariant,
    provider: &P,
) -> Result<SpdiVariant, ConversionError> {
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let (start_tx, end_tx) = if tx_needs_provider(&variant.loc_edit.location) {
        let start = variant.loc_edit.location.start.inner().ok_or_else(|| {
            ConversionError::InvalidPosition {
                description: "cannot convert n. variant with unknown start position".to_string(),
            }
        })?;
        let end = resolve_transcript_end_boundary(&variant.loc_edit.location, "n")?;
        resolve_tx_to_provider_tx(&variant.accession, start, end, provider)?
    } else {
        resolve_tx_exonic_bounded(&variant.accession, &variant.loc_edit.location, provider)?
    };
    emit_spdi_for_edit(
        transcript_axis_sequence(&variant.accession),
        start_tx,
        end_tx,
        edit,
        AlphabetMode::Dna,
        Some(provider),
    )
}

/// Convert an `r.` variant with provider-backed reference fetch. Same
/// coordinate resolution as `n.`; alphabet conversion `u → T` is applied
/// via [`AlphabetMode::Rna`].
fn rna_to_spdi_with_provider<P: ReferenceProvider + ?Sized>(
    variant: &RnaVariant,
    provider: &P,
) -> Result<SpdiVariant, ConversionError> {
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let (start_tx, end_tx) = if rna_needs_provider(&variant.loc_edit.location) {
        let start = variant.loc_edit.location.start.inner().ok_or_else(|| {
            ConversionError::InvalidPosition {
                description: "cannot convert r. variant with unknown start position".to_string(),
            }
        })?;
        let end = resolve_transcript_end_boundary(&variant.loc_edit.location, "r")?;
        resolve_rna_to_provider_tx(&variant.accession, start, end, provider)?
    } else {
        resolve_rna_exonic_bounded(&variant.accession, &variant.loc_edit.location, provider)?
    };
    emit_spdi_for_edit(
        transcript_axis_sequence(&variant.accession),
        start_tx,
        end_tx,
        edit,
        AlphabetMode::Rna,
        Some(provider),
    )
}

// ===========================================================================
// Shared helpers
// ===========================================================================

/// Whether to rewrite RNA bases (`u/U`) to DNA (`T`) in the emitted SPDI
/// deletion/insertion strings.
///
/// `pub(crate)` for [`crate::spdi::apply`], which reads a key's bases back out
/// of the reference window and so has to fold them in the same convention this
/// module's triples already use.
#[derive(Debug, Clone, Copy)]
pub(crate) enum AlphabetMode {
    Dna,
    Rna,
}

/// Unwrap an edit from `Mu`, returning a clear error if the edit is unknown.
fn unwrap_edit<E>(edit: &crate::hgvs::uncertainty::Mu<E>) -> Result<&E, ConversionError> {
    edit.inner()
        .ok_or_else(|| ConversionError::InvalidPosition {
            description: "cannot convert variant with unknown edit".to_string(),
        })
}

/// Resolve the **end** boundary of a `c.`/`n.`/`r.` interval to the single
/// position SPDI needs, declining when the boundary names none (#1804).
///
/// Two shapes reach here with nothing to resolve, and
/// [`crate::hgvs::interval::UncertainBoundary::inner`] answers `None` for both: a **range** boundary
/// (`(20_30)`, the `Range` variant) and an **unknown** one (`?`, `Mu::Unknown`).
/// Each says "the end is somewhere I am not telling you", and SPDI's end is one
/// exact interbase coordinate — there is no notation for either.
///
/// What every one of the seven call sites did instead was substitute the
/// **start**, which is the one answer worse than refusing: it silently
/// re-describes the variant as ending where it begins. Measured on `main` at
/// `439617c2` against `make_intronic_provider`, `c.10_(20_30)delAAAA` and
/// `c.10_?delAAAA` both converted to `NM_INTRON.1:19:AAAA:` — two descriptions
/// that `parse` and `Display` distinctly, sharing one triple with each other and
/// with any third description of that triple. The short forms are worse still:
/// `c.10_(20_30)del` became `NM_INTRON.1:19:A:`, a **one-base** deletion
/// standing in for one whose end is unknown.
///
/// **This is deliberately not the genomic guard again, and the difference is the
/// whole reason #1795 did not sweep these axes in.**
/// [`reject_unresolvable_genomic_position`] walks *both* endpoints and refuses
/// what it finds *inside* a resolvable one — a `+`/`-` offset, a
/// `pter`/`qter`/`cen` landmark — because neither may appear on a genomic
/// position at all. Here an offset is legitimate: `c.(4185+1_4186-1)_(4357+1_4358-1)del`
/// is well-formed HGVS and is the motivating example on
/// [`crate::hgvs::interval::UncertainBoundary`]. So the offset is not the defect
/// and this resolver never looks at one; what it refuses is the *absence* of a
/// single coordinate, which is a property of the boundary rather than of its
/// contents.
///
/// **The verdict is [`ConversionError::InvalidPosition`], and it is taken from
/// this axis rather than inherited from `g.`** — every one of the seven sites
/// already answers the same question for the **start** side, with `ok_or_else`
/// and exactly that variant ("cannot convert c. variant with unknown start
/// position"). The asymmetry between the two sides *is* the defect, so the end
/// gets what the start has always had. The axis' other verdict,
/// `MissingReferenceData`, would be a false statement here: it means "a
/// well-formed position whose SPDI form needs data or a projection I do not
/// have" (see [`resolve_cds_to_tx`], [`resolve_tx_pos`], [`require_simple_tx_pos`]),
/// and it invites the caller to supply a provider. No provider, exon table or
/// genomic projection can resolve `?` — the uncertainty is in the description,
/// not in ferro's reference data.
fn resolve_transcript_end_boundary<'a, T: std::fmt::Display>(
    interval: &'a Interval<T>,
    coord: &str,
) -> Result<&'a T, ConversionError> {
    interval
        .end
        .inner()
        .ok_or_else(|| ConversionError::InvalidPosition {
            description: format!(
                "{coord}. interval end `{}` names no single coordinate: an end that is a range \
                 `(a_b)` or unknown `?` states where the variant ends only approximately, and an \
                 SPDI triple ends at one exact position. Substituting the start would describe a \
                 different variant — one ending where it begins. Give the end a position, or keep \
                 the description in HGVS",
                interval.end
            ),
        })
}

/// Resolve the start position of a `TxInterval` for the simple (no-provider)
/// path. Returns `MissingReferenceData` if the position requires provider
/// data (intronic, downstream `*N`, or non-positive base).
fn tx_pos_for_simple_path(interval: &Interval<TxPos>, coord: &str) -> Result<u64, ConversionError> {
    let start = interval
        .start
        .inner()
        .ok_or_else(|| ConversionError::InvalidPosition {
            description: format!(
                "cannot convert {}. variant with unknown start position",
                coord
            ),
        })?;
    require_simple_tx_pos(start, coord)
}

fn tx_end_for_simple_path(interval: &Interval<TxPos>, coord: &str) -> Result<u64, ConversionError> {
    require_simple_tx_pos(resolve_transcript_end_boundary(interval, coord)?, coord)
}

// No-provider path: with no transcript there is no length to bound against,
// so an exonic position past the 3' end cannot be rejected here (#971). The
// provider-backed path (`hgvs_to_spdi`) does bound; `hgvs_to_spdi_simple`
// callers accept unbounded exonic positions by construction.
fn require_simple_tx_pos(pos: &TxPos, coord: &str) -> Result<u64, ConversionError> {
    if pos.is_intronic() {
        return Err(ConversionError::MissingReferenceData {
            description: format!(
                "intronic {}. position requires reference provider with exon data",
                coord
            ),
        });
    }
    if pos.is_downstream() {
        return Err(ConversionError::MissingReferenceData {
            description: format!(
                "downstream {}. position (*N) requires reference provider with transcript length",
                coord
            ),
        });
    }
    if pos.base < 1 {
        return Err(ConversionError::MissingReferenceData {
            description: format!(
                "non-positive {}. position {} requires reference provider with transcript length",
                coord, pos.base
            ),
        });
    }
    Ok(pos.base as u64)
}

fn rna_pos_for_simple_path(
    interval: &Interval<RnaPos>,
    coord: &str,
) -> Result<u64, ConversionError> {
    let start = interval
        .start
        .inner()
        .ok_or_else(|| ConversionError::InvalidPosition {
            description: format!(
                "cannot convert {}. variant with unknown start position",
                coord
            ),
        })?;
    require_simple_rna_pos(start, coord)
}

fn rna_end_for_simple_path(
    interval: &Interval<RnaPos>,
    coord: &str,
) -> Result<u64, ConversionError> {
    require_simple_rna_pos(resolve_transcript_end_boundary(interval, coord)?, coord)
}

// No-provider path: with no transcript there is no length to bound against,
// so an exonic position past the 3' end cannot be rejected here (#971). The
// provider-backed path (`hgvs_to_spdi`) does bound; `hgvs_to_spdi_simple`
// callers accept unbounded exonic positions by construction.
fn require_simple_rna_pos(pos: &RnaPos, coord: &str) -> Result<u64, ConversionError> {
    if pos.is_intronic() {
        return Err(ConversionError::MissingReferenceData {
            description: format!(
                "intronic {}. position requires reference provider with exon data",
                coord
            ),
        });
    }
    if pos.utr3 {
        return Err(ConversionError::MissingReferenceData {
            description: format!(
                "3' UTR {}. position (*N) requires reference provider with transcript length",
                coord
            ),
        });
    }
    if pos.base < 1 {
        return Err(ConversionError::MissingReferenceData {
            description: format!(
                "non-positive {}. position {} requires reference provider with transcript length",
                coord, pos.base
            ),
        });
    }
    Ok(pos.base as u64)
}

/// True if any endpoint of the interval needs provider data to resolve to a
/// transcript position.
fn tx_needs_provider(interval: &Interval<TxPos>) -> bool {
    let needs = |p: &TxPos| p.is_intronic() || p.is_downstream() || p.base < 1;
    interval.start.inner().is_some_and(needs) || interval.end.inner().is_some_and(needs)
}

fn rna_needs_provider(interval: &Interval<RnaPos>) -> bool {
    let needs = |p: &RnaPos| p.is_intronic() || p.utr3 || p.base < 1;
    interval.start.inner().is_some_and(needs) || interval.end.inner().is_some_and(needs)
}

/// Resolve a CDS-position pair to 1-based transcript positions using the
/// provider's transcript metadata.
///
/// Intronic c. positions (e.g. `c.100+5`) are rejected: SPDI is positional
/// and has no offset notation, so an intronic CDS variant cannot be expressed
/// on the transcript accession without first projecting to genomic coords.
/// That projection is intentionally out of scope for this entry point —
/// callers needing it can use the genomic conversion path explicitly.
fn resolve_cds_to_tx<P: ReferenceProvider + ?Sized>(
    accession: &Accession,
    start: &CdsPos,
    end: &CdsPos,
    provider: &P,
) -> Result<(u64, u64), ConversionError> {
    if start.is_intronic() || end.is_intronic() {
        return Err(ConversionError::MissingReferenceData {
            description: "intronic c. positions cannot be expressed in SPDI without genomic \
                          projection; SPDI is positional and has no offset notation"
                .to_string(),
        });
    }
    let tx_id = accession.transcript_accession();
    let transcript =
        provider
            .get_transcript(&tx_id)
            .map_err(|e| ConversionError::MissingReferenceData {
                description: format!("could not load transcript {}: {}", tx_id, e),
            })?;
    let mapper = CoordinateMapper::new(&transcript);
    let s = mapper
        .cds_to_tx(start)
        .map_err(|e| ConversionError::MissingReferenceData {
            description: format!("could not resolve {} to transcript position: {}", start, e),
        })?;
    let e = mapper
        .cds_to_tx(end)
        .map_err(|e| ConversionError::MissingReferenceData {
            description: format!("could not resolve {} to transcript position: {}", end, e),
        })?;
    // Bound both endpoints against the transcript length: reject a position that
    // maps past the 3' end, not just below base 1. #962 did this for c.*N (3'UTR);
    // #971 extends it to plain exonic c.N, which cds_to_tx can also map off-sequence
    // (e.g. c.99999 -> tx cds_start-1+99999). Keeps c.N / r.N / n.N symmetric.
    let tx_len = transcript.sequence_length();
    let s_u = ensure_tx_in_bounds(s.base, tx_len, "c", start)?;
    let e_u = ensure_tx_in_bounds(e.base, tx_len, "c", end)?;
    Ok((s_u, e_u))
}

/// Resolve an `n.` (TxPos) pair to 1-based transcript positions, including
/// intronic and downstream forms, using the provider.
fn resolve_tx_to_provider_tx<P: ReferenceProvider + ?Sized>(
    accession: &Accession,
    start: &TxPos,
    end: &TxPos,
    provider: &P,
) -> Result<(u64, u64), ConversionError> {
    let tx_id = accession.transcript_accession();
    let transcript =
        provider
            .get_transcript(&tx_id)
            .map_err(|e| ConversionError::MissingReferenceData {
                description: format!("could not load transcript {}: {}", tx_id, e),
            })?;
    let s = resolve_tx_pos(start, &transcript)?;
    let e = resolve_tx_pos(end, &transcript)?;
    Ok((s, e))
}

/// Best-effort upper-bounding of a plain exonic `n.` interval (positive base,
/// non-intronic, non-downstream). When the provider supplies the transcript,
/// each endpoint is bounded against its length via [`resolve_tx_pos`] (#971);
/// if the transcript cannot be loaded, fall back to the unbounded simple-path
/// value so a provider that lacks it still resolves (no regression).
///
/// Only a `MissingReferenceData` error (transcript unavailable) triggers the
/// fallback; an over-length position surfaces as `InvalidPosition` from the
/// bound and is propagated, not swallowed.
fn resolve_tx_exonic_bounded<P: ReferenceProvider + ?Sized>(
    accession: &Accession,
    interval: &Interval<TxPos>,
    provider: &P,
) -> Result<(u64, u64), ConversionError> {
    let start = interval
        .start
        .inner()
        .ok_or_else(|| ConversionError::InvalidPosition {
            description: "cannot convert n. variant with unknown start position".to_string(),
        })?;
    let end = resolve_transcript_end_boundary(interval, "n")?;
    match resolve_tx_to_provider_tx(accession, start, end, provider) {
        Ok(pair) => Ok(pair),
        Err(ConversionError::MissingReferenceData { .. }) => {
            let s = tx_pos_for_simple_path(interval, "n")?;
            let e = tx_end_for_simple_path(interval, "n")?;
            Ok((s, e))
        }
        Err(e) => Err(e),
    }
}

/// Best-effort upper-bounding of a plain exonic `r.` interval, mirroring
/// [`resolve_tx_exonic_bounded`]: bound against the transcript length when the
/// provider supplies it (#971), else fall back to the unbounded simple-path
/// value. Only `MissingReferenceData` (transcript unavailable) triggers the
/// fallback; an over-length position propagates as `InvalidPosition`.
fn resolve_rna_exonic_bounded<P: ReferenceProvider + ?Sized>(
    accession: &Accession,
    interval: &Interval<RnaPos>,
    provider: &P,
) -> Result<(u64, u64), ConversionError> {
    let start = interval
        .start
        .inner()
        .ok_or_else(|| ConversionError::InvalidPosition {
            description: "cannot convert r. variant with unknown start position".to_string(),
        })?;
    let end = resolve_transcript_end_boundary(interval, "r")?;
    match resolve_rna_to_provider_tx(accession, start, end, provider) {
        Ok(pair) => Ok(pair),
        Err(ConversionError::MissingReferenceData { .. }) => {
            let s = rna_pos_for_simple_path(interval, "r")?;
            let e = rna_end_for_simple_path(interval, "r")?;
            Ok((s, e))
        }
        Err(e) => Err(e),
    }
}

fn resolve_rna_to_provider_tx<P: ReferenceProvider + ?Sized>(
    accession: &Accession,
    start: &RnaPos,
    end: &RnaPos,
    provider: &P,
) -> Result<(u64, u64), ConversionError> {
    let tx_id = accession.transcript_accession();
    let transcript =
        provider
            .get_transcript(&tx_id)
            .map_err(|e| ConversionError::MissingReferenceData {
                description: format!("could not load transcript {}: {}", tx_id, e),
            })?;
    let s = resolve_rna_pos(start, &transcript)?;
    let e = resolve_rna_pos(end, &transcript)?;
    Ok((s, e))
}

/// Resolve a single `TxPos` to a 1-based transcript position. Intronic and
/// downstream (`n.*N`) positions are rejected because they have no valid
/// SPDI representation on the transcript accession; an exonic position past
/// the transcript's 3' end is also rejected (#971).
fn resolve_tx_pos(pos: &TxPos, transcript: &Transcript) -> Result<u64, ConversionError> {
    // SPDI is positional and has no offset notation, so intronic n. positions
    // cannot be expressed without genomic projection. Match the sibling
    // `resolve_rna_pos` (and the simple-path helpers) by emitting
    // `MissingReferenceData` here.
    if pos.is_intronic() {
        return Err(ConversionError::MissingReferenceData {
            description: format!(
                "intronic n.{} cannot be expressed in SPDI without genomic projection; \
                 SPDI is positional and has no offset notation",
                pos
            ),
        });
    }
    // n. has no CDS anchor, so `n.*N` is N bases past the transcript end —
    // an off-sequence position on the transcript accession. Reject until
    // genomic projection is wired in; emitting SPDI at `tx_len + N` would
    // produce a coordinate that does not exist on the accession.
    if pos.is_downstream() {
        return Err(ConversionError::InvalidPosition {
            description: format!(
                "downstream n.{} cannot be expressed in SPDI on the transcript accession \
                 without genomic projection",
                pos
            ),
        });
    }
    // Exonic n.N: bound against the transcript length so an over-length position
    // (n.99999 on a 40-base transcript) declines instead of emitting an
    // off-sequence SPDI coordinate (#971), mirroring the c.*N / r.*N bound (#962).
    ensure_tx_in_bounds(pos.base, transcript.sequence_length(), "n", pos)
}

fn resolve_rna_pos(pos: &RnaPos, transcript: &Transcript) -> Result<u64, ConversionError> {
    if pos.is_intronic() {
        return Err(ConversionError::MissingReferenceData {
            description: format!(
                "intronic r.{} cannot be expressed in SPDI without genomic projection; \
                 SPDI is positional and has no offset notation",
                pos
            ),
        });
    }
    if pos.utr3 {
        // r.*N is the Nth base of the 3' UTR — anchored after the stop codon,
        // so r.*1 sits at the base immediately after the last CDS position, not
        // at `tx_len + 1` (one past the entire mRNA). For non-coding transcripts
        // r.*N has no meaning; surface `MissingReferenceData` rather than fall
        // back to `tx_len` (#390 item 2).
        if pos.base < 1 {
            return Err(ConversionError::InvalidPosition {
                description: format!("3' UTR position *{} must be >= 1", pos.base),
            });
        }
        // r.*N and c.*N denote the same transcript position, so resolve r.*N
        // through the same `CoordinateMapper::cds_to_tx` the c.*N path
        // (`resolve_cds_to_tx`) uses, rather than open-coding `cds_end + base`
        // here (#390-follow-up / #944). One conversion, one place: the two
        // arms cannot drift apart, which is the property those issues are
        // about. Since #1619 that conversion is the flat sequence-axis shift —
        // `c.*N` counts the transcript's own bases, including any the exon
        // table does not cover — so the two spellings agree by construction.
        // Requires a CDS end — non-coding transcripts have no 3'UTR anchor.
        if transcript.cds_end.is_none() {
            return Err(ConversionError::MissingReferenceData {
                description: format!(
                    "r.*{} requires a CDS end on the transcript; non-coding \
                     transcripts have no 3'UTR anchor",
                    pos.base
                ),
            });
        }
        let mapper = CoordinateMapper::new(transcript);
        let tx = mapper.cds_to_tx(&CdsPos::utr3(pos.base)).map_err(|e| {
            ConversionError::MissingReferenceData {
                description: format!(
                    "could not resolve r.*{} to transcript position: {}",
                    pos.base, e
                ),
            }
        })?;
        return ensure_tx_in_bounds(tx.base, transcript.sequence_length(), "r", pos);
    }
    if pos.is_5utr() {
        // r.-N is the Nth base of the 5' UTR — numbered upstream from the start
        // codon, exactly mirroring c.-N (numbering.md:58/61: RNA numbering
        // follows the coding DNA reference, which includes c.-N). Resolve it
        // through the same `CoordinateMapper::cds_to_tx` the c.-N path
        // (`resolve_cds_to_tx`) uses — via the equivalent `CdsPos` with base < 1 —
        // so r.-N and c.-N always land at the same transcript/SPDI position.
        // Routing both through one conversion is what guarantees that; since
        // #1619 the conversion itself is the flat `cds_start - N` shift, on the
        // transcript's own bases rather than on the exon table. Requires a CDS
        // start — non-coding transcripts have no 5'UTR anchor.
        if transcript.cds_start.is_none() {
            return Err(ConversionError::MissingReferenceData {
                description: format!(
                    "r.{} requires a CDS start on the transcript; non-coding \
                     transcripts have no 5'UTR anchor",
                    pos.base
                ),
            });
        }
        let mapper = CoordinateMapper::new(transcript);
        let tx = mapper.cds_to_tx(&CdsPos::new(pos.base)).map_err(|e| {
            ConversionError::MissingReferenceData {
                description: format!(
                    "could not resolve r.{} to transcript position: {}",
                    pos.base, e
                ),
            }
        })?;
        return ensure_positive_tx(tx.base, "r", pos);
    }
    // Exonic r.N. On a CODING transcript r. numbering is CDS-relative (#469):
    // r.N denotes the same base as c.N (see `cds_pos_to_rna` in project/rna.rs),
    // so resolve it through the same `CoordinateMapper::cds_to_tx` exactly
    // like c.N and like the r.-N / r.*N branches above — otherwise r.N would be
    // treated as transcript-absolute and disagree with c.N by (cds_start - 1).
    // A NON-coding (NR_) transcript has no CDS anchor, so r.N IS the transcript
    // position directly. Either way, bound against the transcript length (#971).
    if transcript.cds_start.is_some() {
        // Resolving a coding r.N CDS-relative needs a CDS end (`cds_to_tx`
        // requires both anchors). A coding transcript missing its
        // CDS end is malformed/partial annotation: decline with
        // `InvalidPosition` (which propagates past the
        // `resolve_rna_exonic_bounded` `MissingReferenceData` fallback) rather
        // than silently resolving r.N transcript-absolute and unbounded — the
        // same "loaded but inexpressible" treatment `resolve_tx_pos` gives a
        // downstream n.*N. Mirrors the cds_end guard on the r.*N branch above.
        if transcript.cds_end.is_none() {
            return Err(ConversionError::InvalidPosition {
                description: format!(
                    "r.{} on a coding transcript requires a CDS end to resolve \
                     CDS-relative; the transcript has a CDS start but no CDS end",
                    pos.base
                ),
            });
        }
        let mapper = CoordinateMapper::new(transcript);
        let tx = mapper.cds_to_tx(&CdsPos::new(pos.base)).map_err(|e| {
            ConversionError::MissingReferenceData {
                description: format!(
                    "could not resolve r.{} to transcript position: {}",
                    pos.base, e
                ),
            }
        })?;
        return ensure_tx_in_bounds(tx.base, transcript.sequence_length(), "r", pos);
    }
    ensure_tx_in_bounds(pos.base, transcript.sequence_length(), "r", pos)
}

fn ensure_positive_tx<P: std::fmt::Display>(
    base: i64,
    coord: &str,
    pos: P,
) -> Result<u64, ConversionError> {
    if base < 1 {
        return Err(ConversionError::InvalidPosition {
            description: format!(
                "transcript position from {}. coordinate {} resolves to a non-positive base ({})",
                coord, pos, base
            ),
        });
    }
    Ok(base as u64)
}

/// Like [`ensure_positive_tx`], but also rejects a resolved 1-based transcript
/// position that falls *past the transcript's 3' end* (`base > tx_len`).
///
/// A `*N` position numbered beyond the last transcript base — e.g. `r.*99999` /
/// `c.*99999` on a short transcript — maps to a coordinate that does not exist on
/// the accession. `ensure_positive_tx` guards only the lower bound, so without
/// this the resolver would emit an off-sequence SPDI position (#962). `tx_len` is
/// the transcript length in mRNA bases ([`Transcript::sequence_length`]). Applied
/// at the `*N` (3'UTR) resolution sites — `r.*N` in `resolve_rna_pos` and `c.*N`
/// in `resolve_cds_to_tx` — where a position can be numbered past the last
/// transcript base, and at the exonic resolution sites — `n.N` in
/// `resolve_tx_pos`, the exonic tail of `r.N` in `resolve_rna_pos`, and both
/// arms of `c.N` in `resolve_cds_to_tx` — where an over-length position (e.g.
/// `c.99999` on a short transcript) can likewise map past the 3' end (#971).
fn ensure_tx_in_bounds<P: std::fmt::Display + Copy>(
    base: i64,
    tx_len: u64,
    coord: &str,
    pos: P,
) -> Result<u64, ConversionError> {
    let one_based = ensure_positive_tx(base, coord, pos)?;
    if one_based > tx_len {
        return Err(ConversionError::InvalidPosition {
            description: format!(
                "transcript position from {}. coordinate {} resolves to base {} past the \
                 transcript 3' end (length {})",
                coord, pos, one_based, tx_len
            ),
        });
    }
    Ok(one_based)
}

/// Apply edit-specific position arithmetic and emit the SPDI variant.
///
/// `start_one_based` and `end_one_based` are 1-based positions on the SPDI
/// accession (genomic for `g.`/`m.`, transcript for `c.`/`n.`/`r.`).
///
/// # Position convention per edit
///
/// SPDI's `position` field is 0-based. For substitution / deletion /
/// delins (edits that reference specific bases) the SPDI position
/// equals `start_one_based - 1` (computed once at the top of the
/// function as `spdi_pos`). For **insertion** and **duplication** the
/// SPDI position is *interbase* — position N is the boundary AFTER
/// 1-based base N — so the Insertion arm uses `start_one_based`
/// directly and the Duplication arm uses `end_one_based` directly
/// (#390 item 1). The top-of-function `spdi_pos` is unused in those
/// two arms.
///
/// # Provider behavior
///
/// When `provider` is `Some` and the edit is a short-form deletion,
/// duplication, delins, or identity (i.e. lacks an explicit deleted
/// sequence), the reference bases for `[start_one_based, end_one_based]` are
/// fetched via the provider so SPDI's mandatory `del` field can be populated.
/// When `provider` is `None`, those cases return
/// [`ConversionError::MissingReferenceData`].
///
/// A **whole-entity** identity (`g.=`) is the one exception: it names no
/// interval, so it returns [`ConversionError::UnsupportedEditType`] rather
/// than a triple at an arbitrary position.
fn emit_spdi_for_edit<P>(
    sequence: String,
    start_one_based: u64,
    end_one_based: u64,
    edit: &NaEdit,
    alphabet: AlphabetMode,
    provider: Option<&P>,
) -> Result<SpdiVariant, ConversionError>
where
    P: ReferenceProvider + ?Sized,
{
    let hgvs_pos_ob =
        OneBasedPos::try_new(start_one_based).ok_or_else(|| ConversionError::InvalidPosition {
            description: "position 0 is not valid in HGVS".to_string(),
        })?;
    let spdi_pos_zb: ZeroBasedPos = hgvs_pos_ob.to_zero_based();
    let spdi_pos = spdi_pos_zb.value();

    match edit {
        NaEdit::Substitution {
            reference,
            alternative,
        } => Ok(SpdiVariant::new(
            sequence,
            spdi_pos,
            apply_alphabet(&reference.to_string(), alphabet),
            apply_alphabet(&alternative.to_string(), alphabet),
        )),
        NaEdit::Insertion { sequence: inserted } => {
            // An inserted sequence reaches SPDI when its exact bases are known.
            // Three shapes are handled here: a literal; a same-reference range
            // insert (`ins50_57`, `ins50_57inv`), resolved from the reference
            // exactly as `del`/`dup` resolve their omitted bases; and an exact
            // tandem-repeat insert (`insA[10]`, `insAT[3]`), expanded from its
            // spelled unit exactly as the short-form `Repeat` arm expands its
            // own. The fallback rejects everything else, which is now uniformly
            // shapes whose bases are genuinely undetermined: `Count` (`ins10`),
            // `Range` (`ins(10_20)`), `Uncertain` (`ins(?)`), an uncertain or
            // range repeat count (`insA[10_15]`), and named/external references
            // SPDI cannot dereference.
            let ins_str = if let Some(resolved) =
                resolve_position_range_insert(inserted, &sequence, alphabet, provider)?
            {
                resolved
            } else if let Some(resolved) = resolve_exact_repeat_insert(inserted, alphabet)? {
                resolved
            } else if let Some(resolved) =
                resolve_complex_insert(inserted, &sequence, alphabet, provider)?
            {
                resolved
            } else {
                let literal = inserted_sequence_to_string(inserted).ok_or_else(|| {
                    ConversionError::MissingReferenceData {
                        description: "insertion sequence is neither a literal sequence, a \
                                      same-reference position range, nor an exact tandem repeat; \
                                      this shape is not yet encodable as SPDI"
                            .to_string(),
                    }
                })?;
                apply_alphabet(&literal, alphabet)
            };
            // SPDI is 0-based interbase: position N is the boundary
            // AFTER 1-based base N (equivalently between 1-based bases
            // N and N+1). For HGVS `g.{start}_{start+1}ins{seq}` the
            // matching SPDI position is `start_one_based` directly
            // (NOT the -1 conversion used for substitution/deletion,
            // which references a specific base). Closes #390 item 1.
            Ok(SpdiVariant::new(sequence, start_one_based, "", ins_str))
        }
        NaEdit::Duplication {
            sequence: dup_seq, ..
        } => {
            let dup_str = match dup_seq {
                Some(seq) => sequence_to_string(seq),
                None => match provider {
                    Some(p) => fetch_reference_bases(p, &sequence, start_one_based, end_one_based)?,
                    None => {
                        return Err(unspelled_bases_error(
                            UnspelledBases::Duplicated,
                            start_one_based,
                            end_one_based,
                        ));
                    }
                },
            };
            // SPDI encodes a duplication as an insertion immediately
            // after the duplicated region. With the corrected
            // interbase convention (see Insertion arm above), that
            // position is the 1-based end of the duplicated region
            // (`end_one_based`) directly. Closes #390 item 1.
            Ok(SpdiVariant::new(
                sequence,
                end_one_based,
                "",
                apply_alphabet(&dup_str, alphabet),
            ))
        }
        NaEdit::Deletion {
            sequence: del_seq, ..
        } => {
            let del_str = match del_seq {
                Some(seq) => sequence_to_string(seq),
                None => match provider {
                    Some(p) => fetch_reference_bases(p, &sequence, start_one_based, end_one_based)?,
                    None => {
                        return Err(unspelled_bases_error(
                            UnspelledBases::Deleted,
                            start_one_based,
                            end_one_based,
                        ));
                    }
                },
            };
            Ok(SpdiVariant::new(
                sequence,
                spdi_pos,
                apply_alphabet(&del_str, alphabet),
                "",
            ))
        }
        NaEdit::Delins {
            sequence: ins_seq,
            deleted,
            deleted_length: _,
            substitution_reference: None,
        } => {
            // Closes #394 item 3. An inserted sequence reaches SPDI when its
            // exact bases are known. Three shapes are handled: a literal; a
            // same-reference position-range insert (`delins50_57`,
            // `delins50_57inv`), resolved from the reference exactly as
            // `del`/`dup` resolve their omitted bases; and an exact tandem-repeat
            // insert (`delinsA[10]`, `delinsAT[3]`), expanded from its spelled
            // unit exactly as the short-form `Repeat` arm expands its own. The
            // fallback rejects everything else, which is now uniformly shapes
            // whose bases are genuinely undetermined (`Count` `ins10`, `Range`
            // `ins(10_20)`, `Uncertain`, an uncertain or range repeat count, and
            // named/external references SPDI cannot dereference), and declines as
            // `UnsupportedEditType` — matching the sibling arms.
            let ins_str = if let Some(resolved) =
                resolve_position_range_insert(ins_seq, &sequence, alphabet, provider)?
            {
                resolved
            } else if let Some(resolved) = resolve_exact_repeat_insert(ins_seq, alphabet)? {
                resolved
            } else if let Some(resolved) =
                resolve_complex_insert(ins_seq, &sequence, alphabet, provider)?
            {
                resolved
            } else {
                let literal = inserted_sequence_to_string(ins_seq).ok_or_else(|| {
                    ConversionError::UnsupportedEditType {
                        description: "delins inserted sequence is neither a literal sequence, a \
                                      same-reference position range, nor an exact tandem repeat; \
                                      this shape is not yet encodable as SPDI"
                            .to_string(),
                    }
                })?;
                apply_alphabet(&literal, alphabet)
            };
            let del_str = match deleted {
                Some(seq) => sequence_to_string(seq),
                None => match provider {
                    Some(p) => fetch_reference_bases(p, &sequence, start_one_based, end_one_based)?,
                    None => {
                        return Err(unspelled_bases_error(
                            UnspelledBases::DeletedInDelins,
                            start_one_based,
                            end_one_based,
                        ));
                    }
                },
            };
            Ok(SpdiVariant::new(
                sequence,
                spdi_pos,
                apply_alphabet(&del_str, alphabet),
                ins_str,
            ))
        }
        NaEdit::Identity {
            sequence: id_seq,
            whole_entity,
        } => {
            // An identity claims the bases it names, so its triple must span
            // them: `g.263=` is `262:A:A`, not the zero-width `262::`. A
            // zero-width triple is not merely less informative, it *aliases an
            // insertion junction* — put `g.[261_262dup;263=]` through it and
            // both members land on interbase 262, which is indistinguishable
            // from two insertions competing for one interbase.
            //
            // When the input does not spell the base, fetch it, exactly as the
            // Deletion / Delins / Inversion arms do for their own omitted
            // sequences.
            let ref_base = match id_seq {
                Some(seq) => sequence_to_string(seq),
                // `g.=` asserts the whole reference is unchanged and names no
                // interval, so there is no position SPDI could honestly carry.
                // Emitting `0::` looks harmless but round-trips through
                // `spdi_to_hgvs` as `g.1=`, narrowing a statement about the
                // whole sequence to one about base 1. Decline instead, as this
                // module does for every other shape SPDI cannot encode.
                None if *whole_entity => {
                    return Err(ConversionError::UnsupportedEditType {
                        description: "a whole-entity identity (`=`) asserts the entire \
                                      reference is unchanged and names no interval; SPDI \
                                      has no representation for it"
                            .to_string(),
                    });
                }
                None => match provider {
                    Some(p) => fetch_reference_bases(p, &sequence, start_one_based, end_one_based)?,
                    None => {
                        return Err(unspelled_bases_error(
                            UnspelledBases::Unchanged,
                            start_one_based,
                            end_one_based,
                        ));
                    }
                },
            };
            let ref_base = apply_alphabet(&ref_base, alphabet);
            Ok(SpdiVariant::new(
                sequence,
                spdi_pos,
                ref_base.clone(),
                ref_base,
            ))
        }
        NaEdit::Inversion {
            sequence: inv_seq, ..
        } => {
            // SPDI has no native inv; the standard mapping is delins where
            // ins is the reverse-complement of the deleted reference span.
            let del_raw = match inv_seq {
                Some(seq) => sequence_to_string(seq),
                None => match provider {
                    Some(p) => fetch_reference_bases(p, &sequence, start_one_based, end_one_based)?,
                    None => {
                        return Err(unspelled_bases_error(
                            UnspelledBases::Inverted,
                            start_one_based,
                            end_one_based,
                        ));
                    }
                },
            };
            let del_str = apply_alphabet(&del_raw, alphabet);
            let ins_str = reverse_complement(&del_str);
            Ok(SpdiVariant::new(sequence, spdi_pos, del_str, ins_str))
        }
        NaEdit::Repeat {
            sequence: unit_seq,
            count,
            additional_counts,
            trailing,
        } => {
            // SPDI has no native repeat; expand to delins where del is the
            // reference repeat tract and ins is the unit repeated `count`
            // times. Only an explicit unit and exact count are representable.
            if !additional_counts.is_empty() {
                return Err(ConversionError::UnsupportedEditType {
                    description:
                        "genotype-style repeat (multiple counts) cannot be expressed as a single SPDI; emit each allele separately"
                            .to_string(),
                });
            }
            if trailing.is_some() {
                return Err(ConversionError::UnsupportedEditType {
                    description: "repeat with trailing sequence cannot be represented in SPDI"
                        .to_string(),
                });
            }
            let unit = unit_seq
                .as_ref()
                .ok_or_else(|| ConversionError::MissingReferenceData {
                    description:
                        "repeat unit sequence not provided; cannot expand into SPDI delins"
                            .to_string(),
                })?;
            let n_post = match count {
                RepeatCount::Exact(n) => *n as usize,
                // Range, UncertainRange, MinUncertain, MaxUncertain, Unknown:
                // no single expanded ins-string exists, so SPDI cannot encode.
                _ => {
                    return Err(ConversionError::UnsupportedEditType {
                        description:
                            "uncertain or range repeat counts cannot be represented in SPDI"
                                .to_string(),
                    });
                }
            };
            let unit_str = apply_alphabet(&sequence_to_string(unit), alphabet);
            // Bound the expanded ins-string before allocating. The count is
            // user-controlled (u64), so `unit_str.repeat(n_post)` can be
            // forced to allocate gigabytes from a single short input.
            let expansion_bases = unit_str.len().checked_mul(n_post).ok_or_else(|| {
                ConversionError::UnsupportedEditType {
                    description: format!(
                        "repeat expansion {} x {} overflows usize",
                        unit_str.len(),
                        n_post
                    ),
                }
            })?;
            if expansion_bases > MAX_REPEAT_EXPANSION_BASES {
                return Err(ConversionError::UnsupportedEditType {
                    description: format!(
                        "repeat expansion {} bases exceeds SPDI ins-string cap of {} bases",
                        expansion_bases, MAX_REPEAT_EXPANSION_BASES
                    ),
                });
            }
            // A single-position anchor names the tract's *start*, not a
            // one-base span, so widen it to the physical run before fetching
            // (#1431).
            //
            // The spec presents the two spellings as two formats of one
            // variant — "a Community Consultation proposal is being prepared
            // which will suggest to allow only the format where the **entire
            // range** of the repeated sequence is indicated; so
            // `g.123_191CAG[23]`, **not** `g.123CAG[23]`"
            // (`DNA/repeated.md`) — and `123_191` is 69 bases, the whole
            // 23-copy tract. Reading the anchor as a one-base span instead
            // made the two disagree: on a 3-base `A` tract at 263,
            // `g.263A[7]` came out as `262:A:AAAAAAA`, which denotes **nine**
            // `A`s once the two untouched tract bases are counted, while
            // `g.263_265A[7]` correctly denotes seven.
            //
            // The one-base reading also does not survive a multi-base unit at
            // all: `g.263CAG[5]` fetched a single base and died on the
            // divisibility check below ("length 1 is not a multiple of unit
            // length 3"), so the spec's own start-only format was unusable for
            // every unit the spec actually illustrates it with.
            //
            // This is the reading `merge::lowered_repeat` already documents —
            // "a single-position anchor (`e == s`) names only the tract's
            // start and means 'the whole run becomes N', so it absorbs
            // nothing" — so widening here makes the two sites agree rather
            // than introducing a third answer.
            let (del_start, del_end, span_origin) = match provider {
                Some(p) if start_one_based == end_one_based => resolve_repeat_tract_span(
                    p,
                    &sequence,
                    start_one_based,
                    unit_str.as_bytes(),
                    alphabet,
                )?,
                // An explicit range is the caller's own span, so a decline may
                // quote it verbatim.
                _ => (start_one_based, end_one_based, RepeatSpanOrigin::Tract),
            };
            let spdi_pos = del_start - 1;
            let del_str = match provider {
                Some(p) => {
                    fetch_normalized_reference_bases(p, &sequence, del_start, del_end, alphabet)?
                }
                None => {
                    return Err(unspelled_bases_error(
                        UnspelledBases::RepeatTract,
                        start_one_based,
                        end_one_based,
                    ));
                }
            };
            // Two things must hold for the span to be the declared repeat: the
            // HGVS recommendations require an integer number of units
            // (`inversion`-style nonsense otherwise), and a divisible length is
            // necessary but not sufficient — `ATGCAT` is length 6 and is not
            // `AT[3]`.
            //
            // **Both are judged before either is reported**, because the two
            // failures share one diagnostic decision. After a failed tract
            // search `del_start`/`del_end` are this module's unit-wide fallback
            // rather than anything the caller wrote, and *either* check can be
            // the one that trips on it: near the 3' end the fallback is clamped
            // inside the contig, so a multi-base unit leaves a short span that
            // fails divisibility first. Reporting that span names coordinates
            // nobody supplied — the exact defect #1492 is about — so the origin
            // decides the message and the checks only decide *whether* there is
            // one.
            let divisible = !unit_str.is_empty() && del_str.len().is_multiple_of(unit_str.len());
            let matches_unit =
                divisible && del_str == unit_str.repeat(del_str.len() / unit_str.len());
            if !matches_unit {
                return Err(ConversionError::InvalidPosition {
                    description: match span_origin {
                        // The span is the fallback, not the caller's. Name the
                        // anchor and the real fault: no run begins there.
                        RepeatSpanOrigin::NoRunAtAnchor => format!(
                            "no {unit} repeat is anchored at {seq}:{anchor}: the reference \
                             there is not a run of {unit}. A single-position anchor names the \
                             *first base* of the tract, so check the run's start, or spell the \
                             whole range explicitly (`g.<start>_<end>{unit}[n]`).",
                            unit = unit_str,
                            seq = sequence,
                            anchor = start_one_based,
                        ),
                        // The span is one the caller wrote or the search chose,
                        // so quoting it is honest.
                        _ if !divisible => format!(
                            "repeat span {}:{}-{} length {} is not a multiple of unit length {}",
                            sequence,
                            del_start,
                            del_end,
                            del_str.len(),
                            unit_str.len()
                        ),
                        _ => format!(
                            "repeat span {}:{}-{} does not match repeat unit {}",
                            sequence, del_start, del_end, unit_str
                        ),
                    },
                });
            }
            let ins_str = unit_str.repeat(n_post);
            Ok(SpdiVariant::new(sequence, spdi_pos, del_str, ins_str))
        }
        NaEdit::CopyNumber { .. } => Err(ConversionError::UnsupportedEditType {
            description: "copy number variants cannot be represented in SPDI format".to_string(),
        }),
        NaEdit::Conversion { .. } => Err(ConversionError::UnsupportedEditType {
            description: "conversion variants cannot be represented in SPDI format".to_string(),
        }),
        _ => Err(ConversionError::UnsupportedEditType {
            description: format!("unsupported edit type: {:?}", edit),
        }),
    }
}

/// The edit shapes whose SPDI form needs reference bases the description did
/// not spell out. Each names what is unknown and how the caller could spell it.
///
/// Exists so the six arms of [`hgvs_to_spdi`] that hit this wall decline in one
/// voice. They previously carried six hand-written strings, none of which named
/// the interval that could not be resolved or either way out, and one of which
/// (`Delins`) had drifted into a different sentence shape than its siblings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum UnspelledBases {
    Duplicated,
    Deleted,
    DeletedInDelins,
    Unchanged,
    Inverted,
    /// A repeat's pre-expansion tract. The odd one out: see
    /// [`Self::spelled_example`].
    RepeatTract,
}

impl UnspelledBases {
    /// How to describe the missing bases: `"the {} bases at 10..=12"`.
    fn adjective(self) -> &'static str {
        match self {
            Self::Duplicated => "duplicated",
            Self::Deleted | Self::DeletedInDelins => "deleted",
            Self::Unchanged => "unchanged",
            Self::Inverted => "inverted",
            Self::RepeatTract => "pre-expansion repeat-tract",
        }
    }

    /// A spelled-out form of this shape, offered as the example in the message —
    /// or `None` where the description has no way to carry the bases.
    ///
    /// `RepeatTract` is that case, and it is why this returns an `Option`. The
    /// other five arms all read an optional sequence off the edit first and only
    /// consult the provider when it is absent, so "spell it" is genuinely a way
    /// out. The repeat arm has no such field to fill: the notation states the
    /// *unit* (`g.10AC[3]`), never the span the tract currently occupies, so the
    /// provider is the only route and suggesting otherwise would send the caller
    /// after a description they cannot write.
    ///
    /// Every `Some` here **must parse** — an error message that suggests a
    /// description ferro would itself reject is worse than one that suggests
    /// nothing. `spelled_examples_are_parseable` holds them to that.
    fn spelled_example(self) -> Option<&'static str> {
        match self {
            Self::Duplicated => Some("g.10_12dupACG"),
            Self::Deleted => Some("g.10_12delACG"),
            Self::DeletedInDelins => Some("g.10_12delACGinsT"),
            Self::Unchanged => Some("g.10A="),
            Self::Inverted => Some("g.10_12invACG"),
            Self::RepeatTract => None,
        }
    }
}

/// Decline a conversion that needs bases the description left implicit.
///
/// Names the interval and every remedy that actually exists, because a caller who
/// hits this has one or two ways forward — spell the bases in the description
/// where the notation allows it, or pass a provider that can look them up — and
/// neither is guessable from "reference data needed".
fn unspelled_bases_error(
    what: UnspelledBases,
    start_one_based: u64,
    end_one_based: u64,
) -> ConversionError {
    let remedy = match what.spelled_example() {
        Some(example) => format!(
            "Spell them in the description (e.g. `{example}`) or convert with a reference provider."
        ),
        None => "Convert with a reference provider: the repeat notation states the unit, \
                 not the span its tract currently occupies, so the bases cannot be spelled \
                 in the description."
            .to_string(),
    };
    ConversionError::MissingReferenceData {
        description: format!(
            "cannot convert to SPDI: the {} bases at {}..={} are unknown (no reference data). \
             {remedy}",
            what.adjective(),
            start_one_based,
            end_one_based,
        ),
    }
}

/// Why the span [`resolve_repeat_tract_span`] returned is the span it returned.
///
/// Exists so a decline can name the caller's own coordinates (#1492). When the
/// search finds no run, the function falls back to the unit-wide span *at* the
/// anchor — which is not something the caller wrote. Reporting that span in an
/// error told the reader to go look at a window they never named: `g.260CAG[5]`
/// declined with "repeat span 260-262 does not match repeat unit CAG", where
/// `260-262` is this function's invention and the caller's actual mistake — that
/// the run begins at 259, not 260 — went unmentioned. 336 of the 560 rows in
/// #1452's census took that path, so it is the common decline, not a corner.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RepeatSpanOrigin {
    /// A run was found; the span is its full extent. A unit-match failure here
    /// is a real mismatch over coordinates the search itself chose.
    Tract,
    /// The search ran and found no run at the anchor, so the span is the
    /// unit-wide fallback. This is the case that must not quote its span.
    NoRunAtAnchor,
    /// The search never ran — an empty unit, or a provider that serves bases
    /// but not `get_sequence_length`. The span is a placeholder rather than a
    /// judgement, so the generic diagnostic remains the honest one: nothing was
    /// measured that would justify saying "no run begins here".
    NotSearched,
}

/// The physical repeat run a single-position repeat anchor names, as a 1-based
/// inclusive span (#1431).
///
/// A start-only repeat (`g.263A[7]`, `g.123CAG[23]`) addresses the whole tandem
/// run it points into, not the one base or one unit at the anchor — see the
/// call site for the spec text and the two ways the one-base reading went
/// wrong. This finds that run with the same routine the normalizer uses
/// ([`crate::normalize::rules::count_tandem_repeats`]) so the two cannot drift
/// apart.
///
/// **Window growth, and why the cap is not a silent truncation.** The tract's
/// extent is not known before it is read, so the reference is fetched in a
/// window around the anchor and the window is doubled while the run still
/// reaches either edge. A tract that is still growing at
/// [`MAX_REPEAT_SEARCH_BASES`] is *declined* rather than reported at the
/// window's width: returning the clamped span would silently under-count the
/// run and emit a triple denoting the wrong bases, which is the exact defect
/// this function exists to remove.
///
/// Falls back to the unit-wide span at the anchor when no run is found — the
/// unit does not match the reference there — so the caller's existing
/// "does not match repeat unit" diagnostic is what reports it, rather than a
/// second error for one condition.
///
/// **The window is normalised before it is searched** (#1452). `unit` arrives
/// having already been through [`apply_alphabet`], so searching the raw window
/// would compare a soft-masked (lowercase) reference against an uppercase unit
/// and find no run at all. That did not surface as an error: the no-run
/// fallback above returns the *unit-wide* span, whose bases the caller then
/// uppercases, so the unit-match check passed and a truncated triple was
/// emitted. On a 3-copy lowercase `cag` tract, `g.259CAG[5]` came out as
/// `258:CAG:CAGCAGCAGCAGCAG` — seven copies once the two untouched tract
/// units are counted — while the range spelling `g.259_267CAG[5]` correctly
/// gave `258:CAGCAGCAG:CAGCAGCAGCAGCAG`. Both spellings converted; they just
/// denoted different sequences.
///
/// Normalising with [`apply_alphabet`] rather than a bare `to_ascii_uppercase`
/// is what keeps the two comparisons provably identical: on the `r.` axis the
/// unit has had `U` rewritten to `T`, so an uppercase-only window would still
/// miss a `U`-spelled tract and truncate it by the same route.
fn resolve_repeat_tract_span<P>(
    provider: &P,
    accession: &str,
    anchor_one_based: u64,
    unit: &[u8],
    alphabet: AlphabetMode,
) -> Result<(u64, u64, RepeatSpanOrigin), ConversionError>
where
    P: ReferenceProvider + ?Sized,
{
    /// Half-width of the first window. Comfortably covers the tract lengths the
    /// spec illustrates (its `CAG[23]` example is 69 bases) so the common case
    /// costs one fetch.
    const INITIAL_HALF_WIDTH: u64 = 128;

    if unit.is_empty() {
        return Ok((
            anchor_one_based,
            anchor_one_based,
            RepeatSpanOrigin::NotSearched,
        ));
    }
    let unit_len = unit.len() as u64;
    let fallback = (
        anchor_one_based,
        anchor_one_based.saturating_add(unit_len - 1),
    );

    // A failed length lookup must not fail the conversion. `get_sequence_length`
    // carries a trait default that always errors (`provider.rs:424`), so a
    // provider that serves bases perfectly well but does not implement it would
    // otherwise go from converting `g.263A[7]` to refusing it outright — a
    // regression for every out-of-tree `ReferenceProvider`, and one no in-repo
    // test can see because all of ferro's own providers override it.
    //
    // Without a length the window cannot be clamped or its edges judged, so the
    // tract search is not attempted; the unit-wide span is what this returns
    // instead. That leaves such a provider with the pre-#1431 answer rather than
    // the corrected one, which is a known and stated limit, not a silent wrong
    // span — the caller's divisibility and unit-match checks still judge it.
    let Ok(sequence_length) = provider.get_sequence_length(accession) else {
        return Ok((fallback.0, fallback.1, RepeatSpanOrigin::NotSearched));
    };
    if sequence_length == 0 || anchor_one_based > sequence_length {
        return Ok((fallback.0, fallback.1, RepeatSpanOrigin::NotSearched));
    }
    // Now that the length is known, keep the fallback inside the contig: a
    // multi-base unit anchored within `unit_len` of the 3' end would otherwise
    // hand `fetch_reference_bases` a range past the end, replacing the
    // documented "does not match repeat unit" diagnostic with a fetch failure.
    let fallback = (fallback.0, fallback.1.min(sequence_length));

    let mut half_width = INITIAL_HALF_WIDTH;
    loop {
        let window_start = anchor_one_based.saturating_sub(half_width).max(1);
        let window_end = anchor_one_based
            .saturating_add(half_width)
            .min(sequence_length);
        let window = fetch_normalized_reference_bases(
            provider,
            accession,
            window_start,
            window_end,
            alphabet,
        )?;
        let bytes = window.as_bytes();
        let anchor_offset = (anchor_one_based - window_start) as usize;

        let Some((_, tract_start, tract_end)) =
            crate::normalize::rules::count_tandem_repeats(bytes, anchor_offset, unit)
        else {
            return Ok((fallback.0, fallback.1, RepeatSpanOrigin::NoRunAtAnchor));
        };

        // Only grow while the run is still touching an edge the contig has not
        // itself ended at — otherwise the window is not what is bounding it.
        //
        // Both tests are stated in units, not bytes, because
        // `count_tandem_repeats` steps by `unit_len`: a run clipped by the
        // window stops up to `unit_len - 1` bytes short of the edge rather than
        // on it. Testing `tract_end == bytes.len()` therefore only fires when
        // the remaining byte count happens to be a multiple of the unit — true
        // for a 1- or 3-base unit at the initial half-width, false for a 4- or
        // 5-base one — and the clipped span was returned as if it were the whole
        // run. Same on the 5' side: the backward scan stops at
        // `anchor_offset % unit_len`, which need not be 0.
        let open_at_start = tract_start < unit.len() && window_start > 1;
        let open_at_end = tract_end + unit.len() > bytes.len() && window_end < sequence_length;
        if !open_at_start && !open_at_end {
            return Ok((
                window_start + tract_start as u64,
                window_start + tract_end as u64 - 1,
                RepeatSpanOrigin::Tract,
            ));
        }
        // Judged on the whole window, not the half-width, so the figure the
        // message quotes is the window that was actually read.
        if half_width.saturating_mul(2).saturating_add(1) >= MAX_REPEAT_SEARCH_BASES {
            return Err(ConversionError::UnsupportedEditType {
                description: format!(
                    "repeat tract at {accession}:{anchor_one_based} still extends past a \
                     {MAX_REPEAT_SEARCH_BASES}-base search window; spell the whole range \
                     explicitly (`g.<start>_<end>{}[n]`) rather than the start alone",
                    String::from_utf8_lossy(unit)
                ),
            });
        }
        half_width = half_width.saturating_mul(2);
    }
}

/// [`fetch_reference_bases`] with [`apply_alphabet`] already applied, so the
/// caller never holds a window in the reference's own case convention (#1452).
///
/// Reference FASTAs are routinely soft-masked, and the repeat arm compares the
/// fetched bases against a unit that has been through `apply_alphabet` — twice,
/// at two sites. Normalising at the fetch is what stops those two sites from
/// disagreeing: the tract search and the unit-match check now see the same
/// bytes by construction rather than by each remembering to fold. This mirrors
/// how `normalize::merge::canonical_base_byte` centralises the `r.` uracil /
/// thymine equivalence instead of spreading it over its comparison sites.
///
/// The byte offsets a caller computes against the returned string stay valid:
/// `apply_alphabet` rewrites ASCII bytes one-for-one and leaves everything else
/// alone, so the length is preserved. `normalizing_a_window_preserves_its_length`
/// pins that, since a length change would silently mis-map the anchor offset in
/// [`resolve_repeat_tract_span`] rather than fail.
///
/// **Only the repeat arm uses this, and the other arms are not an oversight.**
/// The Duplication / Deletion / Delins / Identity / Inversion arms each fetch
/// through the plain [`fetch_reference_bases`] and fold *after* the `match`,
/// because the bases they fold may instead have come from the description
/// itself (`Some(seq) => sequence_to_string(seq)`, which is a bare
/// `to_string()` and does not fold). Their trailing `apply_alphabet` covers
/// **both** branches. Swapping those fetches to this function would therefore
/// not let the trailing fold be removed — it would only normalize the fetched
/// branch twice, and removing the fold along with it would silently stop
/// normalizing author-spelled bases. This arm is different precisely because it
/// has no author-spelled branch: an unspelled repeat tract is the only way in,
/// so the fetch is the single source and folding there is what makes the tract
/// search and the unit-match check agree by construction.
fn fetch_normalized_reference_bases<P>(
    provider: &P,
    accession: &str,
    start_one_based: u64,
    end_one_based: u64,
    alphabet: AlphabetMode,
) -> Result<String, ConversionError>
where
    P: ReferenceProvider + ?Sized,
{
    let raw = fetch_reference_bases(provider, accession, start_one_based, end_one_based)?;
    let normalized = apply_alphabet(&raw, alphabet);
    debug_assert_eq!(
        normalized.len(),
        raw.len(),
        "alphabet normalization must preserve byte length"
    );
    Ok(normalized)
}

/// Fetch reference bases for a 1-based inclusive interval `[start, end]` on
/// `accession`. Tries [`ReferenceProvider::get_genomic_sequence`] first
/// (correct for `g.`/`m.` and natural for genomic accessions) and falls
/// back to [`ReferenceProvider::get_sequence`] for transcript accessions
/// (`n.`/`r.`/`c.` SPDI emits on the transcript accession per #116).
///
/// The provider takes 0-based half-open coordinates, so the conversion is
/// `[start - 1, end)`.
///
/// Returns [`ConversionError::MissingReferenceData`] when neither call
/// returns data, or when the returned string length does not match the
/// requested interval length (insufficient ref data near a boundary).
fn fetch_reference_bases<P>(
    provider: &P,
    accession: &str,
    start_one_based: u64,
    end_one_based: u64,
) -> Result<String, ConversionError>
where
    P: ReferenceProvider + ?Sized,
{
    if start_one_based < 1 || end_one_based < start_one_based {
        return Err(ConversionError::InvalidPosition {
            description: format!(
                "invalid 1-based interval [{}, {}] for reference fetch",
                start_one_based, end_one_based
            ),
        });
    }
    let zb_start = start_one_based - 1;
    let zb_end = end_one_based;
    let expected_len = (zb_end - zb_start) as usize;

    let bases = match provider.get_genomic_sequence(accession, zb_start, zb_end) {
        Ok(s) => s,
        Err(_) => provider
            .get_sequence(accession, zb_start, zb_end)
            .map_err(|e| ConversionError::MissingReferenceData {
                description: format!(
                    "could not fetch reference for {}:{}-{}: {}",
                    accession, start_one_based, end_one_based, e
                ),
            })?,
    };

    if bases.len() != expected_len {
        return Err(ConversionError::MissingReferenceData {
            description: format!(
                "reference fetch for {}:{}-{} returned {} bases, expected {}",
                accession,
                start_one_based,
                end_one_based,
                bases.len(),
                expected_len
            ),
        });
    }
    Ok(bases)
}

/// Rewrite RNA-alphabet characters to the DNA alphabet for SPDI output.
///
/// SPDI uses the DNA alphabet by convention (RefSeq stores transcript
/// sequences as DNA even on `NR_*` and `NM_*` accessions), so RNA `u`/`U`
/// must become `T`. Other characters are returned unchanged. The output is
/// always uppercase to match SPDI's standard form.
pub(crate) fn apply_alphabet(s: &str, alphabet: AlphabetMode) -> String {
    match alphabet {
        AlphabetMode::Dna => s.to_ascii_uppercase(),
        AlphabetMode::Rna => s
            .chars()
            .map(|c| match c.to_ascii_uppercase() {
                'U' => 'T',
                other => other,
            })
            .collect(),
    }
}

/// Convert an SPDI variant to HGVS genomic format.
///
/// # Arguments
///
/// * `spdi` - The SPDI variant to convert
///
/// # Returns
///
/// * `Ok(HgvsVariant)` - Successfully converted variant
/// * `Err(ConversionError)` - Conversion failed
///
/// # Examples
///
/// ```
/// use ferro_hgvs::spdi::{SpdiVariant, convert::spdi_to_hgvs};
///
/// let spdi = SpdiVariant::new("NC_000001.11", 12344, "A", "G");
/// let hgvs = spdi_to_hgvs(&spdi).unwrap();
/// assert_eq!(hgvs.to_string(), "NC_000001.11:g.12345A>G");
/// ```
pub fn spdi_to_hgvs(spdi: &SpdiVariant) -> Result<HgvsVariant, ConversionError> {
    // Parse the accession using the HGVS parser
    let accession = parse_accession(&spdi.sequence)
        .map(|(_, acc)| acc)
        .map_err(|_| ConversionError::InvalidAccession {
            description: format!("could not parse accession: {}", spdi.sequence),
        })?;

    // Convert 0-based SPDI position to 1-based HGVS position using type-safe conversion
    let spdi_pos_zb = ZeroBasedPos::new(spdi.position);
    let hgvs_pos_ob = spdi_pos_zb.to_one_based();
    let hgvs_pos = hgvs_pos_ob.value();

    // Determine the edit type based on deletion and insertion
    let (interval, edit) = if spdi.is_identity() {
        // A zero-width triple (`NC_000001.11:99::`) names no bases at all, so
        // there is no interval it could be an identity over. Emitting one
        // anyway rendered `g.100=`, an identity asserting a base the triple
        // never claimed — the reverse-direction twin of the invent-a-base
        // error #1362 fixed on the way out. `is_identity()` is `deletion ==
        // insertion`, so a both-empty triple reaches this arm and not the
        // deletion/insertion ones below; refusing here covers it.
        if spdi.deletion.is_empty() {
            return Err(ConversionError::UnsupportedEditType {
                description: "a zero-width SPDI triple names no bases, so it cannot be \
                              converted to an identity over an interval"
                    .to_string(),
            });
        }
        let seq = Some(string_to_sequence(&spdi.deletion)?);
        // An identity claims every base in its triple, so a multi-base one
        // needs a span interval — a point interval would render `g.27GT=`,
        // whose location says one base while its sequence says two. Mirrors
        // the `is_deletion()` arm below.
        let del_len = spdi.deletion.len();
        let interval = if del_len > 1 {
            Interval::new(
                GenomePos::new(hgvs_pos),
                GenomePos::new(hgvs_pos + del_len as u64 - 1),
            )
        } else {
            Interval::point(GenomePos::new(hgvs_pos))
        };
        (
            interval,
            NaEdit::Identity {
                sequence: seq,
                whole_entity: false,
            },
        )
    } else if spdi.deletion.len() == 1 && spdi.insertion.len() == 1 {
        // SNV substitution
        let ref_base = char_to_base(spdi.deletion.chars().next().unwrap())?;
        let alt_base = char_to_base(spdi.insertion.chars().next().unwrap())?;
        (
            Interval::point(GenomePos::new(hgvs_pos)),
            NaEdit::Substitution {
                reference: ref_base,
                alternative: alt_base,
            },
        )
    } else if spdi.is_deletion() {
        // Pure deletion
        let del_len = spdi.deletion.len();
        let del_seq = string_to_sequence(&spdi.deletion)?;
        let interval = if del_len > 1 {
            Interval::new(
                GenomePos::new(hgvs_pos),
                GenomePos::new(hgvs_pos + del_len as u64 - 1),
            )
        } else {
            Interval::point(GenomePos::new(hgvs_pos))
        };
        (
            interval,
            NaEdit::Deletion {
                sequence: Some(del_seq),
                length: None,
            },
        )
    } else if spdi.is_insertion() {
        // Pure insertion. SPDI position N is the 0-based interbase
        // boundary AFTER 1-based base N — i.e. an insertion at SPDI
        // position N corresponds to HGVS `g.N_(N+1)ins{seq}`. Pre-#390
        // this incorrectly used `hgvs_pos = spdi.position + 1`, shifting
        // every emitted ins-form interval by one (`g.(N+1)_(N+2)ins…`).
        // SPDI position 0 represents an insertion before the first
        // base; HGVS has no notation for it, so reject up-front rather
        // than silently emit `g.0_1ins…`.
        if spdi.position == 0 {
            return Err(ConversionError::InvalidPosition {
                description: "SPDI position 0 represents an insertion before the \
                    first base, which has no HGVS notation"
                    .to_string(),
            });
        }
        let ins_seq = string_to_sequence(&spdi.insertion)?;
        (
            Interval::new(
                GenomePos::new(spdi.position),
                GenomePos::new(spdi.position + 1),
            ),
            NaEdit::Insertion {
                sequence: InsertedSequence::Literal(ins_seq),
            },
        )
    } else if !spdi.deletion.is_empty()
        && spdi.deletion.len() == spdi.insertion.len()
        && spdi.deletion.len() >= 2
        && reverse_complement(&spdi.deletion).eq_ignore_ascii_case(&spdi.insertion)
    {
        // Inversion recovery (#270): SPDI delins where the inserted seq
        // is the reverse complement of the deleted seq is canonically a
        // tandem inversion. Length-1 cases are SNVs and handled above;
        // self-RC palindromes of length 2+ (e.g. `AT:AT`, `GC:GC`) only
        // recover when the seqs differ from a plain identity — which
        // they don't, so identity catches them first.
        let inv_seq = string_to_sequence(&spdi.deletion)?;
        let del_len = spdi.deletion.len();
        let interval = Interval::new(
            GenomePos::new(hgvs_pos),
            GenomePos::new(hgvs_pos + del_len as u64 - 1),
        );
        (
            interval,
            NaEdit::Inversion {
                sequence: Some(inv_seq),
                length: None,
            },
        )
    } else {
        // Delins (different lengths or MNV)
        let del_len = spdi.deletion.len();
        let ins_seq = string_to_sequence(&spdi.insertion)?;
        let interval = if del_len > 1 {
            Interval::new(
                GenomePos::new(hgvs_pos),
                GenomePos::new(hgvs_pos + del_len as u64 - 1),
            )
        } else {
            Interval::point(GenomePos::new(hgvs_pos))
        };
        (
            interval,
            NaEdit::Delins {
                sequence: InsertedSequence::Literal(ins_seq),
                deleted: None,
                deleted_length: None,
                substitution_reference: None,
            },
        )
    };

    // m. coord system recovery (#270): NC_012920.* is the canonical
    // human mitochondrial accession; SPDI carries no coord-system tag
    // so the bare reverse path defaults to g. Emit m. when the
    // accession matches a known mitochondrial reference.
    if is_mitochondrial_accession(&accession) {
        return Ok(HgvsVariant::Mt(MtVariant {
            accession,
            gene_symbol: None,
            loc_edit: LocEdit::new(interval, edit),
        }));
    }

    Ok(HgvsVariant::Genome(GenomeVariant {
        accession,
        gene_symbol: None,
        loc_edit: LocEdit::new(interval, edit),
    }))
}

/// Returns true if `accession` is a known mitochondrial reference. Used
/// by `spdi_to_hgvs` to emit `m.` instead of `g.` for these accessions
/// (SPDI carries no coord-system tag). Delegates to
/// [`Accession::is_mitochondrial`], the single source for the accession list.
fn is_mitochondrial_accession(accession: &Accession) -> bool {
    accession.is_mitochondrial()
}

/// Convert an SPDI variant to HGVS, using a reference provider to recover
/// duplication form for SPDI insertions whose inserted sequence equals the
/// immediately-5' reference flank.
///
/// SPDI is a canonical/positional format and represents duplications as
/// insertions of the duplicated sequence. This function performs the inverse
/// recovery for HGVS round-trips: given an SPDI insertion, if the
/// `inserted.len()` reference bases ending at the insertion point equal the
/// inserted sequence (case-insensitive), the change is a tandem duplication
/// per the HGVS spec
/// (`assets/hgvs-nomenclature/docs/recommendations/DNA/duplication.md`) and
/// is rendered as HGVS `dup` rather than `ins`.
///
/// All other SPDI shapes (substitution, deletion, delins, identity) are
/// returned unchanged from [`spdi_to_hgvs`].
///
/// SPDI itself is canonical and places the insertion at its rightmost
/// shiftable offset, so the matched 5' flank is already the most-3' position;
/// no separate 3' shift is performed.
///
/// SPDI→HGVS only produces genomic (`g.`) variants, so dup recovery applies
/// uniformly to genomic and mitochondrial accessions but does not extend to
/// `c.`, `n.`, or `r.` coordinate systems via this function.
///
/// # Arguments
///
/// * `spdi` - The SPDI variant to convert.
/// * `reference` - Reference provider used to fetch the 5'-flanking bases
///   for dup detection. The function tries
///   [`ReferenceProvider::get_genomic_sequence`] first and falls back to
///   [`ReferenceProvider::get_sequence`] for providers that only expose
///   sequences via the transcript path.
///
/// # Errors
///
/// Returns [`ConversionError::MissingReferenceData`] when the reference
/// provider returns an error fetching the 5' flank. Returns the same errors
/// as [`spdi_to_hgvs`] for other failure modes. A short fetch (truncated near
/// a contig boundary) or a non-matching flank silently falls back to the
/// ins-form result, matching the spec recommendation that, absent evidence
/// of tandem flanking, the change is described as an insertion.
///
/// [`ReferenceProvider`]: crate::reference::provider::ReferenceProvider
/// [`ReferenceProvider::get_genomic_sequence`]: crate::reference::provider::ReferenceProvider::get_genomic_sequence
/// [`ReferenceProvider::get_sequence`]: crate::reference::provider::ReferenceProvider::get_sequence
pub fn spdi_to_hgvs_with_ref<R>(
    spdi: &SpdiVariant,
    reference: &R,
) -> Result<HgvsVariant, ConversionError>
where
    R: crate::reference::provider::ReferenceProvider + ?Sized,
{
    // Build the base HGVS variant using the existing reference-free path.
    let base = spdi_to_hgvs(spdi)?;

    // Only insertions are candidates for dup recovery. Everything else
    // (substitution, deletion, delins, identity) passes through unchanged.
    if !spdi.is_insertion() {
        return Ok(base);
    }

    // Try to recover dup form. If the inserted sequence does not match the
    // 5' flank, fall through and return the original `ins`-form variant.
    if let Some(dup_variant) = recover_dup_from_insertion(spdi, reference, &base)? {
        return Ok(dup_variant);
    }
    Ok(base)
}

/// If the SPDI insertion's inserted sequence equals the immediately-5'
/// reference flank (case-insensitive), return an `HgvsVariant` whose edit is
/// `NaEdit::Duplication` over the corresponding 1-based interval. Returns
/// `Ok(None)` when the bases do not match or there are not enough preceding
/// bases. Returns `Err` only when the reference provider returns a hard
/// fetch error (which is propagated as `MissingReferenceData`).
///
/// `base` is the already-built `ins`-form `HgvsVariant`; we reuse its
/// accession and gene_symbol when constructing the dup-form result.
fn recover_dup_from_insertion<R>(
    spdi: &SpdiVariant,
    reference: &R,
    base: &HgvsVariant,
) -> Result<Option<HgvsVariant>, ConversionError>
where
    R: crate::reference::provider::ReferenceProvider + ?Sized,
{
    debug_assert!(spdi.is_insertion());

    let ins = &spdi.insertion;
    let ins_len = ins.len() as u64;
    if ins_len == 0 {
        return Ok(None);
    }

    // Need `ins_len` bases of preceding reference. SPDI position N is
    // the 0-based interbase boundary AFTER 1-based base N, so the
    // 5'-flanking window is the `ins_len` bases ending at HGVS 1-based
    // position `spdi.position`. In the 0-based half-open form used by
    // `get_genomic_sequence` that's `[spdi.position - ins_len, spdi.position)`.
    let flank_end = spdi.position;
    if flank_end < ins_len {
        // Not enough preceding bases (insertion is too close to contig 5' end).
        return Ok(None);
    }
    let flank_start = flank_end - ins_len;

    // Fetch the flanking sequence. Try `get_genomic_sequence` first (the
    // SPDI accession is genomic), then fall back to `get_sequence` for
    // providers that store contigs as transcripts.
    let flank = match reference.get_genomic_sequence(&spdi.sequence, flank_start, flank_end) {
        Ok(s) => s,
        Err(_) => match reference.get_sequence(&spdi.sequence, flank_start, flank_end) {
            Ok(s) => s,
            Err(e) => {
                return Err(ConversionError::MissingReferenceData {
                    description: format!(
                        "could not fetch 5' flank for {}:{}-{}: {}",
                        spdi.sequence, flank_start, flank_end, e
                    ),
                });
            }
        },
    };

    // The fetched window must match `ins` (case-insensitive) and must be
    // exactly `ins_len` bases. A short fetch (e.g., truncated near a contig
    // boundary) means we cannot prove tandem dup → fall back to ins.
    if flank.len() as u64 != ins_len {
        return Ok(None);
    }
    if !flank.eq_ignore_ascii_case(ins) {
        return Ok(None);
    }

    // Build the dup edit. The SPDI 0-based interbase position N and
    // the HGVS 1-based base-N coordinate share the same numeric
    // value (N), even though they describe different things —
    // `spdi.position` (interbase) sits AFTER 1-based base
    // `spdi.position`, which is also the 1-based end of the
    // duplicated region.
    let end_one_based = flank_end; // numerically equal to spdi.position
    let start_one_based = end_one_based + 1 - ins_len;

    let dup_seq = string_to_sequence(ins)?;
    let interval = if ins_len == 1 {
        Interval::point(GenomePos::new(end_one_based))
    } else {
        Interval::new(
            GenomePos::new(start_one_based),
            GenomePos::new(end_one_based),
        )
    };
    let edit = NaEdit::Duplication {
        sequence: Some(dup_seq),
        length: None,
        uncertain_extent: None,
    };

    // Reuse the accession + gene_symbol from the base ins-form variant.
    // SPDI→HGVS produces either Genome or Mt (mitochondrial), depending
    // on the accession; preserve whichever shape the base has.
    match base {
        HgvsVariant::Genome(g) => Ok(Some(HgvsVariant::Genome(GenomeVariant {
            accession: g.accession.clone(),
            gene_symbol: g.gene_symbol.clone(),
            loc_edit: LocEdit::new(interval, edit),
        }))),
        HgvsVariant::Mt(m) => Ok(Some(HgvsVariant::Mt(MtVariant {
            accession: m.accession.clone(),
            gene_symbol: m.gene_symbol.clone(),
            loc_edit: LocEdit::new(interval, edit),
        }))),
        _ => Ok(None),
    }
}

/// Helper to convert a string to a Sequence.
fn string_to_sequence(s: &str) -> Result<Sequence, ConversionError> {
    s.parse().map_err(|_| ConversionError::InvalidPosition {
        description: format!("invalid sequence: {}", s),
    })
}

/// Helper to convert a char to a Base.
fn char_to_base(c: char) -> Result<crate::hgvs::edit::Base, ConversionError> {
    crate::hgvs::edit::Base::from_char(c).ok_or_else(|| ConversionError::InvalidPosition {
        description: format!("invalid base character: {}", c),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;

    /// Build an `n.*N` substitution.
    ///
    /// #1748 refuses that spelling at parse in **every** mode
    /// (`background/numbering.md:52` puts no `*` zone on the non-coding axis),
    /// so a downstream fixture can no longer be parsed. It still has to be
    /// *handled* here: [`TxPos::downstream`](crate::hgvs::location::TxPos::downstream)
    /// is public API, so a caller can hand the conversion path one of these
    /// without the parser ever seeing it. Constructing it is therefore the
    /// faithful fixture, not a workaround.
    fn tx_downstream_sub(accession: &str, base: i64) -> crate::hgvs::HgvsVariant {
        use crate::hgvs::interval::TxInterval;
        use crate::hgvs::location::TxPos;
        let mut variant =
            parse_hgvs(&format!("{accession}:n.{base}A>G")).expect("the template parses");
        let crate::hgvs::HgvsVariant::Tx(tx) = &mut variant else {
            unreachable!("an `n.` description parses as Tx")
        };
        tx.loc_edit.location = TxInterval::point(TxPos::downstream(base));
        variant
    }

    // ------------------------------------------------------------------
    // Unspelled-bases declines (#1388)
    //
    // These pin the message *text*, not just that an error occurs. The arms
    // were already covered for behaviour, which is exactly why a wrapping
    // defect could sit in a branch through a full green suite: nothing read
    // what they said.
    // ------------------------------------------------------------------

    /// Every example offered in a decline must be a description ferro accepts.
    /// Suggesting a form the parser rejects would send the caller somewhere
    /// worse than saying nothing, and the suggestion is a literal that no other
    /// test exercises.
    #[test]
    fn spelled_examples_are_parseable() {
        for what in [
            UnspelledBases::Duplicated,
            UnspelledBases::Deleted,
            UnspelledBases::DeletedInDelins,
            UnspelledBases::Unchanged,
            UnspelledBases::Inverted,
            UnspelledBases::RepeatTract,
        ] {
            let Some(example) = what.spelled_example() else {
                // `RepeatTract` offers no example on purpose — the notation
                // cannot carry the bases. Assert that rather than skipping it,
                // so adding an example here without a parseable form fails.
                assert_eq!(
                    what,
                    UnspelledBases::RepeatTract,
                    "only the repeat tract may decline without offering an example"
                );
                continue;
            };
            let prefixed = format!("NC_000001.11:{example}");
            assert!(
                parse_hgvs(&prefixed).is_ok(),
                "{what:?} suggests `{example}`, which ferro cannot parse as `{prefixed}`"
            );
        }
    }

    /// The repeat arm names the provider as the only route, and does not tell the
    /// caller to spell bases the notation cannot express.
    #[test]
    fn the_repeat_tract_offers_only_the_provider() {
        let message = unspelled_bases_error(UnspelledBases::RepeatTract, 10, 12).to_string();
        assert!(message.contains("10..=12"), "no span named: {message}");
        assert!(
            message.contains("reference provider"),
            "the provider route must be named: {message}"
        );
        assert!(
            !message.contains("Spell them"),
            "a repeat's tract cannot be spelled in the description, so the message \
             must not suggest it: {message}"
        );
    }

    /// The rendered message must carry no accidental whitespace runs. A wrapped
    /// literal without `\` continuations silently embeds the source indentation,
    /// which is invisible in review and obvious to a user.
    #[test]
    fn the_decline_message_has_no_embedded_whitespace_runs() {
        // Every shape, not a representative one: the defect is per-literal, and
        // each shape's message is assembled from a different pair of literals.
        for what in [
            UnspelledBases::Duplicated,
            UnspelledBases::Deleted,
            UnspelledBases::DeletedInDelins,
            UnspelledBases::Unchanged,
            UnspelledBases::Inverted,
            UnspelledBases::RepeatTract,
        ] {
            let message = unspelled_bases_error(what, 10, 12).to_string();
            assert!(
                !message.contains("  "),
                "{what:?} carries a run of consecutive spaces: {message:?}"
            );
            assert!(
                !message.contains('\n') && !message.contains('\t'),
                "{what:?} carries a newline or tab: {message:?}"
            );
        }
    }

    /// The whole point of the message: it names the interval that could not be
    /// resolved, and both ways out.
    #[test]
    fn the_decline_message_names_the_span_and_both_remedies() {
        let message = unspelled_bases_error(UnspelledBases::Unchanged, 10, 12).to_string();
        assert!(message.contains("10..=12"), "no span named: {message}");
        assert!(message.contains("unchanged"), "no shape named: {message}");
        assert!(message.contains("g.10A="), "no example offered: {message}");
        assert!(
            message.contains("Spell them") && message.contains("reference provider"),
            "both remedies must be offered: {message}"
        );
    }

    /// The six shapes must be distinguishable from their messages alone — a
    /// caller reading a log has nothing else to go on.
    #[test]
    fn each_shape_declines_distinguishably() {
        let messages: Vec<String> = [
            UnspelledBases::Duplicated,
            UnspelledBases::Deleted,
            UnspelledBases::DeletedInDelins,
            UnspelledBases::Unchanged,
            UnspelledBases::Inverted,
            UnspelledBases::RepeatTract,
        ]
        .iter()
        .map(|w| unspelled_bases_error(*w, 10, 12).to_string())
        .collect();

        // `Deleted` and `DeletedInDelins` share an adjective by design — both
        // are the deleted span — but differ in the example they offer.
        let unique: std::collections::BTreeSet<&String> = messages.iter().collect();
        assert_eq!(
            unique.len(),
            messages.len(),
            "two shapes decline identically: {messages:#?}"
        );
    }

    /// End-to-end: the message a real provider-less conversion produces, for
    /// each shape that can hit the wall. Pins that the helper is actually wired
    /// in — a unit test on the helper alone would pass with the arms unchanged.
    ///
    /// All **six** shapes, including the repeat tract. It is the one arm whose
    /// remedy differs, so leaving it to the helper-level tests would have left
    /// the only differing branch unproven end-to-end — and it is reachable
    /// without a provider exactly like its siblings.
    #[test]
    fn a_provider_less_conversion_declines_with_the_span() {
        for (descriptor, adjective, span) in [
            ("NC_000001.11:g.10_12dup", "duplicated", "10..=12"),
            ("NC_000001.11:g.10_12del", "deleted", "10..=12"),
            ("NC_000001.11:g.10_12delinsT", "deleted", "10..=12"),
            ("NC_000001.11:g.10_12=", "unchanged", "10..=12"),
            ("NC_000001.11:g.10_12inv", "inverted", "10..=12"),
            (
                "NC_000001.11:g.10_15AC[3]",
                "pre-expansion repeat-tract",
                "10..=15",
            ),
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            let message = hgvs_to_spdi_simple(&variant)
                .expect_err(&format!("`{descriptor}` must decline without a provider"))
                .to_string();
            assert!(
                message.contains(adjective) && message.contains(span),
                "`{descriptor}` declined without naming the {adjective} bases at {span}: {message}"
            );
        }
    }

    // HGVS to SPDI tests

    #[test]
    fn test_hgvs_to_spdi_substitution() {
        let hgvs = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        assert_eq!(spdi.position, 12344);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "G");
    }

    #[test]
    fn test_hgvs_to_spdi_insertion() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_101insATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        // SPDI 0-based interbase: position 100 is the boundary AFTER
        // 1-based base 100, matching HGVS g.100_101ins (closes #390).
        assert_eq!(spdi.position, 100);
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn test_hgvs_to_spdi_deletion_with_seq() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_102delATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn test_hgvs_to_spdi_deletion_without_seq() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_102del").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    #[test]
    fn test_hgvs_to_spdi_delins_without_ref() {
        // Without reference data, delins with unknown deletion sequence returns error
        let hgvs = parse_hgvs("NC_000001.11:g.100_102delinsTTCC").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("unknown"));
    }

    #[test]
    fn test_hgvs_to_spdi_delins_with_explicit_deleted_no_ref() {
        // Issue #120: when the input carries an explicit deleted sequence, the
        // SPDI conversion can succeed without reference data.
        let hgvs = parse_hgvs("NC_000001.11:g.100_102delATGinsTTCC").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs)
            .expect("explicit deleted sequence should not require reference data");
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "TTCC");
    }

    #[test]
    fn test_hgvs_to_spdi_duplication_with_seq() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_102dupATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        // Dup becomes insertion after the duplicated region; SPDI
        // interbase position 102 sits AFTER 1-based base 102 (#390).
        assert_eq!(spdi.position, 102);
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn test_hgvs_to_spdi_identity() {
        let hgvs = parse_hgvs("NC_000001.11:g.100A=").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "A");
    }

    // ------------------------------------------------------------------
    // Reference-window normalization (#1452)
    //
    // The soft-masking behaviour these underwrite is graded end-to-end in
    // `tests/it/issue_1452_soft_masked_repeat_span.rs`. What lives here are
    // the two properties that cannot be isolated from there: the length
    // invariant the anchor offset depends on, and the `r.` axis, which the
    // integration suite reaches only through a full transcript projection.
    // ------------------------------------------------------------------

    /// `resolve_repeat_tract_span` computes `anchor_offset` as a byte index
    /// into the normalized window, so a normalization that changed the byte
    /// length would silently point the tract search at the wrong base rather
    /// than fail. Both `apply_alphabet` arms rewrite ASCII one-for-one and
    /// leave every other byte alone, and this is what holds them to it.
    #[test]
    fn normalizing_a_window_preserves_its_length() {
        for window in ["acgtACGT", "uuuUUU", "acgunACGUN", "", "nnnn"] {
            for alphabet in [AlphabetMode::Dna, AlphabetMode::Rna] {
                assert_eq!(
                    apply_alphabet(window, alphabet).len(),
                    window.len(),
                    "`{window}` changed length under {alphabet:?}"
                );
            }
        }
    }

    /// The window is normalized with `apply_alphabet`, not `to_ascii_uppercase`,
    /// and on the `r.` axis that is the difference between finding a tract and
    /// truncating it.
    ///
    /// A repeat unit reaches `resolve_repeat_tract_span` having already had `U`
    /// rewritten to `T` for SPDI's DNA alphabet, so a `u`-spelled reference must
    /// get the same rewrite or the search compares `T` against `U` and finds no
    /// run — the identical silent-truncation route #1452 fixed for case. The
    /// `Dna` row is the control: there `U` is not `T`, so no tract exists and
    /// the unit-wide fallback is the right answer.
    #[test]
    fn a_uracil_spelled_tract_is_found_on_the_rna_axis() {
        let mut provider = MockProvider::new();
        // A 4-copy `u` tract at 3..=6, lower-case as well, so this covers both
        // normalizations at once.
        provider.add_genomic_sequence("NR_TEST.1", "GGuuuuGG".to_string());

        assert_eq!(
            resolve_repeat_tract_span(&provider, "NR_TEST.1", 3, b"T", AlphabetMode::Rna).unwrap(),
            (3, 6, RepeatSpanOrigin::Tract),
            "the `r.` axis must see the tract its unit was rewritten to match"
        );
        assert_eq!(
            resolve_repeat_tract_span(&provider, "NR_TEST.1", 3, b"T", AlphabetMode::Dna).unwrap(),
            (3, 3, RepeatSpanOrigin::NoRunAtAnchor),
            "on the DNA axis `U` is not `T`, so no tract exists and the \
             unit-wide fallback stands"
        );
    }

    /// A 28-base `ACGT…` contig, so a 1-based position `p` holds
    /// `"ACGT"[(p - 1) % 4]` — base 10 is `C`, and 10..12 is `CGT`.
    fn identity_provider() -> MockProvider {
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_000001.11", "ACGTACGTACGTACGTACGTACGTACGT".to_string());
        provider
    }

    #[test]
    fn an_unspelled_identity_takes_its_base_from_the_reference() {
        // `g.10=` states that base 10 is unchanged without naming it. Every
        // other arm that can omit its sequence (deletion, delins, inversion)
        // fetches it from the provider; this one used to default to the empty
        // string, yielding a zero-width `99::` that claims no base at all.
        let hgvs = parse_hgvs("NC_000001.11:g.10=").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &identity_provider()).unwrap();
        assert_eq!(spdi.position, 9);
        assert_eq!(spdi.deletion, "C");
        assert_eq!(spdi.insertion, "C");
    }

    #[test]
    fn an_unspelled_identity_takes_its_whole_span_from_the_reference() {
        let hgvs = parse_hgvs("NC_000001.11:g.10_12=").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &identity_provider()).unwrap();
        assert_eq!(spdi.position, 9);
        assert_eq!(spdi.deletion, "CGT");
        assert_eq!(spdi.insertion, "CGT");
    }

    /// A same-reference position-range insertion (`ins{start}_{end}`) names its
    /// inserted bases by position rather than spelling them, so `hgvs_to_spdi`
    /// must read them from the reference — exactly as `del`/`dup` read their
    /// omitted bases. On the `ACGT…` contig positions 5..=8 hold `ACGT`.
    #[test]
    fn a_reference_span_insertion_reads_its_bases_from_the_reference() {
        let hgvs = parse_hgvs("NC_000001.11:g.10_11ins5_8").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &identity_provider()).unwrap();
        // SPDI interbase position 10 sits AFTER 1-based base 10 (#390).
        assert_eq!(spdi.position, 10);
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ACGT");
    }

    /// The inverted form (`ins{start}_{end}inv`) inserts the reverse complement
    /// of the named span. Positions 5..=7 hold `ACG`, whose reverse complement
    /// is `CGT` (a non-palindromic span, so the inversion is observable).
    #[test]
    fn an_inverted_reference_span_insertion_reads_the_reverse_complement() {
        let hgvs = parse_hgvs("NC_000001.11:g.10_11ins5_7inv").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &identity_provider()).unwrap();
        assert_eq!(spdi.position, 10);
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "CGT");
    }

    /// A range insert names bases in the reference, so without a provider to
    /// read them the conversion must decline rather than invent them.
    #[test]
    fn a_reference_span_insertion_without_a_provider_is_declined() {
        let hgvs = parse_hgvs("NC_000001.11:g.10_11ins5_8").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(
            matches!(result, Err(ConversionError::MissingReferenceData { .. })),
            "a same-reference range insert cannot be resolved without a provider: {result:?}"
        );
    }

    /// The `delins` arm resolves a same-reference range insert on the same path,
    /// fetching both the deleted span and the inserted span from the reference.
    /// del at 10..=12 is `CGT`; ins 5..=8 is `ACGT`.
    #[test]
    fn a_reference_span_delins_reads_both_spans_from_the_reference() {
        let hgvs = parse_hgvs("NC_000001.11:g.10_12delins5_8").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &identity_provider()).unwrap();
        assert_eq!(spdi.position, 9);
        assert_eq!(spdi.deletion, "CGT");
        assert_eq!(spdi.insertion, "ACGT");
    }

    /// An exact tandem-repeat insertion (`ins{unit}[{n}]`) names its inserted
    /// bases by a spelled unit and an exact copy count. `hgvs_to_spdi` expands
    /// the unit `n` times — no provider is needed, since the bases are named in
    /// the description itself — so the simple path resolves it.
    #[test]
    fn an_exact_repeat_insertion_expands_its_unit() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_101insC[4]").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        // SPDI interbase position 100 sits AFTER 1-based base 100 (#390).
        assert_eq!(spdi.position, 100);
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "CCCC");
    }

    /// A multi-base unit (`ins{unit}[{n}]` with `unit` longer than one base)
    /// repeats the whole unit `n` times.
    #[test]
    fn an_exact_multibase_repeat_insertion_expands_the_whole_unit() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_101insCG[3]").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.position, 100);
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "CGCGCG");
    }

    /// An uncertain or range copy count (`ins{unit}[{min}_{max}]`, `[?]`) names
    /// no single expansion, so it stays undetermined and must decline.
    #[test]
    fn an_uncertain_repeat_count_insertion_is_declined() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_101insC[10_15]").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(
            matches!(result, Err(ConversionError::MissingReferenceData { .. })),
            "a range repeat count names no single expansion: {result:?}"
        );
    }

    /// The `delins` arm expands an exact repeat insert on the same path. del at
    /// 10..=12 is `CGT` (read from the reference); `insCG[3]` expands to
    /// `CGCGCG`.
    #[test]
    fn an_exact_repeat_delins_expands_its_unit() {
        let hgvs = parse_hgvs("NC_000001.11:g.10_12delinsCG[3]").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &identity_provider()).unwrap();
        assert_eq!(spdi.position, 9);
        assert_eq!(spdi.deletion, "CGT");
        assert_eq!(spdi.insertion, "CGCGCG");
    }

    /// A compound insert bracket mixing a literal with a same-reference
    /// coordinate span (`delins[{literal};{start}_{end}]`) is one contiguous
    /// inserted sequence once each part is read out. This is the shape the
    /// original request named — `185_201delins[T;213_271]` — reduced onto the
    /// `ACGT…` contig: del at 10..=12 is `CGT`; the payload is the literal `T`
    /// followed by the span 5..=8 (`ACGT`), so the insertion is `TACGT`.
    #[test]
    fn a_compound_insert_mixing_a_literal_and_a_reference_span_reads_both() {
        let hgvs = parse_hgvs("NC_000001.11:g.10_12delins[T;5_8]").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &identity_provider()).unwrap();
        assert_eq!(spdi.position, 9);
        assert_eq!(spdi.deletion, "CGT");
        assert_eq!(spdi.insertion, "TACGT");
    }

    /// The parts are concatenated in the order written, so swapping them swaps
    /// the resulting bases — `[5_8;T]` is `ACGTT`, not `TACGT`. Pins that the
    /// compound resolver preserves part order rather than sorting or merging.
    #[test]
    fn a_compound_insert_concatenates_its_parts_in_order() {
        let hgvs = parse_hgvs("NC_000001.11:g.10_12delins[5_8;T]").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &identity_provider()).unwrap();
        assert_eq!(spdi.position, 9);
        assert_eq!(spdi.deletion, "CGT");
        assert_eq!(spdi.insertion, "ACGTT");
    }

    /// An inverted span inside a compound bracket reads the reverse complement
    /// of that part only, leaving the surrounding literal untouched. Span 5..=7
    /// is `ACG`, reverse-complemented to `CGT`, so `[A;5_7inv]` is `ACGT`.
    #[test]
    fn a_compound_insert_reverse_complements_an_inverted_span_part() {
        let hgvs = parse_hgvs("NC_000001.11:g.10_12delins[A;5_7inv]").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &identity_provider()).unwrap();
        assert_eq!(spdi.position, 9);
        assert_eq!(spdi.deletion, "CGT");
        assert_eq!(spdi.insertion, "ACGT");
    }

    /// A compound bracket may combine an exact repeat part with a literal; the
    /// repeat expands through the same leaf a bare `insC[3]` uses. `[C[3];A]`
    /// is `CCCA`.
    #[test]
    fn a_compound_insert_expands_an_exact_repeat_part() {
        let hgvs = parse_hgvs("NC_000001.11:g.10_12delins[C[3];A]").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &identity_provider()).unwrap();
        assert_eq!(spdi.position, 9);
        assert_eq!(spdi.deletion, "CGT");
        assert_eq!(spdi.insertion, "CCCA");
    }

    /// A compound insert on the `ins` arm resolves on the same path as `delins`,
    /// with no deletion. `[T;5_8]` at `g.10_11` inserts `TACGT` after base 10.
    #[test]
    fn a_compound_insertion_reads_on_the_ins_arm_too() {
        let hgvs = parse_hgvs("NC_000001.11:g.10_11ins[T;5_8]").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &identity_provider()).unwrap();
        assert_eq!(spdi.position, 10);
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "TACGT");
    }

    /// A compound insert whose only undetermined content is a coordinate span
    /// still needs a provider to read that span, so the provider-less path must
    /// decline rather than emit a triple missing those bases.
    #[test]
    fn a_compound_insert_with_a_span_part_needs_a_provider() {
        let hgvs = parse_hgvs("NC_000001.11:g.10_12delins[T;5_8]").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(
            matches!(result, Err(ConversionError::MissingReferenceData { .. })),
            "a coordinate-span part cannot be resolved without a provider: {result:?}"
        );
    }

    /// A part whose bases are genuinely undetermined declines the whole insert.
    /// A non-exact repeat count (`N[10_15]`) inside the bracket names no single
    /// expansion, so even with a provider the compound insert is unsupported.
    #[test]
    fn a_compound_insert_with_an_undetermined_part_is_declined() {
        let hgvs = parse_hgvs("NC_000001.11:g.10_12delins[N[10_15];A]").unwrap();
        let result = hgvs_to_spdi(&hgvs, &identity_provider());
        assert!(
            matches!(result, Err(ConversionError::UnsupportedEditType { .. })),
            "a non-exact repeat count part names no single expansion: {result:?}"
        );
    }

    #[test]
    fn an_unspelled_identity_fetches_on_the_non_genomic_axes_too() {
        // The arm is shared by g./m./n./r./c., and the r. path additionally
        // runs the fetched bases through `apply_alphabet(_, Rna)`. Each axis is
        // asserted against its own `del` sibling on the same position, since
        // that sibling is the fetch behavior this arm was made to match.
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NR_TEST.1", "ACGTACGTACGTACGTACGTACGTACGT".to_string());
        for (identity, deletion) in [
            ("NR_TEST.1:n.10=", "NR_TEST.1:n.10del"),
            ("NR_TEST.1:r.10=", "NR_TEST.1:r.10del"),
        ] {
            let id_spdi = hgvs_to_spdi(&parse_hgvs(identity).unwrap(), &provider).unwrap();
            let del_spdi = hgvs_to_spdi(&parse_hgvs(deletion).unwrap(), &provider).unwrap();
            assert_eq!(
                (id_spdi.position, id_spdi.deletion.as_str()),
                (del_spdi.position, del_spdi.deletion.as_str()),
                "{identity} must claim the same span `{deletion}` deletes"
            );
            // An identity keeps what it claims; the deletion drops it.
            assert_eq!(
                id_spdi.insertion, id_spdi.deletion,
                "{identity} is an identity"
            );
            assert_eq!(del_spdi.insertion, "", "{deletion} deletes");
        }
    }

    #[test]
    fn an_unspelled_identity_needs_reference_data() {
        // The provider-less path cannot know which bases the span holds, so it
        // reports that rather than emitting a triple whose span is wrong. This
        // matches `g.10del` on the same path.
        let hgvs = parse_hgvs("NC_000001.11:g.10=").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(
            matches!(result, Err(ConversionError::MissingReferenceData { .. })),
            "expected MissingReferenceData, got {result:?}"
        );
    }

    #[test]
    fn an_out_of_range_unspelled_identity_reports_the_fetch_failure() {
        // The error path of the fetch above. An identity past the end of the
        // contig has no bases to claim, so the fetch fails and that failure is
        // propagated — the same answer `g.9999del` already gives, rather than a
        // panic or a triple invented from an empty span. Before this arm
        // consulted the provider it returned `Ok(9998::)` here, so this also
        // pins that an out-of-range identity is no longer silently accepted.
        let hgvs = parse_hgvs("NC_000001.11:g.9999=").unwrap();
        let result = hgvs_to_spdi(&hgvs, &identity_provider());
        assert!(
            matches!(result, Err(ConversionError::MissingReferenceData { .. })),
            "expected MissingReferenceData, got {result:?}"
        );
    }

    #[test]
    fn a_whole_entity_identity_has_no_spdi() {
        // `g.=` asserts the *entire* reference is unchanged and names no
        // interval, so there is no position SPDI could honestly carry. It used
        // to emit `0::`, which `spdi_to_hgvs` reads back as `g.1=` — turning a
        // statement about the whole sequence into one about base 1. Declining
        // is the same answer this module gives every other edit whose shape
        // SPDI cannot encode.
        let hgvs = parse_hgvs("NC_000001.11:g.=").unwrap();
        let result = hgvs_to_spdi(&hgvs, &identity_provider());
        assert!(
            matches!(result, Err(ConversionError::UnsupportedEditType { .. })),
            "expected UnsupportedEditType, got {result:?}"
        );
    }

    #[test]
    fn an_identity_member_does_not_alias_a_sibling_insertion_junction() {
        // The defect this fix exists for. `g.[261_262dup;263=]` is a real
        // normalizer output (5' shuffle, #1321's split spelling). With the
        // identity converting to a zero-width triple, both members landed at
        // interbase 262 — the dup's junction and the identity's empty span —
        // and an applier walking the members cannot tell that apart from two
        // insertions competing for one interbase, so it has to decline the
        // whole description. Giving the identity its real span separates them.
        let provider = identity_provider();
        let variant = parse_hgvs("NC_000001.11:g.[10_11dup;12=]").unwrap();
        let members = match variant {
            HgvsVariant::Allele(allele) => allele.variants,
            other => panic!("expected an allele, got {other}"),
        };
        let triples: Vec<SpdiVariant> = members
            .iter()
            .map(|m| hgvs_to_spdi(m, &provider).expect("member converts"))
            .collect();
        // The dup copies bases 10-11 in at the junction after 11 (interbase 11).
        assert_eq!(triples[0].position, 11);
        assert_eq!(triples[0].deletion, "");
        assert_eq!(triples[0].insertion, "CG");
        // The identity claims base 12 itself: interbase 11..12, not 11..11.
        assert_eq!(triples[1].position, 11);
        assert_eq!(
            triples[1].deletion, "T",
            "the identity must claim its base; a zero-width deletion here is \
             indistinguishable from a second insertion at interbase {}",
            triples[0].position
        );
    }

    #[test]
    fn test_hgvs_to_spdi_simple_cds_requires_provider() {
        // c. variants need transcript metadata to resolve to a transcript
        // position; the simple path therefore returns ProviderRequired.
        let hgvs = parse_hgvs("NM_000088.3:c.100A>G").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::ProviderRequired { .. })
        ));
        let err = result.unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("c."), "message should mention c.: {}", msg);
        assert!(
            msg.contains("provider"),
            "message should mention provider: {}",
            msg
        );
    }

    #[test]
    fn test_hgvs_to_spdi_simple_short_form_inversion_requires_provider() {
        // Short-form inversion (no explicit sequence) cannot determine the
        // reference bases without a provider. Pinned audit: the simple path
        // surfaces MissingReferenceData rather than UnsupportedEditType
        // (the prior, pre-#118 behaviour).
        let hgvs = parse_hgvs("NC_000001.11:g.100_200inv").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    // SPDI to HGVS tests

    #[test]
    fn test_spdi_to_hgvs_substitution() {
        let spdi = SpdiVariant::new("NC_000001.11", 12344, "A", "G");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.12345A>G");
    }

    #[test]
    fn test_spdi_to_hgvs_deletion() {
        let spdi = SpdiVariant::deletion("NC_000001.11", 99, "ATG");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100_102delATG");
    }

    #[test]
    fn test_spdi_to_hgvs_insertion() {
        let spdi = SpdiVariant::insertion("NC_000001.11", 100, "ATG");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        // SPDI 100 = boundary AFTER 1-based 100 = HGVS g.100_101ins (#390).
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100_101insATG");
    }

    #[test]
    fn test_spdi_to_hgvs_delins() {
        let spdi = SpdiVariant::delins("NC_000001.11", 99, "ATG", "TTCC");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100_102delinsTTCC");
    }

    #[test]
    fn test_spdi_to_hgvs_identity() {
        let spdi = SpdiVariant::new("NC_000001.11", 99, "A", "A");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100A=");
    }

    #[test]
    fn test_spdi_to_hgvs_single_del() {
        let spdi = SpdiVariant::deletion("NC_000001.11", 99, "A");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100delA");
    }

    // Roundtrip tests

    #[test]
    fn test_roundtrip_substitution() {
        let original = "NC_000001.11:g.12345A>G";
        let hgvs = parse_hgvs(original).unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(back.to_string(), original);
    }

    #[test]
    fn test_roundtrip_insertion() {
        let original = "NC_000001.11:g.100_101insATG";
        let hgvs = parse_hgvs(original).unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(back.to_string(), original);
    }

    #[test]
    fn test_roundtrip_deletion_with_seq() {
        let original = "NC_000001.11:g.100_102delATG";
        let hgvs = parse_hgvs(original).unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(back.to_string(), original);
    }

    #[test]
    fn test_error_display() {
        let err = ConversionError::UnsupportedVariantType {
            description: "test".to_string(),
        };
        assert!(err.to_string().contains("unsupported variant type"));

        let err = ConversionError::MissingReferenceData {
            description: "test".to_string(),
        };
        assert!(err.to_string().contains("missing reference data"));
    }

    // =========================================================================
    // Issue #117: Reference-aware HGVS→SPDI for del/dup/delins
    // =========================================================================
    //
    // These tests exercise the production `hgvs_to_spdi(variant, provider)`
    // path with a real `MockProvider`. They replace the earlier
    // `MockGenomicRef` prototype, which predated the public provider entry
    // point and is now redundant.

    /// Build a `MockProvider` with a contig where:
    ///   1-based 100..102 = "ATG"
    ///   1-based 200..206 = "GATTACA"
    ///   1-based 1000..1009 = "AAACCCGGGT"
    /// All other positions are filled with 'N'.
    fn make_test_genomic_provider() -> crate::reference::mock::MockProvider {
        let mut p = crate::reference::mock::MockProvider::new();
        let mut contig = String::new();
        contig.push_str(&"N".repeat(99)); // 1-based 1..99 (0-based 0..99)
        contig.push_str("ATG"); // 1-based 100..102
        contig.push_str(&"N".repeat(97)); // pad through 1-based 199
        contig.push_str("GATTACA"); // 1-based 200..206
        contig.push_str(&"N".repeat(793)); // pad through 1-based 999
        contig.push_str("AAACCCGGGT"); // 1-based 1000..1009
        contig.push_str(&"N".repeat(50));
        p.add_genomic_sequence("NC_000001.11", &contig);
        p
    }

    #[test]
    fn fetch_reference_bases_returns_genomic_bases() {
        let provider = make_test_genomic_provider();
        let bases = fetch_reference_bases(&provider, "NC_000001.11", 100, 102).unwrap();
        assert_eq!(bases, "ATG");
    }

    #[test]
    fn fetch_reference_bases_errors_when_provider_lacks_contig() {
        let provider = crate::reference::mock::MockProvider::new();
        let err = fetch_reference_bases(&provider, "NC_000099.99", 100, 102).unwrap_err();
        assert!(matches!(err, ConversionError::MissingReferenceData { .. }));
        let msg = err.to_string();
        assert!(msg.contains("NC_000099.99"));
        assert!(msg.contains("100"));
        assert!(msg.contains("102"));
    }

    #[test]
    fn fetch_reference_bases_errors_on_short_contig() {
        let mut provider = crate::reference::mock::MockProvider::new();
        // Contig is only 3 bases — fetching 1-based 100..102 must fail.
        provider.add_genomic_sequence("NC_000001.11", "ATG");
        let err = fetch_reference_bases(&provider, "NC_000001.11", 100, 102).unwrap_err();
        assert!(matches!(err, ConversionError::MissingReferenceData { .. }));
    }

    #[test]
    fn hgvs_to_spdi_deletion_short_form_with_provider() {
        let provider = make_test_genomic_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102del").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn hgvs_to_spdi_duplication_short_form_with_provider() {
        let provider = make_test_genomic_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102dup").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        // Dup encodes as an SPDI insertion at the 3' end of the duplicated
        // region. SPDI interbase: the position is the 1-based end of
        // the dup region (102), which is the boundary AFTER base 102
        // and matches the equivalent `g.102_103ins…` form (#390).
        assert_eq!(spdi.position, 102);
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn hgvs_to_spdi_single_base_duplication_short_form_with_provider() {
        let provider = make_test_genomic_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100dup").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        // Single-base dup: end_one_based = 100, SPDI position = 100 (#390).
        assert_eq!(spdi.position, 100);
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "A");
    }

    #[test]
    fn hgvs_to_spdi_delins_short_form_with_provider() {
        let provider = make_test_genomic_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102delinsTTCC").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "TTCC");
    }

    #[test]
    fn hgvs_to_spdi_long_deletion_with_provider() {
        let provider = make_test_genomic_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.1000_1009del").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.deletion, "AAACCCGGGT");
        assert_eq!(spdi.deletion.len(), 10);
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn hgvs_to_spdi_explicit_deletion_does_not_consult_provider() {
        // When the user supplied an explicit deleted sequence, ferro emits
        // it as-is. Verified by attaching an empty provider — if we
        // consulted it, the call would fail.
        let provider = crate::reference::mock::MockProvider::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102delATG").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn hgvs_to_spdi_explicit_duplication_does_not_consult_provider() {
        let provider = crate::reference::mock::MockProvider::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102dupATG").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        // SPDI interbase position 102 = boundary AFTER 1-based 102 (#390).
        assert_eq!(spdi.position, 102);
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn hgvs_to_spdi_substitution_unaffected_by_provider() {
        let provider = make_test_genomic_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.to_string(), "NC_000001.11:12344:A:G");
    }

    #[test]
    fn hgvs_to_spdi_mnv_delins_with_provider_round_trips() {
        // Same-length delins should round-trip through SPDI.
        let provider = make_test_genomic_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102delinsGGG").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "GGG");
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(back.to_string(), "NC_000001.11:g.100_102delinsGGG");
    }

    #[test]
    fn hgvs_to_spdi_deletion_round_trip_with_provider() {
        let provider = make_test_genomic_provider();
        let original = parse_hgvs("NC_000001.11:g.100_102del").unwrap();
        let spdi = hgvs_to_spdi(&original, &provider).unwrap();
        let back = spdi_to_hgvs(&spdi).unwrap();
        // SPDI carries the deleted sequence, so the recovered HGVS is the
        // explicit form.
        assert_eq!(back.to_string(), "NC_000001.11:g.100_102delATG");
    }

    #[test]
    fn hgvs_to_spdi_delins_round_trip_with_provider() {
        let provider = make_test_genomic_provider();
        let original = parse_hgvs("NC_000001.11:g.100_102delinsTTCC").unwrap();
        let spdi = hgvs_to_spdi(&original, &provider).unwrap();
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(back.to_string(), "NC_000001.11:g.100_102delinsTTCC");
    }

    #[test]
    fn hgvs_to_spdi_dup_round_trip_emits_ins_form_via_reference_free_path() {
        // HGVS dup → SPDI ins, with the duplicated bases populated from
        // the provider. Without a reference-aware SPDI→HGVS direction the
        // reverse path produces ins form (per the SPDI-as-canonical
        // contract). PR #119 (sibling) adds `spdi_to_hgvs_with_ref` that
        // recovers the dup form when the 5' flank matches.
        let provider = make_test_genomic_provider();
        let original = parse_hgvs("NC_000001.11:g.100_102dup").unwrap();
        let spdi = hgvs_to_spdi(&original, &provider).unwrap();
        // Post-#390: SPDI position = end_one_based (102), matching the
        // equivalent `g.102_103ins` interbase boundary.
        assert_eq!(spdi.position, 102);
        assert_eq!(spdi.insertion, "ATG");
        let recovered = spdi_to_hgvs(&spdi).unwrap();
        // SPDI 102 = boundary AFTER 1-based 102 = HGVS g.102_103ins.
        assert_eq!(recovered.to_string(), "NC_000001.11:g.102_103insATG");
    }

    #[test]
    fn hgvs_to_spdi_mito_short_form_deletion_with_provider() {
        // Mito accession (NC_012920.1) is genomic; verify the same path
        // works there.
        let mut provider = crate::reference::mock::MockProvider::new();
        let mut seq = "N".repeat(16559);
        seq.push_str("GATC"); // 1-based 16560..16563
        seq.push_str(&"N".repeat(20));
        provider.add_genomic_sequence("NC_012920.1", &seq);
        let hgvs = parse_hgvs("NC_012920.1:m.16560_16563del").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.deletion, "GATC");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn hgvs_to_spdi_deletion_with_provider_missing_data() {
        // Provider attached but has no data for the requested contig.
        let provider = crate::reference::mock::MockProvider::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102del").unwrap();
        let result = hgvs_to_spdi(&hgvs, &provider);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("NC_000001.11"));
    }

    #[test]
    fn hgvs_to_spdi_simple_pins_existing_short_form_failures() {
        // The simple (no-provider) entry point continues to err on the
        // three short-form cases, matching the existing audit pin.
        for input in [
            "NC_000001.11:g.100_102del",
            "NC_000001.11:g.100_102dup",
            "NC_000001.11:g.100_102delinsTTCC",
        ] {
            let hgvs = parse_hgvs(input).unwrap();
            let r = hgvs_to_spdi_simple(&hgvs);
            assert!(
                matches!(r, Err(ConversionError::MissingReferenceData { .. })),
                "expected MissingReferenceData for {} (got {:?})",
                input,
                r
            );
        }
    }

    #[test]
    fn hgvs_to_spdi_short_form_is_idempotent() {
        let provider = make_test_genomic_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102del").unwrap();
        let a = hgvs_to_spdi(&hgvs, &provider).unwrap();
        let b = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(a.to_string(), b.to_string());
    }

    // =========================================================================
    // P3: SPDI edge case tests
    // =========================================================================

    #[test]
    fn test_spdi_empty_deletion_insertion() {
        // Pure insertion has empty deletion
        let spdi = SpdiVariant::insertion("NC_000001.11", 100, "ATG");
        assert!(spdi.is_insertion());
        assert!(!spdi.is_deletion());
        assert!(!spdi.is_identity());
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn test_spdi_empty_insertion_deletion() {
        // Pure deletion has empty insertion
        let spdi = SpdiVariant::deletion("NC_000001.11", 100, "ATG");
        assert!(spdi.is_deletion());
        assert!(!spdi.is_insertion());
        assert!(!spdi.is_identity());
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn test_spdi_both_empty_is_identity() {
        // Both empty is identity at position (unusual but valid)
        let spdi = SpdiVariant::new("NC_000001.11", 100, "", "");
        assert!(spdi.is_identity());
        assert!(!spdi.is_insertion());
        assert!(!spdi.is_deletion());
    }

    #[test]
    fn test_spdi_single_base_insertion() {
        let spdi = SpdiVariant::insertion("NC_000001.11", 100, "A");
        assert!(spdi.is_insertion());
        assert_eq!(spdi.insertion.len(), 1);

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("ins"));
    }

    #[test]
    fn test_spdi_single_base_deletion() {
        let spdi = SpdiVariant::deletion("NC_000001.11", 100, "A");
        assert!(spdi.is_deletion());
        assert_eq!(spdi.deletion.len(), 1);

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("del"));
    }

    #[test]
    fn test_spdi_long_insertion_100bp() {
        // Test 100bp insertion (common structural variant size)
        let long_seq = "A".repeat(100);
        let spdi = SpdiVariant::insertion("NC_000001.11", 12345, &long_seq);

        assert!(spdi.is_insertion());
        assert_eq!(spdi.insertion.len(), 100);

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("ins"));
        // Verify the sequence is preserved
        assert!(hgvs.to_string().ends_with(&format!("ins{}", long_seq)));
    }

    #[test]
    fn test_spdi_long_deletion_100bp() {
        // Test 100bp deletion
        let long_seq = "ACGT".repeat(25); // 100bp
        let spdi = SpdiVariant::deletion("NC_000001.11", 12345, &long_seq);

        assert!(spdi.is_deletion());
        assert_eq!(spdi.deletion.len(), 100);

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("del"));
    }

    #[test]
    fn test_spdi_long_indel_asymmetric() {
        // Delete 50bp, insert 150bp (net +100bp)
        let del_seq = "A".repeat(50);
        let ins_seq = "G".repeat(150);
        let spdi = SpdiVariant::delins("NC_000001.11", 12345, &del_seq, &ins_seq);

        assert!(!spdi.is_insertion());
        assert!(!spdi.is_deletion());
        assert_eq!(spdi.deletion.len(), 50);
        assert_eq!(spdi.insertion.len(), 150);

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("delins"));
    }

    #[test]
    fn test_spdi_very_long_insertion_1000bp() {
        // Test 1000bp insertion (larger structural variant)
        let long_seq = "ACGT".repeat(250); // 1000bp
        let spdi = SpdiVariant::insertion("NC_000001.11", 50000, &long_seq);

        assert!(spdi.is_insertion());
        assert_eq!(spdi.insertion.len(), 1000);

        // Should still convert to HGVS
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("ins"));
    }

    #[test]
    fn test_spdi_position_zero() {
        // Position 0 is valid in SPDI (0-based)
        let spdi = SpdiVariant::new("NC_000001.11", 0, "A", "G");
        assert_eq!(spdi.position, 0);

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        // HGVS position should be 1 (1-based)
        assert!(hgvs.to_string().contains("g.1A>G"));
    }

    #[test]
    fn test_spdi_position_max() {
        // Very large position (near chromosome end)
        let spdi = SpdiVariant::new("NC_000001.11", 248956421, "A", "G");

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("248956422")); // 1-based
    }

    #[test]
    fn test_spdi_lowercase_sequence() {
        // SPDI should handle lowercase (though uppercase is standard)
        let spdi = SpdiVariant::new("NC_000001.11", 100, "a", "g");

        // The variant should work even with lowercase
        assert_eq!(spdi.deletion, "a");
        assert_eq!(spdi.insertion, "g");
    }

    #[test]
    fn test_spdi_mixed_case_sequence() {
        // Mixed case sequence
        let spdi = SpdiVariant::new("NC_000001.11", 100, "AtGc", "GcTa");

        assert_eq!(spdi.deletion, "AtGc");
        assert_eq!(spdi.insertion, "GcTa");
    }

    #[test]
    fn test_spdi_n_bases_in_sequence() {
        // N (unknown) bases in sequence
        let spdi = SpdiVariant::new("NC_000001.11", 100, "ANG", "TNC");

        assert_eq!(spdi.deletion, "ANG");
        assert_eq!(spdi.insertion, "TNC");
    }

    #[test]
    fn test_spdi_complex_repeat_sequence() {
        // Repeat sequence (e.g., microsatellite)
        let repeat = "CAG".repeat(30); // 90bp CAG repeat
        let spdi = SpdiVariant::insertion("NC_000004.12", 3074876, &repeat);

        assert!(spdi.is_insertion());
        assert_eq!(spdi.insertion.len(), 90);
    }

    #[test]
    fn test_spdi_to_hgvs_delins_single_base_del() {
        // Single base deletion with multi-base insertion
        let spdi = SpdiVariant::delins("NC_000001.11", 100, "A", "TTTT");

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        // Single base del + ins = delins at single position
        assert!(hgvs.to_string().contains("delinsTTTT"));
    }

    #[test]
    fn test_spdi_to_hgvs_delins_single_base_ins() {
        // Multi-base deletion with single base insertion
        let spdi = SpdiVariant::delins("NC_000001.11", 100, "AAAA", "T");

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("delinsT"));
    }

    #[test]
    fn test_spdi_different_chromosome_formats() {
        // Various accession formats should work
        let test_cases = vec![
            ("NC_000001.11", "NC_000001.11"), // Standard RefSeq
            ("NC_000023.11", "NC_000023.11"), // X chromosome
            ("NC_000024.10", "NC_000024.10"), // Y chromosome
            ("NC_012920.1", "NC_012920.1"),   // Mitochondrial
        ];

        for (input_acc, expected_acc) in test_cases {
            let spdi = SpdiVariant::new(input_acc, 100, "A", "G");
            assert_eq!(spdi.sequence, expected_acc);

            let hgvs = spdi_to_hgvs(&spdi).unwrap();
            assert!(hgvs.to_string().starts_with(expected_acc));
        }
    }

    #[test]
    fn test_spdi_roundtrip_preserves_case_normalized() {
        // Create SPDI substitution (SNV) which roundtrips cleanly
        let spdi = SpdiVariant::new("NC_000001.11", 100, "A", "G");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        let back = hgvs_to_spdi_simple(&hgvs).unwrap();

        // Should preserve the sequence
        assert_eq!(back.deletion, "A");
        assert_eq!(back.insertion, "G");

        // Test multi-base delins - without reference data, the roundtrip
        // cannot reconstruct the deletion sequence and returns an error
        let spdi_delins = SpdiVariant::new("NC_000001.11", 100, "ACGT", "TGCA");
        let hgvs_delins = spdi_to_hgvs(&spdi_delins).unwrap();
        let back_delins = hgvs_to_spdi_simple(&hgvs_delins);

        // Without reference, this should fail since the deletion sequence is unknown
        assert!(back_delins.is_err());
    }

    #[test]
    fn test_spdi_empty_seq_insertion_roundtrip() {
        // Empty deletion (pure insertion) roundtrip
        let spdi = SpdiVariant::insertion("NC_000001.11", 100, "ATG");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();

        // HGVS format for insertion is pos_pos+1insX
        assert!(hgvs.to_string().contains("ins"));

        // Roundtrip back
        let back = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(back.deletion, "");
        assert_eq!(back.insertion, "ATG");
    }

    #[test]
    fn test_hgvs_to_spdi_various_accession_types() {
        // Test conversion from various accession types
        let test_variants = vec![
            "NC_000001.11:g.12345A>G", // Chromosome
            "NC_000023.11:g.12345A>G", // X chromosome
            "NC_012920.1:g.12345A>G",  // Mitochondrial
        ];

        for variant_str in test_variants {
            let hgvs = parse_hgvs(variant_str).unwrap();
            let spdi = hgvs_to_spdi_simple(&hgvs);
            assert!(spdi.is_ok(), "Failed for: {}", variant_str);
        }
    }

    #[test]
    fn test_spdi_display_format() {
        let spdi = SpdiVariant::new("NC_000001.11", 12344, "A", "G");
        assert_eq!(spdi.to_string(), "NC_000001.11:12344:A:G");

        let spdi_del = SpdiVariant::deletion("NC_000001.11", 100, "ATG");
        assert_eq!(spdi_del.to_string(), "NC_000001.11:100:ATG:");

        let spdi_ins = SpdiVariant::insertion("NC_000001.11", 100, "ATG");
        assert_eq!(spdi_ins.to_string(), "NC_000001.11:100::ATG");
    }

    #[test]
    fn test_spdi_identity_various_lengths() {
        // Single base identity
        let spdi1 = SpdiVariant::new("NC_000001.11", 100, "A", "A");
        assert!(spdi1.is_identity());

        // Multi-base identity (unusual but valid)
        let spdi2 = SpdiVariant::new("NC_000001.11", 100, "ATG", "ATG");
        assert!(spdi2.is_identity());

        // Empty identity
        let spdi3 = SpdiVariant::new("NC_000001.11", 100, "", "");
        assert!(spdi3.is_identity());
    }

    // =========================================================================
    // Issue #119: SPDI→HGVS dup recovery
    // =========================================================================

    /// Build a MockProvider with a single genomic sequence registered for
    /// `NC_000001.11`. The string `seq` is the full contig sequence, indexed
    /// from 0-based position 0.
    fn provider_with_genomic(seq: &str) -> crate::reference::mock::MockProvider {
        let mut p = crate::reference::mock::MockProvider::new();
        p.add_genomic_sequence("NC_000001.11", seq);
        p
    }

    #[test]
    fn spdi_to_hgvs_with_ref_recovers_multi_base_dup() {
        // Reference: positions 1-based 100..102 = "ATG"
        // Build a contig where 0-based offsets 99, 100, 101 = 'A', 'T', 'G'
        // We pad with 'N' before and after so the reference fetch works.
        let mut contig = "N".repeat(99);
        contig.push_str("ATG"); // 0-based 99..102 (1-based 100..102)
        contig.push_str(&"N".repeat(50));
        let provider = provider_with_genomic(&contig);

        // SPDI 102::ATG is the canonical SPDI form of g.100_102dupATG
        // under the interbase-correct convention (#390): the insertion
        // sits at the boundary AFTER 1-based base 102.
        let spdi = SpdiVariant::insertion("NC_000001.11", 102, "ATG");

        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100_102dupATG");
    }

    #[test]
    fn spdi_to_hgvs_with_ref_recovers_single_base_dup() {
        // 1-based base 100 = 'A' → SPDI insertion at 100::A is a single-base dup.
        let mut contig = "N".repeat(99);
        contig.push('A'); // 0-based 99 = 1-based 100
        contig.push_str(&"N".repeat(20));
        let provider = provider_with_genomic(&contig);

        // g.100dupA → SPDI 100::A under the interbase-correct
        // convention (#390); SPDI position 100 sits AFTER 1-based 100,
        // and the 5' flank base 1-based 100 = 'A' matches the
        // inserted 'A'.
        let spdi = SpdiVariant::insertion("NC_000001.11", 100, "A");

        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100dupA");
    }

    #[test]
    fn spdi_to_hgvs_with_ref_keeps_ins_when_no_match() {
        // 5' flank of length 3 ending at SPDI position 102 (under the
        // post-#390 interbase convention) = 1-based bases 100..102 =
        // "CCC" — does NOT equal the inserted "ATG".
        let mut contig = "N".repeat(99);
        contig.push_str("CCC"); // 1-based 100..102 = "CCC"
        contig.push_str(&"N".repeat(20));
        let provider = provider_with_genomic(&contig);

        let spdi = SpdiVariant::insertion("NC_000001.11", 102, "ATG");

        // Should fall through to ins-form. spdi_to_hgvs renders an
        // insertion as `g.{pos}_{pos+1}insATG` where `pos` is the
        // SPDI 0-based interbase position itself (#390): SPDI 102 →
        // g.102_103insATG.
        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.102_103insATG");
    }

    #[test]
    fn spdi_to_hgvs_with_ref_rejects_ins_at_contig_start() {
        // SPDI position 0 with a 3-base insertion: there are no
        // preceding bases at all, and HGVS has no standard notation
        // for "insert before the first base" (#390). The conversion
        // surfaces `InvalidPosition` rather than silently emitting
        // `g.1_2insATG` (which mis-represents the actual insertion
        // point).
        let provider = provider_with_genomic("ATGCATGCATGC");

        let spdi = SpdiVariant::insertion("NC_000001.11", 0, "ATG");
        let err = spdi_to_hgvs_with_ref(&spdi, &provider).expect_err(
            "SPDI position 0 (insertion before contig start) must surface InvalidPosition",
        );
        assert!(matches!(err, ConversionError::InvalidPosition { .. }));
    }

    #[test]
    fn spdi_to_hgvs_with_ref_substitution_unchanged() {
        let provider = provider_with_genomic(&"N".repeat(20000));
        let spdi = SpdiVariant::new("NC_000001.11", 12344, "A", "G");
        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.12345A>G");
    }

    #[test]
    fn spdi_to_hgvs_with_ref_deletion_unchanged() {
        let provider = provider_with_genomic(&"N".repeat(2000));
        let spdi = SpdiVariant::deletion("NC_000001.11", 99, "ATG");
        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100_102delATG");
    }

    #[test]
    fn spdi_to_hgvs_with_ref_delins_unchanged() {
        let provider = provider_with_genomic(&"N".repeat(2000));
        let spdi = SpdiVariant::delins("NC_000001.11", 99, "ATG", "TTCC");
        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100_102delinsTTCC");
    }

    #[test]
    fn spdi_to_hgvs_with_ref_identity_unchanged() {
        let provider = provider_with_genomic(&"N".repeat(2000));
        let spdi = SpdiVariant::new("NC_000001.11", 99, "A", "A");
        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100A=");
    }

    /// A zero-width triple is an identity by the predicate (`deletion ==
    /// insertion`) but names no bases, so the conversion arm must refuse it
    /// rather than emit `g.100=` over a base the triple never claimed.
    ///
    /// Asserted on **both** public entry points: `spdi_to_hgvs_with_ref`
    /// delegates to `spdi_to_hgvs` for its base conversion, so the refusal
    /// reaches it too — and a future cut that stopped delegating would silently
    /// lose the guard on the reference-aware path.
    #[test]
    fn zero_width_identity_is_refused_on_both_entry_points() {
        let spdi = SpdiVariant::new("NC_000001.11", 99, "", "");
        assert!(spdi.is_identity(), "both sides empty is an identity");
        assert!(
            matches!(
                spdi_to_hgvs(&spdi),
                Err(ConversionError::UnsupportedEditType { .. })
            ),
            "a triple naming zero bases must not convert to an identity"
        );
        let provider = provider_with_genomic(&"N".repeat(2000));
        assert!(
            matches!(
                spdi_to_hgvs_with_ref(&spdi, &provider),
                Err(ConversionError::UnsupportedEditType { .. })
            ),
            "the reference-aware path must inherit the refusal"
        );
    }

    /// A one-base identity still converts — the guard must be scoped to the
    /// zero-width case and not have cost the arm its ordinary behaviour.
    #[test]
    fn a_one_base_identity_still_converts() {
        let spdi = SpdiVariant::new("NC_000001.11", 99, "A", "A");
        let hgvs = spdi_to_hgvs(&spdi).expect("a one-base identity converts");
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100A=");
    }

    #[test]
    fn spdi_to_hgvs_with_ref_propagates_ref_error() {
        // MockProvider with NO genomic_sequences registered for the
        // requested accession. The fetch returns
        // FerroError::GenomicReferenceNotAvailable for get_genomic_sequence
        // AND no transcript by that name for get_sequence. Both fail, so
        // the helper should return ConversionError::MissingReferenceData.
        let provider = crate::reference::mock::MockProvider::new();
        let spdi = SpdiVariant::insertion("NC_000999.99", 101, "ATG");

        let result = spdi_to_hgvs_with_ref(&spdi, &provider);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    // =========================================================================
    // Issue #116: c./n./r./m. coordinate-system support
    // =========================================================================

    use crate::reference::mock::MockProvider;
    use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand};

    /// Build a small test transcript covering enough cases for c./n./r. tests:
    /// - 5'UTR length 5 (positions 1-5 in tx)
    /// - CDS  length 30 (tx positions 6-35; cds_start=6, cds_end=35)
    /// - 3'UTR length 5 (tx positions 36-40)
    ///
    /// Single exon over the full transcript so the simple no-gap path applies.
    fn make_test_provider() -> MockProvider {
        let tx = Transcript::new(
            "NM_TEST.1".to_string(),
            Some("TEST".to_string()),
            Strand::Plus,
            // 40 bases: 5 (UTR5) + 30 (CDS) + 5 (UTR3)
            "AAAAATGCCCAAAGGGTTTAGGCCCAAAGGGTTATAAA".to_string() + "AA",
            Some(6),
            Some(35),
            vec![Exon::new(1, 1, 40)],
            None,
            None,
            None,
            GenomeBuild::default(),
            ManeStatus::default(),
            None,
            None,
        );
        let mut provider = MockProvider::new();
        provider.add_transcript(tx);
        provider
    }

    /// Build a multi-exon transcript suitable for testing intronic resolution
    /// rejections. Exon 1: tx 1-50, Exon 2: tx 51-100; CDS tx 11-90.
    fn make_intronic_provider() -> MockProvider {
        let tx = Transcript::new(
            "NM_INTRON.1".to_string(),
            Some("INTRON".to_string()),
            Strand::Plus,
            "A".repeat(100),
            Some(11),
            Some(90),
            vec![Exon::new(1, 1, 50), Exon::new(2, 51, 100)],
            None,
            None,
            None,
            GenomeBuild::default(),
            ManeStatus::default(),
            None,
            None,
        );
        let mut provider = MockProvider::new();
        provider.add_transcript(tx);
        provider
    }

    // ----- r.*N / c.*N 3'UTR agreement (#944) --------------------------------

    /// #944: r.*N and c.*N denote the same 3'UTR transcript position, so both
    /// must resolve to the same SPDI position. `resolve_rna_pos` now routes the
    /// 3'UTR case through the same `CoordinateMapper::cds_to_tx` the c.*N path
    /// uses (previously it short-circuited with its own `cds_end + base`), so
    /// one conversion answers both spellings.
    #[test]
    fn r_star_and_c_star_resolve_to_same_spdi_position() {
        let provider = make_test_provider();
        // NM_TEST.1: cds_end = tx 35, so *1 is tx 36 (1-based) → SPDI 0-based 35.
        let c = parse_hgvs("NM_TEST.1:c.*1del").unwrap();
        let r = parse_hgvs("NM_TEST.1:r.*1del").unwrap();
        let c_spdi = hgvs_to_spdi(&c, &provider).unwrap();
        let r_spdi = hgvs_to_spdi(&r, &provider).unwrap();
        assert_eq!(
            c_spdi.position, r_spdi.position,
            "c.*1 and r.*1 must resolve to the same SPDI position"
        );
        assert_eq!(c_spdi.position, 35);
    }

    /// Build a transcript whose 3'UTR lives in a *separate exon across a
    /// tx-coordinate gap*, so the `r.` and `c.` arms only agree if both take
    /// the same route:
    ///
    /// - Exon 1: tx 1..35 (5'UTR tx 1-5, CDS tx 6-35; `cds_end` = 35 is the
    ///   last base of exon 1).
    /// - Gap in the exon table: tx 36..39 are covered by no exon (cdot-style
    ///   alignment gap — `exon1.end + 1 (36) != exon2.start (40)`). The stored
    ///   sequence is 44 bases, so those four bases DO exist on the transcript;
    ///   what they lack is a genomic counterpart.
    /// - Exon 2: tx 40..44.
    ///
    /// **#1619 decided which of those bases `*1` names, and it is tx 36.** The
    /// `c.`/`n.`/`r.` axis counts the transcript's own bases
    /// (`background/numbering.md:21`, `:52`), so `c.*1` is `cds_end + 1` on the
    /// flat sequence. This fixture used to pin the opposite — the exon walk,
    /// which skipped the four unaligned bases and answered tx 40 — and the doc
    /// here called tx 36 "a nonexistent coordinate", which is the mistaken
    /// premise the ruling `c-and-n-positions-are-flat-transcript-offsets`
    /// overturns.
    ///
    /// The fixture stays, because the property it was built for is still the
    /// live one and is still non-vacuous here: `r.*N` and `c.*N` must resolve
    /// identically (#944), and on a single-exon transcript the two routes
    /// coincide so any test there would pass either way.
    fn make_gapped_utr3_provider() -> MockProvider {
        let tx = Transcript::new(
            "NM_GAP.1".to_string(),
            Some("GAP".to_string()),
            Strand::Plus,
            "A".repeat(44),
            Some(6),
            Some(35),
            vec![Exon::new(1, 1, 35), Exon::new(2, 40, 44)],
            None,
            None,
            None,
            GenomeBuild::default(),
            ManeStatus::default(),
            None,
            None,
        );
        let mut provider = MockProvider::new();
        provider.add_transcript(tx);
        provider
    }

    /// #944 (non-vacuous): on a transcript whose 3'UTR sits in a separate exon
    /// across a tx-coordinate gap, r.*1 and c.*1 must STILL resolve to the same
    /// SPDI position — the property #944 is about, and one a single-exon
    /// fixture cannot test because both routes coincide there.
    ///
    /// #1619 fixes WHICH position that is: `*1` is `cds_end + 1` on the flat
    /// transcript (tx 36 → SPDI 35), not the exon walk's tx 40.
    #[test]
    fn r_star_c_star_agree_across_exon_gap() {
        let provider = make_gapped_utr3_provider();
        // cds_end = tx 35, so *1 is tx 36 → SPDI 0-based 35. tx 36-39 are
        // transcript bases the exon table does not cover; `c.*N` counts them.
        let c = parse_hgvs("NM_GAP.1:c.*1A>G").unwrap();
        let r = parse_hgvs("NM_GAP.1:r.*1a>g").unwrap();
        let c_spdi = hgvs_to_spdi(&c, &provider).unwrap();
        let r_spdi = hgvs_to_spdi(&r, &provider).unwrap();
        assert_eq!(
            c_spdi.position, r_spdi.position,
            "c.*1 and r.*1 must resolve to the same SPDI position across the exon gap"
        );
        assert_eq!(
            c_spdi.position, 35,
            "*1 is cds_end + 1 on the flat transcript (tx 36 → SPDI 35); the exon \
             walk this replaced skipped the tx 36-39 bases and answered SPDI 39 (#1619)"
        );
    }

    /// A 3'UTR r.*N on a non-coding transcript (no CDS end) still declines
    /// cleanly rather than mapping through a missing anchor.
    #[test]
    fn r_star_on_non_coding_transcript_declines() {
        let mut provider = MockProvider::new();
        let tx = Transcript::new(
            "NR_TEST.1".to_string(),
            Some("NCTEST".to_string()),
            Strand::Plus,
            "A".repeat(40),
            None,
            None,
            vec![Exon::new(1, 1, 40)],
            None,
            None,
            None,
            GenomeBuild::default(),
            ManeStatus::default(),
            None,
            None,
        );
        provider.add_transcript(tx);
        let r = parse_hgvs("NR_TEST.1:r.*1del").unwrap();
        assert!(hgvs_to_spdi(&r, &provider).is_err());
    }

    // ----- r.-N / c.-N 5'UTR agreement (#960) --------------------------------

    /// #960: r.-N and c.-N denote the same 5'UTR transcript position, so both
    /// must resolve to the same SPDI position. `resolve_rna_pos` now routes the
    /// 5'UTR case through the same `CoordinateMapper::cds_to_tx` the c.-N path
    /// uses (previously it fell through to `ensure_positive_tx` and was
    /// rejected as a non-positive base).
    #[test]
    fn r_minus_and_c_minus_resolve_to_same_spdi_position() {
        let provider = make_test_provider();
        // NM_TEST.1: cds_start = tx 6, so c.-3 / r.-3 is tx 3 (1-based) → SPDI 2.
        let c = parse_hgvs("NM_TEST.1:c.-3del").unwrap();
        let r = parse_hgvs("NM_TEST.1:r.-3del").unwrap();
        let c_spdi = hgvs_to_spdi(&c, &provider).unwrap();
        let r_spdi = hgvs_to_spdi(&r, &provider).unwrap();
        assert_eq!(
            c_spdi.position, r_spdi.position,
            "c.-3 and r.-3 must resolve to the same SPDI position"
        );
        assert_eq!(c_spdi.position, 2);
    }

    /// Build a transcript whose 5'UTR straddles a *tx-coordinate gap*, so the
    /// `r.` and `c.` arms only agree if both take the same route:
    ///
    /// - Exon 1: tx 1..5 — the upstream part of the 5'UTR.
    /// - Gap in the exon table: tx 6..9 are covered by no exon (`exon1.end + 1
    ///   (6) != exon2.start (10)`). The stored sequence is 44 bases, so those
    ///   four are transcript bases with no genomic counterpart.
    /// - Exon 2: tx 10..44 — 5'UTR tx 10-11 then CDS from `cds_start` = tx 12.
    ///
    /// **#1619 decided which base `c.-3` names, and it is tx 9** — three bases
    /// 5' of `cds_start` on the flat transcript. This fixture used to pin the
    /// exon walk's tx 5 and called tx 9 "inside the nonexistent gap"; the
    /// ruling `c-and-n-positions-are-flat-transcript-offsets` overturns that
    /// premise. The fixture stays because the #960 property — `r.-N` and `c.-N`
    /// resolving identically — is still live and still untestable on a
    /// single-exon transcript, where the two routes coincide.
    fn make_gapped_utr5_provider() -> MockProvider {
        let tx = Transcript::new(
            "NM_GAP5.1".to_string(),
            Some("GAP5".to_string()),
            Strand::Plus,
            "A".repeat(44),
            Some(12),
            Some(40),
            vec![Exon::new(1, 1, 5), Exon::new(2, 10, 44)],
            None,
            None,
            None,
            GenomeBuild::default(),
            ManeStatus::default(),
            None,
            None,
        );
        let mut provider = MockProvider::new();
        provider.add_transcript(tx);
        provider
    }

    /// #960 (non-vacuous): on a transcript whose 5'UTR straddles a tx-coordinate
    /// gap, r.-3 and c.-3 must STILL resolve to the same SPDI position — the
    /// property #960 is about, and one a single-exon fixture cannot test.
    ///
    /// #1619 fixes WHICH position that is: `c.-3` is `cds_start - 3` on the
    /// flat transcript (tx 9 → SPDI 8), not the exon walk's tx 5.
    #[test]
    fn r_minus_c_minus_agree_across_exon_gap() {
        let provider = make_gapped_utr5_provider();
        // cds_start = tx 12, so c.-3 / r.-3 is tx 9 → SPDI 0-based 8. tx 6-9 are
        // transcript bases the exon table does not cover; `c.-N` counts them.
        let c = parse_hgvs("NM_GAP5.1:c.-3A>G").unwrap();
        let r = parse_hgvs("NM_GAP5.1:r.-3a>g").unwrap();
        let c_spdi = hgvs_to_spdi(&c, &provider).unwrap();
        let r_spdi = hgvs_to_spdi(&r, &provider).unwrap();
        assert_eq!(
            c_spdi.position, r_spdi.position,
            "c.-3 and r.-3 must resolve to the same SPDI position across the exon gap"
        );
        assert_eq!(
            c_spdi.position, 8,
            "-3 is cds_start - 3 on the flat transcript (tx 9 → SPDI 8); the exon \
             walk this replaced skipped the tx 6-9 bases and answered SPDI 4 (#1619)"
        );
    }

    /// A 5'UTR r.-N on a non-coding transcript (no CDS start) declines cleanly
    /// rather than mapping through a missing anchor.
    #[test]
    fn r_minus_on_non_coding_transcript_declines() {
        let mut provider = MockProvider::new();
        let tx = Transcript::new(
            "NR_TEST5.1".to_string(),
            Some("NCTEST5".to_string()),
            Strand::Plus,
            "A".repeat(40),
            None,
            None,
            vec![Exon::new(1, 1, 40)],
            None,
            None,
            None,
            GenomeBuild::default(),
            ManeStatus::default(),
            None,
            None,
        );
        provider.add_transcript(tx);
        let r = parse_hgvs("NR_TEST5.1:r.-3del").unwrap();
        assert!(hgvs_to_spdi(&r, &provider).is_err());
    }

    // ----- r.*N / c.*N 3' upper bound (#962) ---------------------------------

    /// #962 boundary: the last 3'UTR base still resolves. NM_TEST.1 has tx length
    /// 40 with the 3'UTR at tx 36-40, so `*5` = tx 40 (the final base) — c.*5 and
    /// r.*5 must both succeed and agree at SPDI 0-based 39. Guards the upper-bound
    /// check against off-by-one over-rejection of the last valid position.
    #[test]
    fn r_star_and_c_star_at_last_base_resolve() {
        let provider = make_test_provider();
        let c = parse_hgvs("NM_TEST.1:c.*5del").unwrap();
        let r = parse_hgvs("NM_TEST.1:r.*5del").unwrap();
        let c_spdi = hgvs_to_spdi(&c, &provider).unwrap();
        let r_spdi = hgvs_to_spdi(&r, &provider).unwrap();
        assert_eq!(
            c_spdi.position, r_spdi.position,
            "c.*5 and r.*5 (last transcript base) must resolve to the same SPDI position"
        );
        assert_eq!(c_spdi.position, 39);
    }

    /// #962: a `*N` position past the transcript 3' end is declined (not emitted
    /// as an off-sequence SPDI coordinate), for both c.*N and r.*N, and they
    /// agree in declining. On NM_TEST.1 (tx length 40, 3'UTR ends at *5 = tx 40),
    /// `*6` = tx 41 is one base past the end. Before #962 the resolver only guarded
    /// the lower bound, so `*6` mapped to a nonexistent tx 41 → SPDI 40.
    #[test]
    fn r_star_and_c_star_past_end_both_decline() {
        let provider = make_test_provider();
        let c = parse_hgvs("NM_TEST.1:c.*6del").unwrap();
        let r = parse_hgvs("NM_TEST.1:r.*6del").unwrap();
        assert!(
            hgvs_to_spdi(&c, &provider).is_err(),
            "c.*6 is one base past the transcript 3' end and must decline"
        );
        assert!(
            hgvs_to_spdi(&r, &provider).is_err(),
            "r.*6 is one base past the transcript 3' end and must decline"
        );
        // Far past the end declines just as cleanly.
        let c_far = parse_hgvs("NM_TEST.1:c.*99999del").unwrap();
        let r_far = parse_hgvs("NM_TEST.1:r.*99999del").unwrap();
        assert!(hgvs_to_spdi(&c_far, &provider).is_err());
        assert!(hgvs_to_spdi(&r_far, &provider).is_err());
    }

    // ----- m. (mitochondrial) ------------------------------------------------

    #[test]
    fn test_hgvs_to_spdi_simple_mt_substitution() {
        // m. uses NC_012920.1, which is itself a genomic accession.
        let hgvs = parse_hgvs("NC_012920.1:m.3243A>G").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.sequence, "NC_012920.1");
        assert_eq!(spdi.position, 3242);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "G");
        assert_eq!(spdi.to_string(), "NC_012920.1:3242:A:G");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_mt_insertion() {
        let hgvs = parse_hgvs("NC_012920.1:m.100_101insATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        // SPDI interbase 100 = boundary AFTER 1-based 100 (#390).
        assert_eq!(spdi.position, 100);
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_mt_deletion_with_seq() {
        let hgvs = parse_hgvs("NC_012920.1:m.3243_3245delAGG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.position, 3242);
        assert_eq!(spdi.deletion, "AGG");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_mt_deletion_without_seq_needs_ref() {
        // Without an explicit deleted sequence, the simple path can't supply
        // the SPDI deletion field — same as g.
        let hgvs = parse_hgvs("NC_012920.1:m.3243_3245del").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    #[test]
    fn dup_hgvs_to_spdi_to_hgvs_with_ref_roundtrip_multi_base() {
        // Build a contig with 1-based 100..102 = "ATG"
        let mut contig = "N".repeat(99);
        contig.push_str("ATG");
        contig.push_str(&"N".repeat(20));
        let provider = provider_with_genomic(&contig);

        // Forward: HGVS dup → SPDI ins (interbase position 102 #390)
        let original = "NC_000001.11:g.100_102dupATG";
        let hgvs = parse_hgvs(original).unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.to_string(), "NC_000001.11:102::ATG");

        // Reverse with reference: SPDI ins → HGVS dup
        let recovered = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        assert_eq!(recovered.to_string(), original);
    }

    #[test]
    fn dup_hgvs_to_spdi_to_hgvs_with_ref_roundtrip_single_base() {
        let mut contig = "N".repeat(99);
        contig.push('A');
        contig.push_str(&"N".repeat(20));
        let provider = provider_with_genomic(&contig);

        let original = "NC_000001.11:g.100dupA";
        let hgvs = parse_hgvs(original).unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        // SPDI interbase 100 (#390): boundary AFTER 1-based 100.
        assert_eq!(spdi.to_string(), "NC_000001.11:100::A");

        let recovered = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        assert_eq!(recovered.to_string(), original);
    }

    #[test]
    fn spdi_to_hgvs_with_ref_does_not_false_detect_non_tandem_insertion() {
        // Spec FAQ: ATCGATCGATCG-A-GGGTCCC → ATCGATCGATCG-A-ATCGATCGATCG-GGGTCCC.
        // The 12-base ATCGATCGATCG sequence appears in the reference at
        // 1-based 1..12, but the insertion point (between 1-based 13 and 14,
        // i.e., SPDI position 13) has a 5' flank "TCGATCGATCGA" — NOT
        // matching the inserted "ATCGATCGATCG". Per the spec FAQ this MUST
        // remain ins.
        let contig = "ATCGATCGATCGAGGGTCCC".to_string();
        let provider = provider_with_genomic(&contig);

        // SPDI position 13 = 0-based 13 = inserts between 1-based 13 and 14.
        let spdi = SpdiVariant::insertion("NC_000001.11", 13, "ATCGATCGATCG");
        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        // Expect ins, not dup. Confirm the rendered string contains "ins"
        // and not "dup".
        let s = hgvs.to_string();
        assert!(s.contains("ins"), "expected ins-form, got {}", s);
        assert!(!s.contains("dup"), "expected not dup, got {}", s);
    }

    #[test]
    fn audit_pin_no_ref_spdi_to_hgvs_renders_dup_shape_as_ins() {
        // Pins issue #119 documented behavior: without a reference,
        // spdi_to_hgvs cannot prove tandem dup, so the dup-shaped SPDI
        // 102::ATG (which round-tripped from g.100_102dupATG under the
        // post-#390 interbase-correct convention) is rendered as
        // g.102_103insATG. If a future change attempts to "fix" the
        // round-trip without a reference, this audit pin will fail and
        // demand explicit reconsideration.
        let spdi = SpdiVariant::insertion("NC_000001.11", 102, "ATG");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.102_103insATG");
    }

    #[test]
    fn dup_recovery_is_idempotent_through_two_roundtrips() {
        let mut contig = "N".repeat(99);
        contig.push_str("ATG");
        contig.push_str(&"N".repeat(20));
        let provider = provider_with_genomic(&contig);

        let original = "NC_000001.11:g.100_102dupATG";
        let hgvs1 = parse_hgvs(original).unwrap();
        let spdi1 = hgvs_to_spdi_simple(&hgvs1).unwrap();
        let hgvs2 = spdi_to_hgvs_with_ref(&spdi1, &provider).unwrap();
        let spdi2 = hgvs_to_spdi_simple(&hgvs2).unwrap();
        let hgvs3 = spdi_to_hgvs_with_ref(&spdi2, &provider).unwrap();

        assert_eq!(spdi1, spdi2);
        assert_eq!(hgvs2.to_string(), hgvs3.to_string());
        assert_eq!(hgvs3.to_string(), original);
    }

    #[test]
    fn spdi_to_hgvs_with_ref_recovers_dup_for_mito_accession() {
        // Mock contig for NC_012920.1 with 1-based bases 100..102 = "ATG".
        let mut contig = "N".repeat(99);
        contig.push_str("ATG");
        contig.push_str(&"N".repeat(20));

        let mut provider = crate::reference::mock::MockProvider::new();
        provider.add_genomic_sequence("NC_012920.1", &contig);

        // SPDI 102 (interbase, #390) sits AFTER 1-based 102; the
        // 5' flank "ATG" matches the inserted "ATG" so dup recovers.
        let spdi = SpdiVariant::insertion("NC_012920.1", 102, "ATG");
        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        // After #270, SPDI→HGVS preserves the m. coord system for known
        // mitochondrial accessions (NC_012920.*); dup recovery applies
        // uniformly to genomic and mitochondrial variants.
        assert_eq!(hgvs.to_string(), "NC_012920.1:m.100_102dupATG");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_mt_dup_with_seq() {
        let hgvs = parse_hgvs("NC_012920.1:m.100_102dupATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        // SPDI interbase 102 (#390): boundary AFTER 1-based 102.
        assert_eq!(spdi.position, 102);
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_mt_identity() {
        let hgvs = parse_hgvs("NC_012920.1:m.3243A=").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.position, 3242);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "A");
    }

    // ----- n. (non-coding transcript) ----------------------------------------

    #[test]
    fn test_hgvs_to_spdi_simple_tx_substitution() {
        // n. — SPDI emitted on the transcript accession.
        let hgvs = parse_hgvs("NR_046018.2:n.5C>G").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.sequence, "NR_046018.2");
        assert_eq!(spdi.position, 4); // 5 (1-based) → 4 (0-based)
        assert_eq!(spdi.deletion, "C");
        assert_eq!(spdi.insertion, "G");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_tx_insertion() {
        let hgvs = parse_hgvs("NR_046018.2:n.10_11insATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.sequence, "NR_046018.2");
        // SPDI interbase 10 (#390): boundary AFTER 1-based 10.
        assert_eq!(spdi.position, 10);
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_tx_deletion_with_seq() {
        let hgvs = parse_hgvs("NR_046018.2:n.10_12delATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.position, 9);
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_tx_identity() {
        let hgvs = parse_hgvs("NR_046018.2:n.10A=").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.position, 9);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "A");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_tx_intronic_needs_provider() {
        // n.100+5: intronic offset cannot be expressed as a positional SPDI;
        // the simple path bails with MissingReferenceData (provider needed
        // for genomic projection — out of scope for this PR).
        let hgvs = parse_hgvs("NR_046018.2:n.100+5A>G").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("intronic"), "msg: {}", msg);
    }

    #[test]
    fn test_hgvs_to_spdi_simple_tx_downstream_needs_provider() {
        // n.*5: downstream of transcript end; needs transcript length.
        // Constructed rather than parsed — see `tx_downstream_sub` (#1748).
        let hgvs = tx_downstream_sub("NR_046018.2", 5);
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    #[test]
    fn test_hgvs_to_spdi_simple_tx_negative_base_needs_provider() {
        // n.-3: upstream of transcript start; the simple path has no way to
        // anchor it without knowing the transcript length.
        let hgvs = parse_hgvs("NR_046018.2:n.-3A>G").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    // ----- r. (RNA) ----------------------------------------------------------

    #[test]
    fn test_hgvs_to_spdi_simple_dna_lowercase_uppercased() {
        // Lowercase DNA bases (e.g. `g.100a>g`) must emit uppercase SPDI
        // alleles per the doc on `apply_alphabet`. Regression test: the DNA
        // branch previously preserved input case.
        let hgvs = parse_hgvs("NC_000001.11:g.100a>g").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "G");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_rna_substitution_lowercase() {
        // r.5c>g (lowercase RNA) → SPDI uses uppercase DNA alphabet.
        let hgvs = parse_hgvs("NR_046018.2:r.5c>g").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.sequence, "NR_046018.2");
        assert_eq!(spdi.position, 4);
        assert_eq!(spdi.deletion, "C");
        assert_eq!(spdi.insertion, "G");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_rna_substitution_u_to_t() {
        // r.5u>g: U on the deleted side rewrites to T for SPDI.
        let hgvs = parse_hgvs("NR_046018.2:r.5u>g").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.deletion, "T");
        assert_eq!(spdi.insertion, "G");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_rna_insertion_u_to_t() {
        // r.10_11insauug → SPDI insertion ATTG (a→A, u→T, u→T, g→G).
        let hgvs = parse_hgvs("NR_046018.2:r.10_11insauug").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        // SPDI interbase 10 (#390): boundary AFTER 1-based 10.
        assert_eq!(spdi.position, 10);
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ATTG");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_rna_deletion_with_seq() {
        let hgvs = parse_hgvs("NR_046018.2:r.10_12delauu").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.position, 9);
        assert_eq!(spdi.deletion, "ATT"); // a→A, u→T, u→T
    }

    #[test]
    fn test_hgvs_to_spdi_simple_rna_intronic_needs_provider() {
        let hgvs = parse_hgvs("NR_046018.2:r.10+5a>g").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    // ----- p. (protein) — rejected with helpful error -----------------------

    #[test]
    fn test_hgvs_to_spdi_simple_protein_rejected() {
        let hgvs = parse_hgvs("NP_000079.2:p.Arg600Gln").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::UnsupportedVariantType { .. })
        ));
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("protein") && msg.contains("SPDI"),
            "expected helpful protein-rejection message; got: {}",
            msg
        );
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_protein_rejected() {
        let provider = MockProvider::new();
        let hgvs = parse_hgvs("NP_000079.2:p.Arg600Gln").unwrap();
        let result = hgvs_to_spdi(&hgvs, &provider);
        assert!(matches!(
            result,
            Err(ConversionError::UnsupportedVariantType { .. })
        ));
    }

    // ----- c. (CDS) — provider-aware ----------------------------------------

    #[test]
    fn test_hgvs_to_spdi_with_provider_cds_substitution() {
        let provider = make_test_provider();
        // c.1A>G: cds_start=6, so c.1 → tx 6 → SPDI 5
        let hgvs = parse_hgvs("NM_TEST.1:c.1A>G").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.sequence, "NM_TEST.1");
        assert_eq!(spdi.position, 5);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "G");
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_cds_insertion() {
        let provider = make_test_provider();
        // c.1_2insATG: cds_start=6 → tx start_one_based 6 → SPDI
        // interbase position 6 (boundary AFTER 1-based 6; #390).
        let hgvs = parse_hgvs("NM_TEST.1:c.1_2insATG").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.position, 6);
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_cds_deletion_with_seq() {
        let provider = make_test_provider();
        // c.1_3delATG: tx 6_8 → SPDI 5:ATG:
        let hgvs = parse_hgvs("NM_TEST.1:c.1_3delATG").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.position, 5);
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_cds_5utr() {
        let provider = make_test_provider();
        // c.-3A>G: 3 bases before cds_start (6) → tx 3 → SPDI 2
        let hgvs = parse_hgvs("NM_TEST.1:c.-3A>G").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.position, 2);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "G");
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_cds_3utr() {
        let provider = make_test_provider();
        // c.*2A>G: 2 bases past cds_end (35) → tx 37 → SPDI 36
        let hgvs = parse_hgvs("NM_TEST.1:c.*2A>G").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.position, 36);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "G");
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_cds_exonic_past_3prime_end_declines() {
        // NM_TEST.1 is 40 bases (cds 6-35). c.99999 -> tx 100004, off-sequence.
        // #971: must decline, not emit an off-sequence SPDI coordinate.
        let provider = make_test_provider();
        let hgvs = parse_hgvs("NM_TEST.1:c.99999A>G").unwrap();
        let result = hgvs_to_spdi(&hgvs, &provider);
        assert!(
            matches!(result, Err(ConversionError::InvalidPosition { .. })),
            "over-length exonic c.N must decline with InvalidPosition, got {result:?}"
        );
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_cds_exonic_range_past_3prime_end_declines() {
        // In a range, an over-length end must also decline (#971 acceptance: "in a range").
        let provider = make_test_provider();
        let hgvs = parse_hgvs("NM_TEST.1:c.20_99999del").unwrap();
        let result = hgvs_to_spdi(&hgvs, &provider);
        assert!(
            matches!(result, Err(ConversionError::InvalidPosition { .. })),
            "over-length exonic c. range must decline, got {result:?}"
        );
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_cds_exonic_last_base_ok() {
        // c.30 -> tx 35 (last CDS base, in bounds) must still resolve unchanged.
        let provider = make_test_provider();
        let hgvs = parse_hgvs("NM_TEST.1:c.30A>G").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.position, 34, "c.30 -> tx 35 -> SPDI 34");
    }

    /// r.*N anchors at `cds_end`, not `sequence_length()` — closes
    /// #390 item 2. With the test fixture (cds_end=35, tx_len=40),
    /// r.*2 must resolve to tx position 37 (SPDI 36), matching the
    /// equivalent `c.*2` resolution.
    #[test]
    fn test_hgvs_to_spdi_with_provider_rna_3utr_anchors_at_cds_end() {
        let provider = make_test_provider();
        // For SPDI emission r. lowercase is normalized to upper / U→T
        // by the conversion path; r.*2a>g exercises the 3'UTR anchor.
        let hgvs = parse_hgvs("NM_TEST.1:r.*2a>g").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(
            spdi.position, 36,
            "r.*2 must anchor at cds_end (35) + 2 → tx 37 → SPDI 36; pre-#390 \
             this anchored at sequence_length (40) + 2 → off-sequence SPDI 41"
        );
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "G");
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_rna_exonic_coding_matches_cds() {
        // On a CODING transcript r.N is CDS-relative (#469): r.N and c.N denote the
        // same base, so they must resolve to the same SPDI position. NM_TEST.1
        // cds_start=6, so r.1 == c.1 -> tx 6 -> SPDI 5; r.30 == c.30 -> tx 35 -> SPDI 34.
        let provider = make_test_provider();
        for (r_str, c_str) in [
            ("NM_TEST.1:r.1a>g", "NM_TEST.1:c.1A>G"),
            ("NM_TEST.1:r.30a>g", "NM_TEST.1:c.30A>G"),
        ] {
            let r = hgvs_to_spdi(&parse_hgvs(r_str).unwrap(), &provider).unwrap();
            let c = hgvs_to_spdi(&parse_hgvs(c_str).unwrap(), &provider).unwrap();
            assert_eq!(r.position, c.position, "{r_str} must match {c_str}");
        }
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_rna_exonic_coding_first_base() {
        // Regression on the exact confirmed values: r.1a>g -> SPDI 5 (was 0).
        let provider = make_test_provider();
        let spdi = hgvs_to_spdi(&parse_hgvs("NM_TEST.1:r.1a>g").unwrap(), &provider).unwrap();
        assert_eq!(spdi.position, 5, "coding r.1 -> tx 6 (cds_start) -> SPDI 5");
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_rna_exonic_coding_past_3prime_end_declines() {
        // Now that coding r.N routes through cds_to_tx (#469/#971), an over-length
        // coding r.N must still decline rather than emit an off-sequence SPDI
        // position: cds_to_tx returns an out-of-range base, which
        // ensure_tx_in_bounds then rejects — same as over-length coding c.N.
        let provider = make_test_provider();
        let hgvs = parse_hgvs("NM_TEST.1:r.99999a>g").unwrap();
        let result = hgvs_to_spdi(&hgvs, &provider);
        assert!(
            matches!(result, Err(ConversionError::InvalidPosition { .. })),
            "over-length exonic coding r.N must decline with InvalidPosition, got {result:?}"
        );
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_rna_exonic_coding_range_matches_cds() {
        // The coding r.N == c.N CDS-relative identity must hold across BOTH
        // endpoints of a range, not just single bases: r.1_30del and c.1_30del
        // must produce the same SPDI (position and deleted span).
        let provider = make_test_provider();
        let r = hgvs_to_spdi(&parse_hgvs("NM_TEST.1:r.1_30del").unwrap(), &provider).unwrap();
        let c = hgvs_to_spdi(&parse_hgvs("NM_TEST.1:c.1_30del").unwrap(), &provider).unwrap();
        assert_eq!(
            r.position, c.position,
            "r.1_30del must match c.1_30del position"
        );
        assert_eq!(
            r.deletion.len(),
            c.deletion.len(),
            "r.1_30del must delete the same span as c.1_30del"
        );
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_rna_exonic_coding_missing_cds_end_declines() {
        // A coding transcript (cds_start set) with NO cds_end is malformed. A
        // coding r.N cannot be resolved CDS-relative without a CDS end, so it
        // must DECLINE — not silently fall back to an unbounded,
        // transcript-absolute r.N (off by cds_start-1). The decline surfaces as
        // InvalidPosition, which propagates past the resolve_rna_exonic_bounded
        // MissingReferenceData fallback.
        let tx = Transcript::new(
            "NM_NOCDSEND.1".to_string(),
            Some("TEST".to_string()),
            Strand::Plus,
            "AAAAATGCCCAAAGGGTTTAGGCCCAAAGGGTTATAAA".to_string() + "AA",
            Some(6),
            None,
            vec![Exon::new(1, 1, 40)],
            None,
            None,
            None,
            GenomeBuild::default(),
            ManeStatus::default(),
            None,
            None,
        );
        let mut provider = MockProvider::new();
        provider.add_transcript(tx);
        let result = hgvs_to_spdi(&parse_hgvs("NM_NOCDSEND.1:r.5a>g").unwrap(), &provider);
        assert!(
            matches!(result, Err(ConversionError::InvalidPosition { .. })),
            "coding r.N on a transcript missing cds_end must decline with \
             InvalidPosition, got {result:?}"
        );
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_cds_intronic_rejected() {
        // Intronic c. cannot be expressed as a positional SPDI: SPDI has no
        // offset notation. The provider-aware path surfaces a clear error.
        let provider = make_intronic_provider();
        let hgvs = parse_hgvs("NM_INTRON.1:c.10+5A>G").unwrap();
        let result = hgvs_to_spdi(&hgvs, &provider);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_unknown_transcript() {
        let provider = MockProvider::new();
        let hgvs = parse_hgvs("NM_TEST.1:c.1A>G").unwrap();
        let result = hgvs_to_spdi(&hgvs, &provider);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("transcript") || msg.contains("NM_TEST"));
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_falls_through_to_simple_for_genome() {
        let provider = MockProvider::new();
        let hgvs = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.to_string(), "NC_000001.11:12344:A:G");
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_falls_through_to_simple_for_mt() {
        let provider = MockProvider::new();
        let hgvs = parse_hgvs("NC_012920.1:m.3243A>G").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.to_string(), "NC_012920.1:3242:A:G");
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_falls_through_to_simple_for_exonic_n() {
        // Exonic positive-base n. doesn't actually need provider data; the
        // provider-aware path delegates to the simple path.
        let provider = MockProvider::new();
        let hgvs = parse_hgvs("NR_046018.2:n.5C>G").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.to_string(), "NR_046018.2:4:C:G");
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_tx_exonic_past_3prime_end_declines() {
        // n. treats the base as a direct transcript position; NM_TEST.1 is 40 bases.
        // n.99999 is off-sequence and must decline (#971).
        let provider = make_test_provider();
        let hgvs = parse_hgvs("NM_TEST.1:n.99999C>G").unwrap();
        let result = hgvs_to_spdi(&hgvs, &provider);
        assert!(
            matches!(result, Err(ConversionError::InvalidPosition { .. })),
            "over-length exonic n.N must decline with InvalidPosition, got {result:?}"
        );
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_tx_exonic_range_past_3prime_end_declines() {
        // Over-length end in an n. range must also decline (#971 acceptance: "in a range").
        let provider = make_test_provider();
        let hgvs = parse_hgvs("NM_TEST.1:n.10_99999del").unwrap();
        let result = hgvs_to_spdi(&hgvs, &provider);
        assert!(
            matches!(result, Err(ConversionError::InvalidPosition { .. })),
            "over-length exonic n. range must decline, got {result:?}"
        );
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_tx_exonic_last_base_ok() {
        // n.40 -> tx 40 (last base, in bounds) still resolves.
        let provider = make_test_provider();
        let hgvs = parse_hgvs("NM_TEST.1:n.40C>G").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.position, 39, "n.40 -> tx 40 -> SPDI 39");
    }

    /// Single-exon non-coding transcript: 40 bases, no CDS. Exonic r.N maps the
    /// base directly to the transcript position (NR_ has no CDS anchor), so it is
    /// the clean fixture for r.N 3'-bound tests (avoids the coding-r.N #469 case).
    fn make_noncoding_provider() -> MockProvider {
        let tx = Transcript::new(
            "NR_TEST.1".to_string(),
            Some("NCTEST".to_string()),
            Strand::Plus,
            "A".repeat(40),
            None,
            None,
            vec![Exon::new(1, 1, 40)],
            None,
            None,
            None,
            GenomeBuild::default(),
            ManeStatus::default(),
            None,
            None,
        );
        let mut provider = MockProvider::new();
        provider.add_transcript(tx);
        provider
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_rna_exonic_past_3prime_end_declines() {
        let provider = make_noncoding_provider();
        let hgvs = parse_hgvs("NR_TEST.1:r.99999a>g").unwrap();
        let result = hgvs_to_spdi(&hgvs, &provider);
        assert!(
            matches!(result, Err(ConversionError::InvalidPosition { .. })),
            "over-length exonic r.N must decline with InvalidPosition, got {result:?}"
        );
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_rna_exonic_range_past_3prime_end_declines() {
        let provider = make_noncoding_provider();
        let hgvs = parse_hgvs("NR_TEST.1:r.10_99999del").unwrap();
        let result = hgvs_to_spdi(&hgvs, &provider);
        assert!(
            matches!(result, Err(ConversionError::InvalidPosition { .. })),
            "over-length exonic r. range must decline, got {result:?}"
        );
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_rna_exonic_last_base_ok() {
        // r.40 -> tx 40 (last base, in bounds) still resolves.
        let provider = make_noncoding_provider();
        let hgvs = parse_hgvs("NR_TEST.1:r.40a>g").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.position, 39, "r.40 -> tx 40 -> SPDI 39");
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_falls_through_to_simple_for_exonic_r() {
        // Best-effort (#971): provider lacks the transcript, so the exonic r.
        // conversion falls back to the unbounded simple-path value rather than
        // erroring — mirrors the n. sibling test.
        let provider = MockProvider::new();
        let hgvs = parse_hgvs("NR_046018.2:r.5c>g").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.to_string(), "NR_046018.2:4:C:G");
    }

    #[test]
    fn test_hgvs_to_spdi_with_provider_n_downstream_rejected() {
        // n. has no CDS anchor, so `n.*N` is past the transcript end and
        // cannot be expressed in SPDI on the transcript accession. The
        // provider-aware path must reject rather than silently emit an
        // off-sequence position at `tx_len + N`.
        let tx = Transcript::new(
            "NR_NONCODING.1".to_string(),
            Some("NONCODING".to_string()),
            Strand::Plus,
            "A".repeat(40),
            None,
            None,
            vec![Exon::new(1, 1, 40)],
            None,
            None,
            None,
            GenomeBuild::default(),
            ManeStatus::default(),
            None,
            None,
        );
        let mut provider = MockProvider::new();
        provider.add_transcript(tx);
        // Constructed rather than parsed — see `tx_downstream_sub` (#1748).
        let hgvs = tx_downstream_sub("NR_NONCODING.1", 5);
        let err = hgvs_to_spdi(&hgvs, &provider).unwrap_err();
        assert!(matches!(err, ConversionError::InvalidPosition { .. }));
        let msg = err.to_string();
        assert!(
            msg.contains("downstream n.") && msg.contains("genomic projection"),
            "expected downstream-n rejection, got: {}",
            msg
        );
    }

    // ----- Edit-type passthrough on new coord systems -----------------------

    #[test]
    fn test_hgvs_to_spdi_simple_n_short_form_inversion_requires_provider() {
        // Short-form `n.` inversion needs the provider to fetch reference
        // bases — same shape as `g.`.
        let hgvs = parse_hgvs("NR_046018.2:n.10_20inv").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    #[test]
    fn test_hgvs_to_spdi_simple_m_short_form_inversion_requires_provider() {
        // Short-form `m.` inversion needs the provider to fetch reference
        // bases — same shape as `g.`.
        let hgvs = parse_hgvs("NC_012920.1:m.100_200inv").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    // ---------------------------------------------------------------------
    // Inversion: provider-aware emission (#118)
    // ---------------------------------------------------------------------

    #[test]
    fn test_hgvs_to_spdi_inversion_short_form_with_provider() {
        // 1-based 100..102 is "ATG" → revcomp "CAT"
        let provider = make_test_genomic_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102inv").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "CAT");
    }

    #[test]
    fn test_hgvs_to_spdi_inversion_explicit_sequence_no_provider() {
        // Explicit sequence bypasses the provider entirely — matches the
        // explicit-form policy already used for del/dup.
        let hgvs = parse_hgvs("NC_000001.11:g.100_102invATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "CAT");
    }

    #[test]
    fn test_hgvs_to_spdi_inversion_explicit_sequence_does_not_consult_provider() {
        // Same as above but routed through the provider-aware path with an
        // empty provider — explicit-form must not call the provider.
        let provider = crate::reference::mock::MockProvider::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102invATG").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "CAT");
    }

    #[test]
    fn test_hgvs_to_spdi_inversion_single_base() {
        // Inversion of a single base 'A' → 'T'. `g.100_100inv` cannot be
        // parsed — DNA/inversion.md:16 forbids a one-nucleotide inversion —
        // so the variant is built directly to keep the conversion arm
        // covered for callers that construct an AST themselves.
        let provider = make_test_genomic_provider();
        let hgvs = HgvsVariant::Genome(GenomeVariant {
            accession: parse_accession("NC_000001.11").unwrap().1,
            gene_symbol: None,
            loc_edit: LocEdit::new(
                crate::hgvs::interval::GenomeInterval::point(GenomePos::new(100)),
                NaEdit::Inversion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "T");
    }

    #[test]
    fn test_hgvs_to_spdi_inversion_palindrome_round_trip() {
        // ATAT is its own reverse complement; del == ins.
        let mut provider = crate::reference::mock::MockProvider::new();
        let mut contig = "N".repeat(99);
        contig.push_str("ATAT");
        contig.push_str(&"N".repeat(50));
        provider.add_genomic_sequence("NC_000001.11", &contig);
        let hgvs = parse_hgvs("NC_000001.11:g.100_103inv").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.deletion, "ATAT");
        assert_eq!(spdi.insertion, "ATAT");
    }

    #[test]
    fn test_hgvs_to_spdi_inversion_short_form_missing_provider_data() {
        // Provider has no contig; expect a MissingReferenceData with the
        // accession + 1-based interval in the message (same shape as #117).
        let provider = crate::reference::mock::MockProvider::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102inv").unwrap();
        let err = hgvs_to_spdi(&hgvs, &provider).unwrap_err();
        assert!(matches!(err, ConversionError::MissingReferenceData { .. }));
        let msg = err.to_string();
        assert!(msg.contains("NC_000001.11"));
        assert!(msg.contains("100"));
        assert!(msg.contains("102"));
    }

    #[test]
    fn test_hgvs_to_spdi_inversion_m_short_form_with_provider() {
        // `m.` works the same as `g.` — accession is the mito contig.
        let mut provider = crate::reference::mock::MockProvider::new();
        let mut contig = "N".repeat(99);
        contig.push_str("ATG");
        contig.push_str(&"N".repeat(50));
        provider.add_genomic_sequence("NC_012920.1", &contig);
        let hgvs = parse_hgvs("NC_012920.1:m.100_102inv").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.sequence, "NC_012920.1");
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "CAT");
    }

    #[test]
    fn test_hgvs_to_spdi_inversion_n_short_form_with_provider() {
        // `n.` (non-coding tx): SPDI emits on the transcript accession.
        let mut provider = crate::reference::mock::MockProvider::new();
        let mut contig = "N".repeat(9);
        contig.push_str("ATGCATGC"); // 1-based 10..17
        contig.push_str(&"N".repeat(50));
        provider.add_genomic_sequence("NR_046018.2", &contig);
        let hgvs = parse_hgvs("NR_046018.2:n.10_12inv").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.sequence, "NR_046018.2");
        assert_eq!(spdi.position, 9);
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "CAT");
    }

    #[test]
    fn test_hgvs_to_spdi_inversion_r_short_form_dna_alphabet() {
        // `r.` input — SPDI must come out in DNA alphabet (T, not U).
        // Reference bases on transcripts are stored as DNA, so the
        // del/ins emitted here are uppercase DNA.
        let mut provider = crate::reference::mock::MockProvider::new();
        let mut contig = "N".repeat(9);
        contig.push_str("ATGCATGC");
        contig.push_str(&"N".repeat(50));
        provider.add_genomic_sequence("NR_046018.2", &contig);
        let hgvs = parse_hgvs("NR_046018.2:r.10_12inv").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "CAT");
        // The output must be DNA (no U).
        assert!(!spdi.deletion.contains('U'));
        assert!(!spdi.insertion.contains('U'));
    }

    // ---------------------------------------------------------------------
    // Repeat: provider-aware emission (#118)
    // ---------------------------------------------------------------------

    /// Build a provider with an AT-tandem repeat tract:
    ///   1-based 100..105 = "ATATAT" (3 copies of AT)
    ///   1-based 200..209 = "ATATATATAT" (5 copies of AT)
    fn make_repeat_provider() -> crate::reference::mock::MockProvider {
        let mut p = crate::reference::mock::MockProvider::new();
        let mut contig = "N".repeat(99);
        contig.push_str("ATATAT"); // 1-based 100..105 (6 bases)
        contig.push_str(&"N".repeat(94)); // pad through 1-based 199
        contig.push_str("ATATATATAT"); // 1-based 200..209 (10 bases)
        contig.push_str(&"N".repeat(50));
        p.add_genomic_sequence("NC_000001.11", &contig);
        p
    }

    #[test]
    fn test_hgvs_to_spdi_repeat_expansion_with_provider() {
        // Reference has 3 copies of AT (100..105 = ATATAT); allele has 5.
        let provider = make_repeat_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_105AT[5]").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "ATATAT"); // 3 copies in reference
        assert_eq!(spdi.insertion, "ATATATATAT"); // 5 copies on allele
    }

    #[test]
    fn test_hgvs_to_spdi_repeat_contraction_with_provider() {
        // Reference has 5 copies (200..209); allele has 3.
        let provider = make_repeat_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.200_209AT[3]").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.position, 199);
        assert_eq!(spdi.deletion, "ATATATATAT"); // 5 copies
        assert_eq!(spdi.insertion, "ATATAT"); // 3 copies
    }

    #[test]
    fn test_hgvs_to_spdi_repeat_no_change_with_provider() {
        // Reference and allele both have 3 copies → del == ins (a valid
        // SPDI identity-shape delins).
        let provider = make_repeat_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_105AT[3]").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.deletion, "ATATAT");
        assert_eq!(spdi.insertion, "ATATAT");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_repeat_requires_provider() {
        // No provider → MissingReferenceData (replaces prior
        // UnsupportedEditType behaviour).
        let hgvs = parse_hgvs("NC_000001.11:g.100_105AT[5]").unwrap();
        let err = hgvs_to_spdi_simple(&hgvs).unwrap_err();
        assert!(matches!(err, ConversionError::MissingReferenceData { .. }));
    }

    #[test]
    fn test_hgvs_to_spdi_repeat_missing_provider_data() {
        let provider = crate::reference::mock::MockProvider::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_105AT[5]").unwrap();
        let err = hgvs_to_spdi(&hgvs, &provider).unwrap_err();
        assert!(matches!(err, ConversionError::MissingReferenceData { .. }));
        let msg = err.to_string();
        assert!(msg.contains("NC_000001.11"));
        assert!(msg.contains("100"));
        assert!(msg.contains("105"));
    }

    #[test]
    fn test_hgvs_to_spdi_repeat_uncertain_count_unsupported() {
        // RepeatCount::Range — uncertain — cannot be a single SPDI value.
        let provider = make_repeat_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_105AT[3_5]").unwrap();
        let err = hgvs_to_spdi(&hgvs, &provider).unwrap_err();
        assert!(matches!(err, ConversionError::UnsupportedEditType { .. }));
    }

    #[test]
    fn test_hgvs_to_spdi_repeat_unknown_count_unsupported() {
        let provider = make_repeat_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_105AT[?]").unwrap();
        let err = hgvs_to_spdi(&hgvs, &provider).unwrap_err();
        assert!(matches!(err, ConversionError::UnsupportedEditType { .. }));
    }

    #[test]
    fn test_hgvs_to_spdi_repeat_genotype_unsupported() {
        // Genotype-style additional counts cannot be a single SPDI.
        let provider = make_repeat_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_105AT[3][5]").unwrap();
        let err = hgvs_to_spdi(&hgvs, &provider).unwrap_err();
        assert!(matches!(err, ConversionError::UnsupportedEditType { .. }));
    }

    #[test]
    fn test_hgvs_to_spdi_repeat_no_unit_unsupported() {
        // `g.100_105(5)` parenthesized form: NaEdit::Repeat with no unit.
        // Without a unit we cannot construct the inserted sequence.
        let provider = make_repeat_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_105(5)").unwrap();
        let err = hgvs_to_spdi(&hgvs, &provider).unwrap_err();
        assert!(matches!(err, ConversionError::MissingReferenceData { .. }));
    }

    #[test]
    fn test_hgvs_to_spdi_repeat_span_not_multiple_of_unit() {
        // Span 100..104 is 5 bases; AT unit is 2 bases → not divisible.
        let provider = make_repeat_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_104AT[5]").unwrap();
        let err = hgvs_to_spdi(&hgvs, &provider).unwrap_err();
        assert!(matches!(err, ConversionError::InvalidPosition { .. }));
        let msg = err.to_string();
        assert!(msg.contains("not a multiple"));
    }

    #[test]
    fn test_hgvs_to_spdi_repeat_span_does_not_match_unit() {
        // Reference span is `ATGCAT` (6 bases, multiple of unit length 2)
        // but the locus is NOT an AT tandem repeat. Must reject rather
        // than silently emit a wrong SPDI delins.
        let mut provider = crate::reference::mock::MockProvider::new();
        let mut contig = "N".repeat(99);
        contig.push_str("ATGCAT"); // 1-based 100..105 — not a tandem AT[3]
        contig.push_str(&"N".repeat(50));
        provider.add_genomic_sequence("NC_000001.11", &contig);

        let hgvs = parse_hgvs("NC_000001.11:g.100_105AT[5]").unwrap();
        let err = hgvs_to_spdi(&hgvs, &provider).unwrap_err();
        assert!(matches!(err, ConversionError::InvalidPosition { .. }));
        let msg = err.to_string();
        assert!(
            msg.contains("does not match repeat unit"),
            "expected mismatch message, got: {}",
            msg
        );
    }

    #[test]
    fn test_hgvs_to_spdi_repeat_n_with_provider() {
        // `n.` repeat — SPDI on transcript accession.
        let mut provider = crate::reference::mock::MockProvider::new();
        let mut contig = "N".repeat(9);
        contig.push_str("ATATAT"); // 1-based 10..15
        contig.push_str(&"N".repeat(50));
        provider.add_genomic_sequence("NR_046018.2", &contig);
        let hgvs = parse_hgvs("NR_046018.2:n.10_15AT[5]").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.sequence, "NR_046018.2");
        assert_eq!(spdi.position, 9);
        assert_eq!(spdi.deletion, "ATATAT");
        assert_eq!(spdi.insertion, "ATATATATAT");
    }

    #[test]
    fn test_hgvs_to_spdi_repeat_expansion_too_large() {
        // A huge user-supplied count must be refused before allocating the
        // expanded ins-string. The reference span itself is small (and a
        // valid AT tandem) so this isolates the size guard.
        let provider = make_repeat_provider();
        let huge = MAX_REPEAT_EXPANSION_BASES / 2 + 1; // unit.len() == 2
        let hgvs = parse_hgvs(&format!("NC_000001.11:g.100_105AT[{}]", huge)).unwrap();
        let err = hgvs_to_spdi(&hgvs, &provider).unwrap_err();
        assert!(matches!(err, ConversionError::UnsupportedEditType { .. }));
        let msg = err.to_string();
        assert!(
            msg.contains("exceeds SPDI ins-string cap"),
            "expected size-cap message, got: {}",
            msg
        );
    }

    #[test]
    fn test_hgvs_to_spdi_repeat_m_with_provider() {
        let mut provider = crate::reference::mock::MockProvider::new();
        let mut contig = "N".repeat(99);
        contig.push_str("ATATAT"); // 1-based 100..105
        contig.push_str(&"N".repeat(50));
        provider.add_genomic_sequence("NC_012920.1", &contig);
        let hgvs = parse_hgvs("NC_012920.1:m.100_105AT[5]").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.sequence, "NC_012920.1");
        assert_eq!(spdi.deletion, "ATATAT");
        assert_eq!(spdi.insertion, "ATATATATAT");
    }

    // ------------------------------------------------------------------
    // Genomic offsets (#1628)
    // ------------------------------------------------------------------

    /// A `g.` position carrying a `+`/`-` offset must not be silently
    /// flattened onto its base. `checklist.md:16` prohibits an offset on a
    /// genomic position and SPDI has no offset notation, so the only two
    /// honest answers are "refuse" and "resolve"; there is nothing to
    /// resolve against on a bare genomic accession, so the answer is refuse.
    ///
    /// Dropping the offset is the third, dishonest answer: it makes
    /// `g.266+2del`, `g.266-268del` and `g.266del` — which `normalize`
    /// treats as three different variants — collapse onto one triple.
    ///
    /// All three genomic axes are exercised, on both public entry points, so
    /// each of the six guard call sites is covered: `m.` and `o.` hold the
    /// same `Interval<GenomePos>` as `g.` and are equally able to carry an
    /// offset the parser accepts. The verdict is asserted rather than a bare
    /// `is_err()` — the contract is *refused as an invalid position, naming
    /// the offset*, and "failed somehow" would also be satisfied by an
    /// unrelated error on a path where the guard had been removed. The
    /// message's coordinate prefix is asserted too, since the axis label is
    /// passed in per call site and a copy-paste is otherwise invisible.
    #[test]
    fn a_genomic_offset_is_refused_rather_than_dropped() {
        let provider = identity_provider();
        for (descriptor, coord) in [
            ("NC_000001.11:g.10+2delC", "g"),
            ("NC_000001.11:g.10-2delC", "g"),
            ("NC_000001.11:g.10+2_12delCGT", "g"),
            ("NC_000001.11:g.10_12+2delCGT", "g"),
            ("NC_012920.1:m.10+2delC", "m"),
            ("NC_012920.1:m.10_12+2delCGT", "m"),
            ("NC_001416.1:o.10+2delC", "o"),
            ("NC_001416.1:o.10_12+2delCGT", "o"),
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            for (path, converted) in [
                ("without a provider", hgvs_to_spdi_simple(&variant)),
                ("with a provider", hgvs_to_spdi(&variant, &provider)),
            ] {
                let err = match converted {
                    Err(err) => err,
                    Ok(spdi) => {
                        panic!(
                            "`{descriptor}` converted {path} to `{spdi}`; the offset was dropped"
                        )
                    }
                };
                assert!(
                    matches!(err, ConversionError::InvalidPosition { .. }),
                    "`{descriptor}` {path} must be refused as InvalidPosition, got {err:?}"
                );
                let message = err.to_string();
                assert!(
                    message.contains("carries a +/- offset"),
                    "`{descriptor}` {path} was refused for some other reason: {message}"
                );
                assert!(
                    message.contains(&format!("{coord}. position")),
                    "`{descriptor}` {path} named the wrong coordinate axis: {message}"
                );
            }
        }
    }

    /// The transcript axes are where a `+`/`-` offset is legitimate HGVS, and
    /// the genomic guard must not reach them. Each still declines for its own,
    /// unchanged reason — `MissingReferenceData`, "cannot be expressed in SPDI
    /// without genomic projection" — which is a different verdict from the
    /// genomic `InvalidPosition`, and the difference is the point: `c.10+5` is
    /// a well-formed position SPDI cannot carry, `g.10+2` is not a well-formed
    /// position at all.
    ///
    /// Their exonic siblings on the same providers still convert, so this pins
    /// that the axes are working rather than merely erroring.
    #[test]
    fn a_transcript_axis_offset_declines_for_its_own_reason() {
        let provider = make_intronic_provider();
        for (descriptor, exonic_sibling) in [
            ("NM_INTRON.1:c.10+5A>G", "NM_INTRON.1:c.10A>G"),
            ("NM_INTRON.1:n.10+5A>G", "NM_INTRON.1:n.10A>G"),
            ("NM_INTRON.1:r.10+5a>g", "NM_INTRON.1:r.10a>g"),
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            let err = hgvs_to_spdi(&variant, &provider)
                .expect_err("an intronic transcript position has no SPDI representation");
            assert!(
                matches!(err, ConversionError::MissingReferenceData { .. }),
                "`{descriptor}` must still decline as MissingReferenceData, got {err:?}"
            );
            assert!(
                err.to_string().contains("genomic projection"),
                "`{descriptor}` lost its own reason: {err}"
            );
            let exonic = parse_hgvs(exonic_sibling).expect("fixture must parse");
            assert!(
                hgvs_to_spdi(&exonic, &provider).is_ok(),
                "`{exonic_sibling}` must still convert — the axis is not broken, \
                 only its intronic positions are unrepresentable"
            );
        }
    }

    // ------------------------------------------------------------------
    // Genomic special positions (#1643)
    // ------------------------------------------------------------------

    /// `pter`, `qter` and `cen` name no numeric coordinate, and `GenomePos`
    /// stores `base: 0` for all three. `hgvs_to_spdi` read `base` and never
    /// looked at `special`, so every special position resolved to base 0 — the
    /// sibling of the dropped `offset` of #1628, in the same helper and with
    /// the same consequence: descriptions of different variants collapse onto
    /// one triple.
    ///
    /// **Every descriptor here spells its own bases, and that is the whole
    /// reason the class survived #1641.** A description carrying its bases
    /// needs no provider fetch and no position resolution, so none of the
    /// checks that appear to cover this ever runs. The shapes that do not spell
    /// their bases were refused *incidentally* and for the wrong reason —
    /// `g.pter_10del` by the 1-based check ("position 0 is not valid in HGVS"),
    /// `g.10_qterdel` by the fetch (`invalid 1-based interval [10, 0]`), the
    /// `m.`/`o.` forms by the circular-wraparound check reporting a wraparound
    /// that is not there — which is exactly why reading those refusals as
    /// coverage was wrong.
    ///
    /// The verdict is `InvalidPosition` for the reason #1641 established for
    /// the offset case: a special position cannot be *resolved* on a bare
    /// genomic accession — `pter` and `qter` are landmarks of the assembled
    /// chromosome, and `cen` of its centromere annotation, neither of which a
    /// sequence accession carries — so refusing is the only honest answer, and
    /// it differs from the transcript axes' `MissingReferenceData` because
    /// those decline a *well-formed* position SPDI cannot carry.
    ///
    /// Both entry points and every genomic axis, so all six guard call sites
    /// are covered, and the message is asserted rather than a bare `is_err()`:
    /// "failed somehow" is what every one of those incidental refusals already
    /// satisfied.
    #[test]
    fn a_genomic_special_position_is_refused_rather_than_flattened() {
        let provider = identity_provider();
        for (descriptor, coord) in [
            ("NC_000001.11:g.10_qterdelACGTACGTAC", "g"),
            ("NC_000001.11:g.10_cendelACGTACGTAC", "g"),
            ("NC_000001.11:g.10_qterdupACGTACGTAC", "g"),
            ("NC_000001.11:g.10_qterinvACGTACGTAC", "g"),
            ("NC_000001.11:g.pter_10delACGTACGTAC", "g"),
            ("NC_012920.1:m.10_qterdelACGTACGTAC", "m"),
            ("NC_012920.1:m.pter_10delACGTACGTAC", "m"),
            ("NC_001416.1:o.10_qterdelACGTACGTAC", "o"),
            ("NC_001416.1:o.pter_10delACGTACGTAC", "o"),
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            for (path, converted) in [
                ("without a provider", hgvs_to_spdi_simple(&variant)),
                ("with a provider", hgvs_to_spdi(&variant, &provider)),
            ] {
                let err = match converted {
                    Err(err) => err,
                    Ok(spdi) => panic!(
                        "`{descriptor}` converted {path} to `{spdi}`; the special position \
                         was flattened onto base 0"
                    ),
                };
                assert!(
                    matches!(err, ConversionError::InvalidPosition { .. }),
                    "`{descriptor}` {path} must be refused as InvalidPosition, got {err:?}"
                );
                let message = err.to_string();
                assert!(
                    message.contains("names no numeric coordinate"),
                    "`{descriptor}` {path} was refused for some other reason: {message}"
                );
                assert!(
                    message.contains(&format!("{coord}. position")),
                    "`{descriptor}` {path} named the wrong coordinate axis: {message}"
                );
            }
        }
    }

    /// The confluence failure the flattening produced, measured on `main`
    /// before the guard: a deletion running to the q-arm telomere, one running
    /// to the centromere, and a literal ten-base deletion all converted to
    /// `NC_000001.11:9:ACGTACGTAC:`. Three descriptions of three different
    /// things, one triple.
    ///
    /// Stated as the invariant rather than as three pinned errors, so it keeps
    /// meaning something if a future change makes any of these resolvable: what
    /// must never happen is two of them sharing an answer.
    ///
    /// **The denominator is asserted, because with the guard in place the
    /// `panic!` above is unreachable.** Two of the three descriptors are now
    /// refused and skipped by the `if let Ok`, so `seen` holds one triple and
    /// there is no pair left to collide — exactly the `0 of 0` shape this
    /// repository asserts denominators against. The census below is therefore
    /// what carries the meaning today, and it is written to fail the moment a
    /// special position becomes convertible again: that is precisely when the
    /// collision check goes live and has to be re-read rather than trusted.
    #[test]
    fn special_positions_do_not_collapse_onto_one_another() {
        /// Refused by the #1643 guard today, so they never reach `seen`.
        const CONVERTIBLE_TODAY: usize = 1;

        let provider = identity_provider();
        let descriptors = [
            "NC_000001.11:g.10_qterdelACGTACGTAC",
            "NC_000001.11:g.10_cendelACGTACGTAC",
            "NC_000001.11:g.10_19delACGTACGTAC",
        ];
        let mut seen: Vec<(String, String)> = Vec::new();
        for descriptor in descriptors {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            if let Ok(spdi) = hgvs_to_spdi(&variant, &provider) {
                let triple = spdi.to_string();
                if let Some((other, _)) = seen.iter().find(|(_, t)| *t == triple) {
                    panic!(
                        "`{descriptor}` and `{other}` are different variants but share the \
                         triple `{triple}`"
                    );
                }
                seen.push((descriptor.to_string(), triple));
            }
        }

        assert_eq!(
            seen.len(),
            CONVERTIBLE_TODAY,
            "{} of {} descriptors converted, expected {CONVERTIBLE_TODAY}: {seen:?}. \
             If this grew, a special position is convertible again and the collision check \
             above has just become live — read what it now compares before re-pinning this \
             number. If it shrank, the loop is comparing nothing at all",
            seen.len(),
            descriptors.len()
        );
        assert!(
            seen.len() < descriptors.len(),
            "every descriptor converted, so the #1643 guard is refusing nothing"
        );
    }

    /// The invariant the drop breaks, stated directly: two descriptions that
    /// `to_spdi` maps to the same triple must be the same variant. Here the
    /// offset-free sibling is the one legal spelling, so the offset-carrying
    /// ones must not share its triple.
    #[test]
    fn distinct_genomic_descriptions_do_not_share_one_triple() {
        let provider = identity_provider();
        let plain = hgvs_to_spdi(&parse_hgvs("NC_000001.11:g.10delC").unwrap(), &provider)
            .expect("the offset-free spelling is legal and must convert")
            .to_string();
        for descriptor in ["NC_000001.11:g.10+2delC", "NC_000001.11:g.10-2delC"] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            let got = hgvs_to_spdi(&variant, &provider).map(|s| s.to_string());
            assert_ne!(
                got.as_deref().ok(),
                Some(plain.as_str()),
                "`{descriptor}` and `NC_000001.11:g.10delC` are different variants \
                 to `normalize` but one triple to `to_spdi`"
            );
        }
    }

    // ------------------------------------------------------------------
    // Unresolvable transcript-axis end boundaries (#1804)
    // ------------------------------------------------------------------

    /// Assert one refusal carries #1804's verdict, reason and axis label.
    ///
    /// The verdict is checked rather than a bare `is_err()`: on these axes
    /// "failed somehow" is satisfied by the pre-existing `MissingReferenceData`
    /// declines for an intronic or downstream position, so a bare `is_err()`
    /// would pass on a path where this guard had been removed.
    fn assert_end_boundary_refused(
        descriptor: &str,
        coord: &str,
        path: &str,
        converted: Result<SpdiVariant, ConversionError>,
    ) {
        let err = match converted {
            Err(err) => err,
            Ok(spdi) => panic!(
                "`{descriptor}` converted {path} to `{spdi}`; the unresolvable end was \
                 collapsed onto the start"
            ),
        };
        assert!(
            matches!(err, ConversionError::InvalidPosition { .. }),
            "`{descriptor}` {path} must be refused as InvalidPosition — the same verdict this \
             axis already gives an unresolvable START — got {err:?}"
        );
        let message = err.to_string();
        assert!(
            message.contains("names no single coordinate"),
            "`{descriptor}` {path} was refused for some other reason: {message}"
        );
        assert!(
            message.contains(&format!("{coord}. interval end")),
            "`{descriptor}` {path} named the wrong coordinate axis: {message}"
        );
    }

    /// An end boundary that names no single coordinate — a range `(20_30)` or an
    /// unknown `?` — must be refused, not silently replaced by the start.
    ///
    /// Measured on `main` at `439617c2`: `c.10_(20_30)delAAAA` and
    /// `c.10_?delAAAA` both converted to `NM_INTRON.1:19:AAAA:`, and the
    /// short forms were worse — `c.10_(20_30)del` became `NM_INTRON.1:19:A:`,
    /// a one-base deletion standing in for a deletion of unknown extent.
    ///
    /// All three axes on both public entry points where each is reachable, so
    /// every one of the seven substitution sites is covered. `c.` has no
    /// no-provider path (`hgvs_to_spdi_simple` returns `ProviderRequired` for
    /// it), and `r.*5` is likewise provider-only, so those are listed
    /// separately rather than asserted through a path that declines earlier for
    /// an unrelated reason.
    #[test]
    fn an_unresolvable_transcript_end_is_refused_rather_than_collapsed() {
        let provider = make_intronic_provider();

        // Reachable on both entry points.
        for (descriptor, coord) in [
            ("NM_INTRON.1:n.10_(20_30)delAAAA", "n"),
            ("NM_INTRON.1:n.10_?delAAAA", "n"),
            ("NM_INTRON.1:r.10_(20_30)delaaaa", "r"),
            ("NM_INTRON.1:r.10_?delaaaa", "r"),
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            assert_end_boundary_refused(
                descriptor,
                coord,
                "without a provider",
                hgvs_to_spdi_simple(&variant),
            );
            assert_end_boundary_refused(
                descriptor,
                coord,
                "with a provider",
                hgvs_to_spdi(&variant, &provider),
            );
        }

        // Provider-only: `c.` needs the CDS start, and `r.*5` needs the
        // transcript length, so neither reaches the no-provider path.
        for (descriptor, coord) in [
            ("NM_INTRON.1:c.10_(20_30)delAAAA", "c"),
            ("NM_INTRON.1:c.10_?delAAAA", "c"),
            ("NM_INTRON.1:c.10_(20_30)del", "c"),
            ("NM_INTRON.1:c.10_?del", "c"),
            ("NM_INTRON.1:r.*5_?delaaaa", "r"),
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            assert_end_boundary_refused(
                descriptor,
                coord,
                "with a provider",
                hgvs_to_spdi(&variant, &provider),
            );
        }
    }

    /// **The discriminating case, and the reason this is not #1795's guard
    /// again.** `c.(4185+1_4186-1)_(4357+1_4358-1)del` — the exon-deletion shape
    /// [`crate::hgvs::interval::UncertainBoundary`] names as its motivating
    /// example — carries `+`/`-` offsets *inside* the range boundary, and that
    /// is legitimate HGVS on this axis. The genomic guard refuses an offset
    /// outright; here the offset must be beside the point, so the refusal must
    /// cite the absent coordinate and must not mention it.
    ///
    /// A guard written as "walk both endpoints of the range and refuse an
    /// offset" would pass every assertion in the test above and fail this one,
    /// which is what makes it the discriminating case rather than a restatement.
    #[test]
    fn an_offset_inside_the_range_is_not_what_is_refused() {
        let provider = make_intronic_provider();
        for descriptor in [
            "NM_INTRON.1:c.10_(20+1_21-1)del",
            "NM_INTRON.1:c.10_(20+1_21-1)delAAAA",
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            let err = hgvs_to_spdi(&variant, &provider)
                .expect_err("a range end boundary names no single coordinate");
            let message = err.to_string();
            assert!(
                message.contains("names no single coordinate"),
                "`{descriptor}` must be refused for the absent coordinate: {message}"
            );
            assert!(
                !message.contains("offset"),
                "`{descriptor}` was refused for carrying an offset, but an offset is legitimate \
                 on the `c.` axis — the defect is the absent single coordinate: {message}"
            );
        }
    }

    /// The negative control: a boundary that *does* name a coordinate must keep
    /// converting, on every axis and every entry point that accepted it before.
    ///
    /// `(13)` is the case that discriminates a correct guard from an
    /// over-general one. It is `Mu::Uncertain` — a parenthesised but perfectly
    /// numeric position — so `inner()` answers `Some` and it resolves exactly
    /// like the bare `13` beside it. A guard keyed on "the boundary is
    /// parenthesised", or on `Interval::has_complex_boundaries`, would refuse it
    /// and pass every assertion in the tests above.
    ///
    /// **The denominator is asserted** so a run that converted nothing — the
    /// shape a `0 of 0` pass takes — cannot read as a green negative control.
    #[test]
    fn a_transcript_end_that_names_a_coordinate_still_converts() {
        let provider = make_intronic_provider();
        let both_paths = [
            "NM_INTRON.1:n.10_(13)delAAAA",
            "NM_INTRON.1:n.10_13delAAAA",
            "NM_INTRON.1:n.10delA",
            "NM_INTRON.1:n.10A>G",
            "NM_INTRON.1:r.10_(13)delaaaa",
            "NM_INTRON.1:r.10_13delaaaa",
            "NM_INTRON.1:r.10dela",
        ];
        let provider_only = [
            "NM_INTRON.1:c.10_(13)delAAAA",
            "NM_INTRON.1:c.10_13delAAAA",
            "NM_INTRON.1:c.10delA",
            "NM_INTRON.1:c.10del",
        ];

        let mut converted = 0usize;
        for descriptor in both_paths {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            for (path, result) in [
                ("without a provider", hgvs_to_spdi_simple(&variant)),
                ("with a provider", hgvs_to_spdi(&variant, &provider)),
            ] {
                result.unwrap_or_else(|e| {
                    panic!("`{descriptor}` must still convert {path}, got {e}")
                });
                converted += 1;
            }
        }
        for descriptor in provider_only {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            hgvs_to_spdi(&variant, &provider)
                .unwrap_or_else(|e| panic!("`{descriptor}` must still convert, got {e}"));
            converted += 1;
        }

        assert_eq!(
            converted,
            both_paths.len() * 2 + provider_only.len(),
            "the negative control converted {converted} descriptions; if this shrinks the guard \
             has over-generalised, and if the loop stops running it proves nothing"
        );
    }

    /// The confluence failure the collapse produced, stated as the invariant
    /// rather than as pinned triples: descriptions that `parse` and `Display`
    /// distinctly must not share one SPDI triple.
    ///
    /// On `main` all three of these converted to `NM_INTRON.1:19:AAAA:`.
    ///
    /// **The denominator is asserted, because with the guard in place the
    /// collision `panic!` is unreachable** — two of the three are refused and
    /// skipped, so `seen` holds one triple and there is no pair left to collide.
    /// That is the `0 of 0` shape, so the census below is what carries the
    /// meaning today, and it goes red the moment an unresolvable end becomes
    /// convertible again — which is exactly when the collision check goes live
    /// and has to be re-read rather than trusted.
    #[test]
    fn unresolvable_transcript_ends_do_not_collapse_onto_a_resolvable_one() {
        /// Only the fully-resolvable spelling converts today.
        const CONVERTIBLE_TODAY: usize = 1;

        let provider = make_intronic_provider();
        let descriptors = [
            "NM_INTRON.1:c.10_(20_30)delAAAA",
            "NM_INTRON.1:c.10_?delAAAA",
            "NM_INTRON.1:c.10_13delAAAA",
        ];
        let mut seen: Vec<(&str, String)> = Vec::new();
        for descriptor in descriptors {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            if let Ok(spdi) = hgvs_to_spdi(&variant, &provider) {
                let triple = spdi.to_string();
                if let Some((other, _)) = seen.iter().find(|(_, t)| *t == triple) {
                    panic!(
                        "`{descriptor}` and `{other}` are different descriptions but share the \
                         triple `{triple}`"
                    );
                }
                seen.push((descriptor, triple));
            }
        }

        assert_eq!(
            seen.len(),
            CONVERTIBLE_TODAY,
            "{} of {} descriptors converted, expected {CONVERTIBLE_TODAY}: {seen:?}. If this \
             grew, an unresolvable end is convertible again and the collision check above has \
             just become live — read what it now compares before re-pinning this number. If it \
             shrank, the loop is comparing nothing at all",
            seen.len(),
            descriptors.len()
        );
    }

    /// An intronic start beside an unresolvable end is refused for the **end**,
    /// which is a deliberate change of diagnostic and worth pinning as one.
    ///
    /// On `main`, `n.10+5_?delAAAA` reported the intronic start
    /// (`MissingReferenceData`, "cannot be expressed in SPDI without genomic
    /// projection"); it now reports the end. Both are correct refusals, and the
    /// ordering is structural rather than a preference: `resolve_tx_to_provider_tx`
    /// takes a resolved end **by value**, so there is nothing to hand it until
    /// the boundary yields one. Manufacturing a fake end to satisfy that
    /// signature — the start — is the defect being fixed.
    ///
    /// The intronic decline is unchanged where the end *is* resolvable, which is
    /// the second half of this test and the thing a reader will actually want to
    /// know.
    #[test]
    fn an_intronic_start_with_an_unresolvable_end_reports_the_end() {
        let provider = make_intronic_provider();

        let both_unresolvable = parse_hgvs("NM_INTRON.1:n.10+5_?delAAAA").expect("must parse");
        let err = hgvs_to_spdi(&both_unresolvable, &provider)
            .expect_err("neither endpoint has an SPDI representation");
        assert!(
            matches!(err, ConversionError::InvalidPosition { .. }),
            "expected the end-boundary verdict, got {err:?}"
        );
        assert!(
            err.to_string().contains("names no single coordinate"),
            "expected the end-boundary reason, got {err}"
        );

        let intronic_only = parse_hgvs("NM_INTRON.1:n.10+5_20delAAAA").expect("must parse");
        let err = hgvs_to_spdi(&intronic_only, &provider)
            .expect_err("an intronic transcript position has no SPDI representation");
        assert!(
            matches!(err, ConversionError::MissingReferenceData { .. }),
            "a resolvable end must leave the intronic decline untouched, got {err:?}"
        );
        assert!(
            err.to_string().contains("genomic projection"),
            "the intronic decline lost its own reason: {err}"
        );
    }

    // ------------------------------------------------------------------
    // Complex `(a_b)` boundaries (#1795)
    // ------------------------------------------------------------------

    /// A provider carrying all three genomic axes, so the provider-backed arm
    /// of each test below runs against a populated provider rather than an
    /// empty one.
    ///
    /// **No test here currently reaches a fetch, and the negative control is
    /// the proof**: every one of its rows is asserted to convert identically
    /// through [`hgvs_to_spdi_simple`], which holds no provider at all. That is
    /// because each row spells its own bases (`delACGT`) or its own reference
    /// base (`10C>G`), and [`emit_spdi_for_edit`] consults the provider only
    /// for the *unspelled* arms — a bare `del`/`dup`/`inv`. So the provider is
    /// inert for these rows by construction; it is here so that adding a row
    /// without explicit bases does not also require building a fixture, not
    /// because anything below depends on a fetch succeeding.
    fn three_axis_provider() -> MockProvider {
        let mut provider = MockProvider::new();
        for accession in ["NC_000001.11", "NC_012920.1", "NC_001416.1"] {
            provider.add_genomic_sequence(accession, "ACGTACGTACGTACGTACGTACGTACGT".to_string());
        }
        provider
    }

    /// The #1628 guard reads each boundary through `UncertainBoundary::inner`,
    /// which answers `None` for a `Range` — so a prohibited offset written
    /// inside a complex `(a_b)` boundary was invisible to it and converted.
    ///
    /// This is the same prohibition the simple-position control
    /// (`a_genomic_offset_is_refused_rather_than_dropped`) pins, on the same
    /// three axes and the same two entry points, differing only in that the
    /// offset is written inside a range rather than as the whole boundary.
    /// `checklist.md:16` does not care which of the two it is, and neither
    /// should the guard.
    ///
    /// Both endpoints of the range are exercised — `(20+1_30)` and
    /// `(20_30+1)` — because reading only the near one is the obvious
    /// half-fix, and it passes every assertion that names the first.
    /// A range on the *start* is included too: it was already refused, but for
    /// the wrong reason ("unknown start position"), which diagnoses the
    /// uncertainty and never mentions the prohibited offset that is also there.
    #[test]
    fn a_genomic_offset_inside_a_complex_boundary_is_refused_rather_than_dropped() {
        let provider = three_axis_provider();
        for (descriptor, coord) in [
            ("NC_000001.11:g.10_(20+1_30)delACGT", "g"),
            ("NC_000001.11:g.10_(20-1_30)delACGT", "g"),
            ("NC_000001.11:g.10_(20_30+1)delACGT", "g"),
            ("NC_000001.11:g.10_(20+1_30)insTT", "g"),
            ("NC_000001.11:g.(10+1_11)_20delACGT", "g"),
            ("NC_012920.1:m.10_(20+1_30)delACGT", "m"),
            ("NC_001416.1:o.10_(20+1_30)delACGT", "o"),
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            for (path, converted) in [
                ("without a provider", hgvs_to_spdi_simple(&variant)),
                ("with a provider", hgvs_to_spdi(&variant, &provider)),
            ] {
                let err = match converted {
                    Err(err) => err,
                    Ok(spdi) => panic!(
                        "`{descriptor}` converted {path} to `{spdi}`; the offset inside the \
                         complex boundary was dropped"
                    ),
                };
                assert!(
                    matches!(err, ConversionError::InvalidPosition { .. }),
                    "`{descriptor}` {path} must be refused as InvalidPosition, got {err:?}"
                );
                let message = err.to_string();
                assert!(
                    message.contains("carries a +/- offset"),
                    "`{descriptor}` {path} was refused for some other reason: {message}"
                );
                assert!(
                    message.contains(&format!("{coord}. position")),
                    "`{descriptor}` {path} named the wrong coordinate axis: {message}"
                );
            }
        }
    }

    /// The #1643 half of the same guard is blind in exactly the same way: a
    /// `pter`/`qter`/`cen` written inside a complex boundary was never
    /// inspected. Pinned separately from the offset because the two are
    /// different clauses reached through one walk, and a fix that repairs the
    /// walk for one necessarily repairs it for the other — which is only
    /// evidence if both are asserted.
    #[test]
    fn a_genomic_special_position_inside_a_complex_boundary_is_refused() {
        let provider = three_axis_provider();
        for (descriptor, coord) in [
            ("NC_000001.11:g.10_(20_qter)delACGT", "g"),
            ("NC_000001.11:g.10_(cen_30)delACGT", "g"),
            ("NC_012920.1:m.10_(20_qter)delACGT", "m"),
            ("NC_001416.1:o.10_(20_qter)delACGT", "o"),
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            for (path, converted) in [
                ("without a provider", hgvs_to_spdi_simple(&variant)),
                ("with a provider", hgvs_to_spdi(&variant, &provider)),
            ] {
                let err = match converted {
                    Err(err) => err,
                    Ok(spdi) => panic!(
                        "`{descriptor}` converted {path} to `{spdi}`; the special position \
                         inside the complex boundary was flattened"
                    ),
                };
                assert!(
                    matches!(err, ConversionError::InvalidPosition { .. }),
                    "`{descriptor}` {path} must be refused as InvalidPosition, got {err:?}"
                );
                let message = err.to_string();
                assert!(
                    message.contains("names no numeric coordinate"),
                    "`{descriptor}` {path} was refused for some other reason: {message}"
                );
                assert!(
                    message.contains(&format!("{coord}. position")),
                    "`{descriptor}` {path} named the wrong coordinate axis: {message}"
                );
            }
        }
    }

    /// The enabling half, and the one with no offset anywhere in it:
    /// `get_end_pos(..).unwrap_or(start_pos)` silently collapsed an
    /// unresolvable **end** onto the start, so `g.10_(20_30)delACGT` and
    /// `g.10_?delACGT` both answered as though they had been written
    /// `g.10delACGT`.
    ///
    /// The start side has always refused (`ok_or_else`), and so does the VCF
    /// sibling on both endpoints — the fallback was the asymmetry, not the
    /// refusal. SPDI positions are exact; an uncertain or unknown boundary has
    /// no exact coordinate to offer, so declining is the same answer the guard
    /// above gives for the same reason.
    ///
    /// Both unresolvable shapes are covered — a `Range` boundary and a
    /// `Single(Unknown)` `?` — because `inner()` answers `None` for both and a
    /// fix keyed on `is_range()` would leave `?` collapsing. Edits that consume
    /// the span and edits that do not (`ins`, `dup`) are both included, since
    /// they read the interval through different arms of `emit_spdi_for_edit`.
    #[test]
    fn an_unresolvable_genomic_end_boundary_is_refused_rather_than_collapsed() {
        let provider = three_axis_provider();
        for (descriptor, coord) in [
            ("NC_000001.11:g.10_(20_30)delACGT", "g"),
            ("NC_000001.11:g.10_?delACGT", "g"),
            ("NC_000001.11:g.10_(20_30)insTT", "g"),
            ("NC_000001.11:g.10_(20_30)dupACGT", "g"),
            ("NC_012920.1:m.10_(20_30)delACGT", "m"),
            ("NC_012920.1:m.10_?delACGT", "m"),
            ("NC_001416.1:o.10_(20_30)delACGT", "o"),
            ("NC_001416.1:o.10_?delACGT", "o"),
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            for (path, converted) in [
                ("without a provider", hgvs_to_spdi_simple(&variant)),
                ("with a provider", hgvs_to_spdi(&variant, &provider)),
            ] {
                let err = match converted {
                    Err(err) => err,
                    Ok(spdi) => panic!(
                        "`{descriptor}` converted {path} to `{spdi}`; the end boundary was \
                         collapsed onto the start position"
                    ),
                };
                assert!(
                    matches!(err, ConversionError::InvalidPosition { .. }),
                    "`{descriptor}` {path} must be refused as InvalidPosition, got {err:?}"
                );
                let message = err.to_string();
                assert!(
                    message.contains("names no single coordinate"),
                    "`{descriptor}` {path} was refused for some other reason: {message}"
                );
                assert!(
                    message.contains(&format!("{coord}. end boundary")),
                    "`{descriptor}` {path} named the wrong coordinate axis: {message}"
                );
            }
        }
    }

    /// The discriminating half: an end boundary that *does* name a coordinate
    /// must keep converting, and refusing every parenthesised boundary is the
    /// obvious over-generalisation that would pass all three tests above.
    ///
    /// `(13)` is `Mu::Uncertain` — the author is unsure the position is exactly
    /// 13, but 13 is the coordinate they wrote and `inner()` returns it. That
    /// is a different thing from `(12_14)`, which names a *region* and no
    /// single position, and the two must not be conflated because they share a
    /// pair of parentheses.
    #[test]
    fn a_resolvable_genomic_end_boundary_still_converts() {
        let provider = three_axis_provider();
        for (descriptor, expected) in [
            ("NC_000001.11:g.10_13delACGT", "NC_000001.11:9:ACGT:"),
            ("NC_000001.11:g.10_(13)delACGT", "NC_000001.11:9:ACGT:"),
            ("NC_000001.11:g.(10)_(13)delACGT", "NC_000001.11:9:ACGT:"),
            ("NC_000001.11:g.10C>G", "NC_000001.11:9:C:G"),
            ("NC_012920.1:m.10_13delACGT", "NC_012920.1:9:ACGT:"),
            ("NC_001416.1:o.10_13delACGT", "NC_001416.1:9:ACGT:"),
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            for (path, converted) in [
                ("without a provider", hgvs_to_spdi_simple(&variant)),
                ("with a provider", hgvs_to_spdi(&variant, &provider)),
            ] {
                let spdi = converted
                    .unwrap_or_else(|err| panic!("`{descriptor}` {path} must convert: {err}"));
                assert_eq!(
                    spdi.to_string(),
                    expected,
                    "`{descriptor}` {path} moved to a different triple"
                );
            }
        }
    }

    /// The confluence failure #1628 and #1729 were filed to stop, back on this
    /// shape: five descriptions that `parse` and `Display` keep distinct all
    /// converted to `NC_000001.11:9:ACGT:`. Measured on `origin/main` at
    /// `cc8407bc` before the fix.
    ///
    /// Written as the invariant — no two of them may share a triple — rather
    /// than as five pinned errors, so it keeps meaning something if a future
    /// change makes any of these resolvable. As with
    /// `special_positions_do_not_collapse_onto_one_another`, the denominator is
    /// asserted, because with the fix in place four of the five are refused and
    /// the collision check has nothing left to compare: `0 of 0` must not read
    /// as a pass.
    #[test]
    fn complex_boundaries_do_not_collapse_onto_the_offset_free_spelling() {
        /// Only the fully-determined spelling converts once the fix is in.
        const CONVERTIBLE_TODAY: usize = 1;

        let provider = three_axis_provider();
        let descriptors = [
            "NC_000001.11:g.10_(20+1_30)delACGT",
            "NC_000001.11:g.10_(20+5_30)delACGT",
            "NC_000001.11:g.10_(20-1_30)delACGT",
            "NC_000001.11:g.10_(20_30)delACGT",
            "NC_000001.11:g.10delACGT",
        ];

        // The premise: these are five distinct descriptions, not five
        // spellings the parser has already unified. If this ever fails, the
        // collapse is happening upstream of the conversion and the rest of
        // this test is measuring the wrong component.
        let mut displayed: Vec<String> = descriptors
            .iter()
            .map(|d| parse_hgvs(d).expect("fixture must parse").to_string())
            .collect();
        displayed.sort();
        displayed.dedup();
        assert_eq!(
            displayed.len(),
            descriptors.len(),
            "the descriptions do not survive parse/Display distinctly: {displayed:?}"
        );

        let mut seen: Vec<(&str, String)> = Vec::new();
        for descriptor in descriptors {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            if let Ok(spdi) = hgvs_to_spdi(&variant, &provider) {
                let triple = spdi.to_string();
                if let Some((other, _)) = seen.iter().find(|(_, t)| *t == triple) {
                    panic!(
                        "`{descriptor}` and `{other}` are different descriptions but share \
                         the triple `{triple}`"
                    );
                }
                seen.push((descriptor, triple));
            }
        }

        assert_eq!(
            seen.len(),
            CONVERTIBLE_TODAY,
            "{} of {} descriptors converted, expected {CONVERTIBLE_TODAY}: {seen:?}. \
             If this grew, a complex boundary is convertible again and the collision check \
             above has just become live — read what it now compares before re-pinning this \
             number. If it shrank, `g.10delACGT` stopped converting and the loop is \
             comparing nothing at all",
            seen.len(),
            descriptors.len()
        );
    }
}
