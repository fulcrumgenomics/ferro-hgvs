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
//! | HGVS | SPDI |
//! |------|------|
//! | Substitution `g.12345A>G` | `seq:12344:A:G` |
//! | Deletion `g.100_102del` | `seq:99:NNN:` (requires reference) |
//! | Insertion `g.100_101insATG` | `seq:100::ATG` |
//! | Delins `g.100_102delinsATG` | `seq:99:NNN:ATG` (requires reference) |
//! | Duplication `g.100_102dup` | `seq:102::NNN` (requires reference) |
//!
//! Note: Deletions, delins, and duplications require reference sequence data
//! to determine the deleted sequence.
//!
//! # Coordinate-system support
//!
//! [`hgvs_to_spdi_simple`] accepts coordinate systems whose positions are
//! resolvable without provider data:
//!
//! - `g.` (genomic) — direct.
//! - `m.` (mito) — the mito accession is genomic; same path as `g.`.
//! - `n.` (non-coding tx) — exonic, positive base; SPDI on the transcript
//!   accession.
//! - `r.` (RNA) — exonic, positive base; `u`/`U` rewritten to `T` for
//!   SPDI's DNA alphabet convention.
//!
//! [`hgvs_to_spdi`] additionally handles `c.` (CDS) and intronic / UTR
//! `n.`/`r.` positions by consulting a [`ReferenceProvider`] for transcript
//! metadata. Per the [SPDI spec], the SPDI accession matches the HGVS
//! accession (NCBI Variation Services emits SPDI on transcript accessions
//! the same way).
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
use crate::hgvs::edit::{InsertedSequence, NaEdit, Sequence};
use crate::hgvs::interval::Interval;
use crate::hgvs::location::{CdsPos, GenomePos, RnaPos, TxPos};
use crate::hgvs::parser::accession::parse_accession;
use crate::hgvs::variant::{
    Accession, CdsVariant, GenomeVariant, HgvsVariant, LocEdit, MtVariant, RnaVariant, TxVariant,
};
use crate::reference::provider::ReferenceProvider;
use crate::reference::transcript::Transcript;

/// Error type for conversion failures.
#[derive(Debug, Clone, PartialEq, Eq)]
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

/// Helper to get start position from an interval.
fn get_start_pos(interval: &Interval<GenomePos>) -> Option<u64> {
    interval.start.inner().map(|p| p.base)
}

/// Helper to get end position from an interval.
fn get_end_pos(interval: &Interval<GenomePos>) -> Option<u64> {
    interval.end.inner().map(|p| p.base)
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
/// | `g.` (genomic) | yes | direct 1→0-based conversion |
/// | `m.` (mito) | yes | mito accession is genomic; same as `g.` |
/// | `n.` (non-coding tx) | exonic, positive base | SPDI sits on the transcript accession |
/// | `r.` (RNA) | exonic, positive base | `u`/`U` rewritten to `T`; SPDI uses DNA alphabet |
/// | `c.` (CDS) | NO — needs CDS start | use [`hgvs_to_spdi`] |
/// | `p.` (protein) | NO | not representable in SPDI |
///
/// For `c.`, intronic `n.`/`r.`, UTR `n.`/`r.`, or deletion/dup variants
/// without an explicit deleted sequence, use [`hgvs_to_spdi`] which consults
/// a [`ReferenceProvider`].
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

/// Convert an HGVS variant to SPDI, consulting a reference provider when
/// transcript metadata (CDS start, exon coordinates) is required.
///
/// This is the provider-aware companion to [`hgvs_to_spdi_simple`]. It
/// handles `c.` (CDS) variants, `n.`/`r.` variants with intronic offsets or
/// UTR-style positions, and falls back to the simple path for cases that
/// don't need provider data.
///
/// The resulting SPDI uses the **same accession** as the HGVS variant — for
/// `NM_000088.3:c.1A>G` the SPDI sequence is `NM_000088.3`, not the
/// underlying genomic accession. This matches NCBI Variation Services'
/// behavior: SPDI is positional on whichever accession is provided.
///
/// # Arguments
///
/// * `variant` - The HGVS variant to convert.
/// * `provider` - A reference provider that can return the relevant
///   transcript via `get_transcript`.
///
/// # Errors
///
/// * [`ConversionError::MissingReferenceData`] if the provider does not have
///   the transcript or the position cannot be resolved (e.g. intronic
///   without exon data).
/// * Other `ConversionError` variants for the same reasons as
///   [`hgvs_to_spdi_simple`].
pub fn hgvs_to_spdi<P: ReferenceProvider>(
    variant: &HgvsVariant,
    provider: &P,
) -> Result<SpdiVariant, ConversionError> {
    match variant {
        HgvsVariant::Genome(g) => genome_to_spdi_simple(g),
        HgvsVariant::Mt(m) => mt_to_spdi_simple(m),
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
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_pos = get_start_pos(&variant.loc_edit.location).ok_or_else(|| {
        ConversionError::InvalidPosition {
            description: "cannot convert variant with unknown start position".to_string(),
        }
    })?;
    let end_pos = get_end_pos(&variant.loc_edit.location).unwrap_or(start_pos);
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_pos,
        end_pos,
        edit,
        AlphabetMode::Dna,
    )
}

/// Convert a mitochondrial variant to SPDI (simple conversion).
///
/// Mitochondrial accessions (e.g. `NC_012920.1`) are themselves genomic
/// accessions, so the conversion is identical to the `g.` path with a
/// different coordinate prefix on the HGVS side.
fn mt_to_spdi_simple(variant: &MtVariant) -> Result<SpdiVariant, ConversionError> {
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_pos = get_start_pos(&variant.loc_edit.location).ok_or_else(|| {
        ConversionError::InvalidPosition {
            description: "cannot convert variant with unknown start position".to_string(),
        }
    })?;
    let end_pos = get_end_pos(&variant.loc_edit.location).unwrap_or(start_pos);
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_pos,
        end_pos,
        edit,
        AlphabetMode::Dna,
    )
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
    let end_tx = tx_end_for_simple_path(&variant.loc_edit.location, start_tx, "n")?;
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_tx,
        end_tx,
        edit,
        AlphabetMode::Dna,
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
    let end_pos = rna_end_for_simple_path(&variant.loc_edit.location, start_pos, "r")?;
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_pos,
        end_pos,
        edit,
        AlphabetMode::Rna,
    )
}

/// Convert a CDS (`c.`) variant to SPDI by resolving CDS coordinates to
/// transcript positions through the supplied provider.
///
/// The resulting SPDI uses the variant's transcript accession (e.g.
/// `NM_000088.3`), matching NCBI Variation Services' convention.
fn cds_to_spdi_with_provider<P: ReferenceProvider>(
    variant: &CdsVariant,
    provider: &P,
) -> Result<SpdiVariant, ConversionError> {
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_cds = variant.loc_edit.location.start.inner().ok_or_else(|| {
        ConversionError::InvalidPosition {
            description: "cannot convert c. variant with unknown start position".to_string(),
        }
    })?;
    let end_cds = variant
        .loc_edit
        .location
        .end
        .inner()
        .copied()
        .unwrap_or(*start_cds);
    let (start_tx, end_tx) = resolve_cds_to_tx(&variant.accession, start_cds, &end_cds, provider)?;
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_tx,
        end_tx,
        edit,
        AlphabetMode::Dna,
    )
}

/// Convert an `n.` variant using a provider for cases that the simple path
/// can't handle (intronic offsets, downstream positions, non-positive base).
fn tx_to_spdi_with_provider<P: ReferenceProvider>(
    variant: &TxVariant,
    provider: &P,
) -> Result<SpdiVariant, ConversionError> {
    // Fast path: positions resolvable without provider data
    if !tx_needs_provider(&variant.loc_edit.location) {
        return tx_to_spdi_simple(variant);
    }
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_tx_pos = variant.loc_edit.location.start.inner().ok_or_else(|| {
        ConversionError::InvalidPosition {
            description: "cannot convert n. variant with unknown start position".to_string(),
        }
    })?;
    let end_tx_pos = variant
        .loc_edit
        .location
        .end
        .inner()
        .copied()
        .unwrap_or(*start_tx_pos);
    let (start_tx, end_tx) =
        resolve_tx_to_provider_tx(&variant.accession, start_tx_pos, &end_tx_pos, provider)?;
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_tx,
        end_tx,
        edit,
        AlphabetMode::Dna,
    )
}

/// Convert an `r.` variant using a provider for cases the simple path can't
/// handle. Same coordinate resolution as `n.`; alphabet conversion `u → T`
/// is applied via [`AlphabetMode::Rna`].
fn rna_to_spdi_with_provider<P: ReferenceProvider>(
    variant: &RnaVariant,
    provider: &P,
) -> Result<SpdiVariant, ConversionError> {
    if !rna_needs_provider(&variant.loc_edit.location) {
        return rna_to_spdi_simple(variant);
    }
    let edit = unwrap_edit(&variant.loc_edit.edit)?;
    let start_rna = variant.loc_edit.location.start.inner().ok_or_else(|| {
        ConversionError::InvalidPosition {
            description: "cannot convert r. variant with unknown start position".to_string(),
        }
    })?;
    let end_rna = variant
        .loc_edit
        .location
        .end
        .inner()
        .copied()
        .unwrap_or(*start_rna);
    let (start_tx, end_tx) =
        resolve_rna_to_provider_tx(&variant.accession, start_rna, &end_rna, provider)?;
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_tx,
        end_tx,
        edit,
        AlphabetMode::Rna,
    )
}

// ===========================================================================
// Shared helpers
// ===========================================================================

/// Whether to rewrite RNA bases (`u/U`) to DNA (`T`) in the emitted SPDI
/// deletion/insertion strings.
#[derive(Debug, Clone, Copy)]
enum AlphabetMode {
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

fn tx_end_for_simple_path(
    interval: &Interval<TxPos>,
    fallback: u64,
    coord: &str,
) -> Result<u64, ConversionError> {
    match interval.end.inner() {
        Some(end) => require_simple_tx_pos(end, coord),
        None => Ok(fallback),
    }
}

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
    fallback: u64,
    coord: &str,
) -> Result<u64, ConversionError> {
    match interval.end.inner() {
        Some(end) => require_simple_rna_pos(end, coord),
        None => Ok(fallback),
    }
}

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
fn resolve_cds_to_tx<P: ReferenceProvider>(
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
    let s_u = ensure_positive_tx(s.base, "c", start)?;
    let e_u = ensure_positive_tx(e.base, "c", end)?;
    Ok((s_u, e_u))
}

/// Resolve an `n.` (TxPos) pair to 1-based transcript positions, including
/// intronic and downstream forms, using the provider.
fn resolve_tx_to_provider_tx<P: ReferenceProvider>(
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

fn resolve_rna_to_provider_tx<P: ReferenceProvider>(
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

/// Resolve a single `TxPos` to a 1-based transcript position. For intronic
/// positions, the returned position is the nearest exonic transcript base,
/// because SPDI cannot directly represent intronic offsets. For downstream
/// (*N) positions, the position is offset past the transcript end.
fn resolve_tx_pos(pos: &TxPos, transcript: &Transcript) -> Result<u64, ConversionError> {
    // Reject intronic + downstream-with-offset until #117/#118 wire ref-aware
    // edit material in. SPDI cannot directly carry an intronic offset, so we
    // surface a clear error rather than silently dropping it.
    if pos.is_intronic() {
        return Err(ConversionError::MissingReferenceData {
            description: format!(
                "intronic n.{} cannot be expressed in SPDI without genomic projection; \
                 SPDI is positional and has no offset notation",
                pos
            ),
        });
    }
    let tx_len = transcript.sequence_length();
    if pos.is_downstream() {
        if pos.base < 1 {
            return Err(ConversionError::InvalidPosition {
                description: format!("downstream position *{} must be >= 1", pos.base),
            });
        }
        let value = tx_len.saturating_add(pos.base as u64);
        return Ok(value);
    }
    ensure_positive_tx(pos.base, "n", pos)
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
    let tx_len = transcript.sequence_length();
    if pos.utr3 {
        if pos.base < 1 {
            return Err(ConversionError::InvalidPosition {
                description: format!("3' UTR position *{} must be >= 1", pos.base),
            });
        }
        return Ok(tx_len.saturating_add(pos.base as u64));
    }
    ensure_positive_tx(pos.base, "r", pos)
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

/// Apply edit-specific position arithmetic and emit the SPDI variant.
///
/// `start_one_based` and `end_one_based` are 1-based positions on the SPDI
/// accession (genomic for `g.`/`m.`, transcript for `c.`/`n.`/`r.`).
fn emit_spdi_for_edit(
    sequence: String,
    start_one_based: u64,
    end_one_based: u64,
    edit: &NaEdit,
    alphabet: AlphabetMode,
) -> Result<SpdiVariant, ConversionError> {
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
            let ins_str = inserted_sequence_to_string(inserted).ok_or_else(|| {
                ConversionError::MissingReferenceData {
                    description: "insertion sequence is not a literal sequence".to_string(),
                }
            })?;
            Ok(SpdiVariant::new(
                sequence,
                spdi_pos,
                "",
                apply_alphabet(&ins_str, alphabet),
            ))
        }
        NaEdit::Duplication {
            sequence: dup_seq, ..
        } => {
            if let Some(seq) = dup_seq {
                let end_pos_ob = OneBasedPos::new(end_one_based);
                let spdi_end_zb = end_pos_ob.to_zero_based();
                Ok(SpdiVariant::new(
                    sequence,
                    spdi_end_zb.value(),
                    "",
                    apply_alphabet(&sequence_to_string(seq), alphabet),
                ))
            } else {
                Err(ConversionError::MissingReferenceData {
                    description:
                        "duplication sequence not provided; reference data needed to determine duplicated bases"
                            .to_string(),
                })
            }
        }
        NaEdit::Deletion {
            sequence: del_seq, ..
        } => {
            if let Some(seq) = del_seq {
                Ok(SpdiVariant::new(
                    sequence,
                    spdi_pos,
                    apply_alphabet(&sequence_to_string(seq), alphabet),
                    "",
                ))
            } else {
                Err(ConversionError::MissingReferenceData {
                    description:
                        "deleted sequence not provided; reference data needed to determine deleted bases"
                            .to_string(),
                })
            }
        }
        NaEdit::Delins { sequence: _ins_seq } => {
            let del_len = end_one_based
                .saturating_sub(start_one_based)
                .saturating_add(1) as usize;
            Err(ConversionError::MissingReferenceData {
                description: format!(
                    "Cannot convert delins to SPDI: deleted sequence of length {} is unknown (no reference data)",
                    del_len
                ),
            })
        }
        NaEdit::Identity {
            sequence: id_seq, ..
        } => {
            let ref_base = id_seq
                .as_ref()
                .map(|s| apply_alphabet(&sequence_to_string(s), alphabet))
                .unwrap_or_default();
            Ok(SpdiVariant::new(
                sequence,
                spdi_pos,
                ref_base.clone(),
                ref_base,
            ))
        }
        NaEdit::Inversion { .. } => Err(ConversionError::UnsupportedEditType {
            description: "inversions cannot be represented in SPDI format".to_string(),
        }),
        NaEdit::Repeat { .. } => Err(ConversionError::UnsupportedEditType {
            description: "repeat variants cannot be directly converted to SPDI".to_string(),
        }),
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

/// Rewrite RNA-alphabet characters to the DNA alphabet for SPDI output.
///
/// SPDI uses the DNA alphabet by convention (RefSeq stores transcript
/// sequences as DNA even on `NR_*` and `NM_*` accessions), so RNA `u`/`U`
/// must become `T`. Other characters are returned unchanged. The output is
/// always uppercase to match SPDI's standard form.
fn apply_alphabet(s: &str, alphabet: AlphabetMode) -> String {
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
        // Identity
        let seq = if spdi.deletion.is_empty() {
            None
        } else {
            Some(string_to_sequence(&spdi.deletion)?)
        };
        (
            Interval::point(GenomePos::new(hgvs_pos)),
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
        // Pure insertion
        // SPDI position is where insertion happens (interbase)
        // HGVS: g.pos_pos+1insX means insert between pos and pos+1
        let ins_seq = string_to_sequence(&spdi.insertion)?;
        (
            Interval::new(GenomePos::new(hgvs_pos), GenomePos::new(hgvs_pos + 1)),
            NaEdit::Insertion {
                sequence: InsertedSequence::Literal(ins_seq),
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
            },
        )
    };

    Ok(HgvsVariant::Genome(GenomeVariant {
        accession,
        gene_symbol: None,
        loc_edit: LocEdit::new(interval, edit),
    }))
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

    // Need `ins_len` bases of preceding reference. Under ferro's SPDI
    // insertion convention, `spdi.position` is the 0-based offset of the
    // 1-based base immediately 5' of the insertion point, so the
    // 5'-flanking window is the 0-based half-open interval
    // `[spdi.position + 1 - ins_len, spdi.position + 1)`.
    let flank_end =
        spdi.position
            .checked_add(1)
            .ok_or_else(|| ConversionError::InvalidPosition {
                description: format!("SPDI position {} overflows on +1", spdi.position),
            })?;
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

    // Build the dup edit. 1-based interval: end = spdi.position + 1,
    // start = end + 1 - ins_len.
    let end_one_based = flank_end; // == spdi.position + 1
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

    // Reuse the accession (and gene_symbol, which is None for SPDI inputs)
    // from the base ins-form variant. SPDI→HGVS only produces Genome
    // variants, so we expect HgvsVariant::Genome here.
    let HgvsVariant::Genome(g) = base else {
        // Defensive: SPDI→HGVS should always produce Genome. If it doesn't,
        // do not attempt dup recovery.
        return Ok(None);
    };

    Ok(Some(HgvsVariant::Genome(GenomeVariant {
        accession: g.accession.clone(),
        gene_symbol: g.gene_symbol.clone(),
        loc_edit: LocEdit::new(interval, edit),
    })))
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
        assert_eq!(spdi.position, 99); // 0-based, position before insertion
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
    fn test_hgvs_to_spdi_duplication_with_seq() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_102dupATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        // Dup becomes insertion after the duplicated region
        assert_eq!(spdi.position, 101); // end position, 0-based
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
    fn test_hgvs_to_spdi_unsupported_inversion() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_200inv").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::UnsupportedEditType { .. })
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
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.101_102insATG");
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
    // P1: Reference-dependent conversion tests
    // =========================================================================

    /// Mock genomic reference for testing reference-dependent conversions.
    struct MockGenomicRef {
        /// Map of (accession, start, end) -> sequence
        sequences: std::collections::HashMap<(String, u64, u64), String>,
    }

    impl MockGenomicRef {
        fn new() -> Self {
            let mut sequences = std::collections::HashMap::new();

            // Add test sequences for NC_000001.11
            // Positions are 0-based for internal storage
            // Sequence at g.100_102 (1-based) = indices 99-101 (0-based) = "ATG"
            sequences.insert(("NC_000001.11".to_string(), 99, 102), "ATG".to_string());
            // Sequence at g.200_205 (1-based) = indices 199-204 (0-based) = "GATTACA"
            sequences.insert(
                ("NC_000001.11".to_string(), 199, 206),
                "GATTACA".to_string(),
            );
            // Sequence at g.1000_1009 (1-based) = "AAACCCGGGT"
            sequences.insert(
                ("NC_000001.11".to_string(), 999, 1009),
                "AAACCCGGGT".to_string(),
            );

            Self { sequences }
        }

        fn get_sequence(&self, accession: &str, start: u64, end: u64) -> Option<String> {
            self.sequences
                .get(&(accession.to_string(), start, end))
                .cloned()
        }
    }

    /// Convert deletion HGVS to SPDI using reference data.
    fn hgvs_to_spdi_with_ref(
        variant: &HgvsVariant,
        reference: &MockGenomicRef,
    ) -> Result<SpdiVariant, ConversionError> {
        match variant {
            HgvsVariant::Genome(g) => {
                let sequence = g.accession.to_string();
                let interval = &g.loc_edit.location;

                let edit =
                    g.loc_edit
                        .edit
                        .inner()
                        .ok_or_else(|| ConversionError::InvalidPosition {
                            description: "cannot convert variant with unknown edit".to_string(),
                        })?;

                let start_pos =
                    get_start_pos(interval).ok_or_else(|| ConversionError::InvalidPosition {
                        description: "cannot convert variant with unknown start position"
                            .to_string(),
                    })?;
                let end_pos = get_end_pos(interval).unwrap_or(start_pos);

                // Convert 1-based HGVS position to 0-based SPDI position
                let start_pos_ob = OneBasedPos::new(start_pos);
                let spdi_pos = start_pos_ob.to_zero_based().value();

                match edit {
                    NaEdit::Deletion { sequence: None, .. } => {
                        // Fetch sequence from reference
                        let del_seq = reference
                            .get_sequence(&sequence, spdi_pos, end_pos)
                            .ok_or_else(|| ConversionError::MissingReferenceData {
                                description: format!(
                                    "could not fetch sequence for {}:{}-{}",
                                    sequence, start_pos, end_pos
                                ),
                            })?;
                        Ok(SpdiVariant::new(sequence, spdi_pos, del_seq, ""))
                    }
                    NaEdit::Delins { sequence: ins_seq } => {
                        // Fetch deleted sequence from reference
                        let del_seq = reference
                            .get_sequence(&sequence, spdi_pos, end_pos)
                            .ok_or_else(|| ConversionError::MissingReferenceData {
                                description: format!(
                                    "could not fetch sequence for {}:{}-{}",
                                    sequence, start_pos, end_pos
                                ),
                            })?;
                        let ins_str = inserted_sequence_to_string(ins_seq).ok_or_else(|| {
                            ConversionError::MissingReferenceData {
                                description: "delins inserted sequence is not a literal"
                                    .to_string(),
                            }
                        })?;
                        Ok(SpdiVariant::new(sequence, spdi_pos, del_seq, ins_str))
                    }
                    NaEdit::Duplication { sequence: None, .. } => {
                        // Fetch duplicated sequence from reference
                        let dup_seq = reference
                            .get_sequence(&sequence, spdi_pos, end_pos)
                            .ok_or_else(|| ConversionError::MissingReferenceData {
                                description: format!(
                                    "could not fetch sequence for {}:{}-{}",
                                    sequence, start_pos, end_pos
                                ),
                            })?;
                        // Duplication becomes insertion after the region
                        // Convert 1-based HGVS end position to 0-based SPDI position
                        let end_pos_ob = OneBasedPos::new(end_pos);
                        let spdi_end = end_pos_ob.to_zero_based().value();
                        Ok(SpdiVariant::new(sequence, spdi_end, "", dup_seq))
                    }
                    // Fall back to simple conversion for other types
                    _ => hgvs_to_spdi_simple(variant),
                }
            }
            _ => hgvs_to_spdi_simple(variant),
        }
    }

    #[test]
    fn test_deletion_without_seq_with_reference() {
        let reference = MockGenomicRef::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102del").unwrap();

        // This would fail with simple conversion
        assert!(hgvs_to_spdi_simple(&hgvs).is_err());

        // But works with reference data
        let spdi = hgvs_to_spdi_with_ref(&hgvs, &reference).unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn test_delins_with_reference_provides_actual_deletion() {
        let reference = MockGenomicRef::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102delinsTTCC").unwrap();

        // Simple conversion without reference data returns an error
        let simple = hgvs_to_spdi_simple(&hgvs);
        assert!(simple.is_err());
        assert!(simple.unwrap_err().to_string().contains("unknown"));

        // Reference-based conversion gets actual sequence
        let spdi = hgvs_to_spdi_with_ref(&hgvs, &reference).unwrap();
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "TTCC");
    }

    #[test]
    fn test_delins_roundtrip_with_reference() {
        let reference = MockGenomicRef::new();
        let original_hgvs = parse_hgvs("NC_000001.11:g.100_102delinsTTCC").unwrap();

        // Convert to SPDI with reference
        let spdi = hgvs_to_spdi_with_ref(&original_hgvs, &reference).unwrap();
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "TTCC");

        // Convert back to HGVS
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(back.to_string(), "NC_000001.11:g.100_102delinsTTCC");
    }

    #[test]
    fn test_deletion_roundtrip_with_reference() {
        let reference = MockGenomicRef::new();
        let original_hgvs = parse_hgvs("NC_000001.11:g.100_102del").unwrap();

        // Convert to SPDI with reference
        let spdi = hgvs_to_spdi_with_ref(&original_hgvs, &reference).unwrap();
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "");

        // Convert back to HGVS
        let back = spdi_to_hgvs(&spdi).unwrap();
        // Note: HGVS back includes the sequence
        assert_eq!(back.to_string(), "NC_000001.11:g.100_102delATG");
    }

    #[test]
    fn test_duplication_without_seq_with_reference() {
        let reference = MockGenomicRef::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102dup").unwrap();

        // This would fail with simple conversion
        assert!(hgvs_to_spdi_simple(&hgvs).is_err());

        // But works with reference data
        let spdi = hgvs_to_spdi_with_ref(&hgvs, &reference).unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        // Dup is insertion after the region
        assert_eq!(spdi.position, 101); // end position 0-based
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn test_long_deletion_with_reference() {
        let reference = MockGenomicRef::new();
        let hgvs = parse_hgvs("NC_000001.11:g.1000_1009del").unwrap();

        let spdi = hgvs_to_spdi_with_ref(&hgvs, &reference).unwrap();
        assert_eq!(spdi.deletion, "AAACCCGGGT");
        assert_eq!(spdi.deletion.len(), 10);
    }

    #[test]
    fn test_reference_missing_region() {
        let reference = MockGenomicRef::new();
        // Position not in mock reference
        let hgvs = parse_hgvs("NC_000001.11:g.50000_50002del").unwrap();

        let result = hgvs_to_spdi_with_ref(&hgvs, &reference);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    #[test]
    fn test_substitution_falls_through_to_simple() {
        let reference = MockGenomicRef::new();
        let hgvs = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();

        // Substitution doesn't need reference, falls through
        let spdi = hgvs_to_spdi_with_ref(&hgvs, &reference).unwrap();
        assert_eq!(spdi.to_string(), "NC_000001.11:12344:A:G");
    }

    #[test]
    fn test_mnv_delins_with_reference() {
        // Multi-nucleotide variant (MNV) as delins
        let reference = MockGenomicRef::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102delinsGGG").unwrap();

        let spdi = hgvs_to_spdi_with_ref(&hgvs, &reference).unwrap();
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "GGG");

        // MNV same length should roundtrip
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(back.to_string(), "NC_000001.11:g.100_102delinsGGG");
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

        // SPDI 101::ATG (the canonical SPDI form of g.100_102dupATG per ferro
        // convention; verified by test_hgvs_to_spdi_duplication_with_seq).
        let spdi = SpdiVariant::insertion("NC_000001.11", 101, "ATG");

        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100_102dupATG");
    }

    #[test]
    fn spdi_to_hgvs_with_ref_recovers_single_base_dup() {
        // 1-based base 100 = 'A' → SPDI insertion at 99::A is a single-base dup.
        let mut contig = "N".repeat(99);
        contig.push('A'); // 0-based 99 = 1-based 100
        contig.push_str(&"N".repeat(20));
        let provider = provider_with_genomic(&contig);

        // g.100dupA → SPDI 99::A (5' base is 1-based 100 = 0-based 99).
        let spdi = SpdiVariant::insertion("NC_000001.11", 99, "A");

        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100dupA");
    }

    #[test]
    fn spdi_to_hgvs_with_ref_keeps_ins_when_no_match() {
        // 5' flank of length 3 ending at SPDI position 101 = 1-based bases
        // 100..102 = "CCC" — does NOT equal the inserted "ATG".
        let mut contig = "N".repeat(99);
        contig.push_str("CCC"); // 1-based 100..102 = "CCC"
        contig.push_str(&"N".repeat(20));
        let provider = provider_with_genomic(&contig);

        let spdi = SpdiVariant::insertion("NC_000001.11", 101, "ATG");

        // Should fall through to ins-form. Existing spdi_to_hgvs renders an
        // insertion as `g.{pos}_{pos+1}insATG` where pos = HGVS 1-based of
        // SPDI position (102 here, since SPDI 101 → 1-based 102). The
        // ins-form HGVS is therefore g.102_103insATG.
        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.102_103insATG");
    }

    #[test]
    fn spdi_to_hgvs_with_ref_keeps_ins_at_contig_start() {
        // SPDI position 0 with a 3-base insertion: there are no preceding
        // bases at all (flank_end=1, ins_len=3, flank_end < ins_len).
        let provider = provider_with_genomic("ATGCATGCATGC");

        let spdi = SpdiVariant::insertion("NC_000001.11", 0, "ATG");
        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        // ins-form: SPDI 0 → 1-based 1 → g.1_2insATG
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.1_2insATG");
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
        assert_eq!(spdi.position, 99);
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

        // Forward: HGVS dup → SPDI ins
        let original = "NC_000001.11:g.100_102dupATG";
        let hgvs = parse_hgvs(original).unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.to_string(), "NC_000001.11:101::ATG");

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
        assert_eq!(spdi.to_string(), "NC_000001.11:99::A");

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
        // 101::ATG (which round-tripped from g.100_102dupATG) is rendered
        // as g.102_103insATG. If a future change attempts to "fix" the
        // round-trip without a reference, this audit pin will fail and
        // demand explicit reconsideration.
        let spdi = SpdiVariant::insertion("NC_000001.11", 101, "ATG");
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

        let spdi = SpdiVariant::insertion("NC_012920.1", 101, "ATG");
        let hgvs = spdi_to_hgvs_with_ref(&spdi, &provider).unwrap();
        // Even for mito, the SPDI→HGVS path produces a g. variant; this
        // matches the existing convention in this module.
        assert_eq!(hgvs.to_string(), "NC_012920.1:g.100_102dupATG");
    }

    #[test]
    fn test_hgvs_to_spdi_simple_mt_dup_with_seq() {
        let hgvs = parse_hgvs("NC_012920.1:m.100_102dupATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        // Dup converts to insertion at the end of the duplicated region (0-based).
        assert_eq!(spdi.position, 101);
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
        assert_eq!(spdi.position, 9);
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
        let hgvs = parse_hgvs("NR_046018.2:n.*5A>G").unwrap();
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
        assert_eq!(spdi.position, 9);
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
        // c.1_2insATG: cds_start=6 → tx 6_7 → SPDI position 5 (interbase)
        let hgvs = parse_hgvs("NM_TEST.1:c.1_2insATG").unwrap();
        let spdi = hgvs_to_spdi(&hgvs, &provider).unwrap();
        assert_eq!(spdi.position, 5);
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

    // ----- Edit-type passthrough on new coord systems -----------------------

    #[test]
    fn test_hgvs_to_spdi_simple_n_inversion_unsupported() {
        let hgvs = parse_hgvs("NR_046018.2:n.10_20inv").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::UnsupportedEditType { .. })
        ));
    }

    #[test]
    fn test_hgvs_to_spdi_simple_m_inversion_unsupported() {
        let hgvs = parse_hgvs("NC_012920.1:m.100_200inv").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::UnsupportedEditType { .. })
        ));
    }
}
