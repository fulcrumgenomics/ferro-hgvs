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
//! | Duplication `g.100_102dup` | `seq:101::ATG` | required (short form) |
//! | Duplication `g.100_102dupATG` | `seq:101::ATG` | not required (explicit form) |
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
//! - `g.` (genomic) — direct.
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
use crate::hgvs::interval::Interval;
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
/// [`ConversionError::MissingReferenceData`].
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
/// short-form deletions / delins (and the symmetric `ins` field for
/// short-form duplications) by fetching the reference bases for the
/// variant's interval. Explicit-form input (`g.100_102delATG`,
/// `g.100_102dupATG`, etc.) emits the user-supplied bases as-is and does
/// not consult the provider. Intronic `n.`/`r.` positions remain
/// unsupported (SPDI has no offset notation) and return
/// [`ConversionError::MissingReferenceData`].
///
/// The resulting SPDI uses the **same accession** as the HGVS variant — for
/// `NM_000088.3:c.1A>G` the SPDI sequence is `NM_000088.3`, not the
/// underlying genomic accession. This matches NCBI Variation Services'
/// behavior: SPDI is positional on whichever accession is provided.
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
    let end_pos = get_end_pos(&variant.loc_edit.location).unwrap_or(start_pos);
    emit_spdi_for_edit(
        variant.accession.to_string(),
        start_pos,
        end_pos,
        edit,
        AlphabetMode::Dna,
        None::<&dyn ReferenceProvider>,
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
    let end_pos = rna_end_for_simple_path(&variant.loc_edit.location, start_pos, "r")?;
    emit_spdi_for_edit(
        variant.accession.to_string(),
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
    let end_pos = get_end_pos(&variant.loc_edit.location).unwrap_or(start_pos);
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
    let end_pos = get_end_pos(&variant.loc_edit.location).unwrap_or(start_pos);
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
    let end_pos = get_end_pos(&variant.loc_edit.location).unwrap_or(start_pos);
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
        let end = variant
            .loc_edit
            .location
            .end
            .inner()
            .copied()
            .unwrap_or(*start);
        resolve_tx_to_provider_tx(&variant.accession, start, &end, provider)?
    } else {
        let s = tx_pos_for_simple_path(&variant.loc_edit.location, "n")?;
        let e = tx_end_for_simple_path(&variant.loc_edit.location, s, "n")?;
        (s, e)
    };
    emit_spdi_for_edit(
        variant.accession.to_string(),
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
        let end = variant
            .loc_edit
            .location
            .end
            .inner()
            .copied()
            .unwrap_or(*start);
        resolve_rna_to_provider_tx(&variant.accession, start, &end, provider)?
    } else {
        let s = rna_pos_for_simple_path(&variant.loc_edit.location, "r")?;
        let e = rna_end_for_simple_path(&variant.loc_edit.location, s, "r")?;
        (s, e)
    };
    emit_spdi_for_edit(
        variant.accession.to_string(),
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
    let s_u = ensure_positive_tx(s.base, "c", start)?;
    let e_u = ensure_positive_tx(e.base, "c", end)?;
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
/// SPDI representation on the transcript accession.
fn resolve_tx_pos(pos: &TxPos, _transcript: &Transcript) -> Result<u64, ConversionError> {
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
    if pos.utr3 {
        // r.*N is the Nth base of the 3' UTR — anchored at `cds_end`
        // (the last 1-based CDS position), so r.*1 sits at the base
        // immediately after the stop codon (`cds_end + 1`), not at
        // `tx_len + 1` (which is one past the end of the entire
        // mRNA). For non-coding transcripts r.*N has no meaning;
        // surface `MissingReferenceData` rather than fall back to
        // `tx_len` (#390 item 2).
        //
        // TODO(#390-follow-up): the c.*N path goes through
        // `CoordinateMapper::cds_to_tx` which is exon-aware and
        // honors CIGAR insertions in the 3'UTR; this plain
        // `cds_end + base` short-circuits the same arithmetic. For
        // contiguous single-exon 3'UTRs the two agree, but
        // multi-exon 3'UTRs with cdot tx-coordinate gaps will see
        // r.*N and c.*N land at slightly different transcript
        // positions. Route this through the same exon-aware helper.
        if pos.base < 1 {
            return Err(ConversionError::InvalidPosition {
                description: format!("3' UTR position *{} must be >= 1", pos.base),
            });
        }
        let cds_end = transcript
            .cds_end
            .ok_or_else(|| ConversionError::MissingReferenceData {
                description: format!(
                    "r.*{} requires a CDS end on the transcript; non-coding \
                     transcripts have no 3'UTR anchor",
                    pos.base
                ),
            })?;
        return Ok(cds_end.saturating_add(pos.base as u64));
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
/// duplication, or delins (i.e. lacks an explicit deleted sequence), the
/// reference bases for `[start_one_based, end_one_based]` are fetched via
/// the provider so SPDI's mandatory `del` field can be populated. When
/// `provider` is `None`, those cases return [`ConversionError::MissingReferenceData`].
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
            let ins_str = inserted_sequence_to_string(inserted).ok_or_else(|| {
                ConversionError::MissingReferenceData {
                    description: "insertion sequence is not a literal sequence".to_string(),
                }
            })?;
            // SPDI is 0-based interbase: position N is the boundary
            // AFTER 1-based base N (equivalently between 1-based bases
            // N and N+1). For HGVS `g.{start}_{start+1}ins{seq}` the
            // matching SPDI position is `start_one_based` directly
            // (NOT the -1 conversion used for substitution/deletion,
            // which references a specific base). Closes #390 item 1.
            Ok(SpdiVariant::new(
                sequence,
                start_one_based,
                "",
                apply_alphabet(&ins_str, alphabet),
            ))
        }
        NaEdit::Duplication {
            sequence: dup_seq, ..
        } => {
            let dup_str = match dup_seq {
                Some(seq) => sequence_to_string(seq),
                None => match provider {
                    Some(p) => fetch_reference_bases(p, &sequence, start_one_based, end_one_based)?,
                    None => {
                        return Err(ConversionError::MissingReferenceData {
                            description:
                                "duplication sequence not provided; reference data needed to determine duplicated bases"
                                    .to_string(),
                        });
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
                        return Err(ConversionError::MissingReferenceData {
                            description:
                                "deleted sequence not provided; reference data needed to determine deleted bases"
                                    .to_string(),
                        });
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
        } => {
            // Closes #394 item 3. Non-literal delins inserted sequences
            // (any `InsertedSequence` variant other than `Literal`) are
            // a shape SPDI cannot encode, not a missing-reference
            // problem. Matches the sibling arms for Repeat / Inversion /
            // Conversion which already use `UnsupportedEditType` for the
            // same category.
            let ins_str = inserted_sequence_to_string(ins_seq).ok_or_else(|| {
                ConversionError::UnsupportedEditType {
                    description: "delins inserted sequence is not a literal sequence; \
                                  non-literal inserts cannot be encoded as SPDI"
                        .to_string(),
                }
            })?;
            let del_str = match deleted {
                Some(seq) => sequence_to_string(seq),
                None => match provider {
                    Some(p) => fetch_reference_bases(p, &sequence, start_one_based, end_one_based)?,
                    None => {
                        let del_len = end_one_based
                            .saturating_sub(start_one_based)
                            .saturating_add(1) as usize;
                        return Err(ConversionError::MissingReferenceData {
                            description: format!(
                                "Cannot convert delins to SPDI: deleted sequence of length {} is unknown (no reference data)",
                                del_len
                            ),
                        });
                    }
                },
            };
            Ok(SpdiVariant::new(
                sequence,
                spdi_pos,
                apply_alphabet(&del_str, alphabet),
                apply_alphabet(&ins_str, alphabet),
            ))
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
                        return Err(ConversionError::MissingReferenceData {
                            description:
                                "inversion sequence not provided; reference data needed to determine inverted bases"
                                    .to_string(),
                        });
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
            let del_raw = match provider {
                Some(p) => fetch_reference_bases(p, &sequence, start_one_based, end_one_based)?,
                None => {
                    return Err(ConversionError::MissingReferenceData {
                        description:
                            "repeat reference span not provided; reference data needed to determine pre-expansion bases"
                                .to_string(),
                    });
                }
            };
            let del_str = apply_alphabet(&del_raw, alphabet);
            // The HGVS recommendations require the interval to span an
            // integer number of repeat units. Reject non-divisible spans
            // with a clear message rather than silently emitting nonsense.
            if unit_str.is_empty() || !del_str.len().is_multiple_of(unit_str.len()) {
                return Err(ConversionError::InvalidPosition {
                    description: format!(
                        "repeat span {}:{}-{} length {} is not a multiple of unit length {}",
                        sequence,
                        start_one_based,
                        end_one_based,
                        del_str.len(),
                        unit_str.len()
                    ),
                });
            }
            // Verify the reference tract actually consists of repeated copies
            // of the declared unit. A divisible length is necessary but not
            // sufficient (e.g. `ATGCAT` has length 6 but is not `AT[3]`).
            let pre_count = del_str.len() / unit_str.len();
            if del_str != unit_str.repeat(pre_count) {
                return Err(ConversionError::InvalidPosition {
                    description: format!(
                        "repeat span {}:{}-{} does not match repeat unit {}",
                        sequence, start_one_based, end_one_based, unit_str
                    ),
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
/// (SPDI carries no coord-system tag).
fn is_mitochondrial_accession(accession: &Accession) -> bool {
    let prefix: &str = accession.prefix.as_ref();
    let number: &str = accession.number.as_ref();
    // NC_012920 = human GRCh38 mitochondrion (rCRS).
    // NC_001807 = older rCRS draft (deprecated but still seen).
    matches!((prefix, number), ("NC", "012920") | ("NC", "001807"))
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
        let hgvs = parse_hgvs("NR_NONCODING.1:n.*5A>G").unwrap();
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
        // Inversion of a single base 'A' → 'T'
        let provider = make_test_genomic_provider();
        let hgvs = parse_hgvs("NC_000001.11:g.100_100inv").unwrap();
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
}
