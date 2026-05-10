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
//! [`OneBasedPos`]: crate::coords::OneBasedPos
//! [`ZeroBasedPos`]: crate::coords::ZeroBasedPos

use super::SpdiVariant;
use crate::coords::{OneBasedPos, ZeroBasedPos};
use crate::error::FerroError;
use crate::hgvs::edit::{InsertedSequence, NaEdit, Sequence};
use crate::hgvs::interval::Interval;
use crate::hgvs::location::GenomePos;
use crate::hgvs::parser::accession::parse_accession;
use crate::hgvs::variant::{GenomeVariant, HgvsVariant, LocEdit};

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

/// Convert a genomic HGVS variant to SPDI format.
///
/// This is a "simple" conversion that works for substitutions where the
/// reference and alternate sequences are explicitly stated in the HGVS.
/// For deletions, insertions, and duplications, use the version that takes
/// a reference provider.
///
/// # Arguments
///
/// * `variant` - A genomic HGVS variant (g. coordinate system)
///
/// # Returns
///
/// * `Ok(SpdiVariant)` - Successfully converted variant
/// * `Err(ConversionError)` - Conversion failed
///
/// # Examples
///
/// ```
/// use ferro_hgvs::spdi::convert::hgvs_to_spdi_simple;
/// use ferro_hgvs::parse_hgvs;
///
/// let hgvs = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
/// let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
/// assert_eq!(spdi.to_string(), "NC_000001.11:12344:A:G");
/// ```
pub fn hgvs_to_spdi_simple(variant: &HgvsVariant) -> Result<SpdiVariant, ConversionError> {
    match variant {
        HgvsVariant::Genome(g) => genome_to_spdi_simple(g),
        _ => Err(ConversionError::UnsupportedVariantType {
            description: format!(
                "only genomic (g.) variants can be directly converted to SPDI, got {}",
                variant.variant_type()
            ),
        }),
    }
}

/// Convert a genomic variant to SPDI (simple conversion).
fn genome_to_spdi_simple(variant: &GenomeVariant) -> Result<SpdiVariant, ConversionError> {
    let sequence = variant.accession.to_string();
    let interval = &variant.loc_edit.location;

    // Get the edit (unwrap from Mu)
    let edit = variant
        .loc_edit
        .edit
        .inner()
        .ok_or_else(|| ConversionError::InvalidPosition {
            description: "cannot convert variant with unknown edit".to_string(),
        })?;

    // Get start position
    let start_pos = get_start_pos(interval).ok_or_else(|| ConversionError::InvalidPosition {
        description: "cannot convert variant with unknown start position".to_string(),
    })?;

    // Validate and wrap in type-safe 1-based position
    let hgvs_pos_ob =
        OneBasedPos::try_new(start_pos).ok_or_else(|| ConversionError::InvalidPosition {
            description: "position 0 is not valid in HGVS".to_string(),
        })?;

    // Convert 1-based HGVS position to 0-based SPDI position using type-safe conversion
    let spdi_pos_zb: ZeroBasedPos = hgvs_pos_ob.to_zero_based();
    let spdi_pos = spdi_pos_zb.value();

    match edit {
        NaEdit::Substitution {
            reference,
            alternative,
        } => {
            // Substitution: g.12345A>G -> seq:12344:A:G
            Ok(SpdiVariant::new(
                sequence,
                spdi_pos,
                reference.to_string(),
                alternative.to_string(),
            ))
        }
        NaEdit::Insertion { sequence: inserted } => {
            // Insertion: g.100_101insATG -> seq:100::ATG
            // For insertion, HGVS uses positions flanking the insertion point
            // SPDI position is the 0-based position where insertion happens
            let ins_str = inserted_sequence_to_string(inserted).ok_or_else(|| {
                ConversionError::MissingReferenceData {
                    description: "insertion sequence is not a literal sequence".to_string(),
                }
            })?;
            Ok(SpdiVariant::new(sequence, spdi_pos, "", ins_str))
        }
        NaEdit::Duplication {
            sequence: dup_seq, ..
        } => {
            // Duplication with sequence: g.100_102dupATG -> seq:102::ATG
            // (insertion of the duplicated sequence after the original)
            if let Some(seq) = dup_seq {
                // Position for dup is at the end of the duplicated region
                let end_pos = get_end_pos(interval).unwrap_or(start_pos);
                // Convert 1-based HGVS end position to 0-based SPDI position
                let end_pos_ob = OneBasedPos::new(end_pos);
                let spdi_end_zb = end_pos_ob.to_zero_based();
                Ok(SpdiVariant::new(
                    sequence,
                    spdi_end_zb.value(),
                    "",
                    sequence_to_string(seq),
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
            // Deletion: g.100_102del -> seq:99:ATG:
            // Need the deleted sequence
            if let Some(seq) = del_seq {
                Ok(SpdiVariant::new(
                    sequence,
                    spdi_pos,
                    sequence_to_string(seq),
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
            // Delins: g.100_102delinsATG -> seq:99:XYZ:ATG
            // We need to know the deleted sequence to create a valid SPDI
            // Without reference data, we cannot determine what was deleted
            // Calculate deletion length from interval (saturating to avoid overflow)
            let end_pos = get_end_pos(interval).unwrap_or(start_pos);
            let del_len = (end_pos - start_pos).saturating_add(1) as usize;

            // Without reference data, we can't know the deleted sequence
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
            // Identity: g.100= or g.100A=
            let ref_base = id_seq.as_ref().map(sequence_to_string).unwrap_or_default();
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
    fn test_hgvs_to_spdi_unsupported_coding() {
        let hgvs = parse_hgvs("NM_000088.3:c.100A>G").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::UnsupportedVariantType { .. })
        ));
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
}
