//! Helper functions for Python bindings
//!
//! These functions are separated from the PyO3 code so they can be unit tested
//! without requiring the Python runtime.

use crate::hgvs::edit::{InsertedSequence, NaEdit};
use crate::hgvs::interval::Interval;
use crate::hgvs::location::{CdsPos, GenomePos, RnaPos, TxPos};
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::HgvsVariant;
use crate::normalize::ShuffleDirection;

/// Get the variant type as a string
///
/// Returns a human-readable string describing the type of HGVS variant.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::python_helpers::variant_type_str;
/// use ferro_hgvs::parse_hgvs;
///
/// let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
/// assert_eq!(variant_type_str(&variant), "genomic");
///
/// let variant = parse_hgvs("NM_000088.3:c.100A>G").unwrap();
/// assert_eq!(variant_type_str(&variant), "coding");
/// ```
pub fn variant_type_str(variant: &HgvsVariant) -> &'static str {
    match variant {
        HgvsVariant::Genome(_) => "genomic",
        HgvsVariant::Cds(_) => "coding",
        HgvsVariant::Tx(_) => "non_coding",
        HgvsVariant::Protein(_) => "protein",
        HgvsVariant::Rna(_) => "rna",
        HgvsVariant::Mt(_) => "mitochondrial",
        HgvsVariant::Circular(_) => "circular",
        HgvsVariant::RnaFusion(_) => "rna_fusion",
        HgvsVariant::Allele(_) => "allele",
        HgvsVariant::NullAllele => "null_allele",
        HgvsVariant::UnknownAllele => "unknown_allele",
    }
}

/// Get the edit type as a string from an NaEdit
///
/// Returns a human-readable string describing the type of nucleic acid edit.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::python_helpers::na_edit_type_str;
/// use ferro_hgvs::hgvs::edit::{NaEdit, Base};
///
/// let edit = NaEdit::Substitution { reference: Base::A, alternative: Base::G };
/// assert_eq!(na_edit_type_str(&edit), "substitution");
/// ```
pub fn na_edit_type_str(edit: &NaEdit) -> &'static str {
    match edit {
        NaEdit::Substitution { .. } => "substitution",
        NaEdit::SubstitutionNoRef { .. } => "substitution",
        NaEdit::Deletion { .. } => "deletion",
        NaEdit::Duplication { .. } => "duplication",
        NaEdit::DupIns { .. } => "dupins",
        NaEdit::Insertion { .. } => "insertion",
        NaEdit::Delins { .. } => "delins",
        NaEdit::Inversion { .. } => "inversion",
        NaEdit::Repeat { .. } => "repeat",
        NaEdit::MultiRepeat { .. } => "multi_repeat",
        NaEdit::Identity { .. } => "identity",
        NaEdit::Conversion { .. } => "conversion",
        NaEdit::Unknown { .. } => "unknown",
        NaEdit::Methylation { .. } => "methylation",
        NaEdit::CopyNumber { .. } => "copy_number",
        NaEdit::Splice => "splice",
        NaEdit::NoProduct => "no_product",
        NaEdit::PositionOnly => "position_only",
    }
}

/// Get the edit type from debug string representation
///
/// This is a fallback function that determines the edit type by inspecting
/// the Debug representation of an edit. Use `na_edit_type_str` when possible.
pub fn edit_type_from_debug<T: std::fmt::Debug>(edit: &T) -> &'static str {
    let debug = format!("{:?}", edit);
    if debug.contains("Substitution") {
        "substitution"
    } else if debug.contains("Deletion") {
        "deletion"
    } else if debug.contains("DupIns") {
        "dupins"
    } else if debug.contains("Duplication") {
        "duplication"
    } else if debug.contains("Insertion") {
        "insertion"
    } else if debug.contains("Delins") {
        "delins"
    } else if debug.contains("Inversion") {
        "inversion"
    } else if debug.contains("Repeat") {
        "repeat"
    } else if debug.contains("Identity") {
        "identity"
    } else if debug.contains("Conversion") {
        "conversion"
    } else if debug.contains("Methylation") {
        "methylation"
    } else if debug.contains("CopyNumber") {
        "copy_number"
    } else {
        "unknown"
    }
}

/// Parse a direction string into ShuffleDirection
///
/// Accepts various common formats:
/// - 3prime, 3', 3 -> ThreePrime (default)
/// - 5prime, 5', 5 -> FivePrime
///
/// # Examples
///
/// ```
/// use ferro_hgvs::python_helpers::parse_direction;
/// use ferro_hgvs::ShuffleDirection;
///
/// assert!(matches!(parse_direction("3prime"), ShuffleDirection::ThreePrime));
/// assert!(matches!(parse_direction("5'"), ShuffleDirection::FivePrime));
/// ```
pub fn parse_direction(direction: &str) -> ShuffleDirection {
    match direction.to_lowercase().as_str() {
        "5prime" | "5'" | "5" => ShuffleDirection::FivePrime,
        _ => ShuffleDirection::ThreePrime,
    }
}

/// Get the reference accession from a variant
///
/// Returns the accession string for variants that have one, or an error message
/// for variants that don't support accessions.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::python_helpers::get_variant_reference;
/// use ferro_hgvs::parse_hgvs;
///
/// let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
/// assert_eq!(get_variant_reference(&variant).unwrap(), "NC_000001.11");
/// ```
pub fn get_variant_reference(variant: &HgvsVariant) -> Result<String, &'static str> {
    match variant {
        HgvsVariant::Genome(v) => Ok(v.accession.to_string()),
        HgvsVariant::Cds(v) => Ok(v.accession.to_string()),
        HgvsVariant::Tx(v) => Ok(v.accession.to_string()),
        HgvsVariant::Protein(v) => Ok(v.accession.to_string()),
        HgvsVariant::Rna(v) => Ok(v.accession.to_string()),
        HgvsVariant::Mt(v) => Ok(v.accession.to_string()),
        HgvsVariant::Circular(v) => Ok(v.accession.to_string()),
        HgvsVariant::RnaFusion(v) => Ok(v.five_prime.accession.to_string()),
        HgvsVariant::Allele(a) => {
            if let Some(first) = a.variants.first() {
                get_variant_reference(first)
            } else {
                Err("Empty allele")
            }
        }
        HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => {
            Err("No reference for null/unknown allele")
        }
    }
}

/// Get the edit type from a Mu-wrapped NaEdit
///
/// Handles the uncertainty wrapper, returning "unknown" for Mu::Unknown.
pub fn mu_na_edit_type_str(edit: &Mu<NaEdit>) -> &'static str {
    match edit.inner() {
        Some(e) => na_edit_type_str(e),
        None => "unknown",
    }
}

/// Get the edit type string for any HgvsVariant
///
/// Returns a string describing the type of edit in the variant.
pub fn get_variant_edit_type(variant: &HgvsVariant) -> &'static str {
    match variant {
        HgvsVariant::Genome(v) => mu_na_edit_type_str(&v.loc_edit.edit),
        HgvsVariant::Cds(v) => mu_na_edit_type_str(&v.loc_edit.edit),
        HgvsVariant::Tx(v) => mu_na_edit_type_str(&v.loc_edit.edit),
        HgvsVariant::Rna(v) => mu_na_edit_type_str(&v.loc_edit.edit),
        HgvsVariant::Mt(v) => mu_na_edit_type_str(&v.loc_edit.edit),
        HgvsVariant::Circular(v) => mu_na_edit_type_str(&v.loc_edit.edit),
        HgvsVariant::Protein(_) => "protein",
        HgvsVariant::RnaFusion(_) => "fusion",
        HgvsVariant::Allele(_) => "allele",
        HgvsVariant::NullAllele => "null",
        HgvsVariant::UnknownAllele => "unknown",
    }
}

// ============== Position Accessors ==============

/// Extract the start base position from a genome interval
fn genome_interval_start(interval: &Interval<GenomePos>) -> Option<i64> {
    interval.start.inner().map(|pos| pos.base as i64)
}

/// Extract the end base position from a genome interval
fn genome_interval_end(interval: &Interval<GenomePos>) -> Option<i64> {
    interval.end.inner().map(|pos| pos.base as i64)
}

/// Extract the start base position from a CDS interval
fn cds_interval_start(interval: &Interval<CdsPos>) -> Option<i64> {
    interval.start.inner().map(|pos| pos.base)
}

/// Extract the end base position from a CDS interval
fn cds_interval_end(interval: &Interval<CdsPos>) -> Option<i64> {
    interval.end.inner().map(|pos| pos.base)
}

/// Extract the start base position from a transcript interval
fn tx_interval_start(interval: &Interval<TxPos>) -> Option<i64> {
    interval.start.inner().map(|pos| pos.base)
}

/// Extract the end base position from a transcript interval
fn tx_interval_end(interval: &Interval<TxPos>) -> Option<i64> {
    interval.end.inner().map(|pos| pos.base)
}

/// Extract the start base position from an RNA interval
fn rna_interval_start(interval: &Interval<RnaPos>) -> Option<i64> {
    interval.start.inner().map(|pos| pos.base)
}

/// Extract the end base position from an RNA interval
fn rna_interval_end(interval: &Interval<RnaPos>) -> Option<i64> {
    interval.end.inner().map(|pos| pos.base)
}

/// Get the start position (base only, no offset) for any variant.
///
/// Returns the 1-based start position for genomic, coding, non-coding, RNA, and
/// mitochondrial variants. For single-element alleles, delegates to that
/// sub-variant. Returns None for protein, RNA fusion, null/unknown allele
/// variants, and alleles with multiple sub-variants (whose start is ambiguous).
///
/// Note: 5' UTR (`c.-5A>G`) and 3' UTR (`c.*5A>G`) positions are returned as
/// raw base values and are indistinguishable from CDS positions at the same
/// numeric value. Callers needing to distinguish these cases should inspect
/// the variant string or AST directly.
pub fn get_variant_start(variant: &HgvsVariant) -> Option<i64> {
    match variant {
        HgvsVariant::Genome(v) => genome_interval_start(&v.loc_edit.location),
        HgvsVariant::Cds(v) => cds_interval_start(&v.loc_edit.location),
        HgvsVariant::Tx(v) => tx_interval_start(&v.loc_edit.location),
        HgvsVariant::Rna(v) => rna_interval_start(&v.loc_edit.location),
        HgvsVariant::Mt(v) => genome_interval_start(&v.loc_edit.location),
        HgvsVariant::Circular(v) => genome_interval_start(&v.loc_edit.location),
        HgvsVariant::Allele(a) => match a.variants.as_slice() {
            [single] => get_variant_start(single),
            _ => None,
        },
        HgvsVariant::Protein(_)
        | HgvsVariant::RnaFusion(_)
        | HgvsVariant::NullAllele
        | HgvsVariant::UnknownAllele => None,
    }
}

/// Get the end position (base only, no offset) for any variant.
///
/// Returns the 1-based end position (inclusive). For point variants, end == start.
/// For single-element alleles, delegates to that sub-variant. Returns None for
/// protein, RNA fusion, null/unknown allele variants, and alleles with multiple
/// sub-variants (whose end is ambiguous).
pub fn get_variant_end(variant: &HgvsVariant) -> Option<i64> {
    match variant {
        HgvsVariant::Genome(v) => genome_interval_end(&v.loc_edit.location),
        HgvsVariant::Cds(v) => cds_interval_end(&v.loc_edit.location),
        HgvsVariant::Tx(v) => tx_interval_end(&v.loc_edit.location),
        HgvsVariant::Rna(v) => rna_interval_end(&v.loc_edit.location),
        HgvsVariant::Mt(v) => genome_interval_end(&v.loc_edit.location),
        HgvsVariant::Circular(v) => genome_interval_end(&v.loc_edit.location),
        HgvsVariant::Allele(a) => match a.variants.as_slice() {
            [single] => get_variant_end(single),
            _ => None,
        },
        HgvsVariant::Protein(_)
        | HgvsVariant::RnaFusion(_)
        | HgvsVariant::NullAllele
        | HgvsVariant::UnknownAllele => None,
    }
}

/// Get the intronic offset of the start position for CDS, transcript, and RNA
/// variants.
///
/// For `c.93+1G>T`, returns `Some(1)`. For exonic positions (no offset), returns
/// `None`. Returns `None` for variant types without intronic offsets (genomic,
/// mitochondrial, circular, protein, fusion, allele, null/unknown).
pub fn get_variant_offset(variant: &HgvsVariant) -> Option<i64> {
    match variant {
        HgvsVariant::Cds(v) => v.loc_edit.location.start.inner().and_then(|pos| pos.offset),
        HgvsVariant::Tx(v) => v.loc_edit.location.start.inner().and_then(|pos| pos.offset),
        HgvsVariant::Rna(v) => v.loc_edit.location.start.inner().and_then(|pos| pos.offset),
        _ => None,
    }
}

// ============== Edit Accessors ==============

/// Get the NaEdit from a variant, if it has one.
fn get_na_edit(variant: &HgvsVariant) -> Option<&NaEdit> {
    match variant {
        HgvsVariant::Genome(v) => v.loc_edit.edit.inner(),
        HgvsVariant::Cds(v) => v.loc_edit.edit.inner(),
        HgvsVariant::Tx(v) => v.loc_edit.edit.inner(),
        HgvsVariant::Rna(v) => v.loc_edit.edit.inner(),
        HgvsVariant::Mt(v) => v.loc_edit.edit.inner(),
        HgvsVariant::Circular(v) => v.loc_edit.edit.inner(),
        _ => None,
    }
}

/// Get the substitution reference and alternative bases.
///
/// Returns `Some(('A', 'G'))` for a substitution, `None` for other edit types.
pub fn get_substitution_bases(variant: &HgvsVariant) -> Option<(char, char)> {
    let edit = get_na_edit(variant)?;
    match edit {
        NaEdit::Substitution {
            reference,
            alternative,
        } => {
            let ref_char = (*reference as u8) as char;
            let alt_char = (*alternative as u8) as char;
            Some((ref_char, alt_char))
        }
        _ => None,
    }
}

/// Check if the variant represents an identity (no-change, `=`) edit.
pub fn is_identity(variant: &HgvsVariant) -> bool {
    matches!(get_na_edit(variant), Some(NaEdit::Identity { .. }))
}

/// Compute the span of the affected region from the interval.
///
/// Used for deletion/duplication/delins length calculations; not used for
/// insertions, whose length comes from the inserted sequence rather than the
/// flanking positions.
///
/// Returns `None` for circular (`o.`) variants where `end < start`, which
/// represent origin-crossing intervals whose span depends on the contig length
/// (not available here).
fn compute_span(variant: &HgvsVariant) -> Option<i64> {
    let start = get_variant_start(variant)?;
    let end = get_variant_end(variant)?;
    if matches!(variant, HgvsVariant::Circular(_)) && end < start {
        return None;
    }
    Some((end - start) + 1)
}

/// Compute the net indel length (bases gained or lost) for a variant.
///
/// - Substitution: 0 (same number of bases)
/// - Deletion: -(span)
/// - Insertion: +(inserted length), or None if length is unknowable
/// - Delins: inserted_length - span, or None if inserted length is unknowable
/// - Duplication: +(span)
/// - Inversion: 0
/// - Identity: 0
///
/// Returns None if the indel length cannot be determined (e.g., uncertain inserted
/// sequence length, protein variants, unknown edits).
pub fn get_indel_length(variant: &HgvsVariant) -> Option<i64> {
    let edit = get_na_edit(variant)?;
    match edit {
        NaEdit::Substitution { .. } | NaEdit::SubstitutionNoRef { .. } => Some(0),
        NaEdit::Inversion { .. } => Some(0),
        NaEdit::Identity { .. } => Some(0),
        NaEdit::Deletion { .. } => {
            let span = compute_span(variant)?;
            Some(-span)
        }
        NaEdit::Duplication { .. } => {
            let span = compute_span(variant)?;
            Some(span)
        }
        NaEdit::Insertion { sequence } => {
            let ins_len = inserted_sequence_len(sequence)?;
            Some(ins_len)
        }
        NaEdit::Delins { sequence } => {
            let span = compute_span(variant)?;
            let ins_len = inserted_sequence_len(sequence)?;
            Some(ins_len - span)
        }
        _ => None,
    }
}

/// Get the length of an inserted sequence, if deterministic.
fn inserted_sequence_len(seq: &InsertedSequence) -> Option<i64> {
    seq.len().map(|n| n as i64)
}

/// Check if a variant causes a frameshift (indel_length % 3 != 0).
///
/// Returns false if the indel length is 0 or cannot be determined.
pub fn is_frameshift(variant: &HgvsVariant) -> bool {
    match get_indel_length(variant) {
        Some(len) if len != 0 => len % 3 != 0,
        _ => false,
    }
}

/// Get the number of sub-variants in an allele, or 1 for simple variants.
pub fn get_num_variants(variant: &HgvsVariant) -> usize {
    match variant {
        HgvsVariant::Allele(a) => a.variants.len(),
        _ => 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::edit::{Base, MethylationStatus};
    use crate::parse_hgvs;

    // ===== variant_type_str Tests =====

    #[test]
    fn test_variant_type_str_genomic() {
        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        assert_eq!(variant_type_str(&variant), "genomic");
    }

    #[test]
    fn test_variant_type_str_coding() {
        let variant = parse_hgvs("NM_000088.3:c.100A>G").unwrap();
        assert_eq!(variant_type_str(&variant), "coding");
    }

    #[test]
    fn test_variant_type_str_non_coding() {
        let variant = parse_hgvs("NR_000001.1:n.100A>G").unwrap();
        assert_eq!(variant_type_str(&variant), "non_coding");
    }

    #[test]
    fn test_variant_type_str_protein() {
        let variant = parse_hgvs("NP_000001.1:p.Val100Glu").unwrap();
        assert_eq!(variant_type_str(&variant), "protein");
    }

    #[test]
    fn test_variant_type_str_rna() {
        let variant = parse_hgvs("NM_000088.3:r.100a>g").unwrap();
        assert_eq!(variant_type_str(&variant), "rna");
    }

    #[test]
    fn test_variant_type_str_mitochondrial() {
        let variant = parse_hgvs("NC_012920.1:m.100A>G").unwrap();
        assert_eq!(variant_type_str(&variant), "mitochondrial");
    }

    #[test]
    fn test_variant_type_str_circular() {
        let variant = parse_hgvs("NC_001416.1:o.100A>G").unwrap();
        assert_eq!(variant_type_str(&variant), "circular");
    }

    #[test]
    fn test_variant_type_str_null_allele() {
        let variant = HgvsVariant::NullAllele;
        assert_eq!(variant_type_str(&variant), "null_allele");
    }

    #[test]
    fn test_variant_type_str_unknown_allele() {
        let variant = HgvsVariant::UnknownAllele;
        assert_eq!(variant_type_str(&variant), "unknown_allele");
    }

    // ===== na_edit_type_str Tests =====

    #[test]
    fn test_na_edit_type_str_substitution() {
        let edit = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        };
        assert_eq!(na_edit_type_str(&edit), "substitution");
    }

    #[test]
    fn test_na_edit_type_str_deletion() {
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        assert_eq!(na_edit_type_str(&edit), "deletion");
    }

    #[test]
    fn test_na_edit_type_str_duplication() {
        let edit = NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: None,
        };
        assert_eq!(na_edit_type_str(&edit), "duplication");
    }

    #[test]
    fn test_na_edit_type_str_insertion() {
        use crate::hgvs::edit::{InsertedSequence, Sequence};
        use std::str::FromStr;
        let edit = NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::from_str("ATG").unwrap()),
        };
        assert_eq!(na_edit_type_str(&edit), "insertion");
    }

    #[test]
    fn test_na_edit_type_str_delins() {
        use crate::hgvs::edit::{InsertedSequence, Sequence};
        use std::str::FromStr;
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::from_str("ATG").unwrap()),
        };
        assert_eq!(na_edit_type_str(&edit), "delins");
    }

    #[test]
    fn test_na_edit_type_str_inversion() {
        let edit = NaEdit::Inversion {
            sequence: None,
            length: None,
        };
        assert_eq!(na_edit_type_str(&edit), "inversion");
    }

    #[test]
    fn test_na_edit_type_str_repeat() {
        use crate::hgvs::edit::RepeatCount;
        let edit = NaEdit::Repeat {
            sequence: None,
            count: RepeatCount::Exact(10),
            additional_counts: Vec::new(),
            trailing: None,
        };
        assert_eq!(na_edit_type_str(&edit), "repeat");
    }

    #[test]
    fn test_na_edit_type_str_identity() {
        let edit = NaEdit::Identity {
            sequence: None,
            whole_entity: false,
        };
        assert_eq!(na_edit_type_str(&edit), "identity");
    }

    #[test]
    fn test_na_edit_type_str_unknown() {
        let edit = NaEdit::Unknown {
            whole_entity: false,
        };
        assert_eq!(na_edit_type_str(&edit), "unknown");
    }

    #[test]
    fn test_na_edit_type_str_methylation() {
        let edit = NaEdit::Methylation {
            status: MethylationStatus::GainOfMethylation,
        };
        assert_eq!(na_edit_type_str(&edit), "methylation");
    }

    // ===== edit_type_from_debug Tests =====

    #[test]
    fn test_edit_type_from_debug_substitution() {
        #[derive(Debug)]
        struct TestSubstitution;
        assert_eq!(edit_type_from_debug(&TestSubstitution), "substitution");
    }

    #[test]
    fn test_edit_type_from_debug_deletion() {
        #[derive(Debug)]
        struct TestDeletion;
        assert_eq!(edit_type_from_debug(&TestDeletion), "deletion");
    }

    #[test]
    fn test_edit_type_from_debug_unknown() {
        #[derive(Debug)]
        struct TestOther;
        assert_eq!(edit_type_from_debug(&TestOther), "unknown");
    }

    // ===== parse_direction Tests =====

    #[test]
    fn test_parse_direction_three_prime() {
        assert!(matches!(
            parse_direction("3prime"),
            ShuffleDirection::ThreePrime
        ));
        assert!(matches!(
            parse_direction("3'"),
            ShuffleDirection::ThreePrime
        ));
        assert!(matches!(parse_direction("3"), ShuffleDirection::ThreePrime));
    }

    #[test]
    fn test_parse_direction_five_prime() {
        assert!(matches!(
            parse_direction("5prime"),
            ShuffleDirection::FivePrime
        ));
        assert!(matches!(parse_direction("5'"), ShuffleDirection::FivePrime));
        assert!(matches!(parse_direction("5"), ShuffleDirection::FivePrime));
    }

    #[test]
    fn test_parse_direction_default() {
        assert!(matches!(
            parse_direction("unknown"),
            ShuffleDirection::ThreePrime
        ));
        assert!(matches!(parse_direction(""), ShuffleDirection::ThreePrime));
    }

    #[test]
    fn test_parse_direction_case_insensitive() {
        assert!(matches!(
            parse_direction("5PRIME"),
            ShuffleDirection::FivePrime
        ));
        assert!(matches!(
            parse_direction("5Prime"),
            ShuffleDirection::FivePrime
        ));
    }

    // ===== get_variant_reference Tests =====

    #[test]
    fn test_get_variant_reference_genomic() {
        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        assert_eq!(get_variant_reference(&variant).unwrap(), "NC_000001.11");
    }

    #[test]
    fn test_get_variant_reference_coding() {
        let variant = parse_hgvs("NM_000088.3:c.100A>G").unwrap();
        assert_eq!(get_variant_reference(&variant).unwrap(), "NM_000088.3");
    }

    #[test]
    fn test_get_variant_reference_protein() {
        let variant = parse_hgvs("NP_000001.1:p.Val100Glu").unwrap();
        assert_eq!(get_variant_reference(&variant).unwrap(), "NP_000001.1");
    }

    #[test]
    fn test_get_variant_reference_null_allele() {
        let variant = HgvsVariant::NullAllele;
        assert!(get_variant_reference(&variant).is_err());
    }

    #[test]
    fn test_get_variant_reference_unknown_allele() {
        let variant = HgvsVariant::UnknownAllele;
        assert!(get_variant_reference(&variant).is_err());
    }

    // ===== get_variant_edit_type Tests =====

    #[test]
    fn test_get_variant_edit_type_substitution() {
        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        assert_eq!(get_variant_edit_type(&variant), "substitution");
    }

    #[test]
    fn test_get_variant_edit_type_deletion() {
        let variant = parse_hgvs("NC_000001.11:g.12345del").unwrap();
        assert_eq!(get_variant_edit_type(&variant), "deletion");
    }

    #[test]
    fn test_get_variant_edit_type_duplication() {
        let variant = parse_hgvs("NC_000001.11:g.12345dup").unwrap();
        assert_eq!(get_variant_edit_type(&variant), "duplication");
    }

    #[test]
    fn test_get_variant_edit_type_insertion() {
        let variant = parse_hgvs("NC_000001.11:g.12345_12346insATG").unwrap();
        assert_eq!(get_variant_edit_type(&variant), "insertion");
    }

    #[test]
    fn test_get_variant_edit_type_delins() {
        let variant = parse_hgvs("NC_000001.11:g.12345_12350delinsATG").unwrap();
        assert_eq!(get_variant_edit_type(&variant), "delins");
    }

    #[test]
    fn test_get_variant_edit_type_inversion() {
        let variant = parse_hgvs("NC_000001.11:g.12345_12350inv").unwrap();
        assert_eq!(get_variant_edit_type(&variant), "inversion");
    }

    #[test]
    fn test_get_variant_edit_type_protein() {
        let variant = parse_hgvs("NP_000001.1:p.Val100Glu").unwrap();
        assert_eq!(get_variant_edit_type(&variant), "protein");
    }

    #[test]
    fn test_get_variant_edit_type_null() {
        let variant = HgvsVariant::NullAllele;
        assert_eq!(get_variant_edit_type(&variant), "null");
    }

    #[test]
    fn test_get_variant_edit_type_unknown() {
        let variant = HgvsVariant::UnknownAllele;
        assert_eq!(get_variant_edit_type(&variant), "unknown");
    }

    // ===== Position Accessor Tests =====

    #[test]
    fn test_get_variant_start_genomic() {
        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        assert_eq!(get_variant_start(&variant), Some(12345));
    }

    #[test]
    fn test_get_variant_start_coding() {
        let variant = parse_hgvs("NM_000088.3:c.100A>G").unwrap();
        assert_eq!(get_variant_start(&variant), Some(100));
    }

    #[test]
    fn test_get_variant_start_coding_intronic() {
        let variant = parse_hgvs("NM_000088.3:c.93+1G>T").unwrap();
        assert_eq!(get_variant_start(&variant), Some(93));
    }

    #[test]
    fn test_get_variant_end_range() {
        let variant = parse_hgvs("NC_000001.11:g.12345_12350del").unwrap();
        assert_eq!(get_variant_start(&variant), Some(12345));
        assert_eq!(get_variant_end(&variant), Some(12350));
    }

    #[test]
    fn test_get_variant_end_point() {
        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        assert_eq!(get_variant_end(&variant), Some(12345));
    }

    #[test]
    fn test_get_variant_offset_intronic() {
        let variant = parse_hgvs("NM_000088.3:c.93+1G>T").unwrap();
        assert_eq!(get_variant_offset(&variant), Some(1));
    }

    #[test]
    fn test_get_variant_offset_exonic() {
        let variant = parse_hgvs("NM_000088.3:c.100A>G").unwrap();
        assert_eq!(get_variant_offset(&variant), None);
    }

    #[test]
    fn test_get_variant_offset_genomic() {
        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        assert_eq!(get_variant_offset(&variant), None);
    }

    #[test]
    fn test_get_variant_start_allele_multi() {
        // Multi-sub-variant alleles have ambiguous start; return None.
        let variant = parse_hgvs("NM_000088.3:c.[100A>G;200C>T]").unwrap();
        assert_eq!(get_variant_start(&variant), None);
        assert_eq!(get_variant_end(&variant), None);
    }

    #[test]
    fn test_get_variant_start_allele_single() {
        // Single-sub-variant alleles delegate to that sub-variant.
        let variant = parse_hgvs("NM_000088.3:c.[100A>G]").unwrap();
        assert_eq!(get_variant_start(&variant), Some(100));
        assert_eq!(get_variant_end(&variant), Some(100));
    }

    #[test]
    fn test_get_variant_start_noncoding() {
        let variant = parse_hgvs("NR_046018.2:n.100A>G").unwrap();
        assert_eq!(get_variant_start(&variant), Some(100));
        assert_eq!(get_variant_end(&variant), Some(100));
    }

    #[test]
    fn test_get_variant_offset_noncoding_intronic() {
        let variant = parse_hgvs("NR_046018.2:n.100+5A>G").unwrap();
        assert_eq!(get_variant_start(&variant), Some(100));
        assert_eq!(get_variant_offset(&variant), Some(5));
    }

    #[test]
    fn test_get_variant_start_rna() {
        let variant = parse_hgvs("NM_000088.3:r.100a>g").unwrap();
        assert_eq!(get_variant_start(&variant), Some(100));
        assert_eq!(get_variant_end(&variant), Some(100));
    }

    #[test]
    fn test_get_variant_offset_rna_intronic() {
        let variant = parse_hgvs("NM_000088.3:r.100+5a>g").unwrap();
        assert_eq!(get_variant_offset(&variant), Some(5));
    }

    #[test]
    fn test_get_variant_start_utr3() {
        // c.*5A>G (3' UTR) returns the same base value as c.5A>G — the `*`
        // marker is not exposed by start/end.
        let utr = parse_hgvs("NM_000088.3:c.*5A>G").unwrap();
        let cds = parse_hgvs("NM_000088.3:c.5A>G").unwrap();
        assert_eq!(get_variant_start(&utr), Some(5));
        assert_eq!(get_variant_start(&cds), Some(5));
    }

    #[test]
    fn test_get_variant_start_protein_is_none() {
        let variant = parse_hgvs("NP_000001.1:p.Val100Glu").unwrap();
        assert_eq!(get_variant_start(&variant), None);
        assert_eq!(get_variant_end(&variant), None);
        assert_eq!(get_variant_offset(&variant), None);
        assert_eq!(get_indel_length(&variant), None);
    }

    #[test]
    fn test_get_variant_start_null_allele_is_none() {
        let variant = HgvsVariant::NullAllele;
        assert_eq!(get_variant_start(&variant), None);
        assert_eq!(get_variant_end(&variant), None);
    }

    #[test]
    fn test_get_variant_start_unknown_allele_is_none() {
        let variant = HgvsVariant::UnknownAllele;
        assert_eq!(get_variant_start(&variant), None);
        assert_eq!(get_variant_end(&variant), None);
    }

    // ===== Edit Accessor Tests =====

    #[test]
    fn test_get_substitution_bases() {
        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        assert_eq!(get_substitution_bases(&variant), Some(('A', 'G')));
    }

    #[test]
    fn test_get_substitution_bases_not_sub() {
        let variant = parse_hgvs("NC_000001.11:g.12345del").unwrap();
        assert_eq!(get_substitution_bases(&variant), None);
    }

    #[test]
    fn test_is_identity_true() {
        let variant = parse_hgvs("NM_000088.3:c.100=").unwrap();
        assert!(is_identity(&variant));
    }

    #[test]
    fn test_is_identity_false() {
        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        assert!(!is_identity(&variant));
    }

    // ===== Indel Length Tests =====

    #[test]
    fn test_indel_length_substitution() {
        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        assert_eq!(get_indel_length(&variant), Some(0));
    }

    #[test]
    fn test_indel_length_deletion_point() {
        let variant = parse_hgvs("NC_000001.11:g.12345del").unwrap();
        assert_eq!(get_indel_length(&variant), Some(-1));
    }

    #[test]
    fn test_indel_length_deletion_range() {
        let variant = parse_hgvs("NC_000001.11:g.12345_12350del").unwrap();
        assert_eq!(get_indel_length(&variant), Some(-6));
    }

    #[test]
    fn test_indel_length_insertion() {
        let variant = parse_hgvs("NC_000001.11:g.12345_12346insATG").unwrap();
        assert_eq!(get_indel_length(&variant), Some(3));
    }

    #[test]
    fn test_indel_length_delins() {
        // delins replaces 6 bases with 4 bases → net -2
        let variant = parse_hgvs("NC_000001.11:g.12345_12350delinsATTT").unwrap();
        assert_eq!(get_indel_length(&variant), Some(-2));
    }

    #[test]
    fn test_indel_length_duplication_point() {
        let variant = parse_hgvs("NC_000001.11:g.12345dup").unwrap();
        assert_eq!(get_indel_length(&variant), Some(1));
    }

    #[test]
    fn test_indel_length_duplication_range() {
        let variant = parse_hgvs("NC_000001.11:g.12345_12347dup").unwrap();
        assert_eq!(get_indel_length(&variant), Some(3));
    }

    #[test]
    fn test_indel_length_inversion() {
        let variant = parse_hgvs("NC_000001.11:g.12345_12350inv").unwrap();
        assert_eq!(get_indel_length(&variant), Some(0));
    }

    #[test]
    fn test_indel_length_insertion_count_is_none() {
        // ins10 specifies a count, not a literal sequence — but the count is
        // a known length, so we still report it.
        let variant = parse_hgvs("NC_000001.11:g.12345_12346ins10").unwrap();
        assert_eq!(get_indel_length(&variant), Some(10));
    }

    #[test]
    fn test_indel_length_insertion_range_is_none() {
        // ins(10_20) has unknown exact length — must return None.
        let variant = parse_hgvs("NC_000001.11:g.12345_12346ins(10_20)").unwrap();
        assert_eq!(get_indel_length(&variant), None);
    }

    #[test]
    fn test_indel_length_circular_normal() {
        // Non-wrap-around circular variants compute span normally.
        let variant = parse_hgvs("NC_001416.1:o.100_105del").unwrap();
        assert_eq!(get_indel_length(&variant), Some(-6));
    }

    #[test]
    fn test_indel_length_circular_wraparound_is_none() {
        // Origin-crossing circular intervals (end < start) have an unknowable
        // span without contig length, so indel_length must be None. The parser
        // currently rejects wrap-around at the syntax level, but this guard is
        // defensive for programmatic construction or future parser changes.
        let mut variant = parse_hgvs("NC_001416.1:o.197_4344del").unwrap();
        if let HgvsVariant::Circular(ref mut v) = variant {
            std::mem::swap(&mut v.loc_edit.location.start, &mut v.loc_edit.location.end);
        } else {
            panic!("expected Circular variant");
        }
        assert_eq!(get_variant_start(&variant), Some(4344));
        assert_eq!(get_variant_end(&variant), Some(197));
        assert_eq!(get_indel_length(&variant), None);
        assert!(!is_frameshift(&variant));
    }

    // ===== Frameshift Tests =====

    #[test]
    fn test_is_frameshift_true() {
        // 1bp deletion is frameshift
        let variant = parse_hgvs("NC_000001.11:g.12345del").unwrap();
        assert!(is_frameshift(&variant));
    }

    #[test]
    fn test_is_frameshift_false_in_frame() {
        // 3bp deletion is in-frame
        let variant = parse_hgvs("NC_000001.11:g.12345_12347del").unwrap();
        assert!(!is_frameshift(&variant));
    }

    #[test]
    fn test_is_frameshift_false_substitution() {
        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        assert!(!is_frameshift(&variant));
    }

    // ===== Num Variants Tests =====

    #[test]
    fn test_num_variants_simple() {
        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        assert_eq!(get_num_variants(&variant), 1);
    }

    #[test]
    fn test_num_variants_allele() {
        let variant = parse_hgvs("NM_000088.3:c.[100A>G;200C>T]").unwrap();
        assert_eq!(get_num_variants(&variant), 2);
    }
}
