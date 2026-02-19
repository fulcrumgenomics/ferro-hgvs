//! Helper functions for Python bindings
//!
//! These functions are separated from the PyO3 code so they can be unit tested
//! without requiring the Python runtime.

use crate::hgvs::edit::NaEdit;
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
}
