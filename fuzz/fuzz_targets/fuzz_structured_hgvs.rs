//! Structured fuzz target for HGVS parser
//!
//! Uses the arbitrary crate to generate structured HGVS-like inputs,
//! which is more effective at finding edge cases than pure random bytes.

#![no_main]

use arbitrary::{Arbitrary, Unstructured};
use libfuzzer_sys::fuzz_target;

/// Structured input for generating HGVS-like strings
#[derive(Debug, Arbitrary)]
struct HgvsInput {
    /// Accession type
    accession_type: AccessionType,
    /// Accession number
    accession_num: u32,
    /// Version (0-15)
    version: u8,
    /// Variant type
    variant_type: VariantType,
    /// Position value
    position: i32,
    /// Optional second position for ranges
    position2: Option<i32>,
    /// Edit type
    edit: EditType,
    /// Optional base/amino acid
    base: Option<Base>,
}

#[derive(Debug, Arbitrary)]
enum AccessionType {
    RefSeqNM,
    RefSeqNC,
    RefSeqNR,
    RefSeqNP,
    EnsemblT,
    EnsemblG,
    LRG,
    Custom,
}

#[derive(Debug, Arbitrary)]
enum VariantType {
    Genomic,
    Coding,
    NonCoding,
    Protein,
    Rna,
    Mitochondrial,
}

#[derive(Debug, Arbitrary)]
enum EditType {
    Substitution,
    Deletion,
    Duplication,
    Insertion,
    Delins,
    Repeat,
    Identity,
}

#[derive(Debug, Arbitrary)]
enum Base {
    A,
    T,
    G,
    C,
    U,
    N,
    Ala,
    Arg,
    Asn,
    Asp,
    Cys,
    Gln,
    Glu,
    Gly,
    His,
    Ile,
    Leu,
    Lys,
    Met,
    Phe,
    Pro,
    Ser,
    Thr,
    Trp,
    Tyr,
    Val,
    Ter,
}

impl HgvsInput {
    fn to_hgvs_string(&self) -> String {
        let accession = match self.accession_type {
            AccessionType::RefSeqNM => format!("NM_{:06}.{}", self.accession_num % 1000000, self.version % 16),
            AccessionType::RefSeqNC => format!("NC_{:06}.{}", self.accession_num % 1000000, self.version % 16),
            AccessionType::RefSeqNR => format!("NR_{:06}.{}", self.accession_num % 1000000, self.version % 16),
            AccessionType::RefSeqNP => format!("NP_{:06}.{}", self.accession_num % 1000000, self.version % 16),
            AccessionType::EnsemblT => format!("ENST{:011}", self.accession_num as u64 % 100000000000),
            AccessionType::EnsemblG => format!("ENSG{:011}", self.accession_num as u64 % 100000000000),
            AccessionType::LRG => format!("LRG_{}", self.accession_num % 1000),
            AccessionType::Custom => format!("TEST{}", self.accession_num % 1000),
        };

        let prefix = match self.variant_type {
            VariantType::Genomic => "g",
            VariantType::Coding => "c",
            VariantType::NonCoding => "n",
            VariantType::Protein => "p",
            VariantType::Rna => "r",
            VariantType::Mitochondrial => "m",
        };

        let pos1 = self.position.abs() % 100000;
        let pos_str = if let Some(p2) = self.position2 {
            let pos2 = (p2.abs() % 100000).max(pos1 + 1);
            format!("{}_{}", pos1, pos2)
        } else {
            pos1.to_string()
        };

        let edit_str = match (&self.edit, &self.base, &self.variant_type) {
            (EditType::Substitution, Some(b), VariantType::Protein) => {
                format!("{}>{}", base_to_aa(b), base_to_aa(&Base::Ala))
            }
            (EditType::Substitution, Some(b), _) => {
                format!("{}>{}", base_to_nuc(b), "A")
            }
            (EditType::Substitution, None, _) => "A>G".to_string(),
            (EditType::Deletion, _, _) => "del".to_string(),
            (EditType::Duplication, _, _) => "dup".to_string(),
            (EditType::Insertion, Some(b), VariantType::Protein) => {
                format!("ins{}", base_to_aa(b))
            }
            (EditType::Insertion, _, _) => "insA".to_string(),
            (EditType::Delins, Some(b), VariantType::Protein) => {
                format!("delins{}", base_to_aa(b))
            }
            (EditType::Delins, _, _) => "delinsAA".to_string(),
            (EditType::Repeat, _, _) => "[5]".to_string(),
            (EditType::Identity, _, _) => "=".to_string(),
        };

        format!("{}:{}.{}{}", accession, prefix, pos_str, edit_str)
    }
}

fn base_to_nuc(b: &Base) -> &'static str {
    match b {
        Base::A => "A",
        Base::T => "T",
        Base::G => "G",
        Base::C => "C",
        Base::U => "U",
        Base::N => "N",
        _ => "A",
    }
}

fn base_to_aa(b: &Base) -> &'static str {
    match b {
        Base::Ala => "Ala",
        Base::Arg => "Arg",
        Base::Asn => "Asn",
        Base::Asp => "Asp",
        Base::Cys => "Cys",
        Base::Gln => "Gln",
        Base::Glu => "Glu",
        Base::Gly => "Gly",
        Base::His => "His",
        Base::Ile => "Ile",
        Base::Leu => "Leu",
        Base::Lys => "Lys",
        Base::Met => "Met",
        Base::Phe => "Phe",
        Base::Pro => "Pro",
        Base::Ser => "Ser",
        Base::Thr => "Thr",
        Base::Trp => "Trp",
        Base::Tyr => "Tyr",
        Base::Val => "Val",
        Base::Ter => "Ter",
        _ => "Ala",
    }
}

fuzz_target!(|data: &[u8]| {
    if let Ok(input) = HgvsInput::arbitrary(&mut Unstructured::new(data)) {
        let hgvs_str = input.to_hgvs_string();

        // The parser should never panic
        let _ = ferro_hgvs::parse_hgvs(&hgvs_str);
    }
});
