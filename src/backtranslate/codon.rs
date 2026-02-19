//! Genetic code and codon table.

use crate::hgvs::location::AminoAcid;
use std::collections::HashMap;

/// A single nucleotide base.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Base {
    A,
    T,
    G,
    C,
}

impl Base {
    /// Parse a base from a character.
    pub fn from_char(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            'A' => Some(Base::A),
            'T' | 'U' => Some(Base::T), // U is treated as T
            'G' => Some(Base::G),
            'C' => Some(Base::C),
            _ => None,
        }
    }

    /// Convert to character.
    pub fn to_char(self) -> char {
        match self {
            Base::A => 'A',
            Base::T => 'T',
            Base::G => 'G',
            Base::C => 'C',
        }
    }
}

impl std::fmt::Display for Base {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_char())
    }
}

/// A codon (three nucleotides).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Codon([Base; 3]);

impl Codon {
    /// Create a new codon from three bases.
    pub fn new(b1: Base, b2: Base, b3: Base) -> Self {
        Self([b1, b2, b3])
    }

    /// Parse a codon from a string.
    pub fn parse(s: &str) -> Option<Self> {
        let chars: Vec<char> = s.chars().collect();
        if chars.len() != 3 {
            return None;
        }

        let b1 = Base::from_char(chars[0])?;
        let b2 = Base::from_char(chars[1])?;
        let b3 = Base::from_char(chars[2])?;

        Some(Self([b1, b2, b3]))
    }

    /// Get the three bases.
    pub fn bases(&self) -> &[Base; 3] {
        &self.0
    }
}

impl std::fmt::Display for Codon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}{}", self.0[0], self.0[1], self.0[2])
    }
}

/// Standard genetic code table.
#[derive(Debug, Clone)]
pub struct CodonTable {
    /// Codon to amino acid mapping.
    codon_to_aa: HashMap<Codon, AminoAcid>,
    /// Amino acid to codons mapping.
    aa_to_codons: HashMap<AminoAcid, Vec<Codon>>,
    /// Stop codons.
    stop_codons: Vec<Codon>,
    /// Start codon(s).
    start_codons: Vec<Codon>,
}

impl CodonTable {
    /// Create the standard genetic code.
    pub fn standard() -> Self {
        let mut codon_to_aa = HashMap::new();
        let mut aa_to_codons: HashMap<AminoAcid, Vec<Codon>> = HashMap::new();
        let mut stop_codons = Vec::new();

        // Define the standard genetic code
        let code: &[(&str, Option<AminoAcid>)] = &[
            // Phenylalanine (Phe, F)
            ("TTT", Some(AminoAcid::Phe)),
            ("TTC", Some(AminoAcid::Phe)),
            // Leucine (Leu, L)
            ("TTA", Some(AminoAcid::Leu)),
            ("TTG", Some(AminoAcid::Leu)),
            ("CTT", Some(AminoAcid::Leu)),
            ("CTC", Some(AminoAcid::Leu)),
            ("CTA", Some(AminoAcid::Leu)),
            ("CTG", Some(AminoAcid::Leu)),
            // Isoleucine (Ile, I)
            ("ATT", Some(AminoAcid::Ile)),
            ("ATC", Some(AminoAcid::Ile)),
            ("ATA", Some(AminoAcid::Ile)),
            // Methionine (Met, M) - also start codon
            ("ATG", Some(AminoAcid::Met)),
            // Valine (Val, V)
            ("GTT", Some(AminoAcid::Val)),
            ("GTC", Some(AminoAcid::Val)),
            ("GTA", Some(AminoAcid::Val)),
            ("GTG", Some(AminoAcid::Val)),
            // Serine (Ser, S)
            ("TCT", Some(AminoAcid::Ser)),
            ("TCC", Some(AminoAcid::Ser)),
            ("TCA", Some(AminoAcid::Ser)),
            ("TCG", Some(AminoAcid::Ser)),
            ("AGT", Some(AminoAcid::Ser)),
            ("AGC", Some(AminoAcid::Ser)),
            // Proline (Pro, P)
            ("CCT", Some(AminoAcid::Pro)),
            ("CCC", Some(AminoAcid::Pro)),
            ("CCA", Some(AminoAcid::Pro)),
            ("CCG", Some(AminoAcid::Pro)),
            // Threonine (Thr, T)
            ("ACT", Some(AminoAcid::Thr)),
            ("ACC", Some(AminoAcid::Thr)),
            ("ACA", Some(AminoAcid::Thr)),
            ("ACG", Some(AminoAcid::Thr)),
            // Alanine (Ala, A)
            ("GCT", Some(AminoAcid::Ala)),
            ("GCC", Some(AminoAcid::Ala)),
            ("GCA", Some(AminoAcid::Ala)),
            ("GCG", Some(AminoAcid::Ala)),
            // Tyrosine (Tyr, Y)
            ("TAT", Some(AminoAcid::Tyr)),
            ("TAC", Some(AminoAcid::Tyr)),
            // Stop codons
            ("TAA", None), // Ochre
            ("TAG", None), // Amber
            ("TGA", None), // Opal
            // Histidine (His, H)
            ("CAT", Some(AminoAcid::His)),
            ("CAC", Some(AminoAcid::His)),
            // Glutamine (Gln, Q)
            ("CAA", Some(AminoAcid::Gln)),
            ("CAG", Some(AminoAcid::Gln)),
            // Asparagine (Asn, N)
            ("AAT", Some(AminoAcid::Asn)),
            ("AAC", Some(AminoAcid::Asn)),
            // Lysine (Lys, K)
            ("AAA", Some(AminoAcid::Lys)),
            ("AAG", Some(AminoAcid::Lys)),
            // Aspartic acid (Asp, D)
            ("GAT", Some(AminoAcid::Asp)),
            ("GAC", Some(AminoAcid::Asp)),
            // Glutamic acid (Glu, E)
            ("GAA", Some(AminoAcid::Glu)),
            ("GAG", Some(AminoAcid::Glu)),
            // Cysteine (Cys, C)
            ("TGT", Some(AminoAcid::Cys)),
            ("TGC", Some(AminoAcid::Cys)),
            // Tryptophan (Trp, W)
            ("TGG", Some(AminoAcid::Trp)),
            // Arginine (Arg, R)
            ("CGT", Some(AminoAcid::Arg)),
            ("CGC", Some(AminoAcid::Arg)),
            ("CGA", Some(AminoAcid::Arg)),
            ("CGG", Some(AminoAcid::Arg)),
            ("AGA", Some(AminoAcid::Arg)),
            ("AGG", Some(AminoAcid::Arg)),
            // Glycine (Gly, G)
            ("GGT", Some(AminoAcid::Gly)),
            ("GGC", Some(AminoAcid::Gly)),
            ("GGA", Some(AminoAcid::Gly)),
            ("GGG", Some(AminoAcid::Gly)),
        ];

        for (codon_str, aa_opt) in code {
            let codon = Codon::parse(codon_str).unwrap();

            match aa_opt {
                Some(aa) => {
                    codon_to_aa.insert(codon.clone(), *aa);
                    aa_to_codons.entry(*aa).or_default().push(codon);
                }
                None => {
                    stop_codons.push(codon);
                }
            }
        }

        let start_codons = vec![Codon::parse("ATG").unwrap()];

        Self {
            codon_to_aa,
            aa_to_codons,
            stop_codons,
            start_codons,
        }
    }

    /// Get all codons that encode a given amino acid.
    pub fn codons_for(&self, aa: &AminoAcid) -> &[Codon] {
        self.aa_to_codons
            .get(aa)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Get the amino acid encoded by a codon.
    pub fn amino_acid_for(&self, codon: &Codon) -> Option<&AminoAcid> {
        self.codon_to_aa.get(codon)
    }

    /// Check if a codon is a stop codon.
    pub fn is_stop(&self, codon: &Codon) -> bool {
        self.stop_codons.contains(codon)
    }

    /// Check if a codon is a start codon.
    pub fn is_start(&self, codon: &Codon) -> bool {
        self.start_codons.contains(codon)
    }

    /// Get all stop codons.
    pub fn stop_codons(&self) -> &[Codon] {
        &self.stop_codons
    }

    /// Get all start codons.
    pub fn start_codons(&self) -> &[Codon] {
        &self.start_codons
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_from_char() {
        assert_eq!(Base::from_char('A'), Some(Base::A));
        assert_eq!(Base::from_char('T'), Some(Base::T));
        assert_eq!(Base::from_char('U'), Some(Base::T)); // RNA
        assert_eq!(Base::from_char('G'), Some(Base::G));
        assert_eq!(Base::from_char('C'), Some(Base::C));
        assert_eq!(Base::from_char('N'), None);
    }

    #[test]
    fn test_codon_from_str() {
        let codon = Codon::parse("ATG").unwrap();
        assert_eq!(codon.bases(), &[Base::A, Base::T, Base::G]);
        assert_eq!(codon.to_string(), "ATG");

        // Invalid codons
        assert!(Codon::parse("AT").is_none()); // Too short
        assert!(Codon::parse("ATGC").is_none()); // Too long
        assert!(Codon::parse("ATN").is_none()); // Invalid base
    }

    #[test]
    fn test_standard_code_amino_acids() {
        let table = CodonTable::standard();

        // Check some codons
        assert_eq!(
            table.amino_acid_for(&Codon::parse("ATG").unwrap()),
            Some(&AminoAcid::Met)
        );
        assert_eq!(
            table.amino_acid_for(&Codon::parse("TTT").unwrap()),
            Some(&AminoAcid::Phe)
        );
        assert_eq!(
            table.amino_acid_for(&Codon::parse("GAG").unwrap()),
            Some(&AminoAcid::Glu)
        );
    }

    #[test]
    fn test_standard_code_stop_codons() {
        let table = CodonTable::standard();

        assert!(table.is_stop(&Codon::parse("TAA").unwrap()));
        assert!(table.is_stop(&Codon::parse("TAG").unwrap()));
        assert!(table.is_stop(&Codon::parse("TGA").unwrap()));
        assert!(!table.is_stop(&Codon::parse("ATG").unwrap()));
    }

    #[test]
    fn test_standard_code_start_codons() {
        let table = CodonTable::standard();

        assert!(table.is_start(&Codon::parse("ATG").unwrap()));
        assert!(!table.is_start(&Codon::parse("TTT").unwrap()));
    }

    #[test]
    fn test_codons_for_amino_acid() {
        let table = CodonTable::standard();

        // Phe has 2 codons
        let phe_codons = table.codons_for(&AminoAcid::Phe);
        assert_eq!(phe_codons.len(), 2);

        // Leu has 6 codons
        let leu_codons = table.codons_for(&AminoAcid::Leu);
        assert_eq!(leu_codons.len(), 6);

        // Met has 1 codon
        let met_codons = table.codons_for(&AminoAcid::Met);
        assert_eq!(met_codons.len(), 1);

        // Trp has 1 codon
        let trp_codons = table.codons_for(&AminoAcid::Trp);
        assert_eq!(trp_codons.len(), 1);
    }

    #[test]
    fn test_all_amino_acids_have_codons() {
        let table = CodonTable::standard();

        // Standard 20 amino acids should all have codons
        let amino_acids = [
            AminoAcid::Ala,
            AminoAcid::Arg,
            AminoAcid::Asn,
            AminoAcid::Asp,
            AminoAcid::Cys,
            AminoAcid::Gln,
            AminoAcid::Glu,
            AminoAcid::Gly,
            AminoAcid::His,
            AminoAcid::Ile,
            AminoAcid::Leu,
            AminoAcid::Lys,
            AminoAcid::Met,
            AminoAcid::Phe,
            AminoAcid::Pro,
            AminoAcid::Ser,
            AminoAcid::Thr,
            AminoAcid::Trp,
            AminoAcid::Tyr,
            AminoAcid::Val,
        ];

        for aa in &amino_acids {
            let codons = table.codons_for(aa);
            assert!(
                !codons.is_empty(),
                "Amino acid {:?} should have at least one codon",
                aa
            );
        }
    }

    #[test]
    fn test_codon_count() {
        let table = CodonTable::standard();

        // Total should be 64 codons = 61 sense + 3 stop
        let mut total_sense = 0;
        for aa in &[
            AminoAcid::Ala,
            AminoAcid::Arg,
            AminoAcid::Asn,
            AminoAcid::Asp,
            AminoAcid::Cys,
            AminoAcid::Gln,
            AminoAcid::Glu,
            AminoAcid::Gly,
            AminoAcid::His,
            AminoAcid::Ile,
            AminoAcid::Leu,
            AminoAcid::Lys,
            AminoAcid::Met,
            AminoAcid::Phe,
            AminoAcid::Pro,
            AminoAcid::Ser,
            AminoAcid::Thr,
            AminoAcid::Trp,
            AminoAcid::Tyr,
            AminoAcid::Val,
        ] {
            total_sense += table.codons_for(aa).len();
        }

        assert_eq!(total_sense, 61);
        assert_eq!(table.stop_codons().len(), 3);
    }
}
