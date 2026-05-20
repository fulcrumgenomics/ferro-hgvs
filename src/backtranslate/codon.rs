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

impl Base {
    /// 2-bit packed encoding used by `Codon::index`: A=0, C=1, G=2, T/U=3.
    fn to_index(self) -> u8 {
        match self {
            Base::A => 0,
            Base::C => 1,
            Base::G => 2,
            Base::T => 3,
        }
    }

    fn from_index(i: u8) -> Self {
        match i & 0b11 {
            0 => Base::A,
            1 => Base::C,
            2 => Base::G,
            _ => Base::T,
        }
    }
}

impl Codon {
    /// Create a new codon from three bases.
    pub fn new(b1: Base, b2: Base, b3: Base) -> Self {
        Self([b1, b2, b3])
    }

    /// Parse a codon from a string.
    ///
    /// HGVS codons are always 3 ASCII bases; this fast path avoids the
    /// `Vec<char>` allocation that the previous implementation paid on every
    /// call (1000+ codons per CDS per overlapping transcript inside
    /// `translate_full_cds`).
    pub fn parse(s: &str) -> Option<Self> {
        Self::parse_bytes(s.as_bytes())
    }

    /// Byte-level variant of [`Codon::parse`] for callers that have already
    /// produced a `&[u8]` (e.g. from `chunks_exact(3)` over a CDS sequence).
    /// Skips the `&str` → `&[u8]` indirection so the per-codon loop doesn't
    /// need a UTF-8 validation step that's redundant for ASCII base inputs.
    pub fn parse_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != 3 {
            return None;
        }
        let b1 = Base::from_char(bytes[0] as char)?;
        let b2 = Base::from_char(bytes[1] as char)?;
        let b3 = Base::from_char(bytes[2] as char)?;
        Some(Self([b1, b2, b3]))
    }

    /// Get the three bases.
    pub fn bases(&self) -> &[Base; 3] {
        &self.0
    }

    /// 6-bit packed encoding suitable as an index into a 64-element lookup
    /// table: `b0 << 4 | b1 << 2 | b2` with each base in 2 bits.
    pub fn index(&self) -> u8 {
        (self.0[0].to_index() << 4) | (self.0[1].to_index() << 2) | self.0[2].to_index()
    }

    /// Inverse of [`Codon::index`].
    pub fn from_index(i: u8) -> Self {
        Self([
            Base::from_index(i >> 4),
            Base::from_index(i >> 2),
            Base::from_index(i),
        ])
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
    /// Amino acid -> codons mapping (used by back-translation).
    aa_to_codons: HashMap<AminoAcid, Vec<Codon>>,
    /// Stop codons.
    stop_codons: Vec<Codon>,
    /// Start codon(s).
    start_codons: Vec<Codon>,
    /// 64-entry array indexed by `Codon::index()`. `None` for stop codons.
    /// The hot read path (`amino_acid_for`) hits this directly so codon
    /// translation is a single array load with no hashing.
    aa_lookup: [Option<AminoAcid>; 64],
    /// 64-entry stop mask indexed by `Codon::index()`.
    stop_lookup: [bool; 64],
}

impl CodonTable {
    /// Create the standard genetic code.
    pub fn standard() -> Self {
        let mut aa_to_codons: HashMap<AminoAcid, Vec<Codon>> = HashMap::new();
        let mut stop_codons = Vec::new();
        let mut aa_lookup: [Option<AminoAcid>; 64] = [None; 64];
        let mut stop_lookup: [bool; 64] = [false; 64];

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
            let idx = codon.index() as usize;
            match aa_opt {
                Some(aa) => {
                    aa_lookup[idx] = Some(*aa);
                    aa_to_codons.entry(*aa).or_default().push(codon);
                }
                None => {
                    stop_lookup[idx] = true;
                    stop_codons.push(codon);
                }
            }
        }

        let start_codons = vec![Codon::parse("ATG").unwrap()];

        Self {
            aa_to_codons,
            stop_codons,
            start_codons,
            aa_lookup,
            stop_lookup,
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
    ///
    /// Backed by a 64-element array indexed by [`Codon::index`]; constant time
    /// with no hashing. Returns `None` for stop codons.
    pub fn amino_acid_for(&self, codon: &Codon) -> Option<AminoAcid> {
        self.aa_lookup[codon.index() as usize]
    }

    /// Check if a codon is a stop codon. Constant-time array lookup.
    pub fn is_stop(&self, codon: &Codon) -> bool {
        self.stop_lookup[codon.index() as usize]
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
            Some(AminoAcid::Met)
        );
        assert_eq!(
            table.amino_acid_for(&Codon::parse("TTT").unwrap()),
            Some(AminoAcid::Phe)
        );
        assert_eq!(
            table.amino_acid_for(&Codon::parse("GAG").unwrap()),
            Some(AminoAcid::Glu)
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

    /// Direct-API safety net for the perf rewrite below: every ACGT triple must
    /// resolve to its canonical amino acid via `amino_acid_for`, and only the
    /// three canonical stop codons may report `is_stop = true`.
    #[test]
    fn test_amino_acid_for_all_64_codons() {
        let table = CodonTable::standard();
        let canonical: &[(&str, Option<AminoAcid>)] = &[
            ("TTT", Some(AminoAcid::Phe)),
            ("TTC", Some(AminoAcid::Phe)),
            ("TTA", Some(AminoAcid::Leu)),
            ("TTG", Some(AminoAcid::Leu)),
            ("CTT", Some(AminoAcid::Leu)),
            ("CTC", Some(AminoAcid::Leu)),
            ("CTA", Some(AminoAcid::Leu)),
            ("CTG", Some(AminoAcid::Leu)),
            ("ATT", Some(AminoAcid::Ile)),
            ("ATC", Some(AminoAcid::Ile)),
            ("ATA", Some(AminoAcid::Ile)),
            ("ATG", Some(AminoAcid::Met)),
            ("GTT", Some(AminoAcid::Val)),
            ("GTC", Some(AminoAcid::Val)),
            ("GTA", Some(AminoAcid::Val)),
            ("GTG", Some(AminoAcid::Val)),
            ("TCT", Some(AminoAcid::Ser)),
            ("TCC", Some(AminoAcid::Ser)),
            ("TCA", Some(AminoAcid::Ser)),
            ("TCG", Some(AminoAcid::Ser)),
            ("AGT", Some(AminoAcid::Ser)),
            ("AGC", Some(AminoAcid::Ser)),
            ("CCT", Some(AminoAcid::Pro)),
            ("CCC", Some(AminoAcid::Pro)),
            ("CCA", Some(AminoAcid::Pro)),
            ("CCG", Some(AminoAcid::Pro)),
            ("ACT", Some(AminoAcid::Thr)),
            ("ACC", Some(AminoAcid::Thr)),
            ("ACA", Some(AminoAcid::Thr)),
            ("ACG", Some(AminoAcid::Thr)),
            ("GCT", Some(AminoAcid::Ala)),
            ("GCC", Some(AminoAcid::Ala)),
            ("GCA", Some(AminoAcid::Ala)),
            ("GCG", Some(AminoAcid::Ala)),
            ("TAT", Some(AminoAcid::Tyr)),
            ("TAC", Some(AminoAcid::Tyr)),
            ("TAA", None),
            ("TAG", None),
            ("TGA", None),
            ("CAT", Some(AminoAcid::His)),
            ("CAC", Some(AminoAcid::His)),
            ("CAA", Some(AminoAcid::Gln)),
            ("CAG", Some(AminoAcid::Gln)),
            ("AAT", Some(AminoAcid::Asn)),
            ("AAC", Some(AminoAcid::Asn)),
            ("AAA", Some(AminoAcid::Lys)),
            ("AAG", Some(AminoAcid::Lys)),
            ("GAT", Some(AminoAcid::Asp)),
            ("GAC", Some(AminoAcid::Asp)),
            ("GAA", Some(AminoAcid::Glu)),
            ("GAG", Some(AminoAcid::Glu)),
            ("TGT", Some(AminoAcid::Cys)),
            ("TGC", Some(AminoAcid::Cys)),
            ("TGG", Some(AminoAcid::Trp)),
            ("CGT", Some(AminoAcid::Arg)),
            ("CGC", Some(AminoAcid::Arg)),
            ("CGA", Some(AminoAcid::Arg)),
            ("CGG", Some(AminoAcid::Arg)),
            ("AGA", Some(AminoAcid::Arg)),
            ("AGG", Some(AminoAcid::Arg)),
            ("GGT", Some(AminoAcid::Gly)),
            ("GGC", Some(AminoAcid::Gly)),
            ("GGA", Some(AminoAcid::Gly)),
            ("GGG", Some(AminoAcid::Gly)),
        ];
        assert_eq!(canonical.len(), 64);
        for (codon_str, expected) in canonical {
            let codon = Codon::parse(codon_str).unwrap();
            assert_eq!(
                table.amino_acid_for(&codon),
                *expected,
                "amino_acid_for({codon_str})"
            );
            assert_eq!(
                table.is_stop(&codon),
                expected.is_none(),
                "is_stop({codon_str})"
            );
        }
    }

    /// Cross-check that the array view (`amino_acid_for` / `is_stop`) and the
    /// `aa_to_codons` / `stop_codons` views agree for every 6-bit codon index.
    /// Catches accidental drift if the two are ever populated from different
    /// sources.
    #[test]
    fn test_packed_table_matches_inverse_views() {
        let table = CodonTable::standard();
        for i in 0..64u8 {
            let codon = Codon::from_index(i);
            let via_array = table.amino_acid_for(&codon);
            // Inverse via aa_to_codons: find the (aa, codons) entry that lists
            // this codon; absence ⇒ not a sense codon.
            let via_inverse = table.aa_to_codons.iter().find_map(|(aa, codons)| {
                if codons.contains(&codon) {
                    Some(*aa)
                } else {
                    None
                }
            });
            assert_eq!(via_array, via_inverse, "aa mismatch for {codon}");
            assert_eq!(
                table.is_stop(&codon),
                table.stop_codons.contains(&codon),
                "is_stop mismatch for {codon}"
            );
        }
    }

    #[test]
    fn test_codon_index_roundtrip() {
        for i in 0..64u8 {
            assert_eq!(Codon::from_index(i).index(), i, "index roundtrip {i}");
        }
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
