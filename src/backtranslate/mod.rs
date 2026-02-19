//! Backtranslation from protein to DNA variants.
//!
//! This module provides functionality to reverse-translate protein changes
//! to all possible DNA changes that could produce them.
//!
//! # Example
//!
//! ```
//! use ferro_hgvs::backtranslate::{Backtranslator, CodonTable};
//! use ferro_hgvs::hgvs::location::AminoAcid;
//!
//! let bt = Backtranslator::new(CodonTable::standard());
//!
//! // Find all single-nucleotide changes that could cause Leu -> Phe
//! let leu = AminoAcid::Leu;
//! let phe = AminoAcid::Phe;
//! let changes = bt.backtranslate_substitution(&leu, &phe);
//!
//! for change in changes {
//!     println!("{} -> {}", change.ref_codon, change.alt_codon);
//! }
//! ```

pub mod codon;

pub use codon::{Base, Codon, CodonTable};

use crate::hgvs::location::AminoAcid;

/// A codon change representing a DNA variant.
#[derive(Debug, Clone, PartialEq)]
pub struct CodonChange {
    /// Reference codon (if known).
    pub ref_codon: Codon,
    /// Alternate codon.
    pub alt_codon: Codon,
    /// Position(s) in codon that changed (1-indexed).
    pub changed_positions: Vec<u8>,
    /// The nucleotide changes: (position, ref_base, alt_base).
    pub nucleotide_changes: Vec<(u8, Base, Base)>,
}

impl std::fmt::Display for CodonChange {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} -> {}", self.ref_codon, self.alt_codon)
    }
}

/// Result of backtranslation.
#[derive(Debug, Clone)]
pub struct BacktranslationResult {
    /// Original protein variant description.
    pub protein_variant: String,
    /// All possible DNA changes.
    pub dna_variants: Vec<CodonChange>,
}

/// Backtranslation engine.
#[derive(Debug, Clone)]
pub struct Backtranslator {
    codon_table: CodonTable,
}

impl Backtranslator {
    /// Create a new backtranslator with the given codon table.
    pub fn new(codon_table: CodonTable) -> Self {
        Self { codon_table }
    }

    /// Create a backtranslator with the standard genetic code.
    pub fn standard() -> Self {
        Self::new(CodonTable::standard())
    }

    /// Backtranslate an amino acid substitution.
    ///
    /// Returns all possible codon changes that could cause the given
    /// amino acid substitution. Only returns changes that differ by
    /// a single nucleotide (most common case).
    pub fn backtranslate_substitution(
        &self,
        ref_aa: &AminoAcid,
        alt_aa: &AminoAcid,
    ) -> Vec<CodonChange> {
        let ref_codons = self.codon_table.codons_for(ref_aa);
        let alt_codons = self.codon_table.codons_for(alt_aa);

        let mut results = Vec::new();

        for ref_codon in ref_codons {
            for alt_codon in alt_codons {
                if let Some(change) = self.single_nucleotide_change(ref_codon, alt_codon) {
                    results.push(change);
                }
            }
        }

        results
    }

    /// Backtranslate with known reference codon (context-aware).
    ///
    /// Returns all codons for the alternate amino acid, prioritizing
    /// those that differ by the fewest nucleotides.
    pub fn backtranslate_with_context(
        &self,
        ref_codon: &Codon,
        alt_aa: &AminoAcid,
    ) -> Vec<CodonChange> {
        let alt_codons = self.codon_table.codons_for(alt_aa);

        let mut results: Vec<CodonChange> = alt_codons
            .iter()
            .map(|alt| self.codon_difference(ref_codon, alt))
            .collect();

        // Sort by number of nucleotide changes (fewer is better)
        results.sort_by_key(|c| c.nucleotide_changes.len());

        results
    }

    /// Backtranslate a nonsense mutation (amino acid to stop codon).
    pub fn backtranslate_to_stop(&self, ref_aa: &AminoAcid) -> Vec<CodonChange> {
        let ref_codons = self.codon_table.codons_for(ref_aa);
        let stop_codons = self.codon_table.stop_codons();

        let mut results = Vec::new();

        for ref_codon in ref_codons {
            for stop_codon in stop_codons {
                if let Some(change) = self.single_nucleotide_change(ref_codon, stop_codon) {
                    results.push(change);
                }
            }
        }

        results
    }

    /// Backtranslate a stop loss (stop codon to amino acid).
    pub fn backtranslate_stop_loss(&self, alt_aa: &AminoAcid) -> Vec<CodonChange> {
        let stop_codons = self.codon_table.stop_codons();
        let alt_codons = self.codon_table.codons_for(alt_aa);

        let mut results = Vec::new();

        for stop_codon in stop_codons {
            for alt_codon in alt_codons {
                if let Some(change) = self.single_nucleotide_change(stop_codon, alt_codon) {
                    results.push(change);
                }
            }
        }

        results
    }

    /// Get the codon table.
    pub fn codon_table(&self) -> &CodonTable {
        &self.codon_table
    }

    /// Find single nucleotide change between two codons.
    fn single_nucleotide_change(
        &self,
        ref_codon: &Codon,
        alt_codon: &Codon,
    ) -> Option<CodonChange> {
        let diff = self.codon_difference(ref_codon, alt_codon);
        if diff.nucleotide_changes.len() == 1 {
            Some(diff)
        } else {
            None
        }
    }

    /// Calculate difference between two codons.
    fn codon_difference(&self, ref_codon: &Codon, alt_codon: &Codon) -> CodonChange {
        let ref_bases = ref_codon.bases();
        let alt_bases = alt_codon.bases();

        let mut changed_positions = Vec::new();
        let mut nucleotide_changes = Vec::new();

        for i in 0..3 {
            if ref_bases[i] != alt_bases[i] {
                changed_positions.push((i + 1) as u8);
                nucleotide_changes.push(((i + 1) as u8, ref_bases[i], alt_bases[i]));
            }
        }

        CodonChange {
            ref_codon: ref_codon.clone(),
            alt_codon: alt_codon.clone(),
            changed_positions,
            nucleotide_changes,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_backtranslate_leu_to_phe() {
        let bt = Backtranslator::standard();
        let changes = bt.backtranslate_substitution(&AminoAcid::Leu, &AminoAcid::Phe);

        // Leu has 6 codons, Phe has 2 codons
        // Not all combinations differ by 1 nucleotide
        assert!(!changes.is_empty());

        // All changes should be single nucleotide
        for change in &changes {
            assert_eq!(change.nucleotide_changes.len(), 1);
        }
    }

    #[test]
    fn test_backtranslate_with_context() {
        let bt = Backtranslator::standard();
        let ref_codon = Codon::parse("CTT").unwrap(); // Leu

        let changes = bt.backtranslate_with_context(&ref_codon, &AminoAcid::Phe);

        // Should have results for Phe codons (TTT, TTC)
        assert!(!changes.is_empty());

        // First result should be closest (fewest changes)
        assert!(
            changes[0].nucleotide_changes.len() <= changes.last().unwrap().nucleotide_changes.len()
        );
    }

    #[test]
    fn test_backtranslate_to_stop() {
        let bt = Backtranslator::standard();
        let changes = bt.backtranslate_to_stop(&AminoAcid::Gln);

        // Gln (CAA, CAG) can become stop (TAA, TAG, TGA) with single change
        // CAA -> TAA (C>T at pos 1)
        // CAG -> TAG (C>T at pos 1)
        assert!(changes.len() >= 2);
    }

    #[test]
    fn test_backtranslate_stop_loss() {
        let bt = Backtranslator::standard();
        let changes = bt.backtranslate_stop_loss(&AminoAcid::Gln);

        // Stop codons that can become Gln with single change
        // TAA -> CAA (T>C at pos 1)
        // TAG -> CAG (T>C at pos 1)
        assert!(changes.len() >= 2);
    }

    #[test]
    fn test_val_to_glu_braf() {
        // BRAF V600E: Val -> Glu
        let bt = Backtranslator::standard();
        let changes = bt.backtranslate_substitution(&AminoAcid::Val, &AminoAcid::Glu);

        // Val (GTT, GTC, GTA, GTG) -> Glu (GAA, GAG)
        // GTT -> GAT is Asp, not Glu - need 2 changes for GTT -> GAA
        // But GTG -> GAG (T>A at position 2) is single change
        assert!(!changes.is_empty());

        // Find GTG -> GAG
        let gtg_to_gag = changes
            .iter()
            .find(|c| c.ref_codon.to_string() == "GTG" && c.alt_codon.to_string() == "GAG");
        assert!(gtg_to_gag.is_some());
    }
}
