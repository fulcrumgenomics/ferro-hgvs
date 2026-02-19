//! DNA sequence utilities for CLI operations

/// Reverse complement a DNA sequence
///
/// Reverses the sequence and complements each nucleotide:
/// - A <-> T
/// - G <-> C
/// - Case is preserved
/// - Non-ATGC characters pass through unchanged
///
/// # Examples
///
/// ```
/// use ferro_hgvs::cli::reverse_complement;
///
/// assert_eq!(reverse_complement("ATGC"), "GCAT");
/// assert_eq!(reverse_complement("aattggcc"), "ggccaatt");
/// assert_eq!(reverse_complement("ATGN"), "NCAT");
/// ```
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'a' => 't',
            't' => 'a',
            'g' => 'c',
            'c' => 'g',
            _ => c,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement_basic() {
        assert_eq!(reverse_complement("ATGC"), "GCAT");
    }

    #[test]
    fn test_reverse_complement_empty() {
        assert_eq!(reverse_complement(""), "");
    }

    #[test]
    fn test_reverse_complement_single_base() {
        assert_eq!(reverse_complement("A"), "T");
        assert_eq!(reverse_complement("T"), "A");
        assert_eq!(reverse_complement("G"), "C");
        assert_eq!(reverse_complement("C"), "G");
    }

    #[test]
    fn test_reverse_complement_lowercase() {
        assert_eq!(reverse_complement("atgc"), "gcat");
        assert_eq!(reverse_complement("aattggcc"), "ggccaatt");
    }

    #[test]
    fn test_reverse_complement_mixed_case() {
        assert_eq!(reverse_complement("AtGc"), "gCaT");
        assert_eq!(reverse_complement("AaTtGgCc"), "gGcCaAtT");
    }

    #[test]
    fn test_reverse_complement_with_n() {
        assert_eq!(reverse_complement("ATGN"), "NCAT");
        assert_eq!(reverse_complement("NNNN"), "NNNN");
        assert_eq!(reverse_complement("ANT"), "ANT");
    }

    #[test]
    fn test_reverse_complement_preserves_unknown_chars() {
        assert_eq!(reverse_complement("ATGXYZ"), "ZYXCAT");
        assert_eq!(reverse_complement("A-T-G"), "C-A-T");
    }

    #[test]
    fn test_reverse_complement_palindrome() {
        // ATAT is its own reverse complement
        assert_eq!(reverse_complement("ATAT"), "ATAT");
        // GCGC is its own reverse complement
        assert_eq!(reverse_complement("GCGC"), "GCGC");
    }

    #[test]
    fn test_reverse_complement_long_sequence() {
        let seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC";
        let rc = reverse_complement(seq);
        // Applying reverse complement twice should give original
        assert_eq!(reverse_complement(&rc), seq);
    }

    #[test]
    fn test_reverse_complement_real_codon() {
        // ATG (start codon) -> CAT (reverse complement)
        assert_eq!(reverse_complement("ATG"), "CAT");
        // TAA (stop codon) -> TTA
        assert_eq!(reverse_complement("TAA"), "TTA");
        // TAG (stop codon) -> CTA
        assert_eq!(reverse_complement("TAG"), "CTA");
        // TGA (stop codon) -> TCA
        assert_eq!(reverse_complement("TGA"), "TCA");
    }
}
