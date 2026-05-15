//! DNA sequence utilities shared across the library.

/// Reverse complement a DNA sequence.
///
/// Reverses the sequence and complements each nucleotide, including IUPAC
/// ambiguity codes. Case is preserved. Characters that are not a recognized
/// DNA / IUPAC base pass through unchanged.
///
/// IUPAC complements:
/// - A↔T, G↔C, U→A
/// - R↔Y, S↔S, W↔W, K↔M
/// - B↔V, D↔H, N↔N
///
/// # Examples
///
/// ```
/// use ferro_hgvs::sequence::reverse_complement;
///
/// assert_eq!(reverse_complement("ATGC"), "GCAT");
/// assert_eq!(reverse_complement("aattggcc"), "ggccaatt");
/// assert_eq!(reverse_complement("ATGN"), "NCAT");
/// assert_eq!(reverse_complement("RYK"), "MRY");
/// ```
pub fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(complement_char).collect()
}

fn complement_char(c: char) -> char {
    match c {
        'A' => 'T',
        'T' | 'U' => 'A',
        'G' => 'C',
        'C' => 'G',
        'R' => 'Y',
        'Y' => 'R',
        'S' => 'S',
        'W' => 'W',
        'K' => 'M',
        'M' => 'K',
        'B' => 'V',
        'V' => 'B',
        'D' => 'H',
        'H' => 'D',
        'N' => 'N',
        'a' => 't',
        't' | 'u' => 'a',
        'g' => 'c',
        'c' => 'g',
        'r' => 'y',
        'y' => 'r',
        's' => 's',
        'w' => 'w',
        'k' => 'm',
        'm' => 'k',
        'b' => 'v',
        'v' => 'b',
        'd' => 'h',
        'h' => 'd',
        'n' => 'n',
        other => other,
    }
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
        // X, Z, Q are not IUPAC bases — they pass through unchanged.
        assert_eq!(reverse_complement("ATGXZQ"), "QZXCAT");
        assert_eq!(reverse_complement("A-T-G"), "C-A-T");
    }

    #[test]
    fn test_reverse_complement_iupac_ambiguity_codes() {
        // R (puRine: A/G) ↔ Y (pYrimidine: C/T)
        assert_eq!(reverse_complement("R"), "Y");
        assert_eq!(reverse_complement("Y"), "R");
        // S (G/C, self-complementary), W (A/T, self-complementary)
        assert_eq!(reverse_complement("S"), "S");
        assert_eq!(reverse_complement("W"), "W");
        // K (G/T) ↔ M (A/C)
        assert_eq!(reverse_complement("K"), "M");
        assert_eq!(reverse_complement("M"), "K");
        // B (not A) ↔ V (not T); D (not C) ↔ H (not G)
        assert_eq!(reverse_complement("B"), "V");
        assert_eq!(reverse_complement("V"), "B");
        assert_eq!(reverse_complement("D"), "H");
        assert_eq!(reverse_complement("H"), "D");
        // Combined: RYNK reversed = KNYR, complemented = MNRY
        assert_eq!(reverse_complement("RYNK"), "MNRY");
        // Case preserved for IUPAC codes
        // "ryk" reversed → "kyr"; k→m, y→r, r→y → "mry"
        assert_eq!(reverse_complement("ryk"), "mry");
        // U (RNA) complements to A, just like T. Output uses the DNA letter T
        // is mapped from A, so revcomp of RNA "AUGC" reads back as DNA "GCAT".
        assert_eq!(reverse_complement("AUGC"), "GCAT");
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
