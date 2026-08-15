//! A verbatim slice of GRCh38 that the adjacency and copy-range suites edit.
//!
//! # Why a real slice rather than a synthetic string
//!
//! Every case in those suites is a real edit at real coordinates on GRCh38, so
//! a description this tree pins is one a caller could actually hand to ferro
//! against the assembly. Synthetic cores (`"AT".repeat(..)`) are cheaper, but
//! they make homopolymer and shift behaviour a property of the fixture rather
//! than of the genome — and shift interactions are exactly what these suites
//! are about.
//!
//! # Why it is still a `MockProvider` and not the prepared reference
//!
//! CI never sets `FERRO_MANIFEST`, so a reference-gated test skips green, and a
//! suite whose assertions live only behind that gate reports success while
//! checking nothing. So the bases travel *in this file*: the default lane needs
//! no reference and always runs.
//!
//! # Coordinates
//!
//! The slice is `NC_000001.11:g.1001300..1001899` (1-based, inclusive). The
//! provider serves it under [`LOCAL_CONTIG`] with the slice's first base at
//! local position 1, so local position `p` is assembly position
//! [`hg38_position`]`(p)`.
//!
//! Soft-masking is dropped: the assembly lower-cases repeats in part of this
//! window, HGVS carries no case semantics for DNA, and `MockProvider` compares
//! upper-case bases.

use ferro_hgvs::MockProvider;

/// The assembly contig the slice is taken from (GRCh38 chromosome 1).
pub const HG38_CONTIG: &str = "NC_000001.11";

/// 1-based inclusive assembly position of [`HG38_WINDOW`]'s first base.
pub const HG38_WINDOW_START: u64 = 1_001_300;

/// Contig name the local provider serves the slice under.
///
/// Deliberately *not* [`HG38_CONTIG`]: descriptions against this provider name
/// window-relative positions, and giving them the assembly accession would make
/// them look like assembly descriptions that are wrong by a megabase.
pub const LOCAL_CONTIG: &str = "NC_TESTWIN.1";

/// 600 bases of GRCh38, verbatim and upper-cased.
pub const HG38_WINDOW: &str = concat!(
    "CAAAAGCCTTCTTGACGCCGGGTGGTCCCAAAGGCTTCTGCGGGGTGGGGGGTCCTCAGG",
    "GGGGAAGCCTCAAGGGAGGGCGTGGCATTCCCAGGGTGCGAAGGGGGCGCAGGGACGAGG",
    "GAGGTGGGGAGGGGGAGCTGGGCCAGCGAGAACCGGGAGCTTCTGGTCGGGGAGGGAGTC",
    "GGGGAACTTTTTGGGGAGCTTTTCTGAGCCAGGGAGTCGGCTGATTGGCAGGTTCGCCCC",
    "TGCCCGGGCACCTGGACCCAGGGTTTCTGTGCGGAAGCTTCCCCTCCCCTCGGACCCCAC",
    "GTCTAATCTGGCCCCAAGCAAAGTCCTGCGGCCCACGCGGGAAGGCGCCCTCTTCGCGGC",
    "GCTGACCCCGGCCCTCCGCGGTGCCCCTGAGGCGCCCCCCACACCCCGCCGCTTGCACAG",
    "GGGCGCGGGGGGCTGCGAGGCCGGAGCGGGGGTGGCGCCCTCTCCGCCGAGAGGCTGTCC",
    "GCGCCCCTCGCCGACTGGGGAAAGCCGCGGGGGCTGGGCGGGCGTCTCGGAGGTGGCCCC",
    "GCGAGCACTTAAGCCCCGGCTCTCCTGCCCCGACCTCTCTGCGCGCGCCTCGGCGCTGGA",
);

/// The assembly position of local position `local` (1-based on both sides).
pub fn hg38_position(local: u64) -> u64 {
    assert!(
        local >= 1 && local <= HG38_WINDOW.len() as u64,
        "local position {local} is outside the window"
    );
    HG38_WINDOW_START + local - 1
}

/// A single-contig provider serving [`HG38_WINDOW`] as [`LOCAL_CONTIG`].
pub fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(LOCAL_CONTIG, HG38_WINDOW.to_string());
    provider
}

/// Render a local-coordinate `g.` description against [`LOCAL_CONTIG`].
pub fn local_desc(body: &str) -> String {
    format!("{LOCAL_CONTIG}:g.{body}")
}

/// The base at local position `local`.
pub fn base_at(local: u64) -> u8 {
    HG38_WINDOW.as_bytes()[(local - 1) as usize]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The slice is the length its coordinate arithmetic assumes.
    ///
    /// A truncated or re-wrapped paste would otherwise shift every case in both
    /// suites by a silent offset and still normalize to *something*.
    #[test]
    fn window_is_the_declared_length() {
        assert_eq!(HG38_WINDOW.len(), 600);
    }

    /// Spot-checks against the assembly, so a corrupted paste fails here rather
    /// than as diffuse wrongness in the suites that read it.
    #[test]
    fn window_bases_match_the_assembly() {
        let spots: &[(u64, u8)] = &[
            (1, b'C'),
            (100, b'G'),
            (301, b'G'),
            (450, b'G'),
            (600, b'A'),
        ];
        for &(local, expected) in spots {
            assert_eq!(
                base_at(local),
                expected,
                "local {} (assembly {}) drifted",
                local,
                hg38_position(local)
            );
        }
    }

    /// Only unambiguous bases: an `N` run would make every edit over it vacuous
    /// and every assertion in the suites trivially true.
    #[test]
    fn window_has_no_ambiguous_bases() {
        assert!(
            HG38_WINDOW
                .bytes()
                .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')),
            "window carries a non-ACGT base"
        );
    }

    /// The local/assembly mapping is the identity the suites document.
    #[test]
    fn local_and_assembly_positions_correspond() {
        assert_eq!(hg38_position(1), HG38_WINDOW_START);
        assert_eq!(
            hg38_position(HG38_WINDOW.len() as u64),
            HG38_WINDOW_START + HG38_WINDOW.len() as u64 - 1
        );
    }
}
