//! Genomic (g.) coordinate handling
//!
//! # Coordinate Systems
//!
//! This module handles conversions between coordinate systems:
//!
//! | System | Basis | Notes |
//! |--------|-------|-------|
//! | HGVS genomic (g.) | 1-based | GenomePos.base is 1-based |
//! | Array indexing | 0-based | For accessing sequence bytes |
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::coords::{OneBasedPos, ZeroBasedPos};
use crate::error::FerroError;
use crate::hgvs::location::GenomePos;

/// Validate a genomic position
pub fn validate_genome_pos(pos: &GenomePos, seq_length: u64) -> Result<(), FerroError> {
    if pos.base < 1 {
        return Err(FerroError::InvalidCoordinates {
            msg: "Genomic position must be >= 1".to_string(),
        });
    }

    if pos.base > seq_length {
        return Err(FerroError::InvalidCoordinates {
            msg: format!(
                "Genomic position {} exceeds sequence length {}",
                pos.base, seq_length
            ),
        });
    }

    Ok(())
}

/// Convert a genomic position to 0-based coordinate
pub fn genome_to_zero_based(pos: &GenomePos) -> u64 {
    pos.base - 1
}

/// Convert a 0-based coordinate to genomic position
pub fn zero_based_to_genome(pos: u64) -> GenomePos {
    GenomePos::new(pos + 1)
}

// ============================================================================
// Type-safe conversion functions using newtypes
// ============================================================================

/// Convert a GenomePos to a type-safe 1-based position
///
/// This extracts the base position and wraps it in a OneBasedPos
/// for compile-time safety.
#[inline]
pub fn genome_pos_to_one_based(pos: &GenomePos) -> OneBasedPos {
    // GenomePos.base is documented as 1-based
    OneBasedPos::new(pos.base)
}

/// Convert a type-safe 1-based position to GenomePos
#[inline]
pub fn one_based_to_genome_pos(pos: OneBasedPos) -> GenomePos {
    GenomePos::new(pos.value())
}

/// Convert a GenomePos to a type-safe 0-based position
///
/// This is useful for array indexing operations.
#[inline]
pub fn genome_pos_to_zero_based(pos: &GenomePos) -> ZeroBasedPos {
    ZeroBasedPos::new(pos.base - 1)
}

/// Convert a type-safe 0-based position to GenomePos
#[inline]
pub fn zero_based_to_genome_pos(pos: ZeroBasedPos) -> GenomePos {
    GenomePos::new(pos.value() + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validate_genome_pos() {
        assert!(validate_genome_pos(&GenomePos::new(1), 1000).is_ok());
        assert!(validate_genome_pos(&GenomePos::new(1000), 1000).is_ok());
        assert!(validate_genome_pos(&GenomePos::new(1001), 1000).is_err());
    }

    #[test]
    fn test_coordinate_conversion() {
        let pos = GenomePos::new(100);
        assert_eq!(genome_to_zero_based(&pos), 99);

        let back = zero_based_to_genome(99);
        assert_eq!(back.base, 100);
    }

    #[test]
    fn test_type_safe_one_based_conversion() {
        let pos = GenomePos::new(100);
        let ob = genome_pos_to_one_based(&pos);
        assert_eq!(ob.value(), 100);

        let back = one_based_to_genome_pos(ob);
        assert_eq!(back.base, 100);
    }

    #[test]
    fn test_type_safe_zero_based_conversion() {
        let pos = GenomePos::new(100);
        let zb = genome_pos_to_zero_based(&pos);
        assert_eq!(zb.value(), 99);
        assert_eq!(zb.as_index(), 99);

        let back = zero_based_to_genome_pos(zb);
        assert_eq!(back.base, 100);
    }

    #[test]
    fn test_type_safe_roundtrip() {
        // GenomePos -> OneBasedPos -> ZeroBasedPos -> GenomePos
        let original = GenomePos::new(42);
        let ob = genome_pos_to_one_based(&original);
        let zb = ob.to_zero_based();
        let back = zero_based_to_genome_pos(zb);
        assert_eq!(back.base, original.base);
    }
}
