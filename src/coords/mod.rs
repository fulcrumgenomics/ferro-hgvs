//! Type-safe coordinate system wrappers
//!
//! This module provides newtype wrappers that make the coordinate basis
//! (0-based vs 1-based) explicit at the type level. This prevents bugs
//! from mixing coordinate systems at compile time.
//!
//! # Coordinate Systems in Bioinformatics
//!
//! | Type | Basis | Use Cases |
//! |------|-------|-----------|
//! | [`ZeroBasedPos`] | 0-based | Array indexing, SPDI, cdot genomic, BED |
//! | [`OneBasedPos`] | 1-based | HGVS, VCF, cdot transcript, GFF/GTF, SAM |
//!
//! # Design Principles
//!
//! 1. **No implicit conversion**: Must call explicit methods to convert
//! 2. **Self-documenting**: Type signature shows coordinate system
//! 3. **Zero-cost abstraction**: Compiles to same code as raw integers
//! 4. **Validation**: 1-based positions reject zero at construction
//!
//! # Examples
//!
//! ```
//! use ferro_hgvs::coords::{ZeroBasedPos, OneBasedPos};
//!
//! // Converting between coordinate systems
//! let zb = ZeroBasedPos::new(99);  // index 99
//! let ob = zb.to_one_based();      // position 100
//! assert_eq!(ob.value(), 100);
//!
//! // Back to 0-based
//! let zb_again = ob.to_zero_based();
//! assert_eq!(zb_again.value(), 99);
//!
//! // Use as array index
//! let seq = b"ACGT";
//! let idx = ZeroBasedPos::new(2);
//! assert_eq!(seq[idx.as_index()], b'G');
//! ```

use serde::{Deserialize, Serialize};
use std::fmt;

/// A 0-based position (array-style indexing)
///
/// Used for:
/// - Rust slice/array indexing
/// - SPDI positions
/// - cdot genomic coordinates
/// - BED format
/// - Internal sequence manipulation
///
/// # Examples
///
/// ```
/// use ferro_hgvs::coords::ZeroBasedPos;
///
/// let pos = ZeroBasedPos::new(0);  // First element
/// assert_eq!(pos.as_index(), 0);
/// assert_eq!(pos.to_one_based().value(), 1);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct ZeroBasedPos(u64);

/// A 1-based position (human-readable indexing)
///
/// Used for:
/// - HGVS g., c., n., r., p. positions
/// - VCF POS field
/// - cdot transcript coordinates (tx_start, tx_end)
/// - SAM/BAM format
/// - GFF/GTF format
///
/// # Invariant
///
/// Position must be >= 1. Position 0 is invalid in 1-based systems.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::coords::OneBasedPos;
///
/// let pos = OneBasedPos::new(1);  // First element
/// assert_eq!(pos.value(), 1);
/// assert_eq!(pos.to_zero_based().value(), 0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct OneBasedPos(u64);

impl ZeroBasedPos {
    /// Create a new 0-based position
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::coords::ZeroBasedPos;
    ///
    /// let pos = ZeroBasedPos::new(0);  // Valid: first element
    /// let pos = ZeroBasedPos::new(99); // Valid: 100th element
    /// ```
    #[inline]
    pub const fn new(pos: u64) -> Self {
        Self(pos)
    }

    /// Get the raw value
    #[inline]
    pub const fn value(self) -> u64 {
        self.0
    }

    /// Convert to 1-based position
    ///
    /// This adds 1 to the position value.
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::coords::ZeroBasedPos;
    ///
    /// let zb = ZeroBasedPos::new(99);  // index 99
    /// let ob = zb.to_one_based();      // position 100
    /// assert_eq!(ob.value(), 100);
    /// ```
    #[inline]
    pub const fn to_one_based(self) -> OneBasedPos {
        OneBasedPos(self.0 + 1)
    }

    /// Use as array index (returns usize)
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::coords::ZeroBasedPos;
    ///
    /// let seq = b"ACGT";
    /// let pos = ZeroBasedPos::new(2);
    /// assert_eq!(seq[pos.as_index()], b'G');
    /// ```
    #[inline]
    pub const fn as_index(self) -> usize {
        self.0 as usize
    }

    /// Add an offset to this position
    ///
    /// Returns None if the result would overflow.
    #[inline]
    pub const fn checked_add(self, offset: u64) -> Option<Self> {
        match self.0.checked_add(offset) {
            Some(v) => Some(Self(v)),
            None => None,
        }
    }

    /// Subtract an offset from this position
    ///
    /// Returns None if the result would underflow.
    #[inline]
    pub const fn checked_sub(self, offset: u64) -> Option<Self> {
        match self.0.checked_sub(offset) {
            Some(v) => Some(Self(v)),
            None => None,
        }
    }

    /// Saturating subtraction
    #[inline]
    pub const fn saturating_sub(self, offset: u64) -> Self {
        Self(self.0.saturating_sub(offset))
    }
}

impl OneBasedPos {
    /// Create a new 1-based position
    ///
    /// # Panics
    ///
    /// Panics if pos is 0, which is invalid in 1-based coordinate systems.
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::coords::OneBasedPos;
    ///
    /// let pos = OneBasedPos::new(1);   // Valid: first element
    /// let pos = OneBasedPos::new(100); // Valid: 100th element
    /// ```
    ///
    /// ```should_panic
    /// use ferro_hgvs::coords::OneBasedPos;
    ///
    /// let pos = OneBasedPos::new(0);  // Panics!
    /// ```
    #[inline]
    pub fn new(pos: u64) -> Self {
        assert!(pos > 0, "1-based position cannot be 0");
        Self(pos)
    }

    /// Create without validation (for performance-critical paths)
    ///
    /// # Safety
    ///
    /// Caller must ensure pos > 0. Using this with pos == 0 will result
    /// in incorrect behavior when converting to 0-based positions.
    #[inline]
    pub const fn new_unchecked(pos: u64) -> Self {
        Self(pos)
    }

    /// Try to create a 1-based position, returning None if invalid
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::coords::OneBasedPos;
    ///
    /// assert!(OneBasedPos::try_new(1).is_some());
    /// assert!(OneBasedPos::try_new(0).is_none());
    /// ```
    #[inline]
    pub const fn try_new(pos: u64) -> Option<Self> {
        if pos > 0 {
            Some(Self(pos))
        } else {
            None
        }
    }

    /// Get the raw value
    #[inline]
    pub const fn value(self) -> u64 {
        self.0
    }

    /// Convert to 0-based position
    ///
    /// This subtracts 1 from the position value.
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::coords::OneBasedPos;
    ///
    /// let ob = OneBasedPos::new(100); // position 100
    /// let zb = ob.to_zero_based();    // index 99
    /// assert_eq!(zb.value(), 99);
    /// ```
    #[inline]
    pub const fn to_zero_based(self) -> ZeroBasedPos {
        ZeroBasedPos(self.0 - 1)
    }

    /// Add an offset to this position
    ///
    /// Returns None if the result would overflow.
    #[inline]
    pub const fn checked_add(self, offset: u64) -> Option<Self> {
        match self.0.checked_add(offset) {
            Some(v) => Some(Self(v)),
            None => None,
        }
    }

    /// Subtract an offset from this position
    ///
    /// Returns None if the result would underflow below 1.
    #[inline]
    pub const fn checked_sub(self, offset: u64) -> Option<Self> {
        match self.0.checked_sub(offset) {
            Some(v) if v > 0 => Some(Self(v)),
            _ => None,
        }
    }
}

// Display implementations show the basis clearly in debug output
impl fmt::Display for ZeroBasedPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl fmt::Display for OneBasedPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

// Conversion from raw integers (explicit)
impl From<ZeroBasedPos> for u64 {
    fn from(pos: ZeroBasedPos) -> Self {
        pos.0
    }
}

impl From<OneBasedPos> for u64 {
    fn from(pos: OneBasedPos) -> Self {
        pos.0
    }
}

/// Generic interval type
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Interval<P> {
    pub start: P,
    pub end: P,
}

/// A 0-based half-open interval [start, end)
///
/// Used for: BED format, cdot genomic coordinates, Python slicing conventions.
/// The start is inclusive and the end is exclusive.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::coords::{ZeroBasedPos, ZeroBasedInterval};
///
/// // Interval covering positions 0, 1, 2 (length 3)
/// let interval = ZeroBasedInterval {
///     start: ZeroBasedPos::new(0),
///     end: ZeroBasedPos::new(3),
/// };
/// assert_eq!(interval.len(), 3);
/// ```
pub type ZeroBasedInterval = Interval<ZeroBasedPos>;

/// A 1-based closed interval [start, end]
///
/// Used for: HGVS, VCF, GFF/GTF, SAM.
/// Both start and end are inclusive.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::coords::{OneBasedPos, OneBasedInterval};
///
/// // Interval covering positions 1, 2, 3 (length 3)
/// let interval = OneBasedInterval {
///     start: OneBasedPos::new(1),
///     end: OneBasedPos::new(3),
/// };
/// assert_eq!(interval.len(), 3);
/// ```
pub type OneBasedInterval = Interval<OneBasedPos>;

impl ZeroBasedInterval {
    /// Create a new 0-based half-open interval
    pub const fn new(start: ZeroBasedPos, end: ZeroBasedPos) -> Self {
        Self { start, end }
    }

    /// Convert to 1-based closed interval
    ///
    /// # Panics
    ///
    /// Panics if the interval is empty (end <= start).
    pub fn to_one_based_closed(self) -> OneBasedInterval {
        assert!(
            self.end.value() > self.start.value(),
            "Cannot convert empty interval to 1-based closed"
        );
        Interval {
            start: self.start.to_one_based(),
            // For half-open to closed: end - 1, then convert to 1-based
            end: ZeroBasedPos::new(self.end.value() - 1).to_one_based(),
        }
    }

    /// Length of the interval
    #[inline]
    pub const fn len(&self) -> u64 {
        self.end.value() - self.start.value()
    }

    /// Check if interval is empty
    #[inline]
    pub const fn is_empty(&self) -> bool {
        self.end.value() <= self.start.value()
    }
}

impl OneBasedInterval {
    /// Create a new 1-based closed interval
    pub fn new(start: OneBasedPos, end: OneBasedPos) -> Self {
        Self { start, end }
    }

    /// Convert to 0-based half-open interval
    pub fn to_zero_based_half_open(self) -> ZeroBasedInterval {
        Interval {
            start: self.start.to_zero_based(),
            // For closed to half-open: convert end to 0-based, then add 1
            end: ZeroBasedPos::new(self.end.value()), // end.value() is already +1 from 0-based perspective
        }
    }

    /// Length of the interval (inclusive)
    #[inline]
    pub const fn len(&self) -> u64 {
        self.end.value() - self.start.value() + 1
    }

    /// Check if interval is empty (which shouldn't happen with valid 1-based positions)
    #[inline]
    pub const fn is_empty(&self) -> bool {
        self.end.value() < self.start.value()
    }
}

// ============================================================================
// Convenience conversion functions
// ============================================================================

/// Convert a 1-based HGVS position to a 0-based array index
///
/// This is the most common conversion in bioinformatics code.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::coords::hgvs_pos_to_index;
///
/// assert_eq!(hgvs_pos_to_index(1), 0);   // g.1 -> index 0
/// assert_eq!(hgvs_pos_to_index(100), 99); // g.100 -> index 99
/// ```
#[inline]
pub const fn hgvs_pos_to_index(pos: u64) -> usize {
    (pos - 1) as usize
}

/// Convert a 0-based array index to a 1-based HGVS position
///
/// # Examples
///
/// ```
/// use ferro_hgvs::coords::index_to_hgvs_pos;
///
/// assert_eq!(index_to_hgvs_pos(0), 1);   // index 0 -> g.1
/// assert_eq!(index_to_hgvs_pos(99), 100); // index 99 -> g.100
/// ```
#[inline]
pub const fn index_to_hgvs_pos(idx: usize) -> u64 {
    idx as u64 + 1
}

/// Convert cdot genomic coordinates (0-based half-open) to 1-based closed
///
/// cdot uses 0-based half-open for genomic coordinates.
///
/// # Arguments
///
/// * `start` - 0-based inclusive start
/// * `end` - 0-based exclusive end
///
/// # Returns
///
/// (1-based inclusive start, 1-based inclusive end)
///
/// # Examples
///
/// ```
/// use ferro_hgvs::coords::cdot_genomic_to_closed;
///
/// // cdot: [0, 10) -> 1-based: [1, 10]
/// let (start, end) = cdot_genomic_to_closed(0, 10);
/// assert_eq!(start, 1);
/// assert_eq!(end, 10);
/// ```
#[inline]
pub const fn cdot_genomic_to_closed(start: u64, end: u64) -> (u64, u64) {
    // 0-based start + 1 = 1-based start
    // 0-based exclusive end = 1-based inclusive end (no change needed)
    (start + 1, end)
}

/// Document that cdot transcript coordinates are already 1-based
///
/// This function exists as documentation - cdot tx_start/tx_end are 1-based,
/// NOT 0-based like cdot genomic coordinates. This was the source of bug 944a4e9.
///
/// # Arguments
///
/// * `tx_start` - 1-based transcript start (NOT 0-based!)
/// * `tx_end` - 1-based transcript end (NOT 0-based!)
///
/// # Returns
///
/// The same values, unchanged (this is a documentation function).
#[inline]
pub const fn cdot_tx_coords(tx_start: u64, tx_end: u64) -> (u64, u64) {
    // No conversion needed - already 1-based
    // This function exists to make it explicit that no conversion is needed
    (tx_start, tx_end)
}

/// Convert SPDI 0-based position to HGVS 1-based
///
/// SPDI uses 0-based interbase coordinates.
#[inline]
pub const fn spdi_to_hgvs_pos(spdi_pos: u64) -> u64 {
    spdi_pos + 1
}

/// Convert HGVS 1-based position to SPDI 0-based
///
/// SPDI uses 0-based interbase coordinates.
#[inline]
pub const fn hgvs_to_spdi_pos(hgvs_pos: u64) -> u64 {
    hgvs_pos - 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_based_creation() {
        let pos = ZeroBasedPos::new(0);
        assert_eq!(pos.value(), 0);
        assert_eq!(pos.as_index(), 0);

        let pos = ZeroBasedPos::new(99);
        assert_eq!(pos.value(), 99);
        assert_eq!(pos.as_index(), 99);
    }

    #[test]
    fn test_one_based_creation() {
        let pos = OneBasedPos::new(1);
        assert_eq!(pos.value(), 1);

        let pos = OneBasedPos::new(100);
        assert_eq!(pos.value(), 100);
    }

    #[test]
    #[should_panic(expected = "1-based position cannot be 0")]
    fn test_one_based_rejects_zero() {
        let _ = OneBasedPos::new(0);
    }

    #[test]
    fn test_try_new_one_based() {
        assert!(OneBasedPos::try_new(0).is_none());
        assert!(OneBasedPos::try_new(1).is_some());
    }

    #[test]
    fn test_zero_to_one_based_conversion() {
        // Index 0 (first element) -> position 1
        assert_eq!(ZeroBasedPos::new(0).to_one_based().value(), 1);

        // Index 99 -> position 100
        assert_eq!(ZeroBasedPos::new(99).to_one_based().value(), 100);
    }

    #[test]
    fn test_one_to_zero_based_conversion() {
        // Position 1 -> index 0
        assert_eq!(OneBasedPos::new(1).to_zero_based().value(), 0);

        // Position 100 -> index 99
        assert_eq!(OneBasedPos::new(100).to_zero_based().value(), 99);
    }

    #[test]
    fn test_roundtrip_conversion() {
        // 0-based -> 1-based -> 0-based
        for i in 0..100 {
            let zb = ZeroBasedPos::new(i);
            let ob = zb.to_one_based();
            let zb_again = ob.to_zero_based();
            assert_eq!(zb, zb_again);
        }

        // 1-based -> 0-based -> 1-based
        for i in 1..=100 {
            let ob = OneBasedPos::new(i);
            let zb = ob.to_zero_based();
            let ob_again = zb.to_one_based();
            assert_eq!(ob, ob_again);
        }
    }

    #[test]
    fn test_zero_based_interval_length() {
        // [0, 10) has length 10
        let interval = ZeroBasedInterval::new(ZeroBasedPos::new(0), ZeroBasedPos::new(10));
        assert_eq!(interval.len(), 10);

        // [5, 8) has length 3
        let interval = ZeroBasedInterval::new(ZeroBasedPos::new(5), ZeroBasedPos::new(8));
        assert_eq!(interval.len(), 3);
    }

    #[test]
    fn test_one_based_interval_length() {
        // [1, 10] has length 10
        let interval = OneBasedInterval::new(OneBasedPos::new(1), OneBasedPos::new(10));
        assert_eq!(interval.len(), 10);

        // [5, 7] has length 3
        let interval = OneBasedInterval::new(OneBasedPos::new(5), OneBasedPos::new(7));
        assert_eq!(interval.len(), 3);
    }

    #[test]
    fn test_interval_conversion() {
        // [0, 10) in 0-based half-open = [1, 10] in 1-based closed
        let ho = ZeroBasedInterval::new(ZeroBasedPos::new(0), ZeroBasedPos::new(10));
        let closed = ho.to_one_based_closed();

        assert_eq!(closed.start.value(), 1);
        assert_eq!(closed.end.value(), 10);
        assert_eq!(ho.len(), closed.len()); // Same length

        // [1, 10] in 1-based closed = [0, 10) in 0-based half-open
        let closed = OneBasedInterval::new(OneBasedPos::new(1), OneBasedPos::new(10));
        let ho = closed.to_zero_based_half_open();

        assert_eq!(ho.start.value(), 0);
        assert_eq!(ho.end.value(), 10);
    }

    #[test]
    fn test_hgvs_pos_to_index() {
        assert_eq!(hgvs_pos_to_index(1), 0);
        assert_eq!(hgvs_pos_to_index(100), 99);
    }

    #[test]
    fn test_index_to_hgvs_pos() {
        assert_eq!(index_to_hgvs_pos(0), 1);
        assert_eq!(index_to_hgvs_pos(99), 100);
    }

    #[test]
    fn test_cdot_genomic_to_closed() {
        // [0, 10) -> [1, 10]
        let (start, end) = cdot_genomic_to_closed(0, 10);
        assert_eq!(start, 1);
        assert_eq!(end, 10);

        // [99, 200) -> [100, 200]
        let (start, end) = cdot_genomic_to_closed(99, 200);
        assert_eq!(start, 100);
        assert_eq!(end, 200);
    }

    #[test]
    fn test_cdot_tx_coords_unchanged() {
        // cdot tx coords are already 1-based
        let (start, end) = cdot_tx_coords(1, 100);
        assert_eq!(start, 1);
        assert_eq!(end, 100);
    }

    #[test]
    fn test_spdi_hgvs_conversion() {
        // SPDI 0 = HGVS 1
        assert_eq!(spdi_to_hgvs_pos(0), 1);
        assert_eq!(hgvs_to_spdi_pos(1), 0);

        // SPDI 12344 = HGVS 12345
        assert_eq!(spdi_to_hgvs_pos(12344), 12345);
        assert_eq!(hgvs_to_spdi_pos(12345), 12344);
    }

    #[test]
    fn test_checked_arithmetic() {
        let zb = ZeroBasedPos::new(5);
        assert_eq!(zb.checked_add(3), Some(ZeroBasedPos::new(8)));
        assert_eq!(zb.checked_sub(3), Some(ZeroBasedPos::new(2)));
        assert_eq!(zb.checked_sub(10), None); // Would underflow

        let ob = OneBasedPos::new(5);
        assert_eq!(ob.checked_add(3), Some(OneBasedPos::new(8)));
        assert_eq!(ob.checked_sub(3), Some(OneBasedPos::new(2)));
        assert_eq!(ob.checked_sub(5), None); // Would be 0, invalid for 1-based
    }

    #[test]
    fn test_ordering() {
        let a = ZeroBasedPos::new(5);
        let b = ZeroBasedPos::new(10);
        assert!(a < b);

        let a = OneBasedPos::new(5);
        let b = OneBasedPos::new(10);
        assert!(a < b);
    }
}
