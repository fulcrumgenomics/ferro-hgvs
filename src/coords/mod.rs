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
//! | [`ZeroBasedPos`] | 0-based | Array indexing, SPDI, **raw** cdot genomic, BED |
//! | [`OneBasedPos`] | 1-based | HGVS, VCF, **raw** cdot transcript, GFF/GTF, SAM |
//!
//! # Which cdot layer does this module describe?
//!
//! Every `cdot_*` helper in this module describes the **raw cdot JSON** as it
//! arrives on disk — genomic bounds 0-based half-open, transcript bounds
//! 1-based inclusive. That is *not* the layout ferro holds in memory: on load,
//! `RawCdotTranscript::from_genome_build` (private, in [`crate::data::cdot`])
//! converts both axes (#742), so
//! [`crate::data::cdot::CdotTranscript`] carries genome 1-based `[incl, excl)`
//! and transcript **0-based** half-open. The table at the top of
//! [`crate::data::cdot`] describes that *internal* layout, using the same field
//! names (`tx_start`/`tx_end`) as the raw format.
//!
//! Two tables, two layers, same names: always check which one a statement is
//! about before acting on it.
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
/// # Panics
///
/// Panics if `pos` is `0`, which is not a valid 1-based position. Callers that
/// cannot guarantee a validated position must use
/// [`hgvs_pos_to_index_checked`] instead.
///
/// The check is deliberate and unconditional (#1282). Before it, `pos == 0`
/// underflowed: a debug build panicked from inside the subtraction with no
/// indication of which caller was at fault, and a release build — where
/// `[profile.release]` sets no `overflow-checks` — **wrapped**, silently
/// yielding an index near `usize::MAX`. Failing loudly in both profiles is the
/// lesser of the two, and mirrors [`OneBasedPos::new`], which has always
/// asserted the same invariant.
///
/// The function stays `const`, so the check is stronger than a runtime panic
/// wherever the position is a constant: `const _: usize = hgvs_pos_to_index(0);`
/// is a *compile* error, not a crash at run time.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::coords::hgvs_pos_to_index;
///
/// assert_eq!(hgvs_pos_to_index(1), 0);   // g.1 -> index 0
/// assert_eq!(hgvs_pos_to_index(100), 99); // g.100 -> index 99
/// ```
///
/// ```should_panic
/// use ferro_hgvs::coords::hgvs_pos_to_index;
///
/// let _ = hgvs_pos_to_index(0);  // Panics — there is no HGVS position 0
/// ```
#[inline]
#[track_caller]
pub const fn hgvs_pos_to_index(pos: u64) -> usize {
    assert!(
        pos > 0,
        "1-based HGVS position cannot be 0 (no position 0 exists on any HGVS axis); \
         use hgvs_pos_to_index_checked if the position is not known to be valid"
    );
    (pos - 1) as usize
}

/// Convert a 1-based HGVS position to a 0-based array index, or `None` if the
/// position is not a valid 1-based one.
///
/// The total counterpart of [`hgvs_pos_to_index`], for callers holding a
/// position they have not validated — a coordinate arrived at by arithmetic
/// (a 5'-shift, a clamp, an axis conversion) rather than read from a parsed
/// description. Prefer this at any boundary where "there is no such position"
/// is a real outcome that should become an `Err`, not a panic.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::coords::hgvs_pos_to_index_checked;
///
/// assert_eq!(hgvs_pos_to_index_checked(1), Some(0));
/// assert_eq!(hgvs_pos_to_index_checked(0), None);
/// ```
#[inline]
pub const fn hgvs_pos_to_index_checked(pos: u64) -> Option<usize> {
    if pos > 0 {
        Some((pos - 1) as usize)
    } else {
        None
    }
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

/// Convert **raw cdot JSON** genomic coordinates (0-based half-open) to 1-based closed
///
/// # Layer
///
/// This describes the **raw cdot JSON on disk**, where genomic bounds are
/// 0-based half-open. It does **not** describe ferro's in-memory
/// [`crate::data::cdot::CdotTranscript`], whose `exons[i][0..2]` were already
/// converted on load to genome 1-based `[incl, excl)` (#742) — feeding those
/// through this helper would convert twice.
///
/// Note the *target* here is 1-based **inclusive** on both ends, which differs
/// from the ingestion path's target of 1-based `[incl, excl)`.
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

/// Document that **raw cdot JSON** transcript coordinates are already 1-based
///
/// # Layer
///
/// This is a statement about the **raw cdot JSON on disk**: there, `tx_start` /
/// `tx_end` are 1-based **inclusive**, unlike the same record's genomic bounds,
/// which are 0-based half-open. Mixing the two up is the mistake this
/// no-op helper exists to make visible at the call site.
///
/// It is **not** a statement about ferro's in-memory
/// [`crate::data::cdot::CdotTranscript`]. `RawCdotTranscript::from_genome_build`
/// (private, in [`crate::data::cdot`])
/// converts the raw 1-based-inclusive tx bounds to **0-based half-open** on
/// load (#742), so `CdotTranscript.exons[i][2..4]` are 0-based and must not be
/// passed here. The basis table at the top of [`crate::data::cdot`] describes
/// that internal layout — with the same field names — so read the layer, not
/// just the name.
///
/// # Arguments
///
/// * `tx_start` - raw-cdot 1-based inclusive transcript start (NOT 0-based!)
/// * `tx_end` - raw-cdot 1-based inclusive transcript end (NOT 0-based!)
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
///
/// # Panics
///
/// Panics if `hgvs_pos` is `0`. This is the same defect [`hgvs_pos_to_index`]
/// carries and was found by the sweep #1282 asks for — the bare `hgvs_pos - 1`
/// wrapped to `u64::MAX` in release, producing an SPDI position astronomically
/// past the end of any sequence rather than an error. Use
/// [`hgvs_to_spdi_pos_checked`] where the position is not known to be valid.
///
/// As with [`hgvs_pos_to_index`], the function stays `const`, so a constant `0`
/// is rejected at compile time rather than at run time.
///
/// ```should_panic
/// use ferro_hgvs::coords::hgvs_to_spdi_pos;
///
/// let _ = hgvs_to_spdi_pos(0);  // Panics — there is no HGVS position 0
/// ```
#[inline]
#[track_caller]
pub const fn hgvs_to_spdi_pos(hgvs_pos: u64) -> u64 {
    assert!(
        hgvs_pos > 0,
        "1-based HGVS position cannot be 0 (no position 0 exists on any HGVS axis); \
         use hgvs_to_spdi_pos_checked if the position is not known to be valid"
    );
    hgvs_pos - 1
}

/// Convert HGVS 1-based position to SPDI 0-based, or `None` if the position is
/// not a valid 1-based one.
///
/// The total counterpart of [`hgvs_to_spdi_pos`].
///
/// # Examples
///
/// ```
/// use ferro_hgvs::coords::hgvs_to_spdi_pos_checked;
///
/// assert_eq!(hgvs_to_spdi_pos_checked(1), Some(0));
/// assert_eq!(hgvs_to_spdi_pos_checked(0), None);
/// ```
#[inline]
pub const fn hgvs_to_spdi_pos_checked(hgvs_pos: u64) -> Option<u64> {
    if hgvs_pos > 0 {
        Some(hgvs_pos - 1)
    } else {
        None
    }
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

    /// #1282: the two bare `pos - 1` conversions underflowed on `0`. A debug
    /// build panicked from inside the subtraction, naming neither the function
    /// nor the caller; a release build wrapped, because `[profile.release]`
    /// sets no `overflow-checks`, and handed back an index near `usize::MAX`.
    ///
    /// Both now assert, so the release-build wraparound — the more dangerous
    /// half, since it corrupts silently — cannot happen either.
    #[test]
    #[should_panic(expected = "1-based HGVS position cannot be 0")]
    fn hgvs_pos_to_index_rejects_zero() {
        let _ = hgvs_pos_to_index(0);
    }

    #[test]
    #[should_panic(expected = "1-based HGVS position cannot be 0")]
    fn hgvs_to_spdi_pos_rejects_zero() {
        let _ = hgvs_to_spdi_pos(0);
    }

    /// The checked forms are total: they answer for `0` instead of aborting,
    /// which is what a caller holding an arithmetic result needs.
    #[test]
    fn the_checked_conversions_are_total() {
        assert_eq!(hgvs_pos_to_index_checked(0), None);
        assert_eq!(hgvs_pos_to_index_checked(1), Some(0));
        assert_eq!(hgvs_pos_to_index_checked(100), Some(99));

        assert_eq!(hgvs_to_spdi_pos_checked(0), None);
        assert_eq!(hgvs_to_spdi_pos_checked(1), Some(0));
        assert_eq!(hgvs_to_spdi_pos_checked(100), Some(99));
    }

    /// Both conversions must stay `const`. The `assert!` does not cost that,
    /// and keeping it buys a strictly stronger check: in a const context the
    /// rejection of `0` is a *compile* error rather than a run-time panic.
    /// Dropping `const` would also be a semver-major change to a public API.
    ///
    /// This is a compile-time guard, not a runtime one — it fails the build if
    /// either function loses `const`, which no `#[test]` could catch.
    const _: () = {
        assert!(hgvs_pos_to_index(1) == 0);
        assert!(hgvs_pos_to_index(100) == 99);
        assert!(hgvs_to_spdi_pos(1) == 0);
        assert!(hgvs_to_spdi_pos(100) == 99);
    };

    /// Checked and unchecked must agree everywhere the unchecked one is defined
    /// — otherwise migrating a call site would change behaviour.
    #[test]
    fn checked_and_unchecked_agree_above_zero() {
        for pos in 1..=512u64 {
            assert_eq!(hgvs_pos_to_index_checked(pos), Some(hgvs_pos_to_index(pos)));
            assert_eq!(hgvs_to_spdi_pos_checked(pos), Some(hgvs_to_spdi_pos(pos)));
        }
    }
}
