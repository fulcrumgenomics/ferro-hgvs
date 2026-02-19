//! Interval types for different coordinate systems
//!
//! Intervals represent ranges of positions for multi-position variants
//! like deletions, duplications, and insertions.

use super::location::{CdsPos, GenomePos, ProtPos, RnaPos, TxPos};
use super::uncertainty::Mu;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Represents a boundary position that may be a single position or a range
///
/// In HGVS notation, interval boundaries can be:
/// - A single position: `100` (certain) or `(100)` (uncertain) or `?` (unknown)
/// - A range: `(100_200)` meaning the boundary is somewhere between 100 and 200
///
/// This is used for complex uncertain intervals like:
/// - `c.(4185+1_4186-1)_(4357+1_4358-1)del` - exon deletion with uncertain intronic breakpoints
/// - `g.(?_43044294)_(43125364_?)del` - large deletion with uncertain outer boundaries
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum UncertainBoundary<T> {
    /// A single position (may be certain, uncertain, or unknown)
    Single(Mu<T>),
    /// A range indicating the boundary is somewhere within this range
    /// The parentheses wrap the entire range: (start_end)
    Range {
        /// Start of the range (may be unknown)
        start: Mu<T>,
        /// End of the range (may be unknown)
        end: Mu<T>,
    },
}

impl<T> UncertainBoundary<T> {
    /// Create a certain single boundary
    pub fn certain(pos: T) -> Self {
        Self::Single(Mu::Certain(pos))
    }

    /// Create an uncertain single boundary (position in parentheses)
    pub fn uncertain(pos: T) -> Self {
        Self::Single(Mu::Uncertain(pos))
    }

    /// Create an unknown boundary (?)
    pub fn unknown() -> Self {
        Self::Single(Mu::Unknown)
    }

    /// Create a range boundary (start_end)
    pub fn range(start: Mu<T>, end: Mu<T>) -> Self {
        Self::Range { start, end }
    }

    /// Check if this is a single position boundary
    pub fn is_single(&self) -> bool {
        matches!(self, Self::Single(_))
    }

    /// Check if this is a range boundary
    pub fn is_range(&self) -> bool {
        matches!(self, Self::Range { .. })
    }

    /// Get the inner Mu<T> if this is a single position
    pub fn as_single(&self) -> Option<&Mu<T>> {
        match self {
            Self::Single(mu) => Some(mu),
            Self::Range { .. } => None,
        }
    }

    /// Get the inner position value if this is a single certain or uncertain position
    /// Returns None for unknown positions or range boundaries
    pub fn inner(&self) -> Option<&T> {
        match self {
            Self::Single(mu) => mu.inner(),
            Self::Range { .. } => None,
        }
    }
}

impl<T: Clone> UncertainBoundary<T> {
    /// Convert from a simple Mu<T> to UncertainBoundary
    pub fn from_mu(mu: Mu<T>) -> Self {
        Self::Single(mu)
    }
}

impl<T: fmt::Display> fmt::Display for UncertainBoundary<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Single(mu) => write!(f, "{}", mu),
            Self::Range { start, end } => {
                // Range boundaries are always wrapped in parentheses
                write!(f, "({}_{}", start, end)?;
                write!(f, ")")
            }
        }
    }
}

/// An interval with complex uncertain boundaries
///
/// This extends the basic Interval type to support complex HGVS patterns
/// where boundaries can be ranges rather than single positions.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ComplexInterval<T> {
    pub start: UncertainBoundary<T>,
    pub end: UncertainBoundary<T>,
}

impl<T> ComplexInterval<T> {
    pub fn new(start: UncertainBoundary<T>, end: UncertainBoundary<T>) -> Self {
        Self { start, end }
    }
}

impl<T: fmt::Display> fmt::Display for ComplexInterval<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}_{}", self.start, self.end)
    }
}

/// Type aliases for complex intervals
pub type ComplexCdsInterval = ComplexInterval<CdsPos>;
pub type ComplexGenomeInterval = ComplexInterval<GenomePos>;

/// Generic interval with start and end positions
///
/// Supports both simple intervals (single positions) and complex intervals
/// (range boundaries like `(100+1_101-1)_(200+1_201-1)`).
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Interval<T> {
    pub start: UncertainBoundary<T>,
    pub end: UncertainBoundary<T>,
}

impl<T> Interval<T> {
    pub fn new(start: T, end: T) -> Self {
        Self {
            start: UncertainBoundary::certain(start),
            end: UncertainBoundary::certain(end),
        }
    }

    pub fn point(pos: T) -> Self
    where
        T: Clone,
    {
        Self {
            start: UncertainBoundary::certain(pos.clone()),
            end: UncertainBoundary::certain(pos),
        }
    }

    /// Create an interval with simple uncertainty (Mu<T> for each position)
    pub fn with_uncertainty(start: Mu<T>, end: Mu<T>) -> Self {
        Self {
            start: UncertainBoundary::Single(start),
            end: UncertainBoundary::Single(end),
        }
    }

    /// Create an interval with complex boundaries
    pub fn with_complex_boundaries(start: UncertainBoundary<T>, end: UncertainBoundary<T>) -> Self {
        Self { start, end }
    }

    /// Check if this interval has any complex (range) boundaries
    pub fn has_complex_boundaries(&self) -> bool {
        self.start.is_range() || self.end.is_range()
    }
}

impl<T: fmt::Display + PartialEq> fmt::Display for Interval<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // For simple single-position boundaries, check if it's a point variant
        match (&self.start, &self.end) {
            (UncertainBoundary::Single(start_mu), UncertainBoundary::Single(end_mu)) => {
                // Check if both are unknown (?) - this is a point at unknown position
                if start_mu.is_unknown() && end_mu.is_unknown() {
                    return write!(f, "{}", self.start);
                }
                match (start_mu.inner(), end_mu.inner()) {
                    (Some(start), Some(end))
                        if start == end && start_mu.is_certain() && end_mu.is_certain() =>
                    {
                        // Point variant - show single position
                        write!(f, "{}", self.start)
                    }
                    _ => write!(f, "{}_{}", self.start, self.end),
                }
            }
            // Complex boundaries or mixed - always show both
            _ => write!(f, "{}_{}", self.start, self.end),
        }
    }
}

/// Genomic interval (g. coordinates)
pub type GenomeInterval = Interval<GenomePos>;

/// CDS interval (c. coordinates)
pub type CdsInterval = Interval<CdsPos>;

/// Transcript interval (n. coordinates)
pub type TxInterval = Interval<TxPos>;

/// RNA interval (r. coordinates)
pub type RnaInterval = Interval<RnaPos>;

/// Protein interval (p. coordinates)
pub type ProtInterval = Interval<ProtPos>;

impl GenomeInterval {
    /// Get the length of this interval (returns 0 if either position is unknown)
    pub fn len(&self) -> u64 {
        match (self.start.inner(), self.end.inner()) {
            (Some(start), Some(end)) => {
                if end.base >= start.base {
                    // Use saturating_add to prevent overflow at u64::MAX
                    (end.base - start.base).saturating_add(1)
                } else {
                    0
                }
            }
            _ => 0, // Unknown positions result in length 0
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genome_interval_display() {
        let interval = GenomeInterval::new(GenomePos::new(100), GenomePos::new(200));
        assert_eq!(format!("{}", interval), "100_200");
    }

    #[test]
    fn test_point_interval_display() {
        let interval = GenomeInterval::point(GenomePos::new(100));
        assert_eq!(format!("{}", interval), "100");
    }

    #[test]
    fn test_cds_interval_display() {
        let interval = CdsInterval::new(CdsPos::new(100), CdsPos::with_offset(100, 5));
        assert_eq!(format!("{}", interval), "100_100+5");
    }

    #[test]
    fn test_interval_length() {
        let interval = GenomeInterval::new(GenomePos::new(100), GenomePos::new(105));
        assert_eq!(interval.len(), 6);
    }

    #[test]
    fn test_uncertain_interval() {
        let interval = GenomeInterval::with_uncertainty(
            Mu::Uncertain(GenomePos::new(100)),
            Mu::Uncertain(GenomePos::new(200)),
        );
        assert_eq!(format!("{}", interval), "(100)_(200)");
    }
}
