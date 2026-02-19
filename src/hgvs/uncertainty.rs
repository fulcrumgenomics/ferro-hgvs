//! Uncertainty wrapper for HGVS positions and edits
//!
//! HGVS allows positions and values to be marked as uncertain using parentheses.
//! For example, `c.(100_200)del` indicates uncertain deletion boundaries.
//! Unknown positions are marked with `?` in HGVS notation.

use serde::{Deserialize, Serialize};
use std::fmt;

/// Wrapper type for values that may be certain, uncertain, or unknown
///
/// In HGVS notation:
/// - Uncertain values are wrapped in parentheses: `(100)`
/// - Unknown values are marked with `?`: `c.?_100del`
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Mu<T> {
    /// Value is certain (not wrapped in parentheses)
    Certain(T),
    /// Value is uncertain (wrapped in parentheses in HGVS notation)
    Uncertain(T),
    /// Value is unknown (represented as `?` in HGVS notation)
    Unknown,
}

impl<T> Mu<T> {
    /// Get the inner value if present (None for Unknown)
    pub fn inner(&self) -> Option<&T> {
        match self {
            Mu::Certain(v) | Mu::Uncertain(v) => Some(v),
            Mu::Unknown => None,
        }
    }

    /// Get the inner value, consuming self (None for Unknown)
    pub fn into_inner(self) -> Option<T> {
        match self {
            Mu::Certain(v) | Mu::Uncertain(v) => Some(v),
            Mu::Unknown => None,
        }
    }

    /// Check if the value is certain (known and not in parentheses)
    pub fn is_certain(&self) -> bool {
        matches!(self, Mu::Certain(_))
    }

    /// Check if the value is uncertain (known but in parentheses)
    pub fn is_uncertain(&self) -> bool {
        matches!(self, Mu::Uncertain(_))
    }

    /// Check if the value is unknown (represented as `?`)
    pub fn is_unknown(&self) -> bool {
        matches!(self, Mu::Unknown)
    }

    /// Check if the value is known (either certain or uncertain)
    pub fn is_known(&self) -> bool {
        !self.is_unknown()
    }

    /// Map over the inner value (preserves Unknown)
    pub fn map<U, F: FnOnce(T) -> U>(self, f: F) -> Mu<U> {
        match self {
            Mu::Certain(v) => Mu::Certain(f(v)),
            Mu::Uncertain(v) => Mu::Uncertain(f(v)),
            Mu::Unknown => Mu::Unknown,
        }
    }

    /// Map over a reference to the inner value, preserving uncertainty status
    ///
    /// This is useful when you want to transform the value but don't want to
    /// consume the original Mu.
    pub fn map_ref<U, F: FnOnce(&T) -> U>(&self, f: F) -> Mu<U> {
        match self {
            Mu::Certain(v) => Mu::Certain(f(v)),
            Mu::Uncertain(v) => Mu::Uncertain(f(v)),
            Mu::Unknown => Mu::Unknown,
        }
    }

    /// Create a certain value
    pub fn certain(value: T) -> Self {
        Mu::Certain(value)
    }

    /// Create an uncertain value
    pub fn uncertain(value: T) -> Self {
        Mu::Uncertain(value)
    }

    /// Create an unknown value
    pub fn unknown() -> Self {
        Mu::Unknown
    }
}

impl<T: fmt::Display> fmt::Display for Mu<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Mu::Certain(v) => write!(f, "{}", v),
            Mu::Uncertain(v) => write!(f, "({})", v),
            Mu::Unknown => write!(f, "?"),
        }
    }
}

impl<T: Default> Default for Mu<T> {
    fn default() -> Self {
        Mu::Certain(T::default())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_certain_value() {
        let mu = Mu::Certain(42);
        assert!(mu.is_certain());
        assert!(!mu.is_uncertain());
        assert!(!mu.is_unknown());
        assert!(mu.is_known());
        assert_eq!(*mu.inner().unwrap(), 42);
    }

    #[test]
    fn test_uncertain_value() {
        let mu = Mu::Uncertain(42);
        assert!(!mu.is_certain());
        assert!(mu.is_uncertain());
        assert!(!mu.is_unknown());
        assert!(mu.is_known());
        assert_eq!(*mu.inner().unwrap(), 42);
    }

    #[test]
    fn test_unknown_value() {
        let mu: Mu<i32> = Mu::Unknown;
        assert!(!mu.is_certain());
        assert!(!mu.is_uncertain());
        assert!(mu.is_unknown());
        assert!(!mu.is_known());
        assert!(mu.inner().is_none());
    }

    #[test]
    fn test_display() {
        assert_eq!(format!("{}", Mu::Certain(100)), "100");
        assert_eq!(format!("{}", Mu::Uncertain(100)), "(100)");
        assert_eq!(format!("{}", Mu::<i32>::Unknown), "?");
    }

    #[test]
    fn test_map() {
        let mu = Mu::Certain(10);
        let doubled = mu.map(|x| x * 2);
        assert_eq!(*doubled.inner().unwrap(), 20);
        assert!(doubled.is_certain());
    }

    #[test]
    fn test_map_unknown() {
        let mu: Mu<i32> = Mu::Unknown;
        let mapped = mu.map(|x| x * 2);
        assert!(mapped.is_unknown());
        assert!(mapped.inner().is_none());
    }

    #[test]
    fn test_into_inner() {
        assert_eq!(Mu::Certain(42).into_inner(), Some(42));
        assert_eq!(Mu::Uncertain(42).into_inner(), Some(42));
        assert_eq!(Mu::<i32>::Unknown.into_inner(), None);
    }
}
