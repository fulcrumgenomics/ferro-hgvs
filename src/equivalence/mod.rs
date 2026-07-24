//! Variant equivalence checking.
//!
//! This module provides functionality to determine if two HGVS variants
//! represent the same genomic change, even if expressed differently.
//!
//! # Examples
//!
//! ```
//! use ferro_hgvs::{parse_hgvs, MockProvider};
//! use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};
//!
//! // Create an equivalence checker with test data
//! let provider = MockProvider::with_test_data();
//! let checker = EquivalenceChecker::new(provider);
//!
//! // Check if two variants are equivalent
//! let v1 = parse_hgvs("NM_000088.3:c.10del").unwrap();
//! let v2 = parse_hgvs("NM_000088.3:c.10del").unwrap();
//!
//! let result = checker.check(&v1, &v2).unwrap();
//! assert!(matches!(result.level, EquivalenceLevel::Identical));
//! ```
//!
//! # Equivalence Levels
//!
//! The checker recognizes several levels of equivalence:
//!
//! - **Identical**: Same string representation
//! - **NormalizedMatch**: Same after normalization (e.g., different positions in repeat region)
//! - **SequenceMatch**: Different normalized strings, same resulting sequence when
//!   applied to the reference (e.g. a length-changing `delins` vs a decomposed
//!   cis allele of the same edit)
//! - **AccessionVersionDifference**: Same variant, different accession versions
//! - **NotEquivalent**: Represent different changes
//!
//! The list is open-ended — [`EquivalenceLevel`] is `#[non_exhaustive]`, so
//! recognizing a new class of equivalence adds a level without breaking
//! downstream matches.
//!
//! # References
//!
//! - [HGVS Nomenclature](https://hgvs-nomenclature.org/)
//! - [Variant Normalization](https://www.ncbi.nlm.nih.gov/variation/notation/)

mod checker;

pub use checker::{EquivalenceChecker, EquivalenceLevel, EquivalenceResult};
