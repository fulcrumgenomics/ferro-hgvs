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
//! # Pairwise levels, or a groupable key
//!
//! [`EquivalenceChecker`] answers a rich question about **two** variants, and
//! its answer is a level rather than a value — `AccessionVersionDifference` is
//! decided by comparing *two* accessions' versions, so there is no scalar a
//! bucket could be keyed on. A consumer that wants to *count distinct changes*
//! across a whole call set needs the other shape: one value per variant, equal
//! exactly when the bases are. That is [`SpdiKey`]. Equal keys cover the
//! `Identical`, `NormalizedMatch` and `SequenceMatch` rungs at once — including
//! the case `SequenceMatch` exists for, where two partitions of one edit reach
//! distinct canonical forms, as the spec's own `DNA/delins.md:44-47` pair does.
//! The key does **not** carry `AccessionVersionDifference`, and it has one
//! documented residual where equal bases key unequally; see [`spdi_key`].
//!
//! Note the implication runs one way. Equal keys imply one of those three rungs,
//! but a pair *at* one of them need not key at all: two identical `p.`
//! descriptions are `Identical` and have no key, because [`spdi_key`] refuses
//! rather than guess. [`spdi_key`] enumerates every such class.
//!
//! # References
//!
//! - [HGVS Nomenclature](https://hgvs-nomenclature.org/)
//! - [Variant Normalization](https://www.ncbi.nlm.nih.gov/variation/notation/)

mod checker;
mod key;

pub use checker::{EquivalenceChecker, EquivalenceLevel, EquivalenceResult};
pub use key::{group_by_spdi_key, spdi_key, SpdiKey, SpdiKeyGrouping};
