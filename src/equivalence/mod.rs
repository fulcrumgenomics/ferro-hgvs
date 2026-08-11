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
//! # The relation, stated without circularity
//!
//! `normalize` maps description → description, so it cannot define the
//! equivalence relation it is supposed to respect. The relation is defined by
//! `apply`, which maps description → **sequence**, where equality is byte
//! equality on bases:
//!
//! ```text
//! equivalent(a, b)  ≝  same accession
//!                   ∧  determined(a) = determined(b)
//!                   ∧  ∀ X ∈ determined(a) :  apply_X(a) = apply_X(b)
//! ```
//!
//! An axis is **determined** when the description's coordinates can be carried
//! to it by a mapping that is a function of reference data alone — exon
//! alignments, CDS offsets — and never of normalization. A `c.`/`n.`/`r.`
//! description determines two axes: the transcript, and the genome its exon
//! alignment carries it to. A `g.`/`m.` description determines only the genome.
//! Protein is excluded: translation is many-to-one, and `p.` states a
//! consequence rather than a denotation.
//!
//! # Equivalence Levels
//!
//! The checker recognizes several levels of equivalence:
//!
//! - **Identical**: Same string representation
//! - **NormalizedMatch**: Same after normalization (e.g., different positions in repeat region)
//! - **CrossAxisSequenceMatch**: Same resulting sequence on **every** determined
//!   axis. This is the rung that establishes variant identity, and the one a
//!   confluence gate requires
//! - **SequenceMatch**: Different normalized strings, same resulting sequence
//!   on the description's own axis (e.g. a length-changing `delins` vs a
//!   decomposed cis allele of the same edit). True, and *insufficient for
//!   identity* — `LRG_199t1:c.3921dup` and `c.3922dup` reach it while denoting
//!   genomic changes ~2,790 bp apart, in different exons
//!   (`DNA/duplication.md:148`)
//! - **AccessionVersionDifference**: Same variant, different accession versions
//! - **NotEquivalent**: Represent different changes
//! - **Indeterminate**: Neither — at least one side has no computable
//!   denotation. Not a negative verdict; see [`EquivalenceLevel::is_decided`]
//!
//! [`EquivalenceLevel::is_at_least`] expresses the order over the denotational
//! rungs, and deliberately makes **NormalizedMatch unusable as a gate**: that
//! rung consults the normalizer, so gating on it would restore the circularity
//! above.
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
//! Where the two genuinely disagree is the **axis count**, and this module's
//! own counterexample is the witness: [`spdi_key`] buckets `c.3921dup`
//! **with** `c.3922dup` — they are transcript-apply-equal, so their canonical
//! transcript SPDIs coincide — while
//! [`EquivalenceLevel::CrossAxisSequenceMatch`] separates them. A grouping
//! built on the key answers "how many distinct edits are in this call set?";
//! only the cross-axis rung answers "are these the same variant?".
//!
//! # References
//!
//! - [HGVS Nomenclature](https://hgvs-nomenclature.org/)
//! - [Variant Normalization](https://www.ncbi.nlm.nih.gov/variation/notation/)

mod checker;
mod key;

pub use checker::{EquivalenceChecker, EquivalenceLevel, EquivalenceResult};
pub use key::{group_by_spdi_key, spdi_key, SpdiKey, SpdiKeyGrouping};
