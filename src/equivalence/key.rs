//! A groupable, encoding-invariant key for "do these denote the same bases?".
//!
//! # Why this exists
//!
//! An HGVS descriptor asserts a **partition** of the reference, and this
//! project's normalizer preserves that partition rather than re-deriving a
//! minimal one from the resulting sequence. A deliberate consequence is that two
//! spellings which partition the same bases differently reach *different*
//! canonical forms — `g.8_14delinsGATTA` and `g.[8A>G;9G>A;11C>T;13_14del]` stay
//! apart on purpose.
//!
//! That leaves a consumer who only wants to bucket variants — "how many distinct
//! changes are in this call set?" — with no usable handle. Comparing normalized
//! strings answers a different question, and it couples the consumer's stored
//! buckets to the normalizer, so every future confluence fix churns them.
//!
//! The `canonical-form-choice-when-both-legal` ruling record
//! (`tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`) names this
//! module's job in the same breath as the design it is a consequence of:
//! "Sequence-level equivalence still needs an answer for consumers who dedupe;
//! that is `EquivalenceLevel::SequenceMatch` and a groupable SPDI key, not the
//! canonical string." This is that key.
//!
//! [`crate::spdi::canonical_spdi`] already computes the right thing: apply every
//! member to the reference, then read the difference back out of the resulting
//! *bases*, where partitioning is not represented. This module is the grouping
//! surface over it — a key type that is `Ord` + `Hash` so it can be a map key,
//! a total `Option`-returning accessor, and a bucketing helper.
//!
//! ```
//! use ferro_hgvs::equivalence::{group_by_spdi_key, spdi_key};
//! use ferro_hgvs::{parse_hgvs, MockProvider};
//!
//! let mut provider = MockProvider::new();
//! provider.add_genomic_sequence("NC_KEY.1", "GGATTACAGGCATTAGCCTGA");
//!
//! // One edit, two partitions of it. They are two canonical forms by design …
//! let spanning = parse_hgvs("NC_KEY.1:g.3_7delinsGGCTA").unwrap();
//! let decomposed = parse_hgvs("NC_KEY.1:g.[3A>G;4T>G;5T>C;6A>T;7C>A]").unwrap();
//! assert_ne!(spanning.to_string(), decomposed.to_string());
//!
//! // … and one bucket.
//! assert_eq!(
//!     spdi_key(&spanning, &provider),
//!     spdi_key(&decomposed, &provider),
//! );
//!
//! let grouped = group_by_spdi_key([&spanning, &decomposed], &provider);
//! assert_eq!(grouped.group_count(), 1);
//! assert!(grouped.unkeyable().is_empty());
//! ```
//!
//! # What the key is, and what it is not
//!
//! **It is** equal exactly when two descriptions on one accession produce the
//! same resulting bases, whatever their spelling, member count or member order.
//! It is derived without consulting the normalizer, so it does not move when a
//! canonical form does.
//!
//! **It is not** a canonical SPDI in NCBI's sense, and not a substitute for
//! [`EquivalenceChecker`](crate::equivalence::EquivalenceChecker) — that answers
//! a richer pairwise question (accession-version drift, for instance) which no
//! single scalar key can carry. See [`crate::spdi::apply`]'s module docs for the
//! exact guarantee.
//!
//! # What has no key
//!
//! [`spdi_key`] returns `None` rather than a key that would silently collide.
//! **This is the load-bearing half of the contract**: a key that answered for
//! everything would bucket unrelated variants together, which is worse than
//! refusing. See [`spdi_key`] for the enumerated classes.

use std::borrow::Borrow;
use std::collections::BTreeMap;
use std::fmt;

use crate::hgvs::variant::HgvsVariant;
use crate::reference::ReferenceProvider;
use crate::spdi::{canonical_spdi, SpdiVariant};

/// A groupable key that is equal exactly when two descriptions denote the same
/// bases on the same accession.
///
/// Rendered as an SPDI expression (`accession:position:deletion:insertion`, with
/// a 0-based interbase `position`), which is what a consumer stores. Construct
/// it with [`spdi_key`]; there is deliberately no public constructor from a
/// [`SpdiVariant`], because an arbitrary triple is a *transliteration* of one
/// spelling and carries none of the invariance the type name claims.
///
/// `Ord` is over `(accession, position, deletion, insertion)` — an ordering that
/// exists so the key can index a `BTreeMap` and so a report can be printed
/// deterministically. It sorts positionally within one accession, which is
/// convenient, but nothing about the key's meaning depends on the order.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SpdiKey {
    spdi: SpdiVariant,
}

impl SpdiKey {
    /// The accession the key is on.
    ///
    /// Two keys on different accessions never compare equal, so this is the
    /// axis along which a consumer partitions its buckets before comparing them.
    pub fn accession(&self) -> &str {
        &self.spdi.sequence
    }

    /// The key's structured form, for a caller who wants the coordinates rather
    /// than the string.
    pub fn as_spdi(&self) -> &SpdiVariant {
        &self.spdi
    }

    /// Consume the key and yield its structured form.
    pub fn into_spdi(self) -> SpdiVariant {
        self.spdi
    }
}

impl fmt::Display for SpdiKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.spdi)
    }
}

impl PartialOrd for SpdiKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SpdiKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Field by field rather than through a cloned tuple: a key is a map key,
        // so `cmp` runs on every lookup and must not allocate.
        self.spdi
            .sequence
            .cmp(&other.spdi.sequence)
            .then_with(|| self.spdi.position.cmp(&other.spdi.position))
            .then_with(|| self.spdi.deletion.cmp(&other.spdi.deletion))
            .then_with(|| self.spdi.insertion.cmp(&other.spdi.insertion))
    }
}

/// The grouping key for `variant`, or `None` when no single key describes it.
///
/// Two variants whose keys are `Some` and equal denote the same bases on the
/// same accession. Two whose keys are `Some` and unequal denote different bases.
/// `None` is **not** an equivalence class: two `None`s say nothing about each
/// other, and a caller must not bucket them together.
///
/// # What returns `None`
///
/// Every case below is one where a key would have to be guessed. Refusing keeps
/// the invariant that equal keys mean equal bases.
///
/// - **Protein descriptions** (`p.`). SPDI has no protein axis, and the amino
///   acid consequence of an edit is not recoverable from a nucleotide key. This
///   is a permanent exclusion, not a gap: a consumer bucketing protein variants
///   needs a protein key, which this is not.
/// - **Descriptions naming more than one molecule** — a trans, mosaic, chimeric
///   or unknown-phase allele, and the null (`0`) and unknown (`?`) alleles.
///   There is no single resulting sequence to read a key out of.
/// - **Alleles spanning more than one accession.** A key is per-accession.
/// - **Members that overlap**, or two coincident insertions at one interbase
///   position. Applying them depends on an order the description does not state,
///   so there is no one answer. (HGVS spells a stated order as a single ordered
///   payload, `ins[A;C]`, which *does* key.)
/// - **A stated deletion that disagrees with the reference** — the description
///   and the provider cannot both be right, and picking one would key the
///   variant off bases it does not describe.
/// - **Edits SPDI cannot represent.** Inserted-range payloads
///   (`g.10_11ins20_30`), uncertain or unspecified inserted lengths (`ins(10)`),
///   named elements (`insAluYb8`), and intronic `c.`/`n.`/`r.` offsets, which
///   SPDI has no notation for.
/// - **Reference bases the provider cannot serve** — an absent accession, or a
///   window it declines. A key derived from a truncated window is a *different*
///   key, not a smaller one.
/// - **Spans past the library's bounds** — wider than
///   [`MAX_APPLY_WINDOW`](crate::spdi::MAX_APPLY_WINDOW), or sitting in a repeat
///   tract still running past [`MAX_SHIFT_TRACT`](crate::spdi::MAX_SHIFT_TRACT)
///   bases 3', where how far the change rolls depends on how much reference was
///   read.
///
/// One residual is worth stating because it is a `Some` and not a `None`: a
/// description that changes nothing (`g.3A>A`) keys at its own position, so two
/// no-ops at different positions do not group. That is inherited from
/// [`canonical_spdi`] and is documented there.
///
/// Use [`canonical_spdi`] directly when the *reason* for a refusal is wanted;
/// it returns the same answer as a `Result` carrying a described error.
pub fn spdi_key<P: ReferenceProvider + ?Sized>(
    variant: &HgvsVariant,
    provider: &P,
) -> Option<SpdiKey> {
    canonical_spdi(variant, provider)
        .ok()
        .map(|spdi| SpdiKey { spdi })
}

/// Variants bucketed by [`SpdiKey`], plus the ones that have no key.
///
/// The unkeyable list is a first-class part of the result rather than a silent
/// drop. A consumer reading only `groups()` would under-count its call set and
/// have no way to notice; reading `unkeyable()` tells it exactly how much of the
/// input the buckets do not describe.
#[derive(Debug, Clone)]
pub struct SpdiKeyGrouping<T> {
    groups: BTreeMap<SpdiKey, Vec<T>>,
    unkeyable: Vec<T>,
}

impl<T> SpdiKeyGrouping<T> {
    /// Each distinct key and the variants that share it, in key order.
    pub fn groups(&self) -> impl Iterator<Item = (&SpdiKey, &[T])> {
        self.groups
            .iter()
            .map(|(key, items)| (key, items.as_slice()))
    }

    /// The variants [`spdi_key`] declined, in input order.
    ///
    /// These are **not** one bucket — see [`spdi_key`] for why `None` is not an
    /// equivalence class.
    pub fn unkeyable(&self) -> &[T] {
        &self.unkeyable
    }

    /// How many distinct keys the input reduced to.
    pub fn group_count(&self) -> usize {
        self.groups.len()
    }

    /// The variants sharing `key`, or `None` if no input reached it.
    pub fn get(&self, key: &SpdiKey) -> Option<&[T]> {
        self.groups.get(key).map(Vec::as_slice)
    }

    /// Whether every keyable input landed in a single bucket.
    ///
    /// Vacuously `true` for an input with no keyable members, which is why
    /// [`unkeyable`](Self::unkeyable) has to be read alongside it: "all one
    /// bucket" over zero buckets is not evidence of agreement.
    pub fn is_confluent(&self) -> bool {
        self.groups.len() <= 1
    }

    /// The buckets and the refusals, for a caller that wants to own them.
    pub fn into_parts(self) -> (BTreeMap<SpdiKey, Vec<T>>, Vec<T>) {
        (self.groups, self.unkeyable)
    }
}

/// Bucket `variants` by the bases they denote.
///
/// Accepts anything that borrows an [`HgvsVariant`] — owned variants, `&`
/// variants, or a smart pointer — and hands the same items back inside the
/// buckets, so a caller keeps whatever it put in.
///
/// Each variant is keyed independently; one that cannot be keyed lands in
/// [`SpdiKeyGrouping::unkeyable`] rather than failing the call, because a
/// consumer bucketing a call set wants the buckets it *can* have plus an honest
/// count of what it could not place.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::equivalence::group_by_spdi_key;
/// use ferro_hgvs::{parse_hgvs, MockProvider};
///
/// let mut provider = MockProvider::new();
/// provider.add_genomic_sequence("NC_KEY.1", "GGATTACAGGCATTAGCCTGA");
///
/// let variants = [
///     parse_hgvs("NC_KEY.1:g.13_14insT").unwrap(),
///     parse_hgvs("NC_KEY.1:g.14dup").unwrap(),
///     parse_hgvs("NC_KEY.1:g.3A>G").unwrap(),
/// ];
/// let grouped = group_by_spdi_key(&variants, &provider);
/// assert_eq!(grouped.group_count(), 2);
/// ```
pub fn group_by_spdi_key<T, I, P>(variants: I, provider: &P) -> SpdiKeyGrouping<T>
where
    I: IntoIterator<Item = T>,
    T: Borrow<HgvsVariant>,
    P: ReferenceProvider + ?Sized,
{
    let mut groups: BTreeMap<SpdiKey, Vec<T>> = BTreeMap::new();
    let mut unkeyable = Vec::new();
    for variant in variants {
        match spdi_key(variant.borrow(), provider) {
            Some(key) => groups.entry(key).or_default().push(variant),
            None => unkeyable.push(variant),
        }
    }
    SpdiKeyGrouping { groups, unkeyable }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;
    use crate::reference::MockProvider;

    /// The 40-base contig `spdi::apply`'s own tests use, so a key asserted here
    /// can be checked against the bases by hand.
    ///
    ///  1-based: 1234567890...
    ///           GGATTACAGGCATTAGCCTGAGGATTACAGGCATTAGCCT
    fn provider() -> MockProvider {
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_KEY.1", "GGATTACAGGCATTAGCCTGAGGATTACAGGCATTAGCCT");
        provider.add_genomic_sequence("NC_OTHER.1", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        provider
    }

    fn key(descriptor: &str) -> Option<SpdiKey> {
        let variant = parse_hgvs(descriptor).expect("fixture must parse");
        spdi_key(&variant, &provider())
    }

    /// The property the module exists for: descriptions that assert *different
    /// partitions* of one edit — and so reach different canonical forms by
    /// design — share one key.
    #[test]
    fn different_partitions_of_one_edit_share_a_key() {
        let spanning = key("NC_KEY.1:g.3_7delinsGGCTA").expect("keys");
        let decomposed = key("NC_KEY.1:g.[3A>G;4T>G;5T>C;6A>T;7C>A]").expect("keys");
        assert_eq!(spanning, decomposed);
        assert_eq!(spanning.to_string(), "NC_KEY.1:2:ATTAC:GGCTA");
        assert_eq!(spanning.accession(), "NC_KEY.1");
    }

    /// The control for the test above. Without it, every equality assertion here
    /// would also pass for a key that collapsed everything into one bucket —
    /// which is the failure mode the `None`-rather-than-collide rule is about.
    #[test]
    fn different_edits_do_not_share_a_key() {
        assert_ne!(key("NC_KEY.1:g.3A>G"), key("NC_KEY.1:g.3A>C"));
        assert_ne!(key("NC_KEY.1:g.3_5del"), key("NC_KEY.1:g.3_6del"));
        // Same edit, different accession: a key is per-accession.
        assert_ne!(key("NC_KEY.1:g.3A>G"), key("NC_OTHER.1:g.3T>G"));
    }

    /// The key is not the normalizer's output, so it groups across the axes the
    /// canonical form deliberately keeps apart: member order, and a `dup` versus
    /// the insertion that denotes it.
    #[test]
    fn the_key_ignores_spelling_and_member_order() {
        assert_eq!(key("NC_KEY.1:g.[3A>G;7C>A]"), key("NC_KEY.1:g.[7C>A;3A>G]"));
        assert_eq!(key("NC_KEY.1:g.13_14insT"), key("NC_KEY.1:g.14dup"));
        // The same deletion spelled one base apart inside a tract.
        assert_eq!(key("NC_KEY.1:g.3_5del"), key("NC_KEY.1:g.4_6del"));
    }

    /// Each documented refusal, exercised as the thing the doc claims.
    ///
    /// Written as a table with the reason attached because the value of this
    /// test is the enumeration: a class that quietly started returning a key
    /// would be a silent collision, which is the one failure the type is
    /// supposed to make impossible.
    #[test]
    fn the_documented_classes_have_no_key() {
        for (descriptor, why) in [
            ("NP_000079.2:p.Gly1150Asp", "protein axis — SPDI has none"),
            ("NC_KEY.1:g.[3A>G(;)7C>A]", "trans phase — two molecules"),
            ("NC_KEY.1:g.[3_5del;4T>G]", "overlapping members"),
            (
                "NC_KEY.1:g.[5_6insA;5_6insC]",
                "coincident insertions with no stated order",
            ),
            (
                "[NC_KEY.1:g.3A>G;NC_OTHER.1:g.7T>G]",
                "members on two accessions",
            ),
            (
                "NC_KEY.1:g.3T>G",
                "stated base disagrees with the reference",
            ),
            ("NC_ABSENT.1:g.3A>G", "the provider cannot serve the bases"),
            ("NC_KEY.1:g.10_11ins20_30", "inserted-range payload"),
            ("NC_KEY.1:g.10_11ins(10)", "unspecified inserted length"),
        ] {
            assert_eq!(
                key(descriptor),
                None,
                "`{descriptor}` must have no key ({why}) — a key here is a silent collision"
            );
        }
    }

    /// A member with a stated order in one payload does key, which is what makes
    /// the coincident-insertion refusal a refusal to *guess* rather than a gap.
    #[test]
    fn a_stated_insertion_order_still_keys() {
        assert!(key("NC_KEY.1:g.5_6insAC").is_some());
    }

    #[test]
    fn grouping_separates_buckets_from_refusals() {
        let provider = provider();
        let variants: Vec<HgvsVariant> = [
            "NC_KEY.1:g.3_7delinsGGCTA",
            "NC_KEY.1:g.[3A>G;4T>G;5T>C;6A>T;7C>A]",
            "NC_KEY.1:g.13_14insT",
            "NC_KEY.1:g.[3_5del;4T>G]",
        ]
        .iter()
        .map(|d| parse_hgvs(d).expect("fixture must parse"))
        .collect();

        let grouped = group_by_spdi_key(&variants, &provider);
        assert_eq!(grouped.group_count(), 2);
        assert_eq!(grouped.unkeyable().len(), 1);
        assert!(!grouped.is_confluent());

        let (largest, members) = grouped
            .groups()
            .max_by_key(|(_, members)| members.len())
            .expect("two buckets");
        assert_eq!(members.len(), 2, "the two partitions of one edit");
        assert_eq!(grouped.get(largest).map(<[_]>::len), Some(2));
    }

    /// `into_parts` hands back exactly what went in, so a caller can own the
    /// buckets without re-deriving them.
    #[test]
    fn into_parts_preserves_every_input() {
        let provider = provider();
        let variants: Vec<HgvsVariant> = ["NC_KEY.1:g.3A>G", "NC_ABSENT.1:g.3A>G"]
            .iter()
            .map(|d| parse_hgvs(d).expect("fixture must parse"))
            .collect();
        let (groups, unkeyable) = group_by_spdi_key(variants.clone(), &provider).into_parts();
        let returned = groups.values().map(Vec::len).sum::<usize>() + unkeyable.len();
        assert_eq!(returned, variants.len());
    }

    /// Ordering exists so the key can index a `BTreeMap`; pin that it is
    /// positional within an accession rather than an accident of the string.
    #[test]
    fn keys_order_positionally_within_an_accession() {
        let mut keys = [
            key("NC_KEY.1:g.7C>A").expect("keys"),
            key("NC_KEY.1:g.3A>G").expect("keys"),
        ];
        keys.sort();
        assert!(keys[0].as_spdi().position < keys[1].as_spdi().position);
        assert_eq!(keys[0].clone().into_spdi().position, 2);
    }
}
