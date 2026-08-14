//! A groupable, encoding-invariant key for "do these denote the same bases?".
//!
//! # Why this exists
//!
//! An HGVS descriptor asserts a **partition** of the reference, and two
//! spellings that partition the same bases differently can reach *different*
//! canonical forms. That is not a gap waiting to be closed. The decided ruling
//! `canonical-form-choice-when-both-legal`
//! (`tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`) directs
//! ferro to derive from the resulting sequence "subject to every explicit spec
//! tie-break" — and where two such tie-breaks point opposite ways, neither is
//! selected and both forms survive. The spec supplies its own example:
//! `DNA/delins.md:44` gives a spanning `delins`, `:46` names a four-member split
//! as an "alternative description" of the *same* variant, `:47` recommends the
//! spanning form, and `general.md:34` requires members separated by unchanged
//! nucleotides to stay individual.
//!
//! That leaves a consumer who only wants to bucket variants — "how many distinct
//! changes are in this call set?" — with no usable handle. Comparing normalized
//! strings answers a different question, and it couples the consumer's stored
//! buckets to the normalizer: confluence work moves canonical forms by design,
//! and the project declares those moves as breaking rather than avoiding them,
//! so buckets keyed on the string re-bucket on every such release.
//!
//! This is the handle that does not move: one value per variant, equal exactly
//! when the denoted bases are, derived without consulting the normalizer.
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
//! // One edit, two partitions of it — two distinct descriptions …
//! let spanning = parse_hgvs("NC_KEY.1:g.3_7delinsGGCTA").unwrap();
//! let decomposed = parse_hgvs("NC_KEY.1:g.[3A>G;4T>G;5T>C;6A>T;7C>A]").unwrap();
//! assert_ne!(spanning.to_string(), decomposed.to_string());
//!
//! // … and one bucket. Note the key is read off the denoted bases, so this
//! // holds whether or not the normalizer converges the pair — no `Normalizer`
//! // is involved above.
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
//! Read "the same bases" as case-insensitive, and on the `r.` axis as
//! `U`-insensitive: the key's bases are folded before it is built, so a
//! soft-masked FASTA and an uppercase one give one key for one locus even though
//! their windows differ byte for byte. That fold is what makes the key a
//! property of the sequence rather than of the provider; see
//! [`crate::spdi::apply`]'s module docs for the guarantee and the one boundary it
//! does not reach.
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
/// convenient, but nothing about the key's meaning depends on the order. It is
/// delegated to [`SpdiVariant`]'s derived `Ord`, whose declaration order *is*
/// that tuple, so a field added there joins the comparison rather than being
/// silently dropped from it.
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
        // Delegated to `SpdiVariant`'s *derived* `Ord` rather than written out
        // field by field. The derive is lexicographic in declaration order, which
        // is the `(sequence, position, deletion, insertion)` order documented on
        // the type — so this is the same comparison, and it stays the same
        // comparison when a field is added. Spelling the four fields here type-
        // checked and agreed with `Eq` on every value that existed at the time,
        // and would have gone on doing both after a fifth field appeared while
        // quietly no longer being total: two keys differing only in the new field
        // would compare `Equal` and collapse into one `BTreeMap` bucket while
        // `Eq` called them different. Neither allocates.
        self.spdi.cmp(&other.spdi)
    }
}

/// The grouping key for `variant`, or `None` when no single key describes it.
///
/// Two variants whose keys are `Some` and equal denote the same bases on the
/// same accession. Two whose keys are `Some` and unequal denote different bases —
/// **with one stated exception, the no-op residual below**, which is the only
/// place this direction does not hold. `None` is **not** an equivalence class:
/// two `None`s say nothing about each other, and a caller must not bucket them
/// together.
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
/// - **Edits SPDI cannot represent.** Uncertain or unspecified inserted lengths
///   (`ins(10)`), named elements (`insAluYb8`), and intronic `c.`/`n.`/`r.`
///   offsets, which SPDI has no notation for. (A same-reference position-range
///   insert — `g.10_11ins20_30`, `ins20_30inv` — *does* key: its bases are read
///   from the reference exactly as `del`/`dup` read their omitted bases.)
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
/// no-ops at different positions do not group. They denote the same bases — the
/// unchanged reference — while keying unequally, which is why the second sentence
/// above is qualified: this is the one class where unequal keys do *not* imply
/// different bases. That is inherited from [`canonical_spdi`] and is documented
/// there (`g.3A>A` keys as `3::`, `g.5T>T` as `5::`); giving the two one answer
/// means choosing one, and nothing here needs it.
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
            ("NC_KEY.1:g.10_11ins(10)", "unspecified inserted length"),
        ] {
            assert_eq!(
                key(descriptor),
                None,
                "`{descriptor}` must have no key ({why}) — a key here is a silent collision"
            );
        }
    }

    /// The one documented exception to "unequal keys mean different bases": a
    /// description that changes nothing keys at its own position, so two no-ops
    /// at different positions do not group even though both denote the unchanged
    /// reference.
    ///
    /// Pinned here because [`spdi_key`]'s contract is now stated with this as its
    /// sole carve-out. If the residual is ever given one answer — which
    /// [`canonical_spdi`]'s own docs say is a choice nobody has needed to make —
    /// that sentence has to move with it, and a red test here is what says so.
    #[test]
    fn two_no_ops_at_different_positions_do_not_group() {
        let here = key("NC_KEY.1:g.3A>A").expect("a no-op still keys");
        let there = key("NC_KEY.1:g.5T>T").expect("a no-op still keys");
        assert_ne!(here, there, "the residual is stated on `spdi_key`");
        for no_op in [&here, &there] {
            assert_eq!(no_op.as_spdi().deletion, "", "a no-op deletes nothing");
            assert_eq!(no_op.as_spdi().insertion, "", "a no-op inserts nothing");
        }
    }

    /// A member with a stated order in one payload does key, which is what makes
    /// the coincident-insertion refusal a refusal to *guess* rather than a gap.
    #[test]
    fn a_stated_insertion_order_still_keys() {
        assert!(key("NC_KEY.1:g.5_6insAC").is_some());
    }

    /// A same-reference position-range insert names its bases by position, so
    /// once the provider resolves them (as `del`/`dup` resolve their omitted
    /// bases) it keys off exactly those bases. The 40-base contig holds
    /// `GGCATTAGCCT` at 9..=19, so `g.5_6ins9_19` and the literal spelling of
    /// the same bases must reach one key.
    #[test]
    fn a_same_reference_range_insert_keys_off_its_resolved_bases() {
        let ranged = key("NC_KEY.1:g.5_6ins9_19").expect("a resolvable range insert keys");
        let literal = key("NC_KEY.1:g.5_6insGGCATTAGCCT").expect("the literal spelling keys");
        assert_eq!(
            ranged, literal,
            "a range insert must key off the bases it names, not its spelling"
        );
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

    /// A soft-masked reference must key identically to the uppercase one.
    ///
    /// Real genomic FASTAs lowercase repeat regions, so the *same* accession can
    /// be served soft-masked by one provider and uppercase by another. Case
    /// carries no biological meaning, so the two denote the same bases and the
    /// contract — "equal exactly when the denoted bases are" — requires one key.
    ///
    /// Asserted on `Eq`, `Hash` and `Ord` together rather than on `Eq` alone: a
    /// `BTreeMap` bucket is found by `Ord` and a `HashMap` bucket by `Hash`, so a
    /// key that compared equal while hashing or ordering differently would still
    /// split a consumer's buckets, which is the failure the type exists to
    /// prevent.
    #[test]
    fn a_soft_masked_reference_keys_the_same_as_an_uppercase_one() {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};

        const UPPER: &str = "GGATTACAGGCATTAGCCTGAGGATTACAGGCATTAGCCT";

        let hashed = |key: &SpdiKey| {
            let mut hasher = DefaultHasher::new();
            key.hash(&mut hasher);
            hasher.finish()
        };
        let keyed = |sequence: &str, descriptor: &str| {
            let mut provider = MockProvider::new();
            provider.add_genomic_sequence("NC_KEY.1", sequence.to_string());
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            spdi_key(&variant, &provider)
        };

        // Every shape whose key can carry reference-derived bases: a stated
        // deletion, a deletion the description does not spell, a duplication and
        // an inversion (both of which read their inserted bases *out of* the
        // reference), and a delins mixing both directions.
        for descriptor in [
            "NC_KEY.1:g.3_7delinsGGCTA",
            "NC_KEY.1:g.3_5del",
            "NC_KEY.1:g.3_5dup",
            "NC_KEY.1:g.3_7inv",
            "NC_KEY.1:g.3A>G",
            "NC_KEY.1:g.5_6insAC",
        ] {
            let upper = keyed(UPPER, descriptor).expect("the uppercase reference keys");
            let lower = keyed(&UPPER.to_ascii_lowercase(), descriptor)
                .expect("the soft-masked reference keys");
            assert_eq!(
                upper, lower,
                "`{descriptor}` keys differently on a soft-masked reference: \
                 {upper} vs {lower}"
            );
            assert_eq!(
                hashed(&upper),
                hashed(&lower),
                "`{descriptor}` hashes differently on a soft-masked reference, so a \
                 `HashMap` would split the bucket even if `Eq` agreed"
            );
            assert_eq!(
                upper.cmp(&lower),
                std::cmp::Ordering::Equal,
                "`{descriptor}` orders unequally on a soft-masked reference, so a \
                 `BTreeMap` would split the bucket even if `Eq` agreed"
            );
        }
    }

    /// The delegated `Ord` really is the documented tuple order, and it really
    /// does agree with `Eq`.
    ///
    /// `Ord` is delegated to [`SpdiVariant`]'s derive so a future field cannot
    /// fall out of the comparison. That trades a hand-written `cmp` for a
    /// declaration order, so the thing to pin is that the declaration order is
    /// the documented one: each field is varied in isolation, from the least
    /// significant up, and each must dominate every field below it. Ordering
    /// consistency with `Eq` is asserted alongside, because an `Ord` returning
    /// `Equal` for unequal keys would merge two variants into one `BTreeMap`
    /// bucket — a silent collision, which is the one failure this type exists to
    /// make impossible.
    #[test]
    fn ord_follows_the_documented_field_order_and_agrees_with_eq() {
        use std::cmp::Ordering;

        let of = |spdi: SpdiVariant| SpdiKey { spdi };
        let base = SpdiVariant::new("NC_A.1", 10, "AC", "GG");

        // Each entry raises exactly one field of `base`, most significant last.
        // A greater value in a *less* significant field must never outrank one in
        // a more significant field, which is what "lexicographic" buys.
        let ascending = [
            of(SpdiVariant::new("NC_A.1", 10, "AC", "GG")),
            of(SpdiVariant::new("NC_A.1", 10, "AC", "GT")),
            of(SpdiVariant::new("NC_A.1", 10, "AT", "GG")),
            of(SpdiVariant::new("NC_A.1", 11, "AC", "GG")),
            of(SpdiVariant::new("NC_B.1", 10, "AC", "GG")),
        ];
        for (index, lower) in ascending.iter().enumerate() {
            for higher in &ascending[index + 1..] {
                assert_eq!(
                    lower.cmp(higher),
                    Ordering::Less,
                    "{lower} must sort before {higher}"
                );
            }
        }

        // Consistency with `Eq`, in both directions.
        assert_eq!(of(base.clone()).cmp(&of(base.clone())), Ordering::Equal);
        for other in &ascending[1..] {
            assert_ne!(of(base.clone()), *other);
            assert_ne!(
                of(base.clone()).cmp(other),
                Ordering::Equal,
                "unequal keys must not compare equal — a `BTreeMap` would merge them"
            );
        }
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
