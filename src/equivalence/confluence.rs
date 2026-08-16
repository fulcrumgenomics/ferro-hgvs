//! An opt-in confluence self-check a consumer can run on its own data (#1892).
//!
//! # What confluence is here
//!
//! The property the project states is "`normalize` is constant on each
//! equivalence class" — the ruling
//! `confluence-gate-is-apply-equality-on-every-determined-axis`. The equivalence
//! class is defined by `apply` (description → bases), never by `normalize`, so
//! the statement is not circular: two descriptions that are apply-equal must
//! normalize to one string, and a class that reaches two is a **non-confluence
//! witness**.
//!
//! # Why a consumer needs to run it on their own data
//!
//! The designed cis corpus is enriched for multi-member shapes, so its
//! non-confluence rate is a rate *of that shape family* and tells a consumer
//! nothing about the variants they actually hold. This self-check answers the
//! question they can only answer for themselves: over *my* call set, are any two
//! apply-equal descriptions normalized apart?
//!
//! # What it is, and what it deliberately is not
//!
//! It **composes** two pieces of decided machinery — the equivalence relation
//! ([`EquivalenceChecker`](crate::equivalence::EquivalenceChecker)) and the
//! normalizer — and reports. It groups the inputs into equivalence classes under
//! a caller-chosen [`ConfluenceRelation`], normalizes each member, and records
//! every class that reaches more than one distinct normalized output.
//!
//! It is a **diagnostic**, not an error and not a release gate:
//!
//! * It reports offending groups so a consumer can count them and sample in
//!   production; making it an error would make it unusable in a pipeline.
//! * It does **not** decide *which* relation gates a release — that is the
//!   deferred adjudication in #1890. The relation is the caller's parameter, and
//!   both variants of [`ConfluenceRelation`] name relations that are already
//!   decided and documented in this module. The check emits no pass/fail verdict.
//!
//! # Choosing the relation
//!
//! [`ConfluenceRelation::CrossAxisSequenceMatch`] is the relation the confluence
//! ruling names: apply-equality on **every determined axis**
//! ([`EquivalenceLevel::CrossAxisSequenceMatch`](crate::equivalence::EquivalenceLevel::CrossAxisSequenceMatch)),
//! the rung that establishes variant identity. It is the right default for
//! "are two descriptions of *one variant* normalized apart?".
//!
//! [`ConfluenceRelation::Spdi`] groups by the encoding-invariant
//! [`SpdiKey`](crate::equivalence::SpdiKey) — same bases on the description's own
//! axis. It is **weaker** and *insufficient for variant identity*: it buckets
//! `c.3921dup` with `c.3922dup`, which denote the same transcript bases but
//! project ~2,790 bp apart into different exons. It is offered because a consumer
//! counting distinct edits may want it and because #1892 names it, never as a
//! confluence gate.
//!
//! # Cost: the default relation is quadratic, the SPDI one is not
//!
//! The two relations are not interchangeable on cost, and the difference is a
//! complexity class rather than a constant. [`ConfluenceRelation::Spdi`] keys
//! each input **once** and buckets it, so it is `O(n)` key derivations.
//! [`ConfluenceRelation::CrossAxisSequenceMatch`] has no key — apply-equality is
//! only computable pairwise — so it compares each input against every class
//! representative so far, which is `O(n²)`
//! [`compare_denotations`](crate::equivalence::EquivalenceChecker::compare_denotations)
//! calls in the worst case (a corpus of distinct variants, which is the usual
//! shape of a call set). Each call does its own SPDI conversion and reference
//! fetches.
//!
//! **The ratio is the claim; the milliseconds are not.** Two independent debug
//! builds over `n` distinct genomic substitutions agree on the shape and not on
//! the constants, so read the growth and re-measure the absolutes rather than
//! quoting them. One of the two, `n` = 50 / 100 / 200 / 400: cross-axis
//! 18.7 / 60.0 / 218.9 / 1032.2 ms against SPDI 5.1 / 3.6 / 11.0 / 15.0 ms —
//! that is roughly 3–5x per doubling on the left and no consistent growth at
//! this scale on the right. The corpus is the worst case by construction: every
//! input is a distinct variant, so no class ever absorbs another and every
//! comparison is paid.
//!
//! **There is deliberately no cap and no automatic downgrade.** A cap would have
//! to decide *which* pairs to stop comparing, and every such choice silently
//! shrinks the classes — which is the failure this module exists to prevent (see
//! [`ConfluenceSkip`]): a corpus whose comparisons were truncated reports fewer
//! violations, not more skips. Downgrading to [`ConfluenceRelation::Spdi`] on
//! size would answer a *different question* under the same name, and that
//! relation is documented as insufficient for variant identity. So the cost is
//! stated here and left to the caller, who is the only party that can choose
//! between "the identity relation over a sample" and "a weaker relation over
//! everything". For a call-set-scale run, batch by accession (classes never span
//! accessions under either relation) or sample.
//!
//! ```
//! use ferro_hgvs::equivalence::{ConfluenceRelation, EquivalenceChecker};
//! use ferro_hgvs::{parse_hgvs, MockProvider};
//!
//! let mut provider = MockProvider::new();
//! provider.add_genomic_sequence("NC_KEY.1", "GGATTACAGGCATTAGCCTGA");
//! let checker = EquivalenceChecker::new(provider);
//!
//! // Two spellings of one insertion — apply-equal, and the normalizer converges
//! // them, so the corpus is confluent.
//! let corpus = [
//!     parse_hgvs("NC_KEY.1:g.13_14insT").unwrap(),
//!     parse_hgvs("NC_KEY.1:g.14dup").unwrap(),
//! ];
//! let report = checker.check_confluence(&corpus, ConfluenceRelation::CrossAxisSequenceMatch);
//! assert!(report.is_confluent());
//! ```

use super::checker::TripleDecline;
use std::collections::BTreeSet;

/// The equivalence relation a [`check_confluence`] run groups its inputs by.
///
/// Confluence is "`normalize` is constant on each equivalence class"; this names
/// which relation defines the class. Both variants are **decided, documented**
/// relations already in this crate — the caller chooses; the self-check does not
/// decide which relation gates a release (#1890).
///
/// `#[non_exhaustive]` so a future relation is not a breaking change, mirroring
/// [`EquivalenceLevel`](crate::equivalence::EquivalenceLevel).
///
/// [`check_confluence`]: crate::equivalence::EquivalenceChecker::check_confluence
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ConfluenceRelation {
    /// Apply-equality on **every determined axis**
    /// ([`EquivalenceLevel::CrossAxisSequenceMatch`](crate::equivalence::EquivalenceLevel::CrossAxisSequenceMatch)) —
    /// the relation the confluence ruling names and the rung that establishes
    /// variant identity. Grouped pairwise via
    /// [`EquivalenceChecker::compare_denotations`](crate::equivalence::EquivalenceChecker::compare_denotations),
    /// which never consults the normalizer, so the grouping cannot collapse into
    /// the outputs it is checking.
    CrossAxisSequenceMatch,
    /// Same denoted bases on the description's own axis, via the
    /// encoding-invariant [`SpdiKey`](crate::equivalence::SpdiKey). **Weaker**
    /// than the above and insufficient for variant identity; use it for counting
    /// distinct edits, never as a confluence gate.
    Spdi,
}

impl ConfluenceRelation {
    /// The stable name this relation is spelled with outside Rust — the Python
    /// binding's `relation=` argument and the `relation` field of a rendered
    /// report.
    ///
    /// Exhaustive on purpose: this crate's own matches on a `#[non_exhaustive]`
    /// enum are still checked for exhaustiveness (the attribute only binds
    /// downstream crates), so adding a relation fails the build **here**, next
    /// to [`Self::from_name`], rather than in a binding that silently keeps
    /// answering for the relations it already knows.
    pub fn as_str(self) -> &'static str {
        match self {
            ConfluenceRelation::CrossAxisSequenceMatch => "cross_axis",
            ConfluenceRelation::Spdi => "spdi",
        }
    }

    /// The inverse of [`Self::as_str`], or `None` for an unknown name.
    ///
    /// Deliberately co-located with `as_str`: this direction matches on a
    /// `&str`, so no compiler check can force a new relation to gain an arm.
    /// The pairing is pinned by `every_relation_round_trips_through_its_name`.
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "cross_axis" => Some(ConfluenceRelation::CrossAxisSequenceMatch),
            "spdi" => Some(ConfluenceRelation::Spdi),
            _ => None,
        }
    }
}

impl std::fmt::Display for ConfluenceRelation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// One equivalence class whose members normalized to more than one distinct
/// output — a non-confluence witness.
///
/// A confluent class produces no `ConfluenceGroup`; see
/// [`ConfluenceGroup::from_class`].
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct ConfluenceGroup {
    /// The class members' rendered input descriptions, in input order.
    pub inputs: Vec<String>,
    /// The **distinct** normalized outputs the class reached, sorted. Length is
    /// always at least two — a class with one output is confluent and yields no
    /// group.
    pub outputs: Vec<String>,
}

impl ConfluenceGroup {
    /// Build a violation from one equivalence class's `(input, normalized_output)`
    /// members, or `None` when the class is **confluent** — its members reach at
    /// most one distinct normalized output.
    ///
    /// This is the flagging rule in isolation: it deduplicates and sorts the
    /// outputs and returns `Some` exactly when more than one survives. The inputs
    /// are recorded in the order given so a consumer can see which descriptions
    /// collided.
    pub(crate) fn from_class(members: impl IntoIterator<Item = (String, String)>) -> Option<Self> {
        let mut inputs = Vec::new();
        let mut outputs = BTreeSet::new();
        for (input, output) in members {
            inputs.push(input);
            outputs.insert(output);
        }
        if outputs.len() > 1 {
            Some(ConfluenceGroup {
                inputs,
                outputs: outputs.into_iter().collect(),
            })
        } else {
            None
        }
    }
}

/// Which stage of a [`check_confluence`] run dropped an input, for a consumer
/// that needs to account for coverage rather than read a message.
///
/// The two are **not** interchangeable, and reading them as one is how a
/// coverage figure goes wrong — see [`ConfluenceReport::classes_checked`].
///
/// [`check_confluence`]: crate::equivalence::EquivalenceChecker::check_confluence
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ConfluenceSkipKind {
    /// The chosen relation could not place the input at all, so it is in **no**
    /// class and is not counted in
    /// [`classes_checked`](ConfluenceReport::classes_checked).
    Unplaceable,
    /// The input *was* placed in a class — so its class **is** counted in
    /// [`classes_checked`](ConfluenceReport::classes_checked) — but the
    /// normalizer declined it, so it contributed no output to that class's
    /// output set.
    NormalizationDeclined,
}

impl ConfluenceSkipKind {
    /// The stable name this kind is spelled with outside Rust.
    ///
    /// Exhaustive for the reason given on [`ConfluenceRelation::as_str`].
    pub fn as_str(self) -> &'static str {
        match self {
            ConfluenceSkipKind::Unplaceable => "unplaceable",
            ConfluenceSkipKind::NormalizationDeclined => "normalization_declined",
        }
    }
}

impl std::fmt::Display for ConfluenceSkipKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// An input the confluence analysis could not fully cover, with why.
///
/// A skip is **not** a pass. It is surfaced as a first-class part of the report
/// so a consumer knows exactly how much of its call set the analysis does not
/// cover — a silent drop that read as confluence is the failure this check
/// exists to prevent. Two things land here, kept apart by
/// [`kind`](Self::kind): an input the chosen relation cannot place
/// ([`ConfluenceSkipKind::Unplaceable`] — no computable denotation on this
/// reference, so no [`SpdiKey`](crate::equivalence::SpdiKey) under
/// [`ConfluenceRelation::Spdi`] and no comparable triples under
/// [`ConfluenceRelation::CrossAxisSequenceMatch`]), and an input the normalizer
/// declined ([`ConfluenceSkipKind::NormalizationDeclined`]).
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct ConfluenceSkip {
    /// The rendered input description that was skipped.
    pub input: String,
    /// Which stage dropped it.
    pub kind: ConfluenceSkipKind,
    /// The refusal, **typed**, when one was produced: which class it is in
    /// ([`TripleDecline`]) and which site it came from
    /// ([`DeclineSite`](crate::equivalence::DeclineSite)).
    ///
    /// This is the field to act on. It separates "this description denotes no
    /// sequence to anyone" from "this reference could not serve the bases" —
    /// which read very differently to a consumer with a partly-provisioned
    /// reference — as data rather than as English, so a consumer never has to
    /// sniff [`reason`](Self::reason) for a keyword that a rewording silently
    /// removes.
    ///
    /// `None` only where there is no such refusal to report: a
    /// [`NormalizationDeclined`](ConfluenceSkipKind::NormalizationDeclined) skip,
    /// which is the normalizer's refusal rather than the projection's, and the
    /// [`Spdi`](ConfluenceRelation::Spdi) first-pass mismatch described on
    /// `group_by_spdi`.
    pub decline: Option<TripleDecline>,
    /// Why it could not be analyzed, rendered for a human. A presentation of
    /// [`decline`](Self::decline) where there is one — never the authority on it.
    pub reason: String,
}

/// The result of a [`check_confluence`] run.
///
/// A **diagnostic**: it reports the non-confluence witnesses it found, how many
/// classes it examined, and what it could not place. It carries no pass/fail
/// release verdict; [`is_confluent`](Self::is_confluent) is an observation about
/// *this corpus under this relation*, not a gate.
///
/// [`check_confluence`]: crate::equivalence::EquivalenceChecker::check_confluence
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct ConfluenceReport {
    relation: ConfluenceRelation,
    violations: Vec<ConfluenceGroup>,
    classes_checked: usize,
    skipped: Vec<ConfluenceSkip>,
    undecided_pairs: usize,
}

impl ConfluenceReport {
    /// Construct a report. Crate-internal: a report is only meaningful when its
    /// parts come from one [`check_confluence`] run.
    ///
    /// [`check_confluence`]: crate::equivalence::EquivalenceChecker::check_confluence
    pub(crate) fn new(
        relation: ConfluenceRelation,
        violations: Vec<ConfluenceGroup>,
        classes_checked: usize,
        skipped: Vec<ConfluenceSkip>,
        undecided_pairs: usize,
    ) -> Self {
        Self {
            relation,
            violations,
            classes_checked,
            skipped,
            undecided_pairs,
        }
    }

    /// The relation the run grouped by.
    pub fn relation(&self) -> ConfluenceRelation {
        self.relation
    }

    /// The non-confluence witnesses: every examined class that reached more than
    /// one distinct normalized output.
    pub fn violations(&self) -> &[ConfluenceGroup] {
        &self.violations
    }

    /// How many equivalence classes the run examined, singletons included.
    ///
    /// The denominator a violation count is read against.
    ///
    /// **It is not `inputs - skipped`, and computing coverage as
    /// `classes_checked + skipped.len()` double-counts.** The two
    /// [`ConfluenceSkipKind`]s are excluded differently, which is the whole
    /// reason that enum exists:
    ///
    /// * [`Unplaceable`](ConfluenceSkipKind::Unplaceable) inputs were never
    ///   placed, so they are in no class and contribute nothing here;
    /// * [`NormalizationDeclined`](ConfluenceSkipKind::NormalizationDeclined)
    ///   inputs **were** placed — their class is counted here — and then
    ///   contributed no output to it. A class of two whose only two members both
    ///   decline still counts as one examined class, over an empty output set.
    ///
    /// So the inputs this run actually placed is
    /// `inputs - skipped.filter(kind == Unplaceable).count()`, and the classes
    /// among them is `classes_checked`.
    ///
    /// This paragraph replaces a claim that skips "were never placed in a
    /// class", which was true of one kind and false of the other.
    pub fn classes_checked(&self) -> usize {
        self.classes_checked
    }

    /// The inputs that could not be analyzed, each with a kind and a reason.
    ///
    /// Read this alongside [`is_confluent`](Self::is_confluent): "no violation"
    /// over a corpus that was mostly skipped is not evidence of confluence.
    pub fn skipped(&self) -> &[ConfluenceSkip] {
        &self.skipped
    }

    /// How many pairwise comparisons came back **undecidable** rather than
    /// deciding the pair apart — always `0` under
    /// [`ConfluenceRelation::Spdi`], which compares keys and never a pair.
    ///
    /// This is the *second* way coverage is lost, and it is invisible in
    /// [`skipped`](Self::skipped) because nothing was dropped: both inputs were
    /// placed, they merely failed to merge. Under
    /// [`ConfluenceRelation::CrossAxisSequenceMatch`] a pair whose comparison
    /// declined — the union reference window exceeds the checker's cap, or the
    /// relation errored on a shape it cannot decide — does not merge, exactly as
    /// a pair that was compared and found different does not merge. Grouping
    /// therefore **splits a class it could not examine into singletons**, and
    /// every non-confluence inside that class becomes invisible.
    ///
    /// So a non-zero count means the class structure is a *lower bound* on how
    /// coarse the true classes are, and [`is_confluent`](Self::is_confluent) is
    /// a *lower bound* on the non-confluence present. Zero is what makes
    /// "no violations" a statement about the corpus rather than about the
    /// reference.
    pub fn undecided_pairs(&self) -> usize {
        self.undecided_pairs
    }

    /// Whether **this corpus** showed no non-confluence under **this relation**.
    ///
    /// A diagnostic observation, never a release verdict — and vacuously `true`
    /// for a corpus with no examined classes, which is why
    /// [`skipped`](Self::skipped), [`classes_checked`](Self::classes_checked)
    /// and [`undecided_pairs`](Self::undecided_pairs) must all be read with it.
    /// `true` with a non-empty `skipped` or a non-zero `undecided_pairs` means
    /// "no non-confluence *among the pairs I could examine*", which is a weaker
    /// claim and the one the report is actually entitled to make.
    pub fn is_confluent(&self) -> bool {
        self.violations.is_empty()
    }

    /// Whether the run examined everything it was given: nothing skipped and no
    /// comparison undecidable.
    ///
    /// [`is_confluent`](Self::is_confluent) is only a statement about the whole
    /// corpus when this is `true`. It says nothing about whether violations were
    /// found — a complete run can still report plenty.
    pub fn is_complete(&self) -> bool {
        self.skipped.is_empty() && self.undecided_pairs == 0
    }

    /// The violations and the skips, for a caller that wants to own them.
    pub fn into_parts(self) -> (Vec<ConfluenceGroup>, Vec<ConfluenceSkip>) {
        (self.violations, self.skipped)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn members(pairs: &[(&str, &str)]) -> Vec<(String, String)> {
        pairs
            .iter()
            .map(|(i, o)| (i.to_string(), o.to_string()))
            .collect()
    }

    /// A class whose members all reach one output is confluent, so no group is
    /// produced — whatever the member count.
    #[test]
    fn a_class_with_one_distinct_output_is_confluent() {
        assert_eq!(
            ConfluenceGroup::from_class(members(&[("g.2A>G", "g.2A>G")])),
            None,
            "a singleton is confluent"
        );
        assert_eq!(
            ConfluenceGroup::from_class(members(&[
                ("g.2A>G", "g.2A>G"),
                ("g.2delinsG", "g.2A>G"),
                ("g.2_2delinsG", "g.2A>G"),
            ])),
            None,
            "three spellings, one output"
        );
    }

    /// A class whose members reach two distinct outputs is a violation, and the
    /// outputs are deduplicated and sorted while the inputs keep their order.
    #[test]
    fn a_class_with_two_distinct_outputs_is_a_violation() {
        let group = ConfluenceGroup::from_class(members(&[
            ("input-b", "g.2_7delinsGGCTA"),
            ("input-a", "g.[2A>G;7C>A]"),
            ("input-c", "g.2_7delinsGGCTA"),
        ]))
        .expect("two distinct outputs is non-confluent");

        assert_eq!(
            group.inputs,
            vec!["input-b", "input-a", "input-c"],
            "inputs are recorded in the order given"
        );
        assert_eq!(
            group.outputs,
            // Lexicographic: '2' (0x32) sorts before '[' (0x5B).
            vec!["g.2_7delinsGGCTA", "g.[2A>G;7C>A]"],
            "outputs are the distinct set, sorted"
        );
    }

    /// The boundary: exactly two members, two outputs, is the smallest witness.
    #[test]
    fn two_members_two_outputs_is_the_minimal_witness() {
        let group =
            ConfluenceGroup::from_class(members(&[("a", "x"), ("b", "y")])).expect("two outputs");
        assert_eq!(group.outputs, vec!["x", "y"]);
    }

    /// An empty class produces nothing rather than a spurious violation.
    #[test]
    fn an_empty_class_is_not_a_violation() {
        assert_eq!(ConfluenceGroup::from_class(std::iter::empty()), None);
    }

    /// A report exposes its parts through the accessors, and `into_parts` hands
    /// back exactly the violations and skips it was built from.
    #[test]
    fn a_report_exposes_its_parts() {
        let violation = ConfluenceGroup {
            inputs: vec!["a".to_string(), "b".to_string()],
            outputs: vec!["x".to_string(), "y".to_string()],
        };
        let skip = ConfluenceSkip {
            input: "c".to_string(),
            kind: ConfluenceSkipKind::Unplaceable,
            reason: "no SPDI key".to_string(),
            decline: None,
        };
        let report = ConfluenceReport::new(
            ConfluenceRelation::Spdi,
            vec![violation.clone()],
            3,
            vec![skip.clone()],
            7,
        );

        assert_eq!(report.relation(), ConfluenceRelation::Spdi);
        assert_eq!(report.classes_checked(), 3);
        assert_eq!(report.violations(), std::slice::from_ref(&violation));
        assert_eq!(report.skipped(), std::slice::from_ref(&skip));
        assert_eq!(report.undecided_pairs(), 7);
        assert!(!report.is_confluent(), "one violation means not confluent");
        assert!(!report.is_complete(), "a skip and undecided pairs are gaps");

        let (violations, skipped) = report.into_parts();
        assert_eq!(violations, [violation]);
        assert_eq!(skipped, [skip]);
    }

    /// A run with nothing skipped and nothing undecidable is complete, whatever
    /// it found — completeness is about coverage, not about the verdict.
    #[test]
    fn completeness_is_about_coverage_not_about_the_verdict() {
        let clean = ConfluenceReport::new(
            ConfluenceRelation::CrossAxisSequenceMatch,
            Vec::new(),
            2,
            Vec::new(),
            0,
        );
        assert!(clean.is_complete() && clean.is_confluent());

        let found = ConfluenceReport::new(
            ConfluenceRelation::CrossAxisSequenceMatch,
            vec![ConfluenceGroup {
                inputs: vec!["a".to_string(), "b".to_string()],
                outputs: vec!["x".to_string(), "y".to_string()],
            }],
            2,
            Vec::new(),
            0,
        );
        assert!(
            found.is_complete() && !found.is_confluent(),
            "a complete run can still report violations"
        );

        let blind = ConfluenceReport::new(
            ConfluenceRelation::CrossAxisSequenceMatch,
            Vec::new(),
            2,
            Vec::new(),
            1,
        );
        assert!(
            !blind.is_complete() && blind.is_confluent(),
            "one undecidable comparison and no skip is still an incomplete run"
        );
    }

    /// Every relation's external name round-trips. The rendering direction is
    /// compiler-checked (an exhaustive match on the enum), the parsing
    /// direction cannot be, so the pairing is pinned here — a relation added to
    /// `as_str` without an arm in `from_name` fails this test.
    #[test]
    fn every_relation_round_trips_through_its_name() {
        for relation in [
            ConfluenceRelation::CrossAxisSequenceMatch,
            ConfluenceRelation::Spdi,
        ] {
            let name = relation.as_str();
            assert!(!name.is_empty(), "{relation:?} must have a name");
            assert_eq!(
                ConfluenceRelation::from_name(name),
                Some(relation),
                "{name:?} must parse back to {relation:?}"
            );
            assert_eq!(relation.to_string(), name, "Display matches as_str");
        }
        assert_eq!(ConfluenceRelation::from_name("cross-axis"), None);
        assert_eq!(ConfluenceRelation::from_name(""), None);
    }

    /// The two skip kinds have distinct names, and they round-trip through
    /// `Display` — a consumer reading a rendered report must be able to tell
    /// "never placed" from "placed, then declined".
    #[test]
    fn the_two_skip_kinds_are_distinguishable() {
        let unplaceable = ConfluenceSkipKind::Unplaceable;
        let declined = ConfluenceSkipKind::NormalizationDeclined;
        assert_ne!(unplaceable.as_str(), declined.as_str());
        assert_eq!(unplaceable.to_string(), unplaceable.as_str());
        assert_eq!(declined.to_string(), declined.as_str());
    }
}
