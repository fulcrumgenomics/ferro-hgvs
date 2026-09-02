// Copyright (c) 2024-2025 Fulcrum Genomics LLC
// SPDX-License-Identifier: MIT

//! A **runnable** conformance census over the spec corpus (#1890).
//!
//! # Why this is a library module and not just a test
//!
//! The census — validity, confluence, idempotency, sequence preservation,
//! refusal and the coding-axis instrument, over the spec-derived corpus — used to
//! live only inside `tests/it/spec_conformance_axis.rs`, where its numbers are
//! *pinned*. Pinning is the right thing for a gate, but it means the only way to
//! see a number, diff two revisions, or produce a burn-down was to edit a test.
//!
//! This module is the measurement extracted so it can be *called*. The gate still
//! owns the pins and remains the gate; it consumes [`measure`] with
//! [`Equivalence::ExactString`] so its figures do not move. The three companion
//! pieces designed alongside #1890 — the release-to-release stability harness
//! (#1891) and the consumer-side self-check (#1892) — consume the same
//! measurement, so there is one definition of "conformance" rather than three
//! that drift.
//!
//! # What it deliberately does NOT do
//!
//! - **It states no release verdict.** Which equivalence relation *gates* a
//!   release is an adjudication left open by #1890 on purpose; this module reports
//!   confluence under whichever relation the caller names ([`Equivalence`], no
//!   default at the CLI) and never collapses the census to a single pass/fail or a
//!   single percentage. Rank-1 validity and rank-2 confluence have moved in
//!   opposite directions in one change, so an average would mis-report a
//!   correctness improvement as a regression.
//! - **It refuses a flattering artifact.** A panicking row contributes no output,
//!   so with a seam oracle armed a family silently converges and the census reads
//!   *better* than the truth. [`run_census`] detects `FERRO_ASSERT_*` in the
//!   environment and refuses **itself**, rather than leaving the check to each
//!   caller — a refusal only one of two entry points performs is not a property of
//!   the census.
//! - **It accounts for what it dropped.** The corpus drops designs by
//!   construction, and a real equivalence relation cannot key every output (SPDI
//!   has no offset notation, so intronic outputs have no key). Both populations
//!   are counted — the first through [`CaptureLedger`], under a ceiling that is a
//!   real bound ([`max_dropped_designs`]), the second as
//!   [`CensusReport::confluence_unkeyable_outputs`] — so a partial run cannot read
//!   as a complete one, and every reported zero carries its denominator (the
//!   [`CorpusShape`] travels beside the census).
//! - **It exports what it counted, in full.** The human summary caps the
//!   divergence list for readability and says by how much; `--findings` does not
//!   cap it. A machine-readable export that silently truncates is worse than one
//!   that omits the field, because a consumer diffing two of them sees the same
//!   head of the same list and reads "no change".

use std::collections::{BTreeMap, BTreeSet};

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::completeness::{Allowance, CaptureCounts, CaptureLedger, Shortfall};
use super::spec_corpus::{
    denotation_of, enumerate, CorpusBounds, Denotation, Row, RowKind, SpecCorpus, Strength,
};
use crate::equivalence::spdi_key;
use crate::error_handling::ErrorMode;
use crate::hgvs::interval::{Interval, UncertainBoundary};
use crate::hgvs::location::{CdsPos, TxPos};
use crate::reference::MockProvider;
use crate::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer, ShuffleDirection};

// ---------------------------------------------------------------------------
// The equivalence relation confluence is measured under
// ---------------------------------------------------------------------------

/// The relation under which two normalized outputs count as "the same output"
/// when deciding confluence.
///
/// The corpus builds each family by generating spellings that *re-partition* one
/// variant, so "did every spelling reach ONE output?" depends on what "one" means
/// — which is a different question for a consumer deduping a call set, a
/// maintainer asking whether the normalizer is deterministic, and two tools
/// comparing. The census reports whichever the caller asks for; #1890 leaves
/// which one *gates* a release deliberately unsettled.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Equivalence {
    /// Byte-identical normalized strings. The strictest relation, and the one the
    /// pinned axis test uses so its census does not move; it is **not** offered on
    /// the CLI, whose three relations are the consumer-facing questions below.
    ExactString,
    /// Same denoted bases. Answers "can a consumer dedupe by normalizing?".
    /// Backed by [`spdi_key`] over the whole allele.
    Sequence,
    /// Same denoted bases **and** the same member partition. Answers "is the
    /// normalizer deterministic within the ruling?". The strictest of the three.
    Partition,
    /// The set of per-member SPDI keys matches. Answers "do two tools agree?".
    Spdi,
}

impl Equivalence {
    /// Stable label, used in the artifact and the report.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::ExactString => "exact-string",
            Self::Sequence => "sequence",
            Self::Partition => "partition",
            Self::Spdi => "spdi",
        }
    }

    /// The grouping key for one normalized output under this relation, or `None`
    /// when the relation cannot key it (an intronic offset, an overlapping
    /// member — anything [`spdi_key`] declines).
    ///
    /// [`Self::ExactString`] always keys (the string is its own key), which is
    /// what makes the pinned axis census reproducible: no output is ever
    /// unkeyable, so the only undecidable families are those with fewer than two
    /// normalized spellings, exactly as the gate has always counted them.
    fn key_of(self, output: &str, provider: &MockProvider) -> Option<String> {
        match self {
            Self::ExactString => Some(output.to_string()),
            Self::Sequence => {
                let variant = parse_hgvs(output).ok()?;
                spdi_key(&variant, provider).map(|key| key.to_string())
            }
            Self::Spdi => {
                let variant = parse_hgvs(output).ok()?;
                let mut keys = member_keys(&variant, provider)?;
                keys.dedup(); // a SET of member keys
                Some(keys.join("|"))
            }
            Self::Partition => {
                let variant = parse_hgvs(output).ok()?;
                let whole = spdi_key(&variant, provider)?.to_string();
                let members = member_keys(&variant, provider)?; // a multiset: partition cares about count
                Some(format!("{whole}#{}", members.join("|")))
            }
        }
    }
}

/// The sorted SPDI keys of `variant`'s members, or `None` if any member has no
/// key (which makes the whole output unkeyable under a member-wise relation).
fn member_keys(variant: &HgvsVariant, provider: &MockProvider) -> Option<Vec<String>> {
    let mut keys: Vec<String> = members_of(variant)
        .iter()
        .map(|member| spdi_key(member, provider).map(|key| key.to_string()))
        .collect::<Option<Vec<_>>>()?;
    keys.sort();
    Some(keys)
}

// ---------------------------------------------------------------------------
// The census
// ---------------------------------------------------------------------------

/// The corpus's own shape, independent of any property measured over it.
///
/// Pinned separately from the census so a generator change that shrank the corpus
/// cannot read as an improvement (the #1460 hazard). Stamped into every artifact
/// beside the census so a reported zero always carries its denominator.
#[derive(Debug, Default, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CorpusShape {
    /// Designs the enumeration proposed, before any filtering.
    pub designs_considered: usize,
    /// Rows that survived the ground-truth filter.
    pub rows: usize,
    /// Spellings across every row — the number of normalizations a run performs.
    pub spellings: usize,
    /// Family rows (two or more spellings of one variant).
    pub family_rows: usize,
    /// Single-spelling rows.
    pub single_rows: usize,
    /// Conflict rows (a conflicting allele the normalizer should refuse).
    pub conflict_rows: usize,
    /// Prohibited-shape rows.
    pub prohibited_rows: usize,
    /// Rows with more than one authored member — the priority axis.
    pub multi_member_rows: usize,
    /// Rows carrying a negative guard (the rejected-SVD-WG010 shape).
    pub guarded_rows: usize,
    /// Coding-axis rows separated by two or more unchanged nucleotides — the
    /// coding-axis merge instrument's population.
    pub coding_axis_separation_two_or_more_rows: usize,
}

impl CorpusShape {
    /// Read the shape off a built corpus.
    #[must_use]
    pub fn of(built: &SpecCorpus) -> Self {
        let by_kind = built.by_kind();
        Self {
            designs_considered: built.designs_considered,
            rows: built.rows.len(),
            spellings: built.spellings(),
            family_rows: by_kind.get(&RowKind::Family).copied().unwrap_or(0),
            single_rows: by_kind.get(&RowKind::Single).copied().unwrap_or(0),
            conflict_rows: by_kind.get(&RowKind::Conflict).copied().unwrap_or(0),
            prohibited_rows: by_kind.get(&RowKind::Prohibited).copied().unwrap_or(0),
            multi_member_rows: built.multi_member_rows(),
            guarded_rows: built.by_negative_guard().values().sum(),
            coding_axis_separation_two_or_more_rows: built
                .coding_axis_separation_two_or_more_rows(),
        }
    }
}

/// What one direction's run found. Every field is a count over a partition of the
/// rows, so [`Census::absorb`] over corpus groups is exact.
///
/// `pub` fields, serde-serializable, so the gate can pin it, a run can emit it,
/// and a consumer can diff it.
///
/// # Where the narrative lives
///
/// **Here.** What each figure means, and which way it is allowed to move, is on
/// the fields below — this type moved out of `tests/it/spec_conformance_axis.rs`
/// and the field docs came with it, because a doc comment that describes a field
/// belongs beside the field. The gate keeps what is genuinely its own: the
/// **pinned values**, the row-by-row adjudication of why a particular pin holds
/// the number it holds, and the ruling records those adjudications cite. So a
/// question of the form "what is this counter?" is answered here, and one of the
/// form "why is this counter 3?" is answered there.
#[derive(Debug, Default, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Census {
    /// The error mode the **normalizer** behind every figure below was built with
    /// (#1629). Serialized as a string so the artifact says which mode it was
    /// taken under — a census compared against one taken under a different mode
    /// is a category error.
    ///
    /// **It does not reach the parse half.** Every spelling is read with the bare
    /// `parse_hgvs`, which constructs no `InputPreprocessor` and so applies no
    /// `ErrorConfig` at all — not lenient's repairs, not strict's refusals (#1632,
    /// pinned by `issue_1632_parse_entry_applies_no_mode.rs`). So
    /// [`Self::unparseable_outputs`] and [`Self::outputs`] are mode-*less*
    /// figures, and only the counters downstream of normalization are figures for
    /// this mode. "No mode" is a third thing rather than a synonym for strict,
    /// which is exactly why it has to be said instead of inferred from this field.
    ///
    /// Set from [`measurement_config`], never written by hand, and carried in the
    /// same `assert_eq!(census, pinned)` as every other field — so changing what a
    /// run measures under fails on this field, naming the mode, instead of
    /// silently re-basing every number beside it.
    ///
    /// [`Default`] gives `Strict` (the enum's own default) rather than the lenient
    /// value the measurement actually uses, which is deliberate: an accumulator
    /// that was never stamped must not read as a valid claim.
    #[serde(with = "error_mode_serde")]
    pub measured_under: ErrorMode,
    /// Spellings **attempted**, across every row — incremented before `parse_hgvs`
    /// runs, so it counts every spelling the run reached.
    ///
    /// Deliberately not "spellings normalized": a spelling that fails to parse, or
    /// that normalization declines, is counted here **and** in
    /// [`Self::declined`]. The two are therefore **NOT disjoint buckets**, and
    /// [`report`]'s `"{} outputs, {} declined"` line must not be read as though
    /// they were — subtracting one from the other is wrong the moment `declined`
    /// leaves zero.
    ///
    /// Counting attempts rather than successes is the useful choice: it makes
    /// `outputs == CorpusShape::spellings` a cross-check between the shape and the
    /// census, which a success-only counter could not provide.
    pub outputs: usize,
    /// Spellings normalization declined or panicked on, where an output was
    /// expected. Counted apart from a divergence, so one cannot hide the other,
    /// and **not** disjoint from [`Self::outputs`] — see that field.
    pub declined: usize,
    /// Outputs `parse_hgvs` rejects. Rank-1 invalid.
    pub unparseable_outputs: usize,
    /// Outputs whose members claim overlapping territory. Rank-1 invalid.
    pub outputs_denoting_no_sequence: usize,
    /// Outputs naming an intronic offset on a bare transcript accession
    /// (`checklist.md:20`). Rank-1 invalid.
    pub outputs_leaving_the_transcript: usize,
    /// Outputs naming an intronic offset under the genomic wrapper the clause asks
    /// for. Conformant — not a defect counter.
    pub outputs_intronic_under_a_genomic_wrapper: usize,
    /// Outputs violating an absolute textual prohibition. Rank-1 invalid.
    ///
    /// **Ratchet-pinned by the gate rather than asserted at zero**, even though it
    /// currently reads zero — 32 until #1627, whose 24 × `standards.md:39` rows are
    /// refused rather than re-emitted; 8 until #1628, whose 4 × `checklist.md:16`
    /// `+` offsets and 4 × `checklist.md:45` hyphen ranges are refused too.
    ///
    /// Reaching zero is not the same as the class being unreachable, and this doc
    /// has already gone stale in the other direction once: it said "asserted at
    /// zero" while both censuses pinned **32**, which reads as "ferro emits no
    /// prohibited output" to anyone who does not go and check the pin. Quote the
    /// pin, not this prose, if the two ever disagree again — a count in a doc
    /// comment is the thing that goes stale.
    pub prohibition_violating_outputs: usize,
    /// Families whose spellings all reached one output, under the chosen relation.
    pub converged: usize,
    /// Families reaching exactly two distinct outputs.
    pub split_two: usize,
    /// … exactly three.
    pub split_three: usize,
    /// … four or more.
    pub split_more: usize,
    /// Families confluence could not be decided over — fewer than two spellings
    /// normalized, or an output the relation could not key.
    pub underdetermined: usize,
    /// Outputs that are not their own fixed point.
    pub non_idempotent_outputs: usize,
    /// Outputs whose applied sequence differs from their row's.
    pub sequence_changed: usize,
    /// Conflicting alleles the implementation accepted instead of refusing.
    pub conflicts_accepted: usize,
    /// Shapes the spec prohibits in as many words that were accepted.
    pub prohibited_absolute_accepted: usize,
    /// Shapes a conditional clause prohibits that were accepted.
    pub prohibited_conditional_accepted: usize,
    /// Outputs implementing the rejected SVD-WG010 proposal.
    pub guard_violations: usize,
    /// The denominator [`Self::coding_axis_separation_two_or_more_merges`] is `n
    /// of` — coding-axis rows separated by two or more unchanged nucleotides the
    /// run actually evaluated.
    pub coding_axis_separation_two_or_more_rows: usize,
    /// Coding-axis multi-member alleles that normalized to fewer members than were
    /// authored. An **instrument**, not a verdict — see the gate's docs.
    pub coding_axis_separation_two_or_more_merges: usize,
}

impl Census {
    /// Fold another census in.
    ///
    /// Every field below is a count over a partition of the rows, so summing is
    /// exact rather than approximate — that is what makes [`measure_corpus`]'s
    /// per-group split safe.
    ///
    /// **The destructuring is the point, not style.** It carries no `..`, so adding
    /// a field to [`Census`] makes this function fail to COMPILE until the field is
    /// folded. An earlier version listed the fields by hand and #1710 added two
    /// underneath it; they were silently never summed, and the census reported
    /// `coding_axis_separation_two_or_more_merges` as **0 of 0** — a denominator of
    /// zero, which is the flattering direction the VACUOUS guards exist to catch.
    /// This makes the next one a build error instead. [`compare_reports`] is built
    /// the same way, for the same reason.
    ///
    /// [`Census::measured_under`] is bound and deliberately NOT summed: it is a
    /// mode stamp rather than a counter, and every group is built from the same
    /// [`measurement_config`], so the accumulator's own stamp already carries it.
    /// Binding it keeps it inside the exhaustiveness check.
    fn absorb(&mut self, other: &Census) {
        let Census {
            measured_under: _,
            outputs,
            declined,
            unparseable_outputs,
            outputs_denoting_no_sequence,
            outputs_leaving_the_transcript,
            outputs_intronic_under_a_genomic_wrapper,
            prohibition_violating_outputs,
            converged,
            split_two,
            split_three,
            split_more,
            underdetermined,
            non_idempotent_outputs,
            sequence_changed,
            conflicts_accepted,
            prohibited_absolute_accepted,
            prohibited_conditional_accepted,
            guard_violations,
            coding_axis_separation_two_or_more_rows,
            coding_axis_separation_two_or_more_merges,
        } = other;
        self.outputs += outputs;
        self.declined += declined;
        self.unparseable_outputs += unparseable_outputs;
        self.outputs_denoting_no_sequence += outputs_denoting_no_sequence;
        self.outputs_leaving_the_transcript += outputs_leaving_the_transcript;
        self.outputs_intronic_under_a_genomic_wrapper += outputs_intronic_under_a_genomic_wrapper;
        self.prohibition_violating_outputs += prohibition_violating_outputs;
        self.converged += converged;
        self.split_two += split_two;
        self.split_three += split_three;
        self.split_more += split_more;
        self.underdetermined += underdetermined;
        self.non_idempotent_outputs += non_idempotent_outputs;
        self.sequence_changed += sequence_changed;
        self.conflicts_accepted += conflicts_accepted;
        self.prohibited_absolute_accepted += prohibited_absolute_accepted;
        self.prohibited_conditional_accepted += prohibited_conditional_accepted;
        self.guard_violations += guard_violations;
        self.coding_axis_separation_two_or_more_rows += coding_axis_separation_two_or_more_rows;
        self.coding_axis_separation_two_or_more_merges += coding_axis_separation_two_or_more_merges;
    }
}

/// One divergence: a family whose spellings did not all reach one output.
///
/// Serialized with a `record` discriminator, because the `--findings` export
/// interleaves these with [`Finding`] rows in one JSONL stream and the two shapes
/// are otherwise told apart only by which optional key happens to be present. A
/// consumer key-sniffing a union type is a consumer that breaks the first time
/// either shape gains a field.
#[derive(Debug, Clone, Serialize)]
#[serde(tag = "record", rename = "divergence")]
pub struct Divergence {
    /// The row id.
    pub id: String,
    /// Its distinct normalized outputs.
    pub outputs: Vec<String>,
}

/// A finding worth naming: a validity failure, an accepted prohibition, or a
/// negative-guard violation.
///
/// Carries the same `record` discriminator as [`Divergence`], and for the same
/// reason.
#[derive(Debug, Clone, Serialize)]
#[serde(tag = "record", rename = "finding")]
pub struct Finding {
    /// The row id.
    pub id: String,
    /// What was found.
    pub what: String,
}

/// The result of one direction's run.
#[derive(Debug, Clone)]
pub struct Measured {
    /// The census.
    pub census: Census,
    /// The relation confluence was decided under.
    pub equivalence: Equivalence,
    /// **Every** divergent family, in corpus order.
    ///
    /// Uncapped on purpose. [`MAX_DIVERGENCES`] caps the *human summary* only,
    /// because a panic message naming 721 families is unreadable; it does not cap
    /// this list, and it did until #2063's review found that the machine-readable
    /// `--findings` export inherited the message's cap. The consumer that export
    /// exists for (#1891) diffs two releases to name **which** families moved, and
    /// two truncated files hold the same corpus-order-first 12 whatever moved
    /// behind them — so a release in which 400 families changed read as no change.
    /// A cap on a summary is a courtesy; the same cap on an export is data loss.
    pub divergences: Vec<Divergence>,
    /// Every finding, in corpus order.
    pub findings: Vec<Finding>,
    /// Outputs the relation could not key, dropped from confluence and counted so
    /// a partial confluence decision cannot read as a clean one.
    pub confluence_unkeyable_outputs: usize,
}

/// How many divergent families the **human summary** ([`report`]) names before it
/// says how many more there are.
///
/// It is a readability bound on prose, and nothing else. In particular it is not
/// applied to [`Measured::divergences`], which the `--findings` export writes: see
/// that field.
const MAX_DIVERGENCES: usize = 12;

// ---------------------------------------------------------------------------
// Output predicates (written from the AST, independent of the normalizer)
// ---------------------------------------------------------------------------

/// Whether `output` violates a prohibition the spec states in as many words.
fn violated_prohibition(output: &str) -> Option<&'static str> {
    if output.contains(' ') {
        return Some("general.md:96 — spaces are not permitted in any HGVS description");
    }
    if output.contains('X') || output.contains('x') {
        return Some("standards.md:39 — `X` is an alignment-only symbol, not a nucleotide");
    }
    if let Some(body) = output.split_once(":g.").map(|(_, body)| body) {
        if body.contains('+') || body.contains('*') || body.contains('-') {
            return Some("checklist.md:16 — a `g.` description admits no `+`/`-`/`*` offset");
        }
    }
    None
}

/// Members of a parsed description.
#[must_use]
pub fn members_of(variant: &HgvsVariant) -> Vec<HgvsVariant> {
    match variant {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single.clone()],
    }
}

/// The coding-axis merge instrument's **numerator**: did normalization collapse
/// `row`'s authored members into fewer members?
///
/// # It is `< row.members`, not `< 2`, and the difference is 313 of 997 rows
///
/// The sibling `guard_violations` check may write `< 2` because
/// `is_svd_wg010_shape` pins `kinds.len() == 2`, so for its rows "fewer than two
/// members" and "merged" are the same statement. This population is wider —
/// `MEMBER_COUNTS = [2, 3, 4]` — and there `< 2` means "collapsed **all the way**
/// to one member". A 4-member row merged to 2 is two merges across two or more
/// unchanged nucleotides and would have been counted as **zero**. Measured over
/// the flagged rows: `{2: 684, 3: 159, 4: 154}`, so **313 of 997 (31.4%)** were
/// blind to a partial merge under the narrower test.
///
/// That is the property the counter's doc claims ("the members were merged across
/// the two or more unchanged nucleotides between them"), so measuring anything
/// narrower would be an assertion weaker than the stated contract.
///
/// It still cannot see a merge whose output is unparseable, since the count is
/// taken on a parsed output; [`Census::unparseable_outputs`] is the counter for
/// that.
///
/// Extracted as a named function rather than inlined so the gate's
/// `the_coding_axis_merge_counter_can_observe_a_merge` can exercise it directly —
/// see that test for why a zero pin needs a positive control.
#[must_use]
pub fn coding_axis_merge_observed(row: &Row, output: &HgvsVariant) -> bool {
    members_of(output).len() < row.members
}

/// Whether either endpoint of `interval` names an intronic offset.
fn interval_is_intronic<T>(
    interval: &Interval<T>,
    is_intronic: impl Fn(&T) -> bool + Copy,
) -> bool {
    let boundary = |b: &UncertainBoundary<T>| match b {
        UncertainBoundary::Single(mu) => mu.inner().is_some_and(is_intronic),
        UncertainBoundary::Range { start, end } => {
            start.inner().is_some_and(is_intronic) || end.inner().is_some_and(is_intronic)
        }
    };
    boundary(&interval.start) || boundary(&interval.end)
}

/// Whether `output` names an intronic position under a bare transcript accession —
/// the form `checklist.md:20` refuses. Read from the AST rather than from
/// `src/normalize/`, so the census judges ferro from outside.
fn names_bare_transcript_intronic(output: &HgvsVariant) -> bool {
    members_of(output).iter().any(|member| {
        let (accession, intronic) = match member {
            HgvsVariant::Cds(v) => (
                &v.accession,
                interval_is_intronic(&v.loc_edit.location, CdsPos::is_intronic),
            ),
            HgvsVariant::Tx(v) => (
                &v.accession,
                interval_is_intronic(&v.loc_edit.location, TxPos::is_intronic),
            ),
            _ => return false,
        };
        intronic && accession.genomic_context.is_none()
    })
}

// ---------------------------------------------------------------------------
// Measurement
// ---------------------------------------------------------------------------

/// The built corpus at default bounds.
#[must_use]
pub fn built() -> SpecCorpus {
    super::spec_corpus::corpus(&CorpusBounds::default())
}

/// Rows grouped so every row sharing a synthetic reference is adjacent.
#[must_use]
pub fn grouped(built: &SpecCorpus) -> Vec<&Row> {
    let mut rows: Vec<&Row> = built.rows.iter().collect();
    rows.sort_by(|a, b| {
        a.shape
            .label()
            .cmp(b.shape.label())
            .then_with(|| a.core.cmp(&b.core))
            .then_with(|| a.id.cmp(&b.id))
    });
    rows
}

/// The normalization configuration every figure is measured under. The mode
/// stamped into [`Census::measured_under`] is read off this very config.
fn measurement_config(direction: ShuffleDirection) -> NormalizeConfig {
    NormalizeConfig::default().with_direction(direction)
}

/// The contiguous runs of `rows` that share one synthetic reference.
fn reference_groups<'a>(rows: &'a [&'a Row]) -> Vec<&'a [&'a Row]> {
    let mut groups = Vec::new();
    let mut start = 0usize;
    while start < rows.len() {
        let mut end = start + 1;
        while end < rows.len()
            && rows[end].shape == rows[start].shape
            && rows[end].core == rows[start].core
        {
            end += 1;
        }
        groups.push(&rows[start..end]);
        start = end;
    }
    groups
}

/// Census one reference group. Builds its own frame and normalizer, so it is safe
/// to run alongside its siblings.
fn measure_group(
    group: &[&Row],
    direction: ShuffleDirection,
    equivalence: Equivalence,
) -> Measured {
    let mut census = Census {
        measured_under: measurement_config(direction).error_config.mode,
        ..Census::default()
    };
    let mut divergences: Vec<Divergence> = Vec::new();
    let mut findings: Vec<Finding> = Vec::new();
    let mut confluence_unkeyable_outputs = 0usize;

    let active = group[0].frame();
    let normalizer =
        Normalizer::with_config(active.provider().clone(), measurement_config(direction));
    let (active, normalizer) = (&active, &normalizer);

    for row in group {
        let mut outputs: BTreeSet<String> = BTreeSet::new();
        let mut normalized_spellings = 0usize;
        for spelling in &row.spellings {
            census.outputs += 1;
            let parsed = match parse_hgvs(spelling) {
                Ok(parsed) => parsed,
                Err(_) => {
                    if !matches!(row.kind, RowKind::Conflict | RowKind::Prohibited) {
                        census.declined += 1;
                    }
                    continue;
                }
            };
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                normalizer.normalize(&parsed)
            }));
            let output = match result {
                Ok(Ok(value)) => value.to_string(),
                _ => {
                    if !matches!(row.kind, RowKind::Conflict | RowKind::Prohibited) {
                        census.declined += 1;
                    }
                    continue;
                }
            };
            normalized_spellings += 1;

            match row.kind {
                RowKind::Conflict => {
                    census.conflicts_accepted += 1;
                    findings.push(Finding {
                        id: row.id.clone(),
                        what: format!(
                            "conflicting allele accepted, geometry {}: {spelling} -> {output}",
                            row.geometry.label()
                        ),
                    });
                }
                RowKind::Prohibited => {
                    let (clause, strength) =
                        row.prohibition.expect("a prohibited row names its clause");
                    match strength {
                        Strength::Absolute => census.prohibited_absolute_accepted += 1,
                        Strength::Conditional => census.prohibited_conditional_accepted += 1,
                    }
                    findings.push(Finding {
                        id: row.id.clone(),
                        what: format!(
                            "{} prohibition accepted ({clause}): {spelling} -> {output}",
                            strength.label()
                        ),
                    });
                }
                RowKind::Family | RowKind::Single => {}
            }

            // --- rank 1: validity ---------------------------------------------
            let reparsed = parse_hgvs(&output);
            if reparsed.is_err() {
                census.unparseable_outputs += 1;
                findings.push(Finding {
                    id: row.id.clone(),
                    what: format!("output does not re-parse: {spelling} -> {output}"),
                });
            }
            if let Some(violation) = violated_prohibition(&output) {
                census.prohibition_violating_outputs += 1;
                findings.push(Finding {
                    id: row.id.clone(),
                    what: format!("output violates {violation}: {spelling} -> {output}"),
                });
            }

            // --- sequence preservation, and the overlap half of validity ------
            if let Some(expected) = row.denoted.as_deref() {
                match denotation_of(active.provider(), active.served(), &output) {
                    Denotation::NoSequence => {
                        census.outputs_denoting_no_sequence += 1;
                        findings.push(Finding {
                            id: row.id.clone(),
                            what: format!(
                                "output denotes no sequence (members claim one territory): \
                                 {spelling} -> {output}"
                            ),
                        });
                    }
                    Denotation::Inexpressible => {
                        if reparsed
                            .as_ref()
                            .ok()
                            .is_none_or(names_bare_transcript_intronic)
                        {
                            census.outputs_leaving_the_transcript += 1;
                            findings.push(Finding {
                                id: row.id.clone(),
                                what: format!(
                                    "output left the transcript (an intronic position the input \
                                     did not name, on a BARE transcript accession): \
                                     {spelling} -> {output}"
                                ),
                            });
                        } else {
                            census.outputs_intronic_under_a_genomic_wrapper += 1;
                        }
                    }
                    Denotation::Unparseable => {}
                    Denotation::Sequence(applied) if applied != expected => {
                        census.sequence_changed += 1;
                        findings.push(Finding {
                            id: row.id.clone(),
                            what: format!("output denotes different bases: {spelling} -> {output}"),
                        });
                    }
                    Denotation::Sequence(_) => {}
                }
            }

            // --- idempotency --------------------------------------------------
            if let Ok(parsed_output) = reparsed {
                let again = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    normalizer.normalize(&parsed_output)
                }));
                let stable = matches!(&again, Ok(Ok(value)) if value.to_string() == output);
                if !stable {
                    census.non_idempotent_outputs += 1;
                    findings.push(Finding {
                        id: row.id.clone(),
                        what: format!(
                            "output is not a fixed point: {output} -> {}",
                            match &again {
                                Ok(Ok(value)) => value.to_string(),
                                Ok(Err(error)) => format!("<declined: {error}>"),
                                Err(_) => "<panicked>".to_string(),
                            }
                        ),
                    });
                }
            }

            // --- negative guards ---------------------------------------------
            if !row.negative_guards.is_empty() && spelling == row.authored_spelling() {
                if let Ok(parsed_output) = parse_hgvs(&output) {
                    if members_of(&parsed_output).len() < 2 {
                        census.guard_violations += 1;
                        findings.push(Finding {
                            id: row.id.clone(),
                            what: format!(
                                "merged a frameless separation of one into a single member, \
                                 which is rejected SVD-WG010: {spelling} -> {output}"
                            ),
                        });
                    }
                }
            }

            // --- the coding-axis merge instrument -----------------------------
            if row.coding_axis_separation_two_or_more && spelling == row.authored_spelling() {
                census.coding_axis_separation_two_or_more_rows += 1;
                if let Ok(parsed_output) = parse_hgvs(&output) {
                    if coding_axis_merge_observed(row, &parsed_output) {
                        census.coding_axis_separation_two_or_more_merges += 1;
                        findings.push(Finding {
                            id: row.id.clone(),
                            what: format!(
                                "coding-axis allele merged across {} unchanged nucleotide(s), \
                                 {} authored members -> {} (counted, not adjudicated): \
                                 {spelling} -> {output}",
                                row.separation,
                                row.members,
                                members_of(&parsed_output).len()
                            ),
                        });
                    }
                }
            }

            outputs.insert(output);
        }

        // --- rank 2: confluence, over families only --------------------------
        if row.kind == RowKind::Family {
            if normalized_spellings < 2 {
                census.underdetermined += 1;
                continue;
            }
            let (keyed, unkeyable): (BTreeSet<String>, usize) =
                bucket_outputs(&outputs, equivalence, active.provider());
            confluence_unkeyable_outputs += unkeyable;
            if unkeyable > 0 || keyed.is_empty() {
                // A family the relation could not fully key cannot be certified
                // converged; it is undecided, and the unkeyable outputs are
                // counted above so this does not read as a clean decision.
                census.underdetermined += 1;
                continue;
            }
            match keyed.len() {
                1 => census.converged += 1,
                2 => census.split_two += 1,
                3 => census.split_three += 1,
                _ => census.split_more += 1,
            }
            if keyed.len() > 1 {
                divergences.push(Divergence {
                    id: row.id.clone(),
                    outputs: outputs.iter().cloned().collect(),
                });
            }
        }
    }

    Measured {
        census,
        equivalence,
        divergences,
        findings,
        confluence_unkeyable_outputs,
    }
}

/// Group a family's distinct outputs by the relation, returning the set of
/// distinct group keys and the count of outputs the relation could not key.
fn bucket_outputs(
    outputs: &BTreeSet<String>,
    equivalence: Equivalence,
    provider: &MockProvider,
) -> (BTreeSet<String>, usize) {
    let mut keyed = BTreeSet::new();
    let mut unkeyable = 0usize;
    for output in outputs {
        match equivalence.key_of(output, provider) {
            Some(key) => {
                keyed.insert(key);
            }
            None => unkeyable += 1,
        }
    }
    (keyed, unkeyable)
}

/// Census one direction under one relation, one reference group per rayon task.
///
/// **The result is byte-identical to a serial walk, not merely equivalent**, which
/// is the only basis on which a pinned census may be parallelized at all:
///
/// * every counter is a count over a partition of the rows, so [`Census::absorb`]
///   over the groups is exact and order-free;
/// * `par_iter().collect()` preserves *input* order, so the per-group results fold
///   in corpus order however they finish — which matters for `findings` and
///   `divergences`, uncapped lists whose order the summary reproduces;
/// * both of those lists are appended to in row order within a group, and the
///   groups are contiguous runs of the sorted rows, so the concatenation is
///   exactly the order a serial walk would produce;
/// * nothing is shared across tasks: [`measure_group`] builds its own `Frame` and
///   `Normalizer` (a serial walk rebuilds both per reference too — that is what the
///   grouping is), and the only borrow that crosses is the immutable corpus slice.
#[must_use]
pub fn measure(direction: ShuffleDirection, equivalence: Equivalence) -> Measured {
    measure_corpus(&built(), direction, equivalence)
}

/// Census a specific built corpus, so a caller that also wants the corpus's
/// completeness (via [`CensusReport::build`]) enumerates it exactly once.
#[must_use]
pub fn measure_corpus(
    built: &SpecCorpus,
    direction: ShuffleDirection,
    equivalence: Equivalence,
) -> Measured {
    let rows = grouped(built);
    let groups = reference_groups(&rows);

    let per_group: Vec<Measured> = groups
        .par_iter()
        .map(|group| measure_group(group, direction, equivalence))
        .collect();

    let mut measured = Measured {
        census: Census {
            measured_under: measurement_config(direction).error_config.mode,
            ..Census::default()
        },
        equivalence,
        divergences: Vec::new(),
        findings: Vec::new(),
        confluence_unkeyable_outputs: 0,
    };
    for group in per_group {
        measured.census.absorb(&group.census);
        measured.divergences.extend(group.divergences);
        measured.findings.extend(group.findings);
        measured.confluence_unkeyable_outputs += group.confluence_unkeyable_outputs;
    }
    measured
}

/// A human-readable summary of a run. Never prints a single percentage or a
/// verdict — rank-1 and rank-2 are not commensurable.
#[must_use]
pub fn report(label: &str, measured: &Measured) -> String {
    let census = &measured.census;
    let mut out = format!(
        "spec conformance census ({label}, confluence under `{}`)\n  \
         NORMALIZED UNDER: error mode `{}` — every figure downstream of normalization is a \
         figure for THAT mode; the parse half applies no mode at all (#1632)\n  \
         VALIDITY (rank 1): {} outputs, {} declined, {} unparseable, {} denoting no sequence, \
         {} leaving the transcript on a bare accession, {} intronic under a genomic wrapper \
         (conformant), {} violating an absolute prohibition\n  \
         CONFLUENCE (rank 2): converged {}, split 2 {}, split 3 {}, split 4+ {}, \
         underdetermined {} ({} outputs unkeyable under this relation)\n  \
         IDEMPOTENCY: {} outputs are not their own fixed point\n  \
         SEQUENCE PRESERVATION: {} outputs denote different bases\n  \
         REFUSAL: {} conflicting alleles accepted, {} absolute prohibitions accepted, \
         {} conditional prohibitions accepted\n  \
         NEGATIVE GUARDS: {} outputs implement rejected SVD-WG010\n  \
         CODING-AXIS MERGE (an instrument, not a verdict): {} of {} coding-axis alleles \
         separated by two or more unchanged nucleotides normalized to fewer members than \
         were authored (a full or a partial merge)\n",
        measured.equivalence.label(),
        census.measured_under,
        census.outputs,
        census.declined,
        census.unparseable_outputs,
        census.outputs_denoting_no_sequence,
        census.outputs_leaving_the_transcript,
        census.outputs_intronic_under_a_genomic_wrapper,
        census.prohibition_violating_outputs,
        census.converged,
        census.split_two,
        census.split_three,
        census.split_more,
        census.underdetermined,
        measured.confluence_unkeyable_outputs,
        census.non_idempotent_outputs,
        census.sequence_changed,
        census.conflicts_accepted,
        census.prohibited_absolute_accepted,
        census.prohibited_conditional_accepted,
        census.guard_violations,
        census.coding_axis_separation_two_or_more_merges,
        census.coding_axis_separation_two_or_more_rows,
    );
    for divergence in measured.divergences.iter().take(MAX_DIVERGENCES) {
        out.push_str(&format!(
            "  DIVERGES {} -> {:?}\n",
            divergence.id, divergence.outputs
        ));
    }
    // The summary is capped for readability, so it must say so — a reader who
    // cannot tell a complete list from a truncated one is the failure this whole
    // paragraph of `Measured::divergences` is about, and it applies to prose too.
    let hidden = measured.divergences.len().saturating_sub(MAX_DIVERGENCES);
    if hidden > 0 {
        out.push_str(&format!(
            "  DIVERGES … and {hidden} more (this summary names the first \
             {MAX_DIVERGENCES} in corpus order; `--findings` carries all {} of them)\n",
            measured.divergences.len()
        ));
    }
    let mut grouped: BTreeMap<String, (usize, String)> = BTreeMap::new();
    for finding in &measured.findings {
        let key = finding
            .what
            .split_once(':')
            .map_or_else(|| finding.what.clone(), |(head, _)| head.to_string());
        let entry = grouped
            .entry(key)
            .or_insert_with(|| (0, format!("{} | {}", finding.id, finding.what)));
        entry.0 += 1;
    }
    for (kind, (count, example)) in grouped {
        out.push_str(&format!("  FINDING x{count} [{kind}] e.g. {example}\n"));
    }
    out
}

// ---------------------------------------------------------------------------
// The runnable artifact
// ---------------------------------------------------------------------------

/// Which corpus a run measures.
///
/// Only the hermetic spec corpus is supported today. A real-expression corpus
/// against a prepared reference is the consumer-side self-check's job (#1892),
/// which builds on this instrument.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorpusSource {
    /// The synthetic spec-derived corpus, enumerated hermetically at run time.
    Spec,
}

impl CorpusSource {
    /// Stable label, stamped into the artifact.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Spec => "spec",
        }
    }
}

/// A census run as a machine-readable artifact — the thing `--out` writes and
/// `--compare` reads.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CensusReport {
    /// The corpus measured (`spec`).
    pub corpus: String,
    /// The shuffle direction (`3prime`).
    pub direction: String,
    /// The relation confluence was decided under (`sequence`/`partition`/`spdi`).
    pub equivalence: String,
    /// The corpus's shape, so a reported zero always carries its denominator.
    pub shape: CorpusShape,
    /// The corpus's completeness: `attempted == succeeded + dropped`, with drops
    /// by reason. A partial run cannot read as complete.
    pub capture: CaptureCounts,
    /// The census.
    pub census: Census,
    /// Outputs the relation could not key, dropped from confluence and counted so
    /// a partial confluence decision cannot read as a clean one.
    pub confluence_unkeyable_outputs: usize,
}

/// Why a census run could not produce an artifact.
#[derive(Debug)]
pub enum CensusError {
    /// A `FERRO_ASSERT_*` seam oracle is armed, which would flatter the census.
    /// Constructed by [`run_census`] itself, so every entry point inherits the
    /// refusal — see that function for why it is not the caller's.
    OracleArmed(String),
    /// The corpus dropped more designs than its allowance permits — see
    /// [`max_dropped_designs`] for the ceiling and why it is a real one.
    Incomplete(Shortfall),
}

impl std::fmt::Display for CensusError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::OracleArmed(var) => write!(
                f,
                "refusing to census with a seam oracle armed ({var}): a panicking row contributes \
                 no output, so its family silently converges and the census would read better than \
                 the truth. Unset {var} and re-run."
            ),
            Self::Incomplete(shortfall) => {
                write!(f, "the corpus is incomplete: {shortfall}")
            }
        }
    }
}

impl std::error::Error for CensusError {}

/// A census run: the artifact plus the row-level detail (`findings`,
/// `divergences`) that `--findings` emits and `--out` does not.
#[derive(Debug, Clone)]
pub struct CensusRun {
    /// The machine-readable artifact.
    pub report: CensusReport,
    /// The full measurement, including per-row findings and divergences.
    pub measured: Measured,
}

impl CensusReport {
    /// Run the census and assemble the artifact.
    ///
    /// A thin wrapper over [`run_census`] for callers that only want the artifact.
    pub fn build(
        corpus: CorpusSource,
        direction: ShuffleDirection,
        equivalence: Equivalence,
    ) -> Result<Self, CensusError> {
        Ok(run_census(corpus, direction, equivalence)?.report)
    }
}

/// Run the census, returning both the artifact and the row-level detail.
///
/// Refuses first if a seam oracle is armed, then enumerates the corpus once,
/// routes every design through a [`CaptureLedger`] so drops are accounted, then
/// measures.
///
/// # The oracle check is here, not at the caller
///
/// It used to be the CLI's, and this function's docs said so ("does not read the
/// environment"). That left [`CensusError::OracleArmed`] unconstructible and the
/// protection un-inherited: an in-process consumer (#1891's stability harness,
/// #1892's self-check) calling `run_census` got **no** oracle protection at all
/// while the error type it matched on advertised one. A refusal that only one of
/// two entry points performs is not a property of the census.
///
/// [`oracle_armed`] stays pure over an iterator of names so the prefix rule can be
/// tested without mutating the process environment (which races other tests);
/// [`oracle_armed_in_env`] is the one-line application of it that this function
/// calls, and `a_seam_oracle_in_the_environment_is_refused` exercises *that* path
/// through a subprocess with a controlled environment.
pub fn run_census(
    corpus: CorpusSource,
    direction: ShuffleDirection,
    equivalence: Equivalence,
) -> Result<CensusRun, CensusError> {
    if let Some(var) = oracle_armed_in_env() {
        return Err(CensusError::OracleArmed(var));
    }
    let CorpusSource::Spec = corpus;
    let bounds = CorpusBounds::default();

    // Route the enumeration through the ledger — the last fallible step before the
    // artifact — so a design that produced no row is accounted, not silently
    // absent.
    let attempts = enumerate(&bounds);
    let designs_considered = attempts.len();
    let mut ledger = CaptureLedger::new("spec conformance corpus designs");
    let mut drops: BTreeMap<&'static str, usize> = BTreeMap::new();
    let mut rows = Vec::new();
    for attempt in attempts {
        match attempt {
            Ok(row) => {
                ledger.record_success();
                rows.push(row);
            }
            Err((id, reason)) => {
                ledger.record_drop(id, reason.label());
                *drops.entry(reason.label()).or_default() += 1;
            }
        }
    }
    rows.sort_by(|a, b| a.id.cmp(&b.id));
    let built = SpecCorpus {
        bounds,
        designs_considered,
        drops,
        rows,
    };

    // The corpus drops designs by construction (out of range, denotes the
    // reference, fewer than two spellings). That is enumerated shapes with no valid
    // row, not lost data, so the shortfall is allowed — but named, stamped into the
    // artifact, and BOUNDED.
    let capture = ledger
        .finish_with(Allowance::at_most(
            max_dropped_designs(designs_considered),
            "the spec corpus deliberately drops designs that run out of range, denote the \
             reference, or leave fewer than two spellings; these are enumerated shapes with no \
             valid row. The ceiling is half the enumeration — see `max_dropped_designs`",
        ))
        .map_err(CensusError::Incomplete)?;

    let measured = measure_corpus(&built, direction, equivalence);
    let report = CensusReport {
        corpus: corpus.label().to_string(),
        direction: direction_label(direction).to_string(),
        equivalence: equivalence.label().to_string(),
        shape: CorpusShape::of(&built),
        capture,
        census: measured.census.clone(),
        confluence_unkeyable_outputs: measured.confluence_unkeyable_outputs,
    };
    Ok(CensusRun { report, measured })
}

/// The most designs the enumeration may drop before the census is refused as
/// incomplete: **half** of what it considered.
///
/// # Why a ceiling at all, and why this number
///
/// The first revision passed `designs_considered` itself. With
/// `empty_pass_permitted: false` that refuses exactly one thing — a run that
/// attempted *nothing* — so the CLI's own claim that "a corpus that dropped more
/// than its allowance permits is a refusal" was true and unreachable: every drop
/// count from 0 to 100% passed. An allowance equal to the population is not a
/// bound, it is the shape of one.
///
/// Measured on the corpus this ships against: **3704 of 16650 designs dropped
/// (22.2%)**. Half is therefore ~2.25x the observed rate — loose enough that
/// ordinary generator growth (a new shape family whose members mostly run out of
/// range) does not turn a healthy corpus into a refusal, and tight enough to catch
/// the failure that matters, which is a corpus that mostly failed to *build*
/// reading as a complete census of what it did build. It also gives the census a
/// floor it did not have: with drops capped at half, `succeeded` is at least half
/// the enumeration, so no run can report confluence over a handful of surviving
/// rows.
///
/// It is deliberately a **fraction of the enumeration** rather than a literal: the
/// literal would be a pinned number nothing keeps in step with the corpus, which
/// is the change-detector-not-a-guard failure `CONTRIBUTING.md` names.
fn max_dropped_designs(designs_considered: usize) -> usize {
    designs_considered / 2
}

/// The label a direction is stamped with in the artifact.
fn direction_label(direction: ShuffleDirection) -> &'static str {
    match direction {
        ShuffleDirection::ThreePrime => "3prime",
        ShuffleDirection::FivePrime => "5prime",
    }
}

// ---------------------------------------------------------------------------
// The seam-oracle refusal
// ---------------------------------------------------------------------------

/// The name of the (lexicographically first) armed seam oracle among
/// `env_keys`, if any.
///
/// A `FERRO_ASSERT_*` variable arms an oracle that panics on the very rows the
/// census is trying to count; a panicking row is filed `declined` and never
/// reaches its family's output set, so the census reads better than the truth.
/// The caller refuses rather than emit a flattering artifact — the same reasoning
/// that put `ORACLE_EXCLUDE` in `ci.yml`.
///
/// Pure over an iterator of variable names, so it can be tested without mutating
/// the process environment (which races other tests). [`oracle_armed_in_env`]
/// applies it to the live environment.
pub fn oracle_armed<I>(env_keys: I) -> Option<String>
where
    I: IntoIterator<Item = String>,
{
    env_keys
        .into_iter()
        .filter(|key| key.starts_with("FERRO_ASSERT_"))
        .min()
}

/// [`oracle_armed`] over the live process environment.
#[must_use]
pub fn oracle_armed_in_env() -> Option<String> {
    oracle_armed(std::env::vars().map(|(key, _)| key))
}

// ---------------------------------------------------------------------------
// The burn-down comparison
// ---------------------------------------------------------------------------

/// Which way a counter should move to be virtuous.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DirectionOfVirtue {
    /// A failure counter: falling is an improvement.
    LowerIsBetter,
    /// `converged`: rising is an improvement.
    HigherIsBetter,
    /// A denominator or an instrument: a movement is a population change, not a
    /// verdict.
    Neutral,
}

/// Whether a movement was in the virtuous direction.
///
/// The four variants are four distinct facts, and collapsing the last two is how a
/// burn-down loses a movement: a neutral counter that moved used to be reported as
/// `Unchanged`, which is a statement about the *value* and the value changed. A
/// denominator moving means the corpus changed under the comparison, which is
/// precisely what a reader must see.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Movement {
    /// Moved the virtuous way.
    Improved,
    /// Moved against virtue.
    Worse,
    /// Moved, on a counter that has no direction of virtue — a population change,
    /// reported without a grade.
    Moved,
    /// Did not move at all.
    Unchanged,
}

/// One counter's movement between two runs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CounterMovement {
    /// The census field name.
    pub name: &'static str,
    /// Its value in the base run.
    pub before: usize,
    /// Its value in the candidate run.
    pub after: usize,
    /// Which way it should move to be virtuous.
    pub virtue: DirectionOfVirtue,
    /// Which way it actually moved.
    pub movement: Movement,
}

/// Per-counter movement between two census reports, with the direction of virtue
/// attached. Reports movements only — it never collapses to a single verdict or a
/// percentage, because rank-1 and rank-2 are not commensurable.
///
/// # Every field of both reports is accounted for, and the compiler enforces it
///
/// The first revision built the rows from `before.census.<field>` accesses under a
/// length-annotated array literal (`[…; 20]`) and a comment calling the length "the
/// reminder" that a new field needs a row. **It is not one.** An array's length has
/// nothing to do with a struct's field count, so adding a field to [`Census`]
/// compiles clean and the burn-down silently omits it — which is #1710 exactly, one
/// layer up: two fields added, never summed, and the census reported `0 of 0`.
/// [`Census::absorb`] had already been rewritten to destructure for that reason.
///
/// So both reports are destructured here with **no `..`**, and every row below is
/// built from a binding. Adding a field to [`Census`] or to [`CensusReport`] is a
/// compile error until it is either given a row or explicitly bound to `_` with a
/// reason — which is the point: skipping a counter must be a decision someone
/// wrote down, not a thing that happens.
///
/// # Relation mismatch neutralizes the confluence counters
///
/// Confluence is decided *under a relation*. Comparing a `sequence` run against a
/// `partition` run therefore changes the question, not the behaviour, and grading
/// `converged 11750 -> 11071` as WORSE would be a false verdict on a comparison
/// that has none to give. The caller already prints a NOTE saying so; printing the
/// note and then the verdict is worse than either alone, because the verdict is
/// what a reader carries away. Under mismatched relations the five confluence
/// buckets and `confluence_unkeyable_outputs` are reported [`Neutral`], so they
/// read as `moved` — the movement is still visible, and nothing is graded.
///
/// [`Neutral`]: DirectionOfVirtue::Neutral
#[must_use]
#[allow(clippy::too_many_lines)] // one row per counter; splitting it hides the exhaustiveness
pub fn compare_reports(before: &CensusReport, after: &CensusReport) -> Vec<CounterMovement> {
    use DirectionOfVirtue::{HigherIsBetter, LowerIsBetter, Neutral};

    // `corpus` and `direction` label the run rather than count anything, and
    // `shape` is the corpus's own shape — pinned separately by the gate, and its
    // one figure the census is `n of` (`spellings`) is mirrored by `outputs`
    // below, which is in the burn-down. `equivalence` is read, not compared: it
    // decides whether the confluence rows may be graded at all.
    let CensusReport {
        corpus: _,
        direction: _,
        equivalence: before_relation,
        shape: _,
        capture: before_capture,
        census: before_census,
        confluence_unkeyable_outputs: before_unkeyable,
    } = before;
    let CensusReport {
        corpus: _,
        direction: _,
        equivalence: after_relation,
        shape: _,
        capture: after_capture,
        census: after_census,
        confluence_unkeyable_outputs: after_unkeyable,
    } = after;

    // `measured_under` is a mode stamp rather than a counter, and a census
    // compared against one taken under a different mode is a category error the
    // reader must resolve before reading any row — not a movement.
    let Census {
        measured_under: _,
        outputs: b_outputs,
        declined: b_declined,
        unparseable_outputs: b_unparseable,
        outputs_denoting_no_sequence: b_no_sequence,
        outputs_leaving_the_transcript: b_leaving,
        outputs_intronic_under_a_genomic_wrapper: b_intronic_wrapped,
        prohibition_violating_outputs: b_prohibition_violating,
        converged: b_converged,
        split_two: b_split_two,
        split_three: b_split_three,
        split_more: b_split_more,
        underdetermined: b_underdetermined,
        non_idempotent_outputs: b_non_idempotent,
        sequence_changed: b_sequence_changed,
        conflicts_accepted: b_conflicts,
        prohibited_absolute_accepted: b_absolute,
        prohibited_conditional_accepted: b_conditional,
        guard_violations: b_guards,
        coding_axis_separation_two_or_more_rows: b_coding_rows,
        coding_axis_separation_two_or_more_merges: b_coding_merges,
    } = before_census;
    let Census {
        measured_under: _,
        outputs: a_outputs,
        declined: a_declined,
        unparseable_outputs: a_unparseable,
        outputs_denoting_no_sequence: a_no_sequence,
        outputs_leaving_the_transcript: a_leaving,
        outputs_intronic_under_a_genomic_wrapper: a_intronic_wrapped,
        prohibition_violating_outputs: a_prohibition_violating,
        converged: a_converged,
        split_two: a_split_two,
        split_three: a_split_three,
        split_more: a_split_more,
        underdetermined: a_underdetermined,
        non_idempotent_outputs: a_non_idempotent,
        sequence_changed: a_sequence_changed,
        conflicts_accepted: a_conflicts,
        prohibited_absolute_accepted: a_absolute,
        prohibited_conditional_accepted: a_conditional,
        guard_violations: a_guards,
        coding_axis_separation_two_or_more_rows: a_coding_rows,
        coding_axis_separation_two_or_more_merges: a_coding_merges,
    } = after_census;

    // Confluence is a question asked under a relation; two relations are two
    // questions, and a question change has no direction of virtue.
    let same_relation = before_relation == after_relation;
    let confluence = |virtue: DirectionOfVirtue| {
        if same_relation {
            virtue
        } else {
            Neutral
        }
    };

    let counters: [(&'static str, usize, usize, DirectionOfVirtue); 24] = [
        ("outputs", *b_outputs, *a_outputs, Neutral),
        ("declined", *b_declined, *a_declined, LowerIsBetter),
        (
            "unparseable_outputs",
            *b_unparseable,
            *a_unparseable,
            LowerIsBetter,
        ),
        (
            "outputs_denoting_no_sequence",
            *b_no_sequence,
            *a_no_sequence,
            LowerIsBetter,
        ),
        (
            "outputs_leaving_the_transcript",
            *b_leaving,
            *a_leaving,
            LowerIsBetter,
        ),
        (
            "prohibition_violating_outputs",
            *b_prohibition_violating,
            *a_prohibition_violating,
            LowerIsBetter,
        ),
        (
            "non_idempotent_outputs",
            *b_non_idempotent,
            *a_non_idempotent,
            LowerIsBetter,
        ),
        (
            "sequence_changed",
            *b_sequence_changed,
            *a_sequence_changed,
            LowerIsBetter,
        ),
        (
            "conflicts_accepted",
            *b_conflicts,
            *a_conflicts,
            LowerIsBetter,
        ),
        (
            "prohibited_absolute_accepted",
            *b_absolute,
            *a_absolute,
            LowerIsBetter,
        ),
        (
            "prohibited_conditional_accepted",
            *b_conditional,
            *a_conditional,
            LowerIsBetter,
        ),
        ("guard_violations", *b_guards, *a_guards, LowerIsBetter),
        (
            "split_two",
            *b_split_two,
            *a_split_two,
            confluence(LowerIsBetter),
        ),
        (
            "split_three",
            *b_split_three,
            *a_split_three,
            confluence(LowerIsBetter),
        ),
        (
            "split_more",
            *b_split_more,
            *a_split_more,
            confluence(LowerIsBetter),
        ),
        (
            "underdetermined",
            *b_underdetermined,
            *a_underdetermined,
            confluence(LowerIsBetter),
        ),
        (
            "converged",
            *b_converged,
            *a_converged,
            confluence(HigherIsBetter),
        ),
        // Outputs the relation could not key are dropped from the confluence
        // decision, so a rise makes the decision *less* complete — the same
        // partial-run-reading-as-complete hazard `underdetermined` guards, and it
        // is only visible there indirectly. It is reported in its own right for
        // that reason, and is relation-dependent, so it neutralizes with the rest.
        (
            "confluence_unkeyable_outputs",
            *before_unkeyable,
            *after_unkeyable,
            confluence(LowerIsBetter),
        ),
        (
            "outputs_intronic_under_a_genomic_wrapper",
            *b_intronic_wrapped,
            *a_intronic_wrapped,
            Neutral,
        ),
        (
            "coding_axis_separation_two_or_more_rows",
            *b_coding_rows,
            *a_coding_rows,
            Neutral,
        ),
        (
            "coding_axis_separation_two_or_more_merges",
            *b_coding_merges,
            *a_coding_merges,
            Neutral,
        ),
        // The corpus's own completeness. Populations, not verdicts: a drop is not a
        // defect (the corpus drops by construction), but a run whose enumeration
        // shrank, or whose drop count rose, is not comparable to one whose did not
        // — and without these rows the only way to see that was to open both JSON
        // artifacts by hand.
        (
            "capture_attempted",
            before_capture.attempted,
            after_capture.attempted,
            Neutral,
        ),
        (
            "capture_succeeded",
            before_capture.succeeded,
            after_capture.succeeded,
            Neutral,
        ),
        (
            "capture_dropped",
            before_capture.dropped,
            after_capture.dropped,
            Neutral,
        ),
    ];
    counters
        .into_iter()
        .map(|(name, before, after, virtue)| CounterMovement {
            name,
            before,
            after,
            virtue,
            movement: classify_movement(before, after, virtue),
        })
        .collect()
}

/// Classify a movement given the direction of virtue.
///
/// `pub(crate)` so the CLI's rendering test can go from `(before, after, virtue)`
/// to a rendered line rather than hand-constructing the [`Movement`] it means to
/// be testing the rendering of — a test that supplies its own answer checks only
/// that a `match` arm exists.
pub(crate) fn classify_movement(
    before: usize,
    after: usize,
    virtue: DirectionOfVirtue,
) -> Movement {
    use std::cmp::Ordering::{Equal, Greater, Less};
    match (virtue, after.cmp(&before)) {
        (_, Equal) => Movement::Unchanged,
        (DirectionOfVirtue::Neutral, _) => Movement::Moved,
        (DirectionOfVirtue::LowerIsBetter, Less) | (DirectionOfVirtue::HigherIsBetter, Greater) => {
            Movement::Improved
        }
        (DirectionOfVirtue::LowerIsBetter, Greater) | (DirectionOfVirtue::HigherIsBetter, Less) => {
            Movement::Worse
        }
    }
}

// ---------------------------------------------------------------------------
// serde helper for ErrorMode (which does not derive Serialize itself)
// ---------------------------------------------------------------------------

mod error_mode_serde {
    use super::ErrorMode;
    use serde::{Deserialize, Deserializer, Serializer};

    pub fn serialize<S: Serializer>(mode: &ErrorMode, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&mode.to_string())
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(deserializer: D) -> Result<ErrorMode, D::Error> {
        let raw = String::deserialize(deserializer)?;
        match raw.as_str() {
            "strict" => Ok(ErrorMode::Strict),
            "lenient" => Ok(ErrorMode::Lenient),
            "silent" => Ok(ErrorMode::Silent),
            other => Err(serde::de::Error::custom(format!(
                "unknown error mode `{other}`"
            ))),
        }
    }
}
