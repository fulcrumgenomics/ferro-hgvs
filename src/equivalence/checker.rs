//! Equivalence checker implementation.

use crate::convert::CoordinateMapper;
use crate::error::FerroError;
use crate::hgvs::edit::{InsertedPart, InsertedSequence, NaEdit, RepeatCount};
use crate::hgvs::interval::{GenomeInterval, UncertainBoundary};
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{AllelePhase, HgvsVariant};
use crate::normalize::{NormalizeConfig, Normalizer};
use crate::reference::transcript::Strand;
use crate::reference::ReferenceProvider;
use crate::sequence::reverse_complement;
use crate::spdi::{hgvs_to_spdi, SpdiVariant};

/// Largest reference window (in bases) a sequence-level equivalence rung will
/// reconstruct, bounding the cost of a reference fetch.
///
/// **Declining past this leaves [`EquivalenceLevel::Indeterminate`], not
/// `NotEquivalent`.** This comment used to argue the opposite — that two
/// variants whose edited spans are farther apart than the cap "cannot describe
/// the same local edit", so the decline "only ever leaves the pre-existing
/// `NotEquivalent` verdict in place". That reasoning holds for two *single*
/// edits and fails for the shape this crate constructs on purpose: a cis allele
/// whose members are 150 kb apart has a union window past the cap, and two
/// spellings of it really can be one variant. Reporting a decided negative
/// there is exactly the conflation [`EquivalenceLevel::Indeterminate`] exists
/// to remove, so the wide case is now reported as undecided rather than argued
/// away.
const MAX_SEQUENCE_COMPARE_WINDOW: u64 = 100_000;

/// Level of equivalence between two variants.
///
/// This enum is open-ended: the checker gains rungs as new classes of
/// equivalence become decidable (`SequenceMatch` was added in #1158 for two
/// encodings that produce the same edited sequence). Marked `#[non_exhaustive]`
/// so downstream callers must include a wildcard arm when matching — adding a
/// rung is therefore not a breaking change. Mirrors the same attribute on
/// [`EquivalenceResult`] and on the error/projection API (#1033).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum EquivalenceLevel {
    /// Exactly identical (same string representation).
    Identical,
    /// Equivalent after normalization (e.g., shifted to same position).
    NormalizedMatch,
    /// Equivalent by apply-equality on **every determined axis**: the strongest
    /// denotational rung, and the one that establishes variant identity.
    ///
    /// An axis is *determined* when the description's coordinates can be carried
    /// to it by a mapping that is a function of **reference data alone** — exon
    /// alignments and CDS offsets — and never of normalization. A `c.`/`n.`/`r.`
    /// description on a transcript determines two: the transcript itself, and
    /// the genome its exon alignment carries it to. A `g.`/`m.` description
    /// determines only the genome (a genomic description does not determine a
    /// transcript — many overlap it). Protein is excluded: translation is
    /// many-to-one, and `p.` states a consequence, not a denotation.
    ///
    /// This is the rung a confluence gate should require; see
    /// [`EquivalenceLevel::is_at_least`].
    CrossAxisSequenceMatch,
    /// Equivalent by resulting reference sequence **on the description's own
    /// axis**: the two variants normalize to different HGVS strings but, when
    /// applied to that one reference, produce the same edited sequence (e.g. a
    /// length-changing `delins` vs a decomposed cis allele of the same edit).
    /// See issue #1158.
    ///
    /// **True, and insufficient for identity.** Single-axis apply-equality does
    /// not establish that two descriptions denote the same variant, because a
    /// transcript description also determines a genomic denotation and the two
    /// can disagree. The spec's own worked counterexample:
    /// `LRG_199t1:c.3921dup` and `LRG_199t1:c.3922dup` denote the *same
    /// transcript sequence* — both positions carry `T` — yet project to
    /// `NC_000023.11:g.32441180dup` and `g.32438390dup`, ~2,790 bp apart in
    /// different exons (exon 27 ends at 32441180, exon 28 begins at 32438390 on
    /// the minus-strand `NM_004006.2`). `DNA/duplication.md:148` calls the
    /// second "the wrong nucleotide, in the wrong exon".
    ///
    /// So this rung reports a real relationship and must **not** be used as a
    /// confluence gate. Require [`Self::CrossAxisSequenceMatch`] instead.
    SequenceMatch,
    /// Same variant but different accession versions (e.g., NM_000088.3 vs NM_000088.4).
    ///
    /// **Deliberately off the denotational order** (see
    /// [`EquivalenceLevel::is_at_least`]): it is a claim about two *different*
    /// reference sequences, and the relation's first conjunct is "same
    /// accession". Apply-equality is not defined between two versions — the
    /// bases at a coordinate may themselves differ — so this rung is orthogonal
    /// to the sequence rungs rather than stronger or weaker than them. It stays
    /// a positive, decided verdict; it simply never satisfies a gate.
    AccessionVersionDifference,
    /// Not equivalent - represent different changes.
    ///
    /// A **positive claim**, sound only when some rung above actually examined
    /// the two descriptions. Where nothing could examine them, the verdict is
    /// [`Self::Indeterminate`].
    NotEquivalent,
    /// Neither equivalent nor distinguishable: at least one side has no
    /// computable denotation on some determined axis.
    ///
    /// **Not a negative verdict.** Before this rung, `NotEquivalent` did double
    /// duty — "these denote different sequences" *and* "I could not compute a
    /// denotation" — and [`Self::is_equivalent`] answered `false` to both, so an
    /// undecidable pair was reported as a decided negative. Use
    /// [`Self::is_decided`] to tell the two apart; `is_equivalent` stays `false`
    /// here, because undecidable is not a positive either.
    ///
    /// What lands here is a **closed list**, deliberately, and it is a floor
    /// that grows as shapes are measured rather than a claim of completeness:
    ///
    /// * a non-literal inserted payload (`ins6`, `insN[33]`, `delins(?)`) —
    ///   the description denotes no sequence for `apply` to produce;
    /// * a null or unknown allele (`0`, `?`), and a multi-molecule allele
    ///   (trans / mosaic / chimeric / and-or / products / unknown phase), which
    ///   denotes no *single* sequence;
    /// * a span wider than [`MAX_SEQUENCE_COMPARE_WINDOW`] **on the
    ///   description's own axis**, where the checker declined to reconstruct
    ///   the reference at all. The same cap struck on a *second* determined
    ///   axis does **not** land here: the own-axis comparison ran and agreed,
    ///   so the verdict is [`Self::SequenceMatch`] — see the note on that arm
    ///   in `strengthen_across_axes`.
    ///
    /// **Not** here: two descriptions on different accessions. That fails the
    /// relation's first conjunct outright, which is a decided negative.
    /// **Not yet** here: an uncertain position range (`g.(100_150)del`), which
    /// would need a location accessor generic over every variant kind that this
    /// crate does not currently expose. Add it when a wrong answer is measured,
    /// with the pair that measures it.
    Indeterminate,
}

impl EquivalenceLevel {
    /// Returns true if the variants are considered equivalent.
    ///
    /// [`Self::Indeterminate`] answers `false` — undecidable is not a positive
    /// verdict — so a caller that needs to separate "no" from "cannot tell"
    /// must also consult [`Self::is_decided`].
    pub fn is_equivalent(&self) -> bool {
        !matches!(
            self,
            EquivalenceLevel::NotEquivalent | EquivalenceLevel::Indeterminate
        )
    }

    /// Returns true unless the checker could not decide — i.e. false **only**
    /// for [`Self::Indeterminate`].
    ///
    /// This is the predicate a census wants: a confluence measurement asserts
    /// over the *decided* pairs and reports the undecidable count alongside,
    /// rather than folding it into either side.
    pub fn is_decided(&self) -> bool {
        !matches!(self, EquivalenceLevel::Indeterminate)
    }

    /// Whether this verdict is at least as strong as `floor` on the
    /// **denotational order**.
    ///
    /// The order is
    /// [`Identical`](Self::Identical) ⊐
    /// [`CrossAxisSequenceMatch`](Self::CrossAxisSequenceMatch) ⊐
    /// [`SequenceMatch`](Self::SequenceMatch), and nothing else is on it. A
    /// confluence gate is written `level.is_at_least(CrossAxisSequenceMatch)`.
    ///
    /// **[`NormalizedMatch`](Self::NormalizedMatch) is off the order in both
    /// directions, and that is the point.** That rung consults the normalizer,
    /// so gating on it would define the equivalence relation in terms of the
    /// very function the gate is about — the circularity this whole ladder
    /// exists to remove. Making it unusable as a floor means the circular gate
    /// cannot be written by accident.
    /// [`AccessionVersionDifference`](Self::AccessionVersionDifference) is off
    /// the order for a different reason, given on that variant.
    ///
    /// A `floor` that is not on the order therefore answers `false` for every
    /// verdict, including itself.
    pub fn is_at_least(&self, floor: EquivalenceLevel) -> bool {
        match (self.denotational_rank(), floor.denotational_rank()) {
            (Some(here), Some(want)) => here >= want,
            _ => false,
        }
    }

    /// Position on the denotational order, or `None` for a rung that is not on
    /// it. Private: the order is expressed to callers through
    /// [`Self::is_at_least`], so the numbers never become API.
    fn denotational_rank(&self) -> Option<u8> {
        match self {
            EquivalenceLevel::Identical => Some(3),
            EquivalenceLevel::CrossAxisSequenceMatch => Some(2),
            EquivalenceLevel::SequenceMatch => Some(1),
            EquivalenceLevel::NormalizedMatch
            | EquivalenceLevel::AccessionVersionDifference
            | EquivalenceLevel::NotEquivalent
            | EquivalenceLevel::Indeterminate => None,
        }
    }

    /// Returns a human-readable description of the equivalence level.
    pub fn description(&self) -> &'static str {
        match self {
            EquivalenceLevel::Identical => "Identical representation",
            EquivalenceLevel::NormalizedMatch => "Equivalent after normalization",
            EquivalenceLevel::CrossAxisSequenceMatch => {
                "Equivalent by resulting sequence on every determined axis"
            }
            EquivalenceLevel::SequenceMatch => "Equivalent by resulting sequence",
            EquivalenceLevel::AccessionVersionDifference => {
                "Same variant, different accession versions"
            }
            EquivalenceLevel::NotEquivalent => "Not equivalent",
            EquivalenceLevel::Indeterminate => "Undecidable: no computable denotation",
        }
    }
}

impl std::fmt::Display for EquivalenceLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.description())
    }
}

/// Result of an equivalence check with additional details.
///
/// Marked `#[non_exhaustive]` so a future detail field is not a breaking
/// change. Construct it with [`EquivalenceResult::new`] plus the `with_*`
/// builders rather than a struct literal; the fields stay `pub`, so anything
/// the builders do not cover (setting only one of the two normalized forms,
/// which [`EquivalenceResult::with_normalized`] sets as a pair) is still
/// reachable by assigning to the field directly. Mirrors the same attribute on
/// [`EquivalenceLevel`] and
/// `NormalizeResult` (#1033).
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct EquivalenceResult {
    /// The determined equivalence level.
    pub level: EquivalenceLevel,
    /// The normalized form of the first variant (if normalization was performed).
    pub normalized_first: Option<String>,
    /// The normalized form of the second variant (if normalization was performed).
    pub normalized_second: Option<String>,
    /// Additional notes about the comparison.
    pub notes: Vec<String>,
}

impl EquivalenceResult {
    /// Create a new equivalence result.
    pub fn new(level: EquivalenceLevel) -> Self {
        Self {
            level,
            normalized_first: None,
            normalized_second: None,
            notes: Vec::new(),
        }
    }

    /// Add normalized forms to the result.
    pub fn with_normalized(mut self, first: String, second: String) -> Self {
        self.normalized_first = Some(first);
        self.normalized_second = Some(second);
        self
    }

    /// Add a note to the result.
    pub fn with_note(mut self, note: impl Into<String>) -> Self {
        self.notes.push(note.into());
        self
    }

    /// Returns true if the variants are considered equivalent.
    pub fn is_equivalent(&self) -> bool {
        self.level.is_equivalent()
    }
}

/// Checks equivalence between HGVS variants.
///
/// The checker determines if two variants represent the same genomic change,
/// even if they have different string representations.
pub struct EquivalenceChecker<P: ReferenceProvider> {
    normalizer: Normalizer<P>,
}

impl<P: ReferenceProvider> EquivalenceChecker<P> {
    /// Create a new equivalence checker with a reference provider.
    ///
    /// # Arguments
    ///
    /// * `provider` - Reference sequence provider for normalization
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::MockProvider;
    /// use ferro_hgvs::equivalence::EquivalenceChecker;
    ///
    /// let provider = MockProvider::with_test_data();
    /// let checker = EquivalenceChecker::new(provider);
    /// ```
    pub fn new(provider: P) -> Self {
        Self {
            normalizer: Normalizer::new(provider),
        }
    }

    /// Create a new equivalence checker with a custom normalizer configuration.
    ///
    /// # Arguments
    ///
    /// * `provider` - Reference sequence provider
    /// * `config` - Normalization configuration
    pub fn with_config(provider: P, config: NormalizeConfig) -> Self {
        Self {
            normalizer: Normalizer::with_config(provider, config),
        }
    }

    /// Check if two variants are equivalent.
    ///
    /// This method compares two variants and determines their level of equivalence.
    /// It first checks for string identity, then normalizes both variants and
    /// compares the normalized forms.
    ///
    /// # Arguments
    ///
    /// * `v1` - First variant
    /// * `v2` - Second variant
    ///
    /// # Returns
    ///
    /// * `Ok(EquivalenceResult)` - The equivalence result with details
    /// * `Err(FerroError)` - If normalization fails
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::{parse_hgvs, MockProvider};
    /// use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};
    ///
    /// let provider = MockProvider::with_test_data();
    /// let checker = EquivalenceChecker::new(provider);
    ///
    /// let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
    /// let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
    ///
    /// let result = checker.check(&v1, &v2).unwrap();
    /// assert!(result.is_equivalent());
    /// ```
    pub fn check(
        &self,
        v1: &HgvsVariant,
        v2: &HgvsVariant,
    ) -> Result<EquivalenceResult, FerroError> {
        let str1 = v1.to_string();
        let str2 = v2.to_string();

        // Check for identical string representation
        if str1 == str2 {
            return Ok(EquivalenceResult::new(EquivalenceLevel::Identical)
                .with_normalized(str1.clone(), str2.clone()));
        }

        // Check for accession version differences
        if let Some(result) = self.check_accession_version_difference(v1, v2, &str1, &str2) {
            return Ok(result);
        }

        // Normalize both variants and compare
        let norm1 = self.normalizer.normalize(v1)?;
        let norm2 = self.normalizer.normalize(v2)?;

        let norm_str1 = norm1.to_string();
        let norm_str2 = norm2.to_string();

        if norm_str1 == norm_str2 {
            Ok(EquivalenceResult::new(EquivalenceLevel::NormalizedMatch)
                .with_normalized(norm_str1, norm_str2)
                .with_note("Variants normalize to the same form"))
        } else {
            // Check if they might still be equivalent with version differences
            if let Some(result) =
                self.check_accession_version_difference(&norm1, &norm2, &norm_str1, &norm_str2)
            {
                return Ok(result.with_note("Equivalent after normalization, different versions"));
            }

            // A description that names no definite bases denotes no sequence,
            // so no rung below may report a positive denotational verdict about
            // it. The check runs on the **inputs**, here, and not only on the
            // normalized pair `denotational_verdict` re-checks, because
            // normalization *expands* `insN[10]` into a literal run of `N`s: by
            // the time the sequence rung compares two reconstructed windows, two
            // indeterminate payloads can be byte-equal, and the rung would
            // report agreement on a pair where neither side denotes anything to
            // agree about.
            //
            // No input reaching `SequenceVerdict::Same` this way has been
            // measured — every constructed candidate either converges at
            // `NormalizedMatch` above or is caught by the same check on the
            // normalized pair. This is therefore stated as an
            // invariant rather than as a fix for an observed wrong answer, and
            // `an_indeterminate_input_never_wins_a_decided_denotational_rung`
            // pins it as one: the rule is a property of `check`, not of
            // whichever path a given pair happens to take.
            for side in [v1, v2] {
                if let Some(reason) = indeterminate_reason(side) {
                    return Ok(EquivalenceResult::new(EquivalenceLevel::Indeterminate)
                        .with_normalized(norm_str1, norm_str2)
                        .with_note(reason));
                }
            }

            self.denotational_verdict(&norm1, &norm2, norm_str1, norm_str2)
        }
    }

    /// Decide the **denotational** relation between two descriptions, as
    /// written, without ever asking whether they normalize to one string.
    ///
    /// This is the relation
    /// `rulings[confluence-gate-is-apply-equality-on-every-determined-axis]`
    /// defines: apply-equality on every *determined* axis, a relation over
    /// `apply` whose codomain is bases. [`Self::check`] answers a wider
    /// question — it is the whole ladder, textual rungs included, and it
    /// short-circuits at [`EquivalenceLevel::NormalizedMatch`] the moment the
    /// two normalized forms coincide. That short-circuit is correct for a
    /// caller asking "are these the same?" and useless to one asking "did the
    /// denotational comparison agree?", because a converged pair never reaches
    /// a denotational rung at all — and `NormalizedMatch` is deliberately off
    /// the order, so `is_at_least(CrossAxisSequenceMatch)` *rejects* it.
    ///
    /// So this entry point starts one rung lower and reports what the
    /// comparison found:
    ///
    /// * **Nothing is normalized.** The two variants are compared exactly as
    ///   handed over, so a caller that wants normalized forms compared must
    ///   normalize them itself. That is what makes the answer independent of
    ///   how confluent the normalizer happens to be.
    /// * **There is no `Identical` rung here.** Two identical descriptions are
    ///   re-derived and compared like any other pair, which reports
    ///   [`EquivalenceLevel::CrossAxisSequenceMatch`] when every determined axis
    ///   could be computed, and a lower rung when one could not.
    /// * **`AccessionVersionDifference` is still answered**, because it needs no
    ///   normalizer; without it, one variant spelled on two versions of an
    ///   accession would fall through to a `NotEquivalent` no rung examined.
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};
    /// use ferro_hgvs::{parse_hgvs, MockProvider};
    ///
    /// let checker = EquivalenceChecker::new(MockProvider::with_test_data());
    /// let v = parse_hgvs("NM_000088.3:c.10del").unwrap();
    ///
    /// // `check` stops at the textual rung; the denotational entry point
    /// // re-derives the bases and reports what it compared.
    /// assert_eq!(checker.check(&v, &v).unwrap().level, EquivalenceLevel::Identical);
    /// assert!(checker.compare_denotations(&v, &v).unwrap().level.is_decided());
    /// ```
    pub fn compare_denotations(
        &self,
        v1: &HgvsVariant,
        v2: &HgvsVariant,
    ) -> Result<EquivalenceResult, FerroError> {
        let str1 = v1.to_string();
        let str2 = v2.to_string();
        if let Some(result) = self.check_accession_version_difference(v1, v2, &str1, &str2) {
            return Ok(result);
        }

        // The same hoist [`Self::check`] performs above its own call into
        // `denotational_verdict`, and for the same reason: a description that
        // names no definite bases denotes no sequence, so no rung below may
        // report a positive denotational verdict about it.
        //
        // It cannot be left to the fall-through check inside
        // `denotational_verdict`, because that one runs *after* the sequence
        // rung has had its chance to return. Two byte-equal indeterminate
        // payloads reconstruct to byte-equal windows, so the rung answers
        // `SequenceVerdict::Same` and returns before the fall-through is
        // reached. Measured: `NC_…:g.1001_1002insNNN` against itself reported
        // `CrossAxisSequenceMatch` — a *decided positive* about a pair where
        // neither side denotes anything, and one that
        // `is_at_least(CrossAxisSequenceMatch)` accepts, i.e. exactly what a
        // confluence gate written as the module docs prescribe would swallow.
        // `check` never showed it because its string-identity rung answers
        // `Identical` first; this entry point has no such rung by design.
        //
        // Pinned by
        // `an_indeterminate_input_never_wins_a_decided_denotational_rung`,
        // which drives both entry points — the invariant is a property of the
        // denotational rungs, not of whichever entry point reaches them.
        for side in [v1, v2] {
            if let Some(reason) = indeterminate_reason(side) {
                return Ok(EquivalenceResult::new(EquivalenceLevel::Indeterminate)
                    .with_normalized(str1, str2)
                    .with_note(reason));
            }
        }

        self.denotational_verdict(v1, v2, str1, str2)
    }

    /// The denotational half of the ladder: the sequence rung, the cross-axis
    /// strengthening above it, and the three ways of declining below it.
    ///
    /// Shared verbatim by [`Self::check`] — which reaches it with the
    /// *normalized* pair — and by [`Self::compare_denotations`], which reaches
    /// it with the pair as written. `str1`/`str2` are the two descriptions'
    /// rendered forms, recorded on the result so a caller can see which pair a
    /// verdict was reached over.
    fn denotational_verdict(
        &self,
        v1: &HgvsVariant,
        v2: &HgvsVariant,
        str1: String,
        str2: String,
    ) -> Result<EquivalenceResult, FerroError> {
        // Sequence-level equivalence (issue #1158): two variants may normalize
        // to different HGVS strings yet produce the same edited reference
        // sequence — e.g. a length-changing `delins` vs a decomposed cis allele
        // of the same edit, which ferro deliberately keeps in distinct
        // canonical forms. Reconstruct each edited window from SPDI triples and
        // compare. This is best-effort: any variant that cannot be projected
        // (unsupported edit, missing reference, mixed accessions, or a
        // multi-molecule allele) simply declines, so the rung only ever
        // upgrades a `NotEquivalent` verdict.
        match self.same_resulting_sequence(v1, v2) {
            SequenceVerdict::Same => {
                let (level, note) = self.strengthen_across_axes(v1, v2);
                return Ok(EquivalenceResult::new(level)
                    .with_normalized(str1, str2)
                    .with_note(note));
            }
            // The reference window needed to compare them is wider than the
            // cap, so nothing was reconstructed and nothing was compared.
            SequenceVerdict::WindowTooWide => {
                return Ok(EquivalenceResult::new(EquivalenceLevel::Indeterminate)
                    .with_normalized(str1, str2)
                    .with_note(format!(
                        "The reference window covering both variants exceeds \
                         {MAX_SEQUENCE_COMPARE_WINDOW} bases, so neither resulting \
                         sequence was reconstructed"
                    )));
            }
            // The window is bounded and both sides convert to SPDI on a shared
            // accession, but the provider could not serve its bases, so — as with
            // `WindowTooWide` — nothing was reconstructed and nothing was
            // compared (#1989). This is the "want of reference data" decline, and
            // it must report `Indeterminate` for the same reason its sibling does:
            // falling through to the `NotEquivalent` tail would assert a decided
            // negative about a pair that was never examined. Pinned by
            // `issue_1989_declined_is_indeterminate`.
            SequenceVerdict::ReferenceUnavailable => {
                return Ok(EquivalenceResult::new(EquivalenceLevel::Indeterminate)
                    .with_normalized(str1, str2)
                    .with_note(
                        "The reference window covering both variants could not be read from the \
                         provider, so neither resulting sequence was reconstructed"
                            .to_string(),
                    ));
            }
            // `Different` is a decided negative — both sides were reconstructed
            // and disagree. `Declined` is the catch-all decline. The shapes the
            // checker routes to `NotEquivalent` *deliberately* are an edit SPDI
            // cannot carry (`issue_1578_followup_equivalence_rungs` for a
            // fusion), a cross-accession pair (the mixed-accession arm of
            // `same_resulting_sequence`) and an un-appliable overlapping allele
            // (`issue_1244_equivalence_overlap_panic`), each pinned there.
            //
            // Those are not all of it, and the arm above does not make them so.
            // `edit_triples` discards its conversion error —
            // `hgvs_to_spdi(member, provider).ok()?`, `checker.rs:1114` — so a
            // member whose conversion needed the reference and could not read it
            // also arrives here as `Declined` and is read as a negative:
            // `NC_ABSENT.1:g.2del` vs `g.2delC` answers `NotEquivalent` against
            // a provider serving no bases and `CrossAxisSequenceMatch` against a
            // served contig. So `ReferenceUnavailable` above is the
            // reference-level decline reachable from the WINDOW FETCH, not every
            // reference-level decline; the rest is #2056, a follow-up on
            // `edit_triples`. These two fall through.
            SequenceVerdict::Different | SequenceVerdict::Declined => {}
        }

        // #1578 follow-up: `NotEquivalent` is a *positive claim* that two
        // descriptions denote different variants, and it is only sound if
        // some rung above actually examined them. For a structural
        // rearrangement both deciding rungs are blind at once — see
        // `undecidable_reason` — so falling through here would answer a
        // question nobody asked. Decline instead, which is what the other
        // three modules that hand-roll `Allele` recursion already do with a
        // ring (`spdi::apply`, `vcf::from_hgvs`, `project::projector`).
        //
        // `Indeterminate` below says the same thing as a *verdict* rather
        // than as an error. The two are kept apart deliberately: the error
        // is the shipped contract for a ring and a `sup` marker (pinned by
        // `issue_1578_followup_equivalence_rungs`), and turning it into a
        // verdict would be a silent API change for every caller that
        // matches on the `Err`. A future change may fold them together;
        // until then, read them as one idea with two spellings.
        for side in [v1, v2] {
            if let Some(reason) = self.undecidable_reason(side) {
                return Err(FerroError::UnsupportedVariant {
                    variant_type: reason,
                });
            }
        }

        // Nothing above could compute a denotation for one of the sides, so
        // there is no negative to report either.
        //
        // Only the pair handed to *this* function is examined. On
        // [`Self::check`]'s path that is the **normalized** pair, and its
        // inputs were already cleared above the sequence rung — re-checking
        // them would be dead work, and worse, would read as though this were
        // the only place the rule is applied. The normalized forms still have
        // to be checked on their own account, in the other direction from that
        // hoisted check: normalization can *introduce* an indeterminate payload
        // the input did not name, by resolving a shape into a literal run of
        // `N`s.
        for side in [v1, v2] {
            if let Some(reason) = indeterminate_reason(side) {
                return Ok(EquivalenceResult::new(EquivalenceLevel::Indeterminate)
                    .with_normalized(str1, str2)
                    .with_note(reason));
            }
        }

        // The note says "no rung found them equivalent" rather than the older
        // "do not normalize to the same form", because this function is also
        // reached from [`Self::compare_denotations`], where nothing was
        // normalized and the old wording named a comparison that never ran.
        Ok(EquivalenceResult::new(EquivalenceLevel::NotEquivalent)
            .with_normalized(str1, str2)
            .with_note("The two descriptions differ and no rung found them equivalent"))
    }

    /// Given that the two variants agree on the axis their descriptions are
    /// written in, decide whether they agree on **every** determined axis.
    ///
    /// Returns the rung to report and the note to attach. A transcript
    /// description determines a second axis — the genome its exon alignment
    /// carries it to — so agreement there is what separates
    /// [`EquivalenceLevel::CrossAxisSequenceMatch`] from
    /// [`EquivalenceLevel::SequenceMatch`]. A genomic description determines
    /// only the one axis already compared, so it reaches the stronger rung
    /// directly; without that, a gate stated as "at least
    /// `CrossAxisSequenceMatch`" would reject every `g.` pair.
    ///
    /// **Never consults the normalizer.** The second axis is derived from
    /// [`CoordinateMapper`], which holds nothing but a `&Transcript` and does
    /// exon arithmetic over `exon.genomic_start`/`genomic_end` and
    /// `transcript.strand` — reference data, by construction. The one
    /// whole-variant API that *would* reintroduce the circularity,
    /// `VariantProjector::project_to_genomic_normalized`, is deliberately not
    /// used here.
    fn strengthen_across_axes(
        &self,
        v1: &HgvsVariant,
        v2: &HgvsVariant,
    ) -> (EquivalenceLevel, String) {
        match (self.second_axis(v1), self.second_axis(v2)) {
            (SecondAxis::None, SecondAxis::None) => (
                EquivalenceLevel::CrossAxisSequenceMatch,
                "Variants produce the same resulting sequence on the one axis they determine"
                    .to_string(),
            ),
            (SecondAxis::Genomic(acc1, t1), SecondAxis::Genomic(acc2, t2)) if acc1 == acc2 => {
                match compare_triples(|start, end| self.fetch_window(&acc1, start, end), &t1, &t2) {
                    SequenceVerdict::Same => (
                        EquivalenceLevel::CrossAxisSequenceMatch,
                        "Variants produce the same resulting sequence on every determined axis"
                            .to_string(),
                    ),
                    SequenceVerdict::Different => (
                        EquivalenceLevel::SequenceMatch,
                        format!(
                            "Variants produce the same resulting sequence on their own axis but \
                             not on {acc1}, so they are not the same variant"
                        ),
                    ),
                    // **Not `Indeterminate`, deliberately.** That rung means no
                    // denotation could be computed for a side; here both sides
                    // denote, the own-axis comparison ran, and it agreed. What
                    // is missing is corroboration, not a denotation, and
                    // reporting "cannot tell" would discard a measured fact.
                    // Nothing is over-claimed by staying decided, because a
                    // gate reads the *level*: `is_at_least(
                    // CrossAxisSequenceMatch)` rejects `SequenceMatch` either
                    // way. The distinction between "the second axis disagreed"
                    // and "it could not be asked" survives in the note below and
                    // nowhere else; carrying it in the level needs a rung of its
                    // own rather than a reuse of `Indeterminate`. Pinned by
                    // `a_second_axis_that_cannot_be_computed_stops_at_the_single_axis_rung`.
                    SequenceVerdict::WindowTooWide
                    | SequenceVerdict::ReferenceUnavailable
                    | SequenceVerdict::Declined => (
                        EquivalenceLevel::SequenceMatch,
                        format!(
                            "Variants produce the same resulting sequence on their own axis; the \
                             {acc1} axis they also determine could not be compared"
                        ),
                    ),
                }
            }
            // Two different contigs is a real disagreement, not a decline.
            (SecondAxis::Genomic(acc1, _), SecondAxis::Genomic(acc2, _)) => (
                EquivalenceLevel::SequenceMatch,
                format!(
                    "Variants produce the same resulting sequence on their own axis but sit on \
                     different genomic contigs ({acc1} and {acc2})"
                ),
            ),
            _ => (
                EquivalenceLevel::SequenceMatch,
                "Variants produce the same resulting sequence on their own axis; a determined \
                 axis could not be computed"
                    .to_string(),
            ),
        }
    }

    /// The axis a description determines *in addition to* the one it is written
    /// in.
    ///
    /// Keyed on the variant kind rather than on whether the provider happens to
    /// hold a transcript under the accession, because that is what
    /// "determined" means: the mapping exists as a fact about the reference,
    /// whether or not this provider can serve it.
    fn second_axis(&self, v: &HgvsVariant) -> SecondAxis {
        match v {
            // A genomic description determines only the genome. It does not
            // determine a transcript — many transcripts overlap one locus, and
            // picking one would be a choice, not a mapping.
            HgvsVariant::Genome(_) | HgvsVariant::Mt(_) | HgvsVariant::Circular(_) => {
                SecondAxis::None
            }
            // `p.` never gets here: translation is many-to-one, `p.` states a
            // consequence rather than a denotation, and `hgvs_to_spdi` refuses
            // it — so `edit_triples` already declined upstream.
            HgvsVariant::Cds(_) | HgvsVariant::Tx(_) | HgvsVariant::Rna(_) => {
                match self.edit_triples(v) {
                    Some((accession, triples)) => {
                        match self.transcript_triples_to_genomic(&accession, &triples) {
                            Some((contig, mapped)) => SecondAxis::Genomic(contig, mapped),
                            None => SecondAxis::NotComputed,
                        }
                    }
                    None => SecondAxis::NotComputed,
                }
            }
            // `edit_triples` has already established that a cis allele's
            // members share one accession, so the first member fixes the axis.
            HgvsVariant::Allele(allele) => match allele.variants.first() {
                Some(_) if allele.phase != AllelePhase::Cis => SecondAxis::NotComputed,
                Some(first) => match self.second_axis(first) {
                    SecondAxis::None => SecondAxis::None,
                    SecondAxis::NotComputed => SecondAxis::NotComputed,
                    // Re-derive from the whole allele so every member is mapped
                    // together, rather than reporting only the first member's
                    // triples.
                    SecondAxis::Genomic(..) => match self.edit_triples(v) {
                        Some((accession, triples)) => {
                            match self.transcript_triples_to_genomic(&accession, &triples) {
                                Some((contig, mapped)) => SecondAxis::Genomic(contig, mapped),
                                None => SecondAxis::NotComputed,
                            }
                        }
                        None => SecondAxis::NotComputed,
                    },
                },
                None => SecondAxis::NotComputed,
            },
            _ => SecondAxis::NotComputed,
        }
    }

    /// Re-express transcript-axis SPDI triples on the genomic axis, using the
    /// transcript's exon alignment and nothing else.
    ///
    /// Declines (returns `None`) whenever the genomic counterpart is not one
    /// contiguous triple:
    ///
    /// * a deleted span that crosses an intron — the genomic edit is then two
    ///   or more disjoint pieces, which one SPDI triple cannot hold;
    /// * a pure insertion whose two transcript anchors are not adjacent on the
    ///   genome, i.e. the insertion point *is* an exon junction. This is the
    ///   `c.3921dup` half of the spec's own worked case
    ///   (`DNA/duplication.md:148`): the junction-side spelling has no single
    ///   genomic interbase counterpart, while the far-side spelling has one, so
    ///   the pair cannot reach cross-axis agreement.
    fn transcript_triples_to_genomic(
        &self,
        tx_accession: &str,
        triples: &[SpdiVariant],
    ) -> Option<(String, Vec<SpdiVariant>)> {
        let transcript = self
            .normalizer
            .provider()
            .get_transcript(tx_accession)
            .ok()?;
        let contig = transcript.chromosome.clone()?;
        let strand = transcript.strand;
        let mapper = CoordinateMapper::new(&transcript);
        let to_genomic = |base: u64| -> Option<u64> {
            mapper
                .tx_to_genomic(&crate::hgvs::location::TxPos::new(
                    i64::try_from(base).ok()?,
                ))
                .ok()
                .flatten()
        };

        let mut mapped = Vec::with_capacity(triples.len());
        for triple in triples {
            let (start, deletion, insertion) = if triple.deletion.is_empty() {
                // A pure insertion sits at interbase `position`, between
                // transcript bases `position` and `position + 1` (1-based).
                if triple.position == 0 {
                    return None;
                }
                let left = to_genomic(triple.position)?;
                let right = to_genomic(triple.position + 1)?;
                match strand {
                    Strand::Plus if right == left + 1 => {
                        (left, String::new(), triple.insertion.clone())
                    }
                    Strand::Minus if left == right + 1 => {
                        (right, String::new(), reverse_complement(&triple.insertion))
                    }
                    _ => return None,
                }
            } else {
                let length = triple.deletion.len() as u64;
                let first = to_genomic(triple.position + 1)?;
                let last = to_genomic(triple.position + length)?;
                // `checked_sub` turns the 1-based genomic position into SPDI's
                // 0-based interbase one, and declines rather than underflowing
                // if the span's low end is genomic 0. Nothing between a
                // provider and here asserts `Exon::genomic_start >= 1` — the
                // cdot path takes `e[0]` verbatim and the JSON loader takes
                // whatever the file holds — so a 0-based placement reaches this
                // arithmetic, where on `u64` it panics under debug assertions
                // and wraps to `u64::MAX` in release. There is no interbase
                // point to the left of genomic 0, so declining is the same
                // answer the insertion arm above gives `triple.position == 0`.
                match strand {
                    Strand::Plus if last >= first && last - first + 1 == length => (
                        first.checked_sub(1)?,
                        triple.deletion.clone(),
                        triple.insertion.clone(),
                    ),
                    Strand::Minus if first >= last && first - last + 1 == length => (
                        last.checked_sub(1)?,
                        reverse_complement(&triple.deletion),
                        reverse_complement(&triple.insertion),
                    ),
                    _ => return None,
                }
            };
            mapped.push(SpdiVariant::new(contig.clone(), start, deletion, insertion));
        }
        Some((contig, mapped))
    }

    /// Check if two variants represent the same change but on different accession versions.
    fn check_accession_version_difference(
        &self,
        v1: &HgvsVariant,
        v2: &HgvsVariant,
        str1: &str,
        str2: &str,
    ) -> Option<EquivalenceResult> {
        // This helper only detects whether two variants are the *same change on
        // different accession versions*. `None` here means "not an accession-version
        // difference" — it is not an "unsupported" signal. For compound/null/unknown
        // alleles there is no single accession version to compare (an allele's
        // `accession()` reports only its first member, which would be misleading), so
        // we decline. The caller (`check`) then falls through to the normalizer, which
        // handles alleles directly, so allele equivalence is fully determined — just
        // not via this version-comparison path.
        let acc1 = match v1 {
            HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => return None,
            HgvsVariant::Allele(_) => return None,
            _ => v1.accession()?,
        };
        let acc2 = match v2 {
            HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => return None,
            HgvsVariant::Allele(_) => return None,
            _ => v2.accession()?,
        };

        // Check if base accessions (without version) match but versions differ
        if acc1.prefix == acc2.prefix && acc1.number == acc2.number && acc1.version != acc2.version
        {
            // Extract the part after the accession to compare the variant itself
            let variant_part1 = extract_variant_part(str1);
            let variant_part2 = extract_variant_part(str2);

            if variant_part1 == variant_part2 {
                return Some(
                    EquivalenceResult::new(EquivalenceLevel::AccessionVersionDifference)
                        .with_normalized(str1.to_string(), str2.to_string())
                        .with_note(format!(
                            "Same variant on different versions: {} vs {}",
                            acc1.version
                                .map(|v| v.to_string())
                                .unwrap_or_else(|| "no version".to_string()),
                            acc2.version
                                .map(|v| v.to_string())
                                .unwrap_or_else(|| "no version".to_string())
                        )),
                );
            }
        }

        None
    }

    /// Decide whether two variants, applied to their (shared) reference, produce
    /// the same edited sequence on the axis their descriptions are written in.
    ///
    /// Best-effort and side-effect free. [`SequenceVerdict::Declined`] covers
    /// every variant that cannot be projected to SPDI triples on a single
    /// shared accession. Two declines are split out of it because they must not
    /// be read as a negative: [`SequenceVerdict::WindowTooWide`] when the union
    /// window exceeds [`MAX_SEQUENCE_COMPARE_WINDOW`], and
    /// [`SequenceVerdict::ReferenceUnavailable`] when the provider cannot serve
    /// the window's bases (#1989). Both are resolved to
    /// [`EquivalenceLevel::Indeterminate`] by the caller.
    ///
    /// **The split covers the window fetch, and not every reference failure.**
    /// [`Self::edit_triples`] discards its conversion error —
    /// `hgvs_to_spdi(member, provider).ok()?`, `checker.rs:1114` — so a member
    /// whose SPDI conversion had to read the reference and could not is `None`
    /// there and `Declined` here, indistinguishable from an edit SPDI cannot
    /// carry. Measured on this revision, against a provider serving no bases:
    /// `NC_ABSENT.1:g.2del` vs `g.2delC` — one edit, two spellings — answers
    /// `NotEquivalent`, while the same pair on a served contig answers
    /// `CrossAxisSequenceMatch`. Narrowing that is #2056, a follow-up on
    /// `edit_triples`, not on this function.
    fn same_resulting_sequence(&self, v1: &HgvsVariant, v2: &HgvsVariant) -> SequenceVerdict {
        let Some((acc1, triples1)) = self.edit_triples(v1) else {
            return SequenceVerdict::Declined;
        };
        let Some((acc2, triples2)) = self.edit_triples(v2) else {
            return SequenceVerdict::Declined;
        };
        // A resulting sequence is only comparable on the same reference. Two
        // different accessions fail the relation's first conjunct outright,
        // which is a decided negative rather than something undecidable.
        if acc1 != acc2 {
            return SequenceVerdict::Declined;
        }
        compare_triples(
            |start, end| self.fetch_window(&acc1, start, end),
            &triples1,
            &triples2,
        )
    }

    /// Why no rung can decide equivalence for `variant`, or `None` if some rung
    /// can.
    ///
    /// The checker decides in rungs, and a structural rearrangement defeats both
    /// of the two that compare *meaning* rather than text:
    ///
    /// - [`Normalizer::normalize`] returns a genome ring and a `sup` marker
    ///   unchanged (the pass-through arms in `src/normalize/mod.rs`), so
    ///   comparing normalized forms degenerates into the string comparison the
    ///   `Identical` rung already made; and
    /// - [`hgvs_to_spdi`] cannot represent either, so [`Self::edit_triples`]
    ///   declines and the `SequenceMatch` rung (#1158) cannot fire either.
    ///
    /// With both blind, `NotEquivalent` would be asserted having evaluated
    /// nothing. The two *textual* rungs — `Identical` and
    /// `AccessionVersionDifference` — are unaffected and still answer, because
    /// they need neither normalization nor a provider.
    ///
    /// **The SPDI conjunct is not redundant, and it is what keeps this honest.**
    /// It makes the predicate self-retiring: if `hgvs_to_spdi` ever gains ring
    /// support, `edit_triples` starts succeeding, the `SequenceMatch` rung
    /// begins deciding these pairs, and this decline stops firing on its own
    /// rather than masking the new capability. `Circular` (`o.`) is the live
    /// proof that the conjunct matters: its normalization is a pass-through too,
    /// yet SPDI represents it, so it is decided by the sequence rung and must
    /// never reach this decline.
    /// **Scoped to the two shapes a defect was measured on.**
    /// [`HgvsVariant::RnaFusion`] shares both blindnesses and is deliberately
    /// *not* listed: no pair of fusion spellings denoting one fusion has been
    /// exhibited, so widening the decline to `::` transcript fusions would
    /// refuse pairs the checker answers today on no evidence that any of those
    /// answers is wrong. Add it when a wrong answer is measured, with the pair
    /// that measures it.
    ///
    /// **A third shape now shares the SPDI blindness and is deliberately not
    /// listed either: a genomic `pter`/`qter`/`cen` description (#1643).**
    /// `hgvs_to_spdi` used to flatten those onto `base: 0` and convert; it now
    /// refuses, so [`Self::edit_triples`]'s `.ok()?` swallows the new error and
    /// the `SequenceMatch` rung stops firing for them. The reach is narrower
    /// than it looks, and the reason is worth writing down because it is not
    /// visible from this file: both rungs run on the **normalized** forms, and
    /// normalization *resolves* a marker against the provider
    /// (`g.10_qterdelACGTACGTAC` -> `g.10_28del` on a 28-base contig), so most
    /// of the class never reaches the sequence rung carrying one at all.
    ///
    /// **Which markers survive normalization is a property of the PROVIDER, not
    /// of the marker.** `cen` always does — `resolve_special_genome_pos` returns
    /// `Ok(None)` for it unconditionally, there being no centromere annotation
    /// on a sequence accession to resolve against. But `pter` and `qter` survive
    /// too whenever the length lookup fails: that arm falls through to
    /// `canonicalize_genome_variant`, which rewrites only the edit body and
    /// clones the location, so the marker is echoed. Measured, on a provider
    /// holding no sequence for the accession: `g.10_qterdelACGTACGTAC` ->
    /// `g.10_qterdel` and `g.pter_10delACGTACGTAC` -> `g.pter_10del`. So this is
    /// "`cen` plus any marker whose reference is unservable", not "`cen`".
    ///
    /// **One verdict class moved, and the old answer was wrong.** Both figures
    /// below come from disabling the guard in place and re-running, and they are
    /// **two separate corpora** — do not read either denominator onto the other:
    ///
    /// | corpus | pairs | verdicts that moved |
    /// |---|---|---|
    /// | 22 mixed `cen`/`qter`/`pter` descriptors, on a serving and a non-serving provider | 231 x 2 | **0** |
    /// | 10 descriptors whose 3' anchor is `cen`, serving provider | 45 | **3** |
    ///
    /// The three are what the first corpus could not reach, and they say why: a
    /// `cen` on the END of a *duplication* is what zeroes the triple, because
    /// the `Duplication` arm emits `end_one_based` as the SPDI position. So
    /// `g.10_cendupACGTACGTAC`, `g.12_cendupACGTACGTAC` and
    /// `g.4_cendupACGTACGTAC` all converted to `NC_000001.11:0::ACGTACGTAC`, and
    /// each of the three pairs among them answered `SequenceMatch` — **wrong**,
    /// since they duplicate different spans. All three now answer
    /// `NotEquivalent`. The `del`/`inv` spellings of the same shape never
    /// collided: those arms emit `start_one_based - 1`, which the `cen` does not
    /// touch.
    ///
    /// So the sibling in such a pair has to be another special position, and an
    /// earlier revision of this note gave `g.10_19dupACGTACGTAC` instead. That
    /// is not a case: on a fully numeric anchor the same arm emits
    /// `end_one_based` as a real coordinate, so it converts at **19**, and two
    /// triples at different positions cannot match however blind the rung is.
    ///
    /// **No genuinely-equivalent pair lost an answer**, which is why nothing
    /// here changed. Every equivalent `cen` pair tried (cis decompositions,
    /// member reorderings, delins/substitution respellings) is decided one rung
    /// earlier by `NormalizedMatch`, because normalization's own cis merge
    /// collapses them, and no pair moved in the `NotEquivalent` -> equivalent
    /// direction in either corpus. The residual hazard is the honesty one this
    /// predicate exists for — a `NotEquivalent` reached with the sequence rung
    /// blind — not a wrong answer left standing, so per the standard above it
    /// earns this note rather than a decline.
    fn undecidable_reason(&self, variant: &HgvsVariant) -> Option<String> {
        // SPDI first, at every level: if `edit_triples` succeeds the `SequenceMatch`
        // rung can decide, so there is nothing undecidable to report. Keeping this
        // ahead of the structural test is what makes the predicate self-retiring —
        // ring support in `hgvs_to_spdi` silently switches these pairs back to being
        // answered instead of declined.
        if self.edit_triples(variant).is_some() {
            return None;
        }
        // A cis allele carrying a ring is the same blindness one level down, and it
        // is a shape this codebase constructs deliberately —
        // `issue_1578_followup_ring_declines.rs` pairs a legal genomic member with a
        // ring and requires `spdi::apply`, `vcf::from_hgvs` and the projector to all
        // decline it. Without this recursion the equivalence checker was the one
        // module that answered it: two such alleles compared `NotEquivalent`, which
        // is a positive claim that no rung examined.
        if let HgvsVariant::Allele(allele) = variant {
            return allele
                .variants
                .iter()
                .find_map(|member| self.undecidable_reason(member));
        }
        if !matches!(
            variant,
            HgvsVariant::GenomeRing(_) | HgvsVariant::Supernumerary(_)
        ) {
            return None;
        }
        Some(format!(
            "{variant}: cannot decide equivalence — normalization passes this \
             description through unchanged and SPDI cannot represent it, so \
             neither the normalized-form nor the resulting-sequence comparison \
             examined it"
        ))
    }

    /// Project a variant to the SPDI primitive edits that make up its resulting
    /// sequence, together with the single accession they act on.
    ///
    /// Returns `None` when a single resulting sequence is undefined or cannot be
    /// derived: multi-molecule alleles (trans / mosaic / chimeric / and-or /
    /// products / unknown phase), null/unknown alleles, edits SPDI cannot
    /// represent, or members that span more than one accession.
    fn edit_triples(&self, v: &HgvsVariant) -> Option<(String, Vec<SpdiVariant>)> {
        let members: Vec<&HgvsVariant> = match v {
            HgvsVariant::Allele(allele) => {
                // A single resulting sequence is only well-defined for a cis
                // allele — all edits applied to the same molecule.
                if allele.phase != AllelePhase::Cis {
                    return None;
                }
                allele.variants.iter().collect()
            }
            HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => return None,
            single => vec![single],
        };
        if members.is_empty() {
            return None;
        }

        let mut accession: Option<String> = None;
        let mut triples = Vec::with_capacity(members.len());
        for member in members {
            let spdi = hgvs_to_spdi(member, self.normalizer.provider()).ok()?;
            match &accession {
                None => accession = Some(spdi.sequence.clone()),
                Some(acc) if *acc != spdi.sequence => return None,
                Some(_) => {}
            }
            triples.push(spdi);
        }
        accession.map(|acc| (acc, triples))
    }

    /// Fetch reference bases for the 0-based half-open interval
    /// `[start, end)` on `accession`, trying the genomic provider first and
    /// falling back to the transcript-sequence provider (mirroring the SPDI
    /// converter's own fetch). Returns `None` on any provider error or short
    /// read.
    fn fetch_window(&self, accession: &str, start: u64, end: u64) -> Option<String> {
        let provider = self.normalizer.provider();
        let bases = provider
            .get_genomic_sequence(accession, start, end)
            .or_else(|_| provider.get_sequence(accession, start, end))
            .ok()?;
        if bases.len() as u64 != end - start {
            return None;
        }
        Some(bases)
    }

    /// Check if multiple variants are all equivalent to each other.
    ///
    /// # Arguments
    ///
    /// * `variants` - Slice of variants to compare
    ///
    /// # Returns
    ///
    /// * `true` if all variants are equivalent
    /// * `false` if any pair is not equivalent
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::{parse_hgvs, MockProvider};
    /// use ferro_hgvs::equivalence::EquivalenceChecker;
    ///
    /// let provider = MockProvider::with_test_data();
    /// let checker = EquivalenceChecker::new(provider);
    ///
    /// let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
    /// let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
    /// let v3 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
    ///
    /// assert!(checker.all_equivalent(&[v1, v2, v3]).unwrap());
    /// ```
    pub fn all_equivalent(&self, variants: &[HgvsVariant]) -> Result<bool, FerroError> {
        if variants.len() < 2 {
            return Ok(true);
        }

        let first = &variants[0];
        for variant in &variants[1..] {
            let result = self.check(first, variant)?;
            if !result.is_equivalent() {
                return Ok(false);
            }
        }

        Ok(true)
    }

    /// Group variants by equivalence.
    ///
    /// Returns groups of variants that are equivalent to each other.
    ///
    /// # Arguments
    ///
    /// * `variants` - Slice of variants to group
    ///
    /// # Returns
    ///
    /// Vector of variant groups, where variants in each group are equivalent.
    pub fn group_by_equivalence(
        &self,
        variants: &[HgvsVariant],
    ) -> Result<Vec<Vec<HgvsVariant>>, FerroError> {
        let mut groups: Vec<Vec<HgvsVariant>> = Vec::new();

        for variant in variants {
            let mut found_group = false;

            for group in &mut groups {
                if !group.is_empty() {
                    let result = self.check(&group[0], variant)?;
                    if result.is_equivalent() {
                        group.push(variant.clone());
                        found_group = true;
                        break;
                    }
                }
            }

            if !found_group {
                groups.push(vec![variant.clone()]);
            }
        }

        Ok(groups)
    }
}

/// What a sequence-level comparison concluded.
///
/// Three of the four variants mean "nothing was compared", and they are kept
/// apart because they are not read the same way. `WindowTooWide` and
/// `ReferenceUnavailable` are the two declines [`compare_triples`] raises before
/// it reconstructs anything — the window that would settle the pair was not
/// obtained, once because it exceeds the cap and once because the provider could
/// not serve it — and both must report `Indeterminate`, because arguing either
/// into a decided negative claims a comparison that never ran (#1989).
/// `Declined` is the catch-all; the checker deliberately routes it to
/// `NotEquivalent`, and the shapes it is *meant* to carry — an edit SPDI cannot
/// carry, a cross-accession pair, an un-appliable overlapping allele — are
/// pinned elsewhere. `Different` is the one decided negative that was actually
/// measured: both sides reconstructed, and they disagree.
///
/// **`Declined` is not purely representation-level.** `EquivalenceChecker::edit_triples`
/// discards its conversion error — `hgvs_to_spdi(member, provider).ok()?`,
/// `checker.rs:1114` — so a member whose SPDI conversion had to read the
/// reference and could not is `None` there and `Declined` here, on the same
/// footing as an unrepresentable edit. Measured on this revision:
/// `NC_ABSENT.1:g.2del` vs `g.2delC` (one edit, two spellings) answers
/// `NotEquivalent` against a provider serving no bases and
/// `CrossAxisSequenceMatch` against a served contig; `g.2dup` vs `g.2_3insC`
/// against the unserved provider, and two out-of-bounds substitutions on a
/// *served* 12-base contig, answer `NotEquivalent` by the same route. Each of
/// those declines before the fetch, which is how it is distinguishable from the
/// `ReferenceUnavailable` arm at all. Distinguishing them is #2056, a follow-up
/// on `edit_triples`; the split below is scoped to [`compare_triples`]'s fetch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SequenceVerdict {
    /// Both sides were reconstructed and the resulting sequences are equal.
    Same,
    /// Both sides were reconstructed and the resulting sequences differ.
    Different,
    /// Nothing was reconstructed: the union window exceeds
    /// [`MAX_SEQUENCE_COMPARE_WINDOW`].
    WindowTooWide,
    /// Nothing was reconstructed: the union window is bounded and both sides
    /// convert to SPDI on a shared accession, but the provider could not serve
    /// the reference bases needed to compare them (#1989). The "want of
    /// reference data" decline — the sibling of [`Self::WindowTooWide`], and like
    /// it a decline that must **not** be read as a negative.
    ReferenceUnavailable,
    /// Nothing was reconstructed, for any other reason. The shapes it is meant
    /// to carry are representation-level — an edit SPDI cannot carry, a
    /// cross-accession pair, an overlapping allele that cannot be applied — and
    /// the checker deliberately reads this as a negative. It also carries a
    /// reference failure that happened *before* the fetch, in `edit_triples`'s
    /// discarded conversion error (`checker.rs:1114`) — that is #2056, and the
    /// enum-level note has the measurements.
    Declined,
}

/// The axis a description determines beyond the one it is written in.
#[derive(Debug, Clone, PartialEq, Eq)]
enum SecondAxis {
    /// There is none — the description already names its only axis.
    None,
    /// The genome, as `(contig, triples re-expressed on it)`.
    Genomic(String, Vec<SpdiVariant>),
    /// One is determined but could not be computed from what is available.
    NotComputed,
}

/// Apply two triple sets to a common reference window and compare the results.
///
/// `fetch` reads `[start, end)` on the shared accession, 0-based, returning
/// `None` on any provider error or short read.
fn compare_triples<F>(
    fetch: F,
    triples1: &[SpdiVariant],
    triples2: &[SpdiVariant],
) -> SequenceVerdict
where
    F: FnOnce(u64, u64) -> Option<String>,
{
    // The union window covering every edit of both variants. SPDI positions
    // are 0-based interbase; a triple spans `[position, position + del_len)`.
    let all = triples1.iter().chain(triples2.iter());
    let (mut win_start, mut win_end) = (u64::MAX, u64::MIN);
    for t in all {
        win_start = win_start.min(t.position);
        win_end = win_end.max(t.position + t.deletion.len() as u64);
    }
    if win_start > win_end {
        return SequenceVerdict::Declined;
    }
    if win_end - win_start > MAX_SEQUENCE_COMPARE_WINDOW {
        return SequenceVerdict::WindowTooWide;
    }

    let Some(reference) = fetch(win_start, win_end) else {
        // The window is bounded and both triple sets are in hand, but the
        // provider could not serve the bases. This is "want of reference data",
        // not a representation limit, so it is the `WindowTooWide` sibling and
        // must not be read as a negative (#1989).
        return SequenceVerdict::ReferenceUnavailable;
    };

    match (
        crate::spdi::apply::apply_triples(&reference, win_start, triples1),
        crate::spdi::apply::apply_triples(&reference, win_start, triples2),
    ) {
        // Case-insensitive, like the deletion check in `apply_triples`:
        // reference FASTAs are often soft-masked (repeats lower-cased), so
        // case carries no biological meaning here and must not split two
        // otherwise-identical resulting sequences.
        (Some(a), Some(b)) if a.eq_ignore_ascii_case(&b) => SequenceVerdict::Same,
        (Some(_), Some(_)) => SequenceVerdict::Different,
        _ => SequenceVerdict::Declined,
    }
}

/// Why `variant` has no computable denotation, or `None` when it has one.
///
/// **A closed list, keyed on the description alone** — never on what the
/// provider happens to hold. That distinction is the whole design: a missing
/// transcript is a gap in the *data*, and reporting it as "undecidable" would
/// turn every provider gap into a non-verdict and make the rung useless as a
/// census. What lands here is a description that denotes no sequence to anyone,
/// with any reference.
///
/// The list is a floor and is expected to grow as shapes are measured. It
/// currently covers:
///
/// * a **non-literal inserted payload** — `ins6`, `insN[33]`, `ins(10_20)`,
///   `insAluYb8`, `delins(?)`;
/// * `NaEdit::Unknown` — `c.?`, which asserts a change without naming one;
/// * the **null and unknown alleles** (`0`, `?`), and a **multi-molecule
///   allele** (trans / mosaic / chimeric / and-or / products / unknown phase),
///   which denotes no *single* sequence;
/// * a genomic position carrying an **offset** (`g.10+2del`) or a **special
///   marker** (`pter`, `qter`, `cen`), and an **uncertain or unknown position**
///   (`g.(100_150)del`, `g.?_100del`). The first two are the live holes
///   recorded as #1641 and #1643: `hgvs_to_spdi`'s position helpers read
///   `GenomePos::base` and drop both fields, so `g.10+2delC` and `g.10delC`
///   produce the *same* triple. Catching them here means the checker reports
///   "cannot tell" rather than a confident wrong answer, whichever way that
///   defect is fixed at the conversion.
///
/// Not covered yet: uncertain positions on the `c.`/`n.`/`r.` axes, which need
/// a location accessor generic over every variant kind that this crate does not
/// expose. Add them when a wrong answer is measured, with the pair that
/// measures it.
fn indeterminate_reason(variant: &HgvsVariant) -> Option<String> {
    match variant {
        HgvsVariant::NullAllele => {
            Some("`0` denotes no sequence: it asserts the absence of a product".to_string())
        }
        HgvsVariant::UnknownAllele => {
            Some("`?` denotes no sequence: it asserts that the variant is unknown".to_string())
        }
        HgvsVariant::Allele(allele) => {
            if allele.phase != AllelePhase::Cis {
                return Some(format!(
                    "{variant}: a {} allele names more than one molecule, so it denotes no \
                     single sequence to compare",
                    allele.phase
                ));
            }
            allele.variants.iter().find_map(indeterminate_reason)
        }
        HgvsVariant::Genome(g) => genome_indeterminate_reason(variant, &g.loc_edit.location)
            .or_else(|| edit_indeterminate_reason(variant, &g.loc_edit.edit)),
        HgvsVariant::Mt(m) => genome_indeterminate_reason(variant, &m.loc_edit.location)
            .or_else(|| edit_indeterminate_reason(variant, &m.loc_edit.edit)),
        HgvsVariant::Circular(o) => genome_indeterminate_reason(variant, &o.loc_edit.location)
            .or_else(|| edit_indeterminate_reason(variant, &o.loc_edit.edit)),
        HgvsVariant::Cds(c) => edit_indeterminate_reason(variant, &c.loc_edit.edit),
        HgvsVariant::Tx(n) => edit_indeterminate_reason(variant, &n.loc_edit.edit),
        HgvsVariant::Rna(r) => edit_indeterminate_reason(variant, &r.loc_edit.edit),
        _ => None,
    }
}

/// The genomic-position half of [`indeterminate_reason`].
fn genome_indeterminate_reason(variant: &HgvsVariant, loc: &GenomeInterval) -> Option<String> {
    let describe = |what: &str| format!("{variant}: {what}, so it denotes no definite sequence");
    for boundary in [&loc.start, &loc.end] {
        let mu = match boundary {
            UncertainBoundary::Single(mu) => mu,
            UncertainBoundary::Range { .. } => {
                return Some(describe("a boundary is an uncertain range"))
            }
        };
        {
            match mu {
                Mu::Unknown => return Some(describe("a position is unknown (`?`)")),
                Mu::Uncertain(_) => return Some(describe("a position is uncertain")),
                Mu::Certain(pos) => {
                    if pos.special.is_some() {
                        return Some(describe(
                            "a position is a chromosome landmark (`pter`/`qter`/`cen`) with no \
                             base coordinate",
                        ));
                    }
                    if pos.offset.is_some() {
                        return Some(describe(
                            "a genomic position carries an intronic offset, which has no meaning \
                             on a genomic reference",
                        ));
                    }
                }
            }
        }
    }
    None
}

/// The edit half of [`indeterminate_reason`].
///
/// The edit itself is a [`Mu`]: an uncertain edit (`g.(100del)`) states a
/// *predicted* change rather than an observed one, so it denotes no sequence
/// either.
fn edit_indeterminate_reason(variant: &HgvsVariant, edit: &Mu<NaEdit>) -> Option<String> {
    let edit = match edit {
        Mu::Certain(edit) => edit,
        Mu::Uncertain(_) => {
            return Some(format!(
                "{variant}: the edit is uncertain, so it states a predicted change rather than \
                 one that denotes a sequence"
            ))
        }
        Mu::Unknown => {
            return Some(format!(
                "{variant}: the edit is unknown, so it denotes no sequence"
            ))
        }
    };
    let payload = match edit {
        NaEdit::Unknown { .. } => {
            return Some(format!(
                "{variant}: `?` asserts a change without naming one, so it denotes no sequence"
            ))
        }
        NaEdit::Insertion { sequence }
        | NaEdit::BreakpointInsertion { sequence }
        | NaEdit::Delins { sequence, .. } => sequence,
        _ => return None,
    };
    inserted_is_indeterminate(payload).then(|| {
        format!(
            "{variant}: the inserted payload names no definite bases, so it denotes no sequence"
        )
    })
}

/// Whether an inserted payload names no definite bases.
///
/// Stated as a closed *positive* list of the non-literal shapes rather than as
/// "anything that is not a literal", because several of the remaining shapes —
/// a position range, an external reference — are perfectly resolvable against a
/// provider and must keep their existing verdicts.
fn inserted_is_indeterminate(inserted: &InsertedSequence) -> bool {
    match inserted {
        // A literal built from IUPAC ambiguity codes names a *family* of
        // sequences, not one. It reaches here because normalization expands
        // `insN[10]` into exactly that, so leaving it out would let the
        // expansion launder an indeterminate payload into a definite-looking
        // one.
        InsertedSequence::Literal(sequence) => {
            !sequence.bases().iter().all(|b| base_is_definite(*b))
        }
        InsertedSequence::Count(_)
        | InsertedSequence::Range(_, _)
        | InsertedSequence::Named(_)
        | InsertedSequence::Uncertain => true,
        InsertedSequence::Repeat { base, count } => {
            !matches!(count, RepeatCount::Exact(_)) || !base_is_definite(*base)
        }
        InsertedSequence::SequenceRepeat { count, .. } => !matches!(count, RepeatCount::Exact(_)),
        // A bracket is indeterminate exactly when one of its parts is. The
        // `Literal` arm is not optional: `ins[NN;A]` names the same family of
        // sequences as the bare `insNN` the arm above refuses, so omitting it
        // would let a bracket launder an indeterminate payload into a
        // definite-looking one. The remaining parts — a position range, a CDS
        // position range, an external reference — are resolvable against a
        // provider and keep their existing verdicts.
        InsertedSequence::Complex(parts) => parts.iter().any(|part| match part {
            InsertedPart::Literal(sequence) => {
                !sequence.bases().iter().all(|b| base_is_definite(*b))
            }
            InsertedPart::Repeat { base, count } => {
                !matches!(count, RepeatCount::Exact(_)) || !base_is_definite(*base)
            }
            _ => false,
        }),
        _ => false,
    }
}

/// Whether `base` names exactly one nucleotide. An IUPAC ambiguity code
/// (`N`, `R`, `Y`, …) does not, so a payload built from one denotes a family of
/// sequences rather than a sequence.
fn base_is_definite(base: crate::hgvs::edit::Base) -> bool {
    use crate::hgvs::edit::Base;
    matches!(base, Base::A | Base::C | Base::G | Base::T | Base::U)
}

/// Extract the variant part from an HGVS string (everything after the colon).
fn extract_variant_part(hgvs: &str) -> &str {
    if let Some(pos) = hgvs.find(':') {
        &hgvs[pos..]
    } else {
        hgvs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;
    use crate::reference::MockProvider;

    fn checker() -> EquivalenceChecker<MockProvider> {
        EquivalenceChecker::new(MockProvider::with_test_data())
    }

    #[test]
    fn test_identical_variants() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::Identical);
        assert!(result.is_equivalent());
    }

    #[test]
    fn test_different_variants() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.20A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
        assert!(!result.is_equivalent());
    }

    #[test]
    fn test_identical_compound_alleles_are_equivalent() {
        // Regression: `check_accession_version_difference` declines (returns None) for
        // alleles, and `check()` must fall through to the normalizer, which handles
        // alleles. Two identical compound alleles must compare as Identical, proving
        // there is no silent "unsupported" gap for alleles.
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.[10A>G;20C>T]").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.[10A>G;20C>T]").unwrap();
        assert!(matches!(v1, HgvsVariant::Allele(_)));

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::Identical);
        assert!(result.is_equivalent());
    }

    #[test]
    fn test_differing_alleles_do_not_error() {
        // A pair of non-identical alleles must still flow through `check()` without
        // erroring. The version-difference helper declines for alleles (returns None,
        // not an error), so these resolve via normalization rather than being reported
        // as an accession-version match — confirming the helper's None is
        // "not-applicable", never a silent unsupported failure.
        //
        // Both alleles use `NM_000088.3` (present in the MockProvider fixture) so the
        // test exercises real normalization and its correctness does not hinge on how
        // the normalizer treats an unknown accession. The members differ in
        // coordinates, so the alleles are not identical.
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.[10A>G;20C>T]").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.[11A>G;21C>T]").unwrap();
        assert!(matches!(v1, HgvsVariant::Allele(_)));
        assert!(matches!(v2, HgvsVariant::Allele(_)));

        // The contract under test: `check()` returns Ok for every allele pair.
        let result = checker.check(&v1, &v2).unwrap();
        // The differing member coordinates do not normalize equal, so they compare as
        // NotEquivalent — the point is the Ok, not the level.
        assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
    }

    #[test]
    fn test_equivalence_level_is_equivalent() {
        assert!(EquivalenceLevel::Identical.is_equivalent());
        assert!(EquivalenceLevel::NormalizedMatch.is_equivalent());
        assert!(EquivalenceLevel::SequenceMatch.is_equivalent());
        assert!(EquivalenceLevel::AccessionVersionDifference.is_equivalent());
        assert!(!EquivalenceLevel::NotEquivalent.is_equivalent());
    }

    #[test]
    fn test_equivalence_level_description() {
        assert!(!EquivalenceLevel::Identical.description().is_empty());
        assert!(!EquivalenceLevel::NormalizedMatch.description().is_empty());
        assert!(!EquivalenceLevel::SequenceMatch.description().is_empty());
        assert!(!EquivalenceLevel::AccessionVersionDifference
            .description()
            .is_empty());
        assert!(!EquivalenceLevel::NotEquivalent.description().is_empty());
    }

    #[test]
    fn test_equivalence_level_display() {
        let level = EquivalenceLevel::Identical;
        assert_eq!(level.to_string(), level.description());
    }

    #[test]
    fn test_all_equivalent_empty() {
        let checker = checker();
        assert!(checker.all_equivalent(&[]).unwrap());
    }

    #[test]
    fn test_all_equivalent_single() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        assert!(checker.all_equivalent(&[v1]).unwrap());
    }

    #[test]
    fn test_all_equivalent_same() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v3 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        assert!(checker.all_equivalent(&[v1, v2, v3]).unwrap());
    }

    #[test]
    fn test_all_equivalent_different() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.20A>G").unwrap();
        assert!(!checker.all_equivalent(&[v1, v2]).unwrap());
    }

    #[test]
    fn test_group_by_equivalence() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v3 = parse_hgvs("NM_000088.3:c.20A>G").unwrap();

        let groups = checker.group_by_equivalence(&[v1, v2, v3]).unwrap();
        assert_eq!(groups.len(), 2);
        // First group should have 2 variants (the identical ones)
        assert!(groups.iter().any(|g| g.len() == 2));
        // Second group should have 1 variant
        assert!(groups.iter().any(|g| g.len() == 1));
    }

    #[test]
    fn test_extract_variant_part() {
        assert_eq!(extract_variant_part("NM_000088.3:c.10A>G"), ":c.10A>G");
        assert_eq!(extract_variant_part("no_colon"), "no_colon");
    }

    #[test]
    fn test_equivalence_result_with_note() {
        let result = EquivalenceResult::new(EquivalenceLevel::NormalizedMatch)
            .with_note("Test note 1")
            .with_note("Test note 2");

        assert_eq!(result.notes.len(), 2);
        assert_eq!(result.notes[0], "Test note 1");
        assert_eq!(result.notes[1], "Test note 2");
    }

    #[test]
    fn test_substitution_at_different_positions() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.11A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
    }

    #[test]
    fn test_genomic_variants_identical() {
        let checker = checker();
        let v1 = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        let v2 = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::Identical);
    }

    #[test]
    fn test_genomic_variants_different() {
        // Two genuinely different substitutions must compare `NotEquivalent` —
        // but only once the sequence rung can actually reconstruct their window.
        // The reference is served here for exactly that reason: without it, the
        // comparison never runs and the honest verdict is `Indeterminate`, not a
        // decided negative (#1989). See
        // `issue_1989_declined_is_indeterminate` for the unserved case.
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_000001.11", "ACGTACGTACGT");
        let checker = EquivalenceChecker::new(provider);
        let v1 = parse_hgvs("NC_000001.11:g.1A>G").unwrap();
        let v2 = parse_hgvs("NC_000001.11:g.2C>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
    }

    #[test]
    fn test_protein_variants_identical() {
        let checker = checker();
        let v1 = parse_hgvs("NP_000079.2:p.Val600Glu").unwrap();
        let v2 = parse_hgvs("NP_000079.2:p.Val600Glu").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::Identical);
    }

    #[test]
    fn test_deletion_variants() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10del").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.10del").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::Identical);
    }

    // =========================================================================
    // P2: Cross-accession version tests
    // =========================================================================

    #[test]
    fn test_accession_version_difference_same_variant() {
        let checker = checker();
        // Same position and change, different accession versions
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.4:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
        assert!(result.is_equivalent());
        assert!(result.notes.iter().any(|n| n.contains("version")));
    }

    #[test]
    fn test_accession_version_difference_detected_before_normalize() {
        let checker = checker();
        // Versions differ, should detect before normalization
        let v1 = parse_hgvs("NC_000001.10:g.12345A>G").unwrap();
        let v2 = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
    }

    #[test]
    fn test_accession_version_difference_with_deletion() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10_12del").unwrap();
        let v2 = parse_hgvs("NM_000088.4:c.10_12del").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
    }

    #[test]
    fn test_accession_version_difference_with_insertion() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10_11insATG").unwrap();
        let v2 = parse_hgvs("NM_000088.4:c.10_11insATG").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
    }

    #[test]
    fn test_accession_version_different_variant_not_equivalent() {
        let checker = checker();
        // Different versions AND different variants = not equivalent
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.4:c.20A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
        assert!(!result.is_equivalent());
    }

    #[test]
    fn test_different_accessions_entirely_not_equivalent() {
        let checker = checker();
        // Different base accessions (not just version)
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_001234.1:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
    }

    #[test]
    fn test_protein_accession_version_difference() {
        let checker = checker();
        let v1 = parse_hgvs("NP_000079.1:p.Val600Glu").unwrap();
        let v2 = parse_hgvs("NP_000079.2:p.Val600Glu").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
    }

    #[test]
    fn test_no_version_vs_versioned() {
        let checker = checker();
        // One without version, one with
        let v1 = parse_hgvs("NM_000088:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        // Should still detect as version difference
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
    }

    #[test]
    fn test_version_difference_result_has_normalized_forms() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.4:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert!(result.normalized_first.is_some());
        assert!(result.normalized_second.is_some());
    }

    #[test]
    fn test_genomic_version_difference_grch37_vs_grch38() {
        let checker = checker();
        // GRCh37 (NC_000001.10) vs GRCh38 (NC_000001.11)
        let v1 = parse_hgvs("NC_000001.10:g.100000A>G").unwrap();
        let v2 = parse_hgvs("NC_000001.11:g.100000A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
        // Note: actual position may differ between builds, but this tests the detection
    }

    #[test]
    fn test_version_difference_note_contains_versions() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.5:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        // Notes should mention the specific versions
        let has_version_note = result
            .notes
            .iter()
            .any(|n| n.contains("3") && n.contains("5"));
        assert!(has_version_note, "Notes should contain version numbers");
    }

    #[test]
    fn test_multiple_variants_with_version_mix() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v3 = parse_hgvs("NM_000088.4:c.10A>G").unwrap();

        // all_equivalent should return true since they're all equivalent
        assert!(checker
            .all_equivalent(&[v1.clone(), v2.clone(), v3.clone()])
            .unwrap());

        // group_by_equivalence should put them all in one group
        let groups = checker.group_by_equivalence(&[v1, v2, v3]).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 3);
    }

    #[test]
    fn test_lrg_accession_version_difference() {
        let checker = checker();
        // LRG transcripts
        let v1 = parse_hgvs("LRG_1t1:c.10A>G").unwrap();
        let v2 = parse_hgvs("LRG_1t2:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        // These are different transcripts (t1 vs t2), not versions
        // Should be NotEquivalent
        assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
    }

    #[test]
    fn test_ensembl_accession_version_difference() {
        let checker = checker();
        let v1 = parse_hgvs("ENST00000123456.1:c.10A>G").unwrap();
        let v2 = parse_hgvs("ENST00000123456.2:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
    }

    // =========================================================================
    // P3: Equivalence grouping performance tests
    // =========================================================================

    #[test]
    fn test_grouping_100_identical_variants() {
        let checker = checker();

        // Generate 100 identical variants
        let variants: Vec<HgvsVariant> = (0..100)
            .map(|_| parse_hgvs("NM_000088.3:c.10A>G").unwrap())
            .collect();

        let groups = checker.group_by_equivalence(&variants).unwrap();
        // All should be in one group
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 100);
    }

    #[test]
    fn test_grouping_100_unique_variants() {
        let checker = checker();

        // Generate 100 unique variants (different positions)
        let variants: Vec<HgvsVariant> = (1..=100)
            .map(|i| parse_hgvs(&format!("NM_000088.3:c.{}A>G", i)).unwrap())
            .collect();

        let groups = checker.group_by_equivalence(&variants).unwrap();
        // Each should be in its own group
        assert_eq!(groups.len(), 100);
        for group in &groups {
            assert_eq!(group.len(), 1);
        }
    }

    #[test]
    fn test_grouping_mixed_identical_and_unique() {
        let checker = checker();

        // 50 copies of variant A, 50 copies of variant B
        let mut variants: Vec<HgvsVariant> = Vec::with_capacity(100);

        for _ in 0..50 {
            variants.push(parse_hgvs("NM_000088.3:c.10A>G").unwrap());
        }
        for _ in 0..50 {
            variants.push(parse_hgvs("NM_000088.3:c.20C>T").unwrap());
        }

        let groups = checker.group_by_equivalence(&variants).unwrap();
        assert_eq!(groups.len(), 2);

        // Order may vary, but sizes should be 50 each
        let sizes: Vec<usize> = groups.iter().map(|g| g.len()).collect();
        assert!(sizes.contains(&50));
        assert_eq!(sizes.iter().sum::<usize>(), 100);
    }

    #[test]
    fn test_grouping_interleaved_variants() {
        let checker = checker();

        // Interleave 3 different variants
        let mut variants: Vec<HgvsVariant> = Vec::with_capacity(60);

        for i in 0..20 {
            variants.push(parse_hgvs(&format!("NM_000088.3:c.{}A>G", 10)).unwrap());
            variants.push(parse_hgvs(&format!("NM_000088.3:c.{}A>G", 20)).unwrap());
            variants.push(parse_hgvs(&format!("NM_000088.3:c.{}A>G", 30)).unwrap());
            let _ = i; // silence unused warning
        }

        let groups = checker.group_by_equivalence(&variants).unwrap();
        assert_eq!(groups.len(), 3);

        // Each group should have 20 variants
        for group in &groups {
            assert_eq!(group.len(), 20);
        }
    }

    #[test]
    fn test_grouping_accession_version_differences() {
        let checker = checker();

        // Mix of version 3 and version 4 of same variant
        let mut variants: Vec<HgvsVariant> = Vec::with_capacity(40);

        for _ in 0..20 {
            variants.push(parse_hgvs("NM_000088.3:c.10A>G").unwrap());
        }
        for _ in 0..20 {
            variants.push(parse_hgvs("NM_000088.4:c.10A>G").unwrap());
        }

        let groups = checker.group_by_equivalence(&variants).unwrap();
        // They're considered equivalent (AccessionVersionDifference), so all in one group
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 40);
    }

    #[test]
    fn test_grouping_diverse_variant_types() {
        let checker = checker();

        // Different variant types at different positions
        let variant_strs = [
            "NM_000088.3:c.10A>G",       // substitution
            "NM_000088.3:c.20del",       // deletion
            "NM_000088.3:c.30_31insATG", // insertion
            "NM_000088.3:c.40dup",       // duplication
            "NC_000001.11:g.12345A>G",   // genomic
        ];

        let variants: Vec<HgvsVariant> = variant_strs
            .iter()
            .map(|s| parse_hgvs(s).unwrap())
            .collect();

        let groups = checker.group_by_equivalence(&variants).unwrap();
        // All different, so 5 groups
        assert_eq!(groups.len(), 5);
    }

    #[test]
    fn test_grouping_empty_input() {
        let checker = checker();
        let variants: Vec<HgvsVariant> = vec![];

        let groups = checker.group_by_equivalence(&variants).unwrap();
        assert!(groups.is_empty());
    }

    #[test]
    fn test_grouping_single_variant() {
        let checker = checker();
        let variants = vec![parse_hgvs("NM_000088.3:c.10A>G").unwrap()];

        let groups = checker.group_by_equivalence(&variants).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 1);
    }

    #[test]
    fn test_grouping_performance_200_variants() {
        use std::time::Instant;

        let checker = checker();

        // Generate 200 variants: 10 groups of 20 each
        let mut variants: Vec<HgvsVariant> = Vec::with_capacity(200);

        for pos in (1..=10).map(|x| x * 10) {
            for _ in 0..20 {
                variants.push(parse_hgvs(&format!("NM_000088.3:c.{}A>G", pos)).unwrap());
            }
        }

        let start = Instant::now();
        let groups = checker.group_by_equivalence(&variants).unwrap();
        let duration = start.elapsed();

        assert_eq!(groups.len(), 10);
        for group in &groups {
            assert_eq!(group.len(), 20);
        }

        // Should complete in reasonable time (< 5 seconds)
        assert!(
            duration.as_secs() < 5,
            "Grouping took too long: {:?}",
            duration
        );
    }

    #[test]
    fn test_all_equivalent_performance_100() {
        use std::time::Instant;

        let checker = checker();

        // 100 identical variants
        let variants: Vec<HgvsVariant> = (0..100)
            .map(|_| parse_hgvs("NM_000088.3:c.10A>G").unwrap())
            .collect();

        let start = Instant::now();
        let result = checker.all_equivalent(&variants).unwrap();
        let duration = start.elapsed();

        assert!(result);
        assert!(
            duration.as_secs() < 2,
            "all_equivalent took too long: {:?}",
            duration
        );
    }

    #[test]
    fn test_all_equivalent_early_exit_on_non_equivalent() {
        use std::time::Instant;

        let checker = checker();

        // First 2 are different, rest are same
        let mut variants: Vec<HgvsVariant> = Vec::with_capacity(100);
        variants.push(parse_hgvs("NM_000088.3:c.10A>G").unwrap());
        variants.push(parse_hgvs("NM_000088.3:c.20A>G").unwrap()); // different!

        for _ in 0..98 {
            variants.push(parse_hgvs("NM_000088.3:c.10A>G").unwrap());
        }

        let start = Instant::now();
        let result = checker.all_equivalent(&variants).unwrap();
        let duration = start.elapsed();

        assert!(!result);
        // Should exit early and be very fast
        assert!(
            duration.as_millis() < 1000,
            "all_equivalent should exit early"
        );
    }

    #[test]
    fn test_pairwise_check_performance() {
        use std::time::Instant;

        let checker = checker();

        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.20A>G").unwrap();

        // Run 1000 pairwise checks
        let start = Instant::now();
        for _ in 0..1000 {
            let _ = checker.check(&v1, &v2).unwrap();
        }
        let duration = start.elapsed();

        // 1000 checks should complete in < 2 seconds
        assert!(
            duration.as_secs() < 2,
            "Pairwise checks too slow: {:?}",
            duration
        );
    }

    #[test]
    fn test_grouping_stability_same_results() {
        let checker = checker();

        // Same variants in same order should produce same groups
        let variants: Vec<HgvsVariant> = vec![
            parse_hgvs("NM_000088.3:c.10A>G").unwrap(),
            parse_hgvs("NM_000088.3:c.20A>G").unwrap(),
            parse_hgvs("NM_000088.3:c.10A>G").unwrap(),
            parse_hgvs("NM_000088.3:c.30A>G").unwrap(),
            parse_hgvs("NM_000088.3:c.20A>G").unwrap(),
        ];

        let groups1 = checker.group_by_equivalence(&variants).unwrap();
        let groups2 = checker.group_by_equivalence(&variants).unwrap();

        // Same number of groups
        assert_eq!(groups1.len(), groups2.len());

        // Groups have same sizes
        let mut sizes1: Vec<_> = groups1.iter().map(|g| g.len()).collect();
        let mut sizes2: Vec<_> = groups2.iter().map(|g| g.len()).collect();
        sizes1.sort();
        sizes2.sort();
        assert_eq!(sizes1, sizes2);
    }

    #[test]
    fn test_grouping_preserves_variant_order_in_groups() {
        let checker = checker();

        // Variants added in specific order
        let variants: Vec<HgvsVariant> = vec![
            parse_hgvs("NM_000088.3:c.10A>G").unwrap(), // Group A, first
            parse_hgvs("NM_000088.3:c.10A>G").unwrap(), // Group A, second
            parse_hgvs("NM_000088.3:c.10A>G").unwrap(), // Group A, third
        ];

        let groups = checker.group_by_equivalence(&variants).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 3);

        // All variants in the group should be identical
        for v in &groups[0] {
            assert_eq!(v.to_string(), "NM_000088.3:c.10A>G");
        }
    }

    /// A bracketed payload must not launder what a bare literal refuses.
    /// `insNN` is indeterminate because the literal names a family of sequences
    /// rather than one; `ins[NN;A]` names exactly the same family, and
    /// `InsertedSequence::Complex` must say so.
    ///
    /// Asserted on the **raw** parsed description on purpose. The checker asks
    /// this question of `v1`/`v2` as well as of their normalized forms,
    /// precisely so that a verdict about what the user wrote can see what the
    /// user wrote — and the normalizer happens to flatten a bracket of literals
    /// into one `Literal` (`ins[NN;A]` -> `insNNA`), so the normalized side
    /// masks this at the end-to-end level. The mask is incidental to how
    /// bracketed payloads are expanded today and is not a property this
    /// predicate may lean on.
    #[test]
    fn a_bracketed_payload_cannot_launder_an_ambiguity_code() {
        for description in [
            "NC_000001.11:g.1001_1002ins[NN;A]",
            "NC_000001.11:g.1001_1002ins[AA;N]",
            "NC_000001.11:g.1001_1002ins[R;A]",
            "NC_000001.11:g.1001_1002delins[NN;A]",
        ] {
            let variant = parse_hgvs(description).unwrap();
            assert!(
                indeterminate_reason(&variant).is_some(),
                "{description}: an ambiguity code inside a bracket names no definite bases"
            );
        }

        // The negative control: a bracket of definite literals still denotes
        // exactly one sequence, so it must keep its existing verdict.
        let definite = parse_hgvs("NC_000001.11:g.1001_1002ins[AA;C]").unwrap();
        assert!(indeterminate_reason(&definite).is_none());
    }

    /// `transcript_triples_to_genomic` turns a 1-based genomic position into
    /// SPDI's 0-based interbase one by subtracting 1. Nothing between a
    /// provider and that arithmetic asserts `Exon::genomic_start >= 1` — the
    /// cdot path takes `e[0]` verbatim and the JSON loader takes whatever the
    /// file holds — so a transcript whose exon is placed at genomic 0 reaches
    /// it. On `u64` that subtraction panics under `debug-assertions` and wraps
    /// to `u64::MAX` in a release build, which is the worse of the two. The
    /// deleted-span arm must decline the way the insertion arm above it
    /// declines `triple.position == 0`.
    #[test]
    fn a_transcript_placed_at_genomic_zero_is_declined_rather_than_underflowed() {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Transcript};

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript::new(
            "NM_ZEROBASED.1".to_string(),
            Some("ZERO".to_string()),
            Strand::Plus,
            "ACGTACGTAC".to_string(),
            Some(1),
            Some(9),
            // A 0-based placement: transcript base 1 sits at genomic 0.
            vec![Exon::with_genomic(1, 1, 10, 0, 9)],
            Some("chr_zero".to_string()),
            Some(0),
            Some(9),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        ));
        let checker = EquivalenceChecker::new(provider);

        // A one-base deletion at transcript base 1: `first == last == 0`.
        let triples = vec![SpdiVariant::new(
            "NM_ZEROBASED.1".to_string(),
            0,
            "A".to_string(),
            String::new(),
        )];
        assert!(
            checker
                .transcript_triples_to_genomic("NM_ZEROBASED.1", &triples)
                .is_none(),
            "a genomic position with no interbase point to its left must be declined"
        );
    }

    /// `compare_denotations` deliberately has no `Identical` rung: it answers
    /// what the *comparison* found, so an identical pair is re-derived and
    /// compared like any other. Without this, a caller could not tell "the
    /// bases agree" from "the strings agree", which is the distinction the
    /// entry point exists for.
    #[test]
    fn compare_denotations_re_derives_an_identical_pair_instead_of_short_circuiting() {
        let checker = checker();
        let v = parse_hgvs("NM_000088.3:c.10del").unwrap();

        assert_eq!(
            checker.check(&v, &v).unwrap().level,
            EquivalenceLevel::Identical
        );
        // The exact rung, not merely "not `Identical`". `assert_ne!` alone is
        // satisfied by `NotEquivalent` — a *decided negative* about a
        // description compared with itself — and by `Indeterminate`, which
        // would mean the re-derivation this test is named for never happened.
        //
        // `SequenceMatch` rather than `CrossAxisSequenceMatch` is the measured
        // answer and is correct here: `MockProvider::with_test_data` serves no
        // genomic sequence for `NM_000088.3`, so the second determined axis
        // cannot be computed and the verdict stops one rung short. The note
        // says so in as many words — "a determined axis could not be computed"
        // — so a future fixture that *does* serve that axis will fail here and
        // be re-pinned deliberately rather than drift.
        let denoted = checker.compare_denotations(&v, &v).unwrap();
        assert_eq!(denoted.level, EquivalenceLevel::SequenceMatch);
        assert!(denoted.level.is_decided());
    }

    /// The accession-version rung needs no normalizer, so it is answered here
    /// too. Dropping it would report a decided `NotEquivalent` for one variant
    /// spelled on two versions of its accession — a positive claim no rung
    /// examined, since apply-equality is not defined across two references.
    #[test]
    fn compare_denotations_still_answers_an_accession_version_difference() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.4:c.10A>G").unwrap();

        assert_eq!(
            checker.compare_denotations(&v1, &v2).unwrap().level,
            EquivalenceLevel::AccessionVersionDifference
        );
    }
}
