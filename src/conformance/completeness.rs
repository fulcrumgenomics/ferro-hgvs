//! Completeness accounting for corpus and fixture generators (#1550).
//!
//! Every generator under `examples/` walks a population — corpus rows, spec
//! strings, transcripts, reference accessions — and writes an artifact derived
//! from it. They share one failure mode: **a fallible step whose failure is
//! representable as a legitimate value**. `spans_of(..).unwrap_or_default()`
//! yields an empty vec that classifies as "no gaps"; `else { continue }` skips
//! exactly the record that had a problem; a discarded `Result` reports nothing
//! at all. The dropped population is never counted, so a partial run and a clean
//! run produce indistinguishable output and the only thing separating them is
//! whether a reviewer read the right five lines.
//!
//! A grep cannot find these. Measured on `main`, the idioms above occur at ~145
//! sites across `examples/` and `src/conformance/`, overwhelmingly legitimately;
//! what makes one a defect is not the idiom but that its failure flows into a
//! value that is then *counted*, which is a composition no lint can see. So the
//! mechanical half stays narrow (`tests/it/generator_completeness.rs`) and the
//! semantic half is carried by a type.
//!
//! [`CaptureLedger`] is that type. It is
//! `examples/extract_split_member_separation_windows.rs`'s refusal — *"Refusing
//! to write a fixture built from a partial pass"* — promoted out of one file:
//!
//! ```
//! use ferro_hgvs::conformance::completeness::CaptureLedger;
//!
//! let mut ledger = CaptureLedger::new("transcript windows");
//! for accession in ["NM_000088.3", "NM_003002.2"] {
//!     let resolved: Result<&str, String> = Ok(accession);
//!     ledger.record(accession, resolved);
//! }
//! // `finish` is the gate: it refuses on any shortfall, and the counts it
//! // returns are serializable so the artifact carries its own claim.
//! let counts = ledger.finish().expect("a complete pass");
//! assert_eq!(counts.attempted, 2);
//! assert_eq!(counts.dropped, 0);
//! ```
//!
//! Three properties are deliberate:
//!
//! - **`finish` refuses rather than reports.** A printed count is read by a
//!   human who already believes the run worked. `finish` returns a `Result` the
//!   caller must handle before it can reach its `fs::write`.
//! - **A shortfall is waived only by naming it.** [`Allowance`] takes a cap and
//!   a justification, and is constructed at the call site, so an expected drop
//!   is a declaration in the source rather than an absence.
//! - **An empty pass is a shortfall too.** A generator that attempted nothing
//!   writes an artifact whose emptiness reads as "there was nothing to find" —
//!   the corpus-zero hazard, where a structural blindness is indistinguishable
//!   from a safe result. Waiving it needs
//!   [`Allowance::permitting_an_empty_pass`].
//!
//! # What the ledger accounts for, and what it does not
//!
//! **It accounts for the step you route through it — not for the artifact.**
//! `succeeded` counts records that got past *one* fallible step. It is never
//! compared against the number of rows that actually reach the file, and the
//! type has no way to make that comparison: it never sees the artifact. So a
//! generator can use this exactly as prescribed and still write a partial
//! artifact carrying a spotless claim:
//!
//! ```ignore
//! let mut ledger = CaptureLedger::new("rows");
//! for row in rows {
//!     let Some(v) = ledger.record(row.id(), parse(row)) else { continue };
//!     let Ok(span) = spans_of(v) else { continue };  // second step, unaccounted
//!     artifact.push(span);
//! }
//! let counts = ledger.finish()?;   // {"attempted":10,"succeeded":10,"dropped":0}
//! ```                              // …stamped onto an artifact holding 5 rows
//!
//! That is the *exact* shape #1550 cites for `dump_confluence_divergences.rs`
//! (`spans_of(..)` failing downstream of a successful capture), so it is worth
//! being blunt about: **route the last fallible step before the artifact, not
//! the first.** If a pass has several fallible steps, either route them all
//! through one ledger — `record_drop` at each — or give each its own ledger and
//! `finish` each. A claim scoped to one loop is still worth stamping; it is just
//! a claim about that loop, and an artifact that says so is honest where one
//! that implies more is not.
//!
//! The two mechanical guards in `tests/it/generator_completeness.rs` are
//! narrower still: substring heuristics over generator source text, which catch
//! the common case and prove nothing. See that file's module docs.

use std::collections::BTreeMap;
use std::fmt;

use serde::{Deserialize, Serialize};

/// One dropped record: what was being processed, and why it did not make it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DroppedRecord {
    /// The thing that was being captured — an accession, an input string, a row
    /// id. Named rather than counted so a failure report says *which*.
    pub subject: String,
    /// Why it was dropped. Grouped verbatim in [`CaptureCounts::dropped_by_reason`],
    /// so prefer a small closed vocabulary over an interpolated error string.
    ///
    /// [`CaptureLedger::record`] does not, and cannot — it has only the error to
    /// go on. See its docs for when to reach for
    /// [`CaptureLedger::record_drop`] instead.
    pub reason: String,
}

/// A declared allowance for a shortfall, constructed at the call site.
///
/// The justification is required rather than optional: an allowance without one
/// is indistinguishable from having never considered the question, which is the
/// state this whole module exists to make impossible.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Allowance {
    max_dropped: usize,
    empty_pass_permitted: bool,
    justification: String,
}

impl Allowance {
    /// Permit up to `max_dropped` drops, for the stated reason.
    ///
    /// `Allowance::at_most(0, ..)` is meaningful and is *not* the same as
    /// [`CaptureLedger::finish`]: it still refuses a drop, but it records in the
    /// artifact that the question was asked and answered.
    ///
    /// # Panics
    ///
    /// If `justification` is empty or only whitespace. "Required" was previously
    /// presence rather than content, so `Allowance::at_most(0, "")` was accepted —
    /// an allowance carrying no reason is exactly the "indistinguishable from
    /// having never considered the question" state the type exists to prevent.
    pub fn at_most(max_dropped: usize, justification: impl Into<String>) -> Self {
        let justification = justification.into();
        assert!(
            !justification.trim().is_empty(),
            "an Allowance needs a justification: an empty one is indistinguishable \
             from having never considered the question",
        );
        Self {
            max_dropped,
            empty_pass_permitted: false,
            justification,
        }
    }

    /// Also permit a pass that attempted nothing.
    ///
    /// Needed only by generators whose population is legitimately optional (an
    /// absent input file, a corpus filtered to empty by a CLI flag). For anything
    /// that ships a committed artifact, an empty pass is the corpus-zero hazard
    /// and should stay a refusal.
    #[must_use]
    pub fn permitting_an_empty_pass(mut self) -> Self {
        self.empty_pass_permitted = true;
        self
    }
}

/// The allowance as stamped into an artifact, so a reader of the fixture sees
/// the waiver rather than having to find the generator that granted it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AllowanceNote {
    /// The declared ceiling on drops.
    pub max_dropped: usize,
    /// Whether a pass that attempted nothing was permitted.
    pub empty_pass_permitted: bool,
    /// Why the shortfall is acceptable.
    pub justification: String,
}

/// The completeness claim of one capture pass, serializable so a generator can
/// stamp it into the artifact it writes.
///
/// A fixture that carries its counts can be asserted on by the tests that
/// consume it, which is the difference between a completeness claim and a
/// completeness *hope*: the claim survives the generator run and travels with
/// the data.
///
/// # There is deliberately no `Default`
///
/// `Default` was derived here in the first revision, and
/// `ledger.finish().unwrap_or_default()` therefore compiled clean — fabricating
/// `{"attempted":0,"succeeded":0,"dropped":0}`, which serializes as a *tidier*
/// claim than a real one and reads as clean to anyone diffing a fixture. That is
/// `unwrap_or_default()`, the archetypal defect this module's opening paragraph
/// names, applied to the type meant to kill it. Without the derive, that line no
/// longer compiles and the caller has to say what it wants. Do not re-add it:
/// there is no meaningful "default completeness claim", and the empty one is the
/// most misleading value the type can hold.
///
/// ```compile_fail
/// use ferro_hgvs::conformance::completeness::CaptureLedger;
///
/// let mut ledger = CaptureLedger::new("rows");
/// ledger.record_success();
/// // `CaptureCounts: Default` would make this compile and discard the refusal.
/// let _counts = ledger.finish().unwrap_or_default();
/// ```
///
/// `let _ = ledger.finish();` is still silent — that one is unavoidable in Rust,
/// and `#[must_use]` on `Shortfall`'s `Result` is what argues against it.
///
/// # The arithmetic invariant does not survive deserialization
///
/// `attempted == succeeded + dropped` holds by construction *through
/// [`CaptureLedger`]* — `attempted()` is derived, so no caller can break it. It
/// is **not** re-checked by `Deserialize`, and deserialization is exactly how a
/// consuming test reads a stamped claim. Use [`CaptureCounts::is_self_consistent`]
/// if you are reading counts you did not just produce.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CaptureCounts {
    /// Records the pass tried to capture. Equal to `succeeded + dropped` by
    /// construction — the ledger has no way to increment one without the other.
    pub attempted: usize,
    /// Records that made it into the artifact.
    pub succeeded: usize,
    /// Records that did not.
    pub dropped: usize,
    /// Drops grouped by reason, so an artifact records the *shape* of a waived
    /// shortfall and not merely its size.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub dropped_by_reason: BTreeMap<String, usize>,
    /// The allowance under which a shortfall was accepted, when there was one.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub allowance: Option<AllowanceNote>,
}

impl CaptureCounts {
    /// Whether the pass captured everything it attempted, and attempted
    /// something.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.dropped == 0 && self.attempted > 0
    }

    /// Whether the arithmetic holds: `attempted == succeeded + dropped`, and the
    /// per-reason tally sums to `dropped`.
    ///
    /// Always true of counts produced by [`CaptureLedger`]. Worth asserting on
    /// counts read back out of an artifact, where nothing enforced it — a
    /// hand-edited or truncated claim is otherwise indistinguishable from a real
    /// one.
    ///
    /// Both sums are checked, because the input this method exists for is the
    /// input nothing validated: a claim holding two near-`usize::MAX` counts
    /// panics a debug build on `succeeded + dropped` and wraps a release one,
    /// and a wrapped sum can land back on `attempted` — reporting a claim as
    /// consistent *because* it overflowed. Overflow is treated as inconsistent.
    #[must_use]
    pub fn is_self_consistent(&self) -> bool {
        let total = self.succeeded.checked_add(self.dropped);
        let by_reason = self
            .dropped_by_reason
            .values()
            .try_fold(0usize, |sum, count| sum.checked_add(*count));
        total == Some(self.attempted) && by_reason == Some(self.dropped)
    }
}

/// The detail of a [`Shortfall::Dropped`].
///
/// Split out and boxed because it is the `Err` half of every `finish`, and an
/// inline variant this size makes each of them a large-`Result` clippy denial.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DroppedShortfall {
    /// What was being captured, for the message.
    pub subject: String,
    /// The counts as they stood at `finish`.
    pub counts: CaptureCounts,
    /// Every dropped record, in the order they were recorded.
    pub records: Vec<DroppedRecord>,
    /// The ceiling that was exceeded, when an allowance was declared.
    pub allowed: Option<usize>,
}

/// Why [`CaptureLedger::finish`] refused.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Shortfall {
    /// The pass attempted no records at all.
    NothingAttempted {
        /// What was being captured, for the message.
        subject: String,
    },
    /// The pass dropped records.
    Dropped(Box<DroppedShortfall>),
}

impl Shortfall {
    /// What was being captured.
    #[must_use]
    pub fn subject(&self) -> &str {
        match self {
            Self::NothingAttempted { subject } => subject,
            Self::Dropped(detail) => &detail.subject,
        }
    }

    /// The dropped records, empty for [`Shortfall::NothingAttempted`].
    #[must_use]
    pub fn dropped_records(&self) -> &[DroppedRecord] {
        match self {
            Self::NothingAttempted { .. } => &[],
            Self::Dropped(detail) => &detail.records,
        }
    }
}

/// How many dropped subjects a message names before summarising the rest. Enough
/// to see the shape of a failure without a 20,000-line panic.
const MAX_NAMED_SUBJECTS: usize = 10;

/// How many distinct drop reasons a message lists before summarising the rest.
///
/// This cap is not decoration. [`CaptureLedger::record`] — the documented
/// adoption path — sets `reason` to `error.to_string()`, and for an
/// `io::Error`/`anyhow::Error` that is typically unique per record, so
/// `dropped_by_reason` degenerates to one bucket per drop. Printing that map in
/// full defeated [`MAX_NAMED_SUBJECTS`] entirely: measured at 2000 drops, the
/// refusal message was 2,014 lines and 71 KB — precisely the "20,000-line panic"
/// the subject cap exists to prevent (#1551 review). Both halves are capped now,
/// and the totals are stated exactly either way, because a summary that hid the
/// count would reintroduce the ambiguity this module removes.
const MAX_NAMED_REASONS: usize = 10;

impl fmt::Display for Shortfall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NothingAttempted { subject } => write!(
                f,
                "{subject}: the pass attempted 0 records, so an artifact written from it \
                 would claim an empty population rather than report that it found none. \
                 Refusing to write a fixture built from a partial pass. If an empty pass \
                 is legitimate here, declare it with \
                 `Allowance::at_most(0, ..).permitting_an_empty_pass()`."
            ),
            Self::Dropped(detail) => {
                let DroppedShortfall {
                    subject,
                    counts,
                    records,
                    allowed,
                } = detail.as_ref();
                write!(
                    f,
                    "{subject}: {} of {} records were dropped",
                    counts.dropped, counts.attempted
                )?;
                match allowed {
                    Some(max) => write!(f, ", exceeding the declared allowance of {max}")?,
                    None => write!(f, " and no allowance was declared")?,
                }
                write!(
                    f,
                    ". Refusing to write a fixture built from a partial pass, because a \
                     partial pass and a clean one produce indistinguishable output.\n  \
                     by reason ({} distinct):",
                    counts.dropped_by_reason.len()
                )?;
                for (reason, count) in counts.dropped_by_reason.iter().take(MAX_NAMED_REASONS) {
                    write!(f, "\n    {reason}: {count}")?;
                }
                if counts.dropped_by_reason.len() > MAX_NAMED_REASONS {
                    let hidden: usize = counts
                        .dropped_by_reason
                        .values()
                        .skip(MAX_NAMED_REASONS)
                        .sum();
                    write!(
                        f,
                        "\n    … and {} more reasons covering {hidden} drops",
                        counts.dropped_by_reason.len() - MAX_NAMED_REASONS,
                    )?;
                }
                write!(f, "\n  dropped:")?;
                for record in records.iter().take(MAX_NAMED_SUBJECTS) {
                    write!(f, "\n    {} ({})", record.subject, record.reason)?;
                }
                if records.len() > MAX_NAMED_SUBJECTS {
                    write!(f, "\n    … and {} more", records.len() - MAX_NAMED_SUBJECTS)?;
                }
                Ok(())
            }
        }
    }
}

impl std::error::Error for Shortfall {}

/// Tracks what a capture pass attempted, captured and dropped, and refuses to
/// certify a pass that fell short.
///
/// See the [module docs](self) for why this exists and what it deliberately does
/// not try to do.
#[derive(Debug, Clone)]
pub struct CaptureLedger {
    subject: String,
    succeeded: usize,
    dropped: Vec<DroppedRecord>,
}

impl CaptureLedger {
    /// Start a ledger for `subject` — a short noun phrase naming the population,
    /// used verbatim in failure messages (`"transcript windows"`, `"spec rows"`).
    pub fn new(subject: impl Into<String>) -> Self {
        Self {
            subject: subject.into(),
            succeeded: 0,
            dropped: Vec::new(),
        }
    }

    /// Record a captured record.
    pub fn record_success(&mut self) {
        self.succeeded += 1;
    }

    /// Record a dropped record and why.
    ///
    /// `reason` is grouped verbatim, so keep it a short closed-vocabulary phrase
    /// (`"no stored sequence"`) and put the varying part in `subject`.
    pub fn record_drop(&mut self, subject: impl Into<String>, reason: impl Into<String>) {
        self.dropped.push(DroppedRecord {
            subject: subject.into(),
            reason: reason.into(),
        });
    }

    /// Record the outcome of a fallible step, returning its value on success.
    ///
    /// This is the adoption path for the idioms this module exists to replace:
    /// `let Ok(v) = f() else { continue };` becomes
    /// `let Some(v) = ledger.record(id, f()) else { continue };`, which drops the
    /// same record but can no longer do so silently.
    ///
    /// **It trades reason quality for convenience, on purpose.** `reason` becomes
    /// `error.to_string()`, and for an `io::Error` or an `anyhow::Error` that is
    /// usually unique per record — so
    /// [`dropped_by_reason`](CaptureCounts::dropped_by_reason) degenerates to one
    /// bucket per drop and stops being a grouping. That is an acceptable default
    /// (the detail is the useful thing when a run first refuses), and the
    /// [`Display`](Shortfall) impl caps what it prints so a wide error vocabulary
    /// cannot produce an unreadable panic. When the *shape* of a drop matters —
    /// a shortfall you expect to waive under an [`Allowance`] and read back out
    /// of an artifact — call [`record_drop`](Self::record_drop) with a short
    /// closed-vocabulary reason instead, and keep the varying detail in
    /// `subject`.
    pub fn record<T, E: fmt::Display>(
        &mut self,
        subject: impl Into<String>,
        outcome: Result<T, E>,
    ) -> Option<T> {
        match outcome {
            Ok(value) => {
                self.record_success();
                Some(value)
            }
            Err(error) => {
                self.record_drop(subject, error.to_string());
                None
            }
        }
    }

    /// Records attempted so far.
    #[must_use]
    pub fn attempted(&self) -> usize {
        self.succeeded + self.dropped.len()
    }

    /// Records captured so far.
    #[must_use]
    pub fn succeeded(&self) -> usize {
        self.succeeded
    }

    /// Records dropped so far.
    #[must_use]
    pub fn dropped(&self) -> usize {
        self.dropped.len()
    }

    /// The counts as they stand, without adjudicating them.
    ///
    /// For progress reporting only — reading these instead of calling
    /// [`finish`](Self::finish) is precisely the printed-count pattern that lets
    /// a partial pass reach `fs::write`.
    #[must_use]
    pub fn counts(&self) -> CaptureCounts {
        let mut dropped_by_reason: BTreeMap<String, usize> = BTreeMap::new();
        for record in &self.dropped {
            *dropped_by_reason.entry(record.reason.clone()).or_insert(0) += 1;
        }
        CaptureCounts {
            attempted: self.attempted(),
            succeeded: self.succeeded,
            dropped: self.dropped.len(),
            dropped_by_reason,
            allowance: None,
        }
    }

    /// Close the ledger, refusing on any shortfall.
    ///
    /// # Errors
    ///
    /// [`Shortfall::NothingAttempted`] if the pass attempted no records;
    /// [`Shortfall::Dropped`] if it dropped any.
    pub fn finish(self) -> Result<CaptureCounts, Shortfall> {
        self.close(None)
    }

    /// Close the ledger under a declared [`Allowance`].
    ///
    /// The allowance is stamped into the returned counts, so an artifact records
    /// the waiver alongside the shortfall it waived.
    ///
    /// # Errors
    ///
    /// [`Shortfall::Dropped`] if the drops exceed the allowance's ceiling, and
    /// [`Shortfall::NothingAttempted`] if the pass attempted nothing and the
    /// allowance did not permit that.
    pub fn finish_with(self, allowance: Allowance) -> Result<CaptureCounts, Shortfall> {
        self.close(Some(allowance))
    }

    fn close(self, allowance: Option<Allowance>) -> Result<CaptureCounts, Shortfall> {
        let mut counts = self.counts();
        let empty_permitted = allowance.as_ref().is_some_and(|a| a.empty_pass_permitted);
        if counts.attempted == 0 && !empty_permitted {
            return Err(Shortfall::NothingAttempted {
                subject: self.subject,
            });
        }
        let ceiling = allowance.as_ref().map_or(0, |a| a.max_dropped);
        counts.allowance = allowance.map(|a| AllowanceNote {
            max_dropped: a.max_dropped,
            empty_pass_permitted: a.empty_pass_permitted,
            justification: a.justification,
        });
        if counts.dropped > ceiling {
            let allowed = counts.allowance.as_ref().map(|a| a.max_dropped);
            return Err(Shortfall::Dropped(Box::new(DroppedShortfall {
                subject: self.subject,
                counts,
                records: self.dropped,
                allowed,
            })));
        }
        Ok(counts)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A ledger whose every record succeeded certifies, and reports counts that
    /// add up. `attempted == succeeded + dropped` is an invariant of the type,
    /// not something a caller can get wrong: there is no way to increment
    /// `attempted` on its own.
    #[test]
    fn a_complete_pass_certifies() {
        let mut ledger = CaptureLedger::new("transcript windows");
        ledger.record_success();
        ledger.record_success();

        let counts = ledger.finish().expect("a complete pass certifies");
        assert_eq!(counts.attempted, 2);
        assert_eq!(counts.succeeded, 2);
        assert_eq!(counts.dropped, 0);
        assert!(counts.dropped_by_reason.is_empty());
        assert!(counts.allowance.is_none());
        assert!(counts.is_complete());
    }

    /// The core refusal: one dropped record and `finish` is an `Err`, however
    /// many succeeded. This is what the printed count could not do.
    #[test]
    fn a_single_drop_refuses() {
        let mut ledger = CaptureLedger::new("transcript windows");
        for _ in 0..999 {
            ledger.record_success();
        }
        ledger.record_drop("NM_003002.2", "no stored sequence");

        let error = ledger.finish().expect_err("one drop must refuse");
        match &error {
            Shortfall::Dropped(detail) => {
                assert_eq!(detail.subject, "transcript windows");
                assert_eq!(detail.counts.attempted, 1000);
                assert_eq!(detail.counts.succeeded, 999);
                assert_eq!(detail.counts.dropped, 1);
                assert_eq!(detail.records.len(), 1);
                assert_eq!(detail.records[0].subject, "NM_003002.2");
                assert_eq!(detail.allowed, None);
            }
            other => panic!("expected a Dropped shortfall, got {other:?}"),
        }
    }

    /// The refusal must name the subject, the arithmetic and the reason, and
    /// must say it is refusing — a reader who sees only "3 of 5" cannot tell
    /// whether the generator went on to write.
    #[test]
    fn the_refusal_names_what_was_lost_and_that_it_refused() {
        let mut ledger = CaptureLedger::new("spec rows");
        ledger.record_success();
        ledger.record_drop("LRG_199t1:c.850del", "unresolvable accession");
        ledger.record_drop("NM_007294.3:c.2077delins", "unresolvable accession");

        let message = ledger.finish().expect_err("must refuse").to_string();
        assert!(message.contains("spec rows"), "{message}");
        assert!(message.contains("2 of 3 records were dropped"), "{message}");
        assert!(message.contains("no allowance was declared"), "{message}");
        assert!(
            message.contains("Refusing to write a fixture built from a partial pass"),
            "{message}"
        );
        assert!(message.contains("unresolvable accession: 2"), "{message}");
        assert!(message.contains("LRG_199t1:c.850del"), "{message}");
    }

    /// Long failures are truncated so a panic stays readable, but the total is
    /// still stated exactly — a summary that hid the count would reintroduce the
    /// very ambiguity this module removes.
    #[test]
    fn a_long_refusal_is_truncated_but_still_states_the_total() {
        let mut ledger = CaptureLedger::new("corpus rows");
        for i in 0..(MAX_NAMED_SUBJECTS + 5) {
            ledger.record_drop(format!("row-{i}"), "parse failed");
        }

        let message = ledger.finish().expect_err("must refuse").to_string();
        assert!(
            message.contains("15 of 15 records were dropped"),
            "{message}"
        );
        assert!(message.contains("row-0"), "{message}");
        assert!(!message.contains("row-14"), "{message}");
        assert!(message.contains("and 5 more"), "{message}");
    }

    /// An allowance waives a shortfall up to its ceiling, and the waiver is
    /// carried in the counts so the artifact records who permitted what.
    #[test]
    fn a_declared_allowance_waives_up_to_its_ceiling() {
        let mut ledger = CaptureLedger::new("reference accessions");
        ledger.record_success();
        ledger.record_drop("LRG_199", "unversioned reference");

        let counts = ledger
            .finish_with(Allowance::at_most(
                1,
                "bare LRG_ ids have no versioned index entry",
            ))
            .expect("within the declared allowance");
        assert_eq!(counts.dropped, 1);
        // Borrowed, not moved: the assertion below has to be about the counts
        // the ledger produced. An earlier revision moved the allowance out here
        // and then asserted on a hand-built `CaptureCounts`, which held whatever
        // `finish_with` returned — including a complete claim (#1551 review).
        let note = counts
            .allowance
            .as_ref()
            .expect("the allowance is recorded");
        assert_eq!(note.max_dropped, 1);
        assert_eq!(
            note.justification,
            "bare LRG_ ids have no versioned index entry"
        );
        // Still not a complete pass — an allowance permits a shortfall, it does
        // not turn one into a full capture.
        assert!(!counts.is_complete());
    }

    /// An allowance is a ceiling, not a licence: one drop past it still refuses,
    /// and the message states the ceiling that was exceeded.
    #[test]
    fn an_allowance_is_a_ceiling_not_a_licence() {
        let mut ledger = CaptureLedger::new("reference accessions");
        ledger.record_drop("LRG_199", "unversioned reference");
        ledger.record_drop("LRG_200", "unversioned reference");

        let error = ledger
            .finish_with(Allowance::at_most(1, "at most one bare LRG_ id"))
            .expect_err("two drops exceed a ceiling of one");
        assert!(matches!(&error, Shortfall::Dropped(detail) if detail.allowed == Some(1)));
        assert!(
            error
                .to_string()
                .contains("exceeding the declared allowance of 1"),
            "{error}"
        );
    }

    /// `Allowance::at_most(0, ..)` is meaningful: it records that the question
    /// was asked, and still refuses a drop.
    #[test]
    fn a_zero_allowance_records_the_question_and_still_refuses() {
        let mut clean = CaptureLedger::new("rows");
        clean.record_success();
        let counts = clean
            .finish_with(Allowance::at_most(0, "no drop is acceptable here"))
            .expect("a clean pass passes under a zero allowance");
        assert_eq!(
            counts
                .allowance
                .expect("the allowance is recorded even at zero")
                .max_dropped,
            0
        );

        let mut lossy = CaptureLedger::new("rows");
        lossy.record_drop("row-0", "parse failed");
        assert!(lossy
            .finish_with(Allowance::at_most(0, "no drop is acceptable here"))
            .is_err());
    }

    /// A pass that attempted nothing is a shortfall of its own. An empty
    /// artifact reads as "there was nothing to find", which is exactly how a
    /// structural blindness disguises itself as a safe result.
    #[test]
    fn an_empty_pass_refuses_unless_it_is_permitted() {
        let error = CaptureLedger::new("corpus rows")
            .finish()
            .expect_err("an empty pass must refuse");
        assert!(matches!(error, Shortfall::NothingAttempted { .. }));
        assert!(error.to_string().contains("attempted 0 records"), "{error}");
        assert!(error.dropped_records().is_empty());

        // A bare allowance does not waive it — permitting drops says nothing
        // about permitting an empty population.
        assert!(CaptureLedger::new("corpus rows")
            .finish_with(Allowance::at_most(5, "some rows may fail"))
            .is_err());

        let counts = CaptureLedger::new("corpus rows")
            .finish_with(
                Allowance::at_most(0, "the corpus is optional here").permitting_an_empty_pass(),
            )
            .expect("an explicitly permitted empty pass");
        assert_eq!(counts.attempted, 0);
        assert!(!counts.is_complete());
    }

    /// `record` is the adoption path for `let Ok(v) = .. else { continue }`: it
    /// returns the same `Option` the caller was already branching on, and the
    /// drop is accounted for on the way past.
    #[test]
    fn record_accounts_for_a_fallible_step_and_yields_its_value() {
        let mut ledger = CaptureLedger::new("transcripts");
        let ok: Result<u32, String> = Ok(7);
        let err: Result<u32, String> = Err("no stored sequence".to_string());

        assert_eq!(ledger.record("NM_000088.3", ok), Some(7));
        assert_eq!(ledger.record("NM_003002.2", err), None);
        assert_eq!(ledger.attempted(), 2);
        assert_eq!(ledger.succeeded(), 1);
        assert_eq!(ledger.dropped(), 1);

        let error = ledger.finish().expect_err("the drop must refuse");
        assert_eq!(error.dropped_records()[0].reason, "no stored sequence");
    }

    /// The counts are the artifact's completeness claim, so they must survive a
    /// JSON round trip — a claim a fixture cannot carry is not a claim a
    /// consuming test can assert on.
    #[test]
    fn counts_round_trip_through_json_so_an_artifact_can_carry_them() {
        let mut ledger = CaptureLedger::new("windows");
        ledger.record_success();
        ledger.record_drop("NM_003002.2", "no stored sequence");
        let counts = ledger
            .finish_with(Allowance::at_most(1, "one transcript is not provisioned"))
            .expect("within the allowance");

        let json = serde_json::to_string(&counts).expect("serialize counts");
        let parsed: CaptureCounts = serde_json::from_str(&json).expect("deserialize counts");
        assert_eq!(parsed, counts);
        assert_eq!(parsed.dropped_by_reason["no stored sequence"], 1);
        assert_eq!(
            parsed.allowance.expect("allowance survives").justification,
            "one transcript is not provisioned"
        );
    }

    /// The refusal message stays bounded even when every drop carries its own
    /// reason — which is what [`CaptureLedger::record`] produces for an error
    /// type whose `Display` interpolates the record.
    ///
    /// Before `MAX_NAMED_REASONS`, `MAX_NAMED_SUBJECTS` was defeated by exactly
    /// this: the subject list was capped at 10 while `dropped_by_reason` printed
    /// in full, so 2000 drops rendered a 2,014-line, 71 KB panic (#1551 review).
    #[test]
    fn a_wide_reason_vocabulary_cannot_produce_an_unbounded_refusal() {
        let mut ledger = CaptureLedger::new("corpus rows");
        let drops = 2000;
        for i in 0..drops {
            let outcome: Result<(), String> = Err(format!("row {i} failed to resolve"));
            ledger.record(format!("row-{i}"), outcome);
        }

        let message = ledger.finish().expect_err("must refuse").to_string();
        let lines = message.lines().count();
        assert!(
            lines <= MAX_NAMED_SUBJECTS + MAX_NAMED_REASONS + 8,
            "refusal must stay readable, got {lines} lines:\n{message}"
        );
        // The totals survive the truncation on both axes — a summary that hid
        // the count would reintroduce the ambiguity this module removes.
        assert!(
            message.contains(&format!("{drops} of {drops} records were dropped")),
            "{message}"
        );
        assert!(
            message.contains(&format!("by reason ({drops} distinct)")),
            "{message}"
        );
        assert!(
            message.contains(&format!(
                "and {} more reasons covering {} drops",
                drops - MAX_NAMED_REASONS,
                drops - MAX_NAMED_REASONS
            )),
            "{message}"
        );
    }

    /// An allowance with a blank justification is the state [`Allowance`]'s own
    /// doc calls indistinguishable from never having considered the question, so
    /// "required" has to mean content and not merely presence.
    #[test]
    #[should_panic(expected = "an Allowance needs a justification")]
    fn an_allowance_with_a_blank_justification_is_refused() {
        let _ = Allowance::at_most(0, "   ");
    }

    /// `attempted == succeeded + dropped` holds by construction through the
    /// ledger, but nothing re-checks it on the way back out of an artifact — and
    /// deserializing a stamped claim is exactly how a consuming test reads one.
    #[test]
    fn self_consistency_is_checkable_on_counts_read_back_from_an_artifact() {
        let mut ledger = CaptureLedger::new("windows");
        ledger.record_success();
        ledger.record_drop("NM_003002.2", "no stored sequence");
        let counts = ledger
            .finish_with(Allowance::at_most(1, "one transcript is not provisioned"))
            .expect("within the allowance");
        assert!(counts.is_self_consistent());

        // Arithmetic consistency is not artifact verification, and this is the
        // limitation the module documents rather than a case the checker
        // catches: 100 == 100 + 0 holds, so a claim invented wholesale passes.
        // What `is_self_consistent` rules out is a claim that contradicts
        // *itself*, never one that contradicts the file it is stamped onto.
        let internally_consistent_but_unverifiable: CaptureCounts = serde_json::from_str(
            r#"{"attempted":100,"succeeded":100,"dropped":0,"dropped_by_reason":{}}"#,
        )
        .expect("deserialize");
        assert!(internally_consistent_but_unverifiable.is_self_consistent());

        let inconsistent: CaptureCounts =
            serde_json::from_str(r#"{"attempted":100,"succeeded":10,"dropped":0}"#)
                .expect("deserialize");
        assert!(!inconsistent.is_self_consistent());
    }

    /// The counts this method exists to check are the ones nothing enforced, so
    /// its own arithmetic has to survive whatever a file holds. `succeeded +
    /// dropped` on two near-`usize::MAX` values panics in debug and wraps in
    /// release, and a wrapped sum can land back on `attempted` — reporting a
    /// claim as consistent *because* it overflowed (#1551 review).
    #[test]
    fn self_consistency_treats_overflowing_counts_as_inconsistent() {
        let overflowing: CaptureCounts = serde_json::from_str(&format!(
            r#"{{"attempted":1,"succeeded":{max},"dropped":{max}}}"#,
            max = usize::MAX
        ))
        .expect("deserialize");
        assert!(!overflowing.is_self_consistent());

        let overflowing_reasons: CaptureCounts = serde_json::from_str(&format!(
            r#"{{"attempted":{max},"succeeded":0,"dropped":{max},
                 "dropped_by_reason":{{"a":{max},"b":{max}}}}}"#,
            max = usize::MAX
        ))
        .expect("deserialize");
        assert!(!overflowing_reasons.is_self_consistent());

        // The case the wrap actually costs something on: `MAX + 2` wraps to 1,
        // which *is* `attempted`, and the per-reason tally agrees. A release
        // build of the unchecked version called this claim consistent.
        let wraps_onto_attempted: CaptureCounts = serde_json::from_str(&format!(
            r#"{{"attempted":1,"succeeded":{max},"dropped":2,
                 "dropped_by_reason":{{"unversioned reference":2}}}}"#,
            max = usize::MAX
        ))
        .expect("deserialize");
        assert!(!wraps_onto_attempted.is_self_consistent());
    }

    /// A clean pass's counts stay compact in the artifact: no empty reason map,
    /// no null allowance. A fixture diff should show the claim, not the absence
    /// of one.
    #[test]
    fn a_clean_claim_serializes_without_empty_fields() {
        let mut ledger = CaptureLedger::new("windows");
        ledger.record_success();
        let json = serde_json::to_string(&ledger.finish().expect("clean")).expect("serialize");
        assert_eq!(json, r#"{"attempted":1,"succeeded":1,"dropped":0}"#);
    }
}
