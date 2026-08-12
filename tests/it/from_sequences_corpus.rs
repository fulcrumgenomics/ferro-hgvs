//! `from_sequences` judged against the committed corpora.
//!
//! These are the tests that decide whether the surface is right, and they need
//! **no fixtures of their own**: `to_sequences` turns any HGVS corpus already in
//! the tree into a `from_sequences` corpus, and the two functions check each
//! other. Every description this repository has collected becomes an input.
//!
//! # What each check can and cannot see
//!
//! In increasing strength, and deliberately not collapsed into one assertion —
//! a failure that says *which* property broke is worth more than a boolean:
//!
//! | check | catches |
//! |---|---|
//! | re-parse | a malformed output |
//! | rule-1 prohibitions | a description the spec forbids outright |
//! | **denotation preserved** | a description of the wrong bases |
//! | idempotence | a derivation that moves its own output |
//! | **confluence** | two spellings of one variant reaching two answers |
//!
//! # The corpus is `--axes g,c`; only the `g.` half is reachable
//!
//! Read every number below as a claim about **5 636 genomic classes**, not about
//! the corpus's 11 272. The `c.` classes are drawn against `NM_TEST.1` and
//! `from_sequences` refuses every non-genomic accession by design, so not one of
//! them can reach a comparison — they are excluded *structurally*, before any
//! property is tested.
//!
//! That is stated here, and asserted in [`AxisCensus`], because a corpus zero is
//! a claim about the corpus: "5 636 classes with no divergence" was
//! arithmetically `11 272 / 2` while reading as a result over the whole thing.
//! The `unusable` counter that recorded the refusals was printed and never
//! asserted, and the class floor sat *below* the g-only value, so nothing in the
//! test could have noticed if the `c.` half had gone missing for some other
//! reason entirely.
//!
//! **`EquivalenceChecker` is not among them, and that is deliberate.** It
//! normalizes both sides before comparing, so using it here would be circular —
//! the repository's own guidance says so. Denotation is compared through
//! `to_sequences`, which reaches the bases via `hgvs_to_spdi` and an SPDI
//! splice; that applier is not the normalizer and cannot agree with a derivation
//! merely because the derivation produced it.

use crate::cis_confluence_axis::{corpus, reference_for, Class};
use ferro_hgvs::{from_sequences, FromSequencesOptions, HgvsVariant, MockProvider, Normalizer};
use std::collections::BTreeSet;

/// The pad `to_sequences` is asked for, matching `merge::CANONICAL_PAD` — the
/// room the normalizer's own derivation gives its 3' placement. With `0` every
/// member would sit against a window edge and every answer would be reported as
/// read-dependent, which is true but useless.
const PAD: u64 = 128;

/// One spelling's journey: description -> sequences -> description.
///
/// `None` when the corpus row cannot be put through the pipeline at all — an
/// unparseable spelling, or a variant `to_sequences` declines (overlapping
/// members, an edit SPDI cannot carry). Those are **counted by the caller**
/// rather than silently dropped: a skip that reads as a pass is the failure mode
/// these corpora exist to remove.
fn derive(normalizer: &Normalizer<MockProvider>, spelling: &str) -> Option<HgvsVariant> {
    let variant = ferro_hgvs::parse_hgvs(spelling).ok()?;
    let pair = normalizer.to_sequences(&variant, PAD).ok()?;
    from_sequences(
        &pair.accession,
        pair.position,
        &pair.reference,
        &pair.alternate,
        &FromSequencesOptions::default(),
    )
    .ok()
}

/// Prohibitions that make a description invalid outright, not merely
/// non-preferred.
///
/// A deliberately small set, taken from the clauses `spec_conformance_axis`
/// already reads as absolute. It is **not** a conformance suite: the point here
/// is that a derivation must never produce one of these, which is a much
/// narrower claim than "the output is the recommended form".
fn violated_prohibition(output: &str) -> Option<&'static str> {
    if output.contains(' ') {
        return Some("general.md:96 — spaces are not permitted in any HGVS description");
    }
    if output.contains('X') || output.contains('x') {
        return Some("standards.md:39 — `X` is an alignment-only symbol, not a nucleotide");
    }
    if let Some(body) = output.split_once(":g.").map(|(_, body)| body) {
        if body.contains('+') || body.contains('-') || body.contains('*') {
            return Some("checklist.md:16 — a `g.` description admits no `+`/`-`/`*` offset");
        }
    }
    None
}

/// Group the corpus by the reference its classes are drawn against.
///
/// Classes are sorted by an id beginning with the core index and axis, so every
/// class sharing a reference is contiguous — building one provider per group
/// rather than per class is the difference between seconds and minutes, the
/// same reason `cis_confluence_axis::measure` does it.
fn for_each_group(mut visit: impl FnMut(&Normalizer<MockProvider>, &[Class])) {
    // `chunk_by` rather than an index cursor: the hand-rolled version's
    // `start`/`end` bookkeeping is what the `continue`-vs-`return` note below
    // records having already gone wrong once.
    for group in corpus()
        .classes
        .chunk_by(|a, b| a.core == b.core && a.axis == b.axis)
    {
        let (provider, _reference) = reference_for(&group[0]);
        visit(&Normalizer::new(provider), group);
    }
}

/// **Every spelling of one variant derives to one description.**
///
/// This is the property the whole surface exists for, and the cis confluence
/// corpus is already built to express it: each class is a set of spellings
/// verified — independently of the normalizer — to denote one sequence. Since
/// the derivation is a function of `(reference block, alternate block)` and
/// every spelling in a class produces the same pair, one output per class is
/// **structural** here, not best effort.
///
/// Asserted at zero with the evaluated count asserted non-zero beside it, so a
/// corpus that failed to load cannot pass as a clean result.
#[test]
fn every_spelling_of_a_class_derives_to_one_description() {
    let mut evaluated = 0usize;
    let mut divergent: Vec<String> = Vec::new();
    let mut divergent_count = 0usize;
    let mut unusable = 0usize;
    let mut axis = AxisCensus::default();

    for_each_group(|normalizer, classes| {
        for class in classes {
            let mut outputs = BTreeSet::new();
            let mut seen = 0usize;
            for spelling in &class.spellings {
                match derive(normalizer, spelling) {
                    Some(derived) => {
                        seen += 1;
                        outputs.insert(derived.to_string());
                    }
                    None => unusable += 1,
                }
            }
            axis.record(&class.axis, seen);
            // A class only tests confluence if at least two of its spellings got
            // through; one is a trivial pass and would inflate the denominator.
            //
            // `continue`, not `return`: this closure is called once per *group*
            // of classes sharing a reference, so returning here would abandon
            // every remaining class in the group. That bug read as a corpus that
            // would not load — 32 classes evaluated instead of thousands.
            if seen < 2 {
                continue;
            }
            evaluated += 1;
            if outputs.len() > 1 {
                divergent_count += 1;
                // The sample is capped; the count is not. Reporting `len()` of a
                // capped sample as the result reads as "exactly ten diverged",
                // which is a claim the sample cannot support.
                if divergent.len() < 10 {
                    divergent.push(format!("{}: {:?}", class.id, outputs));
                }
            }
        }
    });

    eprintln!(
        "from_sequences confluence: {evaluated} classes with 2+ usable spellings, \
         {unusable} spellings from_sequences could not derive\n  axis: {axis:?}"
    );
    axis.assert_the_transcript_axis_is_structurally_excluded(unusable);
    // The floor is the **genomic** class count, which is what this test can
    // reach at all — see `AxisCensus`. It was `> 5_000` against a corpus of
    // 11 272 classes of which 5 636 are genomic, so it sat 636 under the only
    // achievable value and could not have noticed the entire `c.` half going
    // missing. Pinned just under the genomic half instead, so a real loss of
    // classes fails rather than passing at half strength.
    assert!(
        evaluated >= axis.genomic_classes,
        "{} of {} genomic classes had fewer than two derivable spellings",
        axis.genomic_classes - evaluated,
        axis.genomic_classes
    );
    assert!(
        evaluated > 5_500,
        "corpus did not load: only {evaluated} classes evaluated"
    );
    assert!(
        divergent_count == 0,
        "{divergent_count} of {evaluated} classes reached more than one description \
         (first {} shown):\n{}",
        divergent.len(),
        divergent.join("\n")
    );
}

/// How many classes each axis contributed, and how many of them this surface
/// could evaluate at all.
///
/// # Why this is here: the headline zero is a claim about half the corpus
///
/// The corpus is generated `--axes g,c`, and its `c.` classes are drawn against
/// `NM_TEST.1`. `from_sequences` refuses every non-genomic accession by design
/// and says so at `from_sequences.rs` — so **no `c.` class can reach `evaluated`
/// at all**. Every one of them lands with `seen == 0`, hits the `seen < 2`
/// `continue`, and vanishes.
///
/// That made "5 636 classes with no divergence" arithmetically `11 272 / 2`: the
/// genomic half, reported as though it were the corpus. `unusable` counted the
/// refusals and was `eprintln!`'d and never asserted, and the `evaluated > 5_000`
/// floor sat *under* the g-only value, so neither could notice.
///
/// This repository's doctrine is that **a corpus zero is a claim about the
/// corpus**, and that a structural zero must be reported as structural. So the
/// exclusion is asserted rather than described: the transcript axis must
/// contribute classes (or the corpus stopped generating them, and the exclusion
/// this documents is no longer the reason for the shape of the result), and none
/// of them may evaluate (or the refusal has changed and the headline now covers
/// more than it says).
#[derive(Debug, Default)]
struct AxisCensus {
    /// Classes on a `g.` accession — the ones this surface can judge.
    genomic_classes: usize,
    /// … of those, the ones with two or more derivable spellings.
    genomic_evaluated: usize,
    /// Classes on a transcript accession (`c.`), which `from_sequences` refuses.
    transcript_classes: usize,
    /// … of those, any that nonetheless produced a derivable spelling. Must be
    /// zero, and is asserted so rather than assumed.
    transcript_derivable: usize,
}

impl AxisCensus {
    fn record(&mut self, axis: &str, seen: usize) {
        if axis == "g" {
            self.genomic_classes += 1;
            if seen >= 2 {
                self.genomic_evaluated += 1;
            }
        } else {
            self.transcript_classes += 1;
            if seen > 0 {
                self.transcript_derivable += 1;
            }
        }
    }

    /// The structural zero, asserted from both sides.
    fn assert_the_transcript_axis_is_structurally_excluded(&self, unusable: usize) {
        assert!(
            self.transcript_classes > 0,
            "the corpus generated no transcript-axis classes, so the exclusion this test \
             documents is not the reason for its shape any more — re-read the census: {self:?}"
        );
        assert_eq!(
            self.transcript_derivable, 0,
            "a transcript-axis class produced a derivable spelling; `from_sequences` is \
             documented as refusing every non-genomic accession, so either the refusal moved \
             or the corpus did: {self:?}"
        );
        // The refused spellings are the ones `unusable` counts, so the counter
        // is tied to something rather than merely printed. It is a floor, not an
        // equality: `to_sequences` also declines a handful of genomic rows.
        assert!(
            unusable >= self.transcript_classes,
            "only {unusable} spellings were unusable, fewer than the {} transcript-axis \
             classes that must all refuse — the accession gate is not firing where this test \
             says it does",
            self.transcript_classes
        );
    }
}

/// **Every derived description denotes the bases it was derived from, parses,
/// violates no absolute prohibition, and is idempotent.**
///
/// Four properties in one sweep over the same corpus, because they share the
/// expensive half — building a provider and applying to the reference — and
/// separating them would triple the runtime to report the same failures.
/// Each is counted and reported apart, so a red run says which one broke.
#[test]
fn every_derived_description_is_valid_and_stable() {
    let mut checked = 0usize;
    let mut unparseable: Vec<String> = Vec::new();
    let mut prohibited: Vec<String> = Vec::new();
    let mut redenoted: Vec<String> = Vec::new();
    let mut redenoted_count = 0usize;
    let mut unstable: Vec<String> = Vec::new();

    for_each_group(|normalizer, classes| {
        for class in classes {
            for spelling in &class.spellings {
                let Ok(variant) = ferro_hgvs::parse_hgvs(spelling) else {
                    continue;
                };
                let Ok(pair) = normalizer.to_sequences(&variant, PAD) else {
                    continue;
                };
                let Ok(derived) = from_sequences(
                    &pair.accession,
                    pair.position,
                    &pair.reference,
                    &pair.alternate,
                    &FromSequencesOptions::default(),
                ) else {
                    continue;
                };
                checked += 1;
                let rendered = derived.to_string();

                if ferro_hgvs::parse_hgvs(&rendered).is_err() && unparseable.len() < 5 {
                    unparseable.push(rendered.clone());
                }
                if let Some(clause) = violated_prohibition(&rendered) {
                    if prohibited.len() < 5 {
                        prohibited.push(format!("{rendered}: {clause}"));
                    }
                }

                // Denotation, compared as an encoding-invariant SPDI key rather
                // than as two windows.
                //
                // Windows are the wrong comparand and comparing them was this
                // test's own first bug: `to_sequences` returns the members' own
                // span, so a derivation that narrows the span narrows the
                // window, and two correct answers compare unequal. `canonical_spdi`
                // is derived from the resulting bases and is independent of both
                // the span and the spelling — which is exactly the question
                // being asked. It also reaches the bases through `hgvs_to_spdi`
                // rather than through the derivation, so it cannot agree by
                // construction.
                match (
                    normalizer.canonical_spdi(&variant),
                    normalizer.canonical_spdi(&derived),
                ) {
                    (Ok(before), Ok(after)) if before != after => {
                        // Counted before sampled, for the reason the confluence
                        // test records: the length of a capped sample is not the
                        // size of the population it was drawn from.
                        redenoted_count += 1;
                        if redenoted.len() < 5 {
                            redenoted
                                .push(format!("{spelling} -> {rendered}: {before} vs {after}"));
                        }
                    }
                    // Either side declining is not evidence about the other:
                    // `canonical_spdi` refuses shapes SPDI cannot carry, and an
                    // input it refuses has no key to compare against.
                    _ => {}
                }

                // Idempotence: deriving from the derived form's own window must
                // reproduce it.
                if let Ok(after) = normalizer.to_sequences(&derived, PAD) {
                    if let Ok(again) = from_sequences(
                        &after.accession,
                        after.position,
                        &after.reference,
                        &after.alternate,
                        &FromSequencesOptions::default(),
                    ) {
                        if again.to_string() != rendered && unstable.len() < 5 {
                            unstable.push(format!("{rendered} -> {again}"));
                        }
                    }
                }
            }
        }
    });

    eprintln!("from_sequences validity: {checked} derived descriptions checked");
    assert!(
        checked > 20_000,
        "corpus did not load: only {checked} checked"
    );
    assert!(
        unparseable.is_empty(),
        "unparseable output:\n{unparseable:#?}"
    );
    assert!(prohibited.is_empty(), "prohibited output:\n{prohibited:#?}");
    assert!(
        redenoted_count == 0,
        "denotation changed on {redenoted_count} of {checked} (first {} shown):\n{redenoted:#?}",
        redenoted.len()
    );
    assert!(unstable.is_empty(), "not idempotent:\n{unstable:#?}");
}
