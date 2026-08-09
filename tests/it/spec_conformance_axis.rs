//! The four properties, over the exhaustive spec-derived conformance corpus.
//!
//! # What this measures, and why it can measure it without an oracle
//!
//! Conformance used to be measured by the generated spec fixture: 934 harvested
//! rows, **52** of them multi-member cis alleles and **40** carrying a
//! spec-stated answer. That ceiling is the spec's, not ours — "matches the form
//! the spec publishes" needs the spec to have published one — so a change to the
//! multi-member partitioner measured `0 gained, 0 lost`. That is measurement
//! blindness, not neutrality.
//!
//! `ferro_hgvs::conformance::spec_corpus` gets past it by generating *many
//! spellings of one variant*: re-partition it, and shift each partition within
//! its ambiguous run. Every spelling in such a family denotes the same sequence
//! **by construction**, which makes four properties checkable with no published
//! answer at all:
//!
//! | rank | property | how it is checked |
//! |---|---|---|
//! | 1 | **validity** | the output re-parses, denotes a sequence, and violates no absolute prohibition |
//! | 2 | **confluence** | every spelling of one variant reaches ONE output |
//! | — | **idempotency** | each output is its own fixed point |
//! | — | **sequence preservation** | each output denotes the input's bases, via `hgvs_to_spdi` |
//!
//! Rank 1 and rank 2 are the operator's precedence order
//! (`rulings[adjudication-precedence-order]`); the other two are invariants
//! rather than ranked choices.
//!
//! # This is a census, not a gate
//!
//! None of these properties holds today and a permanently red suite is worth
//! nothing, so every figure below is a **pinned baseline**: each failure figure
//! may only ever go DOWN and `converged` only ever UP. A change that moves one
//! must re-bless the pin in the same commit and say so in the PR.
//!
//! # The pins are NOT an acceptable state — five counters are rank-1 defects
//!
//! Pinning a non-zero figure records a defect; it does not excuse one. These five
//! are conformance defects that this instrument found on its first run, and none
//! is fixed here (a corpus PR that also changed the normalizer would make its own
//! baseline unmeasurable).
//!
//! **The base they are measured against is `main`.** Every figure below was
//! measured on this branch, which adds no `src/normalize/` change of its own —
//! so the census records `main`'s behaviour, not this PR's. The corpus is new;
//! the defects it counts are not.
//!
//! | counter | 3' | 5' | what it is |
//! |---|---|---|---|
//! | `outputs_leaving_the_transcript` | **371** | 0 | a `c.`/`n.` output naming an INTRONIC position its input did not. `general.md:44` exempts junction-adjacent deletions/duplications from the 3' rule; `checklist.md:20` says a bare `NM_` cannot express an intronic position at all. **Entirely minus-strand**: the same design on a plus-strand transcript clamps at the exon end correctly |
//! | `outputs_denoting_no_sequence` | 10 | 18 | two members of the OUTPUT claim one territory, so it denotes nothing — rank-1 invalid |
//! | `sequence_changed` | 4 | 0 | the output denotes different bases from the input; a member was dropped |
//! | `non_idempotent_outputs` | 7 | 4 | the output is not its own fixed point |
//! | `conflicts_accepted` | 72 | 72 | a conflicting allele was normalized instead of refused — nested, partially overlapping, and two insertions at one interbase |
//! | `prohibited_absolute_accepted` | 32 | 32 | a shape the spec calls "not allowed" was accepted |
//! | `prohibition_violating_outputs` | 32 | 32 | and then EMITTED, so the prohibition is not enforced on output either |
//!
//! `guard_violations` is **0 of 210 guarded rows**, which is a real negative
//! result rather than a vacuous one: ferro does not merge an irreducible frameless
//! separation of one, so it does not implement rejected SVD-WG010. The denominator
//! is asserted non-zero, because `0 of 0` is what a rebuilt #1456 looks like.
//!
//! The 3'/5' asymmetry is itself a finding: the transcript-leaving class is
//! **371 to 0**, and 3' is the default direction.
//!
//! # The honest-zero discipline is enforced, not described
//!
//! Every property whose denominator is zero fails loudly as **VACUOUS** rather
//! than passing. A corpus that stopped generating conflict rows would otherwise
//! report `0 accepted` — indistinguishable from perfect refusal, and exactly what
//! #1456 did three times (`0 of 18,432`, reported as evidence of neutrality when
//! the corpus simply could not build a conflicting allele).
//!
//! # No fixture, deliberately
//!
//! The corpus is enumerated **at test time** from the library. There is no
//! generated artifact to stage, no CI cache entry to add, and no committed data
//! file — so this runs on every `cargo test` and every CI shard with no wiring at
//! all. `examples/generate_spec_conformance_corpus.rs` exists for the other half
//! of the job: checking the committed rule inventory's citations against the spec
//! checkout, which needs the submodule and therefore belongs in CI's
//! `spec-fixtures` job.
//!
//! # Tier 2: a failing row is promoted, not left in the corpus
//!
//! A corpus row can silently stop being generated — that is the #1456/#1460/#1478
//! class exactly. So every row that fails here is promoted to a **named committed
//! test** in `tests/it/spec_corpus_regressions.rs`, with the governing
//! `file.md:line`. Corpus for breadth; named tests for anything that ever went
//! red.

use std::collections::{BTreeMap, BTreeSet};

use ferro_hgvs::conformance::spec_corpus::{
    corpus, denotation_of, denoted_by, CorpusBounds, Denotation, Frame, RefShape, Row, RowKind,
    SpecCorpus, Strength,
};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer, ShuffleDirection};

// ---------------------------------------------------------------------------
// The pinned baselines
// ---------------------------------------------------------------------------

/// The corpus's own shape, pinned separately from the properties measured over
/// it.
///
/// A generator change that halved the corpus would otherwise improve every rate
/// below while measuring less — the failure mode `#1460` is, where a corpus that
/// could not reach `MAX_SPLIT_BLOCK` reported a length-gated change as neutral.
const SHAPE: CorpusShape = CorpusShape {
    designs_considered: 16_650,
    rows: 12_946,
    spellings: 58_552,
    family_rows: 11_882,
    single_rows: 804,
    conflict_rows: 96,
    prohibited_rows: 164,
    multi_member_rows: 12_496,
    guarded_rows: 210,
};

/// The 3'-direction census, pinned. See the module docs for which way each
/// figure is allowed to move.
///
/// # Measured on `main`, and it measures `main`
///
/// This branch adds a corpus, an axis test and a regression suite; it touches no
/// file under `src/normalize/`, `src/hgvs/`, `src/spdi/` or `src/project/`. So
/// every figure here is `main`'s behaviour observed for the first time, not
/// behaviour this change introduced. Re-measured against `main` @ 5b54b2b6, and
/// identical row-for-row to the corpus's first measurement on `main` @ 35de96c8.
///
/// A pin measured somewhere the branch does not reproduce is red on every run
/// and can therefore detect nothing — the instrument would be permanently blind
/// in exactly the way it exists to prevent. That is the whole reason these are
/// pinned at the measured values rather than at zero.
///
/// The confluence figure is the headline and it is **not** good: this corpus is
/// deliberately far harsher than `cis_confluence_axis`'s, because it varies
/// member geometry, transcript geometry and scale rather than holding all three
/// fixed. Read the three asserted-zero counters before the headline, per the
/// module docs.
pub(crate) const THREE_PRIME: Census = Census {
    // -- validity (rank 1) --
    outputs: 58_552,
    declined: 0,
    unparseable_outputs: 0,
    outputs_denoting_no_sequence: 10,
    outputs_leaving_the_transcript: 371,
    prohibition_violating_outputs: 32,
    // -- confluence (rank 2) --
    converged: 9_140,
    split_two: 2_435,
    split_three: 250,
    split_more: 57,
    underdetermined: 0,
    // -- idempotency --
    non_idempotent_outputs: 7,
    // -- sequence preservation --
    sequence_changed: 4,
    // -- refusal --
    conflicts_accepted: 72,
    prohibited_absolute_accepted: 32,
    prohibited_conditional_accepted: 40,
    // -- negative guards --
    guard_violations: 0,
};

/// The 5'-direction census, pinned.
///
/// `--direction 5prime` is a supported public option and confluence is a
/// property of the normalizer rather than of one shuffle direction, so it is
/// measured in full. Two directions landing on materially different numbers
/// would mean a fix was treating a symptom of the shuffle rather than the
/// partitioner.
///
/// Measured on the same base as [`THREE_PRIME`] — `main`, which this branch
/// does not modify.
///
/// The two directions agreeing on every rank-1 counter except
/// `outputs_leaving_the_transcript` is itself the cross-check the two-direction
/// measurement exists for: a figure that appears at 3' and not at 5' is a claim
/// about the code (all three copies of the #670 junction gate are guarded on
/// `ThreePrime` with no 5' mirror), not about one shuffle direction's luck.
pub(crate) const FIVE_PRIME: Census = Census {
    outputs: 58_552,
    declined: 0,
    unparseable_outputs: 0,
    outputs_denoting_no_sequence: 18,
    outputs_leaving_the_transcript: 0,
    prohibition_violating_outputs: 32,
    converged: 8_944,
    split_two: 2_696,
    split_three: 204,
    split_more: 38,
    underdetermined: 0,
    non_idempotent_outputs: 4,
    sequence_changed: 0,
    conflicts_accepted: 72,
    prohibited_absolute_accepted: 32,
    prohibited_conditional_accepted: 40,
    guard_violations: 0,
};

/// The corpus's shape, independent of any property measured over it.
#[derive(Debug, Default, PartialEq, Eq)]
struct CorpusShape {
    designs_considered: usize,
    rows: usize,
    spellings: usize,
    family_rows: usize,
    single_rows: usize,
    conflict_rows: usize,
    prohibited_rows: usize,
    multi_member_rows: usize,
    guarded_rows: usize,
}

impl CorpusShape {
    fn of(built: &SpecCorpus) -> Self {
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
        }
    }
}

/// What one direction's run found. Every field is pinned.
///
/// `pub(crate)`, with `pub(crate)` fields, so sibling modules can reference the
/// pins instead of re-typing them. `corpus_prohibited_inputs.rs` decomposes
/// three of these figures and used to hard-code them; a literal copy of a pin
/// goes stale silently at the next re-bless, and its failure message then
/// accuses the implementation of a defect that is really a stale constant.
#[derive(Debug, Default, PartialEq, Eq)]
pub(crate) struct Census {
    /// Spellings normalized, across every row.
    pub(crate) outputs: usize,
    /// Spellings normalization declined or panicked on, in a stratum where an
    /// output was expected. Counted apart from a divergence, so one cannot hide
    /// the other.
    pub(crate) declined: usize,
    /// Outputs `parse_hgvs` rejects. **Asserted at zero** — rank-1 invalid.
    pub(crate) unparseable_outputs: usize,
    /// Outputs whose members claim overlapping territory, so the description
    /// denotes no sequence. Rank-1 invalid.
    pub(crate) outputs_denoting_no_sequence: usize,
    /// Outputs the oracle cannot express although the input was exonic — which
    /// means the description **left the transcript**: a `c.`/`n.` output naming an
    /// intronic offset position that the input did not name.
    ///
    /// Two clauses bear on it. `general.md:44` exempts "deletions/duplications
    /// around exon/exon junctions using **c.**, **r.** or **n.** reference
    /// sequences" from the 3' rule, so the shift should stop at the junction; and
    /// `checklist.md:20` says an `NM_` "can only be used to describe variants in
    /// introns using a `c.` prefix when a genomic reference sequence is given", so
    /// the output is not even expressible against the accession it carries.
    pub(crate) outputs_leaving_the_transcript: usize,
    /// Outputs violating an absolute textual prohibition.
    ///
    /// **Ratchet-pinned at 32, not asserted at zero** — this doc said the latter
    /// while both censuses pinned 32, which reads as "ferro emits no prohibited
    /// output" to anyone who does not go and check the pin. It emits 32, and
    /// `corpus_prohibited_inputs` decomposes them: 8 × `checklist.md:16` and
    /// 24 × `standards.md:39`, every one a re-emission of an input that already
    /// violated the clause rather than a violation ferro manufactured.
    pub(crate) prohibition_violating_outputs: usize,
    /// Families whose spellings all reached ONE output. The goal is for this to
    /// equal `family_rows`.
    pub(crate) converged: usize,
    /// Families reaching exactly two distinct outputs.
    pub(crate) split_two: usize,
    /// … exactly three.
    pub(crate) split_three: usize,
    /// … four or more.
    pub(crate) split_more: usize,
    /// Families left with fewer than two spellings that normalized at all, so
    /// confluence is not decidable. Counted rather than skipped: a rise here
    /// would hollow the axis out while every other figure improved.
    pub(crate) underdetermined: usize,
    /// Outputs that are not their own fixed point.
    pub(crate) non_idempotent_outputs: usize,
    /// Outputs whose applied sequence differs from their row's.
    ///
    /// **Ratchet-pinned, and it is not zero in both directions**: 4 at 3' and 0
    /// at 5'. The zero is the 5' census alone, so a doc claiming the counter is
    /// asserted at zero describes half the axis. The four are the CDS-end flush
    /// pairs `spec_corpus_regressions` dissects, and they are the most serious
    /// class the axis measures — a normalization that changes the sequence its
    /// input denoted.
    pub(crate) sequence_changed: usize,
    /// Conflicting alleles — two members claiming one territory — the
    /// implementation accepted instead of refusing.
    pub(crate) conflicts_accepted: usize,
    /// Shapes the spec prohibits in as many words ("is not allowed", "is not
    /// correct") that the implementation accepted.
    pub(crate) prohibited_absolute_accepted: usize,
    /// Shapes a conditional clause prohibits that the implementation accepted. A
    /// finding to adjudicate rather than a regression.
    pub(crate) prohibited_conditional_accepted: usize,
    /// Rows whose output implements a REJECTED consultation proposal — a
    /// separation floor of two on a frameless axis, which is SVD-WG010
    /// (`consultation/SVD-WG010.md:8`, "The proposal has been **rejected**").
    pub(crate) guard_violations: usize,
}

// ---------------------------------------------------------------------------
// Absolute prohibitions, checked on the output
// ---------------------------------------------------------------------------

/// Whether `output` violates a prohibition the spec states in as many words, and
/// which.
///
/// Deliberately a **closed, small** list of textual checks rather than a general
/// validator: each entry cites a clause that uses prohibitive words, so a
/// non-zero count is a conformance defect and not an adjudication. The parser
/// already rejects most of these on *input*; the question here is whether the
/// normalizer can emit one.
fn violated_prohibition(output: &str) -> Option<&'static str> {
    // `general.md:96`: "spaces are *not* permitted in any HGVS description".
    if output.contains(' ') {
        return Some("general.md:96 — spaces are not permitted in any HGVS description");
    }
    // `background/standards.md:39` footnotes `X` and `-` as "used in alignment
    // only", so neither is a nucleotide a description may state. `X` cannot
    // appear anywhere in a legal DNA/RNA description.
    if output.contains('X') || output.contains('x') {
        return Some("standards.md:39 — `X` is an alignment-only symbol, not a nucleotide");
    }
    // `checklist.md:16`: a genomic reference "can not have nucleotides with
    // additions like a `+`, `-`, or `*`". Checked only on the description body,
    // since an accession may legitimately hold neither.
    if let Some(body) = output.split_once(":g.").map(|(_, body)| body) {
        if body.contains('+') || body.contains('*') || body.contains('-') {
            return Some("checklist.md:16 — a `g.` description admits no `+`/`-`/`*` offset");
        }
    }
    None
}

/// Members of a parsed description, for the negative guard and the applier.
fn members_of(variant: &HgvsVariant) -> Vec<HgvsVariant> {
    match variant {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single.clone()],
    }
}

// ---------------------------------------------------------------------------
// Measurement
// ---------------------------------------------------------------------------

/// One divergence, for the failure message when a pin moves.
struct Divergence {
    id: String,
    outputs: Vec<String>,
}

/// A finding worth naming in the report: a validity failure, an accepted
/// prohibition, or a negative-guard violation.
struct Finding {
    id: String,
    what: String,
}

struct Measured {
    census: Census,
    divergences: Vec<Divergence>,
    findings: Vec<Finding>,
}

/// The built corpus, shared across both directions.
///
/// Enumerating it costs ~0.3 s, so the two direction tests each build their own
/// rather than sharing a `OnceLock` — which keeps them independent and lets
/// nextest run them in parallel.
fn built() -> SpecCorpus {
    corpus(&CorpusBounds::default())
}

/// Rows in an order that groups every row sharing a synthetic reference, so one
/// provider is built per reference rather than per row.
fn grouped(built: &SpecCorpus) -> Vec<&Row> {
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

fn measure(direction: ShuffleDirection) -> Measured {
    let built = built();
    let mut census = Census::default();
    let mut divergences: Vec<Divergence> = Vec::new();
    let mut findings: Vec<Finding> = Vec::new();

    let rows = grouped(&built);
    let mut frame: Option<(RefShape, String, Frame, Normalizer<MockProvider>)> = None;
    for row in rows {
        let key_matches = frame
            .as_ref()
            .is_some_and(|(shape, core, _, _)| *shape == row.shape && core == &row.core);
        if !key_matches {
            let rebuilt = row.frame();
            let normalizer = Normalizer::with_config(
                rebuilt.provider().clone(),
                NormalizeConfig::default().with_direction(direction),
            );
            frame = Some((row.shape, row.core.clone(), rebuilt, normalizer));
        }
        let (_, _, active, normalizer) = frame.as_ref().expect("a frame was just built");

        let mut outputs: BTreeSet<String> = BTreeSet::new();
        let mut normalized_spellings = 0usize;
        for spelling in &row.spellings {
            census.outputs += 1;
            let parsed = match parse_hgvs(spelling) {
                Ok(parsed) => parsed,
                Err(_) => {
                    // A row whose own spelling does not parse is expected only
                    // where refusal IS the property.
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

            // Refusal is the property for these two kinds, so reaching an output
            // at all is the finding.
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
            //
            // Three outcomes, kept apart deliberately. Collapsing
            // `Inexpressible` into `NoSequence` published a wrong headline once:
            // an intronic output has no SPDI triple, so 381 outputs were reported
            // as "two members claiming one territory" when most had simply left
            // the transcript — a different defect with different clauses.
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
                        census.outputs_leaving_the_transcript += 1;
                        findings.push(Finding {
                            id: row.id.clone(),
                            what: format!(
                                "output left the transcript (an intronic position the input \
                                 did not name): {spelling} -> {output}"
                            ),
                        });
                    }
                    // Already counted as `unparseable_outputs` above.
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
            // Only on the AUTHORED spelling: the spanning-`delins` respelling of
            // a two-member design is a single member already, so evaluating this
            // over every spelling would count the corpus's own candidate as a
            // merge. See `Row::authored_spelling`.
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

            outputs.insert(output);
        }

        // --- rank 2: confluence, over families only --------------------------
        if row.kind == RowKind::Family {
            if normalized_spellings < 2 {
                census.underdetermined += 1;
            } else {
                match outputs.len() {
                    1 => census.converged += 1,
                    2 => census.split_two += 1,
                    3 => census.split_three += 1,
                    _ => census.split_more += 1,
                }
                if outputs.len() > 1 && divergences.len() < 12 {
                    divergences.push(Divergence {
                        id: row.id.clone(),
                        outputs: outputs.into_iter().collect(),
                    });
                }
            }
        }
    }

    Measured {
        census,
        divergences,
        findings,
    }
}

fn report(label: &str, measured: &Measured) -> String {
    let census = &measured.census;
    let mut out = format!(
        "spec conformance axis ({label})\n  \
         VALIDITY (rank 1): {} outputs, {} declined, {} unparseable, {} denoting no sequence, \
         {} leaving the transcript, {} violating an absolute prohibition\n  \
         CONFLUENCE (rank 2): converged {}, split 2 {}, split 3 {}, split 4+ {}, \
         underdetermined {}\n  \
         IDEMPOTENCY: {} outputs are not their own fixed point\n  \
         SEQUENCE PRESERVATION: {} outputs denote different bases\n  \
         REFUSAL: {} conflicting alleles accepted, {} absolute prohibitions accepted, \
         {} conditional prohibitions accepted\n  \
         NEGATIVE GUARDS: {} outputs implement rejected SVD-WG010\n",
        census.outputs,
        census.declined,
        census.unparseable_outputs,
        census.outputs_denoting_no_sequence,
        census.outputs_leaving_the_transcript,
        census.prohibition_violating_outputs,
        census.converged,
        census.split_two,
        census.split_three,
        census.split_more,
        census.underdetermined,
        census.non_idempotent_outputs,
        census.sequence_changed,
        census.conflicts_accepted,
        census.prohibited_absolute_accepted,
        census.prohibited_conditional_accepted,
        census.guard_violations,
    );
    for divergence in &measured.divergences {
        out.push_str(&format!(
            "  DIVERGES {} -> {:?}\n",
            divergence.id, divergence.outputs
        ));
    }
    // Findings are grouped, because an accepted prohibition fires once per row
    // and there are hundreds of rows per clause.
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

/// Assert one direction's census against its pin, printing the measured numbers
/// either way so a moved pin can be re-blessed from the test output.
fn assert_census(direction: ShuffleDirection, label: &str, pinned: &Census) {
    let measured = measure(direction);
    println!("{}", report(label, &measured));
    let census = &measured.census;

    // The honest-zero discipline: a property whose denominator is zero is
    // VACUOUS, not passing.
    assert!(
        census.outputs > 0,
        "{label}: VACUOUS — the corpus produced no outputs at all"
    );

    // The rank-1 counters are pinned rather than asserted at zero, because they
    // are not zero on `main` — see the module docs' table. What is asserted here
    // is that they never RISE, checked ahead of the full census so the message
    // names the property rather than dumping a struct diff.
    for (measured_value, pinned_value, what) in [
        (
            census.unparseable_outputs,
            pinned.unparseable_outputs,
            "outputs that do not re-parse",
        ),
        (
            census.outputs_denoting_no_sequence,
            pinned.outputs_denoting_no_sequence,
            "outputs that denote no sequence (two members claiming one territory)",
        ),
        (
            census.outputs_leaving_the_transcript,
            pinned.outputs_leaving_the_transcript,
            "outputs that left the transcript (an intronic position the input did not name)",
        ),
        (
            census.sequence_changed,
            pinned.sequence_changed,
            "outputs that denote different bases from their input",
        ),
        (
            census.non_idempotent_outputs,
            pinned.non_idempotent_outputs,
            "outputs that are not their own fixed point",
        ),
        (
            census.conflicts_accepted,
            pinned.conflicts_accepted,
            "conflicting alleles accepted instead of refused",
        ),
        (
            census.prohibited_absolute_accepted,
            pinned.prohibited_absolute_accepted,
            "shapes the spec calls \"not allowed\" that were accepted",
        ),
        // The two the loop was missing. Both were still caught by the full
        // census assertion below, but only as a struct diff — which is the one
        // thing this loop exists to avoid, since a reader then has to work out
        // which of twenty fields moved and whether it moved the wrong way.
        (
            census.prohibition_violating_outputs,
            pinned.prohibition_violating_outputs,
            "outputs that violate an absolute textual prohibition",
        ),
        (
            census.prohibited_conditional_accepted,
            pinned.prohibited_conditional_accepted,
            "shapes a conditional clause prohibits that were accepted",
        ),
        (
            census.guard_violations,
            pinned.guard_violations,
            "outputs implementing REJECTED SVD-WG010",
        ),
    ] {
        assert!(
            measured_value <= pinned_value,
            "{label}: {what} rose from {pinned_value} to {measured_value}. This is a rank-1 \
             conformance regression, not a representation choice.\n{}",
            report(label, &measured)
        );
    }

    assert_eq!(
        census,
        pinned,
        "{label}: the conformance census moved. Every failure figure must only ever go DOWN and \
         `converged` only ever UP; if this change lowers one, re-bless the pin in the same commit \
         and say so in the PR, and promote any newly-failing row to a named test in \
         `spec_corpus_regressions.rs`. Measured:\n{}",
        report(label, &measured)
    );
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// The corpus's own shape, before any property is read out of it.
///
/// #1460's lesson mechanized: a generator change that shrank the corpus would
/// otherwise improve every rate below while measuring less.
#[test]
fn the_corpus_has_the_shape_its_censuses_are_measured_over() {
    let built = built();
    let measured = CorpusShape::of(&built);
    println!("corpus shape: {measured:#?}");
    assert_eq!(
        measured, SHAPE,
        "the corpus shape moved; every census below is measured over it, so re-bless both \
         together and say what changed"
    );

    // Each property must have a non-zero denominator, or its census is a claim
    // about the corpus rather than about the implementation.
    assert!(measured.family_rows > 0, "VACUOUS: no family rows");
    assert!(measured.conflict_rows > 0, "VACUOUS: no conflict rows");
    assert!(
        measured.prohibited_rows > 0,
        "VACUOUS: no prohibited-shape rows"
    );
    assert!(
        measured.guarded_rows > 0,
        "VACUOUS: no rows carry a negative guard, so the rejected-SVD-WG010 check measures nothing"
    );
    assert!(
        measured.single_rows > 0,
        "VACUOUS: no single-spelling rows, so intronic and trans/unknown-phase idempotency \
         measures nothing"
    );

    // Multi-member is the priority axis. Real corpora are 592 of 9,949,738.
    let share = measured.multi_member_rows as f64 / measured.rows as f64;
    assert!(
        share > 0.9,
        "multi-member share is {share:.4}; this corpus exists to invert the natural 0.00006"
    );
}

/// Every family really is a family: two or more distinct spellings, each of
/// which denotes the row's sequence when applied independently of the
/// normalizer.
///
/// Without this the confluence figure could be a statement about a generator bug
/// — a "family" whose members denote different variants would diverge for a
/// reason that has nothing to do with the normalizer.
#[test]
fn every_family_denotes_one_sequence_independently_of_the_normalizer() {
    let built = built();
    let mut checked = 0usize;
    // Every family, not a sample. An earlier revision sampled, which is the wrong
    // trade here: the strata are enumerated in order, so any prefix or stride
    // risks skipping `scale` and `repeats` entirely — and those are the two whose
    // ground truth is hardest to get right. Measured at well under a second,
    // because the applier is a splice and a byte compare.
    for row in grouped(&built) {
        if row.kind != RowKind::Family {
            continue;
        }
        let frame = row.frame();
        let expected = row.denoted.as_deref().expect("a family denotes a sequence");
        assert!(row.spellings.len() >= 2, "{} is a singleton", row.id);
        assert_eq!(
            row.spellings.iter().collect::<BTreeSet<_>>().len(),
            row.spellings.len(),
            "{} repeats a spelling",
            row.id
        );
        for spelling in &row.spellings {
            assert_eq!(
                denoted_by(frame.provider(), frame.served(), spelling).as_deref(),
                Some(expected),
                "{}: {spelling} does not denote the row's sequence",
                row.id
            );
        }
        checked += 1;
    }
    assert!(checked > 0, "VACUOUS: no families were checked");
}

#[test]
fn three_prime_conformance_census() {
    assert_census(ShuffleDirection::ThreePrime, "3prime", &THREE_PRIME);
}

#[test]
fn five_prime_conformance_census() {
    assert_census(ShuffleDirection::FivePrime, "5prime", &FIVE_PRIME);
}
