//! The protein axis's first measurement, over the codon-designed protein corpus.
//!
//! # What this measures, and why it can measure it without a published answer
//!
//! `ferro_hgvs::conformance::protein_corpus` generates many spellings of one
//! protein change — a deletion at each equivalent position of a synonymous
//! tandem tract, a substitution beside its `delins` and parenthesised twins, a
//! cis allele in both member orders — and
//! `synthetic_protein::protein_denotation_of` confirms, independently of the
//! normalizer, that every spelling in a family denotes one protein. That makes
//! four properties checkable with no spec-published form at all:
//!
//! | rank | property | how it is checked |
//! |---|---|---|
//! | 1 | **validity** | the output re-parses and denotes a protein |
//! | 2 | **confluence** | every spelling of one protein reaches ONE output |
//! | — | **idempotency** | each output is its own fixed point |
//! | — | **protein preservation** | each output denotes the input's residues, via the oracle |
//!
//! # `Indeterminate` is its own column, on a ruling
//!
//! `rulings[confluence-gate-is-apply-equality-on-every-determined-axis]`
//! (decided 2026-08-10) states it: "Confluence is asserted only over DECIDED
//! pairs; the `Indeterminate` count is reported alongside and never folded into
//! either side." A frameshift, an extension and an `insXaa[n]` payload each fix
//! part of a protein and leave the rest open — the description genuinely does
//! not determine the molecule, so counting such a family as `converged` would
//! flatter the normalizer and counting it as divergent would invent a failure.
//! [`Census::indeterminate_families`] and its two companions therefore sit apart
//! from [`Census::converged`] and the three `split_*` figures, and the two sets
//! are never summed.
//!
//! # The first measurement — read the four findings before quoting a figure
//!
//! Over **540 rows / 1,044 spellings**, in both shuffle directions (which agree
//! byte for byte):
//!
//! | figure | value | reading |
//! |---|---|---|
//! | `unparseable_outputs` | **0** | every output re-parses |
//! | `outputs_denoting_no_protein` | **0** | no output denotes nothing |
//! | `protein_changed` | **0** | no output denotes a different protein from its input |
//! | `non_idempotent_outputs` | **0** | every output is its own fixed point |
//! | `converged` / `split_two` | **198 / 18** | over DETERMINED families only |
//! | `indeterminate_families` | **18**, all converged | its own column, per the ruling |
//! | `distinguished_preserved` | **90 of 90** | no parenthesised twin is collapsed |
//! | `conflicts_accepted` | **72 of 72** | **a rank-1 finding** |
//!
//! Four things to carry away, in descending order of how much they matter:
//!
//! 1. **`conflicts_accepted` is 72 of 72, in strict mode as well as lenient.**
//!    A reversed range, a non-flanking insertion, two allele members claiming
//!    one residue and a position past the terminator are each normalized and
//!    emitted unchanged. `rulings[absolute-prohibition-enforcement-stage]` makes
//!    lenient's behaviour arguably correct and strict's not, and says the ruling
//!    is unimplemented (#1630); `no_conflicting_description_is_refused_in_strict_mode_either`
//!    pins the protein axis's share of that.
//! 2. **The 18 divergent families are one class, and it is a contested
//!    adjudication rather than a defect** — see
//!    `the_repeat_count_and_the_explicit_deletion_of_one_tract_do_not_converge`,
//!    which is `#[ignore]`d with the question written on it.
//! 3. **Transcript geometry is inert for a `p.`-only measurement**, measured
//!    over all three exon layouts and both strands by
//!    `shape_is_inert_for_the_protein_axis`. So the corpus's phase dimension is
//!    a *guard* — it is what would catch the axis becoming geometry-sensitive —
//!    and not, today, a source of signal.
//! 4. **The four zeros above are real but narrow.** They say the normalizer does
//!    not break a protein description it accepts. They say nothing about the
//!    102 validity clauses; see the section below.
//!
//! # This is a census, not a gate
//!
//! Every figure is a **pinned baseline** measured on this branch. A change that
//! moves one must re-bless the pin in the same commit and say which way it
//! moved and why.
//!
//! # What a green run does NOT say — read this before quoting any figure
//!
//! The oracle judges **denotation**, deliberately and by design. A one-for-one
//! `delins`, a `Cys76_Cys76` range and a `p.[Arg76Ser;Cys77Trp]` that
//! `protein/delins.md:21` says should be one `delins` all denote a protein
//! perfectly well, and this census scores each of them as sound.
//!
//! The protein slice is **205 clause units** across nine files — the committed
//! rule inventory's figure, reproduced by re-deriving it from the spec checkout
//! with the inventory's own definition of a clause unit — and the great majority
//! state how a description may be **written**. A green census here is evidence
//! about none of those.
//!
//! (`synthetic_protein.rs` says "149 clauses (102 of them validity
//! requirements)". Neither number reconciles with the inventory and nothing in
//! the repository classifies a clause as a validity requirement, so that figure
//! is not restated here. See `protein_corpus`'s module docs.)
//!
//! So the rule inventory's `not-generatable` classification for
//! `docs/recommendations/protein/*` is **left in place** by this module and its
//! corpus. Lifting it needs a validity judge, which is a separate piece of work
//! and not one a denotation oracle can stand in for.

use std::collections::{BTreeMap, BTreeSet};

use ferro_hgvs::conformance::protein_corpus::{
    corpus, Determinacy, ProteinCorpus, Row, RowKind, Site,
};
use ferro_hgvs::conformance::synthetic_protein::{
    protein_denotation_of, ProteinDenotation, ProteinFrame, ProteinRefShape,
};
use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::normalize::{NormalizeConfig, ShuffleDirection};
use ferro_hgvs::{parse_hgvs, Normalizer};

// ---------------------------------------------------------------------------
// The pinned shape
// ---------------------------------------------------------------------------

/// The corpus's shape, independent of any property measured over it.
///
/// Pinned for the reason `spec_conformance_axis`'s twin is: a generator change
/// that halved the corpus would improve every rate below while measuring less,
/// which is the #1460 failure mode.
///
/// # Every figure here and below is exactly **3x** an independent measurement
///
/// State this whenever one of these numbers is quoted, because it is provable
/// rather than estimated and it changes what they mean. Two committed tests pin
/// it between them: `shape_is_inert_for_the_protein_axis` asserts that every row
/// group is built on **all three** of `ProteinRefShape::all()` — so
/// `rows == 3 * groups` — and then that the three replicas produce byte-identical
/// output. So the replicas are not three measurements; they are one measurement
/// written down three times.
///
/// Read every count accordingly: `rows: 540` is 180 designs on three geometries,
/// `conflicts_accepted: 72` is 24 distinct conflicting descriptions, and the
/// census's 18 divergent families are **6** distinct divergences. Quoting `72 of
/// 72` as seventy-two independent conflicts overstates the evidence threefold.
///
/// This is deliberate and is not a defect to fix here: the geometry dimension is
/// carried as a **guard** against `p.` output becoming splice-structure
/// dependent, not as a source of signal, and the day it stops being inert those
/// tests fail and the replication becomes real. But the multiplier is invisible
/// in the constants themselves, which is why it is written down.
const SHAPE: CorpusShape = CorpusShape {
    designs_considered: 594,
    rows: 540,
    spellings: 1_044,
    family_rows: 234,
    single_rows: 144,
    distinguished_rows: 90,
    conflict_rows: 72,
    determined_rows: 288,
    indeterminate_rows: 144,
    no_protein_rows: 36,
    undenoted_rows: 72,
};

/// The 3'-direction census, pinned.
const THREE_PRIME: Census = Census {
    measured_under: ErrorMode::Lenient,
    outputs: 1_044,
    declined: 0,
    unparseable_outputs: 0,
    outputs_denoting_no_protein: 0,
    protein_changed: 0,
    non_idempotent_outputs: 0,
    // 198 of the 216 determined families reach one output. The 18 that do not
    // are one class, one row per frame — see
    // `the_repeat_count_and_the_explicit_deletion_of_one_tract_do_not_converge`,
    // which names the open question rather than pinning an answer to it.
    converged: 198,
    split_two: 18,
    split_three: 0,
    split_more: 0,
    // The terminator substitution and its `delins` twin: both read past the
    // stop, so the oracle fixes only the prefix. Counted here and in neither
    // figure above.
    indeterminate_families: 18,
    indeterminate_converged: 18,
    indeterminate_split: 0,
    // Every parenthesised twin survives, so no evidential claim is destroyed.
    distinguished_rows: 90,
    distinguished_preserved: 90,
    distinguished_collapsed: 0,
    // 72 of 72 — every conflicting description is normalized rather than
    // refused. A rank-1 finding, and the same figure in strict mode; see
    // `no_conflicting_description_is_refused_in_strict_mode_either`.
    conflicts_accepted: 72,
};

/// The 5'-direction census, pinned.
///
/// Measured in full rather than assumed identical: confluence is a property of
/// the normalizer, not of one shuffle direction, so a figure that appears at 3'
/// and not at 5' is a claim about the code.
///
/// **It is byte-identical to [`THREE_PRIME`], and that is a result rather than
/// a copy.** `general.md:43` makes the 3' rule reach protein descriptions, so a
/// direction-sensitive protein normalizer is conceivable; measured, this one is
/// not. The two constants are written out separately so the day one of them
/// moves, the other does not move with it.
const FIVE_PRIME: Census = Census {
    measured_under: ErrorMode::Lenient,
    outputs: 1_044,
    declined: 0,
    unparseable_outputs: 0,
    outputs_denoting_no_protein: 0,
    protein_changed: 0,
    non_idempotent_outputs: 0,
    // 198 of the 216 determined families reach one output. The 18 that do not
    // are one class, one row per frame — see
    // `the_repeat_count_and_the_explicit_deletion_of_one_tract_do_not_converge`,
    // which names the open question rather than pinning an answer to it.
    converged: 198,
    split_two: 18,
    split_three: 0,
    split_more: 0,
    // The terminator substitution and its `delins` twin: both read past the
    // stop, so the oracle fixes only the prefix. Counted here and in neither
    // figure above.
    indeterminate_families: 18,
    indeterminate_converged: 18,
    indeterminate_split: 0,
    // Every parenthesised twin survives, so no evidential claim is destroyed.
    distinguished_rows: 90,
    distinguished_preserved: 90,
    distinguished_collapsed: 0,
    // 72 of 72 — every conflicting description is normalized rather than
    // refused. A rank-1 finding, and the same figure in strict mode; see
    // `no_conflicting_description_is_refused_in_strict_mode_either`.
    conflicts_accepted: 72,
};

#[derive(Debug, Default, PartialEq, Eq)]
struct CorpusShape {
    designs_considered: usize,
    rows: usize,
    spellings: usize,
    family_rows: usize,
    single_rows: usize,
    distinguished_rows: usize,
    conflict_rows: usize,
    determined_rows: usize,
    indeterminate_rows: usize,
    no_protein_rows: usize,
    undenoted_rows: usize,
}

impl CorpusShape {
    fn of(built: &ProteinCorpus) -> Self {
        let by_kind = built.by_kind();
        let by_determinacy = built.by_determinacy();
        let count = |label: &str| by_determinacy.get(label).copied().unwrap_or(0);
        Self {
            designs_considered: built.designs_considered,
            rows: built.rows.len(),
            spellings: built.spellings(),
            family_rows: by_kind.get(&RowKind::Family).copied().unwrap_or(0),
            single_rows: by_kind.get(&RowKind::Single).copied().unwrap_or(0),
            distinguished_rows: by_kind.get(&RowKind::Distinguished).copied().unwrap_or(0),
            conflict_rows: by_kind.get(&RowKind::Conflict).copied().unwrap_or(0),
            determined_rows: count(Determinacy::Determined.label()),
            indeterminate_rows: count(Determinacy::Indeterminate.label()),
            no_protein_rows: count(Determinacy::NoProtein.label()),
            undenoted_rows: count(Determinacy::Undenoted.label()),
        }
    }
}

// ---------------------------------------------------------------------------
// The census
// ---------------------------------------------------------------------------

/// One direction's measurement.
///
/// The `indeterminate_*` block is deliberately NOT summable with the confluence
/// block above it — see the module docs and the ruling they cite.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct Census {
    /// Which error mode produced these figures. Recorded because refusal
    /// counters are mode-dependent and a census with no mode on it invites
    /// being compared against one measured under another.
    measured_under: ErrorMode,

    // -- validity (rank 1) --
    /// Spellings the census attempted to normalize.
    outputs: usize,
    /// Spellings that did not parse or whose normalization failed, where that
    /// was not the row's property.
    declined: usize,
    /// Outputs `parse_hgvs` cannot read back.
    unparseable_outputs: usize,
    /// Outputs the oracle says denote no single protein — the protein twin of
    /// `outputs_denoting_no_sequence`.
    outputs_denoting_no_protein: usize,

    // -- protein preservation --
    /// Outputs denoting a different protein from their input, counted **per
    /// spelling** over rows whose input denotation is
    /// [`Determinacy::Determined`]. So its denominator is those rows' spellings,
    /// not `determined_rows` — a determined family contributes once per
    /// spelling, and reading this against the row count would understate it.
    protein_changed: usize,

    // -- idempotency --
    /// Outputs that are not their own fixed point.
    non_idempotent_outputs: usize,

    // -- confluence (rank 2), over DETERMINED families only --
    /// Families whose spellings all reached one output.
    converged: usize,
    /// Families reaching exactly two distinct outputs.
    split_two: usize,
    /// Three.
    split_three: usize,
    /// Four or more.
    split_more: usize,

    // -- the Indeterminate column: reported alongside, never folded --
    /// Families whose denotation is [`ProteinDenotation::Indeterminate`].
    indeterminate_families: usize,
    /// …of which reached one output.
    indeterminate_converged: usize,
    /// …of which did not.
    indeterminate_split: usize,

    // -- distinguished pairs: apply-equal, epistemically different --
    /// Rows whose spellings denote one protein and must NOT converge.
    distinguished_rows: usize,
    /// …whose spellings stayed distinct, which is the wanted answer.
    distinguished_preserved: usize,
    /// …that were collapsed onto one output, destroying the evidential claim
    /// `protein/substitution.md:16` attaches to the parentheses.
    distinguished_collapsed: usize,

    // -- refusal --
    /// Conflicting descriptions the normalizer produced an output for.
    conflicts_accepted: usize,
}

/// A family that did not converge, for the report.
#[derive(Debug, Clone)]
struct Divergence {
    id: String,
    determinacy: &'static str,
    outputs: Vec<String>,
}

/// Anything worth naming in the report that is not a counter.
#[derive(Debug, Clone)]
struct Finding {
    id: String,
    what: String,
}

#[derive(Debug, Default)]
struct Measured {
    census: Census,
    divergences: Vec<Divergence>,
    findings: Vec<Finding>,
}

/// How many divergences the report prints before eliding.
const MAX_DIVERGENCES: usize = 12;

impl Census {
    fn absorb(&mut self, other: &Census) {
        self.outputs += other.outputs;
        self.declined += other.declined;
        self.unparseable_outputs += other.unparseable_outputs;
        self.outputs_denoting_no_protein += other.outputs_denoting_no_protein;
        self.protein_changed += other.protein_changed;
        self.non_idempotent_outputs += other.non_idempotent_outputs;
        self.converged += other.converged;
        self.split_two += other.split_two;
        self.split_three += other.split_three;
        self.split_more += other.split_more;
        self.indeterminate_families += other.indeterminate_families;
        self.indeterminate_converged += other.indeterminate_converged;
        self.indeterminate_split += other.indeterminate_split;
        self.distinguished_rows += other.distinguished_rows;
        self.distinguished_preserved += other.distinguished_preserved;
        self.distinguished_collapsed += other.distinguished_collapsed;
        self.conflicts_accepted += other.conflicts_accepted;
    }
}

fn built() -> ProteinCorpus {
    corpus()
}

/// Rows sorted so every row sharing one frame is adjacent.
fn grouped(built: &ProteinCorpus) -> Vec<&Row> {
    let mut rows: Vec<&Row> = built.rows.iter().collect();
    rows.sort_by(|a, b| {
        (a.shape.label(), a.design, &a.id).cmp(&(b.shape.label(), b.design, &b.id))
    });
    rows
}

/// The contiguous runs of `rows` sharing one frame.
fn frame_groups<'a>(rows: &'a [&'a Row]) -> Vec<&'a [&'a Row]> {
    let mut groups = Vec::new();
    let mut start = 0usize;
    while start < rows.len() {
        let mut end = start + 1;
        while end < rows.len() && rows[end].frame_key() == rows[start].frame_key() {
            end += 1;
        }
        groups.push(&rows[start..end]);
        start = end;
    }
    groups
}

fn measurement_config(direction: ShuffleDirection) -> NormalizeConfig {
    NormalizeConfig::default().with_direction(direction)
}

/// Normalize one frame's rows and census them.
fn measure_group(group: &[&Row], direction: ShuffleDirection) -> Measured {
    let mut census = Census {
        measured_under: measurement_config(direction).error_config.mode,
        ..Census::default()
    };
    let mut divergences: Vec<Divergence> = Vec::new();
    let mut findings: Vec<Finding> = Vec::new();

    let frame = group[0].frame();
    let normalizer =
        Normalizer::with_config(frame.provider().clone(), measurement_config(direction));

    for row in group {
        let mut outputs: BTreeSet<String> = BTreeSet::new();
        for spelling in &row.spellings {
            census.outputs += 1;
            let Ok(parsed) = parse_hgvs(spelling) else {
                // A corpus spelling that does not parse is expected only where
                // refusal IS the row's property.
                if row.kind != RowKind::Conflict {
                    census.declined += 1;
                }
                continue;
            };
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                normalizer.normalize(&parsed)
            }));
            let Ok(Ok(value)) = result else {
                if row.kind != RowKind::Conflict {
                    census.declined += 1;
                }
                continue;
            };
            let output = value.to_string();

            if row.kind == RowKind::Conflict {
                census.conflicts_accepted += 1;
                findings.push(Finding {
                    id: row.id.clone(),
                    what: format!("conflicting description accepted: {spelling} -> {output}"),
                });
            }

            // --- rank 1: validity -----------------------------------------
            if parse_hgvs(&output).is_err() {
                census.unparseable_outputs += 1;
                findings.push(Finding {
                    id: row.id.clone(),
                    what: format!("output does not re-parse: {spelling} -> {output}"),
                });
                continue;
            }
            let output_denotation = protein_denotation_of(frame.provider(), &output);
            if matches!(Determinacy::of(&output_denotation), Determinacy::Undenoted)
                && row.kind != RowKind::Conflict
            {
                census.outputs_denoting_no_protein += 1;
                findings.push(Finding {
                    id: row.id.clone(),
                    what: format!(
                        "output denotes no single protein: {spelling} -> {output} \
                         ({output_denotation:?})"
                    ),
                });
            }

            // --- protein preservation --------------------------------------
            if row.determinacy() == Determinacy::Determined
                && output_denotation != row.denoted
                && row.kind != RowKind::Conflict
            {
                census.protein_changed += 1;
                findings.push(Finding {
                    id: row.id.clone(),
                    what: format!("output denotes a different protein: {spelling} -> {output}"),
                });
            }

            // --- idempotency ----------------------------------------------
            if let Ok(reparsed) = parse_hgvs(&output) {
                let again = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    normalizer.normalize(&reparsed)
                }));
                if let Ok(Ok(second)) = again {
                    if second.to_string() != output {
                        census.non_idempotent_outputs += 1;
                        findings.push(Finding {
                            id: row.id.clone(),
                            what: format!("output is not a fixed point: {output} -> {}", second),
                        });
                    }
                }
            }

            outputs.insert(output);
        }

        // --- distinguished pairs -------------------------------------------
        //
        // Measured BEFORE the confluence block and never inside it: the
        // wanted answer here is more than one output, which is the opposite of
        // what `converged` rewards.
        if row.kind == RowKind::Distinguished {
            census.distinguished_rows += 1;
            if outputs.len() >= row.spellings.len() {
                census.distinguished_preserved += 1;
            } else {
                census.distinguished_collapsed += 1;
                findings.push(Finding {
                    id: row.id.clone(),
                    what: format!(
                        "distinguished spellings collapsed: {:?} -> {:?}",
                        row.spellings, outputs
                    ),
                });
            }
            continue;
        }

        // --- rank 2: confluence, over families only ------------------------
        if row.kind != RowKind::Family {
            continue;
        }
        let arity = outputs.len();
        let determinacy = row.determinacy();
        // The ruling's split: DETERMINED families are the confluence
        // denominator; INDETERMINATE ones are counted beside it and never in it.
        match determinacy {
            Determinacy::Indeterminate => {
                census.indeterminate_families += 1;
                if arity <= 1 {
                    census.indeterminate_converged += 1;
                } else {
                    census.indeterminate_split += 1;
                    divergences.push(Divergence {
                        id: row.id.clone(),
                        determinacy: determinacy.label(),
                        outputs: outputs.iter().cloned().collect(),
                    });
                }
            }
            Determinacy::Determined | Determinacy::NoProtein => {
                // Arity 0 — every spelling declined — is folded into
                // `converged` rather than given an arm, and that is only safe
                // because each of those declines already incremented
                // `declined`, which both censuses pin at 0. So a family that
                // converged by producing nothing cannot reach this figure
                // without failing the census first. If `declined` ever stops
                // being pinned at zero, split this arm out.
                match arity {
                    0 | 1 => census.converged += 1,
                    2 => census.split_two += 1,
                    3 => census.split_three += 1,
                    _ => census.split_more += 1,
                }
                if arity > 1 {
                    divergences.push(Divergence {
                        id: row.id.clone(),
                        determinacy: determinacy.label(),
                        outputs: outputs.iter().cloned().collect(),
                    });
                }
            }
            // A family whose authored spelling denotes nothing is a
            // `RowKind::Conflict`, so this arm is unreachable through `corpus()`.
            // It is written out rather than folded into either branch above,
            // because folding it would put an undenoted family into a
            // confluence figure the ruling scopes to decided pairs.
            Determinacy::Undenoted => {}
        }
    }

    Measured {
        census,
        divergences,
        findings,
    }
}

fn measure(direction: ShuffleDirection) -> Measured {
    let built = built();
    let rows = grouped(&built);
    let mut total = Measured {
        census: Census {
            measured_under: measurement_config(direction).error_config.mode,
            ..Census::default()
        },
        ..Measured::default()
    };
    for group in frame_groups(&rows) {
        let measured = measure_group(group, direction);
        total.census.absorb(&measured.census);
        total.divergences.extend(measured.divergences);
        total.findings.extend(measured.findings);
    }
    total.divergences.sort_by(|a, b| a.id.cmp(&b.id));
    total.findings.sort_by(|a, b| a.id.cmp(&b.id));
    total
}

fn report(label: &str, measured: &Measured) -> String {
    use std::fmt::Write;
    let census = &measured.census;
    let mut out = String::new();
    let _ = writeln!(out, "== protein conformance census: {label} ==");
    let _ = writeln!(out, "measured under: {:?}", census.measured_under);
    let _ = writeln!(out, "-- validity --");
    let _ = writeln!(out, "outputs                     {}", census.outputs);
    let _ = writeln!(out, "declined                    {}", census.declined);
    let _ = writeln!(
        out,
        "unparseable_outputs         {}",
        census.unparseable_outputs
    );
    let _ = writeln!(
        out,
        "outputs_denoting_no_protein {}",
        census.outputs_denoting_no_protein
    );
    let _ = writeln!(out, "-- protein preservation --");
    let _ = writeln!(
        out,
        "protein_changed             {}",
        census.protein_changed
    );
    let _ = writeln!(out, "-- idempotency --");
    let _ = writeln!(
        out,
        "non_idempotent_outputs      {}",
        census.non_idempotent_outputs
    );
    let _ = writeln!(out, "-- confluence (DETERMINED families only) --");
    let _ = writeln!(out, "converged                   {}", census.converged);
    let _ = writeln!(out, "split_two                   {}", census.split_two);
    let _ = writeln!(out, "split_three                 {}", census.split_three);
    let _ = writeln!(out, "split_more                  {}", census.split_more);
    let _ = writeln!(
        out,
        "-- INDETERMINATE (its own column; never folded into the block above) --"
    );
    let _ = writeln!(
        out,
        "indeterminate_families      {}",
        census.indeterminate_families
    );
    let _ = writeln!(
        out,
        "indeterminate_converged     {}",
        census.indeterminate_converged
    );
    let _ = writeln!(
        out,
        "indeterminate_split         {}",
        census.indeterminate_split
    );
    let _ = writeln!(
        out,
        "-- DISTINGUISHED (apply-equal, epistemically different; must NOT converge) --"
    );
    let _ = writeln!(
        out,
        "distinguished_rows          {}",
        census.distinguished_rows
    );
    let _ = writeln!(
        out,
        "distinguished_preserved     {}",
        census.distinguished_preserved
    );
    let _ = writeln!(
        out,
        "distinguished_collapsed     {}",
        census.distinguished_collapsed
    );
    let _ = writeln!(out, "-- refusal --");
    let _ = writeln!(
        out,
        "conflicts_accepted          {}",
        census.conflicts_accepted
    );
    if !measured.divergences.is_empty() {
        let _ = writeln!(out, "-- divergences --");
        for divergence in measured.divergences.iter().take(MAX_DIVERGENCES) {
            let _ = writeln!(
                out,
                "  {} [{}] -> {:?}",
                divergence.id, divergence.determinacy, divergence.outputs
            );
        }
        if measured.divergences.len() > MAX_DIVERGENCES {
            let _ = writeln!(
                out,
                "  … and {} more",
                measured.divergences.len() - MAX_DIVERGENCES
            );
        }
    }
    if !measured.findings.is_empty() {
        let _ = writeln!(out, "-- findings --");
        for finding in measured.findings.iter().take(MAX_DIVERGENCES) {
            let _ = writeln!(out, "  {}: {}", finding.id, finding.what);
        }
        if measured.findings.len() > MAX_DIVERGENCES {
            let _ = writeln!(
                out,
                "  … and {} more",
                measured.findings.len() - MAX_DIVERGENCES
            );
        }
    }
    out
}

fn assert_census(direction: ShuffleDirection, label: &str, pinned: &Census) {
    let measured = measure(direction);
    assert_eq!(
        &measured.census,
        pinned,
        "the {label} protein census moved.\n{}",
        report(label, &measured)
    );
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn the_corpus_has_the_shape_its_censuses_are_measured_over() {
    let built = built();
    let shape = CorpusShape::of(&built);
    assert_eq!(
        shape,
        SHAPE,
        "the protein corpus's shape moved, so every figure below is measured over a \
         different population:\n  by kind      {:?}\n  by stratum   {:?}\n  by site      {:?}\n  \
         by edit      {:?}\n  by shape     {:?}\n  by phase     {:?}\n  by determinacy {:?}\n  \
         drops        {:?}",
        built.by_kind(),
        built.by_stratum(),
        built.by_site(),
        built.by_edit(),
        built.by_shape(),
        built.by_junction_phase(),
        built.by_determinacy(),
        built.drops,
    );
    // Naming the drops is part of the shape: a design that stops being
    // generated is invisible in a bare count.
    let dropped: Vec<String> = built
        .dropped
        .iter()
        .map(|(id, reason)| format!("{id} [{}]", reason.label()))
        .collect();
    assert_eq!(
        dropped.len(),
        built.designs_considered - built.rows.len(),
        "the drop list and the drop count disagree: {dropped:?}"
    );
}

/// Every family's spellings denote one protein, established through the oracle
/// and not through the normalizer.
///
/// Without this the confluence figures would be unattributable: a family whose
/// spellings denote different proteins *should* reach different outputs, so a
/// divergence would say nothing about the normalizer.
#[test]
fn every_family_denotes_one_protein_independently_of_the_normalizer() {
    let built = built();
    let mut checked = 0usize;
    for row in &built.rows {
        if row.kind != RowKind::Family {
            continue;
        }
        let frame = row.frame();
        for spelling in &row.spellings {
            assert_eq!(
                protein_denotation_of(frame.provider(), spelling),
                row.denoted,
                "{}: {spelling} does not denote the row's protein",
                row.id
            );
            checked += 1;
        }
    }
    assert!(checked > 0, "no family spellings were checked at all");
}

/// The corpus reaches every determinacy class, so the `Indeterminate` column
/// has a non-empty denominator.
///
/// `0 of 0` is what a census with a vacuous column looks like, and it reads
/// identically to a clean one.
#[test]
fn the_indeterminate_column_has_a_non_empty_denominator() {
    let built = built();
    let by_determinacy = built.by_determinacy();
    let indeterminate = by_determinacy
        .get(Determinacy::Indeterminate.label())
        .copied()
        .unwrap_or(0);
    assert!(
        indeterminate > 0,
        "no indeterminate rows, so the column reports 0 of 0: {by_determinacy:?}"
    );
}

/// **Is transcript geometry inert for a `p.`-only measurement?**
///
/// `tests/it/protein_axis_split_move.rs` argues that "a protein description has
/// no codons of its own", from which it follows that the exon layout under a
/// protein should not change what a `p.` description normalizes to. That is a
/// claim about the axis, and this test measures it rather than assuming it —
/// assuming it would be another entry on `CLAUDE.md`'s structural-blindness
/// inventory, since the corpus varies phase and strand precisely so a
/// junction-sensitive defect *could* show.
///
/// The assertion is on the measured verdict, so if the answer ever changes the
/// test fails and the corpus's geometry dimension stops being decoration.
#[test]
fn shape_is_inert_for_the_protein_axis() {
    let built = built();
    let mut outputs_by_design: BTreeMap<String, BTreeMap<&'static str, Vec<String>>> =
        BTreeMap::new();
    for group in frame_groups(&grouped(&built)) {
        let frame = group[0].frame();
        let normalizer = Normalizer::with_config(
            frame.provider().clone(),
            measurement_config(ShuffleDirection::ThreePrime),
        );
        for row in group {
            // The row id is `<shape>-<rest>`; the rest is the design.
            let (shape_label, rest) = row.id.split_once('-').expect("row ids carry a shape");
            let rendered: Vec<String> = row
                .spellings
                .iter()
                .map(|spelling| match parse_hgvs(spelling) {
                    Ok(parsed) => std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        normalizer.normalize(&parsed)
                    }))
                    .ok()
                    .and_then(Result::ok)
                    .map_or_else(|| "<declined>".to_string(), |value| value.to_string()),
                    Err(_) => "<unparseable>".to_string(),
                })
                .collect();
            let entry = outputs_by_design.entry(rest.to_string()).or_default();
            let label = ProteinRefShape::all()
                .into_iter()
                .find(|shape| shape.label() == shape_label)
                .expect("a row id names a known shape")
                .label();
            entry.insert(label, rendered);
        }
    }

    let mut differing: Vec<&String> = Vec::new();
    for (design, per_shape) in &outputs_by_design {
        let distinct: BTreeSet<&Vec<String>> = per_shape.values().collect();
        if distinct.len() > 1 {
            differing.push(design);
        }
    }
    assert!(
        outputs_by_design.len() > 1,
        "nothing to compare across geometries"
    );
    // …and that is a count of GROUPS, not of comparisons. A group holding one
    // geometry cannot disagree with itself, so it contributes agreement it never
    // demonstrated — the `0 of 0` shape `spec_conformance_axis` already guards
    // its own denominators against. Every row is built on every shape today, so
    // the contract is the full set, not merely "more than one".
    let geometries = ProteinRefShape::all().len();
    let short: Vec<(&String, usize)> = outputs_by_design
        .iter()
        .filter(|(_, per_shape)| per_shape.len() != geometries)
        .map(|(design, per_shape)| (design, per_shape.len()))
        .collect();
    assert!(
        short.is_empty(),
        "{} of {} row groups were measured on fewer than all {geometries} geometries, so those \
         rows were never compared across shapes at all and their inertness is unmeasured. \
         Starting with {:?}. Fix the generator or scope the claim — do not relax this, because a \
         group with one shape reads exactly like a group that agreed.",
        short.len(),
        outputs_by_design.len(),
        short.iter().take(5).collect::<Vec<_>>()
    );
    assert!(
        differing.is_empty(),
        "transcript geometry is NOT inert for the protein axis — {} of {} designs differ across \
         the three exon layouts, starting with {:?}. That is a finding, not a test to relax: it \
         means a `p.` output depends on the splice structure under it, and the corpus's phase \
         dimension has become load-bearing.",
        differing.len(),
        outputs_by_design.len(),
        differing.iter().take(5).collect::<Vec<_>>()
    );
}

/// A conflicting description denotes no protein, so a census that counted zero
/// of them would be measuring nothing.
#[test]
fn the_conflict_stratum_denotes_nothing_before_the_normalizer_is_asked() {
    let built = built();
    let mut conflicts = 0usize;
    for row in &built.rows {
        if row.kind != RowKind::Conflict {
            continue;
        }
        conflicts += 1;
        assert!(
            matches!(row.denoted, ProteinDenotation::NoSingleSequence(_)),
            "{} is filed as a conflict but the oracle says {:?}",
            row.id,
            row.denoted
        );
    }
    assert!(conflicts > 0, "the corpus holds no conflicting rows");
}

/// Every site the corpus declares is reachable on a real frame, including the
/// two a nucleotide corpus has no analogue of.
#[test]
fn the_terminator_and_past_end_sites_resolve_where_the_corpus_says_they_do() {
    let frame = ProteinFrame::from_cds(
        ProteinRefShape::SingleExon,
        &ferro_hgvs::conformance::protein_corpus::codon_core(
            ferro_hgvs::conformance::protein_corpus::all_designs()[0],
        ),
    )
    .expect("the first design is a CDS");
    let residues = frame.residues().len() as u64;
    assert_eq!(Site::Terminator.position(&frame), residues + 1);
    assert_eq!(Site::PastEnd.position(&frame), residues + 2);
    assert_eq!(
        frame.residues_with_stop().len(),
        residues as usize + 1,
        "the terminator is addressable"
    );
}

/// **CONTESTED — an open question, deliberately left failing and `#[ignore]`d.**
///
/// `p.Ala50[2]` and `p.Ala52_Ala53del` denote one protein: the reference tract
/// holds four `Ala`, so stating two copies and deleting two residues leave the
/// same molecule. Ferro preserves both spellings, so the pair is the whole of
/// this census's `split_two`.
///
/// # Why it is not obviously a defect
///
/// `protein/repeated.md:22` states a preference and **conditions it on
/// provenance ferro cannot recover**:
///
/// > "**NOTE**: when the repeat is variable in the population and the reference
/// > sequence has 10 units, the description `p.Ala2[9]` is preferred over
/// > `p.Ala11del`."
///
/// "when the repeat is variable in the population" is a fact about a population,
/// not about the reference or the description, so a normalizer holding only a
/// protein sequence cannot evaluate the antecedent. That is the same shape
/// `rulings[separation-is-a-property-of-the-spelling-not-of-the-variant]`
/// records for `delins.md:79-84`'s provenance rationale, and the shape
/// `rulings[confluence-gate-is-apply-equality-on-every-determined-axis]` warns
/// about in terms: "two descriptions can be apply-equal on every axis and still
/// make different epistemic claims, which is a canonical-form question and must
/// not be encoded as a rung."
///
/// # What is actually open
///
/// Whether `[n]` is ferro's canonical form for a change wholly inside a tandem
/// tract, unconditionally — which would satisfy `:22`'s consequent while
/// ignoring its antecedent — or whether the two spellings stay distinct because
/// the antecedent is unrecoverable. That is an operator decision under the
/// escalation register in `rulings[adjudication-precedence-order]`, and it is
/// not one this test may make. It is written as an assertion so the question is
/// executable the moment it is answered, and `#[ignore]`d so it does not assert
/// an answer nobody has given.
#[test]
#[ignore = "open question: protein/repeated.md:22 prefers the [n] form only under a \
            population condition ferro cannot evaluate; needs an operator ruling"]
fn the_repeat_count_and_the_explicit_deletion_of_one_tract_do_not_converge() {
    let measured = measure(ShuffleDirection::ThreePrime);
    assert_eq!(
        measured.census.split_two,
        0,
        "the repeat-count and explicit-deletion spellings of one tract still \
         diverge:\n{}",
        report("3'", &measured)
    );
}

/// **A rank-1 finding, measured rather than argued.**
///
/// Every conflicting description in the corpus is normalized in lenient mode —
/// 72 of 72 — and the four shapes are a reversed range
/// (`protein/deletion.md:18`), a non-flanking insertion
/// (`protein/insertion.md:17-18`), two allele members claiming one residue, and
/// a position past the terminator.
///
/// `rulings[absolute-prohibition-enforcement-stage]` (decided 2026-08-10) makes
/// the enforcement stage mode-dependent: "strict fails at PARSE (strict
/// validates input conformance, not merely parseability); lenient does not
/// validate input conformance and fails only when it cannot NORMALIZE". So
/// lenient's 72 is arguably the decided behaviour and strict's is not — which
/// makes strict the measurement that matters, and it is taken here rather than
/// inferred.
///
/// The record itself says the ruling is **not implemented** (#1630), so this
/// pins the protein axis's share of that gap. The figure is asserted equal to
/// the lenient one, so the day strict starts refusing, this fails and the gap is
/// closed on the record.
#[test]
fn no_conflicting_description_is_refused_in_strict_mode_either() {
    let built = built();
    let mut conflicts = 0usize;
    let mut accepted = 0usize;
    let mut examples: Vec<String> = Vec::new();
    for group in frame_groups(&grouped(&built)) {
        let frame = group[0].frame();
        let normalizer = Normalizer::with_config(
            frame.provider().clone(),
            NormalizeConfig::default()
                .with_direction(ShuffleDirection::ThreePrime)
                .with_error_mode(ErrorMode::Strict),
        );
        for row in group {
            if row.kind != RowKind::Conflict {
                continue;
            }
            conflicts += 1;
            for spelling in &row.spellings {
                let Ok(parsed) = parse_hgvs(spelling) else {
                    continue;
                };
                if let Ok(value) = normalizer.normalize(&parsed) {
                    accepted += 1;
                    if examples.len() < 4 {
                        examples.push(format!("{spelling} -> {value}"));
                    }
                }
            }
        }
    }
    assert!(conflicts > 0, "no conflicting rows to measure");
    assert_eq!(
        accepted, conflicts,
        "strict mode's refusal rate on the protein axis moved — {accepted} of {conflicts} \
         accepted, examples {examples:?}. If this is now LOWER, #1630 has progressed and the \
         pin should come down with a note saying which shape strict began refusing."
    );
}

#[test]
fn three_prime_conformance_census() {
    assert_census(ShuffleDirection::ThreePrime, "3'", &THREE_PRIME);
}

#[test]
fn five_prime_conformance_census() {
    assert_census(ShuffleDirection::FivePrime, "5'", &FIVE_PRIME);
}
