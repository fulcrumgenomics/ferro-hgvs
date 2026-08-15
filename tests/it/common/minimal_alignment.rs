//! Enumerate every minimal alignment of a `(reference, alternate)` block, and
//! report which reference bases are matched in **all** of them.
//!
//! # What this is for
//!
//! `rulings[unchanged-is-read-over-every-minimal-alignment]` (in
//! `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`) decides
//! that
//!
//! > a reference base is unchanged iff it is matched in EVERY MINIMAL
//! > alignment of the block.
//!
//! Until this module existed, that rule was applied by reading and by hand
//! argument: nothing in the repo could *compute* it. An intersection over a set
//! nobody enumerates is not a checkable claim, and the arithmetic has been got
//! wrong more than once. [`minimal_alignments`] enumerates the set and takes
//! the intersection, so "is this base unchanged?" has an answer a test can
//! assert on.
//!
//! # This is an instrument, not a normalizer input
//!
//! Nothing under `src/` calls this and nothing should. It lives in
//! `tests/it/common/` for three reasons:
//!
//! - **No production consumer.** Normalization decides its own partition
//!   through `src/normalize/seqfirst/align.rs`. Adding a second entry point to
//!   `src/` that only tests call would be crate surface kept alive by an
//!   `#[allow(dead_code)]`.
//! - **A deliberately different notion sits next door.** That module's
//!   `Dominators` asks "is one `(ref, alt)` **cell** a match on every minimal
//!   path?", which is strictly stronger than the record's question and,
//!   per the record, *not a bug* — the two notions differ on purpose. Putting
//!   an exponential column-based enumerator beside the linear cell-based one
//!   would invite exactly the confusion the record was written to end.
//! - **Exponential by construction.** The optimal-alignment count of an
//!   all-mismatch block grows by a factor approaching 5.83 per base, so this is
//!   unfit for a hot path however it were written. See
//!   [`DEFAULT_ALIGNMENT_CAP`].
//!
//! # The closed form this checks against
//!
//! `issue_1539_split_member_separation::forced_unchanged_columns` computes the
//! same set in `Θ(n·m)` without materialising an alignment, and is the detector
//! the record names as where the ruling came from. It is a gate, so it needs to
//! be fast; this is an instrument, so it needs to be *inspectable* — it can say
//! how many alignments there are and what each of them matches, which a
//! closed form cannot. The two are kept as separate implementations on purpose
//! so each can check the other, and
//! `minimal_alignment_enumeration_proptest::the_closed_form_and_the_enumerator_agree`
//! asserts they do.
//!
//! # The record does not state its cost model
//!
//! This is the gap the module makes visible. Two things point at substitution
//! cost 1 — the record's own worked example (`CAG -> AGA`, "position-wise:
//! cost 3" over three substitutions), and the closed-form detector above, which
//! hardcodes it — but neither is a ruling, and the two plausible models
//! **disagree**: on the `1420-v4` block (`ATGC -> GCTG`)
//! [`CostModel::Levenshtein`] calls two reference bases unchanged while
//! [`CostModel::SubstitutionCostsTwo`] calls one. So the model is a required
//! parameter here rather than a default, and
//! `minimal_alignment_enumeration::the_two_cost_models_disagree_on_1420_v4`
//! pins both answers without picking a winner.

/// One step of an alignment: the label on one column.
///
/// `Ord` is derived so a caller can sort alignments to check for duplicates;
/// the order is the declaration order and carries no meaning of its own.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Step {
    /// Reference and alternate bases agree; consumes one of each.
    Match,
    /// Reference and alternate bases differ; consumes one of each.
    Sub,
    /// Consumes a reference base with no alternate base.
    Del,
    /// Consumes an alternate base with no reference base.
    Ins,
}

/// How much each column costs.
///
/// Both models price a match at 0 and each indel at 1; they differ only on a
/// substitution. **Neither is "the" model** — see the module docs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CostModel {
    /// Unit-cost Levenshtein: a substitution costs 1.
    ///
    /// This is what the record's own worked example implies (`CAG -> AGA`,
    /// "position-wise: cost 3" over three substitutions), and it is the model
    /// `src/normalize/seqfirst/align.rs` builds its DAG under.
    Levenshtein,
    /// A substitution costs 2 — exactly what deleting the reference base and
    /// inserting the alternate one costs.
    ///
    /// Under this model a `Sub` column is never *cheaper* than the `Del`+`Ins`
    /// pair that denotes the same thing, only tied. The tie is kept rather than
    /// removed, because dropping the `Sub` step is a third model and the point
    /// of this module is to make the choice visible instead of folding models
    /// together. Keeping it changes the **count** — `ATGC -> GCTG` enumerates 6
    /// optimal alignments with the tied `Sub` step and 4 without it — but
    /// provably not the **unchanged set**: swapping a `Sub` column for the
    /// `Del`+`Ins` pair preserves both the cost and the set of matched
    /// reference offsets, so the intersection is identical either way.
    SubstitutionCostsTwo,
}

impl CostModel {
    /// Every model, so a caller can sweep them rather than naming one.
    pub const ALL: [CostModel; 2] = [CostModel::Levenshtein, CostModel::SubstitutionCostsTwo];

    /// The cost of a column pairing two unequal bases.
    pub fn substitution_cost(self) -> u32 {
        match self {
            CostModel::Levenshtein => 1,
            CostModel::SubstitutionCostsTwo => 2,
        }
    }

    /// The cost of `step`, under this model.
    pub fn step_cost(self, step: Step) -> u32 {
        match step {
            Step::Match => 0,
            Step::Sub => self.substitution_cost(),
            Step::Del | Step::Ins => 1,
        }
    }
}

/// The most optimal alignments [`minimal_alignments`] will enumerate before
/// refusing.
///
/// # Why a cap at all, and why it is an error rather than a truncation
///
/// The unchanged set is an **intersection**, so dropping alignments can only
/// *add* offsets to it. A silently truncated enumeration would therefore be
/// wrong in the unsafe direction — bases would look unchanged that are not,
/// which is precisely the misreading the record exists to prevent. So
/// exceeding the cap returns [`CapExceeded`] and no set at all.
///
/// # Why 65536
///
/// The count grows exponentially, so the cap's job is to fail loudly rather
/// than to be generous — one more base costs roughly six times the budget, and
/// four more costs a thousand times, so no reachable value changes which blocks
/// this can answer for. Measured, on an all-mismatch block `A^n -> C^n` under
/// [`CostModel::SubstitutionCostsTwo`] (the central Delannoy numbers, the worst
/// case at each length):
///
/// | `n` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
/// |---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
/// | optimal alignments | 3 | 13 | 63 | 321 | 1683 | 8989 | 48639 | 265729 | 1462563 |
///
/// 65536 is the largest round value that still fits a 7-base all-mismatch
/// block, and enumerating that many costs a few megabytes and a few
/// milliseconds — small enough that the cap never fires on the blocks this
/// instrument is for. The record's worked examples and the reported
/// `1420`/`1421` blocks are 3-5 bases and enumerate in single digits.
///
/// A caller that needs a different bound names it with
/// [`minimal_alignments_capped`], which is also how the refusal itself is
/// tested.
pub const DEFAULT_ALIGNMENT_CAP: usize = 65_536;

/// Enumeration hit the cap, so no unchanged set can be reported.
///
/// The minimal `cost` is still known — it comes from the dynamic-programming
/// grid, which is polynomial and always completes — so the cost is reported
/// even though the alignment set is not.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CapExceeded {
    /// The model the refused enumeration was run under.
    pub model: CostModel,
    /// The minimal edit cost, which is computable regardless of the cap.
    pub cost: u32,
    /// The bound that was exceeded.
    pub cap: usize,
}

impl std::fmt::Display for CapExceeded {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "more than {} minimal alignments at cost {} under {:?}; \
             no unchanged set reported, because an intersection over a \
             truncated set is too permissive",
            self.cap, self.cost, self.model
        )
    }
}

/// Every minimal alignment of one `(reference, alternate)` block, and what they
/// agree on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MinimalAlignments {
    model: CostModel,
    cost: u32,
    alignments: Vec<Vec<Step>>,
    unchanged: Vec<u32>,
}

impl MinimalAlignments {
    /// The model this was computed under.
    pub fn model(&self) -> CostModel {
        self.model
    }

    /// The minimal edit cost.
    pub fn cost(&self) -> u32 {
        self.cost
    }

    /// How many distinct optimal alignments there are.
    ///
    /// Reported because "unique" and "several that happen to agree" are
    /// different situations: the first pins one reading of the block, the
    /// second means the unchanged set is a genuine intersection and the
    /// alternate coordinates it pairs with are *not* pinned. That distinction
    /// is what separates this notion from the cell-based one — see the module
    /// docs.
    pub fn count(&self) -> usize {
        self.alignments.len()
    }

    /// Whether exactly one alignment attains the minimal cost.
    pub fn is_unique(&self) -> bool {
        self.alignments.len() == 1
    }

    /// Reference offsets matched in **every** minimal alignment — the record's
    /// "unchanged" set. Ascending, deduplicated.
    pub fn unchanged(&self) -> &[u32] {
        &self.unchanged
    }

    /// Every optimal alignment, as its sequence of columns, in a deterministic
    /// order (diagonal before deletion before insertion, depth first).
    pub fn alignments(&self) -> &[Vec<Step>] {
        &self.alignments
    }

    /// Reference offsets that alignment `index` matches, ascending.
    ///
    /// [`unchanged`](Self::unchanged) is the intersection of this over every
    /// index.
    ///
    /// # Panics
    ///
    /// If `index >= self.count()`.
    pub fn matched_reference_offsets(&self, index: usize) -> Vec<u32> {
        matched_reference_offsets(&self.alignments[index])
    }
}

/// Enumerate the minimal alignments of `reference -> alternate` under `model`,
/// capped at [`DEFAULT_ALIGNMENT_CAP`].
pub fn minimal_alignments(
    reference: &[u8],
    alternate: &[u8],
    model: CostModel,
) -> Result<MinimalAlignments, CapExceeded> {
    minimal_alignments_capped(reference, alternate, model, DEFAULT_ALIGNMENT_CAP)
}

/// As [`minimal_alignments`], with a caller-chosen bound on how many optimal
/// alignments will be enumerated before refusing.
///
/// A `cap` of 0 refuses everything, which is a legitimate way to ask only for
/// the cost.
pub fn minimal_alignments_capped(
    reference: &[u8],
    alternate: &[u8],
    model: CostModel,
    cap: usize,
) -> Result<MinimalAlignments, CapExceeded> {
    let grid = SuffixCosts::build(reference, alternate, model);
    let cost = grid.at(0, 0);
    let alignments = enumerate(reference, alternate, model, &grid, cap).ok_or(CapExceeded {
        model,
        cost,
        cap,
    })?;
    debug_assert!(
        !alignments.is_empty(),
        "every block has at least one alignment"
    );
    let unchanged = intersect_matched_reference_offsets(reference.len(), &alignments);
    Ok(MinimalAlignments {
        model,
        cost,
        alignments,
        unchanged,
    })
}

/// Reference offsets where `reference` and `alternate` agree **position-wise**,
/// i.e. under the correspondence `i <-> i`.
///
/// This is the reading the record rules out, kept here so tests can state the
/// difference rather than assert it by hand. Offsets past the shorter of the
/// two are never matched, so the two sequences need not be the same length.
///
/// Neither this set nor [`MinimalAlignments::unchanged`] contains the other in
/// general — see `minimal_alignment_enumeration::neither_notion_contains_the_other`.
pub fn position_wise_matches(reference: &[u8], alternate: &[u8]) -> Vec<u32> {
    reference
        .iter()
        .zip(alternate.iter())
        .enumerate()
        .filter(|(_, (r, a))| r == a)
        .map(|(i, _)| i as u32)
        .collect()
}

/// Reference offsets an alignment matches, ascending.
fn matched_reference_offsets(alignment: &[Step]) -> Vec<u32> {
    let mut matched = Vec::new();
    let mut ref_offset = 0u32;
    for &step in alignment {
        match step {
            Step::Match => {
                matched.push(ref_offset);
                ref_offset += 1;
            }
            Step::Sub | Step::Del => ref_offset += 1,
            Step::Ins => {}
        }
    }
    matched
}

/// The intersection of [`matched_reference_offsets`] over every alignment.
fn intersect_matched_reference_offsets(ref_len: usize, alignments: &[Vec<Step>]) -> Vec<u32> {
    let mut present = vec![true; ref_len];
    for alignment in alignments {
        let mut here = vec![false; ref_len];
        for offset in matched_reference_offsets(alignment) {
            here[offset as usize] = true;
        }
        for (keep, seen) in present.iter_mut().zip(here) {
            *keep &= seen;
        }
    }
    present
        .into_iter()
        .enumerate()
        .filter(|(_, keep)| *keep)
        .map(|(i, _)| i as u32)
        .collect()
}

/// `costs[i][j]` = the minimal cost of aligning `reference[i..]` to
/// `alternate[j..]`.
///
/// Suffix rather than prefix costs, so the enumerator can walk **forwards**
/// from `(0, 0)` taking any edge that keeps the running total minimal. That
/// makes every depth-first branch productive: each one reaches the sink, so the
/// search never explores a dead end and the cap bounds real work rather than
/// wasted work.
struct SuffixCosts {
    width: usize,
    costs: Vec<u32>,
}

impl SuffixCosts {
    fn build(reference: &[u8], alternate: &[u8], model: CostModel) -> Self {
        let n = reference.len();
        let m = alternate.len();
        let width = m + 1;
        let mut costs = vec![0u32; (n + 1) * width];
        // Aligning a suffix against nothing costs one indel per base.
        for i in 0..=n {
            costs[i * width + m] = (n - i) as u32;
        }
        for j in 0..=m {
            costs[n * width + j] = (m - j) as u32;
        }
        for i in (0..n).rev() {
            for j in (0..m).rev() {
                let diagonal = if reference[i] == alternate[j] {
                    0
                } else {
                    model.substitution_cost()
                };
                costs[i * width + j] = (diagonal + costs[(i + 1) * width + j + 1])
                    .min(1 + costs[(i + 1) * width + j])
                    .min(1 + costs[i * width + j + 1]);
            }
        }
        Self { width, costs }
    }

    fn at(&self, i: usize, j: usize) -> u32 {
        self.costs[i * self.width + j]
    }
}

/// One frame of the explicit depth-first stack.
struct Frame {
    i: usize,
    j: usize,
    /// Which of the three out-edges to try next: 0 diagonal, 1 deletion,
    /// 2 insertion.
    next_edge: u8,
}

/// Depth-first enumeration of every optimal alignment; `None` once more than
/// `cap` of them have been found.
///
/// Iterative rather than recursive on purpose: recursion depth would be
/// `reference.len() + alternate.len()`, which a long forced run of matches
/// reaches even though it yields a single alignment, so the stack depth would
/// be bounded by the *inputs* rather than by the cap.
fn enumerate(
    reference: &[u8],
    alternate: &[u8],
    model: CostModel,
    grid: &SuffixCosts,
    cap: usize,
) -> Option<Vec<Vec<Step>>> {
    let n = reference.len();
    let m = alternate.len();
    let mut alignments: Vec<Vec<Step>> = Vec::new();
    let mut path: Vec<Step> = Vec::new();
    let mut stack = vec![Frame {
        i: 0,
        j: 0,
        next_edge: 0,
    }];

    while let Some(&Frame { i, j, next_edge }) = stack.last() {
        if i == n && j == m {
            alignments.push(path.clone());
            if alignments.len() > cap {
                return None;
            }
            stack.pop();
            if !stack.is_empty() {
                path.pop();
            }
            continue;
        }

        let mut descended = None;
        let mut edge = next_edge;
        while edge < 3 {
            let candidate = match edge {
                0 if i < n && j < m => {
                    let step = if reference[i] == alternate[j] {
                        Step::Match
                    } else {
                        Step::Sub
                    };
                    Some((i + 1, j + 1, step))
                }
                1 if i < n => Some((i + 1, j, Step::Del)),
                2 if j < m => Some((i, j + 1, Step::Ins)),
                _ => None,
            };
            edge += 1;
            if let Some((next_i, next_j, step)) = candidate {
                if model.step_cost(step) + grid.at(next_i, next_j) == grid.at(i, j) {
                    descended = Some((next_i, next_j, step));
                    break;
                }
            }
        }

        // `edge` has advanced past every option tried, so this frame resumes
        // after the one it just took.
        stack
            .last_mut()
            .expect("the frame just read is still on the stack")
            .next_edge = edge;

        match descended {
            Some((next_i, next_j, step)) => {
                path.push(step);
                stack.push(Frame {
                    i: next_i,
                    j: next_j,
                    next_edge: 0,
                });
            }
            None => {
                stack.pop();
                if !stack.is_empty() {
                    path.pop();
                }
            }
        }
    }

    Some(alignments)
}
