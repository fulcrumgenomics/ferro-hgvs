//! Adapter rules — Stage-B `merge.rs` passes migrated into
//! [`Rule`](crate::partition::ruled::Rule) impls by **calling the existing pass
//! unchanged** through the [`crate::partition::bridge`] (Fable design §10 step 1,
//! ruling (a)). Contrast [`crate::partition::rules`], whose rules are hand-rewritten
//! NATIVE over a [`Cut`] because their semantics *is* a coarsening that must READ a
//! lock (`SepZero`).
//!
//! # An adapter runs a label-blind pass for GEOMETRY, then types its own output
//!
//! The `merge.rs` passes operate on `Vec<Piece>`, which carries no label and no lock
//! — so a pass can only ever produce geometry, and the run typing happens afterwards.
//! At step 1 the bridge lifted a pass's output entirely `Unlabelled` ([`bridge::lift`]),
//! which is why the module doc used to say an adapter "can only emit unlocked runs" and
//! that a lock-setting rule was thereby forced native. **Step 5 (design §4.8) retires
//! that claim:** an adapter now lifts through [`bridge::lift_relabel`], which re-derives
//! each output run's label from its geometry using the render stage's own recognisers
//! ([`merge::is_inversion`], [`merge::is_tandem_duplication`]) and carries an earlier
//! rule's label forward on runs the pass left untouched. So [`RunInv`] and
//! [`TandemDupRun`] set `Inv`/`Dup` locks while remaining adapters — the pass supplies
//! the geometry, the typer supplies the label. What still forces a rule native is
//! needing to *read* a lock mid-pass, which only [`crate::partition::rules::SepZero`]
//! does. The [`RULE_MIGRATION`] inventory records which each rule is.
//!
//! # The adapter shape (the pattern the adapters copy)
//!
//! ```ignore
//! fn apply(&self, cut: &Cut, ctx: &BlockCtx) -> Option<Cut> {
//!     let mut pieces = cut_to_pieces(cut);
//!     let before = pieces.clone();
//!     merge::<the pass>(&mut pieces, ctx.reference);   // the SHIPPING pass, executed
//!     if pieces == before {
//!         return None;                                 // no-op: rule did not fire
//!     }
//!     // A Coarsen lifts through `lift_coarsen`; a Recut through `lift_recut`.
//!     lift_coarsen(ctx, cut, pieces)                   // re-validate; decline on unsound or lock
//! }
//! ```
//!
//! The `before == after` guard reproduces the pass's own no-op case as a declining
//! rule. **Most** adapters here are **Coarsen**s that only ever *reduce* the run count
//! when they fire, so a firing application strictly decreases `Φ = (unlocked, total)`;
//! they lift through [`bridge::lift_coarsen`], which additionally REFUSES to fire when
//! the pass would absorb a locked `Dup`/`Inv` — the lock guard for the label-blind
//! passes that `SepZero` (native) already applies (design §4.2). The exceptions are
//! [`TandemDupRun`] and [`RunInv`], **Recut**s whose peel/flanking can *increase* the
//! count; each lifts through [`bridge::lift_recut`], which locks the whole decomposition
//! it produced (not only the geometric `Dup`/`Inv`), so `unlocked` strictly decreases
//! and the Recut is `Φ`-valid in the fixed point.
//! Byte-identity to the shipping pipeline is definitional: the pass body executed is
//! the shipping one, and a label rides on top of unchanged geometry, so
//! [`cut_to_pieces`] of the result is unchanged.

use crate::normalize::merge;
use crate::partition::block_ctx::{BlockCtx, FrameContext};
use crate::partition::bridge::{
    cut_to_pieces, lift, lift_coarsen, lift_recut, lift_relabel, preserve_only,
};
use crate::partition::ruled::{
    Authority, Cut, DirectionScope, FrameScope, Label, MoleculeScope, Rule, RuleKind, Scope,
};

/// How a migrated Stage-B pass is realised: an [`Adapter`](Migration::Adapter) that
/// wraps the named `merge.rs` pass through the bridge, or a hand-rewritten
/// [`Native`](Migration::Native) rule. A rule is native only when it must READ a lock
/// mid-pass ([`crate::partition::rules::SepZero`]); a rule that merely SETS a lock
/// stays an adapter, typing its output via [`bridge::lift_relabel`] (module doc). So
/// the classification is fixed by what the rule does with locks, not a free choice.
///
/// The `Migration`/[`RULE_MIGRATION`] table is the design §8 realisation audit,
/// read only by the tests that cross-check it against `pipeline_order()`; nothing
/// on a runtime path reads it, so both are `#[cfg(test)]` (dev-only already, via
/// the `partition` module) — the guardrail is kept without a `dead_code` allow.
#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Migration {
    /// Hand-rewritten over a [`Cut`] — its semantics must *read* a lock (decline to
    /// coarsen across one), which the label-blind pass cannot see.
    Native,
    /// Wraps the named `merge.rs` pass, executed unchanged via the bridge. May still
    /// SET `Inv`/`Dup` locks on its output by typing the pass's geometry
    /// ([`RunInv`], [`TandemDupRun`]).
    Adapter {
        /// The `merge.rs` function this adapter calls, for the byte-identity audit.
        wraps: &'static str,
    },
}

/// One row per migrated rule: its [`Rule::id`] and how it is realised. Enforced-table
/// guardrail (design §8): a rule may not silently change realisation. Grown as rules
/// land; `SepZero` is native (it READS the lock exclusion), every other pass is an
/// adapter — including the two Recuts that SET `Inv`/`Dup` locks on their typed output
/// (step 5; the "an adapter never locks" invariant was retired with [`bridge::lift`]'s
/// monopoly on lifting — see [`bridge::lift_relabel`]).
#[cfg(test)]
pub(crate) const RULE_MIGRATION: &[(&str, Migration)] = &[
    ("sep-zero", Migration::Native),
    (
        "payload-coincidence",
        Migration::Adapter {
            wraps: "coalesce_payload_alignment_split",
        },
    ),
    (
        "compensating-gaps",
        Migration::Adapter {
            wraps: "coalesce_compensating_gap_split",
        },
    ),
    (
        "placed-gap",
        Migration::Adapter {
            wraps: "split_is_a_placed_gap_coincidence",
        },
    ),
    (
        "codon-frame-separation",
        Migration::Adapter {
            wraps: "coalesce_coding_frame_separation",
        },
    ),
    (
        "codon-exception",
        Migration::Adapter {
            wraps: "apply_coding_codon_exception",
        },
    ),
    (
        "tandem-dup-run",
        Migration::Adapter {
            wraps: "coalesce_by_run{peel_tandem_dup_beside_change;coalesce_solid_run}",
        },
    ),
    (
        "run-inv",
        Migration::Adapter {
            wraps: "coalesce_inversion_runs",
        },
    ),
    (
        "split-concealed",
        Migration::Adapter {
            wraps: "split_concealed_separations",
        },
    ),
];

/// `PayloadCoincidence` — the `delins.md:46-47` payload-coincidence merge. Wraps
/// [`merge::coalesce_payload_alignment_split`]: a lone net-deletion split whose
/// payload re-aligns within the span is one spanning `delins`.
///
/// Coarsen. `DnaOnly` (the carve-out is DNA-axis-scoped,
/// `delins-payload-coincidence-carve-out-is-coding-dna-scoped`) + `NetDeletion` (the
/// pass fires only when `payload.len() < span.len()`, which — because unchanged flanks
/// add equally to reference and resulting length — is exactly whole-block net
/// deletion, so the declared direction scope will equal the computed gate at step 3).
/// The molecule limb is now enforced by the engine's `Scope::admits`; deleting the
/// call-site `payload_coalesce_applies` gate and asserting declared == computed is
/// the rest of step 3.
pub struct PayloadCoincidence;

impl Rule for PayloadCoincidence {
    fn id(&self) -> &'static str {
        "payload-coincidence"
    }

    fn kind(&self) -> RuleKind {
        RuleKind::Coarsen
    }

    fn scope(&self) -> Scope {
        Scope {
            molecule: MoleculeScope::DnaOnly,
            frame: FrameScope::Any,
            direction: DirectionScope::NetDeletion,
        }
    }

    fn authority(&self) -> Authority {
        Authority::Ruling("delins-recommendation-reach-when-the-input-arrives-split")
    }

    fn apply(&self, cut: &Cut, ctx: &BlockCtx) -> Option<Cut> {
        let mut pieces = cut_to_pieces(cut);
        let before = pieces.clone();
        merge::coalesce_payload_alignment_split(&mut pieces, ctx.reference);
        if pieces == before {
            // The pass declined this block; the rule does not fire.
            return None;
        }
        // A coarsening preserves the splice, so `Cut::new` (inside the lift)
        // re-validates and returns `Ok`; `lift_coarsen` additionally declines if the
        // pass would absorb a locked `Dup`/`Inv` (the fixed-point lock guard).
        lift_coarsen(ctx, cut, pieces)
    }
}

/// `CompensatingGaps` — the `delins.md:44-47` compensating-gap merge. Wraps
/// [`merge::coalesce_compensating_gap_split`]: a split an alignment manufactured by
/// inserting in one member and deleting in another (a non-monotone offset) is one
/// spanning `delins`.
///
/// Coarsen, `DnaOnly` — the same axis scope as [`PayloadCoincidence`] (its call-site
/// gate `compensating_gap_coalesce_applies` is `cuts_with_canonical && is_dna`). No
/// direction limb: the compensating shape is net-neutral as often as net-deletion.
/// The precise ledger record is owed at step 4 (this pass cites none of its own
/// today — see the migration checklist); the authority here names the clause it
/// implements.
pub struct CompensatingGaps;

impl Rule for CompensatingGaps {
    fn id(&self) -> &'static str {
        "compensating-gaps"
    }

    fn kind(&self) -> RuleKind {
        RuleKind::Coarsen
    }

    fn scope(&self) -> Scope {
        Scope {
            molecule: MoleculeScope::DnaOnly,
            frame: FrameScope::Any,
            direction: DirectionScope::Any,
        }
    }

    fn authority(&self) -> Authority {
        Authority::Clause("DNA/delins.md:44")
    }

    fn apply(&self, cut: &Cut, ctx: &BlockCtx) -> Option<Cut> {
        let mut pieces = cut_to_pieces(cut);
        let before = pieces.clone();
        merge::coalesce_compensating_gap_split(&mut pieces, ctx.reference);
        if pieces == before {
            return None;
        }
        lift_coarsen(ctx, cut, pieces)
    }
}

/// `PlacedGap` — the #1610 placed-gap veto. Wraps
/// [`merge::split_is_a_placed_gap_coincidence`]: a lone, minimal, unequal-length
/// net-deletion split whose payload merely coincides with the reference is not a
/// separation but one whole-block `delins`. Unlike the other adapters this wraps a
/// *predicate* rather than a mutator — when it holds, the caller collapses the
/// trimmed block to a single member spanning `[lo, hi_ref)`, which this adapter
/// reproduces. The ctx is the padded WINDOW (design §12 R1), and the pipeline
/// evaluates the predicate on BLOCK-frame pieces over the trimmed block (inside
/// `partition_block_for_rule`, before any shift), so the adapter re-derives the
/// block with `trim_common_flanks` and moves the runs into its frame first.
///
/// Coarsen, `DnaOnly`/`NetDeletion` — the predicate's own gates are
/// `carve_out.may_disbelieve_a_separation()` (a DNA axis) and `result.len() <
/// reference.len()` with `pieces.len() >= 2`. The engine's `Scope::admits` now gates
/// the RULE on molecule, but the predicate's own signature still takes a
/// `CoincidenceCarveOut`, so the adapter keeps deriving it from `ctx.molecule`
/// (mirroring `bakeoff_arm_override_pieces`) — that local derivation does NOT go away
/// at step 3; only the call-site `_applies` gates do.
pub struct PlacedGap;

impl Rule for PlacedGap {
    fn id(&self) -> &'static str {
        "placed-gap"
    }

    fn kind(&self) -> RuleKind {
        RuleKind::Coarsen
    }

    fn scope(&self) -> Scope {
        Scope {
            molecule: MoleculeScope::DnaOnly,
            frame: FrameScope::Any,
            direction: DirectionScope::NetDeletion,
        }
    }

    fn authority(&self) -> Authority {
        Authority::Ruling("unequal-length-block-a-placed-gap-is-not-a-separation")
    }

    fn apply(&self, cut: &Cut, ctx: &BlockCtx) -> Option<Cut> {
        use crate::normalize::merge::{CoincidenceCarveOut, Piece};
        use crate::partition::block_ctx::Molecule;

        // The pipeline evaluates this predicate inside `partition_block_for_rule`,
        // on BLOCK-frame pieces over the trimmed block, before any offset or shift.
        // The ctx is the padded WINDOW (design §12 R1), so re-derive the block with
        // the same trim the pipeline used and move the pieces into its frame; a
        // piece outside the block cannot occur here (nothing has shifted yet) and
        // declines defensively.
        let (lo, hi_ref, hi_alt) = merge::trim_common_flanks(ctx.reference, ctx.resulting);
        let block_ref = &ctx.reference[lo..hi_ref];
        let block_alt = &ctx.resulting[lo..hi_alt];
        let mut block_pieces = Vec::with_capacity(cut.runs().len());
        for r in cut.runs() {
            if r.ref_start < lo || r.ref_end > hi_ref {
                return None;
            }
            block_pieces.push(Piece {
                ref_start: r.ref_start - lo,
                ref_end: r.ref_end - lo,
                alt: r.alt.clone(),
            });
        }
        // The predicate's own axis gate needs the carve-out; derive it from the
        // block's molecule, the inverse of `bakeoff_arm_override_pieces`'s mapping.
        let carve_out = match ctx.molecule {
            Molecule::Dna => CoincidenceCarveOut::InReach,
            _ => CoincidenceCarveOut::OutOfReach,
        };
        if !merge::split_is_a_placed_gap_coincidence(&block_pieces, block_ref, block_alt, carve_out)
        {
            return None;
        }
        // The caller's collapse, in window coordinates: one piece spanning the
        // block. The predicate requires `pieces.len() >= 2`, so this reduces the
        // run count (Phi strictly decreases); the guard is belt-and-braces.
        let collapsed = vec![Piece {
            ref_start: lo,
            ref_end: hi_ref,
            alt: block_alt.to_vec(),
        }];
        if collapsed == cut_to_pieces(cut) {
            return None;
        }
        lift(ctx, collapsed)
    }
}

/// The pipeline's CDS-axis coordinate of window offset 0, reconstructed from the
/// frame's zones — the `w_lo` both codon adapters pass to their `merge.rs` passes.
///
/// The passes read 1-indexed CDS positions (`same_codon` is `(pos-1)/3` and rejects
/// `pos < 1`), so this is a coordinate, not the raw phase:
///  * `cds_start = Some(s)`: c.1 sits at window offset `s`, so offset 0 is `c.(1-s)`.
///    This also encodes the 5'UTR seam for free — offsets `< s` map to `<= 0`, which
///    `same_codon` rejects (no codon merge into the 5'UTR), exactly as the pipeline's
///    signed coordinate does.
///  * `cds_start = None`: the CDS begins at/before the window, so every offset is
///    `>= c.1`; the absolute number is unknown, but `same_codon` is invariant to a
///    `w_lo` shift of 3, so `grid_phase()+1` (grid_phase `== (w_lo-1) mod 3`, all
///    positions `>= 1`) reproduces the pipeline's coordinate byte-for-byte. The
///    `cds_end_axis` comparison is likewise `w_lo`-invariant (it reduces to a window
///    offset), so `CodonException` stays exact with the synthetic origin too.
///  * `NonCoding`: no reading frame, so a caller's `reading_frame` guard makes the
///    pass a no-op; the returned 0 is unused.
fn cds_axis_origin(frame: &FrameContext) -> i64 {
    match frame {
        FrameContext::Coding {
            cds_start: Some(s), ..
        } => 1 - *s as i64,
        _ => frame.grid_phase().map_or(0, |p| i64::from(p) + 1),
    }
}

/// The inverse of [`cds_axis_origin`]: the `FrameContext` whose origin reproduces
/// the pipeline's `w_lo` (the CDS-axis coordinate of window offset 0) up to a codon,
/// exactly at the 5' seam, and whose `cds_end` preserves the window offset of the
/// last CDS base. The wiring builds the ruled `BlockCtx` with it.
///
/// * `w_lo <= 1`: c.1 sits at window offset `1 - w_lo`, so `cds_start = Some` pins
///   the grid and the origin is exact (the 5'UTR seam is then encoded).
/// * `w_lo > 1`: the CDS start is off-window; only the phase is known, and
///   `same_codon` is invariant to a shift of 3.
/// * `cds_end_axis = Some(e)`: the last CDS base is window offset `e - w_lo`, so one
///   past it is `e - w_lo + 1`, clamped at 0 — a window wholly past the CDS end has
///   every offset in the 3'UTR, which is how the pipeline's `within_cds` reads it.
pub(crate) fn frame_context_for_axis_origin(
    reading_frame: bool,
    w_lo: i64,
    cds_end_axis: Option<i64>,
) -> FrameContext {
    if !reading_frame {
        return FrameContext::NonCoding;
    }
    let cds_end = cds_end_axis
        .map(|e| usize::try_from((e - w_lo + 1).max(0)).expect("clamped to non-negative"));
    if w_lo <= 1 {
        let s = usize::try_from(1 - w_lo).expect("w_lo <= 1");
        FrameContext::coding(((3 - s % 3) % 3) as u8, Some(s), cds_end)
    } else {
        FrameContext::coding(((w_lo - 1) % 3) as u8, None, cds_end)
    }
}

/// `CodonFrameSeparation` — the early half of the `delins.md:18` codon exception.
/// Wraps [`merge::coalesce_coding_frame_separation`]: on a coding, length-changing
/// block, two members one base apart that fall in one codon are merged (the unchanged
/// base between them becomes interior). This pass runs EARLY in the pipeline (before
/// the payload/compensating passes), so it is a rule of its own — NOT bundled with the
/// late `apply_coding_codon_exception`, which is a separate rule (`CodonException`).
///
/// Coarsen, `NeedsCodingFrame`. The pass's `w_lo` is the 1-indexed CDS-axis coordinate
/// of window offset 0 (`same_codon` reads `(pos-1)/3` and rejects `pos<1`), which the
/// adapter reconstructs from the frame's CDS zones — see `apply` for the derivation and
/// why it reproduces both the codon phase and the 5'UTR seam byte-for-byte.
/// `reading_frame` is whether the frame is coding and `length_changing` is
/// `reference.len() != resulting.len()`. The pass self-gates on `reading_frame`, so it
/// is a no-op off a coding frame even before the engine's frame limb is wired (step 3).
pub struct CodonFrameSeparation;

impl Rule for CodonFrameSeparation {
    fn id(&self) -> &'static str {
        "codon-frame-separation"
    }

    fn kind(&self) -> RuleKind {
        RuleKind::Coarsen
    }

    fn scope(&self) -> Scope {
        Scope {
            molecule: MoleculeScope::Any,
            frame: FrameScope::NeedsCodingFrame,
            direction: DirectionScope::Any,
        }
    }

    fn authority(&self) -> Authority {
        Authority::Ruling("delins-codon-carve-out-gap-one")
    }

    fn apply(&self, cut: &Cut, ctx: &BlockCtx) -> Option<Cut> {
        let mut pieces = cut_to_pieces(cut);
        let before = pieces.clone();
        let reading_frame = ctx.frame.grid_phase().is_some();
        let length_changing = ctx.reference.len() != ctx.resulting.len();
        merge::coalesce_coding_frame_separation(
            &mut pieces,
            reading_frame,
            length_changing,
            cds_axis_origin(ctx.frame),
            ctx.reference,
        );
        if pieces == before {
            return None;
        }
        lift_relabel(ctx, cut, pieces, preserve_only)
    }
}

/// `CodonException` — the late half of the `delins.md:18` codon exception. Wraps
/// [`merge::apply_coding_codon_exception`]: a `[Sub; unchanged; Sub]` triplet whose
/// two substitutions fall in one codon (and within the CDS) is merged into one
/// three-base delins. Runs LATE in the pipeline (after the payload/compensating/
/// placed-gap passes), which is why it is a rule of its own rather than bundled with
/// the early `CodonFrameSeparation`.
///
/// Coarsen, `NeedsCodingFrame`. Shares `CodonFrameSeparation`'s `w_lo` reconstruction
/// ([`cds_axis_origin`]) and additionally passes `cds_end_axis` — the 3'UTR seam — as
/// `w_lo + cds_end - 1` (the axis coordinate of the last CDS base). The pass's
/// `within_cds` test only ever compares `w_lo + offset <= cds_end_axis`, which reduces
/// to `offset < cds_end`, so this is exact regardless of `w_lo`'s absolute value. Off a
/// coding frame the pass self-gates to a no-op.
pub struct CodonException;

impl Rule for CodonException {
    fn id(&self) -> &'static str {
        "codon-exception"
    }

    fn kind(&self) -> RuleKind {
        RuleKind::Coarsen
    }

    fn scope(&self) -> Scope {
        Scope {
            molecule: MoleculeScope::Any,
            frame: FrameScope::NeedsCodingFrame,
            direction: DirectionScope::Any,
        }
    }

    fn authority(&self) -> Authority {
        Authority::Ruling("delins-codon-carve-out-gap-one")
    }

    fn apply(&self, cut: &Cut, ctx: &BlockCtx) -> Option<Cut> {
        let mut pieces = cut_to_pieces(cut);
        let before = pieces.clone();
        let reading_frame = ctx.frame.grid_phase().is_some();
        let w_lo = cds_axis_origin(ctx.frame);
        // The CDS end on the same axis as `w_lo + offset`: the coordinate of the last
        // CDS base. `within_cds` only compares it against `w_lo + offset`, so the check
        // is `w_lo`-invariant (it reduces to `offset < cds_end`).
        let cds_end_axis = match ctx.frame {
            FrameContext::Coding {
                cds_end: Some(e), ..
            } => Some(w_lo + *e as i64 - 1),
            _ => None,
        };
        merge::apply_coding_codon_exception(
            &mut pieces,
            reading_frame,
            w_lo,
            ctx.reference,
            cds_end_axis,
        );
        if pieces == before {
            return None;
        }
        lift_coarsen(ctx, cut, pieces)
    }
}

/// `TandemDupRun` — the #2175/#2174 per-run tandem-dup peel plus solid-run collapse.
/// Wraps the pipeline's combined `coalesce_by_run` closure: within each FrameWall-
/// bounded run, peel a tandem dup abutting a change (so the dup survives), then
/// collapse the residual solid run into one delins.
///
/// **This is ONE rule, not two**, even though it runs two `merge.rs` passes: the
/// pipeline applies `{ peel; coalesce_solid_run }` as a single per-run closure inside
/// ONE `coalesce_by_run`, so splitting it into two rules (two separate groupings)
/// would regroup differently — `coalesce_by_run`'s wall placement and overlap
/// reassembly depend on the whole closure — and break byte-identity.
///
/// Recut, `DnaOnly` (call-site gate `cuts_with_canonical() && is_dna`). The peel
/// EXPOSES a dup, so this can INCREASE the run count — it is NOT a Coarsen. As of
/// step 5 it LOCKS the peeled dup, typing it `Dup` via the render stage's
/// [`merge::is_tandem_duplication`] recogniser ([`bridge::lift_relabel`]); the lock
/// is what makes `unlocked` strictly decrease so the Recut is `Φ`-valid under the
/// fixed point, and what stops `SepZero` re-absorbing the dup. It stays an adapter —
/// the pass supplies the geometry, the typer the label. Its label authority is
/// `duplication-must-ranks-the-label-not-the-partition` (see [`Self::authority`]).
pub struct TandemDupRun;

impl Rule for TandemDupRun {
    fn id(&self) -> &'static str {
        "tandem-dup-run"
    }

    fn kind(&self) -> RuleKind {
        RuleKind::Recut
    }

    fn scope(&self) -> Scope {
        Scope {
            molecule: MoleculeScope::DnaOnly,
            frame: FrameScope::Any,
            direction: DirectionScope::Any,
        }
    }

    fn authority(&self) -> Authority {
        // The peel EXPOSES a tandem dup by re-deriving from the resulting sequence
        // (`canonical-form-choice-when-both-legal`, rule 6) and KEEPS it. The label
        // is governed by `duplication-must-ranks-the-label-not-the-partition`
        // (2026-08-13), which reads `duplication.md:17`/`:18` as ranking the LABEL
        // for a change re-derived from the resulting sequence, not as requiring a
        // partition that exposes one — so a `dup`-shaped piece must be a dup. Where
        // that keep-the-dup choice deviates from a spec worked example at separation
        // zero (W43), `separation-zero-dup-member-is-preserved-not-merged` records
        // the deviation. This cites the general label ruling, not the bare `:18`
        // clause the doc above marked as owed.
        Authority::Ruling("duplication-must-ranks-the-label-not-the-partition")
    }

    fn apply(&self, cut: &Cut, ctx: &BlockCtx) -> Option<Cut> {
        let mut pieces = cut_to_pieces(cut);
        let before = pieces.clone();
        merge::coalesce_by_run(&mut pieces, ctx.reference, |run, reference| {
            merge::peel_tandem_dup_beside_change(run, reference, merge::PeelReach::PlusK);
            merge::coalesce_solid_run(run, reference);
        });
        if pieces == before {
            return None;
        }
        // The peel's per-run source-corruption gate (`peel_tandem_dup_beside_change`
        // "Gate A") sees only the run it cut, so under `PlusK` it can mint a dup
        // whose source copy a SIBLING piece overwrote. The render seam
        // (`duplication_anchor`, `source_start < previous_ref_end`) is then
        // guaranteed to refuse it, and `canonicalize_from_sequence_with_rule`
        // silently falls back to the per-member pipeline — which types each input
        // spelling's own members, so two spellings of one variant leak apart
        // (the `s05-c-m4-sep1-p1-all-ins` confluence loss). Decline the recut
        // instead: leaving the pre-peel cut converges every spelling on the same
        // canonical form. The predicate is the renderer's own, so this fires iff
        // the render would refuse — it never drops a peel that would render (a
        // disjoint-source `[dup;sub]` like #2175 is untouched).
        if merge::dup_source_overlaps_prior_piece(&pieces, ctx.reference) {
            return None;
        }
        // Codon-frame precedence — restore the shipped coding-axis merge (#1744).
        // On a reading-frame axis, `CodonFrameSeparation`
        // (`delins-codon-carve-out-gap-one`) has already coarsened a within-one-
        // codon net-frameshift block to ONE delins before this rule runs (or, on
        // a block this tight, the DAG seed itself derives the identical merged
        // form directly — the pass then has nothing left to do, but the shape is
        // still the codon-frame merge's output SIGNATURE and must be protected
        // the same way) — the form `origin/main` (`CanonicalCoalesced`,
        // `PeelReach::TractOnly`) ships and `coding_frame_merge_axis_asymmetry`
        // pins, and the deliberately-kept status quo of GitHub #1744 (a rule-2
        // preference departure per #1725). The step-8 `PlusK` peel would
        // re-split it into `[sub; dup]`; decline so the merge stands, consistent
        // with `codon-exception-vs-coincidence-carve-out-precedence`'s
        // codon-first ruling.
        //
        // Shape gate: EXACTLY the codon-frame merge's output signature — a lone,
        // single-reference-base (`ref_end - ref_start == 1`), length-changing
        // delins on a reading-frame axis whose alt carries the ONE retained
        // reference base as an interior byte. `coalesce_coding_frame_separation`
        // merges ANY pair one unchanged base apart within one codon, not just two
        // insertions (`codon-carve-out-shape-restriction`: WIDEN, edit-type-
        // independent), so the retained base's ALT POSITION varies by which side
        // of the pair is the pure insertion:
        //
        //   H1 `[ins; unchanged; ins]` — both sides are zero-width, so the
        //   merge's OWN one-base span *is* the retained base: it can sit
        //   anywhere in the alt's interior (both insertions non-empty).
        //   H2 `[sub; unchanged; ins]` — the substitution is exactly one alt
        //   byte (HGVS substitutions are always single-base), so the retained
        //   base is `ref[ref_end]` (one past the piece's own span) and sits at
        //   `alt[1]`, position-exact.
        //   H3 `[ins; unchanged; sub]` (the mirror) — the substitution is again
        //   exactly one alt byte, this time the LAST one, so the retained base
        //   is `ref[ref_start - 1]` (one before the span) and sits at
        //   `alt[alt.len() - 2]`, position-exact.
        //
        // H2/H3 additionally require `same_codon` over the two positions the
        // retained base bridges — not a defensive add-on but the pass's own
        // `one_amino_acid` precondition (`delins.md:18`'s "together affecting
        // one amino acid"), re-checked here because a single collapsed piece no
        // longer carries the two-piece geometry that precondition was tested
        // against. Without it, W43 (`c.[13G>A;13_14insTA]` -> delins `ATA`) is a
        // BYTE-LEVEL look-alike of H3: it is a SUB immediately followed by an
        // insertion at separation ZERO (not one), so `ref[ref_start - 1]`
        // (`c.12` in that fixture) is an arbitrary, unrelated base that happens
        // to equal `alt[alt.len() - 2]` in that fixture's sequence — but `c.12`
        // and `c.13` straddle a codon boundary (codon 4 / codon 5), so
        // `same_codon` is false and the shape gate correctly declines to match
        // it. Genomic/frameless dup-beside-change exposure (the R7 disclosure)
        // and multi-codon blocks are untouched — `ctx.frame.grid_phase()` gates
        // on a reading frame at all, and `same_codon` gates the width within it.
        //
        // Independent-dup veto (measured, not guessed): a shape-gate match is
        // NECESSARY but not SUFFICIENT. `ATGCCTGAAACCACGTACGTACGT:c.[10A>C;
        // 11_12insCC]` matches H2 byte-for-byte (`alt[1]=='C'==ref[ref_end]`,
        // same codon) yet `origin/main` KEEPS the split
        // (`c.[10A>C;11_12dup]`): the payload `CC` is *itself*, independent of
        // the substitution, a verbatim tandem duplication of the two reference
        // bases it abuts (`duplication.md:18`'s MUST), so `PeelReach::TractOnly`
        // — the shipped canonical reach, no PlusK needed — finds and keeps that
        // same dup on `before` alone. Re-running the narrower TractOnly peel
        // here and comparing to `before` is the canonical arm's own oracle for
        // "is this dup mandatory", not a heuristic: when TractOnly agrees there
        // is nothing to peel (as for this rule's actual repro and for W43 — the
        // dup is real in both, but is reachable only via `PlusK`'s wider reach,
        // or not at all in W43's separation-zero case), the veto does not fire
        // and the shape-gate match stands.
        if ctx.frame.grid_phase().is_some() && before.len() == 1 {
            let d = &before[0];
            let w_lo = cds_axis_origin(ctx.frame);
            let span_one = d.ref_end - d.ref_start == 1;
            let long_enough = d.alt.len() >= 3;
            let h1_ins_unchanged_ins = span_one
                && long_enough
                && d.alt[1..d.alt.len() - 1].contains(&ctx.reference[d.ref_start]);
            let h2_sub_unchanged_ins = span_one
                && long_enough
                && ctx.reference.get(d.ref_end).is_some_and(|&b| d.alt[1] == b)
                && merge::same_codon(w_lo + d.ref_start as i64, w_lo + d.ref_end as i64);
            let h3_ins_unchanged_sub = span_one
                && long_enough
                && d.ref_start >= 1
                && ctx
                    .reference
                    .get(d.ref_start - 1)
                    .is_some_and(|&b| d.alt[d.alt.len() - 2] == b)
                && merge::same_codon(w_lo + d.ref_start as i64 - 1, w_lo + d.ref_start as i64);
            if h1_ins_unchanged_ins || h2_sub_unchanged_ins || h3_ins_unchanged_sub {
                let mut tract_only = before.clone();
                merge::coalesce_by_run(&mut tract_only, ctx.reference, |run, reference| {
                    merge::peel_tandem_dup_beside_change(
                        run,
                        reference,
                        merge::PeelReach::TractOnly,
                    );
                    merge::coalesce_solid_run(run, reference);
                });
                if tract_only == before {
                    return None;
                }
            }
        }
        // Type a peeled dup with the render stage's own recogniser and LOCK it
        // (design §4.2/§4.8). A Recut lift (`lift_recut`) locks EVERY run this pass
        // produced — the dup and the residual it was peeled from — so the whole
        // committed decomposition survives the fixed point and `Φ`'s unlocked count
        // strictly decreases. Runs the pass left untouched keep their incoming labels.
        lift_recut(ctx, cut, pieces, |p, r| {
            merge::is_tandem_duplication(p, r).then_some(Label::Dup)
        })
    }
}

/// `RunInv` — inversion typing over runs of consecutive pieces. Wraps
/// [`merge::coalesce_inversion_runs`]: it recognises a whole-block, whole-span, or
/// interior/flanked reverse complement and rewrites the affected run into the inv
/// geometry (whole-span → one member; flanked → `[del; inv; del]`).
///
/// **This is a geometry adapter, NOT a native rule — resolved 2026-09-06 (the
/// migration's OPEN QUESTION).** `coalesce_inversion_runs` is a pure `&mut Vec<Piece>`
/// mutator that materialises the reverse complement into a member's `alt` and leaves
/// TYPING to the render stage ([`merge::is_inversion`]/`anchor_for_piece`). Step 5
/// applies that same recogniser to the pass's output via [`bridge::lift_relabel`], so
/// this adapter now SETS the `Inv` label + lock on the whole-span inversion — it did
/// NOT have to go native (the Fable framing that the inversions "cannot be adapters"
/// held only while [`bridge::lift`] was the only lift; the hybrid retires it).
///
/// Recut, `Any/Any/Any` — the pipeline calls it UNGATED (`merge.rs`'s pass call has
/// no `is_dna`/frame guard), and it can INCREASE the run count (the flanked route
/// invents `[del; inv; del]` from one piece). It lifts through
/// [`bridge::lift_recut`], which locks the WHOLE decomposition, not only the geometric
/// `inv`: the flanked route's flanking `del`s are `Unlabelled` but locked, so
/// `unlocked` falls `1 → 0` and the Recut strictly decreases `Φ` even though `total`
/// rises (design §4.5 — a re-cut turns one unlocked run into `n` locked runs). Whole-
/// span inversions (`1 → 1` locked) are `Φ`-safe the same way. The ctx is the padded
/// canonical WINDOW (design §12 R1), so the
/// adapter re-derives the trimmed block with the pipeline's own `trim_common_flanks`
/// and passes `(pieces, ctx.reference, lo, &ctx.reference[lo..hi_ref],
/// &ctx.resulting[lo..hi_alt])` — exactly what the pipeline hands
/// `coalesce_inversion_runs`. On a fully-trimmed block (`lo == 0`) this is the
/// whole-window call `merge.rs`'s own unit tests exercise.
pub struct RunInv;

impl Rule for RunInv {
    fn id(&self) -> &'static str {
        "run-inv"
    }

    fn kind(&self) -> RuleKind {
        RuleKind::Recut
    }

    fn scope(&self) -> Scope {
        Scope {
            molecule: MoleculeScope::Any,
            frame: FrameScope::Any,
            direction: DirectionScope::Any,
        }
    }

    fn authority(&self) -> Authority {
        Authority::Ruling("whole-span-reverse-complement-types-as-inv")
    }

    fn apply(&self, cut: &Cut, ctx: &BlockCtx) -> Option<Cut> {
        let mut pieces = cut_to_pieces(cut);
        let before = pieces.clone();
        // The ctx is the padded canonical WINDOW (design §12 R1), so the block the
        // pass types is re-derived with the pipeline's own trim. The pipeline calls
        // it `(pieces, &ref_bytes, lo, &ref_bytes[lo..hi_ref], &result[lo..hi_alt])`;
        // Route 0 reads the block (is it a reverse complement?), so handing it the
        // whole window would never fire on a padded one.
        let (lo, hi_ref, hi_alt) = merge::trim_common_flanks(ctx.reference, ctx.resulting);
        merge::coalesce_inversion_runs(
            &mut pieces,
            ctx.reference,
            lo,
            &ctx.reference[lo..hi_ref],
            &ctx.resulting[lo..hi_alt],
        );
        if pieces == before {
            return None;
        }
        // Type each inversion run with the render stage's own recogniser and LOCK
        // it (design §4.2/§4.8). Runs the pass left untouched keep their incoming
        // labels.
        let result = lift_recut(ctx, cut, pieces, |p, r| {
            merge::is_inversion(p, r).then_some(Label::Inv)
        })?;
        // F4: retire the FLANKED route. `coalesce_inversion_runs`'s flanked route
        // ([del; inv; del]) produces an `inv` at separation zero from a
        // reference-consuming run, which `delins-adjacent-members-when-both-consume-
        // reference` merges into one delins — and `SepZero` now does exactly that
        // (F4). Emitting the split here would CYCLE with `SepZero` (the fixed point
        // oscillates and `Φ` rises, panicking `run_engine`'s progress assert), and
        // the merged delins is the decided form anyway. So decline when the result
        // carries a sep-zero `Inv` beside a reference-consuming run: a whole-span
        // inversion is seeded as one locked `Inv` (F1) and never reaches here, and a
        // genuinely separated (exact-window) inversion has an unchanged base beside
        // it, so both are left untouched.
        let flanked = result.runs().windows(2).any(|w| {
            w[0].ref_end == w[1].ref_start
                && ((matches!(w[0].label, Label::Inv) && w[1].ref_end > w[1].ref_start)
                    || (matches!(w[1].label, Label::Inv) && w[0].ref_end > w[0].ref_start))
        });
        if flanked {
            return None;
        }
        Some(result)
    }
}

/// `SplitConcealed` — the #1539 concealed-separation audit. Wraps
/// [`merge::split_concealed_separations`]: a member a split emitted that itself hides
/// a `general.md:34` separation (a base unchanged in every minimal alignment,
/// flanked by changes that fall in different codons) is cut apart into its own
/// members.
///
/// **A geometry adapter, and — unlike [`RunInv`] and [`TandemDupRun`] — it stays one
/// even at step 2.** It SPLITS a piece into more pieces, all of them plain
/// `Unlabelled` geometry; it sets no `Inv`/`Dup` label and takes no lock, so the
/// bridge's structural guardrail is not even in tension. It is the cleanest member of
/// the "lock-touching" batch the handoff named: it does not touch a lock at all.
///
/// Recut, `Any/Any/Any` — it can INCREASE the piece count (that is its whole job), so
/// it is a Recut and asserts no `Φ` decrease; the pipeline call is UNGATED across
/// axes. Args mirror [`CodonFrameSeparation`]'s translation: `reading_frame =
/// ctx.frame.grid_phase().is_some()` (`merge.rs`'s pipeline passes
/// `frame.carries_translated_frame()`, the same predicate the codon adapters already
/// map this way), `length_changing = ctx.reference.len() != ctx.resulting.len()`
/// (`hi_ref != hi_alt`), and `w_lo = cds_axis_origin(ctx.frame)` — the same CDS-axis
/// origin the codon adapters reconstruct. The pass takes no `cds_end`; it self-gates
/// to a no-op on a length-neutral block (`!length_changing`), so a substitution-only
/// equal-length block declines even before the engine's direction limb (step 3).
pub struct SplitConcealed;

impl Rule for SplitConcealed {
    fn id(&self) -> &'static str {
        "split-concealed"
    }

    fn kind(&self) -> RuleKind {
        RuleKind::Recut
    }

    fn scope(&self) -> Scope {
        Scope {
            molecule: MoleculeScope::Any,
            frame: FrameScope::Any,
            direction: DirectionScope::Any,
        }
    }

    fn authority(&self) -> Authority {
        Authority::Clause("general.md:34")
    }

    fn apply(&self, cut: &Cut, ctx: &BlockCtx) -> Option<Cut> {
        let mut pieces = cut_to_pieces(cut);
        let before = pieces.clone();
        let reading_frame = ctx.frame.grid_phase().is_some();
        let length_changing = ctx.reference.len() != ctx.resulting.len();
        merge::split_concealed_separations(
            &mut pieces,
            reading_frame,
            length_changing,
            cds_axis_origin(ctx.frame),
            ctx.reference,
        );
        if pieces == before {
            return None;
        }
        lift_relabel(ctx, cut, pieces, preserve_only)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normalize::merge::Piece;
    use crate::partition::block_ctx::{FrameContext, Molecule, Provenance};
    use crate::partition::ruled::{Label, Run};

    /// Whether `run` is geometrically a tandem duplication / inversion, by the render
    /// stage's own recogniser — so a test can assert an adapter typed the run the
    /// renderer would.
    fn as_piece(run: &Run) -> Piece {
        Piece {
            ref_start: run.ref_start,
            ref_end: run.ref_end,
            alt: run.alt.clone(),
        }
    }

    /// Splice `runs` back over `reference` to the sequence they denote — so a test
    /// can build a guaranteed-sound seed from geometry alone.
    fn denote(reference: &[u8], runs: &[(usize, usize, &[u8])]) -> Vec<u8> {
        let mut out = Vec::new();
        let mut cursor = 0usize;
        for &(s, e, alt) in runs {
            out.extend_from_slice(&reference[cursor..s]);
            out.extend_from_slice(alt);
            cursor = e;
        }
        out.extend_from_slice(&reference[cursor..]);
        out
    }

    /// Build a sound seed from `(reference, runs)`, run `rule` and the shipping `raw`
    /// pass on the same pieces, and assert they agree exactly: the adapter returns
    /// `Some` iff the pass changed the pieces, to the same geometry, strictly
    /// decreasing `Φ` and never locking a run. Returns whether the rule fired, so the
    /// caller can guard non-vacuity. This is the byte-identity check at the unit
    /// level — the adapter is compared against the very pass it wraps.
    fn adapter_agrees_with_raw<R: Rule>(
        reference: &[u8],
        runs: &[(usize, usize, &[u8])],
        rule: &R,
        raw: impl Fn(&mut Vec<Piece>, &[u8]),
    ) -> bool {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        let resulting = denote(reference, runs);
        let c = ctx(reference, &resulting, &frame, &prov);
        let seed = Cut::new(
            &c,
            runs.iter()
                .map(|&(s, e, a)| Run::unlabelled(s, e, a.to_vec()))
                .collect(),
        )
        .expect("sound seed");

        let mut raw_pieces = cut_to_pieces(&seed);
        raw(&mut raw_pieces, reference);
        let raw_changed = raw_pieces != cut_to_pieces(&seed);

        match rule.apply(&seed, &c) {
            Some(out) => {
                assert!(raw_changed, "adapter fired but the raw pass was a no-op");
                assert_eq!(
                    cut_to_pieces(&out),
                    raw_pieces,
                    "adapter geometry != raw pass"
                );
                assert!(
                    out.phi() < seed.phi(),
                    "a fired Coarsen must strictly decrease Phi"
                );
                assert!(
                    out.runs().iter().all(|r| !r.locked),
                    "an adapter never returns a locked run",
                );
                true
            }
            None => {
                assert!(
                    !raw_changed,
                    "adapter declined but the raw pass changed the pieces"
                );
                false
            }
        }
    }

    #[test]
    fn compensating_gaps_agrees_with_the_raw_pass_and_fires() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // Ground-truth firing input: merge.rs's own `a_manufactured_split_is_merged_
        // and_a_genuine_one_is_not` row "#1040". The canonical partition of this block
        // is 4 members whose cumulative offset is non-monotone (-1,-1,-1,0 -- deleted
        // then inserted back), the manufactured compensating split the pass merges
        // into one whole-block delins. Built from `partition_block_canonical` rather
        // than hand-typed because the compensating shape needs >= 3 members
        // (`MIN_MANUFACTURED_SPLIT_MEMBERS`) that only the real partitioner produces.
        let reference = b"CAGTGACTAG";
        let alt = b"TGTCACGACT";
        let pieces = merge::partition_block_canonical(reference, alt).expect("canonical partition");
        assert_eq!(
            pieces.len(),
            4,
            "precondition: the canonical partition is 4 members"
        );
        let c = ctx(reference, alt, &frame, &prov);
        let seed = Cut::new(
            &c,
            pieces
                .iter()
                .map(|p| Run::unlabelled(p.ref_start, p.ref_end, p.alt.clone()))
                .collect(),
        )
        .expect("sound seed");

        // The raw shipping pass on the same pieces, for the byte-identity comparison.
        let mut raw_pieces = cut_to_pieces(&seed);
        merge::coalesce_compensating_gap_split(&mut raw_pieces, reference);
        assert_ne!(
            raw_pieces,
            cut_to_pieces(&seed),
            "precondition: the raw pass fires"
        );

        let out = CompensatingGaps
            .apply(&seed, &c)
            .expect("the compensating split merges");
        assert_eq!(
            cut_to_pieces(&out),
            raw_pieces,
            "adapter geometry != raw pass"
        );
        assert!(
            out.phi() < seed.phi(),
            "a fired Coarsen strictly decreases Phi"
        );
        assert!(
            out.runs().iter().all(|r| !r.locked),
            "an adapter never returns a locked run"
        );

        // A plain two-substitution split carries no compensating gap, so it declines.
        let declined = !adapter_agrees_with_raw(
            b"AC",
            &[(0, 1, b"T"), (1, 2, b"G")],
            &CompensatingGaps,
            merge::coalesce_compensating_gap_split,
        );
        assert!(declined, "the substitution split must not fire");
    }

    fn ctx<'a>(
        reference: &'a [u8],
        resulting: &'a [u8],
        frame: &'a FrameContext,
        prov: &'a Provenance,
    ) -> BlockCtx<'a> {
        BlockCtx {
            reference,
            resulting,
            frame,
            molecule: Molecule::Dna,
            provenance: prov,
        }
    }

    #[test]
    fn payload_coincidence_merges_a_net_deletion_coincident_split() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // reference ACGT -> CG: delins AC>C at [0,2) — a *gap-bearing* member, since
        // its alt is one base but it consumes two (`split_carries_a_gap_bearing_insert`
        // requires `alt.len() != consumed`, which a substitution never is) — then del T
        // at [3,4). payload "CG" is an ordered subsequence of "ACGT" and shorter than
        // the [0,4) hull (net deletion), so the pass merges to one [0,4)delinsCG.
        let c = ctx(b"ACGT", b"CG", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 2, b"C".to_vec()),
                Run::unlabelled(3, 4, Vec::new()),
            ],
        )
        .expect("sound seed");

        let out = PayloadCoincidence.apply(&seed, &c).expect("merged");
        assert_eq!(
            out.runs().len(),
            1,
            "the split collapses to one spanning delins"
        );
        assert_eq!((out.runs()[0].ref_start, out.runs()[0].ref_end), (0, 4));
        assert_eq!(out.runs()[0].alt, b"CG");
        // Coarsen obligation: Phi strictly decreased (total 2 -> 1).
        assert!(out.phi() < seed.phi());
        // The adapter guardrail, demonstrated: an adapter can never emit a locked run
        // (the bridge has nowhere to carry a lock).
        assert!(
            out.runs().iter().all(|r| !r.locked),
            "an adapter never returns a locked run"
        );
    }

    #[test]
    fn payload_coincidence_declines_a_length_neutral_split() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // reference AC -> TG: two adjacent substitutions. Neither is gap-bearing
        // (a substitution's alt.len() equals what it consumes) and the block is
        // length-neutral, so the pass declines and the rule does not apply.
        let c = ctx(b"AC", b"TG", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(1, 2, b"G".to_vec()),
            ],
        )
        .expect("sound seed");
        assert!(PayloadCoincidence.apply(&seed, &c).is_none());
    }

    #[test]
    fn placed_gap_collapses_a_coincident_net_deletion_split() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // The #1610 shape (merge.rs's `the_placed_gap_predicate_discriminates`): TGCA
        // -> AAC, cut at the payload's C into [{0,2}AA; {3,4}del]. Net deletion,
        // single-base separation, DNA axis -> the split is a placed-gap coincidence
        // and collapses to one whole-block [0,4)delinsAAC.
        let c = ctx(b"TGCA", b"AAC", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 2, b"AA".to_vec()),
                Run::unlabelled(3, 4, Vec::new()),
            ],
        )
        .expect("sound seed");
        let out = PlacedGap.apply(&seed, &c).expect("collapsed");
        assert_eq!(out.runs().len(), 1);
        assert_eq!((out.runs()[0].ref_start, out.runs()[0].ref_end), (0, 4));
        assert_eq!(out.runs()[0].alt, b"AAC");
        assert!(out.phi() < seed.phi());
        assert!(
            out.runs().iter().all(|r| !r.locked),
            "an adapter never returns a locked run"
        );
    }

    #[test]
    fn placed_gap_declines_a_two_base_separation() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // TGCAT -> AACA: [{0,2}AA; {4,5}del], two unchanged bases (CA) between the
        // members. `general.md:34` keeps them individual, so the predicate declines.
        let c = ctx(b"TGCAT", b"AACA", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 2, b"AA".to_vec()),
                Run::unlabelled(4, 5, Vec::new()),
            ],
        )
        .expect("sound seed");
        assert!(PlacedGap.apply(&seed, &c).is_none());
    }

    #[test]
    fn placed_gap_declines_on_the_rna_axis() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // The SAME collapsing block, but on an RNA molecule: `delins.md:47` does not
        // reach `r.`, so the adapter derives an OutOfReach carve-out from ctx.molecule
        // and the split stands. This calls apply directly (not through the engine), so
        // it exercises the adapter's own carve-out derivation; the engine's now-wired
        // `Scope::admits` DnaOnly limb is the separate, rule-level gate.
        let c = BlockCtx {
            reference: b"TGCA",
            resulting: b"AAC",
            frame: &frame,
            molecule: Molecule::Rna,
            provenance: &prov,
        };
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 2, b"AA".to_vec()),
                Run::unlabelled(3, 4, Vec::new()),
            ],
        )
        .expect("sound seed");
        assert!(PlacedGap.apply(&seed, &c).is_none());
    }

    #[test]
    fn codon_frame_separation_merges_two_members_in_one_codon() {
        // Coding frame, grid phase 0 (codons [0,3),[3,6),...). reference ACGTAC ->
        // CTTAC: del A at [0,1) and sub G>T at [2,3), one unchanged base (C at 1)
        // between them, both inside codon 0. Net length change (-1), so the coding
        // exception merges them into one [0,3)delinsCT, the unchanged C carried
        // interior.
        let frame = FrameContext::coding(0, Some(0), None);
        let prov = Provenance::none();
        let c = ctx(b"ACGTAC", b"CTTAC", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, Vec::new()),
                Run::unlabelled(2, 3, b"T".to_vec()),
            ],
        )
        .expect("sound seed");
        let out = CodonFrameSeparation.apply(&seed, &c).expect("merged");
        assert_eq!(out.runs().len(), 1);
        assert_eq!((out.runs()[0].ref_start, out.runs()[0].ref_end), (0, 3));
        assert_eq!(out.runs()[0].alt, b"CT");
        assert!(out.phi() < seed.phi());
        assert!(
            out.runs().iter().all(|r| !r.locked),
            "an adapter never returns a locked run"
        );
    }

    #[test]
    fn codon_frame_separation_declines_off_a_coding_frame() {
        // The SAME block on a non-coding frame: no reading frame, so the codon
        // exception does not apply and the two members stay split.
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        let c = ctx(b"ACGTAC", b"CTTAC", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, Vec::new()),
                Run::unlabelled(2, 3, b"T".to_vec()),
            ],
        )
        .expect("sound seed");
        assert!(CodonFrameSeparation.apply(&seed, &c).is_none());
    }

    #[test]
    fn codon_frame_separation_respects_the_five_prime_utr_seam() {
        // The SAME two-member block, but with c.1 at offset 3 (cds_start=Some(3)), so
        // offsets 0 and 2 are in the 5'UTR (c.-3, c.-1). The codon exception is a CDS
        // rule, so the merge must NOT happen. This distinguishes the zone-aware w_lo
        // from a naive phase-only one: grid_phase()+1 == 1 would wrongly merge, while
        // the reconstructed w_lo (1 - 3 = -2) makes same_codon reject the sub-1
        // positions, exactly as the pipeline's signed coordinate does.
        let frame = FrameContext::coding(0, Some(3), None);
        let prov = Provenance::none();
        let c = ctx(b"ACGTAC", b"CTTAC", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, Vec::new()),
                Run::unlabelled(2, 3, b"T".to_vec()),
            ],
        )
        .expect("sound seed");
        assert!(
            CodonFrameSeparation.apply(&seed, &c).is_none(),
            "a 5'UTR block is not a codon merge",
        );
    }

    #[test]
    fn codon_exception_merges_a_two_substitution_triplet_in_one_codon() {
        // Coding frame, grid phase 0. reference ACGTAC -> TCATAC: sub A>T at c.1
        // (offset 0) and sub G>A at c.3 (offset 2), one unchanged base (C at c.2)
        // between them; the triplet c.1_3 is one codon. The exception merges them into
        // c.1_3delinsTCA (the unchanged C carried interior).
        let frame = FrameContext::coding(0, Some(0), None);
        let prov = Provenance::none();
        let c = ctx(b"ACGTAC", b"TCATAC", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(2, 3, b"A".to_vec()),
            ],
        )
        .expect("sound seed");
        let out = CodonException.apply(&seed, &c).expect("merged");
        assert_eq!(out.runs().len(), 1);
        assert_eq!((out.runs()[0].ref_start, out.runs()[0].ref_end), (0, 3));
        assert_eq!(out.runs()[0].alt, b"TCA");
        assert!(out.phi() < seed.phi());
        assert!(
            out.runs().iter().all(|r| !r.locked),
            "an adapter never returns a locked run"
        );
    }

    #[test]
    fn codon_exception_declines_off_a_coding_frame() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        let c = ctx(b"ACGTAC", b"TCATAC", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(2, 3, b"A".to_vec()),
            ],
        )
        .expect("sound seed");
        assert!(CodonException.apply(&seed, &c).is_none());
    }

    #[test]
    fn codon_exception_respects_the_three_prime_utr_seam() {
        // The SAME triplet, but the CDS ends at offset 2 (cds_end=Some(2)), so the
        // triplet's third base sits past the last CDS base — in the 3'UTR. "Together
        // affecting one amino acid" is unstatable past the stop codon, so the merge
        // must NOT happen. This locks cds_end_axis: w_lo+2-1 = 2 makes within_cds(3)
        // false; an impl ignoring cds_end would wrongly merge.
        let frame = FrameContext::coding(0, Some(0), Some(2));
        let prov = Provenance::none();
        let c = ctx(b"ACGTAC", b"TCATAC", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(2, 3, b"A".to_vec()),
            ],
        )
        .expect("sound seed");
        assert!(
            CodonException.apply(&seed, &c).is_none(),
            "a triplet reaching past the stop codon stays two members",
        );
    }

    #[test]
    fn tandem_dup_run_agrees_with_the_raw_pass_and_fires() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // #2175: a CA tandem expansion (CACA -> CACACA) abutting a substitution
        // (g.28A>C). The minimal-edit derivation spreads the tandem copy across
        // single-base insertions and the dup vanishes; coalesce_by_run's
        // {peel; solid-run} restores g.[26_27dup;28A>C]. From
        // tests/it/issue_2175_dup_abutting_change.rs.
        let reference = b"TAGTAAACCATTTTACGGAGGATCACAAATTCCTCCTTAT";
        let alt = b"TAGTAAACCATTTTACGGAGGATCACACACATTCCTCCTTAT";
        let pieces = merge::partition_block_canonical(reference, alt).expect("canonical partition");
        let c = ctx(reference, alt, &frame, &prov);
        let seed = Cut::new(
            &c,
            pieces
                .iter()
                .map(|p| Run::unlabelled(p.ref_start, p.ref_end, p.alt.clone()))
                .collect(),
        )
        .expect("sound seed");

        // The raw combined per-run pass, exactly as the pipeline applies it.
        let mut raw = cut_to_pieces(&seed);
        merge::coalesce_by_run(&mut raw, reference, |run, reference| {
            merge::peel_tandem_dup_beside_change(run, reference, merge::PeelReach::TractOnly);
            merge::coalesce_solid_run(run, reference);
        });
        assert_ne!(
            raw,
            cut_to_pieces(&seed),
            "precondition: the combined pass fires"
        );

        let out = TandemDupRun
            .apply(&seed, &c)
            .expect("the peel/solid-run pass fires");
        assert_eq!(cut_to_pieces(&out), raw, "adapter geometry != raw pass");
        // Step 5: the peeled dup is TYPED `Dup` and LOCKED so the fixed point cannot
        // re-absorb it.
        let dup_runs: Vec<_> = out
            .runs()
            .iter()
            .filter(|r| r.label == Label::Dup)
            .collect();
        assert_eq!(dup_runs.len(), 1, "exactly one peeled dup is typed");
        assert!(dup_runs[0].locked, "the peeled dup is locked");
        assert!(
            merge::is_tandem_duplication(&as_piece(dup_runs[0]), reference),
            "the typed run is geometrically a tandem duplication",
        );
        // The Recut lift locks the whole peeled decomposition, so `Φ`'s unlocked count
        // strictly decreases even though the peel can raise the run count — which is
        // what makes this Recut usable inside the fixed-point engine.
        assert!(
            out.phi() < seed.phi(),
            "lift_recut's locks make Phi decrease: {:?} -> {:?}",
            seed.phi(),
            out.phi(),
        );

        // Two substitutions separated by an unchanged base: no tandem dup to peel and
        // no solid run to collapse (the gap breaks the run), so the pass declines.
        // (Adjacent subs would NOT decline — they are a solid run coalesce_solid_run
        // merges, which is the pass firing, not a no-op.)
        let plain = ctx(b"ACGT", b"TCAT", &frame, &prov);
        let plain_seed = Cut::new(
            &plain,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(2, 3, b"A".to_vec()),
            ],
        )
        .expect("sound seed");
        assert!(
            TandemDupRun.apply(&plain_seed, &plain).is_none(),
            "no tandem dup/solid run -> declines"
        );
    }

    #[test]
    fn run_inv_types_a_whole_span_reverse_complement() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // merge.rs's own `a_whole_span_revcomp_of_two_substitutions_types_as_inv`:
        // reference AAGCTA -> TAGCTT is an exact reverse complement, cut as two lone
        // substitutions [{0,1}T; {5,6}T] (the unchanged AGCT interior between them).
        // `whole-span-reverse-complement-types-as-inv` types the whole [0,6) span as
        // one inv, so the pass coalesces the two members into one spanning [0,6)member.
        let reference = b"AAGCTA";
        let resulting = denote(reference, &[(0, 1, b"T"), (5, 6, b"T")]);
        let c = ctx(reference, &resulting, &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(5, 6, b"T".to_vec()),
            ],
        )
        .expect("sound seed");

        // The raw pass in the exact block-relative frame the adapter uses.
        let mut raw = cut_to_pieces(&seed);
        merge::coalesce_inversion_runs(&mut raw, reference, 0, reference, &resulting);
        assert_ne!(
            raw,
            cut_to_pieces(&seed),
            "precondition: the raw pass fires"
        );

        let out = RunInv
            .apply(&seed, &c)
            .expect("the whole-span revcomp types as inv");
        assert_eq!(cut_to_pieces(&out), raw, "adapter geometry != raw pass");
        assert_eq!(
            out.runs().len(),
            1,
            "the two members coalesce to one inv span"
        );
        assert_eq!((out.runs()[0].ref_start, out.runs()[0].ref_end), (0, 6));
        // Step 5: the whole-span inversion is now TYPED `Inv` and LOCKED (via the
        // render stage's own `is_inversion` recogniser), which both makes the Recut
        // Phi-valid (1 unlocked -> 0 unlocked) and stops SepZero re-absorbing it.
        assert_eq!(out.runs()[0].label, Label::Inv, "the span is typed Inv");
        assert!(out.runs()[0].locked, "the inv is locked");
        assert!(
            merge::is_inversion(&as_piece(&out.runs()[0]), reference),
            "the typed run is geometrically an inversion",
        );
        assert_eq!(
            out.phi(),
            (0, 1),
            "one locked run: unlocked decreased 1 -> 0"
        );
    }

    /// #2161: a SUB-span reverse complement beside another change, handed to the
    /// adapter in the 3'-most placement the retired mirror never re-ran. The window
    /// `{insG@42, del@[47,48)}` has hull `[42,48)`; its difference hull `[42,47)`
    /// is `revcomp(TTAAC)`, so it types as one locked `Inv` and the substitution
    /// two unchanged bases away is carried forward untouched — separation two, so
    /// the F4 sep-zero decline does not fire.
    #[test]
    fn run_inv_types_a_sub_span_reverse_complement_off_the_three_prime_placement() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        let reference = b"AACCGGTTAATCGATCGATTGCACGTACGTGCAATCGATCGATTAACCGGTTAACCGGTTAACCGG";
        let resulting = denote(reference, &[(42, 42, b"G"), (47, 48, b""), (49, 50, b"T")]);
        let c = ctx(reference, &resulting, &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(42, 42, b"G".to_vec()),
                Run::unlabelled(47, 48, Vec::new()),
                Run::unlabelled(49, 50, b"T".to_vec()),
            ],
        )
        .expect("sound seed");
        let out = RunInv
            .apply(&seed, &c)
            .expect("the sub-span revcomp types off the 3' placement");
        assert_eq!(out.runs().len(), 2, "[inv; sub]");
        let inv = &out.runs()[0];
        assert_eq!((inv.ref_start, inv.ref_end), (42, 47));
        assert_eq!(&reference[42..47], b"TTAAC");
        assert_eq!(inv.alt, b"GTTAA", "revcomp(TTAAC)");
        assert_eq!(inv.label, Label::Inv);
        assert!(inv.locked);
        assert!(merge::is_inversion(&as_piece(inv), reference));
        let sub = &out.runs()[1];
        assert_eq!((sub.ref_start, sub.ref_end), (49, 50));
        assert_eq!(
            sub.label,
            Label::Unlabelled,
            "the untouched sub keeps its label"
        );
        assert!(!sub.locked);
        assert!(out.phi() < seed.phi(), "{:?} < {:?}", out.phi(), seed.phi());
        // And the recut is final: re-applying over its own locked output declines.
        assert!(RunInv.apply(&out, &c).is_none());
    }

    #[test]
    fn the_codon_exception_does_not_diverge_under_the_fixed_point_r2() {
        use crate::partition::ruled::run_engine;
        use crate::partition::rules::SepZero;
        // Design §12 R2: `CodonGapOne` is history-dependent by clause, so a fixed point
        // could in principle reach a different codon accumulation than the single
        // left-to-right pass on a chain of three codon-adjacent members. Structurally it
        // cannot: the codon exception merges a `[Sub; unchanged; Sub]` triplet whose two
        // subs lie in ONE codon, and two such triplets sharing a member would need three
        // GAPPED positions inside one 3-wide codon — impossible. So the exception cannot
        // chain, is idempotent on its own output (a merged triplet is a delins, not a
        // Sub), and the fixed point reaches the same accumulation as one pass. This pins
        // that, so the corpus's inability to build coding chains (R2/#1478) is not the
        // only evidence for the coding axis.
        let frame = FrameContext::coding(0, Some(0), None);
        let prov = Provenance::none();
        // Coding phase 0 (codons c.1_3, c.4_6). Three subs at c.1, c.3, c.5: {c.1,c.3}
        // share codon 1 (a codon exception, gap at c.2); c.5 is alone in codon 2.
        let reference = b"ACGTACGTAC";
        let resulting = denote(reference, &[(0, 1, b"T"), (2, 3, b"T"), (4, 5, b"C")]);
        let c = ctx(reference, &resulting, &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(2, 3, b"T".to_vec()),
                Run::unlabelled(4, 5, b"C".to_vec()),
            ],
        )
        .expect("sound");

        // Single pass: SepZero then CodonException, once each in the driver's order.
        let single = {
            let mut cut = seed.clone();
            for rule in [&SepZero as &dyn Rule, &CodonException] {
                if rule.scope().admits(&c) {
                    cut = rule.apply(&cut, &c).unwrap_or(cut);
                }
            }
            cut
        };
        // Fixed point over the same rule set.
        let rules: [&dyn Rule; 2] = [&SepZero, &CodonException];
        let fixed = run_engine(seed, &c, &rules);

        // Non-vacuity: the exception actually fired (merged the c.1_3 triplet).
        assert!(
            single.runs().len() < 3,
            "precondition: the codon exception fired"
        );
        // The merge is exactly the c.1_3 triplet; c.5 stays its own sub.
        assert_eq!(
            cut_to_pieces(&single),
            vec![
                Piece {
                    ref_start: 0,
                    ref_end: 3,
                    alt: b"TCT".to_vec()
                },
                Piece {
                    ref_start: 4,
                    ref_end: 5,
                    alt: b"C".to_vec()
                },
            ]
        );
        assert_eq!(
            cut_to_pieces(&fixed),
            cut_to_pieces(&single),
            "the fixed point reaches the same codon accumulation as one pass",
        );
    }

    #[test]
    fn a_run_inv_lock_survives_the_fixed_point_engine() {
        use crate::partition::ruled::run_engine;
        use crate::partition::rules::SepZero;
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // Whole-span inv (AAGCTA -> TAGCTT), seeded as two subs. Under the fixed-point
        // engine RunInv types+locks the span; on restart RunInv is a no-op on its own
        // output (Phi-safe: 1 locked run), and SepZero cannot re-absorb a locked run —
        // so the engine terminates with the `Inv` label and lock intact. This is the
        // step-5 property the whole label machinery exists for: the lock reaches
        // SepZero's guard through the restart-from-top loop.
        let reference = b"AAGCTA";
        let resulting = denote(reference, &[(0, 1, b"T"), (5, 6, b"T")]);
        let c = ctx(reference, &resulting, &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(5, 6, b"T".to_vec()),
            ],
        )
        .expect("sound");
        let rules: [&dyn Rule; 2] = [&RunInv, &SepZero];
        let normal = run_engine(seed, &c, &rules);
        assert_eq!(normal.runs().len(), 1, "the span is one inv member");
        assert_eq!(normal.runs()[0].label, Label::Inv);
        assert!(
            normal.runs()[0].locked,
            "the Inv lock survives to the engine's normal form"
        );
    }

    #[test]
    fn run_inv_declines_the_flanked_route_that_sep_zero_owns() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // A TRIMMED flanked inversion: reference GG|ACGTAA|GG whose payload is the
        // reverse complement of the interior [2,8) alone (`ACGTAA` -> `TTACGT`), the
        // flanks deleted. G != T at both ends, so trim_common_flanks leaves the whole
        // window as the block (lo = 0) — the shape merge.rs's own
        // `a_flanked_inversion_falls_out_as_del_inv_del` exercises, re-flanked so the
        // block is what the pipeline would actually hand the pass (the old TTACGTAAGG
        // fixture trims to a plain AAGG deletion). The raw pass INVENTS members
        // (1 -> 3, an adjacent [del; inv; del]) — but F4 retires that route: it places
        // an `inv` at separation zero from a reference-consuming `del`, which
        // `delins-adjacent-members-when-both-consume-reference` (SepZero) merges into
        // one delins. Emitting the split here would CYCLE with SepZero — the fixed
        // point oscillates and Φ rises, panicking `run_engine`'s progress assert — so
        // RunInv declines and leaves the merged delins to the rule that owns it.
        let reference = b"GGACGTAAGG";
        let core = b"TTACGT"; // rc(reference[2..8]) == rc("ACGTAA")
        let resulting = denote(reference, &[(0, 10, core)]);
        assert_eq!(
            resulting, core,
            "the payload is the inverted interior alone"
        );
        let c = ctx(reference, &resulting, &frame, &prov);
        let seed = Cut::new(&c, vec![Run::unlabelled(0, 10, core.to_vec())]).expect("sound seed");

        // Precondition: the underlying pass fires and produces the adjacent
        // [del; inv; del] (1 -> 3) — the geometry F4 declines to emit. Sabotage: delete
        // the flanked-route decline in `RunInv::apply` and this test's `is_none()`
        // becomes a `Some` carrying exactly these three members.
        let mut raw = cut_to_pieces(&seed);
        merge::coalesce_inversion_runs(&mut raw, reference, 0, reference, &resulting);
        assert_ne!(
            raw,
            cut_to_pieces(&seed),
            "precondition: the raw pass fires"
        );
        assert_eq!(
            raw.len(),
            3,
            "the raw pass invents an adjacent [del; inv; del]"
        );
        assert!(
            raw.windows(2).all(|w| w[0].ref_end == w[1].ref_start),
            "the raw members are adjacent (separation zero), which is SepZero's to merge",
        );

        // F4: RunInv declines the flanked route rather than cycle with SepZero.
        assert!(
            RunInv.apply(&seed, &c).is_none(),
            "a flanked [del; inv; del] is SepZero's delins, so RunInv declines it",
        );
    }

    #[test]
    fn run_inv_reads_the_trimmed_block_inside_a_padded_window() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // A whole-span reverse complement (AAGCTA -> TAGCTT, the block of
        // `run_inv_types_a_whole_span_reverse_complement`) padded CC..CC on both sides:
        // the pipeline calls the pass with the WINDOW, block_lo = 2, and the block
        // slices. The adapter must derive (lo, hi_ref, hi_alt) = (2, 8, 8) itself and
        // pass exactly that. A whole-span revcomp coalesces to ONE inv run — no
        // flanking del — so F4's flanked-route decline does not apply and RunInv fires
        // positively. This is the padded-window POSITIVE case (its unpadded sibling is
        // `run_inv_types_a_whole_span_reverse_complement`).
        let window = b"CCAAGCTACC";
        let resulting = denote(window, &[(2, 3, b"T"), (7, 8, b"T")]);
        assert_eq!(resulting, b"CCTAGCTTCC");
        let c = ctx(window, &resulting, &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(2, 3, b"T".to_vec()),
                Run::unlabelled(7, 8, b"T".to_vec()),
            ],
        )
        .expect("sound");

        // In the block frame the pass fires (AAGCTA -> TAGCTT is a revcomp), coalescing
        // the two subs into one inv [2,8).
        let mut raw = cut_to_pieces(&seed);
        merge::coalesce_inversion_runs(&mut raw, window, 2, &window[2..8], &resulting[2..8]);
        assert_ne!(
            raw,
            cut_to_pieces(&seed),
            "precondition: the raw pass fires in the block frame"
        );
        // The whole-window call — what the adapter would do if it did NOT trim — is a
        // no-op: the whole padded window is not a reverse complement, so route 0 never
        // fires. Had the adapter passed the whole window, RunInv would DECLINE (nothing
        // changed), so the positive fire below is itself proof that the adapter trimmed.
        let mut whole = cut_to_pieces(&seed);
        merge::coalesce_inversion_runs(&mut whole, window, 0, window, &resulting);
        assert_eq!(
            whole,
            cut_to_pieces(&seed),
            "the whole-window (untrimmed) call is a no-op — it is not a revcomp",
        );

        let out = RunInv
            .apply(&seed, &c)
            .expect("the whole-span revcomp in a padded window fires");
        assert_eq!(
            cut_to_pieces(&out),
            raw,
            "adapter must pass (window, lo, block_ref, block_alt)"
        );
        assert_eq!(
            out.runs().len(),
            1,
            "the whole-span revcomp coalesces to one inv"
        );
        assert_eq!((out.runs()[0].ref_start, out.runs()[0].ref_end), (2, 8));
        assert_eq!(out.runs()[0].label, Label::Inv, "typed Inv");
        assert!(out.runs()[0].locked, "the inv is locked");
    }

    #[test]
    fn placed_gap_collapses_to_the_block_span_inside_a_padded_window() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // #1610's TGCA -> AAC, padded GG..GG: trim gives (2, 6, 5). The pipeline
        // evaluates the predicate on BLOCK-frame pieces over the block and collapses
        // to one whole-block piece, which in window coordinates is [2,6) -> "AAC" —
        // never [0, window.len()).
        let window = b"GGTGCAGG";
        let resulting = b"GGAACGG";
        let c = ctx(window, resulting, &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(2, 4, b"AA".to_vec()),
                Run::unlabelled(5, 6, Vec::new()),
            ],
        )
        .expect("sound seed");
        let out = PlacedGap.apply(&seed, &c).expect("collapsed");
        assert_eq!(out.runs().len(), 1);
        assert_eq!((out.runs()[0].ref_start, out.runs()[0].ref_end), (2, 6));
        assert_eq!(out.runs()[0].alt, b"AAC");
    }

    #[test]
    fn run_inv_declines_an_ambiguous_flank_split() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // merge.rs's `an_ambiguous_flank_split_is_refused`: `TTTT` reverse-complements
        // to `AAAA`, so a `TTTT` payload over `AAAAAA` matches revcomp at three offsets.
        // The pass refuses rather than pick one, so the rule does not fire.
        let reference = b"AAAAAA";
        let resulting = denote(reference, &[(0, 6, b"TTTT")]);
        let c = ctx(reference, &resulting, &frame, &prov);
        let seed = Cut::new(&c, vec![Run::unlabelled(0, 6, b"TTTT".to_vec())]).expect("sound seed");
        assert!(
            RunInv.apply(&seed, &c).is_none(),
            "three equally-compliant flank splits -> the pass refuses, so the rule declines",
        );
    }

    #[test]
    fn split_concealed_cuts_a_member_hiding_a_forced_unchanged_base() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // merge.rs's `a_member_concealing_a_forced_unchanged_base_is_cut_apart`
        // (#1539): reference GGTCGGTCC -> ACGTGG, a net deletion partitioned as
        // [{0,5}AC; {7,9}GG]. The first member `GGTCG -> AC` hides the C at offset 3
        // (unchanged in every minimal alignment), so the audit cuts it into
        // [{0,3}A; {4,5}del], giving [{0,3}A; {4,5}del; {7,9}GG] (2 -> 3 pieces).
        // Frameless here (NonCoding), where `general.md:34` stands with no codon
        // exception to spare the base — merge.rs's frameless arm cuts it just the same.
        let reference = b"GGTCGGTCC";
        let resulting = denote(reference, &[(0, 5, b"AC"), (7, 9, b"GG")]);
        assert_eq!(resulting, b"ACGTGG", "the seed denotes the #1539 block");
        let c = ctx(reference, &resulting, &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 5, b"AC".to_vec()),
                Run::unlabelled(7, 9, b"GG".to_vec()),
            ],
        )
        .expect("sound seed");

        // The raw pass with the adapter's own derived args: reading_frame=false
        // (NonCoding), length_changing=true (9 != 6), w_lo=0 (cds_axis_origin of a
        // non-coding frame).
        let mut raw = cut_to_pieces(&seed);
        merge::split_concealed_separations(&mut raw, false, true, 0, reference);
        assert_ne!(
            raw,
            cut_to_pieces(&seed),
            "precondition: the raw pass fires"
        );

        let out = SplitConcealed
            .apply(&seed, &c)
            .expect("the concealed separation is cut apart");
        assert_eq!(cut_to_pieces(&out), raw, "adapter geometry != raw pass");
        assert_eq!(
            out.runs().len(),
            3,
            "the hidden separation becomes a gap (2 -> 3)"
        );
        assert!(
            out.phi() > seed.phi(),
            "a split increases the run count, so Phi rises (Recut)"
        );
        assert!(
            out.runs().iter().all(|r| !r.locked),
            "a split emits plain geometry: no label, no lock",
        );
    }

    #[test]
    fn split_concealed_declines_a_length_neutral_block() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // Two substitutions on an equal-length block (ACGT -> TCAT). The pass audits
        // only length-changing partitions (`!length_changing` early-returns), and the
        // adapter derives length_changing = ref.len() != resulting.len() = false, so
        // the rule declines — matching merge.rs's own equal-length exclusion.
        let reference = b"ACGT";
        let resulting = denote(reference, &[(0, 1, b"T"), (2, 3, b"A")]);
        assert_eq!(resulting, b"TCAT");
        assert_eq!(
            reference.len(),
            resulting.len(),
            "the block is length-neutral"
        );
        let c = ctx(reference, &resulting, &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(2, 3, b"A".to_vec()),
            ],
        )
        .expect("sound seed");
        assert!(
            SplitConcealed.apply(&seed, &c).is_none(),
            "a length-neutral block is not audited, so the rule declines",
        );
    }

    #[test]
    fn frame_context_inverts_cds_axis_origin_up_to_a_codon() {
        // `same_codon` is invariant to a w_lo shift of 3 and rejects pos < 1, so the
        // inverse must reproduce w_lo modulo 3 with every origin >= 1 — and exactly at
        // the 5' seam (w_lo <= 1), where cds_start pins the grid.
        for w_lo in 1..=12 {
            let f = frame_context_for_axis_origin(true, w_lo, None);
            let origin = cds_axis_origin(&f);
            assert!(origin >= 1, "w_lo={w_lo} gave origin {origin}");
            assert_eq!((origin - w_lo).rem_euclid(3), 0, "w_lo={w_lo}");
        }
        assert_eq!(
            cds_axis_origin(&frame_context_for_axis_origin(true, 1, None)),
            1
        );
        assert_eq!(
            cds_axis_origin(&frame_context_for_axis_origin(true, 0, None)),
            0
        );
        assert_eq!(
            cds_axis_origin(&frame_context_for_axis_origin(true, -2, None)),
            -2
        );
    }

    #[test]
    fn frame_context_preserves_the_window_offset_of_the_cds_end() {
        // CodonException derives cds_end_axis' = origin + cds_end - 1; `within_cds`
        // compares origin + offset <= cds_end_axis', i.e. offset <= cds_end - 1. The
        // pipeline's test is w_lo + offset <= cds_end_axis, i.e. offset <= e - w_lo.
        // So cds_end - 1 must equal (e - w_lo), clamped at -1 when the window lies
        // wholly past the CDS end (every offset is then 3'UTR).
        for (w_lo, e) in [(1i64, 12i64), (5, 12), (7, 7), (9, 4)] {
            let f = frame_context_for_axis_origin(true, w_lo, Some(e));
            let FrameContext::Coding {
                cds_end: Some(end), ..
            } = f
            else {
                panic!("coding with a 3' seam");
            };
            assert_eq!(end as i64 - 1, (e - w_lo).max(-1), "w_lo={w_lo} e={e}");
        }
    }

    #[test]
    fn frame_context_is_non_coding_without_a_reading_frame() {
        assert_eq!(
            frame_context_for_axis_origin(false, 5, Some(9)),
            FrameContext::NonCoding
        );
    }

    #[test]
    fn rule_migration_inventory_classifies_the_landed_rules() {
        // The exemplar adapter is recorded as an adapter over its shipping pass, and
        // the odd-one-out native rule as native — the two shapes step 1 is allowed.
        let by_id: std::collections::HashMap<_, _> = RULE_MIGRATION.iter().copied().collect();
        assert_eq!(
            by_id.get(PayloadCoincidence.id()),
            Some(&Migration::Adapter {
                wraps: "coalesce_payload_alignment_split"
            }),
        );
        assert_eq!(by_id.get("sep-zero"), Some(&Migration::Native));
        // The two lock-touching passes resolved as geometry adapters at step 1.
        assert_eq!(
            by_id.get(RunInv.id()),
            Some(&Migration::Adapter {
                wraps: "coalesce_inversion_runs"
            }),
        );
        assert_eq!(
            by_id.get(SplitConcealed.id()),
            Some(&Migration::Adapter {
                wraps: "split_concealed_separations"
            }),
        );
    }
}
