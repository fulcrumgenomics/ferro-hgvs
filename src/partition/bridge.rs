//! The `Piece`↔`Cut` bridge — the enabling primitive for the **adapter-first**
//! step-1 migration (Fable design §10, ruling (a)).
//!
//! An *adapter* rule ([`crate::partition::rules::adapters`]) reproduces a Stage-B
//! `merge.rs` pass byte-for-byte by *calling that pass unchanged*, not by
//! transcribing it: it maps the [`Cut`] it is handed down to the `Vec<Piece>` the
//! legacy pass mutates ([`cut_to_pieces`]), runs the pass, and maps the result back
//! up to a validated [`Cut`] ([`pieces_to_cut`]). Because the pass is the shipping
//! code executed, the adapter is byte-identical to baseline *definitionally* rather
//! than empirically — which is the whole reason (a) beats a hand-rewrite for the
//! label-blind passes.
//!
//! # This is a confined, deliberate breach of `partition`'s independence
//!
//! [`crate::partition`]'s module doc states that only [`crate::partition::registry`]
//! reaches into `normalize::merge` internals. This module is the **second** such
//! confined reacher, and it exists for the same reason: the migration cannot wrap a
//! `merge`-private pass without touching `merge`'s [`Piece`]. Like `registry`, the
//! coupling is one narrow, dev-scoped seam ([`Piece`] plus the handful of
//! `pub(crate)` label-blind passes the adapters call), not a general dependency — and
//! it is temporary: step 3 retires the adapters as the passes move native.
//!
//! # Two lifts: the label-blind one, and the typed one (step 5)
//!
//! A [`Piece`] carries no label and no lock bit, so the *pass* is label-blind: what a
//! `merge.rs` pass produces is pure geometry. [`pieces_to_cut`] and [`lift`] reflect
//! that — they mint [`Run::unlabelled`], unlocked runs, and `bridge_loses_labels_and_locks`
//! pins that loss as intended for the label-blind path. Through step 4 that was the
//! only lift, so an adapter could not return a locked `Cut` and any lock-setting rule
//! was forced native.
//!
//! Step 5 adds [`lift_relabel`], which types the pass's OUTPUT: it re-derives each
//! run's label from its geometry with the render stage's own recognisers
//! ([`merge::is_inversion`], [`merge::is_tandem_duplication`]) and carries an earlier
//! rule's label forward on unchanged runs. So an adapter CAN now set an `Inv`/`Dup`
//! lock — the pass supplies the geometry, the typer the label. What still forces a
//! rule native is needing to *read* a lock mid-pass, which only
//! [`crate::partition::rules::SepZero`] does; setting one no longer does.

use crate::normalize::merge::{self, Piece};
use crate::normalize::ShuffleDirection;
use crate::partition::block_ctx::BlockCtx;
use crate::partition::partitioner::PartitionError;
use crate::partition::ruled::{Cut, Label, Run};
use std::cell::Cell;

/// Project a [`Cut`] down to the `Vec<Piece>` a legacy `merge.rs` pass mutates.
///
/// Labels and lock bits are **dropped** — a [`Piece`] cannot carry them. The runs
/// are already ascending and disjoint (a [`Cut`] guarantees it), which is the order
/// every pass expects, so no re-sort is needed.
pub(crate) fn cut_to_pieces(cut: &Cut) -> Vec<Piece> {
    cut.runs()
        .iter()
        .map(|r| Piece {
            ref_start: r.ref_start,
            ref_end: r.ref_end,
            alt: r.alt.clone(),
        })
        .collect()
}

/// Whether `pieces` still equals `cut_to_pieces(cut)` — i.e. a `merge.rs` pass left
/// the cut's geometry untouched — computed without materialising the comparison
/// `Vec`. Exactly equivalent to `pieces == cut_to_pieces(cut)`: a [`Piece`] carries
/// only `(ref_start, ref_end, alt)` and derives `PartialEq` structurally, and
/// [`cut_to_pieces`] copies each of those three from the corresponding [`Run`] and
/// nothing else — so comparing `pieces` against `cut.runs()` field by field answers
/// the same question and allocates nothing. Adapters use it as the no-op guard in
/// place of a `let before = pieces.clone()` snapshot.
pub(crate) fn pieces_match_cut(pieces: &[Piece], cut: &Cut) -> bool {
    let runs = cut.runs();
    pieces.len() == runs.len()
        && pieces
            .iter()
            .zip(runs)
            .all(|(p, r)| p.ref_start == r.ref_start && p.ref_end == r.ref_end && p.alt == r.alt)
}

/// Lift a `Vec<Piece>` a legacy pass produced back to a validated [`Cut`].
///
/// Every piece becomes an [`Run::unlabelled`] run — so the returned cut is entirely
/// `Unlabelled` and unlocked, which is the bridge's structural label-lock guardrail
/// (see the module doc). [`Cut::new`] re-checks soundness and renderability against
/// `ctx`, so a pass that somehow broke the splice is caught here rather than shipped;
/// an adapter propagates that [`PartitionError`] by declining (`.ok()`), never panics.
pub(crate) fn pieces_to_cut(ctx: &BlockCtx, pieces: Vec<Piece>) -> Result<Cut, PartitionError> {
    let runs: Vec<Run> = pieces
        .into_iter()
        .map(|p| Run::unlabelled(p.ref_start, p.ref_end, p.alt))
        .collect();
    Cut::new(ctx, runs)
}

thread_local! {
    /// Piece lists this thread could not lift to a `Cut` (design step-1 seam).
    /// Per-thread because a driver reads it as a delta around ONE block.
    static LIFT_REFUSALS: Cell<u64> = const { Cell::new(0) };
}

/// `pieces_to_cut` that **records** a refusal instead of swallowing it.
///
/// The coincident-insertion case this seam was built for is **retired** as of
/// design §10 step 2: `Cut::new` now folds two content-bearing insertions at one
/// interbase into the single `Ins` they denote, so `pieces_to_cut` no longer
/// refuses that shape — the peel's `#2201` intermediate and the flush-shift pair
/// both lift and fold rather than declining. What remains is the residual guard: a
/// piece list a pass leaves genuinely `Unsound` (it does not splice back to
/// `ctx.resulting`) still refuses here, is counted, and the driver declines the
/// block (`Err(Unrenderable)`) so the wiring falls back to the shipping chain —
/// byte-identical by construction. Measured `0` over the corpus either way.
pub(crate) fn lift(ctx: &BlockCtx, pieces: Vec<Piece>) -> Option<Cut> {
    match pieces_to_cut(ctx, pieces) {
        Ok(cut) => Some(cut),
        Err(_) => {
            LIFT_REFUSALS.with(|c| c.set(c.get() + 1));
            None
        }
    }
}

/// Lift refusals recorded on this thread so far (monotone). Read before and after
/// a chain; a difference means some intermediate was unrepresentable.
pub(crate) fn lift_refusals() -> u64 {
    LIFT_REFUSALS.with(Cell::get)
}

/// `lift`, but each output piece is TYPED — so a label a rule set survives the
/// pass's `cut_to_pieces`/mutate/lift round-trip, which plain [`lift`] erases (it
/// mints every run `Unlabelled`). This is the step-5 machinery that lets `Dup`/`Inv`
/// locks reach [`crate::partition::rules::SepZero`]'s lock-guard under the fixed
/// point (design §4.2). Each piece's label is resolved, in priority order:
///
/// 1. `typer(piece, ctx.reference)` — a geometry recogniser (`is_inversion`,
///    `is_tandem_duplication`). A recognised run is typed and **locked iff the label
///    locks** ([`Label::locks`]), which is how `RunInv`/`TandemDupRun` set their
///    `Inv`/`Dup` locks over exactly the render stage's own geometry (§4.8).
/// 2. an incoming run with identical geometry (span and alt bytes) — its label and
///    lock are carried forward, so a coarsening that leaves a run untouched preserves
///    an earlier rule's typing.
/// 3. otherwise `Unlabelled`/unlocked — a freshly-produced run, exactly `lift`'s mint.
///
/// A [`Cut`]'s runs are disjoint and coincident insertions are folded, so the
/// geometry match in (2) is unambiguous. Records a refusal like [`lift`] on an
/// unsound intermediate.
pub(crate) fn lift_relabel(
    ctx: &BlockCtx,
    incoming: &Cut,
    pieces: Vec<Piece>,
    typer: impl Fn(&Piece, &[u8]) -> Option<crate::partition::ruled::Label>,
) -> Option<Cut> {
    lift_typed(ctx, incoming, pieces, typer, false)
}

/// A Recut's lift: like [`lift_relabel`], but a NEWLY-produced run — one matching no
/// incoming run's geometry and not recognised by the typer — is `Unlabelled` yet
/// **locked**, not unlocked. This is what makes a Recut strictly decrease `Φ`'s
/// unlocked count under the fixed point (design §4.5: a re-cut turns one unlocked run
/// into `n` LOCKED runs), including the flanked-inversion route `[del; inv; del]`
/// whose flanking `del`s are the Recut's committed output and must not be re-coarsened
/// even though they are not themselves inversions.
///
/// Byte-identical to [`lift_relabel`] in the single-pass driver: nothing reads a lock
/// there except `SepZero`, which runs before the Recuts, and the render re-derives
/// kinds from geometry — so the extra locks change no output and matter only under the
/// fixed point, where they both keep `Φ` decreasing and stop `SepZero` re-absorbing
/// the Recut's members.
pub(crate) fn lift_recut(
    ctx: &BlockCtx,
    incoming: &Cut,
    pieces: Vec<Piece>,
    typer: impl Fn(&Piece, &[u8]) -> Option<crate::partition::ruled::Label>,
) -> Option<Cut> {
    // A Recut is final once it locks its output, so it declines to re-cut a locked run —
    // on restart-from-top the pass re-runs over its own locked output, and without this
    // it would re-peel/re-type it every iteration (`Φ` stalls, the debug_assert fires).
    // On the first, legitimate application the input is unlocked, so nothing is absorbed.
    if absorbed_a_lock(incoming, &pieces) {
        return None;
    }
    lift_typed(ctx, incoming, pieces, typer, true)
}

/// The shared body of [`lift_relabel`] and [`lift_recut`]. Per piece: a typer-recognised
/// run is typed and locked iff its label locks; an unchanged incoming run keeps its
/// label and lock; an unmatched (new) run is `Unlabelled` and locked iff `lock_new`.
fn lift_typed(
    ctx: &BlockCtx,
    incoming: &Cut,
    pieces: Vec<Piece>,
    typer: impl Fn(&Piece, &[u8]) -> Option<crate::partition::ruled::Label>,
    lock_new: bool,
) -> Option<Cut> {
    let runs: Vec<Run> =
        pieces
            .into_iter()
            .map(|p| {
                if let Some(label) = typer(&p, ctx.reference) {
                    return Run {
                        ref_start: p.ref_start,
                        ref_end: p.ref_end,
                        alt: p.alt,
                        locked: label.locks(),
                        label,
                    };
                }
                match incoming.runs().iter().find(|r| {
                    r.ref_start == p.ref_start && r.ref_end == p.ref_end && r.alt == p.alt
                }) {
                    Some(r) => Run {
                        ref_start: p.ref_start,
                        ref_end: p.ref_end,
                        alt: p.alt,
                        label: r.label,
                        locked: r.locked,
                    },
                    None => Run {
                        ref_start: p.ref_start,
                        ref_end: p.ref_end,
                        alt: p.alt,
                        label: crate::partition::ruled::Label::Unlabelled,
                        locked: lock_new,
                    },
                }
            })
            .collect();
    match Cut::new(ctx, runs) {
        Ok(cut) => Some(cut),
        Err(_) => {
            LIFT_REFUSALS.with(|c| c.set(c.get() + 1));
            None
        }
    }
}

/// The label-preserving typer that recognises no new geometry — carries incoming
/// labels forward and mints the rest `Unlabelled`. The coarsening adapters lift
/// through this so a `Dup`/`Inv` label set by an earlier rule is not erased when a
/// later coarsening leaves that run untouched.
pub(crate) fn preserve_only(_p: &Piece, _r: &[u8]) -> Option<crate::partition::ruled::Label> {
    None
}

/// A coarsening's lift: [`lift_relabel`] with [`preserve_only`], but it DECLINES
/// (returns `None`) if the pass ABSORBED a locked incoming run.
///
/// A locked `Dup`/`Inv` outranks a coarsening (design §4.2), but a `merge.rs` coarsen
/// pass operates on `Piece`s and is blind to the lock — so the adapter enforces the
/// lock here, after the fact: a locked run is "absorbed" when its exact geometry (span
/// **and** alt bytes) is not among the output pieces, and if any is, the whole
/// coarsening is refused so the run survives. `SepZero` (native) reads the lock
/// directly; this is the same guard for the label-blind adapters, which is what stops
/// the fixed point re-coarsening an `Inv` a `RunInv` just typed (the
/// `delins_hiding_an_inversion` reversion) or a peeled `Dup`.
///
/// Byte-identical in the single-pass driver: the payload/compensating coarsens run
/// before any Recut locks a run, and the one coarsen that runs after them
/// (`CodonException`) merges a `[Sub; unchanged; Sub]` triplet, none of whose members
/// is a zero-width `Dup` or a `≥2`-base `Inv`, so it can never absorb a lock. The guard
/// therefore only ever fires under the fixed point.
pub(crate) fn lift_coarsen(ctx: &BlockCtx, incoming: &Cut, pieces: Vec<Piece>) -> Option<Cut> {
    if absorbed_a_lock(incoming, &pieces) {
        // Not a refusal (the intermediate is representable) — the coarsening simply
        // may not fire across a lock, so the rule declines.
        return None;
    }
    lift_relabel(ctx, incoming, pieces, preserve_only)
}

/// Whether a pass's output `pieces` dropped a LOCKED incoming run — i.e. a locked
/// `Dup`/`Inv` whose exact geometry (span + alt bytes) is not among the pieces. A
/// locked run is final (design §4.2), so any adapter whose pass would absorb one must
/// decline; both [`lift_coarsen`] and [`lift_recut`] gate on this. For a Recut it is
/// what makes the pass fire at most once per run: on restart the pass re-runs over its
/// own now-LOCKED output and, if it would re-cut it, this refuses — so `Φ` cannot stall
/// on a rule endlessly re-cutting what it already locked.
fn absorbed_a_lock(incoming: &Cut, pieces: &[Piece]) -> bool {
    incoming.runs().iter().filter(|r| r.locked).any(|locked| {
        !pieces.iter().any(|p| {
            p.ref_start == locked.ref_start && p.ref_end == locked.ref_end && p.alt == locked.alt
        })
    })
}

/// `merge::shift_pieces` over a `Cut`: the 3'/5' placement of every pure indel
/// within the window, bounded by its neighbours. Not a rule (Φ-neutral,
/// direction-parameterised); the driver's placement step. Total: on a lift refusal
/// the input is returned and the refusal stands recorded.
pub(crate) fn shift(cut: Cut, ctx: &BlockCtx, direction: ShuffleDirection) -> Cut {
    let mut pieces = cut_to_pieces(&cut);
    merge::shift_pieces(&mut pieces, ctx.reference, direction);
    lift(ctx, pieces).unwrap_or(cut)
}

/// Place every run in `direction` — the final, per-run 3'/5' placement (design
/// §4.7's step-7 placement) — **preserving each run's label and lock**.
///
/// Unlike [`shift`], which drops labels through the label-blind `Piece` bridge,
/// this re-attaches them by index: [`merge::shift_pieces`] moves each run within
/// its own walls without changing the run count, the run order, or any payload, so
/// the *i*-th placed piece is the *i*-th input run with its coordinates moved. Two
/// properties follow, and they are what let step 7 retire the direction mirror:
/// placement is **label-preserving** (a locked `Inv`/`Dup` stays locked and typed
/// through it) and **member-count-preserving** (a placement can never change how
/// many members there are). So the 5' arm is the 3' partition re-placed 5', and the
/// #1542 direction-symmetry of the member count holds by construction rather than
/// by mirroring.
pub(crate) fn place_only(cut: Cut, ctx: &BlockCtx, direction: ShuffleDirection) -> Cut {
    let input_runs = cut.runs().to_vec();
    let labels: Vec<_> = input_runs.iter().map(|r| (r.label, r.locked)).collect();
    let mut pieces = cut_to_pieces(&cut);
    merge::shift_pieces(&mut pieces, ctx.reference, direction);
    debug_assert_eq!(
        pieces.len(),
        labels.len(),
        "a placement shift never changes the run count",
    );
    let mut runs: Vec<Run> = pieces
        .into_iter()
        .zip(labels)
        .map(|(p, (label, locked))| Run {
            ref_start: p.ref_start,
            ref_end: p.ref_end,
            alt: p.alt,
            label,
            locked,
        })
        .collect();
    // F2b (#1542): a `Dup` is placed by WHOLE UNITS, not the base-by-base insertion
    // shuffle `shift_pieces` applied. `shift_pieces` treats a dup as a generic
    // insertion and can roll it mid-tract, off any valid dup slot, so it renders as
    // a plain `ins` — which the downstream `merge_consecutive_edits` then folds into
    // a `delins`, changing the MEMBER COUNT with direction (`5' g.[21del;24delinsAAC]`
    // vs `3' g.[21del;24G>A;25_26dup]`). A dup must stay a dup: re-place it from its
    // canonical (3'-most, the input) insertion interbase — 3' keeps that slot; a
    // LONE dup (no member abutting its 3' side) walks 5' to the 5'-most copy of the
    // unit in the reference tandem tract. The payload is the reference unit verbatim
    // (unit steps never rotate it), so the resulting sequence, the label and the
    // count are all preserved by construction;
    // `separation-is-a-property-of-the-spelling-not-of-the-variant` requires exactly
    // this direction-independence of the count.
    //
    // A dup that was PEELED beside a change (`peel_tandem_dup_beside_change`, #2175)
    // is a different case: the peel's descending scan places it 3'-most, ABUTTING the
    // change, and takes no direction — the peeled shape is direction-independent
    // (`issue_2175_dup_abutting_change`, `fragmentation_corpus`). Walking such a dup
    // 5' would detach it from the change and give the 5' arm a different member
    // spelling than the 3' arm. So a dup that abuts a 3' sibling (a following run
    // starting at the dup's interbase) is frozen at its 3'-most slot; only a lone dup
    // shuffles.
    // F2b sibling clamp: track the highest reference offset a preceding run has
    // claimed, exactly as the render seam's `dup_source_overlaps_prior_piece`
    // does (`previous_ref_end`), so the 5' whole-unit walk below cannot step a
    // dup onto a source a sibling overwrote. `runs` is ascending and
    // non-overlapping (a placement preserves both), so this reassigned-each-run
    // value is the running maximum the renderer sees.
    let mut previous_ref_end = 0usize;
    for i in 0..runs.len() {
        if matches!(runs[i].label, Label::Dup) {
            let input = &input_runs[i];
            let k = input.alt.len();
            let mut p = input.ref_start; // == ref_end: a dup is a zero-width insertion
            debug_assert_eq!(input.ref_start, input.ref_end, "a Dup run is zero-width");
            // A dup peeled beside a change abuts a following run that starts at the
            // dup's interbase; such a dup is frozen at its 3'-most slot (see above).
            let abuts_three_prime_sibling = input_runs
                .get(i + 1)
                .is_some_and(|next| next.ref_start == input.ref_start);
            if k > 0 && direction == ShuffleDirection::FivePrime && !abuts_three_prime_sibling {
                // Step 5' by whole units while another copy of the unit sits just
                // before the current source `[p-k, p)` AND the stepped-to source
                // `[p-2k, p-k)` does not reach into a sibling run's territory
                // (`previous_ref_end`). Without the clamp the dup can walk onto a
                // source a preceding edit overwrote, which the render seam then
                // refuses (`dup_source_overlaps_prior_piece`) and falls back to
                // the per-member pipeline — the s02 5' F2b defect. The 3'-most
                // input slot is sound by construction, and stepping 5' only
                // lowers `p`, so the clamp keeps every placed dup renderable.
                while p >= 2 * k
                    && ctx.reference[p - 2 * k..p - k] == input.alt[..]
                    && p - 2 * k >= previous_ref_end
                {
                    p -= k;
                }
            }
            runs[i].ref_start = p;
            runs[i].ref_end = p;
            runs[i].alt = input.alt.clone();
        }
        previous_ref_end = runs[i].ref_end;
    }
    // A placement preserves the splice, so lifting back cannot fail; the fallback
    // to the input is belt-and-braces and never a silent change of geometry.
    Cut::new(ctx, runs).unwrap_or(cut)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::block_ctx::{FrameContext, Molecule, Provenance};
    use crate::partition::ruled::Label;

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
    fn pieces_match_cut_agrees_with_a_materialised_compare() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // A sound single-substitution cut: ACGT -> AGGT.
        let c = ctx(b"ACGT", b"AGGT", &frame, &prov);
        let cut = Cut::new(&c, vec![Run::unlabelled(1, 2, b"G".to_vec())]).expect("sound sub");

        let identical = cut_to_pieces(&cut);
        let mut changed_alt = identical.clone();
        changed_alt[0].alt = b"T".to_vec();
        let mut changed_span = identical.clone();
        changed_span[0].ref_end = 3;
        let longer = {
            let mut v = identical.clone();
            v.push(Piece {
                ref_start: 3,
                ref_end: 3,
                alt: b"A".to_vec(),
            });
            v
        };
        let empty: Vec<Piece> = vec![];

        for candidate in [&identical, &changed_alt, &changed_span, &longer, &empty] {
            assert_eq!(
                pieces_match_cut(candidate, &cut),
                *candidate == cut_to_pieces(&cut),
                "pieces_match_cut must equal a materialised compare for {candidate:?}",
            );
        }
        assert!(
            pieces_match_cut(&identical, &cut),
            "the unchanged projection matches the cut",
        );
        assert!(
            !pieces_match_cut(&changed_alt, &cut),
            "a changed payload does not match, even at equal length",
        );

        // A two-run cut, with a perturbation on the SECOND piece only — guards the
        // zip against an ordering or early-exit slip that a one-run cut cannot see.
        let c2 = ctx(b"ACGTACGT", b"AXGTAYGT", &frame, &prov);
        let cut2 = Cut::new(
            &c2,
            vec![
                Run::unlabelled(1, 2, b"X".to_vec()),
                Run::unlabelled(5, 6, b"Y".to_vec()),
            ],
        )
        .expect("sound two-sub cut");
        let base = cut_to_pieces(&cut2);
        let mut second_differs = base.clone();
        second_differs[1].alt = b"Z".to_vec();
        assert!(
            pieces_match_cut(&base, &cut2),
            "the unchanged two-run projection matches",
        );
        assert_eq!(
            pieces_match_cut(&second_differs, &cut2),
            second_differs == cut_to_pieces(&cut2),
            "a change to the second piece is caught, matching a materialised compare",
        );
        assert!(!pieces_match_cut(&second_differs, &cut2));
    }

    #[test]
    fn place_only_keeps_a_dup_on_its_slot_in_both_directions() {
        // F2b (#1542): a `Dup` is placed by whole units, so it stays a dup in both
        // directions instead of `shift_pieces` rolling it off its tandem slot into
        // a plain insertion (which the render then folds, changing the member count
        // with direction). reference TACAC -> TACACAC: an "AC" dup whose tandem tract
        // is AC at [1,3) and [3,5); the 3'-most dup inserts at interbase 5, the
        // 5'-most at interbase 3, both denoting TACACAC.
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        let c = ctx(b"TACAC", b"TACACAC", &frame, &prov);
        let dup = Run {
            ref_start: 5,
            ref_end: 5,
            alt: b"AC".to_vec(),
            label: Label::Dup,
            locked: true,
        };
        let cut = Cut::new(&c, vec![dup]).expect("sound 3'-most dup");
        // 3' keeps the canonical (3'-most) slot.
        let three = place_only(cut.clone(), &c, ShuffleDirection::ThreePrime);
        assert_eq!((three.runs()[0].ref_start, three.runs()[0].ref_end), (5, 5));
        assert!(
            matches!(three.runs()[0].label, Label::Dup),
            "stays a dup at 3'"
        );
        // 5' walks it one whole unit to interbase 3, still a valid dup.
        let five = place_only(cut, &c, ShuffleDirection::FivePrime);
        let r = &five.runs()[0];
        assert_eq!(
            (r.ref_start, r.ref_end),
            (3, 3),
            "the dup steps 5' by one unit, not base-by-base",
        );
        assert!(matches!(r.label, Label::Dup), "and stays a dup at 5'");
        assert_eq!(
            &c.reference[r.ref_start - r.alt.len()..r.ref_start],
            r.alt.as_slice(),
            "the 5' slot's payload is still the reference unit immediately before it",
        );
    }

    /// F2b sibling clamp (regression guard): `place_only` steps a locked `Dup`
    /// 5' by whole units (reading only the reference) and, before the clamp,
    /// could walk it onto a source a SIBLING run overwrote. The s02 5' arm was
    /// its natural witness, but the `pure-tandem-expansion-insertion-is-not-
    /// repartitioned` guard now stops that dup ever being minted, so this pins
    /// the fix by direct construction.
    ///
    /// The s02 3' ruled cut is `[8_9 T>A (sub); 11_11 AT dup]` over
    /// `TAAAATTATATTTATTATTT`. Placed 5', the unclamped walk stepped the `AT` dup
    /// from interbase 11 to 9, whose source `[7,9)` lies under the substitution at
    /// `[8,9)`, so `dup_source_overlaps_prior_piece` held — the render seam would
    /// refuse the dup and fall back to the per-member pipeline (the s05 leak
    /// mechanism). With the sibling clamp the dup stays at interbase 11 (source
    /// `[9,11)`, disjoint from the sub), so the overlap predicate is false and the
    /// dup renders.
    #[test]
    fn place_only_clamps_a_locked_dup_off_a_sibling_clobbered_source() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        let reference = b"TAAAATTATATTTATTATTT";
        let c = ctx(reference, b"TAAAATTAAATATTTATTATTT", &frame, &prov);
        let cut = Cut::new(
            &c,
            vec![
                Run::unlabelled(8, 9, b"A".to_vec()),
                Run {
                    ref_start: 11,
                    ref_end: 11,
                    alt: b"AT".to_vec(),
                    label: Label::Dup,
                    locked: true,
                },
            ],
        )
        .expect("sound s02 3' cut");
        let placed = place_only(cut, &c, ShuffleDirection::FivePrime);
        assert!(
            !merge::dup_source_overlaps_prior_piece(&cut_to_pieces(&placed), reference),
            "the sibling clamp must keep the placed dup off a source a sibling overwrote, so the \
             render seam accepts it (F2b fix)",
        );
        // Concretely: the dup stays at its 3'-most slot (interbase 11), one whole
        // unit shy of the 5' walk that would have collided with the sub at [8,9).
        let dup = placed
            .runs()
            .iter()
            .find(|r| matches!(r.label, Label::Dup))
            .expect("the dup survives placement");
        assert_eq!(
            (dup.ref_start, dup.ref_end),
            (11, 11),
            "the clamp stops the dup one unit shy of the sibling-clobbered source",
        );
    }

    #[test]
    fn an_unlabelled_cut_round_trips_through_the_bridge_unchanged() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // reference ACGT -> TCTT: sub A>T at [0,1), sub G>T at [2,3), one unchanged
        // base (C at 1) between them. splice: T + C + T + T = TCTT.
        let c = ctx(b"ACGT", b"TCTT", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(2, 3, b"T".to_vec()),
            ],
        )
        .expect("sound seed");

        let pieces = cut_to_pieces(&seed);
        assert_eq!(pieces.len(), 2);
        assert_eq!(
            (
                pieces[0].ref_start,
                pieces[0].ref_end,
                pieces[0].alt.as_slice()
            ),
            (0, 1, b"T".as_slice())
        );
        assert_eq!(
            (
                pieces[1].ref_start,
                pieces[1].ref_end,
                pieces[1].alt.as_slice()
            ),
            (2, 3, b"T".as_slice())
        );

        // An unlabelled cut is exactly what the bridge preserves: geometry and alt
        // bytes survive, and there were no labels/locks to lose, so it is identity.
        let back = pieces_to_cut(&c, pieces).expect("sound round-trip");
        assert_eq!(back, seed);
    }

    #[test]
    fn bridge_loses_labels_and_locks_by_construction() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // reference ACGT -> ACGTACGT: a tandem dup of [0,4). Label the run `Dup` and
        // lock it; a Piece has nowhere to carry either, so the bridge must strip both.
        let c = ctx(b"ACGT", b"ACGTACGT", &frame, &prov);
        let locked_dup = Run {
            ref_start: 4,
            ref_end: 4,
            alt: b"ACGT".to_vec(),
            label: Label::Dup,
            locked: true,
        };
        let seed = Cut::new(&c, vec![locked_dup]).expect("sound seed");
        assert!(
            seed.runs()[0].locked,
            "precondition: the seed run is locked"
        );

        let back = pieces_to_cut(&c, cut_to_pieces(&seed)).expect("sound round-trip");
        // The geometry and bytes survive; the label and lock do NOT. This is the
        // guardrail: nothing that flows through the bridge can come back locked.
        assert_eq!(back.runs().len(), 1);
        assert_eq!((back.runs()[0].ref_start, back.runs()[0].ref_end), (4, 4));
        assert_eq!(back.runs()[0].alt, b"ACGT");
        assert_eq!(back.runs()[0].label, Label::Unlabelled);
        assert!(!back.runs()[0].locked, "the bridge cannot restore a lock");
        assert_ne!(
            back, seed,
            "label loss makes this deliberately non-identity"
        );
    }

    #[test]
    fn lift_relabel_carries_an_incoming_label_forward_and_types_new_geometry() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // Incoming cut: a locked `Inv` over [0,6) (AAGCTA -> TAGCTT). A coarsening that
        // left this run untouched hands `lift_relabel` the same geometry back; with the
        // `preserve_only` typer the label and lock must survive — this is what lets a
        // `Dup`/`Inv` lock reach `SepZero` across a later adapter's cut_to_pieces round
        // trip.
        let c = ctx(b"AAGCTA", b"TAGCTT", &frame, &prov);
        let incoming = Cut::new(
            &c,
            vec![Run {
                ref_start: 0,
                ref_end: 6,
                alt: b"TAGCTT".to_vec(),
                label: Label::Inv,
                locked: true,
            }],
        )
        .expect("sound");
        let same = cut_to_pieces(&incoming);
        let back = lift_relabel(&c, &incoming, same, preserve_only).expect("sound");
        assert_eq!(
            back, incoming,
            "unchanged geometry keeps its label and lock"
        );

        // A DIFFERENT geometry over the same block (two subs) has no incoming match,
        // so `preserve_only` mints it `Unlabelled`/unlocked — a freshly-produced run.
        let fresh = vec![
            Piece {
                ref_start: 0,
                ref_end: 1,
                alt: b"T".to_vec(),
            },
            Piece {
                ref_start: 5,
                ref_end: 6,
                alt: b"T".to_vec(),
            },
        ];
        let minted = lift_relabel(&c, &incoming, fresh, preserve_only).expect("sound");
        assert!(
            minted
                .runs()
                .iter()
                .all(|r| r.label == Label::Unlabelled && !r.locked),
            "unmatched pieces are minted unlabelled/unlocked"
        );
    }

    #[test]
    fn lift_coarsen_declines_when_it_would_absorb_a_locked_run() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // Incoming carries a locked `Inv` over [0,6) (AAGCTA -> TAGCTT). A coarsening
        // whose output no longer contains that exact geometry has absorbed the lock;
        // `lift_coarsen` must decline so the inv survives the fixed point (the
        // `delins_hiding_an_inversion` reversion this guard exists to stop).
        let c = ctx(b"AAGCTA", b"TAGCTT", &frame, &prov);
        let incoming = Cut::new(
            &c,
            vec![Run {
                ref_start: 0,
                ref_end: 6,
                alt: b"TAGCTT".to_vec(),
                label: Label::Inv,
                locked: true,
            }],
        )
        .expect("sound");

        // An output that no longer carries the [0,6) inv run (same bases, re-cut into
        // two pieces) — the lock was absorbed.
        let absorbing = vec![
            Piece {
                ref_start: 0,
                ref_end: 3,
                alt: b"TAG".to_vec(),
            },
            Piece {
                ref_start: 3,
                ref_end: 6,
                alt: b"CTT".to_vec(),
            },
        ];
        assert!(
            lift_coarsen(&c, &incoming, absorbing).is_none(),
            "a coarsening may not absorb a locked run",
        );

        // Control: an output that preserves the locked run's exact geometry lifts and
        // keeps the label and lock.
        let preserving = cut_to_pieces(&incoming);
        let out = lift_coarsen(&c, &incoming, preserving).expect("the lock is preserved");
        assert_eq!(out, incoming, "the locked inv survives untouched");
    }

    #[test]
    fn lift_recut_declines_to_recut_a_locked_run() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // A Recut re-running over its OWN locked output must decline, or it re-cuts and
        // re-locks endlessly (Φ stalls). Model it: incoming carries a locked `Dup`
        // (zero-width "AC" at interbase 2 on AC -> ACAC), and the pass "re-cut" it into
        // a spanning geometry that no longer contains the [2,2)dup — a lock absorbed.
        let c = ctx(b"AC", b"ACAC", &frame, &prov);
        let incoming = Cut::new(
            &c,
            vec![Run {
                ref_start: 2,
                ref_end: 2,
                alt: b"AC".to_vec(),
                label: Label::Dup,
                locked: true,
            }],
        )
        .expect("sound");
        let recut = vec![Piece {
            ref_start: 0,
            ref_end: 2,
            alt: b"ACAC".to_vec(),
        }];
        assert!(
            lift_recut(&c, &incoming, recut, |_, _| None).is_none(),
            "a Recut may not re-cut a locked run",
        );

        // Control: a FIRST application over UNLOCKED input, producing NEW geometry,
        // locks that new output as a Recut must (Φ decrease). ACGT -> TCTT seeded as two
        // subs; the "recut" merges them into one new [0,3) delins.
        let c2 = ctx(b"ACGT", b"TCTT", &frame, &prov);
        let unlocked = Cut::new(
            &c2,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(2, 3, b"T".to_vec()),
            ],
        )
        .expect("sound");
        let merged = vec![Piece {
            ref_start: 0,
            ref_end: 3,
            alt: b"TCT".to_vec(),
        }];
        let out = lift_recut(&c2, &unlocked, merged, |_, _| None).expect("unlocked input lifts");
        assert_eq!(out.runs().len(), 1);
        assert!(out.runs()[0].locked, "new geometry from a Recut is locked");
    }

    #[test]
    fn lift_folds_a_coincident_insertion_pair_now_that_the_seam_is_retired() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // Two content-bearing insertions at interbase 2 (ACGT -> ACACGT). Before
        // step 2 this was the lift-refusal seam — `Cut::new` rejected it as
        // Unrenderable and `lift` counted a refusal. Now `Cut::new` FOLDS coincident
        // insertions, so `lift` succeeds with the one `Ins` they denote and records
        // no refusal.
        let c = ctx(b"ACGT", b"ACACGT", &frame, &prov);
        let before = lift_refusals();
        let colliding = vec![
            Piece {
                ref_start: 2,
                ref_end: 2,
                alt: b"A".to_vec(),
            },
            Piece {
                ref_start: 2,
                ref_end: 2,
                alt: b"C".to_vec(),
            },
        ];
        let cut = lift(&c, colliding).expect("the coincident pair now folds");
        assert_eq!(
            cut_to_pieces(&cut),
            vec![Piece {
                ref_start: 2,
                ref_end: 2,
                alt: b"AC".to_vec(),
            }]
        );
        assert_eq!(
            lift_refusals(),
            before,
            "no refusal: the fold retired the seam"
        );
    }

    #[test]
    fn lift_still_records_a_refusal_for_a_genuinely_unsound_intermediate() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // A piece list that does not splice back to `resulting` is Unsound, not a
        // renderability question the fold can repair — `Cut::new` still rejects it,
        // so `lift` returns None and counts the refusal. This is the seam's residual
        // purpose after coincident insertions are folded.
        let c = ctx(b"ACGT", b"AGGT", &frame, &prov);
        let before = lift_refusals();
        // Claims ref[1..2] (`C`) -> `T`, which splices to `ATGT`, not `AGGT`.
        let unsound = vec![Piece {
            ref_start: 1,
            ref_end: 2,
            alt: b"T".to_vec(),
        }];
        assert!(lift(&c, unsound).is_none());
        assert_eq!(lift_refusals(), before + 1, "the refusal is recorded");
    }

    #[test]
    fn shift_agrees_with_shift_pieces_in_both_directions() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // GAAAAT -> GAAAT: one A deleted from the run [1,5). Seeded 3'-most at [4,5);
        // the 5' shift walks it to [1,2), the 3' shift leaves it.
        let c = ctx(b"GAAAAT", b"GAAAT", &frame, &prov);
        let seed = Cut::new(&c, vec![Run::unlabelled(4, 5, Vec::new())]).expect("sound");
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            let mut raw = cut_to_pieces(&seed);
            merge::shift_pieces(&mut raw, c.reference, direction);
            let out = shift(seed.clone(), &c, direction);
            assert_eq!(cut_to_pieces(&out), raw, "{direction:?}");
        }
        let five = shift(seed.clone(), &c, ShuffleDirection::FivePrime);
        assert_eq!((five.runs()[0].ref_start, five.runs()[0].ref_end), (1, 2));
        let three = shift(seed, &c, ShuffleDirection::ThreePrime);
        assert_eq!((three.runs()[0].ref_start, three.runs()[0].ref_end), (4, 5));
    }
}
