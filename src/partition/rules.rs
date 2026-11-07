//! Concrete ruled-cut rules — the Stage-B `merge.rs` passes migrated one at a time
//! into [`Rule`](crate::partition::ruled::Rule) impls (design §10 step 1). Each rule
//! reproduces the semantics of the `merge.rs` function it replaces; the migration is
//! byte-identity-gated at each step. These rules are on the shipping path:
//! `PartitionRule::Ruled` is the default rule, and [`crate::partition::driver`] runs
//! them on every block it partitions.

use crate::partition::block_ctx::BlockCtx;
use crate::partition::ruled::{Authority, Cut, Label, Rule, RuleKind, Scope};

/// `SepZero` — merge two runs with no unchanged reference base between them into one
/// (separation zero). Ports `merge.rs::coalesce_adjacent_pieces`: at separation zero
/// the split spelling is `class="invalid"` (`substitution.md:32`), so two adjacent
/// members are one `delins`.
///
/// Coarsen, [`Scope::ANY`]. Whole-merge in one `apply`, matching the original pass's
/// single-call semantics (the step-1 engine is single-pass). The merged run is
/// re-`Unlabelled` — the new geometry has not been typed.
///
/// **Lock guard (design §4.2, scoped by F4).** Two unlocked adjacent runs always
/// merge. A *locked* run is absorbed only when both runs consume reference bases,
/// which is the reach of the ruling this rule cites
/// (`delins-adjacent-members-when-both-consume-reference`):
///
/// * a locked `Dup` is a zero-width insertion, so it consumes no reference and is
///   **never** absorbed. Once a peeled `Dup` is locked, no order of rule application
///   can coarsen it away (#2175, the peeled-dup-beside-a-change invariant). This is
///   the data invariant that replaces the shipping "do not re-run
///   `coalesce_adjacent_pieces` after the peel" call-order hack (see
///   `merge_coincident_insertions` in `src/normalize/merge.rs`);
/// * a locked `Inv`, or a Recut's locked `Unlabelled` output, **is** absorbed when
///   the run beside it also consumes reference (a `sub`/`del`/`delins`/`inv`).
///
/// The merged run is locked only when every run it absorbed was locked. Unlocking
/// it would raise `Φ = (unlocked, total)` when both inputs are locked — `(0, 2)`
/// becomes `(1, 1)`, which trips `run_engine`'s progress assert — and would let a
/// Recut re-split the run SepZero just merged. A merge with any unlocked input stays
/// unlocked and open to later rules.
///
/// The lock guard matters at the step-5 fixed point, where `SepZero` re-fires after
/// the Recuts have set locks. Before it, in the pre-shrink application, the only
/// lock present is a whole-span `Inv` the seed emits as a single run.
pub struct SepZero;

impl Rule for SepZero {
    fn id(&self) -> &'static str {
        "sep-zero"
    }

    fn kind(&self) -> RuleKind {
        RuleKind::Coarsen
    }

    fn scope(&self) -> Scope {
        Scope::ANY
    }

    fn authority(&self) -> Authority {
        Authority::Ruling("delins-adjacent-members-when-both-consume-reference")
    }

    fn apply(&self, cut: &Cut, ctx: &BlockCtx) -> Option<Cut> {
        let mut runs = cut.runs().to_vec();
        let mut merged = false;
        let mut i = 1;
        while i < runs.len() {
            debug_assert!(
                runs[i - 1].ref_end <= runs[i].ref_start,
                "strict overlap between adjacent runs: the 3'-shuffle cannot produce this",
            );
            // Lock guard, scoped by the ruling this rule cites (F4). A sep-zero
            // merge absorbs a locked run only when
            // `delins-adjacent-members-when-both-consume-reference` reaches it:
            // both members consume reference bases (a locked `Inv` beside a
            // sub/del/delins/inv). A locked `Dup` is a zero-width insertion — it
            // consumes no reference, so `both_consume_ref` is false and it is never
            // absorbed (#2175, the peeled-dup-beside-a-change invariant). Two
            // UNLOCKED runs always merge.
            let l = &runs[i - 1];
            let r = &runs[i];
            let both_consume_ref = l.ref_end > l.ref_start && r.ref_end > r.ref_start;
            let lock_ok = (!l.locked && !r.locked) || both_consume_ref;
            if l.ref_end == r.ref_start && lock_ok {
                let next = runs.remove(i);
                let prev = &mut runs[i - 1];
                prev.ref_end = next.ref_end;
                prev.alt.extend(next.alt);
                // The merged geometry is a new, untyped run. It stays locked only if
                // both sides were locked, so `Φ` strictly decreases either way
                // (see the type doc).
                prev.label = Label::Unlabelled;
                prev.locked = prev.locked && next.locked;
                merged = true;
            } else {
                i += 1;
            }
        }
        if !merged {
            return None;
        }
        // A coarsening preserves the splice, so the result round-trips whenever the
        // input did; `Cut::new` re-checks it (and renderability). A defensive `.ok()`
        // rather than `.expect()` — a rule declines rather than panics.
        Cut::new(ctx, runs).ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::block_ctx::{FrameContext, Molecule, Provenance};
    use crate::partition::ruled::Run;

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
    fn two_adjacent_runs_merge_at_separation_zero() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // reference ACGT; run0 [0,1)->T (A>T), run1 [1,2)->T (C>T); zero separation.
        // splice: T + T + reference[2..4]=GT = TTGT.
        let c = ctx(b"ACGT", b"TTGT", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(1, 2, b"T".to_vec()),
            ],
        )
        .expect("sound seed");
        let out = SepZero.apply(&seed, &c).expect("merged");
        assert_eq!(out.runs().len(), 1);
        assert_eq!((out.runs()[0].ref_start, out.runs()[0].ref_end), (0, 2));
        assert_eq!(out.runs()[0].alt, b"TT");
        // total dropped 2 -> 1, so Phi strictly decreased (the engine's obligation).
        assert!(out.phi() < seed.phi());
    }

    #[test]
    fn a_gap_between_runs_is_not_merged() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // reference ACGT; run0 [0,1)->T (A>T), run1 [2,3)->T (G>T); one unchanged base
        // (C at 1) between them. splice: T + C + T + T = TCTT.
        let c = ctx(b"ACGT", b"TCTT", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(2, 3, b"T".to_vec()),
            ],
        )
        .expect("sound seed");
        // Nothing is adjacent, so the rule declines.
        assert!(SepZero.apply(&seed, &c).is_none());
    }

    #[test]
    fn sep_zero_does_not_absorb_a_locked_run() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // reference ACGT -> ACCTT: a locked `Dup` (zero-width insert "C" at interbase
        // 2) immediately followed by a Sub G>T at [2,3). They are adjacent at ref
        // offset 2, so an unlocked SepZero merges them; the lock forbids it. This is
        // the peeled-dup-beside-a-change invariant (#2175) as data, not call order.
        let c = ctx(b"ACGT", b"ACCTT", &frame, &prov);
        let locked_dup = Run {
            ref_start: 2,
            ref_end: 2,
            alt: b"C".to_vec(),
            label: Label::Dup,
            locked: true,
        };
        let sub = Run::unlabelled(2, 3, b"T".to_vec());
        let seed = Cut::new(&c, vec![locked_dup.clone(), sub.clone()]).expect("sound");
        assert!(
            SepZero.apply(&seed, &c).is_none(),
            "a locked Dup is never absorbed by a sep-zero merge"
        );

        // Control: the SAME geometry with the dup UNLOCKED merges into one run — so
        // the decline above is the lock's doing, not the geometry's.
        let unlocked_dup = Run {
            label: Label::Unlabelled,
            locked: false,
            ..locked_dup
        };
        let seed2 = Cut::new(&c, vec![unlocked_dup, sub]).expect("sound");
        let merged = SepZero.apply(&seed2, &c).expect("the unlocked pair merges");
        assert_eq!(merged.runs().len(), 1);
        assert_eq!(
            (merged.runs()[0].ref_start, merged.runs()[0].ref_end),
            (2, 3)
        );
        assert_eq!(merged.runs()[0].alt, b"CT");
    }

    #[test]
    fn sep_zero_never_raises_phi_when_both_adjacent_runs_are_locked() {
        use crate::partition::ruled::run_engine;
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // reference ACGT -> GTAC: two adjacent locked `Inv` runs, [0,2) AC->GT and
        // [2,4) GT->AC. Both consume reference, so
        // `delins-adjacent-members-when-both-consume-reference` merges them. With both
        // inputs locked, Phi = (0 unlocked, 2 total); an UNLOCKED merged run would be
        // (1, 1), which is lexicographically HIGHER and trips `run_engine`'s progress
        // assert. The merge must keep Phi strictly decreasing.
        let c = ctx(b"ACGT", b"GTAC", &frame, &prov);
        let inv = |ref_start, ref_end, alt: &[u8]| Run {
            ref_start,
            ref_end,
            alt: alt.to_vec(),
            label: Label::Inv,
            locked: true,
        };
        let seed = Cut::new(&c, vec![inv(0, 2, b"GT"), inv(2, 4, b"AC")]).expect("sound");
        assert_eq!(seed.phi(), (0, 2));
        let out = SepZero
            .apply(&seed, &c)
            .expect("both consume reference, so they merge");
        assert!(
            out.phi() < seed.phi(),
            "Phi must strictly decrease: {:?} -> {:?}",
            seed.phi(),
            out.phi()
        );
        assert_eq!(out.runs().len(), 1);
        assert_eq!((out.runs()[0].ref_start, out.runs()[0].ref_end), (0, 4));
        assert_eq!(out.runs()[0].alt, b"GTAC");
        assert!(
            out.runs()[0].locked,
            "a merge of two locked runs stays locked"
        );

        // Through the engine's restart-from-top loop, which debug-asserts progress.
        let rules: [&dyn Rule; 1] = [&SepZero];
        let normal = run_engine(seed, &c, &rules);
        assert_eq!(normal.runs().len(), 1);
    }

    #[test]
    fn sep_zero_merge_of_a_locked_and_an_unlocked_run_is_unlocked() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // reference ACGT -> GTTT: a locked `Inv` [0,2) AC->GT beside an unlocked sub
        // [2,3) G->T. Phi (1, 2) -> (1, 1) with an unlocked result, so the merged
        // geometry stays open to later rules exactly as before.
        let c = ctx(b"ACGT", b"GTTT", &frame, &prov);
        let locked_inv = Run {
            ref_start: 0,
            ref_end: 2,
            alt: b"GT".to_vec(),
            label: Label::Inv,
            locked: true,
        };
        let seed =
            Cut::new(&c, vec![locked_inv, Run::unlabelled(2, 3, b"T".to_vec())]).expect("sound");
        let out = SepZero.apply(&seed, &c).expect("merged");
        assert!(out.phi() < seed.phi());
        assert_eq!(out.runs().len(), 1);
        assert!(!out.runs()[0].locked);
        assert_eq!(out.runs()[0].label, Label::Unlabelled);
    }

    #[test]
    fn three_adjacent_runs_collapse_to_one() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // reference ACGT -> TTTT via three adjacent single-base subs, zero separation.
        let c = ctx(b"ACGT", b"TTTT", &frame, &prov);
        let seed = Cut::new(
            &c,
            vec![
                Run::unlabelled(0, 1, b"T".to_vec()),
                Run::unlabelled(1, 2, b"T".to_vec()),
                Run::unlabelled(2, 3, b"T".to_vec()),
            ],
        )
        .expect("sound seed");
        let out = SepZero.apply(&seed, &c).expect("merged");
        assert_eq!(out.runs().len(), 1);
        assert_eq!((out.runs()[0].ref_start, out.runs()[0].ref_end), (0, 3));
        assert_eq!(out.runs()[0].alt, b"TTT");
    }
}
