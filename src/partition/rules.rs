//! Concrete ruled-cut rules — the Stage-B `merge.rs` passes migrated one at a time
//! into [`Rule`](crate::partition::ruled::Rule) impls (design §10 step 1). Each rule
//! reproduces the semantics of the `merge.rs` function it replaces; the migration is
//! byte-identity-gated at each step. Additive: nothing here is wired into the shipping
//! path until `PartitionRule::Ruled` is turned on.

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
/// **Lock exclusion (design §4.2, step 2).** A separation-zero merge does not absorb
/// a *locked* run — the two labels that lock are `Dup` and `Inv`, the ones that
/// outrank a merge. This is the data invariant that replaces the shipping
/// "do not re-run `coalesce_adjacent_pieces` after the peel" call-order hack
/// (`merge.rs:11766`): once a peeled `Dup` is locked, no order of rule application
/// can coarsen it away. It is byte-identical in the single-pass driver — `SepZero`
/// runs before any `Label` rule sets a lock, so no locked run is ever present when
/// it fires — and becomes load-bearing at the step-5 fixed point, where the mirror
/// could otherwise re-run `SepZero` after the peel.
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
            // UNLOCKED runs merge as before (byte-identical to the single-pass
            // driver, where nothing is locked when `SepZero` runs).
            let l = &runs[i - 1];
            let r = &runs[i];
            let both_consume_ref = l.ref_end > l.ref_start && r.ref_end > r.ref_start;
            let lock_ok = (!l.locked && !r.locked) || both_consume_ref;
            if l.ref_end == r.ref_start && lock_ok {
                let next = runs.remove(i);
                let prev = &mut runs[i - 1];
                prev.ref_end = next.ref_end;
                prev.alt.extend(next.alt);
                // The merged geometry is a new, untyped run.
                prev.label = Label::Unlabelled;
                prev.locked = false;
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
