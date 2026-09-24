//! The driver of the ruled cut over `Seed::canonical`, in the shipping pipeline's
//! order (design §10).
//!
//! [`partition_ruled`] runs three stages:
//!
//! 1. the pre-placement rules, once each, in the shipping order;
//! 2. the 3'-most shift, the pre-shrink `SepZero` merge and the shrink, then the
//!    post-shrink rules ([`fixpoint_rules`]) run to a NORMAL FORM under
//!    [`run_engine`]'s restart-from-top (step 5) — a rule re-fires after another
//!    rule changed the cut, which a single pass could not reach;
//! 3. final placement (step 7): the direction mirror is retired, so the partition
//!    is derived once on the 3' placement and a 5' request re-places it with
//!    `bridge::place_only`, which preserves labels and member count.
//!
//! This module reaches NO `merge` internal: the seed, the bridge and the rules do. Order
//! is data (`pipeline_order`), pinned by a test against the `RULE_MIGRATION` inventory so
//! a rule cannot be registered and not driven.

use crate::normalize::ShuffleDirection;
use crate::partition::adapters::{
    CodonException, CodonFrameSeparation, CompensatingGaps, PayloadCoincidence, PlacedGap, RunInv,
    SplitConcealed, TandemDupRun,
};
use crate::partition::block_ctx::BlockCtx;
use crate::partition::bridge;
use crate::partition::partitioner::PartitionError;
use crate::partition::ruled::{run_engine, Cut, Rule};
use crate::partition::rules::SepZero;
use crate::partition::seed::{self, Seed};

/// Rules applied before placement, in the shipping order: the placed-gap
/// collapse (inside `partition_block_for_rule`), the early codon merge, the
/// concealed-separation audit.
static PRE_PLACEMENT_RULES: [&(dyn Rule + Sync); 3] =
    [&PlacedGap, &CodonFrameSeparation, &SplitConcealed];

/// Rules applied inside the direction mirror AFTER `shift` -> `SepZero` ->
/// shrink, in the shipping order.
static POST_SHRINK_RULES: [&(dyn Rule + Sync); 5] = [
    &PayloadCoincidence,
    &CompensatingGaps,
    &TandemDupRun,
    &RunInv,
    &CodonException,
];

/// The rules the post-shrink fixed point runs, in registration (priority) order:
/// `SepZero` leads so it re-fires (its lock-guard protects a `Dup`/`Inv`) whenever
/// a coarsen or recut exposes a fresh sep-0 adjacency, then `POST_SHRINK_RULES`
/// verbatim, indexed so the fixpoint set cannot drift from [`pipeline_order`].
///
/// Public so an out-of-crate instrument (the critical-pair census) sweeps exactly
/// the system [`partition_ruled`] runs, rather than a copy that can go stale.
pub fn fixpoint_rules() -> [&'static dyn Rule; 6] {
    [
        &SepZero,
        POST_SHRINK_RULES[0],
        POST_SHRINK_RULES[1],
        POST_SHRINK_RULES[2],
        POST_SHRINK_RULES[3],
        POST_SHRINK_RULES[4],
    ]
}

/// The nine rule ids in execution order — the pipeline's order, stated once.
pub fn pipeline_order() -> [&'static str; 9] {
    let mut out = [""; 9];
    let mut i = 0;
    for rule in PRE_PLACEMENT_RULES {
        out[i] = rule.id();
        i += 1;
    }
    out[i] = SepZero.id();
    i += 1;
    for rule in POST_SHRINK_RULES {
        out[i] = rule.id();
        i += 1;
    }
    out
}

/// Apply `rule` once if its scope admits `ctx`; a declining rule leaves the cut
/// as it was. Single pass, no fixed point, no Φ assertion (Recuts may raise it).
fn apply_once(rule: &dyn Rule, cut: Cut, ctx: &BlockCtx) -> Cut {
    if !rule.scope().admits(ctx) {
        return cut;
    }
    rule.apply(&cut, ctx).unwrap_or(cut)
}

/// Partition `ctx` (a canonical window and its resulting sequence) with the
/// ruled cut: the pre-placement rules once each, the post-shrink rules to a fixed
/// point, then final placement in `direction` (see the module doc).
///
/// `Err(GridTooLarge)` when the seed declines; `Err(Unrenderable)` when some
/// intermediate could not be represented as a `Cut` (the bridge's lift-refusal
/// seam). Either way the caller takes the shipping chain for the block.
pub fn partition_ruled(ctx: &BlockCtx, direction: ShuffleDirection) -> Result<Cut, PartitionError> {
    let refusals_before = bridge::lift_refusals();
    let mut cut = Seed::canonical(ctx)?;
    for rule in PRE_PLACEMENT_RULES {
        cut = apply_once(rule, cut, ctx);
    }
    // Placement out (design §10 step 7): the direction mirror is retired. The rules
    // derive THE canonical partition once, on the seed's own 3'-most placement (the
    // 3' direction), and placement is then a final, per-run, label- and
    // member-count-preserving re-placement in the requested direction. So the 5'
    // arm is the 3' partition re-placed 5', and the member count is
    // direction-independent by construction — not by running everything twice and
    // adopting the fewer-member direction.
    let cut = {
        let cut = bridge::shift(cut, ctx, ShuffleDirection::ThreePrime);
        // The pre-shrink sep-0 merge, kept exactly where the single-pass driver had
        // it (the byte-identical entry into the fixed point).
        let cut = apply_once(&SepZero, cut, ctx);
        let cut = seed::shrink_to_differences(cut, ctx);
        // Fixed point (step 5): run the post-shrink rules to a normal form under
        // restart-from-top; see `fixpoint_rules` for the set and its order.
        run_engine(cut, ctx, &fixpoint_rules())
    };
    // Final placement. The 3' arm is already at its 3'-most placement; the 5' arm is
    // the same partition re-placed 5' (label- and count-preserving).
    let cut = match direction {
        ShuffleDirection::ThreePrime => cut,
        ShuffleDirection::FivePrime => bridge::place_only(cut, ctx, ShuffleDirection::FivePrime),
    };
    if bridge::lift_refusals() != refusals_before {
        return Err(PartitionError::Unrenderable);
    }
    Ok(cut)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normalize::merge;
    use crate::partition::adapters::RULE_MIGRATION;
    use crate::partition::block_ctx::{FrameContext, Molecule, Provenance};
    use crate::partition::bridge::{cut_to_pieces, place_only};

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
    fn the_pipeline_order_is_the_shipping_order_and_drives_every_migrated_rule() {
        assert_eq!(
            pipeline_order(),
            [
                "placed-gap",
                "codon-frame-separation",
                "split-concealed",
                "sep-zero",
                "payload-coincidence",
                "compensating-gaps",
                "tandem-dup-run",
                "run-inv",
                "codon-exception",
            ]
        );
        let mut driven: Vec<&str> = pipeline_order().to_vec();
        driven.sort_unstable();
        let mut inventory: Vec<&str> = RULE_MIGRATION.iter().map(|(id, _)| *id).collect();
        inventory.sort_unstable();
        assert_eq!(
            driven, inventory,
            "every RULE_MIGRATION row is driven exactly once"
        );
    }

    #[test]
    fn a_declined_seed_propagates_grid_too_large() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        let reference = vec![b'A'; 4200];
        let resulting = vec![b'C'; 4200];
        let c = ctx(&reference, &resulting, &frame, &prov);
        assert!(matches!(
            partition_ruled(&c, ShuffleDirection::ThreePrime),
            Err(PartitionError::GridTooLarge {
                ref_len: 4200,
                alt_len: 4200
            })
        ));
    }

    /// The driver reproduces the shipping Stage-B chain on a NonCoding DNA
    /// window. The replica below IS the chain's body
    /// (`canonicalize_from_sequence_with_rule`, re-grep before trusting its
    /// order), transcribed over pieces; the definitive gate is the wiring-level
    /// comparison in merge.rs, this pins the driver in isolation.
    fn shipping_chain(
        reference: &[u8],
        resulting: &[u8],
        direction: ShuffleDirection,
    ) -> Vec<merge::Piece> {
        let (lo, hi_ref, hi_alt) = merge::trim_common_flanks(reference, resulting);
        let block_ref = &reference[lo..hi_ref];
        let block_alt = &resulting[lo..hi_alt];
        let mut pieces = merge::partition_block_canonical(block_ref, block_alt).expect("DAG");
        if merge::split_is_a_placed_gap_coincidence(
            &pieces,
            block_ref,
            block_alt,
            merge::CoincidenceCarveOut::InReach,
        ) {
            pieces = vec![merge::Piece {
                ref_start: 0,
                ref_end: block_ref.len(),
                alt: block_alt.to_vec(),
            }];
        }
        for p in &mut pieces {
            p.ref_start += lo;
            p.ref_end += lo;
        }
        let length_changing = hi_ref != hi_alt;
        merge::coalesce_coding_frame_separation(&mut pieces, false, length_changing, 0, reference);
        merge::split_concealed_separations(&mut pieces, false, length_changing, 0, reference);
        merge::place_direction_symmetrically(&mut pieces, direction, |pieces, dir| {
            merge::shift_pieces(pieces, reference, dir);
            merge::coalesce_adjacent_pieces(pieces);
            merge::shrink_pieces_to_differences(pieces, reference);
            merge::coalesce_payload_alignment_split(pieces, reference);
            merge::coalesce_compensating_gap_split(pieces, reference);
            merge::coalesce_by_run(pieces, reference, |run, reference| {
                merge::peel_tandem_dup_beside_change(run, reference, merge::PeelReach::TractOnly);
                merge::coalesce_solid_run(run, reference);
            });
            merge::coalesce_inversion_runs(pieces, reference, lo, block_ref, block_alt);
            merge::apply_coding_codon_exception(pieces, false, 0, reference, None);
        });
        pieces
    }

    #[test]
    fn the_driver_reproduces_the_shipping_chain() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        let cases: [(&[u8], &[u8]); 8] = [
            // #2161: `43_47inv` beside `50G>T` on the inverted-repeat contig. The
            // DAG seed is `[insG; del C; G>T]` and the inv is an exact window hull
            // only under the 5' placement of the del — the shipping mirror found it
            // by running 5'; the ruled arm must find it from the 3' placement alone.
            (
                b"AACCGGTTAATCGATCGATTGCACGTACGTGCAATCGATCGATTAACCGGTTAACCGGTTAACCGG",
                b"AACCGGTTAATCGATCGATTGCACGTACGTGCAATCGATCGAGTTAACGTTTAACCGGTTAACCGG",
            ),
            // #2175 tandem dup beside a change (peel + solid run)
            (
                b"TAGTAAACCATTTTACGGAGGATCACAAATTCCTCCTTAT",
                b"TAGTAAACCATTTTACGGAGGATCACACACATTCCTCCTTAT",
            ),
            // #1040 compensating split, padded so the seed has flanks to roll into
            (b"ACGTACGTCAGTGACTAGACGTACGT", b"ACGTACGTTGTCACGACTACGTACGT"),
            // a whole-span reverse complement of two substitutions, padded
            (b"GGGGAAGCTAGGGG", b"GGGGTAGCTTGGGG"),
            // net deletion: one A of the run lost and the C after it substituted
            // (the DAG seeds this as one delins; a length-changing control that
            // exercises SplitConcealed/PayloadCoincidence gating)
            (b"CCCGAAAACCCC", b"CCCGAAATCCC"),
            // The three critical-pair witnesses the step-6 census surfaced
            // (`partition_critical_pair_census`): whole-span reverse complements
            // where RunInv competes with a coincidence coarsen / a tandem-dup run.
            (b"GACA", b"TGTC"), // compensating-gaps X run-inv
            (b"CACA", b"GCGC"), // compensating-gaps X tandem-dup-run
            (b"GCA", b"TGC"),   // run-inv X tandem-dup-run
        ];
        let mut fired = 0usize;
        let mut diverged_from_shipping: Vec<String> = Vec::new();
        for (reference, resulting) in cases {
            let c = ctx(reference, resulting, &frame, &prov);
            let p3 = partition_ruled(&c, ShuffleDirection::ThreePrime).expect("3' answers");
            let p5 = partition_ruled(&c, ShuffleDirection::FivePrime).expect("5' answers");

            // Step-7 contract: the 5' arm IS the 3' partition re-placed 5', so the
            // member count is direction-independent by construction — the property
            // the retired mirror enforced by running everything twice.
            assert_eq!(
                cut_to_pieces(&p5),
                cut_to_pieces(&place_only(p3.clone(), &c, ShuffleDirection::FivePrime)),
                "5' arm is the 3' partition re-placed on {}",
                String::from_utf8_lossy(reference)
            );
            assert_eq!(
                p3.runs().len(),
                p5.runs().len(),
                "member count is direction-independent on {}",
                String::from_utf8_lossy(reference)
            );

            for (direction, got) in [
                (ShuffleDirection::ThreePrime, &p3),
                (ShuffleDirection::FivePrime, &p5),
            ] {
                let want = shipping_chain(reference, resulting, direction);
                if cut_to_pieces(got) != want {
                    diverged_from_shipping.push(format!(
                        "{direction:?} on {}",
                        String::from_utf8_lossy(reference)
                    ));
                }
                let seed = Seed::canonical(&c).expect("seed");
                if cut_to_pieces(&seed) != want {
                    fired += 1;
                }
            }
        }
        // Step 7 retires `place_direction_symmetrically` and replaces it with
        // `place_only`. The one case where that once diverged from the shipping
        // mirror — the #2175 tandem-dup block's 5' arm, where `place_only` walked
        // the inserted `CA` to its 5'-most slot while the mirror kept it 3'-most —
        // is closed by the F2b sibling clamp: a dup peeled BESIDE a change abuts a
        // 3' sibling, so `place_only` now freezes it at the 3'-most (abutting) slot
        // in both directions, exactly as the shipping chain places it
        // (`issue_2175_dup_abutting_change`, `fragmentation_corpus` pin the same
        // 3'-most form for 5' and 3'). So the ruled driver now reproduces the
        // shipping chain on the whole corpus.
        assert!(
            diverged_from_shipping.is_empty(),
            "ruled must reproduce the shipping chain on every case; diverged on {diverged_from_shipping:?}"
        );
        assert!(fired > 0, "non-vacuity: some rule fired on some case");
    }
}
