//! The seed of the ruled cut (design §4.3): the DAG canonical alignment of the
//! trimmed block, 3'-most tie-break included, offset back into the window and
//! shrunk so every run is its own difference hull.
//!
//! This is the third confined reacher into `normalize::merge` (after `registry`
//! and `bridge`): the seed IS `partition_block_canonical`, so it calls it.
//! `ruled.rs` stays clean of `merge` internals; only this module may hold the
//! shrink, because the shrink is the seed's invariant, not the engine's.
//!
//! # The context is the WINDOW (design §12 R1, stated)
//!
//! `ctx.reference`/`ctx.resulting` are the padded canonical window and its
//! resulting sequence — exactly what every Stage-B pass in
//! `canonicalize_from_sequence_with_rule` reads — and the runs are in window
//! coordinates. The seed trims to the changed block itself (the same
//! `trim_common_flanks` call the pipeline makes) and offsets the DAG's pieces by
//! `lo`, so a rule that reaches 5' of the block (the tandem-tract scan, the
//! peel's source, the 5' mirror) sees the pad the pipeline gives it.

use crate::normalize::merge;
use crate::partition::block_ctx::BlockCtx;
use crate::partition::bridge::{cut_to_pieces, lift, pieces_to_cut};
use crate::partition::partitioner::PartitionError;
use crate::partition::ruled::{Cut, Label, Run};

/// The ruled cut's starting point. Stateless: [`Seed::canonical`] builds the seed
/// [`Cut`] for a block context, and the ruled driver's rules rewrite from there.
pub struct Seed;

impl Seed {
    /// The canonical seed of `ctx`: trim the window to its changed block, cut
    /// the block with the DAG canonical alignment (fewest members, 3'-most —
    /// `merge::partition_block_canonical`), offset the pieces back into window
    /// coordinates, and shrink each to its difference hull.
    ///
    /// `Err(GridTooLarge)` when the DAG declines the block (grid-cell or span
    /// cap); the caller then takes the shipping per-member pipeline — the same
    /// exit the window cap already forces (design §1.1). An identical pair
    /// seeds the empty cut, which is sound by construction.
    pub fn canonical(ctx: &BlockCtx) -> Result<Cut, PartitionError> {
        let (lo, hi_ref, hi_alt) = merge::trim_common_flanks(ctx.reference, ctx.resulting);
        if lo == hi_ref && lo == hi_alt {
            return Cut::new(ctx, Vec::new());
        }
        let block_ref = &ctx.reference[lo..hi_ref];
        let block_alt = &ctx.resulting[lo..hi_alt];

        // F1 (`whole-span-reverse-complement-types-as-inv`): a span whose whole
        // content is its exact reverse complement is one `inv`, and the ruling
        // settles that typing *before any cut is considered*. Deciding it here —
        // as one locked `Inv` run, skipping the DAG — keeps a clean whole-span
        // inversion out of every later rule's reach (no coarsen or peel can
        // absorb a locked `Inv`), instead of letting the DAG's minimal alignment
        // shred it into indel-bearing members and relying on `RunInv` — registered
        // last in the fixpoint — to reassemble what the earlier rules have locked.
        // Equal-length and >= 2 bases by the helper's own gate, so `shift_pieces`
        // never moves it and the predicate is a pure function of the block.
        if let Some(pieces) = merge::whole_span_reverse_complement(block_ref, block_alt, lo) {
            debug_assert_eq!(
                pieces.len(),
                1,
                "a whole-span reverse complement is one run"
            );
            let piece = &pieces[0];
            let run = Run {
                ref_start: piece.ref_start,
                ref_end: piece.ref_end,
                alt: piece.alt.clone(),
                label: Label::Inv,
                locked: true,
            };
            return Cut::new(ctx, vec![run]);
        }

        let mut pieces = merge::partition_block_canonical(block_ref, block_alt).ok_or(
            PartitionError::GridTooLarge {
                ref_len: block_ref.len(),
                alt_len: block_alt.len(),
            },
        )?;
        for piece in &mut pieces {
            piece.ref_start += lo;
            piece.ref_end += lo;
        }
        // A run is always its own difference hull (design §4.3). A no-op on a
        // DAG seed — `the_seed_is_a_fixed_point_of_shrink` — but stated here
        // so the invariant is the seed's, not an accident of the aligner.
        merge::shrink_pieces_to_differences(&mut pieces, ctx.reference);
        pieces_to_cut(ctx, pieces)
    }
}

/// Shrink every run of `cut` to the hull of its own differences —
/// `merge::shrink_pieces_to_differences` over a `Cut`. Φ-neutral (the run count
/// is unchanged), so it is not a rule; the driver applies it where the shipping
/// chain does (after `SepZero`, before the payload passes). A shrink preserves
/// the splice, so lifting back cannot fail; the fallback to the input is
/// belt-and-braces and never a silent change of geometry.
pub(crate) fn shrink_to_differences(cut: Cut, ctx: &BlockCtx) -> Cut {
    let mut pieces = cut_to_pieces(&cut);
    merge::shrink_pieces_to_differences(&mut pieces, ctx.reference);
    lift(ctx, pieces).unwrap_or(cut)
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

    fn geometry(cut: &Cut) -> Vec<(usize, usize, Vec<u8>)> {
        cut.runs()
            .iter()
            .map(|r| (r.ref_start, r.ref_end, r.alt.clone()))
            .collect()
    }

    #[test]
    fn the_seed_is_the_dag_partition_of_the_trimmed_block_offset_into_the_window() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // Window GGACGTGG -> GGTTGG. trim_common_flanks eats the GG prefix and,
        // on the 3' side, the shared TGG suffix (ref[5]=T matches result[3]=T),
        // so the block is ref[2..5]="ACG" -> result[2..3]="T", i.e. (2, 5, 3).
        // The DAG cuts "ACG" -> "T" as one net-deletion delins [0,3)->"T", which
        // offset by lo is window [2,5)->"T".
        let c = ctx(b"GGACGTGG", b"GGTTGG", &frame, &prov);
        let seed = Seed::canonical(&c).expect("the DAG accepts a 3x1 block");

        // Definitional check: the raw DAG on the trimmed block, offset by lo.
        let (lo, hi_ref, hi_alt) = merge::trim_common_flanks(c.reference, c.resulting);
        assert_eq!((lo, hi_ref, hi_alt), (2, 5, 3));
        let raw =
            merge::partition_block_canonical(&c.reference[lo..hi_ref], &c.resulting[lo..hi_alt])
                .expect("raw DAG");
        let raw_offset: Vec<(usize, usize, Vec<u8>)> = raw
            .iter()
            .map(|p| (p.ref_start + lo, p.ref_end + lo, p.alt.clone()))
            .collect();
        assert_eq!(geometry(&seed), raw_offset, "seed != DAG offset by lo");
        // And the concrete geometry, so a DAG change is noticed here too.
        assert_eq!(geometry(&seed), vec![(2, 5, b"T".to_vec())]);
        assert!(seed.runs().iter().all(|r| !r.locked), "a seed is unlocked");
    }

    #[test]
    fn a_whole_span_reverse_complement_seeds_as_one_locked_inv() {
        // F1 (`whole-span-reverse-complement-types-as-inv`): a block whose whole
        // content is its exact reverse complement seeds as ONE locked `Inv` run,
        // not the DAG's indel-bearing minimal alignment. Sabotage: delete the F1
        // branch in `Seed::canonical` -> the DAG cuts this block into several runs
        // (or a delins), and both the geometry and the label assertions fail.
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // Window GGACGG -> GGGTGG. trim eats the GG flanks; block ref[2..4]="AC"
        // -> result[2..4]="GT" = revcomp("AC").
        let c = ctx(b"GGACGG", b"GGGTGG", &frame, &prov);
        let seed = Seed::canonical(&c).expect("seed");
        assert_eq!(geometry(&seed), vec![(2, 4, b"GT".to_vec())]);
        let run = &seed.runs()[0];
        assert!(
            matches!(run.label, Label::Inv),
            "a whole-span revcomp seeds as Inv"
        );
        assert!(run.locked, "the whole-span inv is locked");
        // A 16-base coding-frame inversion (the c.612_627inv shape) likewise.
        let c2 = ctx(
            b"GGATATCAATTAGGCATAGG",
            b"GGTATGCCTAATTGATATGG",
            &frame,
            &prov,
        );
        let seed2 = Seed::canonical(&c2).expect("seed");
        assert_eq!(seed2.runs().len(), 1, "one inv run, not the DAG's shred");
        assert!(matches!(seed2.runs()[0].label, Label::Inv));
    }

    #[test]
    fn the_seed_declines_a_grid_past_the_budget_with_grid_too_large() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // 4200 x 4200 with no common flank: (4201)^2 cells > MAX_SEQFIRST_GRID_CELLS
        // (4097^2). Refused before the grid is allocated, so this is cheap.
        let reference = vec![b'A'; 4200];
        let resulting = vec![b'C'; 4200];
        let c = ctx(&reference, &resulting, &frame, &prov);
        assert_eq!(
            Seed::canonical(&c),
            Err(PartitionError::GridTooLarge {
                ref_len: 4200,
                alt_len: 4200
            })
        );
    }

    #[test]
    fn the_seed_of_identical_sequences_is_the_empty_cut() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        let c = ctx(b"ACGT", b"ACGT", &frame, &prov);
        let seed = Seed::canonical(&c).expect("an empty cut is sound");
        assert!(seed.runs().is_empty());
    }

    #[test]
    fn shrink_narrows_a_widened_run_to_its_difference_hull() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // Window GATTG -> GATG. A run [1,4)->"AT" spells `ATT -> AT`, whose
        // difference hull is the deletion of the last T: [3,4)->"".
        let c = ctx(b"GATTG", b"GATG", &frame, &prov);
        let wide = Cut::new(&c, vec![Run::unlabelled(1, 4, b"AT".to_vec())]).expect("sound");
        let narrow = shrink_to_differences(wide, &c);
        assert_eq!(geometry(&narrow), vec![(3, 4, Vec::new())]);
    }

    #[test]
    fn shrink_leaves_a_no_change_run_whole() {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        // A run that trims away entirely spells no change; the pass leaves it
        // whole rather than emptying it (merge.rs's own rule).
        let c = ctx(b"GATTG", b"GATTG", &frame, &prov);
        let identity = Cut::new(&c, vec![Run::unlabelled(1, 3, b"AT".to_vec())]).expect("sound");
        let same = shrink_to_differences(identity.clone(), &c);
        assert_eq!(same, identity);
    }

    #[test]
    fn the_seed_is_a_fixed_point_of_shrink() {
        // The shrink inside `Seed::canonical` is a no-op on a DAG seed — a run
        // of a minimal alignment is already its own difference hull. Pinned on
        // the two real blocks the adapters already use.
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        for (reference, resulting) in [
            (
                b"TAGTAAACCATTTTACGGAGGATCACAAATTCCTCCTTAT".as_slice(),
                b"TAGTAAACCATTTTACGGAGGATCACACACATTCCTCCTTAT".as_slice(),
            ),
            (b"CAGTGACTAG".as_slice(), b"TGTCACGACT".as_slice()),
        ] {
            let c = ctx(reference, resulting, &frame, &prov);
            let seed = Seed::canonical(&c).expect("seed");
            assert_eq!(shrink_to_differences(seed.clone(), &c), seed);
        }
    }
}
