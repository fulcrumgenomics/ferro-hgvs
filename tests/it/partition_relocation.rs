//! Behavior-equivalence pin for the arm relocation (plan Task 1).
//!
//! The move of the arm machinery from `bakeoff/` to `partition/` re-signatures the
//! composition layer from `&Case` to `&BlockCtx`. "Behavior-preserving" is a claim
//! that needs evidence, so this test pins it: the four R2 winner arms
//! (`trim_only` / `axis_peel` / `peel_comp` / `peel_strict`), over a fixed set of
//! `(reference, resulting)` blocks — a substitution, a deletion, an insertion, and
//! a net-insertion/duplication-shaped block — must produce byte-identical
//! `Vec<Member>` via the NEW `BlockCtx`-based `Strategy::run`.
//!
//! The `expected` members below were captured from the PRE-MOVE build, running
//! each winner arm through the OLD `Strategy::run(&Case)` API over these exact
//! blocks (a throwaway dumper, since removed). That the same bytes fall out of the
//! new `&BlockCtx` API IS the behavior-preservation proof: the two APIs bracket
//! the refactor, and equality across them shows the plumbing change moved no
//! output. (All four arms happen to agree on these blocks — they share one L2
//! typer and Dial-B config and differ only in L1 walls, which coincide here — so
//! each row's `expected` is one vector asserted against all four arms.)

use ferro_hgvs::partition::arms::r2_winner;
use ferro_hgvs::partition::block_ctx::{BlockCtx, FrameContext, Molecule, Provenance};
use ferro_hgvs::partition::metrics::{round_trips, RefApplier};
use ferro_hgvs::partition::output::{EditKind, Member};
use ferro_hgvs::partition::Strategy;

fn member(kind: EditKind, ref_start: usize, ref_end: usize, inserted: &[u8]) -> Member {
    Member {
        kind,
        ref_start,
        ref_end,
        inserted: inserted.to_vec(),
    }
}

/// `(block name, reference, resulting, expected members)` — one row per block.
type Block = (&'static str, &'static [u8], &'static [u8], Vec<Member>);

/// The four required block shapes, with members captured from the pre-move build.
fn blocks() -> Vec<Block> {
    vec![
        (
            "sub",
            b"AAACAAAGAAA",
            b"AAATAAAGAAA",
            vec![member(EditKind::Sub, 3, 4, b"T")],
        ),
        (
            "del",
            b"ACGTTTGCA",
            b"ACGTGCA",
            vec![member(EditKind::Del, 4, 6, b"")],
        ),
        (
            "ins",
            b"ACGTACGT",
            b"ACGTGGACGT",
            vec![member(EditKind::Ins, 4, 4, b"GG")],
        ),
        (
            "net_ins_dup",
            b"ACGTTCAGGTCACAATT",
            b"ACGTTCAGGTCACACACT",
            vec![
                member(EditKind::Dup, 14, 14, b"CA"),
                member(EditKind::Delins, 14, 16, b"C"),
            ],
        ),
    ]
}

fn winners() -> Vec<(&'static str, Strategy)> {
    vec![
        ("trim_only", r2_winner::trim_only()),
        ("axis_peel", r2_winner::axis_peel()),
        ("peel_comp", r2_winner::peel_comp()),
        ("peel_strict", r2_winner::peel_strict()),
    ]
}

#[test]
fn winner_arms_produce_identical_partitions_after_the_move() {
    for (block_name, reference, resulting, expected) in blocks() {
        let frame = FrameContext::NonCoding;
        let provenance = Provenance::none();
        let ctx = BlockCtx {
            reference,
            resulting,
            frame: &frame,
            molecule: Molecule::from_axis("g"),
            provenance: &provenance,
        };
        for (arm_name, arm) in winners() {
            // `run` now takes &BlockCtx; byte-identity vs the pre-move capture is
            // the binding proof that the plumbing change moved no output.
            let got = arm.run(&ctx);
            assert_eq!(
                got.members, expected,
                "arm {arm_name} on block {block_name}: members drifted vs the pre-move capture",
            );
            assert!(
                !got.members.is_empty(),
                "arm {arm_name} on block {block_name}: empty partition",
            );
            assert!(
                round_trips(&RefApplier, &ctx, &got),
                "arm {arm_name} on block {block_name}: partition does not round-trip",
            );
        }
    }
}
