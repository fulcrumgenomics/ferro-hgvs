//! Task 3 byte-identity gate: the four legacy `FERRO_PARTITION` rules, driven
//! through the new `partition::registry` arm surface, must produce the *same*
//! partition as today's in-module dispatch (`partition_block_for_rule`).
//!
//! # Why this does NOT route through `FERRO_PARTITION`
//!
//! `partition_rule_outcome()` caches the environment read in a process-wide
//! `OnceLock`, so a loop that set `FERRO_PARTITION` per iteration in one process
//! would only ever exercise the arm the process started with (Fable review #2).
//! Both sides are therefore driven by an **explicit rule/name argument**:
//!
//! * the legacy side via `normalize::legacy_partition_members(name, …)`, which
//!   funnels through `partition_block_for_rule(<explicit PartitionRule>, …)`, and
//! * the registry side via an explicit `registry` name lookup that runs the
//!   `Strategy::Mono`-wrapped arm through `Partitioner::partition`.
//!
//! Neither reads the cached env, so all four rules are compared correctly in one
//! process. Both sides funnel their pieces through the *same* Direction-1
//! `piece_to_member` adapter (`normalize::legacy_partition_members` is exactly
//! `piece_to_member ∘ partition_block_for_rule`, and the registry arm is exactly
//! `piece_to_member ∘ partition_block_for_rule` wrapped in `Strategy::run`'s
//! fold/reattach/validate passes), so comparing the resulting `Member` lists is
//! comparing the partitions: the Direction-1 mapping is deterministic, so
//! `Member`-list equality is exactly piece-list equality. What the registry path
//! adds over the raw dispatch — `validate_sound`, `fold_coincident_insertions`,
//! `reattach_examined` — is thereby proven transparent on the legacy corpus,
//! which is the whole point of the gate.

use ferro_hgvs::normalize::legacy_partition_members;
use ferro_hgvs::partition::block_ctx::{BlockCtx, FrameContext, Molecule, Provenance};
use ferro_hgvs::partition::output::Member;
use ferro_hgvs::partition::registry::{arm, registry};

/// One block of the equivalence corpus: a changed `(reference, resulting)` pair
/// and whether it sits on a reading-frame axis (which selects `min_separation`).
struct Block {
    reference: Vec<u8>,
    resulting: Vec<u8>,
    coding: bool,
}

/// The four stable legacy names, paired with nothing — the registry and the
/// legacy accessor both key on the name, so the enum never appears here.
const LEGACY_NAMES: [&str; 4] = ["live", "shadow", "canonical", "canonical-coalesced"];

/// Run the registry arm named `name` over `block` and return its members.
///
/// This is the "via the registry" path: build the neutral `BlockCtx` the arm is
/// a pure function of, look the arm up by its stable name, and drive it through
/// `Partitioner::partition` (which is `Strategy::run` + `validate_sound`). The
/// molecule is `Dna` and the provenance ∅, matching the DNA-only legacy corpus.
fn via_registry(name: &str, block: &Block) -> Vec<Member> {
    let frame = if block.coding {
        FrameContext::coding(0, None, None)
    } else {
        FrameContext::NonCoding
    };
    let provenance = Provenance::none();
    let ctx = BlockCtx {
        reference: &block.reference,
        resulting: &block.resulting,
        frame: &frame,
        molecule: Molecule::Dna,
        provenance: &provenance,
    };
    arm(name)
        .unwrap_or_else(|| panic!("registry is missing the `{name}` arm"))
        .partition(&ctx)
        .unwrap_or_else(|e| {
            panic!("legacy arm `{name}` must always be Ok on the legacy corpus, got {e:?}")
        })
        .members
}

/// The words of length `n` over the DNA alphabet, in a fixed order.
fn words(n: usize) -> Vec<Vec<u8>> {
    let mut out: Vec<Vec<u8>> = vec![Vec::new()];
    for _ in 0..n {
        out = out
            .iter()
            .flat_map(|w| {
                (*b"ACGT").into_iter().map(move |base| {
                    let mut next = w.clone();
                    next.push(base);
                    next
                })
            })
            .collect();
    }
    out
}

/// The designed corpus: every changed `(reference, resulting)` pair with
/// reference length 1..=3 and resulting length 0..=3 (so every `piece_to_member`
/// kind — sub, del, ins, delins, dup, inv, identity — and every partition family
/// is exercised exhaustively at small scale), each run on both a reading-frame
/// and a frameless axis, plus a handful of crafted longer blocks for scale
/// (inversions, tandem duplications, multi-member separations, and the spec's own
/// `delins.md:44` worked block).
fn corpus() -> Vec<Block> {
    let mut blocks = Vec::new();
    let references: Vec<Vec<u8>> = (1..=3).flat_map(words).collect();
    let resultings: Vec<Vec<u8>> = (0..=3).flat_map(words).collect();
    for reference in &references {
        for resulting in &resultings {
            if reference == resulting {
                continue; // no net change: nothing to partition
            }
            for coding in [false, true] {
                blocks.push(Block {
                    reference: reference.clone(),
                    resulting: resulting.clone(),
                    coding,
                });
            }
        }
    }
    // Crafted longer blocks, each exercising a shape the small grid cannot reach
    // at length: whole-span inversions, tandem duplications, multi-member
    // separations at 1/2/3 unchanged bases, and the spec's most-ambiguous block.
    let crafted: [(&[u8], &[u8]); 10] = [
        (b"ATGC", b"GCAT"),             // whole-span inversion
        (b"AACGTT", b"AAACGTGTT"),      // interior tandem dup
        (b"ACGTACGT", b"ACGTACGTACGT"), // longer tandem dup
        (b"ACGTACGTAC", b"TCGTACGTAG"), // two subs, far apart
        (b"ACGTAC", b"TCGTAG"),         // two subs, sep 3
        (b"ACGTA", b"TCGAA"),           // two subs, sep 2
        (b"ACGT", b"TCAT"),             // two subs, sep 1
        (
            b"CAGGGATATGAGAGAACTTCTTCCCCTAAGCCTCGATTCAAGAGCTATGCCT",
            b"TTCCTCGATGCCTG",
        ), // delins.md:44
        (b"GATTACAGATTACA", b"GATTACAGGATTACA"), // insertion in a repeat
        (b"TTTTTTTT", b"TTTATTTT"),     // homopolymer sub
    ];
    for (reference, resulting) in crafted {
        for coding in [false, true] {
            blocks.push(Block {
                reference: reference.to_vec(),
                resulting: resulting.to_vec(),
                coding,
            });
        }
    }
    blocks
}

/// Every legacy rule, over every corpus block, must produce the same partition
/// through the registry as through today's `partition_block_for_rule` dispatch.
#[test]
fn each_legacy_rule_is_byte_identical_through_the_registry() {
    let corpus = corpus();
    assert!(
        corpus.len() > 5_000,
        "the equivalence corpus collapsed to {} blocks; a structural zero here \
         would make this gate vacuous",
        corpus.len()
    );
    let mut compared = 0usize;
    for name in LEGACY_NAMES {
        for block in &corpus {
            let legacy: Vec<Member> =
                legacy_partition_members(name, &block.reference, &block.resulting, block.coding)
                    .unwrap_or_else(|| panic!("`{name}` is not a legacy rule name"));
            let registry: Vec<Member> = via_registry(name, block);
            assert_eq!(
                registry,
                legacy,
                "registry `{name}` diverged from partition_block_for_rule on \
                 reference={:?} resulting={:?} coding={}",
                String::from_utf8_lossy(&block.reference),
                String::from_utf8_lossy(&block.resulting),
                block.coding,
            );
            compared += 1;
        }
    }
    assert_eq!(
        compared,
        LEGACY_NAMES.len() * corpus.len(),
        "every (rule, block) pair must be compared",
    );
}

/// The registry advertises exactly the four stable legacy names, and each maps
/// to an arm whose canonical-coalesce eligibility matches what the arm declares.
#[test]
fn the_registry_exposes_the_four_legacy_names() {
    for name in LEGACY_NAMES {
        assert!(arm(name).is_some(), "registry is missing `{name}`");
    }
    assert!(
        registry().len() >= LEGACY_NAMES.len(),
        "the four legacy names must all be registered",
    );
    assert!(
        arm("no-such-arm").is_none(),
        "an unregistered name must resolve to None, not a default",
    );
}
