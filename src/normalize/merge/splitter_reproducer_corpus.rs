//! The harvested reproducer corpus, as executable rows.
//!
//! Each row is one reproducer reduced to the only thing [`partition_block`]
//! can see: the maximally affix-trimmed `(ref_block, alt_block)` pair, plus the
//! member count recorded for that case. The splitter is a pure function of that
//! pair, so this table is the measurement a future splitter change is scored
//! against — it says what today's splitter does on every block the corpus
//! reaches, not what any candidate design ought to do.
//!
//! # Why this lives here and not in `tests/it/`
//!
//! [`partition_block`] and [`Piece`](super::Piece) are private to `merge.rs`.
//! Reaching them from an integration test would mean widening their visibility
//! for test convenience, so the corpus is a `#[cfg(test)]` child module of
//! `merge` instead — child modules see their ancestors' private items, and
//! nothing in the shipping API surface moves.
//!
//! # `recorded_members` vs. what the splitter returns
//!
//! `recorded_members` is the member count recorded for the case's *normalized
//! output*. For most rows the block alone determines it and the assertion is
//! direct. Nine rows carry a `splitter_returns` override: on those, the block
//! does **not** determine the recorded count. They fall in three groups, and the
//! third is a finding rather than a known design seam:
//!
//! 1. **Two-member spellings** (`#422-two-member`, `#999-neg-split`,
//!    `#1260a-split`, `#1262a-split`, `#1262b-split`). Each shares its block
//!    with a spanning sibling row in this table that the splitter answers
//!    identically. Their recorded two-member outputs therefore cannot come from
//!    the splitter; see the collision census below.
//!
//!    `#1260b` is the exception, and it is the override *inverted*: the two-gap
//!    alignment (#1260, PR #1285) made the splitter answer the **split** count,
//!    so `#1260b-spanning` is the row that now carries the override.
//! 2. **Coding-axis codon merges** (`#1235-c-codon-exception`, `#1235-r-coding`).
//!    The block splits into two; the one-member coding answer is produced after
//!    partitioning, so it is not a block-level fact.
//! 3. **`#1262-window`** — `("AAAAAAA", "ACAAAA")`, the untrimmed poly-A window
//!    pinned at one member for the sequence-first partitioner
//!    (`seqfirst::partition`'s `calibration_corpus_partitions_as_measured`).
//!    `partition_block` returns **two** pieces for the same input. The two
//!    splitters disagree on this window. Its `#1260` sibling
//!    (`("AAAAAAA", "AACACAAAA")`) agrees at one member, so the disagreement is
//!    specific to this pair, not to untrimmed windows in general.
//!
//! # Collision census
//!
//! Six block pairs each carry two different recorded member counts. Since a
//! pure function of the pair returns one answer, at most one of each pair's two
//! recorded counts can be a splitter decision. [`COLLISIONS`] pins all six:
//! two are recorded with **both** answers asserted correct in the repo, three
//! remain **known defects**, and one — `#1260b` — is **resolved at the
//! splitter**, which now answers its split count rather than its spanning one.
//!
//! # Recorded but deliberately not asserted here
//!
//! - **`#1157` case B, template instance** (`TEMPLATE:g.10_20delinsTC` /
//!   `g.12_20del`, and the `g.19_29delinsTG` second instance) — STALE. The
//!   idempotence violation the issue records no longer reproduces. Its two
//!   VERIFIED siblings, `#1157-B-del` and `#1157-B-dup`, are rows below.
//! - **`#1158`** — STALE. The reported `NotEquivalent` no longer reproduces, and
//!   the committed test deliberately asserts only `level.is_equivalent()` rather
//!   than a rung. It is an equivalence-level case and derives no block pair of
//!   its own; its blocks are `#1157`'s.
//! - **`LRG_199t1:c.850_901` end to end** — UNVERIFIABLE in the harvest: the
//!   ferro-prepared reference volumes were offline, so the accession could not
//!   be resolved. Only its block-level form is asserted (`LRG_199` below), from
//!   the literal already committed in `seqfirst::partition` and
//!   `seqfirst::align`. No bases were reconstructed from memory.
//! - **`#1234`'s distant-sibling negative control** (`g.[4_6del;60T>C]`) — its
//!   two members are 54 nt apart, so the splitter never receives that span as
//!   one block.
//! - **`#1235`'s mixed-type and overlap-refusal cases** — the harvest derives no
//!   block pair for them; they are member-typing and refusal behaviour.
//!
//! The five `#1286`/`#1297`/`#1316`/`#1321`/`#1323` rows are recorded as
//! producing no `partition_block` call at all on either spelling — they converge
//! through the repeat/`dup` path. They are kept as rows because their trimmed
//! blocks pin what the splitter *would* have no choice about, not because the
//! splitter decides them.

use super::{partition_block, Piece};

/// The `LRG_199t1:c.850_901` reference block, 52 nt of real sequence, copied
/// from the literal already committed in `seqfirst::partition` (`LRG199_REF`)
/// and `seqfirst::align`'s `lrg199_has_exactly_one_dominator_match`.
const LRG199_REF: &str = "CAGGGATATGAGAGAACTTCTTCCCCTAAGCCTCGATTCAAGAGCTATGCCT";

/// The `LRG_199` payload, 14 nt, from the same committed literal.
const LRG199_ALT: &str = "TTCCTCGATGCCTG";

/// One harvested reproducer, reduced to the pair the splitter sees.
struct BlockRow {
    /// Issue or case identifier, as the harvest names it.
    id: &'static str,
    /// The maximally affix-trimmed reference block.
    ref_block: &'static str,
    /// The trimmed replacement block.
    alt_block: &'static str,
    /// Member count recorded for this case's normalized output.
    recorded_members: usize,
    /// Set only where `partition_block` does not return `recorded_members`:
    /// the count it does return, and where the difference comes from. The row
    /// is asserted at this value; `recorded_members` stays as recorded.
    splitter_returns: Option<(usize, &'static str)>,
    /// One line about the block itself.
    note: &'static str,
}

impl BlockRow {
    /// The member count this row asserts against the live splitter.
    fn expected(&self) -> usize {
        self.splitter_returns
            .map_or(self.recorded_members, |(n, _)| n)
    }
}

/// Shorthand for a row whose block determines the recorded member count.
const fn row(
    id: &'static str,
    ref_block: &'static str,
    alt_block: &'static str,
    recorded_members: usize,
    note: &'static str,
) -> BlockRow {
    BlockRow {
        id,
        ref_block,
        alt_block,
        recorded_members,
        splitter_returns: None,
        note,
    }
}

/// Shorthand for a row where the block does not determine the recorded count.
const fn diverging(
    id: &'static str,
    ref_block: &'static str,
    alt_block: &'static str,
    recorded_members: usize,
    splitter_returns: (usize, &'static str),
    note: &'static str,
) -> BlockRow {
    BlockRow {
        id,
        ref_block,
        alt_block,
        recorded_members,
        splitter_returns: Some(splitter_returns),
        note,
    }
}

/// Reason string shared by the six two-member-spelling rows.
const SPELLING: &str = "recorded count is the two-member spelling's output; the block \
                        is shared with this row's spanning sibling";

/// The corpus.
///
/// Ordered as the harvest orders it: cluster A (contiguous replacements that
/// must stay one member), cluster B (blocks that must split), cluster C, then
/// the calibration set.
#[rustfmt::skip]
const CORPUS: &[BlockRow] = &[
    // --- cluster A: contiguous replacements that must stay ONE member -------
    row("#1034", "CTG", "ACA", 1,
        "revcomp(CTG) = CAG != ACA, so the 2 nt sub-run TG->CA is not carved out"),
    row("#1034b", "TCC", "GAG", 1,
        "the 2 nt prefix TC->GA is a revcomp; the run still stays whole"),
    row("#1034-guard-3nt", "GCT", "AGC", 1,
        "whole-run revcomp: one piece, which is what lets it render as inv"),
    row("#1034-guard-2nt", "GG", "CC", 1, "whole-run revcomp, 2 nt"),
    row("#1040", "CAGTGACTAG", "TGTCACGACT", 1,
        "10 nt run whose interior GTGA->TCAC is a revcomp; not carved"),
    row("#1040-affix-trim", "GA", "TC", 1,
        "the shared AC.. / ..GT affixes are already trimmed off here"),
    row("#1040-length-changing", "GCT", "AG", 1,
        "a length-changing block can never be an inv"),
    row("#1040-two-invs", "TGACA", "CAATG", 2,
        "the unchanged A at offset 2 separates two revcomp pairs"),
    row("#1041", "GCA", "TGC", 1, "whole-run revcomp: revcomp(GCA) = TGC"),
    row("#160-1092", "GG", "CC", 1, "the only part of #160 that survived #1034"),
    row("#160-1150", "TCC", "GAG", 1, "same block as #1034b, reached from #160"),
    row("#160-1099", "TGCTC", "AAGCT", 1,
        "interior 3 nt GCT->AGC is a revcomp; the whole run is not"),
    row("#182-a", "GATCC", "TCTAA", 2, "unchanged base at offset 2"),
    row("#182-b", "CCTGA", "AATTC", 2, "unchanged base at offset 2"),
    row("#182-c", "CCTGATCC", "AATTCTAA", 3, "unchanged bases at offsets 2 and 5"),
    row("#182-d", "CTGA", "ATTC", 2, "a single-base run stays its own piece"),
    row("#182-e", "GATCTC", "TCTATA", 3, "unchanged bases at offsets 2 and 4"),
    row("#422-spanning", "CTATAG", "AAACCCC", 1,
        "6 nt replaced by 7 nt with one coincidentally preserved interior base"),
    diverging("#422-two-member", "CTATAG", "AAACCCC", 2, (1, SPELLING),
        "same block as #422-spanning: collision, both counts pinned correct"),

    // --- cluster B: blocks that must SPLIT ----------------------------------
    row("#1232", "CAATT", "TA", 2, "5 nt -> 2 nt; the minimal edit set is not unique"),
    row("#1231", "AAT", "CAA", 2, "offsets 0 and 2 change, offset 1 does not"),
    row("#1233", "AAT", "TAA", 2, "offsets 0 and 2 change, offset 1 does not"),
    row("#1157-A-padded", "AGTCAGT", "GATTA", 3,
        "7 nt -> 5 nt splitting at two unchanged bases: the cluster-B calibration"),
    row("#1157-A-template", "AGTCACC", "GATTA", 3,
        "the issue body's own instance of the same shape, one base different"),
    row("#1157-B-del", "AA", "", 1, "pure 2 nt deletion inside an A-tract"),
    row("#1157-B-dup", "", "G", 1, "pure 1 nt insertion inside a G-tract"),
    row("#1229", "CAC", "AGT", 1,
        "three consecutive changed nt are one piece; revcomp(CAC) = GTG != AGT"),
    row("#1230", "GATG", "CATC", 2,
        "only offsets 0 and 3 change; the interior AT is untouched"),
    row("#1235-n", "CAA", "TAG", 2, "no reading frame, so the separation stands"),
    diverging("#1235-c-codon-exception", "CAA", "TAG", 1,
        (2, "the coding-axis one-amino-acid merge runs after partitioning"),
        "the same block as #1235-n and #1241, on a coding axis"),
    row("#1235-c-codon-boundary", "ATT", "GTA", 2,
        "separation 1 on a coding axis still splits when the codon differs"),
    diverging("#1235-r-coding", "CAA", "TAG", 1,
        (2, "the coding-axis one-amino-acid merge runs after partitioning"),
        "the r. sibling of #1235-c-codon-exception, same block"),
    row("#1241", "CAA", "TAG", 2,
        "the non-coding r. half: no CDS, so no codon frame, same block"),

    // --- cluster C ----------------------------------------------------------
    row("#1234", "ACCA", "T", 1, "4 nt -> 1 nt with no preserved interior base"),
    row("#1234-lone-del", "CCA", "", 1, "the lone-deletion negative control"),

    // --- calibration --------------------------------------------------------
    row("LRG_199", LRG199_REF, LRG199_ALT, 1,
        "the spec's own worked example: 52 nt -> 14 nt must stay one member"),
    row("LRG_199-standin", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", LRG199_ALT, 1,
        "the periodic synthetic stand-in for the row above, same net-deletion regime"),
    row("long-delins-40nt",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        "CATGCATGCATGCATGCATTCATGCATGCATGCATGCATG", 2,
        "equal length, so every column is compared in place; only offset 19 matches"),
    diverging("#1260a-spanning", "A", "CAC", 1,
        (2, "the two-gap alignment keeps the reference A between the two insertions"),
        "net +2; the splitter now answers the SPLIT count"),
    row("#1260a-split", "A", "CAC", 2,
        "same block as #1260a-spanning: the collision the two-gap alignment resolved"),
    row("#1260-window", "AAAAAAA", "AACACAAAA", 1,
        "the UNTRIMMED 7 nt poly-A window, 9 nt on the alt side, so above \
         MAX_SINGLE_BASE_SEPARATION_BLOCK and held whole. The pipeline never \
         sees this: trim_common_flanks reduces it to #1260a's (A, CAC), which \
         is 3 nt and does split -- which is how #1260 converges"),
    diverging("#1260b-spanning", "AA", "CAAC", 1,
        (2, "the two-gap alignment keeps the whole reference between the two insertions"),
        "the separation-2 instance of #1260; the splitter now answers the SPLIT count"),
    row("#1260b-split", "AA", "CAAC", 2,
        "same block as #1260b-spanning: the collision the two-gap alignment resolved"),
    row("#1262a-spanning", "AA", "C", 1, "net -1 with no preserved reference base"),
    diverging("#1262a-split", "AA", "C", 2, (1, SPELLING),
        "same block as #1262a-spanning: collision, recorded as a known defect"),
    diverging("#1262-window", "AAAAAAA", "ACAAAA", 1,
        (2, "the sequence-first partitioner pins one member for this same input"),
        "the untrimmed poly-A window where the two splitters disagree"),
    diverging("#1262b-spanning", "AA", "CAC", 1,
        (2, "ins + sub outranks a spanning delins, so the split stands"),
        "the sub-plus-insertion instance of #1262; the splitter now answers the SPLIT count"),
    row("#1262b-split", "AA", "CAC", 2,
        "same block as #1262b-spanning: the collision the separation threshold resolved"),
    row("#999-positive", "G", "CT", 1,
        "an insertion shifted flush against a substitution leaves no gap"),
    diverging("#999-neg-spanning", "GC", "CGA", 1,
        (2, "ins + sub outranks a spanning delins, so the required split stands"),
        "net +1 with exactly one preserved interior base, like #422-spanning"),
    row("#999-neg-split", "GC", "CGA", 2,
        "same block as #999-neg-spanning: the splitter now derives the REQUIRED \
         two-member count from sequence alone, with no help from the spelling"),
    diverging("#1000", "C", "GCA", 1,
        (2, "base for base #1260's `A -> CAC`: two insertions either side of the \
             retained C, re-blessed with it"),
        "a 2 nt insertion shifted flush against a substitution"),

    // --- spelling-confluence gap: trimmed blocks are pure insertions --------
    row("#1287", "", "GAAA", 1, "pure insertion: an empty ref block cannot be partitioned"),
    row("#1290", "", "CA", 1, "pure insertion"),
    row("#1296", "C", "A", 1,
        "the only one of the thirteen whose trimmed block is a 1 nt substitution"),
    row("#1301", "", "CAAA", 1, "pure insertion"),
    row("#1304", "", "GA", 1, "pure insertion"),
    row("#1308", "", "TTGG", 1, "pure insertion"),
    row("#1312", "", "AACC", 1, "pure insertion"),
    row("#1320", "", "CAAAAA", 1, "pure insertion"),
    row("#1286", "", "AA", 1, "pure insertion; recorded as reaching no partition_block call"),
    row("#1297", "", "AA", 1, "pure insertion; recorded as reaching no partition_block call"),
    row("#1316", "", "AAAA", 1, "pure insertion; recorded as reaching no partition_block call"),
    row("#1321", "", "GA", 1, "pure insertion; recorded as reaching no partition_block call"),
    row("#1323", "", "AGA", 1, "pure insertion; recorded as reaching no partition_block call"),
];

/// How the repo records the two answers a colliding block pair carries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Pinning {
    /// Both recorded answers are asserted correct somewhere in the repo.
    BothCorrect,
    /// The divergence is recorded as a defect that must stay visible.
    KnownDefect,
    /// The defect is **fixed at the splitter**: the splitter now answers the
    /// *split* row's count rather than the spanning one, so the two spellings
    /// no longer disagree about how to divide the block. Kept as a row rather
    /// than deleted, because the pair is what makes the resolution legible —
    /// and because the two recorded counts still differ until the rest of the
    /// pipeline follows the splitter.
    ResolvedToSplit,
}

/// One collision: two corpus rows that share a block pair but carry different
/// recorded member counts.
struct Collision {
    /// The row whose recorded count the splitter reproduces.
    spanning: &'static str,
    /// The row whose recorded count the splitter does not reproduce.
    split: &'static str,
    pinning: Pinning,
}

/// The collision census: six pairs, two recorded correct, four recorded as
/// known defects.
const COLLISIONS: &[Collision] = &[
    Collision {
        spanning: "#422-spanning",
        split: "#422-two-member",
        pinning: Pinning::BothCorrect,
    },
    Collision {
        spanning: "#999-neg-spanning",
        split: "#999-neg-split",
        pinning: Pinning::ResolvedToSplit,
    },
    Collision {
        spanning: "#1260a-spanning",
        split: "#1260a-split",
        pinning: Pinning::ResolvedToSplit,
    },
    Collision {
        spanning: "#1260b-spanning",
        split: "#1260b-split",
        pinning: Pinning::ResolvedToSplit,
    },
    Collision {
        spanning: "#1262a-spanning",
        split: "#1262a-split",
        pinning: Pinning::KnownDefect,
    },
    Collision {
        spanning: "#1262b-spanning",
        split: "#1262b-split",
        pinning: Pinning::ResolvedToSplit,
    },
];

/// Look up a row by id, panicking with the id when it is missing so a renamed
/// row fails loudly rather than silently dropping a collision.
fn find(id: &str) -> &'static BlockRow {
    CORPUS
        .iter()
        .find(|r| r.id == id)
        .unwrap_or_else(|| panic!("no corpus row with id {id}"))
}

/// Substitute `pieces` into `reference` and return the sequence they denote.
///
/// This is the splitter's own contract stated independently: the pieces must be
/// disjoint, ascending, and reconstruct the alt block exactly.
fn apply(reference: &[u8], pieces: &[Piece]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut cursor = 0;
    for piece in pieces {
        assert!(
            piece.ref_start >= cursor,
            "pieces must be disjoint and ascending; got {pieces:?}"
        );
        assert!(
            piece.ref_end <= reference.len(),
            "piece runs past the reference; got {pieces:?}"
        );
        out.extend_from_slice(&reference[cursor..piece.ref_start]);
        out.extend_from_slice(&piece.alt);
        cursor = piece.ref_end;
    }
    out.extend_from_slice(&reference[cursor..]);
    out
}

/// Every block in the corpus, against the live splitter.
#[test]
fn corpus_blocks_partition_as_the_live_splitter_does() {
    let mut wrong = Vec::new();
    for row in CORPUS {
        let pieces = partition_block(row.ref_block.as_bytes(), row.alt_block.as_bytes());
        if pieces.len() != row.expected() {
            wrong.push(format!(
                "  {}: ({}, {}) expected {} pieces, got {} {:?}\n    recorded as: {}",
                row.id,
                row.ref_block,
                row.alt_block,
                row.expected(),
                pieces.len(),
                pieces,
                row.note
            ));
        }
    }
    assert!(wrong.is_empty(), "corpus rows moved:\n{}", wrong.join("\n"));
}

/// The pieces of every corpus block reconstruct that block's alt exactly.
///
/// A member count alone would be satisfied by a partition that denotes the
/// wrong sequence, so the count assertions above are only worth as much as this.
#[test]
fn corpus_pieces_reconstruct_every_alt_block() {
    for row in CORPUS {
        let reference = row.ref_block.as_bytes();
        let pieces = partition_block(reference, row.alt_block.as_bytes());
        assert_eq!(
            String::from_utf8_lossy(&apply(reference, &pieces)),
            row.alt_block,
            "{} pieces do not denote the alt block: {pieces:?}",
            row.id
        );
    }
}

/// The rows where the block does not determine the recorded member count.
///
/// Pinned as an exact set: a change that closes one of these — or opens a new
/// one — must say so here rather than move a count quietly.
#[test]
fn recorded_and_splitter_counts_diverge_on_exactly_these_rows() {
    let diverging: Vec<&str> = CORPUS
        .iter()
        .filter(|r| r.splitter_returns.is_some())
        .map(|r| r.id)
        .collect();
    assert_eq!(
        diverging,
        vec![
            "#422-two-member",
            "#1235-c-codon-exception",
            "#1235-r-coding",
            "#1260a-spanning",
            "#1260b-spanning",
            "#1262a-split",
            "#1262-window",
            "#1262b-spanning",
            "#999-neg-spanning",
            "#1000",
        ]
    );
    for row in CORPUS.iter().filter(|r| r.splitter_returns.is_some()) {
        assert_ne!(
            row.recorded_members,
            row.expected(),
            "{} carries an override that matches the recorded count",
            row.id
        );
    }
}

/// The collision census, executable.
///
/// For each pair: the two rows really are the same splitter input, they carry
/// different recorded member counts, and the splitter returns one answer for
/// both. That is what makes the pair unsatisfiable for a splitter that sees
/// only the block — at most one of the two recorded counts can be its decision.
#[test]
fn colliding_block_pairs_are_one_splitter_input_with_two_recorded_answers() {
    for collision in COLLISIONS {
        let spanning = find(collision.spanning);
        let split = find(collision.split);

        assert_eq!(
            (spanning.ref_block, spanning.alt_block),
            (split.ref_block, split.alt_block),
            "{} and {} are recorded as a collision but are different blocks",
            spanning.id,
            split.id
        );
        assert_ne!(
            spanning.recorded_members, split.recorded_members,
            "{} and {} record the same member count, so they do not collide",
            spanning.id, split.id
        );

        // Which side the splitter reproduces is the census's whole content, so
        // it is asserted per-pinning rather than assumed to be the spanning one.
        let pieces = partition_block(spanning.ref_block.as_bytes(), spanning.alt_block.as_bytes());
        let (reproduced, not_reproduced) = match collision.pinning {
            Pinning::ResolvedToSplit => (split, spanning),
            Pinning::BothCorrect | Pinning::KnownDefect => (spanning, split),
        };
        assert_eq!(
            pieces.len(),
            reproduced.recorded_members,
            "{} is the recorded count the splitter reproduces",
            reproduced.id
        );
        assert_ne!(
            pieces.len(),
            not_reproduced.recorded_members,
            "{} is the recorded count the splitter does not reproduce",
            not_reproduced.id
        );
    }
}

/// One collision is recorded with both answers correct (`#422`, the deliberate
/// pin), one remains a known defect (`#1262a`), and **four are now resolved at
/// the splitter** — the two-gap alignment plus a separation threshold of one
/// unchanged base made the splitter derive the split count from sequence alone.
/// That split is the census's whole point, so it is pinned.
#[test]
fn the_census_is_one_correct_pin_one_known_defect_and_four_resolved() {
    let correct: Vec<&str> = COLLISIONS
        .iter()
        .filter(|c| c.pinning == Pinning::BothCorrect)
        .map(|c| c.split)
        .collect();
    assert_eq!(correct, vec!["#422-two-member"]);

    let defects: Vec<&str> = COLLISIONS
        .iter()
        .filter(|c| c.pinning == Pinning::KnownDefect)
        .map(|c| c.split)
        .collect();
    assert_eq!(defects, vec!["#1262a-split"]);

    let resolved: Vec<&str> = COLLISIONS
        .iter()
        .filter(|c| c.pinning == Pinning::ResolvedToSplit)
        .map(|c| c.split)
        .collect();
    assert_eq!(
        resolved,
        vec![
            "#999-neg-split",
            "#1260a-split",
            "#1260b-split",
            "#1262b-split"
        ]
    );
}

/// The corpus is a table, so a duplicated or renamed id would silently shadow a
/// case. Ids are unique; shared *blocks* are expected and are not deduplicated.
#[test]
fn corpus_ids_are_unique() {
    let mut ids: Vec<&str> = CORPUS.iter().map(|r| r.id).collect();
    let total = ids.len();
    ids.sort_unstable();
    ids.dedup();
    assert_eq!(ids.len(), total, "duplicate corpus id");
    assert_eq!(total, 65, "corpus size changed; update the module doc");
}
