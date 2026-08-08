//! Confluence over designed multi-member cis alleles (#1443).
//!
//! # The gap
//!
//! `canonicalize_from_sequence` is gated on `members.len() > 1`, so only a
//! multi-member cis allele reaches the partitioner at all. Part 1 of #1443
//! harvested every such input from the bulk corpora and found **650 rows out of
//! 9 949 738** — real corpora are observed variants, where two changes rarely
//! land close enough to share a block. Consumers submitting *designed* alleles
//! have the opposite distribution, so the partitioner is effectively untested
//! and no real corpus can fix that.
//!
//! `examples/generate_cis_confluence_corpus.rs` fills it. Each class it emits is
//! one synthetic reference, one denoted sequence, and N ≥ 2 distinct spellings
//! verified — by an applier independent of the normalizer — to denote that
//! sequence.
//!
//! # READ THIS FIRST: what these numbers measure changed, and the target did too
//!
//! This module used to ask one question — **does a class reach exactly one
//! output?** — and treat the answer as a defect count with a target of 100%.
//! That framing belonged to the re-derivation partitioner, under which a
//! description's members were only a hint: `canonicalize_from_sequence` threw
//! the input's own partition away, minimised edit distance over the resulting
//! *sequence*, and so answered every spelling of one sequence identically by
//! construction. Sequence-confluence was the right target because the model
//! made the description's own assertion irrelevant.
//!
//! `FERRO_PARTITION` now defaults to `preserve`, and under it **a description
//! asserts a partition of the reference into changed blocks and unchanged runs**.
//! `g.19_33delinsCGG` asserts one changed block; `g.[19_23del;27_33del]` asserts
//! two. Those are different assertions about the same resulting sequence, and
//! `partition_block_preserving` keeps each, moving a boundary only where
//! `general.md:34-35` / `DNA/delins.md:44-47` license it (merge two members
//! closer than the axis floor) or where `general.md:34` requires it (split an
//! **equal-length** member whose interior holds an unchanged run reaching that
//! floor).
//!
//! So **two spellings asserting different partitions are no longer required to
//! converge, and a fall in sequence-confluence is not by itself a defect.** The
//! `converged` figures below fell when the default flipped — 3': 6 629 ->
//! 5 611; 5': 6 626 -> 5 684, A/B-measured against `FERRO_PARTITION=live`, which
//! reproduces the old numbers exactly — and essentially all of that fall is the
//! lone-`delins` spelling of each class no longer being re-derived into the
//! multi-member form (see the corpus's shape, below).
//!
//! ## The three questions, and which of them is a gate
//!
//! | figure | what it measures | is it a target? |
//! |---|---|---|
//! | `sequence_changed`, `not_fixed_point` | correctness: an output denotes its class's sequence, and re-reading it re-asserts the same partition | **yes, absolute zero** |
//! | `same_partition_*` | confluence *within* one asserted partition — an upper bound, see [`input_arity`] | a watched residual, not a ratchet |
//! | `converged`, `split_*`, `cross_partition_divergence` | sequence-confluence across *different* assertions | **no.** A change detector whose movement must be argued, in either direction |
//!
//! Only the first row is a property the partition model promises. The second is
//! the property that *replaces* the old target, but the grouping key available
//! to a test is a proxy, so it bounds the residual rather than counting it. The
//! third is retained because it is what a downstream consumer experiences —
//! a re-normalization — and it must never move silently.
//!
//! ## Why the corpus makes the fall almost arithmetic
//!
//! Every one of the 11 272 classes is exactly **one single-member spelling plus
//! a set of spellings that all carry the same member count** (2, 3 or 4);
//! measured, and pinned by
//! [`the_corpus_is_a_dense_set_of_real_confluence_classes`]. So each class poses
//! exactly one cross-partition question (lone `delins` against the multi-member
//! group) and at most one same-partition question (the multi-member group
//! against itself). `cross_partition_divergence` is 5 281 of 5 661 diverging
//! classes at 3': the lone spelling disagreeing with the group it never asserted
//! the same partition as.
//!
//! Read the pins together with `CLAUDE.md`'s note on representation stability.
//! The bar there is **confluence plus disclosure, not stability** — but note
//! which confluence: the confluence that blocks the downstream consumer is
//! between two spellings *of one variant as they store it*, and where those two
//! spellings assert different partitions this model answers them differently on
//! purpose. That is a representation change and it owes the release a
//! `dump_normalized_corpus` measurement and a `Representation-Change:` trailer,
//! whichever direction the strings moved.
//!
//! # Sequence preservation is the harder half, and it survived the model change
//!
//! A class that converged on the *wrong* sequence would be worse than one that
//! diverged, so every output is also applied back to the reference — through
//! `hgvs_to_spdi`, the same way `tests/it/common/synthetic.rs`'s
//! `assert_padded_preserving` does it — and compared with the class's denoted
//! sequence. That count is asserted at zero, and is still zero under
//! `preserve`. Its sibling `not_fixed_point` — added with the partition model,
//! and also zero — asks the dual question: re-reading ferro's own answer must
//! re-assert the same partition. Those two zeros are what make the fall in
//! `converged` a change of representation rather than of meaning.
//!
//! # Corpus
//!
//! The corpus is a generated artifact (gitignored, like the spec fixtures) and
//! is regenerated on demand through `common::fixture_gen`, so a fresh checkout
//! just works. Its parameters are what the pins are measured over: regenerating
//! with different `--seeds` re-rolls every number here.
//!
//! # HISTORY — the notes below describe `FERRO_PARTITION=live`
//!
//! Everything from here down was written under the re-derivation partitioner
//! and is kept because `live` is still selectable and still A/B'd against. Do
//! not read these numbers as the current pins; the constants above are the
//! current pins, and each carries its own `live` baseline for the comparison.
//!
//! # The pins were once per oracle configuration, and no longer are (#1454)
//!
//! Each direction used to carry two pinned censuses, selected by an
//! `expected_census` helper: one for a plain run and one for
//! `FERRO_ASSERT_IDEMPOTENT=1`. Two classes normalized to a non-fixed-point
//! output (#1454), so the oracle turned them into panics that were counted
//! `declined` rather than two distinct outputs, and one pin could only ever be
//! green in one of CI's two jobs.
//!
//! #1454 is fixed, `declined` is 0 in both configurations, and the two censuses
//! are now byte-identical — so there is one constant per direction again.
//!
//! **The collapse did not land where this module predicted, and the difference
//! is the useful part.** The old text expected `converged` to rise by two. It
//! rises by *one* per direction (3': 6632 -> 6633; 5': 6627 -> 6628), and it
//! *falls* by one against the oracle pins (6634 -> 6633; 6629 -> 6628).
//!
//! Measured per class rather than inferred from the net figure — a `+1` is also
//! what "converged three, broke two" looks like, and the two must be told apart.
//! Diffing every class's output set across the fix: **exactly one class changes
//! convergence status, and none regresses.**
//!
//! ```text
//! s04-c-m2-sep3-p4-rot4   the pure #1454 class
//!   before  [13_15delinsTAA;16_19delinsTAAT] | [13_15delinsTAA;16_19inv]
//!   after   [13_15delinsTAA;16_19inv]                            CONVERGED
//!
//! s00-c-m2-sep5-p1-rot4   the issue's own reproducer
//!   before  [10_12delinsTAA;13_15delinsTAT]  | [9dup;15del]
//!   after   [10_12delinsTAA;13_15inv]        | [9dup;15del]      still splits
//! ```
//!
//! The second class had the #1454 defect too, and the fix corrects it — that
//! member is a fixed point now. The class nonetheless still splits, because its
//! competing output is a different **partition** of the same variant rather than
//! a different typing of one member, and no amount of inversion typing merges
//! `[9dup;15del]` with a two-member delins. That residual belongs to the
//! #1235 / #1419-#1421 partition-non-uniqueness family and is deliberately left
//! alone here: choosing between those two partitions moves shipped
//! representations, which is its own measurement.
//!
//! So the `+2` was never reachable by fixing #1454. It was an artifact: under
//! the oracle, `s00-c-m2-sep5-p1-rot4`'s delins spelling *panicked*, leaving
//! `[9dup;15del]` as the only surviving outcome, and a class with one outcome
//! scores as `converged`. The oracle was crediting the defect with a convergence
//! that was really a suppressed output.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::PathBuf;
use std::sync::OnceLock;

use serde::Deserialize;

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

use crate::common::cis_apply_oracle::apply_with;
use crate::common::fixture_gen;

const CORPUS_RELATIVE_PATH: &str = "tests/fixtures/cis/cis_confluence_corpus.json";

/// Kept in step with `examples/generate_cis_confluence_corpus.rs`. The test
/// rebuilds the providers rather than the corpus carrying 512 bases of padding
/// per class, so the two must agree on the padding, the accessions and the CDS.
const PAD_OFFSET: usize = 256;
const GENOMIC_CONTIG: &str = "NC_TEST.1";
const TX_ACCESSION: &str = "NM_TEST.1";
const TX_CONTIG: &str = "chr_synth";
const CDS_START: u64 = 1;
const CDS_END: u64 = 63;

// ---------------------------------------------------------------------------
// The pinned baseline
// ---------------------------------------------------------------------------

/// The 3'-direction census, pinned.
///
/// # Sequence-confluence: **49.8%** (5 611 of 11 272), down from 58.8%
///
/// Read the module docs before reading this as a regression. Under
/// `FERRO_PARTITION=live` — the re-derivation partitioner, still selectable and
/// A/B'd for exactly this — the same corpus gives `converged` 6 629, `split 2`
/// 4 509, `split 3` 115, `split 4+` 19. Flipping the default to `preserve`
/// moves it to 5 611 / 5 524 / 135 / 2, a fall of **1 018 classes**.
///
/// Of the 5 661 classes that now diverge, **5 281 diverge only across different
/// asserted partitions** (`cross_partition_divergence`) — the lone `delins`
/// spelling answers one question and the multi-member group answers another.
/// `live` converged those by discarding the lone spelling's assertion and
/// re-deriving the group's; `preserve` cannot, and is not meant to. This is the
/// single largest representation change the partition rule makes and it is the
/// one to declare in the changelog.
///
/// # The three correctness figures are still hard zeros
///
/// `declined` 0, `sequence_changed` 0 and `not_fixed_point` 0, over 47 392
/// spellings. So every one of the 5 661 divergences is a pair of well-formed,
/// sequence-preserving, self-stable descriptions — nothing here is a parse
/// failure, a corrupted sequence, or an unstable answer wearing a divergence's
/// clothes. `not_fixed_point` is new and is the strongest single result in this
/// file: the partition rule's central claim is that reading ferro's own answer
/// back re-asserts the same partition, and at corpus scale it does, exactly.
///
/// # The same-partition residual went UP, and "210 new defects" is the wrong reading
///
/// `same_partition_divergence` 205 -> 380 (`same_partition_converged` 11 067 ->
/// 10 892). That direction deserves suspicion and got it: the classes this test
/// prints were read end to end, and they are the **proxy** rather than defects.
/// The grouping key is the input's member count ([`input_arity`]), and under
/// `live` that proxy was tight for a reason unrelated to correctness —
/// re-derivation collapsed *every* same-arity spelling onto one string, so a
/// group could hardly split. Under `preserve` two spellings can carry the same
/// member count while asserting members on different territory, and then they
/// legitimately answer differently. Worked instance, straight out of this test's
/// own failure report:
///
/// ```text
/// s00-c-m2-sep0-p4-rot3      both inputs are two-member
///   c.[9T>A;10_13dup]     -> c.[8_9insAAA;11dup]    members at 9 and gap 13: separated, kept
///   c.[9T>A;9_10insAATA]  -> c.9delinsAAATA         members at 9 and gap 9: touching, merged
/// ```
///
/// Same arity, different territory, therefore different partitions, therefore
/// two answers — the merge in the second row is the licensed move firing, not a
/// disagreement. So this figure bounds the residual from above and must be read
/// as a bound. Sharpening it needs a shift-invariant territory signature, which
/// a test at this level cannot compute without re-implementing the shuffle.
///
/// # What the pre-`preserve` history of this constant recorded
///
/// Kept because the reasoning is still sound about `live`, which is still
/// selectable: `converged` was 6 632 before #1454 and 6 633 after, then #1524
/// took it to 6 629 deliberately (`s00-c-m3-sep0-p8-rot4`,
/// `s00-c-m4-sep0-p8-rot4`, `s02-c-m4-sep1-p4-rot1`, `s03-c-m3-sep1-p8-rot2`)
/// by removing a `build_split_variants` step that re-split a member into two
/// *touching* members which the legacy per-member merge then re-merged wider.
/// Note where that argument now lands: the four classes it cost were exactly
/// lone-`delins`-against-multi-member disagreements, which is the whole
/// `cross_partition_divergence` bucket above. The #1524 note called the
/// multi-member form "the *better* of the two forms" on `general.md:34`
/// grounds; the partition rule reaches the same conclusion structurally instead
/// of by cost.
/// # RE-MEASURED 2026-08-08, after the partition rule was made to govern every
/// block. Read this block before reading anything above as current.
///
/// Everything above describes a **blend**: `partition_block_preserving` used to
/// decline on 33.4% of invocations and fall back to `partition_block`, so a third
/// of the figures were the re-derivation partitioner's. The decline rate is now
/// 0.0% in both directions, and the census moved again as a result. Measured on
/// this tree, both arms, same corpus and same instrument:
///
/// ```text
///                                  live (re-derive)   preserve (this pin)
///   converged                                 8 026                3 416
///   split 2                                   3 112                7 552
///   split 3                                     115                  283
///   split 4+                                     19                   21
///   same-partition converged                 11 116               10 459
///   cross-partition divergence                3 090                7 043
///   same-partition divergence                   156                  813
///   sequence changed                              0                    0
///   not a fixed point                             0                    0
/// ```
///
/// So the representation change to declare is **4 610 classes losing
/// sequence-confluence at 3'** (8 026 -> 3 416) and 4 518 at 5' (8 021 ->
/// 3 503). The previous pins (5 611 / 5 684) were neither arm — they were the
/// blend, so the branch-to-branch delta of 2 195 understates the change against
/// the release by half.
///
/// **The three correctness figures are still hard zeros**, over 47 392 spellings:
/// `declined` 0, `sequence_changed` 0, `not_fixed_point` 0. That is what makes
/// the fall above a representation change rather than a defect: every one of the
/// 7 856 divergences is a pair of well-formed, sequence-preserving, self-stable
/// descriptions.
///
/// ## AND A ZERO HERE IS A CLAIM ABOUT THIS CORPUS, NOT ABOUT THE TREE
///
/// This must be said next to the zeros, because it was nearly missed.
/// `sequence_changed: 0` and `not_fixed_point: 0` hold on the designed cis corpus
/// while `cis_junction_crossing_shift`'s exhaustive sweeps — on the same tree, on
/// the same day — report **187 776** two- and three-member cases whose output
/// overlaps or denotes no single sequence, every one of them a repeat or `dup`
/// grown over a sibling's junction (`g.[2_3dup;5_6insA]` -> `g.[1_9T[11];5_6insA]`).
/// The corpus cannot build that shape: `generate_cis_confluence_corpus` draws its
/// members from a fixed rotation set over a 63-nt core and never places a
/// junction-only member strictly inside a homopolymer tract that another member
/// grows. So these two zeros are real and narrow. They are not evidence that the
/// tree preserves sequence.
const THREE_PRIME: Census = Census {
    classes: 11_272,
    spellings: 47_392,
    declined: 0,
    converged: 3_416,
    split_two: 7_552,
    split_three: 283,
    split_more: 21,
    underdetermined: 0,
    sequence_changed: 0,
    same_partition_groups: 11_272,
    same_partition_converged: 10_459,
    cross_partition_divergence: 7_043,
    same_partition_divergence: 813,
    not_fixed_point: 0,
};

/// The 5'-direction census, pinned. `--direction 5prime` is a supported public
/// option, and both confluence and partition preservation are properties of the
/// normalizer rather than of one shuffle direction, so it is measured in full
/// rather than spot-checked.
///
/// It tracks the 3' figure closely under `preserve` (5 684 against 5 611) as it
/// did under `live` (6 626 against 6 629), which is the reading to keep: what
/// changed is which *question* the axis answers, not which end of an ambiguous
/// run the shuffle walks to. A rule change that moved only one of the two
/// numbers would be treating a symptom.
///
/// The `live` baseline for the A/B is `converged` 6 626, `split 2` 4 532,
/// `split 3` 105, `split 4+` 9, `same_partition_converged` 11 079,
/// `same_partition_divergence` 193.
///
/// Two differences from 3' worth naming: 5' converges slightly *more* here
/// (5 684 against 5 611), and `same_partition_divergence` is slightly *higher*
/// (389 against 380). Both follow from the cause the residual note on
/// [`THREE_PRIME`] sets out — which anchor a member shuffles to decides whether
/// two same-arity spellings land on the same territory — so neither is a
/// direction-specific defect.
/// RE-MEASURED 2026-08-08 with [`THREE_PRIME`]; the A/B table and the
/// corpus-blindness caveat there apply verbatim. The `live` baseline for this
/// direction is `converged` 8 021, `split 2` 3 137, `split 3` 105, `split 4+` 9,
/// `same_partition_converged` 11 132, `cross_partition_divergence` 3 111,
/// `same_partition_divergence` 140.
///
/// It still tracks 3' closely (3 503 against 3 416), which remains the reading to
/// keep: what the partition rule changed is which *question* the axis answers,
/// not which end of an ambiguous run the shuffle walks to.
const FIVE_PRIME: Census = Census {
    classes: 11_272,
    spellings: 47_392,
    declined: 0,
    converged: 3_503,
    split_two: 7_562,
    split_three: 194,
    split_more: 13,
    underdetermined: 0,
    sequence_changed: 0,
    same_partition_groups: 11_272,
    same_partition_converged: 10_540,
    cross_partition_divergence: 7_037,
    same_partition_divergence: 732,
    not_fixed_point: 0,
};

// ---------------------------------------------------------------------------
// Corpus
// ---------------------------------------------------------------------------

/// One confluence class. A structural mirror of the generator's own `Class`;
/// an integration test cannot import a type out of an example, so the two are
/// kept in step by this file failing to deserialize if they drift.
#[derive(Deserialize)]
struct Class {
    id: String,
    axis: String,
    core: String,
    denoted: String,
    members: usize,
    separation: usize,
    spellings: Vec<String>,
}

#[derive(Deserialize)]
struct Corpus {
    classes: Vec<Class>,
}

fn corpus_path() -> PathBuf {
    fixture_gen::fixture_path(CORPUS_RELATIVE_PATH)
}

fn corpus() -> &'static Corpus {
    static CORPUS: OnceLock<Corpus> = OnceLock::new();
    CORPUS.get_or_init(|| {
        let path = corpus_path();
        fixture_gen::ensure_generated_example_fixture(
            &path,
            "generate_cis_confluence_corpus",
            ".cis_confluence_corpus",
            "cis confluence corpus",
            || {},
        );
        let text = fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
        serde_json::from_str(&text)
            .unwrap_or_else(|e| panic!("failed to parse {}: {e}", path.display()))
    })
}

// ---------------------------------------------------------------------------
// Providers — the same construction the generator used
// ---------------------------------------------------------------------------

fn padded(core: &str) -> String {
    let pad = "ACGT".repeat(PAD_OFFSET / 4);
    format!("{pad}{core}{pad}")
}

/// The provider a class is drawn against, and the full sequence its coordinates
/// address.
fn reference_for(class: &Class) -> (MockProvider, String) {
    if class.axis == "g" {
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence(GENOMIC_CONTIG, padded(&class.core));
        return (provider, padded(&class.core));
    }
    let mut provider = MockProvider::new();
    let tx_len = class.core.len() as u64;
    let g_start = PAD_OFFSET as u64 + 1;
    let g_end = PAD_OFFSET as u64 + tx_len;
    let exon = Exon::with_genomic(1, 1, tx_len, g_start, g_end);
    let transcript = Transcript::new(
        TX_ACCESSION.to_string(),
        Some("SYNTH".to_string()),
        Strand::Plus,
        class.core.clone(),
        Some(CDS_START),
        Some(CDS_END),
        vec![exon],
        Some(TX_CONTIG.to_string()),
        Some(g_start),
        Some(g_end),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );
    provider.add_genomic_sequence(TX_CONTIG, padded(&class.core));
    provider.add_transcript(transcript);
    (provider, class.core.clone())
}

/// The sequence a class denotes, on the axis's own reference.
fn expected_sequence(class: &Class) -> String {
    if class.axis == "g" {
        padded(&class.denoted)
    } else {
        class.denoted.clone()
    }
}

// ---------------------------------------------------------------------------
// Measurement
// ---------------------------------------------------------------------------

/// What one direction's run found. Every field is pinned; see the module docs
/// for which way each is allowed to move.
#[derive(Debug, Default, PartialEq, Eq)]
struct Census {
    /// Classes measured.
    classes: usize,
    /// Spellings normalized across all classes.
    spellings: usize,
    /// Spellings normalization declined (or panicked on). Not folded into the
    /// distinct-output count: a decline is a different failure from a
    /// divergence and conflating them would let one hide the other.
    declined: usize,
    /// Classes whose surviving spellings all reached **one** output. The goal
    /// is for this to equal `classes`.
    converged: usize,
    /// Classes that split into exactly two distinct outputs.
    split_two: usize,
    /// … exactly three.
    split_three: usize,
    /// … four or more.
    split_more: usize,
    /// Classes left with fewer than two spellings that normalized at all, so
    /// confluence is not decidable for them. Counted rather than silently
    /// skipped: a rise here would hollow the axis out while every other number
    /// improved.
    underdetermined: usize,
    /// Outputs whose applied sequence differs from the class's. Asserted at
    /// zero — a class that converges on the wrong sequence is worse than one
    /// that diverges.
    sequence_changed: usize,
    /// Same-partition groups with at least two spellings that normalized: the
    /// denominator of the measurement this module now leads with. See
    /// [`input_arity`] for why the input's member count is the grouping key.
    same_partition_groups: usize,
    /// … of those, the groups that reached **one** output. Under the partition
    /// model this is the figure that must go up, because every spelling in a
    /// group asserts the same number of changed blocks over the same variant.
    same_partition_converged: usize,
    /// Diverging classes where every same-partition group nevertheless
    /// converged — so the class's two outputs are two *different assertions*
    /// about where the changed blocks are, not two answers to one question.
    /// Expected to be large under `preserve` and **not** a defect count.
    cross_partition_divergence: usize,
    /// Diverging classes where some same-partition group itself split — an
    /// **upper bound** on same-partition non-confluence, not a defect count.
    /// See [`input_arity`] for why the bound is loose and the constants for
    /// what the bound measured before and after.
    same_partition_divergence: usize,
    /// Outputs that are not fixed points: `norm(norm(x)) != norm(x)`. Asserted
    /// at zero, absolutely rather than as a baseline. This is the one property
    /// the partition model makes a *hard* prediction about at corpus scale — a
    /// description asserts a partition, and re-reading ferro's own answer must
    /// re-assert the same one — so unlike every figure above it there is no
    /// legitimate non-zero value. #1454 was an instance and is the reason the
    /// pins were once per oracle configuration.
    not_fixed_point: usize,
}

/// A worst-offender line, for the failure message when a pin moves.
///
/// Same-partition divergences are listed first and are what a reader should
/// look at: a cross-partition split is two answers to two different questions
/// and tells nobody anything, while a same-partition split is the residual this
/// axis is now for.
struct Divergence {
    id: String,
    same_partition: bool,
    outputs: Vec<String>,
    /// `spelling -> output`, one line per spelling. Only worth printing for a
    /// same-partition divergence, where the question is which two spellings
    /// that asserted the same number of blocks disagreed.
    detail: Vec<String>,
}

/// How many changed blocks a description **asserts** — its member count.
///
/// This is the grouping key for the same-partition census, and it is a proxy
/// rather than the partition itself, deliberately stated as such. Two
/// descriptions with different member counts *cannot* assert the same partition
/// (a partition with `k` changed blocks is not one with `j != k`), so grouping
/// by arity never puts two different partitions' answers into one group's
/// numerator by mistake in the direction that matters — it can only be too
/// *coarse*, folding two same-arity-but-different-territory assertions together
/// and so **over**-reporting `same_partition_divergence`. That bias is the safe
/// one: the residual it reports is an upper bound on the defects, never a floor
/// that hides one.
///
/// The generated corpus makes the grouping exact in practice: every class is one
/// single-member spelling plus a set of spellings that all carry the same
/// member count (measured — `the_corpus_is_a_dense_set_of_real_confluence_classes`
/// pins it), so each class contributes at most one decidable group.
///
/// Read off the string rather than the parsed variant: what a partition model
/// preserves is the *description's* assertion, and re-deriving it from the code
/// under test would check that code against itself. No corpus spelling nests a
/// bracket (a repeat unit) inside the allele brackets, which the same corpus
/// contract pins, so counting separators is exact here.
fn input_arity(spelling: &str) -> usize {
    let Some((_, body)) = spelling.split_once(':') else {
        return 1;
    };
    let Some(coordinates) = body.get(2..) else {
        return 1;
    };
    if coordinates.starts_with('[') {
        coordinates.matches(';').count() + 1
    } else {
        1
    }
}

fn measure(direction: ShuffleDirection) -> (Census, Vec<Divergence>) {
    let classes = &corpus().classes;
    let mut census = Census::default();
    let mut worst: Vec<Divergence> = Vec::new();

    // Classes are sorted by an id that starts with the core index and the axis,
    // so every class sharing a reference is contiguous. Building one provider
    // and one normalizer per group instead of per class is the difference
    // between this running in seconds and in minutes.
    let mut start = 0usize;
    while start < classes.len() {
        let mut end = start + 1;
        while end < classes.len()
            && classes[end].axis == classes[start].axis
            && classes[end].core == classes[start].core
        {
            end += 1;
        }
        let (provider, reference) = reference_for(&classes[start]);
        let normalizer = Normalizer::with_config(
            provider.clone(),
            NormalizeConfig::default().with_direction(direction),
        );

        for class in &classes[start..end] {
            census.classes += 1;
            let expected = expected_sequence(class);
            let mut outputs: BTreeSet<String> = BTreeSet::new();
            // Outputs bucketed by the member count the *input* asserted, so a
            // class's spellings can be asked the partition-model question as
            // well as the sequence-confluence one. The `usize` counts the
            // spellings that reached the bucket, which a `BTreeSet` of outputs
            // cannot: two agreeing spellings collapse to one output.
            let mut by_partition: BTreeMap<usize, (usize, BTreeSet<String>)> = BTreeMap::new();
            let mut detail: Vec<String> = Vec::new();
            let mut normalized_spellings = 0usize;
            for spelling in &class.spellings {
                census.spellings += 1;
                let Ok(variant) = parse_hgvs(spelling) else {
                    census.declined += 1;
                    continue;
                };
                let normalized = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    normalizer.normalize(&variant)
                }));
                let output = match normalized {
                    Ok(Ok(value)) => value.to_string(),
                    _ => {
                        census.declined += 1;
                        continue;
                    }
                };
                normalized_spellings += 1;
                if apply_with(&provider, &reference, &output).as_deref() != Some(expected.as_str())
                {
                    census.sequence_changed += 1;
                }
                // Re-normalizing ferro's own answer must be a no-op. Measured
                // here rather than left to `FERRO_ASSERT_IDEMPOTENT`, which is
                // compiled out in release builds — and this axis runs in
                // release.
                let again = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    parse_hgvs(&output).ok().map(|v| normalizer.normalize(&v))
                }));
                match again {
                    Ok(Some(Ok(second))) if second.to_string() == output => {}
                    _ => census.not_fixed_point += 1,
                }
                let arity = input_arity(spelling);
                detail.push(format!("m{arity} {spelling} -> {output}"));
                let bucket = by_partition.entry(arity).or_default();
                bucket.0 += 1;
                bucket.1.insert(output.clone());
                outputs.insert(output);
            }

            // The same-partition census. A group of one spelling decides
            // nothing and is not counted either way, exactly as a class of one
            // is `underdetermined` above.
            let mut every_group_converged = true;
            for (spellings_in_group, group) in by_partition.values() {
                if *spellings_in_group < 2 {
                    continue;
                }
                census.same_partition_groups += 1;
                if group.len() == 1 {
                    census.same_partition_converged += 1;
                } else {
                    every_group_converged = false;
                }
            }

            if normalized_spellings < 2 {
                census.underdetermined += 1;
                continue;
            }
            match outputs.len() {
                1 => census.converged += 1,
                2 => census.split_two += 1,
                3 => census.split_three += 1,
                _ => census.split_more += 1,
            }
            if outputs.len() > 1 {
                if every_group_converged {
                    census.cross_partition_divergence += 1;
                } else {
                    census.same_partition_divergence += 1;
                }
                worst.push(Divergence {
                    id: class.id.clone(),
                    same_partition: !every_group_converged,
                    outputs: outputs.into_iter().collect(),
                    detail,
                });
            }
        }
        start = end;
    }

    (census, worst)
}

fn report(direction: &str, census: &Census, worst: &[Divergence]) -> String {
    let mut out = format!(
        "cis confluence ({direction}): {} classes, {} spellings, {} declined\n  \
         converged: {}\n  split 2: {}\n  split 3: {}\n  split 4+: {}\n  \
         underdetermined: {}\n  sequence changed: {}\n  \
         same-partition groups: {}\n  same-partition converged: {}\n  \
         cross-partition divergence: {}\n  same-partition divergence: {}\n  \
         not a fixed point: {}\n",
        census.classes,
        census.spellings,
        census.declined,
        census.converged,
        census.split_two,
        census.split_three,
        census.split_more,
        census.underdetermined,
        census.sequence_changed,
        census.same_partition_groups,
        census.same_partition_converged,
        census.cross_partition_divergence,
        census.same_partition_divergence,
        census.not_fixed_point,
    );
    // Same-partition first, and capped: the whole list is thousands of classes,
    // almost all of them cross-partition and so uninformative.
    let mut listed = 0usize;
    for same_partition in [true, false] {
        for divergence in worst.iter().filter(|d| d.same_partition == same_partition) {
            if listed >= 10 {
                break;
            }
            listed += 1;
            out.push_str(&format!(
                "  [{}] {} -> {:?}\n",
                if same_partition { "same" } else { "cross" },
                divergence.id,
                divergence.outputs
            ));
            if same_partition {
                for line in &divergence.detail {
                    out.push_str(&format!("      {line}\n"));
                }
            }
        }
    }
    out
}

/// Assert one direction's census against its pin, printing the measured numbers
/// either way so a moved pin can be re-blessed from the test output.
///
/// The two zeros are asserted first and separately, because they are the only
/// figures here with a *correct* value rather than a measured one: a class that
/// converges on the wrong sequence, or an answer that will not survive being
/// read back, is a defect whatever the confluence figures say, and reporting it
/// as "the census moved" would bury it.
fn assert_census(direction: ShuffleDirection, label: &str, pinned: &Census) {
    let (measured, worst) = measure(direction);
    println!("{}", report(label, &measured, &worst));
    assert_eq!(
        measured.sequence_changed, 0,
        "{label}: {} normalized outputs no longer denote their class's sequence — a class that \
         converges on the wrong sequence is worse than one that diverges",
        measured.sequence_changed
    );
    assert_eq!(
        measured.not_fixed_point, 0,
        "{label}: {} normalized outputs are not fixed points — reading ferro's own answer back \
         re-asserts a different partition, which is the one thing the partition model promises \
         outright",
        measured.not_fixed_point
    );
    assert_eq!(
        &measured,
        pinned,
        "{label}: the census moved. Read the module docs before re-blessing: these figures do \
         NOT all have a good direction any more.\n  \
         * `sequence_changed` and `not_fixed_point` are absolute zeros — non-zero is a defect, \
         never a re-bless.\n  \
         * `same_partition_converged` up / `same_partition_divergence` down is progress, but the \
         grouping key is a proxy (see `input_arity`), so movement needs classes read, not just \
         counted.\n  \
         * `converged` and the `split_*` figures measure sequence-confluence ACROSS different \
         asserted partitions. Under `FERRO_PARTITION=preserve` two such spellings are two \
         different assertions and are not required to agree, so a fall here is a representation \
         change to declare, not automatically a regression. A/B it with `FERRO_PARTITION=live` \
         and say in the PR which forms moved and how many.\n\
         Measured:\n{}",
        report(label, &measured, &worst)
    );
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// The corpus is only worth measuring if it is actually dense in the shape it
/// claims to cover, so its own contract is checked before anything is read out
/// of it. A generator change that quietly emitted a hundred classes would
/// otherwise show up as an improved convergence rate.
#[test]
fn the_corpus_is_a_dense_set_of_real_confluence_classes() {
    let corpus = corpus();
    assert!(
        corpus.classes.len() >= 5_000,
        "the corpus collapsed to {} classes; it is meant to be thousands",
        corpus.classes.len()
    );
    for class in &corpus.classes {
        assert!(
            class.spellings.len() >= 2,
            "{} is a singleton, which is not a confluence class",
            class.id
        );
        assert_eq!(
            class.spellings.iter().collect::<BTreeSet<_>>().len(),
            class.spellings.len(),
            "{} repeats a spelling",
            class.id
        );
        assert_ne!(
            class.denoted, class.core,
            "{} denotes its own reference, so it describes no variant",
            class.id
        );
    }
    // Both axes, every member count and every separation the generator
    // enumerates must be present — a parameter contributing nothing would make
    // its share of the census meaningless.
    for axis in ["g", "c"] {
        assert!(corpus.classes.iter().any(|c| c.axis == axis));
    }
    for members in [2, 3, 4] {
        assert!(corpus.classes.iter().any(|c| c.members == members));
    }
    for separation in [0, 1, 2, 3, 5, 8] {
        assert!(corpus.classes.iter().any(|c| c.separation == separation));
    }

    // The shape the same-partition census rests on: every class is one
    // single-member spelling plus a set that all carry one member count. Pinned
    // rather than assumed, because it is what makes each class pose exactly one
    // cross-partition question and at most one same-partition question — and a
    // generator change that emitted, say, both two- and three-member spellings
    // per class would silently re-mean `same_partition_groups` while every
    // other number here still looked plausible.
    let mut single_member_spellings = 0usize;
    for class in &corpus.classes {
        let arities: BTreeSet<usize> = class.spellings.iter().map(|s| input_arity(s)).collect();
        assert_eq!(
            arities.len(),
            2,
            "{} carries member counts {arities:?}; a class is meant to be one single-member \
             spelling against one multi-member arity",
            class.id
        );
        assert!(
            arities.contains(&1),
            "{} has no single-member spelling, so it poses no cross-partition question",
            class.id
        );
        single_member_spellings += class
            .spellings
            .iter()
            .filter(|s| input_arity(s) == 1)
            .count();
        assert!(
            class
                .spellings
                .iter()
                .filter(|s| input_arity(s) > 1)
                .count()
                >= 2,
            "{} has fewer than two multi-member spellings, so its same-partition group decides \
             nothing",
            class.id
        );
    }
    assert_eq!(
        single_member_spellings,
        corpus.classes.len(),
        "exactly one single-member spelling per class"
    );
}

#[test]
fn three_prime_confluence_census() {
    assert_census(ShuffleDirection::ThreePrime, "3prime", &THREE_PRIME);
}

#[test]
fn five_prime_confluence_census() {
    assert_census(ShuffleDirection::FivePrime, "5prime", &FIVE_PRIME);
}
