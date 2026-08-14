//! Confluence over designed multi-member cis alleles on the `n.` and `r.` axes.
//!
//! # The gap this closes
//!
//! [`super::cis_confluence_axis`] censuses the same generator's `g.` and `c.`
//! output. Those were the two axes #1443 shipped, and they leave the other half
//! of `Normalizer::is_splittable_single_member`'s gate unmeasured: that gate
//! admits `HV::Genome` and `HV::Mt` only, so widening it to `c.`/`n.`/`r.` moves
//! the `n.` and `r.` axes with no committed census to move *against*. The same
//! is true of `FERRO_PARTITION`: the partitioner is axis-blind, but "axis-blind
//! in the source" and "measured on that axis" are different claims, and only the
//! second is evidence.
//!
//! So this module is the `n.`/`r.` half, deliberately kept as a separate corpus
//! and a separate set of pins rather than folded into the `g,c` ones. Folding
//! them in would have re-rolled every number `cis_confluence_axis` pins, which
//! is a change to an existing baseline dressed up as new coverage.
//!
//! # Everything else is the `g,c` module's contract
//!
//! Same generator, same class shape, same census fields, same
//! must-only-move-one-way discipline, same zero-tolerance on
//! `sequence_changed`. Read [`super::cis_confluence_axis`]'s module docs for
//! all of it; only what differs is restated here.
//!
//! # What differs
//!
//! * **The reference.** Both axes are drawn against one `NR_` transcript with
//!   **no CDS**, where `g.`/`c.` use a padded contig and a CDS-bearing `NM_`.
//! * **The alphabet.** An `r.` spelling is lowercase over `acgu`. It is only a
//!   rendering difference: `hgvs_to_spdi` rewrites `u` to `T` before any
//!   sequence comparison, so the stored sequence is DNA on both axes.
//! * **The census is per axis.** `g,c` pins one figure across both of its axes;
//!   here each axis is pinned separately, because the question the axis gate
//!   asks is whether `n.` and `r.` behave *like each other* and like `c.`. A
//!   combined number could hide one axis regressing while the other improved.
//! * **What the reference cannot vary.** Stated rather than left to be
//!   rediscovered, because a census is a claim about its corpus first and about
//!   the normalizer second. The `NR_` transcript is built by
//!   [`reference_for`] as a **single** exon spanning its whole length
//!   (`Exon::with_genomic(1, 1, tx_len, ..)`), so no class here carries an
//!   intronic offset or straddles a junction between two exons, and the
//!   `--axes n,r` generator has no shape that could emit one. What *is*
//!   reachable is the exon's own first and last base, so a boundary clamp there
//!   is exercised. Note the `g,c` corpus is blind in exactly the same way —
//!   `coding_provider` builds one exon too — so this is a property of the shared
//!   generator and not a respect in which the transcript-typed axes are weaker.
//!   Worth stating because the sibling census reads as though it were not:
//!   `cis_confluence_axis::FIVE_PRIME` attributes part of its own #1484 move to
//!   an `enclosing_exon` off-by-one "that let a member sitting on an exon's
//!   first base escape the window clamp and shuffle across the junction". The
//!   first-base half of that is what either corpus can build; the crossing is a
//!   description of the code path that was fixed, not of a shape any class here
//!   or there enumerates. The 64-base core the generator draws
//!   (`generate_cis_confluence_corpus`'s `CORE_LEN`) is likewise far under
//!   `merge::MAX_SPLIT_BLOCK` (1024), so nothing here reaches a length-gated
//!   path — both limits are inherited from the shared generator and hold for
//!   the `g,c` corpus too.
//! * **One pin per axis and direction, not two.** `g,c` carries a second census
//!   per direction for `FERRO_ASSERT_IDEMPOTENT=1`, because two of its `c.`
//!   classes normalize to a non-fixed-point output (#1454) and the oracle turns
//!   those into `declined` panics. Neither axis here has such a class: measured
//!   under all three oracles the four censuses below are unchanged, so there is
//!   no second configuration to pin and no `expected_census` selector. If a
//!   future change gives `n.` or `r.` a non-fixed-point output, these pins go
//!   red in CI's `test-oracle` job while staying green in `Test` — which is the
//!   signal to add the split rather than to re-bless one number.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::PathBuf;
use std::sync::OnceLock;

use rayon::prelude::*;
use serde::Deserialize;

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

use crate::common::cis_apply_oracle::apply_with;
use crate::common::fixture_gen;

const CORPUS_RELATIVE_PATH: &str = "tests/fixtures/cis/cis_confluence_nr_corpus.json";

/// Kept in step with `examples/generate_cis_confluence_corpus.rs`, exactly as
/// [`super::cis_confluence_axis`] keeps its own copies: the test rebuilds the
/// provider rather than the corpus carrying 512 bases of padding per class.
const PAD_OFFSET: usize = 256;
const NR_ACCESSION: &str = "NR_TEST.1";
const TX_CONTIG: &str = "chr_synth";

// ---------------------------------------------------------------------------
// The pinned baseline
// ---------------------------------------------------------------------------

/// The `n.` 3'-direction census, pinned. See the module docs: every divergence
/// figure is a baseline that must only ever go down, and `converged` a floor
/// that must only ever go up.
///
/// **71.2% of designed `n.` classes converge** (4 012 of 5 636); `r.` is
/// 71.1% (4 010, see [`R_THREE_PRIME`]). `cis_confluence_axis` pins **8 026 of
/// 11 272 — also 71.2%** across its combined `g,c` corpus, so the
/// transcript-typed axes now converge at the same rate as the two originals.
///
/// Take that combined figure from `cis_confluence_axis::THREE_PRIME` and not
/// from that module's prose: its headline sentence still reads "58.8% ...
/// (6 633 of 11 272)", which its own later notes have since moved past. That
/// drift is on `main` and predates this file; it is not this module's to
/// correct, but it must not be cited from here either.
///
/// Re-blessed when this branch was rebased onto `main`, as all four pins here
/// were: `converged` 2_628 -> 4_012, `split_two` 2_941 -> 1_558, `split_three`
/// 58 -> 57. `split_more`, `classes` and `spellings` are unchanged and
/// `sequence_changed` stayed 0, so this is the same corpus measured against a
/// moved normalizer rather than a different corpus. The mover is #1484; see
/// [`R_THREE_PRIME`] for the `r.` figures and the same argument in full.
///
/// # What these constants refuted
///
/// Recorded because the belief was load-bearing and will recur. Before the
/// re-bless below, all four censuses here sat at 46.6% (2 624–2 628 of 5 636)
/// against a `g,c` corpus already pinned near 71%, and the obvious reading was
/// that the gap *was* the axis gate:
/// `Normalizer::is_splittable_single_member` admits `HV::Genome` and `HV::Mt`
/// only, the corpus carries a lone spanning `delins` spelling for every class,
/// and on `g.` that spelling is re-derived to agree with the bracketed ones
/// while on a transcript-typed axis it stays whole and diverges from them.
///
/// #1484 closed the gap **without changing that gate's axis list**, so the
/// difference between the axes was never a measurement of what the gate is
/// worth. The mechanism above may still be real; this census no longer
/// evidences it, and nothing here quantifies the gate. Note also that no
/// per-axis `g.`-only or `c.`-only convergence figure is pinned anywhere in
/// the repo — `cis_confluence_axis` pins the two axes combined — so any such
/// split quoted in a comment or a PR body is unmeasured and should not be
/// re-cited.
///
/// Note what the remainder is *not*: `declined` and `sequence_changed` are both
/// zero, so every divergence is two well-formed, sequence-preserving outputs for
/// one variant — the confluence defect itself, not a parse failure or a
/// corrupted sequence wearing its clothes.
/// **Re-blessed again by #1649's two-deletion alignment: 4 012 -> 4 179.**
/// `merge`'s payload splitter gains the *deletion, retained reference, deletion*
/// shape it was missing — the mirror of the two-insertion shape it already had —
/// so a two-gap alignment no longer has to be flattened onto a coarser partition
/// chosen by the input's spelling. `split_two` falls 138, `split_three` 24 and
/// `split_more` 5, summing exactly to the 167 gained, so every moved class
/// converged outright rather than dropping an arity. No divergence figure rises.
///
/// **That this fires on `n.` at all is the finding.** The pass is not
/// coding-axis machinery, and the four censuses in this file move by 167/167/172/172
/// against `cis_confluence_axis`'s 334/344 — the `n.` and `r.` axes gain in the
/// same proportion as the `c.` one, which is what a partitioner change should
/// look like and what an axis-scoped carve-out would not.
/// # Re-blessed by the partition default flip (#1835) — and the licence here is
/// NOT the one the `c.` axis uses
///
/// `converged` **4 179 -> 5 635**, `split_two` 1 420 -> 1, `split_three` and
/// `split_more` to 0. `declined`, `underdetermined` and `sequence_changed` stay
/// **0**, so the ratchet is satisfied and no class converged onto other bases.
///
/// **Read which half of the new default is doing the work, because it is easy to
/// cite the wrong record here.** `CanonicalCoalesced` is two things: the
/// canonical partitioner, and the `delins.md:44-47` coalescing pass on top of
/// it. The pass is **`c.`-only** — `delins-payload-coincidence-carve-out-is-coding-dna-scoped`
/// is decided and puts `g.`/`m.`/`o.`/`n.`/`r.` outside `:47` entirely, so on the
/// two axes this file measures the coalescing arm is *identical to* `Canonical`
/// and contributes nothing. Every class gained here is gained by deriving the
/// partition from the resulting sequence instead of from the input's spelling,
/// which is `canonical-form-choice-when-both-legal` and
/// `derivation-may-not-be-bounded-by-the-inputs-spelling` — both decided, both
/// axis-neutral.
///
/// So a reader must not cite `:47` for this movement, and must not read the
/// movement as evidence that the axis scope leaks. The scope is what makes the
/// residues here differ from the `c.` axis's, not what closed the gap.
///
/// The chronology entry above — "that this fires on `n.` at all is the finding",
/// with the four censuses moving in proportion to the `c.` one — is the same
/// observation one change earlier, and it holds again: `n.` and `r.` gain
/// 1 456 and 1 456 against `cis_confluence_axis`'s 2 910, which is what an
/// axis-neutral partitioner change looks like on a corpus half the size.
/// # #1878 — `converged` FALLS ON ALL FOUR, AND THE CAUSE IS #1440'S
///
/// `n.` 3': 5 635 -> 5 386. `n.` 5': 5 636 -> 5 387. `r.` 3': 5 633 -> 5 384.
/// `r.` 5': 5 634 -> 5 385. `declined` and `sequence_changed` stay 0 on all
/// four, so the whole movement is between `converged` and `split_two`.
///
/// This file's ratchet says `converged` may only go up, so the fall is argued
/// rather than re-blessed. The argument is `cis_confluence_axis`'s, in full and
/// not repeated here: `rulings[equal-length-block-column-correspondence-is-unique]`
/// (decided) moves the OUTPUT, and the input-relative weight bound
/// (`rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]`, decided,
/// tracked as #1440) is what SPLITS the classes. Measured with that bound
/// disabled, all four return to their previous values.
///
/// Note the four move by 249 apiece against the `c.` axis's 498 — the same
/// axis-neutral proportion this section's paragraphs above record for the
/// gains, read in the other direction.
const N_THREE_PRIME: Census = Census {
    classes: 5_636,
    spellings: 23_696,
    declined: 0,
    // #1835 raised these; #1878 lowers them — see the section above.
    converged: 5_386,
    split_two: 250,
    split_three: 0,
    split_more: 0,
    underdetermined: 0,
    sequence_changed: 0,
};

/// The `r.` 3'-direction census, pinned. See [`N_THREE_PRIME`].
///
/// Re-blessed when this branch was rebased onto `main`: `converged` 2_627 ->
/// 4_010, `split_two` 2_941 -> 1_560, `split_three` 59 -> 57. Every figure moved
/// in the direction [`assert_census`] demands, `sequence_changed` stayed 0, and
/// `classes`/`spellings` are unchanged, so the corpus is the same one. The mover
/// is #1484, which changed the `n.`/`r.` axes with no test covering them — this
/// file is the coverage, and the re-bless is the first measurement of what that
/// PR did here.
/// **Re-blessed again by #1649: 4 010 -> 4 177**, an identical +167 to its `n.`
/// twin with identical `split_*` deltas. The `r.` and `n.` axes tracking each
/// other exactly is the cross-check this pair of pins exists for: the alphabet
/// differs and the partition does not.
/// **Re-blessed by the partition default flip (#1835): 4 177 -> 5 633**, with
/// `split_two` 1 422 -> 3 and the two higher arities to 0. See [`N_THREE_PRIME`]
/// for the licensing, and in particular for why `delins.md:47` must **not** be
/// cited on this axis.
///
/// **The `n.`/`r.` twins stop tracking each other exactly, by two classes**, and
/// that is worth flagging because their agreement has been offered as a
/// cross-check in every entry above. `n.` keeps one divergent class here and `r.`
/// keeps three; at 5' the same pair reads 0 and 2.
///
/// **The cause is not measured, and is deliberately not guessed at here.** The
/// obvious hypothesis is the alphabet — an `r.` payload is spelled in lower-case
/// `acgu`, so a coincidence test could in principle land differently — but
/// nothing in this change measured that, and the two censuses have tracked each
/// other exactly through four previous re-blesses, so the divergence is new
/// information rather than a known asymmetry. What *is* asserted is what the
/// ratchet asserts: both figures moved down, `sequence_changed` is 0 on both, and
/// the residue is two classes rather than a family. If a later change needs the
/// mechanism, dump the divergent class ids on the two axes and diff them; do not
/// cite this comment as though it had.
const R_THREE_PRIME: Census = Census {
    classes: 5_636,
    spellings: 23_696,
    declined: 0,
    // #1835 raised these; #1878 lowers them — see `N_THREE_PRIME`.
    converged: 5_384,
    split_two: 252,
    split_three: 0,
    split_more: 0,
    underdetermined: 0,
    sequence_changed: 0,
};

/// The `n.` 5'-direction census, pinned. Confluence is a property of the
/// normalizer rather than of one shuffle direction, so it is measured in full
/// rather than spot-checked.
///
/// The 5' direction is **no longer a public option** — see `README.md` rule 6
/// and `tests/it/five_prime_public_surface_removed.rs`. That changes the reason
/// this census exists, not its scope: the 5' arm is ferro's differential oracle over its own 3' output — the instrument that found #1542, where 7 of 8 `FERRO_PARTITION` x direction configurations agreed and only the shipped `live`/3' arm diverged. An arm that is only spot-checked
/// cannot serve as that oracle.
///
/// Re-blessed alongside [`N_THREE_PRIME`]: `converged` 2_625 -> 4_011,
/// `split_two` 2_954 -> 1_569, `split_three` 53 -> 52, every other field
/// unchanged. Its `r.` twin is [`R_FIVE_PRIME`].
/// **Re-blessed again by #1649: 4 011 -> 4 183** (+172), `split_two` -154,
/// `split_three` -16, `split_more` -2. The 5' direction gains five more classes
/// than the 3' does, the same asymmetry [`N_THREE_PRIME`]'s history records: the
/// finer partition is chosen before the shift, so the members are shuffled
/// independently afterwards and only some of those landings coincide.
/// **Re-blessed by the partition default flip (#1835): 4 183 -> 5 636**, i.e.
/// every class, with all three divergence figures at 0. See [`N_THREE_PRIME`]
/// for the licensing.
///
/// This is the second of the four censuses to reach total convergence — the
/// other is `cis_confluence_axis::FIVE_PRIME`, on the `g,c` corpus. Both are 5'.
/// The 3' direction keeps a residue on both corpora, which is the same
/// direction-asymmetry every entry above records, read at its limit: the
/// partition is chosen before the shift, so a class whose two spellings land
/// apart can only do so in the direction that walks toward the ambiguity.
const N_FIVE_PRIME: Census = Census {
    classes: 5_636,
    spellings: 23_696,
    declined: 0,
    // #1835 raised these; #1878 lowers them — see `N_THREE_PRIME`.
    converged: 5_387,
    split_two: 249,
    split_three: 0,
    split_more: 0,
    underdetermined: 0,
    sequence_changed: 0,
};

/// The `r.` 5'-direction census, pinned. See [`N_FIVE_PRIME`].
///
/// Re-blessed alongside [`R_THREE_PRIME`], and by almost exactly the same
/// margin: `converged` 2_624 -> 4_009, `split_two` 2_954 -> 1_571,
/// `split_three` 54 -> 52. That the two directions moved together (+1_383 and
/// +1_385) is itself the evidence the change is a property of the normalizer
/// rather than of one shuffle direction.
/// **Re-blessed again by #1649: 4 009 -> 4 181**, again matching its `n.` twin
/// exactly at +172. All four censuses in this file move monotonically — every
/// `converged` up, every divergence figure down — so the ratchet is satisfied in
/// the direction it was written for on both axes and both directions.
/// **Re-blessed by the partition default flip (#1835): 4 181 -> 5 634**, with
/// `split_two` 1 417 -> 2 and the higher arities to 0. See [`N_THREE_PRIME`] for
/// the licensing and [`R_THREE_PRIME`] for the `n.`/`r.` residue difference,
/// which appears at 5' as 0 against 2 and is likewise not explained here.
///
/// All four censuses in this file still move monotonically — every `converged`
/// up, every divergence figure down — so the ratchet is satisfied in the
/// direction it was written for, on both axes and both directions, by the
/// largest margin it has ever taken.
const R_FIVE_PRIME: Census = Census {
    classes: 5_636,
    spellings: 23_696,
    declined: 0,
    // #1835 raised these; #1878 lowers them — see `N_THREE_PRIME`.
    converged: 5_385,
    split_two: 251,
    split_three: 0,
    split_more: 0,
    underdetermined: 0,
    sequence_changed: 0,
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
            &["--axes", "n,r"],
            ".cis_confluence_nr_corpus",
            "cis confluence n./r. corpus",
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
///
/// One construction serves both axes: an `NR_` record stores DNA whichever way
/// a spelling renders it, so an `r.` class differs from its `n.` twin only in
/// the spellings, never in the reference. Mirrors the generator's
/// `noncoding_provider`.
fn reference_for(class: &Class) -> (MockProvider, String) {
    let mut provider = MockProvider::new();
    let tx_len = class.core.len() as u64;
    let g_start = PAD_OFFSET as u64 + 1;
    let g_end = PAD_OFFSET as u64 + tx_len;
    let exon = Exon::with_genomic(1, 1, tx_len, g_start, g_end);
    let transcript = Transcript::new(
        NR_ACCESSION.to_string(),
        Some("SYNTH_NR".to_string()),
        Strand::Plus,
        class.core.clone(),
        None,
        None,
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

// ---------------------------------------------------------------------------
// Measurement
// ---------------------------------------------------------------------------

/// What one axis in one direction found. Field-for-field the same shape as
/// [`super::cis_confluence_axis`]'s `Census`, so the two sets of numbers are
/// directly comparable; see that module for what each field means and which way
/// it is allowed to move.
#[derive(Debug, Default, PartialEq, Eq)]
struct Census {
    classes: usize,
    spellings: usize,
    declined: usize,
    converged: usize,
    split_two: usize,
    split_three: usize,
    split_more: usize,
    underdetermined: usize,
    sequence_changed: usize,
}

impl Census {
    /// Fold another census in. Every field is a count over a partition of the
    /// classes, so summing is exact rather than approximate — which is what
    /// makes [`measure`]'s per-group split safe.
    fn absorb(&mut self, other: &Census) {
        self.classes += other.classes;
        self.spellings += other.spellings;
        self.declined += other.declined;
        self.converged += other.converged;
        self.split_two += other.split_two;
        self.split_three += other.split_three;
        self.split_more += other.split_more;
        self.underdetermined += other.underdetermined;
        self.sequence_changed += other.sequence_changed;
    }
}

/// A worst-offender line, for the failure message when a pin moves.
struct Divergence {
    id: String,
    outputs: Vec<String>,
}

/// How many divergent classes the failure message names. Applied per group and
/// again to the concatenation, which is what keeps the parallel result
/// byte-identical to a serial one — see [`measure`].
const WORST_LIMIT: usize = 10;

/// The contiguous runs of `classes` that share one reference and sit on `axis`.
///
/// Classes are sorted by an id that starts with the core index and the axis, so
/// every class sharing a reference is contiguous. Building one provider and one
/// normalizer per group instead of per class is the difference between this
/// running in seconds and in minutes; [`measure`] additionally runs the groups
/// concurrently, which it can do only because each one owns its provider.
fn reference_groups<'a>(classes: &'a [Class], axis: &str) -> Vec<&'a [Class]> {
    let mut groups = Vec::new();
    let mut start = 0usize;
    while start < classes.len() {
        let mut end = start + 1;
        while end < classes.len()
            && classes[end].axis == classes[start].axis
            && classes[end].core == classes[start].core
        {
            end += 1;
        }
        if classes[start].axis == axis {
            groups.push(&classes[start..end]);
        }
        start = end;
    }
    groups
}

/// Census one reference group. Owns its provider and normalizer, reads nothing
/// outside its own slice, and so is safe to run alongside its siblings.
fn measure_group(group: &[Class], direction: ShuffleDirection) -> (Census, Vec<Divergence>) {
    let mut census = Census::default();
    let mut worst: Vec<Divergence> = Vec::new();

    let (provider, reference) = reference_for(&group[0]);
    let normalizer = Normalizer::with_config(
        provider.clone(),
        NormalizeConfig::default().with_direction(direction),
    );

    for class in group {
        census.classes += 1;
        let expected = class.denoted.clone();
        let mut outputs: BTreeSet<String> = BTreeSet::new();
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
            if apply_with(&provider, &reference, &output).as_deref() != Some(expected.as_str()) {
                census.sequence_changed += 1;
            }
            outputs.insert(output);
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
        if outputs.len() > 1 && worst.len() < WORST_LIMIT {
            worst.push(Divergence {
                id: class.id.clone(),
                outputs: outputs.into_iter().collect(),
            });
        }
    }

    (census, worst)
}

/// Census one axis in one direction, one reference group per rayon task.
///
/// **The result is byte-identical to the serial walk, not merely equivalent**,
/// which is the only basis on which a pinned census may be parallelized at all:
///
/// * every `Census` field is a count over a partition of the classes, so
///   [`Census::absorb`] over the groups is exact and order-free;
/// * `par_iter().collect()` preserves *input* order, so the per-group results
///   are folded in corpus order however they finish;
/// * `worst` is "the first [`WORST_LIMIT`] divergent classes in corpus order".
///   Each group keeps its own first `WORST_LIMIT`, and concatenating the groups
///   in corpus order and truncating again yields exactly that set — the serial
///   loop's global cap can only have kept the same ones.
///
/// Nothing is shared across tasks: each group builds its own `MockProvider` and
/// `Normalizer` (it already did — that is the per-group loop's whole point), and
/// the only borrow that crosses is the immutable corpus slice. The normalization
/// oracles' recursion guard is a thread-local, so it is per-task too.
fn measure(direction: ShuffleDirection, axis: &str) -> (Census, Vec<Divergence>) {
    let groups = reference_groups(&corpus().classes, axis);

    let per_group: Vec<(Census, Vec<Divergence>)> = groups
        .par_iter()
        .map(|group| measure_group(group, direction))
        .collect();

    let mut census = Census::default();
    let mut worst: Vec<Divergence> = Vec::new();
    for (group_census, group_worst) in per_group {
        census.absorb(&group_census);
        worst.extend(group_worst);
    }
    worst.truncate(WORST_LIMIT);

    (census, worst)
}

fn report(label: &str, census: &Census, worst: &[Divergence]) -> String {
    let mut out = format!(
        "cis confluence ({label}): {} classes, {} spellings, {} declined\n  \
         converged: {}\n  split 2: {}\n  split 3: {}\n  split 4+: {}\n  \
         underdetermined: {}\n  sequence changed: {}\n",
        census.classes,
        census.spellings,
        census.declined,
        census.converged,
        census.split_two,
        census.split_three,
        census.split_more,
        census.underdetermined,
        census.sequence_changed,
    );
    for divergence in worst {
        out.push_str(&format!(
            "  {} -> {:?}\n",
            divergence.id, divergence.outputs
        ));
    }
    out
}

/// Assert one axis-and-direction census against its pin, printing the measured
/// numbers either way so a moved pin can be re-blessed from the test output.
fn assert_census(direction: ShuffleDirection, axis: &str, label: &str, pinned: &Census) {
    let (measured, worst) = measure(direction, axis);
    println!("{}", report(label, &measured, &worst));
    assert_eq!(
        measured.sequence_changed, 0,
        "{label}: {} normalized outputs no longer denote their class's sequence — a class that \
         converges on the wrong sequence is worse than one that diverges",
        measured.sequence_changed
    );
    assert_eq!(
        &measured,
        pinned,
        "{label}: the confluence census moved. Every divergence figure must only ever go DOWN \
         and `converged` only ever UP; if this change lowers one, re-bless the pin in the same \
         commit and say so in the PR. Measured:\n{}",
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
fn the_nr_corpus_is_a_dense_set_of_real_confluence_classes() {
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
    // This corpus must be the `n,r` one and nothing else. A stale `g,c` artifact
    // installed under this path would make every census below measure zero
    // classes, which `assert_census` does catch — it compares the whole `Census`
    // for equality rather than only checking the one-way movement its failure
    // message describes, so `classes: 0` cannot pass as "went down". What this
    // guard buys is therefore one failure that names the wrong artifact, instead
    // of four opaque census mismatches to work backwards from.
    assert!(
        corpus
            .classes
            .iter()
            .all(|c| c.axis == "n" || c.axis == "r"),
        "the corpus carries an axis other than n./r.; it is the wrong artifact"
    );
    for axis in ["n", "r"] {
        assert!(corpus.classes.iter().any(|c| c.axis == axis));
    }
    for members in [2, 3, 4] {
        assert!(corpus.classes.iter().any(|c| c.members == members));
    }
    for separation in [0, 1, 2, 3, 5, 8] {
        assert!(corpus.classes.iter().any(|c| c.separation == separation));
    }
}

/// `r.` spellings must actually be RNA-rendered, or the axis is `n.` twice over
/// and its census is a duplicate rather than a measurement.
///
/// The two per-spelling checks at the bottom are both **negative** — no
/// uppercase, no `t` — and negatives alone cannot tell an RNA rendering from a
/// lowercasing that never mapped `t` to `u`: that carries no `t` and no capital
/// either, and it is precisely the "`n.` twice over" outcome above. Nor does
/// looking for a `u` fix it, because `dup` supplies one in every spelling that
/// names a duplication.
///
/// So the positive half is asserted against the `n.` twin, which is the actual
/// contract: an `r.` class differs from the `n.` class with the same id only in
/// how its bases are spelled, never in its positions, its members or its
/// reference. `mapped_a_base` then pins that the `t` -> `u` arm was really
/// taken, so a corpus whose payloads happened to be `t`-free could not pass this
/// vacuously.
#[test]
fn the_rna_axis_is_spelled_in_the_rna_alphabet() {
    let corpus = corpus();
    let rna: Vec<&Class> = corpus.classes.iter().filter(|c| c.axis == "r").collect();
    assert!(!rna.is_empty(), "no r. classes to check");

    let by_id: BTreeMap<&str, &Class> = corpus
        .classes
        .iter()
        .filter(|c| c.axis == "n")
        .map(|c| (c.id.as_str(), c))
        .collect();
    let mut mapped_a_base = false;
    for class in &rna {
        let twin_id = class.id.replacen("-r-", "-n-", 1);
        let twin = by_id.get(twin_id.as_str()).unwrap_or_else(|| {
            panic!(
                "{} has no n. twin at {twin_id}, so the two axes are not the same corpus",
                class.id
            )
        });
        assert_eq!(
            class.spellings.len(),
            twin.spellings.len(),
            "{} and {twin_id} are one design on two axes but carry different \
             spelling counts",
            class.id
        );
        for (rendered, dna) in class.spellings.iter().zip(&twin.spellings) {
            let (_, dna_body) = dna
                .split_once(":n.")
                .unwrap_or_else(|| panic!("{dna} is not an n. description"));
            if dna_body.contains('T') {
                mapped_a_base = true;
            }
            // Safe to fold the whole body: positions are digits and every edit
            // keyword (`del`, `dup`, `ins`, `delins`) is already lowercase and
            // `t`-free, so only bases are touched.
            let expected: String = dna_body
                .to_ascii_lowercase()
                .chars()
                .map(|base| if base == 't' { 'u' } else { base })
                .collect();
            let (_, rna_body) = rendered
                .split_once(":r.")
                .unwrap_or_else(|| panic!("{rendered} is not an r. description"));
            assert_eq!(
                rna_body, expected,
                "{rendered} is not {dna} with its bases spelled in acgu"
            );
        }
    }
    assert!(
        mapped_a_base,
        "no n. spelling carried a `T`, so the t -> u mapping was never exercised \
         and these r. classes could be their n. twins with only the case changed"
    );

    for class in rna {
        for spelling in &class.spellings {
            let (_, description) = spelling
                .split_once(":r.")
                .unwrap_or_else(|| panic!("{spelling} is not an r. description"));
            assert!(
                !description.chars().any(|c| c.is_ascii_uppercase()),
                "{spelling} carries an uppercase base, so it is not RNA-rendered"
            );
            assert!(
                !description.contains('t'),
                "{spelling} carries a `t`, which RNA spells `u`"
            );
        }
    }
}

/// [`measure`]'s byte-identity argument rests on [`reference_groups`] returning
/// a **partition** of the axis's classes, in corpus order. Nothing else observes
/// that: the censuses would still sum to something if a class were dropped or
/// counted twice, and `worst` only ever surfaces inside a failure message — so a
/// broken partition is invisible on a green run and would move the pins
/// silently, which is the one way a parallelized census can go wrong without
/// saying so.
#[test]
fn the_reference_groups_partition_each_axis_in_corpus_order() {
    let classes = &corpus().classes;
    for axis in ["n", "r"] {
        let groups = reference_groups(classes, axis);
        assert!(
            !groups.is_empty(),
            "{axis}: no groups, so nothing is measured"
        );

        // Every group is non-empty and shares one reference, which is what lets
        // `measure_group` build a single provider from `group[0]`.
        for group in &groups {
            assert!(
                !group.is_empty(),
                "{axis}: an empty group would panic on group[0]"
            );
            assert!(
                group
                    .iter()
                    .all(|c| c.axis == group[0].axis && c.core == group[0].core),
                "{axis}: a group spans two references, so one provider cannot serve it"
            );
        }

        // Concatenating the groups reproduces the axis's classes exactly, in
        // corpus order — no class dropped, none counted twice, none reordered.
        let flattened: Vec<&str> = groups
            .iter()
            .flat_map(|group| group.iter())
            .map(|c| c.id.as_str())
            .collect();
        let expected: Vec<&str> = classes
            .iter()
            .filter(|c| c.axis == axis)
            .map(|c| c.id.as_str())
            .collect();
        assert_eq!(
            flattened, expected,
            "{axis}: the groups are not a partition of the axis in corpus order"
        );
    }
}

#[test]
fn noncoding_three_prime_confluence_census() {
    assert_census(
        ShuffleDirection::ThreePrime,
        "n",
        "n. 3prime",
        &N_THREE_PRIME,
    );
}

#[test]
fn rna_three_prime_confluence_census() {
    assert_census(
        ShuffleDirection::ThreePrime,
        "r",
        "r. 3prime",
        &R_THREE_PRIME,
    );
}

#[test]
fn noncoding_five_prime_confluence_census() {
    assert_census(ShuffleDirection::FivePrime, "n", "n. 5prime", &N_FIVE_PRIME);
}

#[test]
fn rna_five_prime_confluence_census() {
    assert_census(ShuffleDirection::FivePrime, "r", "r. 5prime", &R_FIVE_PRIME);
}
