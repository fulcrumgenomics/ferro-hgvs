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
//! sequence. This module normalizes every spelling in every class and asks the
//! confluence question: **does a class reach exactly one output?**
//!
//! # This is a census, not a pass/fail gate
//!
//! It does not converge today, and a permanently red suite would be worth
//! nothing. So the numbers below are **pinned baselines**: the target is total
//! convergence, and each pinned figure must only ever move in the direction of
//! it. A change that raises a divergence count is a regression and must be
//! explained; a change that lowers one should re-bless the number in the same
//! commit.
//!
//! Read the pins together with `CLAUDE.md`'s note on representation stability.
//! Confluence and stability are different properties, and a fix for the first
//! moves shipped strings — so a commit that lowers a divergence count here still
//! owes the release its `dump_normalized_corpus` measurement.
//!
//! # Sequence preservation is the harder half
//!
//! A class that converged on the *wrong* sequence would be worse than one that
//! diverged, so every output is also applied back to the reference — through
//! `hgvs_to_spdi`, the same way `tests/it/common/synthetic.rs`'s
//! `assert_padded_preserving` does it — and compared with the class's denoted
//! sequence. That count is asserted at zero.
//!
//! # Corpus
//!
//! The corpus is a generated artifact (gitignored, like the spec fixtures) and
//! is regenerated on demand through `common::fixture_gen`, so a fresh checkout
//! just works. Its parameters are what the pins are measured over: regenerating
//! with different `--seeds` re-rolls every number here.
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

use std::collections::BTreeSet;
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

/// The 3'-direction census, pinned. See the module docs: every divergence
/// figure is a baseline that must only ever go down, and `converged` a floor
/// that must only ever go up.
///
/// **100.0% of designed classes converge** (11 271 of 11 272). The 1 that does
/// not is the measurement this module exists to make: nothing in the repo could
/// state that number before, because only 650 real rows reach the partitioner at
/// all — and `multi_member_cis_axis` measures convergence over the respellable
/// ones of those, which is what a corpus of far-apart members looks like.
///
/// **Read the percentage as rounded, not as exact.** 11 271 of 11 272 is
/// 99.991 %, and the guard below rounds to the nearest tenth, so the headline
/// says 100.0 while one class still diverges. The residual is named in the
/// re-bless section below; a reader who needs "all of them" should read the
/// count, which is the figure the guard actually asserts.
///
/// # This headline was wrong for four movements of the pin
///
/// It read **"58.8% … (6 633 of 11 272). The 4 639 that do not"** until
/// 2026-08-11. That is not the figure as *first* landed: it landed in #1458
/// reading 6 632 / 4 640 and was carried forward exactly once, by #1493 — the
/// change that closed #1454 — and never again. So #1454 is the one movement it
/// did track, and the ones it missed are #1524, the axis gate, #1539 and #1599.
/// Each of those four is recorded in the chronology below, which terminates
/// correctly at the pinned figure; only the headline was left behind. The one
/// number a reader takes away was therefore the one number in the file that was
/// wrong, and it was wrong by 1 394 classes — larger than any single movement the
/// chronology records.
///
/// Both counts here were re-derived from the pin rather than carried across the
/// rebase onto #1599: the movement to 8 027 landed while this correction was in
/// review, which turned "three movements" into four and 1 393 into 1 394. That is
/// the drift this section is about, arriving once more during the change that
/// documents it — and it is why the numbers below are asserted rather than
/// written down.
///
/// It arrived a **third** time in #1649, which moved the pin to 8 361 while that
/// very paragraph was on the branch. The headline moved with it because
/// [`the_headline_matches_the_pin`] made it impossible not to — which is the
/// first movement of this pin the prose has tracked without a human noticing it
/// had to, and the argument for the guard stated as an outcome rather than as an
/// intention.
///
/// Restating a pinned constant in prose is the drift this repository names as
/// its recurring failure mode, and prose cannot be derived from a `const`. The
/// restatement is therefore *checked* instead:
/// [`the_headline_matches_the_pin`] reads this file's own source and fails when
/// the three figures here disagree with [`THREE_PRIME`]. A future re-bless that
/// forgets the prose now reddens rather than silently ageing.
///
/// Note the entries below are **historical and correctly labelled** — in
/// particular the axis gate's "6 633 -> 8 006 … by 1 373", which the three
/// sections after it correct to 8 003, then 8 026, then 8 027. They are a
/// changelog, not competing claims, and an earlier review of this file mistook
/// one of them for a second stale figure.
///
/// Note what the remainder is *not*: `declined` and `sequence_changed` are both
/// zero, so every divergence is two well-formed, sequence-preserving outputs
/// for one variant. That is the confluence defect itself, not a parse failure
/// or a corrupted sequence wearing its clothes.
///
/// `converged` was 6 632 before #1454; see the module docs for which single
/// class moved and why the other #1454 class could not.
///
/// # `converged` went **down** by 4 in #1524, deliberately
///
/// 6 633 -> 6 629, `split 2` 4 505 -> 4 509. This is the direction this pin
/// exists to catch, so it is stated rather than absorbed. The four classes are
/// `s00-c-m3-sep0-p8-rot4`, `s00-c-m4-sep0-p8-rot4`, `s02-c-m4-sep1-p4-rot1`
/// and `s03-c-m3-sep1-p8-rot2` (5' loses two of the same family).
///
/// In each, one multi-member spelling used to reach the lone-`delins` form only
/// **via** the defect #1524 removes: `Normalizer::build_split_variants` re-split
/// an intermediate member into two *touching* members, the legacy per-member
/// merge then re-merged them into a wider `delins`, and that wider form is what
/// the other spellings had settled on. Take the illegal intermediate away and
/// the chain stops one step earlier. Traced end to end on
/// `s00-c-m3-sep0-p8-rot4`: `c.[7_8insTTAATATA;17_24del;24_25insCCAACCCC]` now
/// gives `c.[18_20delinsAAT;24_25insCCAACCCC]` where the other three spellings
/// give `c.18_24delinsAATATATCCAACCCC`.
///
/// Two things make this an acceptable price rather than a hidden regression.
/// `sequence_changed` is still 0, so no class converged on the wrong sequence;
/// and the newly-divergent output is the *better* of the two forms — the
/// single `delins` spans three unchanged nucleotides (`c.21_23`) that
/// `general.md:34` says to describe individually, while the two-member form
/// does not. What is lost is agreement, not correctness, and the disagreement
/// is between a lone `delins` and a multi-member allele — the exact pair
/// #1235's axis-gate widening is about, since `is_splittable_single_member`
/// still refuses to re-derive a lone `c.` delins from its sequence.
///
/// Measured against 58 corpus rows fixed. The alternative shape of the fix —
/// joining the codon triplet to a touching run on both edges instead of
/// declining the exception on the left — cost **44** classes here rather than
/// 4, and was rejected for it.
/// **Re-blessed when the axis gate opened to `c.`/`n.`/`r.`.** Widening
/// `is_splittable_single_member` lets a lone transcript-axis member reach
/// `sequence_first_pass`, which is the entire point of that change: converged
/// rises **6 633 -> 8 006** and `split_two` falls by the same **1 373**, while
/// `split_three` (115) and `split_more` (19) are **unchanged**. Every
/// divergence figure moved down or stayed flat, which is the direction this pin
/// demands.
///
/// The figure is now identical with and without `FERRO_ASSERT_IDEMPOTENT`
/// (measured in both, `declined: 0` in both). That is #1493's doing rather than
/// this branch's — it closed #1454, so the two classes that used to panic the
/// oracle no longer do, and the per-configuration pins this file once carried
/// were already collapsed on `main` before this change landed.
///
/// **The pinned figure is 8 003, not 8 006.** Both numbers above are true of
/// the change that produced them, and neither is the number to pin: the axis
/// gate raises converged by 1 373 against a `main` that already carries
/// #1537's deliberate reduction of 4 (3') and 2 (5'). The composition was
/// measured rather than arithmetic — 8 006 is what this branch scored against
/// a `main` predating #1537, and re-running after the rebase gives 8 003 in
/// **both** directions. That the two directions land on the same number, as
/// they did before the rebase, is the reading to trust.
///
/// # Raised again by the #1539 member audit: 8 003 -> 8 026
///
/// `split_two` falls by the same 23 and `split_three` (115) and `split_more`
/// (19) are unchanged, so every class that moved moved from two outputs to one.
/// `split_concealed_separations` cuts a member that conceals a separation
/// `general.md:34` requires, and the classes this converges are ones where the
/// lone-`delins` spelling and the multi-member spelling previously reached the
/// concealed form and the individual form respectively.
///
/// `declined`, `underdetermined` and `sequence_changed` are all still 0, and
/// that is the load-bearing part rather than a formality: an earlier revision of
/// that pass dropped payload bases when it cut a run of two matched columns with
/// an insertion between them, and this census is what reported it — 45 declined
/// and 10 underdetermined, with `converged` *rising* at the same time. A pass
/// that corrupts a sequence can raise the convergence figure, so the three zeros
/// have to be read before the headline.
///
/// # Raised again by the amino-acid precondition: 8 026 -> 8 027
///
/// `coalesce_coding_frame_separation` now tests the second conjunct of
/// `DNA/delins.md:18` — "together affecting one amino acid" — instead of the
/// `length_changing` proxy, so it declines a merge whose span crosses a codon
/// boundary. `split_two` falls by the same 1 and `split_three`/`split_more` are
/// unchanged, so the one class that moved moved from two outputs to one, and the
/// three zeros are still zero.
///
/// One is a small number and worth saying plainly: the change's value is
/// conformance, not this census.
///
/// **And the reading it invites is wrong — scoped, 2026-08-11.** This paragraph
/// used to end "What the census rules out is the opposite reading — that
/// restoring the precondition *costs* confluence. It does not, in either
/// direction." That is true **of this corpus** and false in general, so it must
/// not be read as a claim about the change. `spec_conformance_axis`'s corpus is
/// deliberately harsher — it varies member geometry, transcript geometry and
/// scale rather than holding them fixed — and there the 3' direction **loses**
/// one: `converged` 9,140 -> 9,139, six classes losing convergence against five
/// gaining it. Both censuses are correct about their own corpus; neither
/// generalises to the other, which is the whole reason two of them exist. The
/// promoted rows are in
/// `spec_corpus_regressions::the_codon_gate_splits_a_spanning_delins_its_own_members_do_not`.
/// **Raised again by #1649's two-deletion alignment: 8 027 -> 8 361**, the
/// largest movement this pin has taken since the axis gate. `merge`'s payload
/// splitter could express *insertion, retained reference, insertion* but not the
/// mirrored *deletion, retained reference, deletion*, so a payload whose
/// alignment against its span carries two gaps was forced onto a coarser
/// partition than the sequence supports — and which coarser partition it landed
/// on depended on the spelling it started from, which is what a confluence
/// failure is. `split_two` falls by 276, `split_three` by 48 and `split_more` by
/// 10; the three sum to the 334 `converged` gains, so **every** moved class went
/// straight to convergence and none merely dropped an arity. No divergence
/// figure rises, so the ratchet is satisfied in the direction it was written for.
///
/// # Re-blessed by #1716, and this one costs a class rather than buying one
///
/// `merge_consecutive_edits`' per-member codon-frame predicate tested its left
/// anchor's **right edge** for codon membership while authorising a `delins` over
/// `prev_a.start ..= next.end`, so a left member wider than one base merged
/// across a codon boundary — `general.md:35`'s "together affecting one amino
/// acid" cannot hold over four or more positions. Testing the span instead
/// declines those merges.
///
/// `converged` is **unchanged at 8 361 / 8 367**, and exactly one class moves, in
/// both directions: `s03-c-m3-sep1-p8-rot2` goes from two outputs to three, so
/// `split_two` falls by one and `split_three` rises by one. The full divergence
/// list was dumped either side (the `worst` cap raised locally, then restored)
/// and diffed by class id: **no class loses convergence, none gains it, and no
/// other class changes arity** in either direction.
///
/// **`split_three` rising is a divergence figure going UP**, which this pin's
/// ratchet forbids without a re-bless — hence this section. The price is the same
/// one #1240 records two paragraphs above, and for the same reason: the
/// newly-divergent third output is `c.[9_16delinsATGAGTAA;18delinsGTCCTGGCA]`,
/// and what it stopped agreeing with is `c.9_18delinsATGAGTAAAGTCCTGGCA` — a
/// merge across `c.9..c.18`, i.e. codons 3 to 6, which `general.md:34` says to
/// describe individually. Agreement is lost; correctness is not. `declined`,
/// `underdetermined` and `sequence_changed` are all still 0.
/// # Re-blessed by the partition default flip (#1835), and this is the largest
/// move this pin has ever taken
///
/// `converged` **8 361 -> 11 271** and every divergence figure falls to its
/// floor: `split_two` 2 834 -> 1, `split_three` 68 -> 0, `split_more` 9 -> 0.
/// Divergence is 2 911 classes under `live` and **1** under the new default, so
/// the ratchet this pin is written for — divergence only ever down, `converged`
/// only ever up — is satisfied in the strongest direction it admits. `declined`,
/// `underdetermined` and `sequence_changed` are all still **0**, which is the
/// triple read first: no class converged onto a different sequence, so this is
/// agreement being gained and not correctness being spent.
///
/// **The figure is not derived here — it is the one a decided operator ruling
/// already published.** `delins-recommendation-reach-when-the-input-arrives-split`
/// records, in its own CONFLUENCE IS NOT A GROUND paragraph, that over this very
/// 11 272-class corpus "divergent classes fall from 2 911 under `live` to 1
/// under `canonical-coalesced`, the spanning-vs-split family (1 756 classes)
/// closes entirely, and cross-class disagreements go 16 -> 0", and it attributes
/// those figures to "the PARTITION DEFAULT, not to this ruling". This change is
/// that default. So the re-bless is the ledger's own measurement arriving, not a
/// new claim.
///
/// **Why converging thousands of classes is not itself the argument.** That same
/// paragraph is explicit that confluence is *not* a ground for the ruling:
/// `coalesce_payload_alignment_split` runs last on pieces already derived from
/// the sequence, so each equivalence class reaches one answer whichever way the
/// arm is set. What the arm decides is **which form** that answer takes, and the
/// forms are licensed one family at a time by the records this PR cites — the
/// `c.`-axis payload-coincidence scope, the gap-bearing-insert scope, and
/// `duplication-must-ranks-the-label-not-the-partition`. Reading the convergence
/// gain as the licence would be reading a count as a property, which is the one
/// move this file's own module docs forbid.
///
/// **The single surviving 3' divergence is `s06-g-m4-sep1-p1-all-dup`**, which
/// reaches `g.[265dup;267_271A[7];271_272insA]` and `g.[265dup;267_271A[8]]` —
/// a repeat-notation disagreement over whether a trailing copy is spelled into
/// the tract count or beside it. Neither spelling is a partition this change
/// touches; `delins-adjacent-members-when-both-consume-reference` already
/// records the adjacent repeat-expansion shape as outside its reach, and the
/// same holds here. The 5' direction has no residue at all.
const THREE_PRIME: Census = Census {
    classes: 11_272,
    spellings: 47_392,
    declined: 0,
    // #1835: the partition default flip — see the section above.
    converged: 11_271,
    split_two: 1,
    split_three: 0,
    split_more: 0,
    underdetermined: 0,
    sequence_changed: 0,
};

/// The 5'-direction census, pinned. Confluence is a property of the normalizer
/// rather than of one shuffle direction, so it is measured in full rather than
/// spot-checked.
///
/// The 5' direction is **no longer a public option** — see `README.md` rule 6
/// and `tests/it/five_prime_public_surface_removed.rs`. That changes the reason
/// this census exists, not its scope: the 5' arm is ferro's differential oracle over its own 3' output — the instrument that found #1542, where 7 of 8 `FERRO_PARTITION` x direction configurations agreed and only the shipped `live`/3' arm diverged. An arm that is only spot-checked
/// cannot serve as that oracle.
///
/// It landed within five classes of the 3' figure (6 628 against 6 633), which
/// is worth reading as evidence: the divergence is a property of how the
/// partitioner splits a block, not of which end of an ambiguous run the shuffle
/// walks to. A fix that moved only one of these two numbers would be treating a
/// symptom.
///
/// Lowered by **two** in #1524 (6 628 -> 6 626, `split 2` 4 530 -> 4 532), for
/// the reason `THREE_PRIME` records at length. That the two directions move by
/// different amounts is itself the expected reading of the note above: the
/// affected classes are ones whose chain to a common form ran through an
/// intermediate that only some shuffle directions produce.
/// **Re-blessed with the axis gate**, and the two directions now agree exactly
/// (both 8 006) rather than within five. The 5' figure moves by **1 378**
/// against the 3' direction's 1 373 — five more, because the same change fixes
/// an off-by-one in `enclosing_exon` that let a member sitting on an exon's
/// first base escape the window clamp and shuffle across the junction. That is
/// a 5'-only symptom, since the 3' walk moves away from the exon start. Both
/// directions ending on the same number is the reading to trust: the residual
/// is a property of the partitioner, not of a shuffle direction.
///
/// **Raised again by the #1539 member audit: 8 003 -> 8 021**, with `split_two`
/// falling by the same 18 and `split_three`/`split_more` unchanged. The 3'
/// direction gains 23 rather than 18, which is the expected asymmetry: the audit
/// runs before the shift, so a member it cuts is then 3'- or 5'-shifted on its
/// own, and only some of those landings coincide with the other spelling's.
///
/// **Raised again by the amino-acid precondition: 8 021 -> 8 023**, with
/// `split_two` falling by the same 2 and `split_three`/`split_more` unchanged.
/// The 3' direction gains 1 rather than 2 — the same asymmetry the paragraph
/// above records, and for the same reason: the precondition is tested before the
/// shift, so a pair it leaves split is then shuffled independently and only some
/// of those landings coincide with the other spelling's.
///
/// See `THREE_PRIME` for why the three zeros are read first — and for why the
/// entries above, including "both 8 006" and the 8 003 and 8 021 the last two
/// start from, are historical and correctly labelled rather than stale figures.
/// An earlier review of this file mistook the 3' side's equivalent entry for a
/// second stale number; this side carries the same shape.
/// **Raised again by #1649's two-deletion alignment: 8 023 -> 8 367**, with
/// `split_two` down 309, `split_three` down 31 and `split_more` down 4 — again
/// summing exactly to the 344 `converged` gains, so every moved class converged
/// outright. The 5' direction gains 344 against the 3' direction's 334, the same
/// ten-class asymmetry the entries above record and for the same reason: the
/// finer partition is chosen before the shift, so the two members are then
/// shuffled independently and only some of those landings coincide with the
/// other spelling's.
///
/// **Re-blessed by #1716 with `converged` unchanged at 8 367**, `split_two` down
/// one and `split_three` up one — the same single class (`s03-c-m3-sep1-p8-rot2`)
/// as at 3', and no asymmetry this time because the merge the fix declines is
/// decided before either shuffle. See `THREE_PRIME` for the class's two forms and
/// for why losing this agreement is not losing correctness.
/// # Re-blessed by the partition default flip (#1835), to total convergence
///
/// `converged` **8 367 -> 11 272**, i.e. every class, with all three divergence
/// figures at **0**. See `THREE_PRIME` for the licensing, which is shared: the
/// ledger record `delins-recommendation-reach-when-the-input-arrives-split`
/// published these figures against the partition default before this change
/// existed, and the per-family licences are the records this PR cites rather
/// than the convergence gain itself.
///
/// **The asymmetry finally disappears, and that is the reading to trust.** Every
/// entry above records the two directions moving by slightly different amounts,
/// always for the same reason: the partition is chosen before the shift, so two
/// members are then shuffled independently and only some landings coincide. Here
/// the 5' direction converges completely while the 3' keeps one class. That is
/// the same mechanism read from the other end — the residue is
/// `s06-g-m4-sep1-p1-all-dup`, a repeat-notation class whose two spellings differ
/// in which end of an ambiguous tract a trailing copy lands on, so it is visible
/// only to the direction that walks toward it.
///
/// `sequence_changed` remains **0** in both directions, so nothing converged onto
/// different bases.
const FIVE_PRIME: Census = Census {
    classes: 11_272,
    spellings: 47_392,
    declined: 0,
    // #1835: the partition default flip — see `THREE_PRIME`.
    converged: 11_272,
    split_two: 0,
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
pub(crate) struct Class {
    pub(crate) id: String,
    pub(crate) axis: String,
    pub(crate) core: String,
    denoted: String,
    members: usize,
    separation: usize,
    pub(crate) spellings: Vec<String>,
}

#[derive(Deserialize)]
pub(crate) struct Corpus {
    pub(crate) classes: Vec<Class>,
}

fn corpus_path() -> PathBuf {
    fixture_gen::fixture_path(CORPUS_RELATIVE_PATH)
}

pub(crate) fn corpus() -> &'static Corpus {
    static CORPUS: OnceLock<Corpus> = OnceLock::new();
    CORPUS.get_or_init(|| {
        let path = corpus_path();
        fixture_gen::ensure_generated_example_fixture(
            &path,
            "generate_cis_confluence_corpus",
            // The generator's own default is `g,c`, but naming it here keeps the
            // corpus these pins are measured over independent of that default:
            // widening it would otherwise silently re-roll every number below.
            &["--axes", "g,c"],
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
pub(crate) fn reference_for(class: &Class) -> (MockProvider, String) {
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
}

/// A worst-offender line, for the failure message when a pin moves.
struct Divergence {
    id: String,
    outputs: Vec<String>,
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

/// How many divergent classes the failure message names. Applied per group and
/// again to the concatenation, which is what keeps the parallel result
/// byte-identical to a serial one — see [`measure`].
const WORST_LIMIT: usize = 10;

/// The contiguous runs of `classes` that share one reference.
///
/// Classes are sorted by an id that starts with the core index and the axis, so
/// every class sharing a reference is contiguous. Building one provider and one
/// normalizer per group instead of per class is the difference between this
/// running in seconds and in minutes; [`measure`] additionally runs the groups
/// concurrently, which it can do only because each one owns its provider.
fn reference_groups(classes: &[Class]) -> Vec<&[Class]> {
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
        groups.push(&classes[start..end]);
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
        let expected = expected_sequence(class);
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

/// Census one direction over both axes, one reference group per rayon task.
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
fn measure(direction: ShuffleDirection) -> (Census, Vec<Divergence>) {
    let groups = reference_groups(&corpus().classes);

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

fn report(direction: &str, census: &Census, worst: &[Divergence]) -> String {
    let mut out = format!(
        "cis confluence ({direction}): {} classes, {} spellings, {} declined\n  \
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

/// Assert one direction's census against its pin, printing the measured numbers
/// either way so a moved pin can be re-blessed from the test output.
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
}

#[test]
fn three_prime_confluence_census() {
    assert_census(ShuffleDirection::ThreePrime, "3prime", &THREE_PRIME);
}

#[test]
fn five_prime_confluence_census() {
    assert_census(ShuffleDirection::FivePrime, "5prime", &FIVE_PRIME);
}

/// [`THREE_PRIME`]'s headline restates three figures that live in the `const`
/// itself, and prose cannot be derived from a constant — so it is checked here.
///
/// This exists because the restatement **did** drift, silently, through four
/// movements of the pin (#1524, the axis gate, #1539, #1599). Every one of those
/// was recorded in the chronology directly below the headline, and none of them
/// updated the headline, so the file's most-read sentence was 1 394 classes out
/// while everything under it was right. It was carried forward exactly once, by
/// #1493, and never again. A census re-bless is exactly the moment this is
/// easiest to forget — #1599 moved the pin again while this very correction was
/// in review — which is why the guard belongs next to the pin rather than in a
/// docs-audit test somewhere else.
///
/// It guards the headline and nothing else, and that is the whole class rather
/// than a first instalment: every other figure in this file's prose sits inside a
/// dated chronology entry and is a claim about *then*, so it stays true as the
/// pin moves. That was established by reading them, not assumed — see the note on
/// [`THREE_PRIME`] about the entry an earlier review misread as a second stale
/// figure.
///
/// It reads this file's own source rather than a doc-comment reflection API,
/// because none exists: `include_str!` on the module's own path is the only way
/// to get the rendered prose at runtime. That makes the test sensitive to the
/// headline's *wording* — deliberately. If someone rewrites the sentence, they
/// should have to look at this test and confirm the numbers are still asserted,
/// rather than have a silent regex quietly stop matching and pass.
#[test]
fn the_headline_matches_the_pin() {
    const MARKER: &str = "% of designed classes converge**";
    let source = include_str!("cis_confluence_axis.rs");

    let headline = source
        .lines()
        // Doc-comment lines only. `MARKER` is a literal in this very function, so
        // an unrestricted search always matches *something* — the line above —
        // which makes the panic below unreachable and turns a reworded headline
        // into `got [] from "const MARKER: …"` on the count assertion instead of
        // the instruction written for exactly that case.
        .find(|line| line.trim_start().starts_with("///") && line.contains(MARKER))
        .unwrap_or_else(|| {
            panic!(
                "no line contains `{MARKER}` — the headline was reworded. Re-point this \
                 test at the new wording and keep asserting the three figures against \
                 THREE_PRIME; do not delete it, or the prose goes back to ageing silently"
            )
        });

    // The headline spells thousands with ASCII spaces (`8 027`), so digits are
    // gathered across a single space and joined before parsing.
    // `usize`, matching the `Census` fields these are compared against.
    let numbers: Vec<usize> = {
        let flattened = headline.replace(' ', "");
        let mut found = Vec::new();
        let mut digits = String::new();
        for character in flattened.chars() {
            if character.is_ascii_digit() {
                digits.push(character);
            } else if !digits.is_empty() {
                // `71.2` contributes `71` and `2`; the fractional part is
                // re-joined below rather than dropped, so a wrong tenth fails.
                found.push(digits.parse::<usize>().expect("digits parse"));
                digits.clear();
            }
        }
        if !digits.is_empty() {
            found.push(digits.parse::<usize>().expect("digits parse"));
        }
        found
    };

    // `71`, `2`, `8027`, `11272`, `3245` — the percentage's two parts, then the
    // three counts, in the order the sentence states them.
    assert_eq!(
        numbers.len(),
        5,
        "the headline should state a percentage and three counts, got {numbers:?} from \
         {headline:?}"
    );

    let tenths = numbers[0] * 10 + numbers[1];
    // Rounded to the nearest tenth rather than truncated, because the headline is
    // a percentage a human reads and a reader rounds. Truncation makes no
    // difference at the current pin (8 027 / 11 272 = 71.212 %), but at some
    // future pin whose fraction lands at or above .x5 it would oblige the prose
    // to state a tenth that is not the correctly-rounded one — the guard
    // dictating a wrong-looking number, which is the opposite of its job.
    let expected_tenths =
        (THREE_PRIME.converged * 1000 + THREE_PRIME.classes / 2) / THREE_PRIME.classes;
    assert_eq!(
        tenths,
        expected_tenths,
        "the headline says {}.{}% but THREE_PRIME is {} of {} = {}.{}%",
        numbers[0],
        numbers[1],
        THREE_PRIME.converged,
        THREE_PRIME.classes,
        expected_tenths / 10,
        expected_tenths % 10,
    );
    assert_eq!(
        numbers[2], THREE_PRIME.converged,
        "the headline's converged count disagrees with the pin"
    );
    assert_eq!(
        numbers[3], THREE_PRIME.classes,
        "the headline's class count disagrees with the pin"
    );
    assert_eq!(
        numbers[4],
        THREE_PRIME.classes - THREE_PRIME.converged,
        "the headline's non-converging count is not `classes - converged`"
    );
}
