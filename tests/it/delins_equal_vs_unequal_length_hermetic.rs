//! The equal- vs unequal-length `delins` discriminator, on a synthetic contig,
//! so it runs on every PR.
//!
//! # Read this together with its manifest-backed twin
//!
//! This file is the hermetic companion to
//! `tests/it/delins_equal_vs_unequal_length_discriminator.rs`
//! (`it::delins_equal_vs_unequal_length_discriminator::equal_and_unequal_length_delins_get_opposite_verdicts`).
//! That file pins two **real** loci — BRCA2 `NC_000013.11:g.32340866_32340868delinsATC`
//! and MSH2 `NC_000002.11:g.47639670_47639673delinsTT` — and carries the full
//! spec argument for both. Read it first; this module does not repeat the
//! argument, only the mechanism.
//!
//! The two files are deliberately not redundant:
//!
//! * the twin needs `FERRO_MANIFEST`, so it is the **nightly's reality check**
//!   that the synthetic shape below still matches biology;
//! * this file needs nothing, so it is the **PR gate**.
//!
//! # Why a hermetic copy exists at all
//!
//! PR CI provisions no prepared reference — `.github/workflows/ci.yml` says so
//! and emits a `::notice` about it — so the twin's single test skips green on
//! every pull request. The question it guards is the most re-litigated one in
//! this project (`delins.md:17` against `:47`, argued at least four times), and
//! a guard that cannot fire on the change that would break it is not a guard.
//! Skip-on-absent is one of the mechanisms this campaign exists to remove, so
//! adding one while removing others would be a poor trade. Hence: same
//! discriminator, no reference.
//!
//! # The construction
//!
//! One synthetic contig, one span, two payloads. Reference `A,T,G` at
//! 264..266 (see [`SPAN_BASES`]); the flanking bases are `G` at 263 and `C` at
//! 267, so neither arm's members can 3'-shift and what is measured is the
//! partition rather than the shuffle.
//!
//! ```text
//! EQUAL    NC_TEST.1:g.264_266delinsCTC  -> NC_TEST.1:g.[264A>C;266G>C]
//! UNEQUAL  NC_TEST.1:g.264_266delinsTC   -> NC_TEST.1:g.[264del;266G>C]
//! ```
//!
//! Both strings above are copy-pasted from this base's output.
//!
//! A third, separate arm — [`MERGE_INPUT`], on its own contig — demonstrates
//! the thing neither of the two above can: since #2155, `:47`'s carve-out
//! reaches this axis too, and a block whose derived split is a placed gap plus
//! a `delins`-rendering member (not [`UNEQUAL_INPUT`]'s plain substitution) now
//! collapses to the spanning form on `g.`, exactly as it already did on `c.`.
//! See [`MERGE_INPUT`]'s doc comment for the construction and the module
//! section "The unequal arm still does not merge" for why the two arms above
//! land differently even though both are net-deletion, single-unchanged-base
//! blocks.
//!
//! ## The two arms differ only in payload length. Asserted, not asserted-in-prose.
//!
//! This is the entire point of the file, and it is the thing the two real loci
//! cannot show: they sit on different chromosomes, at different spans, with
//! different reference bases, so a reader can always tell themselves the
//! difference in verdict came from something else.
//!
//! Here it cannot have. The two inputs share accession, span, and therefore
//! reference bases; the unequal payload is the equal payload **with its first
//! base dropped**. Net length goes `0` to `-1` and nothing else moves.
//! [`the_two_arms_differ_only_in_payload_length`] asserts exactly that, by
//! reconstructing one input from the other rather than by comparing two
//! hand-written strings.
//!
//! ## Why the verdicts differ, in one paragraph each
//!
//! **Equal (3 -> 3): ferro follows `delins.md:17`.** The description denotes a
//! column-for-column reading, `A,T,G` becomes `C,T,C`, so denotation *forces*
//! the interior `T` to be unchanged. The changed columns are 264 and 266,
//! separated by one unchanged column, and `delins.md:17` — "two variants
//! separated by one or more nucleotides should be described individually and
//! **not** as a `delins`" — recommends the individual description. It is
//! lowercase prose and so requires nothing read strictly (the repository
//! `docs/READING_THE_SPEC.md` has the RFC 2119 census); what makes it decisive is that it is the
//! only clause whose preconditions this shape meets. `delins.md:16` cannot apply
//! (it needs consecutive changed positions) and the codon exception at
//! `delins.md:18` / `general.md:35` cannot reach a `g.` description at all, since
//! a genomic description has no reading frame. The decided ruling below does not
//! reach this arm either: at equal length the interior column is retained by
//! denotation, not by the payload/reference coincidence that ruling is scoped to.
//!
//! **Unequal (3 -> 2): ferro's split is still CONFORMANT, but no longer because
//! `:47` is off this axis — it isn't, any more.** There is no column-for-column
//! reading, so nothing in the input forces the `T` at 265 to be retained. It
//! survives into ferro's output *only* because the payload happens to open with
//! a `T` that an aligner can pair it with. That coincidence is the construction
//! `delins.md:46` builds — "parts of the inserted sequence *align* with the
//! reference sequence, giving an alternative description" — and `delins.md:47`
//! answers it: **"The `delins` format is recommended"**.
//!
//! Read the section below ("The unequal arm still does not merge") before
//! reading the record quotations that follow as an indictment of this arm — the
//! axis half of that indictment was retired by #2155, and what is left is a
//! *type* argument, not an axis one.
//!
//! The governing record is the ruling `delins-merge-vs-individual-gap-two-or-more`
//! in `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`, status
//! **decided** (operator ruling, 2026-08-07). Its scope paragraph is what this
//! arm falls under:
//!
//! > This record settles ONLY the case `:44-47` describes: a MINIMAL single
//! > `delins` that would be split because payload bases coincide with reference
//! > bases — the alignment-driven split `:46` constructs and `:47` advises
//! > against. It is NOT a general licence to merge changes separated by two or
//! > more nucleotides.
//!
//! ## Minimality is asserted, because the ruling is scoped to a minimal delins
//!
//! "Minimal" here means the reference span and the payload share neither a
//! prefix nor a suffix — otherwise the description could be trimmed and the
//! ruling would be about a different string. Both arms are minimal (`A` vs `C`
//! and `A` vs `T` at the front; `G` vs `C` at the back), and
//! [`both_arms_are_minimal`] checks it against the provider's own bases rather
//! than trusting this paragraph.
//!
//! ## One boundary caveat, carried over from the twin
//!
//! The ruling's **id and question** are framed at a separation of "two or more"
//! nucleotides, while this arm's separation is **one** (the single unchanged
//! column 265). So the row matches the ruling's *scope paragraph* — an
//! alignment-driven split off a minimal `delins` — while sitting outside its
//! *title*. The scope paragraph is the operative text; anyone acting on this row
//! should say which reading they are using.
//!
//! # The assert-then-flip is WITHDRAWN — the unequal arm was never going to merge
//!
//! This section used to say the unequal arm's expectation was ferro's "current,
//! **wrong**" answer, and that when `delins-merge-vs-individual-gap-two-or-more`
//! was implemented the row would become the spanning
//! `NC_TEST.1:g.264_266delinsTC`. **That instruction is wrong and must not be
//! followed.** It was written before the axis scope was decided.
//!
//! `delins-payload-coincidence-carve-out-is-coding-dna-scoped` (decided) scopes
//! `DNA/delins.md:47` — the clause that recommends the spanning form over a
//! payload-coincidence split — to the **coding DNA axis**. Both arms here are
//! `NC_TEST.1:g.`, a genomic reference. So `:47` does not reach them, and
//! `general.md:34` governs unqualified: two variants separated by an unchanged
//! nucleotide are described individually. The unequal arm's split is therefore
//! **conformant**, and it always was; what was wrong was the note, not the
//! output.
//!
//! Note the direction scope points the same way about *which* record applies but
//! not about the outcome: this row is a net deletion (3 nt span, 2 nt payload),
//! so `delins-merge-vs-individual-gap-two-or-more` would reach it by direction —
//! the AXIS scope is what excludes it, and the axis scope is the later and
//! narrower of the two. Read them in that order.
//!
//! The sibling file `delins_equal_vs_unequal_length_discriminator.rs` withdrew
//! the identical instruction for its own `g.` row under #1835, on this reasoning.
//! The withdrawal was applied to one of the two siblings and not the other; this
//! is the other.
//!
//! # CORRECTED 2026-08-17 (#2155): the axis half of the paragraph above is now
//! # wrong. The unequal arm still does not merge — read on for why.
//!
//! `delins-payload-coincidence-carve-out-is-coding-dna-scoped` was **superseded
//! to all DNA axes** on 2026-08-17 (#2155): `:47`'s carve-out now reaches
//! `c./g./m./n.`, and `NC_TEST.1:g.` is one of them. So the sentence above
//! reading "`:47` does not reach them" is no longer true, and if this file is
//! read only that far it now looks like a defect that this arm still splits.
//! It is not one. [`MERGE_INPUT`] below shows the axis genuinely reaches this
//! contig now: same axis, same net-deletion, same single-unchanged-base shape,
//! and it merges.
//!
//! What still keeps *this* arm split is the other half of
//! [`UNEQUAL_SPANNING_FORM`]'s doc comment, which was already on record before
//! this correction and needed no change: the derived split
//! `[264del;266G>C]` has a lone **substitution** as its second member, and
//! `merge.rs`'s `piece_renders_as_delins` excludes exactly that shape (a
//! single reference base replaced by a single payload base is a real,
//! rank-1 type the split buys — see `split_buys_no_higher_priority_type`'s doc
//! comment — not an alignment artefact). The carve-out's own machinery
//! (`split_is_a_placed_gap_coincidence`) requires every non-gap member to
//! render as a `delins`; a plain substitution never does, on any axis. So this
//! arm now declines the merge on **type**, where it used to decline on
//! **axis** — the outcome is unchanged, the reason is not, and a reader citing
//! the axis for it after today would be citing a superseded reason for a
//! still-correct output.
//!
//! # Both arms are also checked against the bases they denote
//!
//! Each arm now asserts that the pinned output and its input denote the same
//! contig sequence, through `cis_apply_oracle` — `hgvs_to_spdi` plus an SPDI
//! splice, with no normalization in the path, so it cannot agree with an output
//! merely because normalization produced it (#1626).
//!
//! Without it a string equality cannot tell a **re-spelling** from a
//! **corruption**: an output describing different bases that happened to render
//! as the pinned string would satisfy it. Both arms pass, so both of today's
//! answers are genuine re-spellings.
//!
//! It matters most on the unequal arm, where the two candidate descriptions
//! differ in **partition** rather than in the bases they denote. Asserting the
//! denotation separately is what lets this file say that the split and the
//! spanning form are two spellings of one variant, so the choice between them is
//! a representation question and never a correctness one.
//!
//! Fully hermetic: a `JsonProvider` (the type the sibling hermetic tests still
//! name by its `MockProvider` alias), no `FERRO_MANIFEST`, no fixtures, no
//! out-of-tree corpora.

use ferro_hgvs::reference::JsonProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

use crate::common::cis_apply_oracle::apply_reason;

/// Padding on each side of [`CORE`], so the measured span sits far from either
/// contig end and no boundary clamp is in play.
///
/// Matches `coding_frame_merge_axis_asymmetry.rs` and `cis_confluence_axis.rs`
/// so the geometry is the one the rest of the suite reasons about.
const PAD_OFFSET: usize = 256;

const CONTIG: &str = "NC_TEST.1";

/// The core the span is drawn from. Shared verbatim with
/// `coding_frame_merge_axis_asymmetry.rs`, chosen there for low local
/// repetitiveness: a homopolymer run would let the 3'-shift dominate and would
/// measure shuffling rather than partitioning.
const CORE: &str = "ACGTACGATGCAAGTCCTGAGGTCAATCGGATCCTAGGACTTCAGGTACCATGAGTCATGCAA";

/// 1-based inclusive genomic span both arms address. `CORE[7..10] == "ATG"`.
const SPAN: (u64, u64) = (264, 266);

/// The reference bases at [`SPAN`]. Asserted against the provider, not quoted.
const SPAN_BASES: &str = "ATG";

/// The bases immediately outside [`SPAN`], asserted so that "neither arm can
/// 3'-shift" stays true if `CORE` is ever edited.
const FLANKS: (char, char) = ('G', 'C');

/// EQUAL arm — 3 nt span, 3 nt payload. `delins.md:17` is the only clause that
/// reaches this shape, and ferro emits what it recommends.
const EQUAL_INPUT: &str = "NC_TEST.1:g.264_266delinsCTC";

/// Measured on this base.
const EQUAL_EXPECTED: &str = "NC_TEST.1:g.[264A>C;266G>C]";

/// UNEQUAL arm — the same span, the same reference bases, one fewer payload
/// base. `delins.md:47` recommends the spanning form for this shape, and since
/// #2155 its carve-out reaches this axis (`g.`) — but the derived split still
/// declines to merge, on **type**: see [`UNEQUAL_SPANNING_FORM`]'s doc comment.
const UNEQUAL_INPUT: &str = "NC_TEST.1:g.264_266delinsTC";

/// Measured on this base, and **conformant** — see the module docs' 2026-08-17
/// correction. This is not an assert-then-flip; there is no flip owed, before
/// or after #2155.
const UNEQUAL_EXPECTED: &str = "NC_TEST.1:g.[264del;266G>C]";

/// The spanning form, named only so the guard below can assert ferro does **not**
/// emit it here.
///
/// It was previously called `UNEQUAL_RULING_CONFORMANT` and documented as the
/// answer a future implementation should reach. That was wrong — the rename is
/// deliberate, so a reader grepping the old name from the ruling record lands on
/// a correction rather than on an instruction to flip.
///
/// # #1610 still does not reach this arm — but only for ONE of its two reasons now
///
/// #1610 implements a neighbouring rule in the partitioner — see
/// `unequal-length-block-a-placed-gap-is-not-a-separation` — so the obvious
/// inference is that this arm should have moved with it. It has not, but this
/// doc comment used to give **two** independent reasons and now only one of
/// them still holds.
///
/// **Still true: the record's fourth condition.** This block's derived split is
/// `[264del;266G>C]`, whose second member is a lone **substitution**. A
/// substitution is a rank-1 type the split genuinely buys
/// (`piece_renders_as_delins` excludes exactly this shape), so
/// `split_is_a_placed_gap_coincidence` declines and the members stay
/// individual. This reason does not depend on the axis at all — it would hold
/// identically on `c.`.
///
/// **No longer true, corrected 2026-08-17 (#2155): the axis half.** This doc
/// comment used to add a second, independent reason — "#1610's rule is gated on
/// `CoincidenceCarveOut::may_disbelieve_a_separation`, so it does not run on a
/// `g.` description at all." `delins-payload-coincidence-carve-out-is-coding-
/// dna-scoped` was superseded to all DNA axes that day, so
/// `may_disbelieve_a_separation` now answers `true` for `g.` too and #1610's
/// rule DOES run on this description. It just still declines, on the type
/// reason alone. See [`MERGE_INPUT`] for a shape where the axis reason used to
/// be doing real work and the merge now goes through.
const UNEQUAL_SPANNING_FORM: &str = "NC_TEST.1:g.264_266delinsTC";

/// A second, independent contig demonstrating what [`UNEQUAL_INPUT`] cannot:
/// the #2155 axis widening actually reaching a `g.` block and merging it.
///
/// Net deletion (4 nt span, 2 nt payload), one unchanged interior base, same
/// shape class as [`UNEQUAL_INPUT`] — but the derived split here is a placed
/// gap (`AG` -> `T`, a `delins`-rendering piece) plus a pure deletion (`G` ->
/// nothing), **not** a plain substitution. Both derived pieces satisfy
/// `split_is_a_placed_gap_coincidence`'s "placed gap or renders-as-`delins`"
/// requirement, so nothing blocks the merge on type, and since #2155 nothing
/// blocks it on axis either. This is the same shape as the real MSH2 locus in
/// `delins_equal_vs_unequal_length_discriminator.rs`
/// (`NC_000002.11:g.47639670_47639673delinsTT`, `AGTG` -> `TT`), reproduced
/// here hermetically.
const MERGE_CONTIG: &str = "NC_TEST2.1";

/// Flanking `C` on each side, chosen only so neither the first nor the last
/// reference base of [`MERGE_SPAN_BASES`] repeats into its flank — `A`/`C` and
/// `G`/`C` both differ, so nothing here can 3'-shift and what is measured is
/// the partition rather than the shuffle. Asserted, not assumed —
/// [`the_merge_fixture_is_the_one_it_was_measured_on`].
const MERGE_PAD: usize = 64;

/// 1-based inclusive genomic span, on [`MERGE_CONTIG`].
const MERGE_SPAN: (u64, u64) = (MERGE_PAD as u64 + 1, MERGE_PAD as u64 + 4);

/// The reference bases at [`MERGE_SPAN`]. Asserted against the provider, not
/// quoted.
const MERGE_SPAN_BASES: &str = "AGTG";

/// Minimal against [`MERGE_SPAN_BASES`]: `A` vs `T` at the front, `G` vs `T` at
/// the back — neither end shares a base, so the description cannot be trimmed
/// and the ruling's minimality precondition holds.
const MERGE_INPUT: &str = "NC_TEST2.1:g.65_68delinsTT";

/// Measured on this base. Identical to [`MERGE_INPUT`] itself: the spanning
/// form of a `delins` whose own span and payload are what it already names, so
/// a normalizer that keeps the block whole reports back the same string it was
/// given (`changed: false` in the CLI's own vocabulary — see
/// `MEMBERS_COALESCED_FROM_REPORTED_FORM`, which fires only when the *input*
/// named more than one member and this one names one throughout).
const MERGE_EXPECTED: &str = "NC_TEST2.1:g.65_68delinsTT";

fn padded_contig() -> String {
    let pad = "ACGT".repeat(PAD_OFFSET / 4);
    format!("{pad}{CORE}{pad}")
}

/// [`MERGE_CONTIG`]'s sequence: [`MERGE_SPAN_BASES`] flanked by non-repeating
/// `C` padding, per [`MERGE_PAD`].
fn merge_padded_contig() -> String {
    let pad = "C".repeat(MERGE_PAD);
    format!("{pad}{MERGE_SPAN_BASES}{pad}")
}

fn provider() -> JsonProvider {
    let mut provider = JsonProvider::new();
    provider.add_genomic_sequence(CONTIG, padded_contig());
    provider.add_genomic_sequence(MERGE_CONTIG, merge_padded_contig());
    provider
}

/// The bases of `padded_contig()` over a 1-based inclusive span.
fn bases(start: u64, end: u64) -> String {
    padded_contig()[(start as usize - 1)..(end as usize)].to_string()
}

/// The bases of `merge_padded_contig()` over a 1-based inclusive span.
fn merge_bases(start: u64, end: u64) -> String {
    merge_padded_contig()[(start as usize - 1)..(end as usize)].to_string()
}

fn normalize(input: &str) -> String {
    let normalizer = Normalizer::with_config(provider(), NormalizeConfig::default());
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    normalizer
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("normalize {input}: {e}"))
        .to_string()
}

/// The contig sequence `description` denotes, reached **without** the
/// normalizer (#1626).
///
/// `cis_apply_oracle` converts each member through `hgvs_to_spdi` and splices
/// the reference, so nothing in this path can agree with an output merely
/// because normalization produced it. A hand-rolled applier here would be a
/// second thing to drift; this is the one the cis sweeps already rest on.
///
/// It panics rather than returning an `Option`. Every description this module
/// feeds it is a two-member allele on disjoint columns of a synthetic contig,
/// so a decline is a defect in the description and not a limit of the oracle —
/// and a silently skipped comparison is the failure mode the whole check exists
/// to remove.
fn denotes(description: &str) -> String {
    apply_reason(&provider(), &padded_contig(), description)
        .unwrap_or_else(|why| panic!("{description} denotes no single sequence: {why:?}"))
}

/// [`denotes`], against [`MERGE_CONTIG`]'s sequence instead of [`CONTIG`]'s.
fn merge_denotes(description: &str) -> String {
    apply_reason(&provider(), &merge_padded_contig(), description)
        .unwrap_or_else(|why| panic!("{description} denotes no single sequence: {why:?}"))
}

/// Assert that `input` normalizes to `expected` **and** that the two denote the
/// same contig sequence.
///
/// The second half is the one that is new. Both arms below are re-spellings —
/// the pinned output is never the input — and a string equality cannot tell a
/// re-spelling from a corruption: an output describing different bases that
/// happened to render as the pinned string would satisfy it. That is not a
/// hypothetical class in this repository; it is what #1615 exists for, and
/// #1592 and #1600 are live instances where a well-formed, in-bounds,
/// idempotent, re-parseable output denotes different bases.
///
/// It matters most on the unequal arm, where the two candidate descriptions
/// differ in partition rather than in the bases they denote, so the choice
/// between them is a representation question and never a correctness one.
fn assert_normalizes_preserving(input: &str, expected: &str, whichever: &str) {
    let actual = normalize(input);
    assert_eq!(actual, expected, "{whichever}");
    assert_eq!(
        denotes(&actual),
        denotes(input),
        "{input} -> {actual} is not a re-spelling: it denotes different bases"
    );
}

/// The payload of a `<accession>:g.<start>_<end>delins<payload>` input.
fn payload_of(input: &str) -> &str {
    input
        .rsplit_once("delins")
        .unwrap_or_else(|| panic!("{input} is not a `delins` description"))
        .1
}

/// The fixture is the one the pinned strings were measured on.
///
/// Without this, an edit to `CORE` or `PAD_OFFSET` silently re-points both arms
/// at different bases and the module docs stop describing what runs.
#[test]
fn the_fixture_is_the_one_the_arms_were_measured_on() {
    assert_eq!(PAD_OFFSET % 4, 0, "the ACGT pad must tile exactly");
    let (start, end) = SPAN;
    assert_eq!(
        bases(start, end),
        SPAN_BASES,
        "{CONTIG}:{start}_{end} must read {SPAN_BASES} for this module's clause arguments to hold"
    );
    assert_eq!(
        (
            bases(start - 1, start - 1).chars().next().unwrap(),
            bases(end + 1, end + 1).chars().next().unwrap(),
        ),
        FLANKS,
        "the flanking bases pin that neither arm's members can 3'-shift; a change here means \
         the outputs below are measuring the shuffle rather than the partition"
    );
    // The interior column is the one whose fate the whole file is about.
    assert_eq!(
        SPAN_BASES.as_bytes()[1] as char,
        'T',
        "the interior base must be the one both payloads carry"
    );
}

/// The fixture [`MERGE_INPUT`]/[`MERGE_EXPECTED`] were measured on.
///
/// Without this, an edit to [`MERGE_SPAN_BASES`] or [`MERGE_PAD`] silently
/// re-points the arm at different bases and the module docs stop describing
/// what runs — the same discipline
/// [`the_fixture_is_the_one_the_arms_were_measured_on`] applies to the other
/// two arms, on the second contig.
#[test]
fn the_merge_fixture_is_the_one_it_was_measured_on() {
    let (start, end) = MERGE_SPAN;
    assert_eq!(
        merge_bases(start, end),
        MERGE_SPAN_BASES,
        "{MERGE_CONTIG}:{start}_{end} must read {MERGE_SPAN_BASES} for this arm's argument to hold"
    );
    assert_ne!(
        merge_bases(start - 1, start - 1).chars().next().unwrap(),
        MERGE_SPAN_BASES.chars().next().unwrap(),
        "the 5' flank must differ from the span's first base, or the block could 3'-shift and \
         this would measure the shuffle rather than the partition"
    );
    assert_ne!(
        merge_bases(end + 1, end + 1).chars().next().unwrap(),
        MERGE_SPAN_BASES.chars().last().unwrap(),
        "the 3' flank must differ from the span's last base, or the block could 3'-shift and \
         this would measure the shuffle rather than the partition"
    );
}

/// The discriminator is **net length and nothing else**.
///
/// The unequal input is reconstructed from the equal one by dropping a single
/// payload base, so the claim is checked rather than asserted in prose. The two
/// real loci in the manifest-backed twin cannot support this claim — different
/// chromosomes, spans and bases — which is the specific gap this file closes.
#[test]
fn the_two_arms_differ_only_in_payload_length() {
    let equal_payload = payload_of(EQUAL_INPUT);
    let unequal_payload = payload_of(UNEQUAL_INPUT);

    // Same accession, same span: everything left of the payload is identical.
    let equal_prefix = &EQUAL_INPUT[..EQUAL_INPUT.len() - equal_payload.len()];
    let unequal_prefix = &UNEQUAL_INPUT[..UNEQUAL_INPUT.len() - unequal_payload.len()];
    assert_eq!(
        equal_prefix, unequal_prefix,
        "the two arms must address the same accession and span, or the verdict difference is \
         attributable to something other than net length"
    );

    // And the payloads differ by exactly one dropped base.
    assert_eq!(
        &equal_payload[1..],
        unequal_payload,
        "the unequal payload must be the equal payload with its first base dropped"
    );

    // Net length: 0 for the equal arm, -1 for the unequal one.
    let span_len = (SPAN.1 - SPAN.0 + 1) as usize;
    assert_eq!(
        equal_payload.len(),
        span_len,
        "equal arm is {span_len}->{span_len}"
    );
    assert_eq!(
        unequal_payload.len(),
        span_len - 1,
        "unequal arm is {span_len}->{}",
        span_len - 1
    );
}

/// Both arms are **minimal**: no shared prefix or suffix between the reference
/// span and the payload.
///
/// The decided ruling is scoped to a *minimal* single `delins`, so an arm that
/// could be trimmed would be a different question. Checked against the
/// provider's bases rather than against [`SPAN_BASES`] alone.
#[test]
fn both_arms_are_minimal() {
    let (start, end) = SPAN;
    let reference = bases(start, end);
    for input in [EQUAL_INPUT, UNEQUAL_INPUT] {
        let payload = payload_of(input);
        assert_ne!(
            reference.as_bytes().first(),
            payload.as_bytes().first(),
            "{input}: reference {reference} and payload {payload} share a prefix, so the \
             description is not minimal and the ruling's scope does not reach it"
        );
        assert_ne!(
            reference.as_bytes().last(),
            payload.as_bytes().last(),
            "{input}: reference {reference} and payload {payload} share a suffix, so the \
             description is not minimal and the ruling's scope does not reach it"
        );
    }
}

/// The discriminator itself: one span, two payloads, opposite verdicts.
///
/// One test rather than two because the point is the *contrast* — read apart,
/// each arm looks like an arbitrary call, and that is how this pair has already
/// been re-argued three times.
///
/// The manifest-backed twin is
/// `it::delins_equal_vs_unequal_length_discriminator::equal_and_unequal_length_delins_get_opposite_verdicts`.
#[test]
fn equal_and_unequal_length_delins_get_opposite_verdicts_hermetically() {
    // EQUAL — ferro follows `delins.md:17`, the only clause whose preconditions
    // this shape meets: the equal-length span denotes an unchanged interior
    // column, and `delins.md:18`'s codon exception cannot reach a `g.`
    // description.
    assert_normalizes_preserving(
        EQUAL_INPUT,
        EQUAL_EXPECTED,
        &format!(
            "equal-length {EQUAL_INPUT} is described individually, which is what \
             `delins.md:17` recommends and the only clause of the three that reaches this shape"
        ),
    );

    // UNEQUAL — still conformant, but see the module docs' 2026-08-17
    // correction before reading this as "off axis": since #2155 the axis
    // reaches `g.` too. What keeps this specific arm split is the derived
    // split's TYPE, not its axis — see [`UNEQUAL_SPANNING_FORM`].
    assert_normalizes_preserving(
        UNEQUAL_INPUT,
        UNEQUAL_EXPECTED,
        &format!(
            "unequal-length {UNEQUAL_INPUT} is split on a payload/reference coincidence at 265, \
             which is the alignment-driven split `delins.md:46` constructs and `delins.md:47` \
             advises against. Since #2155 `:47`'s carve-out reaches this axis (`g.`), so the \
             axis is not what keeps this split; what does is that the derived split's second \
             member, 266G>C, is a lone substitution — a real type the split buys, not payload \
             coincidence — so `split_is_a_placed_gap_coincidence` declines regardless of axis. \
             If this assertion just failed with {UNEQUAL_SPANNING_FORM}, that would mean the \
             carve-out started merging across a real substitution, which is a correctness \
             regression, not a fix: read the module docs before re-pinning"
        ),
    );

    // The split and the spanning form are distinct strings, so the assertion
    // above is testing something. Kept as a guard against the two constants
    // being edited into agreement, which would make this file vacuous.
    assert_ne!(
        UNEQUAL_EXPECTED, UNEQUAL_SPANNING_FORM,
        "the pinned output equals the spanning form, so this file no longer distinguishes the \
         two descriptions it exists to keep apart"
    );
}

/// The #2155 positive control: a placed-gap-plus-`delins` block, which used to
/// stay split on `g.` for axis reasons alone, now collapses to the spanning
/// form — the same collapse [`UNEQUAL_INPUT`] cannot show, because its derived
/// split carries a plain substitution instead. See [`MERGE_INPUT`]'s doc
/// comment for the construction.
#[test]
fn the_widened_carve_out_now_merges_a_placed_gap_shape_on_g() {
    let actual = normalize(MERGE_INPUT);
    assert_eq!(
        actual, MERGE_EXPECTED,
        "{MERGE_INPUT} must collapse to the spanning form: `AG`->`T` renders as a `delins` and \
         `G`->nothing is a placed gap, so every derived piece satisfies \
         `split_is_a_placed_gap_coincidence`'s \"placed gap or renders-as-`delins`\" test, and \
         #2155 widened the axis this carve-out reaches to include `g.`. If this assertion just \
         failed with the split form, the widening has regressed"
    );
    assert_eq!(
        merge_denotes(&actual),
        merge_denotes(MERGE_INPUT),
        "{MERGE_INPUT} -> {actual} is not a re-spelling: it denotes different bases"
    );
}
