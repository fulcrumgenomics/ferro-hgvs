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
//! EQUAL    NC_TEST.1:g.264_266delinsCTC -> NC_TEST.1:g.[264A>C;266G>C]
//! UNEQUAL  NC_TEST.1:g.264_266delinsTC  -> NC_TEST.1:g.[264del;266G>C]
//! ```
//!
//! Both strings above are copy-pasted from this base's output.
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
//! `CLAUDE.md` has the RFC 2119 census); what makes it decisive is that it is the
//! only clause whose preconditions this shape meets. `delins.md:16` cannot apply
//! (it needs consecutive changed positions) and the codon exception at
//! `delins.md:18` / `general.md:35` cannot reach a `g.` description at all, since
//! a genomic description has no reading frame. The decided ruling below does not
//! reach this arm either: at equal length the interior column is retained by
//! denotation, not by the payload/reference coincidence that ruling is scoped to.
//!
//! **Unequal (3 -> 2): ferro VIOLATES a decided ruling.** There is no
//! column-for-column reading, so nothing in the input forces the `T` at 265 to
//! be retained. It survives into ferro's output *only* because the payload
//! happens to open with a `T` that an aligner can pair it with. That coincidence
//! is the construction `delins.md:46` builds — "parts of the inserted sequence
//! *align* with the reference sequence, giving an alternative description" — and
//! `delins.md:47` answers it: **"The `delins` format is recommended"**.
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
//! # Assert-then-flip
//!
//! The unequal arm's expectation is ferro's **current, wrong** answer. When the
//! ruling is implemented it becomes [`UNEQUAL_RULING_CONFORMANT`] — the input
//! string unchanged, i.e. the spanning `delins`. It is asserted rather than
//! `#[ignore]`d so that the day the behaviour moves, this file has to be read.
//! The flip is an output-moving change and owes the release a
//! `Representation-Change:` trailer.
//!
//! Fully hermetic: a `JsonProvider` (the type the sibling hermetic tests still
//! name by its `MockProvider` alias), no `FERRO_MANIFEST`, no fixtures, no LFS.

use ferro_hgvs::reference::JsonProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

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
/// base. The decided ruling on `delins.md:47` governs, and ferro does not obey.
const UNEQUAL_INPUT: &str = "NC_TEST.1:g.264_266delinsTC";

/// Measured on this base. **Assert-then-flip**: this is the wrong answer, pinned
/// deliberately. See [`UNEQUAL_RULING_CONFORMANT`].
const UNEQUAL_EXPECTED: &str = "NC_TEST.1:g.[264del;266G>C]";

/// The right answer for the unequal arm once
/// `delins-merge-vs-individual-gap-two-or-more` is implemented: the spanning
/// `delins`, i.e. the input unchanged.
///
/// Named rather than left in prose so the flip is a two-line edit and so this
/// string is greppable from the ruling record.
const UNEQUAL_RULING_CONFORMANT: &str = "NC_TEST.1:g.264_266delinsTC";

fn padded_contig() -> String {
    let pad = "ACGT".repeat(PAD_OFFSET / 4);
    format!("{pad}{CORE}{pad}")
}

fn provider() -> JsonProvider {
    let mut provider = JsonProvider::new();
    provider.add_genomic_sequence(CONTIG, padded_contig());
    provider
}

/// The bases of `padded_contig()` over a 1-based inclusive span.
fn bases(start: u64, end: u64) -> String {
    padded_contig()[(start as usize - 1)..(end as usize)].to_string()
}

fn normalize(input: &str) -> String {
    let normalizer = Normalizer::with_config(provider(), NormalizeConfig::default());
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    normalizer
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("normalize {input}: {e}"))
        .to_string()
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
    assert_eq!(
        normalize(EQUAL_INPUT),
        EQUAL_EXPECTED,
        "equal-length {EQUAL_INPUT} is described individually, which is what `delins.md:17` \
         recommends and the only clause of the three that reaches this shape"
    );

    // UNEQUAL — WRONG, pinned deliberately. See "Assert-then-flip" above.
    assert_eq!(
        normalize(UNEQUAL_INPUT),
        UNEQUAL_EXPECTED,
        "unequal-length {UNEQUAL_INPUT} is currently split on a payload/reference coincidence at \
         265, which is the alignment-driven split `delins.md:46` constructs and `delins.md:47` \
         advises against. The decided ruling `delins-merge-vs-individual-gap-two-or-more` says \
         the spanning form wins, i.e. {UNEQUAL_RULING_CONFORMANT}. If this assertion just failed \
         with that string, the defect is FIXED — flip the expectation and declare the \
         representation change"
    );

    // Guard the flip against being confused with a different move: state
    // explicitly that today's output is not yet the conformant one.
    assert_ne!(
        UNEQUAL_EXPECTED, UNEQUAL_RULING_CONFORMANT,
        "the pinned output equals the ruling-conformant form — the assert-then-flip note above \
         is stale and this test now pins correct behaviour"
    );
}
