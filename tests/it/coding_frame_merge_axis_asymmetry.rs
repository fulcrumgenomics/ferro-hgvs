//! One variant, three axes, two answers: the `c.` axis merges what `n.` and
//! `g.` split.
//!
//! # The measurement
//!
//! Take a single reference `T` and replace it with `ATC`. That is a net **+2**
//! edit, and it has two obvious spellings which denote the same sequence:
//!
//! ```text
//! X:{axis}.{p}delinsATC        and        X:{axis}.[{p-1}_{p}insA;{p}_{p+1}insC]
//! ```
//!
//! Run both spellings, at every position in a 63 nt synthetic transcript whose
//! reference base is a `T`, on three axes over the **same 63 bases** — a coding
//! transcript (`c.`), a non-coding transcript (`n.`), and the genomic contig the
//! transcript sits on (`g.`). Fourteen positions, three axes, two spellings:
//! 84 normalizations. The result is unanimous and it is two facts, not one.
//!
//! **Fact 1 — every axis is self-confluent.** In all 42 pairs the two spellings
//! reach the same output. That is worth a committed guard: it is what PR #1484
//! bought when it widened `Normalizer::is_splittable_single_member` to admit
//! `c.`/`n.`/`r.`, and before that widening a lone transcript-axis `delins`
//! could not be re-derived from its sequence at all.
//!
//! **Fact 2 — the axes disagree with each other, 14 times out of 14.** `n.` and
//! `g.` always split into two members, and always identically (the `g.` output
//! is the `n.` output with every coordinate shifted by the 256 nt of padding).
//! `c.` always keeps the lone `delins`. Same bases, same edit, different string.
//!
//! ```text
//! p=9   c.9delinsATC      n.[8dup;9_10insC]        g.[264dup;265_266insC]
//! p=4   c.4delinsATC      n.[3_4insA;4_5insC]      g.[259_260insA;260_261insC]
//! ```
//!
//! # Where it comes from, verified by counterfactual rather than by reading
//!
//! `merge::canonicalize_from_sequence` partitions the block **identically on all
//! three axes** — two pure-insertion pieces separated by the one unchanged `T`.
//! The single divergence is `merge::coalesce_coding_frame_separation`, whose
//! first act is:
//!
//! ```text
//! if !reading_frame || !length_changing || !pieces.windows(2).any(joins) {
//!     return None;
//! }
//! ```
//!
//! `reading_frame` comes from `merge::axis_frame`, which returns `true` for
//! `CisKind::Cds` and `false` for `CisKind::Genome`, `CisKind::Mt` and
//! `CisKind::Tx`. So on `c.` the pass re-absorbs the intervening reference base
//! into the first piece's payload and rebuilds the lone `delins` **before the
//! 3'-shift ever runs**; on `n.`/`g.` it declines on its first line, both pieces
//! survive, and the leading `insA` shifts into a `dup` wherever the preceding
//! base allows.
//!
//! Confirmed by forcing that guard to `return None` unconditionally in a scratch
//! copy of the tree: the `c.` axis then produces `c.[3_4insA;4_5insC]`,
//! `c.[8dup;9_10insC]` and `c.[51dup;52_53insC]` — the `n.` forms exactly, at
//! every position tested. Restoring the guard restores the lone `delins`.
//!
//! # The `c.` merge is unlicensed by the clause it cites
//!
//! `coalesce_coding_frame_separation` is the live implementation of
//! `general.md:35` / `delins.md:18`:
//!
//! > **exception**: two variants separated by one nucleotide, **together
//! > affecting one amino acid**, should be described as a "delins"
//!
//! It implements the *distance* half of that clause and **not** the amino-acid
//! half — it has no codon test at all. (The codon-precise pass,
//! `merge::apply_coding_codon_exception`, is guarded to equal-length blocks and
//! so cannot reach this shape.) That omission is load-bearing here, because
//! these rows are a net **+2** edit — a **frameshift**. A frameshifting pair
//! cannot satisfy "together affecting one amino acid", so the pass fires on a
//! block where its own cited precondition is unmet.
//!
//! That is a finding about the **justification**, and it is as far as the
//! evidence reaches. It is deliberately *not* a finding that the `c.` output is
//! wrong, for two separate reasons:
//!
//! * **Neither clause is normative.** `general.md:35` is lowercase prose, and so
//!   is the `general.md:34` separation rule the `n.`/`g.` axes follow — a census
//!   of the pinned spec checkout finds an uppercase RFC 2119 keyword in exactly
//!   one place outside `style.md` (see the repository `CLAUDE.md`). So neither
//!   side *requires* anything, and which recommendation ferro honours on a
//!   frameshifting block is a house-style question the spec leaves open.
//! * **Nothing has ruled on it.** The nearest record is
//!   `codon-carve-out-shape-restriction`, which is `undecided` and whose
//!   rationale says as much in terms: "WHETHER THAT NARROWING IS WRONG HAS NOT
//!   BEEN DEMONSTRATED". It does supply the constraint this file rests on —
//!
//!   > `general.md:35`'s "together affecting one amino acid" cannot cover a
//!   > frameshift pair, because a frameshift does not affect one amino acid.
//!
//!   — but a constraint recorded *inside* an unsettled question is not a
//!   verdict, and reading it as one is how a record that decides nothing
//!   smuggles in a decision nobody made.
//!
//! So the table below pins **observed behaviour on this base**, together with
//! the two checks that rule out the reading which would rescue the `c.` output —
//! "the two changes share a codon, so the exception applies". At **p=4** the two
//! changed columns straddle the codon 1 / codon 2 boundary, and at **p=52** they
//! straddle codon 17 / codon 18. The `c.` axis merges those two exactly as it
//! merges the twelve same-codon rows. The pass never asks the question, so
//! sharing a codon is not what is driving it.
//!
//! # Assert-then-flip
//!
//! The `c.` expectations below are ferro's current answers. The **candidate**
//! replacement is the `n.` form on `c.` coordinates — `c.[8dup;9_10insC]` for
//! p=9 — reached either by restoring `general.md:35`'s amino-acid precondition
//! to `coalesce_coding_frame_separation` (so it declines on frameshifting
//! blocks), or by an explicit ruling that ferro merges here as house style, in
//! which case the pass should stop citing a clause whose precondition it does
//! not check. Which of those two ships is exactly what is unsettled above; this
//! file does not answer it and must not be read as directing the flip.
//!
//! Asserted rather than `#[ignore]`d so the day this moves, this file is one of
//! the things that has to be read. It is an output-moving change on the `c.`
//! axis — PR #1484 measured the neighbouring widening at 1 957 real corpus rows
//! — so whichever way it is settled, the move owes the release a
//! `Representation-Change:` trailer.
//!
//! # What this file is NOT
//!
//! It is **not** the non-confluence it was originally chased as. A prior note
//! held that `c.9delinsATC` and `c.[8_9insA;9_10insC]` froze as *two* outputs
//! while the genomic pair converged, and attributed it to
//! `is_splittable_single_member` matching only `HV::Genome`/`HV::Mt`. Neither
//! half survives measurement on this base: that function admits all five axis
//! variants as of #1484, and the `c.` pair converges. What is left is a
//! cross-axis representation difference, which is a different (and less severe)
//! defect than a confluence failure — both spellings of the variant do reach one
//! answer on every axis. Fact 1 below is what stops that regressing.
//!
//! # Every pinned output is also checked against the bases it denotes
//!
//! The strings above are change detectors, and a string equality cannot tell a
//! **re-spelling** from a **corruption**: an output describing different bases
//! that happened to render as the pinned string would satisfy it. That is not a
//! hypothetical class here — it is what #1615 exists for, and #1592 and #1600
//! are live instances on `main` where a well-formed, in-bounds, idempotent,
//! re-parseable output denotes different bases.
//!
//! It bites hardest on this table because the outputs are `dup`s and `ins`s. A
//! `dup` names no bases in its rendering, so a composition that duplicated the
//! *wrong* base still prints the identical string.
//!
//! So each of the 42 comparisons — 14 rows x 3 axes — now also asserts that the
//! output and the input denote the same sequence, through `cis_apply_oracle` —
//! `hgvs_to_spdi` plus an SPDI splice, with no normalization in the path, so it
//! cannot agree with an output merely because normalization produced it (#1626).
//! All 42 pass.
//!
//! **28 of the 42 can fail; the `c.` 14 cannot, and saying which is which is the
//! point.** The `n.` and `g.` answers are genuine re-spellings — a lone `delins`
//! in, an insertion pair out — so a corrupted composition moves the applied
//! sequence and the comparison catches it. The `c.` answer is the lone input
//! *unchanged* (that non-split is the finding; [`ROWS`] records that the `c.`
//! output is `c.{position}delinsATC` at every one of the fourteen, which is
//! exactly the string [`spellings`] feeds in), so there `actual == input` and the
//! comparison is an identity. That leaves the `c.` column pinned but not freshly:
//! `CDS_START == 1`, so `c.p` and `n.p` are one position and the two axes' inputs
//! denote the same bases of [`CORE`] — the `n.` row's assertion therefore already
//! fixes what the `c.` answer denotes, transitively. What carries the `c.` row on
//! its own is the string equality.
//!
//! Note what this does **not** buy. It cannot say the `c.` merge is the right
//! form — the whole point of the assert-then-flip note above is that nothing has
//! ruled on that — only that today's answer and tomorrow's flipped one are
//! answers about the same bases.
//!
//! Fully hermetic: a `MockProvider`, no `FERRO_MANIFEST`, no fixtures.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

use crate::common::cis_apply_oracle::apply_reason;

/// Padding on each side of the core on the genomic contig, so that transcript
/// position `p` is genomic position `PAD_OFFSET + p`.
///
/// The **offset** matches `cis_confluence_axis.rs` and
/// `examples/generate_cis_confluence_corpus.rs` so the geometry is the one the
/// rest of the suite reasons about. The transcript's parent **contig name**
/// deliberately does not: those files serve the transcript's genome under a
/// separate `chr_synth`, which is fine there because they never address that
/// contig on the `g.` axis. This file does, so it uses one name for both — see
/// [`transcript_provider`].
const PAD_OFFSET: usize = 256;

const GENOMIC_CONTIG: &str = "NC_TEST.1";
const CODING_TX: &str = "NM_TEST.1";
const NONCODING_TX: &str = "NR_TEST.1";

/// The 63 nt transcript. `CDS_START = 1` on the coding transcript, so `c.p` and
/// `n.p` address the same base and the three axes are directly comparable.
///
/// Chosen for low local repetitiveness rather than at random: a homopolymer run
/// around `p` would let the 3'-shift dominate and would measure shuffling rather
/// than partitioning.
const CORE: &str = "ACGTACGATGCAAGTCCTGAGGTCAATCGGATCCTAGGACTTCAGGTACCATGAGTCATGCAA";
const CDS_START: u64 = 1;
const CDS_END: u64 = 63;

/// One measured position: `(position, same_codon, noncoding, genomic)`.
///
/// * `position` — 1-based transcript position of the reference `T` replaced.
/// * `same_codon` — whether the two changed columns lie in one codon
///   (`CDS_START == 1`). `false` at p=4 and p=52, the rows that refute the
///   shared-codon reading of the `c.` behaviour.
/// * `noncoding`, `genomic` — the outputs measured on this base.
///
/// A tuple rather than a named struct so [`ROWS`] fits one line per position:
/// the point of the table is that fourteen rows say the same thing, and a
/// six-line-per-row expansion hides exactly that.
///
/// The `c.` output is not stored. It is `c.{position}delinsATC` at every one of
/// the fourteen, and that uniformity is half the finding, so it is computed and
/// the computation is the assertion.
type Row = (u64, bool, &'static str, &'static str);

/// Every position in [`CORE`] carrying a `T`, with room for a neighbour on each
/// side. Measured on this base; see the module docs.
///
/// **These are change detectors, not conformance expectations.** No ruling
/// blesses any string in this table, and the repository `CLAUDE.md` is explicit
/// that pinning today's output is not an adjudication. What the table is *for* is
/// the cross-axis disagreement it makes visible — 14 of 14, unanimous — and the
/// two `same_codon == false` rows, which refute the shared-codon explanation of
/// the `c.` column. Read a failure here as "the partitioner moved", then decide
/// whether that was intended; do not read the strings as the answers the spec
/// requires.
#[rustfmt::skip]
const ROWS: &[Row] = &[
    (4,  false, "NR_TEST.1:n.[3_4insA;4_5insC]",     "NC_TEST.1:g.[259_260insA;260_261insC]"),
    (9,  true,  "NR_TEST.1:n.[8dup;9_10insC]",       "NC_TEST.1:g.[264dup;265_266insC]"),
    (15, true,  "NR_TEST.1:n.[14_15insA;17dup]",     "NC_TEST.1:g.[270_271insA;273dup]"),
    (18, true,  "NR_TEST.1:n.[17_18insA;18_19insC]", "NC_TEST.1:g.[273_274insA;274_275insC]"),
    (23, true,  "NR_TEST.1:n.[22_23insA;24dup]",     "NC_TEST.1:g.[278_279insA;280dup]"),
    (27, true,  "NR_TEST.1:n.[26dup;28dup]",         "NC_TEST.1:g.[282dup;284dup]"),
    (32, true,  "NR_TEST.1:n.[31dup;34dup]",         "NC_TEST.1:g.[287dup;290dup]"),
    (35, true,  "NR_TEST.1:n.[34_35insA;35_36insC]", "NC_TEST.1:g.[290_291insA;291_292insC]"),
    (41, true,  "NR_TEST.1:n.[40_41insA;41_42insC]", "NC_TEST.1:g.[296_297insA;297_298insC]"),
    (42, true,  "NR_TEST.1:n.[41_42insA;43dup]",     "NC_TEST.1:g.[297_298insA;299dup]"),
    (47, true,  "NR_TEST.1:n.[46_47insA;47_48insC]", "NC_TEST.1:g.[302_303insA;303_304insC]"),
    (52, false, "NR_TEST.1:n.[51dup;52_53insC]",     "NC_TEST.1:g.[307dup;308_309insC]"),
    (56, true,  "NR_TEST.1:n.[55_56insA;57dup]",     "NC_TEST.1:g.[311_312insA;313dup]"),
    (59, true,  "NR_TEST.1:n.[58dup;59_60insC]",     "NC_TEST.1:g.[314dup;315_316insC]"),
];

fn padded() -> String {
    let pad = "ACGT".repeat(PAD_OFFSET / 4);
    format!("{pad}{CORE}{pad}")
}

fn genomic_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(GENOMIC_CONTIG, padded());
    provider
}

/// A single-exon transcript over [`CORE`]. `cds` present makes it a `c.` axis
/// reference; absent makes it `n.`. Nothing else differs between the two.
///
/// The transcript's genomic parent is [`GENOMIC_CONTIG`] — the same contig the
/// `g.` arm addresses — so "three axes over the same 63 bases" holds by
/// construction rather than by two providers happening to be handed the same
/// padded string under different contig names.
fn transcript_provider(accession: &str, cds: Option<(u64, u64)>) -> MockProvider {
    let mut provider = MockProvider::new();
    let tx_len = CORE.len() as u64;
    let g_start = PAD_OFFSET as u64 + 1;
    let g_end = PAD_OFFSET as u64 + tx_len;
    let transcript = Transcript::new(
        accession.to_string(),
        Some("SYNTH".to_string()),
        Strand::Plus,
        CORE.to_string(),
        cds.map(|c| c.0),
        cds.map(|c| c.1),
        vec![Exon::with_genomic(1, 1, tx_len, g_start, g_end)],
        Some(GENOMIC_CONTIG.to_string()),
        Some(g_start),
        Some(g_end),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );
    provider.add_genomic_sequence(GENOMIC_CONTIG, padded());
    provider.add_transcript(transcript);
    provider
}

fn normalize(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::default());
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    normalizer
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("normalize {input}: {e}"))
        .to_string()
}

/// The sequence `description` denotes, reached **without** the normalizer
/// (#1626).
///
/// `cis_apply_oracle` converts each member through `hgvs_to_spdi` and splices
/// the reference, so nothing in this path can agree with an output merely
/// because normalization produced it. Reusing it rather than hand-rolling a
/// second applier is the point: a per-file copy is a per-file way to drift.
///
/// It panics rather than returning an `Option`. Every description this file
/// feeds it is a lone `delins` or a two-member insertion pair on a synthetic
/// reference, so a decline is a defect and not a limit of the oracle — and a
/// comparison that silently skips is the failure mode this check exists to
/// remove.
fn denotes(provider: &MockProvider, reference: &str, description: &str) -> String {
    apply_reason(provider, reference, description)
        .unwrap_or_else(|why| panic!("{description} denotes no single sequence: {why:?}"))
}

/// Assert that `input` normalizes to `expected` on `provider` **and** that the
/// two denote the same bases of `reference`.
///
/// The second half is what is new (#1626). On `n.` and `g.` the answer is a
/// **re-spelling** — a one-base span with a three-base payload comes back as an
/// insertion pair — and a string equality cannot tell a re-spelling from a
/// corruption. It is especially blind on those outputs, since a `dup` names no
/// bases in its rendering: a composition that duplicated the *wrong* base still
/// prints the identical string.
///
/// On `c.` the answer is the lone input *unchanged*, so `actual == input` and the
/// second assertion is an identity there — it cannot fail. It is applied
/// uniformly rather than special-cased: the row is carried by the string
/// equality, and what the `c.` answer denotes is fixed transitively by the `n.`
/// row, whose input denotes the same bases (`CDS_START == 1`). See the module
/// docs.
///
/// Nor does the sibling [`every_axis_is_self_confluent_across_both_spellings`]
/// close it. Agreement between two spellings cannot distinguish "both correct"
/// from "both wrong" — the argument `cis_allele_confluence_proptest.rs` makes
/// about itself.
fn assert_normalizes_preserving(
    provider: &MockProvider,
    reference: &str,
    input: &str,
    expected: &str,
    context: &str,
) {
    let actual = normalize(provider.clone(), input);
    assert_eq!(actual, expected, "{context}");
    assert_eq!(
        denotes(provider, reference, &actual),
        denotes(provider, reference, input),
        "{input} -> {actual} is not a re-spelling: it denotes different bases"
    );
}

/// The two spellings of "replace the `T` at `position` with `ATC`", on one axis.
fn spellings(accession: &str, axis: char, position: u64) -> (String, String) {
    (
        format!("{accession}:{axis}.{position}delinsATC"),
        format!(
            "{accession}:{axis}.[{}_{position}insA;{position}_{}insC]",
            position - 1,
            position + 1
        ),
    )
}

/// Guards the fixture the two tests below rest on. A `CORE` edit that moved a
/// `T`, or a padding change, would otherwise silently re-point every row.
#[test]
fn the_fixture_is_the_one_the_rows_were_measured_on() {
    assert_eq!(CORE.len() as u64, CDS_END, "CDS spans the whole transcript");
    assert_eq!(PAD_OFFSET % 4, 0, "the ACGT pad must tile exactly");
    for &(position, same_codon, _, _) in ROWS {
        let base = CORE.as_bytes()[position as usize - 1];
        assert_eq!(
            base as char, 'T',
            "c.{position} must be a T for `delinsATC` to retain it"
        );
        // `CDS_START == 1`, so codon k spans positions 3k-2 ..= 3k.
        let expected_same_codon = position.div_ceil(3) == (position - 1).div_ceil(3);
        assert_eq!(
            same_codon, expected_same_codon,
            "c.{position}: recorded same_codon does not match the frame"
        );
    }
    assert_eq!(
        ROWS.iter()
            .filter(|(_, same_codon, _, _)| !same_codon)
            .count(),
        2,
        "p=4 and p=52 are the codon-straddling rows that refute the shared-codon reading; \
         losing them would hollow out the argument in this file's docs"
    );
}

/// **Fact 1.** Both spellings reach one output, on every axis, at every
/// position. This is what #1484's axis widening bought; a regression here is a
/// confluence regression and is more serious than the disagreement Fact 2 pins.
#[test]
fn every_axis_is_self_confluent_across_both_spellings() {
    for &(p, _, _, _) in ROWS {
        for (accession, axis, provider) in [
            (
                CODING_TX,
                'c',
                transcript_provider(CODING_TX, Some((CDS_START, CDS_END))),
            ),
            (NONCODING_TX, 'n', transcript_provider(NONCODING_TX, None)),
            (GENOMIC_CONTIG, 'g', genomic_provider()),
        ] {
            let position = if axis == 'g' {
                p + PAD_OFFSET as u64
            } else {
                p
            };
            let (lone, allele) = spellings(accession, axis, position);
            let from_lone = normalize(provider.clone(), &lone);
            let from_allele = normalize(provider, &allele);
            assert_eq!(
                from_lone, from_allele,
                "{axis}. axis at position {position}: the lone `delins` spelling and the \
                 two-member spelling denote one variant and must normalize alike"
            );
        }
    }
}

/// **Fact 2.** The three axes do not agree with each other: `n.` and `g.` split,
/// `c.` does not.
///
/// The `c.` expectation is **assert-then-flip** — see the module docs. When
/// `coalesce_coding_frame_separation` stops firing on frameshifting blocks, the
/// `c.` column becomes the `n.` column with `n.` rewritten to `c.`.
#[test]
fn the_coding_axis_merges_what_the_other_axes_split() {
    for &(p, same_codon, noncoding, genomic) in ROWS {
        // `n.` and `g.` — the same two-member answer, 256 apart.
        let (n_lone, _) = spellings(NONCODING_TX, 'n', p);
        assert_normalizes_preserving(
            &transcript_provider(NONCODING_TX, None),
            CORE,
            &n_lone,
            noncoding,
            &format!("n. axis at {p}"),
        );
        let (g_lone, _) = spellings(GENOMIC_CONTIG, 'g', p + PAD_OFFSET as u64);
        assert_normalizes_preserving(
            &genomic_provider(),
            &padded(),
            &g_lone,
            genomic,
            &format!("g. axis at {p}"),
        );

        // `c.` — the lone `delins`, at every position, in and out of frame.
        let (c_lone, _) = spellings(CODING_TX, 'c', p);
        let coding_n_parity = noncoding.replace("NR_TEST.1:n.", "NM_TEST.1:c.");
        assert_normalizes_preserving(
            &transcript_provider(CODING_TX, Some((CDS_START, CDS_END))),
            CORE,
            &c_lone,
            &format!("NM_TEST.1:c.{p}delinsATC"),
            &format!(
            "c. axis at {p} (same_codon={same_codon}): `coalesce_coding_frame_separation` re-merges the \
             two pieces because `axis_frame` gives `CisKind::Cds` a reading frame, without \
             checking `general.md:35`'s \"together affecting one amino acid\" precondition — \
             which a net +2 frameshift cannot meet. The `n.`-parity form is {coding_n_parity}. \
             If this assertion just failed with that string, the pass has stopped firing on this \
             shape: that is the candidate flip in this file's docs, so check that the house-style \
             question it names has actually been settled, then flip the expectation and declare \
             the representation change"
            ),
        );

        // Stated explicitly so the flip cannot be mistaken for a no-op.
        assert_ne!(
            format!("NM_TEST.1:c.{p}delinsATC"),
            coding_n_parity,
            "at {p} the pinned `c.` output already equals the `n.`-parity form — the \
             assert-then-flip note in this file's docs is stale"
        );
    }
}
