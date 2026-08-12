//! The mechanism behind the corpus's largest class: a `c.`/`n.` output naming an
//! intronic position its input did not.
//!
//! # What this file pins, and what it deliberately does not
//!
//! `spec_conformance_axis` measured `outputs_leaving_the_transcript` at **371**
//! on the 3' direction and **0** on the 5'; `spec_corpus_regressions`'s
//! `a_minus_strand_junction_shift_leaves_the_transcript` names one row of it. Both
//! describe the *symptom*. These tests pin the **mechanism**, so that a fix can be
//! judged against the cause rather than against one string.
//!
//! # FIXED by #1704 — and the mechanism below is unchanged, which is the point
//!
//! The file kept its name and every one of its measurements, because none of them
//! was wrong: the trigger is still the landing on the exon's last base, the
//! traveller is still the intron's own bases rather than the strand, the exon→exon
//! case is still correctly refused, and the 5' direction still has no mirror of the
//! gate. **What changed is one thing only — the accession the answer is rendered
//! against.** `checklist.md:20` needs a genomic reference for an intronic position,
//! the genomic re-shuffle already resolved one to compute the crossing at all, and
//! `Normalizer::reparent_junction_exit` now renders the result as
//! `NC_SYNTH.1(NM_TEST.1):c.10+2del`. Every expectation below moved by exactly that
//! prefix and by nothing else; a diff of this file against its pre-#1704 form is
//! the cleanest statement available of what the fix did and did not do.
//!
//! The two candidates named at the end of this doc are therefore settled:
//! re-parenting, not refusing. `Normalizer::reparent_junction_exit` carries the
//! reasoning, including why refusing could not have closed the class (the census is
//! a *lenient*-mode figure, so a strict-mode refusal leaves all 371 in place).
//!
//! # The fixtures are built here rather than drawn from `spec_corpus`
//!
//! Deliberate, for two reasons. A named regression test must survive an edit to
//! the generator (the #1456/#1460/#1478 lesson that
//! `spec_corpus_regressions.rs` opens with). And — decisively for *this* defect —
//! the corpus's synthetic intron is a **confound**: `transcript_provider` lays a
//! literal `"GATTACA"…` intron into the contig and reverse-complements only the
//! exon blocks, so the intron read in *transcript* direction is `GATT…` on the
//! plus strand and `AATC…` on the minus. Against the corpus's `AT`-alphabet cores
//! that single fact decides which strand leaks. [`junction_provider`] below takes
//! the intron as a parameter and lays it out in **transcript** direction on both
//! strands, which is what makes the strand a controlled variable instead of a
//! hidden one.
//!
//! # The three findings, in order of how badly they were mis-stated before
//!
//! 1. **It is not minus-strand-only.** The plus strand leaves the exon under an
//!    intron that continues the run, and the minus strand does not under one that
//!    does not. See [`the_exit_follows_the_intron_bases_not_the_strand`].
//! 2. **It is not a missing clamp.** The exon-confined shuffle is correctly
//!    bounded (`normalize/boundary.rs:142-160`). What leaves the exon is a
//!    deliberate second pass, `#670` at `src/normalize/mod.rs:4768-4834`, which
//!    fires **only** when the confined shuffle lands exactly on the exon's last
//!    base, re-runs the shuffle in genomic space via
//!    `normalize_boundary_spanning_cds`, and adopts the answer **only if it
//!    crossed into the intron**. The intronic output is that branch's success
//!    condition. See [`the_exit_needs_the_shuffle_to_land_on_the_exons_last_base`].
//! 3. **The 3'/0 asymmetry is one `if`.** All three copies of that gate —
//!    `normalize_cds` (`mod.rs:4785`), `normalize_tx` (`:5407`), `normalize_rna`
//!    (`:6917`) — are guarded on `shuffle_direction == ThreePrime`, and there is
//!    no 5' mirror anywhere. So the corpus's `0` at 5' is a claim about the
//!    **code**, not about the corpus: the property is varied here across four
//!    intron alphabets and both strands and still measures zero. See
//!    [`the_five_prime_direction_never_leaves_the_exon`].
//!
//! # The clause the existing tests cite says the opposite of what they say it does
//!
//! `spec_conformance_axis:45` and `spec_corpus_regressions:70-109` both read
//! `general.md:44` as exempting a junction-adjacent deletion from the 3' rule, and
//! pin `NM_TEST.1:c.20del` as the correct answer. `general.md:44` points at
//! `background/numbering.md#DNAc`, and that passage is narrower in one direction
//! and explicit in the other:
//!
//! - `background/numbering.md:23` — "the 3' rule is not applied when there is a
//!   deletion/duplication around **exon/exon** junctions with identical
//!   nucleotides flanking the junction, where shifting the variant 3' would place
//!   it **in the next exon**."
//! - `background/numbering.md:26` — "**NOTE**: this exception **does not apply**
//!   to a deletion/duplication around exon/**intron** and intron/exon junctions
//!   with identical nucleotides flanking the junction".
//!
//! So the exception covers the exon→next-exon case, which ferro already blocks in
//! two independent places (the exon-confined shuffle bound, and
//! `merge.rs:2336`/`:2350`'s `crosses_exon_junction` / `enclosing_exon`), and
//! `:26` explicitly withholds it from the exon→intron case that `#670`
//! implements. **The shift is right.** What is wrong is the accession it is
//! rendered against: `checklist.md:20` — an `NM_` "can only be used to describe
//! variants in introns using a `c.` prefix when a genomic reference sequence is
//! given", e.g. `NC_000023.10(NM_004006.2):c.94-2A>G`.
//!
//! That relocates the defect from the coordinate to the rendering, and it rules
//! one candidate fix out: converging on `c.20del` would implement the exception
//! `numbering.md:26` withholds, and would silently revert `#670`, whose own
//! committed fixtures (`src/normalize/mod.rs:13329-13405`) pin the intronic answer
//! as correct. The two live candidates are re-parenting the output onto
//! `NC_SYNTH.1(NM_TEST.1)` or refusing it.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

// ---------------------------------------------------------------------------
// A three-exon fixture whose intron is a parameter
// ---------------------------------------------------------------------------

pub(crate) const CODING: &str = "NM_TEST.1";
pub(crate) const NONCODING: &str = "NR_TEST.1";
pub(crate) const CONTIG: &str = "NC_SYNTH.1";

/// Bases of `ACGT` padding either side of the transcript's genomic span.
///
/// **160 is a floor, not a decoration.** `normalize_boundary_spanning_cds`
/// fetches `g +/- NormalizeConfig::window_size` (100 by default) and
/// `MockProvider::get_genomic_sequence` *errors* rather than truncating when that
/// runs off the contig — and the `#670` call site swallows the error
/// (`if let Ok(..)`), falling through to the exon-confined answer. A short pad
/// therefore turns every test below silently green for the wrong reason. Measured:
/// at `PAD = 64` the minus-strand fixture reported `c.10del` for a neighbourhood
/// that yields `c.10+2del`, because exon 1 sits at the *high* end of a
/// minus-strand contig and `g + 100` overran it.
const PAD: usize = 160;

/// Exon block length. Three of them, so the transcript is 60 bases.
const EXON_LEN: usize = 20;

/// 1-based inclusive CDS bounds. `cds_start = 11` puts `c.1` at tx 11, which
/// makes **`c.10` the last base of exon 1** — the junction every test below
/// works against. `cds_end = 52` leaves an eight-base 3'UTR.
const CDS: (u64, u64) = (11, 52);

/// Reverse complement, uppercase DNA.
fn revcomp(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// A provider holding one three-exon transcript and the contig it sits on.
///
/// `exons` are the three 20-base blocks in **transcript** order and `intron` is
/// the intronic block, also written in **transcript** direction — the contig is
/// reverse-complemented wholesale for a minus-strand transcript, so a caller
/// sees the identical transcript-direction neighbourhood on both strands. That
/// is the property the corpus's own fixture lacks, and it is what turns the
/// strand into a controlled variable here.
///
/// `coding` selects the accession and whether [`CDS`] is applied: an `NR_`
/// record carries no CDS, which is what puts the `n.` axis on the same footing.
pub(crate) fn junction_provider(
    strand: Strand,
    exons: [&str; 3],
    intron: &str,
    coding: bool,
) -> MockProvider {
    junction_provider_on_chromosome(strand, exons, intron, coding, Some(CONTIG))
}

pub(crate) fn junction_provider_on_chromosome(
    strand: Strand,
    exons: [&str; 3],
    intron: &str,
    coding: bool,
    chromosome: Option<&str>,
) -> MockProvider {
    junction_provider_named(
        strand,
        exons,
        intron,
        coding,
        chromosome,
        if coding { CODING } else { NONCODING },
    )
}

/// As [`junction_provider_on_chromosome`], but the transcript's own accession is
/// a parameter too.
///
/// Only [`a_self_referential_wrapper_is_declined`] needs this: every other test
/// works on `NM_`/`NR_`, where the accession is decided by `coding`. `coding`
/// here selects only whether [`CDS`] is applied, so an LRG record can be served
/// with coding coordinates.
pub(crate) fn junction_provider_named(
    strand: Strand,
    exons: [&str; 3],
    intron: &str,
    coding: bool,
    chromosome: Option<&str>,
    accession: &str,
) -> MockProvider {
    for block in exons {
        assert_eq!(block.len(), EXON_LEN, "each exon block is {EXON_LEN} bases");
    }
    let tx: String = exons.concat();
    // The unspliced transcript-direction sequence, and each exon's 0-based
    // offset within it.
    let mut unspliced = String::new();
    let mut offsets = Vec::with_capacity(3);
    for (index, block) in exons.iter().enumerate() {
        if index > 0 {
            unspliced.push_str(intron);
        }
        offsets.push(unspliced.len());
        unspliced.push_str(block);
    }

    let pad = "ACGT".repeat(PAD / 4);
    let oriented = if strand == Strand::Plus {
        unspliced.clone()
    } else {
        revcomp(&unspliced)
    };
    let contig = format!("{pad}{oriented}{pad}");

    // 1-based genomic bounds of each exon, in transcript order. On the minus
    // strand a block at transcript-direction offset `o` occupies the mirrored
    // window, so exon 1 sits at the HIGHEST genomic coordinates.
    let span = unspliced.len();
    let exon_records: Vec<Exon> = offsets
        .iter()
        .enumerate()
        .map(|(index, &offset)| {
            let (g_start, g_end) = if strand == Strand::Plus {
                (PAD + offset + 1, PAD + offset + EXON_LEN)
            } else {
                (PAD + span - offset - EXON_LEN + 1, PAD + span - offset)
            };
            Exon::with_genomic(
                index as u32 + 1,
                (index * EXON_LEN + 1) as u64,
                ((index + 1) * EXON_LEN) as u64,
                g_start as u64,
                g_end as u64,
            )
        })
        .collect();

    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(CONTIG, contig.clone());
    // Serve the contig under whatever the transcript NAMES as its chromosome, so
    // the crossing itself succeeds and the rendering guards are the only thing
    // left to judge. Without this, renaming the chromosome also breaks the
    // genomic fetch the `#670` gate needs, the answer falls back to the
    // exon-confined `c.10del`, and a guard test would pass for the wrong reason.
    if let Some(name) = chromosome.filter(|name| *name != CONTIG) {
        provider.add_genomic_sequence(name, contig);
    }
    provider.add_transcript(Transcript::new(
        accession.to_string(),
        Some("SYNTH".to_string()),
        strand,
        tx,
        coding.then_some(CDS.0),
        coding.then_some(CDS.1),
        exon_records,
        chromosome.map(str::to_string),
        Some(PAD as u64 + 1),
        Some((PAD + span) as u64),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

/// Exon blocks whose only ambiguity is a four-base `run` of `base` flush against
/// **exon 1's 3' end** (tx 17..=20, i.e. `c.7`..`c.10`).
///
/// Everything else is `C`/`G` filler that cannot extend an `A` or `T` run, so
/// the shift has exactly one place to go and the tests read as an A/B on the
/// intron.
pub(crate) fn exons_with_run_at_exon1_end(base: char) -> [String; 3] {
    let filler = if base == 'C' { 'G' } else { 'C' };
    let mut exon1: String = std::iter::repeat_n(filler, EXON_LEN).collect();
    exon1.replace_range(16..20, &std::iter::repeat_n(base, 4).collect::<String>());
    let rest: String = std::iter::repeat_n(filler, EXON_LEN).collect();
    [exon1, rest.clone(), rest]
}

/// The mirror: a four-base run flush against **exon 2's 5' end** (tx 21..=24,
/// i.e. `c.11`..`c.14`), for the 5'-direction tests.
fn exons_with_run_at_exon2_start(base: char) -> [String; 3] {
    let filler = if base == 'C' { 'G' } else { 'C' };
    let mut exon2: String = std::iter::repeat_n(filler, EXON_LEN).collect();
    exon2.replace_range(0..4, &std::iter::repeat_n(base, 4).collect::<String>());
    let rest: String = std::iter::repeat_n(filler, EXON_LEN).collect();
    [rest.clone(), exon2, rest]
}

/// A 30-base intron whose first bases (in transcript direction) are `lead`.
///
/// The tail alternates between two bases **chosen to differ from `lead`'s last
/// character**, so `lead` alone decides how far a shift can travel past the
/// junction. Deriving the tail rather than fixing it at `CGCG…` is not
/// fastidiousness: with a fixed tail, `intron_leading_with("CC")` yields
/// `CCC…` and silently lengthens the run by one, which showed up as a
/// `c.10+3del` where `c.10+2del` was pinned.
pub(crate) fn intron_leading_with(lead: &str) -> String {
    let last = lead.chars().next_back().expect("a non-empty lead");
    let tail: Vec<char> = "ACGT"
        .chars()
        .filter(|&base| base != last)
        .take(2)
        .collect();
    let mut intron = lead.to_string();
    while intron.len() < 30 {
        intron.push(tail[intron.len() % 2]);
    }
    assert_eq!(intron.len(), 30, "the intron is a fixed 30 bases");
    assert_ne!(
        intron.as_bytes()[lead.len()],
        last as u8,
        "the tail must not extend the lead, or the run length is not `lead`"
    );
    intron
}

fn normalize(provider: &MockProvider, input: &str, direction: ShuffleDirection) -> String {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));
    let normalizer = Normalizer::with_config(
        provider.clone(),
        NormalizeConfig::default().with_direction(direction),
    );
    normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("{input} must normalize: {e}"))
        .to_string()
}

fn normalize_3prime(provider: &MockProvider, input: &str) -> String {
    normalize(provider, input, ShuffleDirection::ThreePrime)
}

// ---------------------------------------------------------------------------
// 1. The strand is a confound
// ---------------------------------------------------------------------------

/// **Question.** The corpus records this class as "**entirely** minus-strand",
/// and `a_minus_strand_junction_shift_leaves_the_transcript` pins the plus strand
/// as the control that "clamps correctly". Is the strand the variable?
///
/// **No — the intron's bases are.** Holding the strand fixed and varying only the
/// two intronic bases immediately 3' of the junction flips the answer on *both*
/// strands, and holding the intron fixed makes the two strands agree exactly.
///
/// Where the corpus's appearance of strand-specificity comes from is arithmetic
/// rather than behaviour: `spec_corpus::transcript_provider` writes a literal
/// `"GATTACA"…` intron into the contig and reverse-complements only the exon
/// blocks, so in transcript direction the intron reads `GATT…` on the plus strand
/// and `AATC…` on the minus — and the corpus's cores for this stratum are drawn
/// from an `AT` alphabet. `G` never continues an `A`/`T` run; `AA` always does.
///
/// **This matters for the fix.** A fix keyed on `Strand::Minus`, or a regression
/// test that treats the plus strand as a passing control, would be keyed on the
/// fixture's intron.
#[test]
fn the_exit_follows_the_intron_bases_not_the_strand() {
    for base in ['A', 'C', 'G', 'T'] {
        let exons = exons_with_run_at_exon1_end(base);
        let blocks = [exons[0].as_str(), exons[1].as_str(), exons[2].as_str()];

        // An intron that continues the run: the shift leaves the exon by exactly
        // the two bases the intron extends it, on BOTH strands.
        let continuing = intron_leading_with(&format!("{base}{base}"));
        for strand in [Strand::Plus, Strand::Minus] {
            let provider = junction_provider(strand, blocks, &continuing, true);
            assert_eq!(
                normalize_3prime(&provider, &format!("{CODING}:c.7del")),
                format!("{CONTIG}({CODING}):c.10+2del"),
                "base {base} on the {strand:?} strand, intron {continuing}. The SHIFT is \
                 correct (numbering.md:26 withholds the exon/exon exception from exon/intron \
                 junctions) and the ACCESSION now carries the genomic wrapper checklist.md:20 \
                 requires (#1704)."
            );
        }

        // An intron that does not: the shift stops at the exon's last base, on
        // both strands. Same code path, same strands, opposite answer.
        let blocking = intron_leading_with(if base == 'C' { "GG" } else { "CC" });
        for strand in [Strand::Plus, Strand::Minus] {
            let provider = junction_provider(strand, blocks, &blocking, true);
            assert_eq!(
                normalize_3prime(&provider, &format!("{CODING}:c.7del")),
                format!("{CODING}:c.10del"),
                "base {base} on the {strand:?} strand with a non-continuing intron must stop \
                 at the junction — this is the control that makes the row above a defect"
            );
        }
    }
}

/// **Question.** How far does the description travel once it leaves the exon?
///
/// **Exactly as far as the intron continues the run**, which pins that this is an
/// ordinary sequence-following 3' shift over the *unspliced* neighbourhood rather
/// than a fixed one-base overrun. It also shows the exit is not a `+1` off-by-one:
/// a three-base continuation yields `c.10+3`.
#[test]
fn the_distance_travelled_is_the_length_of_the_matching_intron_prefix() {
    let exons = exons_with_run_at_exon1_end('A');
    let blocks = [exons[0].as_str(), exons[1].as_str(), exons[2].as_str()];
    // The `CC` row stays inside the exon, so it keeps the bare accession; every
    // row that travels into the intron carries the genomic wrapper (#1704).
    for (lead, expected) in [
        ("CC", format!("{CODING}:c.10del")),
        ("AC", format!("{CONTIG}({CODING}):c.10+1del")),
        ("AAC", format!("{CONTIG}({CODING}):c.10+2del")),
        ("AAAC", format!("{CONTIG}({CODING}):c.10+3del")),
    ] {
        let intron = intron_leading_with(lead);
        for strand in [Strand::Plus, Strand::Minus] {
            let provider = junction_provider(strand, blocks, &intron, true);
            assert_eq!(
                normalize_3prime(&provider, &format!("{CODING}:c.7del")),
                expected,
                "intron prefix {lead} on the {strand:?} strand"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// 2. It is not a missing clamp — the trigger is a landing, and it is narrow
// ---------------------------------------------------------------------------

/// **Question.** Is the exon-junction clamp missing, or is it the wrong
/// coordinate space?
///
/// **Neither.** The exon-confined shuffle is present and correct; what leaves the
/// exon is a *second*, deliberate pass. `#670` at `src/normalize/mod.rs:4768`
/// fires only when the confined shuffle's result ends exactly on the exon's last
/// base (`new_tx_end == exon_only.right`), then re-runs the shuffle in genomic
/// space and adopts it only if it crossed into the intron.
///
/// The gate's narrowness is the observable signature, and it is what this pins:
/// move the same run **one base 5'** so the confined shuffle lands at `c.9`
/// instead of `c.10`, and the identical continuing intron produces no exit at
/// all. A genuinely missing clamp would leak in both placements; a wrong
/// coordinate space would leak whenever the genomic neighbourhood matched.
#[test]
fn the_exit_needs_the_shuffle_to_land_on_the_exons_last_base() {
    let intron = intron_leading_with("AA");

    // Flush against the junction: tx 17..=20 == c.7..=c.10. Exits.
    let flush = exons_with_run_at_exon1_end('A');
    // One base short: tx 16..=19 == c.6..=c.9, with tx 20 broken. Does not.
    let mut short = flush.clone();
    short[0].replace_range(15..20, "AAAAC");

    for strand in [Strand::Plus, Strand::Minus] {
        let flush_blocks = [flush[0].as_str(), flush[1].as_str(), flush[2].as_str()];
        assert_eq!(
            normalize_3prime(
                &junction_provider(strand, flush_blocks, &intron, true),
                &format!("{CODING}:c.7del")
            ),
            format!("{CONTIG}({CODING}):c.10+2del"),
            "the {strand:?}-strand run landing ON the junction exits"
        );

        let short_blocks = [short[0].as_str(), short[1].as_str(), short[2].as_str()];
        assert_eq!(
            normalize_3prime(
                &junction_provider(strand, short_blocks, &intron, true),
                &format!("{CODING}:c.6del")
            ),
            format!("{CODING}:c.9del"),
            "a run ending ONE BASE short of the junction never reaches the #670 gate, under \
             the identical intron and strand — so the exit is a landing trigger, not a \
             missing clamp ({strand:?})"
        );
    }
}

/// **Question.** Is the exon→**exon** case — the one `numbering.md:23` really
/// does exempt — also broken?
///
/// **No, and that is worth pinning as a negative result.** With a run flush
/// against exon 1's end *and* continuing into exon 2's first bases, the shift
/// stops at the junction: it does not roll into the next exon. So the clause
/// ferro is accused of violating is in fact the one it honours, and the clause it
/// honours (`numbering.md:26`, the 3' rule applying across an exon/intron border)
/// is the one that produces the 371 rows.
///
/// Recorded per `CLAUDE.md`'s "record what was refuted, not only what was
/// decided": the belief that the 371 class is a junction-clamp failure is the
/// belief this measurement kills.
#[test]
fn the_shift_still_refuses_to_cross_into_the_next_exon() {
    // A run spanning the junction in TRANSCRIPT space: exon 1 ends `AAAA`,
    // exon 2 begins `AAAA`. The intron between them blocks, so the only way to
    // travel is through the splice.
    let filler: String = std::iter::repeat_n('C', EXON_LEN).collect();
    let mut exon1 = filler.clone();
    exon1.replace_range(16..20, "AAAA");
    let mut exon2 = filler.clone();
    exon2.replace_range(0..4, "AAAA");
    let intron = intron_leading_with("CC");

    for strand in [Strand::Plus, Strand::Minus] {
        let provider = junction_provider(
            strand,
            [exon1.as_str(), exon2.as_str(), filler.as_str()],
            &intron,
            true,
        );
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.7del")),
            format!("{CODING}:c.10del"),
            "numbering.md:23 exempts a junction-adjacent deletion from the 3' rule where \
             shifting would place it IN THE NEXT EXON; ferro honours that ({strand:?})"
        );
    }
}

// ---------------------------------------------------------------------------
// 3. The 3'/5' asymmetry is one `if`
// ---------------------------------------------------------------------------

/// **Question.** The corpus measures `outputs_leaving_the_transcript` at 371 on
/// the 3' direction and **0** on the 5'. Is that zero a claim about the corpus?
///
/// **No — it is a claim about the code, and this varies the property to show it.**
/// All three copies of the junction-crossing continuation are gated on
/// `shuffle_direction == ThreePrime` (`normalize_cds` `mod.rs:4785`,
/// `normalize_tx` `:5407`, `normalize_rna` `:6917`) and there is no 5' mirror
/// anywhere in the crate.
///
/// Swept here across four run alphabets and both strands, with an intron whose
/// **transcript-direction 3' end** continues the run — i.e. the exact mirror of
/// the geometry that makes the 3' direction leak — a 5' shuffle stops at exon 2's
/// first base every time. Eight placements (four run alphabets x two strands),
/// zero exits.
///
/// **This is a structural zero, and it is reported as one rather than suppressed**
/// (the honest form PR #1480 established). It is not evidence that the 5'
/// direction is correct: it is evidence that the 5' direction cannot express the
/// answer the 3' rule's mirror would give, which makes the two directions
/// non-confluent across a junction.
#[test]
fn the_five_prime_direction_never_leaves_the_exon() {
    let mut exits = Vec::new();
    // `(base, strand, output)` for every placement, collected across the WHOLE
    // sweep before anything is asserted. Pinning inside the loop is what made
    // `exits` dead: the exact-string assert fires on the first offending row, so
    // the sweep aborts and the `is_empty()` check below never sees a non-empty
    // vector — it could only ever report a property it had not finished
    // measuring. Ordering the push before the assert *within* an iteration does
    // not fix that; the two assertions have to be separated by the whole sweep.
    let mut placements = Vec::new();
    for base in ['A', 'C', 'G', 'T'] {
        let exons = exons_with_run_at_exon2_start(base);
        let blocks = [exons[0].as_str(), exons[1].as_str(), exons[2].as_str()];
        // An intron ENDING in a run of EXACTLY TWO of the run's base, so a 5'
        // shift over the unspliced neighbourhood would travel `c.11-1`,
        // `c.11-2` and no further.
        //
        // Three positions are rewritten, not two, and the third is the point.
        // `intron_leading_with("CG")` fills its tail by alternating `A`/`C`, so
        // index 27 is already `C`: overwriting only 28..30 left the `C` row
        // ending `…A C C C`, a run of three, while the other three bases got the
        // intended two. The zero this sweep reports is unaffected — a longer run
        // only makes an exit easier, so the row still refutes nothing it claimed
        // to — but the controlled variable the sweep advertises was not held on
        // one of its four rows, which is the kind of quiet non-uniformity that
        // makes a later reader trust a comparison that was never made. Pinning
        // index 27 to something the run cannot extend makes all four rows the
        // same experiment.
        let mut intron = intron_leading_with("CG");
        let separator = "ACGT"
            .chars()
            .find(|&candidate| candidate != base)
            .expect("three of the four bases differ from `base`");
        intron.replace_range(27..30, &format!("{separator}{base}{base}"));
        debug_assert_eq!(
            intron.chars().rev().take_while(|&c| c == base).count(),
            2,
            "the intron must end in exactly two `{base}`, or the 5' shift has further to travel \
             on this row than on its siblings"
        );

        for strand in [Strand::Plus, Strand::Minus] {
            let provider = junction_provider(strand, blocks, &intron, true);
            let output = normalize(
                &provider,
                &format!("{CODING}:c.14del"),
                ShuffleDirection::FivePrime,
            );
            if output.contains('-') && !output.contains(":c.-") {
                exits.push(format!("{base}/{strand:?}: {output}"));
            }
            placements.push((base, strand, output));
        }
    }

    // The exit set first, over the completed sweep — this is the property the
    // module doc's "eight placements, zero exits" cites.
    assert!(
        exits.is_empty(),
        "the 5' direction left the transcript: {exits:?}"
    );
    // Non-vacuity: eight is four run alphabets times two strands. Without this
    // the `is_empty()` above would pass just as happily over a sweep that
    // generated nothing at all — and adding it is what caught this test's doc
    // claiming SIXTEEN, a figure copied from the 3' sibling above, which carries
    // a third loop this one does not.
    assert_eq!(
        placements.len(),
        8,
        "expected eight placements (4 run alphabets x 2 strands), got {}",
        placements.len()
    );
    // Only then the per-placement strings.
    for (base, strand, output) in &placements {
        assert_eq!(
            *output,
            format!("{CODING}:c.11del"),
            "5' shuffle, base {base} on the {strand:?} strand, intron ending {base}{base}"
        );
    }
}

// ---------------------------------------------------------------------------
// Reach: which axes and which edit types
// ---------------------------------------------------------------------------

/// **Question.** Is the exit confined to `c.` deletions at the first junction?
///
/// **No.** It reaches every junction that has a following intron, both `del` and
/// `dup`, and the `n.` axis as well as `c.` — `normalize_tx` carries the same
/// gate at `mod.rs:5407` (`#704` sub-problem A) and `normalize_rna` a third at
/// `:6917`.
///
/// Pinned as one test because the fix has to cover all three call sites: a repair
/// applied to `normalize_cds` alone would leave the `n.` half of the corpus's 371
/// exactly where it is.
#[test]
fn the_exit_reaches_the_second_junction_the_dup_spelling_and_the_n_axis() {
    let intron = intron_leading_with("AA");
    let filler: String = std::iter::repeat_n('C', EXON_LEN).collect();

    // A run flush against EXON 2's end: tx 37..=40, which is `c.27`..`c.30`.
    let mut exon2 = filler.clone();
    exon2.replace_range(16..20, "AAAA");
    let second_junction = [filler.as_str(), exon2.as_str(), filler.as_str()];

    // A run flush against exon 1's end, for the `n.` axis: tx 17..=20.
    let mut exon1 = filler.clone();
    exon1.replace_range(16..20, "AAAA");
    let first_junction = [exon1.as_str(), filler.as_str(), filler.as_str()];

    for strand in [Strand::Plus, Strand::Minus] {
        let provider = junction_provider(strand, second_junction, &intron, true);
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.27del")),
            format!("{CONTIG}({CODING}):c.30+2del"),
            "the exon2/exon3 junction crosses identically ({strand:?})"
        );
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.27dup")),
            format!("{CONTIG}({CODING}):c.30+2dup"),
            "a duplication crosses identically to a deletion ({strand:?})"
        );

        // `n.` — no CDS at all, so the axis is the plain transcript numbering and
        // tx 20 is `n.20`.
        let noncoding = junction_provider(strand, first_junction, &intron, false);
        assert_eq!(
            normalize_3prime(&noncoding, &format!("{NONCODING}:n.17del")),
            format!("{CONTIG}({NONCODING}):n.20+2del"),
            "the `n.` axis carries its own copy of the gate, and its own wrapper ({strand:?})"
        );
    }
}

/// **Question.** The corpus exemplar is an *allele*,
/// `s01-c3m-junction-1-del-del-p1-sep0`: `NM_TEST.1:c.[18del;21del]` reaching an
/// output whose first member is intronic. Does the allele add anything to the
/// lone-member mechanism?
///
/// **No — the allele's members each take the path independently.** The member
/// flush against the junction exits; its sibling, which lands mid-exon, does not.
/// So the 371 figure is a count of *members* in this shape, not of a multi-member
/// interaction, and the `member independence` family that
/// `the_cds_end_flush_pair_is_its_two_members_normalized_separately` documents is
/// not implicated here.
///
/// Pinned because it rules out a whole class of fix: nothing about the partition,
/// the merge or the sibling geometry is involved.
#[test]
fn an_allele_exits_member_by_member() {
    let intron = intron_leading_with("AA");
    let filler: String = std::iter::repeat_n('C', EXON_LEN).collect();
    // A run flush against exon 1's end (tx 17..=20 == c.7..=c.10) and a second,
    // independent run wholly inside exon 2 (tx 25..=28 == c.15..=c.18).
    let mut exon1 = filler.clone();
    exon1.replace_range(16..20, "AAAA");
    let mut exon2 = filler.clone();
    exon2.replace_range(4..8, "AAAA");
    let blocks = [exon1.as_str(), exon2.as_str(), filler.as_str()];

    for strand in [Strand::Plus, Strand::Minus] {
        let provider = junction_provider(strand, blocks, &intron, true);
        // Each member on its own.
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.7del")),
            format!("{CONTIG}({CODING}):c.10+2del")
        );
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.15del")),
            format!("{CODING}:c.18del"),
            "the mid-exon member never crosses, so on its own it keeps the bare accession"
        );
        // And the allele: precisely those two answers, side by side — under ONE
        // accession. The wrapper is lifted to the whole description rather than
        // applied per member, because `AlleleVariant`'s compact form
        // (`ACC:c.[a;b]`) requires the members to share an accession: wrapping
        // only the crossing member would drop this into the expanded
        // `[NC_SYNTH.1(NM_TEST.1):c.10+2del;NM_TEST.1:c.18del]`, a much larger
        // representation change than the defect it repairs. A genomic wrapper on
        // an exonic position is unremarkable, so lifting it is the cheap side.
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.[7del;15del]")),
            format!("{CONTIG}({CODING}):c.[10+2del;18del]"),
            "the allele is its two members' independent answers under one accession, so the \
             crossing is a per-member property and the wrapper a whole-description one \
             ({strand:?})"
        );
    }
}

/// **Question.** Is the wrapped output a fixed point, and does ferro's own parser
/// read it back?
///
/// **Both, and both are load-bearing.** `FERRO_ASSERT_IDEMPOTENT` and
/// `FERRO_ASSERT_REPARSE` run over the whole suite in CI, so a wrapper that
/// normalized away on the second pass — or that `parse_hgvs` refused — would be a
/// regression wearing a fix's clothes, and it is the second pass that is the real
/// question: the re-parented output *is* an intronic position on a transcript
/// axis, so it re-enters the very path that produced it.
///
/// **Why this test still exists after the fix.** It was written to record that
/// **no seam oracle covered the corpus's largest defect class** — the output was
/// a fixed point, so `FERRO_ASSERT_IDEMPOTENT` was blind; it re-parsed, so
/// `FERRO_ASSERT_REPARSE` was blind (`parse_hgvs` accepts an intronic position on
/// a bare `NM_`, pinned by
/// `spec_corpus_regressions::a_bare_transcript_accession_accepts_an_intronic_position`);
/// and `FERRO_ASSERT_IN_BOUNDS` has its own blind spot for it
/// (`in_bounds_oracle_scope.rs`). That is still true of the *shape*, which is why
/// the class had to be caught by a census rather than by a seam oracle — and why
/// this file, not an oracle, is where the fix is guarded.
#[test]
fn the_wrapped_output_is_a_fixed_point_and_re_parses() {
    let intron = intron_leading_with("AA");
    let exons = exons_with_run_at_exon1_end('A');
    let blocks = [exons[0].as_str(), exons[1].as_str(), exons[2].as_str()];

    for strand in [Strand::Plus, Strand::Minus] {
        let provider = junction_provider(strand, blocks, &intron, true);
        let once = normalize_3prime(&provider, &format!("{CODING}:c.7del"));
        assert_eq!(once, format!("{CONTIG}({CODING}):c.10+2del"));
        assert_eq!(
            normalize_3prime(&provider, &once),
            once,
            "the wrapped output is its own fixed point on the {strand:?} strand — the wrapper \
             survives a second pass rather than being re-derived or stripped"
        );
        assert!(
            parse_hgvs(&once).is_ok(),
            "and ferro's own parser reads the compound reference back"
        );
    }
}

// ---------------------------------------------------------------------------
// #1704: the rendering, which is the half that was actually wrong
// ---------------------------------------------------------------------------

/// **Question.** `checklist.md:20` says a bare `NM_`/`NR_` "can only be used to
/// describe variants in introns using a `c.` prefix when a genomic reference
/// sequence is given". The junction crossing above is legal
/// (`background/numbering.md:26` withholds `general.md:44`'s exon/exon exception
/// from an exon/intron junction), so the coordinate stays — what must change is
/// the **accession** it is rendered against. Does it?
///
/// **It must, and this is the pin.** The crossing is computed in genomic space
/// against a named contig, so the genomic reference `checklist.md:20` asks for is
/// already in hand at the moment the answer is adopted; the output is rendered as
/// `NC_SYNTH.1(NM_TEST.1):c.10+2del`, the exact compound form the clause names.
///
/// Both axes and both strands, because the fix has to cover all three copies of
/// the `#670` gate.
#[test]
fn a_junction_exit_is_rendered_against_the_genomic_wrapper() {
    let intron = intron_leading_with("AA");
    let exons = exons_with_run_at_exon1_end('A');
    let blocks = [exons[0].as_str(), exons[1].as_str(), exons[2].as_str()];

    for strand in [Strand::Plus, Strand::Minus] {
        assert_eq!(
            normalize_3prime(
                &junction_provider(strand, blocks, &intron, true),
                &format!("{CODING}:c.7del")
            ),
            format!("{CONTIG}({CODING}):c.10+2del"),
            "checklist.md:20 — an intronic position needs the genomic wrapper, and the \
             crossing resolved that contig itself ({strand:?})"
        );
        assert_eq!(
            normalize_3prime(
                &junction_provider(strand, blocks, &intron, false),
                &format!("{NONCODING}:n.17del")
            ),
            format!("{CONTIG}({NONCODING}):n.20+2del"),
            "the same on the n. axis, whose gate is a separate copy ({strand:?})"
        );
    }
}

/// **Question.** The fix adds the wrapper where ferro's own shuffle manufactured
/// the offset. What about an input that **already** names an intronic position on
/// a bare transcript?
///
/// **Left exactly as authored.** That class is settled by the
/// `bare-transcript-intronic-position` ruling — strict input hygiene refuses it
/// (`IntronicOnBareTranscript` / W4007), lenient accepts it as written — and
/// re-spelling it here would overturn one decided record as a side effect of
/// implementing another. It is also what keeps the two halves of `checklist.md:20`
/// answered by one predicate: the same
/// `intronic_on_bare_transcript_warning` that decides the strict refusal decides
/// whether an output needs the wrapper.
#[test]
fn an_authored_intronic_position_is_left_as_authored() {
    let intron = intron_leading_with("AA");
    let exons = exons_with_run_at_exon1_end('A');
    let blocks = [exons[0].as_str(), exons[1].as_str(), exons[2].as_str()];

    for strand in [Strand::Plus, Strand::Minus] {
        let provider = junction_provider(strand, blocks, &intron, true);
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.10+2del")),
            format!("{CODING}:c.10+2del"),
            "an authored intronic position is neither re-parented nor moved ({strand:?})"
        );
    }
}

// ---------------------------------------------------------------------------
// The rendering guards, and the residue they leave
// ---------------------------------------------------------------------------

/// **Question.** `reparent_junction_exit` refuses to build a wrapper from three
/// conditions — trailing input after the accession parse, a
/// `is_valid_compound_outer` rejection, and an outer equal to the bare accession.
/// Does any of them do anything?
///
/// **All three do, but only once the named chromosome also serves sequence** —
/// and that precondition is the whole reason they looked like dead code. Renaming
/// the chromosome alone does not reach them: the `#670` gate fetches genomic
/// sequence *by that name*, so a name the provider cannot serve breaks the
/// crossing itself and the answer falls back to the exon-confined `c.10del`,
/// never reaching the rendering stage. A guard test built that way passes while
/// asserting nothing about the guard, which is why
/// [`junction_provider_on_chromosome`] serves the contig under whatever the
/// transcript names.
///
/// Each case below therefore has a crossing to render (`c.10+2del` is reached)
/// and is declined at the rendering step, returning the pre-`#1704` bare string.
/// That is the documented safe direction, not a repair — see
/// [`a_mixed_allele_still_ships_a_manufactured_offset_bare`] for what it costs.
#[test]
fn each_rendering_guard_declines_rather_than_emitting_a_bad_reference() {
    let intron = intron_leading_with("AA");
    let exons = exons_with_run_at_exon1_end('A');
    let blocks = [exons[0].as_str(), exons[1].as_str(), exons[2].as_str()];
    let declined = format!("{CODING}:c.10+2del");

    // Which conjunct actually catches each case is MEASURED, not assumed —
    // dropping `outer == bare` alone changes nothing *here*, because `bare` is a
    // RefSeq transcript and `is_valid_compound_outer` rejects a transcript outer
    // one conjunct earlier. So the third case below is the pairing rule firing
    // again, and reading it as the self-reference check would credit that
    // conjunct with a decline it does not make on this accession.
    //
    // It is not dead code, though — see
    // [`a_self_referential_wrapper_is_declined`], which reaches it on the one
    // accession class where the pairing rule lets a self-reference through.
    for (chromosome, conjunct) in [
        ("NC_SYNTH.1junk", "trailing input after the accession parse"),
        (
            "NM_OTHER.1",
            "the pairing rule — a known transcript outer, backwards",
        ),
        (
            CODING,
            "the pairing rule again, NOT the `outer == bare` conjunct",
        ),
    ] {
        let provider =
            junction_provider_on_chromosome(Strand::Plus, blocks, &intron, true, Some(chromosome));
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.7del")),
            declined,
            "chromosome {chromosome:?} must be declined by {conjunct}, leaving the bare \
             pre-#1704 string rather than a reference ferro cannot justify"
        );
    }

    // The control: the one usable chromosome renders the wrapper. Without this
    // the three assertions above would all pass on a build that never wraps.
    let provider =
        junction_provider_on_chromosome(Strand::Plus, blocks, &intron, true, Some(CONTIG));
    assert_eq!(
        normalize_3prime(&provider, &format!("{CODING}:c.7del")),
        format!("{CONTIG}({CODING}):c.10+2del"),
        "and a usable genomic accession still wraps, so the three declines above are the \
         guards firing rather than the wrapper being unreachable"
    );
}

/// **Question.** A SAM-style `chr17` must not be emitted as `chr17(NM_…)` —
/// `refseq.md` admits no such reference. What actually stops it?
///
/// **The bare accession parse, NOT the compound-outer rule.** This is worth
/// pinning because the obvious reading is the other way round, and the code
/// comment asserted the other way round until this test was written.
/// `is_valid_compound_outer` deliberately *admits* an unclassifiable custom
/// accession (`inferred_variant_type().is_none()`, #1146) so that a custom or
/// assembly reference can still carry a specification — and the assertion below
/// proves it, because ferro's own parser reads `chr17(NM_TEST.1):c.10+2del` back
/// happily. What excludes it is that `parse_accession("chr17")` fails on the
/// *bare* string, so the guard never gets as far as the pairing rule.
///
/// `GRCh38` is the same story from the other direction: `is_assembly_ref` would
/// admit it, and the bare parse is again what declines.
#[test]
fn a_sam_style_chromosome_is_excluded_by_the_bare_accession_parse() {
    let intron = intron_leading_with("AA");
    let exons = exons_with_run_at_exon1_end('A');
    let blocks = [exons[0].as_str(), exons[1].as_str(), exons[2].as_str()];

    for chromosome in ["chr17", "GRCh38"] {
        assert!(
            parse_hgvs(&format!("{chromosome}({CODING}):c.10+2del")).is_ok(),
            "the compound form is one ferro's own parser accepts, so the pairing rule is NOT \
             what refuses {chromosome:?}"
        );
        let provider =
            junction_provider_on_chromosome(Strand::Plus, blocks, &intron, true, Some(chromosome));
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.7del")),
            format!("{CODING}:c.10+2del"),
            "and ferro still declines to emit it — the bare parse of {chromosome:?} is the \
             mechanism, which is the half a reader would otherwise attribute wrongly"
        );
    }
}

/// **Question.** Is the `outer == bare` conjunct reachable, or is it dead code?
///
/// **Reachable, on exactly one accession class — a bare LRG record.** This is
/// the case the sibling test above cannot make, and the reason that conjunct is
/// kept rather than deleted.
///
/// For every RefSeq accession the guard can be handed, the conjunct is
/// unreachable and the sibling test measures that: `bare` comes from
/// [`bare_transcript_intronic_accession`], `is_valid_compound_outer` rejects a
/// known transcript outer, so the pairing rule fires first and `NM_TEST.1` as its
/// own chromosome is declined one conjunct earlier.
///
/// **LRG is the exception, and it is an accident of two rules meeting.** The
/// `checklist.md:20` predicate admits an LRG reference by *prefix*
/// (`Accession::is_lrg`), which covers the bare genomic record `LRG_<N>` as well
/// as the transcript `LRG_<N>t<M>` — while `is_valid_compound_outer` keys off
/// `inferred_variant_type`, which reads a bare `LRG_<N>` as **genomic** and so
/// admits it as an outer. So for `bare == LRG_5` both earlier conjuncts pass, and
/// only `outer == bare` stands between ferro and `LRG_5(LRG_5):c.10+2del`.
///
/// That string is not caught downstream either: `parse_hgvs` reads it back
/// happily (asserted below), so `FERRO_ASSERT_REPARSE` would never see it. A
/// self-referential wrapper is nonsense that no oracle in this repo detects.
///
/// The provider shape is not contrived. An LRG record *is* its own genomic
/// reference — there is no separate contig for it — so a provider that fills
/// `Transcript::chromosome` with the LRG accession is a reasonable thing to
/// build, which is what makes this worth guarding rather than deleting.
#[test]
fn a_self_referential_wrapper_is_declined() {
    let intron = intron_leading_with("AA");
    let exons = exons_with_run_at_exon1_end('A');
    let blocks = [exons[0].as_str(), exons[1].as_str(), exons[2].as_str()];
    let lrg = "LRG_5";

    // Both earlier conjuncts pass for this accession, so a decline here can only
    // be the third. Stated as assertions rather than as prose, because the whole
    // point of the test is that the first two do NOT fire.
    assert!(
        parse_hgvs(&format!("{lrg}({lrg}):c.10+2del")).is_ok(),
        "ferro's own parser accepts the self-referential wrapper, so neither the bare parse nor \
         the pairing rule refuses it — and no re-parse oracle would catch it either"
    );

    let provider = junction_provider_named(Strand::Plus, blocks, &intron, true, Some(lrg), lrg);
    assert_eq!(
        normalize_3prime(&provider, &format!("{lrg}:c.7del")),
        format!("{lrg}:c.10+2del"),
        "an LRG record naming itself as its chromosome must be declined by `outer == bare`, \
         leaving the bare pre-#1704 string rather than the meaningless `{lrg}({lrg})`"
    );

    // The control: the same LRG record on a real contig still wraps, so the
    // decline above is the self-reference and not LRG being unsupported.
    let hosted = junction_provider_named(Strand::Plus, blocks, &intron, true, Some(CONTIG), lrg);
    assert_eq!(
        normalize_3prime(&hosted, &format!("{lrg}:c.7del")),
        format!("{CONTIG}({lrg}):c.10+2del"),
        "an LRG record hosted on a genomic contig wraps normally, so the decline above is \
         specifically the self-reference"
    );
}

/// **Question.** A transcript with no chromosome at all — same path?
///
/// **No, and the difference matters when reading a failure.** With no name there
/// is nothing to fetch, so the `#670` gate never resolves a crossing and the
/// answer is the exon-confined `c.10del` — the rendering guards are never
/// reached. Every case in
/// [`each_rendering_guard_declines_rather_than_emitting_a_bad_reference`] reaches
/// `c.10+2del` first. Two distinct mechanisms with two distinct outputs, pinned
/// apart so a regression in either cannot be read as the other.
#[test]
fn an_unnamed_chromosome_never_reaches_the_rendering_guards() {
    let intron = intron_leading_with("AA");
    let exons = exons_with_run_at_exon1_end('A');
    let blocks = [exons[0].as_str(), exons[1].as_str(), exons[2].as_str()];
    let provider = junction_provider_on_chromosome(Strand::Plus, blocks, &intron, true, None);
    assert_eq!(
        normalize_3prime(&provider, &format!("{CODING}:c.7del")),
        format!("{CODING}:c.10del"),
        "no chromosome means no genomic fetch, so the crossing never happens and there is no \
         intronic offset to render — a different failure from a declined wrapper"
    );
}

/// **Question — and this one records a DEFECT, not a decision.** An allele mixes
/// a member whose intronic offset ferro *manufactured* with a member the author
/// spelled intronic. What ships?
///
/// **The manufactured offset ships on a bare accession** — the exact class
/// `#1704` exists to close, surviving in the one shape its guard cannot see.
///
/// The cause is a granularity mismatch. Provenance ("did ferro manufacture this
/// offset?") is a property of each **leaf**, but it is asked here at
/// whole-description granularity: `names_bare_transcript_intronic(input)` is an
/// ANY-leaf existence test, so one authored member vetoes the wrapper for every
/// other member, including one ferro moved itself. Compare the lone cases pinned
/// in [`a_junction_exit_is_rendered_against_the_genomic_wrapper`] and
/// [`an_authored_intronic_position_is_left_as_authored`]: each is correct alone,
/// and only the mixture is wrong.
///
/// **This is pinned as it is, deliberately.** Closing it needs per-leaf
/// provenance carried from the site that manufactures the offset, because it
/// cannot be recovered by comparing input and output after the fact —
/// normalization reorders, merges and splits members, so no identity map from
/// output leaf back to input leaf survives. And the rendering question that then
/// arises (the compact form `ACC:c.[a;b]` requires members to share an accession,
/// so a mixed description must either lift the wrapper to the whole description
/// or expand to per-member accessions) is a representation-policy choice that is
/// not this PR's to make. Pinning the current string keeps the residue from
/// changing unnoticed in either direction.
#[test]
fn a_mixed_allele_still_ships_a_manufactured_offset_bare() {
    let intron = intron_leading_with("AA");
    let filler: String = std::iter::repeat_n('C', EXON_LEN).collect();
    // A run flush against exon 1's end (c.7..=c.10), and a second run wholly
    // inside exon 2, so the two members take visibly different paths.
    let mut exon1 = filler.clone();
    exon1.replace_range(16..20, "AAAA");
    let mut exon2 = filler.clone();
    exon2.replace_range(4..8, "AAAA");
    let blocks = [exon1.as_str(), exon2.as_str(), filler.as_str()];

    for strand in [Strand::Plus, Strand::Minus] {
        let provider = junction_provider(strand, blocks, &intron, true);
        // The authored member sits at the SECOND junction, so it cannot collide
        // with the offset the first member manufactures.
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.[7del;30+5del]")),
            format!("{CODING}:c.[10+2del;30+5del]"),
            "RESIDUE: `10+2` is ferro's own, and it ships bare because the authored `30+5` \
             vetoed the wrapper for the whole description ({strand:?})"
        );
        // The same defect on the other side of the same intron, so the pin is on
        // the mixture rather than on one offset's spelling.
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.[7del;11-5del]")),
            format!("{CODING}:c.[10+2del;11-5del]"),
            "and with the authored member spelled from the 3' side of intron 1 ({strand:?})"
        );
        // The negative control: with NO manufactured offset the description is
        // correctly left alone, so the veto is not simply "any allele is skipped".
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.[15del;30+5del]")),
            format!("{CODING}:c.[18del;30+5del]"),
            "no member crosses, so nothing is manufactured and leaving it bare is correct \
             ({strand:?})"
        );
    }
}
