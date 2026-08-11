//! The mechanism behind the corpus's largest class: a `c.`/`n.` output naming an
//! intronic position its input did not.
//!
//! # What this file pins, and what it deliberately does not
//!
//! `spec_conformance_axis` measures `outputs_leaving_the_transcript` at **371**
//! on the 3' direction and **0** on the 5'; `spec_corpus_regressions`'s
//! `a_minus_strand_junction_shift_leaves_the_transcript` names one row of it. Both
//! describe the *symptom*. These tests pin the **mechanism**, so that a fix can be
//! judged against the cause rather than against one string.
//!
//! Nothing here is fixed. The fix locus is `src/normalize/`, and this is a
//! characterization PR.
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
    provider.add_genomic_sequence(CONTIG, contig);
    provider.add_transcript(Transcript::new(
        if coding { CODING } else { NONCODING }.to_string(),
        Some("SYNTH".to_string()),
        strand,
        tx,
        coding.then_some(CDS.0),
        coding.then_some(CDS.1),
        exon_records,
        Some(CONTIG.to_string()),
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
                format!("{CODING}:c.10+2del"),
                "PINNED DEFECT — base {base} on the {strand:?} strand, intron {continuing}. \
                 checklist.md:20 makes an intronic position inexpressible against a bare NM_. \
                 The SHIFT is correct (numbering.md:26 withholds the exon/exon exception from \
                 exon/intron junctions); the ACCESSION is not."
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
    for (lead, expected) in [
        ("CC", "c.10del"),
        ("AC", "c.10+1del"),
        ("AAC", "c.10+2del"),
        ("AAAC", "c.10+3del"),
    ] {
        let intron = intron_leading_with(lead);
        for strand in [Strand::Plus, Strand::Minus] {
            let provider = junction_provider(strand, blocks, &intron, true);
            assert_eq!(
                normalize_3prime(&provider, &format!("{CODING}:c.7del")),
                format!("{CODING}:{expected}"),
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
            format!("{CODING}:c.10+2del"),
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
            format!("{CODING}:c.30+2del"),
            "the exon2/exon3 junction leaks identically ({strand:?})"
        );
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.27dup")),
            format!("{CODING}:c.30+2dup"),
            "a duplication leaks identically to a deletion ({strand:?})"
        );

        // `n.` — no CDS at all, so the axis is the plain transcript numbering and
        // tx 20 is `n.20`.
        let noncoding = junction_provider(strand, first_junction, &intron, false);
        assert_eq!(
            normalize_3prime(&noncoding, &format!("{NONCODING}:n.17del")),
            format!("{NONCODING}:n.20+2del"),
            "the `n.` axis carries its own copy of the gate ({strand:?})"
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
            format!("{CODING}:c.10+2del")
        );
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.15del")),
            format!("{CODING}:c.18del")
        );
        // And the allele: precisely those two answers, side by side.
        assert_eq!(
            normalize_3prime(&provider, &format!("{CODING}:c.[7del;15del]")),
            format!("{CODING}:c.[10+2del;18del]"),
            "PINNED DEFECT — the allele is its two members' independent answers, so the exit \
             is a per-member property ({strand:?})"
        );
    }
}

/// **Question.** Is the intronic output at least a fixed point, so that the
/// defect is one description rather than an oscillation?
///
/// **Yes, and that is the bad news.** Re-normalizing `c.10+2del` returns it
/// unchanged, which means `FERRO_ASSERT_IDEMPOTENT` cannot see this class either —
/// the same blind spot `FERRO_ASSERT_IN_BOUNDS` has for it
/// (`in_bounds_oracle_scope.rs`). Two of the three seam oracles are silent here
/// for two independent reasons, and the third (`FERRO_ASSERT_REPARSE`) is silent
/// because `parse_hgvs` accepts an intronic position on a bare `NM_` (pinned by
/// `spec_corpus_regressions::a_bare_transcript_accession_accepts_an_intronic_position`).
///
/// So **no seam oracle covers the corpus's largest defect class**, and this test
/// is where that is stated as a measurement rather than as an inference.
#[test]
fn the_intronic_output_is_a_fixed_point_so_no_seam_oracle_sees_it() {
    let intron = intron_leading_with("AA");
    let exons = exons_with_run_at_exon1_end('A');
    let blocks = [exons[0].as_str(), exons[1].as_str(), exons[2].as_str()];

    for strand in [Strand::Plus, Strand::Minus] {
        let provider = junction_provider(strand, blocks, &intron, true);
        let once = normalize_3prime(&provider, &format!("{CODING}:c.7del"));
        assert_eq!(once, format!("{CODING}:c.10+2del"));
        assert_eq!(
            normalize_3prime(&provider, &once),
            once,
            "the intronic output is its own fixed point on the {strand:?} strand, so the \
             idempotency oracle is blind to it too"
        );
        assert!(
            parse_hgvs(&once).is_ok(),
            "and it re-parses, so the re-parse oracle is blind to it as well"
        );
    }
}
