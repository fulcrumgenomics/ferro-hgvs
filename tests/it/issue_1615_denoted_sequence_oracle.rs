//! The denoted-sequence seam oracle, checked against the defects it exists for.
//!
//! `FERRO_ASSERT_SEQUENCE` is the fourth oracle at `normalize_core_checked`'s
//! exit (#1615). The three before it — `FERRO_ASSERT_IDEMPOTENT`,
//! `FERRO_ASSERT_REPARSE`, `FERRO_ASSERT_IN_BOUNDS` — are all *form* questions,
//! and a description denoting entirely different bases can satisfy every one of
//! them at once. Eight defects did exactly that.
//!
//! This file is the oracle's own regression guard. It does **not** normalize:
//! every row pins a `(reference, input, wrong output)` triple recorded on the
//! issue that reported it, and asserts that
//! [`compare_denoted_sequences`](ferro_hgvs::spdi::compare_denoted_sequences) —
//! the predicate `Normalizer::assert_denoted_sequence` wraps — reports the pair
//! as a fault. Pinning the recorded strings rather than re-normalizing is what
//! keeps the guard meaningful after each defect is fixed: a test that ran the
//! normalizer would go green on the fix and stop saying anything about the
//! oracle.
//!
//! The two halves matter equally. The rows below must **fire**; the rows in
//! `respellings_that_denote_the_same_bases` must **not**, because an oracle that
//! fires on ordinary correct normalization is one nobody can leave enabled.

use crate::common::cis_apply_oracle::provider as template_provider;
use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::spdi::{compare_denoted_sequences, DenotedSequenceComparison, NotComparable};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Which reference a recorded defect was found on.
enum Reference {
    /// A bare contig under the accession `TEMPLATE`, as the cis-allele tests
    /// build it.
    Template(&'static str),
    /// A core wrapped in 256 bases of `ACGT` padding under `NC_TEST.1`, so core
    /// base 1 sits at `g.257`. The shape every `SyntheticBuilder`-backed issue
    /// reproducer uses.
    PaddedCore(&'static str),
}

impl Reference {
    fn build(&self) -> MockProvider {
        match self {
            Self::Template(sequence) => template_provider(sequence),
            Self::PaddedCore(core) => SyntheticBuilder::genomic(core).build(),
        }
    }
}

/// A defect this oracle is answerable for: the reference it was found on, the
/// input, and the description the normalizer emitted for it at the time.
struct Recorded {
    issue: &'static str,
    reference: Reference,
    input: &'static str,
    /// The output as recorded on the issue. **Not** what `main` emits today for
    /// the six that are fixed — see the module docs for why that is the point.
    wrong_output: &'static str,
}

/// The eight defects named on #1615, each with the wrong output its own issue
/// records.
///
/// Six were found by hand, one shape at a time, and fixed and regression-tested
/// separately (#1254, #1281, #1290, #1304, #1308, #1312). #1592 and #1600 were
/// found by a proptest run past its default case budget and are open at the time
/// of writing. All eight are the same class, and none of the three older oracles
/// sees any of them.
const RECORDED_DEFECTS: &[Recorded] = &[
    Recorded {
        issue: "#1254",
        reference: Reference::Template("TAATATATATAATATATATT"),
        input: "TEMPLATE:g.[3_4del;9del]",
        // The deletion 3'-shifted straight past its sibling, the two merged,
        // and the merged member shifted again.
        wrong_output: "TEMPLATE:g.12_14del",
    },
    Recorded {
        issue: "#1281",
        reference: Reference::Template("TTTTTTTTTAATATATTTTAATAT"),
        input: "TEMPLATE:g.[1del;9_10delinsA]",
        // Two members claiming base 1. This one denotes no sequence at all,
        // which is a *different* and more severe verdict — see
        // `an_output_that_contradicts_itself_is_not_a_skip`.
        wrong_output: "TEMPLATE:g.[1del;1del]",
    },
    Recorded {
        issue: "#1290",
        reference: Reference::PaddedCore("ATACAGAAAATCAGGGCATA"),
        input: "NC_TEST.1:g.[263_264insA;265_266insC]",
        wrong_output: "NC_TEST.1:g.[265_266insC;266dup]",
    },
    Recorded {
        issue: "#1304",
        reference: Reference::PaddedCore("GCATGAAAAT"),
        input: "NC_TEST.1:g.[260_261insGA;261_262insA;264del]",
        wrong_output: "NC_TEST.1:g.[261_262insA;261_262dup;265del]",
    },
    Recorded {
        issue: "#1308",
        reference: Reference::PaddedCore("CAGAAGATGAATAA"),
        input: "NC_TEST.1:g.[263_264insTG;264_265insTG]",
        wrong_output: "NC_TEST.1:g.[264_265insTG;264_265dup]",
    },
    Recorded {
        issue: "#1312",
        reference: Reference::PaddedCore("TAAAACCA"),
        input: "NC_TEST.1:g.[260_261insAC;261_262insAC]",
        wrong_output: "NC_TEST.1:g.261_262insACCA",
    },
    Recorded {
        issue: "#1592",
        reference: Reference::PaddedCore("TCCAGCAGATAT"),
        input: "NC_TEST.1:g.[262_263insAG;264G>T;266T>A]",
        // The repeat names the lone `A` at 265 rather than the one at 267.
        wrong_output: "NC_TEST.1:g.265A[3]",
    },
    Recorded {
        issue: "#1600",
        reference: Reference::PaddedCore("GAACAGCAGAAGCGA"),
        input: "NC_TEST.1:g.261_267delinsGCAGAAAC",
        wrong_output: "NC_TEST.1:g.[261del;266_267insCA]",
    },
];

/// The oracle's acceptance criterion: it fires on every defect #1615 names.
///
/// "Fires" is either verdict the seam panics on — a changed sequence, or an
/// output that denotes none while the input does. A `NotComparable` or an
/// `Agree` here would be a miss, and the message says which so a regression
/// reads as itself rather than as a generic failure.
#[test]
fn the_oracle_fires_on_every_recorded_defect() {
    let mut missed = Vec::new();
    for defect in RECORDED_DEFECTS {
        let provider = defect.reference.build();
        let input = parse_hgvs(defect.input).expect("input must parse");
        let output = parse_hgvs(defect.wrong_output).expect("recorded output must parse");
        match compare_denoted_sequences(&input, &output, &provider) {
            DenotedSequenceComparison::Differ { .. }
            | DenotedSequenceComparison::OutputContradictsItself => {}
            DenotedSequenceComparison::Agree => missed.push(format!(
                "{}: the oracle sees no difference between `{}` and `{}`",
                defect.issue, defect.input, defect.wrong_output
            )),
            DenotedSequenceComparison::NotComparable(reason) => missed.push(format!(
                "{}: the oracle declined to compare `{}` against `{}` — {reason}",
                defect.issue, defect.input, defect.wrong_output
            )),
        }
    }
    assert!(
        missed.is_empty(),
        "the denoted-sequence oracle must fire on every defect #1615 names, but missed \
         {}/{}:\n  {}",
        missed.len(),
        RECORDED_DEFECTS.len(),
        missed.join("\n  "),
    );
}

/// #1281's verdict is specifically the severe one, not merely "different".
///
/// `g.[1del;1del]` has no resulting sequence — two members claim base 1, so
/// applying them depends on an order the description does not state. Reporting
/// that as a skip would file the worse defect under the milder one's exemption,
/// which is the failure mode a sequence oracle exists to remove.
#[test]
fn an_output_that_contradicts_itself_is_not_a_skip() {
    let provider = template_provider("TTTTTTTTTAATATATTTTAATAT");
    let input = parse_hgvs("TEMPLATE:g.[1del;9_10delinsA]").expect("parse");
    let output = parse_hgvs("TEMPLATE:g.[1del;1del]").expect("parse");
    assert_eq!(
        compare_denoted_sequences(&input, &output, &provider),
        DenotedSequenceComparison::OutputContradictsItself,
    );
}

/// The negative control at unit scale: legitimate re-spellings must not fire.
///
/// Every row here is a pair the normalizer is *allowed* to move between — a
/// 3'-shift through a tract, a merge of adjacent members, a decomposition, a
/// repeat spelling of a duplication. The union window is what makes the first
/// kind work: `g.3_4del` and `g.7_8del` denote the same bases over any window
/// containing both, and a per-description window would give each its own frame
/// and report a difference that is not there.
#[test]
fn respellings_that_denote_the_same_bases_do_not_fire() {
    let rows: &[(Reference, &str, &str)] = &[
        // #1254's correct answer, against its own input: a clamped 3'-shift plus
        // a merge of the two now-adjacent deletions.
        (
            Reference::Template("TAATATATATAATATATATT"),
            "TEMPLATE:g.[3_4del;9del]",
            "TEMPLATE:g.7_9del",
        ),
        // A pure 3'-shift with no merge, across a window neither description
        // covers alone.
        (
            Reference::Template("TAATATATATAATATATATT"),
            "TEMPLATE:g.3_4del",
            "TEMPLATE:g.7_8del",
        ),
        // A spanning delins and its full decomposition — #1159's shape.
        (
            Reference::PaddedCore("GCATGAAAAT"),
            "NC_TEST.1:g.257_258delinsAA",
            "NC_TEST.1:g.[257G>A;258C>A]",
        ),
        // A duplication and the insertion that denotes it.
        (
            Reference::PaddedCore("GCATGAAAAT"),
            "NC_TEST.1:g.261_262insA",
            "NC_TEST.1:g.262dup",
        ),
        // #1290's correct answer: two insertions that merge into one.
        (
            Reference::PaddedCore("ATACAGAAAATCAGGGCATA"),
            "NC_TEST.1:g.[263_264insA;265_266insC]",
            "NC_TEST.1:g.266_267insCA",
        ),
        // #1600's correct answer, which the bracketed spelling already reaches.
        (
            Reference::PaddedCore("GAACAGCAGAAGCGA"),
            "NC_TEST.1:g.261_267delinsGCAGAAAC",
            "NC_TEST.1:g.[261del;267_268insAC]",
        ),
    ];
    for (reference, input, output) in rows {
        let provider = reference.build();
        let input_variant = parse_hgvs(input).expect("parse");
        let output_variant = parse_hgvs(output).expect("parse");
        assert_eq!(
            compare_denoted_sequences(&input_variant, &output_variant, &provider),
            DenotedSequenceComparison::Agree,
            "`{input}` and `{output}` denote the same bases; the oracle must stay silent",
        );
    }
}

/// An input that denotes no single sequence is a *counted skip*, not a pass.
///
/// There is no baseline to compare an output against, so the comparison
/// declines — the same discipline `assert_reparseable` and `assert_in_bounds`
/// apply to inputs that were already broken. What matters is that the reason is
/// reported rather than collapsed into silence, so a run can say how much it
/// actually compared.
#[test]
fn an_input_that_denotes_no_sequence_is_reported_as_such() {
    let provider = template_provider("TAATATATATAATATATATT");
    // A trans allele describes two molecules, so there is no one resulting
    // sequence for either side.
    let input = parse_hgvs("TEMPLATE:g.[3_4del];[9del]").expect("parse");
    let output = parse_hgvs("TEMPLATE:g.[7_8del];[9del]").expect("parse");
    assert_eq!(
        compare_denoted_sequences(&input, &output, &provider),
        DenotedSequenceComparison::NotComparable(NotComparable::InputDenotesNoSequence),
    );
}

/// Two descriptions on different accessions are not comparable at all.
///
/// Normalization can legitimately substitute an accession (#785's transcript
/// version selection), and comparing bases across two sequences would report a
/// difference that says nothing about the normalization.
#[test]
fn a_changed_accession_is_reported_rather_than_compared() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_ONE.1", "GGATTACAGGCATTAGCCTGA");
    provider.add_genomic_sequence("NC_TWO.1", "TTTTTTTTTTTTTTTTTTTTT");
    let input = parse_hgvs("NC_ONE.1:g.5del").expect("parse");
    let output = parse_hgvs("NC_TWO.1:g.5del").expect("parse");
    assert_eq!(
        compare_denoted_sequences(&input, &output, &provider),
        DenotedSequenceComparison::NotComparable(NotComparable::AccessionChanged),
    );
}

/// The comparison reports the bases, not just a verdict.
///
/// The seam's panic message quotes all three strings, because "the sequence
/// changed" is not actionable on its own — #1592 and #1600 were both diagnosed
/// from the reference/input/output triple, and a bare boolean would have made
/// each of them a fresh investigation.
#[test]
fn a_difference_carries_the_three_sequences() {
    let provider = template_provider("TAATATATATAATATATATT");
    let input = parse_hgvs("TEMPLATE:g.[3_4del;9del]").expect("parse");
    let output = parse_hgvs("TEMPLATE:g.12_14del").expect("parse");
    let DenotedSequenceComparison::Differ {
        accession,
        start,
        reference,
        from_input,
        from_output,
    } = compare_denoted_sequences(&input, &output, &provider)
    else {
        panic!("#1254's recorded pair must differ");
    };
    assert_eq!(accession, "TEMPLATE");
    // The window is the union of both descriptions' spans: interbase 2 (the
    // input's `3_4del`) through 14 (the output's `12_14del`).
    assert_eq!(start, 2);
    assert_eq!(reference, "ATATATATAATA");
    assert_eq!(from_input, "ATATTAATA");
    assert_eq!(from_output, "ATATATATA");
    assert_ne!(from_input, from_output);
}

/// The measured false-positive classes, each pinned to the verdict that closed
/// it.
///
/// This is the half of the oracle that took the work. Its first run over the
/// suite raised **344** fires, of which 16 were real: an oracle that cannot tell
/// "the applier could not derive this" from "the normalizer got this wrong" is
/// not usable, and every row below is a class — not a case — that was
/// misclassified until the verdict beside it existed.
///
/// The parenthesised counts below are each measured on the run that isolated that
/// class, after the earlier ones were closed — they overlap and sum well above
/// 344, so they are not a partition of it. `docs/ORACLES.md`'s table names the classes
/// without their counts; the caveat lives here.
///
/// **Which of the seven measured classes are here, and which are not** — the
/// name says "the measured false-positive classes", so say where each one
/// actually lives rather than letting the name imply all seven:
///
/// | class | covered by |
/// |---|---|
/// | untransliterable output (328) | two rows below |
/// | insertion flush against a deletion (233) | `an_insertion_flush_against_a_deletion_is_not_an_overlap` |
/// | overlap-conflicting input (~12) | a row below |
/// | `pter`/`qter` (6) | a row below |
/// | uncertain `[(…)]` allele (1) | a row below |
/// | `r.` payload vs a DNA reference (1) | `an_rna_description_over_a_dna_reference_does_not_fire` |
/// | corrected `REFSEQ_MISMATCH` (5) | **not covered here, and not coverable here** — it is skipped in `Normalizer::assert_denoted_sequence`, which reads the warning list, not in this predicate, so no input to `compare_denoted_sequences` can reach it |
///
/// The unifying mistake was reading a *derivation failure* as a *fault*. The two
/// sides of the comparison do not need the reference equally: a description that
/// states its own bases converts with no provider, and the canonical form of the
/// same variant usually does not.
#[test]
fn the_measured_false_positive_classes_stay_silent() {
    // A provider holding no bases at all, which is what most of the suite
    // normalizes against.
    let bare = MockProvider::new();
    let tract = template_provider("TAATATATATAATATATATT");

    let rows: &[(&MockProvider, &str, &str, NotComparable)] = &[
        // 328 fires. The input states its deleted bases and needs no reference;
        // the output does not state them and cannot be derived without one.
        (
            &bare,
            "NC_000001.11:g.[1000G>A;1001A>C]",
            "NC_000001.11:g.1000_1001delinsAC",
            NotComparable::OutputUntransliterable,
        ),
        (
            &bare,
            "NC_000007.13:g.21641025delT",
            "NC_000007.13:g.21641025del",
            NotComparable::OutputUntransliterable,
        ),
        // 6 fires. `pter` carries no numeric coordinate, and `hgvs_to_spdi`
        // reads the `base` field it leaves at 0.
        (
            &tract,
            "TEMPLATE:g.pter_5del",
            "TEMPLATE:g.1_5del",
            NotComparable::UnresolvableSpecialPosition,
        ),
        // 1 fire. The members of a predicted allele are uncertain by
        // construction, so the normalizer does not clamp them and the output may
        // legitimately still overlap.
        (
            &tract,
            "TEMPLATE:g.[(3_4del;9del)]",
            "TEMPLATE:g.[(7_9del;9del)]",
            NotComparable::UncertainAllele,
        ),
        // ~12 fires. An insertion *interior* to a deletion: the input itself
        // denotes nothing, so there is no baseline. This is the
        // overlap-conflicting allele the repository declines to canonicalize.
        (
            &tract,
            "TEMPLATE:g.[2_4del;2_3insAA]",
            "TEMPLATE:g.2_4delinsAA",
            NotComparable::InputDenotesNoSequence,
        ),
    ];

    for (provider, input, output, expected) in rows {
        let input_variant = parse_hgvs(input).expect("parse");
        let output_variant = parse_hgvs(output).expect("parse");
        assert_eq!(
            compare_denoted_sequences(&input_variant, &output_variant, *provider),
            DenotedSequenceComparison::NotComparable(*expected),
            "`{input}` -> `{output}` must decline as {expected}",
        );
    }
}

/// A pure insertion flush against the 5' end of a deletion is **not** an
/// overlap, and the pair denotes one well-defined sequence.
///
/// 233 of the first run's fires were this shape, `g.[…;insA;A[5]]`-flavoured, and
/// they came from `apply_triples`' disjointness guard rather than from anything
/// the descriptions did. `splice_denoted_sequence` uses `cis_apply_oracle`'s
/// walk instead: longer deletions first among members at one position, so the
/// span-claiming member is applied before the zero-length one flush against it.
///
/// The companion case — an insertion *interior* to a deletion — is a genuine
/// overlap and is covered above; keeping both pinned is what stops a fix to
/// either from quietly becoming a fix to the other.
/// Two assertions, because they fail for different reasons and only the second
/// can see the walk order. Against itself the pair only has to *apply at all* —
/// that is the 233-fire regression, since a disjointness refusal makes both sides
/// denote nothing and the verdict becomes `InputDenotesNoSequence` rather than
/// `Agree`. But identical strings cannot distinguish "the members were applied in
/// the right order" from "both sides made the same mistake": a walk that dropped
/// the insertion, or applied it after the deletion had consumed base 4, would
/// still agree with itself. The re-spelling closes that — `g.4_5delinsCC` states
/// the resulting bases in one member, with no ordering left to get wrong, so it
/// agrees with the two-member form only if the walk really produced `TAA` + `CC` +
/// base 6 onwards.
#[test]
fn an_insertion_flush_against_a_deletion_is_not_an_overlap() {
    let provider = template_provider("TAATATATATAATATATATT");
    // `3_4insCC` sits at the interbase 5' of `4_5del`'s first base.
    let flush_pair = parse_hgvs("TEMPLATE:g.[3_4insCC;4_5del]").expect("parse");
    assert_eq!(
        compare_denoted_sequences(&flush_pair, &flush_pair, &provider),
        DenotedSequenceComparison::Agree,
        "the flush pair must apply at all; a disjointness refusal reports it as \
         denoting no sequence, which is the 233-fire class",
    );

    // Bases 1-5 are `TAATA`, so deleting `4_5` and inserting `CC` leaves
    // `TAA` + `CC` + base 6 onwards — the same bases the two-member spelling
    // denotes, in a single member that fixes no order.
    let respelled = parse_hgvs("TEMPLATE:g.4_5delinsCC").expect("parse");
    assert_eq!(
        compare_denoted_sequences(&flush_pair, &respelled, &provider),
        DenotedSequenceComparison::Agree,
        "`g.[3_4insCC;4_5del]` and `g.4_5delinsCC` denote the same bases; a walk \
         that mis-orders the flush members disagrees with the single-member form",
    );
}

/// The seventh measured class: an `r.` description against a DNA reference.
///
/// One fire, and the only thing that closes it is the `U`-to-`T` fold in
/// `canonical_base`. `hgvs_to_spdi` transliterates the `r.` axis to RNA while the
/// window it is spliced into is served as DNA, so every `r.` comparison rests on
/// that fold — without it these pairs read as `UGC` against `TGC` and the
/// comparison declines rather than agreeing.
///
/// `SyntheticBuilder::noncoding` is what creates the asymmetry: it stores the
/// transcript's bases as **DNA** under `NR_TEST.1`, which is what a real prepared
/// reference does too.
#[test]
fn an_rna_description_over_a_dna_reference_does_not_fire() {
    // Transcript positions 1..20; `r.6_8` is `TGC` in the served DNA, `ugc` on
    // the RNA axis. Positions 6..12 carry no tract, so neither spelling shifts.
    let provider = SyntheticBuilder::noncoding("ACGTGTGCAAGTACGTTTTT", Strand::Plus).build();
    let rows: &[(&str, &str)] = &[
        // The pair `canonical_base`'s doc comment names: a stated RNA payload
        // against the same duplication read off the DNA reference.
        ("NR_TEST.1:r.6_8dupugc", "NR_TEST.1:r.6_8dup"),
        // And the same question with no stated payload at all, so the fold is
        // reached through the deleted-bases check rather than the final compare.
        ("NR_TEST.1:r.6_8delugc", "NR_TEST.1:r.6_8del"),
    ];
    for (input, output) in rows {
        let input_variant = parse_hgvs(input).expect("parse");
        let output_variant = parse_hgvs(output).expect("parse");
        assert_eq!(
            compare_denoted_sequences(&input_variant, &output_variant, &provider),
            DenotedSequenceComparison::Agree,
            "`{input}` and `{output}` denote the same bases in two alphabets; the oracle \
             must stay silent",
        );
    }
}

/// The seam is wired, and a run with the flag on really compares something.
///
/// Everything else in this file exercises `compare_denoted_sequences` directly,
/// so **deleting the `assert_denoted_sequence` call from
/// `Normalizer::assert_seam_oracles` would leave this whole suite green** — and
/// every gating job that arms the flag green with zero comparisons made. That
/// is the failure mode the counters exist for; this is the assertion that reads
/// them, rather than leaving it to a human to remember to.
///
/// Those jobs are `sweeps`, which appends this module to its selection for
/// exactly this reason, and — as of #1815 — `test-oracle`, whose own selection
/// already contains this module (nothing in its `-E` negates it). So the counter
/// guard is live in two gating jobs rather than one; read that as defence in
/// depth, not as a second copy to remove. This paragraph used to call `sweeps`
/// "the flag's only gating home", which was true when it was written.
///
/// `compared` is monotonic and this test contributes to it itself, so a `> 0`
/// assertion is safe under `cargo test`'s shared process as well as under
/// nextest's process-per-test — no ordering between tests is assumed.
///
/// Gated on the flag because the oracle is: with `FERRO_ASSERT_SEQUENCE` unset,
/// `assert_denoted_sequence` returns before touching either counter, so there is
/// nothing to assert. The skip says so out loud rather than passing quietly.
#[cfg(debug_assertions)]
#[test]
fn the_seam_compares_when_the_flag_is_set() {
    let enabled = std::env::var("FERRO_ASSERT_SEQUENCE")
        .map(|v| !v.is_empty() && v != "0")
        .unwrap_or(false);
    if !enabled {
        eprintln!(
            "skipping: FERRO_ASSERT_SEQUENCE is unset, so the seam oracle is off and its \
             counters stay at zero by design. CI's `sweeps` job sets it."
        );
        return;
    }

    let provider = template_provider("TAATATATATAATATATATT");
    let variant = parse_hgvs("TEMPLATE:g.3_4del").expect("parse");
    Normalizer::new(provider)
        .normalize(&variant)
        .expect("normalize");

    let (compared, skipped) = ferro_hgvs::normalize::denoted_sequence_oracle_counts();
    assert!(
        compared > 0,
        "FERRO_ASSERT_SEQUENCE is set and a normalization just ran, but the denoted-sequence \
         oracle reports {compared} comparisons ({skipped} skipped). Either the call to \
         `assert_denoted_sequence` is no longer on `assert_seam_oracles`' path, or every \
         normalization is declining — both make a green run meaningless."
    );
}
