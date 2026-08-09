//! Two validity items: a conflict lenient mode drops, and a projection split.
//!
//! Both are descriptions ferro should refuse and does not, and both are rank-1
//! validity failures — "the output re-parses, denotes a sequence, and violates no
//! absolute prohibition" is rank 1 of the operator's precedence order
//! (`rulings[adjudication-precedence-order]`, tabulated at
//! `tests/it/spec_conformance_axis.rs:20`). They are kept in one module because
//! the *shape* of each defect is the same: a stage detects something and the
//! finding does not survive to the caller.
//!
//! # Item 1 — the conflict is detected and then dropped at the API seam
//!
//! An allele whose members claim the same reference territory denotes no
//! sequence. Ferro knows this: `detect_overlap_conflicts` fires on the raw
//! members, `ErrorMode::Strict` promotes the W5002 warning to an error, and
//! lenient mode hands the description back exactly as authored (#395/#486/#1004,
//! and #1406's gate at `src/normalize/mod.rs:3450`). All of that is correct and
//! is pinned below as such.
//!
//! What is **not** correct is that the finding never reaches the caller.
//! `Normalizer::normalize` returns `Result<HgvsVariant, _>` — one value, no
//! warning channel (`src/normalize/mod.rs:1764`) — so a lenient caller receives
//! `Ok(<a description denoting no sequence>)` and nothing else. Measured through
//! the CLI, `--error-mode lenient` and `--error-mode silent` are byte-identical:
//! the description on stdout, empty stderr, exit 0. Yet `ErrorMode` itself says
//! they differ — `emits_warnings()` is `true` for `Lenient` and `false` for
//! `Silent` (`src/error_handling/types.rs:46-48`), and `Lenient`'s own doc says
//! "warnings will be generated" (`:20-24`).
//!
//! So the laundering is not in the normalizer's algorithm. It is at the seam:
//! **the mode that promises a diagnostic has no channel to deliver one through
//! the entry point that returns the output.**
//!
//! ## The tension, stated rather than glossed
//!
//! Validity is rank 1 and is never traded. Lenient mode exists precisely to keep
//! going on bad input. Refusing honours the first and abolishes the second;
//! silently accepting honours the second and abolishes the first.
//!
//! **The adjudication: pass the input through, and report.** The precedence order
//! ranks candidate outputs a normalizer *chose*; preservation is a decline to
//! choose, not a choice. Handing back an invalid description is not ferro
//! asserting the description is valid — it becomes exactly that once the
//! invalidity goes unreported, because then no caller can tell the two apart.
//! The diagnostic is therefore not a nicety attached to the pass-through; it is
//! what makes the pass-through a decline rather than an endorsement. Refusal
//! stays strict mode's job, and strict mode already does it.
//!
//! # Item 2 — a single-member description projects into several members
//!
//! Same shape one stage over: a rule fires on one axis and does not survive to
//! the others. `LRG_199t1:c.145_147delinsTGG` — the spec's own worked example at
//! `DNA/delins.md:37` — projects as one member on `c.`, `n.`, `r.` and `p.`, and
//! as `LRG_199:g.[494841C>T;494843C>G]` on the genomic axis. Those two genomic
//! positions are `c.145` and `c.147`, so that is exactly the description
//! `DNA/delins.md:42` renders `class="invalid"` and calls "not correct".
//!
//! The cause is that the codon-frame exception (`DNA/delins.md:18`, decided for
//! ferro in `rulings[delins-codon-carve-out-gap-one]`) is evaluated against the
//! axis being rendered rather than against the variant. The `g.` axis has no
//! reading frame, so `general.md:34`'s separation floor of one applies unopposed.
//! None of the cis-allele partitioning work is implicated: the input has no
//! members to partition.
//!
//! The repository already counts the *shape* — `Status::ProjectionSplitsSingleMember`,
//! budgeted at 9 — but the shape is not by itself a defect, and the section
//! header below says which of the nine rows are and which are not.

use ferro_hgvs::conformance::spec_corpus::{
    corpus_cores, denotation_of, Denotation, Frame, RefShape, DENSE_CORE_LEN,
};
use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

/// The `AT`-alphabet core the corpus draws first, as `spec_corpus_regressions`
/// uses it. Long homopolymer runs, which is where an ambiguous shift lives.
fn at_core() -> String {
    corpus_cores(1, DENSE_CORE_LEN).remove(0)
}

/// The three conflicting geometries `spec_corpus_regressions::
/// conflicting_member_geometries_are_normalized_instead_of_refused` names, each
/// of which the corpus verifies denotes no sequence.
///
/// Repeated here rather than imported so this module states its own premise: a
/// row that quietly stopped being a conflict would otherwise satisfy every
/// assertion below while pinning nothing.
const CONFLICTING_ALLELES: [(&str, &str); 3] = [
    ("nested", "NM_TEST.1:c.[9_12del;10_11del]"),
    ("overlapping", "NM_TEST.1:c.[9_12del;11_14del]"),
    ("coincident insertions", "NM_TEST.1:c.[8_9insAC;8_9insCC]"),
];

/// A normalizer over `frame` in the given error mode.
fn normalizer(frame: &Frame, mode: ErrorMode) -> Normalizer<MockProvider> {
    Normalizer::with_config(
        frame.provider().clone(),
        NormalizeConfig::default().with_error_mode(mode),
    )
}

/// Normalize through the plain entry point — the one that returns only a value.
fn normalize_in(frame: &Frame, mode: ErrorMode, input: &str) -> Result<String, String> {
    normalizer(frame, mode)
        .normalize(&parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}")))
        .map(|v| v.to_string())
        .map_err(|e| e.to_string())
}

/// The W5002 diagnostics `input` produces in `mode`, via the entry point that
/// has a channel for them.
fn conflict_diagnostics(frame: &Frame, mode: ErrorMode, input: &str) -> usize {
    normalizer(frame, mode)
        .normalize_with_diagnostics(
            &parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}")),
        )
        .expect("a non-strict mode must not reject")
        .warnings
        .iter()
        .filter(|w| {
            matches!(
                w,
                ferro_hgvs::normalize::NormalizationWarning::OverlapConflict { .. }
            )
        })
        .count()
}

// ---------------------------------------------------------------------------
// Item 1: what lenient mode does with an allele that denotes no sequence
// ---------------------------------------------------------------------------

/// **Question.** What does lenient mode emit for an allele whose members claim
/// one another's territory, and does the output denote anything?
///
/// **It emits the input, verbatim, and the output denotes nothing** — the same
/// nothing the input denoted. Measured on all three geometries the corpus builds:
/// input and output are byte-identical and both classify `Denotation::NoSequence`.
///
/// The relevant clause is `general.md:58` — "descriptions removing part of a
/// reference sequence and replacing it with part of the same sequence are not
/// allowed (e.g., `NM_004006.2:c.[762_768del;767_774dup]`)" — and it is worth
/// being exact about how far it reaches: it names the `del`+`dup` spelling, which
/// ferro already rejects at parse. It does not name `del`+`del` or two insertions
/// at one interbase. The plainer statement, the one that does reach all three, is
/// that a description must denote one sequence and these denote none.
///
/// **ADJUDICATED CORRECT, in part.** Emitting the input rather than a repaired
/// form is ferro's standing answer since #395/#486/#1004 and is right: the
/// independent applier declines an overlapping description outright, so there are
/// no "input's bases" a repair could preserve — any output would come from
/// *choosing* an apply order, which is what W5002 refuses to do. The half that is
/// not right is pinned in the next test.
#[test]
fn lenient_hands_back_a_conflicting_allele_that_denotes_no_sequence() {
    let core = at_core();
    let frame = Frame::build(RefShape::CodingSingleExon, &core);
    let served = frame.served();

    for (geometry, input) in CONFLICTING_ALLELES {
        // The premise. Without it this test asserts a preference, not a defect.
        assert_eq!(
            denotation_of(frame.provider(), served, input),
            Denotation::NoSequence,
            "{geometry}: `{input}` must denote no sequence for this test to mean anything"
        );

        let output = normalize_in(&frame, ErrorMode::Lenient, input)
            .unwrap_or_else(|e| panic!("{geometry}: lenient must not reject `{input}`: {e}"));
        assert_eq!(
            output, input,
            "{geometry}: the conflicting allele must come back as authored"
        );
        assert_eq!(
            denotation_of(frame.provider(), served, &output),
            Denotation::NoSequence,
            "{geometry}: and the output denotes no sequence — a rank-1 validity failure \
             that lenient mode returns as `Ok`"
        );
        // Stable in one pass, so preserving is an answer rather than a pause.
        assert_eq!(
            normalize_in(&frame, ErrorMode::Lenient, &output).expect("lenient"),
            output,
            "{geometry}: the preserved form must be a fixed point"
        );
    }
}

/// **Question.** Is the conflict reported to a lenient caller?
///
/// **Only to one that asks a different question.** The W5002 `OverlapConflict`
/// warning is produced — `normalize_with_diagnostics` carries exactly one per
/// input — but `Normalizer::normalize`, the entry point the CLI and every
/// in-repo caller use, returns `Result<HgvsVariant, _>` and has nowhere to put
/// it (`src/normalize/mod.rs:1764`).
///
/// **PINNED DEFECT.** The consequence, asserted rather than described: through
/// `normalize`, `ErrorMode::Lenient` and `ErrorMode::Silent` are
/// indistinguishable — same `Ok`, same string — while `ErrorMode` itself
/// declares they differ (`emits_warnings()`, `src/error_handling/types.rs:46-48`)
/// and `Lenient`'s doc promises "warnings will be generated" (`:20-24`).
/// Confirmed end-to-end at the CLI: `ferro normalize --error-mode lenient` and
/// `--error-mode silent` both print the description on stdout with empty stderr
/// and exit 0.
///
/// **Correct behaviour:** a lenient normalization that declined to canonicalize
/// must say so through the same call that returned the description — the input
/// still passes through unchanged, refusal stays strict mode's job. Until then a
/// caller cannot distinguish "ferro normalized this" from "ferro refused to, and
/// this denotes no sequence".
#[test]
fn the_conflict_is_dropped_by_the_entry_point_that_returns_the_output() {
    let core = at_core();
    let frame = Frame::build(RefShape::CodingSingleExon, &core);

    for (geometry, input) in CONFLICTING_ALLELES {
        // The finding exists. Asserted first, so the pin below is about the
        // channel rather than about detection having quietly stopped.
        assert_eq!(
            conflict_diagnostics(&frame, ErrorMode::Lenient, input),
            1,
            "{geometry}: the conflict must be detected for this test to be about reporting"
        );
        // And strict mode reaches it, which is what makes the drop a seam
        // problem rather than a detection one.
        let strict = normalize_in(&frame, ErrorMode::Strict, input)
            .expect_err("strict mode must reject a conflicting allele");
        assert!(
            strict.contains("W5002") || strict.contains("OverlapConflictingEdits"),
            "{geometry}: strict must reject `{input}` AS a conflict; got: {strict}"
        );

        // The defect. `ErrorMode` says these two modes differ...
        assert!(
            ErrorMode::Lenient.emits_warnings() && !ErrorMode::Silent.emits_warnings(),
            "the modes must claim to differ for the pin below to be a contradiction"
        );
        // ...and through `normalize` they do not.
        assert_eq!(
            normalize_in(&frame, ErrorMode::Lenient, input),
            normalize_in(&frame, ErrorMode::Silent, input),
            "PINNED DEFECT ({geometry}) — lenient and silent are indistinguishable through \
             `Normalizer::normalize`, which has no warning channel. Correct behaviour: a \
             lenient caller must learn from this call that `{input}` was declined and denotes \
             no sequence."
        );
    }
}

/// **Question.** Does the pass-through at least keep the conflict *visible in the
/// description*, so that re-reading the output reaches the same verdict?
///
/// **Yes, and that is the property to protect.** #1406 row 3 was the case where
/// it did not: lenient mode repaired the members one at a time until the output
/// no longer looked conflicting, so strict mode accepted a description it had
/// just rejected. The gate at `src/normalize/mod.rs:3450-3481` closed that, and
/// `issue_1406_lenient_output_keeps_the_conflict.rs` pins it for the
/// `[inv;ins]` family.
///
/// **ADJUDICATED CORRECT**, asserted here for the three *corpus* geometries,
/// which that file does not cover. It is the reason refusal is not needed in
/// lenient mode: the invalidity is still recoverable from the output by anyone
/// who re-reads it strictly. What it does not fix is that the caller who got the
/// output was never told to.
#[test]
fn the_lenient_output_is_still_refused_by_strict() {
    let core = at_core();
    let frame = Frame::build(RefShape::CodingSingleExon, &core);

    for (geometry, input) in CONFLICTING_ALLELES {
        let output = normalize_in(&frame, ErrorMode::Lenient, input).expect("lenient");
        let strict = normalize_in(&frame, ErrorMode::Strict, &output)
            .expect_err("the lenient output must still be a conflict");
        assert!(
            strict.contains("W5002") || strict.contains("OverlapConflictingEdits"),
            "{geometry}: lenient turned `{input}` into `{output}`, which strict does not \
             reject as a conflict — the conflict was laundered out of the description"
        );
    }
}

// ---------------------------------------------------------------------------
// Item 2: a single-member description projecting into a bracketed allele
// ---------------------------------------------------------------------------
//
// The repository already has a name and a counter for the shape — the spec
// enumeration's `Status::ProjectionSplitsSingleMember`, budgeted at 9
// (`tests/it/spec_enumeration_tests.rs:356`) — but no test names the rows and
// no issue has ever been filed on it. It has only ever been mentioned inside
// other work: created by #1272, moved 10 -> 9 by #1271, measured 9 -> 13 -> 9
// across #1235's branch.
//
// **The shape alone is not the defect, and saying otherwise is the trap.** Of
// the nine rows, seven are inputs whose members end up three or more bases
// apart, where `general.md:34` — "two variants separated by one or more
// nucleotides should be described individually and **not** as a 'delins'" —
// plainly asks for the split. The enumeration's own comment says as much: the
// rows "enter by *shape* … not by the divergence `delins.md:42` names".
//
// The two rows that ARE a violation are the ones pinned below, and what makes
// them one is not the split but *which axis* splits.

use ferro_hgvs::conformance::spec_projection::load_slice;
use ferro_hgvs::data::{CdotMapper, CdotTranscript, Projector};
use ferro_hgvs::project::VariantProjector;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::VariantProjection;

/// The committed hermetic reference slice the spec enumeration projects on —
/// transcripts, genomic windows and `NG_`/`LRG_` placements, regenerated by
/// `cargo run --features dev --example extract_spec_enumeration_windows`.
///
/// Committed, unlike the enumeration artifact it feeds, so these tests need no
/// prepared reference and no `FERRO_MANIFEST`.
const SLICE: &str = "tests/fixtures/grammar/spec_enumeration_windows.json";

/// A projector over the committed slice, built exactly as
/// `generate_spec_enumeration` builds its own.
fn slice_projector() -> VariantProjector<ferro_hgvs::conformance::reference_window::WindowProvider>
{
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(SLICE);
    let windows = load_slice(&path).unwrap_or_else(|e| panic!("{} must load: {e}", path.display()));
    let cdot = CdotMapper::from_transcripts(windows.transcripts.iter());
    VariantProjector::new(Projector::new(cdot), windows.to_provider())
}

/// Every rendered axis of `projection`, as `(axis letter, description)`.
fn rendered_axes(projection: &VariantProjection) -> Vec<(char, String)> {
    [
        ('g', &projection.genomic),
        ('c', &projection.coding),
        ('n', &projection.noncoding),
        ('r', &projection.rna),
        ('p', &projection.protein),
    ]
    .into_iter()
    .filter_map(|(code, axis)| axis.as_ref().map(|v| (code, v.to_string())))
    .collect()
}

/// Does `rendered` carry a top-level allele bracket?
///
/// The enumeration's own test (`renders_multiple_members`,
/// `examples/generate_spec_enumeration.rs:1858`) asks the same question the same
/// way — a `[` immediately after the axis letter's `.` — because a bracket can
/// also appear inside an inserted-sequence payload (`ins[A;100_110]`), which is
/// one member, not two.
fn splits_into_members(rendered: &str) -> bool {
    let body = rendered.rsplit_once(':').map_or(rendered, |(_, tail)| tail);
    body.split_once('.')
        .is_some_and(|(_, edit)| edit.starts_with('[') || edit.starts_with("(["))
}

/// **Question.** May a projection of a description that names ONE variant render
/// that variant as an allele of two?
///
/// **Not this one, and the spec says so about this exact string.**
/// `LRG_199t1:c.145_147delinsTGG` is the spec's own worked example
/// (`DNA/delins.md:37`), and the NOTE two lines below it, `DNA/delins.md:42`,
/// reads: "two variants separated by one nucleotide, together affecting one
/// amino acid, should be described as a "delins", so the description
/// `c.[145C>T;147C>G]` is not correct" — with that split spelling rendered
/// `class="invalid"` in the published page.
///
/// **PINNED DEFECT.** Ferro applies that exception on four axes and not on the
/// fifth. One input, one projection:
///
/// ```text
/// c   LRG_199t1:c.145_147delinsTGG      1 member   correct
/// n   LRG_199t1:n.389_391delinsTGG      1 member
/// r   LRG_199t1:r.(145_147delinsugg)    1 member
/// p   LRG_199p1:p.(Arg49Trp)            1 member
/// g   LRG_199:g.[494841C>T;494843C>G]   2 members  <- the invalid form
/// ```
///
/// `g.494841` and `g.494843` are `c.145` and `c.147`, so the genomic axis emits
/// precisely the description `:42` calls not correct, re-coordinated. **Correct
/// output: `LRG_199:g.494841_494843delinsTGG`** — one member, because the
/// exception at `DNA/delins.md:18` ("**exception**: two variants separated by one
/// nucleotide, together affecting one amino acid, should be described as a
/// 'delins'") is a property of the variant, and this variant affects one amino
/// acid whichever coordinate system it is spelled in. That reading is already
/// settled for ferro: `rulings[delins-codon-carve-out-gap-one]` is **decided**,
/// `delins.md:18` governing.
///
/// The mechanism is that the exception is evaluated against the axis being
/// rendered rather than against the variant. A `g.` axis has no reading frame of
/// its own, so `general.md:34`'s floor of one applies unopposed and the pair
/// separates. Nothing about the cis-allele partitioner is implicated — the input
/// has no members to partition.
#[test]
fn the_genomic_axis_alone_splits_the_specs_own_codon_delins() {
    let projector = slice_projector();

    for (input, transcript, genomic) in [
        (
            "LRG_199t1:c.145_147delinsTGG",
            "LRG_199t1",
            "LRG_199:g.[494841C>T;494843C>G]",
        ),
        // `delins.md:19`'s example, whose stated purpose is to stop tools
        // predicting `p.[Lys79*;Lys79Asn]` for one Lys79Tyr change. Ferro's `p.`
        // axis gets that right — `p.(Lys79Tyr)` — while its `g.` axis publishes
        // the split the note exists to prevent.
        (
            "LRG_199t1:c.235_237delinsTAT",
            "LRG_199t1",
            "LRG_199:g.[499798A>T;499800G>T]",
        ),
    ] {
        let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));
        let projection = projector
            .project_variant(&variant, transcript)
            .unwrap_or_else(|e| panic!("{input} must project: {e}"));
        let axes = rendered_axes(&projection);

        // The premise: this input names exactly one variant. Asserted, because a
        // row that quietly acquired a bracket in its *input* would satisfy
        // everything below while pinning nothing.
        assert!(
            !splits_into_members(input),
            "{input} must be a single-member description for this test to mean anything"
        );
        assert!(
            axes.len() >= 4,
            "{input} rendered only {axes:?}; the slice must serve every axis for the \
             four-against-one comparison below to be a comparison"
        );

        let split: Vec<&(char, String)> = axes
            .iter()
            .filter(|(_, rendered)| splits_into_members(rendered))
            .collect();
        assert_eq!(
            split.len(),
            1,
            "PINNED DEFECT — exactly one axis of `{input}` must currently split, and it \
             must be the genomic one. Rendered axes: {axes:?}"
        );
        assert_eq!(
            (split[0].0, split[0].1.as_str()),
            ('g', genomic),
            "PINNED DEFECT — the genomic axis emits the form `DNA/delins.md:42` calls \
             \"not correct\". Correct output: one member, the delins re-coordinated. \
             When this is fixed, assert the single-member genomic form here and drop \
             the split assertion above."
        );
    }
}

/// **Question.** Is the split a property of the variant or of how the caller
/// spelled it?
///
/// **Of the spelling**, which is the sharper statement of the defect and needs
/// no reference at all. On a synthetic single-exon transcript placed so that
/// `g.N == c.N`, one variant — replacing codon 2's `CTG` with `ATA` — comes back
/// as one member when authored on `c.` and as two when authored on `g.`, and in
/// the `g.` case the split is carried onto the coding axis as well:
///
/// ```text
/// tx:c.4_6delinsATA    ->  c. c.4_6delinsATA        n. r.  all ONE member
/// REF:g.4_6delinsATA   ->  g. g.[4C>A;6G>A]
///                          c. c.[4C>A;6G>A]         n. r.  all TWO members
/// ```
///
/// Both spellings denote the same bases and the same `p.(Leu2Ile)`.
///
/// **Be exact about which line is the defect here, because one of them is not.**
/// `REF:g.[4C>A;6G>A]` as the answer to a *bare* `g.` input is #79's deliberate
/// scope — that issue's own table lists "`g.[145G>A;147C>T]` (no codon frame in
/// `g.`) | unchanged" — and a bare genomic description names no transcript, so
/// there is no frame to consult. **PINNED, not a defect.**
///
/// The defect is the line below it. `REF(tx):c.[4C>A;6G>A]` is a **coding** axis,
/// rendered against a named transcript whose frame is known and whose codon 2 is
/// exactly what the pair affects — and it is the shape `DNA/delins.md:42` calls
/// "not correct". The split was decided once, on the axis the caller happened to
/// author on, and every derived axis inherited it.
///
/// So the two tests state one rule between them: **the exception is evaluated
/// once per rendered axis, against the transcript the projection is against —
/// never inherited from the authoring axis.** The reference-backed test above is
/// that rule failing in the `c.` -> `g.` direction; this one is the same rule
/// failing in the `g.` -> `c.` direction, and it runs everywhere.
///
/// **Correct output:** `REF(tx):c.4_6delinsATA` and `REF(tx):n.4_6delinsATA`,
/// matching the transcript spelling. The `g.` line stays as it is.
#[test]
fn the_same_codon_delins_splits_when_authored_genomically_and_not_when_authored_on_the_transcript()
{
    /// Met-Leu-Arg-Trp-Ala-Stop. Codon 2 is `CTG`; replacing it with `ATA`
    /// changes bases 1 and 3 and leaves base 2, so the pair is separated by
    /// exactly one nucleotide and affects exactly one amino acid — the shape
    /// `DNA/delins.md:18` carves out.
    const CDS: &str = "ATGCTGCGTTGGGCCTAA";

    let projector = synthetic_projector(CDS);
    let coding = axis_map(&projector, "tx:c.4_6delinsATA");
    let genomic = axis_map(&projector, "REF:g.4_6delinsATA");

    // The transcript spelling: the carve-out applies, so nothing splits.
    assert_eq!(
        coding,
        vec![
            ('c', "tx:c.4_6delinsATA".to_string()),
            ('n', "tx:n.4_6delinsATA".to_string()),
            ('r', "tx:r.(4_6delinsaua)".to_string()),
            ('p', "txp:p.(Leu2Ile)".to_string()),
        ],
        "the transcript spelling is the control and must stay merged"
    );

    // The genomic spelling of the same variant. The `g.` line is #79's scope and
    // is fine; the `c.`/`n.` lines below it are the defect, because those axes
    // are rendered against a transcript whose frame is known.
    assert_eq!(
        genomic,
        vec![
            ('g', "REF:g.[4C>A;6G>A]".to_string()),
            ('c', "REF(tx):c.[4C>A;6G>A]".to_string()),
            ('n', "REF(tx):n.[4C>A;6G>A]".to_string()),
            ('r', "tx:r.[(4c>a;6g>a)]".to_string()),
            ('p', "txp:p.(Leu2Ile)".to_string()),
        ],
        "PINNED DEFECT (the transcript axes only) — `REF(tx):c.[4C>A;6G>A]` is the \
         form `DNA/delins.md:42` calls \"not correct\", emitted on an axis whose \
         reading frame ferro has in hand. Correct output: `REF(tx):c.4_6delinsATA` \
         and `REF(tx):n.4_6delinsATA`, matching the transcript spelling above; the \
         bare `g.` line stays as it is."
    );

    // Stated as the property too, so a future fix that moves the strings but
    // keeps the disagreement cannot satisfy this test by re-blessing.
    //
    // Counted over the TRANSCRIPT axes only, which is what makes the two columns
    // comparable: the `g.` axis has no frame to consult when the caller wrote a
    // bare `g.` description, so including it would compare an axis that has the
    // frame against one that does not and the difference would prove nothing.
    let transcript_splits = |axes: &[(char, String)]| {
        axes.iter()
            .filter(|(code, rendered)| *code != 'g' && splits_into_members(rendered))
            .count()
    };
    assert!(
        transcript_splits(&coding) == 0 && transcript_splits(&genomic) > 0,
        "PINNED DEFECT — one variant, two spellings, and its transcript axes disagree \
         about member count ({} split from the `c.` spelling, {} from the `g.` one). \
         Correct behaviour: both zero, because both are rendered against the same \
         transcript and the same codon.",
        transcript_splits(&coding),
        transcript_splits(&genomic)
    );
}

/// Every rendered axis of `input`, projected on `tx`.
fn axis_map(projector: &VariantProjector<MockProvider>, input: &str) -> Vec<(char, String)> {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));
    let projection = projector
        .project_variant(&variant, "tx")
        .unwrap_or_else(|e| panic!("{input} must project: {e}"));
    rendered_axes(&projection)
}

/// A single-exon plus-strand coding transcript `tx` on contig `REF`, placed so
/// that `g.N == c.N` — which is what makes the two spellings above an A/B on the
/// axis alone rather than on the coordinates.
fn synthetic_projector(cds: &str) -> VariantProjector<MockProvider> {
    let len = cds.len() as u64;
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "tx".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("GENE".to_string()),
            contig: "REF".to_string(),
            strand: Strand::Plus,
            exons: vec![[1, len + 1, 0, len]],
            cds_start: Some(0),
            cds_end: Some(len),
            gene_id: None,
            protein: Some("txp".to_string()),
            exon_cigars: Vec::new(),
        },
    );

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "tx".to_string(),
        Some("GENE".to_string()),
        TxStrand::Plus,
        cds.to_string(),
        Some(1),
        Some(len),
        vec![Exon::new(1, 1, len)],
        Some("REF".to_string()),
        Some(1),
        Some(len),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    provider.add_genomic_sequence("REF", format!("{cds}NN"));

    VariantProjector::new(Projector::new(cdot), provider)
}
